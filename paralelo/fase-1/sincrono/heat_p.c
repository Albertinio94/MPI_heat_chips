/* heat_p.c

	 Difusion del calor en 2 dimensiones      Version en paralelo

	 Se analizan las posiciones de los chips en una tarjeta, para conseguir la temperatura minima
	 de la tarjeta. Se utiliza el metodo de tipo Poisson, y la tarjeta se discretiza en una rejilla
	 de puntos 2D.

	 Entrada: card > la definicion de la tarjeta y las configuraciones a simular
	 Salida: la mejor configuracion y la temperatura media
		  card_s.chips: situacion termica inicial
		card_s.res: la situacion termica final

	 defines.h: definiciones de ciertas variables y estructuras de datos

	 Compilar con estos dos ficheros:
	   diffusion.c: insertar calor, difundir y calcular la temperatura media hasta que se estabilice
	   faux.c: ciertas funciones auxiliares

************************************************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <time.h>

#include "defines.h"
#include "faux_p.h"
#include "diffusion_p.h"

/************************************************************************************/
void init_grid_chips(int conf, struct info_param param, struct info_chips *chips, int **chip_coord, float *grid_chips) {
	int i, j, n;

	for (i = 0; i < NROW; i++)
		for (j = 0; j < NCOL; j++)
			grid_chips[i * NCOL + j] = param.t_ext;

	for (n = 0; n < param.nchip; n++)
		for (i = chip_coord[conf][2 * n] * param.scale; i < (chip_coord[conf][2 * n] + chips[n].h) * param.scale; i++)
			for (j = chip_coord[conf][2 * n + 1] * param.scale; j < (chip_coord[conf][2 * n + 1] + chips[n].w) * param.scale; j++)
				grid_chips[(i + 1) * NCOL + (j + 1)] = chips[n].tchip;
}

/************************************************************************************/
void init_grids(struct info_param param, float *grid, float *grid_aux, int tam) {
	int i, j;

	for (i = 0; i < tam; i++)
		for (j = 0; j < NCOL; j++)
			grid[i * NCOL + j] = param.t_ext;
}

/************************************************************************************/
/************************************************************************************/
int main(int argc, char *argv[])
{
	struct info_param param;
	struct info_chips *chips;
	int **chip_coord;

	float *grid, *final_grid, *grid_chips_root, *grid_chips_local, *grid_aux;
	struct info_results BT;

	int conf, i;
	double t0, t1;
	double tej, Tmean;

	int pid, npr, size;
	int quotient, remainder, *tam, *dist;

	struct MPI_info_param mpi_param;

	MPI_Datatype row;
	
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &npr);

	// Check if card description file has been given
	if (argc != 2) {
		if (pid == ROOT) printf("\n\nERROR: needs a card description file\n\n");
		exit(-1);
	}

	if (pid == ROOT) {
		// reading initial data file
		read_data(argv[1], &param, &chips, &chip_coord);

		printf("\n  ===================================================================");
		printf("\n    Thermal diffusion - PARALLEL version ");
		printf("\n    %d x %d points, %d chips", RSIZE * param.scale, CSIZE * param.scale, param.nchip);
		printf("\n    T_ext = %1.1f, Tmax_chip = %1.1f, T_delta: %1.3f, Max_iter: %d", param.t_ext, param.tmax_chip, param.t_delta, param.max_iter);
		printf("\n  ===================================================================\n\n");

		t0 = MPI_Wtime();
	}

	// Create struct in order to send simulation parameters to every process
	create_mpi_params(&param, &mpi_param);

	// Send simulation parameters to every process (param.scale needed in order to calculate NROW)
	MPI_Bcast(&param, 1, mpi_param.param, ROOT, MPI_COMM_WORLD);

	// Create row type
	MPI_Type_vector(NCOL, 1, 1, MPI_FLOAT, &row);
	MPI_Type_commit(&row);

	// Calculate piece of data to send to each process
	size = npr * sizeof(int);	
	// vector with sizes of each piece
	tam = (int *) malloc(size);
	// distances from the beggining of each piece
	dist = (int *) malloc(size);

	quotient = (NROW) / npr;
	remainder = (NROW) % npr;

	// Calculate the distribution (remainder distributed one by one among first pieces)
	for (i = 0; i < npr; i++) {
		tam[i] = quotient;
		if (i < remainder) tam[i]++;
		if (i == 0) dist[i] = 0;
		else dist[i] = dist[i-1] + tam[i-1];
	}

	// each one takes its pieces + space for neighbor rows
	grid = (float *) malloc((tam[pid]+2)*NCOL * sizeof(float));
	grid_aux = (float *) malloc(tam[pid]*NCOL * sizeof(float));
	grid_chips_local = (float *) malloc(tam[pid]*NCOL * sizeof(float));

	if (pid == ROOT) {
		final_grid = (float *) malloc(NROW*NCOL * sizeof(float));
		grid_chips_root = (float *) malloc(NROW*NCOL * sizeof(float));
		BT.bgrid = (float *) malloc(NROW*NCOL * sizeof(float));
		BT.cgrid = (float *) malloc(NROW*NCOL * sizeof(float));
		BT.Tmean = MAXDOUBLE;
	}

	// loop to process chip configurations
	for (conf = 0; conf < param.nconf; conf++) {
		// initial values for grids
		if (pid == ROOT) init_grid_chips(conf, param, chips, chip_coord, grid_chips_root);
		init_grids(param, grid, grid_aux, tam[pid]+2);
		MPI_Scatterv(grid_chips_root, tam, dist, row, grid_chips_local, tam[pid], row, ROOT, MPI_COMM_WORLD);

		// main loop: thermal injection/disipation until convergence (t_delta or max_iter)
		Tmean = calculate_Tmean(param, grid, grid_chips_local, grid_aux, pid, npr, tam[pid], row);

		// gather calculated grids
		MPI_Gatherv(&grid[NCOL], tam[pid], row, final_grid, tam, dist, row, ROOT, MPI_COMM_WORLD);

		if (pid == ROOT) {
			printf("  Config: %2d    Tmean: %1.2f\n", conf + 1, Tmean);
			// processing configuration results
			results_conf(conf, Tmean, param, final_grid, grid_chips_root, &BT);
		}
	}
	
	if (pid == ROOT) {
		t1 = MPI_Wtime();
		tej = t1 - t0;
		printf("\n\n >>> Best configuration: %2d    Tmean: %1.2f\n", BT.conf + 1, BT.Tmean);
		printf("   > Time (parallel): %1.3f s \n\n", tej);
		// writing best configuration results
		results(param, &BT, argv[1]);
	}

	free(tam);
	free(dist);
	free(grid);
	free(grid_aux);
	// free(grid_chips_local); // descomentarlo provoca segmentation fault
	if (pid == ROOT) {
		for (i = 0; i < param.nconf; i++) free(chip_coord[i]);
		free(chip_coord);
		free(chips);
		free(BT.bgrid);
		free(BT.cgrid);
		// free(final_grid); // // descomentarlo provoca segmentation fault
		// free(grid_chips_root); // descomentarlo provoca segmentation fault
	}

	MPI_Type_free(&row);
	MPI_Finalize();

	return (0);
}