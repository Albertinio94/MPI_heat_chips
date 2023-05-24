/* heat_p.c

	 Difusion del calor en 2 dimensiones      Version en paralelo

	 Se analizan las posiciones de los chips en una tarjeta, para conseguir la temperatura minima
	 de la tarjeta. Se utiliza el metodo de tipo Poisson, y la tarjeta se discretiza en una rejilla
	 de puntos 2D.

	 Entrada: card > la definicion de la tarjeta y las configuraciones a simular
	 Salida: la mejor configuracion y la temperatura media
		card_par.chips: situacion termica inicial
		card_par.res: la situacion termica final

	 defines.h: definiciones de ciertas variables y estructuras de datos

	 Compilar con estos dos ficheros:
	   diffusion_p.c: insertar calor, difundir y calcular la temperatura media hasta que se estabilice
	   faux_p.c: ciertas funciones auxiliares

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
//Function-like macro to avoid having to pass a lot of parameters
#define process_request() {\
	MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);\
	if (status.MPI_TAG == CONF_RESULT) {\
		MPI_Recv(buf_resultado_conf, TAM_PACK, MPI_PACKED, status.MPI_SOURCE, CONF_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);\
		/* extract Tmean and final grid*/\
		pos = 0;\
		MPI_Unpack(buf_resultado_conf, TAM_PACK, &pos, &Tmean, 1, MPI_DOUBLE, MPI_COMM_WORLD);\
		MPI_Unpack(buf_resultado_conf, TAM_PACK, &pos, final_grid, NROW*NCOL, MPI_FLOAT, MPI_COMM_WORLD);\
		/*Save configuration*/\
		printf("  Config: %2d    Tmean: %1.2f\n", last_conf_sent_to_groups[status.MPI_SOURCE/P], Tmean);\
		results_conf(last_conf_sent_to_groups[status.MPI_SOURCE/P], Tmean, param, final_grid, &BT);\
	} else { /* prevent buffer block*/\
		MPI_Recv(NULL, 0, MPI_INT, status.MPI_SOURCE, CONF_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);\
	}\
}

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

	float *grid, *final_grid, *grid_chips_global, *grid_chips_local, *grid_aux;
	float *current_matrix, *not_current_matrix;
	struct info_results BT;

	int conf = 0, i, j;
	double t0, t1;
	double tej, Tmean;

	int pid, pid_worker, group_leader, npr, size;
	int quotient, remainder, *tam, *dist;

	MPI_Comm comm_worker;
	int worker_group_key, *last_conf_sent_to_groups, pos, tag;
	char *buf_resultado_conf; // Tmean + new grid
	MPI_Status status;

	struct MPI_info_param mpi_param;

	MPI_Datatype row;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &npr);

	// Check if card description file has been given
	if (argc != 2) {
		if (IS_MANAGER) printf("\n\nERROR: needs a card description file.\n\n");
		exit(-1);
	} else if (npr == 1 || npr != 1 + P * ((npr - 1) / P)) {
		if (IS_MANAGER) printf("\n\nERROR: number of processors must be 1 + k*%d, with k > 0.\n\n", P);
		exit(-1);
	}

	if (IS_MANAGER) {
		// reading initial data file
		read_data(argv[1], &param, &chips, &chip_coord);

		printf("\n  ===================================================================");
		printf("\n    Thermal diffusion - PARALLEL version (manager/worker)");
		printf("\n    %d x %d points, %d chips", RSIZE * param.scale, CSIZE * param.scale, param.nchip);
		printf("\n    T_ext = %1.1f, Tmax_chip = %1.1f, T_delta: %1.3f, Max_iter: %d", param.t_ext, param.tmax_chip, param.t_delta, param.max_iter);
		printf("\n  ===================================================================\n\n");

		t0 = MPI_Wtime();
	}

	// Create struct in order to send simulation parameters to every process
	create_mpi_params(&param, &mpi_param);

	// Send simulation parameters to every process (param.scale needed in order to calculate NROW)
	MPI_Bcast(&param, 1, mpi_param.param, MANAGER, MPI_COMM_WORLD);

	// Create row type
	MPI_Type_vector(NCOL, 1, 1, MPI_FLOAT, &row);
	MPI_Type_commit(&row);

	// Create communicators for groups of workers
	worker_group_key = pid/P; // manager = npr - 1
	MPI_Comm_split(MPI_COMM_WORLD, worker_group_key, pid, &comm_worker);
	MPI_Comm_rank (comm_worker, &pid_worker);
	group_leader = worker_group_key * P;

	// // Calculate piece of data to send to each process
	size = P * sizeof(int);	
	// vector with sizes of each piece
	tam = (int *) malloc(size);
	// distances from the beggining of each piece
	dist = (int *) malloc(size);

	quotient = NROW / P;
	remainder = NROW % P;

	// Calculate the distribution (remainder distributed one by one among first pieces)
	for (i = 0; i < P; i++) {
		tam[i] = quotient;
		if (i < remainder) tam[i]++;
		if (i == 0) dist[i] = 0;
		else dist[i] = dist[i-1] + tam[i-1];
	}
	
	if (IS_MANAGER) {
		grid_chips_global = (float *) malloc(NROW*NCOL * sizeof(float));
		final_grid = (float *) malloc(NROW*NCOL * sizeof(float));	
		BT.bgrid = (float *) malloc(NROW*NCOL * sizeof(float));
		BT.cgrid = (float *) malloc(NROW*NCOL * sizeof(float));
		BT.Tmean = MAXDOUBLE;
		buf_resultado_conf = (char *) malloc(TAM_PACK);
		last_conf_sent_to_groups = (int *) malloc(((npr-1)/P) * sizeof(int));
	} else if (IS_LEADER) {
		grid_chips_global = (float *) malloc(NROW*NCOL * sizeof(float));
		grid_chips_local = (float *) malloc(tam[pid_worker]*NCOL * sizeof(float));
		final_grid = (float *) malloc(NROW*NCOL * sizeof(float));	
		grid = (float *) malloc((tam[pid_worker]+2)*NCOL * sizeof(float));
		grid_aux = (float *) malloc(tam[pid_worker]*NCOL * sizeof(float));
		buf_resultado_conf = (char *) malloc(TAM_PACK);
	} else {
		grid = (float *) malloc((tam[pid_worker]+2)*NCOL * sizeof(float));
		grid_aux = (float *) malloc(tam[pid_worker]*NCOL * sizeof(float));
		grid_chips_local = (float *) malloc(tam[pid_worker]*NCOL * sizeof(float));
	}

	// loop to process chip configurations
	if (IS_MANAGER) {
		for (conf = 0; conf < param.nconf; conf++) {
			init_grid_chips(conf, param, chips, chip_coord, grid_chips_global);

			process_request();
			// Load new grid chips and send it to leader
			last_conf_sent_to_groups[status.MPI_SOURCE/P] = conf+1;
			MPI_Send(grid_chips_global, NROW*NCOL, MPI_FLOAT, status.MPI_SOURCE, CONF, MPI_COMM_WORLD);
		}

		// attend last request of each group
		for (i = 0; i < (npr-1)/P; i++) {
			process_request();

			MPI_Send(NULL, 0, MPI_INT, status.MPI_SOURCE, FINALIZE, MPI_COMM_WORLD);
		}

		// Get best initial grid
		init_grid_chips(BT.conf, param, chips, chip_coord, grid_chips_global);
		for (i = 1; i < NROW - 1; i++)
			for (j = 1; j < NCOL - 1; j++)
				BT.cgrid[i * NCOL + j] = grid_chips_global[i * NCOL + j];
		// Compute time
		t1 = MPI_Wtime();
		tej = t1 - t0;
		printf("\n\n >>> Best configuration: %2d    Tmean: %1.2f\n", BT.conf, BT.Tmean);
		printf("   > Time (parallel): %1.3f s \n\n", tej);
		// Write best configuration results
		results(param, &BT, argv[1]);

	} else {
		if (IS_LEADER) {
			MPI_Send(NULL, 0, MPI_INT, MANAGER, CONF_REQUEST, MPI_COMM_WORLD);
		}

		while (1) { // recibir configuraciones
			if (IS_LEADER) {
				MPI_Probe(MANAGER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				tag = status.MPI_TAG;
				if (tag == CONF) MPI_Recv(grid_chips_global, NROW*NCOL, MPI_FLOAT, MANAGER, CONF, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				else MPI_Recv(NULL, 0, MPI_INT, MANAGER, FINALIZE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			// enviar tags al resto de workers (0 ==> seguimos, 1 ==> mpi_finalize)
			if (!IS_LEADER) tag = 5;
			MPI_Bcast(&tag, 1, MPI_INT, LEADER, comm_worker);

			if (tag == CONF) {
				init_grids(param, grid, grid_aux, tam[pid_worker]+2);
				MPI_Barrier(comm_worker);
				MPI_Scatterv(grid_chips_global, tam, dist, row, grid_chips_local, tam[pid_worker], row, LEADER, comm_worker);

				// main loop: thermal injection/disipation until convergence (t_delta or max_iter)
				Tmean = calculate_Tmean(param, grid, grid_chips_local, grid_aux, pid_worker, P, tam[pid_worker], row, comm_worker);

				// gather calculated grids
				MPI_Gatherv(&grid[NCOL], tam[pid_worker], row, final_grid, tam, dist, row, LEADER, comm_worker);

				if (IS_LEADER) {
					pos = 0;
					MPI_Pack(&Tmean, 1, MPI_DOUBLE, buf_resultado_conf, TAM_PACK, &pos, MPI_COMM_WORLD);
					MPI_Pack(final_grid, NROW*NCOL, MPI_FLOAT, buf_resultado_conf, TAM_PACK, &pos, MPI_COMM_WORLD);
					MPI_Send(buf_resultado_conf, TAM_PACK, MPI_PACKED, MANAGER, CONF_RESULT, MPI_COMM_WORLD);
				}
			} else { // tag = FINALIZE
				break;
			}
		}
		
	}

	// free(tam);
	// free(dist);
	// free(grid);
	// free(grid_aux);
	// // free(grid_chips_local); // uncommenting it causes segmentation fault
	// if (pid == MANAGER) {
	// 	for (i = 0; i < param.nconf; i++) free(chip_coord[i]);
	// 	free(chip_coord);
	// 	free(chips);
	// 	free(BT.bgrid);
	// 	free(BT.cgrid);
	// 	// free(final_grid); // // uncommenting it causes segmentation fault
	// 	// free(grid_chips_root); // uncommenting it causes segmentation fault
	// }


	// MPI_Type_free(&row);
	// MPI_Comm_free(&comm_worker);
	MPI_Finalize();
	return (0);
}