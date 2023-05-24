/* File: diffusion_s.c */ 

#include <mpi.h>
#include <stdio.h>
#include <math.h>

#include "defines.h"



/************************************************************************************/
void thermal_update(struct info_param param, float *grid, float *grid_chips, int pid, int npr, int tam) {
	int i, j, a, b;
	// prevent thermal diffusion on borders
	int first_row = (pid == 0) ? 2 : 1;
	int last_row = (pid == npr-1) ? tam : tam+1;

	// heat injection at chip positions
	for (i = first_row; i < last_row; i++)
		for (j = 1; j < NCOL - 1; j++)
			if (grid_chips[(i-1) * NCOL + j] > grid[i * NCOL + j])
				grid[i * NCOL + j] += 0.05 * (grid_chips[(i-1) * NCOL + j] - grid[i * NCOL + j]);

	// air cooling at the middle of the card
	a = 0.44 * (NCOL - 2) + 1;
	b = 0.56 * (NCOL - 2) + 1;

	for (i = first_row; i < last_row; i++)
		for (j = a; j < b; j++)
			grid[i * NCOL + j] -= 0.01 * (grid[i * NCOL + j] - param.t_ext);
}

/************************************************************************************/
double thermal_diffusion(struct info_param param, float *grid, float *grid_aux, int pid, int npr, int tam) {
	int i, j;
	float T;
	double Tfull = 0.0;
	// prevent thermal diffusion on borders
	int first_row = (pid == 0) ? 2 : 1;
	int last_row = (pid == npr-1) ? tam : tam+1;

	for (i = first_row; i < last_row; i++)
		for (j = 1; j < NCOL - 1; j++) {
			T = grid[i * NCOL + j] +
				0.10 * (grid[(i + 1) * NCOL + j] + grid[(i - 1) * NCOL + j] + grid[i * NCOL + (j + 1)] + grid[i * NCOL + (j - 1)] +
						grid[(i + 1) * NCOL + j + 1] + grid[(i - 1) * NCOL + j + 1] + grid[(i + 1) * NCOL + (j - 1)] + grid[(i - 1) * NCOL + (j - 1)] - 8 * grid[i * NCOL + j]);

			grid_aux[(i-1) * NCOL + j] = T;
			Tfull += T;
		}

	// new values for the grid
	for (i = first_row; i < last_row; i++)
		for (j = 1; j < NCOL - 1; j++)
			grid[i * NCOL + j] = grid_aux[(i-1) * NCOL + j];

	return (Tfull);
}

/************************************************************************************/
double calculate_Tmean(struct info_param param, float *grid, float *grid_chips, float *grid_aux, int pid, int npr, int tam, MPI_Datatype row, MPI_Comm comm_worker) {
	int i, j, end, niter;
	float Tfull;
	double t_local, t_global, t_mean, t_mean_0 = param.t_ext;
	MPI_Request reqs[2];

	end = 0;
	niter = 0;

	while (end == 0) {
		niter++;
		t_mean = 0.0;

		// heat injection and air cooling
		thermal_update(param, grid, grid_chips, pid, npr, tam);

		// Send border rows to neighbors
		if (pid < npr-1) MPI_Isend(&grid[tam*NCOL], 1, row, pid+1, 0, comm_worker, &reqs[0]);
		if (pid > 0) MPI_Isend(&grid[NCOL], 1, row, pid-1, 0, comm_worker, &reqs[1]);

		if (pid < npr-1) MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
		if (pid > 0) MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);

		if (pid < npr-1) MPI_Irecv(&grid[(tam+1)*NCOL], 1, row, pid+1, 0, comm_worker, &reqs[1]);
		if (pid > 0) MPI_Irecv(&grid[0], 1, row, pid-1, 0, comm_worker, &reqs[0]);

		if (pid < npr-1) MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
		if (pid > 0) MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
		
		// thermal diffusion
		t_local = thermal_diffusion(param, grid, grid_aux, pid, npr, tam);

		// convergence every 10 iterations
		if (niter % 10 == 0) {
			MPI_Allreduce(&t_local, &t_global, 1, MPI_DOUBLE, MPI_SUM, comm_worker);
			t_mean = t_global / ((NCOL - 2) * (NROW - 2));

			if ((fabs(t_mean - t_mean_0) < param.t_delta) || (niter > param.max_iter)) end = 1; 
			else t_mean_0 = t_mean;
		}
	} // end while
	
	return (t_mean);
}