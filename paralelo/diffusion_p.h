/* File: diffusion.h */ 
#include <mpi.h>

double calculate_Tmean(struct info_param param, float *grid, float *grid_chips, float *grid_aux, int pid, int npr, int tam, MPI_Datatype row);