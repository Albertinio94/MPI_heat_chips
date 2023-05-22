/**********************************************************
 * defines.h
 * definitions for the thermal simulation of the card
 ***********************************************************/

#include <mpi.h>

// minimal card and maximum size
#define RSIZE 200
#define CSIZE 100

#define NROW (RSIZE * param.scale + 2) // extended row number
#define NCOL (CSIZE * param.scale + 2) // extended column number

#define N_PARAMS 7 // number of general parameters
#define ROOT 0

struct info_param {
	// num. of configurations to test, num. of chips, max. num. iterations, card size scale
	int nconf, nchip, max_iter, scale; 
	// external temp., max. temp. of a chip, temp. incr. for convergence
	float t_ext, tmax_chip, t_delta;
};

struct info_chips {
	int h, w; // size (h, w)
	float tchip; // temperature
};

struct info_results {
	double Tmean; // mean temp.
	int conf; // conf number
	float *bgrid; // final grid
	float *cgrid; // initial grid (chips)
};

struct MPI_info_param {
	MPI_Datatype param;
	int sizes[N_PARAMS];
	MPI_Aint distances[N_PARAMS];
	MPI_Datatype types[N_PARAMS];
};
