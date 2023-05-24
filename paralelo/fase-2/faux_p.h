/* File: faux.h */ 

extern void read_data (char *, struct info_param *, struct info_chips **, int ***);
extern void results_conf (int, double, struct info_param, float *, struct info_results *);
extern void results (struct info_param, struct info_results *, char *);
extern void create_mpi_params(struct info_param *param, struct MPI_info_param *mpi_param);