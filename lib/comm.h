#ifndef COMM_H
#define COMM_H 1

void comm_init(int* p_argc, char*** p_argv);
void comm_finalise();
int comm_this_rank();
int comm_n_nodes();
void comm_barrier();
int comm_bcast_int(const int val);
char* comm_bcast_char(char const * const src);
//double* comm_bcast_double(double const * const double_src, const int n_src);
void comm_mpi_bcast_double(double* const p, const int count);

int comm_allreduce_min_int(int val);
void comm_abort();
#endif
