#ifdef WITHMPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <cstring>
#include <cassert>
#include "msg.h"
#include "comm.h"

static int this_rank= -1;
static int n_nodes= 0;

#ifdef WITHMPI
void comm_init(int* p_argc, char*** p_argv)
{
  int ret= MPI_Init(p_argc, p_argv); assert(ret == MPI_SUCCESS);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
}


void comm_finalise()
{
  MPI_Barrier(MPI_COMM_WORLD);
  msg_printf(msg_verbose, "Finishing program with MPI_Finalize()\n");
  MPI_Finalize();
}

void comm_barrier()
{
  MPI_Barrier(MPI_COMM_WORLD);
}

int comm_bcast_int(const int val)
{
  int val_comm= val;
  MPI_Bcast(&val_comm, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return val_comm;
}

char* comm_bcast_char(char const * const src)
{
  int n=0;
  if(this_rank == 0)
    n= strlen(src);
  n= comm_bcast_int(n);

  char* const str= (char*) malloc(n+1); assert(str);
  strncpy(str, src, n+1);
  
  int ret= MPI_Bcast(str, n, MPI_CHAR, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  return str;
}

int comm_allreduce_min_int(int val)
{
  int val_reduced;
  MPI_Allreduce(&val, &val_reduced, 1, MPI_INT, MPI_MIN,  MPI_COMM_WORLD);
  return val_reduced;
}

/*
double* comm_bcast_double(double* const double_src, const int n)
{
  double* const double_array= (double*) malloc(sizeof(double)*n);
  assert(double_array);
  for(int i=0; i<n; ++i)
    double_array[i]= double_src[i];
  
  MPI_Bcast(double_array, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  free(double_src);
  
  return double_array;
}
*/

void comm_mpi_bcast_double(double* const p, const int count)
{
  int ret= MPI_Bcast(p, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);
}

void comm_abort()
{
  MPI_Abort(MPI_COMM_WORLD, 1);
}

#else
//
// Serial: dummy communication functions
//
void comm_init(int* p_argc, char*** p_argv)
{
  this_rank= 0;
  n_nodes= 1;

  msg_printf(msg_verbose, "Serial mode\n");
}
void comm_finalise() {}
void comm_barrier() {}

char* comm_bcast_char(char const * const src)
{
  const size_t n= strlen(src);
  char* const str= (char*) malloc(n+1); assert(str);
  strncpy(str, src, n+1);

  return str;
}

/*
double* comm_bcast_double(double * const double_src, const int n_src)
{
  return double_src;
}
*/

 //int comm_mpi_bcast_double(double* const p, const int count){};
int comm_bcast_int(const int val)
{
  return val;
}

void comm_mpi_bcast_double(double* const p, const int count)
{
}

int comm_allreduce_min_int(int val)
{
  return val;
}


void comm_abort()
{
  abort();
}


#endif


int comm_this_rank()
{
  return this_rank;
}


int comm_n_nodes()
{
  return n_nodes;
}


