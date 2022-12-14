
#ifdef USE_MPI
#include "mpi.h"
#include "communicators.h"
#endif
#include <stdlib.h>
#include <stdio.h>

void Reduce_Times(double time, double *time_array){
  int nprocs, my_rank; 
#ifdef USE_MPI
  MPI_Comm_size(enzo_comm, &nprocs);
  MPI_Comm_rank(enzo_comm, &my_rank);

  MPI_Gather(&time, 1, MPI_DOUBLE, 
	     time_array, 1, MPI_DOUBLE, 
	     0, enzo_comm);
#else
  time_array[0] = time;
#endif

  return;
}
