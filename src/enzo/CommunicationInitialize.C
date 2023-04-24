/***********************************************************************
/
/  COMMUNICATION ROUTINE: INITIALIZE
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#include <stdlib.h>
#endif /* USE_MPI */

#include <stdio.h>
#include <string.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "communication.h"
#ifdef NBODY
#include "communicators.h"
#endif
 
/* function prototypes */
void my_exit(int exit_status);

#ifdef USE_MPI
void CommunicationErrorHandlerFn(MPI_Comm *comm, MPI_Arg *err, ...);
#ifdef NBODY
#define NumberOfNbodyProcessors 1
	MPI_Comm nbody_comm;
	MPI_Comm inter_comm;
	MPI_Comm enzo_comm;
#endif
#endif



int CommunicationInitialize(Eint32 *argc, char **argv[])
{
 
#ifdef USE_MPI
 
  /* Initialize MPI and get info. */

  MPI_Arg mpi_rank;
  MPI_Arg mpi_size;
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_create_errhandler(CommunicationErrorHandlerFn, &CommunicationErrorHandler);
  MPI_Comm_set_errhandler(comm, CommunicationErrorHandler);

#ifdef NBODY
	// by YS, N body function will be activated only if mpi_size > 12
	if (mpi_size >= 2) {
		// by YS, create the group of processes in MPI_COMM_WORLD
		MPI_Group world_group;
		MPI_Comm_group(MPI_COMM_WORLD, &world_group);

		const int mpi_nbody_size = NumberOfNbodyProcessors;
		const int mpi_enzo_size  = mpi_size-mpi_nbody_size;
		int ranks_nbody[mpi_nbody_size];
		for (int i=0;i<mpi_nbody_size;i++) {
			ranks_nbody[i] = i + mpi_enzo_size;
		}
		int ranks_enzo[mpi_enzo_size];
		for (int i=0;i<mpi_enzo_size;i++) {
			ranks_enzo[i] = i;
		}
		int tag = mpi_rank/(mpi_enzo_size);



		// Construct a enzo group and corresponding communicator
		MPI_Group enzo_group;
		MPI_Group_incl(world_group, mpi_enzo_size, ranks_enzo, &enzo_group);
		//MPI_Comm_create_group(MPI_COMM_WORLD, enzo_group, tag, &enzo_comm);
		MPI_Comm_create(MPI_COMM_WORLD, enzo_group, &enzo_comm);

		// Construct a nbody group and corresponding communicator
		MPI_Group nbody_group;
		MPI_Group_incl(world_group, mpi_nbody_size, ranks_nbody, &nbody_group);
		MPI_Comm_create(MPI_COMM_WORLD, nbody_group, &nbody_comm);


		int new_rank, new_size;
		//if (mpi_rank >= (mpi_size-mpi_nbody_size)) {
		if (MPI_COMM_NULL != nbody_comm) {
			MPI_Comm_rank(nbody_comm, &new_rank);
			MPI_Comm_size(nbody_comm, &new_size);
			fprintf(stderr,"nbody: (%d, %d)\n", new_rank,new_size);
			NumberOfProcessors = mpi_nbody_size;
			MPI_Intercomm_create(nbody_comm, 0, MPI_COMM_WORLD, 0, 
					 10, &inter_comm); 
		}
		if (MPI_COMM_NULL != enzo_comm) {
			MPI_Comm_rank(enzo_comm, &new_rank);
			MPI_Comm_size(enzo_comm, &new_size);
			fprintf(stderr,"enzo: (%d, %d)\n", new_rank,new_size);
			NumberOfProcessors = mpi_size-mpi_nbody_size;
			MPI_Intercomm_create(enzo_comm, 0, MPI_COMM_WORLD, mpi_enzo_size, 
					 10, &inter_comm); 
		}
		TotalNumberOfProcessors = mpi_size;
	}
#else
  NumberOfProcessors = mpi_size;
#endif
  MyProcessorNumber = mpi_rank;
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("MPI_Init: NumberOfProcessors = %"ISYM"\n", NumberOfProcessors);
    printf("MPI_Init: TotalNumberOfProcessors = %"ISYM"\n", TotalNumberOfProcessors);
	}
 
#else /* USE_MPI */
 
  MyProcessorNumber  = 0;
  NumberOfProcessors = 1;
 
#endif /* USE_MPI */
 
  CommunicationTime = 0;
 
  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  return SUCCESS;
}
 
#ifdef USE_MPI
void CommunicationErrorHandlerFn(MPI_Comm *comm, MPI_Arg *err, ...)
{
  char error_string[1024];
  MPI_Arg length, error_class;
  if (*err != MPI_ERR_OTHER) {
      MPI_Error_class(*err, &error_class);
      MPI_Error_string(error_class, error_string, &length);
      fprintf(stderr, "P%"ISYM": %s\n", MyProcessorNumber, error_string);
      MPI_Error_string(*err, error_string, &length);
      fprintf(stderr, "P%"ISYM": %s\n", MyProcessorNumber, error_string);
      ENZO_FAIL("MPI communication error.");
  } // ENDIF MPI_ERROR
  return;
}
#endif /* USE_MPI */
 
int CommunicationFinalize()
{
 
#ifdef USE_MPI
#ifdef NBODY
	MPI_Comm_free(&enzo_comm);
	MPI_Comm_free(&nbody_comm);
	MPI_Comm_free(&inter_comm);
#endif
  MPI_Errhandler_free(&CommunicationErrorHandler);
  MPI_Finalize();
#endif /* USE_MPI */
 
  return SUCCESS;
}

void CommunicationAbort(int status)
{

#ifdef USE_MPI
  MPI_Abort(MPI_COMM_WORLD,status);
#else
  //  my_exit(status);
#endif

  return;
}
