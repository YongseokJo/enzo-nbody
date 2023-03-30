/***********************************************************************
	/
	/  FIND ALL NBODY PARTICLES OVER ALL PROCESSORS
	/
	/  written by: Yongseok Jo
	/  date:       November, 2022
	/
	/  PURPOSE: First synchronizes particle information in the normal and 
	/           nbody particles.  Then we make a global particle list, which
	/           simplifies Nbody calculations.
	/
 ************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "NbodyRoutines.h"  //added  
#include "phys_constants.h"


int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		HierarchyEntry **Grids[]);
int GetUnits(float *DensityUnits, float *LengthUnits,
		float *TemperatureUnits, float *TimeUnits,
		float *VelocityUnits, double *MassUnits, FLOAT Time);

#ifdef NBODY
int FinalizeNbodyComputation(LevelHierarchyEntry *LevelArray[], int level)
{

	int i, GridNum, LocalNumberOfNbodyParticles;
	LevelHierarchyEntry *Temp;
	int start_index;

	if (level == MaximumRefinementLevel) {

		LocalNumberOfNbodyParticles = FindTotalNumberOfNbodyParticles(LevelArray);

		int *NbodyParticleIDTemp;
		float *NbodyParticleMassTemp;
		float *NbodyParticlePositionTemp[MAX_DIMENSION];
		float *NbodyParticleVelocityTemp[MAX_DIMENSION];
		float *NbodyParticleAccelerationNoStarTemp[MAX_DIMENSION];

		NbodyParticleIDTemp   = new int[LocalNumberOfNbodyParticles];
		NbodyParticleMassTemp = new float[LocalNumberOfNbodyParticles];
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			NbodyParticlePositionTemp[dim]            = new float[LocalNumberOfNbodyParticles];
			NbodyParticleVelocityTemp[dim]            = new float[LocalNumberOfNbodyParticles];
			NbodyParticleAccelerationNoStarTemp[dim]  = new float[LocalNumberOfNbodyParticles];
		}


		/* Find the index of the array */
		start_index = FindStartIndex(&LocalNumberOfNbodyParticles);

		fprintf(stderr,"Done?1\n");

#ifdef USE_MPI
		if (MyProcessorNumber == ROOT_PROCESSOR) {

			/* Receiving Index, NumberOfParticles, NbodyArrays from other processs */
			int* start_index_all;
			int* LocalNumberAll;
			start_index_all = new int[NumberOfProcessors];
			LocalNumberAll = new int[NumberOfProcessors];
			MPI_Request request;
			MPI_Status status;


			MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, LocalNumberAll, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
			MPI_Gather(&start_index, 1, IntDataType, start_index_all, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);

			CommunicationBarrier();

			for (i=0;i<NumberOfProcessors;i++) {
				fprintf(stderr,"proc=%d: (%d, %d)\n",i, start_index_all[i],LocalNumberAll[i]);
			}

			fprintf(stderr,"Done?2\n");



			  /*-----------------------------------------------*/
			 /******** Recv Arrays to Fortran Nbody6++    *****/
			/*-----------------------------------------------*/
			MPI_Recv(&NumberOfNbodyParticles, 1, MPI_INT, NumberOfProcessors, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(NbodyParticleMass, NumberOfNbodyParticles, MPI_DOUBLE, NumberOfProcessors, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Recv(NbodyParticlePosition[dim], NumberOfNbodyParticles, MPI_DOUBLE, NumberOfProcessors, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(NbodyParticleVelocity[dim], NumberOfNbodyParticles, MPI_DOUBLE, NumberOfProcessors, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(NbodyParticleAccelerationNoStar[dim], NumberOfNbodyParticles, MPI_DOUBLE, NumberOfProcessors, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int order=0; order<HERMITE_ORDER;order++)
					MPI_Recv(NbodyParticleAcceleration[dim][order], NumberOfNbodyParticles, MPI_DOUBLE, NumberOfProcessors, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			fprintf(stderr,"NumberOfParticles after NBODY=%d\n",NumberOfNbodyParticles);

			  /*-----------------------------------------------*/
			 /******** Copy ID and Acc to Old arrays! *********/
			/*-----------------------------------------------*/
			CopyNbodyArrayToOld();
			fprintf(stderr,"Done?2\n");


			CommunicationBarrier();

			/* Sending Index, NumberOfParticles, NbodyArrays to other processs */
			for (i=0;i<NumberOfProcessors;i++) {
				fprintf(stderr,"(%d, %d)",start_index_all[i],LocalNumberAll[i]);
			}

			MPI_Iscatterv(NbodyParticleMass, LocalNumberAll, start_index_all, MPI_DOUBLE,
					NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
			MPI_Iscatterv(NbodyParticleID, LocalNumberAll, start_index_all, IntDataType,
					NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Iscatterv(NbodyParticlePosition[dim], LocalNumberAll, start_index_all, MPI_DOUBLE,
						NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
						&request);
				MPI_Wait(&request, &status);
				MPI_Iscatterv(NbodyParticleVelocity[dim], LocalNumberAll, start_index_all, MPI_DOUBLE,
						NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
						&request);
				MPI_Wait(&request, &status);
			}

			CommunicationBarrier();
			if (start_index_all != NULL)
				delete [] start_index_all;
			start_index_all = NULL;
			if (LocalNumberAll != NULL)
				delete [] LocalNumberAll;
			LocalNumberAll = NULL;
			DeleteNbodyArrays();

		} else {

			MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
			MPI_Gather(&start_index, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);

			CommunicationBarrier();

			MPI_Request request;
			MPI_Status status;

			CommunicationBarrier();

			/* Receiving Index, NumberOfParticles, NbodyArrays from the root processs */
			MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
					NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
					&request);
			MPI_Wait(&request, &status);
			MPI_Iscatterv(NULL, NULL, NULL, IntDataType,
					NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType, ROOT_PROCESSOR, enzo_comm,
					&request);
			MPI_Wait(&request, &status);

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
						NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
						&request);
				MPI_Wait(&request, &status);
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
						NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
						&request);
				MPI_Wait(&request, &status);
			}
			CommunicationBarrier();
		} // end else
#endif


		/* Update Particle Velocity and Position Back to Grids */
		int count = 0;
		for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++)
			for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel)
				if (Temp->GridData->UpdateNbodyParticles(&count, NbodyParticleIDTemp, NbodyParticleMassTemp,
							NbodyParticlePositionTemp, NbodyParticleVelocityTemp) == FAIL) {
					ENZO_FAIL("Error in grid::CopyNbodyParticles.");
				}

		fprintf(stderr,"Done?4\n");



		/* Destruct Arrays*/
		if (NbodyParticleIDTemp != NULL)
			delete [] NbodyParticleIDTemp;
		NbodyParticleIDTemp = NULL;
		fprintf(stderr,"Done?5\n");

		if (NbodyParticleMassTemp != NULL)
			delete [] NbodyParticleMassTemp;
		NbodyParticleMassTemp = NULL;
		fprintf(stderr,"Done?6\n");

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			if (NbodyParticlePositionTemp[dim] != NULL)
				delete [] NbodyParticlePositionTemp[dim];
			NbodyParticlePositionTemp[dim] = NULL;
			fprintf(stderr,"Done?7\n");

			if (NbodyParticleVelocityTemp[dim] != NULL)
				delete [] NbodyParticleVelocityTemp[dim];
			NbodyParticleVelocityTemp[dim] = NULL;
			fprintf(stderr,"Done?8\n");

			if (NbodyParticleAccelerationNoStarTemp[dim] != NULL)
				delete [] NbodyParticleAccelerationNoStarTemp[dim];
			NbodyParticleAccelerationNoStarTemp[dim] = NULL;
			fprintf(stderr,"Done?9\n");
		}
			CommunicationBarrier();
		fprintf(stderr,"Done?10\n");

		} // ENDIF level
		return SUCCESS;
	}
#endif
