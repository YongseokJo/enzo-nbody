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


void InitializeNbodyArrays(int);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		HierarchyEntry **Grids[]);
int GetUnits(float *DensityUnits, float *LengthUnits,
		float *TemperatureUnits, float *TimeUnits,
		float *VelocityUnits, double *MassUnits, FLOAT Time);

#ifdef NBODY
int FinalizeNbodyComputation(LevelHierarchyEntry *LevelArray[], int level)
{


	if (level == MaximumRefinementLevel) {
		int i, GridNum, LocalNumberOfNbodyParticles;
		LevelHierarchyEntry *Temp;
		int start_index;

		LocalNumberOfNbodyParticles = FindTotalNumberOfNbodyParticles(LevelArray);

		int *NbodyParticleIDTemp;
		float *NbodyParticlePositionTemp[MAX_DIMENSION];
		float *NbodyParticleVelocityTemp[MAX_DIMENSION];

		NbodyParticleIDTemp   = new int[LocalNumberOfNbodyParticles];
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			NbodyParticlePositionTemp[dim]            = new float[LocalNumberOfNbodyParticles];
			NbodyParticleVelocityTemp[dim]            = new float[LocalNumberOfNbodyParticles];
		}


		fprintf(stderr,"FNC in\n");
		/* Find the index of the array */
		start_index = FindStartIndex(&LocalNumberOfNbodyParticles);
		fprintf(stderr,"Fianl Done?1\n");

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


			for (i=0;i<NumberOfProcessors;i++) {
				fprintf(stderr,"proc=%d: (%d, %d)\n",i, start_index_all[i],LocalNumberAll[i]);
			}

			fprintf(stderr,"Done?2\n");



			  /*-----------------------------------------------*/
			 /******** Recv Arrays to Fortran Nbody6++    *****/
			/*-----------------------------------------------*/
			InitializeNbodyArrays(1);
			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Recv(NbodyParticlePosition[dim], NumberOfNbodyParticles, MPI_DOUBLE, NumberOfProcessors, 300, MPI_COMM_WORLD, &status);
				MPI_Recv(NbodyParticleVelocity[dim], NumberOfNbodyParticles, MPI_DOUBLE, NumberOfProcessors, 400, MPI_COMM_WORLD, &status);
			}

			fprintf(stderr,"NumberOfParticles after NBODY=%d\n",NumberOfNbodyParticles);
			fprintf(stderr,"enzo: X=%f, V=%f\n ",NbodyParticlePosition[0], NbodyParticleVelocity[0]);


			/* Sending Index, NumberOfParticles, NbodyArrays to other processs */
			/*
			MPI_Iscatterv(NbodyParticleID, LocalNumberAll, start_index_all, IntDataType,
					NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
			*/

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
		fprintf(stderr,"Done?2-5\n");

			if (start_index_all != NULL)
				delete [] start_index_all;
			start_index_all = NULL;
			if (LocalNumberAll != NULL)
				delete [] LocalNumberAll;
			LocalNumberAll = NULL;
			DeleteNbodyArrays();

		}  // ENDIF: Root processor
		else {

			MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
			MPI_Gather(&start_index, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);


			MPI_Request request;
			MPI_Status status;


		fprintf(stderr,"Done?2-5\n");
			/* Receiving Index, NumberOfParticles, NbodyArrays from the root processs */
			/*
			MPI_Iscatterv(NULL, NULL, NULL, IntDataType,
					NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType, ROOT_PROCESSOR, enzo_comm,
					&request);
			MPI_Wait(&request, &status);
			*/

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
		} // end else
#endif

		//CommunicationBarrier();
		fprintf(stderr,"Done?3\n");

		/* Update Particle Velocity and Position Back to Grids */
		int count = 0;
		for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++)
			for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel)
				if (Temp->GridData->UpdateNbodyParticles(&count, NbodyParticleIDTemp,
							NbodyParticlePositionTemp, NbodyParticleVelocityTemp) == FAIL) {
					ENZO_FAIL("Error in grid::CopyNbodyParticles.");
				}

		fprintf(stderr,"Done?4\n");



		/* Destruct Arrays*/
		if (NbodyParticleIDTemp != NULL)
			delete [] NbodyParticleIDTemp;
		NbodyParticleIDTemp = NULL;

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			if (NbodyParticlePositionTemp[dim] != NULL)
				delete [] NbodyParticlePositionTemp[dim];
			NbodyParticlePositionTemp[dim] = NULL;

			if (NbodyParticleVelocityTemp[dim] != NULL)
				delete [] NbodyParticleVelocityTemp[dim];
			NbodyParticleVelocityTemp[dim] = NULL;
		}
		fprintf(stderr,"Done?10\n");

		} // ENDIF level
		return SUCCESS;
	}
#endif

