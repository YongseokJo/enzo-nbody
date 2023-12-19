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


	if (LevelArray[level+1] == NULL) {
		int i, GridNum, LocalNumberOfNbodyParticles, NewLocalNumberOfNbodyParticles;
		LevelHierarchyEntry *Temp;
		int start_index, start_index_new;

		FindTotalNumberOfNbodyParticles(LevelArray, &LocalNumberOfNbodyParticles, &NewLocalNumberOfNbodyParticles);

		if (NumberOfNbodyParticles == 0)
			return SUCCESS;

		float *NbodyParticlePositionTemp[MAX_DIMENSION];
		float *NbodyParticleVelocityTemp[MAX_DIMENSION];

		float *NewNbodyParticlePositionTemp[MAX_DIMENSION];
		float *NewNbodyParticleVelocityTemp[MAX_DIMENSION];

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			NbodyParticlePositionTemp[dim]            = new float[LocalNumberOfNbodyParticles];
			NbodyParticleVelocityTemp[dim]            = new float[LocalNumberOfNbodyParticles];
		}

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			NewNbodyParticlePositionTemp[dim]            = new float[NewLocalNumberOfNbodyParticles];
			NewNbodyParticleVelocityTemp[dim]            = new float[NewLocalNumberOfNbodyParticles];
		}


		/* Find the index of the array */
		start_index = FindStartIndex(&LocalNumberOfNbodyParticles);
		start_index_new = FindStartIndex(&NewLocalNumberOfNbodyParticles);

#ifdef USE_MPI
		if (MyProcessorNumber == ROOT_PROCESSOR) {

			//float *NewNbodyParticleMass;
			float *NewNbodyParticleVelocity[MAX_DIMENSION]; // feedback can affect velocity
			float *NewNbodyParticlePosition[MAX_DIMENSION]; // feedback can affect velocity

			//NewNbodyParticleMass = new float[NumberOfNewNbodyParticles];
			for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
				NewNbodyParticlePosition[dim]            = new float[NumberOfNewNbodyParticles];
				NewNbodyParticleVelocity[dim]            = new float[NumberOfNewNbodyParticles];
			}

			/* Receiving Index, NumberOfParticles, NbodyArrays from other processs */
			int* start_index_all;
			int* LocalNumberAll;
			start_index_all = new int[NumberOfProcessors];
			LocalNumberAll = new int[NumberOfProcessors];

			int* start_index_all_new;
			int* NewLocalNumberAll;
			start_index_all_new = new int[NumberOfProcessors];
			NewLocalNumberAll = new int[NumberOfProcessors];

			MPI_Request request;
			MPI_Status status;
			int ierr;
			int errclass,resultlen;
			char err_buffer[MPI_MAX_ERROR_STRING];

			MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, LocalNumberAll, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
			MPI_Gather(&start_index, 1, IntDataType, start_index_all, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);


			if (NumberOfNewNbodyParticles > 0)  {
				MPI_Gather(&NewLocalNumberOfNbodyParticles, 1, IntDataType, NewLocalNumberAll, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index_new, 1, IntDataType, start_index_all_new, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
			}



			  /*-----------------------------------------------*/
			 /******** Recv Arrays to Fortran Nbody6++    *****/
			/*-----------------------------------------------*/
			InitializeNbodyArrays(1);
			CommunicationInterBarrier();
			fprintf(stderr,"NumberOfParticles=%d\n",NumberOfNbodyParticles);
			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				ierr = MPI_Recv(NbodyParticlePosition[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 300, inter_comm, &status);
				ierr = MPI_Recv(NbodyParticleVelocity[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 400, inter_comm, &status);
			}
			fprintf(stderr,"NewNumberOfParticles=%d\n",NumberOfNewNbodyParticles);
			if (NumberOfNewNbodyParticles > 0) 
				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					ierr = MPI_Recv(NewNbodyParticlePosition[dim], NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 500, inter_comm, &status);
					ierr = MPI_Recv(NewNbodyParticleVelocity[dim], NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 600, inter_comm, &status);
				}
			CommunicationInterBarrier();

			fprintf(stderr,"NumberOfParticles after NBODY=%d\n",NumberOfNbodyParticles);
			fprintf(stderr,"enzo: X=%e, V=%e\n ",NbodyParticlePosition[0][0], NbodyParticleVelocity[0][0]);


			/* Sending Index, NumberOfParticles, NbodyArrays to other processs */
			/*
			MPI_Iscatterv(NbodyParticleID, LocalNumberAll, start_index_all, IntDataType,
					NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
			*/

			//CommunicationBarrier();
			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Iscatterv(NbodyParticlePosition[dim], LocalNumberAll, start_index_all, MPI_DOUBLE,
						NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
						&request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NbodyParticleVelocity[dim], LocalNumberAll, start_index_all, MPI_DOUBLE,
						NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
						&request);
				ierr  = MPI_Wait(&request, &status);
			}
			if (NumberOfNewNbodyParticles > 0) 
				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					MPI_Iscatterv(NewNbodyParticlePosition[dim], NewLocalNumberAll, start_index_all_new, MPI_DOUBLE,
							NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							&request);
					ierr = MPI_Wait(&request, &status);
					MPI_Iscatterv(NewNbodyParticleVelocity[dim], NewLocalNumberAll, start_index_all_new, MPI_DOUBLE,
							NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							&request);
					ierr  = MPI_Wait(&request, &status);
				}



				fprintf(stderr,"Root:Done?2-5\n");

		if (start_index_all != NULL)
			delete [] start_index_all;
		start_index_all = NULL;
		if (start_index_all != NULL)
			delete [] start_index_all_new;
		start_index_all_new = NULL;
		if (LocalNumberAll != NULL)
			delete [] LocalNumberAll;
		LocalNumberAll = NULL;
		if (NewLocalNumberAll != NULL)
			delete [] NewLocalNumberAll;
		NewLocalNumberAll = NULL;

		/*
		if (NewNbodyParticleMass != NULL)
			delete [] NewNbodyParticleMass;
		NewNbodyParticleMass = NULL;
		*/


		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			if (NewNbodyParticlePosition[dim] != NULL)
				delete [] NewNbodyParticlePosition[dim];
			NewNbodyParticlePosition[dim] = NULL;

			if (NewNbodyParticleVelocity[dim] != NULL)
				delete [] NewNbodyParticleVelocity[dim];
			NewNbodyParticleVelocity[dim] = NULL;
		}


			DeleteNbodyArrays();

		}  // ENDIF: Root processor
		else {

			MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
			MPI_Gather(&start_index, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);

			if (NumberOfNewNbodyParticles > 0)  {
				MPI_Gather(&NewLocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index_new, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
			}

			MPI_Request request;
			MPI_Status status;
			int ierr;


			/* Receiving Index, NumberOfParticles, NbodyArrays from the root processs */
			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
						NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
						&request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
						NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
						&request);
				ierr = MPI_Wait(&request, &status);
			}
			if (NumberOfNewNbodyParticles > 0) 
				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							NewNbodyParticlePositionTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							&request);
					ierr = MPI_Wait(&request, &status);
					MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							NewNbodyParticleVelocityTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							&request);
					ierr = MPI_Wait(&request, &status);
				}

		} // end else
#endif


		/* Update Particle Velocity and Position Back to Grids */
		int count = 0;
		for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++)
			for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel)
				if (Temp->GridData->UpdateNbodyParticles(&count, 
							LocalNumberOfNbodyParticles, NbodyParticleIDTemp, 
							NbodyParticlePositionTemp, NbodyParticleVelocityTemp,
							NewLocalNumberOfNbodyParticles, NewNbodyParticleIDTemp, 
							NewNbodyParticlePositionTemp, NewNbodyParticleVelocityTemp
							) == FAIL) {
					ENZO_FAIL("Error in grid::CopyNbodyParticles.");
				}

		fprintf(stderr,"Proc %d, # of Nbody = %d, count = %d\n", 
				MyProcessorNumber, LocalNumberOfNbodyParticles, count);



		/* Destruct Arrays*/
		if (NbodyParticleIDTemp != NULL)
			delete [] NbodyParticleIDTemp;
		NbodyParticleIDTemp = NULL;

		if (NewNbodyParticleIDTemp != NULL)
			delete [] NewNbodyParticleIDTemp;
		NewNbodyParticleIDTemp = NULL;

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			if (NbodyParticlePositionTemp[dim] != NULL)
				delete [] NbodyParticlePositionTemp[dim];
			NbodyParticlePositionTemp[dim] = NULL;

			if (NbodyParticleVelocityTemp[dim] != NULL)
				delete [] NbodyParticleVelocityTemp[dim];
			NbodyParticleVelocityTemp[dim] = NULL;
		}

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			if (NewNbodyParticlePositionTemp[dim] != NULL)
				delete [] NewNbodyParticlePositionTemp[dim];
			NewNbodyParticlePositionTemp[dim] = NULL;

			if (NewNbodyParticleVelocityTemp[dim] != NULL)
				delete [] NewNbodyParticleVelocityTemp[dim];
			NewNbodyParticleVelocityTemp[dim] = NULL;
		}
		NumberOfNewNbodyParticles = 0;

		} // ENDIF level
		return SUCCESS;
	}
#endif

