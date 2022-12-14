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



extern "C" void FORTRAN_NAME(nbody6)();
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		HierarchyEntry **Grids[]);
//int FindTotalNumberOfNbodyParticles(LevelHierarchyEntry *LevelArray[]);
//int FindStartIndex(int* LocalNumberOfNbodyParticles);
//void InitializeNbodyArrays(void);
//void DeleteNbodyArrays(void);
//void CopyNbodyArrayToOld(void);
//void MatchAccelerationWithIndex(void);


#define NBODY
#ifdef NBODY
int PrepareNbodyComputation(LevelHierarchyEntry *LevelArray[], int level)
{

	int i, GridNum, LocalNumberOfNbodyParticles;
	int SavedP3IMFCalls;
	Star *LocalStars = NULL, *GridStars = NULL, *cstar = NULL, *lstar = NULL;
	LevelHierarchyEntry *Temp;
	int *NumberOfStarsInGrids;
	int start_index;

	if (MyProcessorNumber == ROOT_PROCESSOR && NbodyFirst) {
		NbodyParticleIDOld =  NULL;
		NumberOfNbodyParticlesOld = 0;
		for (int dim=0; dim<MAX_DIMENSION; dim++) 
			for (int order=0; order<HERMITE_ORDER; order++)
				NbodyParticleAccelerationOld[dim][order] = NULL;
		NbodyFirst=FALSE;
		fprintf(stderr,"Proc %d: All Good 0, NbodyFirst=%d\n", MyProcessorNumber,NbodyFirst);
		fprintf(stderr,"ID Old = %d\n",NbodyParticleIDOld);
		fprintf(stderr,"Num Old = %d\n",NumberOfNbodyParticlesOld);
	}

	fprintf(stderr,"All Good 1\n");

	/* Store AccelerationNoStar of particles to sustainable ParticleAttribute */
	for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
		Temp->GridData->CopyAccelerationToAttribute();

	fprintf(stderr,"All Good 2\n");


	if (level == MaximumRefinementLevel) {

			fprintf(stderr,"ID Old = %d\n",NbodyParticleIDOld);
			fprintf(stderr,"Num Old = %d\n",NumberOfNbodyParticlesOld);

		LocalNumberOfNbodyParticles = FindTotalNumberOfNbodyParticles(LevelArray);
		fprintf(stderr,"Local Number Of Nbody Particles = %d\n", LocalNumberOfNbodyParticles);
		fprintf(stderr,"Total Number Of Nbody Particles = %d\n", NumberOfNbodyParticles);

			fprintf(stderr,"ID Old = %d\n",NbodyParticleIDOld);
			fprintf(stderr,"Num Old = %d\n",NumberOfNbodyParticlesOld);

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
			NbodyParticleAccelerationNoStarTemp[dim] = new float[LocalNumberOfNbodyParticles];
		}



		/* Get particle information from Grids */
		int count = 0;
		for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++)
			for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel) 
				if (Temp->GridData->CopyNbodyParticles(&count, NbodyParticleIDTemp, NbodyParticleMassTemp,
							NbodyParticlePositionTemp, NbodyParticleVelocityTemp, NbodyParticleAccelerationNoStarTemp) == FAIL) {
					//NbodyParticleAccelerationTemp
					ENZO_FAIL("Error in grid::CopyNbodyParticles.");
				}

		if (count != LocalNumberOfNbodyParticles)
			ENZO_FAIL("Error in grid::CopyNbodyParticles.");

		/* Find the index of the array */
		start_index = FindStartIndex(&LocalNumberOfNbodyParticles);

		fprintf(stderr,"Done?1\n");

		/* We need to gather all the information from other processes */

		if (MyProcessorNumber == ROOT_PROCESSOR) {

			/* Initialize the nbody array used for direct Nbody calculation*/
			InitializeNbodyArrays();

			fprintf(stderr,"ID Old = %d\n", NbodyParticleIDOld);
			fprintf(stderr,"Num Old = %d\n",NumberOfNbodyParticlesOld);

			/* Receiving Index, NumberOfParticles, NbodyArrays from other processs */
#ifdef USE_MPI
			int* start_index_all;
			int* LocalNumberAll;
			start_index_all = new int[NumberOfProcessors];
			LocalNumberAll = new int[NumberOfProcessors];

			MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, LocalNumberAll, 1, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);
			MPI_Gather(&start_index, 1, IntDataType, start_index_all, 1, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);
			for (i=0;i<NumberOfProcessors;i++) {
				fprintf(stderr,"(%d, %d)\n",start_index_all[i],LocalNumberAll[i]);
			}
			fprintf(stderr,"ID Old = %d\n",NbodyParticleIDOld);
			fprintf(stderr,"Num Old = %d\n",NumberOfNbodyParticlesOld);

			MPI_Gatherv(NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_FLOAT,
					NbodyParticleMass, LocalNumberAll, start_index_all, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
			MPI_Gatherv(NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType,
					NbodyParticleID, LocalNumberAll, start_index_all, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);

			fprintf(stderr,"ID Old = %d\n",NbodyParticleIDOld);
			fprintf(stderr,"Num Old = %d\n",NumberOfNbodyParticlesOld);
			fprintf(stderr,"Done?1-1\n");

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Gatherv(NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT,
						NbodyParticlePosition[dim], LocalNumberAll, start_index_all, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);

				MPI_Gatherv(NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT,
						NbodyParticleVelocity[dim], LocalNumberAll, start_index_all, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);

				MPI_Gatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT,
						NbodyParticleAccelerationNoStar[dim], LocalNumberAll, start_index_all, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
			}
			fprintf(stderr,"Done?1-2\n");


			fprintf(stderr,"ID Old = %d\n",NbodyParticleIDOld);
			fprintf(stderr,"Num Old = %d\n",NumberOfNbodyParticlesOld);
			//if (NbodyParticleIDOld != NULL) {
			if (NumberOfNbodyParticlesOld != 0) {
				/* Now Assign Acceleration History According to ParticleID */
				fprintf(stderr,"Done?1-2-1\n");
				MatchAccelerationWithIndex();
			}

			fprintf(stderr,"Done?1-2-2\n");

			/* Do direct calculation!*/
			nbody6();


			/* Copy ID and Acc to Old arrays!*/
			CopyNbodyArrayToOld();
			fprintf(stderr,"Done?2\n");


			CommunicationBarrier();

			/* Sending Index, NumberOfParticles, NbodyArrays to other processs */
			for (i=0;i<NumberOfProcessors;i++) {
				fprintf(stderr,"(%d, %d)",start_index_all[i],LocalNumberAll[i]);
			}

			MPI_Scatterv(NbodyParticleMass, LocalNumberAll, start_index_all, MPI_FLOAT,
					NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
			MPI_Scatterv(NbodyParticleID, LocalNumberAll, start_index_all, IntDataType,
					NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Scatterv(NbodyParticlePosition[dim], LocalNumberAll, start_index_all, MPI_FLOAT,
						NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
				MPI_Scatterv(NbodyParticleVelocity[dim], LocalNumberAll, start_index_all, MPI_FLOAT,
						NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
			}

			delete [] start_index_all;
			delete [] LocalNumberAll;
			DeleteNbodyArrays();


		} else {
			/* Sending Index, NumberOfParticles, NbodyArrays to the root processs */

			MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);
			MPI_Gather(&start_index, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);

			MPI_Gatherv(NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_FLOAT,
					NULL, NULL, NULL, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
			MPI_Gatherv(NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType,
					NULL, NULL, NULL, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Gatherv(NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT,
						NULL, NULL, NULL, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
				MPI_Gatherv(NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT,
						NULL, NULL, NULL, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
				MPI_Gatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT,
						NULL, NULL, NULL, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
			}

			CommunicationBarrier();

			/* Receiving Index, NumberOfParticles, NbodyArrays from the root processs */
			MPI_Scatterv(NULL, NULL, NULL, MPI_FLOAT,
					NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
			MPI_Scatterv(NULL, NULL, NULL, IntDataType,
					NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Scatterv(NULL, NULL, NULL, MPI_FLOAT,
						NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
				MPI_Scatterv(NULL, NULL, NULL, MPI_FLOAT,
						NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_FLOAT, ROOT_PROCESSOR, MPI_COMM_WORLD);
			}
#endif
		}

		fprintf(stderr,"Done?3\n");
		/* Update Particle Velocity and Position Back to Grids */
		count = 0;
		for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++)
			for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel)
				if (Temp->GridData->UpdateNbodyParticles(&count, NbodyParticleIDTemp, NbodyParticleMassTemp,
							NbodyParticlePositionTemp, NbodyParticleVelocityTemp) == FAIL) {
					ENZO_FAIL("Error in grid::CopyNbodyParticles.");
				}

		fprintf(stderr,"Done?4\n");



		/* Destruct Arrays*/
		delete [] NbodyParticleIDTemp;
		delete [] NbodyParticleMassTemp;
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			delete [] NbodyParticlePositionTemp[dim];
			delete [] NbodyParticleVelocityTemp[dim];
			delete [] NbodyParticleAccelerationNoStarTemp[dim];
		}

		/* somewhere here I should include direct nbody codes */

		/* update back the particle positions and velocities */

	} // ENDIF level


	return SUCCESS;

}
#endif

