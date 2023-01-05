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



int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		HierarchyEntry **Grids[]);
int FindTotalNumberOfNbodyParticles(LevelHierarchyEntry *LevelArray[]);
int FindStartIndex(int LocalNumberOfNbodyParticles);
void InitializeNbodyArrays(void);


#define NBODY
#ifdef NBODY
int NbodyParticleFindAll(LevelHierarchyEntry *LevelArray[], int true_level)
{

	int i, level, GridNum, LocalNumberOfNbodyParticles;
	int SavedP3IMFCalls;
	Star *LocalStars = NULL, *GridStars = NULL, *cstar = NULL, *lstar = NULL;
	LevelHierarchyEntry *Temp;
	int *NumberOfStarsInGrids;
	int start_index;

	LocalNumberOfNbodyParticles = FindTotalNumberOfNbodyParticles(LevelArray);
	fprintf(stderr,"Local Number Of Nbody Particles = %d\n", LocalNumberOfNbodyParticles);
	fprintf(stderr,"Total Number Of Nbody Particles = %d\n", NumberOfNbodyParticles);


	if (NumberOfNbodyParticles == 0) return SUCCESS;

	// Initialization!
	float *NbodyParticleMassTemp;
	float *NbodyParticlePositionTemp[MAX_DIMENSION];
	float *NbodyParticleVelocityTemp[MAX_DIMENSION];
	float *NbodyParticleAccelerationTemp[MAX_DIMENSION][HERMITE_ORDER];
	float *NbodyParticleAccelerationNoStarTemp[MAX_DIMENSION];
	int *NbodyParticleIDTemp;

	fprintf(stderr,"All Good 1\n");

	if (NbodyFirst) {
		InitializeNbodyArrays();
		NbodyFirst = FALSE;
	}

	fprintf(stderr,"All Good 2\n");

	if (LocalNumberOfNbodyParticles == 0)  return SUCCESS;

	NbodyParticleMassTemp = new float[LocalNumberOfNbodyParticles];
	NbodyParticleIDTemp   = new int[LocalNumberOfNbodyParticles];

	for (int dim=0; dim<MAX_DIMENSION; dim++) {

		NbodyParticlePositionTemp[dim] = new float[LocalNumberOfNbodyParticles];
		NbodyParticleVelocityTemp[dim] = new float[LocalNumberOfNbodyParticles];
		NbodyParticleAccelerationNoStarTemp[dim] = new float[LocalNumberOfNbodyParticles];

		for (int i=0; i<HERMITE_ORDER; i++) 
			NbodyParticleAccelerationTemp[dim][i] = new float[NumberOfNbodyParticles];

	}


	fprintf(stderr,"All Good 3\n");

	int* count{ new int (0) };
	// Going through grids, copy nbody particles to temp arrays
	//for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
	for (Temp = LevelArray[true_level]; Temp; Temp = Temp->NextGridThisLevel) {
		fprintf(stderr,"All Good 3-1\n");
		if (Temp->GridData->CopyNbodyParticles(count, NbodyParticleIDTemp, NbodyParticleMassTemp,
					NbodyParticlePositionTemp, NbodyParticleVelocityTemp, NbodyParticleAccelerationNoStarTemp,
					NbodyParticleAccelerationTemp) == FAIL) {
			ENZO_FAIL("Error in grid::CopyNbodyParticles."); 
		}
	}// ENDFOR grid
	//} // ENDFOR level 


	fprintf(stderr,"All Good 4\n");
	//if (LocalNumberOfNbodyParticles != *count) {
	//		fprintf(stderr, "LNNP = %d, cout = %d\n", LocalNumberOfNbodyParticles, *count);
	//	}

	if (true_level == MaximumRefinementLevel) {
		//  Back to Nbody Array
		start_index = FindStartIndex(LocalNumberOfNbodyParticles);
		std::copy(NbodyParticleMassTemp,NbodyParticleMassTemp+LocalNumberOfNbodyParticles,NbodyParticleMass+start_index);
		std::copy(NbodyParticleIDTemp,NbodyParticleIDTemp+LocalNumberOfNbodyParticles,NbodyParticleID+start_index);

		// This can be controlled by address
		for (int dim = 0; dim<MAX_DIMENSION; dim++) {

			std::copy(NbodyParticlePositionTemp[dim],NbodyParticlePositionTemp[dim]+LocalNumberOfNbodyParticles,
					NbodyParticlePosition[dim]+start_index);
			std::copy(NbodyParticleVelocityTemp[dim],NbodyParticleVelocityTemp[dim]+LocalNumberOfNbodyParticles,
					NbodyParticleVelocity[dim]+start_index);
			std::copy(NbodyParticleAccelerationNoStarTemp[dim],NbodyParticleAccelerationNoStarTemp[dim]+LocalNumberOfNbodyParticles,
					NbodyParticleAccelerationNoStar[dim]+start_index);

			for (int order; order<HERMITE_ORDER; order++) {
				std::copy(NbodyParticleAccelerationTemp[dim][order],NbodyParticleAccelerationTemp[dim][order]+LocalNumberOfNbodyParticles,
						NbodyParticleAcceleration[dim][order]+start_index);
			} // ENDFOR order

		} // ENFOR dim
	}

	/* Synchronize number of stars across processors */
#ifdef What
	CommunicationAllReduceValues(NumberOfStarsInGrids, NumberOfGrids, MPI_MAX);
#endif

	delete count;


	return SUCCESS;

}
#endif

