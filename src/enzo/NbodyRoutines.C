/***********************************************************************
	/
	/  FIND ALL STAR PARTICLES OVER ALL PROCESSORS
	/
	/  written by: Yongseok Jo
	/  date:       November, 2022
	/
 ************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm> 
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

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		HierarchyEntry **Grids[]);


#ifdef NBODY
#ifdef USE_MPI
void scan(int *in, int *inout, int *len, MPI_Datatype *dptr);
#endif

void InitializeNbodyArrays(bool NbodyFirst);
void InitializeNbodyArrays(void);
void InitializeNbodyArrays(int); 



void InitializeNbodyArrays(bool NbodyFirst) {

	if (NbodyParticleMass != NULL)
		delete [] NbodyParticleMass;
	NbodyParticleMass = new float[NumberOfNbodyParticles];

	if (NbodyParticleID != NULL)
		delete [] NbodyParticleID;
	NbodyParticleID = new int[NumberOfNbodyParticles];


	for (int dim=0; dim<MAX_DIMENSION; dim++) {

		if (NbodyParticlePosition[dim] != NULL) delete [] NbodyParticlePosition[dim];
		NbodyParticlePosition[dim] = new float[NumberOfNbodyParticles];

		if (NbodyParticleVelocity[dim] != NULL)
			delete [] NbodyParticleVelocity[dim];
		NbodyParticleVelocity[dim] = new float[NumberOfNbodyParticles];
	}

	for (int dim=0; dim<MAX_DIMENSION; dim++) {

		if (NbodyParticleAccelerationNoStar[dim] != NULL)
			delete [] NbodyParticleAccelerationNoStar[dim];
		NbodyParticleAccelerationNoStar[dim] = new float[NumberOfNbodyParticles];

	}
}


void InitializeNbodyArrays(void) {

	if (NbodyParticleMass != NULL)
		delete [] NbodyParticleMass;
	NbodyParticleMass = new float[NumberOfNbodyParticles];

	if (NbodyParticleID != NULL)
		delete [] NbodyParticleID;
	NbodyParticleID = new int[NumberOfNbodyParticles];


	for (int dim=0; dim<MAX_DIMENSION; dim++) {

		if (NbodyParticleVelocity[dim] != NULL)
			delete [] NbodyParticleVelocity[dim];
		NbodyParticleVelocity[dim] = new float[NumberOfNbodyParticles];
	}

	for (int dim=0; dim<MAX_DIMENSION; dim++) { 

		if (NbodyParticleAccelerationNoStar[dim] != NULL)
			delete [] NbodyParticleAccelerationNoStar[dim];
		NbodyParticleAccelerationNoStar[dim] = new float[NumberOfNbodyParticles];

	}
}

void InitializeNbodyArrays(int) {


	if (NbodyParticleID != NULL)
		delete [] NbodyParticleID;
	NbodyParticleID = new int[NumberOfNbodyParticles];


	for (int dim=0; dim<MAX_DIMENSION; dim++) {

		if (NbodyParticlePosition[dim] != NULL) delete [] NbodyParticlePosition[dim];
		NbodyParticlePosition[dim] = new float[NumberOfNbodyParticles];

		if (NbodyParticleVelocity[dim] != NULL)
			delete [] NbodyParticleVelocity[dim];
		NbodyParticleVelocity[dim] = new float[NumberOfNbodyParticles];

	}
}



void DeleteNbodyArrays(void) {

	if (NbodyParticleMass != NULL) {
		delete [] NbodyParticleMass;
	NbodyParticleMass = NULL;
	}
	if (NbodyParticleID != NULL) {
		delete [] NbodyParticleID;
	NbodyParticleID = NULL;
	}

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		if (NbodyParticlePosition[dim] != NULL) {
			delete [] NbodyParticlePosition[dim];
			NbodyParticlePosition[dim] = NULL;
		}

		if (NbodyParticleVelocity[dim] != NULL) {
			delete [] NbodyParticleVelocity[dim];
			NbodyParticleVelocity[dim] = NULL;
		}
	}

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		if (NbodyParticleAccelerationNoStar[dim] != NULL) {
			delete [] NbodyParticleAccelerationNoStar[dim];
			NbodyParticleAccelerationNoStar[dim] = NULL;
		}
	}

		/*
		for (int i=0; i<HERMITE_ORDER; i++) {
			delete [] NbodyParticleAcceleration[dim][i];
			NbodyParticleAcceleration[dim][i] = NULL;
		}
		*/
}


void CopyNbodyArrayToOld(void) {

	/* Create New Array */
	NbodyParticleIDOld =  new int[NumberOfNbodyParticles]{0};
	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		for (int order=0; order<HERMITE_ORDER; order++) {
			NbodyParticleAccelerationOld[dim][order] = new float[NumberOfNbodyParticles]{0};
		}
	}

	/* Assgin */
	NumberOfNbodyParticlesOld = NumberOfNbodyParticles;
	for (int i=0; i<NumberOfNbodyParticles; i++) {
		NbodyParticleIDOld[i] = NbodyParticleID[i];
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			for (int order=0; order<HERMITE_ORDER; order++) {
				NbodyParticleAccelerationOld[dim][order][i] = NbodyParticleAcceleration[dim][order][i];
			}
		}
	}
}


void MatchAccelerationWithIndex(void) {

	/* Match!*/
	for (int i=0; i<NumberOfNbodyParticlesOld; i++) {
		for (int j=0; j<NumberOfNbodyParticlesOld; j++) {
			if (NbodyParticleIDOld[i] == NbodyParticleID[j]) {
				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					for (int order=0; order<HERMITE_ORDER; order++){
						NbodyParticleAcceleration[dim][order][j] = NbodyParticleAccelerationOld[dim][order][i];
					} //ENDFOR order
				} //ENDFOR dim
			} // ENDIF they match
		} // ENFOR old particles
	} //ENDFOR particles

	/* Delete Old Array */
	delete [] NbodyParticleIDOld;
	NbodyParticleIDOld = NULL;
	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		for (int order=0; order<HERMITE_ORDER; order++) {
			delete [] NbodyParticleAccelerationOld[dim][order];
			NbodyParticleAccelerationOld[dim][order] = NULL;
		}
	}
} // MatchAccelerationWithIndex




void FindTotalNumberOfNbodyParticles(LevelHierarchyEntry *LevelArray[], int *LocalNumberOfNbodyParticles) {

	*LocalNumberOfNbodyParticles=0;
	int level;
  LevelHierarchyEntry *Temp;
	NumberOfNbodyParticles = 0;

	//fprintf(stderr,"In the FindTotalNbody\n");

	for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
		for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
			Temp->GridData->SetNumberOfNbodyParticles();
			*LocalNumberOfNbodyParticles += Temp->GridData->ReturnNumberOfNbodyParticles();
		}
	}

#ifdef USE_MPI
	MPI_Allreduce(LocalNumberOfNbodyParticles, &NumberOfNbodyParticles, 1,
			IntDataType, MPI_SUM, enzo_comm);
#else
	NumberOfNbodyParticles = *LocalNumberOfNbodyParticles;
#endif

}



void FindTotalNumberOfNbodyParticles(LevelHierarchyEntry *LevelArray[],
		int *LocalNumberOfNbodyParticles, int *NewLocalNumberOfNbodyParticles) {

	*LocalNumberOfNbodyParticles=0;
	*NewLocalNumberOfNbodyParticles=0;
	int level;
  LevelHierarchyEntry *Temp;
	NumberOfNbodyParticles = 0;

	//fprintf(stderr,"In the FindTotalNbody\n");

	for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
		for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
			Temp->GridData->SetNumberOfNbodyParticles();
			*LocalNumberOfNbodyParticles += Temp->GridData->ReturnNumberOfNbodyParticles();
			*NewLocalNumberOfNbodyParticles += Temp->GridData->ReturnNumberOfNewNbodyParticles();
		}
	}

#ifdef USE_MPI
	//MPI_Allgather(&LocalNumberOfNbodyParticles, 1, MPI_INT, &NumberOfNbodyParticles, 1, MPI_INT, enzo_comm);
	MPI_Allreduce(LocalNumberOfNbodyParticles, &NumberOfNbodyParticles, 1,
			IntDataType, MPI_SUM, enzo_comm);
	MPI_Allreduce(NewLocalNumberOfNbodyParticles, &NumberOfNewNbodyParticles, 1,
			IntDataType, MPI_SUM, enzo_comm);
#else
	NumberOfNbodyParticles = *LocalNumberOfNbodyParticles;
	NumberOfNewNbodyParticles = *NewLocalNumberOfNbodyParticles;
#endif

}



int FindStartIndex(int* LocalNumberOfNbodyParticles) {

	int start_index = 0;
#ifdef USE_MPI
	//MPI_Op  myOp;
	//MPI_Op_create((MPI_User_function *)scan, 1, &myOp);
	//MPI_Scan(LocalNumberOfNbodyParticles, &start_index, 1, MPI_INT, myOp, enzo_comm);

	MPI_Scan(LocalNumberOfNbodyParticles, &start_index, 1, MPI_INT, MPI_SUM, enzo_comm);
	start_index -= *LocalNumberOfNbodyParticles;

	fprintf(stderr, "Proc: %d, LocalNumberOfNbodyParticles:%d\n",MyProcessorNumber, *LocalNumberOfNbodyParticles);
	fprintf(stderr, "Proc: %d, Start Index:%d\n",MyProcessorNumber, start_index);

	//MPI_Op_free(&myOp);
#endif
	return start_index; //return start index of Nbody arrays for each processor
}




#ifdef nouse
void scan(int *in, int *inout, int *len, MPI_Datatype *dptr)
{
	fprintf(stderr,"in=%d, inout=%d, len=%d", in[0], inout[0], *len);

	for (int i = 0; i < *len-1; ++i) {
		inout[i] += in[i];
	}
}
#endif



#endif


