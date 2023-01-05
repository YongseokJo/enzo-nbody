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
void scan(int *in, int *inout, int *len, MPI_Datatype *dptr);




void InitializeNbodyArrays(void) {

	if (NbodyParticleMass != NULL)
		delete [] NbodyParticleMass;
	NbodyParticleMass = new float[NumberOfNbodyParticles];

	if (NbodyParticleID != NULL)
		delete [] NbodyParticleID;
	NbodyParticleID   = new int[NumberOfNbodyParticles];

	for (int dim=0; dim<MAX_DIMENSION; dim++) {

		if (NbodyParticlePosition[dim] != NULL)
			delete [] NbodyParticlePosition[dim];
		NbodyParticlePosition[dim] = new float[NumberOfNbodyParticles];

		if (NbodyParticleVelocity[dim] != NULL)
			delete [] NbodyParticleVelocity[dim];
		NbodyParticleVelocity[dim] = new float[NumberOfNbodyParticles];

		for (int i=0; i<HERMITE_ORDER; i++) {

			if (NbodyParticleVelocity[dim] != NULL)
				delete [] NbodyParticleAcceleration[dim][i];
			NbodyParticleAcceleration[dim][i] = new float[NumberOfNbodyParticles];

		}
	}
}





void DeleteNbodyArrays(void) {

	delete [] NbodyParticleMass;
	delete [] NbodyParticleID;
	NbodyParticleMass = NULL;
	NbodyParticleID = NULL;

	for (int dim=0; dim<MAX_DIMENSION; dim++) {

		delete [] NbodyParticlePosition[dim];
		delete [] NbodyParticleVelocity[dim];
		NbodyParticlePosition[dim] = NULL;
		NbodyParticleVelocity[dim] = NULL;

		for (int i=0; i<HERMITE_ORDER; i++) {
			delete [] NbodyParticleAcceleration[dim][i];
			NbodyParticleAcceleration[dim][i] = NULL;
		}
	}
}


/***
void CopyNbodyArrays(float *pos[], float *vel[], float *mass, float *id, float *acc[][]) {

}
***/






int FindTotalNumberOfNbodyParticles(LevelHierarchyEntry *LevelArray[]) {

	int LocalNumberOfNbodyParticles=0;
	int level;
	int num_tmp; // delete
  LevelHierarchyEntry *Temp;
	NumberOfNbodyParticles = 0;

	//fprintf(stderr,"In the FindTotalNbody\n");

	for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
		num_tmp = 0;
		for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
			num_tmp += Temp->GridData->ReturnNumberOfNbodyParticles();
			LocalNumberOfNbodyParticles += Temp->GridData->ReturnNumberOfNbodyParticles();
		}
		fprintf(stderr,"level=%d, LocalNumber=%d",level, num_tmp);
	}

#ifdef USE_MPI
	//MPI_Allgather(&LocalNumberOfNbodyParticles, 1, MPI_INT, &NumberOfNbodyParticles, 1, MPI_INT, MPI_COMM_WORLD);
	MPI_Allreduce(&LocalNumberOfNbodyParticles, &NumberOfNbodyParticles, 1,
			IntDataType, MPI_SUM, MPI_COMM_WORLD);
#endif

	return LocalNumberOfNbodyParticles;
}





int FindStartIndex(int LocalNumberOfNbodyParticles) {

	int *start_index = 0;
#ifdef USE_MPI
	MPI_Op  myOp;
	MPI_Op_create((MPI_User_function *)scan, 0, &myOp);
	MPI_Scan(&LocalNumberOfNbodyParticles, &start_index, 1, MPI_INT, myOp, MPI_COMM_WORLD);

	fprintf(stderr, "Proc: %d, LocalNumberOfNbodyParticles:%d\n",MyProcessorNumber, LocalNumberOfNbodyParticles);
	fprintf(stderr, "Proc: %d, Start Index:%d\n",MyProcessorNumber,start_index);

#endif
	return *start_index; //return start index of Nbody arrays for each processor
}





void scan(int *in, int *inout, int *len, MPI_Datatype *dptr)
{
	int i;

	for (i = 1; i < *len; ++i) {
		*inout = *in + *inout;
		in++;
		inout++;
	}
}



#endif


