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

	*LocalNumberOfNbodyParticles    = 0;
	*NewLocalNumberOfNbodyParticles = 0;
	int level;
  LevelHierarchyEntry *Temp;
	NumberOfNbodyParticles    = 0;
	NumberOfNewNbodyParticles = 0;

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

	//fprintf(stderr, "Proc: %d, LocalNumberOfNbodyParticles:%d\n",MyProcessorNumber, *LocalNumberOfNbodyParticles);
	//fprintf(stderr, "Proc: %d, Start Index:%d\n",MyProcessorNumber, start_index);

	//MPI_Op_free(&myOp);
#endif
	return start_index; //return start index of Nbody arrays for each processor
}


void IdentifyNbodyParticlesEvolveLevel(LevelHierarchyEntry *LevelArray[], int level) {

	if (!isNbodyParticleIdentification)
		return;

  LevelHierarchyEntry *Temp;

	if (NbodyClusterPosition[0] == -1 || NbodyClusterPosition[0] == 0) {
		double TotalMass=0;

		NbodyClusterPosition[0] = 0.;
		NbodyClusterPosition[1] = 0.;
		NbodyClusterPosition[2] = 0.;

		for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++) {
			for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
				Temp->GridData->GetNbodyCenterOfMass(TotalMass);
			}
		}

		CommunicationAllReduceValues(NbodyClusterPosition, 3, MPI_SUM);
		CommunicationAllReduceValues(&TotalMass, 1, MPI_SUM);

		if (TotalMass==0)
			return;

		NbodyClusterPosition[0] /= TotalMass;
		NbodyClusterPosition[1] /= TotalMass;
		NbodyClusterPosition[2] /= TotalMass;

		//fprintf(stderr, "GetCenterOfMass\n");
		//fprintf(stderr, "NbodyClusterPosition = (%.e, %.e, %.e), R2=%.3e\n", NbodyClusterPosition[0], NbodyClusterPosition[1], NbodyClusterPosition[2], NbodyClusterPosition[3]);
	}

	for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++) {
		for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel) {
			Temp->GridData->IdentifyNbodyParticles();
		}
	}
}

void GetCenterOfMass(double *mass, double *x[MAX_DIMENSION], double *v[MAX_DIMENSION], double x_com[], double v_com[], int N) {
	double total_mass=0.;
	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		x_com[dim] = 0.;
		v_com[dim] = 0.;
	}

	for (int i=0; i<N; i++) {
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			x_com[dim] += mass[i]*x[dim][i];
			v_com[dim] += mass[i]*v[dim][i];
		}
		total_mass += mass[i];
	}

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		x_com[dim] /= total_mass;
		v_com[dim] /= total_mass;
	}
}





void grid::GetNbodyCenterOfMass(double &TotalMass) {

	if (MyProcessorNumber != ProcessorNumber) return;


	double Mass;

	float dv = CellWidth[0][0]*CellWidth[0][0]*CellWidth[0][0];

	for (int i=0; i < NumberOfParticles; i++) {
		if (ParticleType[i] == PARTICLE_TYPE_NBODY) {
			Mass = ParticleMass[i]*dv;
			NbodyClusterPosition[0] += ParticlePosition[0][i]*Mass;
			NbodyClusterPosition[1] += ParticlePosition[1][i]*Mass;
			NbodyClusterPosition[2] += ParticlePosition[2][i]*Mass;
			TotalMass += Mass;
		}
	}
}


#ifdef no_use
void scan(int *in, int *inout, int *len, MPI_Datatype *dptr)
{
	fprintf(stderr,"in=%d, inout=%d, len=%d", in[0], inout[0], *len);

	for (int i = 0; i < *len-1; ++i) {
		inout[i] += in[i];
	}
}
#endif



#endif


