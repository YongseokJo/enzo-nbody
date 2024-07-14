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
#include "CosmologyParameters.h"

void InitializeNbodyArrays(bool NbodyFirst);
void InitializeNbodyArrays(void);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		HierarchyEntry **Grids[]);
int GetUnits(float *DensityUnits, float *LengthUnits,
		float *TemperatureUnits, float *TimeUnits,
		float *VelocityUnits, double *MassUnits, FLOAT Time);
int SendToNbodyFirst(LevelHierarchyEntry *LevelArray[], int level);
int SendToNbody(LevelHierarchyEntry *LevelArray[], int level);

#ifdef NBODY
int PrepareNbodyComputation(LevelHierarchyEntry *LevelArray[], int level)
{

	LevelHierarchyEntry *Temp;

	/* Store AccelerationNoStar of particles to sustainable ParticleAttribute */
	for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
		Temp->GridData->CopyAccelerationToAttribute();



	if (LevelArray[level+1] == NULL) {
		if (NbodyFirst) {
			if (SendToNbodyFirst(LevelArray, level))
				NbodyFirst=FALSE;
		} else {
			SendToNbody(LevelArray, level);	
		}
	} // ENDIF level
	return SUCCESS;
}








int SendToNbodyFirst(LevelHierarchyEntry *LevelArray[], int level) {

	int i, GridNum, LocalNumberOfNbodyParticles=0, NewLocalNumberOfNbodyParticles=0;
	LevelHierarchyEntry *Temp;
	int start_index;

	NumberOfNewNbodyParticles = 0;
	fprintf(stdout, "ENZO: Entering SendToNbodyFirst ...\n");
	//FindTotalNumberOfNbodyParticles(LevelArray, &LocalNumberOfNbodyParticles);
	FindTotalNumberOfNbodyParticles(LevelArray, &LocalNumberOfNbodyParticles, &NewLocalNumberOfNbodyParticles);
	fprintf(stderr, "ENZO: Entering SendToNbodyFirst ...\n");
	//fprintf(stderr,"NewNumberOfParticles=%d\n",NumberOfNewNbodyParticles);
	//fprintf(stderr,"NumberOfParticles=%d\n",NumberOfNbodyParticles);


	//if (NumberOfNbodyParticles == 0)
		//return SUCCESS;
	/* Find the index of the array */
	start_index = FindStartIndex(&LocalNumberOfNbodyParticles);

	fprintf(stdout, "ENZO: Finding star index finished, LocalNumberOfNbodyParticles=%d.\n", LocalNumberOfNbodyParticles);

	/* Do direct calculation!*/
	float dt = 1e-3, scale_factor=1.0;
	float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
				TemperatureUnits=1;
	double MassUnits=1;
	float Time, TimeStep;
	Time = LevelArray[level]->GridData->ReturnTime(); // Not sure for cosmology?
	TimeStep = LevelArray[level]->GridData->ReturnTimeStep(); // Not sure for cosmology?
	if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
				&TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
		ENZO_FAIL("Error in GetUnits.");
	}


	/* At first, all the Nbody info should be updated to fortran codes only once*/
	NbodyParticleIDOld =  NULL;
	NumberOfNbodyParticlesOld = 0;
	float *NbodyParticleMassTemp;
	float *NbodyParticleCreationTimeTemp;
	float *NbodyParticleDynamicalTimeTemp;
	float *NbodyParticlePositionTemp[MAX_DIMENSION];
	float *NbodyParticleVelocityTemp[MAX_DIMENSION];
	float *NbodyParticleAccelerationNoStarTemp[MAX_DIMENSION];

	NbodyParticleIDTemp            = new int[LocalNumberOfNbodyParticles];
	NbodyParticleMassTemp          = new float[LocalNumberOfNbodyParticles];
	NbodyParticleCreationTimeTemp  = new float[LocalNumberOfNbodyParticles];
	NbodyParticleDynamicalTimeTemp = new float[LocalNumberOfNbodyParticles];

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		NbodyParticlePositionTemp[dim]            = new float[LocalNumberOfNbodyParticles];
		NbodyParticleVelocityTemp[dim]            = new float[LocalNumberOfNbodyParticles];
		NbodyParticleAccelerationNoStarTemp[dim]  = new float[LocalNumberOfNbodyParticles];
	}

	fprintf(stdout, "ENZO: Allocate variables.\n");

	/* Get particle information from Grids */
	int count = 0;
	if (LocalNumberOfNbodyParticles > 0) {
		for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++) {
			for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel) {
				if (Temp->GridData->CopyNbodyParticlesFirst(&count, NbodyParticleIDTemp, NbodyParticleMassTemp,
							NbodyParticlePositionTemp, NbodyParticleVelocityTemp, NbodyParticleAccelerationNoStarTemp,
							NbodyParticleCreationTimeTemp, NbodyParticleDynamicalTimeTemp) == FAIL) {
					ENZO_FAIL("Error in grid::CopyNbodyParticlesFirst.");
				}
			}
		}
		if (count != LocalNumberOfNbodyParticles) {
			ENZO_FAIL("Error in grid::CopyNbodyParticlesFirst.");
		}
	}

	fprintf(stdout, "ENZO: 1-1\n");

#ifdef USE_MPI
	if (MyProcessorNumber == ROOT_PROCESSOR) {
		/* Receiving Index, NumberOfParticles, NbodyArrays from other processs */
		int* start_index_all;
		int* LocalNumberAll;
		start_index_all = new int[NumberOfProcessors];
		LocalNumberAll = new int[NumberOfProcessors];
		MPI_Request request;
		MPI_Status status;
		int ierr;
		int errclass,resultlen;
		char err_buffer[MPI_MAX_ERROR_STRING];
		float *NbodyParticleCreationTime;
		float *NbodyParticleDynamicalTime;

		NbodyParticleCreationTime  = new float[NumberOfNbodyParticles];
		NbodyParticleDynamicalTime = new float[NumberOfNbodyParticles];


		MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, LocalNumberAll, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
		MPI_Gather(&start_index, 1, IntDataType, start_index_all, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);

		fprintf(stdout, "ENZO: 1-2\n");
		/* Initialize the nbody array used for direct Nbody calculation*/
		InitializeNbodyArrays(NbodyFirst);

		fprintf(stdout, "ENZO: 1-3\n");


		/*-----------------------------------------------*/
		/******  Gather Arrays from other processes  *****/
		/*-----------------------------------------------*/
		MPI_Igatherv(NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
				NbodyParticleMass, LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm, &request);
		MPI_Wait(&request, &status);
		MPI_Igatherv(NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType,
				NbodyParticleID, LocalNumberAll, start_index_all, IntDataType, ROOT_PROCESSOR, enzo_comm,&request);
		MPI_Wait(&request, &status);
		MPI_Igatherv(NbodyParticleCreationTimeTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
				NbodyParticleCreationTime, LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
		MPI_Wait(&request, &status);
		MPI_Igatherv(NbodyParticleDynamicalTimeTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
				NbodyParticleDynamicalTime, LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
		MPI_Wait(&request, &status);

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			MPI_Igatherv(NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
					NbodyParticlePosition[dim], LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);

			MPI_Igatherv(NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
					NbodyParticleVelocity[dim], LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
		}

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			MPI_Igatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
					NbodyParticleAccelerationNoStar[dim], LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);

		}
		//fprintf(stderr,"Done?1-2\n");


		fprintf(stdout, "First TimeStep: %e, %f\n", TimeStep, TimeUnits);



		/*----------------------------------------------------------*/
		/******** Send Arrays to Fortran Nbody6++  First Time  *****/
		/*--------------------------------------------------------*/
		fprintf(stdout, "ENZO: Waiting for NBODY+ to send data (first) \n");
		fprintf(stderr, "ENZO: Waiting for NBODY+ to send data (first) \n");
		//fprintf(stderr,"NewNumberOfParticles=%d\n",NumberOfNewNbodyParticles);
		CommunicationInterBarrier();
		MPI_Send(&NumberOfNbodyParticles, 1, MPI_INT, 1, 100, inter_comm);
		if (NumberOfNbodyParticles != 0) {
			MPI_Send(NbodyParticleID           , NumberOfNbodyParticles, MPI_INT   , 1, 200, inter_comm);
			MPI_Send(NbodyParticleMass         , NumberOfNbodyParticles, MPI_DOUBLE, 1, 201, inter_comm);
			MPI_Send(NbodyParticleCreationTime , NumberOfNbodyParticles, MPI_DOUBLE, 1, 202, inter_comm);
			MPI_Send(NbodyParticleDynamicalTime, NumberOfNbodyParticles, MPI_DOUBLE, 1, 203, inter_comm);

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Send(NbodyParticlePosition[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 300, inter_comm);
				MPI_Send(NbodyParticleVelocity[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 400, inter_comm);
			}
			// the fourth component of acceleration carries potential
			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Send(NbodyParticleAccelerationNoStar[dim], NumberOfNbodyParticles, MPI_DOUBLE,
					 	1, 500, inter_comm);
			}
		}
		fprintf(stderr, "ENZO: data sending! \n");
		MPI_Send(&TimeStep,                    1, MPI_DOUBLE, 1,  600, inter_comm);
		MPI_Send(&TimeUnits,                   1, MPI_DOUBLE, 1,  700, inter_comm);
		MPI_Send(&LengthUnits,                 1, MPI_DOUBLE, 1,  800, inter_comm);
		MPI_Send(&DensityUnits,                1, MPI_DOUBLE, 1,  900, inter_comm);
		MPI_Send(&VelocityUnits,               1, MPI_DOUBLE, 1, 1000, inter_comm);
		MPI_Send(&StarParticleFeedback       , 1, MPI_INT   , 1, 1100, inter_comm);
		MPI_Send(&StarMassEjectionFraction   , 1, MPI_DOUBLE, 1, 1200, inter_comm);
		MPI_Send(&Time									     , 1, MPI_DOUBLE, 1, 1300, inter_comm);
		MPI_Send(&NbodySmoothingLength       , 1, MPI_DOUBLE, 1, 1400, inter_comm);
		MPI_Send(&NbodyTimeStepConstant	     , 1, MPI_DOUBLE, 1, 1500, inter_comm);
		MPI_Send(&NbodyNeighborRadius		     , 1, MPI_DOUBLE, 1, 1600, inter_comm);
		MPI_Send(&NbodyClusterPosition[3]    , 1, MPI_DOUBLE, 1, 1700, inter_comm);
		MPI_Send(&isNbodyParticleIdentification, 1, MPI_INT   , 1, 1750, inter_comm);
		MPI_Send(&isIdentificationOnTheFly     , 1, MPI_INT   , 1, 1775, inter_comm);
		MPI_Send(&NbodyFixNumNeighbor        , 1, MPI_INT   , 1, 1800, inter_comm);
		MPI_Send(&NbodyBinaryRegularization  , 1, MPI_INT   , 1, 1900, inter_comm);
		MPI_Send(&NbodyBinaryDistance        , 1, MPI_DOUBLE, 1, 2000, inter_comm);
		MPI_Send(&NbodyBinaryTimeStep        , 1, MPI_DOUBLE, 1, 2100, inter_comm);
		//MPI_Send(&HydroMethod         , 1, MPI_INT   , 1, 1200, inter_comm);

		fprintf(stderr, "ENZO: ComovingCoordinates=%d\n",ComovingCoordinates);
		MPI_Send(&ComovingCoordinates        , 1, MPI_INT   , 1, 3000, inter_comm);
		if (ComovingCoordinates) {
			fprintf(stderr, "ENZO: cosmo data sending! \n");
			MPI_Send(&HubbleConstantNow        , 1, MPI_DOUBLE, 1, 3010, inter_comm);
			MPI_Send(&OmegaMatterNow           , 1, MPI_DOUBLE, 1, 3020, inter_comm);
			MPI_Send(&OmegaDarkMatterNow       , 1, MPI_DOUBLE, 1, 3030, inter_comm);
			MPI_Send(&OmegaLambdaNow           , 1, MPI_DOUBLE, 1, 3040, inter_comm);
			MPI_Send(&OmegaRadiationNow        , 1, MPI_DOUBLE, 1, 3050, inter_comm);
			MPI_Send(&ComovingBoxSize          , 1, MPI_DOUBLE, 1, 3060, inter_comm);
			MPI_Send(&MaxExpansionRate         , 1, MPI_DOUBLE, 1, 3070, inter_comm);
			MPI_Send(&InitialTimeInCodeUnits   , 1, MPI_DOUBLE, 1, 3080, inter_comm);
			MPI_Send(&InitialRedshift          , 1, MPI_DOUBLE, 1, 3090, inter_comm);
			MPI_Send(&FinalRedshift            , 1, MPI_DOUBLE, 1, 3100, inter_comm);
			MPI_Send(&CosmologyTableNumberOfBins, 1, MPI_INT  , 1, 3110, inter_comm);
			MPI_Send(&CosmologyTableLogtIndex   , 1, MPI_INT  , 1, 3120, inter_comm);
			MPI_Send(&CosmologyTableLogaInitial , 1, MPI_DOUBLE, 1, 3130, inter_comm);
			MPI_Send(&CosmologyTableLogaFinal   , 1, MPI_DOUBLE, 1, 3140, inter_comm);
			MPI_Send(CosmologyTableLoga, CosmologyTableNumberOfBins, MPI_DOUBLE, 1, 3150, inter_comm);
			MPI_Send(CosmologyTableLogt, CosmologyTableNumberOfBins, MPI_DOUBLE, 1, 3160, inter_comm);
		}
		CommunicationInterBarrier();
		fprintf(stderr, "ENZO: data sent! \n");

		if (start_index_all != NULL)
			delete [] start_index_all;
		start_index_all = NULL;
		if (LocalNumberAll != NULL)
			delete [] LocalNumberAll;
		LocalNumberAll = NULL;
		DeleteNbodyArrays();
		if (NbodyParticleCreationTime != NULL)
			delete [] NbodyParticleCreationTime;
		NbodyParticleCreationTime = NULL;
		if (NbodyParticleDynamicalTime != NULL)
			delete [] NbodyParticleDynamicalTime;
		NbodyParticleDynamicalTime = NULL;
	} // endif : root processor
	else {
		/* Sending Index, NumberOfParticles, NbodyArrays to the root processs */
		MPI_Request request;
		MPI_Status status;

		MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
		MPI_Gather(&start_index, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);


		MPI_Igatherv(NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
				NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm, &request);
		MPI_Wait(&request, &status);
		MPI_Igatherv(NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType,
				NULL, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm, &request);
		MPI_Wait(&request, &status);
		MPI_Igatherv(NbodyParticleCreationTimeTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
				NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
		MPI_Wait(&request, &status);
		MPI_Igatherv(NbodyParticleDynamicalTimeTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
				NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
		MPI_Wait(&request, &status);

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			MPI_Igatherv(NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
					NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);

			MPI_Igatherv(NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
					NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
		}

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			MPI_Igatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
					NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);

		}
	} // endelse: other processor
#endif

	/* Destruct Arrays*/
	if (NbodyParticleMassTemp != NULL)
		delete [] NbodyParticleMassTemp;
	NbodyParticleMassTemp = NULL;

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		if (NbodyParticlePositionTemp[dim] != NULL)
			delete [] NbodyParticlePositionTemp[dim];
		NbodyParticlePositionTemp[dim] = NULL;

		if (NbodyParticleVelocityTemp[dim] != NULL)
			delete [] NbodyParticleVelocityTemp[dim];
		NbodyParticleVelocityTemp[dim] = NULL;
	}

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		if (NbodyParticleAccelerationNoStarTemp[dim] != NULL)
			delete [] NbodyParticleAccelerationNoStarTemp[dim];
		NbodyParticleAccelerationNoStarTemp[dim] = NULL;
	}

	if (NbodyParticleCreationTimeTemp != NULL)
		delete [] NbodyParticleCreationTimeTemp;
	NbodyParticleCreationTimeTemp = NULL;

	if (NbodyParticleDynamicalTimeTemp != NULL)
		delete [] NbodyParticleDynamicalTimeTemp;
	NbodyParticleDynamicalTimeTemp = NULL;

	return SUCCESS;
}










int SendToNbody(LevelHierarchyEntry *LevelArray[], int level) {

	int i, LocalNumberOfNbodyParticles, NewLocalNumberOfNbodyParticles;
	LevelHierarchyEntry *Temp;
	int start_index, start_index_new;

	FindTotalNumberOfNbodyParticles(LevelArray, &LocalNumberOfNbodyParticles, &NewLocalNumberOfNbodyParticles);

	/* Find the index of the array */
	start_index     = FindStartIndex(&LocalNumberOfNbodyParticles);
	start_index_new = FindStartIndex(&NewLocalNumberOfNbodyParticles);


	/* Do direct calculation!*/
	float dt = 1e-3, scale_factor=1.0;
	float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
				TemperatureUnits=1;
	double MassUnits=1;
	float Time, TimeStep;
	Time = LevelArray[level]->GridData->ReturnTime(); // Not sure ?
	TimeStep = LevelArray[level]->GridData->ReturnTimeStep(); // Not sure ?
	if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
				&TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
		ENZO_FAIL("Error in GetUnits.");
	}


	//int *NbodyParticleIDTemp;
	//int *NewNbodyParticleIDTemp;
	float *NbodyParticleMassTemp;
	float *NewNbodyParticleMassTemp;
	float *NewNbodyParticleCreationTimeTemp;
	float *NewNbodyParticleDynamicalTimeTemp;
	float *NewNbodyParticlePositionTemp[MAX_DIMENSION]; // feedback can affect velocity
	float *NewNbodyParticleVelocityTemp[MAX_DIMENSION]; // feedback can affect velocity
	float *NewNbodyParticleAccelerationNoStarTemp[MAX_DIMENSION];
	float *NbodyParticleAccelerationNoStarTemp[MAX_DIMENSION];

	NbodyParticleIDTemp               = new int[LocalNumberOfNbodyParticles];
	NbodyParticleMassTemp             = new double[LocalNumberOfNbodyParticles];
	NewNbodyParticleIDTemp            = new int[NewLocalNumberOfNbodyParticles];
	NewNbodyParticleMassTemp          = new float[NewLocalNumberOfNbodyParticles];
	NewNbodyParticleCreationTimeTemp  = new float[NewLocalNumberOfNbodyParticles];
	NewNbodyParticleDynamicalTimeTemp = new float[NewLocalNumberOfNbodyParticles];
	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		NewNbodyParticlePositionTemp[dim]            = new float[NewLocalNumberOfNbodyParticles];
		NewNbodyParticleVelocityTemp[dim]            = new float[NewLocalNumberOfNbodyParticles];
		NewNbodyParticleAccelerationNoStarTemp[dim]  = new float[NewLocalNumberOfNbodyParticles];
		NbodyParticleAccelerationNoStarTemp[dim]     = new float[LocalNumberOfNbodyParticles];
	}

	/* Get particle information from Grids */
	int count = 0;
	int count_new = 0;
	if (LocalNumberOfNbodyParticles > 0 || NewLocalNumberOfNbodyParticles > 0) {
		for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++)
			for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel) 
				if (Temp->GridData->CopyNbodyParticles(&count, NbodyParticleIDTemp, NbodyParticleMassTemp, NbodyParticleAccelerationNoStarTemp,
							&count_new, NewNbodyParticleIDTemp, NewNbodyParticleMassTemp,
							NewNbodyParticlePositionTemp, NewNbodyParticleVelocityTemp, NewNbodyParticleAccelerationNoStarTemp,
							NewNbodyParticleCreationTimeTemp, NewNbodyParticleDynamicalTimeTemp
							) == FAIL) {

					//NbodyParticleAccelerationTemp
					ENZO_FAIL("Error in grid::CopyNbodyParticles.");
				}

		if (count != LocalNumberOfNbodyParticles)
			ENZO_FAIL("Error in grid::CopyNbodyParticles.");

		if (count_new != NewLocalNumberOfNbodyParticles)
			ENZO_FAIL("Error in grid::CopyNbodyParticles.");
	}

#ifdef USE_MPI
	if (MyProcessorNumber == ROOT_PROCESSOR) {

	int *NewNbodyParticleID;
	float *NbodyParticleMass;
	float *NewNbodyParticleMass;
	float *NewNbodyParticleCreationTime;
	float *NewNbodyParticleDynamicalTime;
	float *NewNbodyParticleVelocity[MAX_DIMENSION]; // feedback can affect velocity
	float *NewNbodyParticlePosition[MAX_DIMENSION]; // feedback can affect velocity
	float *NewNbodyParticleAccelerationNoStar[MAX_DIMENSION];

	NbodyParticleMass             = new float[NumberOfNbodyParticles];
	NewNbodyParticleID            = new int[NumberOfNewNbodyParticles];
	NewNbodyParticleMass          = new float[NumberOfNewNbodyParticles];
	NewNbodyParticleCreationTime  = new float[NumberOfNewNbodyParticles];
	NewNbodyParticleDynamicalTime = new float[NumberOfNewNbodyParticles];

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		NewNbodyParticlePosition[dim]            = new float[NumberOfNewNbodyParticles];
		NewNbodyParticleVelocity[dim]            = new float[NumberOfNewNbodyParticles];
		NewNbodyParticleAccelerationNoStar[dim]  = new float[NumberOfNewNbodyParticles];
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


		MPI_Gather(&NewLocalNumberOfNbodyParticles, 1, IntDataType, NewLocalNumberAll, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
		MPI_Gather(&start_index_new, 1, IntDataType, start_index_all_new, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);

		/* Initialize the nbody array used for direct Nbody calculation*/
		InitializeNbodyArrays();


		/*-----------------------------------------------*/
		/******  Gather Arrays from other processes  *****/
		/*-----------------------------------------------*/
		MPI_Igatherv(NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType,
				NbodyParticleID, LocalNumberAll, start_index_all, IntDataType, ROOT_PROCESSOR, enzo_comm,&request);
		MPI_Wait(&request, &status);
		MPI_Igatherv(NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
				NbodyParticleMass, LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
		MPI_Wait(&request, &status);
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			MPI_Igatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
					NbodyParticleAccelerationNoStar[dim], LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
		}

		if (NumberOfNewNbodyParticles != 0) {
			MPI_Igatherv(NewNbodyParticleMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
					NewNbodyParticleMass, NewLocalNumberAll, start_index_all_new, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm, &request);
			MPI_Wait(&request, &status);
			MPI_Igatherv(NewNbodyParticleIDTemp, NewLocalNumberOfNbodyParticles, IntDataType,
					NewNbodyParticleID, NewLocalNumberAll, start_index_all_new, IntDataType, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
			MPI_Igatherv(NewNbodyParticleCreationTimeTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
					NewNbodyParticleCreationTime, NewLocalNumberAll, start_index_all_new, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
			MPI_Igatherv(NewNbodyParticleDynamicalTimeTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
					NewNbodyParticleDynamicalTime, NewLocalNumberAll, start_index_all_new, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Igatherv(NewNbodyParticlePositionTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
						NewNbodyParticlePosition[dim], NewLocalNumberAll, start_index_all_new, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
				MPI_Wait(&request, &status);
				MPI_Igatherv(NewNbodyParticleVelocityTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
						NewNbodyParticleVelocity[dim], NewLocalNumberAll, start_index_all_new, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
				MPI_Wait(&request, &status);
			}

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Igatherv(NewNbodyParticleAccelerationNoStarTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
						NewNbodyParticleAccelerationNoStar[dim], NewLocalNumberAll, start_index_all_new, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
				MPI_Wait(&request, &status);
			}
		}

		//fprintf(stderr,"Done?1-2\n");



		fprintf(stdout, "TimeStep: %e, %f\n", TimeStep, TimeUnits);
		//MPI_Ssend(&NumberOfNbodyParticles, 1, MPI_INT, NumberOfProcessors, 10, inter_comm);
		/*-----------------------------------------------*/
		/******** Send Arrays to Fortran Nbody6++    *****/
		/*-----------------------------------------------*/

		fprintf(stdout, "ENZO: Waiting for NBODY+ to send data \n");
		//MPI_Ssend(&NumberOfNbodyParticles, 1, MPI_INT, NumberOfProcessors, 10, inter_comm);
		CommunicationInterBarrier();

		if (NumberOfNbodyParticles != 0) {
		for (int j=0; j<NumberOfNbodyParticles; j++)
			fprintf(stdout, "ENZO: PID in ENZO =%d in PNC\n", NbodyParticleID[j]);
		//MPI_Ssend(&NumberOfNbodyParticles, 1, MPI_INT, NumberOfProcessors, 10, inter_comm);
			MPI_Send(NbodyParticleID  , NumberOfNbodyParticles, MPI_INT   , 1, 25, inter_comm);
			MPI_Send(NbodyParticleMass, NumberOfNbodyParticles, MPI_DOUBLE, 1, 40, inter_comm);
			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Send(NbodyParticleAccelerationNoStar[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 50, inter_comm);
			}
		}


		//fprintf(stderr, "NumberOfNewNbodyParticles=%d in PNC\n", NumberOfNewNbodyParticles);
		MPI_Send(&NumberOfNewNbodyParticles, 1, MPI_INT, 1, 100, inter_comm);
		if (NumberOfNewNbodyParticles != 0) {
			/*
			for (int i = 0; i < NumberOfNewNbodyParticles; i++) {
				fprintf(stderr, "Mass Of NewNbodyParticles=%.3e\n", NewNbodyParticleMass[i]);
				fprintf(stderr, "Vel  Of NewNbodyParticles=(%.3e, %.3e, %.3e)\n", 
						NewNbodyParticleVelocity[0][i], NewNbodyParticleVelocity[1][i], NewNbodyParticleVelocity[2][i]);
				fprintf(stderr, "Pos  Of NewNbodyParticles=(%.3e, %.3e, %.3e)\n",
						NewNbodyParticlePosition[0][i], NewNbodyParticlePosition[1][i], NewNbodyParticlePosition[2][i]);
				fprintf(stderr, "Acc  Of NewNbodyParticles=(%.3e, %.3e, %.3e)\n",
						NewNbodyParticleAccelerationNoStar[0][i], NewNbodyParticleAccelerationNoStar[1][i], NewNbodyParticleAccelerationNoStar[2][i]);
			}
			*/
			MPI_Send(NewNbodyParticleID           , NumberOfNewNbodyParticles, MPI_INT   , 1, 200, inter_comm);
			MPI_Send(NewNbodyParticleMass         , NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 201, inter_comm);
			MPI_Send(NewNbodyParticleCreationTime , NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 202, inter_comm);
			MPI_Send(NewNbodyParticleDynamicalTime, NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 203, inter_comm);

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Send(NewNbodyParticlePosition[dim]          , NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 300, inter_comm);
				MPI_Send(NewNbodyParticleVelocity[dim]          , NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 400, inter_comm);
				MPI_Send(NewNbodyParticleAccelerationNoStar[dim], NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 500, inter_comm);
			}
		}

		ierr = MPI_Send(&TimeStep, 1, MPI_DOUBLE, 1, 600, inter_comm);
		ierr = MPI_Send(&Time    , 1, MPI_DOUBLE, 1, 700, inter_comm);
		//ierr = MPI_Send(&TimeUnits, 1, MPI_DOUBLE, 1, 700, inter_comm);
		ierr = CommunicationInterBarrier();
		fprintf(stdout, "ENZO: Data sent.\n");


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
		DeleteNbodyArrays();


		fprintf(stdout, "ENZO: 1\n");

//#Merge  part



		if (NbodyParticleMass != NULL)
			delete [] NbodyParticleMass;
		NbodyParticleMass = NULL;

		if (NewNbodyParticleMass != NULL)
			delete [] NewNbodyParticleMass;
		NewNbodyParticleMass = NULL;

		if (NewNbodyParticleID != NULL)
			delete [] NewNbodyParticleID;
		NewNbodyParticleID = NULL;

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			if (NewNbodyParticlePosition[dim] != NULL)
				delete [] NewNbodyParticlePosition[dim];
			NewNbodyParticlePosition[dim] = NULL;

			if (NewNbodyParticleVelocity[dim] != NULL)
				delete [] NewNbodyParticleVelocity[dim];
			NewNbodyParticleVelocity[dim] = NULL;
		}

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			if (NewNbodyParticleAccelerationNoStar[dim] != NULL)
				delete [] NewNbodyParticleAccelerationNoStar[dim];
			NewNbodyParticleAccelerationNoStar[dim] = NULL;
		}

		if (NewNbodyParticleCreationTime != NULL)
			delete [] NewNbodyParticleCreationTime;
		NewNbodyParticleCreationTime = NULL;

		if (NewNbodyParticleDynamicalTime != NULL)
			delete [] NewNbodyParticleDynamicalTime;
		NewNbodyParticleDynamicalTime = NULL;

		fprintf(stdout, "ENZO: 2\n");

	} // endif : root processor
	else {
		/* Sending Index, NumberOfParticles, NbodyArrays to the root processs */

		MPI_Request request;
		MPI_Status status;

		MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
		MPI_Gather(&start_index, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);


		MPI_Gather(&NewLocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
		MPI_Gather(&start_index_new, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);



		MPI_Igatherv(NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType,
				NULL, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm, &request);
		MPI_Wait(&request, &status);
		MPI_Igatherv(NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
				NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm, &request);
		MPI_Wait(&request, &status);
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			//MPI_Igatherv(NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
			//			NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			//	MPI_Wait(&request, &status);
			MPI_Igatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
					NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
		}



		if (NumberOfNewNbodyParticles != 0) {
			MPI_Igatherv(NewNbodyParticleMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
					NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm, &request);
			MPI_Wait(&request, &status);
			MPI_Igatherv(NewNbodyParticleIDTemp, NewLocalNumberOfNbodyParticles, IntDataType,
					NULL, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm, &request);
			MPI_Wait(&request, &status);
			MPI_Igatherv(NewNbodyParticleCreationTimeTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
					NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);
			MPI_Igatherv(NewNbodyParticleDynamicalTimeTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
					NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
			MPI_Wait(&request, &status);

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Igatherv(NewNbodyParticlePositionTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
						NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
				MPI_Wait(&request, &status);
				MPI_Igatherv(NewNbodyParticleVelocityTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
						NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
				MPI_Wait(&request, &status);
			}

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				MPI_Igatherv(NewNbodyParticleAccelerationNoStarTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE,
						NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
				MPI_Wait(&request, &status);
			}
		}



	} // endelse: other processor
#endif

	/* Destruct Arrays*/
	/*
		 if (NbodyParticleIDTemp != NULL)
		 delete [] NbodyParticleIDTemp;
		 NbodyParticleIDTemp = NULL;
		 */

	fprintf(stdout, "ENZO: 3\n");

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		if (NbodyParticleAccelerationNoStarTemp[dim] != NULL)
			delete [] NbodyParticleAccelerationNoStarTemp[dim];
		NbodyParticleAccelerationNoStarTemp[dim] = NULL;
	}
	if (NbodyParticleMassTemp != NULL)
		delete [] NbodyParticleMassTemp;
	NbodyParticleMassTemp = NULL;
	if (NewNbodyParticleMassTemp != NULL)
		delete [] NewNbodyParticleMassTemp;
	NewNbodyParticleMassTemp = NULL;
	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		if (NewNbodyParticlePositionTemp[dim] != NULL)
			delete [] NewNbodyParticlePositionTemp[dim];
		NewNbodyParticlePositionTemp[dim] = NULL;
		if (NewNbodyParticleVelocityTemp[dim] != NULL)
			delete [] NewNbodyParticleVelocityTemp[dim];
		NewNbodyParticleVelocityTemp[dim] = NULL;
	}
	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		if (NewNbodyParticleAccelerationNoStarTemp[dim] != NULL)
			delete [] NewNbodyParticleAccelerationNoStarTemp[dim];
		NewNbodyParticleAccelerationNoStarTemp[dim] = NULL;
	}
	if (NewNbodyParticleCreationTimeTemp != NULL)
		delete [] NewNbodyParticleCreationTimeTemp;
	NewNbodyParticleCreationTimeTemp = NULL;
	if (NewNbodyParticleDynamicalTimeTemp != NULL)
		delete [] NewNbodyParticleDynamicalTimeTemp;
	NewNbodyParticleDynamicalTimeTemp = NULL;

	fprintf(stdout, "ENZO: 4\n");
	return SUCCESS;
}

#endif


