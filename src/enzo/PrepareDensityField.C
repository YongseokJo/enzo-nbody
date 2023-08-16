/***********************************************************************
	/
	/  PREPARE DENSITY FIELD (CALLED BY EVOLVE LEVEL)
	/
	/  written by: Greg Bryan
	/  date:       June, 1999
	/  modifiedN:  Robert Harkness
	/  date:       February, 2008
	/
	/ ======================================================================= 
	/ This routine prepares the density field for all the grids on this level,
	/ both particle and baryonic densities.  It also calculates the potential
	/ field if this is level 0 (since this involves communication). 
	/
	/   This is part of a collection of routines called by EvolveLevel.
	/   These have been optimized for enhanced message passing
	/   performance by performing two passes -- one which generates
	/   sends and the second which receives them.
	/
	/  modified: Robert Harkness, December 2007
	/
 ************************************************************************/

#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */

#include <stdio.h>
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "performance.h"
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
#include "communication.h"
#include "CommunicationUtilities.h"
#include "phys_constants.h"
#include "ActiveParticle.h"

/* function prototypes */

int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time, bool NoStar);
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0, bool NoStar=FALSE);

int CommunicationBufferPurge(void);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[],
		int NumberOfSubgrids[],
		int FluxFlag,
		TopGridData* MetaData, bool NoStar);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
		int NumberOfSubgrids[] = NULL,
		int FluxFlag = FALSE,
		TopGridData* MetaData = NULL, bool NoStar=FALSE);

int ActiveParticleDepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
		int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
		int level);

int PrepareGravitatingMassField1(HierarchyEntry *Grid);
#ifdef FAST_SIB
int PrepareGravitatingMassField2a(HierarchyEntry *Grid, int grid1,
		SiblingGridList SiblingList[],
		TopGridData *MetaData, int level,
		FLOAT When);
#else
int PrepareGravitatingMassField2a(HierarchyEntry *Grid, TopGridData *MetaData,
		LevelHierarchyEntry *LevelArray[], int level,
		FLOAT When);
#endif

int PrepareGravitatingMassField2b(HierarchyEntry *Grid, int level);



#ifdef NBODY
int PrepareGravitatingMassFieldNoStar1(HierarchyEntry *Grid);
#ifdef FAST_SIB
int PrepareGravitatingMassFieldNoStar2a(HierarchyEntry *Grid, int grid1,
		SiblingGridList SiblingList[],
		TopGridData *MetaData, int level,
		FLOAT When);
#else
int PrepareGravitatingMassFieldNoStar2a(HierarchyEntry *Grid, TopGridData *MetaData,
		LevelHierarchyEntry *LevelArray[], int level,
		FLOAT When);
#endif

int PrepareGravitatingMassFieldNoStar2b(HierarchyEntry *Grid, int level);


#endif

#ifdef FAST_SIB
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
		SiblingGridList SiblingList[],
		HierarchyEntry *Grids[], int NumberOfGrids);
#else
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
		HierarchyEntry *Grids[], int NumberOfGrids);
#endif

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		HierarchyEntry **Grids[]);



extern int CopyPotentialFieldAverage;

#define GRIDS_PER_LOOP 100000



#ifdef FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
		int level, TopGridData *MetaData, FLOAT When,SiblingGridList **SiblingGridListStorage)
#else   // !FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
		int level, TopGridData *MetaData, FLOAT When)
#endif  // end FAST_SIB
{

	/* Return if this does not concern us */
	if (!SelfGravity) return SUCCESS;

	LCAPERF_START("PrepareDensityField");

	int grid1, grid2, StartGrid, EndGrid;

	/* Set the time for evaluation of the fields, etc. */

	FLOAT EvaluateTime = LevelArray[level]->GridData->ReturnTime() +
		When*LevelArray[level]->GridData->ReturnTimeStep();

	/* If level is above MaximumGravityRefinementLevel, then just
		 update the gravity at the MaximumGravityRefinementLevel. */

	int reallevel = level;
	level = min(level, MaximumGravityRefinementLevel);

	/* Create an array (Grids) of all the grids. */

	typedef HierarchyEntry* HierarchyEntryPointer;
	HierarchyEntry **Grids;
	int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
	SiblingGridList *SiblingList = SiblingGridListStorage[level];

	/************************************************************************/
	/* Grids: Deposit particles in their GravitatingMassFieldParticles.
		 (Do a batch of grids at a time; this is a loop over the batches)
		 */

	if (traceMPI) 
		fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField (Send)\n");

#ifdef FORCE_MSG_PROGRESS 
	CommunicationBarrier();
#endif

	TIME_MSG("Depositing particle mass field");
	LCAPERF_START("DepositParticleMassField");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		/* First, generate the receive calls. */

		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
		CommunicationDirection = COMMUNICATION_POST_RECEIVE;
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			DepositParticleMassField(Grids[grid1], EvaluateTime, FALSE);

#ifdef FORCE_MSG_PROGRESS 
		CommunicationBarrier();
#endif

		if (traceMPI) 
			fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField"
					" (Receive)\n");

		/* Next, send data and process grids on the same processor. */

		fprintf(stdout,"4-1\n"); // by YS
		CommunicationDirection = COMMUNICATION_SEND;
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			DepositParticleMassField(Grids[grid1], EvaluateTime, FALSE);

		/* Finally, receive the data and process it. */

		fprintf(stdout,"4-2\n"); // by YS

	} // ENDFOR grid batches
	LCAPERF_STOP("DepositParticleMassField");
	fprintf(stdout,"\nProc:%d 4-2 ends.\n", MyProcessorNumber); //by YS


#ifdef FORCE_BUFFER_PURGE
	CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
	CommunicationBarrier();
#endif


	////////////// GravitatingMassFieldParticles for No Star
#ifdef NBODY
	TIME_MSG("Depositing particle mass field");
	LCAPERF_START("DepositParticleMassField");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		/* First, generate the receive calls. */

		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
		CommunicationDirection = COMMUNICATION_POST_RECEIVE;
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			DepositParticleMassField(Grids[grid1], EvaluateTime, TRUE);

#ifdef FORCE_MSG_PROGRESS 
		CommunicationBarrier();
#endif

		if (traceMPI) 
			fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField"
					" (Receive)\n");

		/* Next, send data and process grids on the same processor. */

		fprintf(stdout,"\nProc:%d 4-10\n", MyProcessorNumber); // by YS
		CommunicationDirection = COMMUNICATION_SEND;
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			DepositParticleMassField(Grids[grid1], EvaluateTime, TRUE);

		/* Finally, receive the data and process it. */

		fprintf(stdout,"\nProc:%d 4-20\n", MyProcessorNumber); // by YS
		CommunicationReceiveHandler(NULL,NULL,FALSE,NULL,TRUE);

	} // ENDFOR grid batches
	LCAPERF_STOP("DepositParticleMassField");


#endif
	////////////// GravitatingMassFieldParticles for No Star Ends
	///
	///
	///

// by YS debug
#ifdef NBODY 

	int ndiff=0;
	float diff=0;
	float org,nostar;
	int size=1;
	for (StartGrid = 0; StartGrid< NumberOfGrids; StartGrid+=GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		ndiff=0;
		diff=0;
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
			for (int dim=0; dim<Grids[grid1]->GridData->GetGridRank(); dim++) {
				size *= Grids[grid1]->GridData->ReturnGravitatingMassFieldDimension(dim);
			}
			/*
			for (int i=0; i<size; i++) {
				org = Grids[grid1]->GridData->GetGravitatingMassFieldParticles()[0][i];
				nostar = Grids[grid1]->GridData->GetGravitatingMassFieldParticles()[1][i];
				diff = abs(org - nostar);
				if (diff>1e-5) {
					fprintf(stderr, "In EvolveLevel, diff = %f\n", diff);
					ndiff++;
				} 
			}*/
		}
	}
	fprintf(stderr, "In EvolveLevel, ndiff = (%d/%d)\n", ndiff,size);
#endif







	/******************************************************************/
	/* Grids: compute the GravitatingMassField (baryons & particles). */
	/*   This is now split into two section. */

	if (traceMPI) 
		fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF1 (send)\n", 
				MyProcessorNumber);

	TIME_MSG("PrepareGravitatingMassField1");
	LCAPERF_START("PrepareGravitatingMassField1");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		/* ----- section 1 ---- */
		/* First, generate the receive calls. */

		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
		CommunicationDirection = COMMUNICATION_POST_RECEIVE;

		fprintf(stdout,"Proc:%d, 4-30\n",MyProcessorNumber); // by YS
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassField1(Grids[grid1]);

		fprintf(stdout,"Proc:%d, 4-31\n",MyProcessorNumber); // by YS
		/* Next, send data and process grids on the same processor. */

		CommunicationDirection = COMMUNICATION_SEND;
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassField1(Grids[grid1]);

		/* Finally, receive the data and process it. */
		fprintf(stdout,"Proc:%d, 4-32\n", MyProcessorNumber); // by YS

		CommunicationReceiveHandler();

	} // ENDFOR grid batches
	LCAPERF_STOP("PrepareGravitatingMassField1");
	fprintf(stdout,"Proc:%d, 4-4\n", MyProcessorNumber); // by YS


#ifdef FORCE_MSG_PROGRESS 
	CommunicationBarrier();
#endif

	if (traceMPI) 
		fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF2 (receive)\n", 
				MyProcessorNumber);

	TIME_MSG("PrepareGravitatingMassField2");
	LCAPERF_START("PrepareGravitatingMassField2a");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		/* ----- section 2 ---- */
		/* First, generate the receive calls. */

		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
#ifdef BITWISE_IDENTICALITY
		CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
#else
		CommunicationDirection = COMMUNICATION_POST_RECEIVE;
#endif

#ifdef FAST_SIB
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassField2a(Grids[grid1], grid1, SiblingList,
					MetaData, level, When);
#else
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassField2a(Grids[grid1], MetaData, LevelArray,
					level, When);
#endif
		fprintf(stdout,"4-5"); // by YS

#ifndef BITWISE_IDENTICALITY
		/* Next, send data and process grids on the same processor. */

		CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassField2a(Grids[grid1], grid1, SiblingList,
					MetaData, level, When);
#else
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassField2a(Grids[grid1], MetaData, LevelArray,
					level, When);
#endif

		fprintf(stdout,"4-6"); // by YS
		CommunicationReceiveHandler();
#endif /* BITWISE_IDENTICALITY */
		fprintf(stdout,"Proc:%d, 4-7\n", MyProcessorNumber); // by YS

	} // ENDFOR grid batches
	LCAPERF_STOP("PrepareGravitatingMassField2a");

#ifdef FORCE_BUFFER_PURGE
	CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
	CommunicationBarrier();
#endif

	/************************************************************************/
	LCAPERF_START("PrepareGravitatingMassField2b");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		/* ----- section 2 ---- */
		/* First, generate the receive calls. */

		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
		CommunicationDirection = COMMUNICATION_POST_RECEIVE;

		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassField2b(Grids[grid1], level);

		/* Next, send data and process grids on the same processor. */

		CommunicationDirection = COMMUNICATION_SEND;
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassField2b(Grids[grid1], level);

		CommunicationReceiveHandler();

	} // ENDFOR grid batches
	LCAPERF_STOP("PrepareGravitatingMassField2b");

	fprintf(stdout,"Proc:%d, 4-8\n", MyProcessorNumber); // by YS






#ifdef NBODY
	/*******************************************************************************/
	// No Star Starts
	/******************************************************************/
	/* Grids: compute the GravitatingMassField (baryons & particles). */
	/*   This is now split into two section. */

	if (traceMPI) 
		fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF1 (send)\n", 
				MyProcessorNumber);

	TIME_MSG("PrepareGravitatingMassField1");
	LCAPERF_START("PrepareGravitatingMassField1");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		/* ----- section 1 ---- */
		/* First, generate the receive calls. */

		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
		CommunicationDirection = COMMUNICATION_POST_RECEIVE;

		fprintf(stdout,"Proc:%d, 4-30\n",MyProcessorNumber); // by YS
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassFieldNoStar1(Grids[grid1]);

		fprintf(stdout,"Proc:%d, 4-31\n",MyProcessorNumber); // by YS
		/* Next, send data and process grids on the same processor. */

		CommunicationDirection = COMMUNICATION_SEND;
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassFieldNoStar1(Grids[grid1]);

		/* Finally, receive the data and process it. */
		fprintf(stdout,"Proc:%d, 4-32\n", MyProcessorNumber); // by YS

		CommunicationReceiveHandler(NULL,NULL,FALSE,NULL,TRUE);

	} // ENDFOR grid batches
	LCAPERF_STOP("PrepareGravitatingMassField1");
	fprintf(stdout,"Proc:%d, 4-4\n", MyProcessorNumber); // by YS


#ifdef FORCE_MSG_PROGRESS 
	CommunicationBarrier();
#endif

	if (traceMPI) 
		fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF2 (receive)\n", 
				MyProcessorNumber);

	TIME_MSG("PrepareGravitatingMassField2");
	LCAPERF_START("PrepareGravitatingMassField2a");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		/* ----- section 2 ---- */
		/* First, generate the receive calls. */

		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
#ifdef BITWISE_IDENTICALITY
		CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
#else
		CommunicationDirection = COMMUNICATION_POST_RECEIVE;
#endif

#ifdef FAST_SIB
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassFieldNoStar2a(Grids[grid1], grid1, SiblingList,
					MetaData, level, When);
#else
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassFieldNoStar2a(Grids[grid1], MetaData, LevelArray,
					level, When);
#endif
		fprintf(stdout,"4-5"); // by YS

#ifndef BITWISE_IDENTICALITY
		/* Next, send data and process grids on the same processor. */

		CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassFieldNoStar2a(Grids[grid1], grid1, SiblingList,
					MetaData, level, When);
#else
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassFieldNoStar2a(Grids[grid1], MetaData, LevelArray,
					level, When);
#endif

		fprintf(stdout,"4-6"); // by YS
		CommunicationReceiveHandler(NULL,NULL,FALSE,NULL,TRUE);
#endif /* BITWISE_IDENTICALITY */
		fprintf(stdout,"Proc:%d, 4-7\n", MyProcessorNumber); // by YS

	} // ENDFOR grid batches
	LCAPERF_STOP("PrepareGravitatingMassField2a");

#ifdef FORCE_BUFFER_PURGE
	CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
	CommunicationBarrier(); //by YS
#endif

	/************************************************************************/
	LCAPERF_START("PrepareGravitatingMassField2b");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		/* ----- section 2 ---- */
		/* First, generate the receive calls. */

		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
		CommunicationDirection = COMMUNICATION_POST_RECEIVE;

		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassFieldNoStar2b(Grids[grid1], level);

		/* Next, send data and process grids on the same processor. */

		CommunicationDirection = COMMUNICATION_SEND;
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			PrepareGravitatingMassFieldNoStar2b(Grids[grid1], level);

		CommunicationReceiveHandler(NULL,NULL,FALSE,NULL,TRUE);

	} // ENDFOR grid batches
	LCAPERF_STOP("PrepareGravitatingMassField2b");

	fprintf(stdout,"Proc:%d, 4-8\n", MyProcessorNumber); // by YS

	/*******************************************************************************/
	// No Star Done
	/*******************************************************************************/
#endif


	CommunicationBarrier(); //by YS





	/************************************************************************/
	/* Copy overlapping mass fields to ensure consistency and B.C.'s. */

	//  if (level > 0)

	if (traceMPI) 
		fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): COMF1 (send)\n", 
				MyProcessorNumber);

	TIME_MSG("CopyOverlappingMassField");
	LCAPERF_START("CopyOverlappingMassField");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		CommunicationDirection = COMMUNICATION_POST_RECEIVE;
		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;

		fprintf(stdout,"4-8'\n"); // by YS
#ifdef FAST_SIB
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
				Grids[grid1]->GridData->
					CheckForOverlap(SiblingList[grid1].GridList[grid2],

							MetaData->LeftFaceBoundaryCondition,
							MetaData->RightFaceBoundaryCondition,
							&grid::CopyOverlappingMassField);
#else
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
				Grids[grid1]->GridData->
					CheckForOverlap(Grids[grid2]->GridData,
							MetaData->LeftFaceBoundaryCondition,
							MetaData->RightFaceBoundaryCondition,
							&grid::CopyOverlappingMassField);
#endif
		fprintf(stdout,"4-9\n"); // by YS

		CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
				Grids[grid1]->GridData->
					CheckForOverlap(SiblingList[grid1].GridList[grid2],
							MetaData->LeftFaceBoundaryCondition,
							MetaData->RightFaceBoundaryCondition,
							&grid::CopyOverlappingMassField);
#else
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
				Grids[grid1]->GridData->
					CheckForOverlap(Grids[grid2]->GridData,
							MetaData->LeftFaceBoundaryCondition,
							MetaData->RightFaceBoundaryCondition,
							&grid::CopyOverlappingMassField);
#endif

		CommunicationReceiveHandler();

	} // ENDFOR grid batches
	LCAPERF_STOP("CopyOverlappingMassField");
	fprintf(stdout,"4-10\n"); // by YS


#ifdef NBODY
	/************************************************************************/
	// No Star Starts
	/************************************************************************/
	/* Copy overlapping mass fields to ensure consistency and B.C.'s. */

	//  if (level > 0)

	if (traceMPI) 
		fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): COMF1 (send)\n", 
				MyProcessorNumber);

	TIME_MSG("CopyOverlappingMassField");
	LCAPERF_START("CopyOverlappingMassField");
	for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

		CommunicationDirection = COMMUNICATION_POST_RECEIVE;
		CommunicationReceiveIndex = 0;
		CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;

#ifdef FAST_SIB
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
				Grids[grid1]->GridData->
					CheckForOverlap(SiblingList[grid1].GridList[grid2],

							MetaData->LeftFaceBoundaryCondition,
							MetaData->RightFaceBoundaryCondition,
							&grid::CopyOverlappingMassFieldNoStar);
#else
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
				Grids[grid1]->GridData->
					CheckForOverlap(Grids[grid2]->GridData,
							MetaData->LeftFaceBoundaryCondition,
							MetaData->RightFaceBoundaryCondition,
							&grid::CopyOverlappingMassFieldNoStar);
#endif
		fprintf(stdout,"4-9\n"); // by YS

		CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
				Grids[grid1]->GridData->
					CheckForOverlap(SiblingList[grid1].GridList[grid2],
							MetaData->LeftFaceBoundaryCondition,
							MetaData->RightFaceBoundaryCondition,
							&grid::CopyOverlappingMassFieldNoStar);
#else
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
			for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
				Grids[grid1]->GridData->
					CheckForOverlap(Grids[grid2]->GridData,
							MetaData->LeftFaceBoundaryCondition,
							MetaData->RightFaceBoundaryCondition,
							&grid::CopyOverlappingMassFieldNoStar);
#endif

		CommunicationReceiveHandler(NULL,NULL,FALSE,NULL,TRUE);

	} // ENDFOR grid batches
	LCAPERF_STOP("CopyOverlappingMassField");
	fprintf(stdout,"4-10?\n"); // by YS

	/************************************************************************/
	// No Star Ends
	/************************************************************************/
#endif
















#ifdef FORCE_BUFFER_PURGE
	CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
	CommunicationBarrier();
#endif

	CommunicationDirection = COMMUNICATION_SEND_RECEIVE;




	/************************************************************************/
	/* Here Active particles have the opportunity to deposit mass into the
		 temporary mass buffers.
		 */

	ActiveParticleDepositMass(Grids, MetaData, NumberOfGrids, LevelArray,
			level);


	/************************************************************************/
	/* Compute the potential for the top grid. */

	fprintf(stdout,"Proc:%d, Pre barrier\n",MyProcessorNumber); //by YS 
	//CommunicationBarrier(); //by YS 
	//fprintf(stdout,"Proc:%d, Post barrier\n",MyProcessorNumber);//by YS 

	if (level == 0) {
		TIME_MSG("ComputePotentialFieldLevelZero");
		LCAPERF_START("ComputePotentialFieldLevelZero");
		TIMER_START("ComputePotentialFieldLevelZero");
		if (traceMPI) 
			fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): CPFLZero "
					"(send-receive)\n", MyProcessorNumber);
#ifdef FAST_SIB
		ComputePotentialFieldLevelZero(MetaData, SiblingList,
				Grids, NumberOfGrids);
#else
		ComputePotentialFieldLevelZero(MetaData, Grids, NumberOfGrids);
#endif
		TIMER_STOP("ComputePotentialFieldLevelZero");
		LCAPERF_STOP("ComputePotentialFieldLevelZero");
	}

	//CommunicationBarrier();// by Jo
	fprintf(stdout,"4-11"); // by YS



	/*
	
	for (StartGrid = 0; StartGrid< NumberOfGrids; StartGrid+=GRIDS_PER_LOOP) {
		EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);
		for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
			if (Grids[grid1]->GridData->GetPotentialField()[0]==NULL) 
				fprintf(stdout,"\nProc:%d, Grav Field is Null.\n", MyProcessorNumber); //by YS
			else
				fprintf(stdout,"\nProc:%d, Grav Field: %e\n", MyProcessorNumber, Grids[grid1]->GridData->GetPotentialField()[0][0]); //by YS
			if (Grids[grid1]->GridData->GetPotentialField()[1]==NULL) 
				fprintf(stdout,"\nProc:%d, Grav Field No Star is Null.\n", MyProcessorNumber); //by YS
			else
				fprintf(stdout,"Proc:%d, Grav Field NoStar: %e\n", MyProcessorNumber, Grids[grid1]->GridData->GetPotentialField()[1][0]);//by YS
		}
	}
*/





	/************************************************************************/
	/* Compute a first iteration of the potential and share BV's. */
	int iterate;
	if (level > 0) {
		LCAPERF_START("SolveForPotential");
		TIMER_START("SolveForPotential");
		CopyPotentialFieldAverage = 1;
		for (iterate = 0; iterate < PotentialIterations; iterate++) {

			if (iterate > 0)
				CopyPotentialFieldAverage = 2;


			for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
				Grids[grid1]->GridData->SolveForPotential(level, EvaluateTime);
				if (CopyGravPotential)
					Grids[grid1]->GridData->CopyPotentialToBaryonField();
			}
			fprintf(stdout,"4-12"); // by YS

			if (traceMPI) fprintf(tracePtr, "ITPOT post-recv\n");

#ifdef FORCE_MSG_PROGRESS 
			CommunicationBarrier();
#endif




			TIME_MSG("CopyPotentialField");
			for (StartGrid = 0; StartGrid < NumberOfGrids; 
					StartGrid += GRIDS_PER_LOOP) {
				EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

#ifdef BITWISE_IDENTICALITY
				CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
#else
				CommunicationDirection = COMMUNICATION_POST_RECEIVE;
#endif
				CommunicationReceiveIndex = 0;
				CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
#ifdef FAST_SIB
				for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

					//fprintf(stderr, "#SIBSend on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);

					// for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
					for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
						Grids[grid1]->GridData->
							CheckForOverlap(SiblingList[grid1].GridList[grid2],
									MetaData->LeftFaceBoundaryCondition,
									MetaData->RightFaceBoundaryCondition,
									&grid::CopyPotentialField);

					grid2 = grid1;
					Grids[grid1]->GridData->
						CheckForOverlap(Grids[grid2]->GridData,
								MetaData->LeftFaceBoundaryCondition,
								MetaData->RightFaceBoundaryCondition,
								&grid::CopyPotentialField);

				} // ENDFOR grid1
#else
				for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
					for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
						Grids[grid1]->GridData->
							CheckForOverlap(Grids[grid2]->GridData,
									MetaData->LeftFaceBoundaryCondition,
									MetaData->RightFaceBoundaryCondition,
									&grid::CopyPotentialField);
#endif

				fprintf(stdout,"4-13"); // by YS
#ifndef BITWISE_IDENTICALITY
#ifdef FORCE_MSG_PROGRESS 
				CommunicationBarrier();
#endif

				if (traceMPI) fprintf(tracePtr, "ITPOT send\n");

				CommunicationDirection = COMMUNICATION_SEND;


#ifdef FAST_SIB
				for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

					//fprintf(stderr, "#SIBRecv on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);

					// for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
					for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
						Grids[grid1]->GridData->
							CheckForOverlap(SiblingList[grid1].GridList[grid2],
									MetaData->LeftFaceBoundaryCondition,
									MetaData->RightFaceBoundaryCondition,
									&grid::CopyPotentialField);

					grid2 = grid1;
					Grids[grid1]->GridData->
						CheckForOverlap(Grids[grid2]->GridData,
								MetaData->LeftFaceBoundaryCondition,
								MetaData->RightFaceBoundaryCondition,
								&grid::CopyPotentialField);

				} // ENDFOR grid1
#else
				for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
					for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
						Grids[grid1]->GridData->
							CheckForOverlap(Grids[grid2]->GridData,
									MetaData->LeftFaceBoundaryCondition,
									MetaData->RightFaceBoundaryCondition,
									&grid::CopyPotentialField);
#endif

				CommunicationReceiveHandler();
#endif





#ifdef NBODY
/******************************************************************************************/
				//No Star Starts
/******************************************************************************************/
#ifdef BITWISE_IDENTICALITY
				CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
#else
				CommunicationDirection = COMMUNICATION_POST_RECEIVE;
#endif
				CommunicationReceiveIndex = 0;
				CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
#ifdef FAST_SIB
				for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

					//fprintf(stderr, "#SIBSend on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);

					// for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
					for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
						Grids[grid1]->GridData->
							CheckForOverlap(SiblingList[grid1].GridList[grid2],
									MetaData->LeftFaceBoundaryCondition,
									MetaData->RightFaceBoundaryCondition,
									&grid::CopyPotentialFieldNoStar);

					grid2 = grid1;
					Grids[grid1]->GridData->
						CheckForOverlap(Grids[grid2]->GridData,
								MetaData->LeftFaceBoundaryCondition,
								MetaData->RightFaceBoundaryCondition,
								&grid::CopyPotentialFieldNoStar);

				} // ENDFOR grid1
#else
				for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
					for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
						Grids[grid1]->GridData->
							CheckForOverlap(Grids[grid2]->GridData,
									MetaData->LeftFaceBoundaryCondition,
									MetaData->RightFaceBoundaryCondition,
									&grid::CopyPotentialFieldNoStar);
#endif

				fprintf(stdout,"4-13"); // by YS
#ifndef BITWISE_IDENTICALITY
#ifdef FORCE_MSG_PROGRESS 
				CommunicationBarrier();
#endif

				if (traceMPI) fprintf(tracePtr, "ITPOT send\n");

				CommunicationDirection = COMMUNICATION_SEND;


#ifdef FAST_SIB
				for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

					//fprintf(stderr, "#SIBRecv on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);

					// for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
					for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
						Grids[grid1]->GridData->
							CheckForOverlap(SiblingList[grid1].GridList[grid2],
									MetaData->LeftFaceBoundaryCondition,
									MetaData->RightFaceBoundaryCondition,
									&grid::CopyPotentialFieldNoStar);

					grid2 = grid1;
					Grids[grid1]->GridData->
						CheckForOverlap(Grids[grid2]->GridData,
								MetaData->LeftFaceBoundaryCondition,
								MetaData->RightFaceBoundaryCondition,
								&grid::CopyPotentialFieldNoStar);

				} // ENDFOR grid1
#else
				for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
					for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
						Grids[grid1]->GridData->
							CheckForOverlap(Grids[grid2]->GridData,
									MetaData->LeftFaceBoundaryCondition,
									MetaData->RightFaceBoundaryCondition,
									&grid::CopyPotentialFieldNoStar);
#endif

				CommunicationReceiveHandler();
#endif



/******************************************************************************************/
				//No Star Ends
/******************************************************************************************/
#endif // ENDIF nbody
			} // ENDFOR grid batches
		} // ENDFOR iterations
		CopyPotentialFieldAverage = 0;
		TIMER_STOP("SolveForPotential");
		LCAPERF_STOP("SolveForPotential");
	} // ENDIF level > 0













	fprintf(stdout,"4-14"); // by YS
	/* if level > MaximumGravityRefinementLevel, then do final potential
		 solve (and acceleration interpolation) here rather than in the main
		 EvolveLevel since it involves communications. */

	if (reallevel > MaximumGravityRefinementLevel) {

		/* compute potential and acceleration on coarser level [LOCAL]
			 (but only if there is at least a subgrid -- it should be only
			 if there is a subgrrid on reallevel, but this is ok). */

		for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
			if (Grids[grid1]->NextGridNextLevel != NULL) {
				Grids[grid1]->GridData->SolveForPotential(MaximumGravityRefinementLevel);
				if (CopyGravPotential)
					Grids[grid1]->GridData->CopyPotentialToBaryonField();
				else
					Grids[grid1]->GridData->ComputeAccelerationField
						((HydroMethod == Zeus_Hydro) ? DIFFERENCE_TYPE_STAGGERED : 
						 DIFFERENCE_TYPE_NORMAL, MaximumGravityRefinementLevel);
			}

		fprintf(stdout,"4-15"); // by YS
		/* Interpolate potential for reallevel grids from coarser grids. */

		if (!CopyGravPotential) {

			int Dummy, GridCount;
			LevelHierarchyEntry *Temp, *LastTemp;
			HierarchyEntry *Temp3;
			LevelHierarchyEntry *FirstTemp = LevelArray[reallevel];

#ifdef FORCE_MSG_PROGRESS 
			CommunicationBarrier();
#endif

			do {

				GridCount = 0;
				CommunicationDirection = COMMUNICATION_POST_RECEIVE;
				CommunicationReceiveIndex = 0;
				CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
				Temp = FirstTemp;
				while (Temp != NULL && GridCount++ < GRIDS_PER_LOOP) {
					Temp3 = Temp->GridHierarchyEntry;
					for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
						Temp3 = Temp3->ParentGrid;
					Temp->GridData->InterpolateAccelerations(Temp3->GridData);
					Temp = Temp->NextGridThisLevel;
				} // ENDWHILE
				LastTemp = Temp;

				CommunicationDirection = COMMUNICATION_SEND;
				Temp = FirstTemp;
				while (Temp != LastTemp) {
					Temp3 = Temp->GridHierarchyEntry;
					for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
						Temp3 = Temp3->ParentGrid;
					Temp->GridData->InterpolateAccelerations(Temp3->GridData);
					Temp = Temp->NextGridThisLevel;
				}
				FirstTemp = LastTemp;

				CommunicationReceiveHandler();

			} while (LastTemp != NULL);

		} // end:  if (!CopyGravPotential)

	} // end: if (reallevel > MaximumGravityRefinementLevel)

	//CommunicationBarrier();// by YS Jo
	fprintf(stdout,"4-16"); // by YS


	// --------------------------------------------------
	// MEMORY LEAK FIX
	//
	// valgrind error: "1,388,304 (67,352 direct, 1,320,952 indirect)
	// bytes in 130 blocks are definitely lost in loss record 22 of 46"
	//
	// Adding missing delete [] () for Grids[] allocated in
	// GenerateGridArray()
	// --------------------------------------------------
	delete [] Grids;

	// --------------------------------------------------

	LCAPERF_STOP("PrepareDensityField");
	return SUCCESS;

}

