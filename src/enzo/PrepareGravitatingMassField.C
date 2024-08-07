/***********************************************************************
/
/  PREPARE THE GRAVITATING MASS FIELD
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "communication.h"

/* function prototypes */

#ifndef FAST_SIB
int CopyOverlappingParticleMassFields(grid* CurrentGrid,
				      TopGridData *MetaData,
				      LevelHierarchyEntry *LevelArray[],
				      int level);
int CopyOverlappingParticleMassFieldsNoStar(grid* CurrentGrid,
				      TopGridData *MetaData,
				      LevelHierarchyEntry *LevelArray[],
				      int level);
#endif
int DepositBaryons(HierarchyEntry *Grid, FLOAT When, bool NoStar);
 
/* EvolveHierarchy function */
 

int PrepareGravitatingMassField1(HierarchyEntry *Grid)
{

  /* declarations */

  int RefinementFactor = RefineBy;
  grid *CurrentGrid = Grid->GridData;

  /* Gravity: initialize and clear the gravitating mass field. */
		fprintf(stdout,"Proc:%d, 4-3-1\n",MyProcessorNumber); // by YS

  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
		fprintf(stdout,"Proc:%d, 4-3-2\n",MyProcessorNumber); // by YS
    if (CurrentGrid->InitializeGravitatingMassField(RefinementFactor) == FAIL){
      ENZO_FAIL("Error in grid->InitializeGravitatingMassField.\n");
    }
    CurrentGrid->ClearGravitatingMassField();
    //CurrentGrid->ClearGravitatingMassFieldNoStar();
		fprintf(stdout,"Proc:%d, 4-3-3\n",MyProcessorNumber); // by YS
	}

  /* Baryons: copy parent density (no interpolation) to regions in
     GravitatingMassField which are beyond the boundary of the current grid. */

  int CommunicationReceiveIndexLast = CommunicationReceiveIndex;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  if (Grid->ParentGrid != NULL)
   if (CurrentGrid->CopyParentToGravitatingFieldBoundary(
				         Grid->ParentGrid->GridData) == FAIL) {
     ENZO_FAIL("Error in grid->CopyParentToGravitatingFieldBoundary.\n");
   } 
  //  if (CommunicationReceiveIndex != CommunicationReceiveIndexLast)
  //    CommunicationReceiveCurrentDependsOn = CommunicationReceiveIndex-1;

  return SUCCESS;
}

/************************************************************************/

#ifdef FAST_SIB
int PrepareGravitatingMassField2a(HierarchyEntry *Grid, int grid1,
				 SiblingGridList SiblingList[],
				 TopGridData *MetaData, int level,
				 FLOAT When)
#else
int PrepareGravitatingMassField2a(HierarchyEntry *Grid, TopGridData *MetaData,
				 LevelHierarchyEntry *LevelArray[], int level,
				 FLOAT When)
#endif
{
 
  /* declarations */
 
  int grid2;
  grid *CurrentGrid = Grid->GridData;
 
  /* Baryons: deposit mass into GravitatingMassField. */
 
		fprintf(stdout,"Proc:%d, 4-4-1\n",MyProcessorNumber); // by YS
  // IF STATEMENT HERE TO MAKE IT SO NO GAS CONTRIBUTES TO GRAVITY
  if(!SelfGravityGasOff){
    if (DepositBaryons(Grid, When, FALSE) == FAIL) {
      ENZO_FAIL("Error in DepositBaryons\n");
      printf("      Potential calculated for the gas\n");
    }
  }
 
		fprintf(stdout,"Proc:%d, 4-4-2\n",MyProcessorNumber); // by YS
  /* Particles: go through all the other grids on this level and add all
     their overlapping GravitatingMassFieldParticles to this grid's
     GravitatingMassField.  Handle periodicity properly. */
 
//  fprintf(stderr, "  PGMF - CopyOverlappingParticleMassField\n");
 
#ifdef FAST_SIB
  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
    if (CurrentGrid->CheckForOverlap(SiblingList[grid1].GridList[grid2],
                                     MetaData->LeftFaceBoundaryCondition,
                                     MetaData->RightFaceBoundaryCondition,
                                     &grid::AddOverlappingParticleMassField)
        == FAIL) {
      fprintf(stderr, "Error in grid->AddOverlappingParticleMassFields.\n");
    }
  //  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
  //fprintf(stderr,"grid %i %i\n", grid1, grid2);

		fprintf(stdout,"Proc:%d, 4-4-3\n",MyProcessorNumber); // by YS
#else
  if (CopyOverlappingParticleMassFields(CurrentGrid, MetaData,
                                        LevelArray, level) == FAIL) {
    ENZO_FAIL("Error in CopyOverlappingParticleMassFields.\n");
  }
#endif

		fprintf(stdout,"Proc:%d, 4-4-4\n",MyProcessorNumber); // by YS
  /* If we are using comoving coordinates, we must adjust the source term. */
 
  if (CommunicationDirection == COMMUNICATION_SEND ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {

    if (ComovingCoordinates)
      if (CurrentGrid->ComovingGravitySourceTerm() == FAIL) {
	ENZO_FAIL("Error in grid->ComovingGravitySourceTerm.\n");
      }
 
  } // end: if (CommunicationDirection != COMMUNICATION_SEND)
 
		fprintf(stdout,"Proc:%d, 4-4-5\n",MyProcessorNumber); // by YS
  return SUCCESS;
}

/************************************************************************/

int PrepareGravitatingMassField2b(HierarchyEntry *Grid, int level)
{
 
  /* declarations */
 
  grid *CurrentGrid = Grid->GridData;

  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  if (level > 0)

    CurrentGrid->PreparePotentialField(Grid->ParentGrid->GridData);
 
  return SUCCESS;
}





















#ifdef NBODY
int PrepareGravitatingMassFieldNoStar1(HierarchyEntry *Grid)
{

  /* declarations */

  int RefinementFactor = RefineBy;
  grid *CurrentGrid = Grid->GridData;

  /* Gravity: initialize and clear the gravitating mass field. */

  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {

		// by YS debug
    if (CurrentGrid->InitializeGravitatingMassField(RefinementFactor) == FAIL){
      ENZO_FAIL("Error in grid->InitializeGravitatingMassField.\n");
    }

    CurrentGrid->ClearGravitatingMassFieldNoStar();
  }

  /* Baryons: copy parent density (no interpolation) to regions in
     GravitatingMassField which are beyond the boundary of the current grid. */

  int CommunicationReceiveIndexLast = CommunicationReceiveIndex;
	CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  if (Grid->ParentGrid != NULL)
   if (CurrentGrid->CopyParentToGravitatingFieldBoundaryNoStar(
				         Grid->ParentGrid->GridData) == FAIL) {
     ENZO_FAIL("Error in grid->CopyParentToGravitatingFieldBoundary.\n");
   }
  //  if (CommunicationReceiveIndex != CommunicationReceiveIndexLast)
  //    CommunicationReceiveCurrentDependsOn = CommunicationReceiveIndex-1;

  return SUCCESS;
}

/************************************************************************/

#ifdef FAST_SIB
int PrepareGravitatingMassFieldNoStar2a(HierarchyEntry *Grid, int grid1,
				 SiblingGridList SiblingList[],
				 TopGridData *MetaData, int level,
				 FLOAT When)
#else
int PrepareGravitatingMassFieldNoStar2a(HierarchyEntry *Grid, TopGridData *MetaData,
				 LevelHierarchyEntry *LevelArray[], int level,
				 FLOAT When)
#endif
{
 
  /* declarations */
 
  int grid2;
  grid *CurrentGrid = Grid->GridData;
 
  /* Baryons: deposit mass into GravitatingMassField. */
	if(!SelfGravityGasOff){
		if (DepositBaryons(Grid, When, TRUE) == FAIL) {
			ENZO_FAIL("Error in DepositBaryons\n");
			printf("      Potential calculated for the gas\n");
		}
  }
 
  /* Particles: go through all the other grids on this level and add all
     their overlapping GravitatingMassFieldParticles to this grid's
     GravitatingMassField.  Handle periodicity properly. */
 
//  fprintf(stderr, "  PGMF - CopyOverlappingParticleMassField\n");
 
#ifdef FAST_SIB
  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
    if (CurrentGrid->CheckForOverlap(SiblingList[grid1].GridList[grid2],
                                     MetaData->LeftFaceBoundaryCondition,
                                     MetaData->RightFaceBoundaryCondition,
                                     &grid::AddOverlappingParticleMassFieldNoStar)
        == FAIL) {
      fprintf(stderr, "Error in grid->AddOverlappingParticleMassFields.\n");
    }
  //  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
  //fprintf(stderr,"grid %i %i\n", grid1, grid2);

#else
  if (CopyOverlappingParticleMassFieldsNoStar(CurrentGrid, MetaData,
                                        LevelArray, level) == FAIL) { // by YS did not touch
    ENZO_FAIL("Error in CopyOverlappingParticleMassFields.\n");
  }
#endif

  /* If we are using comoving coordinates, we must adjust the source term. */
 
  if (CommunicationDirection == COMMUNICATION_SEND ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {

    if (ComovingCoordinates)
      if (CurrentGrid->ComovingGravitySourceTermNoStar() == FAIL) {
	ENZO_FAIL("Error in grid->ComovingGravitySourceTerm.\n");
      }
 
  } // end: if (CommunicationDirection != COMMUNICATION_SEND)
 
  return SUCCESS;
}

/************************************************************************/

int PrepareGravitatingMassFieldNoStar2b(HierarchyEntry *Grid, int level)
{
 
  /* declarations */
 
  grid *CurrentGrid = Grid->GridData;

  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  if (level > 0)

    CurrentGrid->PreparePotentialFieldNoStar(Grid->ParentGrid->GridData);
 
  return SUCCESS;
}
#endif





