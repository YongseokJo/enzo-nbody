/***********************************************************************
/
/  DEPOSIT PARTICLES INTO PARTICLEMASSFIELD IN THIS GRID AND ALL SUBGRIDS
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
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
#include "TopGridData.h"
#include "LevelHierarchy.h"
 
/* function prototypes */
 
#ifdef NBODY
int DepositParticleMassFieldChildren(HierarchyEntry *DepositGrid,
				     HierarchyEntry *Grid, FLOAT Time, bool NoStar);
#else
int DepositParticleMassFieldChildren(HierarchyEntry *DepositGrid,
				     HierarchyEntry *Grid, FLOAT Time);
#endif 
#ifdef NBODY
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT TimeMidStep, bool NoStar)
#else
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT TimeMidStep)
#endif
{
 
  /* Get the time and dt for this grid.  Compute time+1/2 dt. */
 
  if (TimeMidStep < 0)
    TimeMidStep =     Grid->GridData->ReturnTime() +
                  0.5*Grid->GridData->ReturnTimeStep();
 
  /* Initialize the gravitating mass field only if in send-receive mode
     (i.e. this routine is called only once) or if in the first of the
     three communication modes (post-receive). */

  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
			CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {

		/* Initialize the gravitating mass field parameters (if necessary). */

		if (Grid->GridData->InitializeGravitatingMassFieldParticles(RefineBy)
				== FAIL) {
			ENZO_FAIL("Error in grid->InitializeGravitatingMassFieldParticles.\n");
		}

		/* Clear the GravitatingMassFieldParticles. */

#ifdef NBODY
		if (!NoStar) {
#endif
			if (Grid->GridData->ClearGravitatingMassFieldParticles() == FAIL) {
				ENZO_FAIL("Error in grid->ClearGravitatingMassFieldParticles.\n");
			}
#ifdef NBODY
		} else {
			if (Grid->GridData->ClearGravitatingMassFieldParticlesNoStar() == FAIL) {
				ENZO_FAIL("Error in grid->ClearGravitatingMassFieldParticles.\n");
			}
		}
#endif
		fprintf(stdout,"4-1-1\n"); // by YS
 
//  fprintf(stderr, "--DepositParticleMassField (Send) Initialize & Clear\n");
 
  } // end: if (CommunicationDirection != COMMUNICATION_SEND)
 
  /* Deposit particles to GravitatingMassFieldParticles in this grid. */
 
//  fprintf(stderr, "--DepositParticleMassField Call DepositParticlePositions\n");
 
#ifdef NBODY
  if (Grid->GridData->DepositParticlePositions(Grid->GridData, TimeMidStep,
				 GRAVITATING_MASS_FIELD_PARTICLES, NoStar) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }
#else
  if (Grid->GridData->DepositParticlePositions(Grid->GridData, TimeMidStep,
				 GRAVITATING_MASS_FIELD_PARTICLES) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }
#endif

		fprintf(stdout,"4-1-2\n"); // by YS
 
  /* Recursively deposit particles in children (at TimeMidStep). */
  if (Grid->NextGridNextLevel != NULL) {
#ifdef NBODY 
		if (DepositParticleMassFieldChildren(Grid, Grid->NextGridNextLevel,
					TimeMidStep, NoStar)
				== FAIL) {
			ENZO_FAIL("Error in DepositParticleMassFieldChildren.\n");
		}
#else
		if (DepositParticleMassFieldChildren(Grid, Grid->NextGridNextLevel,
					TimeMidStep)
				== FAIL) {
			ENZO_FAIL("Error in DepositParticleMassFieldChildren.\n");
		}
#endif
	}
 
  return SUCCESS;
}
 
 
 
#ifdef NBODY 
int DepositParticleMassFieldChildren(HierarchyEntry *DepositGrid,
				     HierarchyEntry *Grid, FLOAT DepositTime, bool NoStar)
#else
int DepositParticleMassFieldChildren(HierarchyEntry *DepositGrid,
				     HierarchyEntry *Grid, FLOAT DepositTime)
#endif
{
 
  /* Deposit particles in Grid into DepositGrid at the given time. */
 
  if (Grid->GridData->DepositParticlePositions(DepositGrid->GridData,
		     DepositTime, GRAVITATING_MASS_FIELD_PARTICLES, NoStar) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }
 
  /* Next grid on this level. */
 
  if (Grid->NextGridThisLevel != NULL)
    if (DepositParticleMassFieldChildren(DepositGrid, Grid->NextGridThisLevel,
					 DepositTime, NoStar) == FAIL) {
      ENZO_FAIL("Error in DepositParticleMassFieldChildren(1).\n");
    }
 
  /* Recursively deposit particles in children. */
 
  if (Grid->NextGridNextLevel != NULL)
    if (DepositParticleMassFieldChildren(DepositGrid, Grid->NextGridNextLevel,
					 DepositTime,NoStar) == FAIL) {
      ENZO_FAIL("Error in DepositParticleMassFieldChildren(2).\n");

    }
 
 
  return SUCCESS;
}



