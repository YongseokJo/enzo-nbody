/***********************************************************************
/
/  GRID CLASS (ALLOCATE AND CLEAR THE GRAVITATING MASS FIELD PARTICLES)
/
/  written by: Greg Bryan
/  date:       June, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
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

#define NBODY
 
/* function prototypes */
 
int grid::ClearGravitatingMassFieldParticles()
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Error check. */
 
  if (GravitatingMassFieldParticlesCellSize == FLOAT_UNDEFINED) {
    ENZO_FAIL("GravitatingMassFieldParticles uninitialized.\n");
  }
 
  /* Compute size of the gravitating mass field. */
 
  int dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldParticlesDimension[dim];
 
  /* Allocate and clear the field. */
 
//  if (GravitatingMassFieldParticles != NULL)
//    fprintf(stderr, "ClearGravitatingMassField: Warning! Field not NULL.\n");
 
  if (GravitatingMassFieldParticles[0] == NULL) {
		/* by YS Jo, 0 for the original field; 1 for the gravity with stars */
#ifdef NBODY
		//GravitatingMassFieldParticles = new float*[2];
		GravitatingMassFieldParticles[0] = new float[size];
		//GravitatingMassFieldParticles[1] = new float[size];
#else
		GravitatingMassFieldParticles = new float[size];
#endif
	}
	if (GravitatingMassFieldParticles[0] == NULL) {
		ENZO_FAIL("malloc error (out of memory?)\n");
  }
 
  for (int i = 0; i < size; i++) {
#ifdef NBODY
    //GravitatingMassFieldParticles[1][i] = 0.0;
    GravitatingMassFieldParticles[0][i] = 0.0;
#else
    GravitatingMassFieldParticles[i] = 0.0;
#endif
	}
 
  return SUCCESS;
}

#ifdef NBODY 
int grid::ClearGravitatingMassFieldParticlesNoStar()
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Error check. */
 
  if (GravitatingMassFieldParticlesCellSize == FLOAT_UNDEFINED) {
    ENZO_FAIL("GravitatingMassFieldParticles uninitialized.\n");
  }
 
  /* Compute size of the gravitating mass field. */
 
  int dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldParticlesDimension[dim];
 
  /* Allocate and clear the field. */
 
//  if (GravitatingMassFieldParticles != NULL)
//    fprintf(stderr, "ClearGravitatingMassField: Warning! Field not NULL.\n");
 
	/*In principle, this field should be initialized in the above function.*/
  if (GravitatingMassFieldParticles[1] == NULL) {
		/* by YS Jo, 0 for the original field; 1 for the gravity with stars */
    //GravitatingMassFieldParticles = new float*[2];
    //GravitatingMassFieldParticles[0] = new float[size];
    GravitatingMassFieldParticles[1] = new float[size];
	}
  if (GravitatingMassFieldParticles[1] == NULL) 
    ENZO_FAIL("malloc error (out of memory?)\n");

  for (int i = 0; i < size; i++) 
    GravitatingMassFieldParticles[1][i] = 0.0;

  return SUCCESS;
}
#endif
