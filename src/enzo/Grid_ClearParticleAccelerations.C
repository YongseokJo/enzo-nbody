/***********************************************************************
/
/  GRID CLASS (ALLOCATES AND CLEARS THE ACCELERATION FIELD FOR PARTICLESS)
/
/  written by: Greg Bryan
/  date:       March, 1995
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
 
 
 
int grid::ClearParticleAccelerations()
{
 
  int i, dim;
 
  if (NumberOfParticles > 0)
 
    /* Loop over active dimension */
 
    for (dim = 0; dim < GridRank+ComputePotential; dim++) {
 
      /* Error check. */
 
      if (ParticleAcceleration[dim][0] != NULL)
	fprintf(stderr, "ClearParticleAccelerations: Field not NULL.\n");
 
      /* Allocate accleration field. */
#ifdef NBODY 
      //ParticleAcceleration[dim] = new float*[2];
      ParticleAcceleration[dim][0] = new float[NumberOfParticles];
      //ParticleAcceleration[dim][1] = new float[NumberOfParticles];
#else
      ParticleAcceleration[dim] = new float[NumberOfParticles];
#endif
 
      /* Clear it. */
 
      for (i = 0; i < NumberOfParticles; i++) {
#ifdef NBODY
				ParticleAcceleration[dim][1][i] = 0.0;
				//ParticleAcceleration[dim][0][i] = 0.0;
#else
				ParticleAcceleration[dim][i] = 0.0;
#endif
			}
 
    }

  if (NumberOfActiveParticles > 0)

    /* Loop over active dimension */
    
    for (dim = 0; dim < GridRank+ComputePotential; dim++) {
      
      /* Error check. */
      
      if (ActiveParticleAcceleration[dim] != NULL)
        fprintf(stderr, "ClearParticleAccelerations: Field not NULL.\n");
      
      /* Allocate accleration field. */
      
      ActiveParticleAcceleration[dim] = new float[NumberOfActiveParticles];
      
      /* Clear it. */
      
      for (i = 0; i < NumberOfActiveParticles; i++)
        ActiveParticleAcceleration[dim][i] = 0.0;
      
    }
 
  return SUCCESS;
}

#ifdef NBODY
 
int grid::ClearParticleAccelerationsNoStar()
{
 
  int i, dim;
 
  if (NumberOfParticles > 0)
 
    /* Loop over active dimension */
 
    for (dim = 0; dim < GridRank+ComputePotential; dim++) {
 
      /* Error check. */
 
      if (ParticleAcceleration[dim][1] != NULL)
	fprintf(stderr, "ClearParticleAccelerations: Field not NULL.\n");
 
      /* Allocate accleration field. */
#ifdef NBODY 
      //ParticleAcceleration[dim] = new float*[2];
      ParticleAcceleration[dim][1] = new float[NumberOfParticles];
      //ParticleAcceleration[dim][1] = new float[NumberOfParticles];
#else
      ParticleAcceleration[dim] = new float[NumberOfParticles];
#endif
 
      /* Clear it. */
 
      for (i = 0; i < NumberOfParticles; i++) {
#ifdef NBODY
				ParticleAcceleration[dim][1][i] = 0.0;
				//ParticleAcceleration[dim][0][i] = 0.0;
#else
				ParticleAcceleration[dim][i] = 0.0;
#endif
			}
 
    }

  if (NumberOfActiveParticles > 0)

    /* Loop over active dimension */
    
    for (dim = 0; dim < GridRank+ComputePotential; dim++) {
      
      /* Error check. */
      
      if (ActiveParticleAcceleration[dim] != NULL)
        fprintf(stderr, "ClearParticleAccelerations: Field not NULL.\n");
      
      /* Allocate accleration field. */
      
      ActiveParticleAcceleration[dim] = new float[NumberOfActiveParticles];
      
      /* Clear it. */
      
      for (i = 0; i < NumberOfActiveParticles; i++)
        ActiveParticleAcceleration[dim][i] = 0.0;
      
    }
 
  return SUCCESS;
}
#endif
