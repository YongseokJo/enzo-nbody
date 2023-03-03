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
 
  int i, j, dim;
 
  if (NumberOfParticles > 0)
 
    /* Loop over active dimension */
 
    for (dim = 0; dim < GridRank+ComputePotential; dim++) {
 
      /* Error check. */
 
      if (ParticleAcceleration[dim] != NULL)
				fprintf(stderr, "ClearParticleAccelerations: Field not NULL.\n");


			/* Allocate accleration field. */
			ParticleAcceleration[dim] = new float[NumberOfParticles];

#ifdef NBODY 
      if (ParticleAccelerationNoStar[dim] != NULL)
				fprintf(stderr, "ClearStarBackGroundAcceleration: Field not NULL.\n");

			/* Allocate accleration field. */
			ParticleAccelerationNoStar[dim] = new float[NumberOfParticles]; //by YS, we have to change this to
#endif

			/* Clear it. */
			//for (i = 0; i < NumberOfStars; i++) {
				//number of nbody stars
			//}

			for (i = 0; i < NumberOfParticles; i++) {
				ParticleAcceleration[dim][i] = 0.0;
#ifdef NBODY
				ParticleAccelerationNoStar[dim][i] = 0.0; //by YS, we have to change this to
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
 
	int i, j, dim;

	if (NumberOfParticles > 0)

		/* Loop over active dimension */

		for (dim = 0; dim < GridRank+ComputePotential; dim++) {

			/* Error check. */

			if (ParticleAcceleration[dim] != NULL)
				fprintf(stderr, "ClearParticleAccelerations: Field not NULL.\n");

			/* Allocate accleration field. */
			ParticleAccelerationNoStar[dim] = new float[NumberOfParticles]; //by YS, we have to change this to
																																			 //number of nbody stars

			/* Clear it. */
			for (i = 0; i < NumberOfStars; i++) {
					ParticleAccelerationNoStar[dim][i] = 0;  //by YS, we have to change this to
																										 //number of nbody stars
			}
		}

	return SUCCESS;
}
#endif
