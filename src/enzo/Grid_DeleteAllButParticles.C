/***********************************************************************
/
/  GRID CLASS (REMOVE ALL FIELDS BUT LEAVE THE PARTICLES ALONE)
/
/  written by: Greg Bryan
/  date:       April, 1996
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
void grid::DeleteAllButParticles()
{
 
  int i, j;
 
  //  this->DeleteParticles();
 
  for (i = 0; i < MAX_DIMENSION; i++) {
#ifdef NBODY
    delete [] AccelerationFieldNoStar[i];
		if (ParticleAccelerationNoStar[i] != NULL) {
			delete [] ParticleAccelerationNoStar[i];
			ParticleAccelerationNoStar[i] = NULL;
		}
    AccelerationFieldNoStar[i]         = NULL;
#endif
    delete [] ParticleAcceleration[i];
    delete [] AccelerationField[i];
    AccelerationField[i]         = NULL;
    ParticleAcceleration[i]      = NULL;
  }
#ifdef NBODY
	delete [] ParticleAccelerationNoStar[MAX_DIMENSION];
	ParticleAccelerationNoStar[MAX_DIMENSION] = NULL;
#endif
  delete [] ParticleAcceleration[MAX_DIMENSION];
  ParticleAcceleration[MAX_DIMENSION] = NULL;
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    delete [] BaryonField[i];
    delete [] OldBaryonField[i];
    BaryonField[i]    = NULL;
    OldBaryonField[i] = NULL;
  }

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++)
    if (OldAccelerationField[i] != NULL) {
      delete [] OldAccelerationField[i];
      OldAccelerationField[i] = NULL;
    }
#endif
 


#ifdef NBODY
  delete [] PotentialFieldNoStar;
  delete [] GravitatingMassFieldNoStar;
  delete [] GravitatingMassFieldParticlesNoStar;
 
  PotentialFieldNoStar                = NULL;
  GravitatingMassFieldNoStar          = NULL;
  GravitatingMassFieldParticlesNoStar = NULL;
#endif

  delete [] PotentialField;
  delete [] GravitatingMassField;
  delete [] GravitatingMassFieldParticles;
 
  PotentialField                = NULL;
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;
 
}
