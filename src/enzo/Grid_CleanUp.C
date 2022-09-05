/***********************************************************************
/
/  GRID CLASS (BEFORE REBUILDING, REMOVED UNNEEDED ARRAYS)
/
/  written by: Greg Bryan
/  date:       June, 1995
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
 
 
void grid::CleanUp()
{
 
  int i;
 
  for (i = 0; i < MAX_DIMENSION; i++) {
#ifdef NBODY
    delete [] ParticleAcceleration[i][1];
    delete [] ParticleAcceleration[i][0];
    ParticleAcceleration[i][0]      = NULL;
    ParticleAcceleration[i][1]      = NULL;
#else
    delete [] ParticleAcceleration[i];
    ParticleAcceleration[i]      = NULL;
#endif
//    delete [] AccelerationField[i];
 
//    AccelerationField[i]         = NULL;
  }
#ifdef NBODY
  delete [] ParticleAcceleration[MAX_DIMENSION][0];
  delete [] ParticleAcceleration[MAX_DIMENSION][1];
  ParticleAcceleration[MAX_DIMENSION][0] = NULL;
  ParticleAcceleration[MAX_DIMENSION][1] = NULL;
#else
  delete [] ParticleAcceleration[MAX_DIMENSION];
  ParticleAcceleration[MAX_DIMENSION] = NULL;
#endif
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    delete [] OldBaryonField[i];
    OldBaryonField[i] = NULL;
  }
 
  delete [] GravitatingMassField;
  delete [] GravitatingMassFieldParticles;

#ifdef NBODY	
  GravitatingMassField[0]          = NULL;
  GravitatingMassField[1]          = NULL;
  GravitatingMassFieldParticles[0] = NULL;
  GravitatingMassFieldParticles[1] = NULL;
#else
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;
#endif

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++)
    if (OldAccelerationField[i] != NULL) {
      delete [] OldAccelerationField[i];
      OldAccelerationField[i] = NULL;
    }
#endif

  if( UseMHDCT ){
    for(i=0;i<3;i++){
      if( OldMagneticField[i] != NULL ){
	delete [] OldMagneticField[i];
	OldMagneticField[i] = NULL;
      }

      
    }
  }

}
