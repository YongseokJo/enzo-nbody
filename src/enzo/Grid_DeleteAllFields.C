/***********************************************************************
/
/  GRID CLASS (REMOVE ALL FIELDS)
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
 
void grid::DeleteAllFields()
{
 
  int i, j;
 
  this->DeleteParticles();
 
  for (i = 0; i < MAX_DIMENSION; i++) {
#ifdef NBODY
    delete [] AccelerationField[i][0];
    delete [] AccelerationField[i][1];
    delete [] ParticleAcceleration[i];
		for (j = 0; j < HERMITE_ORDER; j++) {
			delete [] StarBackGroundAcceleration[i][j];
			StarBackGroundAcceleration[i][j] = NULL;
		}
    ParticleAcceleration[i]         = NULL;
    AccelerationField[i][0]         = NULL;
    AccelerationField[i][1]         = NULL;
#else
    delete [] ParticleAcceleration[i];
    ParticleAcceleration[i]      = NULL;
    delete [] AccelerationField[i];
    AccelerationField[i]         = NULL;
#endif
  }
#ifdef NBODY
  delete [] ParticleAcceleration[MAX_DIMENSION];
  ParticleAcceleration[MAX_DIMENSION] = NULL;
	for (j = 0; j < HERMITE_ORDER; j++) {
		delete [] StarBackGroundAcceleration[MAX_DIMENSION][j];
		StarBackGroundAcceleration[MAX_DIMENSION][j] = NULL;
	}
#else
  delete [] ParticleAcceleration[MAX_DIMENSION];
  ParticleAcceleration[MAX_DIMENSION] = NULL;
#endif

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
 
  for(i=0;i<3;i++){
    if(MagneticField[i] != NULL){
      delete [] MagneticField[i];
      MagneticField[i] = NULL;
    }
    if( ElectricField[i] != NULL ){
      delete [] ElectricField[i];
      ElectricField[i] = NULL;
    }
    if(OldMagneticField[i] != NULL){
      delete [] OldMagneticField[i];
      OldMagneticField[i] = NULL;
    }
    if(OldElectricField[i] != NULL){
      delete [] OldElectricField[i];
      OldElectricField[i] = NULL;
    }

    if( AvgElectricField[i] != NULL ){
      delete[] AvgElectricField[i];
      AvgElectricField[i] = NULL;
    }
  }

#ifdef NBODY
  delete [] PotentialField[0];
  delete [] PotentialField[1];
  delete [] GravitatingMassField[0];
  delete [] GravitatingMassField[1];
  delete [] GravitatingMassFieldParticles[0];
  delete [] GravitatingMassFieldParticles[1];
 
  PotentialField[0]                = NULL;
  PotentialField[1]                = NULL;
  GravitatingMassField[0]          = NULL;
  GravitatingMassField[1]          = NULL;
  GravitatingMassFieldParticles[0] = NULL;
  GravitatingMassFieldParticles[1] = NULL;
#else
  delete [] PotentialField;
  delete [] GravitatingMassField;
  delete [] GravitatingMassFieldParticles;
 
  PotentialField                = NULL;
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;
#endif
 
}
