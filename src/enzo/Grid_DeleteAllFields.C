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
    delete [] ParticleAcceleration[i];
    delete [] AccelerationField[i];
    ParticleAcceleration[i]      = NULL;
    AccelerationField[i]         = NULL;
  }

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

  delete [] PotentialField;
  delete [] GravitatingMassField;
  delete [] GravitatingMassFieldParticles;
 
  PotentialField                = NULL;
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;
 
}


#ifdef NBODY
void grid::DeleteAllFieldsNoStar()
{

  int i, j;

  for (i = 0; i < MAX_DIMENSION; i++) {

    delete [] AccelerationFieldNoStar[i];
		if (ParticleAccelerationNoStar[i] != NULL) {
		delete [] ParticleAccelerationNoStar[i];
		ParticleAccelerationNoStar[i] = NULL;
		}
    AccelerationFieldNoStar[i]         = NULL;

  }

	delete [] ParticleAccelerationNoStar[MAX_DIMENSION];
	ParticleAccelerationNoStar[MAX_DIMENSION] = NULL;


  delete [] PotentialFieldNoStar;
  delete [] GravitatingMassFieldNoStar;
  delete [] GravitatingMassFieldParticlesNoStar;

  PotentialFieldNoStar                = NULL;
  GravitatingMassFieldNoStar          = NULL;
  GravitatingMassFieldParticlesNoStar = NULL;

}
#endif
