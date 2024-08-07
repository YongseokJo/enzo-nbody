/***********************************************************************
/
/  GRID CLASS (SEARCH FOR ALL STAR PARTICLES AND RETURN HOW MANY)
/
/  written by: Greg Bryan
/  date:       September, 2000
/  modified1:  JHK & JHW (2009)
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/

#include <stdlib.h> 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ActiveParticle.h"

void grid::SetNewParticleIndex(int &NumberCount1, PINT &NumberCount2)
{
  int n, abstype;
  for (n = 0; n < NumberOfParticles; n++) 
    if (ParticleNumber[n] == INT_UNDEFINED) {
      abstype = ABS(ParticleType[n]);
#ifdef NBODY
      if (abstype == PARTICLE_TYPE_STAR || abstype == PARTICLE_TYPE_NBODY_NEW ||
#else
      if (abstype == PARTICLE_TYPE_STAR ||
#endif
	  (abstype >= PARTICLE_TYPE_MUST_REFINE &&
	   abstype != PARTICLE_TYPE_MBH))
	ParticleNumber[n] = NumberCount1++ + NumberCount2;
      else 
	ParticleNumber[n] = NumberCount1 + NumberCount2++;
//      printf("New star particle index = %d (%d %d)\n",
//	     ParticleNumber[n], NumberCount1, NumberCount2);
    }
  return;
}

#define NO_DEBUG

void grid::SetNewActiveParticleIndex(PINT &next_id)
{

  int n, abstype;
  int ori_count = next_id;

  for (n = 0; n < NumberOfActiveParticles; n++)
    if (ActiveParticles[n]->Identifier == INT_UNDEFINED) {
      ActiveParticles[n]->Identifier = next_id++;
#ifdef DEBUG
      std::cout << "SNPI[" << MyProcessorNumber << "] " << "GridID: "
        << this->ID << " APID: " << ActiveParticles[n]->Identifier
        << std::endl;
#endif
    }

  return;
}
