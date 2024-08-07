/***********************************************************************
/
/  GRID CLASS (DEPOSIT LOCAL PARTICLE POSITIONS ONTO THE SPECIFIED FIELD)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  May, 2009 by John Wise
/                Only consider particles local on the processor.  This
/                happens in FindSubgrids (RebuildHierarchy) when we
/                keep the particles on the subgrid's local processor.
/
/  PURPOSE:
/     This routine deposits the particle living in this grid into either
/       the GravitatingMassField or the GravitatingMassFieldParticles of
/       itself, depending on the value of DepositField.
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* This controls the maximum particle mass which will be deposited in
   the MASS_FLAGGING_FIELD.  Only set in Grid_SetFlaggingField. */
 
extern float DepositParticleMaximumParticleMass;
 
int grid::DepositParticlePositionsLocal(FLOAT DepositTime, int DepositField,
					bool BothFlags)
{
 
  /* Declarations. */
 
  int dim, i, size;
  float MassFactor = 1.0, *ParticleMassTemp, *ParticleMassPointer;
//#ifdef NBODY
//  float *ParticleMassTempNoStar, *ParticleMassPointerNoStar;
//#endif
 
  /* If there are no particles, don't deposit anything. */
 
  if (NumberOfParticles == 0 && NumberOfActiveParticles == 0)
    return SUCCESS;
 
  /* If the Target is this grid and the DepositField is MassFlaggingField,
     then multiply the Particle density by the volume to get the mass. */
 
  if (DepositField == MASS_FLAGGING_FIELD || 
      DepositField == PARTICLE_MASS_FLAGGING_FIELD)
    for (dim = 0; dim < GridRank; dim++)
      MassFactor *= CellWidth[dim][0];
 
  /* If required, Change the mass of particles in this grid. */
 
	if (MassFactor != 1.0) {
		ParticleMassTemp = new float[NumberOfParticles];
		//#ifdef NBODY
		//		ParticleMassTempNoStar = new float[NumberOfParticles];
		//#endif
		for (i = 0; i < NumberOfParticles; i++) {
			ParticleMassTemp[i] = ParticleMass[i]*MassFactor;
			/*#ifdef NBODY
				if ( ParticleType == PARTICLE_TYPE_STAR)
				ParticleMassTempNoStar[i] = 0;
				else
				ParticleMassTempNoStar[i] = ParticleMass[i]*MassFactor;
				#endif */
		}
		ParticleMassPointer = ParticleMassTemp;
		//#ifdef NBODY
		//		ParticleMassPointerNoStar = ParticleMassTempNoStar;
		//#endif

	} else
		/*#ifdef NBODY
			{
			ParticleMassPointer = ParticleMass;
			for (i = 0; i < NumberOfParticles; i++) {
			if ( ParticleType == PARTICLE_TYPE_STAR)
			ParticleMassTempNoStar[i] = 0;
			else
			ParticleMassTempNoStar[i] = ParticleMass[i];
			}
			}
			#else*/
		ParticleMassPointer = ParticleMass;
	//#endif

  /* Allocate and fill the ActiveParticleMassPointer, obtain
     ActiveParticlePosition from the grid object */
  
  float* ActiveParticleMassPointer = new float[NumberOfActiveParticles];
  for (i = 0; i < NumberOfActiveParticles; i++)
    ActiveParticleMassPointer[i] = ActiveParticles[i]->ReturnMass()*MassFactor;

  FLOAT** ActiveParticlePosition = new FLOAT*[GridRank];
  for (dim = 0; dim < GridRank; dim++)
    ActiveParticlePosition[dim] = new FLOAT[NumberOfActiveParticles];

  this->GetActiveParticlePosition(ActiveParticlePosition);
 
  /* If the target field is MASS_FLAGGING_FIELD, then set masses of
     particles which are too large to zero (to prevent run-away refinement). */
 
  if ((DepositField == MASS_FLAGGING_FIELD ||
       DepositField == PARTICLE_MASS_FLAGGING_FIELD) &&
      DepositParticleMaximumParticleMass > 0 && MassFactor != 1.0) {
    for (i = 0; i < NumberOfParticles; i++)
      ParticleMassPointer[i] = min(DepositParticleMaximumParticleMass,
                                   ParticleMassPointer[i]);
    for (i = 0; i < NumberOfActiveParticles; i++)
      ActiveParticleMassPointer[i] = min(DepositParticleMaximumParticleMass,
                                         ActiveParticleMassPointer[i]);
  }
 
  /* Compute difference between current time and DepositTime. */
 
  float TimeDifference = DepositTime - Time;
 
  /* Move particles to positions at Time + TimeDifference. */
 
  this->UpdateParticlePosition(TimeDifference);
 
  /* Deposit particles. */
 
//  fprintf(stderr, "----DPP Call this->DepositPositions with NP = %"ISYM"\n", NumberOfParticles);
/*#ifdef NBODY 
  if (this->DepositPositions(ParticlePosition, ParticleMassPointer, ParticleMassPointerNoStar,
			     NumberOfParticles, DepositField) == FAIL) {
    ENZO_FAIL("Error in grid->DepositPositions\n");
  }
#else*/
  if (this->DepositPositions(ParticlePosition, ParticleMassPointer,
			     NumberOfParticles, DepositField) == FAIL) {
    ENZO_FAIL("Error in grid->DepositPositions\n");
  }
//#endif

  if (this->DepositPositions(ActiveParticlePosition, ActiveParticleMassPointer,
                 NumberOfActiveParticles, DepositField) == FAIL) {
    ENZO_FAIL("Error in grid->DepositPositions\n");
  }

  /* If requested, only consider cells that have already been flagged. */
  if (BothFlags) {

    float *DepositFieldPointer;
    switch (DepositField) {
    case GRAVITATING_MASS_FIELD:
      DepositFieldPointer = GravitatingMassField;
      break;
    case GRAVITATING_MASS_FIELD_PARTICLES:
      DepositFieldPointer = GravitatingMassFieldParticles;
      break;
    case MASS_FLAGGING_FIELD:
      DepositFieldPointer = MassFlaggingField;
      break;
    case PARTICLE_MASS_FLAGGING_FIELD:
      DepositFieldPointer = ParticleMassFlaggingField;
      break;
    } // ENDSWITCH
    
    for (dim = 0, size = 1; dim < GridRank; dim++)
      size *= GridDimension[dim];
    for (i = 0; i < size; i++)
      if (FlaggingField[i] == 0)
	DepositFieldPointer[i] = 0.0;
  }
  
  /* Delete the ActiveParticlePosition array */

  for (dim = 0; dim < GridRank; dim++)
    delete [] ActiveParticlePosition[dim];
  delete [] ActiveParticlePosition;
  delete [] ActiveParticleMassPointer;
 
  /* If necessary, delete the particle mass temporary. */
 
  if (MassFactor != 1.0)
    delete [] ParticleMassTemp;

  if (BothFlags) {
    delete [] FlaggingField;
    FlaggingField = NULL;
  }

  /* Return particles to positions at Time. */
 
  this->UpdateParticlePosition(-TimeDifference);
 
  return SUCCESS;
}
