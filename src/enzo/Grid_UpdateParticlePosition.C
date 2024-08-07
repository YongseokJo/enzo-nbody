/***********************************************************************
	/
	/  GRID CLASS (UPDATE PARTICLE POSITION FOR VELOCITY)
	/
	/  written by: Greg Bryan
	/  date:       May, 1995
	/  modified1:  David Collins
	/  date:       July, 2009
	/              Added OffProcessorUpdate.  This is to fix an error in 
	/              DepositParticlePositions wherein particles need to be updated
	/              before being interpolated into the Parent grid's
	/              GravitatingMassFieldParticles, and the Parent may
	/              be on a different process.
	/           
	/
	/  PURPOSE:
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

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


int grid::UpdateParticlePosition(float TimeStep, int OffProcessorUpdate)
{

	/* Return if this doesn't concern us. */
	/* OffProcessorUpdate defaults to FALSE */
	if (ProcessorNumber != MyProcessorNumber && OffProcessorUpdate == FALSE)
		return SUCCESS;

	if (NumberOfParticles == 0 && NumberOfActiveParticles == 0) return SUCCESS;

	FLOAT a = 1.0, dadt;
	int i, dim;

	//  if (debug)
	//    printf("UpdateParticlePosition: moving %"ISYM" particles forward by %"FSYM".\n",
	//	   NumberOfParticles, TimeStep);

	/* If using comoving coordinates, divide the acceleration by a(t) first.
		 (We use abs(TimeStep) because this routine is occasionally used to
		 move particles forward dt and then immediately reverse this (-dt);
		 using abs(dt) keeps things consistent). */

	if (ComovingCoordinates)
		if (CosmologyComputeExpansionFactor(Time + 0.5*fabs(TimeStep), &a, &dadt)
				== FAIL) {
			ENZO_FAIL("Error in CsomologyComputeExpansionFactors.");
		}

	float Coefficient = TimeStep/a;
	/* Loop over dimensions. */

	if (NumberOfParticles > 0) {
		for (dim = 0; dim < GridRank; dim++) {

			/* Error check. */

			if (ParticleVelocity[dim] == NULL) {
				ENZO_FAIL("No ParticleVelocity present.");
			}

			/* update velocities. */


			for (i = 0; i < NumberOfParticles; i++)
				ParticlePosition[dim][i] += Coefficient*ParticleVelocity[dim][i];
			// by YS Jo
			//if ( ParticleType[i] == PARTICLE_TYPE_NBODY ) continue;
			//&& GridLevel != MaximumRefinementLevel )


			/* wrap particle positions for periodic case.
				 (now done in CommunicationTransferParticles) */

		}
	}

	if (NumberOfActiveParticles > 0) {
		for (i = 0; i < NumberOfActiveParticles; i++) {

			FLOAT* appos;
			float* apvel;
			appos = ActiveParticles[i]->ReturnPosition();
			apvel = ActiveParticles[i]->ReturnVelocity();

			for (dim = 0; dim < GridRank; dim++)
				appos[dim] += Coefficient*apvel[dim];

			ActiveParticles[i]->SetPosition(appos);

			FLOAT period[3];
			for (dim = 0; dim < 3; dim++) {
				period[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
			}
			ActiveParticles[i]->SetPositionPeriod(period);
		}
	}
	return SUCCESS;
}

#ifdef NBODY
int grid::UpdateParticlePositionNoStar(float TimeStep, int OffProcessorUpdate)
{

	/* Return if this doesn't concern us. */
	/* OffProcessorUpdate defaults to FALSE */
	if (ProcessorNumber != MyProcessorNumber && OffProcessorUpdate == FALSE)
		return SUCCESS;

	if (NumberOfParticles == 0 && NumberOfActiveParticles == 0) return SUCCESS;

	FLOAT a = 1.0, dadt;
	int i, dim;

	//  if (debug)
	//    printf("UpdateParticlePosition: moving %"ISYM" particles forward by %"FSYM".\n",
	//	   NumberOfParticles, TimeStep);

	/* If using comoving coordinates, divide the acceleration by a(t) first.
		 (We use abs(TimeStep) because this routine is occasionally used to
		 move particles forward dt and then immediately reverse this (-dt);
		 using abs(dt) keeps things consistent). */

	if (ComovingCoordinates)
		if (CosmologyComputeExpansionFactor(Time + 0.5*fabs(TimeStep), &a, &dadt)
				== FAIL) {
			ENZO_FAIL("Error in CsomologyComputeExpansionFactors.");
		}

	float Coefficient = TimeStep/a;
	/* Loop over dimensions. */

	if (NumberOfParticles > 0) {
		for (dim = 0; dim < GridRank; dim++) {

			/* Error check. */

			if (ParticleVelocity[dim] == NULL) {
				ENZO_FAIL("No ParticleVelocity present.");
			}

			/* update velocities. */


			for (i = 0; i < NumberOfParticles; i++) {
				/* Only particle not nbody will be updated */
				//&& GridLevel != MaximumRefinementLevel )
				if ( ParticleType[i] != PARTICLE_TYPE_NBODY && ParticleType[i] != PARTICLE_TYPE_NBODY_NEW )  {
					ParticlePosition[dim][i] += Coefficient*ParticleVelocity[dim][i];
				}
			} // ENDFOR particles
			/* wrap particle positions for periodic case.
				 (now done in CommunicationTransferParticles) */
		}
	}

	if (NumberOfActiveParticles > 0) {
		for (i = 0; i < NumberOfActiveParticles; i++) {

			FLOAT* appos;
			float* apvel;
			appos = ActiveParticles[i]->ReturnPosition();
			apvel = ActiveParticles[i]->ReturnVelocity();

			for (dim = 0; dim < GridRank; dim++)
				appos[dim] += Coefficient*apvel[dim];

			ActiveParticles[i]->SetPosition(appos);

			FLOAT period[3];
			for (dim = 0; dim < 3; dim++) {
				period[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
			}
			ActiveParticles[i]->SetPositionPeriod(period);
		}
	}
	return SUCCESS;
}
#endif
