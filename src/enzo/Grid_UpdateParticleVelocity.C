/***********************************************************************
	/
	/  GRID CLASS (UPDATE PARTICLE VELOCITY FROM ACCELERATIONS)
	/
	/  written by: Greg Bryan
	/  date:       May, 1995
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

/* function prototypes */

#define VELOCITY_METHOD3

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


int grid::UpdateParticleVelocity(float TimeStep)
{

	/* Return if this doesn't concern us. */

	if (ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	if ((NumberOfParticles == 0 && NumberOfActiveParticles == 0) || ParticleAcceleration[0] == NULL)
		return SUCCESS;

	FLOAT a = 1.0, dadt;
#if defined(VELOCITY_METHOD1) || defined(VELOCITY_METHOD2)
	float VelocityMidStep;
#endif
	int i, dim , dim1;

	FLOAT coef, coef1, coef2;

	/* If using comoving coordinates, divide by a(t) first. */

	if (ComovingCoordinates)
		if (CosmologyComputeExpansionFactor(Time + TimeStep, &a, &dadt)
				== FAIL) {
			ENZO_FAIL("Error in CsomologyComputeExpansionFactors.");
		}

	/* Loop over dimensions. */

	for (int dim = 0; dim < GridRank; dim++) {

		/* Error check. */

		if (ParticleAcceleration[dim] == NULL) {
			ENZO_FAIL("No ParticleAcceleration present.");
		}

		/* Update velocities.  */

		if (ComovingCoordinates) {

			coef = 0.5*dadt/a*TimeStep;
			coef1 = 1.0 - coef;
			coef2 = 1.0 / (1.0 + coef);

			/* If using comoving coordinates, subtract the (time-centered)
				 drag-like term and add the acceleration. The acceleration has
				 already been divided by a(t). */

			for (i = 0; i < NumberOfParticles; i++) {
#ifdef NBODY
				// by YS Jo
				if ( ParticleType[i] == PARTICLE_TYPE_NBODY || ParticleType[i] == PARTICLE_TYPE_NBODY_NEW ) continue;
				//&& GridLevel != MaximumRefinementLevel )
#endif				

#ifdef VELOCITY_METHOD1

				/* i) partially time-centered. */

				VelocityMidStep = ParticleVelocity[dim][i] +
					ParticleAcceleration[dim][i]*0.5*TimeStep;

#ifdef NBODY
				ParticleVelocity[dim][i] +=
					(-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;
#endif

#endif /* VELOCITY_METHOD1 */

#ifdef VELOCITY_METHOD2

				/* ii) partially backward. */

				VelocityMidStep = ParticleVelocity[dim][i] ;

#ifdef NBODY
				ParticleVelocity[dim][i] +=
					(-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;
#endif
				//ParticleVelocity[dim][i] +=
				//  (-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;

#endif /* VELOCITY_METHOD2 */

#ifdef VELOCITY_METHOD3

				/* iii) Semi-implicit way */
#ifdef NBODY 
				ParticleVelocity[dim][i] = (coef1*ParticleVelocity[dim][i] +
						ParticleAcceleration[dim][i]*TimeStep)*coef2;
#endif


#endif /* VELOCITY_METHOD3 */

			}
		}
		else {

			/* Otherwise, just add the acceleration. */

			for (i = 0; i < NumberOfParticles; i++) {
#ifdef NBODY
				if ( ParticleType[i] != PARTICLE_TYPE_NBODY  && ParticleType[i] != PARTICLE_TYPE_NBODY_NEW ) 
#endif
					ParticleVelocity[dim][i] += ParticleAcceleration[dim][i] * TimeStep;
			}
		}

	}


	if (ProblemType == 29)
		for (i = 0; i < NumberOfParticles; i++)
			printf("id=%"PISYM"  %"PSYM" %"PSYM" %"PSYM"  %"ESYM" %"ESYM" %"ESYM" \n", ParticleNumber[i],
					ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i],
					ParticleVelocity[0][i], ParticleVelocity[1][i], ParticleVelocity[2][i]);

	if (NumberOfActiveParticles > 0) {
		if (ComovingCoordinates) {
			if (CosmologyComputeExpansionFactor(Time + TimeStep, &a, &dadt)
					== FAIL) {
				ENZO_FAIL("Error in CsomologyComputeExpansionFactors.");
			}

			for (i = 0; i < NumberOfActiveParticles; i++) {

				float* apvel = ActiveParticles[i]->ReturnVelocity();

				for (dim = 0; dim < GridRank; dim++) {

#ifdef VELOCITY_METHOD1

					/* i) partially time-centered. */

					VelocityMidStep = apvel[dim] +
						ActiveParticleAcceleration[dim][i]*0.5*TimeStep;

					apvel[dim] +=
						(-VelocityMidStep*dadt/a + ActiveParticleAcceleration[dim][i]) * TimeStep;

#endif /* VELOCITY_METHOD1 */

#ifdef VELOCITY_METHOD2

					/* ii) partially backward. */

					VelocityMidStep = apvel[dim];

					apvel[dim] +=
						(-VelocityMidStep*dadt/a + ActiveParticleAcceleration[dim][i]) * TimeStep;

#endif /* VELOCITY_METHOD2 */

#ifdef VELOCITY_METHOD3

					/* iii) Semi-implicit way */

					apvel[dim] = (coef1*apvel[dim] + ActiveParticleAcceleration[dim][i]*TimeStep)*coef2;


#endif /* VELOCITY_METHOD3 */

				}
				ActiveParticles[i]->SetVelocity(apvel);
			}
		}
		else {

			/* Otherwise, just add the acceleration. */

			for (i = 0; i < NumberOfActiveParticles; i++) {
				float* apvel = ActiveParticles[i]->ReturnVelocity();

				for (dim = 0; dim < GridRank; dim++)
					apvel[dim] += ActiveParticleAcceleration[dim][i] * TimeStep;

				ActiveParticles[i]->SetVelocity(apvel);
			}
		}
	}

	return SUCCESS;
}
