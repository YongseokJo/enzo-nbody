/***********************************************************************
/
/  GRID CLASS (SOLVE POISSON ON POTENTIAL FIELD)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int MultigridSolver(float *RHS, float *Solution, int Rank, int TopDims[],
		    float &norm, float &mean, int start_depth,
		    float tolerance, int max_iter);
extern "C" void FORTRAN_NAME(smooth2)(float *source, float *dest, int *ndim,
                                   int *sdim1, int *sdim2, int *sdim3);
 
#define TOLERANCE 2.0e-6
#define MAX_ITERATION 20
 
int grid::SolveForPotential(int level, FLOAT PotentialTime)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
#ifdef NBODY
  if (GravitatingMassField[0] == NULL)  // if this is not set we have nothing to do.
		return SUCCESS;
  if (GravitatingMassField[1] == NULL)  // if this is not set we have nothing to do.
		return SUCCESS;
#else
  if (GravitatingMassField == NULL)  // if this is not set we have nothing to do.
		return SUCCESS;
#endif

  LCAPERF_START("grid_SolveForPotential");
 
  /* declarations */
 
  int dim, size = 1, i;
  float tol_dim = TOLERANCE * POW(0.1, 3-GridRank);
  //  if (GridRank == 3)
  //    tol_dim = 1.0e-5;
 
  /* Compute adot/a at time = t+1/2dt (time-centered). */
 
  if (PotentialTime < 0)
    PotentialTime = Time + 0.5*dtFixed;
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(PotentialTime, &a, &dadt) == FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.\n");
    }
 
  /* Compute right hand side. */
 
  float InverseVolumeElement = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GravitatingMassFieldDimension[dim];
    InverseVolumeElement *= (GravitatingMassFieldDimension[dim]-1);
  }
  tol_dim = max(sqrt(float(size))*1e-6, tol_dim);
 
#ifdef NBODY
	float *rhs[2];
	rhs[0]	= new float[size];
	rhs[1]	= new float[size];
#else
  float *rhs = new float[size];
#endif
 
  float Constant = GravitationalConstant * InverseVolumeElement *
                   POW(GravitatingMassFieldCellSize, 2) / a;
 
#define NO_SMOOTH_SOURCE
#ifdef SMOOTH_SOURCE
 
#ifdef NBODY
	FORTRAN_NAME(smooth2)(GravitatingMassField[0], rhs[0], &GridRank,
			GravitatingMassFieldDimension,
			GravitatingMassFieldDimension+1,
			GravitatingMassFieldDimension+2);
	FORTRAN_NAME(smooth2)(GravitatingMassField[1], rhs[1], &GridRank,
			GravitatingMassFieldDimension,
			GravitatingMassFieldDimension+1,
			GravitatingMassFieldDimension+2);
#else
	FORTRAN_NAME(smooth2)(GravitatingMassField, rhs, &GridRank,
			GravitatingMassFieldDimension,
			GravitatingMassFieldDimension+1,
			GravitatingMassFieldDimension+2);
#endif
#if 0
  FORTRAN_NAME(smooth2)(rhs, GravitatingMassField, &GridRank,
			GravitatingMassFieldDimension,
			GravitatingMassFieldDimension+1,
			GravitatingMassFieldDimension+2);
  FORTRAN_NAME(smooth2)(GravitatingMassField, rhs, &GridRank,
			GravitatingMassFieldDimension,
			GravitatingMassFieldDimension+1,
			GravitatingMassFieldDimension+2);
#endif
  for (i = 0; i < size; i++)
    rhs[i] *= Constant;
 
#else /* SMOOTH_SOURCE */
 


		//by YS debug
	int ndiff=0;
	int npdiff=0;
	float diff=0;
	float pdiff=0;
  for (i = 0; i < size; i++) {
#ifdef NBODY
    rhs[0][i] = GravitatingMassField[0][i] * Constant;
    rhs[1][i] = GravitatingMassField[1][i] * Constant;

		diff = abs(GravitatingMassField[0][i] - GravitatingMassField[1][i]);
		pdiff = abs(GravitatingMassFieldParticles[0][i] - GravitatingMassFieldParticles[1][i]);
		if (diff>1e-5) {
			fprintf(stderr, "In SolveForPotential, diff = %f\n", diff);
			ndiff++;
		}
		if (pdiff>1e-5) {
			fprintf(stderr, "pdiff = %f\n", pdiff);
			npdiff++;
		}
#else
    rhs[i] = GravitatingMassField[i] * Constant;
#endif

	}
	/*
	fprintf(stderr, "ndiff = %d\n", ndiff);
	fprintf(stderr, "npdiff = %d\n", npdiff);
	*/

 
#endif /* SMOOTH_SOURCE */
 
  /* Restrict fields to lower resolution if desired. */
 
  int GravitySmooth = max(level - MaximumGravityRefinementLevel, 0);
  GravitySmooth = 0;
 
  /* Iterate with multigrid. */
 
  float norm = huge_number, mean = norm;
#ifdef UNUSED
  int iteration = 0;
#endif /* UNUSED */
#ifdef NBODY 
  if (MultigridSolver(rhs[0], PotentialField[0], GridRank,
		      GravitatingMassFieldDimension, norm, mean,
		      GravitySmooth, tol_dim, MAX_ITERATION) == FAIL) {
    ENZO_FAIL("Error in MultigridDriver.\n");
  }

  if (MultigridSolver(rhs[1], PotentialField[1], GridRank,
		      GravitatingMassFieldDimension, norm, mean,
		      GravitySmooth, tol_dim, MAX_ITERATION) == FAIL) {
    ENZO_FAIL("Error in MultigridDriver.\n");
  }
#else
  if (MultigridSolver(rhs, PotentialField, GridRank,
		      GravitatingMassFieldDimension, norm, mean,
		      GravitySmooth, tol_dim, MAX_ITERATION) == FAIL) {
    ENZO_FAIL("Error in MultigridDriver.\n");
  }
#endif
 
#ifdef UNUSED
  while (norm/mean > tol_dim) {
    if (MultigridSolver(rhs, PotentialField, GridRank,
			GravitatingMassFieldDimension, norm, mean,
			GravitySmooth) == FAIL) {
      ENZO_FAIL("Error in MultigridDriver.\n");
    }
    printf("%"ISYM" %"GSYM"\n", iteration, norm/mean);
    if (iteration++ > MAX_ITERATION) {
      ENZO_VFAIL("exceeding iteration count (%"ISYM")\n", iteration)
    }
  }
#endif /* UNUSED */
 
  /* Clean up. */
 
#ifdef NBODY
  delete [] rhs[0];
  delete [] rhs[1];
#else
  delete [] rhs;
#endif

#define NO_POTENTIALDEBUGOUTPUT
//#define POTENTIALDEBUGOUTPUT
#ifdef POTENTIALDEBUGOUTPUT
  for (int i=0;i<GridDimension[0]; i++) {
    int igrid = GRIDINDEX_NOGHOST(i,(GridEndIndex[0]+GridStartIndex[0])/2,(GridEndIndex[0]+GridStartIndex[0])/2);
#ifdef NBODY
		printf("i: %i \t SolvedSub %g\n", i, PotentialField[0][igrid]);
		printf("i: %i \t SolvedSub %g\n", i, PotentialField[1][igrid]);
#else
		printf("i: %i \t SolvedSub %g\n", i, PotentialField[igrid]);
#endif //NBODY
	}
  float maxPot=-1e30, minPot=1e30;    
  float maxGM=-1e30, minGM=1e30;
  for (int i=0;i<size; i++) {
#ifdef NBODY
		maxPot = max(maxPot,PotentialField[0][i]);
		minPot = min(minPot,PotentialField[0][i]);
    maxGM = max(maxGM,GravitatingMassField[0][i]);
    minGM = min(minGM,GravitatingMassField[0][i]);
#else
		maxPot = max(maxPot,PotentialField[i]);
		minPot = min(minPot,PotentialField[i]);
    maxGM = max(maxGM,GravitatingMassField[i]);
    minGM = min(minGM,GravitatingMassField[i]);
#endif // NBODY
  }
  if (debug1) printf("SolvedPotential: Potential minimum: %g \t maximum: %g\n", minPot, maxPot);
  if (debug1) printf("SolvedPotential: GM minimum: %g \t maximum: %g\n", minGM, maxGM);

#endif // POTENTIALDEBUGOUTPUT
 
  LCAPERF_STOP("grid_SolveForPotential");
  return SUCCESS;
}
 
