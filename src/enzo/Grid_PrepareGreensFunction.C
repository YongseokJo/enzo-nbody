/***********************************************************************
/
/  GRID CLASS (PREPARES THE GREENS FUNCTION IN REAL SPACE)
/
/  written by: Greg Bryan
/  date:       April, 1998
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
//  Compute derived quantites
//
 
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
#include "phys_constants.h"
 
int grid::PrepareGreensFunction()
{
 
  /* Declarations. */
 
  int i, j, k, dim;
 
  /* Error check. */
 

  if (PotentialField != NULL) {
    ENZO_FAIL("Potential field not null.\n");
  }
 
  if (GravityBoundaryType != TopGridPeriodic) {
    ENZO_VFAIL("GravityBoundaryType %"ISYM" not supported.\n",
	    GravityBoundaryType)
  }
 
  /* Compute size and allocate field with size of GravitatingMassField. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldDimension[dim];
 
#ifdef NBODY
  PotentialFieldNoStar = new float[size];
#endif
  PotentialField = new float[size];
 
  /* Set the constant to be used. */
 
  float GravConst_factor;
  if (GridRank == 3)
    GravConst_factor = -GravitationalConstant/(4.0*pi);
  if (GridRank == 2)
    GravConst_factor = -GravitationalConstant*0.5/pi;
  if (GridRank == 1)
    GravConst_factor = -GravitationalConstant*0.5;
 
  /* Set Greens' function. */
 
  int n = 0;
  float xpos, ypos, zpos, r;
  for (k = 0; k < GravitatingMassFieldDimension[2]; k++) {
    zpos = GravitatingMassFieldLeftEdge[2] +
            float(k)*GravitatingMassFieldCellSize;
    for (j = 0; j < GravitatingMassFieldDimension[1]; j++) {
      ypos = GravitatingMassFieldLeftEdge[1] +
	     float(j)*GravitatingMassFieldCellSize;
      for (i = 0; i < GravitatingMassFieldDimension[0]; i++, n++) {
	xpos = GravitatingMassFieldLeftEdge[0] +
	       float(i)*GravitatingMassFieldCellSize;
	r = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
	r = max(r, GravitatingMassFieldCellSize);
	r *= GravitatingMassFieldCellSize;
#ifdef NBODY
	if (GridRank == 3) {
	  PotentialFieldNoStar[n] = GravConst_factor/r;
	}
	if (GridRank == 2) {
	  PotentialFieldNoStar[n] = GravConst_factor*log(r);
	}
	if (GridRank == 1) {
	  PotentialFieldNoStar[n] = GravConst_factor*r;
	}
#endif
	if (GridRank == 3)
	  PotentialField[n] = GravConst_factor/r;
	if (GridRank == 2)
	  PotentialField[n] = GravConst_factor*log(r);
	if (GridRank == 1)
	  PotentialField[n] = GravConst_factor*r;

 
      }
    }
  }
 
  return SUCCESS;
}
 
