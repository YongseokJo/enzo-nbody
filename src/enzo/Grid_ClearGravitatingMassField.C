/***********************************************************************
/
/  GRID CLASS (ALLOCATE AND CLEAR THE GRAVITATING MASS FIELD)
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
 
 
int grid::ClearGravitatingMassField()
{
 
  /* Return is this is not the right processor. */
	fprintf(stdout,"Proc:%d, In Clear\n", MyProcessorNumber); // by YS
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* Error check. */
 
  if (GravitatingMassFieldCellSize == FLOAT_UNDEFINED) {
    ENZO_FAIL("GravitatingMassField uninitialized.\n");
  }
 
  /* compute size of the gravitating mass field */
 
  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldDimension[dim];
 
	//fprintf(stdout,"Proc:%d, Size: %d\n", MyProcessorNumber, size); // by YS
  /* allocate and clear the field */
 
  //  if (GravitatingMassField != NULL)
  //    fprintf(stderr, "ClearGravitatingMassField: Warning! Field not NULL.\n");
 
		//fprintf(stdout,"Proc:%d, is this NULL\n", MyProcessorNumber); // by YS
#ifdef NBODY
  if (GravitatingMassField[0] == NULL) {
     GravitatingMassField[0] = new float[size];
#else
    GravitatingMassField = new float[size];
#endif
		}
  if (GravitatingMassField[0] == NULL) {
    ENZO_FAIL("malloc error (out of memory?)\n");
  }
 
  for (i = 0; i < size; i++) {
#ifdef NBODY
    GravitatingMassField[0][i] = 0.0;
    //GravitatingMassField[1][i] = 0.0;
#else
    GravitatingMassField[i] = 0.0;
#endif
	}
  return SUCCESS;
}

#ifdef NBODY
int grid::ClearGravitatingMassFieldNoStar()
{
 
  /* Return is this is not the right processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* Error check. */
 
  if (GravitatingMassFieldCellSize == FLOAT_UNDEFINED) {
    ENZO_FAIL("GravitatingMassField uninitialized.\n");
  }
 
  /* compute size of the gravitating mass field */
 
  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldDimension[dim];
 
  /* allocate and clear the field */
 
  //if (GravitatingMassField[1] != NULL)
		//fprintf(stderr, "ClearGravitatingMassField: Warning! Field not NULL.\n");
 
  if (GravitatingMassField[1] == NULL) 
     GravitatingMassField[1] = new float[size];

  if (GravitatingMassField[1] == NULL) 
    ENZO_FAIL("malloc error (out of memory?)\n");
 
  for (i = 0; i < size; i++) 
    GravitatingMassField[1][i] = 0.0;
 
  return SUCCESS;
}
#endif
