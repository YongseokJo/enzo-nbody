/***********************************************************************
/
/  LLF WITH ZERO ORDER RECONSTRUCTION
/
/  written by: Peng Wang
/  date:       May, 2007
/  modified1: J. S. Oishi (for MHD)
/  date:       Apr, 2011 
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "EOS.h"
#include "ReconstructionRoutines.h"

int llf_mhd(float **FluxLine, float **priml, float **primr, float **prim, int ActiveSize);
int LLF_Zero_MHD(float **prim, float **priml, float **primr,
	    float **species, float **colors,  float **FluxLine, int ActiveSize,
	    char direc, int ij, int ik)
{

  // compute priml and primr                                                                                
  int iprim;
  int idual = (DualEnergyFormalism) ? 1 : 0;
  for (int field = 0; field < NEQ_MHD-idual; field++) {
    for (int i = 0; i < ActiveSize+1; i++) {
      //      printf("NEQ_MHD = %d; NumberOfGhostZones =  %d; ActiveSize = %d \n", NEQ_MHD, NumberOfGhostZones, ActiveSize);  
      iprim = i + NumberOfGhostZones - 1;
      priml[field][i] = prim[field][iprim];
      primr[field][i] = prim[field][iprim+1];
    }
  }

  // compute FluxLine
  if (llf_mhd(FluxLine, priml, primr, prim, ActiveSize)==FAIL) {
    return FAIL;
  }

  float sum;
  if (NSpecies > 0) {
    for (int i = 0; i < ActiveSize+1; i++) {
      iprim = i+NumberOfGhostZones-1;
      for (int field = 0; field < NSpecies; field++) {
        if (FluxLine[iD][i] >= 0) {
          species[field][i] = prim[field+NEQ_MHD-idual][iprim  ];
        } else {
          species[field][i] = prim[field+NEQ_MHD-idual][iprim+1];
        }
      }
      sum = 0;
      for (int field = 0; field < NSpecies; field++) {
        sum += species[field][i];
      }
      for (int field = 0; field < NSpecies; field++) {
        species[field][i] /= sum;
      }
      for (int field = 0; field < NSpecies; field++) {
        FluxLine[field+NEQ_MHD][i] = FluxLine[iD][i]*species[field][i];
      }
    }
  }

  if (NColor > 0) {
    for (int i = 0; i < ActiveSize+1; i++) {
      iprim = i+NumberOfGhostZones-1;
      for (int field = 0; field < NColor; field++) {
        if (FluxLine[iD][i] >= 0) {
          colors[field][i] = prim[field+NEQ_MHD+NSpecies-idual][iprim  ];
        } else {
          colors[field][i] = prim[field+NEQ_MHD+NSpecies-idual][iprim+1];
        }
      }
      for (int field = 0; field < NColor; field++) {
        FluxLine[field+NEQ_MHD+NSpecies][i] = FluxLine[iD][i]*colors[field][i];
      }
    }
  }
  
  return SUCCESS;
}
