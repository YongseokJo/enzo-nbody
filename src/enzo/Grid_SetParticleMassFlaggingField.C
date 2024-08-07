/***********************************************************************
/
/  GRID CLASS (SET THE MASS FLAGGING FIELD FOR PARTICLES)
/
/  written by: John Wise
/  date:       May, 2009
/  modified1:  
/
/  PURPOSE:    This routine sums the particle mass flagging field in a 
/              non-blocking fashion.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
 
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
#include "communication.h"

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */
int Return_MPI_Tag(int grid_num, int proc);

/* The following is defined in Grid_DepositParticlePositions.C. */
 
extern float DepositParticleMaximumParticleMass;
 
int grid::SetParticleMassFlaggingField(int StartProc, int EndProc, int level, 
				       int ParticleMassMethod, int MustRefineMethod,
				       int *SendProcs, int NumberOfSends)
{

  //printf("grid::SetParticleMassFlaggingField called \n");
  int i, irecv, dim, size, proc, MPI_Tag;
  bool and_flag;

  /* Return if we're not needed here */

//  if (MyProcessorNumber == ProcessorNumber &&
//      CommunicationDirection == COMMUNICATION_SEND)
//    return SUCCESS;

  //printf("    grid::SetParticleMassFlaggingField   \n");

  if (MyProcessorNumber != ProcessorNumber &&
      (CommunicationDirection == COMMUNICATION_RECEIVE ||
       CommunicationDirection == COMMUNICATION_POST_RECEIVE))
    return SUCCESS;

//  printf("--> SetPMFlag[P%"ISYM"/%"ISYM"]: level %"ISYM", grid %"ISYM", "
//	 "comm_dir = %"ISYM", npart = %"ISYM"\n", 
//	 MyProcessorNumber, ProcessorNumber, level, GridNum, 
//	 CommunicationDirection, NumberOfParticles);

#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
  MPI_Arg Source;
  float *buffer = NULL;
#endif /* USE_MPI */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (CommunicationDirection == COMMUNICATION_SEND) {

    int method, NumberOfFlaggedCells;
    bool KeepFlaggingField;

    /* Calculate the flagging field only if 
       1) this grid belongs to this grid and it's the first pass, or
       2) this grid isn't local and this processor is between StartProc 
          and EndProc
    */

    // If local, this is already calculated on the first pass
    if (MyProcessorNumber == ProcessorNumber && StartProc > 0)
      return SUCCESS;

    // If not local, only calculate if the processor is within the
    // given range and has particles.
    if (MyProcessorNumber != ProcessorNumber &&
	(MyProcessorNumber < StartProc || MyProcessorNumber >= EndProc ||
	 NumberOfParticles == 0 && NumberOfActiveParticles == 0))
      return SUCCESS;

    /***********************************************************************/
    /* beginning of Cell flagging criterion routine                        */
    /*
       NOTE (JHW 05/2014): If using MustRefineParticlesCreateParticles
       > 0, then only flag cells for particle overdensity refinement
       if they are flagged by the must-refine particles on level ==
       MustRefineParticlesRefineToLevel ONLY so that the entire region
       occupied by the high-resolution region is not refined.  Thus,
       must-refine flagging must be done first.
    */
  
    /* Allocate and clear mass flagging field. */
 
    this->ClearParticleMassFlaggingField();

    /* ==== METHOD 8: BY POSITION OF MUST-REFINE PARTICLES  ==== */
 
    if (MustRefineMethod >= 0 &&
	level <= MustRefineParticlesRefineToLevel) {

      KeepFlaggingField = (level == MustRefineParticlesRefineToLevel);
      NumberOfFlaggedCells = 
	this->DepositMustRefineParticles(ParticleMassMethod, level,
					 KeepFlaggingField);

      if (NumberOfFlaggedCells < 0) {
	ENZO_FAIL("Error in grid->DepositMustRefineParticles.\n");
      }

    } // ENDIF MustRefineMethod
    
    /* ==== METHOD 4: BY PARTICLE MASS ==== */

    if (ParticleMassMethod >= 0 &&
	(level >= MustRefineParticlesRefineToLevel ||
	 MustRefineParticlesCreateParticles == 0)) {

      and_flag = (level == MustRefineParticlesRefineToLevel &&
		  MustRefineParticlesCreateParticles > 0);
      
      /* Set the maximum particle mass to be deposited (cleared below). */
      
      DepositParticleMaximumParticleMass =
	0.99999*MinimumMassForRefinement[ParticleMassMethod]*POW(RefineBy,
 	 level*MinimumMassForRefinementLevelExponent[ParticleMassMethod]);
 
      /* Deposit particles in this grid to MassFlaggingField. */
 
      this->DepositParticlePositionsLocal(this->ReturnTime(),
					  PARTICLE_MASS_FLAGGING_FIELD,
					  and_flag);
      
      DepositParticleMaximumParticleMass = 0;

    } // ENDIF ParticleMassMethod

    /* End of Cell flagging criterion routine                              */
    /***********************************************************************/

    /* Send the particle mass flagging field to the grid's processor */

#ifdef USE_MPI
    if (MyProcessorNumber != ProcessorNumber) {
      //MPI_Tag = Return_MPI_Tag(GridNum, MyProcessorNumber);
//      printf("----> SetPMFlag[P%"ISYM"/%"ISYM"]: sending %"ISYM" floats.\n",
//	     MyProcessorNumber, ProcessorNumber, size);
      CommunicationBufferedSend(ParticleMassFlaggingField, size, DataType,
				ProcessorNumber, MPI_SENDPMFLAG_TAG, 
				enzo_comm, size*sizeof(float));
      delete [] ParticleMassFlaggingField;
      ParticleMassFlaggingField = NULL;
    }
#endif /* USE_MPI */

    return SUCCESS;

  } // ENDIF: COMMUNICATION_SEND

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_RECEIVE) {

    /* Receive the data and sum it */

    buffer = CommunicationReceiveBuffer[CommunicationReceiveIndex];
    for (i = 0; i < size; i++)
      ParticleMassFlaggingField[i] += buffer[i];
    delete [] buffer;

  } // ENDIF receive

  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {

    /* Post receive requests */

    Count = size;
    for (proc = 0; proc < NumberOfSends; proc++) {
      Source = SendProcs[proc];
//      printf("----> SetPMFlag[P%"ISYM"/%"ISYM"]: "
//	     "posting receive for %"ISYM" floats, coming from P%"ISYM".\n",
//	     MyProcessorNumber, ProcessorNumber, size, Source);

      if (Source >= StartProc && Source < EndProc) {
	buffer = new float[size];
	MPI_Irecv(buffer, Count, DataType, Source, MPI_SENDPMFLAG_TAG, enzo_comm, 
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);

	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = NULL;
	CommunicationReceiveCallType[CommunicationReceiveIndex] = 13;
	CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer;
//	CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = StartProc;
//	CommunicationReceiveArgumentInt[1][CommunicationReceiveIndex] = EndProc;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;
      } // ENDIF inside processor range

    } // ENDFOR processors

  } // ENDIF post receive

#endif /* USE_MPI */
   
  return SUCCESS;
 
}

/************************************************************************/

int Return_MPI_Tag(int grid_num, int proc)
{
  // Return a somewhat-unique MPI tag for communication.  The factors
  // are prime.
  return 6373*MPI_SENDPMFLAG_TAG + 4041*grid_num + 1973*proc;
}

#ifdef UNUSED
void InitializeParticleMassFlaggingFieldCommunication(void)
{
#ifdef USE_MPI
  int i;
  irq = 0;
  for (i = 0; i < MAX_REQUEST_HANDLES; i++) {
    FlagRequestHandle[i] = NULL;
    if (FlagBuffer[i] != NULL)

      delete [] FlagBuffer[i];
    FlagBuffer[i] = NULL;
  }
#endif USE_MPI
  return;
}
#endif
