/***********************************************************************
	/
	/  GRID CLASS (SEND FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
	/
	/  written by: Greg Bryan
	/  date:       December, 1997
	/  modified1:
	/
	/  PURPOSE:
	/
	/  INPUTS:
	/
 ************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
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
#include "CommunicationUtilities.h"

// function prototypes

extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
		int *sdim1, int *sdim2, int *sdim3,
		int *ddim1, int *ddim2, int *ddim3,
		int *sstart1, int *sstart2, int *sstart3,
		int *dstart1, int *dstart2, int *dststart3);

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
		int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */


int grid::CommunicationSendRegion(grid *ToGrid, int ToProcessor,int SendField,
		int NewOrOld, int RegionStart[], int RegionDim[])
{
#ifdef USE_MPI 
	MPI_Request  RequestHandle;
	MPI_Status Status;
	MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
	MPI_Arg Count;
	MPI_Arg Source;

	/* Return if not on processor. */

	if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
		return SUCCESS;

	//  if (MyProcessorNumber != ProcessorNumber && 
	//      MyProcessorNumber != ToProcessor)
	//    return SUCCESS;

	int index, field, dim, Zero[] = {0, 0, 0};

	// Compute size of region to transfer

	int NumberOfFields = ((SendField == ALL_FIELDS)? NumberOfBaryonFields : 1) *
		((NewOrOld == NEW_AND_OLD)? 2 : 1);

	if (SendField == ACCELERATION_FIELDS)
		NumberOfFields = GridRank;

	int RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];
	int TransferSize = RegionSize * NumberOfFields;

	/* MHD Dimension stuff */

	int MHDRegionDim[3][3], MHDRegionSize[3]={1,1,1};
	int MHDeRegionDim[3][3], MHDeRegionSize[3]={1,1,1};

	if( UseMHDCT ){
		//Account for face centered field.  Note that I don't want to communicate the OldCenteredField
		TransferSize += ((SendField == ALL_FIELDS)? 3*RegionSize : 0 )*
			((NewOrOld == NEW_AND_OLD)? 2 : 1);

		for(field =0; field<3;field++){
			for(dim=0;dim<3;dim++){
				MHDRegionDim[field][dim] = RegionDim[dim]+MHDAdd[field][dim];
				MHDRegionSize[field] *= MHDRegionDim[field][dim];
				MHDeRegionDim[field][dim]= RegionDim[dim]+( (field==dim)?0:1);
				MHDeRegionSize[field] *=MHDeRegionDim[field][dim];
			}

			TransferSize += ((SendField == ALL_FIELDS)? MHDRegionSize[field]: 0 )*
				((NewOrOld == NEW_AND_OLD)? 2 : 1);

		}//field

	}//if(UseMHDCT)

	// Allocate buffer
#define NBODY
#ifdef NBODY 
	float *buffer1 = NULL;
	if (CommunicationDirection == COMMUNICATION_RECEIVE)
		buffer1 = CommunicationReceiveBuffer[CommunicationReceiveIndex];
	else	   
		buffer1 = new float[TransferSize];
	float *buffer2 = NULL;
	if (CommunicationDirection == COMMUNICATION_RECEIVE)
		buffer2 = CommunicationReceiveBuffer1[CommunicationReceiveIndex];
	else	   
		buffer2 = new float[TransferSize];
#endif
	float *buffer = NULL;
	if (CommunicationDirection == COMMUNICATION_RECEIVE)
		buffer = CommunicationReceiveBuffer[CommunicationReceiveIndex];
	else	   
		buffer = new float[TransferSize];

	// If this is the from processor, pack fields

	if (MyProcessorNumber == ProcessorNumber) {

		//  printf("SendRegion: RegionStart = %"ISYM" %"ISYM" %"ISYM"\n", RegionStart[0], RegionStart[1], RegionStart[2]);

		index = 0;

		if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
			for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
				if (field == SendField || SendField == ALL_FIELDS) {
					FORTRAN_NAME(copy3d)(BaryonField[field], &buffer[index],
							GridDimension, GridDimension+1, GridDimension+2,
							RegionDim, RegionDim+1, RegionDim+2,
							Zero, Zero+1, Zero+2,
							RegionStart, RegionStart+1, RegionStart+2);
					index += RegionSize;
				}

		if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
			for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
				if (field == SendField || SendField == ALL_FIELDS) {
					FORTRAN_NAME(copy3d)(OldBaryonField[field], &buffer[index],
							GridDimension, GridDimension+1, GridDimension+2,
							RegionDim, RegionDim+1, RegionDim+2,
							Zero, Zero+1, Zero+2,
							RegionStart, RegionStart+1, RegionStart+2);
					index += RegionSize;
				}

		if(UseMHDCT && SendField == ALL_FIELDS){


			if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY ){

				for(field=0;field<3;field++){
					FORTRAN_NAME(copy3d)(MagneticField[field], &buffer[index],
							&MagneticDims[field][0], 
							&MagneticDims[field][1], 
							&MagneticDims[field][2], 
							&MHDRegionDim[field][0],
							&MHDRegionDim[field][1],
							&MHDRegionDim[field][2],
							Zero, Zero+1, Zero+2,
							RegionStart, RegionStart+1, RegionStart+2);
					index += MHDRegionSize[field];
				}
			}

			if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY){
				for(field=0;field<3;field++){
					FORTRAN_NAME(copy3d)(OldMagneticField[field], &buffer[index],
							&MagneticDims[field][0],
							&MagneticDims[field][1],
							&MagneticDims[field][2],
							&MHDRegionDim[field][0],
							&MHDRegionDim[field][1],
							&MHDRegionDim[field][2],
							Zero, Zero+1, Zero+2,
							RegionStart, RegionStart+1, RegionStart+2);
					index += MHDRegionSize[field];
				}
			}//new and old or old only
		} // if(UseMHDCT && SendField == ALL_FIELDS)

		if (SendField == GRAVITATING_MASS_FIELD_PARTICLES) {
#ifdef NBODY
			FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles[0], buffer1,
					GravitatingMassFieldParticlesDimension,
					GravitatingMassFieldParticlesDimension+1,
					GravitatingMassFieldParticlesDimension+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					RegionStart, RegionStart+1, RegionStart+2);
		FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles[1], buffer2,
				GravitatingMassFieldParticlesDimension,
				GravitatingMassFieldParticlesDimension+1,
				GravitatingMassFieldParticlesDimension+2,
				RegionDim, RegionDim+1, RegionDim+2,
				Zero, Zero+1, Zero+2,
				RegionStart, RegionStart+1, RegionStart+2);
#else
		FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles, buffer,
				GravitatingMassFieldParticlesDimension,
				GravitatingMassFieldParticlesDimension+1,
				GravitatingMassFieldParticlesDimension+2,
				RegionDim, RegionDim+1, RegionDim+2,
				Zero, Zero+1, Zero+2,
				RegionStart, RegionStart+1, RegionStart+2);
#endif
		}

		if (SendField == GRAVITATING_MASS_FIELD) {
#ifdef NBODY
			FORTRAN_NAME(copy3d)(GravitatingMassField[0], buffer1,
					GravitatingMassFieldDimension,
					GravitatingMassFieldDimension+1,
					GravitatingMassFieldDimension+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					RegionStart, RegionStart+1, RegionStart+2);
		FORTRAN_NAME(copy3d)(GravitatingMassField[1], buffer2,
				GravitatingMassFieldDimension,
				GravitatingMassFieldDimension+1,
				GravitatingMassFieldDimension+2,
				RegionDim, RegionDim+1, RegionDim+2,
				Zero, Zero+1, Zero+2,
				RegionStart, RegionStart+1, RegionStart+2);
#else
		FORTRAN_NAME(copy3d)(GravitatingMassField, buffer,
				GravitatingMassFieldDimension,
				GravitatingMassFieldDimension+1,
				GravitatingMassFieldDimension+2,
				RegionDim, RegionDim+1, RegionDim+2,
				Zero, Zero+1, Zero+2,
				RegionStart, RegionStart+1, RegionStart+2);
#endif
		}

		if (SendField == POTENTIAL_FIELD) {
#ifdef NBODY
			FORTRAN_NAME(copy3d)(PotentialField[0], buffer1,
					GravitatingMassFieldDimension,
					GravitatingMassFieldDimension+1,
					GravitatingMassFieldDimension+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					RegionStart, RegionStart+1, RegionStart+2);
			FORTRAN_NAME(copy3d)(PotentialField[1], buffer2,
					GravitatingMassFieldDimension,
					GravitatingMassFieldDimension+1,
					GravitatingMassFieldDimension+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					RegionStart, RegionStart+1, RegionStart+2);
#else
			FORTRAN_NAME(copy3d)(PotentialField, buffer,
					GravitatingMassFieldDimension,
					GravitatingMassFieldDimension+1,
					GravitatingMassFieldDimension+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					RegionStart, RegionStart+1, RegionStart+2);
#endif
		}

		if (SendField == ACCELERATION_FIELDS)
			for (dim = 0; dim < GridRank; dim++) {
#ifdef NBODY
				FORTRAN_NAME(copy3d)(AccelerationField[dim][0], &buffer1[index],
						GridDimension, GridDimension+1, GridDimension+2,
						RegionDim, RegionDim+1, RegionDim+2,
						Zero, Zero+1, Zero+2,
						RegionStart, RegionStart+1, RegionStart+2);
				FORTRAN_NAME(copy3d)(AccelerationField[dim][1], &buffer2[index],
						GridDimension, GridDimension+1, GridDimension+2,
						RegionDim, RegionDim+1, RegionDim+2,
						Zero, Zero+1, Zero+2,
						RegionStart, RegionStart+1, RegionStart+2);
#else
				FORTRAN_NAME(copy3d)(AccelerationField[dim], &buffer[index],
						GridDimension, GridDimension+1, GridDimension+2,
						RegionDim, RegionDim+1, RegionDim+2,
						Zero, Zero+1, Zero+2,
						RegionStart, RegionStart+1, RegionStart+2);
#endif
				index += RegionSize;
			}
	}

	/* Send buffer */

	/* Only send if processor numbers are not identical */

	if (ProcessorNumber != ToProcessor) {

#ifdef MPI_INSTRUMENTATION
		starttime = MPI_Wtime();
#endif

		//    fprintf(stderr, "P(%d) communication for %d floats from %d to %d (phase %d)\n",
		//    	    MyProcessorNumber, TransferSize, ProcessorNumber,
		//    	    ToProcessor, CommunicationDirection);

		/* Send the data if on send processor, but leave buffer until the data
			 has been transfered out. */

		if (MyProcessorNumber == ProcessorNumber) {
#ifdef MPI_INSTRUMENTATION
			if (traceMPI) 
				fprintf(tracePtr, "CSR Sending %"ISYM" floats from %"ISYM" to %"ISYM"\n", 
						TransferSize, MyProcessorNumber, ToProcessor);
#endif
#ifdef NBODY
			if (SendField == ACCELERATION_FIELDS || SendField == GRAVITATING_MASS_FIELD
					|| SendField == GRAVITATING_MASS_FIELD_PARTICLES) {
				CommunicationBufferedSend(buffer1, TransferSize, DataType, ToProcessor, 
						MPI_SENDREGION_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);
				CommunicationBufferedSend(buffer2, TransferSize, DataType, ToProcessor, 
						MPI_SENDREGION_NOSTAR_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);
			} else {
				CommunicationBufferedSend(buffer, TransferSize, DataType, ToProcessor, 
						MPI_SENDREGION_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);
			}
#else
			CommunicationBufferedSend(buffer, TransferSize, DataType, ToProcessor, 
					MPI_SENDREGION_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);
#endif
		}

		if (MyProcessorNumber == ToProcessor) {

			//      fprintf(stderr, "Waiting for %d floats at %d from %d\n", TransferSize, 
			//	      MyProcessorNumber, ProcessorNumber);

			/* Post the receive message without waiting for the message to
				 be received.  When the data arrives, this will be called again
				 in (the real) receive mode. */

			if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {

				//	printf("Posting receive from P%"ISYM" for %"ISYM" floats in "
				//	       "comm index %"ISYM"\n", ProcessorNumber, TransferSize, 
				//	       CommunicationReceiveIndex);

#ifdef NBODY
				if (SendField == ACCELERATION_FIELDS || SendField == GRAVITATING_MASS_FIELD
						|| SendField == GRAVITATING_MASS_FIELD_PARTICLES) {
					MPI_Irecv(buffer1, TransferSize, DataType, ProcessorNumber, 
							MPI_SENDREGION_TAG, MPI_COMM_WORLD, 
							CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
					CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer1;
					CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
						MPI_Irecv(buffer2, TransferSize, DataType, ProcessorNumber, 
								MPI_SENDREGION_NOSTAR_TAG, MPI_COMM_WORLD, 
								CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
					CommunicationReceiveBuffer1[CommunicationReceiveIndex] = buffer2;
					CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
						CommunicationReceiveCurrentDependsOn;
					CommunicationReceiveIndex++;
				} else {
					MPI_Irecv(buffer, TransferSize, DataType, ProcessorNumber, 
							MPI_SENDREGION_TAG, MPI_COMM_WORLD, 
							CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
					CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer;
					CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
						CommunicationReceiveCurrentDependsOn;
					CommunicationReceiveIndex++;		
				}

#else
				MPI_Irecv(buffer, TransferSize, DataType, ProcessorNumber, 
						MPI_SENDREGION_TAG, MPI_COMM_WORLD, 
						CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
				CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer;
				CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
					CommunicationReceiveCurrentDependsOn;
				CommunicationReceiveIndex++;
#endif
			}

			/* If in send-receive mode, then wait for the message now. */

			if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
#ifdef NBODY
				if (SendField == ACCELERATION_FIELDS || SendField == GRAVITATING_MASS_FIELD
						|| SendField == GRAVITATING_MASS_FIELD_PARTICLES) {
					MPI_Recv(buffer1, TransferSize, DataType, ProcessorNumber, 
							MPI_SENDREGION_TAG, MPI_COMM_WORLD, &Status);
					MPI_Recv(buffer2, TransferSize, DataType, ProcessorNumber, 
							MPI_SENDREGION_NOSTAR_TAG, MPI_COMM_WORLD, &Status);
				} else {
					MPI_Recv(buffer, TransferSize, DataType, ProcessorNumber, 
							MPI_SENDREGION_TAG, MPI_COMM_WORLD, &Status);
				}
#else
				MPI_Recv(buffer, TransferSize, DataType, ProcessorNumber, 
						MPI_SENDREGION_TAG, MPI_COMM_WORLD, &Status);
#endif
			}
		} // ENDIF ToProcessor


#ifdef MPI_INSTRUMENTATION
		endtime = MPI_Wtime();
		timer[5] += endtime-starttime;
		counter[5] ++;
		timer[6] += double(TransferSize);
		RecvComm += endtime-starttime;
		CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */


	} // ENDIF different processors

	/* If this is the to processor, and we're either in send-receive mode
		 or receive mode, then unpack the data. */

	if (MyProcessorNumber == ToProcessor &&
			(CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
			 CommunicationDirection == COMMUNICATION_RECEIVE)) {

		//    if (ToProcessor != ProcessorNumber)
		//      fprintf(stderr, "Received %d floats at %d from %d\n", TransferSize, 
		//	      MyProcessorNumber, ProcessorNumber);

		index = 0;

		if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
			for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
				if (field == SendField || SendField == ALL_FIELDS) {
					delete ToGrid->BaryonField[field];
					ToGrid->BaryonField[field] = new float[RegionSize];
					FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->BaryonField[field],
							RegionDim, RegionDim+1, RegionDim+2,
							RegionDim, RegionDim+1, RegionDim+2,
							Zero, Zero+1, Zero+2,
							Zero, Zero+1, Zero+2);
					index += RegionSize;
				}

		if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
			for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
				if (field == SendField || SendField == ALL_FIELDS) {
					delete ToGrid->OldBaryonField[field];
					ToGrid->OldBaryonField[field] = new float[RegionSize];
					FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->OldBaryonField[field],
							RegionDim, RegionDim+1, RegionDim+2,
							RegionDim, RegionDim+1, RegionDim+2,
							Zero, Zero+1, Zero+2,
							Zero, Zero+1, Zero+2);
					index += RegionSize;
				}

		if( UseMHDCT && SendField == ALL_FIELDS ){

			/* send Bf */
			if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY){
				for(field=0;field<3;field++){

					delete ToGrid->MagneticField[field];
					ToGrid->MagneticField[field] = new float[MHDRegionSize[field] ];
					FORTRAN_NAME(copy3d)(&buffer[index], MagneticField[field],
							&MHDRegionDim[field][0],
							&MHDRegionDim[field][1],
							&MHDRegionDim[field][2],
							&MHDRegionDim[field][0],
							&MHDRegionDim[field][1],
							&MHDRegionDim[field][2],
							Zero, Zero+1, Zero+2,
							Zero, Zero+1, Zero+2);
					index += MHDRegionSize[field];
				}
			}

			if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
				for(field=0;field<3;field++){
					if( OldMagneticField[field] != NULL )
						fprintf(stderr,"Warning: CommSendRegion: OldMagneticField != NULL" );
					delete ToGrid->OldMagneticField[field];
					ToGrid->OldMagneticField[field] = new float[MHDRegionSize[field]];
					FORTRAN_NAME(copy3d)(&buffer[index], OldMagneticField[field],
							&MHDRegionDim[field][0],
							&MHDRegionDim[field][1],
							&MHDRegionDim[field][2],
							&MHDRegionDim[field][0],
							&MHDRegionDim[field][1],
							&MHDRegionDim[field][2],
							Zero, Zero+1, Zero+2,
							Zero, Zero+1, Zero+2);
					index += MHDRegionSize[field];
				}
		} // if( UseMHDCT && SendField == ALL_FIELDS )

		if (SendField == GRAVITATING_MASS_FIELD_PARTICLES) {
#ifdef NBODY
			delete ToGrid->GravitatingMassFieldParticles[1];
			delete ToGrid->GravitatingMassFieldParticles[0];
			//ToGrid->GravitatingMassFieldParticles = new float*[2];
			ToGrid->GravitatingMassFieldParticles[0] = new float[RegionSize];
			ToGrid->GravitatingMassFieldParticles[1] = new float[RegionSize];
			FORTRAN_NAME(copy3d)(buffer1, ToGrid->GravitatingMassFieldParticles[0],
					RegionDim, RegionDim+1, RegionDim+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					Zero, Zero+1, Zero+2);
			FORTRAN_NAME(copy3d)(buffer2, ToGrid->GravitatingMassFieldParticles[1],
					RegionDim, RegionDim+1, RegionDim+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					Zero, Zero+1, Zero+2);
#else
			delete ToGrid->GravitatingMassFieldParticles;
			ToGrid->GravitatingMassFieldParticles = new float[RegionSize];
			FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassFieldParticles,
					RegionDim, RegionDim+1, RegionDim+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					Zero, Zero+1, Zero+2);
#endif
		}

		if (SendField == GRAVITATING_MASS_FIELD) {
#ifdef NBODY
			delete ToGrid->GravitatingMassField[0];
			delete ToGrid->GravitatingMassField[1];
			//ToGrid->GravitatingMassField = new float*[2];
			ToGrid->GravitatingMassField[0] = new float[RegionSize];
			ToGrid->GravitatingMassField[1] = new float[RegionSize];
			FORTRAN_NAME(copy3d)(buffer1, ToGrid->GravitatingMassField[1],
					RegionDim, RegionDim+1, RegionDim+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					Zero, Zero+1, Zero+2);
			FORTRAN_NAME(copy3d)(buffer2, ToGrid->GravitatingMassField[1],
					RegionDim, RegionDim+1, RegionDim+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					Zero, Zero+1, Zero+2);

#else
			delete ToGrid->GravitatingMassField;
			ToGrid->GravitatingMassField = new float[RegionSize];
			FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassField,
					RegionDim, RegionDim+1, RegionDim+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					Zero, Zero+1, Zero+2);
#endif
		}

		if (SendField == POTENTIAL_FIELD) {
#ifdef NBODY
			delete ToGrid->PotentialField[1];
			delete ToGrid->PotentialField[0];
			//ToGrid->PotentialField = new float*[2];
			ToGrid->PotentialField[0] = new float[RegionSize];
			ToGrid->PotentialField[1] = new float[RegionSize];
			FORTRAN_NAME(copy3d)(buffer1, ToGrid->PotentialField[0],
					RegionDim, RegionDim+1, RegionDim+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					Zero, Zero+1, Zero+2);
			FORTRAN_NAME(copy3d)(buffer2, ToGrid->PotentialField[1],
					RegionDim, RegionDim+1, RegionDim+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					Zero, Zero+1, Zero+2);
#else
			delete ToGrid->PotentialField;
			ToGrid->PotentialField = new float[RegionSize];
			FORTRAN_NAME(copy3d)(buffer, ToGrid->PotentialField,
					RegionDim, RegionDim+1, RegionDim+2,
					RegionDim, RegionDim+1, RegionDim+2,
					Zero, Zero+1, Zero+2,
					Zero, Zero+1, Zero+2);
#endif
		}

		if (SendField == ACCELERATION_FIELDS)
			for (dim = 0; dim < GridRank; dim++) {
#ifdef NBODY
				delete ToGrid->AccelerationField[dim][0];
				delete ToGrid->AccelerationField[dim][1];
				//ToGrid->AccelerationField[dim] = new float*[2];
				ToGrid->AccelerationField[dim][0] = new float[RegionSize];
				ToGrid->AccelerationField[dim][1] = new float[RegionSize];
				FORTRAN_NAME(copy3d)(&buffer1[index], ToGrid->AccelerationField[dim][0],
						RegionDim, RegionDim+1, RegionDim+2,
						RegionDim, RegionDim+1, RegionDim+2,
						Zero, Zero+1, Zero+2,
						Zero, Zero+1, Zero+2);
				FORTRAN_NAME(copy3d)(&buffer2[index], ToGrid->AccelerationField[dim][1],
						RegionDim, RegionDim+1, RegionDim+2,
						RegionDim, RegionDim+1, RegionDim+2,
						Zero, Zero+1, Zero+2,
						Zero, Zero+1, Zero+2);
#else
				delete ToGrid->AccelerationField[dim];
				ToGrid->AccelerationField[dim] = new float[RegionSize];
				FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->AccelerationField[dim],
						RegionDim, RegionDim+1, RegionDim+2,
						RegionDim, RegionDim+1, RegionDim+2,
						Zero, Zero+1, Zero+2,
						Zero, Zero+1, Zero+2);
#endif
				index += RegionSize;
			}

		/* Only delete the buffer if we're in receive mode (in send mode
			 it will be deleted by CommunicationBufferedSend and if we're in
			 post-receive mode then it will be deleted when we get to
			 receive-mode). */

		delete [] buffer;

	} // ENDIF unpack

#endif /* USE_MPI */ 

	return SUCCESS;
}
