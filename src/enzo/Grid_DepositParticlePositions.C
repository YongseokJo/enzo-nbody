/***********************************************************************
	/
	/  GRID CLASS (DEPOSIT PARTICLE POSITIONS ONTO THE SPECIFIED FIELD)
	/
	/  written by: Greg Bryan
	/  date:       May, 1995
	/  modified1:
	/
	/  PURPOSE:
	/     This routine deposits the particle living in this grid into either
	/       the GravitatingMassField or the GravitatingMassFieldParticles of
	/       the TargetGrid, depending on the value of DepositField.
	/     It also moves the particles in this grid so they are at the correct
	/       time and adjusts their 'mass' to be appropriate for TargetGrid
	/       (since particle 'mass' changed from level to level).
	/
	/  NOTE:
	/
 ************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

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
#include "ActiveParticle.h"
#include "communication.h"

/* function prototypes */

extern "C" void PFORTRAN_NAME(cic_deposit)(FLOAT *posx, FLOAT *posy,
		FLOAT *posz, int *ndim, int *npositions,
		float *densfield, float *field, FLOAT *leftedge,
		int *dim1, int *dim2, int *dim3, float *cellsize,
		float *cloudsize);
extern "C" void PFORTRAN_NAME(ngp_deposit)(FLOAT *posx, FLOAT *posy,
		FLOAT *posz, int *ndim, int *npositions,
		float *densfield, float *field, FLOAT *leftedge,
		int *dim1, int *dim2, int *dim3, float *cellsize);

extern "C" void PFORTRAN_NAME(smooth_deposit)(FLOAT *posx, FLOAT *posy,
		FLOAT *posz, int *ndim, int *npositions,
		float *densfield, float *field, FLOAT *leftedge,
		int *dim1, int *dim2, int *dim3, float *cellsize,
		float *rsmooth);

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
		int Target, int Tag, MPI_Comm CommWorld, 
		int BufferSize);
#endif /* USE_MPI */
double ReturnWallTime(void);

/* This controls the maximum particle mass which will be deposited in
	 the MASS_FLAGGING_FIELD.  Only set in Grid_SetFlaggingField. */

float DepositParticleMaximumParticleMass = 0;
#ifdef NBODY
int grid::DepositParticlePositions(grid *TargetGrid, FLOAT DepositTime,
		int DepositField, bool NoStar)
#else
int grid::DepositParticlePositions(grid *TargetGrid, FLOAT DepositTime,
		int DepositField)
#endif
{

	/* Return if this doesn't concern us. */

	if (TargetGrid->CommunicationMethodShouldExit(this) ||
			(NumberOfParticles == 0 && NumberOfActiveParticles == 0))
		return SUCCESS;

	//  fprintf(stderr, "----DPP: MyPN = %"ISYM", PN = %"ISYM", TGPN = %"ISYM", DIR (R=1,S=2) = %"ISYM", NP = %"ISYM"\n",
	//	  MyProcessorNumber, ProcessorNumber, TargetGrid->ProcessorNumber, CommunicationDirection, NumberOfParticles);

	/* Declarations. */

	int dim, i, j, k, size, index1, index2;
	int Dimension[] = {1,1,1};
	int OriginalDimension[] = {1,1,1};
	long_int Offset[] = {0,0,0};
	float MassFactor = 1.0, *ParticleMassTemp, *ParticleMassPointer; 
	FLOAT CellSize, CloudSize;
	float *ParticleMassPointerSink;    
	float TimeDifference = 0;
	FLOAT LeftEdge[MAX_DIMENSION], OriginalLeftEdge[MAX_DIMENSION];
	float *DepositFieldPointer, *OriginalDepositFieldPointer;


	/* 1) GravitatingMassField. */

	if (DepositField == GRAVITATING_MASS_FIELD) {
		if (TargetGrid->GravitatingMassFieldCellSize <= 0)
			TargetGrid->InitializeGravitatingMassField(RefineBy);
		/* by YS Jo, 0 for the original field; 1 for the gravity with stars */
#ifdef NBODY
		if (NoStar)
			DepositFieldPointer = TargetGrid->GravitatingMassFieldNoStar;
		else
			DepositFieldPointer = TargetGrid->GravitatingMassField;
#else
		DepositFieldPointer = TargetGrid->GravitatingMassField;
#endif
		CellSize            = TargetGrid->GravitatingMassFieldCellSize;
		CloudSize           = CellWidth[0][0];
		for (dim = 0; dim < GridRank; dim++) {
			LeftEdge[dim]  = TargetGrid->GravitatingMassFieldLeftEdge[dim];
			Dimension[dim] = TargetGrid->GravitatingMassFieldDimension[dim];
		}
	}

	/* 2) GravitatingMassFieldParticles. */

	else if (DepositField == GRAVITATING_MASS_FIELD_PARTICLES) {
		if (TargetGrid->GravitatingMassFieldParticlesCellSize <= 0)
			TargetGrid->InitializeGravitatingMassFieldParticles(RefineBy);
		/* by YS Jo, 0 for the original field; 1 for the gravity with stars */
#ifdef NBODY
		if (NoStar)
			DepositFieldPointer = TargetGrid->GravitatingMassFieldParticlesNoStar;
		else
			DepositFieldPointer = TargetGrid->GravitatingMassFieldParticles;
#else
		DepositFieldPointer = TargetGrid->GravitatingMassFieldParticles;
#endif
		CellSize            = TargetGrid->CellWidth[0][0];
		CloudSize            = CellWidth[0][0];
		for (dim = 0; dim < GridRank; dim++) {
			LeftEdge[dim]  = TargetGrid->GravitatingMassFieldParticlesLeftEdge[dim];
			Dimension[dim] = TargetGrid->GravitatingMassFieldParticlesDimension[dim];
		}
	}

	/* 3) MassFlaggingField */

	else if (DepositField == MASS_FLAGGING_FIELD) {
		DepositFieldPointer = TargetGrid->MassFlaggingField;
		CellSize            = TargetGrid->CellWidth[0][0];
		CloudSize           = CellWidth[0][0];
		for (dim = 0; dim < GridRank; dim++) {
			LeftEdge[dim]  = TargetGrid->CellLeftEdge[dim][0];
			Dimension[dim] = TargetGrid->GridDimension[dim];
		}
	}


	/* 4) ParticleMassFlaggingField */

	//  else if (DepositField == PARTICLE_MASS_FLAGGING_FIELD) {
	//    DepositFieldPointer = TargetGrid->ParticleMassFlaggingField;
	//    CellSize            = float(TargetGrid->CellWidth[0][0]);
	//    for (dim = 0; dim < GridRank; dim++) {
	//      LeftEdge[dim]  = TargetGrid->CellLeftEdge[dim][0];
	//      Dimension[dim] = TargetGrid->GridDimension[dim];
	//    }
	//  }

	/* 5) error */

	else {
		ENZO_VFAIL("DepositField = %"ISYM" not recognized.\n", DepositField)
	}  

	/* If on different processors, generate a temporary field to hold
		 the density. */

	if (ProcessorNumber != TargetGrid->ProcessorNumber) {
		/* If this is the target grid processor, then record the orginal
			 field characteristics so we can add it in when the data arrives. */

		for (dim = 0; dim < GridRank; dim++) {
			OriginalLeftEdge[dim] = LeftEdge[dim];
			OriginalDimension[dim] = Dimension[dim];
		}
		OriginalDepositFieldPointer = DepositFieldPointer;


		/* Resize the deposit region so it is just big enough to contain the
			 grid where the particles reside. */

		size = 1;
		for (dim = 0; dim < GridRank; dim++) {
			LeftEdge[dim] = (long_int((FLOAT)GridLeftEdge[dim]/CellSize)-2)*CellSize;
			Offset[dim] = nlongint((LeftEdge[dim] - OriginalLeftEdge[dim])/CellSize);
			if (Offset[dim] < 0) {
				fprintf(stderr, "P(%d)(1): dx=%"GOUTSYM"/%"GOUTSYM" = %"GOUTSYM"\n",
						MyProcessorNumber, CellSize, CellWidth[0][0], 
						CellSize/CellWidth[0][0]);
				fprintf(stderr, "P(%d)(2): %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
						MyProcessorNumber, OriginalLeftEdge[0], OriginalLeftEdge[1], 
						OriginalLeftEdge[2]);
				fprintf(stderr, "P(%d)(3): %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
						MyProcessorNumber, GridLeftEdge[0], GridLeftEdge[1], 
						GridLeftEdge[2]);
				fprintf(stderr, "P(%d)(4): %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
						MyProcessorNumber, GridRightEdge[0], GridRightEdge[1], 
						GridRightEdge[2]);
				fprintf(stderr, "P(%d)(5): %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
						MyProcessorNumber, LeftEdge[0], LeftEdge[1], LeftEdge[2]);
				fprintf(stderr, "P(%d)(6): %ld %ld %ld - %ld %ld %ld\n",
						MyProcessorNumber, Offset[0], Offset[1], Offset[2],
						Dimension[0], Dimension[1], Dimension[2]);
				fprintf(stderr, "P(%d)(7): %"GOUTSYM" %ld\n",
						MyProcessorNumber, (int(GridLeftEdge[dim]/CellSize)-2)*CellSize, 
						int(GridLeftEdge[dim]/CellSize));

				ENZO_VFAIL("Offset[%d] = %d < 0\n", dim, Offset[dim])
			}
			Dimension[dim] = int((GridRightEdge[dim] - LeftEdge[dim])/CellSize) + 3;
			size *= Dimension[dim];
		}

		/* Allocate buffer to communicate deposit region (unless in receive-mode
			 in which case the buffer was already allocated in post-receive mode). */

#ifdef USE_MPI
		if (CommunicationDirection == COMMUNICATION_RECEIVE) {
			DepositFieldPointer = 
				CommunicationReceiveBuffer[CommunicationReceiveIndex];
		} else {
			DepositFieldPointer = new float[size];
			if (MyProcessorNumber == ProcessorNumber) {
				for (i = 0; i < size; i++) {
					DepositFieldPointer[i] = 0.0;
				}
			}
		}
#endif /* USE_MPI */

	} // ENDIF different processors

	fprintf(stdout,"\nProc:%d 4-2-1\n", MyProcessorNumber); // by YS
	if (MyProcessorNumber == ProcessorNumber) {

		/* If using CIC-mode deposit, then set cloudsize equal to cellsize. */

		if (ParticleSubgridDepositMode == CIC_DEPOSIT)
			CloudSize = CellSize;

		/* If the Target is this grid and the DepositField is MassFlaggingField,
			 then multiply the Particle density by the volume to get the mass. */

		if (this == TargetGrid && DepositField == MASS_FLAGGING_FIELD)
			for (dim = 0; dim < GridRank; dim++)
				MassFactor *= CellWidth[dim][0];

		/* If the DepositGrid and this grid are not the same, we must adjust the
			 particle 'mass'. */

		if (this != TargetGrid) {

			/* Find the difference in resolution between this grid and TargetGrid. */

			float RefinementFactors[MAX_DIMENSION];
			this->ComputeRefinementFactorsFloat(TargetGrid, RefinementFactors);

			/* Compute the implied difference in 'mass' between particles in this
				 grid and those in TargetGrid. */

			for (dim = 0; dim < GridRank; dim++)
				MassFactor *= RefinementFactors[dim];

		} // ENDIF (this != TargetGrid)

		/* Check if we are smoothing. */

		int SmoothField = (DepositPositionsParticleSmoothRadius <= CellSize) ? FALSE : TRUE;

		/* If required, Change the mass of particles in this grid. */

#ifdef NBODY
		if (NoStar) {
			ParticleMassTemp = new float[NumberOfParticles];
			float MassFactorTemp = 1.;

			if (MassFactor != 1.0 ||
					((StarParticleCreation == (1 << SINK_PARTICLE)) &&
					 SmoothField == TRUE)) 
				MassFactorTemp = MassFactor;

			for (i = 0; i < NumberOfParticles; i++) {
				if ((ParticleType[i] == PARTICLE_TYPE_NBODY) ||
						(ParticleType[i] == PARTICLE_TYPE_NBODY_NEW))
					ParticleMassTemp[i] = 0;
				else
					ParticleMassTemp[i] = ParticleMass[i]*MassFactorTemp;
			}
			ParticleMassPointer = ParticleMassTemp;
		} else { // NoStar
			if (MassFactor != 1.0 ||
					((StarParticleCreation == (1 << SINK_PARTICLE)) &&
					 SmoothField == TRUE)) {
				ParticleMassTemp = new float[NumberOfParticles];

				for (i = 0; i < NumberOfParticles; i++)
					ParticleMassTemp[i] = ParticleMass[i]*MassFactor;
				ParticleMassPointer = ParticleMassTemp;
			}
			else
				ParticleMassPointer = ParticleMass;
		} // endif NoStar

		/*
		// by YS debug
		ParticleMassTemp = new float[NumberOfParticles];

		for (i = 0; i < NumberOfParticles; i++)
			ParticleMassTemp[i] = ParticleMass[i]*MassFactor;
		ParticleMassPointer = ParticleMassTemp;
		*/
#else
		if (MassFactor != 1.0 ||
				((StarParticleCreation == (1 << SINK_PARTICLE)) &&
				 SmoothField == TRUE)) {
			ParticleMassTemp = new float[NumberOfParticles];

			for (i = 0; i < NumberOfParticles; i++)
				ParticleMassTemp[i] = ParticleMass[i]*MassFactor;
			ParticleMassPointer = ParticleMassTemp;
		}
		else
			ParticleMassPointer = ParticleMass;
#endif



		/* If the target field is MASS_FLAGGING_FIELD, then set masses of
			 particles which are too large to zero (to prevent run-away
			 refinement). */

		if (DepositField == MASS_FLAGGING_FIELD &&
				DepositParticleMaximumParticleMass > 0 && MassFactor != 1.0)
			for (i = 0; i < NumberOfParticles; i++)
				ParticleMassPointer[i] = min(DepositParticleMaximumParticleMass,
						ParticleMassPointer[i]);

		/* Compute difference between current time and DepositTime. */

		TimeDifference = DepositTime - Time;

		/* Move particles to positions at Time + TimeDifference. */

		//The second argument forces the update even if 
		//MyProcessor == Target->ProcessorNumber != this->ProcessorNumber
		this->UpdateParticlePosition(TimeDifference, TRUE);

		/*
			 Right now all active particles are unsmoothed.
			 Will need to come back to this to optionally add smoothing.
			 */
		float FCellSize = (float)CellSize;
		float FCloudSize = (float)CloudSize;
		/* If using sink particles, then create a second field of unsmoothed sink particles
			 (since we don't want sink particles smoothed -- they are stellar sized). */

		/* Note that several types of particles may be appropriate for this,
			 but they will have to be added if needed. */
		if ((this->ReturnNumberOfStarParticles() > 0) && 
				(StarParticleCreation == (1 << SINK_PARTICLE)) && SmoothField == TRUE) {
			ParticleMassPointerSink = new float[NumberOfParticles];
			for (i = 0; i < NumberOfParticles; i++) {
				if (ParticleType[i] == PARTICLE_TYPE_STAR) {
					ParticleMassPointerSink[i] = ParticleMassPointer[i];
					ParticleMassPointer[i] = 0;
				} else {
					ParticleMassPointerSink[i] = 0;
				}
			}

			/* Deposit sink particles (only) to field using CIC or NGP. 
				 (only use NGP if cellsize > cloudsize - i.e. source is subgrid) */

			if (ParticleSubgridDepositMode == NGP_DEPOSIT && CellSize > 1.5*CloudSize) {
				PFORTRAN_NAME(ngp_deposit)(
						ParticlePosition[0], ParticlePosition[1], ParticlePosition[2], 
						&GridRank, &NumberOfParticles, ParticleMassPointerSink, DepositFieldPointer, 
						LeftEdge, Dimension, Dimension+1, Dimension+2, &FCellSize);
			} else {
				PFORTRAN_NAME(cic_deposit)(
						ParticlePosition[0], ParticlePosition[1], ParticlePosition[2], 
						&GridRank, &NumberOfParticles, ParticleMassPointerSink, DepositFieldPointer, 
						LeftEdge, Dimension, Dimension+1, Dimension+2, &FCellSize, &FCloudSize);
			}

			delete [] ParticleMassPointerSink;

		}

		/* Deposit particles. */

		if (SmoothField == FALSE) {

			//  fprintf(stderr, "------DP Call Fortran cic_deposit with CellSize = %"GSYM"\n", CellSize);

			/* Deposit sink particles (only) to field using CIC or NGP. 
				 (only use NGP if cellsize > cloudsize - i.e. source is subgrid) */

			if (ParticleSubgridDepositMode == NGP_DEPOSIT && CellSize > 1.5*CloudSize) {
				PFORTRAN_NAME(ngp_deposit)
					(ParticlePosition[0], ParticlePosition[1], ParticlePosition[2], 
					 &GridRank, &NumberOfParticles, ParticleMassPointer, DepositFieldPointer, 
					 LeftEdge, Dimension, Dimension+1, Dimension+2, &FCellSize);
			} else {
				PFORTRAN_NAME(cic_deposit)
					(ParticlePosition[0], ParticlePosition[1], ParticlePosition[2], 
					 &GridRank, &NumberOfParticles, ParticleMassPointer, DepositFieldPointer, 
					 LeftEdge, Dimension, Dimension+1, Dimension+2, &FCellSize, &FCloudSize);
			}

		} else {

			/* Deposit to field using large-spherical CIC, with radius of
				 DepositPositionsParticleSmoothRadius */

			//  fprintf(stderr, "------DP Call Fortran smooth_deposit with DPPSmoothRadius = %"GSYM"\n", DepositPositionsParticleSmoothRadius);

			PFORTRAN_NAME(smooth_deposit)
				(ParticlePosition[0], ParticlePosition[1], ParticlePosition[2], &GridRank,
				 &NumberOfParticles, ParticleMassPointer, DepositFieldPointer, LeftEdge, 
				 Dimension, Dimension+1, Dimension+2, &FCellSize, 
				 &DepositPositionsParticleSmoothRadius);
		}


		if ((this->ReturnNumberOfStarParticles() > 0) && 
				(StarParticleCreation == (1 << SINK_PARTICLE)) && SmoothField == TRUE) {
			for (i = 0; i < NumberOfParticles; i++) {
				if (ParticleType[i] == PARTICLE_TYPE_STAR) {
					ParticleMassPointer[i] = ParticleMassPointerSink[i];
				}
			}
			delete [] ParticleMassPointerSink;
		}

		if (NumberOfActiveParticles > 0) {
			FLOAT** ActiveParticlePosition = new FLOAT*[GridRank];
			for (dim = 0; dim < GridRank; dim++)
				ActiveParticlePosition[dim] = new FLOAT[NumberOfActiveParticles];
			this->GetActiveParticlePosition(ActiveParticlePosition);

			float* ActiveParticleMassPointer = new float[NumberOfActiveParticles];
			for (i = 0; i < NumberOfActiveParticles; i++) {
				if ((MassFactor != 1.0) || (SmoothField == TRUE))
					ActiveParticleMassPointer[i] = ActiveParticles[i]->ReturnMass()*MassFactor;
				else
					ActiveParticleMassPointer[i] = ActiveParticles[i]->ReturnMass();

				if (DepositField == MASS_FLAGGING_FIELD &&
						DepositParticleMaximumParticleMass > 0 && MassFactor != 1.0)
					ActiveParticleMassPointer[i] = min(DepositParticleMaximumParticleMass,
							ParticleMassPointer[i]);
			}

			if (SmoothField == FALSE) {

				PFORTRAN_NAME(cic_deposit)(
						ActiveParticlePosition[0], ActiveParticlePosition[1], ActiveParticlePosition[2],
						&GridRank, &NumberOfActiveParticles, ActiveParticleMassPointer, DepositFieldPointer,
						LeftEdge, Dimension, Dimension+1, Dimension+2, &FCellSize, &FCloudSize);

			}
			else {

				PFORTRAN_NAME(smooth_deposit)(
						ActiveParticlePosition[0], ActiveParticlePosition[1], ActiveParticlePosition[2],
						&GridRank, &NumberOfActiveParticles, ActiveParticleMassPointer, DepositFieldPointer,
						LeftEdge, Dimension, Dimension+1, Dimension+2, &FCellSize,
						&DepositPositionsParticleSmoothRadius);
			}

			for (dim = 0; dim < GridRank; dim++)
				delete [] ActiveParticlePosition[dim];
			delete [] ActiveParticlePosition;
			delete [] ActiveParticleMassPointer;
		}

	} // ENDIF this processor


		fprintf(stdout,"\nProc:%d 4-2-4\n", MyProcessorNumber); // by YS

	/* If any girds are on different processors, copy deposited field back to the
		 target grid and add to the correct field. */
	if (ProcessorNumber != TargetGrid->ProcessorNumber) {

#ifdef USE_MPI

		/* If posting a receive, then record details of call. */

		if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
			CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
			CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = TargetGrid;
			CommunicationReceiveCallType[CommunicationReceiveIndex] = 3;
			CommunicationReceiveArgument[0][CommunicationReceiveIndex] = DepositTime;
			CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] =
				DepositField;
		}

		MPI_Status status;
		MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
		MPI_Arg Count = size;
		MPI_Arg Source = ProcessorNumber;
		MPI_Arg Dest = TargetGrid->ProcessorNumber;

		double time1 = ReturnWallTime();

		if (MyProcessorNumber == ProcessorNumber)  {
			CommunicationBufferedSend(DepositFieldPointer, Count, DataType, 
					Dest, MPI_SENDREGION_TAG, 
					enzo_comm, BUFFER_IN_PLACE);
		}
		if (MyProcessorNumber == TargetGrid->ProcessorNumber &&
				CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
			MPI_Recv(DepositFieldPointer, Count, DataType, Source, 
					MPI_SENDREGION_TAG, enzo_comm, &status);
		}

		if (MyProcessorNumber == TargetGrid->ProcessorNumber &&
				CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
			MPI_Irecv(DepositFieldPointer, Count, DataType, Source, 
					MPI_SENDREGION_TAG, enzo_comm, 
					CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
			CommunicationReceiveBuffer[CommunicationReceiveIndex] = 
				DepositFieldPointer;
			CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
				CommunicationReceiveCurrentDependsOn;
			CommunicationReceiveIndex++;
		}


		CommunicationTime += ReturnWallTime() - time1;

#endif /* USE_MPI */

		if (MyProcessorNumber == TargetGrid->ProcessorNumber &&
				CommunicationDirection != COMMUNICATION_POST_RECEIVE) {
			index1 = 0;
			for (k = 0; k < Dimension[2]; k++)
				for (j = 0; j < Dimension[1]; j++) {
					index2 = ((k+Offset[2])*OriginalDimension[1] + j + Offset[1])*
						OriginalDimension[0] + Offset[0];
					for (i = 0; i < Dimension[0]; i++) {
						OriginalDepositFieldPointer[index2++] +=
							DepositFieldPointer[index1++];
					}
				}

			delete [] DepositFieldPointer;
		} // end: if (MyProcessorNumber == TargetGrid->ProcessorNumber)

	} // end: If (ProcessorNumber != TargetGrid->ProcessorNumber)

	if (MyProcessorNumber == ProcessorNumber) {

		/* If necessary, delete the particle mass temporary. */

		if (MassFactor != 1.0)

			delete [] ParticleMassTemp;

		/* Return particles to positions at Time. */

		this->UpdateParticlePosition(-TimeDifference);

	}

	return SUCCESS;
}
