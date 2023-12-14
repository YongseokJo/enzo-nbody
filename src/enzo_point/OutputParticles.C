/***********************************************************************
	/
	/  FIND ALL NBODY PARTICLES OVER ALL PROCESSORS
	/
	/  written by: Yongseok Jo
	/  date:       November, 2022
	/
	/  PURPOSE: First synchronizes particle information in the normal and 
	/           nbody particles.  Then we make a global particle list, which
	/           simplifies Nbody calculations.
	/
 ************************************************************************/

#include "mpi.h"
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "NbodyRoutines.h"  //added  
#include "phys_constants.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		HierarchyEntry **Grids[]);


#define GravTest
#ifdef GravTest
int OutputParitcles(LevelHierarchyEntry *LevelArray[], int level)
{

	int i, GridNum, LocalNumberOfNbodyParticles;
	LevelHierarchyEntry *Temp;
	int start_index;
	int length = 100000;
	int numpart = 0;
	float *pos[MAX_DIMENSION];
	float *acc[MAX_DIMENSION];

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		pos[dim] = new float[numpart];
		acc[dim] = new float[numpart];
	}
	//fprintf(stderr,"All Good 1\n");

	/* Store AccelerationNoStar of particles to sustainable ParticleAttribute */
	if (LevelArray[level+1] == NULL) {
		for (int l=1; l<level+1; level++)
			for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
				Temp->GridData->ReturnParticles();
	}


	numall = FindTotalNumberOfNbodyParticles(LevelArray);
	start_index = FindStartIndex(&numpart);

	float *PosAll[MAX_DIMENSION];
	float *AccAll[MAX_DIMENSION];

	for (int dim=0; dim<MAX_DIMENSION; dim++) {
		PosAll[dim] = new float[numpart];
		AccAll[dim] = new float[numpart];
	}

	if (MyProcessorNumber == ROOT_PROCESSOR) {
				/* Receiving Index, NumberOfParticles, NbodyArrays from other processs */
				int* start_index_all;
				int* localnumber;
				start_index_all = new int[NumberOfProcessors];
				localnumber = new int[NumberOfProcessors];
				MPI_Request request;
				MPI_Status status;
				int ierr;
				int errclass,resultlen;
				char err_buffer[MPI_MAX_ERROR_STRING];

				MPI_Gather(&numpart, 1, IntDataType, NumberAll, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index, 1, IntDataType, start_index_all, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);

				/*-----------------------------------------------*/
				/******  Gather Arrays from other processes  *****/
				/*-----------------------------------------------*/

				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					MPI_Igatherv(pos[dim], numpar, MPI_DOUBLE,
							PosAll[dim], NumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
					MPI_Wait(&request, &status);

					MPI_Igatherv(acc[dim], numpart, MPI_DOUBLE,
							NbodyParticleVelocity[dim], LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
					MPI_Wait(&request, &status);
				}

				for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
					MPI_Igatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
							NbodyParticleAccelerationNoStar[dim], LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
					MPI_Wait(&request, &status);
				}
				fprintf(stderr,"Done?1-2\n");



				fprintf(stdout, "First TimeStep: %e, %f\n", TimeStep, TimeUnits);

				fprintf(stderr,"Potential = %f,", NbodyParticleAccelerationNoStar[MAX_DIMENSION][0]);
				fprintf(stderr,"Potential = %f", NbodyParticleAccelerationNoStar[MAX_DIMENSION][1]);
				//fprintf(stderr, "mass:%e \n", NbodyParticleMass[0]);
				//fprintf(stderr, "x:%e \n", NbodyParticlePosition[0][0]);
				//fprintf(stderr, "vel:%e \n", NbodyParticleVelocity[0][0]);

				/*----------------------------------------------------------*/
				/******** Send Arrays to Fortran Nbody6++  First Time  *****/
				/*--------------------------------------------------------*/
				MPI_Send(&NumberOfNbodyParticles, 1, MPI_INT, 1, 100, inter_comm);
				MPI_Send(NbodyParticleMass, NumberOfNbodyParticles, MPI_DOUBLE, 1, 200, inter_comm);
				MPI_Send(NbodyParticleID, NumberOfNbodyParticles, MPI_INT, 1, 250, inter_comm);
				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					MPI_Send(NbodyParticlePosition[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 300, inter_comm);
					MPI_Send(NbodyParticleVelocity[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 400, inter_comm);
				}
				for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
					MPI_Send(NbodyParticleAccelerationNoStar[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 500, inter_comm);
				}
				fprintf(stderr,"Potential = %f,", NbodyParticleAccelerationNoStar[MAX_DIMENSION][0]);
				fprintf(stderr,"Potential = %f", NbodyParticleAccelerationNoStar[MAX_DIMENSION][1]);
				MPI_Send(&TimeStep, 1, MPI_DOUBLE, 1, 600, inter_comm);
				MPI_Send(&TimeUnits, 1, MPI_DOUBLE, 1, 700, inter_comm);
				MPI_Send(&LengthUnits, 1, MPI_DOUBLE, 1, 800, inter_comm);
				MPI_Send(&MassUnits, 1, MPI_DOUBLE, 1, 900, inter_comm);
				MPI_Send(&VelocityUnits, 1, MPI_DOUBLE, 1, 1000, inter_comm);

				if (start_index_all != NULL)
					delete [] start_index_all;
				start_index_all = NULL;
				if (LocalNumberAll != NULL)
					delete [] LocalNumberAll;
				LocalNumberAll = NULL;
				DeleteNbodyArrays();

			} // endif : root processor
			else {
				/* Sending Index, NumberOfParticles, NbodyArrays to the root processs */

				MPI_Request request;
				MPI_Status status;

				MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);


				MPI_Igatherv(NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
						NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm, &request);
				MPI_Wait(&request, &status);

				MPI_Igatherv(NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType,
						NULL, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm, &request);
				MPI_Wait(&request, &status);

				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					MPI_Igatherv(NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
							NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
					MPI_Wait(&request, &status);

					MPI_Igatherv(NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
							NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
					MPI_Wait(&request, &status);
				}

				for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
					MPI_Igatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
							NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
					MPI_Wait(&request, &status);
				}


			} // endelse: other processor
#endif

			NbodyFirst=FALSE;
			fprintf(stderr,"Proc %d: All Good 0, NbodyFirst=%d\n", MyProcessorNumber,NbodyFirst);
			fprintf(stderr,"Done?3\n");
			/* Destruct Arrays*/
			/*
			if (NbodyParticleIDTemp != NULL)
				delete [] NbodyParticleIDTemp;
			NbodyParticleIDTemp = NULL;
			*/

			if (NbodyParticleMassTemp != NULL)
				delete [] NbodyParticleMassTemp;
			NbodyParticleMassTemp = NULL;

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				if (NbodyParticlePositionTemp[dim] != NULL)
					delete [] NbodyParticlePositionTemp[dim];
				NbodyParticlePositionTemp[dim] = NULL;

				if (NbodyParticleVelocityTemp[dim] != NULL)
					delete [] NbodyParticleVelocityTemp[dim];
				NbodyParticleVelocityTemp[dim] = NULL;
			}

			for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
				if (NbodyParticleAccelerationNoStarTemp[dim] != NULL)
					delete [] NbodyParticleAccelerationNoStarTemp[dim];
				NbodyParticleAccelerationNoStarTemp[dim] = NULL;
			}
		}// endif: nbody first 
		else {
			//int *NbodyParticleIDTemp;
			//float *NbodyParticleMassTemp;
			//float *NbodyParticleVelocityTemp[MAX_DIMENSION]; // feedback can affect velocity
			float *NbodyParticleAccelerationNoStarTemp[MAX_DIMENSION+1];

			NbodyParticleIDTemp   = new int[LocalNumberOfNbodyParticles];
			//NbodyParticleMassTemp = new float[LocalNumberOfNbodyParticles];
			for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
				//NbodyParticleVelocityTemp[dim]            = new float[LocalNumberOfNbodyParticles];
				NbodyParticleAccelerationNoStarTemp[dim]  = new float[LocalNumberOfNbodyParticles];
			}

			/* Get particle information from Grids */
			int count = 0;
			for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++)
				for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel) 
					if (Temp->GridData->CopyNbodyParticles(&count, NbodyParticleIDTemp, NbodyParticleAccelerationNoStarTemp) == FAIL) {
						//NbodyParticleAccelerationTemp
						ENZO_FAIL("Error in grid::CopyNbodyParticles.");
					}

			if (count != LocalNumberOfNbodyParticles)
				ENZO_FAIL("Error in grid::CopyNbodyParticles.");

#ifdef USE_MPI
			if (MyProcessorNumber == ROOT_PROCESSOR) {
				/* Receiving Index, NumberOfParticles, NbodyArrays from other processs */
				int* start_index_all;
				int* LocalNumberAll;
				start_index_all = new int[NumberOfProcessors];
				LocalNumberAll = new int[NumberOfProcessors];
				MPI_Request request;
				MPI_Status status;
				int ierr;
				int errclass,resultlen;
				char err_buffer[MPI_MAX_ERROR_STRING];

				MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, LocalNumberAll, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index, 1, IntDataType, start_index_all, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);

				/* Initialize the nbody array used for direct Nbody calculation*/
				InitializeNbodyArrays();



				/*-----------------------------------------------*/
				/******  Gather Arrays from other processes  *****/
				/*-----------------------------------------------*/
				MPI_Igatherv(NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType,
						NbodyParticleID, LocalNumberAll, start_index_all, IntDataType, ROOT_PROCESSOR, enzo_comm,&request);
				MPI_Wait(&request, &status);
				for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
					MPI_Igatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
							NbodyParticleAccelerationNoStar[dim], LocalNumberAll, start_index_all, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
					MPI_Wait(&request, &status);
				}
				fprintf(stderr,"Done?1-2\n");



				fprintf(stdout, "TimeStep: %e, %f\n", TimeStep, TimeUnits);
				/*-----------------------------------------------*/
				/******** Send Arrays to Fortran Nbody6++    *****/
				/*-----------------------------------------------*/
				//MPI_Ssend(&NumberOfNbodyParticles, 1, MPI_INT, NumberOfProcessors, 100, inter_comm);
				CommunicationInterBarrier();
				MPI_Send(NbodyParticleID, NumberOfNbodyParticles, MPI_INT, 1, 250, inter_comm);
				for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
					ierr = MPI_Send(NbodyParticleAccelerationNoStar[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 500, inter_comm);
				}
				ierr = MPI_Send(&TimeStep, 1, MPI_DOUBLE, 1, 600, inter_comm);
				ierr = MPI_Send(&TimeUnits, 1, MPI_DOUBLE, 1, 700, inter_comm);
				ierr = CommunicationInterBarrier();

				/**
				if (ierr == MPI_SUCCESS) {
					fprintf(stderr,"Receiving from process %d with tag %d\n", status.MPI_SOURCE, status.MPI_TAG);
				} else {
					fprintf(stderr,"Error receiving message from process %d with tag %d\n", status.MPI_SOURCE, status.MPI_TAG);
					MPI_Error_class(ierr,&errclass);
					if (errclass== MPI_ERR_RANK) {
						fprintf(stderr,"Invalid rank used in MPI send call\n");
						MPI_Error_string(ierr,err_buffer,&resultlen);
						fprintf(stderr,err_buffer);
						MPI_Finalize();
					}
				}
				**/

				if (start_index_all != NULL)
					delete [] start_index_all;
				start_index_all = NULL;
				if (LocalNumberAll != NULL)
					delete [] LocalNumberAll;
				LocalNumberAll = NULL;
				DeleteNbodyArrays();
			} // endif : root processor
			else {
				/* Sending Index, NumberOfParticles, NbodyArrays to the root processs */

				MPI_Request request;
				MPI_Status status;

				MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);


				//MPI_Igatherv(NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE,
				//		NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm, &request);
				//MPI_Wait(&request, &status);
				MPI_Igatherv(NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType,
						NULL, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm, &request);
				MPI_Wait(&request, &status);

				for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
					//MPI_Igatherv(NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
					//			NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
					//	MPI_Wait(&request, &status);
					MPI_Igatherv(NbodyParticleAccelerationNoStarTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE,
							NULL, NULL, NULL, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,&request);
					MPI_Wait(&request, &status);
				}

			} // endelse: other processor
#endif

			/* Destruct Arrays*/
			/*
			if (NbodyParticleIDTemp != NULL)
				delete [] NbodyParticleIDTemp;
			NbodyParticleIDTemp = NULL;
			*/

			for (int dim=0; dim<MAX_DIMENSION+1; dim++) {
				if (NbodyParticleAccelerationNoStarTemp[dim] != NULL)
					delete [] NbodyParticleAccelerationNoStarTemp[dim];
				NbodyParticleAccelerationNoStarTemp[dim] = NULL;
			}
		} // endelse: not first
	} // ENDIF level
	return SUCCESS;

				int rows = 3;
				int cols = 3;
				const char *filename = "output.csv";

				FILE *csvFile = fopen(filename, "w");

				if (csvFile == NULL) {
					printf("Failed to open the CSV file.\n");
					return;
				}

				// Write the 2D array elements to the CSV file
				for (int i = 0; i < rows; i++) {
					for (int j = 0; j < cols; j++) {
						fprintf(csvFile, "%d", arr[i][j]);

						// Add a comma after each element except the last one in a row
						if (j < cols - 1) {
							fprintf(csvFile, ",");
						}
					}
					//                                                                                                                             // Add a newline character at the end of each row
					fprintf(csvFile, "\n");
				}

				// Close the CSV file
				fclose(csvFile);
			}


			write2DArrayToCSV(arr, rows, cols, filename);

			printf("2D array has been written to %s\n", filename);

}
#endif

