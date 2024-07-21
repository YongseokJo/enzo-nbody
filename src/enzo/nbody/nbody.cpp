#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
//#include "defs.h"
#include "global.h"
#include "nbody.h"
#include "cuda/cuda_functions.h"



//Global Variables
int NNB, newNNB; double global_time; //bool debug;
std::vector<Particle*> RegularList;
std::vector<Binary*> BinaryList;
std::vector<Particle*> BinaryCandidateList;
ULL NextRegTimeBlock=0.;
//MPI_Comm  comm, inter_comm, nbody_comm;
double EnzoTimeStep;
FILE* binout;
FILE* nbpout;
FILE* gpuout;

int nbody(int MyProcessorNumber) {
	std::vector<Particle*> particle{};
	global_time = 0.;
	int irank = 0;


	/*********************************************************************
	 *  Configuration for outputing log files
	 *********************************************************************/
	binout = fopen("binary_output.txt", "w");
	nbpout = fopen("nbody_output.txt", "w");
	gpuout = fopen("cuda_output.txt", "w");
	//pfmout = fopen("performace_output.txt", "w");
	fprintf(nbpout, "Nbody Output Starts!\n");
	fprintf(binout, "Binary Output Starts!\n");
	fprintf(gpuout, "CUDA Output Starts!\n");
	//fprintf(pfmout, "Performance log!\n");
	fflush(nbpout);
	fflush(binout);
	fflush(gpuout);
	//fflush(pfmout);





	/*********************************************************************
	 *  Nbody Calculations
	 *********************************************************************/
	std::cout << "Staring Nbody+ ..." << endl;
	InitializeDevice(&irank);
	InitialCommunication(particle);
	InitializeParticle(particle);
	Evolve(particle);

	// Particle should be deleted at some point

	/*********************************************************************
	 *  Exit
	 *********************************************************************/
	fclose(binout);
	fclose(nbpout);
	fclose(gpuout);
	//fclose(pfmout);

	return true;
}
