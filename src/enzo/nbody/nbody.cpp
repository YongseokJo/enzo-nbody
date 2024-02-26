#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
//#include "defs.h"
#include "global.h"
#include "nbody.h"
#include "cuda/cuda_functions.h"


using namespace std;

//Global Variables
int NNB, newNNB; double global_time; //bool debug;
int NumNeighborMax = 100;
std::vector<int> RegIndexList; 
//MPI_Comm  comm, inter_comm, nbody_comm;
double EnzoTimeStep;
const double dt_min = 0.03125;
const int dt_level_min = -5;

int nbody(int MyProcessorNumber) {
	std::vector<Particle*> particle{};
	global_time = 0.;
	int irank = 0;

	// Generate a unique filename for each process
	std::ostringstream filenameStream;
	filenameStream << "nbody_output";
	std::string filename = filenameStream.str();

	// Open a file for each process
	std::ofstream outFile(filename);

	// Redirect cout to the output file
	std::streambuf* coutBuffer = std::cout.rdbuf();
	std::cout.rdbuf(outFile.rdbuf());


	std::cout << "Staring Nbody+ ..." << endl;

	InitializeDevice(&irank);
	InitialCommunication(particle);
	InitializeParticle(particle);
	Evolve(particle);

	// Particle should be deleted at some point

	// Close the file and restore cout
	outFile.close();
	std::cout.rdbuf(coutBuffer);

	return true;
}
