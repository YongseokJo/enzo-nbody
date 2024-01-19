#include <mpi.h>
#include <iostream>
//#include "defs.h"
#include "global.h"
#include "nbody.h"


using namespace std;

//Global Variables
int NNB; double global_time; //bool debug;
std::vector<int> LevelList;
//MPI_Comm  comm, inter_comm, nbody_comm;
double EnzoTimeStep;
const double dt_min = 0.03125;
const int dt_level_min = -5;

int nbody(int MyProcessorNumber) {
	cout << "Staring Nbody+ ..." << endl;
	std::vector<Particle*> particle{};

	global_time = 0.;
	//debug = true;

	InitialCommunication(particle);
	//Parser(argc, argv);

	//if (readData(particle) == FAIL)
		//fprintf(stderr, "Read Data Failed!\n");

	InitializeParticle(particle);

	Evolve(particle);

	// Particle should be deleted at some point

	return 0;
}
