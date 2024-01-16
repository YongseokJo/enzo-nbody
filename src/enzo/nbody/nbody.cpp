#include <mpi.h>
#include <iostream>
//#include "defs.h"
#include "global.h"
#include "nbody.h"


using namespace std;

//Global Variables
int NNB; double global_time; //bool debug;
double dt_min=1e-20;
std::vector<int> LevelList;
//MPI_Comm  comm, inter_comm, nbody_comm;

int nbody(int MyProcessorNumber) {
	cout << "Staring Nbody+ ..." << endl;
	std::vector<Particle*> particle{};

	//comm        = com;
	//inter_comm  = inter_com;
	//nbody_comm  = nbody_com;

	global_time = 0.;
	//debug = true;

	InitialCommunication(particle);
	//Parser(argc, argv);

	//if (readData(particle) == FAIL)
		//fprintf(stderr, "Read Data Failed!\n");

	/***
		for (Particle* elem: particle) {
		std::cout << elem->Position[0] <<" ";
		}
		std::cout << std::endl;
	 ***/

	initializeParticle(particle);


	//createComputationChain(particle);

	for (Particle* elem: particle) {
		std::cout << elem->TimeStepIrr <<" ";
	}

	/*
	for (Particle* elem: particle) {
		fprintf(stdout, "PID=%d, TReg=%e, TIrr=%e\n", elem->getPID(),elem->TimeStepReg, elem->TimeStepIrr);
		fprintf(stdout, "%e, %e, %e, %e, %e, %e\n",
				elem->Force[0], elem->Force[1], elem->Force[2], elem->ForceDot[0], elem->ForceDot[1], elem->ForceDot[2]);
		fprintf(stdout, "%e, %e, %e, %e, %e, %e\n",
				elem->dFReg[0][0], elem->dFReg[0][1], elem->dFReg[0][2], elem->dFIrr[0][0], elem->dFIrr[0][1], elem->dFIrr[0][3]);
		fprintf(stdout, "%e, %lf, %lf, %lf, %e, %e, %e\n\n",
				elem->Mass, elem->Position[0], elem->Position[1], elem->Position[2], elem->Velocity[0], elem->Velocity[1], elem->Velocity[2]);

	}
	std::cout << std::endl;
	*/



	Evolve(particle);

	// Particle should be deleted at some point

	return 0;
}
