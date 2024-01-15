#include <mpi.h>
#include <iostream>
#include "global.h"
#include "defs.h"
#include "Particle/Particle.h"


int InitialCommunication(MPI_Comm &inter_comm, std::vector<Particle*> &particle) {

	int *PID;
	double *Mass, *Position[Dim], *Velocity[Dim], *BackgroundAcceleration[Dim];
	double TimeStep, TimeUnits, LengthUnits, VelocityUnits;

	MPI_Request request;
	MPI_Status status;

	fprintf(stdout, "NBODY+: Receiving data from Enzo...\n");
  MPI_Recv(&NNB, 1, MPI_INT, 0, 100, inter_comm, &status);
	PID = new int[NNB];
  MPI_Recv(&PID, NNB, MPI_INT, 0, 250, inter_comm, &status);
	Mass = new double[NNB];
  MPI_Recv(&Mass, NNB, MPI_DOUBLE, 0, 200, inter_comm, &status);
	for (int dim=0; dim<Dim; dim++) {

		Position[dim] = new double[NNB];
		Velocity[dim] = new double[NNB];
		MPI_Recv(&Position, NNB, MPI_DOUBLE, 0, 300, inter_comm, &status);
		MPI_Recv(&Velocity, NNB, MPI_DOUBLE, 0, 400, inter_comm, &status);
	}

	for (int dim=0; dim<Dim; dim++) {
		BackgroundAcceleration[dim] = new double[NNB];
		MPI_Recv(&BackgroundAcceleration, NNB, MPI_DOUBLE, 0, 500, inter_comm, &status);
	}
  MPI_Recv(&TimeStep, 1, MPI_DOUBLE, 0, 600, inter_comm, &status);
  MPI_Recv(&TimeUnits, 1, MPI_DOUBLE, 0, 700, inter_comm, &status);
  MPI_Recv(&LengthUnits, 1, MPI_DOUBLE, 0, 800, inter_comm, &status);
  MPI_Recv(&VelocityUnits, 1, MPI_DOUBLE, 0, 900, inter_comm, &status);


	// Unit conversion needed
	for (int i=0; i<NNB; i++)
		particle[i]->setParticleInfo(PID, Mass, Position, Velocity, BackgroundAcceleration, i);

	std::cout << "Particle loaded." << std::endl;

}
