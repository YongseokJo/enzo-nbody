#include <mpi.h>
#include <iostream>
#include "global.h"
#include "defs.h"
#include "Particle/Particle.h"

double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
double EnzoTimeStep;

int CommunicationInterBarrier();

int InitialCommunication(std::vector<Particle*> &particle) {

	int *PID;
	double *Mass, *Position[Dim], *Velocity[Dim], *BackgroundAcceleration[Dim];
	double TimeStep, TimeUnits, LengthUnits, VelocityUnits, MassUnits;

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
		MPI_Recv(&Position[dim], NNB, MPI_DOUBLE, 0, 300, inter_comm, &status);
		MPI_Recv(&Velocity[dim], NNB, MPI_DOUBLE, 0, 400, inter_comm, &status);
	}

	for (int dim=0; dim<Dim; dim++) {
		BackgroundAcceleration[dim] = new double[NNB];
		MPI_Recv(&BackgroundAcceleration[dim], NNB, MPI_DOUBLE, 0, 500, inter_comm, &status);
	}
  MPI_Recv(&TimeStep,      1, MPI_DOUBLE, 0,  600, inter_comm, &status);
  MPI_Recv(&TimeUnits,     1, MPI_DOUBLE, 0,  700, inter_comm, &status);
  MPI_Recv(&LengthUnits,   1, MPI_DOUBLE, 0,  800, inter_comm, &status);
  MPI_Recv(&MassUnits,     1, MPI_DOUBLE, 0,  900, inter_comm, &status);
  MPI_Recv(&VelocityUnits, 1, MPI_DOUBLE, 0, 1000, inter_comm, &status);

	// Enzo to Nbody unit convertors
	EnzoMass         = MassUnits/Msun*mass_unit;
	EnzoLength       = LengthUnits/pc*position_unit;
	EnzoVelocity     = VelocityUnits/pc*yr*position_unit/time_unit;
	EnzoTime         = TimeUnits/yr*time_unit;
	EnzoAcceleration = LengthUnits*LengthUnits/TimeUnits/pc/pc*yr*position_unit*position_unit/time_unit;

	EnzoTimeStep     = TimeStep*EnzoTime;

	for (int i=0; i<NNB; i++)
		particle[i]->setParticleInfo(PID, Mass, Position, Velocity, BackgroundAcceleration, i);

	std::cout << "Particle loaded." << std::endl;

}


int ReceiveFromEzno(std::vector<Particle*> &particle) {

	int *PID, *newPID, newNNB;
	double *BackgroundAcceleration[Dim];
	double *newMass, *newPosition[Dim], *newVelocity[Dim], *newBackgroundAcceleration[Dim];
	double TimeStep;

	MPI_Request request;
	MPI_Status status;
	fprintf(stdout, "NBODY+: Receiving data from Enzo...\n");

	// Existing Particle Information
	PID = new int[NNB];
	CommunicationInterBarrier();
	MPI_Recv(&PID, NNB, MPI_DOUBLE, 0, 25, inter_comm, &status);
	for (int dim=0; dim<Dim; dim++) {
		BackgroundAcceleration[dim] = new double[NNB];
		MPI_Recv(&BackgroundAcceleration, NNB, MPI_DOUBLE, 0, 50, inter_comm, &status);
	}


	// New Particle Information
  MPI_Recv(&newNNB, 1, MPI_INT, 0, 100, inter_comm, &status);
	newPID = new int[newNNB];
  MPI_Recv(&newPID, NNB, MPI_INT, 0, 250, inter_comm, &status);
	newMass = new double[newNNB];
  MPI_Recv(&newMass, NNB, MPI_DOUBLE, 0, 200, inter_comm, &status);

	for (int dim=0; dim<Dim; dim++) {
		newPosition[dim] = new double[newNNB];
		newVelocity[dim] = new double[newNNB];
		MPI_Recv(&newPosition[dim], newNNB, MPI_DOUBLE, 0, 300, inter_comm, &status);
		MPI_Recv(&newVelocity[dim], newNNB, MPI_DOUBLE, 0, 400, inter_comm, &status);
	}
	for (int dim=0; dim<Dim; dim++) {
		newBackgroundAcceleration[dim] = new double[newNNB];
		MPI_Recv(&newBackgroundAcceleration[dim], newNNB, MPI_DOUBLE, 0, 500, inter_comm, &status);
	}

	MPI_Recv(&TimeStep,      1, MPI_DOUBLE, 0,  600, inter_comm, &status);
	CommunicationInterBarrier();


	EnzoTimeStep     = TimeStep*EnzoTime;

	// Update Existing Particles
	for (int i=0; i<NNB; i++)
		particle[i]->setParticleInfo(PID, BackgroundAcceleration, i);

	Particle *ptcl = new Particle[newNNB];
	// Update New Particles
	for (int i=0; i<newNNB; i++) {
		ptcl[i].setParticleInfo(newPID, newMass, newPosition, newVelocity, newBackgroundAcceleration, -1, ++i);
		particle.push_back(&ptcl[i]);
	}

	// Adjust accelerations according to new particles


}



int SendToEzno(std::vector<Particle*> &particle) {

	MPI_Request request;
	MPI_Status status;

	double *Position[Dim], *Velocity[Dim], *newPosition[Dim], *newVelocity[Dim];

	// I have to maintain the order
	//for 
	CommunicationInterBarrier();
	//fprintf(stderr,"NumberOfParticles=%d\n",NumberOfNbodyParticles);
	for (int dim=0; dim<Dim; dim++) {
		MPI_Send(Position[dim], NNB, MPI_DOUBLE, 0, 300, inter_comm);
		MPI_Send(Velocity[dim], NNB, MPI_DOUBLE, 0, 400, inter_comm);
	}


	//fprintf(stderr,"NewNumberOfParticles=%d\n",NumberOfNewNbodyParticles);
	if (NNB > 0) 
		for (int dim=0; dim<Dim; dim++) {
			MPI_Send(newPosition[dim], NNB, MPI_DOUBLE, 0, 500, inter_comm);
			MPI_Send(newVelocity[dim], NNB, MPI_DOUBLE, 0, 600, inter_comm);
		}
	CommunicationInterBarrier();

}




