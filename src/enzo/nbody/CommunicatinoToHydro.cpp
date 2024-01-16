#include <mpi.h>
#include <iostream>
#include "global.h"
#include "defs.h"

Particle* FirstParticleInEnzo;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;

void InitializeParticle(int newNNB, Particle* newParticle, std::vector<Particle*> &particle);
int CommunicationInterBarrier();
int newNNB = 0;


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

	Particle* ptclPtr;
	Particle* ptcl = new Particle[NNB];
	for (int i=0; i<NNB; i++) {
		if (i<(NNB)-1)
			ptclPtr = &ptcl[i+1];
		else
			ptclPtr = nullptr;
		ptcl[i].setParticleInfo(PID, Mass, Position, Velocity, BackgroundAcceleration, i, ptclPtr);
		particle.push_back(&ptcl[i]);
	}
	std::cout << "Particle loaded." << std::endl;
	delete PID, Mass, Position, Velocity, BackgroundAcceleration;
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


	EnzoTimeStep = TimeStep*EnzoTime;

	// Update Existing Particles
	// need to update if the ids match between Enzo and Nbody
	// set the last next pointer array to null

	Particle* updatedNextPtr;
	updatedNextPtr = nullptr;

	// loop for PID, going backwards to update the NextParticle
	for (int i=NNB-1; i>=0; i--) {
		for (int j=0; j<NNB; j++) {
			if (PID[i] == particle[j]->PID) {
				particle[i]->setParticleInfo(PID, BackgroundAcceleration, i, updatedNextPtr);
				updatedNextPtr = particle[i];
				if (i==0)
					FirstParticleInEnzo = particle[j];
				continue;
			}
		}
	}


	Particle* ptclPtr;
	Particle *ptcl = new Particle[newNNB];
	// Update New Particles
	for (int i=0; i<newNNB; i++) {
		if (i<(newNNB-1)) {
			ptclPtr = &ptcl[(i+1)];
		} else {
			ptclPtr = nullptr;
		}
		ptcl[i].setParticleInfo(newPID, newMass, newPosition, newVelocity, newBackgroundAcceleration, -1, i, ptclPtr);
		//initialize current time
		particle.push_back(&ptcl[i]);
	}

	// This includes modification of regular force and irregular force
	InitializeParticle(newNNB, ptcl, particle); 

	delete PID, BackgroundAcceleration, newPID, newMass, newPosition, newVelocity, newBackgroundAcceleration;

}



int SendToEzno(std::vector<Particle*> &particle) {

	MPI_Request request;
	MPI_Status status;

	double *Position[Dim], *Velocity[Dim], *newPosition[Dim], *newVelocity[Dim];

	int i;
	Particle *ptcl;
	ptcl = FirstParticleInEnzo;
	
	for (int dim=0; dim<Dim; dim++) {
		Position[dim]    = new double[NNB-newNNB];
		Velocity[dim]    = new double[NNB-newNNB];
		newPosition[dim] = new double[newNNB];
		newVelocity[dim] = new double[newNNB];
	}
	// Construct arrays
	int count=0;
	while (ptcl) {
		if (count < NNB-newNNB)
			for (int dim=0; dim<Dim; dim++) {
				Position[count][dim] = ptcl->Position[dim];
				Velocity[count][dim] = ptcl->Velocity[dim];
			}
		else 
			for (int dim=0; dim<Dim; dim++) {
				newPosition[count-NNB][dim] = ptcl->Position[dim];
				newVelocity[count-NNB][dim] = ptcl->Velocity[dim];
			}
		count++;
		ptcl = ptcl->NextParticleInEnzo;
	}


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

	for (int dim=0; dim<Dim; dim++) {
		delete Position[dim], Velocity[dim], newPosition[dim], newVelocity[dim];
	}
}




