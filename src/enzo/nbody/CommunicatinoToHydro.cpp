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

	fprintf(stdout, "NBODY+: Waiting for Enzo to receive data...\n");

	CommunicationInterBarrier();
	fprintf(stdout, "NBODY+: Receiving data from Enzo...\n");
  MPI_Recv(&NNB, 1, MPI_INT, 0, 100, inter_comm, &status);
	if (NNB != 0) {
		PID = new int[NNB];
		MPI_Recv(PID, NNB, MPI_INT, 0, 250, inter_comm, &status);
		Mass = new double[NNB];
		MPI_Recv(Mass, NNB, MPI_DOUBLE, 0, 200, inter_comm, &status);

		for (int dim=0; dim<Dim; dim++) {
			Position[dim] = new double[NNB];
			Velocity[dim] = new double[NNB];
			MPI_Recv(Position[dim], NNB, MPI_DOUBLE, 0, 300, inter_comm, &status);
			MPI_Recv(Velocity[dim], NNB, MPI_DOUBLE, 0, 400, inter_comm, &status);
		}

		for (int dim=0; dim<Dim; dim++) {
			BackgroundAcceleration[dim] = new double[NNB];
			MPI_Recv(BackgroundAcceleration[dim], NNB, MPI_DOUBLE, 0, 500, inter_comm, &status);
		}
	}
  MPI_Recv(&TimeStep,      1, MPI_DOUBLE, 0,  600, inter_comm, &status);
  MPI_Recv(&TimeUnits,     1, MPI_DOUBLE, 0,  700, inter_comm, &status);
  MPI_Recv(&LengthUnits,   1, MPI_DOUBLE, 0,  800, inter_comm, &status);
  MPI_Recv(&MassUnits,     1, MPI_DOUBLE, 0,  900, inter_comm, &status);
  MPI_Recv(&VelocityUnits, 1, MPI_DOUBLE, 0, 1000, inter_comm, &status);
	CommunicationInterBarrier();
	std::cout << "Data received!\n" << std::endl;

	// Enzo to Nbody unit convertors
	EnzoMass         = MassUnits/Msun*mass_unit;
	EnzoLength       = LengthUnits/pc*position_unit;
	EnzoVelocity     = VelocityUnits/pc*yr*position_unit/time_unit;
	EnzoTime         = TimeUnits/yr*time_unit;
	EnzoAcceleration = LengthUnits*LengthUnits/TimeUnits/pc/pc*yr*position_unit*position_unit/time_unit;

	EnzoTimeStep     = TimeStep*EnzoTime;

	if (NNB != 0) {
		Particle* ptclPtr;
		Particle* ptcl = new Particle[NNB];
		for (int i=0; i<NNB; i++) {
			if (i == NNB-1)
				ptclPtr = nullptr;
			else
				ptclPtr = &ptcl[i+1];
			ptcl[i].setParticleInfo(PID, Mass, Position, Velocity, BackgroundAcceleration, ptclPtr, i);
			particle.push_back(&ptcl[i]);
		}
		std::cout << NNB << " particles loaded!" << std::endl;
		delete PID, Mass, Position, Velocity, BackgroundAcceleration;
	}
}


int ReceiveFromEzno(std::vector<Particle*> &particle) {

	int *PID, *newPID, newNNB;
	double *BackgroundAcceleration[Dim];
	double *newMass, *newPosition[Dim], *newVelocity[Dim], *newBackgroundAcceleration[Dim];
	double TimeStep;

	MPI_Request request;
	MPI_Status status;
	fprintf(stdout, "NBODY+: Waiting for Enzo to receive data...\n");

	CommunicationInterBarrier();
	fprintf(stdout, "NBODY+: Receiving data from Enzo...\n");
	// Existing Particle Information
	if (NNB != 0)	{
		PID = new int[NNB];
		MPI_Recv(PID, NNB, MPI_DOUBLE, 0, 25, inter_comm, &status);
		for (int dim=0; dim<Dim; dim++) {
			BackgroundAcceleration[dim] = new double[NNB];
			MPI_Recv(BackgroundAcceleration[dim], NNB, MPI_DOUBLE, 0, 50, inter_comm, &status);
		}
	}


	// New Particle Information
  MPI_Recv(&newNNB, 1, MPI_INT, 0, 100, inter_comm, &status);
	if (newNNB > 0) {
		newPID = new int[newNNB];
		MPI_Recv(newPID, NNB, MPI_INT, 0, 250, inter_comm, &status);
		newMass = new double[newNNB];
		MPI_Recv(newMass, NNB, MPI_DOUBLE, 0, 200, inter_comm, &status);

		for (int dim=0; dim<Dim; dim++) {
			newPosition[dim] = new double[newNNB];
			newVelocity[dim] = new double[newNNB];
			MPI_Recv(newPosition[dim], newNNB, MPI_DOUBLE, 0, 300, inter_comm, &status);
			MPI_Recv(newVelocity[dim], newNNB, MPI_DOUBLE, 0, 400, inter_comm, &status);
		}
		for (int dim=0; dim<Dim; dim++) {
			newBackgroundAcceleration[dim] = new double[newNNB];
			MPI_Recv(newBackgroundAcceleration[dim], newNNB, MPI_DOUBLE, 0, 500, inter_comm, &status);
		}
	}

	// Timestep
	MPI_Recv(&TimeStep,      1, MPI_DOUBLE, 0,  600, inter_comm, &status);
	CommunicationInterBarrier();

	fprintf(stdout, "NBODY+: Data trnsferred!\n");

	EnzoTimeStep = TimeStep*EnzoTime;

	// Update Existing Particles
	// need to update if the ids match between Enzo and Nbody
	// set the last next pointer array to null

	Particle* updatedNextPtr;
	updatedNextPtr = nullptr;

	if (NNB != 0) {
		// loop for PID, going backwards to update the NextParticle
		for (int i=NNB-1; i>=0; i--) {
			for (int j=0; j<NNB; j++) {
				if (PID[i] == particle[j]->PID) {
					particle[i]->setParticleInfo(PID, BackgroundAcceleration, updatedNextPtr, i);
					updatedNextPtr = particle[i];
					if (i==0)
						FirstParticleInEnzo = particle[j];
					continue;
				} //endif PID
			} //endfor j
		} //endfor i
	} //endif nnb


	if (newNNB > 0) {
		Particle* ptclPtr;
		Particle *ptcl = new Particle[newNNB];
		// Update New Particles
		for (int i=0; i<newNNB; i++) {
			if (i<(newNNB-1)) {
				ptclPtr = &ptcl[(i+1)];
			} else {
				ptclPtr = nullptr;
			}
			ptcl[i].setParticleInfo(newPID, newMass, newPosition, newVelocity,
					newBackgroundAcceleration, NormalStar+SingleParticle+NewParticle, ptclPtr, i);
			//initialize current time
			particle.push_back(&ptcl[i]);
		}

		NNB += newNNB;
		// This includes modification of regular force and irregular force
		InitializeParticle(newNNB, ptcl, particle);
	}

	if (NNB != 0)
		delete PID, BackgroundAcceleration; 
	if (newNNB != 0)
		delete newPID, newMass, newPosition, newVelocity, newBackgroundAcceleration;

}



int SendToEzno(std::vector<Particle*> &particle) {

	if (NNB == 0)
		return DONE;
	MPI_Request request;
	MPI_Status status;

	double *Position[Dim], *Velocity[Dim], *newPosition[Dim], *newVelocity[Dim];

	int i;
	Particle *ptcl;
	ptcl = FirstParticleInEnzo;
	
	for (int dim=0; dim<Dim; dim++) {
		Position[dim]    = new double[NNB-newNNB];
		Velocity[dim]    = new double[NNB-newNNB];

		if (newNNB > 0) {
			newPosition[dim] = new double[newNNB];
			newVelocity[dim] = new double[newNNB];
		}
	}
	// Construct arrays
	int count=0;
	while (ptcl != nullptr) {
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


	fprintf(stdout, "NBODY+: Waiting for Enzo to sent data...\n");
	CommunicationInterBarrier();
	fprintf(stdout, "NBODY+: Sending data to Enzo...\n");
	//fprintf(stderr,"NumberOfParticles=%d\n",NumberOfNbodyParticles);
	for (int dim=0; dim<Dim; dim++) {
		MPI_Send(Position[dim], NNB, MPI_DOUBLE, 0, 300, inter_comm);
		MPI_Send(Velocity[dim], NNB, MPI_DOUBLE, 0, 400, inter_comm);
	}


	//fprintf(stderr,"NewNumberOfParticles=%d\n",NumberOfNewNbodyParticles);
	if (newNNB > 0) 
		for (int dim=0; dim<Dim; dim++) {
			MPI_Send(newPosition[dim], NNB, MPI_DOUBLE, 0, 500, inter_comm);
			MPI_Send(newVelocity[dim], NNB, MPI_DOUBLE, 0, 600, inter_comm);
		}
	CommunicationInterBarrier();
	fprintf(stdout, "NBODY+: Data sent!\n");

	for (int dim=0; dim<Dim; dim++) {
		delete Position[dim], Velocity[dim];
		if (newNNB > 0) 
			delete	newPosition[dim], newVelocity[dim];
	}
}




