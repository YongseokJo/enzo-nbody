#include <mpi.h>
#include <iostream>
#include "global.h"
#include "defs.h"

Particle* FirstParticleInEnzo = nullptr;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
double EnzoCurrentTime, ClusterRadius2;

void InitializeParticle(Particle* newParticle, std::vector<Particle*> &particle);
int CommunicationInterBarrier();



int InitialCommunication(std::vector<Particle*> &particle) {

	int *PID;
	double *Mass, *Position[Dim], *Velocity[Dim], *BackgroundAcceleration[Dim];
	double *CreationTime, *DynamicalTime;
	double TimeStep, TimeUnits, LengthUnits, VelocityUnits, MassUnits;


	MPI_Request request;
	MPI_Status status;

	std::cout << "NBODY+: Waiting for Enzo to receive data...\n" << std::flush;

	CommunicationInterBarrier();
	std::cout << "NBODY+: Receiving data from Enzo...\n" << std::flush;
	MPI_Recv(&NNB, 1, MPI_INT, 0, 100, inter_comm, &status);
	std::cout << "NBODY+: NNB="<< NNB << std::endl;
	if (NNB != 0) {
		PID           = new int[NNB];
		Mass          = new double[NNB];
		CreationTime  = new double[NNB];
		DynamicalTime = new double[NNB];
		MPI_Recv(PID          , NNB, MPI_INT   , 0, 200, inter_comm, &status);
		MPI_Recv(Mass         , NNB, MPI_DOUBLE, 0, 201, inter_comm, &status);
		MPI_Recv(CreationTime , NNB, MPI_DOUBLE, 0, 202, inter_comm, &status);
		MPI_Recv(DynamicalTime, NNB, MPI_DOUBLE, 0, 203, inter_comm, &status);

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
	MPI_Recv(&TimeStep,                 1, MPI_DOUBLE, 0,  600, inter_comm, &status);
	MPI_Recv(&TimeUnits,                1, MPI_DOUBLE, 0,  700, inter_comm, &status);
	MPI_Recv(&LengthUnits,              1, MPI_DOUBLE, 0,  800, inter_comm, &status);
	MPI_Recv(&MassUnits,                1, MPI_DOUBLE, 0,  900, inter_comm, &status);
	MPI_Recv(&VelocityUnits,            1, MPI_DOUBLE, 0, 1000, inter_comm, &status);
	MPI_Recv(&StarParticleFeedback    , 1, MPI_INT   , 0, 1100, inter_comm, &status);
	MPI_Recv(&StarMassEjectionFraction, 1, MPI_DOUBLE, 0, 1200, inter_comm, &status);
	MPI_Recv(&EnzoCurrentTime         , 1, MPI_DOUBLE, 0, 1300, inter_comm, &status);
	MPI_Recv(&ClusterRadius2          , 1, MPI_DOUBLE, 0, 1400, inter_comm, &status);
	//MPI_Recv(&HydroMethod         , 1, MPI_INT   , 0, 1200, inter_comm, &status);
	CommunicationInterBarrier();
	std::cout << "Data received!\n" << std::endl;


	// Enzo to Nbody unit convertors
	EnzoMass         = MassUnits/Msun/mass_unit;
	EnzoLength       = LengthUnits/pc/position_unit;
	EnzoVelocity     = VelocityUnits/pc*yr/velocity_unit;
	EnzoTime         = TimeUnits/yr/time_unit;
	EnzoAcceleration = LengthUnits/TimeUnits/TimeUnits/pc*yr*yr/position_unit*time_unit*time_unit;

	EnzoTimeStep     = TimeStep*EnzoTime;
	EnzoCurrentTime  = EnzoCurrentTime*EnzoTime;
	ClusterRadius2  *= EnzoLength*EnzoLength;

	std::cout << "enzo Time:" << TimeStep << std::endl;
	std::cout << "nbody Time:" << EnzoTimeStep << std::endl;

	Particle* ptclPtr;
	Particle* ptcl = new Particle[NNB];
	if (NNB != 0) {
		for (int i=0; i<NNB; i++) {
			if (i == NNB-1)
				ptclPtr = nullptr;
			else
				ptclPtr = &ptcl[i+1];
			ptcl[i].setParticleInfo(PID, Mass, CreationTime, DynamicalTime, Position, Velocity,
					BackgroundAcceleration, ptclPtr, i);
			particle.push_back(&ptcl[i]);
		}

		FirstParticleInEnzo = &ptcl[0];
		std::cout << "enzo Pos     :" << Position[0][0] << ", " << Position[0][1] << std::endl;
		std::cout << "nbody Pos    :" << Position[0][0]*EnzoLength << ", " << Position[0][1]*EnzoLength << std::endl;
		std::cout << "enzo Vel     :" << Velocity[1][0] << ", " << Velocity[1][1] << std::endl;
		std::cout << "nbody Vel    :" << Velocity[1][0]*EnzoVelocity << ", " << Velocity[1][1]*EnzoVelocity << std::endl;
		std::cout << "km/s Vel     :" << Velocity[1][0]*VelocityUnits/1e5 << ", " << Velocity[1][1]*VelocityUnits/1e5 << std::endl;
		std::cout << "enzo  Mass   :" << Mass[0] << std::endl;
		std::cout << "nbody Mass   :" << Mass[0]*EnzoMass << std::endl;
		std::cout << "CreationTime :" << CreationTime[0] << std::endl;
		std::cout << "DynamicalTime:" << DynamicalTime[0] << std::endl;
		std::cout << "StarMassEjectionFraction:" << StarMassEjectionFraction << std::endl;
		std::cout << "StarParticleFeedback:" << StarParticleFeedback << std::endl;
		std::cout << "NBODY+       : "  << NNB << " particles loaded!" << std::endl;

		delete [] PID;
		delete [] Mass;
		delete [] CreationTime;
		delete [] DynamicalTime;
		for (int dim=0; dim<Dim; dim++) {
			delete[] Position[dim];
			delete[] Velocity[dim];
			delete[] BackgroundAcceleration[dim];
		}
	}
	return true;
}


int ReceiveFromEzno(std::vector<Particle*> &particle) {

	int *PID, *newPID;
	double *BackgroundAcceleration[Dim];
	double *Mass, *newMass, *newPosition[Dim], *newVelocity[Dim], *newBackgroundAcceleration[Dim];
	double *newCreationTime, *newDynamicalTime;
	double TimeStep;

	MPI_Request request;
	MPI_Status status;
	std::cout << "NBODY+: Waiting for Enzo to receive data..." << std::endl;

	fprintf(stdout, "NBODY+: Waiting for Enzo to receive data...\n");
	CommunicationInterBarrier();
	std::cout << "NBODY+: Receiving data from Enzo..." << std::endl;
	// Existing Particle Information
	if (NNB != 0)	{
		PID = new int[NNB];
		Mass = new double[NNB];
		MPI_Recv(PID , NNB, MPI_INT   , 0, 25, inter_comm, &status);
		MPI_Recv(Mass, NNB, MPI_DOUBLE, 0, 40, inter_comm, &status);
		for (int dim=0; dim<Dim; dim++) {
			BackgroundAcceleration[dim] = new double[NNB];
			MPI_Recv(BackgroundAcceleration[dim], NNB, MPI_DOUBLE, 0, 50, inter_comm, &status);
		}
	}


	// New Particle Information
	MPI_Recv(&newNNB, 1, MPI_INT, 0, 100, inter_comm, &status);
	if (newNNB > 0)
	{
		newPID           = new int[newNNB];
		newMass          = new double[newNNB];
		newCreationTime  = new double[newNNB];
		newDynamicalTime = new double[newNNB];
		MPI_Recv(newPID          , newNNB, MPI_INT   , 0, 200, inter_comm, &status);
		MPI_Recv(newMass         , newNNB, MPI_DOUBLE, 0, 201, inter_comm, &status);
		MPI_Recv(newCreationTime , newNNB, MPI_DOUBLE, 0, 202, inter_comm, &status);
		MPI_Recv(newDynamicalTime, newNNB, MPI_DOUBLE, 0, 203, inter_comm, &status);

		for (int dim=0; dim<Dim; dim++) {
			newPosition[dim]               = new double[newNNB];
			newVelocity[dim]               = new double[newNNB];
			newBackgroundAcceleration[dim] = new double[newNNB];
			MPI_Recv(newPosition[dim]              , newNNB, MPI_DOUBLE, 0, 300, inter_comm, &status);
			MPI_Recv(newVelocity[dim]              , newNNB, MPI_DOUBLE, 0, 400, inter_comm, &status);
			MPI_Recv(newBackgroundAcceleration[dim], newNNB, MPI_DOUBLE, 0, 500, inter_comm, &status);
		}
	}

	/***************************************************
	 * at some point we have to reconstruct *particle* vector because there is redundance due to particle removal.
	***************************************************/

	// Timestep
	MPI_Recv(&TimeStep       , 1, MPI_DOUBLE, 0, 600, inter_comm, &status);
	double OldEnzoCurrentTime = EnzoCurrentTime;
	MPI_Recv(&EnzoCurrentTime, 1, MPI_DOUBLE, 0, 700, inter_comm, &status);
	CommunicationInterBarrier();

	std::cout << "Enzo  Time    :" << EnzoCurrentTime << std::endl;
	//std::cout << "Nbody Time    :" << OldEnzoCurrentTime+particle[0]->CurrentTimeReg*EnzoTimeStep << std::endl;
	EnzoTimeStep    = TimeStep*EnzoTime;

	std::cout << "NBODY+: Data trnsferred!" << std::endl;
	fprintf(stdout, "NBODY+: Data trnsferred!\n");
	EnzoCurrentTime = EnzoCurrentTime*EnzoTime;
	std::cout << "Enzo  Time    :" << EnzoCurrentTime << std::endl;
	std::cout << "Enzo  TimeStep:" << TimeStep << std::endl;
	std::cout << "Enzo  TimeStep:" << EnzoTimeStep << std::endl;
	//std::cout << "Nbody Mass    :" << particle[0]->Mass << std::endl;
	//std::cout << "Enzo  Mass    :" << Mass[0]*EnzoMass << std::endl;
	//std::cout << "enzo Time :" << TimeStep << std::endl;
	//std::cout << "nbody Time:" << EnzoTimeStep << std::endl;

	// Update Existing Particles
	// need to update if the ids match between Enzo and Nbody
	// set the last next pointer array to null

	Particle *NextPtr, *LastPtr;
	NextPtr = nullptr;
	LastPtr = nullptr;

	std::cout << "NBODY+: 1" << std::endl;

	if (NNB != 0) {
		std::cout << "NBODY+: 2"  << std::endl;
		// loop for PID, going backwards to update the NextParticle
		for (int i=NNB-1; i>=0; i--) {
			for (int j=0; j<NNB; j++) {
				// PID of the received particle matches the PID of the existing particle
				if (PID[i] == particle[j]->PID) {
					particle[j]->setParticleInfo(Mass, BackgroundAcceleration, NextPtr, i);
					NextPtr = particle[j];
					fprintf(stdout, "NBODY+: pid= %d, x=%e\n",particle[j]->PID,particle[j]->Position[0]);
					if (i == 0) {
						FirstParticleInEnzo = particle[j];
					}
					if (i == NNB-1) {
						LastPtr = particle[j];
					}
					break;
				} //endif PID
			} //endfor j
		} //endfor i
	} //endif nnb

	//std::cout << "NBODY+: FirstParticleInEnzo PID= " << \
	FirstParticleInEnzo->PID << "in ReceiveFromEzno" << std::endl;

	if (newNNB > 0) {
		Particle* ptclPtr;
		Particle *newPtcl = new Particle[newNNB]; // we shouldn't delete it, maybe make it to vector

		// Update New Particles
		for (int i=0; i<newNNB; i++) {
			if (i<(newNNB-1)) {
				ptclPtr = &newPtcl[(i+1)];
			} else {
				ptclPtr = nullptr;
			}
			newPtcl[i].setParticleInfo(newPID, newMass, newCreationTime, newDynamicalTime, newPosition, newVelocity,
					newBackgroundAcceleration, NormalStar+SingleParticle+NewParticle, ptclPtr, i);
			//initialize current time
			particle.push_back(&newPtcl[i]);
		}

		if (FirstParticleInEnzo == nullptr) {
			FirstParticleInEnzo = particle[0];
		}
		if (LastPtr != nullptr) {
			LastPtr->NextParticleInEnzo = particle[0];
		}

		std::cout << "enzo Pos  :" << newPosition[0][0] << ", " << newPosition[0][1] << std::endl;
		std::cout << "nbody Pos :" << newPosition[0][0]*EnzoLength << ", " << newPosition[0][1]*EnzoLength << std::endl;
		std::cout << "enzo Vel  :" << newVelocity[1][0] << ", " << newVelocity[1][1] << std::endl;
		std::cout << "nbody Vel :" << newVelocity[1][0]*EnzoVelocity << ", " << newVelocity[1][1]*EnzoVelocity << std::endl;
		std::cout << "km/s Vel  :" << newVelocity[1][0]*velocity_unit/1e5/yr*pc << ", " << newVelocity[1][1]*velocity_unit/1e5/yr*pc << std::endl;
		std::cout << "enzo  Mass:" << newMass[0] << std::endl;
		std::cout << "nbody Mass:" << newMass[0]*EnzoMass << std::endl;
		std::cout << "NBODY+    : "  << newNNB << " new particles loaded!" << std::endl;

		// This includes modification of regular force and irregular force
		InitializeParticle(newPtcl, particle);
	}

	if (NNB != 0) {
		delete[] PID;
		delete[] Mass;
		for (int dim=0; dim<Dim; dim++) {
			delete[] BackgroundAcceleration[dim];
		}
	}
	if (newNNB != 0) {
		delete[] newPID;
		delete[] newMass;
		delete[] newDynamicalTime;
		delete[] newCreationTime;
		for (int dim=0; dim<Dim; dim++) {
			delete[] newBackgroundAcceleration[dim];
			delete[] newPosition[dim];
			delete[] newVelocity[dim];
		}
	}
	NNB += newNNB;
	return true;
}



int SendToEzno(std::vector<Particle*> &particle) {

	std::cout << "NBODY+: Entering SentToEnzo..." << std::endl;
	if (NNB == 0 && newNNB == 0) {
		std::cout << "NBODY+: Skipping SentToEnzo..." << std::endl;
		return DONE;
	}
	MPI_Request request;
	MPI_Status status;

	double TimeStep;
	double *Position[Dim], *Velocity[Dim], *newPosition[Dim], *newVelocity[Dim];

	Particle *ptcl;

	for (int dim=0; dim<Dim; dim++) {
		Position[dim]    = new double[NNB-newNNB];
		Velocity[dim]    = new double[NNB-newNNB];

		if (newNNB > 0) {
			newPosition[dim] = new double[newNNB];
			newVelocity[dim] = new double[newNNB];
		}
	}

	std::cout << "NBODY+: NNB=" << NNB << ", newNNB=" << newNNB << std::endl;
	std::cout << "NBODY+: size=" << particle.size() << std::endl;
	std::cout << "NBODY+: FirstParticleInEnzo PID=" << \
		FirstParticleInEnzo->PID << " in SendToEzno" << std::endl;

	// Construct arrays
	ptcl = FirstParticleInEnzo;
	if (ptcl == nullptr) {
		std::cout << "NBODY+: Warning! FirstParticleInEnzo is Null!" << std::endl;
	}

	int EscapeParticleNum = 0;
	double r2;
	int index;
	std::vector<Particle*> EscapeList;

	if (NNB - newNNB > 0) {
		for (int i=0; i<NNB-newNNB; i++) {
			r2 = 0;
			for (int dim=0; dim<Dim; dim++) {
				Position[dim][i] = ptcl->Position[dim]/EnzoLength;
				Velocity[dim][i] = ptcl->Velocity[dim]/EnzoVelocity;
				r2              += Position[dim][i]*Position[dim][i];
			}
			if (r2 > ClusterRadius2) { // in Enzo Unit
				Position[0][i] += 222;
				EscapeList.push_back(ptcl);
				EscapeParticleNum++;
			}
			//fprintf(stdout, "NBODY+: pid= %d, x=%e\n",ptcl->PID,Position[0][i]);
			ptcl = ptcl->NextParticleInEnzo;

			if ((ptcl == nullptr) && (i != NNB-1))
			{
				std::cout << "NBODY+: Warning! FirstParticleInEnzo is Null! No!" << std::endl;
			}
		}
	}

	if (newNNB > 0) {
		for (int i=0; i<newNNB; i++) {
			r2 = 0;
			for (int dim=0; dim<Dim; dim++) {
				newPosition[dim][i] = ptcl->Position[dim]/EnzoLength;
				newVelocity[dim][i] = ptcl->Velocity[dim]/EnzoVelocity;
				r2                 += newPosition[dim][i]*newPosition[dim][i];
			}
			if (r2 > ClusterRadius2) {
				newPosition[0][i] += 222;
				EscapeList.push_back(ptcl);
				EscapeParticleNum++;
			}
			ptcl = ptcl->NextParticleInEnzo;
		}
	}

	if (ptcl != nullptr) {
		std::cout << "NBODY+: Warrning! NextParticleInEnzo does not match!" << std::endl;
	}
	std::cout << "NBODY+: Waiting for Enzo to sent data..." << std::endl;
	fprintf(stdout, "NBODY+: Waiting for Enzo to sent data...\n");
	CommunicationInterBarrier();
	std::cout << "NBODY+: Sending data to Enzo..." << std::endl;
	if (NNB-newNNB != 0)
	{
		for (int dim = 0; dim < Dim; dim++)
		{
			MPI_Send(Position[dim], NNB - newNNB, MPI_DOUBLE, 0, 300, inter_comm);
			MPI_Send(Velocity[dim], NNB - newNNB, MPI_DOUBLE, 0, 400, inter_comm);
		}
	}

	//fprintf(stderr,"NewNumberOfParticles=%d\n",NumberOfNewNbodyParticles);
	if (newNNB > 0) {
		for (int dim=0; dim<Dim; dim++) {
			MPI_Send(newPosition[dim], newNNB, MPI_DOUBLE, 0, 500, inter_comm);
			MPI_Send(newVelocity[dim], newNNB, MPI_DOUBLE, 0, 600, inter_comm);
		}
	}
	CommunicationInterBarrier();
	std::cout << "NBODY+: Data sent!" << std::endl;
	fprintf(stdout, "NBODY+: Data sent!\n");


	for (Particle* elem2:EscapeList) {
		int i=0;
		for (Particle* elem1:particle) {
			if (elem1 == elem2) {
				particle.erase(particle.begin() + i);
				break;
			}
		}
	}


	// Initialization
	int i=0;
	RegIndexList.clear();
	for (Particle* elem:particle) {
		for (int dim=0; dim<Dim; dim++) {
			elem->PredPosition[dim] =  elem->Position[dim];
			elem->PredVelocity[dim] =  elem->Velocity[dim];
		}
		elem->isRegular = true;
		RegIndexList.push_back(i++);
		elem->CurrentTimeIrr = 0.;
		elem->CurrentTimeReg = 0.;
	}


	for (int dim = 0; dim < Dim; dim++)
	{
		delete[] Position[dim];
		delete[] Velocity[dim];
		if (newNNB > 0)
		{
			delete[] newPosition[dim];
			delete[] newVelocity[dim];
		}
	}

	std::cout << "NBODY+: Sending data finished." << std::endl;

	NNB -= EscapeParticleNum;

	return true;
}




