#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "global.h"
#include "defs.h"

Particle* FirstParticleInEnzo = nullptr;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
double EnzoCurrentTime, ClusterRadius2;
double ClusterAcceleration[Dim], ClusterPosition[Dim], ClusterVelocity[Dim];
double EPS2, eta, InitialRadiusOfAC;
int FixNumNeighbor, FixNumNeighbor0;
int BinaryRegularization;
double KSTime;
double KSDistance;

void InitializeNewParticle(std::vector<Particle*> &particle, int offset);
void GetCenterOfMass(double *mass, double *x[Dim], double *v[Dim], double x_com[], double v_com[], int N);
void GetNewCenterOfMass(std::vector<Particle*> &particle, double *mass2, double *x2[Dim], double *v2[Dim], int n2, double x_X[], double v_X[]);
void UpdateNextRegTime(std::vector<Particle*> &particle);
int CommunicationInterBarrier();


using namespace std;
const int width = 18;

int InitialCommunication(std::vector<Particle*> &particle) {

	int *PID;
	double *Mass, *Position[Dim], *Velocity[Dim], *BackgroundAcceleration[Dim];
	double *CreationTime, *DynamicalTime;
	double TimeStep, TimeUnits, LengthUnits, VelocityUnits, MassUnits;

	MPI_Request request;
	MPI_Status status;

	fprintf(nbpout, "NBODY+: First Waiting for Enzo to receive data...\n");

	CommunicationInterBarrier();
	fprintf(nbpout, "NBODY+: Frist Receiving data from Enzo...\n");
	MPI_Recv(&NNB, 1, MPI_INT, 0, 100, inter_comm, &status);
	fprintf(nbpout, "NBODY+: NNB=%d\n", NNB);
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
	MPI_Recv(&EPS2                    , 1, MPI_DOUBLE, 0, 1400, inter_comm, &status);
	MPI_Recv(&eta                     , 1, MPI_DOUBLE, 0, 1500, inter_comm, &status);
	MPI_Recv(&InitialRadiusOfAC       , 1, MPI_DOUBLE, 0, 1600, inter_comm, &status);
	MPI_Recv(&ClusterRadius2          , 1, MPI_DOUBLE, 0, 1700, inter_comm, &status);
	MPI_Recv(&FixNumNeighbor          , 1, MPI_INT   , 0, 1800, inter_comm, &status);
	MPI_Recv(&BinaryRegularization    , 1, MPI_INT   , 0, 1900, inter_comm, &status);
	MPI_Recv(&KSDistance              , 1, MPI_DOUBLE, 0, 2000, inter_comm, &status);
	MPI_Recv(&KSTime                  , 1, MPI_DOUBLE, 0, 2100, inter_comm, &status);
	//MPI_Recv(&HydroMethod         , 1, MPI_INT   , 0, 1200, inter_comm, &status);
	CommunicationInterBarrier();
	fprintf(nbpout, "Data received!\n");


	// Enzo to Nbody unit convertors
	EnzoMass         = MassUnits/Msun/mass_unit;
	EnzoLength       = LengthUnits/pc/position_unit;
	EnzoVelocity     = VelocityUnits/pc*yr/velocity_unit;
	EnzoTime         = TimeUnits/yr/time_unit;
	//EnzoAcceleration = LengthUnits/TimeUnits/TimeUnits/pc*yr*yr/position_unit*time_unit*time_unit;
	EnzoAcceleration = EnzoLength/EnzoTime/EnzoTime;

	// Unit conversion
	EnzoTimeStep       = TimeStep*EnzoTime;
	EnzoCurrentTime   *= EnzoTime;
	if (EPS2 < 0)
		EPS2 = -1;
	else {
		EPS2 *= EnzoLength;
		EPS2 *= EPS2;
	}

	InitialRadiusOfAC *= EnzoLength;
	FixNumNeighbor0    = FixNumNeighbor;

	fprintf(nbpout, "Enzo Time                = %lf\n", TimeStep);
	fprintf(nbpout, "Nbody Time               = %lf\n", EnzoTimeStep);
	fprintf(nbpout, "EPS2                     = %lf\n", EPS2);
	fprintf(nbpout, "InitialRadiusOfAC        = %.2e\n", InitialRadiusOfAC);
	fprintf(nbpout, "eta                      = %lf\n", eta);
	fprintf(nbpout, "ClusterRadius2           = %.2e\n", ClusterRadius2);
	fprintf(nbpout, "StarMassEjectionFraction = %lf\n", StarMassEjectionFraction);
	fprintf(nbpout, "StarParticleFeedback     = %d\n", StarParticleFeedback);
	fprintf(nbpout, "FixNumNeighbor           = %d\n", FixNumNeighbor);
	fprintf(nbpout, "BinaryRegularization     = %d\n", BinaryRegularization);
	fprintf(nbpout, "KSTime                   = %lf\n", KSTime);
	fprintf(nbpout, "KSDistance               = %lf\n", KSDistance);


	if (NNB != 0) {
		// set COM for background acc
		GetCenterOfMass(Mass, Position, Velocity, ClusterPosition, ClusterVelocity, NNB);

		ClusterAcceleration[0] = 0.;
		ClusterAcceleration[1] = 0.;
		ClusterAcceleration[2] = 0.;
		double total_mass=0.;
		for (int i=0; i<NNB; i++) {
			for (int dim=0; dim<Dim; dim++) {
				ClusterAcceleration[dim] += Mass[i]*BackgroundAcceleration[dim][i];
			}
			total_mass += Mass[i];
		}
		
		for (int dim=0; dim<Dim; dim++) 
			ClusterAcceleration[dim] /= total_mass;

		Particle* ptclPtr;

		for (int i=0; i<NNB; i++) {
			for (int dim=0; dim<Dim; dim++) {
				BackgroundAcceleration[dim][i] -= ClusterAcceleration[dim];
				Position[dim][i]               -= ClusterPosition[dim];
				Velocity[dim][i]               -= ClusterVelocity[dim];
			}

			particle.push_back(new Particle(PID, Mass, CreationTime, DynamicalTime, Position, Velocity,
						BackgroundAcceleration, ptclPtr, i));
		}

		fprintf(stderr, "Initial PID order= ");
		for (int i=0; i<NNB; i++) {
			if (i == NNB-1)
				ptclPtr = nullptr;
			else
				ptclPtr = particle[i+1];
			particle[i]->NextParticleInEnzo = ptclPtr;
			fprintf(stderr, "%d, ", particle[i]->PID);
		}
		FirstParticleInEnzo = particle[0];
		fprintf(stderr, "\n");

		/*
		std::cout << "enzo Pos     :" << Position[0][0] << ", " << Position[0][1] << std::endl;
		std::cout << "nbody Pos    :" << Position[0][0]*EnzoLength << ", " << Position[0][1]*EnzoLength << std::endl;
		std::cout << "enzo Vel     :" << Velocity[1][0] << ", " << Velocity[1][1] << std::endl;
		std::cout << "nbody Vel    :" << Velocity[1][0]*EnzoVelocity << ", " << Velocity[1][1]*EnzoVelocity << std::endl;
		std::cout << "km/s Vel     :" << Velocity[1][0]*VelocityUnits/1e5 << ", " << Velocity[1][1]*VelocityUnits/1e5 << std::endl;
		std::cout << "enzo  Mass   :" << Mass[0] << std::endl;
		std::cout << "nbody Mass   :" << Mass[0]*EnzoMass << std::endl;
		std::cout << "CreationTime :" << CreationTime[0] << std::endl;
		std::cout << "DynamicalTime:" << DynamicalTime[0] << std::endl;
		*/

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
	fprintf(nbpout, "NBODY+: %d particles loaded!\n", NNB);
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
	fprintf(nbpout, "NBODY+: Waiting for Enzo to receive data...\n");


	fprintf(stdout, "NBODY+: Waiting for Enzo to receive data...\n");
	CommunicationInterBarrier();
	fprintf(nbpout, "NBODY+: Receiving data from Enzo...\n");

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
	std::cout << "Enzo Time    :" << EnzoCurrentTime << std::endl;
	std::cout << "Enzo TimeStep:" << TimeStep << std::endl;
	std::cout << "EnzoTimeStep :" << EnzoTimeStep << std::endl;
	//std::cout << "Nbody Mass    :" << particle[0]->Mass << std::endl;
	//std::cout << "Enzo  Mass    :" << Mass[0]*EnzoMass << std::endl;
	//std::cout << "enzo Time :" << TimeStep << std::endl;
	//std::cout << "nbody Time:" << EnzoTimeStep << std::endl;



	// COM conversion
	// later on we might need to take mass weight into account 
	// 1. F=ma; 2. F -> F_com; 3. F_com -> a_com
	ClusterAcceleration[0] = 0.0;
	ClusterAcceleration[1] = 0.0;
	ClusterAcceleration[2] = 0.0;

	double total_mass = 0.;
	if (NNB != 0) {
		for (int i=0; i<NNB; i++) {
			for (int dim=0; dim<Dim; dim++) {
				ClusterAcceleration[dim] += Mass[i]*BackgroundAcceleration[dim][i];
			}
			total_mass += Mass[i];
		}
		std::cout << 
			std::setw(width)  << "NBODY : ClusterAcceleration = " << 
			std::setw(width)  << ClusterAcceleration[0] << 
			std::setw(width)  << ClusterAcceleration[1] << 
			std::setw(width)  << ClusterAcceleration[2] << 
			std::endl;

		std::cout << 
			std::setw(width)  << "NBODY : ClusterPosition = " << 
			std::setw(width)  << ClusterPosition[0] << 
			std::setw(width)  << ClusterPosition[1] << 
			std::setw(width)  << ClusterPosition[2] << 
			std::endl;	 	
		std::cout << 
			std::setw(width)  << "NBODY : ClusterVelocity = " << 
			std::setw(width)  << ClusterVelocity[0] << 
			std::setw(width)  << ClusterVelocity[1] << 
			std::setw(width)  << ClusterVelocity[2] << 
			std::endl;
	}



	if (newNNB != 0) {
		// we need to make adjustment to COM
		GetNewCenterOfMass(particle, newMass, newPosition, newVelocity, newNNB, 
				ClusterPosition, ClusterVelocity);

		for (int i=0; i<newNNB; i++) {
			for (int dim=0; dim<Dim; dim++) {
				ClusterAcceleration[dim] += newMass[i]*newBackgroundAcceleration[dim][i];
			}
			total_mass += newMass[i];
		}
	}

	if (total_mass != 0.) {
		for (int dim=0; dim<Dim; dim++) {
			ClusterAcceleration[dim] /= total_mass;
		}
	}



	// Update Existing Particles
	// need to update if the ids match between Enzo and Nbody
	// set the last next pointer array to null

	Particle *NextPtr, *LastPtr;
	NextPtr = nullptr;
	LastPtr = nullptr;
	if (NNB != 0) {
		// loop for PID, going backwards to update the NextParticle
		for (int i=NNB-1; i>=0; i--) {
			for (int dim=0; dim<Dim; dim++) {
				BackgroundAcceleration[dim][i] -= ClusterAcceleration[dim];
			}
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


	fprintf(stderr, "Receiving PID order= ");
	NextPtr = FirstParticleInEnzo;
	for (int i=0; i<NNB; i++) {
		fprintf(stderr, "%d, ", NextPtr->PID);
		NextPtr = NextPtr->NextParticleInEnzo;
	}
	fprintf(stderr, "\n");


	//std::cout << "NBODY+: FirstParticleInEnzo PID= " << \
	FirstParticleInEnzo->PID << "in ReceiveFromEzno" << std::endl;

	Particle* ptclPtr;
	//Particle *newPtcl = new Particle[newNNB]; // we shouldn't delete it, maybe make it to vector
	
	if (newNNB > 0) {
		/*
		for (int i = 0; i < newNNB; i++) {
			fprintf(stderr, "NBODY: Mass Of NewNbodyParticles=%.3e\n", newMass[i]);
			fprintf(stderr, "NBODY: Vel  Of news=(%.3e, %.3e, %.3e)\n", 
					newVelocity[0][i], newVelocity[1][i], newVelocity[2][i]);
			fprintf(stderr, "NBODY: Pos  Of news=(%.3e, %.3e, %.3e)\n",
					newPosition[0][i], newPosition[1][i], newPosition[2][i]);
			fprintf(stderr, "NBODY: Acc  Of news=(%.3e, %.3e, %.3e)\n",
					newBackgroundAcceleration[0][i], newBackgroundAcceleration[1][i], newBackgroundAcceleration[2][i]);
		}
		*/


		// Update New Particles
		for (int i=0; i<newNNB; i++) {
			for (int dim=0; dim<Dim; dim++) {
				newBackgroundAcceleration[dim][i] -= ClusterAcceleration[dim];
				newPosition[dim][i]               -= ClusterPosition[dim];
				newVelocity[dim][i]               -= ClusterVelocity[dim];
			}

			//initialize current time
			particle.push_back(new Particle(newPID, newMass, newCreationTime, newDynamicalTime, newPosition, newVelocity,
					newBackgroundAcceleration, NormalStar+SingleParticle+NewParticle, ptclPtr, i));
		}

		// build paticle chain
		for (int i=0; i<newNNB; i++) {
			if (i == newNNB-1)
				ptclPtr = nullptr;
			else 
				ptclPtr = particle[(i+1+NNB)];
			particle[NNB+i]->NextParticleInEnzo = ptclPtr;
		}

		// what are these?
		if (FirstParticleInEnzo == nullptr || NNB == 0) {
			FirstParticleInEnzo = particle[0];
		}
		// this connects nnb to newnnb
		if (LastPtr != nullptr) {
			LastPtr->NextParticleInEnzo = particle[NNB];
		}

		/*
		std::cout << "enzo Pos  :" << newPosition[0][0] << ", " << newPosition[0][1] << std::endl;
		std::cout << "nbody Pos :" << newPosition[0][0]*EnzoLength << ", " << newPosition[0][1]*EnzoLength << std::endl;
		std::cout << "enzo Vel  :" << newVelocity[1][0] << ", " << newVelocity[1][1] << std::endl;
		std::cout << "nbody Vel :" << newVelocity[1][0]*EnzoVelocity << ", " << newVelocity[1][1]*EnzoVelocity << std::endl;
		std::cout << "km/s Vel  :" << newVelocity[1][0]*velocity_unit/1e5/yr*pc << ", " << newVelocity[1][1]*velocity_unit/1e5/yr*pc << std::endl;
		std::cout << "enzo  Mass:" << newMass[0] << std::endl;
		std::cout << "nbody Mass:" << newMass[0]*EnzoMass << std::endl;
		*/
		std::cout << "NBODY+    : "  << newNNB << " new particles loaded!" << std::endl;

		// This includes modification of regular force and irregular force
		InitializeNewParticle(particle, NNB);
	}

	/*
	for (Particle* ptcl:particle) {
		fprintf(stderr, "NBODY2: Mass Of NewNbodyParticles=%.3e\n", ptcl->Mass);
		fprintf(stderr, "NBODY2: Vel  Of news=(%.3e, %.3e, %.3e)\n", 
				ptcl->Velocity[0], ptcl->Velocity[1], ptcl->Velocity[2]);
		fprintf(stderr, "NBODY2: Pos  Of news=(%.3e, %.3e, %.3e)\n",
				ptcl->Position[0], ptcl->Position[1], ptcl->Position[2]);
		fprintf(stderr, "NBODY2: Acc  Of news=(%.3e, %.3e, %.3e)\n",
				ptcl->a_tot[0][0], ptcl->a_tot[1][0], ptcl->a_tot[2][0]);
	}
	*/

	// Initialization
	int i=0;
	//RegularList.clear();
	for (Particle* elem:particle) {
		for (int dim=0; dim<Dim; dim++) {
			elem->PredPosition[dim] =  elem->Position[dim];
			elem->PredVelocity[dim] =  elem->Velocity[dim];
		}
		//RegularList.push_back(elem);
		elem->CurrentTimeIrr  = 0.;
		elem->CurrentBlockIrr = 0.;
		elem->CurrentTimeReg  = 0.;
		elem->CurrentBlockReg = 0.;
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

	FixNumNeighbor = std::min((int) std::floor(NNB/2), FixNumNeighbor0);


	fprintf(nbpout, "NBODY+    : In ReceiveFromEzno (after new particle might be added): \n");
	fprintf(nbpout, "NBODY+    : original NNB      = %d (+%d)\n", NNB-newNNB, newNNB);
	fprintf(nbpout, "NBODY+    : newly updated NNB = %d\n", NNB, newNNB);
	fprintf(nbpout, "NBODY+    : Particle size     = %d\n", particle.size());
	fprintf(nbpout, "NBODY+    : RegularList size = %d\n", RegularList.size());


	fprintf(stdout, "NBODY+    : In ReceiveFromEzno (after new particle might be added): \n");
	fprintf(stdout, "NBODY+    : original NNB      = %d (+%d)\n", NNB-newNNB, newNNB);
	fprintf(stdout, "NBODY+    : newly updated NNB = %d\n", NNB, newNNB);
	fprintf(stdout, "NBODY+    : Particle size     = %d\n", particle.size());
	fprintf(stdout, "NBODY+    : RegularList size = %d\n", RegularList.size());
	return true;
}



int SendToEzno(std::vector<Particle*> &particle) {

	std::cout << "NBODY+: Entering SendToEnzo..." << std::endl;
	if (NNB == 0 && newNNB == 0) {
		std::cout << "NBODY+: Skipping SendToEnzo..." << std::endl;
		return DONE;
	}
	MPI_Request request;
	MPI_Status status;

	double *Position[Dim], *Velocity[Dim], *newPosition[Dim], *newVelocity[Dim];

	Particle *ptcl;

	for (int dim=0; dim<Dim; dim++) {
		if (NNB-newNNB != 0) {
			Position[dim]    = new double[NNB-newNNB];
			Velocity[dim]    = new double[NNB-newNNB];
		}

		if (newNNB > 0) {
			newPosition[dim] = new double[newNNB];
			newVelocity[dim] = new double[newNNB];
		}
	}

	std::cout << "NBODY+: NNB=" << NNB << ", newNNB=" << newNNB << std::endl;
	std::cout << "NBODY+: size=" << particle.size() << std::endl;
	std::cout << "NBODY+: FirstParticleInEnzo PID=" << \
		FirstParticleInEnzo->PID << " in SendToEzno" << std::endl;



	int EscapeParticleNum = 0;
	double r2;
	int index;
	double TimeStep=EnzoTimeStep/EnzoTime;
	std::vector<Particle*> EscapeList;



	// COM evolution	
	for (int dim=0; dim<Dim; dim++) {	
		ClusterPosition[dim] += ClusterVelocity[dim]*TimeStep;
		ClusterPosition[dim] += ClusterAcceleration[dim]*TimeStep*TimeStep/2;
		ClusterVelocity[dim] += ClusterAcceleration[dim]*TimeStep;
	}

	ptcl = FirstParticleInEnzo;
	if (ptcl == nullptr) {
		std::cout << "NBODY+: Warning! FirstParticleInEnzo is Null!" << std::endl;
	}
	if (NNB - newNNB > 0) {
		fprintf(stderr, "Sending PID order= ");
		for (int i=0; i<NNB-newNNB; i++) {
			fprintf(stderr, "%d, ",ptcl->PID);
			r2 = 0;
			for (int dim=0; dim<Dim; dim++) {
				Position[dim][i]  = ptcl->Position[dim]/EnzoLength;
				Velocity[dim][i]  = ptcl->Velocity[dim]/EnzoVelocity;
				r2               += Position[dim][i]*Position[dim][i];

				// COM correction
				Position[dim][i] += ClusterPosition[dim];
				Velocity[dim][i] += ClusterVelocity[dim];
			}
			if (ClusterRadius2 > 0 && r2 > ClusterRadius2) { // in Enzo Unit
				Position[0][i] += 1e10;
				EscapeList.push_back(ptcl);
				ptcl->isErase = true;
				EscapeParticleNum++;
			}
			//fprintf(stdout, "NBODY+: pid= %d, x=%e\n",ptcl->PID,Position[0][i]);
			ptcl = ptcl->NextParticleInEnzo;

			if ((ptcl == nullptr) && (i != NNB-newNNB-1))
			{
				std::cout << "NBODY+: Warning! FirstParticleInEnzo is Null! No!" << std::endl;
			}
		}
		fprintf(stderr, "\n");
	}

	if (newNNB > 0) {
		for (int i=0; i<newNNB; i++) {
			r2 = 0;
			for (int dim=0; dim<Dim; dim++) {
				newPosition[dim][i]  = ptcl->Position[dim]/EnzoLength;
				newVelocity[dim][i]  = ptcl->Velocity[dim]/EnzoVelocity;
				r2                  += newPosition[dim][i]*newPosition[dim][i];

				// COM correction
				newPosition[dim][i] += ClusterPosition[dim];
				newVelocity[dim][i] += ClusterVelocity[dim];
			}
			if (ClusterRadius2 > 0 && r2 > ClusterRadius2) {
				newPosition[0][i] += 1e10;
				EscapeList.push_back(ptcl);
				ptcl->isErase = true;
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
	std::cout << "NBODY+: Escape particles=" << EscapeParticleNum << std::endl;

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

	fprintf(stderr, "NBODY: PID=");
	for (Particle* ptcl:particle) {
		fprintf(stderr, "%d, ", ptcl->PID);
		//fprintf(stderr, "NBODY: PID=%d, X=(%.3e, %.3e, %.3e)\n", 
				//ptcl->PID,ptcl->Position[0], ptcl->Position[1], ptcl->Position[2]);
		//fprintf(stderr, "NBODY: Pos  Of news=(%.3e, %.3e, %.3e)\n",
				//ptcl->Position[0], ptcl->Position[1], ptcl->Position[2]);
		//fprintf(stderr, "NBODY: Mass Of NewNbodyParticles=%.3e\n", ptcl->Mass);
		//fprintf(stderr, "NBODY: Vel  Of news=(%.3e, %.3e, %.3e)\n", 
		//	  ptcl->Velocity[0], ptcl->Velocity[1], ptcl->Velocity[2]);
		//fprintf(stderr, "NBODY: Acc  Of news=(%.3e, %.3e, %.3e)\n",
		//		ptcl->a_tot[0][0], ptcl->a_tot[1][0], ptcl->a_tot[2][0]);
	}
	fprintf(stderr, "\n");

	/*
	for (Particle* ptcl:particle) {
		fprintf(stderr, "NBODY:(particle before) PID=%d\n", ptcl->PID);
	}
	*/


		fprintf(stderr, "NBODY:(Escape) PID=\n");
	for (Particle* ptcl:EscapeList) {
		fprintf(stderr, "%d, ", ptcl->PID);
		/*
		fprintf(stderr, "NBODY3-1: Vel  Of news=(%.3e, %.3e, %.3e)\n", 
				ptcl->Velocity[0], ptcl->Velocity[1], ptcl->Velocity[2]);
		fprintf(stderr, "NBODY3-1: Pos  Of news=(%.3e, %.3e, %.3e)\n",
				ptcl->Position[0], ptcl->Position[1], ptcl->Position[2]);
		fprintf(stderr, "NBODY3-1: Acc  Of news=(%.3e, %.3e, %.3e)\n",
				ptcl->a_tot[0][0], ptcl->a_tot[1][0], ptcl->a_tot[2][0]);
				*/
	}
		fprintf(stderr, "\n ");

	// erase from other particles' neighbor
	for (Particle* ptcl: particle) {
		ptcl->ACList.erase(
				std::remove_if(ptcl->ACList.begin(), ptcl->ACList.end(),
					[](Particle* p) {
					bool to_remove = p->isErase;
					return to_remove;
					}),
				ptcl->ACList.end());
	}

	// erase from particle vector
	particle.erase(
			std::remove_if(particle.begin(), particle.end(),
				[](Particle* p) {
				bool to_remove = p->isErase;
				if (to_remove) delete p;
				return to_remove;
				}),
			particle.end());

	EscapeList.clear();


	fprintf(stderr, "NBODY:(particle after) PID=\n");
	for (Particle* ptcl:particle) {
		fprintf(stderr, "%d, ", ptcl->PID);
		/*
		fprintf(stderr, "NBODY3-2: Mass Of NewNbodyParticles=%.3e\n", ptcl->Mass);
		fprintf(stderr, "NBODY3-2: Vel  Of news=(%.3e, %.3e, %.3e)\n", 
				ptcl->Velocity[0], ptcl->Velocity[1], ptcl->Velocity[2]);
		fprintf(stderr, "NBODY3-2: Pos  Of news=(%.3e, %.3e, %.3e)\n",
				ptcl->Position[0], ptcl->Position[1], ptcl->Position[2]);
		fprintf(stderr, "NBODY3-2: Acc  Of news=(%.3e, %.3e, %.3e)\n",
				ptcl->a_tot[0][0], ptcl->a_tot[1][0], ptcl->a_tot[2][0]);
				*/
	}
		fprintf(stderr, "\n ");



	for (int dim = 0; dim < Dim; dim++)
	{
		if (NNB-newNNB != 0) {
			delete[] Position[dim];
			delete[] Velocity[dim];
		}
		if (newNNB > 0)
		{
			delete[] newPosition[dim];
			delete[] newVelocity[dim];
		}
	}

	NNB -= EscapeParticleNum;
	std::cout << "NBODY+: Sending data finished." << std::endl;


	fprintf(stdout, "NBODY+    : In SendToEzno (after particle might be escaped): \n");
	fprintf(stdout, "NBODY+    : original NNB      = %d (-%d)\n", NNB+EscapeParticleNum, EscapeParticleNum);
	fprintf(stdout, "NBODY+    : newly updated NNB = %d\n", NNB);
	fprintf(stdout, "NBODY+    : Particle size     = %d\n", particle.size());

	return true;
}

/* Adjust COM according to new particles
X = (a1+b1+c1)/M
Y = (a2+b2)/N
Z = (a1+b1+c1+a2+b2)/(M+N) = (X*M+Y*N)/(M+N) -> New COM
Z-X = N*(Y-X)/(M+N) -> This should be applied to particles
*/
// this should be improved by using iterative loop
void GetNewCenterOfMass(std::vector<Particle*> &particle, double *mass2, double *x2[Dim], double *v2[Dim], int n2, double x_X[], double v_X[]) {

	double M=0., N=0.;
	double x_Y[Dim], v_Y[Dim], x_Z[Dim], v_Z[Dim];

	for (int dim=0; dim<Dim; dim++) {
		x_Y[dim] = 0;
		x_Z[dim] = 0;
		v_Y[dim] = 0;
		v_Z[dim] = 0;
	}

	for (Particle* ptcl:particle)
		M += ptcl->Mass/EnzoMass;


	for (int i=0; i<n2; i++) {
		for (int dim=0; dim<Dim; dim++) {
			x_Y[dim] += mass2[i]*x2[dim][i];
			v_Y[dim] += mass2[i]*v2[dim][i];
		}
		N += mass2[i];
	}
	for (int dim=0; dim<Dim; dim++) {
		x_Y[dim] = x_Y[dim]/N;
		v_Y[dim] = v_Y[dim]/N;
		x_Z[dim] = (x_X[dim]*M+x_Y[dim]*N)/(M+N);
		v_Z[dim] = (v_X[dim]*M+v_Y[dim]*N)/(M+N);
	}


	// Adjustment to particles
	for (Particle* ptcl:particle) {
		for (int dim=0; dim<Dim;dim++) {
			ptcl->Position[dim] - (x_Z[dim] - x_X[dim]) * EnzoLength;
			ptcl->Velocity[dim] - (v_Z[dim] - v_X[dim]) * EnzoVelocity;
		}
	}

	for (int dim=0; dim<Dim; dim++) {
		x_X[dim] = x_Z[dim];
		v_X[dim] = v_Z[dim];
	}
}


