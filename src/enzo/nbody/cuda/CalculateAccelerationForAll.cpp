#include <vector>
#include <iostream>
#include "../global.h"
#include <cmath>
#include "defs.h"
#include "cuda_functions.h"

/*
 *  Purporse: calculate acceleration and neighbors of all the particles on GPU.
 *  -> calculates based on their current positions
 *
 *  Date    : 2024.01.17  by Seoyoung Kim
 *
 */

void CalculateSingleAcceleration(Particle *ptcl1, Particle *ptcl2, double (&a)[3], double (&adot)[3], int sign);
void SendAllParticlesToGPU(double time, std::vector <Particle*> &particle);
void UpdateNextRegTime(std::vector<Particle*> &particle);

void CalculateAllAccelerationOnGPU(std::vector<Particle*> &particle){




	//int NumGpuCal;

	// variables for saving variables to send to GPU
	// only regular particle informations are stored here

	const int mpi_rank  = 0; // not effective for now
	int NeighborIndex; // this size should coincide with number of threads
	int ListSize = particle.size();
	int *IndexList = new int[ListSize];
	
	double (*AccRegReceive)[Dim];
	double (*AccRegDotReceive)[Dim];
	double (*AccIrr)[Dim];
	double (*AccIrrDot)[Dim];

	//double* PotSend;
	int **ACListReceive;
	int *NumNeighborReceive;
	int MassFlag;


	double a_tmp[Dim]{0}, adot_tmp[Dim]{0};
	double da, dadot;

	Particle * ptcl;

	// need to make array to send to GPU
	// allocate memory to the temporary variables

	AccRegReceive    = new double[ListSize][Dim];
	AccRegDotReceive = new double[ListSize][Dim];
	AccIrr           = new double[ListSize][Dim];
	AccIrrDot        = new double[ListSize][Dim];

	NumNeighborReceive  = new int[ListSize];
	ACListReceive      = new int*[ListSize];

	for (int i=0; i<ListSize; i++) {
		ACListReceive[i] = new int[NumNeighborMax];
		for (int dim=0; dim<Dim; dim++) {
			AccIrr[i][dim]				   = 0;
			AccIrrDot[i][dim]			 	 = 0;
			AccRegReceive[i][dim]    = 0;
			AccRegDotReceive[i][dim] = 0;
		}
	}


	for (int i=0; i<ListSize; i++) {
		IndexList[i] = i;
	} // endfor copy info

	// send the arrays to GPU
	//SendToDevice(&NNB, MassSend, PositionSend, VelocitySend, MdotSend, &FixNumNeighbor);


	// Particles have been already at T_new through irregular time step
	SendAllParticlesToGPU(0., particle);  // needs to be updated
	

	 // endfor copy info

	// calculate the force by sending the particles to GPU
	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive, AccRegDotReceive, NumNeighborReceive, ACListReceive);

	// Calculate the irregular acceleration components based on neighbors of current regular time.
	for (int i=0; i<ListSize; i++) {

		ptcl = particle[i];  // regular particle in particle list

		for (int dim=0; dim<Dim; dim++) {
			a_tmp[dim]    = 0.;
			adot_tmp[dim] = 0.;
		}

		/*******************************************************
		 * Acceleartion correction according to current neighbor
		 ********************************************************/
		for (int j=0;  j<NumNeighborReceive[i]; j++) {
			NeighborIndex = ACListReceive[i][j];  // gained neighbor particle (in next time list)
			CalculateSingleAcceleration(ptcl, particle[NeighborIndex], a_tmp, adot_tmp, 1);
		} // endfor j1, over neighbor at current time

		// update force
		for (int dim=0; dim<Dim; dim++) {
			ptcl->a_reg[dim][0] = AccRegReceive[i][dim];
			ptcl->a_reg[dim][1] = AccRegDotReceive[i][dim];
			ptcl->a_irr[dim][0] = a_tmp[dim];    //AccIrr[i][dim];
			ptcl->a_irr[dim][1] = adot_tmp[dim]; //AccIrrDot[i][dim];
			ptcl->a_tot[dim][0] = ptcl->a_reg[dim][0] + ptcl->a_irr[dim][0]; //+ ptcl->BackgroundAcceleration[dim];
			ptcl->a_tot[dim][1] = ptcl->a_reg[dim][1] + ptcl->a_irr[dim][1];
			// in case
		}

		ptcl->ACList.clear();
		ptcl->NumberOfAC = NumNeighborReceive[i];
		for (int j=0; j<ptcl->NumberOfAC;j++) {
			NeighborIndex = ACListReceive[i][j];  // gained neighbor particle (in next time list)
			ptcl->ACList.push_back(particle[NeighborIndex]);
		}
		ptcl->UpdateRadius();
	} 


	// free all temporary variables
	delete[] NumNeighborReceive;
	delete[] AccRegReceive;
	delete[] AccRegDotReceive;
	delete[] AccIrr;
	delete[] AccIrrDot;
	for (int i=0; i<ListSize; i++) {
		delete[] ACListReceive[i];
	}
	delete[] ACListReceive;
	
} // calculate 0th, 1st derivative of force + neighbors on GPU ends




void SendAllParticlesToGPU(double time, std::vector <Particle*> &particle) {

	// variables for saving variables to send to GPU
	double * Mass;
	double * Mdot;
	double * Radius2;
	double(*Position)[Dim];
	double(*Velocity)[Dim];
	int size = (int) particle.size();

	// allocate memory to the temporary variables
	Mass     = new double[size];
	Mdot     = new double[size];
	Radius2  = new double[size];
	Position = new double[size][Dim];
	Velocity = new double[size][Dim];


	// copy the data of particles to the arrays to be sent
	for (int i=0; i<size; i++) {
		Mass[i]    = particle[i]->Mass;
		Mdot[i]    = 0; //particle[i]->Mass;
		Radius2[i] = particle[i]->RadiusOfAC*particle[i]->RadiusOfAC; // mass wieght?
		if (particle[i]->NumberOfAC == 0)
			particle[i]->predictParticleSecondOrder(time);
		else
			particle[i]->predictParticleSecondOrderIrr(time);

		for (int dim=0; dim<Dim; dim++) {
			Position[i][dim] = particle[i]->PredPosition[dim];
			Velocity[i][dim] = particle[i]->PredVelocity[dim];
		}
	}

	//fprintf(stdout, "Sending particles to GPU...\n");
	//fflush(stdout);
	// send the arrays to GPU
	SendToDevice(&size, Mass, Position, Velocity, Radius2, Mdot);

	//fprintf(stdout, "Done.\n");
	//fflush(stdout);
	// free the temporary variables
	delete[] Mass;
	delete[] Mdot;
	delete[] Radius2;
	delete[] Position;
	delete[] Velocity;
}

