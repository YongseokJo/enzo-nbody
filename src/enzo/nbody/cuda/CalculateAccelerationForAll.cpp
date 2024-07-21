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
 *  Date    : 2024.06.28  by Yongseok Jo
 *
 */

void CalculateSingleAcceleration(Particle *ptcl1, Particle *ptcl2, double (&a)[3], double (&adot)[3], int sign);
void SendAllParticlesToGPU(double time, std::vector <Particle*> &particle);
void UpdateNextRegTime(std::vector<Particle*> &particle);


void CalculateAllAccelerationOnGPU(std::vector<Particle*> &particle){



	const int mpi_rank  = 0; // not effective for now
	int NeighborIndex; // this size should coincide with number of threads

	//int NumGpuCal;

	// variables for saving variables to send to GPU
	// only regular particle informations are stored here
	double *MassSend;
	double *MdotSend;
	double(*PositionSend)[Dim];
	double(*VelocitySend)[Dim];

	double* RadiusOfAC2Send;
	double* TimeStepRegSend;

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


	// need to make array to send to GPU
	// allocate memory to the temporary variables
	MassSend        = new double[NNB];
	MdotSend        = new double[NNB];
	PositionSend    = new double[NNB][Dim];
	VelocitySend    = new double[NNB][Dim];

	RadiusOfAC2Send = new double[NNB];
	TimeStepRegSend = new double[NNB];
	//PotSend         = new double[NNB];

	AccRegReceive    = new double[NNB][Dim];
	AccRegDotReceive = new double[NNB][Dim];
	AccIrr           = new double[NNB][Dim];
	AccIrrDot        = new double[NNB][Dim];

	NumNeighborReceive  = new int[NNB];
	ACListReceive      = new int*[NNB];

	for (int i=0; i<NNB; i++) {
		ACListReceive[i] = new int[NumNeighborMax];
		for (int dim=0; dim<Dim; dim++) {
			AccIrr[i][dim]    = 0;
			AccIrrDot[i][dim] = 0;
		}
	}


	// send the arrays to GPU
	//SendToDevice(&NNB, MassSend, PositionSend, VelocitySend, MdotSend, &FixNumNeighbor);


	// Particles have been already at T_new through irregular time step
	SendAllParticlesToGPU(0., particle);  // needs to be updated
	
	
	Particle * ptcl;

	for (int i=0; i<particle.size(); i++) {
		ptcl = particle[i];
		MassSend[i]   = ptcl->Mass;
		MdotSend[i]   = 0.; // I will do it later
		RadiusOfAC2Send[i] = ptcl->RadiusOfAC*ptcl->RadiusOfAC/ptcl->Mass; // mass weighted

		for (int dim=0; dim<Dim; dim++) {
			PositionSend[i][dim] = ptcl->Position[dim];
			VelocitySend[i][dim] = ptcl->Velocity[dim];
		}
	} // endfor copy info

	// calculate the force by sending the particles to GPU
	CalculateAccelerationOnDevice(&NNB, PositionSend, VelocitySend, AccRegReceive, AccRegDotReceive,
			MdotSend, RadiusOfAC2Send, NumNeighborReceive, ACListReceive, 0.);

	// Calculate the irregular acceleration components based on neighbors of current regular time.
	for (int i=0; i<particle.size(); i++) {

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
	delete[] MassSend;
	delete[] MdotSend;
	delete[] PositionSend;
	delete[] VelocitySend;
	delete[] RadiusOfAC2Send;
	delete[] TimeStepRegSend;
	//delete[] PotSend;

	delete[] NumNeighborReceive;
	delete[] AccRegReceive;
	delete[] AccRegDotReceive;
	delete[] AccIrr;
	delete[] AccIrrDot;
	for (int i=0; i<NNB; i++) {
		delete[] ACListReceive[i];
	}
	delete[] ACListReceive;
	
} // calculate 0th, 1st derivative of force + neighbors on GPU ends



/*
 *  Purporse: calculate acceleration and neighbors of some particles on GPU.
 *  the particles subject for calculation of regular force are given by std::vector<Particle*> list
 *  predicted positions and velocities are used for calculation on GPU for now
 * 
 *  -> calculates based on their current positions
 *
 *  Date    : 2024.01.17  by Seoyoung Kim
 *
 */
void CalculateListAccelerationOnGPU(std::vector<int> &IndexList, std::vector<Particle*> &particle){

	// variables for opening GPU
	int ListSize   = IndexList.size();
	int numGpuSend = 1024;
	int mpi_rank   = 0;
	int numGpuCal;

	// variables for GPU
	double *MassSend;
	double *MdotSend;
	double(*PositionSend)[Dim];
	double(*VelocitySend)[Dim];

	double *r2OfACSend;
	double *TimeStepRegSend;

	double(*AccSend)[Dim];
	double(*AccDotSend)[Dim];
	double *PotSend;
	int **AClistGpu;
	int massFlag;

	// temporary variables for calculating the irregular force

	double dx[Dim], dv[Dim];
	double rij2, dr2i, dr3i, drdv;

	// extra variables
	int ACnumi2, ACjid, i2reg;


	// need to make array to send to GPU
	// allocate memory to the temporary variables
	MassSend     = new double[ListSize];
	MdotSend     = new double[ListSize];
	PositionSend = new double[ListSize][Dim];
	VelocitySend = new double[ListSize][Dim];

	r2OfACSend   = new double[ListSize];
	TimeStepRegSend  = new double[ListSize];

	AccSend       = new double[ListSize][Dim];
	AccDotSend    = new double[ListSize][Dim];
	PotSend       = new double[ListSize];

	AClistGpu = new int*[ListSize];
	for (int i=0; i<(NNB); i++) {
		AClistGpu[i] = new int[FixNumNeighbor];
	}


	// copy the data of particles to the arrays to be sent
	// predicted positions and velocities are sent

	for (int i=0; i<ListSize; i++) {
		rij2          = 0;
		MassSend[i]   = particle[(IndexList[i])]->Mass;
		MdotSend[i]   = particle[(IndexList[i])]->Mass;
		r2OfACSend[i] = particle[(IndexList[i])]->RadiusOfAC * particle[(IndexList[i])]->RadiusOfAC;

		for (int dim=0; dim<Dim; dim++) {
			PositionSend[i][dim] = particle[(IndexList[i])]->PredPosition[dim];
			VelocitySend[i][dim] = particle[(IndexList[i])]->PredVelocity[dim];
		}
		TimeStepRegSend[i] = particle[(IndexList[i])]->TimeStepReg;
	}

	// send the arrays to GPU
	//SendAllParticlesToGPU(particle);

	// calculate the force by sending the particles to GPU in multiples of 1024
	for (int i=0; i<ListSize; i+=numGpuSend) {
		numGpuCal = std::min(1024,(NNB-i));

		//CalculateAccelerationOnDevice(&NNB,&i, PositionSend, VelocitySend,
		//			AccSend, AccDotSend, MdotSend, r2OfACSend);

		// copy the values of regular forces and neighbors obtained in GPU to particles
		// and also calculate the irregular forces in the process
		for (int i2=i; i2<(i+1024); i2++){
			i2reg = (IndexList[i2]);

			// need to clear the AClist and the number of neighbors before copying the new values
			particle[i2reg]->NumberOfAC = 0;
			particle[i2reg]->ACList.clear();

			for (int dim=0; dim<Dim; dim++) {
				particle[i2reg]->a_reg[dim][0] += AccSend[i2][dim];
				particle[i2reg]->a_reg[dim][1] += AccDotSend[i2][dim];
			}

			ACnumi2 = AClistGpu[i][0];

			// if there are neighbors, save them to neighbor list and calculate the irregular force
			if (ACnumi2 > 0){
				particle[i2reg]->NumberOfAC = ACnumi2;

				for (int j=1; j<(ACnumi2+1); j++) {
					ACjid = AClistGpu[i][j];
					particle[i2reg]->ACList.push_back(particle[ACjid]);
					rij2 = 0.0;
					drdv = 0.0;

					for (int dim=0; dim<Dim; dim++) {
						dx[dim] = particle[ACjid]->PredPosition - particle[i2reg]->PredPosition;
						dv[dim] = particle[ACjid]->PredVelocity - particle[i2reg]->PredVelocity;
						rij2   += dx[dim]*dx[dim];
						drdv   += dx[dim]*dv[dim];
					}
					dr2i = 1.0/rij2;
					dr3i = particle[ACjid]->Mass*dr2i*sqrt(dr2i);

					for (int dim=0; dim<Dim; dim++) {
						particle[i2reg]->a_irr[dim][0] += dx[dim]*dr3i;
						particle[i2reg]->a_irr[dim][1] += (dv[dim]-dx[dim]*drdv)*dr3i;
					}
				} // loop on neighbors end
			} else { // in case of no neighbors, just set the neighbor number to 0 just in case. 
				particle[i2]->NumberOfAC = 0;
			} // if statement on neighbor number ends

			for (int dim=0; dim<3; dim++){ // copy the values to other values as well
				particle[i2reg]->a_tot[dim][0] = particle[i2reg]->a_reg[dim][0] + particle[i2reg]->a_irr[dim][0];
				particle[i2reg]->a_tot[dim][1] = particle[i2reg]->a_reg[dim][1] + particle[i2reg]->a_irr[dim][1];
			}
		} // saving the 1024 results from GPU to local class ends
	} // loop of total calculation ends


	// free all temporary variables
	delete[] MassSend;
	delete[] PositionSend;
	delete[] VelocitySend;

	delete[] r2OfACSend;
	delete[] TimeStepRegSend;

	delete[] AccSend;
	delete[] AccDotSend;
	delete[] PotSend;
	delete[] AClistGpu;

} // calculate 0th, 1st derivative of force + neighbors on GPU ends



/*
 *  Purporse: send the information of all particles to GPU in regular integration steps
 *  send the predicted positions and velocities (consistent prediction must be performed before sending)
 *
 *  -> calculates based on their current positions
 *
 *  Date    : 2024.01.17  by Seoyoung Kim
 *
 */


void SendAllParticlesToGPU(double time, std::vector <Particle*> &particle) {

	// variables for saving variables to send to GPU
	double * Mass;
	double * Mdot;
	double(*Position)[Dim];
	double(*Velocity)[Dim];
	int size = (int) particle.size();
	int num=100; // this is garbage

	// allocate memory to the temporary variables
	Mass     = new double[size];
	Mdot     = new double[size];
	Position = new double[size][Dim];
	Velocity = new double[size][Dim];

	// copy the data of particles to the arrays to be sent
	for (int i=0; i<size; i++) {
		Mass[i] = particle[i]->PredMass;
		Mdot[i] = 0; //particle[i]->Mass;
		if (particle[i]->NumberOfAC == 0)	
			particle[i]->predictParticleSecondOrder(time);
		else
			particle[i]->predictParticleSecondOrderIrr(time);

		for (int dim=0; dim<Dim; dim++) {
			Position[i][dim] = particle[i]->PredPosition[dim];
			Velocity[i][dim] = particle[i]->PredVelocity[dim];
		}
	}

	// send the arrays to GPU
	SendToDevice(&size, Mass, Position, Velocity, Mdot, &num);

	// free the temporary variables
	delete[] Mass;
	delete[] Mdot;
	delete[] Position;
	delete[] Velocity;
}

