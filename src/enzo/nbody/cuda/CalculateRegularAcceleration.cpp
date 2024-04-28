#include <vector>
#include <iostream>
#include "../global.h"
#include <cmath>
#include "defs.h"
#include "cuda_functions.h"


void SendAllParticlesToGPU(double current_time, double next_time, std::vector <Particle*> &particle);
void CalculateSingleAcceleration(Particle *ptcl1, Particle *ptcl2, double (&a)[3], double (&adot)[3]);


/*
 *  Purporse: calculate acceleration and neighbors of regular particles by sending them to GPU
 *
 *  Date    : 2024.01.18  by Seoyoung Kim
 *
 */
void CalculateRegAccelerationOnGPU(std::vector<int> IndexList, std::vector<Particle*> &particle){
	// regIds are the list of positions of particles subject to regular force calculation in std::vector list particle

	// variables for opening GPU
	// const int buffer = 10;
	// int numGpuOpen = NNB+buffer;
	//const int NumPtclPerEachCalMax = 2048; // this also caps the number of particles computed each iteration
	const int NumPtclPerEachCalMax = 512; // this also caps the number of particles computed each iteration
	const int mpi_rank  = 0; // not effective for now
	int NumPtclPerEachCal; // this size should coincide with number of threads 
	int NeighborIndex; // this size should coincide with number of threads 
	int ListSize = IndexList.size();

	//int NumGpuCal;

	// variables for saving variables to send to GPU
	// only regular particle informations are stored here

	double *MassSend;
	double *MdotSend;
	double(*PositionSend)[Dim];
	double(*VelocitySend)[Dim];
	double(*AccReceive)[Dim];
	double(*AccDotReceive)[Dim];

	double* RadiusOfAC2Send;
	double* TimeStepRegSend;

	double(*AccRegReceive[2])[Dim];
	double(*AccRegDotReceive[2])[Dim];
	double(*AccIrr[2])[Dim];
	double(*AccIrrDot[2])[Dim];

	//double* PotSend;
	int **ACListReceive[2];
	int *NumNeighborReceive[2] = {0};
	int MassFlag;

	// the current and next times to calculate the regular force on
	double treg[2];  // 0 is current time and 1 is next time

	// temporary variables saving the calculated acceleration values
	double a_tmp[Dim], adot_tmp[Dim];

	// variables for calculating the higher order acceleration derivatives
	double a2, a3, da_dt2, adot_dt, dt, dt2, dt3, dt4;

	Particle *ptcl;
	int j2;



	// need to make array to send to GPU
	// allocate memory to the temporary variables
	MassSend        = new double[ListSize];
	MdotSend        = new double[ListSize];
	PositionSend    = new double[ListSize][Dim];
	VelocitySend    = new double[ListSize][Dim];
	AccReceive      = new double[ListSize][Dim];
	AccDotReceive   = new double[ListSize][Dim];

	RadiusOfAC2Send = new double[ListSize];
	TimeStepRegSend = new double[ListSize];
	//PotSend         = new double[ListSize];

	for (int p=0; p<2; p++) {
		AccRegReceive[p]    = new double[ListSize][Dim];
		AccRegDotReceive[p] = new double[ListSize][Dim];
		AccIrr[p]           = new double[ListSize][Dim];
		AccIrrDot[p]        = new double[ListSize][Dim];

		NumNeighborReceive[p] = new int[ListSize];
		ACListReceive[p]      = new int*[ListSize];
		for (int i=0; i<ListSize; i++) {
			ACListReceive[p][i] = new int[NumNeighborMax];
			for (int dim=0; dim<Dim; dim++) {
				AccIrr[p][i][dim]    = 0;
				AccIrrDot[p][i][dim] = 0;
			}
		}
	}

	// perform the loop twice in order to obtain
	// both the current and predicted acceleration
	// set the current time to 0 and next time to 1
	// and set the time step dt to regular time step
	treg[0] = particle[IndexList[0]]->CurrentTimeReg;  // current regular time
	dt      = particle[IndexList[0]]->TimeStepReg;
	treg[1] = treg[0] + dt;  // next regular time
	if (treg[1] != NextRegTime) {
		fprintf(stderr, "Something wrong! NextRegTime does not match! :CalculateAcceleration.C:105\n");
		fprintf(stderr, "NextRegTime=%.3e, treg[1]=%.3e\n", NextRegTime, treg[1]);
	}
	dt     *= EnzoTimeStep;  // unit conversion


	// first let's open the GPU
	//OpenDevice(&mpi_rank);

	//std::cout <<  "Starting Calculation On Device ..." << std::endl;
	for (int p=0; p<2; p++) {
		// send information of all the particles to GPU
		// includes prediction
#ifdef time_trace
		_time.reg_sendall.markStart();
#endif

		SendAllParticlesToGPU(particle[IndexList[0]]->CurrentTimeReg, treg[p], particle);

#ifdef time_trace
		_time.reg_sendall.markEnd();
		_time.reg_sendall.getDuration();
#endif

		// copy the data of regular particles to the arrays to be sent
		// predicted positions and velocities should be sent
		// but predictions are already done when sending all particles, so no need for duplicated calculation
		for (int i=0; i<ListSize; i++) {
			MassSend[i]   = particle[IndexList[i]]->Mass;
			MdotSend[i]   = 0.; // I will do it later
			RadiusOfAC2Send[i] = particle[IndexList[i]]->RadiusOfAC*particle[IndexList[i]]->RadiusOfAC;

			for (int dim=0; dim<Dim; dim++) {
				PositionSend[i][dim] = particle[IndexList[i]]->PredPosition[dim];
				VelocitySend[i][dim] = particle[IndexList[i]]->PredVelocity[dim];
			}
			//TimeStepRegSend[i] = particle[IndexList[i]]->TimeStepReg;
		}

		// calculate the force by sending the particles to GPU
#ifdef time_trace
		_time.reg_gpu.markStart();
#endif
		CalculateAccelerationOnDevice(&ListSize, PositionSend, VelocitySend, AccRegReceive[p], AccRegDotReceive[p],
			 	MdotSend, RadiusOfAC2Send, NumNeighborReceive[p], ACListReceive[p]);
#ifdef time_trace
		_time.reg_gpu.markEnd();
		_time.reg_gpu.getDuration();

		_time.reg_cpu1.markStart();
#endif

		// Calculate the irregular acceleration components based on neighbors of current regular time.
		for (int i2=0; i2<ListSize; i2++) {
			ptcl = particle[IndexList[i2]];  // regular particle in particle list

			//std::cout <<  "MyIndex=" << IndexList[i2];
			//std::cout <<  "(" << NumNeighborReceive[0][i2] << ")" << std::endl;
			//std::cout <<  "NeighborIndex = ";
			for (int j1=0;  j1<NumNeighborReceive[0][i2]; j1++) {
				NeighborIndex = ACListReceive[0][i2][j1];  // gained neighbor particle (in next time list)
				//std::cout <<  NeighborIndex << "  ";
				CalculateSingleAcceleration(ptcl, particle[NeighborIndex], a_tmp, adot_tmp);

				for (int dim=0; dim<Dim; dim++) {
					AccIrr[p][i2][dim]           += a_tmp[dim];
					AccIrrDot[p][i2][dim]        += adot_tmp[dim];
					AccRegReceive[p][i2][dim]    -= a_tmp[dim];
					AccRegDotReceive[p][i2][dim] -= adot_tmp[dim];
				} // endfor dim
			} // endfor j1, over neighbor at current time
			//std::cout << std::endl;

			if (p == 0) {
				for (int dim=0; dim<Dim; dim++) {
					ptcl->a_reg[dim][0] = AccRegReceive[0][i2][dim];
					ptcl->a_reg[dim][1] = AccRegDotReceive[0][i2][dim];
					ptcl->a_irr[dim][0] = AccIrr[0][i2][dim];
					ptcl->a_irr[dim][1] = AccIrrDot[0][i2][dim];
					ptcl->a_tot[dim][0] = ptcl->a_reg[dim][0] + ptcl->a_irr[dim][0];
					ptcl->a_tot[dim][1] = ptcl->a_reg[dim][1] + ptcl->a_irr[dim][1];
				}
			} // current time update ptcl acc
		} // endfor i2, over regular particles
#ifdef time_trace
		_time.reg_cpu1.markEnd();
		_time.reg_cpu1.getDuration();
#endif
	} // endfor p,

	//std::cout <<  "Calculation On Device Done ..." << std::endl;

	for (int i=0; i<ListSize; i++) {
#ifdef time_trace
		_time.reg_cpu2.markStart();
#endif
		ptcl = particle[IndexList[i]];  // regular particle in particle list
		// CALCULATE AND SAVE THE 0TH, 1ST, 2ND, 3RD DERIVATIVES OF REGULAR AND IRREGULAR FORCE
		dt2 = dt*dt;
		dt3 = dt2*dt;
		dt4 = dt3*dt;

		// calculated the final corrected forces
		for (int dim=0; dim<Dim; dim++) {
			da_dt2  = (   AccRegReceive[0][i][dim] - AccRegReceive[1][i][dim]   ) / dt2;
			adot_dt = (AccRegDotReceive[0][i][dim] + AccRegDotReceive[1][i][dim]) / dt;

			a2 =  -6*da_dt2 - 2*adot_dt - 2*AccRegDotReceive[0][i][dim]/dt;
			a3 = (12*da_dt2 + 6*adot_dt)/dt;

			ptcl->a_reg[dim][2] = a2;
			ptcl->a_reg[dim][3] = a3;
		}

		for (int dim=0; dim<Dim; dim++) {
			da_dt2  = (   AccIrr[0][i][dim] - AccIrr[1][i][dim]   ) / dt2;
			adot_dt = (AccIrrDot[0][i][dim] + AccIrrDot[1][i][dim]) / dt;

			a2 =  -6*da_dt2  - 2*adot_dt - 2*AccIrrDot[0][i][dim]/dt;
			a3 = (12*da_dt2 + 6*adot_dt)/dt;

			ptcl->a_irr[dim][2] = a2;
			ptcl->a_irr[dim][3] = a3;

			ptcl->a_tot[dim][2] = ptcl->a_reg[dim][2] + ptcl->a_irr[dim][2];
			ptcl->a_tot[dim][3] = ptcl->a_reg[dim][3] + ptcl->a_irr[dim][3];
		}

#ifdef time_trace
		_time.reg_cpu2.markEnd();
		_time.reg_cpu2.getDuration();

		_time.reg_cpu3.markStart();
#endif
		// Neighbor update
		ptcl->ACList.clear();
		ptcl->NumberOfAC = NumNeighborReceive[0][i];
		for (int j=0; j<ptcl->NumberOfAC;j++) {
			NeighborIndex = ACListReceive[0][i][j];  // gained neighbor particle (in next time list)
			ptcl->ACList.push_back(particle[NeighborIndex]);
		}

#ifdef time_trace
		_time.reg_cpu3.markEnd();
		_time.reg_cpu3.getDuration();

		_time.reg_cpu4.markStart();
#endif
		// Particle Update
		if (ptcl->NumberOfAC == 0) {
			//ptcl->updateParticle(ptcl->CurrentTimeReg, treg[1], ptcl->a_tot);
			ptcl->correctParticleFourthOrder(ptcl->CurrentTimeReg, treg[1], ptcl->a_tot);
			ptcl->CurrentTimeReg += ptcl->TimeStepReg;
			ptcl->CurrentTimeIrr  = ptcl->CurrentTimeReg;
		}
		else {
			//ptcl->CurrentTimeReg += ptcl->TimeStepReg;
			ptcl->correctParticleFourthOrder(ptcl->CurrentTimeReg, treg[1], ptcl->a_reg);
			ptcl->CurrentTimeReg = ptcl->CurrentTimeIrr;
		}
		ptcl->calculateTimeStepReg(ptcl->a_reg, ptcl->a_reg); // open question
#ifdef time_trace
		_time.reg_cpu4.markEnd();
		_time.reg_cpu4.getDuration();

		_time.reg_cpu5.markStart();
#endif
		ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
		ptcl->isRegular = false;
#ifdef time_trace
		_time.reg_cpu5.markEnd();
		_time.reg_cpu5.getDuration();
#endif
	} // correction and calculation of higher orders finished





	/*
	std::cout <<  "3. a_tot= "<< particle[0]->a_tot[0][0]<< ',' << particle[0]->a_tot[1][0]\
		<< ',' << particle[0]->a_tot[2][0] << std::endl;
	std::cout <<  "4. a_tot= "<< particle[1]->a_tot[0][0]<< ',' << particle[1]->a_tot[1][0]\
		<< ',' << particle[1]->a_tot[2][0] << std::endl;

	std::cout <<  "3. a_irr= "<< particle[0]->a_irr[0][0]<< ',' << particle[0]->a_irr[1][0]\
		<< ',' << particle[0]->a_irr[2][0] << std::endl;
	std::cout <<  "4. a_irr= "<< particle[1]->a_irr[0][0]<< ',' << particle[1]->a_irr[1][0]\
		<< ',' << particle[1]->a_irr[2][0] << std::endl;
		*/


	// free all temporary variables
	delete[] MassSend;
	delete[] MdotSend;
	delete[] PositionSend;
	delete[] VelocitySend;
	delete[] RadiusOfAC2Send;
	delete[] TimeStepRegSend;
	//delete[] PotSend;

	for (int p=0; p<2; p++) {
		delete[] NumNeighborReceive[p];
		delete[] AccRegReceive[p];
		delete[] AccRegDotReceive[p];
		delete[] AccIrr[p];
		delete[] AccIrrDot[p];
		delete[] ACListReceive[p];
	}

	//CloseDevice();
} // calculate 0th, 1st derivative of force + neighbors on GPU ends




/*
 *  Purporse: send the information of all particles to GPU in regular integration steps
 *  send the predicted positions and velocities (consistent prediction must be performed before sending)
 *
 *  -> calculates based on their current positions
 *
 *  Date    : 2024.01.17  by Seoyoung Kim
 *  Date    : 2024.02.07  by Yongseok Jo
 *
 */

void CalculateSingleAcceleration(Particle *ptcl1, Particle *ptcl2, double (&a)[3], double (&adot)[3]) {
	double dx[Dim], dv[Dim];
	double dr2;
	double dxdv;
	double m_r3;

	if (ptcl1 == ptcl2) {
		for (int dim=0; dim<Dim; dim++){
			a[dim]    = 0;
			adot[dim] = 0;
		}
		return;
	}

	dr2  = 0.0;
	dxdv = 0.0;

	for (int dim=0; dim<Dim; dim++) {
		dx[dim] = ptcl2->PredPosition[dim] - ptcl1->PredPosition[dim];
		dv[dim] = ptcl2->PredVelocity[dim] - ptcl1->PredVelocity[dim];
		dr2    += dx[dim]*dx[dim];
		dxdv   += dx[dim]*dv[dim];
	}

	/*
	if (dr2 < EPS2) {
		dr2 += EPS2;
	}
	*/

	m_r3 = ptcl2->Mass/dr2/sqrt(dr2);

	for (int dim=0; dim<Dim; dim++){
		a[dim]    = m_r3*dx[dim];
		adot[dim] = m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
	}
}
