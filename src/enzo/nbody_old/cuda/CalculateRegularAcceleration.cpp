#include <vector>
#include <iostream>
#include "../global.h"
#include <cmath>
#include "defs.h"
#include "cuda_functions.h"


void SendAllParticlesToGPU(double time, std::vector <Particle*> &particle);
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
	const int mpi_rank  = 0; // not effective for now
	int NeighborIndex; // this size should coincide with number of threads
	int ListSize = IndexList.size();

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
	double new_time;  // 0 is current time and 1 is next time
	double a2, a3, da_dt2, adot_dt, dt, dt2, dt3, dt4, dt5;


	double DFR, FRD, SUM, AT3, BT2;
	double DTR, DTSQ, DT2, DT6,DTSQ12, DTR13;

	Particle *ptcl;



	// need to make array to send to GPU
	// allocate memory to the temporary variables
	MassSend        = new double[ListSize];
	MdotSend        = new double[ListSize];
	PositionSend    = new double[ListSize][Dim];
	VelocitySend    = new double[ListSize][Dim];

	RadiusOfAC2Send = new double[ListSize];
	TimeStepRegSend = new double[ListSize];
	//PotSend         = new double[ListSize];

	AccRegReceive    = new double[ListSize][Dim];
	AccRegDotReceive = new double[ListSize][Dim];
	AccIrr           = new double[ListSize][Dim];
	AccIrrDot        = new double[ListSize][Dim];

	NumNeighborReceive  = new int[ListSize];
	ACListReceive      = new int*[ListSize];
	for (int i=0; i<ListSize; i++) {
		ACListReceive[i] = new int[NumNeighborMax];
		for (int dim=0; dim<Dim; dim++) {
			AccIrr[i][dim]    = 0;
			AccIrrDot[i][dim] = 0;
		}
	}

	// perform the loop twice in order to obtain
	// both the current and predicted acceleration
	// set the current time to 0 and next time to 1
	// and set the time step dt to regular time step

	dt    = particle[IndexList[0]]->TimeStepReg;
	new_time = particle[IndexList[0]]->CurrentTimeReg + dt;  // next regular time
	if (new_time != NextRegTime) {
		if (NextRegTime == 0) {
			fprintf(stderr, "First RegularCalculation skips! :CalculateAcceleration.C:105\n");
		}
		else{
			fprintf(stderr, "Something wrong! NextRegTime does not match! :CalculateAcceleration.C:105\n");
			fprintf(stderr, "NextRegTime=%.3e, treg[1]=%.3e\n", NextRegTime, new_time);
		}
		return;
	}


	dt *= EnzoTimeStep;  // unit conversion
	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;

	/*
	DTR = dt;
	DTSQ = DTR*DTR;
	DT6 = 6.0/(DTR*DTSQ);
	DT2 = 2.0/DTSQ;
	DTSQ12 = DTSQ/12;
	DTR13 = DTR/3;
	*/

	//std::cout <<  "Starting Calculation On Device ..." << std::endl;
	// send information of all the particles to GPU
	// includes prediction
#ifdef time_trace
		_time.reg_sendall.markStart();
#endif

		// Particles have been already at T_new through irregular time step
		SendAllParticlesToGPU(new_time, particle);  // needs to be updated

#ifdef time_trace
		_time.reg_sendall.markEnd();
		_time.reg_sendall.getDuration();
#endif

		// copy the data of regular particles to the arrays to be sent
		// predicted positions and velocities should be sent
		// but predictions are already done when sending all particles, so no need for duplicated calculation
		for (int i=0; i<ListSize; i++) {
			ptcl = particle[IndexList[i]];
			MassSend[i]   = ptcl->Mass;
			MdotSend[i]   = 0.; // I will do it later
			RadiusOfAC2Send[i] = ptcl->RadiusOfAC*ptcl->RadiusOfAC/ptcl->Mass; // mass weighted

			for (int dim=0; dim<Dim; dim++) {
				PositionSend[i][dim] = ptcl->PredPosition[dim];
				VelocitySend[i][dim] = ptcl->PredVelocity[dim];
			}
		} // endfor copy info

		// calculate the force by sending the particles to GPU
#ifdef time_trace
		_time.reg_gpu.markStart();
#endif
		CalculateAccelerationOnDevice(&ListSize, PositionSend, VelocitySend, AccRegReceive, AccRegDotReceive,
				MdotSend, RadiusOfAC2Send, NumNeighborReceive, ACListReceive, dt);
#ifdef time_trace
		_time.reg_gpu.markEnd();
		_time.reg_gpu.getDuration();

		_time.reg_cpu1.markStart();
#endif

		// Calculate the irregular acceleration components based on neighbors of current regular time.
		for (int i=0; i<ListSize; i++) {

			ptcl = particle[IndexList[i]];  // regular particle in particle list

			for (int dim=0; dim<Dim; dim++) {
				a_tmp[dim]    = 0.;
				adot_tmp[dim] = 0.;
			}


			/*******************************************************
			* Acceleartion correction according to past neighbor
			********************************************************/
			for (Particle *neighbor:ptcl->ACList) {
				CalculateSingleAcceleration(ptcl, neighbor, a_tmp, adot_tmp);
			}

			/*******************************************************
			* Position and velocity correction due to 4th order correction
			********************************************************/
			for (int dim=0; dim<Dim; dim++) {
				
				/*
				a_tmp[dim]    =  AccRegReceive[i][dim]    - a_tmp[dim] ;
				adot_tmp[dim] =  AccRegDotReceive[i][dim] - adot_tmp[dim];
				DFR = (ptcl->a_reg[dim][0] - a_tmp[dim]);
				FRD = adot_tmp[dim];
				SUM = ptcl->a_reg[dim][1] + adot_tmp[dim];
				// do the higher order correcteion
				
				AT3 = 2.0*DFR + DTR*SUM;
				BT2 = -3.0*DFR - DTR*(SUM + FRD);
				
				ptcl->NewPosition[dim] = ptcl->PredPosition[dim] + (0.6*AT3 + BT2)*DTSQ12;
				ptcl->NewVelocity[dim] = ptcl->PredVelocity[dim] + (0.75*AT3 + BT2)*DTR13;

				ptcl->a_reg[dim][2] = (3.0*AT3 + BT2)*DT2;
				ptcl->a_reg[dim][3] = AT3*DT6;
				*/

				da_dt2  = (ptcl->a_reg[dim][0] - a_tmp[dim]   ) / dt2;
				adot_dt = (ptcl->a_reg[dim][1] + adot_tmp[dim]) / dt;

				a2 =  -6*da_dt2 - 2*adot_dt - 2*adot_tmp[dim]/dt;
				a3 = (12*da_dt2 + 6*adot_dt)/dt;
				// note that these higher order terms and lowers have different neighbors

				// 4th order correction
				// save the values in the temporary variables
				ptcl->NewPosition[dim] = ptcl->PredPosition[dim] + a2*dt4/24 + a3*dt5/120;
				ptcl->NewVelocity[dim] = ptcl->PredVelocity[dim] + a2*dt3/6  + a3*dt4/24;

				ptcl->a_reg[dim][2] = a2;
				ptcl->a_reg[dim][3] = a3;
				// reset for future use
				a_tmp[dim]    = 0.;
				adot_tmp[dim] = 0.;
			}


			/*******************************************************
			* Acceleartion correction according to current neighbor
			********************************************************/
			std::cout <<  "MyIndex=" << IndexList[i];
			std::cout <<  "(" << NumNeighborReceive[i] << ")" << std::endl;
			std::cout <<  "NeighborIndex = ";
			for (int j=0;  j<NumNeighborReceive[i]; j++) {
				NeighborIndex = ACListReceive[i][j];  // gained neighbor particle (in next time list)
																								std::cout <<  NeighborIndex << "  ";
				CalculateSingleAcceleration(ptcl, particle[NeighborIndex], a_tmp, adot_tmp);
			} // endfor j1, over neighbor at current time
			std::cout << std::endl;
			for (int dim=0; dim<Dim; dim++) {
				AccIrr[i][dim]           += a_tmp[dim];
				AccIrrDot[i][dim]        += adot_tmp[dim];
				AccRegReceive[i][dim]    -= a_tmp[dim];
				AccRegDotReceive[i][dim] -= adot_tmp[dim];
			} // endfor dim

			// update force
			for (int dim=0; dim<Dim; dim++) {
				ptcl->a_reg[dim][0] = AccRegReceive[i][dim];
				ptcl->a_reg[dim][1] = AccRegDotReceive[i][dim];
				ptcl->a_irr[dim][0] = AccIrr[i][dim];
				ptcl->a_irr[dim][1] = AccIrrDot[i][dim];
				ptcl->a_tot[dim][0] = ptcl->a_reg[dim][0] + ptcl->a_irr[dim][0];
				ptcl->a_tot[dim][1] = ptcl->a_reg[dim][1] + ptcl->a_irr[dim][1];
				// in case
				if (ptcl->NumberOfAC == 0) {
					ptcl->a_tot[dim][2] = ptcl->a_reg[dim][2];
					ptcl->a_tot[dim][3] = ptcl->a_reg[dim][3];
				}
			}

			ptcl->ACList.clear();
			ptcl->NumberOfAC = NumNeighborReceive[i];
			for (int j=0; j<ptcl->NumberOfAC;j++) {
				NeighborIndex = ACListReceive[i][j];  // gained neighbor particle (in next time list)
				ptcl->ACList.push_back(particle[NeighborIndex]);
			}
		} // endfor i, over regular particles
#ifdef time_trace
		_time.reg_cpu1.markEnd();
		_time.reg_cpu1.getDuration();
#endif

		/*******************************************************
		 * Finally update particles
		 ********************************************************/
		for (int i=0; i<ListSize; i++) {
			ptcl = particle[IndexList[i]];
			ptcl->updateParticle();
			ptcl->CurrentTimeReg = NextRegTime;
			ptcl->calculateTimeStepReg();
			if (ptcl->NumberOfAC == 0)
				ptcl->CurrentTimeIrr = NextRegTime;
		}
	//std::cout <<  "Calculation On Device Done ..." << std::endl;







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

	delete[] NumNeighborReceive;
	delete[] AccRegReceive;
	delete[] AccRegDotReceive;
	delete[] AccIrr;
	delete[] AccIrrDot;
	delete[] ACListReceive;

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

	if (dr2 < EPS2) {
		dr2 = EPS2;
	}

	m_r3 = ptcl2->Mass/dr2/sqrt(dr2);

	for (int dim=0; dim<Dim; dim++){
		a[dim]    += m_r3*dx[dim];
		adot[dim] += m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
	}
}
