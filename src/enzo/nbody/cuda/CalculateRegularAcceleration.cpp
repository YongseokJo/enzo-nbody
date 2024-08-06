#include <vector>
#include <iostream>
#include "../global.h"
#include <cmath>
#include "defs.h"
#include <cassert>
#include "cuda_functions.h"


void UpdateNextRegTime(std::vector<Particle*> &particle);
void SendAllParticlesToGPU(double time, std::vector <Particle*> &particle);
void CalculateSingleAcceleration(Particle *ptcl1, Particle *ptcl2, double (&a)[3], double (&adot)[3], int sign);


/*
 *  Purporse: calculate acceleration and neighbors of regular particles by sending them to GPU
 *
 *  Date    : 2024.01.18  by Seoyoung Kim
 *
 */
void CalculateRegAccelerationOnGPU(std::vector<Particle*> RegularList, std::vector<Particle*> &particle){



	// regIds are the list of positions of particles subject to regular force calculation in std::vector list particle

	// variables for opening GPU
	// const int buffer = 10;
	// int numGpuOpen = NNB+buffer;
	//const int NumPtclPerEachCalMax = 2048; // this also caps the number of particles computed each iteration
	const int mpi_rank  = 0; // not effective for now
	int NeighborIndex; // this size should coincide with number of threads
	int ListSize = RegularList.size();
	int *IndexList = new int[ListSize];

	//int NumGpuCal;

	// variables for saving variables to send to GPU
	// only regular particle informations are stored here
	double (*AccRegReceive)[Dim];
	double (*AccRegDotReceive)[Dim];
	double (*AccIrr)[Dim];
	double (*AccIrrDot)[Dim];
	//int (*ACListReceive)[NumNeighborMax];

	//double* PotSend;
	int **ACListReceive;
	int *NumNeighborReceive;
	int MassFlag;


	double a_tmp[Dim]{0}, adot_tmp[Dim]{0};
	double da, dadot;
	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;


	Particle *ptcl;


	double dt       = RegularList[0]->TimeStepReg;
	double new_time = RegularList[0]->CurrentTimeReg + dt;  // next regular time
	ULL new_block = RegularList[0]->CurrentBlockReg + RegularList[0]->TimeBlockReg;  // next regular time
	
	if (new_block != NextRegTimeBlock) {
		if (NextRegTimeBlock == 0) {
			UpdateNextRegTime(particle);
			fprintf(stderr, "First RegularCalculation skips! :CalculateAcceleration.C:105\n");
		}
		else{
			fprintf(stderr, "Something wrong! NextRegTime does not match! :CalculateAcceleration.C:105\n");
			fprintf(stderr, "NextRegTime=%llu, treg[1]=%llu\n", NextRegTimeBlock, new_block);
		}
		return;
	}


	// need to make array to send to GPU
	// allocate memory to the temporary variables
	//PotSend         = new double[ListSize];

	AccRegReceive    = new double[ListSize][Dim];
	AccRegDotReceive = new double[ListSize][Dim];
	AccIrr           = new double[ListSize][Dim];
	AccIrrDot        = new double[ListSize][Dim];

	NumNeighborReceive  = new int[ListSize];
	ACListReceive      = new int*[ListSize];
	for (int i=0; i<ListSize; i++) {
		ACListReceive[i] = new int[MaxNumNeighbor];
		for (int dim=0; dim<Dim; dim++) {
			AccRegReceive[i][dim]    = 0;
			AccRegDotReceive[i][dim] = 0;
			AccIrr[i][dim]           = 0;
			AccIrrDot[i][dim]        = 0;
		}
	}

	// perform the loop twice in order to obtain
	// both the current and predicted acceleration
	// set the current time to 0 and next time to 1
	// and set the time step dt to regular time step




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
		IndexList[i] = RegularList[i]->ParticleOrder;
	} // endfor copy info

	// calculate the force by sending the particles to GPU
#ifdef time_trace
	_time.reg_gpu.markStart();
#endif
	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive, AccRegDotReceive, NumNeighborReceive, ACListReceive);
#ifdef time_trace
	_time.reg_gpu.markEnd();
	_time.reg_gpu.getDuration();

	_time.reg_cpu1.markStart();
#endif

	// Calculate the irregular acceleration components based on neighbors of current regular time.
	for (int i=0; i<ListSize; i++) {

		ptcl = RegularList[i];  // regular particle in particle list

		for (int dim=0; dim<Dim; dim++) {
			a_tmp[dim]    = 0.;
			adot_tmp[dim] = 0.;
		}


		/*******************************************************
		 * Acceleartion correction according to past neighbor
		 ********************************************************/


		/*
		std::cout <<  "MyPID=" <<  ptcl->PID;
		std::cout <<  "(" << NumNeighborReceive[i] << ", ";
		std::cout <<  "(" << ptcl->RadiusOfAC << ")" << std::endl;
		//std::cout <<  "NeighborIndex = ";
		for (int j=0;  j<NumNeighborReceive[i]; j++) {
			NeighborIndex = ACListReceive[i][j];  // gained neighbor particle (in next time list)
																						//std::cout <<  NeighborIndex << "  (" << particle[NeighborIndex]->PID << "), ";
			std::cout <<  particle[NeighborIndex]->PID << ", ";
		}
		std::cout << std::endl;
		*/

		//fprintf(stderr,"%d Neighbor Correction new=%d, old=%d\n", ptcl->PID, NumNeighborReceive[i], ptcl->NumberOfAC);
		int sign = 1;

		/*
		if (NumNeighborReceive[i]>NumNeighborMax) {
			std::cerr <<  "MyPID=" <<  ptcl->PID << ", NN=" << NumNeighborReceive[i] << std::endl;
		}
		*/

		for (int j=0;  j<NumNeighborReceive[i]; j++) {
			NeighborIndex = ACListReceive[i][j];  // gained neighbor particle (in next time list)
			for (auto it = ptcl->ACList.begin(); it != ptcl->ACList.end(); ) {
				//fprintf(stderr,"New PID = %d, Old PID = %d\n", particle[NeighborIndex]->PID, (*it)->getPID());
				if ((*it)->getPID() == particle[NeighborIndex]->PID) {
					it = ptcl->ACList.erase(it);  // Erase the element and update the iterator
					sign = -1;
					break;
				}
				++it;
			}

			if (sign == -1) {
				sign = 1;	
				continue;
			}

			// here, particles are only in new but not in old neighbors
			CalculateSingleAcceleration(ptcl, particle[NeighborIndex], a_tmp, adot_tmp, sign); // so we have to add it 
		}

		// These particles are in the old but not in new, so should be removed for correction.
		for (Particle *neighbor: ptcl->ACList)
			CalculateSingleAcceleration(ptcl, neighbor, a_tmp, adot_tmp, 0);


		/*******************************************************
		 * Position and velocity correction due to 4th order correction
		 ********************************************************/
		dt  = ptcl->TimeStepReg*EnzoTimeStep;  // unit conversion
		dt2 = dt*dt;
		dt3 = dt2*dt;
		dt4 = dt3*dt;
		dt5 = dt4*dt;

		//fprintf(stdout, "PID=%d\n", ptcl->PID);
		for (int dim=0; dim<Dim; dim++) {
			//fprintf(stdout, "a0   =%.3e, a   =%.3e\n", ptcl->a_reg[dim][0], (AccRegReceive[i][dim] + a_tmp[dim]));
			//fprintf(stdout, "aodt0=%.3e, adot=%.3e\n", ptcl->a_reg[dim][1], (AccRegDotReceive[i][dim] + adot_tmp[dim]));
			da_dt2  = (ptcl->a_reg[dim][0] - AccRegReceive[i][dim] - a_tmp[dim]   ) / dt2;
			adot_dt = (ptcl->a_reg[dim][1] + AccRegDotReceive[i][dim] + adot_tmp[dim]) / dt;


			a2 =  -6*da_dt2 - 2*adot_dt - 2*ptcl->a_reg[dim][1]/dt;
			a3 = (12*da_dt2 + 6*adot_dt)/dt;
			// note that these higher order terms and lowers have different neighbors

			//fprintf(stdout, "da_dt2 =%.3e, adot_dt =%.3e, dt=%.3e\n", da_dt2, adot_dt, dt);
			//fprintf(stdout, "a2     =%.3e, a3      =%.3e\n", a2, a3);
			/*
			if (ptcl->PID == 753) {
				fprintf(stderr, "dim=%d, a2=%.3e, a3=%.3e/a0=%.3e, atot=%.3e, a_tmp=%.3e, adot_tmp=%.3e, dt=%.3e\n", 
						dim, a2,a3,ptcl->a_reg[dim][0],AccRegReceive[i][dim],a_tmp[dim],adot_tmp[dim],dt*1e10/1e6);
				fprintf(stderr, "dim=%d, da_dt2=%.3e, adot_dt=%.3e\n", 
						dim, da_dt2, adot_dt);
			}
			*/

			// 4th order correction
			// save the values in the temporary variables
			ptcl->NewPosition[dim] = ptcl->PredPosition[dim] + a2*dt4/24 + a3*dt5/120;
			ptcl->NewVelocity[dim] = ptcl->PredVelocity[dim] + a2*dt3/6  + a3*dt4/24;

			//ptcl->NewPosition[dim] = ptcl->PredPosition[dim];
			//ptcl->NewVelocity[dim] = ptcl->PredVelocity[dim];


			ptcl->a_reg[dim][2] = a2;
			ptcl->a_reg[dim][3] = a3;
			// reset for future use
			a_tmp[dim]    = 0.;
			adot_tmp[dim] = 0.;
		}


		/*******************************************************
		 * Acceleartion correction according to current neighbor
		 ********************************************************/
		//std::cout <<  "MyIndex=" <<  ;

		/*
		if (dt*1e4<1e-8) {
			//std::cout <<  "MyPID=" <<  ptcl->PID;
			//std::cout <<  "(" << NumNeighborReceive[i] << ")" << std::endl;
			//std::cout <<  "NeighborIndex = ";
			for (int j=0;  j<NumNeighborReceive[i]; j++) {
				NeighborIndex = ACListReceive[i][j];  // gained neighbor particle (in next time list)
																							//std::cout <<  NeighborIndex << "  (" << particle[NeighborIndex]->PID << "), ";
				//std::cout <<  particle[NeighborIndex]->PID << ", ";
			}
			//std::cout << std::endl;
		}
		*/

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
			ptcl->a_tot[dim][0] = ptcl->a_reg[dim][0] + ptcl->a_irr[dim][0];
			ptcl->a_tot[dim][1] = ptcl->a_reg[dim][1] + ptcl->a_irr[dim][1];
			// in case
			if (ptcl->NumberOfAC == 0) {
				ptcl->a_tot[dim][2] = ptcl->a_reg[dim][2];
				ptcl->a_tot[dim][3] = ptcl->a_reg[dim][3];
			}
		}

		/*
		ptcl->ACList.clear();
		ptcl->NumberOfAC = NumNeighborReceive[i];
		for (int j=0; j<ptcl->NumberOfAC;j++) {
			NeighborIndex = ACListReceive[i][j];  // gained neighbor particle (in next time list)
			ptcl->ACList.push_back(particle[NeighborIndex]);
		}
		*/
	} // endfor i, over regular particles
#ifdef time_trace
	_time.reg_cpu1.markEnd();
	_time.reg_cpu1.getDuration();

	_time.reg_cpu2.markStart();
#endif

	/*******************************************************
	 * Finally update particles
	 ********************************************************/
	//for (Particle* ptcl: RegularList) {
	for (int i=0; i<ListSize; i++) {
		ptcl = RegularList[i];  // regular particle in particle list
		ptcl->CurrentBlockReg = NextRegTimeBlock;
		ptcl->CurrentTimeReg  = NextRegTimeBlock*time_step;
		ptcl->calculateTimeStepReg();
		ptcl->calculateTimeStepIrr(ptcl->a_tot,ptcl->a_irr);

		/*
		if (ptcl->TimeLevelReg > ptcl->TimeLevelIrr+3 || 
				mag0(ptcl->a_irr)>mag0(ptcl->a_reg)*1e3) {
			//ptcl->updateParticle();
			ptcl->calculateTimeStepIrr(ptcl->a_tot,ptcl->a_irr);
			continue;
		}
		else {
			ptcl->updateParticle();
		}
		*/
		ptcl->updateParticle();
		//ptcl->calculateTimeStepReg();
		//ptcl->calculateTimeStepIrr(ptcl->a_tot,ptcl->a_irr);
		if (ptcl->NumberOfAC == 0) {
			ptcl->CurrentBlockIrr = NextRegTimeBlock;
			ptcl->CurrentTimeIrr = NextRegTimeBlock*time_step;
		}
		/*
		if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0]) {
			fprintf(stdout, "after, myself = %d\n", ptcl->PID);
			fprintf(stdout, "x[0]=%e, a[0]=%e\n", ptcl->Position[0], ptcl->a_tot[0][0]);
			fflush(stdout);
			throw std::runtime_error("CalculateRegAcceleration.cpp::357\n");
			//assert(ptcl->Position[0] ==  ptcl->Position[0]);
		}
		*/
		ptcl->ACList.clear();
		ptcl->NumberOfAC = NumNeighborReceive[i];
		for (int j=0; j<ptcl->NumberOfAC;j++) {
			NeighborIndex = ACListReceive[i][j];  // gained neighbor particle (in next time list)
			ptcl->ACList.push_back(particle[NeighborIndex]);
		}
		ptcl->UpdateRadius();
	}
#ifdef time_trace
	_time.reg_cpu2.markEnd();
	_time.reg_cpu2.getDuration();

	_time.reg_cpu3.markStart();
#endif
	//UpdateNextRegTime(particle);

	//for (Particle* ptcl: RegularList)
		//ptcl->calculateTimeStepIrr(ptcl->a_tot,ptcl->a_irr);
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
	//delete[] PotSend;

	delete[] NumNeighborReceive;
	delete[] AccRegReceive;
	delete[] AccRegDotReceive;
	delete[] AccIrr;
	delete[] AccIrrDot;
	for (int i=0; i<ListSize; i++) {
		delete[] ACListReceive[i];
	}
	delete[] ACListReceive;


#ifdef time_trace
	_time.reg_cpu3.markEnd();
	_time.reg_cpu3.getDuration();
#endif
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

void CalculateSingleAcceleration(Particle *ptcl1, Particle *ptcl2, double (&a)[3], double (&adot)[3], int sign) {
	double dx[Dim], dv[Dim];
	double dr2;
	double dxdv;
	double m_r3;

	if (ptcl1->PID == ptcl2->PID) {
		return;
	}

	dr2  = 0.0;
	dxdv = 0.0;

	/*
	if (ptcl2->PredPosition[0] !=  ptcl2->PredPosition[0]) {
		fprintf(stdout, "target = %d, neighborhood = %d\n", ptcl1->PID, ptcl2->PID);
		fprintf(stdout, "x[0] = %e\n", ptcl2->PredPosition[0]);
		//fflush(stdout);
		assert(ptcl2->PredPosition[0] ==  ptcl2->PredPosition[0]);
	}
	*/

	for (int dim=0; dim<Dim; dim++) {
		dx[dim] = ptcl2->PredPosition[dim] - ptcl1->PredPosition[dim];
		dv[dim] = ptcl2->PredVelocity[dim] - ptcl1->PredVelocity[dim];
		dr2    += dx[dim]*dx[dim];
		dxdv   += dx[dim]*dv[dim];
	}

	/*
	if (dr2 < EPS2) {
		dr2 = EPS2;
	}
	*/

	m_r3 = ptcl2->Mass/dr2/sqrt(dr2);
	if (sign == 0)
		m_r3 *= -1;

	for (int dim=0; dim<Dim; dim++){
		a[dim]    += m_r3*dx[dim];
		adot[dim] += m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
	}
}
