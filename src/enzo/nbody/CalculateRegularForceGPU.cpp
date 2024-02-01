#include <vector>
#include <iostream>
#include "global.h"
#include <cmath>
#include "defs.h"
#include "cuda_functions.h"


void sendAllParticlesToGpu(double sendTime, std::vector <Particle*> &particle);
void calculatePredAccelerationIJ(Particle *ptcl1, Particle *ptcl2, double (&a)[3], double (&adot)[3]);


/*
 *  Purporse: calculate acceleration and neighbors of regular particles by sending them to GPU
 *
 *  Date    : 2024.01.18  by Seoyoung Kim
 *
 */

void calRegAccelerationOnGPU(std::vector<int> &regIds, std::vector<Particle*> &particle){ 
	// regIds are the list of positions of particles subject to regular force calculation in std::vector list particle

	// variables for opening GPU

	int regSize;
	int numGpuSend = 1024;
	int numGpuCal;
	int mpiRank = 0;

	// variables for saving variables to send to GPU
	// only regular particle informations are stored here

	double* massSend;
    double(*positionSend)[Dim];
	double(*velocitySend)[Dim];

	double* r2OfACSend;
	double* stepRegSend;

	double(*accRegTmp[2])[Dim];
	double(*accRegDotTmp[2])[Dim];
	double(*accIrrTmp[2])[Dim];
	double(*accIrrDotTmp[2])[Dim];

	double* potGpu;
	int** AClistGpu[2];
	int massFlag;

	// the current and next times to calculate the regular force on

	double treg[2];  // 0 is current time and 1 is next time

	// temporary variables saving the calculated acceleration values

	double a_tmp[Dim], adot_tmp[Dim];

	// variables for calculating the higher order acceleration derivatives

	double a2, a3, da_dt2, adot_dt, dt, dt2, dt3, dt4;

	// extra variables

	int ACnumi1, ACnumi2;

	Particle *regPtcl;
	int ACjid;
	int j2;


	// need to make array to send to GPU
    // first check the size of the array to be made

	regSize = regIds.size();
    

	// allocate memory to the temporary variables

	massSend = new double[regSize];
	positionSend = new double[regSize][Dim];
	velocitySend = new double[regSize][Dim];

	r2OfACSend = new double[regSize];
	stepRegSend = new double[regSize];
	potGpu = new double[regSize];

	for (int p=0; p<2; i++) {

		accRegTmp[p] = new double[regSize][Dim];
		accRegDotTmp[p] = new double[regSize][Dim];
		accIrrTmp[p] = new double[regSize][Dim];
		accIrrDotTmp[p] = new double[regSize][Dim];

		AClistGpu[p] = new int*[regSize];

		for (int i=0; i<(NNB); i++) {
			AClistGpu[i] = new int[maxNeighborNum];
		}
	}

	// perform the loop twice in order to obtain
	// both the current and predicted values of acceleration

	// set the values of current time to 0 and next time to 1
	// and set the time step dt to regular time step



	treg[0] = particle[regIds[0]]->CurrentTimeReg;  // current regular time
	dt = particle[regIds[0]]->TimeStepReg; 
	treg[1] = treg[0] + dt;  // next regular time

	std::cout <<  "SENDING REGULAR PARTICLES TO GPU" << std::flush;

	for (int p=0; p<2; p++) {

		// send information of all the particles to GPU

		sendAllParticlesToGpu(treg[p], particle);

		// copy the data of regular particles to the arrays to be sent
		// predicted positions and velocities should be sent
		// but predictions are already done when sending all particles, so no need for duplicated calculation

		for (int i=0; i<regSize; i++) {

			massSend[i] = particle[(regIds[i])]->Mass;
			r2OfACSend[i] = (particle[(regIds[i])]->RadiusOfAC)*(particle[(regIds[i])]->RadiusOfAC);

			for (int dim=0; dim<Dim; dim++) {
				positionSend[i][dim] = particle[(regIds[i])]->PredPosition[dim];
				velocitySend[i][dim] = particle[(regIds[i])]->PredVelocity[dim];
			}

			stepRegSend[i] = particle[(regIds[i])]->TimeStepReg;

		}


		// calculate the force by sending the particles to GPU in multiples of 1024

		for (int i=0; i<regSize; i+=numGpuSend) {
			
			numGpuCal = std::min(1024,(NNB-i));

			gpunb_regf_(&numGpuCal,&r2OfACSend[i],&stepRegSend[i],&positionSend[i],&velocitySend[i],
						&accRegTmp[p][i],&accRegDotTmp[p][i],&potGpu[i],&maxNeighborNum,&maxNeighborNum,AClistGpu[p][i],&massFlag);


			// CALCULATE THE IRREGULAR ACCELERATION COMPONENTS
			// calculate the irregular acceleration components based on neighbors of current regular time.

			for (int i2=i; i2<(i+numGpuCal); i2++) {

				regPtcl = particle[(regIds[i])];  // regular particle in particle list

				for (int j1=1;  j1<(ACnumi1+1); j1++) {

					ACjid = AClistGpu[0][i2][j1];  // gained neighbor particle (in next time list)

					calculatePredAccelerationIJ(regPtcl, particle[ACjid], a_tmp, adot_tmp);
					
					for (int dim=0; dim<Dim; dim++) {
						accRegTmp[p][i2][dim] += a_tmp[dim];
						accRegDotTmp[p][i2][dim] += adot_tmp[dim];						
					}

				}
			}

		} // loop of total calculation ends

	} // 0th, 1st derivative acceleration calculation + neighbors search of regular particles ends



	
	// now we should calculate the correction terms using the information we have

	// we need to 1. make corrections in regular force for neighbor changes on current and next time steps 2. calculate the 2nd and 3rd derivatives using the hermite method
	// we will use the neighbors of the current steps to obtain the 2nd, 3rd derivatives of regular force

	for (int i=0; i<regSize; i++) {

		regPtcl = particle[(regIds[i])];  // regular particle in particle list


		// CORRECTIONS PERFORMED DUE TO NEIGHBOR CHANGES

		// if the neighbor numbers are different in current and next regular times, then correct the calculated regular forces
		// since the last predictions performed are on the next time steps, there is no need to predict the positions again for calculating the correction terms

		ACnumi1 = AClistGpu[0][i][0];  // neighbor number of current time
		ACnumi2 = AClistGpu[1][i][0];  // neighbor number of the following time

		if (ACnumi1 > ACnumi2) { // neighbor particles are lost at next time

			j2 = 1;

			for (int j1=1;  j1<(ACnumi2+1); j1++) { // loop over the neighbor list of next regular time

				while (AClistGpu[1][i][j1] != AClistGpu[0][i][j2]) {  // if a member not on list2 is found ...

					// since neighbor particles are lost in next time, their force contributions are included in the regular force
					// therefore, we need to subtract the forces of neighbors from regular force

					ACjid = AClistGpu[0][i][j2]; // lost neighbor particle (in current time list)

					calculatePredAccelerationIJ(regPtcl, particle[ACjid], a_tmp, adot_tmp);
					
					if (j2>=j1){
						for (int dim=0; dim<Dim; dim++) {
							accRegTmp[1][i][dim] -= a_tmp[dim];
							accRegDotTmp[1][i][dim] -= adot_tmp[dim];						
						}
					}

					j2++; // search the current neighbor list (with j2 index) until next member is found

					if (j2>(ACnumi1+1)) {
						std::cout <<  "WARNING: there are unrecognized gained neighbor particles in regular acceleration routines" << std::flush;
						j2 = 1;
					}
				}
				j2++;
			} // loop for neighbor list of next regular time ends

			// additional loop in case that neighbors lost are last particles in the particle list

			for (int j1=j2; j1<(ACnumi1+1); j1++) {
				
				ACjid = AClistGpu[0][i][j1];  // lost neighbor particle (in current time list)

				calculatePredAccelerationIJ(regPtcl, particle[ACjid], a_tmp, adot_tmp);
				
				for (int dim=0; dim<Dim; dim++) {
					accRegTmp[1][i][dim] -= a_tmp[dim];
					accRegDotTmp[1][i][dim] -= adot_tmp[dim];						
				}
			} // end of additional loop 

		} else if (ACnumi1 < ACnumi2) { // new neighbor particles are added at next time
		
			j2 = 1;

			for (int j1=1;  j1<(ACnumi1+1); j1++) { // loop over the neighbor list of current time

				while (AClistGpu[0][i][j1] != AClistGpu[1][i][j2]) {  // if a member not on list1 is found ...

					// since neighbor particles are gained in next time, their force contributions are left out in the regular force
					// therefore, we need to add the forces of neighbors to regular force

					ACjid = AClistGpu[1][i][j2];  // gained neighbor particle (in next time list)

					calculatePredAccelerationIJ(regPtcl, particle[ACjid], a_tmp, adot_tmp);
					
					if (j2>=j1){
						for (int dim=0; dim<Dim; dim++) {
							accRegTmp[1][i][dim] += a_tmp[dim];
							accRegDotTmp[1][i][dim] += adot_tmp[dim];						
						}
					}

					j2++; // search the current neighbor list (with j2 index) until next member is found

					if (j2>(ACnumi2+1)) {
						std::cout <<  "WARNING: there are unrecognized lost neighbor particles in regular acceleration routines" << std::flush;
						j2 = 1;
					}
				}
				j2++;

			} // loop for neighbor list of next regular time ends

			// additional loop in case that neighbors gained are last particles in the particle list

			for (int j1=j2; j1<(ACnumi2+1); j1++) {
				
				ACjid = AClistGpu[1][i][j1];  // gained neighbor particle (in next time list)

				calculatePredAccelerationIJ(regPtcl, particle[ACjid], a_tmp, adot_tmp);
				
				for (int dim=0; dim<Dim; dim++) {
					accRegTmp[1][i][dim] += a_tmp[dim];
					accRegDotTmp[1][i][dim] += adot_tmp[dim];						
				}
			} // end of additional loop 
		} // corrections on change of neighbors are finished


		// CALCULATE AND SAVE THE 0TH, 1ST, 2ND, 3RD DERIVATIVES OF REGULAR AND IRREGULAR FORCE

		dt2 = dt*dt;
		dt3 = dt2*dt;
		dt4 = dt3*dt;


		// calculated the final corrected forces
		for (int dim=0; dim<Dim; dim++) {
			da_dt2  = (   accRegTmp[0][i][dim] - accRegTmp[1][i][dim]   ) / dt2;
			adot_dt = (accRegDotTmp[0][i][dim] + accRegDotTmp[1][i][dim]) / dt;

			a2 =  -6*da_dt2 - 2*adot_dt - 2*accRegDotTmp[0][i][dim]/dt;
			a3 = (12*da_dt2 + 6*adot_dt)/dt;
			//Position[dim] += a2*dt4/24 + a3*dt4*dt/120;
			//Velocity[dim] += a2*dt3/6  + a3*dt4/24;

			regPtcl->a_reg[dim][0] = accRegTmp[0][i][dim];
			regPtcl->a_reg[dim][1] = accRegDotTmp[0][i][dim];
			regPtcl->a_reg[dim][2] = a2;
			regPtcl->a_reg[dim][3] = a3;
		}

		for (int dim=0; dim<Dim; dim++) {
			da_dt2  = (   accIrrTmp[0][i][dim] - accIrrTmp[1][i][dim]   ) / dt2;
			adot_dt = (accIrrDotTmp[0][i][dim] + accIrrDotTmp[1][i][dim]) / dt;

			a2 =  -6*da_dt2  - 2*adot_dt - 2*accIrrDotTmp[0][i][dim]/dt;
			a3 = (12*da_dt2 + 6*adot_dt)/dt;
			//Position[dim] += a2*dt4/24 + a3*dt4*dt/120;
			//Velocity[dim] += a2*dt3/6  + a3*dt4/24;

			regPtcl->a_irr[dim][0] = accIrrTmp[0][i][dim];
			regPtcl->a_irr[dim][1] = accIrrDotTmp[0][i][dim];
			regPtcl->a_irr[dim][2] = a2;
			regPtcl->a_irr[dim][3] = a3;

			regPtcl->a_tot[dim][0] = regPtcl->a_reg[dim][0] + regPtcl->a_irr[dim][0];
			regPtcl->a_tot[dim][1] = regPtcl->a_reg[dim][1] + regPtcl->a_irr[dim][1];
			regPtcl->a_tot[dim][2] = regPtcl->a_reg[dim][2] + regPtcl->a_irr[dim][2];
			regPtcl->a_tot[dim][3] = regPtcl->a_reg[dim][3] + regPtcl->a_irr[dim][3];
		}

	} // correction and calculation of higher orders finished



	// free all temporary variables

	delete[] massSend;
    delete[] positionSend;
	delete[] velocitySend;

	delete[] r2OfACSend;
	delete[] stepRegSend;

	delete[] accRegTmp;
	delete[] accRegDotTmp;

	delete[] accIrrTmp;
	delete[] accIrrDotTmp;

	delete[] potGpu;
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


void sendAllParticlesToGpu(double sendTime, std::vector <Particle*> &particle) {

	// variables for saving variables to send to GPU

	double* massSendAll;
    double(*positionSendAll)[Dim];
	double(*velocitySendAll)[Dim];

    
    // allocate memory to the temporary variables

	massSendAll = new double[NNB];
	positionSendAll = new double[NNB][Dim];
	velocitySendAll = new double[NNB][Dim];


	// copy the data of particles to the arrays to be sent

	for (int i=0; i<NNB; i++) {

		massSendAll[i] = particle[i]->Mass;

		// before sending them, predict their positions and velocities again just in case

		particle[i]->predictParticleSecondOrder(sendTime);

		for (int dim=0; dim<Dim; dim++) {
			positionSendAll[i][dim] = particle[i]->PredPosition[dim];
			velocitySendAll[i][dim] = particle[i]->PredVelocity[dim];
		}
	}

	// send the arrays to GPU

	gpunb_send_(&NNB, massSendAll, positionSendAll, velocitySendAll);


    // free the temporary variables
    
	delete[] massSendAll;
    delete[] positionSendAll;
	delete[] velocitySendAll;

}


void calculatePredAccelerationIJ(Particle *ptcl1, Particle *ptcl2, double (&a)[3], double (&adot)[3]) {

	double dx[Dim], dv[Dim];
	double rij2;
	double dxdv;

	double m_r3;

	rij2 = 0.0;
	dxdv = 0.0;

	for (int dim=0; dim<Dim; dim++) {
		dx[dim] = ptcl2->PredPosition - ptcl1->PredPosition;
		dv[dim] = ptcl2->PredVelocity - ptcl1->PredVelocity;
		rij2 += dx[dim]*dx[dim];
		dxdv += dx[dim]*dv[dim];
	}

	m_r3 = ptcl2->Mass/rij2/sqrt(rij2);

	for (int dim=0; dim<Dim; dim++){
		a[dim]    = m_r3*dx[dim];
		adot[dim] = m_r3*(dv[dim] - 3*dx[dim]*dxdv/rij2);
	}

}