#include <vector>
#include <iostream>
#include "global.h"
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

void sendAllParticlesToGpu(std::vector <Particle*> &particle);

void calAllAccelerationOnGPU(std::vector<Particle*> &particle){

	// variables for opening GPU

	int numGpuOpen = NNB+10;
	int numGpuSend = 1024;
	int numGpuCal;
	int mpiRank = 0;

	// variables for saving variables to send to GPU

	double* massSend;
    double(*positionSend)[Dim];
	double(*velocitySend)[Dim];

	double* r2OfACSend;
	double stepRegTmp, ri2;
	double* stepRegSend;

	double(*accGpu)[Dim];
	double(*accDotGpu)[Dim];
	double* potGpu;
	int** AClistGpu;
	int massFlag;

	// temporary variables for calculating the irregular force

	double dx[Dim];
	double dv[Dim];
	double rij2,dr2i,dr3i,drdv;

	// extra variables

	bool neighborOK;
	int ACnumi2;
	int ACjid;



	// first let's open the GPU

	gpunb_open_(&numGpuOpen, &mpiRank);


	// need to make array to send to GPU

	// allocate memory to the temporary variables

	massSend = new double[NNB];
	positionSend = new double[NNB][Dim];
	velocitySend = new double[NNB][Dim];

	r2OfACSend = new double[NNB];
	stepRegSend = new double[NNB];

	accGpu = new double[NNB][Dim];
	accDotGpu = new double[NNB][Dim];
	potGpu = new double[NNB];


	AClistGpu = new int*[NNB];
	for (int i=0; i<(NNB); i++) {
		AClistGpu[i] = new int[maxNeighborNum];
	}


	// copy the data of particles to the arrays to be sent

	for (int i=0; i<NNB; i++) {

		ri2 = 0;

		massSend[i] = particle[i]->Mass;
		r2OfACSend[i] = (particle[i]->RadiusOfAC)*(particle[i]->RadiusOfAC);

		for (int dim=0; dim<Dim; dim++) {
			positionSend[i][dim] = particle[i]->Position[dim];
			velocitySend[i][dim] = particle[i]->Velocity[dim];
			ri2 += particle[i]->Position[dim]*particle[i]->Position[dim];
		}

		stepRegTmp = 1.0/8.0*sqrt(1.0 + ri2);
		stepRegSend[i] = std::min(stepRegTmp,1.0);

	}

	// send the arrays to GPU

	gpunb_send_(&NNB, massSend, positionSend, velocitySend);


	// calculate the force by sending the particles to GPU in multiples of 1024

	for (int i=0; i<NNB; i+=numGpuSend) {
		
		numGpuCal = std::min(1024,(NNB-i));

		gpunb_regf_(&numGpuCal,&r2OfACSend[i],&stepRegSend[i],&positionSend[i],&velocitySend[i],
					&accGpu[i],&accDotGpu[i],&potGpu[i],&maxNeighborNum,&maxNeighborNum,AClistGpu[i],&massFlag);


		// copy the values of regular forces and neighbors obtained in GPU to particles
		// and also calculate the irregular forces in the process

		for (int i2=i; i2<(i+1024); i2++){

			for (int dim=0; dim<Dim; dim++) {
				
				particle[i2]->a_reg[dim][0] += accGpu[i2][dim];
				particle[i2]->a_reg[dim][1] += accDotGpu[i2][dim];

			}

            ACnumi2 = AClistGpu[i][0];

            // if there are neighbors, save them to neighbor list and calculate the irregular force
			if (ACnumi2 > 0){

                particle[i2]->NumberOfAC = ACnumi2;

                for (int j=1; j<(ACnumi2+1); j++) {

                    ACjid = AClistGpu[i][j];
                    particle[i2]->ACList.push_back(particle[ACjid]);

                    rij2 = 0.0;
				    drdv = 0.0;

				    for (int dim=0; dim<Dim; dim++) {
				    	dx[dim] = particle[ACjid]->Position - particle[i2]->Position;
			    		dv[dim] = particle[ACjid]->Velocity - particle[i2]->Velocity;
			    		rij2 += dx[dim]*dx[dim];
			    		drdv += dx[dim]*dv[dim];
			    	}

			    	dr2i = 1.0/rij2;
			    	dr3i = particle[ACjid]->Mass*dr2i*sqrt(dr2i);
				
			    	for (int dim=0; dim<Dim; dim++) {
			    		particle[i2]->a_irr[dim][0] += dx[dim]*dr3i;
			    		particle[i2]->a_irr[dim][1] += (dv[dim]-dx[dim]*drdv)*dr3i;
			    	}
                } // loop on neighbors end

            // in case of no neighbors, just set the neighbor number to 0 just in case. 

            } else {
                particle[i2]->NumberOfAC = 0;
            } // if statement on neighbor number ends


			// copy the values to other values as well

			for (int dim=0; dim<3; dim++){
			
				particle[i2]->a_tot[dim][0] = particle[i2]->a_reg[dim][0] + particle[i2]->a_irr[dim][0];
				particle[i2]->a_tot[dim][1] = particle[i2]->a_reg[dim][1] + particle[i2]->a_irr[dim][1];

			}
		} // saving the 1024 results from GPU to local class ends

	} // loop of total calculation ends


	// free all temporary variables

	delete[] massSend;
    delete[] positionSend;
	delete[] velocitySend;

	delete[] r2OfACSend;
	delete[] stepRegSend;

	delete[] accGpu;
	delete[] accDotGpu;
	delete[] potGpu;
	delete[] AClistGpu;

	// close GPU

	gpunb_close_();

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


void calListAccelerationOnGPU(std::vector<int> &regIds, std::vector<Particle*> &particle){

	// variables for opening GPU

	int regSize;
	int numGpuSend = 1024;
	int numGpuCal;
	int mpiRank = 0;

	// variables for saving variables to send to GPU

	double* massSend;
    double(*positionSend)[Dim];
	double(*velocitySend)[Dim];

	double* r2OfACSend;
	double* stepRegSend;

	double(*accGpu)[Dim];
	double(*accDotGpu)[Dim];
	double* potGpu;
	int** AClistGpu;
	int massFlag;

	// temporary variables for calculating the irregular force

	double dx[Dim];
	double dv[Dim];
	double rij2,dr2i,dr3i,drdv;

	// extra variables

	int ACnumi2;
	int ACjid;

	int i2reg;



	// need to make array to send to GPU

    // first check the size of the array to be made

	regSize = regIds.size();
    

	// allocate memory to the temporary variables

	massSend = new double[regSize];
	positionSend = new double[regSize][Dim];
	velocitySend = new double[regSize][Dim];

	r2OfACSend = new double[regSize];
	stepRegSend = new double[regSize];

	accGpu = new double[regSize][Dim];
	accDotGpu = new double[regSize][Dim];
	potGpu = new double[regSize];


	AClistGpu = new int*[regSize];
	for (int i=0; i<(NNB); i++) {
		AClistGpu[i] = new int[maxNeighborNum];
	}


	// copy the data of particles to the arrays to be sent
	// predicted positions and velocities are sent

	for (int i=0; i<regSize; i++) {

		rij2 = 0;

		massSend[i] = particle[(regIds[i])]->Mass;
		r2OfACSend[i] = (particle[(regIds[i])]->RadiusOfAC)*(particle[(regIds[i])]->RadiusOfAC);

		for (int dim=0; dim<Dim; dim++) {
			positionSend[i][dim] = particle[(regIds[i])]->PredPosition[dim];
			velocitySend[i][dim] = particle[(regIds[i])]->PredVelocity[dim];
		}

		stepRegSend[i] = particle[(regIds[i])]->TimeStepReg;

	}

	// send the arrays to GPU

    sendAllParticlesToGpu(particle);

	// calculate the force by sending the particles to GPU in multiples of 1024

	for (int i=0; i<regSize; i+=numGpuSend) {
		
		numGpuCal = std::min(1024,(NNB-i));

		gpunb_regf_(&numGpuCal,&r2OfACSend[i],&stepRegSend[i],&positionSend[i],&velocitySend[i],
					&accGpu[i],&accDotGpu[i],&potGpu[i],&maxNeighborNum,&maxNeighborNum,AClistGpu[i],&massFlag);


		// copy the values of regular forces and neighbors obtained in GPU to particles
		// and also calculate the irregular forces in the process

		for (int i2=i; i2<(i+1024); i2++){

			i2reg = (regIds[i2]);

			// need to clear the AClist and the number of neighbors before copying the new values

			particle[i2reg]->NumberOfAC = 0;
			particle[i2reg]->ACList.clear();

			for (int dim=0; dim<Dim; dim++) {
				
				particle[i2reg]->a_reg[dim][0] += accGpu[i2][dim];
				particle[i2reg]->a_reg[dim][1] += accDotGpu[i2][dim];

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
			    		rij2 += dx[dim]*dx[dim];
			    		drdv += dx[dim]*dv[dim];
			    	}

			    	dr2i = 1.0/rij2;
			    	dr3i = particle[ACjid]->Mass*dr2i*sqrt(dr2i);
				
			    	for (int dim=0; dim<Dim; dim++) {
			    		particle[i2reg]->a_irr[dim][0] += dx[dim]*dr3i;
			    		particle[i2reg]->[dim][1] += (dv[dim]-dx[dim]*drdv)*dr3i;
			    	}
                } // loop on neighbors end

            // in case of no neighbors, just set the neighbor number to 0 just in case. 

            } else {
                particle[i2]->NumberOfAC = 0;
            } // if statement on neighbor number ends


			// copy the values to other values as well

			for (int dim=0; dim<3; dim++){
			
				particle[i2reg]->a_tot[dim][0] = particle[i2reg]->a_reg[dim][0] + particle[i2reg]->a_irr[dim][0];
				particle[i2reg]->a_tot[dim][1] = particle[i2reg]->a_reg[dim][1] + particle[i2reg]->a_irr[dim][1];

			}
		} // saving the 1024 results from GPU to local class ends

	} // loop of total calculation ends


	// free all temporary variables

	delete[] massSend;
    delete[] positionSend;
	delete[] velocitySend;

	delete[] r2OfACSend;
	delete[] stepRegSend;

	delete[] accGpu;
	delete[] accDotGpu;
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


void sendAllParticlesToGpu(std::vector <Particle*> &particle) {

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
