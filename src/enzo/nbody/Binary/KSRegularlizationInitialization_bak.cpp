#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"

void direct_sum(double *x, double *v, double r2, double vx,
		        double mass, double (&a)[3], double (&adot)[3]);



//       /////////////////////        //
//       /////////////////////        //
//       CALCULATEACCELERATION        //
//       /////////////////////        //
//       /////////////////////        //


void CalculateKSAcceleration(Particle* ptclI, Particle* ptclJ, Particle* ptclCM, std::vector<Particle*> &particle, double current_time) {

	int j=0;
	Particle* ptcl1;
	double x[Dim], v[Dim], a[Dim], adot[Dim];
	double r2 = 0;
	double vx = 0;
	double v2 = 0;
	double m_r3, vx_r2, v2x2_r4,v2_r2__ax_r2__v2x2_r4, a2dot, a3dot;
	double A,B;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]    = 0.;
		v[dim]    = 0.;
		a[dim]    = 0.;
		adot[dim] = 0.;
	}

	//std::cout << "nbody+: Entering CalculateInitialAcceleration  ..." << std::endl;


	for (int i=0; i<2; i++) {
	if (i==0) {
		ptcl1 = ptclI;
	} else {
		ptcl1 = ptclJ;
	}	

        // updated the predicted positions and velocities just in case
            
        ptcl1->predictParticleSecondOrder(current_time);

        for (Particle *ptcl2: particle) {

            r2 = 0;
            vx = 0;
            v2 = 0;

	    if (ptcl1 == ptcl2) {
		continue;
	    }
            // if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
            ptcl2->predictParticleSecondOrder(current_time);

            for (int dim=0; dim<Dim; dim++) {
                x[dim] = ptcl1->PredPosition[dim] - ptcl2->PredPosition[dim];
                v[dim] = ptcl1->PredVelocity[dim] - ptcl2->PredVelocity[dim];
                r2    += x[dim]*x[dim];
                vx    += v[dim]*x[dim];
                v2    += v[dim]*v[dim];
            }

            //r2  += EPS2;
            m_r3 = ptcl2->Mass/r2/sqrt(r2); 

            for (int dim=0; dim<Dim; dim++) {
                // Calculate 0th and 1st derivatives of acceleration
                if ((ptcl1->NumberOfAC==0) || (ptcl2 != ptcl1->ACList[j])) {
                    ptcl1->a_reg[dim][0] += m_r3*x[dim];
                    ptcl1->a_reg[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
                }
                else {
                    ptcl1->a_irr[dim][0] += m_r3*x[dim];
                    ptcl1->a_irr[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
                    j++;
                }
            }
    
        } // end of loop for ptcl2 (full particle)


        // update total acceleration as well

        for (int dim=0; dim<Dim; dim++)	 {
            for (int order=0; order<HERMITE_ORDER; order++) {
                ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order]; 
            }
	}

    }// end of loop for pair particle, ptclI and ptclJ


    // copy the calculated values to CM particle
    // and initialize the 3rd and 4th order derivatives just in case

    for (int dim=0; dim<Dim; dim++)	 {
		for (int order=0; order<HERMITE_ORDER; order++) {
            if (order<2) {
		ptclCM->a_reg[dim][order] = (ptclI->a_reg[dim][order]*ptclI->Mass + ptclJ->a_reg[dim][order]*ptclJ->Mass)/(ptclCM->Mass); 
                ptclCM->a_irr[dim][order] = (ptclI->a_irr[dim][order]*ptclI->Mass + ptclJ->a_irr[dim][order]*ptclJ->Mass)/(ptclCM->Mass); 
                ptclCM->a_tot[dim][order] = (ptclI->a_tot[dim][order]*ptclI->Mass + ptclJ->a_tot[dim][order]*ptclJ->Mass)/(ptclCM->Mass); 
            } else {
                ptclCM->a_reg[dim][order] = 0.0; 
                ptclCM->a_irr[dim][order] = 0.0;
                ptclCM->a_tot[dim][order] = 0.0;
            }
	}
    }


    // updated the predicted positions and velocities just in case
            
    //ptcl1 = ptclCM;
    //ptcl1->predictParticleSecondOrder(current_time);


    for (Particle *ptcl2: particle) {

            r2 = 0;
            vx = 0;
            v2 = 0;

	    if ((ptcl2 == ptclI)||(ptcl2==ptclJ)) {
		continue;
	    }

            // if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
            ptcl2->predictParticleSecondOrder(current_time);

            for (int dim=0; dim<Dim; dim++) {
                x[dim] = ptclCM->PredPosition[dim] - ptcl2->PredPosition[dim];
                v[dim] = ptclCM->PredVelocity[dim] - ptcl2->PredVelocity[dim];
                r2    += x[dim]*x[dim];
                vx    += v[dim]*x[dim];
                v2    += v[dim]*v[dim];
            }

            //r2  += EPS2;
            m_r3 = ptcl2->Mass/r2/sqrt(r2);

            vx_r2   = vx/r2;
            v2x2_r4 = vx_r2*vx_r2;
            v2_r2__ax_r2__v2x2_r4 = (v2+x[0]*a[0]+x[1]*a[1]+x[2]*a[2])/r2+v2x2_r4;
            A = (9*(v[0]*a[0]+v[1]*a[1]+v[2]*a[2]) + 3*(x[0]*adot[0]+x[1]*adot[1]+x[2]*adot[2]))/r2\
                    +3*vx_r2*(3*v2_r2__ax_r2__v2x2_r4 - 4*v2x2_r4);

            for (int dim=0; dim<Dim; dim++) {
                B     = v[dim] - 3*x[dim]*vx_r2;
                a2dot = (v[dim] - 6*B*vx_r2                 - 3*v2_r2__ax_r2__v2x2_r4*x[dim])*m_r3;
                a3dot = (a[dim] - 9*B*v2_r2__ax_r2__v2x2_r4 - A*x[dim]                      )*m_r3\
                                - 9*vx_r2*a2dot;
                if ((ptclCM->NumberOfAC==0) || (ptcl2 != ptcl1->ACList[j])) {
                    ptclCM->a_reg[dim][2] += a2dot;
                    ptclCM->a_reg[dim][3] += a3dot;
                }
                else {
                    ptclCM->a_irr[dim][2] += a2dot;
                    ptclCM->a_irr[dim][3] += a3dot;
                    j++;
                }
            } // endfor dim

    } // end of loop for ptcl2 (full particle)

    for (int dim=0; dim<Dim; dim++)	 {
		for (int order=2; order<HERMITE_ORDER; order++) {
	            ptclCM->a_tot[dim][order] = (ptclCM->a_reg[dim][order] + ptclCM->a_reg[dim][order]*ptclJ->Mass); 
		}
	}

    // save the regular and irregular time steps as well
    




}



//        ///////////////////        //
//        ///////////////////        //
//           ISKSCANDIDATE           //
//        ///////////////////        //
//        ///////////////////        //



// made 2024.02.19 by Seoyoung Kim

// from particles with shortest time steps...
// see if they meet the KS regularlization conditions
// reference - search.f

void Particle::isKSCandidate(double next_time) {

    // temporary calculation variables

    double x[Dim],v[Dim];
    double r2;

    double rmin;
    Particle* minPtcl;

    // t
    int numberOfPairCandidate;



    // initialize variables and count numbers just in case

    r2 = 0.0;
    rmin = 1e8;
    numberOfPairCandidate = 0;


    // predict the particle position to obtain more information
    // particle regularlized if the conditions are satisfied at a future time


    this->predictParticleSecondOrder(next_time);

    // need to consider case when c.m particles are the cause of the small steps
    // need to be added later - CHECK

    for (Particle* ptcl: ACList) {

        // if particle time step is too large, skip
        // if the neighbor step is larger then 8 times of the candidate particle, then skip

	r2 = 0.0;
        
        if ((ptcl->TimeStepIrr) > (8*TimeStepIrr)) {
            continue;
        }


        // find out what the paired particle is

        ptcl->predictParticleSecondOrder(next_time);

        for (int dim=0; dim<Dim; dim++) {
				
        	// calculate position and velocity differences
		x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
		v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];
		
		// calculate the square of radius and inner product of r and v for each case
		r2 += x[dim]*x[dim];
        }

        // find out the close particles

        if (r2<KSDistance) {

            numberOfPairCandidate += 1;

            if (r2<rmin) {
                rmin = r2;
                minPtcl = ptcl;
            }
        }
    }


    // save the KS pair information

    if (numberOfPairCandidate>0) {
        isBinary = true;
        BinaryPairParticle = minPtcl;
    }

        
}




//        ///////////////////        //
//        ///////////////////        //
//           STARTNEWKSREG           //
//        ///////////////////        //
//        ///////////////////        //



// start new KS regularlization



//        ///////////////////        //
//        ///////////////////        //
//        NEWKSINITIALIZATION        //
//        ///////////////////        //
//        ///////////////////        //



// initialize conditions for new KS regularlization
// reference: ksinit.F

void NewKSInitialization(Particle* ptclI, std::vector<Particle*> &particle, double current_time) {

    // basic variables for calculation

    Particle *ptclJ;
    Particle *ptclCM;
    Binary *ptclBin;
    std::vector<Particle*> KSNeighborList;

    int ptclIIndex;
    int ptclJIndex;

    std::cout <<"Starting Routine NewKSInitialization" << std::endl;

    // BASIC DEFINITIONS

    // define the pair particle for particle I

    ptclJ = ptclI->BinaryPairParticle;

    //std::cout << "Predicting Particle Positions" << std::endl;

    ptclI->predictParticleSecondOrder(current_time);
    ptclJ->predictParticleSecondOrder(current_time);


    // need to put option if there aren't any close neighbors


    // define the new center of mass particle

    //std::cout << "Make new particle and binary information" << std::endl;

    ptclCM = new Particle;
    ptclBin = new Binary;

    // calculate the values of the center of mass particle
    // and save it to the new center of mass particle

    ptclCM->Mass = ptclI->Mass + ptclJ->Mass;

    for (int dim=0; dim<Dim; dim++) {
        ptclCM->Position[dim] = (ptclI->PredPosition[dim]*ptclI->Mass + ptclJ->PredPosition[dim]*ptclJ->Mass)/ptclCM->Mass;
        ptclCM->Velocity[dim] = (ptclI->PredVelocity[dim]*ptclI->Mass + ptclJ->PredVelocity[dim]*ptclJ->Mass)/ptclCM->Mass;
        ptclCM->PredPosition[dim] = ptclCM->Position[dim];
        ptclCM->PredVelocity[dim] = ptclCM->Velocity[dim];

        ptclBin->Position[dim] = (ptclI->PredPosition[dim]*ptclI->Mass + ptclJ->PredPosition[dim]*ptclJ->Mass)/ptclCM->Mass;
        ptclBin->Velocity[dim] = (ptclI->PredVelocity[dim]*ptclI->Mass + ptclJ->PredVelocity[dim]*ptclJ->Mass)/ptclCM->Mass;
        ptclBin->PredPosition[dim] = ptclCM->Position[dim];
        ptclBin->PredVelocity[dim] = ptclCM->Velocity[dim];

    }

    ptclCM->CurrentTimeReg = ptclI->CurrentTimeReg;
    ptclCM->CurrentTimeIrr = ptclI->CurrentTimeIrr;
    ptclCM->TimeStepIrr = ptclI->TimeStepIrr;
    ptclCM->TimeStepReg = ptclI->TimeStepReg;
    ptclCM->TimeLevelIrr = ptclI->TimeLevelIrr;
    ptclCM->TimeLevelReg = ptclI->TimeLevelReg;
    ptclCM->PredTime = current_time;

    ptclCM->BinaryParticleI = ptclI;
    ptclCM->BinaryParticleJ = ptclJ;
    ptclCM->BinaryInfo = ptclBin;

    ptclCM->isCMptcl = true;
    ptclJ->isBinary = true;
    ptclBin->ptclCM = ptclCM;



    // add it to binary list
    BinaryList.push_back(ptclBin);


    // copy the neighbor list for c.m particle

    ptclCM->RadiusOfAC = ptclI->RadiusOfAC;

    for (Particle* ptcl: ptclI->ACList) {

        if (ptcl == ptclJ) {
            continue;
        }

        ptclCM->ACList.push_back(ptcl);
        ptclCM->NumberOfAC += 1;

    }


    // predict coordinates of neighbor particles

    for (Particle* ptcl: ptclI->ACList) {

        ptcl->predictParticleSecondOrder(current_time);

        // need to make routine for resolving KS components if a KS pair is in the neighbor

    }


    // calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately for the binary pair particle and the cm particle
    //std::cout << "CalculateKSAcceleration starts" << std::endl;


    CalculateKSAcceleration(ptclI,ptclJ,ptclCM,particle,current_time);

    //std::cout << "Calculating Time steps" << std::endl;

    ptclCM->calculateTimeStepIrr(ptclCM->a_tot, ptclCM->a_irr); // calculate irregular time step based on total force
    ptclCM->calculateTimeStepReg(ptclCM->a_reg, ptclCM->a_reg); // calculate regular time step based on total force


    fprintf(binout, "KSRegularlizationInitialization.cpp: result of CM particle value calculation\n");
    fprintf(binout, "from function NewKSInitialization\n");

    fprintf(binout, "Position - x:%f, y:%f, z:%f, \n", ptclCM->Position[0], ptclCM->Position[1], ptclCM->Position[2]);
    fprintf(binout, "Velocity - vx:%f, vy:%f, vz:%f, \n", ptclCM->Velocity[0], ptclCM->Velocity[1], ptclCM->Velocity[2]);

    fprintf(binout, "Total Acceleration - ax:%f, ay:%f, az:%f, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
    fprintf(binout, "Total Acceleration - axdot:%f, aydot:%f, azdot:%f, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
    fprintf(binout, "Total Acceleration - ax2dot:%f, ay2dot:%f, az2dot:%f, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
    fprintf(binout, "Total Acceleration - ax3dot:%f, ay3dot:%f, az3dot:%f, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
    fprintf(binout, "Time Steps - irregular:%f, regular:%f/n", ptclCM->TimeStepIrr, ptclCM->TimeStepReg);




    // calculate the initial values of relevant variables
    std::cout << "Initializing Binary Information" << std::endl;

    ptclBin->InitializeBinary(current_time);

    fprintf(binout, "KS coordinates - u1:%f, u2:%f, u3:%f, u4:%f/n", ptclBin->u[0], ptclBin->u[1], ptclBin->u[2], ptclBin->u[3]);
    fprintf(binout, "KS coordinates - udot1:%f, udot2:%f, udot3:%f, udot4:%f/n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
    fprintf(binout, "KS coordinates - udot1:%f, udot2:%f, udot3:%f, udot4:%f/n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
    fprintf(binout, "Other important KS variables - r:%f, h:%f, gamma: %f, tau:%f, step:%f/n", ptclBin->r, ptclBin->h, ptclBin->gamma, ptclBin->dTau, ptclBin->TimeStep);



    // delete the individual particles and add the new CM particle

    // delete the original components from the list

    //std::cout << "Removing Binary Particles from list" << std::endl;

    ptclIIndex = -1;
    ptclJIndex = -1;

    for (Particle* ptcl : particle) {

        ptclIIndex += 1;

        if (ptcl == ptclI) {
            break;
        }
    }

    particle.erase(particle.begin() + ptclIIndex);


    for (Particle* ptcl : particle) {

        ptclJIndex += 1;

        if (ptcl == ptclJ) {
            break;
        }
    }

    particle.erase(particle.begin() + ptclJIndex);

    // add the original particles

    //std::cout << "add the CM particle to list" << std::endl;
    particle.push_back(ptclCM);


    // we also need to change the neighbor list of Particles
    // assuming that all the neighbors are bidirectional
    // may need to update later if the radius for neighbor differs depending on the particle

    // first for particle I
    std::cout << "Changing the Neighbor List of Particles" << std::endl;
//  previous code that was unsuccessful
//    for (Particle* ptcl: ptclI->ACList) {

        // change the ptclI to ptclCM
//        for (int i = 1; i<ptcl->NumberOfAC; i++) {
//            if (ptcl->ACList[i] == ptclI) {
//                ptcl->ACList[i] = ptclCM;
//                auto it = std::find(ptclJ->ACList.begin(), ptclJ->ACList.end(), ptcl);
//                if (it != ptclJ->ACList.end()) {
//                    ptcl->ACList.erase(it);
//                }
//            }        

//        }
//    }

    // repeat the same loop for particle J with reduced ptcl list

//    for (Particle* ptcl: ptclJ->ACList) {
        // change the ptclI to ptclCM
//        for (int i = 0; i<ptcl->NumberOfAC; i++) {
//            if (ptcl->ACList[i] == ptclJ) {
//                ptcl->ACList[i] = ptclCM;
//            }        
//        }
//    }

    for (Particle* ptcl: particle){

	ptclIIndex = -1;
	ptclJIndex = -1;

	auto itI = std::find(ptcl->ACList.begin(), ptcl->ACList.end(), ptclI);

	if (itI != ptcl->ACList.end()) {
		ptclIIndex = std::distance(ptcl->ACList.begin(),itI);
		ptcl->ACList[ptclIIndex] = ptclCM;
	}

        auto itJ = std::find(ptcl->ACList.begin(), ptcl->ACList.end(), ptclJ);

        if (itJ != ptcl->ACList.end()) {
                ptclJIndex = std::distance(ptcl->ACList.begin(),itJ);
		if ( ptclIIndex == -1 ) {
			ptcl->ACList[ptclJIndex] = ptclCM;
		} else {
			ptcl->ACList.erase(itJ); 
		}
        }

    }




    // Add the binary to binary integration list
    //std::cout << "Add the ptclBin to Binary List" << std::endl;
    BinaryList.push_back(ptclBin);

}

