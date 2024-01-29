#include <vector>
#include <iostream>
#include <cmath>
#include "global.h"
#include "defs.h"


void FindNeighbor(Particle* ptcl1, std::vector<Particle*> &particle);
void CalculateInitialAcceleration(Particle* ptcl1, std::vector<Particle*> &particle);
void direct_sum(double *x, double *v, double r2, double vx,
		double mass, double a[3], double adot[3]);

int InitializeTimeStep(Particle* particle, int size);
int InitializeTimeStep(std::vector<Particle*> &particle);

/*
 *  Purporse: Initialize particles
 *
 *  Date    : 2024.01.02  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *
 */

void InitializeParticle(std::vector<Particle*> &particle) {

	if (NNB == 0) {
		std::cout << "Initialization skips." << std::endl;
		return;
	}
	std::cout << "Initialization starts." << std::endl;

	// loop over particles to initialize their values
	for (Particle* ptcl:particle) {
		FindNeighbor(ptcl, particle);
		CalculateInitialAcceleration(ptcl, particle);
	}

	std::cout << "Timestep initializing..." << std::endl;
	InitializeTimeStep(particle);
	std::cout << "Timestep finished." << std::endl;
	for (Particle* elem:particle) {
		for (int dim=0; dim<Dim; dim++) {
			elem->PredPosition[dim] =  elem->Position[dim];
			elem->PredVelocity[dim] =  elem->Velocity[dim];
		}
	}
	std::cout << "Initialization finished." << std::endl;
}


/*
 *  Purporse: Initialize of new particles
 *
 *  Date    : 2024.01.16  by Seoyoung Kim
 *
 */
void InitializeParticle(Particle* newParticle, std::vector<Particle*> &particle) {

	std::cout << "Initialization of " << newNNB << " New Particles starts." << std::endl;

	// loop over particles to initialize acceleration
	for (int i=0; i<newNNB; i++) {
		FindNeighbor(&newParticle[i], particle);
		CalculateInitialAcceleration(&newParticle[i], particle);

		for (int dim=0; dim<Dim; dim++) {
			newParticle[i].PredPosition[dim] =  newParticle[i].Position[dim];
			newParticle[i].PredVelocity[dim] =  newParticle[i].Velocity[dim];
		}
	}

	std::cout << "Timestep initializing..." << std::endl;
	InitializeTimeStep(newParticle, newNNB);
	std::cout << "Timestep finished." << std::endl;

	std::cout << "Initialization of New Particles finished." << std::endl;
}



/*
 *  Purporse: find neighbors
 *
 *  Date    : 2024.01.02  by Seoyoung Kim
 *  Modified: 2024.01.16  by Seoyoung Kim
 *
 */


// recieve the time we want to calculate the neighbor information, 
// (which would be EnzoTimeMark in the case of New Particles)
// the particle we want to find the neighbors of and the full particle list
void FindNeighbor(Particle* ptcl1, std::vector<Particle*> &particle) {

	// No need to find neighbors if the total number of particles is less than 100
	if (NNB<=100) return;

	double dx;
	double r2;

	//ptcl1->predictParticleSecondOrder(newTime);
	// search for neighbors for ptcl
	for (Particle *ptcl2:particle) {

		if  (ptcl1 == ptcl2) {
			continue;
		}

		r2 = 0.0;

		//ptcl2->predictParticleSecondOrder(newTime);
		for (int dim=0; dim<Dim; dim++) {
			dx = ptcl2->Position[dim] - ptcl1->Position[dim];
			r2 += dx*dx;
		}

		if (sqrt(r2) < ptcl1->RadiusOfAC) {
			ptcl1->ACList.push_back(ptcl2);
			ptcl1->NumberOfAC++;
		}
	} // endfor 
	
	std::cout << "# of Neighbors = " << ptcl1->NumberOfAC << std::endl;
	return ;

}




/*
 *  Purporse: calculate 0th, 1st, 2nd, 3rd Derivative of Force
 *
 *  Date    : 2024.01.02  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *  Modified: 2024.01.16  by Seoyoung Kim
 *
 */


// recieve the time we want to calculate the acceleration information, 
// (which would be EnzoTimeMark in the case of New Particles)
// the particle we want to calculate and the full particle list
void CalculateInitialAcceleration(Particle* ptcl1, std::vector<Particle*> &particle) {

	int j=0;
	double x[Dim], v[Dim], a[Dim], adot[Dim];
	double vx_r2, m_r3, v2x2_r4,v2_r2__ax_r2__v2x2_r4, a2dot, a3dot;
	double A, B, v2;
	double r2 = 0;
	double vx = 0;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]    = 0.;
		v[dim]    = 0.;
		a[dim]    = 0.;
		adot[dim] = 0.;
	}

	//std::cout << "nbody+: Entering CalculateInitialAcceleration  ..." << std::endl;

	//ptcl1->predictParticleSecondOrder(newTime);
	for (Particle *ptcl2:particle) {
		r2 = 0;
		vx = 0;
		v2 = 0;

		if (ptcl1 == ptcl2) {
			continue;
		}

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		//ptcl2->predictParticleSecondOrder(newTime);
		for (int dim=0; dim<Dim; dim++) {
			x[dim] = ptcl2->Position[dim] - ptcl1->Position[dim];
			v[dim] = ptcl2->Velocity[dim] - ptcl1->Velocity[dim];
			r2    += x[dim]*x[dim];
			vx    += v[dim]*x[dim];
			v2    += v[dim]*v[dim];
		}

		r2  += EPS2;
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
			}
		} // endfor ptcl

		// Calculate 2nd and 3rd derivatives of acceleration
		if (restart) {
			;
			/*
				 for (int dim=0; dim<Dim; dim++) {
					 }*/
		}

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
			if ((ptcl1->NumberOfAC==0) || (ptcl2 != ptcl1->ACList[j])) {
				ptcl1->a_reg[dim][2] += a2dot;
				ptcl1->a_reg[dim][3] += a3dot;
			}
			else {
				ptcl1->a_irr[dim][2] += a2dot;
				ptcl1->a_irr[dim][3] += a3dot;
				j++;
			}
		} // endfor dim


		// if the ptcl1 is new particle and ptcl2 is old particle
		// then we need to change the values of ptcl2 as well. 

		// add the regular and irregualar forces and change the neighbor list of the other particle


		// needed to update
		if ( (ptcl1->ParticleType == NewParticle+NormalStar+SingleParticle)
			 	&& (ptcl2->ParticleType == SingleParticle+NormalStar)) {
			if (sqrt(r2)<ptcl2->RadiusOfAC) {
				for (int dim=0; dim<Dim; dim++) {
					ptcl2->a_irr[dim][0] += -m_r3*x[dim];
					ptcl2->a_irr[dim][1] += -m_r3*(v[dim] - 3*x[dim]*vx/r2);
					ptcl2->a_irr[dim][2] += -a2dot;
					ptcl2->a_irr[dim][3] += -a3dot;
				}
				ptcl2->ACList.push_back(ptcl2);
				ptcl2->NumberOfAC++;
			} else {
				for (int dim=0; dim<Dim; dim++) {
					ptcl2->a_reg[dim][0] += -m_r3*x[dim];
					ptcl2->a_reg[dim][1] += -m_r3*(v[dim] - 3*x[dim]*vx/r2);
					ptcl2->a_reg[dim][2] += -a2dot;
					ptcl2->a_reg[dim][3] += -a3dot;
				}
			}
		}

	} // endfor ptcl2


	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=0; order<HERMITE_ORDER; order++) {
			ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order]; 
		}
	}

	/*
	std::cout << "NBODY+: total acceleartion\n" << std::flush;
	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=0; order<HERMITE_ORDER; order++) {
			std::cout << ptcl1->a_reg[dim][order]*position_unit/time_unit/time_unit << " ";
		}
	} // endfor dim
	std::cout << std::endl;
	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=0; order<HERMITE_ORDER; order++) {
			std::cout << ptcl1->a_irr[dim][order]*position_unit/time_unit/time_unit << " ";
		}
	} // endfor dim
	std::cout << '\n' << std::endl;
	*/
}


