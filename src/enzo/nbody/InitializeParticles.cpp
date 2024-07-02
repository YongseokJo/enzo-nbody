#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "global.h"
#include "defs.h"



void FindNeighbor(Particle* ptcl1, std::vector<Particle*> &particle);
void FindNewNeighbor(Particle* ptcl1, std::vector<Particle*> &particle);
void CalculateAcceleration01(Particle* ptcl1, std::vector<Particle*> &particle);
void CalculateAcceleration23(Particle* ptcl1, std::vector<Particle*> &particle);
void CalculateAcceleration23_new(Particle* ptcl1, std::vector<Particle*> &particle);
void direct_sum(double *x, double *v, double r2, double vx,
		double mass, double a[3], double adot[3]);
int InitializeTimeStep(std::vector<Particle*> &particle, int size);
int InitializeTimeStep(std::vector<Particle*> &particle);
void UpdateNextRegTime(std::vector<Particle*> &particle);

/*
 *  Purporse: Initialize particles
 *
 *  Date    : 2024.01.02  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *  Modified: 2024.06.28  by Yongseok Jo
 *
 */

void InitializeParticle(std::vector<Particle*> &particle) {

	if (NNB == 0) {
		fprintf(nbpout, "Initialization skips.\n");
		return;
	}
	fprintf(nbpout, "Initialization starts.\n");

	// loop over particles to initialize their values
	int j = 0;

	for (Particle* ptcl:particle) {
		//std::cout <<  j << ": ";
		//std::cout << std::flush;
		FindNeighbor(ptcl, particle);
		CalculateAcceleration01(ptcl, particle);
	}
	for (Particle* ptcl:particle) {
		CalculateAcceleration23(ptcl, particle);
		std::cout << std::flush;
	}

	fprintf(nbpout, "Timestep initializing.\n");
	InitializeTimeStep(particle);
	fprintf(nbpout, "Timestep finished.\n");

	UpdateNextRegTime(particle);

	/*
	RegularList.clear();
	for (Particle* ptcl:particle) {
		for (int dim=0; dim<Dim; dim++) {
			ptcl->PredPosition[dim] =  ptcl->Position[dim];
			ptcl->PredVelocity[dim] =  ptcl->Velocity[dim];
		}
		RegularList.push_back(ptcl);
	}
	*/
	fprintf(nbpout, "Initialization finished.\n");
}


/*
 *  Purporse: Initialize of new particles
 *
 *  Date    : 2024.01.16  by Seoyoung Kim
 *  Modified: 2024.06.28  by Yongseok Jo
 *
 */
void InitializeNewParticle(std::vector<Particle*> &particle, int offset) {

	std::cout << "Initialization of " << newNNB << " New Particles starts." << std::endl;

	// loop over particles to initialize acceleration
	for (int i=0; i<newNNB; i++) {
		FindNewNeighbor(particle[offset+i], particle);
		CalculateAcceleration01(particle[offset+i], particle);
	}

	for (int i=0; i<newNNB; i++) {
		CalculateAcceleration23(particle[offset+i], particle);

		for (int dim=0; dim<Dim; dim++) {
			particle[i]->PredPosition[dim] =  particle[i]->Position[dim];
			particle[i]->PredVelocity[dim] =  particle[i]->Velocity[dim];
		}
	}

	std::cout << "Timestep initializing..." << std::endl;
	InitializeTimeStep(particle, offset);
	std::cout << "Timestep finished." << std::endl;

	std::cout << "Initialization of New Particles finished." << std::endl;
}



/*
 *  Purporse: find neighbors
 *
 *  Date    : 2024.01.02  by Seoyoung Kim
 *  Modified: 2024.01.16  by Seoyoung Kim
 *  Modified: 2024.06.28  by Yongseok Jo
 *
 */
// recieve the time we want to calculate the neighbor information, 
// (which would be EnzoTimeMark in the case of New Particles)
// the particle we want to find the neighbors of and the full particle list
void FindNeighbor(Particle* ptcl1, std::vector<Particle*> &particle) {

	// No need to find neighbors if the total number of particles is less than 100
	//if (NNB<=100) return;

	double dx;
	double r2;
	double r_max=0;
	double r_nb[FixNumNeighbor];
	int index_max, i=0; // nb_index[FixNumNeighbor], 
	std::vector<int> nb_index(FixNumNeighbor); 

	// search for neighbors for ptcl
	for (Particle *ptcl2:particle) {

		if  (ptcl1 == ptcl2) {
			i++;
			continue;
		}

		r2 = 0.0;
		for (int dim=0; dim<Dim; dim++) {
			dx = ptcl2->Position[dim] - ptcl1->Position[dim];
			r2 += dx*dx;
		}

		if (ptcl1->NumberOfAC < FixNumNeighbor) {
			if (r2 > r_max) {
				r_max     = r2;
				index_max = ptcl1->NumberOfAC;
			}
			nb_index[ptcl1->NumberOfAC] = i;
			r_nb[ptcl1->NumberOfAC]   = r2;
			ptcl1->NumberOfAC++;
		}
		else {
			if ( r2 < r_max) {
				r_nb[index_max]     = r2;
				nb_index[index_max] = i;
				// update new r_max
				r_max = r2;
				for (int k=0; k<FixNumNeighbor; k++) {
					if (r_nb[k] > r_max) {
						r_max     = r_nb[k];
						index_max = k;
					}
				}
			}
		}
		i++;
	} // endfor


	std::sort(nb_index.begin(), nb_index.end());
	for (int j:nb_index) {
		ptcl1->ACList.push_back(particle[j]);
		//std::cout << j << ", ";
	}


	fprintf(stdout, "Neighbors of %d : \n", ptcl1->PID);
	for (Particle* ptcl:ptcl1->ACList) {
		std::cout << ptcl->PID << ", ";
	}
	std::cout << std::endl;

	//std::sort(ptcl1->ACList.begin(),ptcl1->ACList.end());
	//std::cout << "# of Neighbors = " << ptcl1->NumberOfAC << std::endl;
	return ;
}


void FindNewNeighbor(Particle* newPtcl, std::vector<Particle*> &particle) {

	// No need to find neighbors if the total number of particles is less than 100
	//if (NNB<=100) return;
	
	double r0, r1;
	double r_max;
	int index, max_index;

	FindNeighbor(newPtcl, particle);

	for (Particle* ptcl: newPtcl->ACList) {
		index = 0;
		r_max = 0;
		max_index = -1;

		r0 = dist2(ptcl->Position, newPtcl->Position);

		for (Particle* neighbor:ptcl->ACList) {
			r1 = dist2(ptcl->Position, neighbor->Position);

			if (r1 > r_max) {
				r_max = r1;
				max_index = index;
			}

			index++;
		}

		if (r0 < r_max) {
			if (max_index == -1 || r_max == 0) {
				fprintf(stderr, "Max Index = %d, r_max = %e\n", max_index, r_max);
				throw std::runtime_error("Fatal Error in FindNewNeighbor.cpp\n");	
			}
			ptcl->ACList.erase(ptcl->ACList.begin() + max_index);
			ptcl->ACList.push_back(newPtcl);
		}
	}

	return ;
}


void ReInitializeKSParticle(Particle* KSParticle, std::vector<Particle*> &particle) {

	std::cout << "Re-Initialization of KS Pair Particle starts.\n" << std::endl;

	std::cout << "Finding Neighbors... \n" << std::endl;	
	KSParticle->ACList.clear();
	FindNeighbor(KSParticle, particle);

	std::cout << "Calculating Acceleration... \n" << std::endl;	

	CalculateAcceleration01(KSParticle, particle);
	CalculateAcceleration23(KSParticle, particle);

	std::cout << "copying the position and velocities to predictions... \n" << std::endl;	

	for (int dim=0; dim<Dim; dim++) {
		KSParticle->PredPosition[dim] =  KSParticle->Position[dim];
		KSParticle->PredVelocity[dim] =  KSParticle->Velocity[dim];
		KSParticle->NewPosition[dim]  =  KSParticle->Position[dim];
		KSParticle->NewVelocity[dim]  =  KSParticle->Velocity[dim];
	}

	std::cout << "Timestep calculation...\n" << std::endl;


	// advance here because cm particle didn't advance and was removed.
	// ComputationList.push_back(KSParticle);

	std::cout << "Timestep finished.\n" << std::endl;
	std::cout << "Initialization of KS Pair Particle finished.\n" << std::endl;
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
void CalculateAcceleration01(Particle* ptcl1, std::vector<Particle*> &particle) {

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

		m_r3 = ptcl2->Mass/r2/sqrt(r2); 

		if ((ptcl1->NumberOfAC==0) || (ptcl2 != ptcl1->ACList[j])) {
			for (int dim=0; dim<Dim; dim++) {
				// Calculate 0th and 1st derivatives of acceleration
				ptcl1->a_reg[dim][0] += m_r3*x[dim];
				ptcl1->a_reg[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++) {
				ptcl1->a_irr[dim][0] += m_r3*x[dim];
				ptcl1->a_irr[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
			j++;
		} // endfor dim
	} // endfor ptcl2
	
	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=0; order<2; order++) {
			ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order] + ptcl1->BackgroundAcceleration[dim]; 
		}
	}
	return;
}


	/*
		// Calculate 2nd and 3rd derivatives of acceleration
		if (restart) {
			;
		}
	*/

void CalculateAcceleration23(Particle* ptcl1, std::vector<Particle*> &particle) {

	int j=0;
	double x[Dim], v[Dim], a21[Dim], a21dot[Dim], a1[Dim], a2[Dim], a1dot[Dim], a2dot[Dim];
	double a, b, c;
	double rdf_r2, vdf_r2, rdfdot_r2, v2, r2, r3, vr, m_r3;
	double adot2, adot3;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]      = 0.;
		v[dim]      = 0.;
		a21[dim]    = 0.;
		a21dot[dim] = 0.;
		a1[dim]     = ptcl1->a_tot[dim][0];
		a1dot[dim]  = ptcl1->a_tot[dim][1];
	}

	for (Particle *ptcl2:particle) {
		r2 = 0;
		r3 = 0;
		v2 = 0;
		vr = 0;
		rdf_r2 = 0;
		vdf_r2 = 0;
		rdfdot_r2 = 0;

		if (ptcl1 == ptcl2) {
			continue;
		}

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		for (int dim=0; dim<Dim; dim++) {
			a2[dim]    = ptcl2->a_tot[dim][0];
			a2dot[dim] = ptcl2->a_tot[dim][1];
			x[dim]     = ptcl2->Position[dim] - ptcl1->Position[dim];
			v[dim]     = ptcl2->Velocity[dim] - ptcl1->Velocity[dim];
			r2        += x[dim]*x[dim];
			vr        += v[dim]*x[dim];
			v2        += v[dim]*v[dim];
		}

		r3   = r2*sqrt(r2);
		m_r3 = ptcl2->Mass/r3; 

		for (int dim=0; dim<Dim; dim++) {
			a21[dim]    = m_r3*x[dim];
			a21dot[dim] = m_r3*(v[dim] - 3*x[dim]*vr/r2);
			rdf_r2     += x[dim]*(a1[dim]-a2[dim])/r2;
			vdf_r2     += v[dim]*(a1[dim]-a2[dim])/r2;
			rdfdot_r2  += x[dim]*(a1dot[dim]-a2dot[dim])/r2;
		}

		a = vr/r2;
		b = v2/r2 + rdf_r2 + a*a;
		c = 3*vdf_r2 + rdfdot_r2 + a*(3*b-4*a*a);


		if ((ptcl1->NumberOfAC==0) || (ptcl2 != ptcl1->ACList[j])) {
			for (int dim=0; dim<Dim; dim++) {
				adot2 = -ptcl2->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
				adot3 = -ptcl2->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
				ptcl1->a_reg[dim][2] += adot2;
				ptcl1->a_reg[dim][3] += adot3;
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++) {
				adot2 = -ptcl2->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
				adot3 = -ptcl2->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
				ptcl1->a_irr[dim][2] += adot2;
				ptcl1->a_irr[dim][3] += adot3;
			}
			j++;
		} // endfor if
	} //endfor ptcl2
	
	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=2; order<4; order++) {
			ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order]; 
		}
	}

		// if the ptcl1 is new particle and ptcl2 is old particle
		// then we need to change the values of ptcl2 as well. 

		// add the regular and irregualar forces and change the neighbor list of the other particle


		/*
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
					ptcl2->a_reg[dim][0] += - m_r3*x[dim];
					ptcl2->a_reg[dim][1] += - m_r3*(v[dim] - 3*x[dim]*vx/r2);
					ptcl2->a_reg[dim][2] += - a2dot;
					ptcl2->a_reg[dim][3] += - a3dot;
				}
			}
		}
		*/



	/*
	std::cout << "NBODY+: total acceleartion = " << std::flush;
	for (int order=0; order<HERMITE_ORDER; order++) {
		std::cout << "a" << order << "(";
		for (int dim=0; dim<Dim; dim++)	 {
			std::cout << std::scientific << std::setprecision(2) << ptcl1->a_tot[dim][order] << ", ";
		}
		std::cout << "), ";
	} // endfor dim
	std::cout << std::endl;

	//std::cout << "AC=" << ptcl1->ACList[0] << std::endl;

	for (int order=0; order<HERMITE_ORDER; order++) {
		std::cout << "a_reg" << order << "(";
		for (int dim=0; dim<Dim; dim++)	 {
			std::cout << std::scientific << std::setprecision(2) << ptcl1->a_reg[dim][order] << ", "; //position_unit/time_unit/time_unit  
		}
		std::cout << "), ";
	} // endfor dim
	std::cout << std::endl;

	for (int order=0; order<HERMITE_ORDER; order++) {
		std::cout << "a_irr" << order << "(";
		for (int dim=0; dim<Dim; dim++)	{
			std::cout << std::scientific << std::setprecision(2) << ptcl1->a_irr[dim][order] << ", "; //position_unit/time_unit/time_unit  
		}
		std::cout << "), ";
	} // endfor dim
	std::cout << std::endl;
	*/
}


void CalculateAcceleration23_new(Particle* ptcl1, std::vector<Particle*> &particle) {

	int j=0;
	double a2dot[Dim], a3dot[Dim], r2, a1dotk;
	double a[12];
	double a13, a14, a15, a16, a17, a18, a19, a20, a21, a22;

	for (Particle *ptcl2:particle) {
		r2 = 0;

		if (ptcl1 == ptcl2) {
			continue;
		}

		for (int dim=0; dim<Dim; dim++) {
			a[dim]   = ptcl2->Position[dim] - ptcl1->Position[dim];
			a[3+dim] = ptcl2->Velocity[dim] - ptcl1->Velocity[dim];
			a[6+dim] = ptcl2->a_tot[dim][0] - ptcl1->a_tot[dim][0];
			a[9+dim] = ptcl2->a_tot[dim][1] - ptcl1->a_tot[dim][1];
			r2 += a[dim]*a[dim];
		}
		/*
		if (r2 > 0.2 || (ptcl2 != ptcl1->ACList[j]))
			continue;
			*/

		a13 = 1/r2;
		a14 = ptcl2->Mass*a13*sqrt(a13);
		a15 = (a[0]*a[3]+a[1]*a[4]+a[2]*a[5])*a13;
		a16 = a15*a15;
		a17 = 3*a15;
		a18 = 6*a15;
		a19 = 9*a15;
		a20 = (a[3]*a[3]+a[4]*a[4]+a[5]*a[5]+a[0]*a[6]+a[1]*a[7]+a[2]*a[8])*a13 + a16;
		a21 = 9.0*a20;
		a20 = 3.0*a20;
		a22 = (9.0*(a[3]*a[6] + a[4]*a[7] + a[5]*a[8]) +\
				3.0*(a[0]*a[9] + a[1]*a[10] + a[2]*a[11]))*a13 +\
					a17*(a20 - 4.0*a16);


		for (int dim=0; dim<Dim; dim++) {
			a1dotk     = a[dim+3] - a17*a[dim];
			a2dot[dim] = (a[dim+6] - a18*a1dotk - a20*a[dim])*a14;
			a3dot[dim] = (a[dim+9] - a21*a1dotk - a22*a[dim])*a14\
									 - a19*a2dot[dim];
		
		}

		if ((ptcl1->NumberOfAC==0) || (ptcl2 != ptcl1->ACList[j])) {
			for (int dim=0; dim<Dim; dim++) {
				ptcl1->a_reg[dim][2] += a2dot[dim];
				ptcl1->a_reg[dim][3] += a3dot[dim];
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++) {
				ptcl1->a_irr[dim][2] += a2dot[dim];
				ptcl1->a_irr[dim][3] += a3dot[dim];
			}
			j++;
		} // endfor if
	} //endfor ptcl2
	
	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=2; order<4; order++) {
			ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order]; 
		}
	}
}



