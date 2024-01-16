#include <vector>
#include <iostream>
#include <cmath>
#include "global.h"
#include "defs.h"




void FindNeighbor(std::vector<Particle*> &particle);
void calculateAcceleration(std::vector<Particle*> &particle);
void direct_sum(double *x, double *v, double r2, double vx,
		double mass, double a[3], double adot[3]);

/*
 *  Purporse: Initialize particles
 *
 *  Date    : 2024.01.02  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *
 */

void initializeParticle(std::vector<Particle*> &particle) {

	std::cout << "Initialization starts." << std::endl;
	FindNeighbor(particle);

	calculateAcceleration(particle);
	std::cout << "Timestep initializing..." << std::endl;
	for (Particle* elem:particle) {
		elem->initializeTimeStep();

		for (int dim=0; dim<Dim; dim++) {
			elem->PredPosition[dim] =  elem->Position[dim];
		}
	}

	std::cout << "Initialization finished." << std::endl;
}

/*
 *  Purporse: find neighbors
 *
 *  Date    : 2024.01.02  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *
 */
void FindNeighbor(std::vector<Particle*> &particle) {

	// No need to find neighbors if the total number of particles is less than 100
	if (NNB<=100) return;


	double r;
	Particle* sbj_ptcl;

	std::cout << "Finding neighbors ..." << std::endl;

	// search for neighbors for ptcl
	for (Particle *ptcl: particle) {
		for (int i=0; i<NNB; i++) {
			sbj_ptcl = particle[i];

			r = ptcl->getDistanceTo(sbj_ptcl);

			if (r < ptcl->RadiusOfAC && ptcl != sbj_ptcl) {
				ptcl->ACList.push_back(sbj_ptcl);
				ptcl->NumberOfAC++;
			}
		}
	}
	/***
		if (debug) {
		std::cout << particle[2]->NumberOfAC << std::endl;
		for (Particle* element : particle[2]->ACList)
		std::cout << element->getPID() << " ";
		}
		for (int i=0; i<NNB; i++) {
		std::cout << particle[i]->NumberOfAC << ' ';
		}
	 ***/
	std::cout << "Finding neighbors finished" << std::endl;
}


void calculateAcceleration(std::vector<Particle*> &particle) {

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

	std::cout << "Entering calculateAcceleration  ..." << std::endl;

	for (Particle *ptcl1:particle) {
		j=0;
		for (Particle *ptcl2:particle) {
			r2 = 0;
			vx = 0;
			v2 = 0;
			if (ptcl1 != ptcl2) {
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
					if ((ptcl1->NumberOfAC==0)||(ptcl2 != ptcl1->ACList[j])) {
						ptcl1->a_reg[dim][0] = m_r3*x[dim];
						ptcl1->a_reg[dim][1] = m_r3*(v[dim] - 3*x[dim]*vx/r2);
					}
					else {
						ptcl1->a_irr[dim][0] = m_r3*x[dim];
						ptcl1->a_irr[dim][1] = m_r3*(v[dim] - 3*x[dim]*vx/r2);
					}
				}


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
					if ((ptcl1->NumberOfAC==0)||(ptcl2 != ptcl1->ACList[j])) {
						ptcl1->a_reg[dim][2] = a2dot;
						ptcl1->a_reg[dim][3] = a3dot;
					}
					else {
						ptcl1->a_irr[dim][2] = a2dot;
						ptcl1->a_irr[dim][3] = a3dot;
						j++;
					}
				} // endfor dim
			} // endif ptcl1 == ptcl2
		} // endfor ptcl2
		std::cout << "\ntotal acceleartion\n" << std::flush;
		for (int dim=0; dim<Dim; dim++)	 {
			for (int order=0; order<HERMITE_ORDER; order++) {
				ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order];
				std::cout << ptcl1->a_tot[dim][order]*position_unit/time_unit/time_unit << " ";
			}
			std::cout << "\n" << std::endl;
		} // endfor dim
	} // endfor ptcl1
}



