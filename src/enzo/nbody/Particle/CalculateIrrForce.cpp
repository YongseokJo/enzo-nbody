#include <vector>
#include <iostream>
#include <cmath>
#include "../global.h"


void direct_sum(double *x, double *v, double r2, double vx,
	 	        double mass, double mdot, double (&a)[2][3], double (&adot)[2][3], int p) {
	double _r3;

	if (r2 < EPS2)
		r2 = EPS2;  // add softening length


	_r3 = 1/r2/sqrt(r2);

	for (int dim=0; dim<Dim; dim++){
		a[p][dim]    += mass*_r3*x[dim];
		adot[p][dim] += mass*_r3*(v[dim] - 3*x[dim]*vx/r2)+mdot*x[dim]*_r3;
	}
}

/*
 *  Purporse: Calculate the and update the irregular acceleration of particles
 *            using a_0, d(a_0)/dt, a_p, d(a_p)/dt; Nbody6++ manual Eq. 8.9 and 8.10
 * 			  Note that a_p, d(a_p)/dt are calculated from predicted positions and velocities
 *
 *  Date    : 2024.01.09  by Seoyoung Kim
 *  Modified: 2024.01.11  by Seoyoung Kim
 *
 */


void Particle::calculateIrrForce() {

	if (this->NumberOfAC == 0) {
		//CurrentTimeIrr += TimeStepIrr;
		std::cout << "Error: No neighbor in Irregular force!!" << std::endl;
		return;
	}

	double dt, mdot, epsilon=1e-6;
	double tirr[2]; // 0 for current and 1 for advanced times

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a0_irr[2][Dim], a0dot_irr[2][Dim]; // 0 for current and 1 for predicted accelerations


	dt      = TimeStepIrr*EnzoTimeStep; // interval of time step
	tirr[0] = CurrentTimeIrr; // current time of particle
	tirr[1] = CurrentTimeIrr + TimeStepIrr; // the time to be advanced to

	// initialize irregular force terms for ith particle just in case
	for (int p=0; p<2; p++) {
		for (int dim=0; dim<Dim; dim++){
			a0_irr[p][dim]    = 0.0;
			a0dot_irr[p][dim] = 0.0;
		}
	}

	// scan over the neighbor lists to find the irregular acceleration components
	//std::cout <<  "Looping single particles to calculate irregular acceleration...\n" << std::flush;

	for (int p=0; p<2; p++) {
		for (Particle* ptcl: ACList) {

			// reset temporary variables at the start of a new calculation
			r2 = 0.0;
			vx = 0.0;

			ptcl->predictParticleSecondOrder(tirr[p]);
			this->predictParticleSecondOrder(tirr[p]);
			for (int dim=0; dim<Dim; dim++) {
				// calculate position and velocity differences for current time
				x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
				v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];

				// calculate the square of radius and inner product of r and v for each case
				r2 += x[dim]*x[dim];
				vx += v[dim]*x[dim];
			}


			mdot = ptcl->evolveStarMass(tirr[0],
					tirr[0]+TimeStepIrr*1.01)/TimeStepIrr*1e-2; // derivative can be improved
																													 //
			// add the contribution of jth particle to acceleration of current and predicted times
			direct_sum(x ,v, r2, vx, ptcl->Mass, mdot, a0_irr, a0dot_irr, p);
			if (p == 0) {
				for (int dim=0; dim<Dim; dim++) {
					a_tot[dim][0] = a_reg[dim][0] + a0_irr[0][dim];
					a_tot[dim][1] = a_reg[dim][1] + a0dot_irr[0][dim];
				}
			}
		} // endfor ptcl
	} //endfor i, 0 for current time and 1 for provisional


	// correct the force using the 4th order hermite method
	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4;

	// set the values of calculation variables
	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;

	// calculate the correction terms and ...
	// calculate the final acceleration using correction terms

	for (int dim=0; dim<Dim; dim++) {
		da_dt2  = (   a0_irr[0][dim] - a0_irr[1][dim]   ) / dt2;
		adot_dt = (a0dot_irr[0][dim] + a0dot_irr[1][dim]) / dt;

		a2 =  -6*da_dt2 - 2*adot_dt - 2*a0dot_irr[0][dim]/dt;
		a3 = (12*da_dt2 + 6*adot_dt)/dt;

		a_irr[dim][0] = a0_irr[0][dim];
		a_irr[dim][1] = a0dot_irr[0][dim];
		a_irr[dim][2] = a2;
		a_irr[dim][3] = a3;

		a_tot[dim][0] = a_reg[dim][0] + a_irr[dim][0];
		a_tot[dim][1] = a_reg[dim][1] + a_irr[dim][1];
		a_tot[dim][2] = a_reg[dim][2] + a_irr[dim][2];
		a_tot[dim][3] = a_reg[dim][3] + a_irr[dim][3];
	}

	/*
	std::cout << "\ntotal acceleartion\n" << std::flush;
	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=0; order<HERMITE_ORDER; order++) {
			a_tot[dim][order] = a_reg[dim][order] + a_irr[dim][order];
			std::cout << a_tot[dim][order]*position_unit/time_unit/time_unit << " ";
		}
		std::cout << "\n" << std::endl;
		} // endfor dim
		*/
	
	std::cout << "\nIrregular Calculation\n" << std::flush;
	std::cout <<  "3. a_irr= "<< a_irr[0][0]<< ',' << a_irr[1][0]\
		<< ',' << a_irr[2][0] << std::endl;
	std::cout <<  "4. a_irr= "<< a_irr[0][0]<< ',' << a_irr[1][0]\
		<< ',' << a_irr[2][0] << std::endl;
	std::cout <<  "5. a_irr= "<< a_irr[0][0]<< ',' << a_irr[1][0]\
		<< ',' << a_irr[2][0] << std::endl;

	// update the current irregular time and irregular time steps
	//this->updateParticle((CurrentTimeIrr+TimeStepIrr)*EnzoTimeStep, a_irr);
	this->updateParticle(CurrentTimeIrr+TimeStepIrr, a_tot);
	CurrentTimeIrr += TimeStepIrr;
	this->calculateTimeStepIrr(a_tot, a_irr); // calculate irregular time step based on total force
}


