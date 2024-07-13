#include <vector>
#include <iostream>
#include <cmath>
#include "../global.h"
#include <cassert>


void direct_sum(double *x, double *v, double r2, double vx,
	 	        double mass, double mdot, double a[3], double adot[3]) {
	double _r3;

	if (r2 < EPS2)
		r2 = EPS2;  // add softening length

	_r3 = 1/r2/sqrt(r2);

	for (int dim=0; dim<Dim; dim++){
		a[dim]    += mass*_r3*x[dim];
		adot[dim] += mass*_r3*(v[dim] - 3*x[dim]*vx/r2); //+mdot*x[dim]*_r3;
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
	double new_time; // 0 for current and 1 for advanced times

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a_tmp[Dim], adot_tmp[Dim]; // 0 for current and 1 for predicted accelerations

	double m_r3;
	new_time = CurrentTimeIrr + TimeStepIrr; // the time to be advanced to
	dt       = TimeStepIrr*EnzoTimeStep; // interval of time step


	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<Dim; dim++){
		a_tmp[dim]    = 0.0;
		adot_tmp[dim] = 0.0;
	}



	/*******************************************************
	 * Irregular Acceleartion Calculation
	 ********************************************************/
	this->predictParticleSecondOrderIrr(new_time);

	

	for (Particle* ptcl: ACList) {
		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predictParticleSecondOrderIrr(new_time);
		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
			v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}
		if (r2 == 0)  {
			fprintf(stderr, "r2 is zero (%d and %d)\n",PID, ptcl->PID);
			//fflush(stderr);
			continue;
		}

		/*
		if (StarParticleFeedback != 0)
			mdot = ptcl->evolveStarMass(CurrentTimeIrr,
					CurrentTimeIrr+TimeStepIrr*1.01)/TimeStepIrr*1e-2; // derivative can be improved
																													   // add the contribution of jth particle to acceleration of current and predicted times
		*/

		if (r2 < EPS2)
			r2 = EPS2;  // add softening length


		m_r3 = ptcl->PredMass/r2/sqrt(r2);


		for (int dim=0; dim<Dim; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}

		//direct_sum(x ,v, r2, vx, ptcl->Mass, mdot, a_tmp, adot_tmp);
	} // endfor ptcl


	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;
	double dt_ex = (new_time - this->CurrentTimeReg)*EnzoTimeStep;

	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;

	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	for (int dim=0; dim<Dim; dim++) {


		// do the higher order correcteion
		da_dt2  = (a_irr[dim][0] - a_tmp[dim]) / dt2; //  - a_reg[dim][1]*dt_ex
		adot_dt = (a_irr[dim][1] + adot_tmp[dim]) / dt;
		a2 =  -6*da_dt2  - 2*adot_dt - 2*a_irr[dim][1]/dt;
		a3 =  (12*da_dt2 + 6*adot_dt)/dt;

		//fprintf(stderr, "da_dt2=%.2e, adot_dt=%.2e, a2=%.2e, a3=%.2e\n", da_dt2, adot_dt, a2, a3);

		// 4th order correction
		// save the values in the temporary variables
		NewPosition[dim] = PredPosition[dim] + a2*dt4/24 + a3*dt5/120;
		NewVelocity[dim] = PredVelocity[dim] + a2*dt3/6  + a3*dt4/24;

		//NewPosition[dim] = PredPosition[dim];// + a2*dt4/24 + a3*dt5/120;
		//NewVelocity[dim] = PredVelocity[dim];// + a2*dt3/6  + a3*dt4/24;

		// note that these higher order terms and lowers have different neighbors
		a_irr[dim][0] = a_tmp[dim];
		a_irr[dim][1] = adot_tmp[dim];
		a_irr[dim][2] = a2;
		a_irr[dim][3] = a3;
	}


	for (int dim=0; dim<Dim; dim++) {
		a_tot[dim][0] = a_reg[dim][0] + a_irr[dim][0] + a_reg[dim][1]*dt_ex; //+BackgroundAcceleration[dim]; // affect the next
		a_tot[dim][1] = a_reg[dim][1] + a_irr[dim][1];
		a_tot[dim][2] = a_reg[dim][2] + a_irr[dim][2];
		a_tot[dim][3] = a_reg[dim][3] + a_irr[dim][3];
	}

}


