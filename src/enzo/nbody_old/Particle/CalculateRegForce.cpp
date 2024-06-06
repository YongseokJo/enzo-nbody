#include <vector>
#include <iostream>
#include <cmath>
#include "../global.h"
#include "../defs.h"


/*
 *  Purporse: Particle Class
 *
 *  Date    : 2024.01.02  by Yongseok Jo
 *  Modified: 2024.01.30  by Yongseok Jo - Mass correction added
 *
 */

void direct_sum(double *x, double *v, double r2, double vx,
		double mass, double mdot, double (&a)[3], double (&adot)[3]) {

	double _r3;

	r2 += EPS2;  // add softening length
	_r3 = 1/r2/sqrt(r2);

	for (int dim=0; dim<Dim; dim++) {
		a[dim]    += mass*_r3*x[dim];
		adot[dim] += mass*_r3*(v[dim] - 3*x[dim]*vx/r2)+mdot*_r3*x[dim];
	}
}



/*
 *  Purporse: Calculate regular forces up to 2nd order
 *
 *  Date    : 2024.01.09  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *  Modified: 2024.01.25  by Yongseok Jo
 *  Modified: 2024.01.30  by Yongseok Jo - Mass correction added
 *
 */
void Particle::calculateRegAccelerationSecondOrder(std::vector<Particle*> &particle) {

	int maxNeighborNum = 250;
	double rs2 = 0.1;  // scale radius...? need to check

	double dt, r, r2, vx;
	double a0_reg[Dim], a0dot_reg[Dim], x[Dim], v[Dim];
	double a0_irr[Dim], a0dot_irr[Dim];
	double treg, mdot;

	for (int i=0; i<Dim; i++) {
		x[i]            = 0.;
		v[i]            = 0.;
		a0_reg[i]    = 0.;
		a0_irr[i]    = 0.;
		a0dot_reg[i] = 0.;
		a0dot_irr[i] = 0.;
	}

	dt      = TimeStepReg*EnzoTimeStep;
	ACList.clear();
	NumberOfAC = 0;

	// scan over all single particles to find the regular force components
	for (Particle* ptcl: particle) {
		if (ptcl == this)
			continue;

		r2 = 0;
		vx = 0;

		ptcl->predictParticleSecondOrder(CurrentTimeReg); // this takes nbody unit
		this->predictParticleSecondOrder(CurrentTimeReg);
		for (int dim=0; dim<Dim; dim++) {
			// When particles are not at the current time, extrapolate up to 2nd order
			x[dim] = ptcl->PredPosition[dim] - Position[dim];
			v[dim] = ptcl->PredVelocity[dim] - Velocity[dim];

			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		r    = sqrt(r2);
		mdot = ptcl->evolveStarMass(CurrentTimeReg,
			 	CurrentTimeReg+TimeStepReg*1e-2)/TimeStepReg*1e-2; // derivative can be improved
																																									//
		//std::cout << "Mass=" << ptcl->Mass << ", mdot=" << mdot << std::endl;

		if ((NNB > 100) && ((r < this->RadiusOfAC) || (r < ptcl->RadiusOfAC))) {
			NumberOfAC++;
			this->ACList.push_back(ptcl);
			// irregular acceleration
			direct_sum(x ,v, r2, vx, ptcl->Mass, mdot, a0_irr, a0dot_irr);
		}
		else {
			// regular acceleration
			//fprintf(stdout, "r2=%lf, vx=%lf\n", r2, vx);
			direct_sum(x ,v, r2, vx, ptcl->Mass, mdot, a0_reg, a0dot_reg);
			//fprintf(stdout, "a0_reg=%lf\n", a0_reg[0][0]);
		}
	} // endfor ptcl

	for (int dim=0; dim<Dim; dim++) {
		a_reg[dim][0] = a0_reg[dim];
		a_reg[dim][1] = a0dot_reg[dim];
		a_irr[dim][0] = a0_irr[dim];
		a_irr[dim][1] = a0dot_irr[dim];
		a_tot[dim][0] = a_reg[dim][0] + a_irr[dim][0];
		a_tot[dim][1] = a_reg[dim][1] + a_irr[dim][1];
	}
}

/*
 *  Purporse: Calculate regular forces up to 4th order
 *
 *  Modified: 2024.01.25  by Yongseok Jo
 *  Modified: 2024.01.30  by Yongseok Jo - Mass correction added
 *
 */
void Particle::calculateRegAccelerationFourthOrder(std::vector<Particle*> &particle) {

	double dt, r, r2, vx, mdot, epsilon = 1e-6;
	double a0_reg[Dim], a0dot_reg[Dim], x[Dim], v[Dim];
	double a0_irr[Dim], a0dot_irr[Dim];
	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4;

	for (int i=0; i<Dim; i++) {
		x[i]         = 0.;
		v[i]         = 0.;
		a0_reg[i]    = 0.;
		a0_irr[i]    = 0.;
		a0dot_reg[i] = 0.;
		a0dot_irr[i] = 0.;
	}

	dt = TimeStepReg*EnzoTimeStep;

	// scan over all single particles to find the regular force components
	//std::cout <<  "Looping single particles for force...\n" << std::flush;
	for (Particle* ptcl: particle) {
		if (ptcl == this)
			continue;

		r2 = 0;
		vx = 0;

		ptcl->predictParticleSecondOrder(CurrentTimeReg+TimeStepReg); // this takes nbody unit
		this->predictParticleSecondOrder(CurrentTimeReg + TimeStepReg);
		for (int dim=0; dim<Dim; dim++) {
			// When particles are not at the current time, extrapolate up to 2nd order
			x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
			v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];

			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}
		mdot = ptcl->evolveStarMass(CurrentTimeReg+TimeStepReg,
			 	CurrentTimeReg+TimeStepReg*1.01)/TimeStepReg*1e-2; // derivative can be improved

		// regular acceleration
		direct_sum(x ,v, r2, vx, ptcl->Mass, mdot, a0_reg, a0dot_reg);
	} // endfor ptcl

	// scan over neighors for irregular forces
	for (Particle* ptcl: ACList) {

		r2 = 0;
		vx = 0;

		ptcl->predictParticleSecondOrder(CurrentTimeReg); // this takes nbody unit
		for (int dim=0; dim<Dim; dim++) {
			// When particles are not at the current time, extrapolate up to 2nd order
			x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
			v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];

			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}
		mdot = ptcl->evolveStarMass(CurrentTimeReg+TimeStepReg,
			 	CurrentTimeReg+TimeStepReg*1.01)/TimeStepReg*1e-2; // derivative can be improved

		direct_sum(x ,v, r2, vx, ptcl->Mass, mdot, a0_irr, a0dot_irr);
		// subtract the irr acc that had been added to reg in the first loop
		direct_sum(x ,v, r2, vx, -ptcl->Mass, mdot, a0_reg, a0dot_reg);
	}


	dt2 = dt*dt;
	dt3 = dt*dt2;
	dt4 = dt*dt3;
	// calculated the final corrected forces
	for (int dim=0; dim<Dim; dim++) {
		da_dt2  = (   a_reg[dim][0] - a0_reg[dim]   ) / dt2;
		adot_dt = (   a_reg[dim][1] + a0dot_reg[dim]) / dt;

		a2 =  -6*da_dt2 - 2*adot_dt - 2*a_reg[dim][1]/dt;
		a3 = (12*da_dt2 + 6*adot_dt)/dt;

		a_reg[dim][2] = a2;
		a_reg[dim][3] = a3;

		da_dt2  = (   a_irr[dim][0] - a0_irr[dim]   ) / dt2;
		adot_dt = (   a_irr[dim][1] + a0dot_irr[dim]) / dt;

		a2 =  -6*da_dt2  - 2*adot_dt - 2*a_irr[dim][1]/dt;
		a3 = (12*da_dt2 + 6*adot_dt)/dt;

		a_irr[dim][2] = a2;
		a_irr[dim][3] = a3;

		a_tot[dim][2] = a_reg[dim][2] + a_irr[dim][2];
		a_tot[dim][3] = a_reg[dim][3] + a_irr[dim][3];
	}



	/*
	double ATOT=0;
	for (int dim=0; dim<Dim; dim++)	 {
		ATOT += a_tot[dim][0]*a_tot[dim][0];
	}
	ATOT = std::sqrt(ATOT)*position_unit/time_unit/time_unit;
	*/
	//std::cout << ATOT << std::endl;
	//fprintf(stdout, "NBODY+: a_tot = %e\n", ATOT);
}




