#include "../global.h"
#include <vector>
#include <iostream>
#include <cmath>





/*
 *  Purporse: Predict particle positions and velocities up to second order 
 *            using a_0 and d(a_0)/dt; refer to Nbody6++ manual Eq. 8.3 and 8.4
 *
 *  Date    : 2024.04.29  by Seoyoung Kim
 *  Date    : 2024.01.09  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *
 */
void Particle::predictParticleSecondOrder(double time) {

	// Doubling check
	// temporary variables for calculation
	if (this->CurrentTimeReg == time) {
		for (int dim=0; dim<Dim; dim++) {
			PredPosition[dim] = Position[dim];
			PredVelocity[dim] = Velocity[dim];
		}
		return;
	}

	double dt;
	//dt = (time - this->CurrentTimeIrr)*EnzoTimeStep;
	dt = (time - this->CurrentTimeReg)*EnzoTimeStep;

	// only predict the positions if necessary
	// how about using polynomial correction here?
	for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] = ((a_tot[dim][1]*dt/3 + a_tot[dim][0])*dt/2 + Velocity[dim])*dt + Position[dim];
		PredVelocity[dim] =  (a_tot[dim][1]*dt/2 + a_tot[dim][0])*dt   + Velocity[dim];
	}

	return;
}



void Particle::predictParticleSecondOrderIrr(double time) {

	// Doubling check
	// temporary variables for calculation
	if (this->CurrentTimeIrr == time) {
		for (int dim=0; dim<Dim; dim++) {
			PredPosition[dim] = Position[dim];
			PredVelocity[dim] = Velocity[dim];
		}
		return;
	}

	double dt;
	//dt = (time - this->CurrentTimeIrr)*EnzoTimeStep;
	dt = (time - this->CurrentTimeIrr)*EnzoTimeStep;

	// only predict the positions if necessary
	// how about using polynomial correction here?
	for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] = ((a_tot[dim][1]*dt/3 + a_tot[dim][0])*dt/2 + Velocity[dim])*dt + Position[dim];
		PredVelocity[dim] =  (a_tot[dim][1]*dt/2 + a_tot[dim][0])*dt   + Velocity[dim];
	}

	return;
}

/*
 *  Purporse: Correct particle positions and velocities up to fourth order
 *            using a_p and d(a_p)/dt; refer to Nbody6++ manual Eq. 8.9 and 8.10
 *
 *  Date    : 2024.01.10  by Yongseok Jo
 *  Modified: 2024.01.11  by Seoyoung Kim
 *
 */
void Particle::correctParticleFourthOrder(double current_time, double next_time, double a[3][4]) {

	if (current_time == next_time) {
		for (int dim=0; dim<Dim; dim++) {
			PredPosition[dim] = Position[dim];
			PredVelocity[dim] = Velocity[dim];
		}
		return;
	}

	double dt;
	double dt3,dt4,dt5;

	dt = (next_time - current_time)*EnzoTimeStep;

	dt3 = dt*dt*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;

	// correct the predicted values positions and velocities at next_time
	// and save the corrected values to particle positions and velocities
	// the latest values of a2dot and a3dots (obtained from hermite method) are used
	for (int dim=0; dim<Dim; dim++) {
		Position[dim] = PredPosition[dim]+ a[dim][2]*dt4/24 + a[dim][3]*dt5/120;
		Velocity[dim] = PredVelocity[dim]+ a[dim][2]*dt3/6  + a[dim][3]*dt4/24;
	}
}



void Particle::updateParticle() {
	for (int dim=0; dim<Dim; dim++) {
		Position[dim] = NewPosition[dim];
		Velocity[dim] = NewVelocity[dim];
	}
}

/*
void Particle::updateParticle(double current_time, double next_time, double a[3][4]) {

	if (current_time == next_time) {
		for (int dim=0; dim<Dim; dim++) {
			PredPosition[dim] = Position[dim];
			PredVelocity[dim] = Velocity[dim];
		}
		return;
	}
	predictParticleSecondOrder(current_time, next_time, a);
	correctParticleFourthOrder(current_time, next_time, a);
	Mass += evolveStarMass(current_time, next_time); // CurrentTimeIrr
}*/





