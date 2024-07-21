#include "../global.h"
#include "../defs.h"
#include <algorithm>
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
		PredMass = Mass;
		return;
	}

	double dt;
	//dt = (time - this->CurrentTimeIrr)*EnzoTimeStep;
	dt = (time - this->CurrentTimeReg)*EnzoTimeStep;

	// only predict the positions if necessary
	// how about using polynomial correction here?
	for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] = ((a_tot[dim][1]*dt/3 + a_tot[dim][0] + BackgroundAcceleration[dim])*dt/2 + Velocity[dim])*dt + Position[dim];
		PredVelocity[dim] =  (a_tot[dim][1]*dt/2 + a_tot[dim][0] + BackgroundAcceleration[dim])*dt   + Velocity[dim];
	}

	if (StarParticleFeedback != 0) {
		PredMass = Mass + evolveStarMass(CurrentTimeIrr, CurrentTimeIrr+TimeStepIrr*1.01);
		Mdot     = (PredMass - Mass)/TimeStepIrr*1e-2;
	}
	else {
		PredMass = Mass;
		Mdot     = 0;
	}

			////TimeStepIrr*1e-2; // derivative can be improved

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
		PredMass = Mass;
		return;
	}

	double dt;
	//dt = (time - this->CurrentTimeIrr)*EnzoTimeStep;
	dt = (time - this->CurrentTimeIrr)*EnzoTimeStep;

	// only predict the positions if necessary
	// how about using polynomial correction here?
	for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] = ((a_tot[dim][1]*dt/3 + (a_tot[dim][0]+BackgroundAcceleration[dim]))*dt/2 + Velocity[dim])*dt + Position[dim];
		PredVelocity[dim] =  (a_tot[dim][1]*dt/2 + (a_tot[dim][0]+BackgroundAcceleration[dim]))*dt   + Velocity[dim];
	}

	if (StarParticleFeedback != 0) {
		PredMass = Mass + evolveStarMass(CurrentTimeIrr, CurrentTimeIrr+TimeStepIrr*1.01);
		Mdot     = (PredMass - Mass)/TimeStepIrr*1e-2;
	}
	else {
		PredMass = Mass;
		Mdot     = 0;
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


void Particle::polynomialPrediction(double current_time) {

}

void Particle::updateParticle() {
	for (int dim=0; dim<Dim; dim++) {
		Position[dim] = NewPosition[dim];
		Velocity[dim] = NewVelocity[dim];
	}
	Mass = PredMass;
}




void Particle::UpdateRadius() {

	if (LocalDensity == 0) {
		//RadiusOfAC = 0.11;
		RadiusOfAC = InitialRadiusOfAC;
		//RadiusOfAC = 1.00;
		LocalDensity = 10;
	}
	else {
		/* exponential (aggressive) */
		/*
			 const double c = 0.5;
			 const double b = std::log(2) / (NumNeighborMax);  // ln(2) / 40
			 double exp = a * (std::exp(b * NumberOfAC) - 1);
			 */

		/* n=2 polynomial (mild) as n increases it grows mild */
		const int n = 3;
		const double c = (NumNeighborMax-FixNumNeighbor);
		const double b = 0.5 / pow(c,n);  // ln(2) / 40
		double x = NumberOfAC-FixNumNeighbor;
		//double a = n%2==0 ? b*std::abs(x)*pow(x,n-1) : b*pow(x,n);
		double a = n%2==0 ?  b*pow(x,n) : b*std::abs(x)*pow(x,n-1);
		RadiusOfAC *= (1-a);
	}
	/*
	if (NumberOfAC > FixNumNeighbor) {
		if (NumberOfAC > 2*FixNumNeighbor) {
			RadiusOfAC *= 0.90;
		}
		else if (NumberOfAC > 3*FixNumNeighbor) {
			RadiusOfAC *= 0.80;
		}
		else {
			RadiusOfAC *= 0.95;
		}
	}
	if (NumberOfAC < FixNumNeighbor)
		RadiusOfAC *= 1.05;
		*/

	/*
	double MeanRadius=0, TotalMass=0, LocalDensity0=0;
	for (Particle *neighbor:ACList) {
		MeanRadius += neighbor->Mass*dist(Position,neighbor->Position);
		TotalMass  += neighbor->Mass;
	}
	MeanRadius        /= TotalMass/1.2;  // this is the mean radius with some factor
	LocalDensity0      = TotalMass/std::pow(MeanRadius,3.0);

	if (LocalDensity == 0) {
		this->RadiusOfAC   = std::pow(LocalDensity0/(this->Mass*FixNumNeighbor), -1.0/3.0); 

		if (LocalDensity0 > 1.3*LocalDensity) 
			this->RadiusOfAC *= 0.9;
		else if (LocalDensity0 < 0.7*LocalDensity) 
			this->RadiusOfAC *= 1.1;

		if (NumberOfAC < FixNumNeighbor/2)
			this->RadiusOfAC *= 1.2;

		if (NumberOfAC > FixNumNeighbor*2)
			this->RadiusOfAC *= 0.8;
	}
	*/


	//this->LocalDensity = LocalDensity0;

	//this->RadiusOfAC   = std::pow(LocalDensity0/(this->Mass*FixNumNeighbor), -1.0/3.0); 
	//fprintf(stdout, "PID %d : TWR = %.3e, TM=%.3e, LD0=%.3e, RAC=%.3e, mass=%.3e\n", 
			//PID, MeanRadius, TotalMass, LocalDensity0, RadiusOfAC, Mass);
	//fflush(stdout); 
}


void Particle::UpdateNeighbor(std::vector<Particle*> &particle) {
	bool isExcess=false;
	ACList.clear();
	NumberOfAC = 0;
	for (Particle* ptcl:particle) {

		if (ptcl->PID==PID)
			continue;

		if (dist(Position, ptcl->Position)<this->RadiusOfAC) {
			ACList.push_back(ptcl);	
			NumberOfAC++;
			if (NumberOfAC >= NumNeighborMax) {
				isExcess = true; 
				break;
			}
		}
	}

	if (isExcess) {
		std::cerr << "Number of neighbors exceeded." << std::endl;
		this->RadiusOfAC *= 0.8;
		UpdateNeighbor(particle);	
	}
	else {
		std::sort(ACList.begin(),ACList.end(),
				[](Particle* p1, Particle* p2) { return p1->ParticleOrder < p2->ParticleOrder;});
	}
}



void Particle::NeighborCorrection(double r0, Particle* newPtcl, std::vector<Particle*> &particle) {
	double r1, r_max;
	int index, max_index;
	bool isExcess;
	double a[Dim], adot[Dim];

	if (newPtcl->NumberOfAC < NumNeighborMax) {
		ACList.push_back(newPtcl);
		NumberOfAC++;
		this->ComputeAcceleration(newPtcl,a,adot);
		for (int dim=0; dim<Dim; dim++) {
			this->a_irr[dim][0] += a[dim];
			this->a_irr[dim][1] += adot[dim];
		}
	}
	else {
		index = 0;
		r_max = 0;
		max_index = -1;
		RadiusOfAC *= 0.9;
		for (Particle* neighbor:ACList) {
			if (this->PID == neighbor->PID)
				continue;

			r1 = dist(Position, neighbor->Position);

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
			//fprintf(stderr, "particle removed=%d\n", ptcl->ACList[max_index]->PID);
			this->ComputeAcceleration(ACList[max_index],a,adot);
			for (int dim=0; dim<Dim; dim++) {
				this->a_irr[dim][0] -= a[dim];
				this->a_irr[dim][1] -= adot[dim];
				this->a_reg[dim][0] += a[dim];
				this->a_reg[dim][1] += adot[dim];
			}
			this->ComputeAcceleration(newPtcl,a,adot);
			for (int dim=0; dim<Dim; dim++) {
				this->a_irr[dim][0] += a[dim];
				this->a_irr[dim][1] += adot[dim];
				this->a_reg[dim][0] -= a[dim];
				this->a_reg[dim][1] -= adot[dim];
			}
			ACList.erase(ACList.begin() + max_index);
			ACList.push_back(newPtcl);
		}
	}

}




