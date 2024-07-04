#include <iostream>
#include <vector>
#include <cmath>
#include "../global.h"
#include "../../macros_and_parameters.h"

double star_feedback7(double t1, double t2, Particle *self);

double Particle::evolveStarMass(double t1, double t2) {

	if ((isStarEvolution == false) || (t2 == 0)) {
		return 0;
	}

	double dm;

	switch (StarParticleFeedback) {
		case 0:
			dm = 0.;
			isStarEvolution = false;
			break;

		case 128:
			dm = star_feedback7(t1*EnzoTimeStep, t2*EnzoTimeStep, this);
			break;

		default:
			std::cout << "EvolveStarMass: Invalid choice" << std::endl;
	}

	return dm;
}

// Everything here is in enzo unit
double star_feedback7(double t1, double t2, Particle *self) {
	double dm, tau1, tau2, m1, m2;

	t1 += EnzoCurrentTime;
	t2 += EnzoCurrentTime;

	tau1 = (t1 - self->CreationTime)/self->DynamicalTime;
	tau2 = (t2 - self->CreationTime)/self->DynamicalTime;

	dm = self->InitialMass*((1+tau1)*std::exp(-tau1)-(1+tau2)*std::exp(-tau2));
	dm = -max(min(dm, self->Mass), 0.)*StarMassEjectionFraction;

	if (tau1 > 12.0) {
		return 0.;
	}

	return dm;
}

