#include "Particle.h"


bool Particle::checkNeighborForEvolution() {
	int count=0;
	for (Particle* ptcl:ACList) {
		if ((this->CurrentTimeIrr+this->TimeStepIrr - mytolerance) <=ptcl->CurrentTimeIrr+ptcl->TimeStepIrr) {
			count++;
		}
	}
	if (count == NumberOfAC) return true;
	else return false;
}
