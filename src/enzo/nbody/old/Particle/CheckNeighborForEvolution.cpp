#include "../global.h"


bool Particle::checkNeighborForEvolution() {
	if (NumberOfAC == 0) return true;
	int count=0;
	for (Particle* ptcl:ACList) {
		if ((this->CurrentTimeIrr+this->TimeStepIrr) <= (ptcl->CurrentTimeIrr+ptcl->TimeStepIrr)) {
			count++;
		}
	}
	if (count == NumberOfAC) return true;
	else return false;
}
