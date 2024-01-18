#include "../global.h"


void Particle::updateEvolveParticle(std::vector<Particle*> &list, double MinRegTime) {

	double next_time;

	next_time = CurrentTimeIrr + TimeStepIrr;

	if (MinRegTime >= next_time 
			&& (this->checkNeighborForEvolution())) {
		return;
	} else {
		this->isEvolve = 0;
		int i=0;
		for (Particle* ptcl:list) {
			if (ptcl == this) {
				list.erase(list.begin()+i);
			}
			//	break;
			i++;
		}
	}
}




