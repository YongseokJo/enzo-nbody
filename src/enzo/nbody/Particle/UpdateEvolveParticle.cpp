#include "Particle.h"


void Particle::updateEvolveParticle(std::vector<Particle*> &list, double MinRegTime) {
	double cond;
	cond	= MinRegTime-(CurrentTimeIrr+TimeStepIrr);
	//if (PID == 227) {
	fprintf(stderr,"UpdateEVolve, PID=%d\n",PID);
	fprintf(stderr, "RegTime=%e, CurrentTime=%e, TimeStep=%e\n",MinRegTime, CurrentTimeIrr,TimeStepIrr);
	std::cerr << cond <<std::flush;
	fprintf(stderr,"\n");
	//}
	if (MinRegTime >= (this->CurrentTimeIrr+this->TimeStepIrr - mytolerance)
			&& (this->checkNeighborForEvolution())) {
		return;
	} else {
		this->isEvolve = 0;
		int i=0;
		for (Particle* ptcl:list) {
			if (PID == 227) 
				std::cerr << ptcl->getPID() << ' ';
			if (ptcl == this) {
				list.erase(list.begin()+i);
			}
			//	break;
			i++;
		}
		//std::cerr << i <<std::flush;
		fprintf(stderr,"\n");
	}
}




