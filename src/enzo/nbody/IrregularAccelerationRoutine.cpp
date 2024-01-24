#include <iostream>
#include "global.h"


void UpdateEvolveParticle(std::vector<Particle*> &particle, std::vector<Particle*> &list);

// 1. sort particles according to next irregular time
// 2. the sorted list should be update regularly
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle)
{

    std::cout << "Creating a computation chain ...\n";
    if (CreateComputationChain(particle) == false) {
        std::cout << "No irregular particle to update ...\n";
        return true;
    }


    std::cout << "Calculating irregular force ...\n"
              << std::flush;
    for (Particle *ptcl : ComputationChain) 
    {
        fprintf(stdout, "PID=%d, NextRegTime= %e, NextIrrTime = %e\n",
                ptcl->getPID(), MinRegTime, ptcl->CurrentTimeIrr + ptcl->TimeStepIrr);
        fprintf(stdout, "CurrentTimeIrr = %e, TimeStepIrr = %e, CurrentTimeReg=%e, TimeStepReg=%e\n",
                ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->CurrentTimeReg, ptcl->TimeStepReg);

        ptcl->calculateIrrForce(); // this includes particle position and time evolution.
        SortComputationChain(ComputationChain);
    }
}




void UpdateEvolveParticle(std::vector<Particle*> &particle, std::vector<Particle*> &list) {
	std::cout << "nbody+: Updating EvolveParticle..." << std::endl;
	double time, next_time;
	list.clear();
	for (Particle* ptcl: particle) {
		time      = ptcl->CurrentTimeIrr;
		next_time = ptcl->CurrentTimeIrr + ptcl->TimeStepIrr;

		if ((MinRegTime >= next_time) // Regular timestep
				&& ptcl->checkNeighborForEvolution()) { // neighbor
			ptcl->isEvolve = 1;
			list.push_back(ptcl);
		}

		//set global time as the time of a particle that has most advanced
		if (time > global_time)
			global_time = time;
	}
	std::cout << "nbody+: EvolveParticle: ";
	std::cout << list.size();
	std::cout << '\n' << std::flush;
}