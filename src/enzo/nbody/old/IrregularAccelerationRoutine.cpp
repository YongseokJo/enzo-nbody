#include <iostream>
#include "global.h"


bool CreateComputationChain(std::vector<Particle*> &particle);
bool SortComputationChain(std::vector<Particle*> ComputationChain);
Particle *SortComputationChain(Particle* ptcl);

Particle *FirstComputation;
// 1. sort particles according to next irregular time
// 2. the sorted list should be update regularly
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle)
{

		Particle *ParticleForComputation;
    std::cout << "Creating a computation chain ...\n";
    if (CreateComputationChain(particle) == false) {
        std::cout << "No irregular particle to update ...\n";
        return true;
    }

    std::cout << "Calculating irregular force ...\n" << std::endl;

		// Caculating irregular acceleration
		ParticleForComputation = FirstComputation;
		while (ParticleForComputation != nullptr) {

			fprintf(stdout, "PID=%d, NextRegTime= %e, NextIrrTime = %e\n",
					ParticleForComputation->getPID(), NextRegTime, ParticleForComputation->CurrentTimeIrr + ParticleForComputation->TimeStepIrr);
			fprintf(stdout, "CurrentTimeIrr = %e, TimeStepIrr = %e, CurrentTimeReg=%e, TimeStepReg=%e\n",
					ParticleForComputation->CurrentTimeIrr, ParticleForComputation->TimeStepIrr, ParticleForComputation->CurrentTimeReg, ParticleForComputation->TimeStepReg);


			ParticleForComputation->calculateIrrForce(); // this includes particle position and time evolution.
			//ParticleForComputation = ParticleForComputation->NextParticleForComputation;
			ParticleForComputation = SortComputationChain(ParticleForComputation);
		}

		/*
		for (Particle *ptcl : ComputationChain)
		{
			fprintf(stdout, "PID=%d, NextRegTime= %e, NextIrrTime = %e\n",
					ptcl->getPID(), NextRegTime, ptcl->CurrentTimeIrr + ptcl->TimeStepIrr);
			fprintf(stdout, "CurrentTimeIrr = %e, TimeStepIrr = %e, CurrentTimeReg=%e, TimeStepReg=%e\n",
					ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->CurrentTimeReg, ptcl->TimeStepReg);

			ptcl->calculateIrrForce(); // this includes particle position and time evolution.
			SortComputationChain(ComputationChain);
		}
		*/

    std::cout << "Finishing irregular force ...\n" << std::endl;
		return true;
}



