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

#ifdef time_trace
		_time.irr_chain.markStart();
#endif
		Particle *ParticleForComputation;
    //std::cout << "Creating a computation chain ...\n";
    if (CreateComputationChain(particle) == false) {
        std::cout << "No irregular particle to update ...\n";
        return true;
    }

    //std::cout << "Calculating irregular force ...\n" << std::endl;
#ifdef time_trace
		_time.irr_chain.markEnd();
		_time.irr_chain.getDuration();
#endif
		// Caculating irregular acceleration
		ParticleForComputation = FirstComputation;
		while (ParticleForComputation != nullptr) {

			fprintf(stdout, "PID=%d, NextRegTime= %e Myr, NextIrrTime = %e Myr\n",
					ParticleForComputation->getPID(), NextRegTime*EnzoTimeStep*1e10/1e6, (ParticleForComputation->CurrentTimeIrr + ParticleForComputation->TimeStepIrr)*EnzoTimeStep*1e10/1e6);
			fprintf(stdout, "CurrentTimeIrr = %e Myr, TimeStepIrr = %e Myr, CurrentTimeReg=%e Myr, TimeStepReg=%e Myr\n",
					ParticleForComputation->CurrentTimeIrr*EnzoTimeStep*1e10/1e6, ParticleForComputation->TimeStepIrr*EnzoTimeStep*1e10/1e6, ParticleForComputation->CurrentTimeReg*EnzoTimeStep*1e10/1e6, ParticleForComputation->TimeStepReg*EnzoTimeStep*1e10/1e6);


#ifdef time_trace
		_time.irr_force.markStart();
#endif
			ParticleForComputation->calculateIrrForce(); // this includes particle position and time evolution.
#ifdef time_trace
		_time.irr_force.markEnd();
		_time.irr_force.getDuration();
#endif
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

    //std::cout << "Finishing irregular force ...\n" << std::endl;
		return true;
}



