#include <iostream>
#include "global.h"


std::vector<Particle*> ComputationList{};
int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int outputNum);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool SortComputationChain(std::vector<Particle*> ComputationChain);
Particle *SortComputationChain(Particle* ptcl);
bool CreateComputationList(Particle* ptcl);

Particle *FirstComputation;
// 1. sort particles according to next irregular time
// 2. the sorted list should be update regularly
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle)
{

#ifdef time_trace
		_time.irr_chain.markStart();
#endif
		std::cout << "Create Chain\n" << std::flush;
		Particle *ParticleForComputation;
    //std::cout << "Creating a computation chain ...\n";
    if (CreateComputationChain(particle) == false) {
        std::cout << "No irregular particle to update ...\n";
        return true;
    }
		std::cout << "Create List\n" << std::flush;
		ParticleForComputation = FirstComputation;

		if (CreateComputationList(ParticleForComputation) == false) {
			std::cout << "Something's wrong\n";
		}

		std::cout << "List size=" << ComputationList.size() << std::endl;

    //std::cout << "Calculating irregular force ...\n" << std::endl;
#ifdef time_trace
		_time.irr_chain.markEnd();
		_time.irr_chain.getDuration();
#endif
		// Caculating irregular acceleration

		computation:
		std::cout << "Start IRR\n" << std::flush;
		//while (ParticleForComputation != nullptr) {
		for (Particle* ParticleForComputation:ComputationList) {

			/*
			fprintf(stdout, "PID=%d, NextRegTimeBlock= %llu, NextIrrBlock = %llu (%.2e Myr)\n",
					ParticleForComputation->getPID(), NextRegTimeBlock, ParticleForComputation->CurrentBlockIrr+ParticleForComputation->TimeBlockIrr, (ParticleForComputation->CurrentTimeIrr + ParticleForComputation->TimeStepIrr)*EnzoTimeStep*1e10/1e6);
			fprintf(stdout, "CurrentTimeIrr = %.2e Myr, TimeStepIrr = %.2e Myr, CurrentBlockReg=%llu, TimeBlockReg=%llu\n",
					ParticleForComputation->CurrentTimeIrr*EnzoTimeStep*1e10/1e6, ParticleForComputation->TimeStepIrr*EnzoTimeStep*1e10/1e6, ParticleForComputation->CurrentBlockReg, ParticleForComputation->TimeBlockReg);
					*/

			/*

			fprintf(stdout, "a_tot = (%.2e,%.2e,%.2e), a_irr = (%.2e,%.2e,%.2e), a1_irr = (%.2e,%.2e,%.2e), a2_irr = (%.2e,%.2e,%.2e), a3_irr = (%.2e,%.2e,%.2e), n_n=%d\n",
					ParticleForComputation->a_tot[0][0],
					ParticleForComputation->a_tot[1][0],
					ParticleForComputation->a_tot[2][0],
					ParticleForComputation->a_irr[0][0],
					ParticleForComputation->a_irr[1][0],
					ParticleForComputation->a_irr[2][0],
					ParticleForComputation->a_irr[0][1],
					ParticleForComputation->a_irr[1][1],
					ParticleForComputation->a_irr[2][1],
					ParticleForComputation->a_irr[0][2],
					ParticleForComputation->a_irr[1][2],
					ParticleForComputation->a_irr[2][2],
					ParticleForComputation->a_irr[0][3],
					ParticleForComputation->a_irr[1][3],
					ParticleForComputation->a_irr[2][3],
					ParticleForComputation->NumberOfAC);
					*/
					

#ifdef time_trace
		_time.irr_force.markStart();
#endif
			ParticleForComputation->calculateIrrForce(); // this includes particle position


#ifdef time_trace
		_time.irr_force.markEnd();
		_time.irr_force.getDuration();
#endif


			//ParticleForComputation = ParticleForComputation->NextParticleForComputation;// this includes time evolution.
			//ParticleForComputation = SortComputationChain(ParticleForComputation);



#define no_IRR_TEST
#ifdef IRR_TEST
			// create output at appropriate time intervals just for irr
			if (outputTime <= particle[0]->CurrentTimeIrr ) {
				writeParticle(particle, particle[0]->CurrentTimeIrr, outNum++);
				outputTime += outputTimeStep;
			}
#endif

		}

		std::cout << "Update and Sort\n" << std::flush;
		// update particles and chain
		for (Particle* ptcl:ComputationList) {
			ptcl->updateParticle();
			ptcl->CurrentBlockIrr += ptcl->TimeBlockIrr;
			ptcl->CurrentTimeIrr   = ptcl->CurrentBlockIrr*time_step;
			ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
			ParticleForComputation = SortComputationChain(ptcl);
		}
#ifdef time_trace
		_time.irr_sort.markStart();
#endif
		std::cout << "CreateComputationList\n" << std::flush;
		if (CreateComputationList(ParticleForComputation) && ComputationList.size() != 0) {
			goto computation;
		}
#ifdef time_trace
		_time.irr_sort.markEnd();
		_time.irr_sort.getDuration();
#endif
    std::cout << "Finishing irregular force ...\n" << std::endl;


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



