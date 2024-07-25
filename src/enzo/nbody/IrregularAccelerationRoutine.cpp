#include <iostream>
#include "global.h"
#include <cassert>


std::vector<Particle*> ComputationList{};
int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int outputNum);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool UpdateComputationChain(Particle* ptcl);
bool CreateComputationList(Particle* ptcl);
bool AddNewBinariesToList(std::vector<Particle*> &particle);
void BinaryAccelerationRoutine(double next_time, std::vector<Particle*> &particle);
void KSTermination(Particle* ptclCM, std::vector<Particle*> &particle, double current_time, ULL current_block);



Particle *FirstComputation;
// 1. sort particles according to next irregular time
// 2. the sorted list should be update regularly
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle)
{
	bool bin_termination;
	fprintf(nbpout, "Irregular force starts...\n");
	//fflush(nbpout);
#ifdef time_trace
	_time.irr_chain.markStart();
#endif
	Particle *ParticleForComputation;

	if (CreateComputationChain(particle) == false) {
		fprintf(nbpout, "No irregular particle to update ...\n");
		return true;
	}
	fprintf(nbpout, "ComputationChain done ...\n");
	fflush(nbpout);

#ifdef time_trace
	_time.irr_chain.markEnd();
	_time.irr_chain.getDuration();
#endif
	// Caculating irregular acceleration


	while (CreateComputationList(FirstComputation) && ComputationList.size() != 0) {

		if (BinaryRegularization) {
#ifdef time_trace
	_time.irr_bin.markStart();
#endif
			if (AddNewBinariesToList(particle) && ComputationList.size() == 0) {
				fprintf(nbpout, "No irregular particle to update afte binary formation.\n");
				break;
			}

			if ((BinaryList.size()>0)&(binary_time_prev != binary_time)) {
				fprintf(binout, "-------------------------------------\n");
				fprintf(binout, "irr_time = %e \n",
						ComputationList[0]->CurrentTimeIrr*1e10/1e6*EnzoTimeStep);
				fprintf(binout, "binary_time = %e \n",
						binary_time*1e10/1e6*EnzoTimeStep);
				fprintf(binout, "Evolve.cpp: integrating binaries\n");
				fprintf(binout, "# of binaries = %d \n",int(BinaryList.size()));
				BinaryAccelerationRoutine(
						ComputationList[0]->CurrentTimeIrr+ComputationList[0]->TimeStepIrr,
				 	particle);
			}
#ifdef time_trace
	_time.irr_bin.markEnd();
	_time.irr_bin.getDuration();
#endif
		}

		for (Particle* ptcl:ComputationList) {

			/*
			if (ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr || ptcl->CurrentBlockReg+ptcl->TimeBlockReg < ptcl->CurrentBlockIrr)
				fprintf(nbpout,"--------------------error--------------------------------------------------------------------\n");

			fprintf(nbpout, "PID=%d, NextRegTimeBlock= %llu, NextIrrBlock = %llu (%.2e Myr)\n",
					ptcl->getPID(), NextRegTimeBlock,
					ptcl->CurrentBlockIrr+ptcl->TimeBlockIrr,
					(ptcl->CurrentTimeIrr + ptcl->TimeStepIrr)*EnzoTimeStep*1e10/1e6);
			fprintf(nbpout, "CurrentTimeIrr = %.2e Myr (%llu), CurrentTimeReg = %.2e Myr (%llu), TimeStepIrr = %.2e Myr (%llu), TimeStepReg= %.2e Myr (%llu)\n",
					ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6, ptcl->CurrentBlockIrr,
					ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6, ptcl->CurrentBlockReg,
					ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6, ptcl->TimeBlockIrr, 
					ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6, ptcl->TimeBlockReg
					); 

			if (ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr || ptcl->CurrentBlockReg+ptcl->TimeBlockReg < ptcl->CurrentBlockIrr)
				fprintf(nbpout,"----------------------------------------------------------------------------------------\n");
			//fflush(nbpout);


			fprintf(nbpout, "a_tot = (%.2e,%.2e,%.2e), a_irr = (%.2e,%.2e,%.2e), a1_irr = (%.2e,%.2e,%.2e), a2_irr = (%.2e,%.2e,%.2e), a3_irr = (%.2e,%.2e,%.2e), n_n=%d\n",
					ptcl->a_tot[0][0],
					ptcl->a_tot[1][0],
					ptcl->a_tot[2][0],
					ptcl->a_irr[0][0],
					ptcl->a_irr[1][0],
					ptcl->a_irr[2][0],
					ptcl->a_irr[0][1],
					ptcl->a_irr[1][1],
					ptcl->a_irr[2][1],
					ptcl->a_irr[0][2],
					ptcl->a_irr[1][2],
					ptcl->a_irr[2][2],
					ptcl->a_irr[0][3],
					ptcl->a_irr[1][3],
					ptcl->a_irr[2][3],
					ptcl->NumberOfAC);
			fprintf(nbpout, "x = (%.2e,%.2e,%.2e), v = (%.2e,%.2e,%.2e)\n",
					ptcl->Position[0],ptcl->Position[1],ptcl->Position[2],
					ptcl->Velocity[0],ptcl->Velocity[1],ptcl->Velocity[2]
					);
			//fflush(nbpout);
			*/

			/*
				 fprintf(stdout, "in irr, PID=%d : ", ptcl->PID);
				 for (Particle* nn:ptcl->ACList) {
				 fprintf(stdout,"%d, ",nn->PID);	
				 }
				 fprintf(stdout,"\n");	
				 */

			if (mag0(ptcl->a_irr) == 0. && ptcl->NumberOfAC != 0) {
				if (ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr || ptcl->CurrentBlockReg+ptcl->TimeBlockReg < ptcl->CurrentBlockIrr)
					fprintf(nbpout,"--------------------error--------------------------------------------------------------------\n");

				fprintf(nbpout, "PID=%d, NextRegTimeBlock= %llu, NextIrrBlock = %llu (%.2e Myr)\n",
						ptcl->getPID(), NextRegTimeBlock,
						ptcl->CurrentBlockIrr+ptcl->TimeBlockIrr,
						(ptcl->CurrentTimeIrr + ptcl->TimeStepIrr)*EnzoTimeStep*1e10/1e6);
				fprintf(nbpout, "CurrentTimeIrr = %.2e Myr (%llu), CurrentTimeReg = %.2e Myr (%llu), TimeStepIrr = %.2e Myr (%llu), TimeStepReg= %.2e Myr (%llu)\n",
						ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6, ptcl->CurrentBlockIrr,
						ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6, ptcl->CurrentBlockReg,
						ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6, ptcl->TimeBlockIrr, 
						ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6, ptcl->TimeBlockReg
						); 

				if (ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr || ptcl->CurrentBlockReg+ptcl->TimeBlockReg < ptcl->CurrentBlockIrr)
					fprintf(nbpout,"----------------------------------------------------------------------------------------\n");
				//fflush(nbpout);


				fprintf(nbpout, "a_tot = (%.2e,%.2e,%.2e), a_irr = (%.2e,%.2e,%.2e), a1_irr = (%.2e,%.2e,%.2e), a2_irr = (%.2e,%.2e,%.2e), a3_irr = (%.2e,%.2e,%.2e), n_n=%d\n",
						ptcl->a_tot[0][0],
						ptcl->a_tot[1][0],
						ptcl->a_tot[2][0],
						ptcl->a_irr[0][0],
						ptcl->a_irr[1][0],
						ptcl->a_irr[2][0],
						ptcl->a_irr[0][1],
						ptcl->a_irr[1][1],
						ptcl->a_irr[2][1],
						ptcl->a_irr[0][2],
						ptcl->a_irr[1][2],
						ptcl->a_irr[2][2],
						ptcl->a_irr[0][3],
						ptcl->a_irr[1][3],
						ptcl->a_irr[2][3],
						ptcl->NumberOfAC);
				fprintf(nbpout, "x = (%.2e,%.2e,%.2e), v = (%.2e,%.2e,%.2e)\n",
						ptcl->Position[0],ptcl->Position[1],ptcl->Position[2],
						ptcl->Velocity[0],ptcl->Velocity[1],ptcl->Velocity[2]
						);
				//fflush(nbpout);
				fprintf(stderr, "my PID = %d : ",ptcl->PID);
				fprintf(nbpout, "my PID = %d : ",ptcl->PID);
				for (Particle* nn:ptcl->ACList) {
					fprintf(stderr, "(%d, %e, %e), ",nn->PID, nn->Mass, nn->a_irr[0][0]);
					fprintf(nbpout, "%d ",nn->PID);
				}
				fprintf(stderr, "\n");
				fprintf(nbpout, "\n");
				fflush(stderr);
				fflush(nbpout);
				throw std::runtime_error("IrregularAccelerationRoutine.cpp");
			}

#ifdef time_trace
			_time.irr_force.markStart();
#endif
			ptcl->calculateIrrForce(); // this includes particle position

#ifdef time_trace
			_time.irr_force.markEnd();
			_time.irr_force.getDuration();
#endif

		}




#ifdef time_trace
		_time.irr_sort.markStart();
#endif
		//std::cout << "Update and Sort\n" << std::flush;

		// update particles and chain
		// The next particle of the particle calculated lastly should be the start of the next iteration.
		if (BinaryRegularization) {
			binary_time_prev = binary_time;
			binary_block     = ComputationList[0]->CurrentBlockIrr+ComputationList[0]->TimeBlockIrr;
			binary_time      = binary_block*time_step;
			bin_termination  = false;	
		}
		global_time_irr    = ComputationList[0]->CurrentBlockIrr+ComputationList[0]->TimeBlockIrr;
		for (Particle* ptcl:ComputationList) {
			//std::cout << ptcl->PID << " " ;
			ptcl->updateParticle();
			ptcl->CurrentBlockIrr += ptcl->TimeBlockIrr;
			ptcl->CurrentTimeIrr   = ptcl->CurrentBlockIrr*time_step;
			//std::cout << "before TimeStepCal\n" << std::flush;
			ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
			//std::cout << "after TimeStepCal\n" << std::flush;
			if (ptcl->isCMptcl) {
				if (BinaryRegularization && 
						(ptcl->BinaryInfo->r > ptcl->BinaryInfo->r0*2.0
						|| ptcl->BinaryInfo->TimeStep > 2.0*KSTime)) {
					fprintf(binout, "Terminating Binary at time : %e \n", binary_time);
					KSTermination(ptcl, particle, binary_time, binary_block);
					bin_termination=true;
					continue;
				}
			}
			UpdateComputationChain(ptcl);
		}

		// because NextRegTime might have been changed.
		if (BinaryRegularization && bin_termination)
			CreateComputationChain(particle);

#ifdef time_trace
		_time.irr_sort.markEnd();
		_time.irr_sort.getDuration();
#endif
		//fflush(stdout);
	}
	fprintf(nbpout, "Finishing irregular force ...\n");
	//fflush(nbpout);


	return true;
}



