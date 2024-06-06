#include <iostream>
#include <cmath>
#include "global.h"

void UpdateNextRegTime(std::vector<Particle*> &particle);
void CalculateRegAccelerationOnGPU(std::vector<int> IndexList, std::vector<Particle*> &particle);

bool RegularAccelerationRoutine(std::vector<Particle*> &particle)
{
	std::cout << "Calculating regular force ..." << std::endl;

	// Calulating regular acceleration of the particles
#define GPU
#ifdef GPU
	if (RegIndexList.size() > 0) {
		for (int i: RegIndexList) {
			fprintf(stdout, "PID=%d, CurrentTimeIrr = %.2e Myr, CurrentTimeReg= %.2e Myr, NextRegTime=%llu Myr,\n"\
					"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n",
					particle[i]->PID,
					particle[i]->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
					particle[i]->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
					NextRegTimeBlock,
					particle[i]->TimeStepIrr*EnzoTimeStep*1e10/1e6,
					particle[i]->TimeStepReg*EnzoTimeStep*1e10/1e6,
					particle[i]->TimeBlockIrr,
					particle[i]->TimeLevelIrr,
					particle[i]->TimeBlockReg,
					particle[i]->TimeLevelReg
					);

					/*
					 a_tot = (%.2e,%.2e,%.2e), a_reg = (%.2e,%.2e,%.2e), a_irr = (%.2e,%.2e,%.2e), n_n=%d\n\
					a1_reg = (%.2e,%.2e,%.2e), a2_reg = (%.2e,%.2e,%.2e), a3_reg = (%.2e,%.2e,%.2e)\n\
					a1_irr = (%.2e,%.2e,%.2e), a2_irr = (%.2e,%.2e,%.2e), a3_irr = (%.2e,%.2e,%.2e)\n, 
					particle[i]->a_tot[0][0],
					particle[i]->a_tot[1][0],
					particle[i]->a_tot[2][0],
					particle[i]->a_reg[0][0],
					particle[i]->a_reg[1][0],
					particle[i]->a_reg[2][0],
					particle[i]->a_irr[0][0],
					particle[i]->a_irr[1][0],
					particle[i]->a_irr[2][0],
					particle[i]->NumberOfAC,
					particle[i]->a_reg[0][1],
					particle[i]->a_reg[1][1],
					particle[i]->a_reg[2][1],
					particle[i]->a_reg[0][2],
					particle[i]->a_reg[1][2],
					particle[i]->a_reg[2][2],
					particle[i]->a_reg[0][3],
					particle[i]->a_reg[1][3],
					particle[i]->a_reg[2][3],
					particle[i]->a_irr[0][1],
					particle[i]->a_irr[1][1],
					particle[i]->a_irr[2][1],
					particle[i]->a_irr[0][2],
					particle[i]->a_irr[1][2],
					particle[i]->a_irr[2][2],
					particle[i]->a_irr[0][3],
					particle[i]->a_irr[1][3],
					particle[i]->a_irr[2][3]
					);
					*/
		}
		fflush(stdout);
		CalculateRegAccelerationOnGPU(RegIndexList, particle);
	}
#else
	for (Particle *ptcl : particle) {
		if (ptcl->isRegular)
		{
			fprintf(stdout, "Particle ID=%d, Time = %.4e, dtIrr = %.4e, dtReg = %.4e\n",
				 	ptcl->getPID(), ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->TimeStepReg);
			// This only computes accelerations without updating particles.
			ptcl->calculateRegAccelerationSecondOrder(particle);
		}
	}

	for (Particle *ptcl : particle) {
		if (ptcl->isRegular)
		{
			ptcl->calculateRegAccelerationFourthOrder(particle);
		}
	}

	// Update particles
	for (Particle *ptcl : particle) {
		// update the regular time step
		if (ptcl->isRegular) {
			if ((ptcl->NumberOfAC == 0) && (NextRegTimeBlock == ptcl->CurrentBlockReg + ptcl->TimeBlockReg))
			{
				ptcl->updateParticle(NextRegTimeBlock, ptcl->a_tot);
				ptcl->CurrentTimeReg += ptcl->TimeStepReg;
				ptcl->CurrentTimeIrr  = ptcl->CurrentTimeReg;
			}
			else
			{ // Not sure about it
				ptcl->CurrentTimeReg = ptcl->CurrentTimeIrr;
			}

			ptcl->calculateTimeStepReg(ptcl->a_reg, ptcl->a_reg);
			ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
			//std::cout << "NBODY+: time step = " <<  ptcl->TimeStepReg*EnzoTimeStep << std::endl;
			//std::cout << "NBODY+: time step = " <<  ptcl->TimeStepReg*EnzoTimeStep*time_unit << "yr" << std::endl;
			ptcl->isRegular = false;
		}
	}
#endif



	global_time = NextRegTimeBlock*time_step;
	// update the next regular time step
	UpdateNextRegTime(particle);
	//std::cout << "Finishing regular force ...\n" << std::flush;
	return true;
}



void UpdateNextRegTime(std::vector<Particle*> &particle) {

	ULL time_tmp, time=block_max;

	for (Particle *ptcl : particle)
	{
		// Next regular time step
		time_tmp = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;

		// Find the minum regular time step
		if (time_tmp < time)
			time = time_tmp;
	}

	NextRegTimeBlock = std::min(time, block_max);

	// Set isRegular of the particles that will be updated next to 1
	//std::cerr << "Regular: ";
	RegIndexList.clear();
	int i = 0;
	for (Particle* ptcl: particle) {
		time_tmp = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;
		if (NextRegTimeBlock == time_tmp) {
			//std::cerr << ptcl->getPID() << ' ';
			ptcl->isRegular = true;
			RegIndexList.push_back(i);
		}
		i++;
	}
	//std::cerr << '\n' << std::flush;
}

