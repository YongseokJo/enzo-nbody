#include "global.h"
#include <iostream>
#include <cmath>


int time_block = -30;
//int time_block = 0;
ULL block_max = static_cast<ULL>(pow(2., -time_block));
//double time_step = std::pow(2.,time_block);
double time_step = -1;
double getNewTimeStep(double f[3][4], double df[3][4]);
void getBlockTimeStep(double dt, int& TimeLevel, ULL &TimeBlock, double &TimeStep);

/*
 *  Purporse: Initialize timsteps 
 *
 *  Modified: 2024.01.16  by Yongseok Jo
 *  Modified: 2024.05.29  by Yongseok Jo
 *
 */
int InitializeTimeStep(std::vector<Particle*> &particle) {

	std::cout << "Initializing timesteps ..." << std::endl;

	int min_time_level=0;
	double dtIrr, dtReg;

	for (Particle* ptcl: particle) {
		if (ptcl->NumberOfAC != particle.size()-1) {
			dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
			getBlockTimeStep(dtReg, ptcl->TimeLevelReg, ptcl->TimeBlockReg, ptcl->TimeStepReg);

			ptcl->TimeStepReg  = std::min(1.,ptcl->TimeStepReg);
			ptcl->TimeBlockReg = std::min(block_max, ptcl->TimeBlockReg);
			ptcl->TimeLevelReg = std::min(0, ptcl->TimeLevelReg);
		}

		if (ptcl->NumberOfAC != 0) {
			dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
			getBlockTimeStep(dtIrr, ptcl->TimeLevelIrr, ptcl->TimeBlockIrr, ptcl->TimeStepIrr);

			ptcl->TimeStepIrr  = std::min(1.,ptcl->TimeStepIrr);
			ptcl->TimeBlockIrr = std::min(block_max, ptcl->TimeBlockIrr);
			ptcl->TimeLevelIrr = std::min(0, ptcl->TimeLevelIrr);

			if (ptcl->NumberOfAC == particle.size()-1) {
				ptcl->TimeStepReg  = ptcl->TimeStepIrr;
				ptcl->TimeBlockReg = ptcl->TimeBlockIrr;
				ptcl->TimeLevelReg = ptcl->TimeLevelIrr;
			}
		}
		else {
			ptcl->TimeBlockIrr = ptcl->TimeBlockReg;
			ptcl->TimeLevelIrr = ptcl->TimeLevelReg;
			ptcl->TimeStepIrr  = ptcl->TimeStepReg;
		}

		ptcl->CurrentTimeIrr  = 0;
		ptcl->CurrentTimeReg  = 0;
		ptcl->CurrentBlockIrr = 0;
		ptcl->CurrentBlockReg = 0;
	}



	// Irregular Time Step Correction
	for (Particle* ptcl: particle) {
		if (ptcl->NumberOfAC != 0) {
			while (ptcl->TimeLevelIrr >= ptcl->TimeLevelReg 
					&& ptcl->NumberOfAC != particle.size()-1) {
				ptcl->TimeStepIrr *= 0.5;
				ptcl->TimeBlockIrr *= 0.5;
				ptcl->TimeLevelIrr--;
			}
			if (ptcl->TimeLevelIrr < min_time_level) {
				min_time_level = ptcl->TimeLevelIrr;
			}
		}
	}


	// resetting time_block based on the system
	time_block = std::max(-60, min_time_level-MIN_LEVEL_BUFFER);
	block_max = static_cast<ULL>(pow(2., -time_block));
	time_step = std::pow(2.,time_block);

	for (Particle* ptcl: particle) {
		ptcl->TimeBlockIrr = static_cast<ULL>(pow(2., ptcl->TimeLevelIrr-time_block));
		ptcl->TimeBlockReg = static_cast<ULL>(pow(2., ptcl->TimeLevelReg-time_block));
		if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 < KSTime)
			BinaryCandidateList.push_back(ptcl);
	}


	fprintf(stdout, "nbody+:time_block = %d, EnzoTimeStep=%e\n", time_block, EnzoTimeStep);
	return true;
}



int InitializeTimeStep(std::vector<Particle*> &particle, int offset, int newSize) {
	fprintf(nbpout, "Initializing timesteps ...\n");
	int min_time_level=0;
	double dtIrr, dtReg;
	Particle *ptcl;

	for (int i=offset; i<offset+newSize; i++){
		ptcl = particle[i];
		if (ptcl->NumberOfAC != particle.size()-1) {
			dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
			getBlockTimeStep(dtReg, ptcl->TimeLevelReg, ptcl->TimeBlockReg, ptcl->TimeStepReg);

			ptcl->TimeStepReg  = std::min(1.,ptcl->TimeStepReg);
			ptcl->TimeBlockReg = std::min(block_max, ptcl->TimeBlockReg);
			ptcl->TimeLevelReg = std::min(0, ptcl->TimeLevelReg);
		}

		if (ptcl->NumberOfAC != 0) {
			dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
			getBlockTimeStep(dtIrr, ptcl->TimeLevelIrr, ptcl->TimeBlockIrr, ptcl->TimeStepIrr);

			ptcl->TimeStepIrr  = std::min(1.,ptcl->TimeStepIrr);
			ptcl->TimeBlockIrr = std::min(block_max, ptcl->TimeBlockIrr);
			ptcl->TimeLevelIrr = std::min(0, ptcl->TimeLevelIrr);

			if (ptcl->NumberOfAC == particle.size()-1) {
				ptcl->TimeStepReg  = ptcl->TimeStepIrr;
				ptcl->TimeBlockReg = ptcl->TimeBlockIrr;
				ptcl->TimeLevelReg = ptcl->TimeLevelIrr;
			}
		}
		else {
			ptcl->TimeBlockIrr = ptcl->TimeBlockReg;
			ptcl->TimeLevelIrr = ptcl->TimeLevelReg;
			ptcl->TimeStepIrr  = ptcl->TimeStepReg;
		}

		ptcl->CurrentTimeIrr  = 0;
		ptcl->CurrentTimeReg  = 0;
		ptcl->CurrentBlockIrr = 0;
		ptcl->CurrentBlockReg = 0;
	} // endfor size


	// Irregular Time Step Correction
	for (int i=offset; i<offset+newSize; i++){
		ptcl = particle[i];
		if (ptcl->NumberOfAC != 0) {
			while (ptcl->TimeLevelIrr >= ptcl->TimeLevelReg 
					&& ptcl->NumberOfAC != particle.size()-1) {
				ptcl->TimeStepIrr *= 0.5;
				ptcl->TimeBlockIrr *= 0.5;
				ptcl->TimeLevelIrr--;
			}
		}
	}

	if (time_step == -1) {
		for (int i=offset; i<offset+newSize; i++){
			ptcl = particle[i];
			if (ptcl->NumberOfAC != 0) {
				if (ptcl->TimeLevelIrr < min_time_level) {
					min_time_level = ptcl->TimeLevelIrr;
				}
			}
		}
		// resetting time_block based on the system
		time_block = std::max(-60, min_time_level-MIN_LEVEL_BUFFER);
		block_max = static_cast<ULL>(pow(2., -time_block));
		time_step = std::pow(2.,time_block);
	}

	for (int i=offset; i<offset+newSize; i++) {
		ptcl = particle[i];
		ptcl->TimeBlockIrr = static_cast<ULL>(pow(2., ptcl->TimeLevelIrr-time_block));
		ptcl->TimeBlockReg = static_cast<ULL>(pow(2., ptcl->TimeLevelReg-time_block));
		if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 < KSTime)
			BinaryCandidateList.push_back(ptcl);
	}

	return true;
}


