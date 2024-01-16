#include "global.h"
#include <iostream>
#include <cmath>


double dt_block;
int dt_block_level;

double getNewTimeStep(double f[3][4], double df[3][4]);
double getBlockTimeStep(double dt, int &TimeLevel, double &TimeStep);

/*
 *  Purporse: Initialize timsteps 
 *
 *  Modified: 2024.01.16  by Yongseok Jo
 *
 */
int InitializeTimeStep(std::vector<Particle*> &particle) {
	std::cout << "Initializing timesteps ..." << std::endl;
	double timestep_min=1e30;
	double dtIrr, dtReg;
	for (Particle* ptcl: particle) {
		dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
		dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
		getBlockTimeStep(dtIrr, ptcl->TimeLevelIrr, ptcl->TimeStepIrr);
		getBlockTimeStep(dtReg, ptcl->TimeLevelReg, ptcl->TimeStepReg);

		ptcl->CurrentTimeIrr = 0;
		ptcl->CurrentTimeReg = 0;


		if (ptcl->TimeStepIrr < timestep_min) { 
			dt_block       = ptcl->TimeStepIrr;
			dt_block_level = ptcl->TimeLevelIrr;
		}
	}

	for (Particle* ptcl: particle) {
		ptcl->TimeStepReg = std::min(1.,ptcl->TimeStepReg);
		ptcl->TimeLevelReg = std::min(0,ptcl->TimeLevelReg);

		while (ptcl->TimeStepIrr >= ptcl->TimeStepReg) {
			ptcl->TimeStepIrr *= 0.5; 
			ptcl->TimeLevelIrr--; 
		}

		if (ptcl->TimeLevelIrr < dt_block_level + dt_level_min ) {
			std::cerr << "Timestep is too small" << std::endl;
			ptcl->TimeStepIrr  = std::max(dt_block*dt_min, ptcl->TimeStepIrr);
			ptcl->TimeLevelIrr = std::max(dt_block_level+dt_level_min, ptcl->TimeLevelIrr);
		}
	}
}


int InitializeTimeStep(Particle* particle, int size) {
	std::cout << "Initializing timesteps ..." << std::endl;
	double timestep_min=1e30;
	double dtIrr, dtReg;
	Particle *ptcl;

	for (int i=0; i<size; i++){
		ptcl = &particle[i];
		dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
		dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
		getBlockTimeStep(dtIrr, ptcl->TimeLevelIrr, ptcl->TimeStepIrr);
		getBlockTimeStep(dtReg, ptcl->TimeLevelReg, ptcl->TimeStepReg);

		ptcl->CurrentTimeIrr = 0;
		ptcl->CurrentTimeReg = 0;


		if (ptcl->TimeStepIrr < timestep_min) { 
			dt_block       = ptcl->TimeStepIrr;
			dt_block_level = ptcl->TimeLevelIrr;
		}
	}

	for (int i=0; i<size; i++){
		ptcl->TimeStepReg = std::min(1.,ptcl->TimeStepReg);
		ptcl->TimeLevelReg = std::min(0,ptcl->TimeLevelReg);

		while (ptcl->TimeStepIrr >= ptcl->TimeStepReg) {
			ptcl->TimeStepIrr *= 0.5; 
			ptcl->TimeLevelIrr--; 
		}

		if (ptcl->TimeLevelIrr < dt_block_level + dt_level_min ) {
			std::cerr << "Timestep is too small" << std::endl;
			ptcl->TimeStepIrr  = std::max(dt_block*dt_min, ptcl->TimeStepIrr);
			ptcl->TimeLevelIrr = std::max(dt_block_level+dt_level_min, ptcl->TimeLevelIrr);
		}
	}
}
