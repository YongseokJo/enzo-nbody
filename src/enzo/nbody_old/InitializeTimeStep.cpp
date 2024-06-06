#include "global.h"
#include <iostream>
#include <ostream>
#include <cmath>


double dt_block;
int dt_block_level;

double getNewTimeStep(double f[3][4], double df[3][4]);
void getBlockTimeStep(double dt, int &TimeLevel, double &TimeStep);

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
		dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
		//std::cout << "dtReg=" << dtReg << std::endl;
		getBlockTimeStep(dtReg, ptcl->TimeLevelReg, ptcl->TimeStepReg);

		if (ptcl->NumberOfAC != 0) {
			dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
			getBlockTimeStep(dtIrr, ptcl->TimeLevelIrr, ptcl->TimeStepIrr);
		}
		else {
			ptcl->TimeLevelIrr = ptcl->TimeLevelReg;
			ptcl->TimeStepIrr  = ptcl->TimeStepReg;
		}

		ptcl->TimeStepReg = std::min(1.,ptcl->TimeStepReg);
		ptcl->TimeLevelReg = std::min(0,ptcl->TimeLevelReg);

		ptcl->CurrentTimeIrr = 0;
		ptcl->CurrentTimeReg = 0;

#define no_IRR_TEST
#ifdef IRR_TEST
		ptcl->TimeStepReg = 1;
#endif
	}



	// Irregular Time Step Correction
	for (Particle* ptcl: particle) {
		if (ptcl->NumberOfAC != 0) {
			while (ptcl->TimeStepIrr >= ptcl->TimeStepReg) {
				ptcl->TimeStepIrr *= 0.5; 
				ptcl->TimeLevelIrr--; 
			}
		}
		if (ptcl->TimeStepIrr < timestep_min) { 
			dt_block       = ptcl->TimeStepIrr;
			dt_block_level = ptcl->TimeLevelIrr;
		}
	}

	for (Particle* ptcl: particle) {
		if (ptcl->TimeLevelIrr < dt_block_level + dt_level_min ) {
			std::cerr << "Timestep is too small" << std::endl;
			std::cout << "In Initialization, timestep is too small." \
				<< " TimeStepIrr=" << ptcl->TimeStepIrr	<< std::endl;
			ptcl->TimeStepIrr  = std::max(dt_block*dt_min, ptcl->TimeStepIrr);
			ptcl->TimeLevelIrr = std::max(dt_block_level+dt_level_min, ptcl->TimeLevelIrr);
		}
	}

	fprintf(stdout, "nbody+: dtIrr = %e, dtReg= %e, TimeStepIrr=%e, TimeStepReg=%e\n",
			dtIrr, dtReg, particle[0]->TimeStepIrr, particle[0]->TimeStepReg);
	fprintf(stdout, "nbody+:dt_block = %e, dt_block_level= %d, EnzoTimeStep=%e\n",
			dt_block, dt_level_min, EnzoTimeStep);
	return true;
}


int InitializeTimeStep(std::vector<Particle*> &particle, int offset) {
	std::cout << "Initializing timesteps ..." << std::endl;
	double timestep_min=1e30;
	double dtIrr, dtReg;
	Particle *ptcl;

	for (int i=offset; i<offset+newNNB; i++){
		ptcl = particle[i];
		dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
		getBlockTimeStep(dtReg, ptcl->TimeLevelReg, ptcl->TimeStepReg);

		if (ptcl->NumberOfAC != 0) {
			dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
			getBlockTimeStep(dtIrr, ptcl->TimeLevelIrr, ptcl->TimeStepIrr);
		}
		else {
			ptcl->TimeLevelIrr = ptcl->TimeLevelReg;
			ptcl->TimeStepIrr  = ptcl->TimeStepReg;
		}

		ptcl->TimeStepReg = std::min(1.,ptcl->TimeStepReg);
		ptcl->TimeLevelReg = std::min(0,ptcl->TimeLevelReg);

		ptcl->CurrentTimeIrr = 0;
		ptcl->CurrentTimeReg = 0;


		if (ptcl->TimeStepIrr < timestep_min) { 
			dt_block       = ptcl->TimeStepIrr;
			dt_block_level = ptcl->TimeLevelIrr;
		}
	} // endfor size

	for (int i=offset; i<offset+NNB; i++){
		ptcl = particle[i];

		if (ptcl->NumberOfAC != 0) {
			if (ptcl->TimeStepIrr >= ptcl->TimeStepReg) {
				ptcl->TimeStepIrr = 0.5*ptcl->TimeStepReg; 
				ptcl->TimeLevelIrr = ptcl->TimeLevelReg-1; 
			}
		}

		if (ptcl->TimeLevelIrr < dt_block_level + dt_level_min ) {
			std::cerr << "Timestep is too small" << std::endl;
			std::cout << "In Initialization, timestep is too small." \
				<< " TimeStepIrr=" << ptcl->TimeStepIrr	<< std::endl;
			ptcl->TimeStepIrr  = std::max(dt_block*dt_min, ptcl->TimeStepIrr);
			ptcl->TimeLevelIrr = std::max(dt_block_level+dt_level_min, ptcl->TimeLevelIrr);
		}
	} //endfor size
	return 1;
}



