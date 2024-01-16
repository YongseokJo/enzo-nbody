#include "global.h"
#include <iostream>
#include <cmath>



double getNewTimeStep(double f[3][4], double df[3][4]);
double getBlockTimeStep(double dt);

/*
 *  Purporse: Initialize timsteps 
 *
 *  Modified: 2024.01.16  by Yongseok Jo
 *
 */
int InitializeTimeStep(std::vector<Particle*> &particle) {
	std::cout << "Initializing timesteps ..." << std::endl;
	double timestep_min=1e30;
	for (Particle* ptcl: particle) {
		dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
		if (dtIrr < timestep_min) 
			timestep_min = dtIrr;
	}
	
	// dt_block_level = log2(dt_enzo/dt_irr) 
	// dt_block       = dt_enzo/2^dt_block_level
	dt_block_level = static_cast<int>(std::ceil(log((EnzoTimeStep/timestep_min)/log(2.0))));
	dt_block       = std::pow(2, -dt_block_level)*EnzoTimeStep;

	for (Particle* ptcl: particle) {
		dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
		dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
		TimeStepIrr = getBlockTimeStep(dtIrr);
		TimeStepReg = getBlockTimeStep(dtReg);
	}
}
