#include <iostream>
#include <cmath>
#include "../global.h"


double getNewTimeStep(double f[3][4], double df[3][4]);
void getBlockTimeStep(double dt, int& TimeLevel, ULL &TimeBlock, double &TimeStep);

// Update TimeStepIrr
void Particle::calculateTimeStepIrr(double f[3][4],double df[3][4]) {
	double TimeStepTmp;
	int TimeLevelTmp;
	ULL TimeBlockTmp;

	if (this->NumberOfAC == 0)
		return;

	getBlockTimeStep(getNewTimeStep(a_tot, a_irr), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);


	//std::cout << "TimeStepIrrTmp=" << TimeStepIrrTmp << std::endl;
	while ((CurrentBlockIrr+TimeBlockTmp > CurrentBlockReg+TimeBlockReg) || (TimeLevelTmp >= TimeLevelReg)) {
		TimeLevelTmp--;
		TimeBlockTmp *= 0.5;
	}

	if (TimeLevelTmp > TimeLevelIrr+1) {
		if (fmod(CurrentBlockIrr, 2*TimeBlockIrr)==0) {
			TimeLevelTmp++;
		}
		else {
			TimeLevelTmp   = TimeLevelIrr;
		}
	}
	else if (TimeLevelTmp < TimeLevelIrr) {
		if (TimeLevelTmp < TimeLevelIrr-1) {
			TimeLevelTmp -= 2;
		}
		else {
			TimeLevelTmp--;
		}
	} else {
		TimeLevelTmp = TimeLevelIrr;
	}

	TimeLevelIrr = TimeLevelTmp;

	if (TimeLevelIrr < time_block) {
		//std::cerr << "Timestep is too small" << std::endl;
		TimeLevelIrr = std::max(time_block, TimeLevelIrr);
	}

	TimeStepIrr = static_cast<double>(pow(2, TimeLevelIrr));
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-time_block));
	//std::cout << "TimeStepIrr=" << TimeStepIrr << std::endl;
}



// Update TimeStepReg // need to review
void Particle::calculateTimeStepReg() {
	//fprintf(stdout, "Number of AC=%d\n", NumberOfAC);
	//std::cout << NumberOfAC << std::flush;
	double TimeStepTmp;
	ULL TimeBlockTmp;
	int TimeLevelTmp, TimeLevelTmp0;

	getBlockTimeStep(getNewTimeStep(a_reg, a_reg), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);

	//fprintf(stderr, "in CalReg, raw time step=%.2eMyr, ", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);


	//std::cout << "NBODY+: TimeStepRegTmp = " << TimeStepRegTmp << std::endl;

	TimeLevelTmp0 = TimeLevelTmp;

	if (TimeLevelTmp > TimeLevelReg+1) {
		if (fmod(CurrentBlockReg, 2*TimeBlockReg)==0 \
				&& CurrentTimeReg != 0) {
			TimeLevelTmp   = TimeLevelReg+1;
			TimeBlockTmp   = TimeBlockReg*2;

			while ((TimeLevelTmp0 > TimeLevelTmp) \
					&& (fmod(CurrentBlockReg, 2*TimeBlockTmp) == 0)) {
				TimeBlockTmp = 2*TimeBlockTmp;
				TimeLevelTmp   = TimeLevelTmp + 1;
			}
		}
		else {
			TimeLevelTmp   = TimeLevelReg;
		}
	}
	else if (TimeLevelTmp < TimeLevelReg) {
		TimeLevelTmp--;
		if (TimeLevelTmp > TimeLevelTmp0)
			TimeLevelTmp--;
	}
	else {
		TimeLevelTmp = TimeLevelReg;
	}


	// update needed. regcor_gpu.for:725 (Makino, ApJ, 369)
	/*
	if (TimeStepRegTmp > 0.1 && TimeStepRegTmp > TimeStepReg) {
		double v2 = 0., a2=0., dt;
		for (int dim=0; dim<Dim; dim++) {
			v2 += (PredVelocity[dim]-NewVelocity[dim])*(PredVelocity[dim]-NewVelocity[dim]);
			a2 += a_reg[dim][0]*a_reg[dim][0];
		}
		dt = TimeStepReg*std::pow((1e-4*TimeStepReg*TimeStepReg*a2/v2),0.1);
		if (dt < TimeStepRegTmp) {
			TimeStepRegTmp = TimeStepReg;
		}	
	}
	*/

	//fprintf(stderr, " final time step=%.2eMyr\n", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);

	TimeLevelReg = std::max(time_block,TimeLevelTmp);
	//TimeLevelReg = std::max(time_block, TimeLevelReg);

	if (this->NumberOfAC == 0) {
		TimeLevelIrr = TimeLevelReg;
	}

	while (TimeLevelReg < TimeLevelIrr) {
		TimeLevelIrr--;
	}

	TimeStepIrr = static_cast<double>(pow(2, TimeLevelIrr));
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-time_block));

	if (CurrentTimeReg+TimeStepReg > 1 && CurrentTimeReg != 1.0) {
		TimeStepReg = 1 - CurrentTimeReg;
		TimeBlockReg = block_max-CurrentBlockReg;
	}

	//std::cout << "NBODY+: TimeStepReg = " << TimeStepReg << std::endl;
}
