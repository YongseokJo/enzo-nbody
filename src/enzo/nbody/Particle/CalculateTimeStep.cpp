#include <iostream>
#include <cmath>
#include "../global.h"


double getNewTimeStep(double f[3][4], double df[3][4]);
double getBlockTimeStep(double dt, int& TimeLevel, double &TimeStep);

// Update TimeStepIrr
void Particle::calculateTimeStepIrr(double f[3][4],double df[3][4]) {
	double TimeStepIrrTmp;
	int TimeLevelTmp;

	getBlockTimeStep(getNewTimeStep(f, df), TimeLevelTmp, TimeStepIrrTmp);

	while (CurrentTimeIrr+TimeStepIrrTmp > CurrentTimeReg+TimeStepReg) {
		TimeStepIrrTmp *= 0.5;
		TimeLevelTmp--;
	}

	if (TimeStepIrrTmp > 2*TimeStepIrr) {
		if (fmod(CurrentTimeIrr, 2*TimeStepIrr)==0) {
			TimeStepIrrTmp = 2*TimeStepIrr;
			TimeLevelTmp++;
		}
		else {
			TimeStepIrrTmp = TimeStepIrr;
			TimeLevelTmp   = TimeLevelIrr;
		}
	}
	else if (TimeStepIrrTmp < TimeStepIrr) {
		if (0.5*TimeStepIrr > TimeStepIrrTmp) {
			TimeStepIrrTmp = TimeStepIrr/4;
			TimeLevelTmp -= 2;
		}
		else {
			TimeStepIrrTmp = TimeStepIrr/2;
			TimeLevelTmp--;
		}
	}
	else {
		TimeStepIrr = TimeStepIrrTmp;
		TimeLevelIrr = TimeLevelTmp;
	}

	if (TimeLevelIrr < dt_block_level+dt_level_min) {
		std::cerr << "Timestep is too small" << std::endl;
		TimeStepIrr  = std::max(dt_block*dt_min,             TimeStepIrr);
		TimeLevelIrr = std::max(dt_block_level+dt_level_min, TimeLevelIrr);
	}
}

// Update TimeStepReg
void Particle::calculateTimeStepReg(double f[3][4], double df[3][4]) {
	fprintf(stdout, "Number of AC=%d\n", NumberOfAC);
	//std::cout << NumberOfAC << std::flush;
	double TimeStepRegTmp;
	int TimeLevelTmp;
	getBlockTimeStep(getNewTimeStep(f, df), TimeLevelTmp, TimeStepRegTmp);


	if (TimeStepRegTmp > 2*TimeStepReg) {
		if (fmod(CurrentTimeReg, 2*TimeStepReg)==0) {
			TimeStepRegTmp = 2*TimeStepReg;
			TimeLevelTmp++;
			while (fmod(CurrentTimeReg, 2*TimeStepRegTmp)==0) {
				TimeStepRegTmp *= 2;
				TimeLevelTmp++;
			}
		}
		else {
			TimeStepRegTmp = TimeStepReg;
			TimeLevelTmp   = TimeLevelReg;
		}
	}
	else if (TimeStepRegTmp < TimeStepReg) {
		if (0.5*TimeStepReg > TimeStepRegTmp) {
			TimeStepRegTmp = TimeStepReg/4;
			TimeLevelTmp -= 2;
		}
		else {
			TimeStepRegTmp = TimeStepReg/2;
			TimeLevelTmp--;
		}
	}
	else {
		TimeStepReg  = TimeStepRegTmp;
		TimeLevelReg = TimeLevelTmp;
	}

	if (NumberOfAC == 0) {
		TimeStepIrr    = TimeStepReg;
		TimeLevelIrr   = TimeLevelReg;
		CurrentTimeIrr = CurrentTimeReg;
	}

	TimeStepReg  = std::min(1.,TimeStepReg);
	TimeLevelReg = std::min(0,TimeLevelReg);

	while (TimeStepIrr >= TimeStepReg) {
		TimeStepIrr *= 0.5;
		TimeLevelIrr--;
	}
}
