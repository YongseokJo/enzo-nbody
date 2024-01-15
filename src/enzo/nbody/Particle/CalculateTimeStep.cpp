#include "Particle.h"
#include <cmath>
#include <iostream>


double getNewTimeStep(double f[3][4], double df[3][4]);
double getBlockTimeStep(double dt);

// Update TimeStepIrr
void Particle::calculateTimeStepIrr(double f[3][4],double df[3][4]) {
	double TimeStepIrrTmp;
	TimeStepIrrTmp = getBlockTimeStep(getNewTimeStep(f, df));
	while (CurrentTimeIrr+TimeStepIrrTmp > CurrentTimeReg+TimeStepReg) {
		TimeStepIrrTmp *= 0.5;
		//fprintf(stdout, "CTS, time irr =%e\n", TimeStepIrr);
		//std::cout << std::flush;
	}
	if (TimeStepIrrTmp > 2*TimeStepIrr) {
		if (fmod(CurrentTimeIrr, 2*TimeStepIrr)==0)
			TimeStepIrrTmp = 2*TimeStepIrr;
		else
			TimeStepIrrTmp = TimeStepIrr;
	}
	else if (TimeStepIrrTmp < TimeStepIrr) {
		if (0.5*TimeStepIrr > TimeStepIrrTmp)
			TimeStepIrrTmp = TimeStepIrr/4;
		else
			TimeStepIrrTmp = TimeStepIrr/2;
	}
	else
		TimeStepIrr = TimeStepIrrTmp;
}

// Update TimeStepReg
void Particle::calculateTimeStepReg(double f[3][4], double df[3][4]) {
	fprintf(stdout, "Number of AC=%d\n", NumberOfAC);
	//std::cout << NumberOfAC << std::flush;
	double TimeStepRegTmp;
	TimeStepRegTmp = getBlockTimeStep(getNewTimeStep(f, df));
	if (TimeStepIrr > TimeStepReg)
		;

	if (TimeStepRegTmp > 2*TimeStepReg) {
		if (fmod(CurrentTimeReg, 2*TimeStepReg)==0) {
			TimeStepRegTmp = 2*TimeStepReg;
			while (fmod(CurrentTimeReg, 2*TimeStepRegTmp)==0)
				TimeStepRegTmp *= 2;
		}
		else {
			TimeStepRegTmp = TimeStepReg;
		}
	}
	else if (TimeStepRegTmp < TimeStepReg) {
		if (0.5*TimeStepReg > TimeStepRegTmp)
			TimeStepRegTmp = TimeStepReg/4;
		else
			TimeStepRegTmp = TimeStepReg/2;
	}
	else
		TimeStepReg = TimeStepRegTmp;
	if (NumberOfAC == 0) {
		TimeStepIrr = TimeStepReg;
		CurrentTimeIrr = CurrentTimeReg;
	}

}
