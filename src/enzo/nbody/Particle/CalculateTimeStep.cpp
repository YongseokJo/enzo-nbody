#include <iostream>
#include <cmath>
#include "../global.h"


double getNewTimeStep(double f[3][4], double df[3][4], double dt);
void getBlockTimeStep(double dt, int& TimeLevel, double &TimeStep);

// Update TimeStepIrr
void Particle::calculateTimeStepIrr(double f[3][4],double df[3][4]) {
	double TimeStepIrrTmp;
	int TimeLevelTmp;

	if (this->NumberOfAC == 0)
		return;

	getBlockTimeStep(getNewTimeStep(f, df, TimeStepIrr), TimeLevelTmp, TimeStepIrrTmp);

	while ((CurrentTimeIrr+TimeStepIrrTmp > CurrentTimeReg+TimeStepReg) || (TimeStepIrrTmp >= TimeStepReg)) {
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
		if (TimeStepIrrTmp < 0.5*TimeStepIrr) {
			TimeStepIrrTmp = TimeStepIrr/4;
			TimeLevelTmp -= 2;
		}
		else {
			TimeStepIrrTmp = TimeStepIrr/2;
			TimeLevelTmp--;
		}
	} else {
		TimeStepIrrTmp = TimeStepIrr;
		TimeLevelTmp = TimeLevelIrr;
	}

	TimeStepIrr = TimeStepIrrTmp;
	TimeLevelIrr = TimeLevelTmp;

	if (TimeLevelIrr < dt_block_level+dt_level_min) {
		std::cerr << "Timestep is too small" << std::endl;
		std::cout << "In CalTimeStep, timestep is too small."\
		 	<< " TimeStepIrr=" << TimeStepIrr	<< std::endl;
		TimeStepIrr  = std::max(dt_block*dt_min,             TimeStepIrr);
		TimeLevelIrr = std::max(dt_block_level+dt_level_min, TimeLevelIrr);
	}
	TimeStepIrr  = std::max(dt_min/2,             TimeStepIrr);
	TimeLevelIrr = std::max(dt_level_min-1, TimeLevelIrr);
}

// Update TimeStepReg // need to review
void Particle::calculateTimeStepReg(double f[3][4], double df[3][4]) {
	//fprintf(stdout, "Number of AC=%d\n", NumberOfAC);
	//std::cout << NumberOfAC << std::flush;
	double TimeStepRegTmp;
	int TimeLevelTmp;

#ifdef time_trace
		_time.reg_dt1.markStart();
#endif
	getBlockTimeStep(getNewTimeStep(f, df, TimeStepReg), TimeLevelTmp, TimeStepRegTmp);

#ifdef time_trace
		_time.reg_dt1.markEnd();
		_time.reg_dt1.getDuration();

		_time.reg_dt2.markStart();
#endif

	//std::cout << "NBODY+: TimeStepRegTmp = " << TimeStepRegTmp << std::endl;

	if (TimeStepRegTmp > 2*TimeStepReg) {
		if (fmod(CurrentTimeReg, 2*TimeStepReg)==0 && CurrentTimeReg != 0) {
			TimeStepRegTmp = 2*TimeStepReg;
			TimeLevelTmp   = TimeLevelReg + 1;
		}
		else {
			TimeStepRegTmp = TimeStepReg;
			TimeLevelTmp   = TimeLevelReg;
		}
	}
	else if (TimeStepRegTmp < TimeStepReg) {
		if (TimeStepRegTmp < 0.5*TimeStepReg ) {
			TimeStepRegTmp = TimeStepReg/4;
			TimeLevelTmp -= 2;
		}
		else {
			TimeStepRegTmp = TimeStepReg/2;
			TimeLevelTmp--;
		}
	}
	else {
		TimeStepRegTmp  = TimeStepReg;
		TimeLevelTmp = TimeLevelReg;
	}


#ifdef time_trace
		_time.reg_dt2.markEnd();
		_time.reg_dt2.getDuration();

		_time.reg_dt3.markStart();
#endif

	TimeStepReg  = std::min(1.,TimeStepRegTmp);
	TimeLevelReg = std::min(0,TimeLevelTmp);

	if (CurrentTimeReg+TimeStepReg > 1 && CurrentTimeReg != 1.0) {
		TimeStepReg = 1 - CurrentTimeReg;
	}


	if (this->NumberOfAC == 0) {
		TimeStepIrr  = TimeStepReg;
		TimeLevelIrr = TimeLevelReg;
	}

	TimeStepReg  = std::max(dt_min,       TimeStepReg);
	TimeLevelReg = std::max(dt_level_min, TimeLevelReg);

#ifdef time_trace
		_time.reg_dt3.markEnd();
		_time.reg_dt3.getDuration();
#endif

	//std::cout << "NBODY+: TimeStepReg = " << TimeStepReg << std::endl;
}
