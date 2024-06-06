#include <iostream>
#include <cmath>
#include "../global.h"


double getNewTimeStep(double f[3][4], double df[3][4]);
void getBlockTimeStep(double dt, int& TimeLevel, double &TimeStep);

// Update TimeStepIrr
void Particle::calculateTimeStepIrr(double f[3][4],double df[3][4]) {
	double TimeStepIrrTmp;
	int TimeLevelTmp;

	if (this->NumberOfAC == 0)
		return;

	getBlockTimeStep(getNewTimeStep(a_tot, a_irr), TimeLevelTmp, TimeStepIrrTmp);

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
		//std::cerr << "Timestep is too small" << std::endl;
		//std::cout << "In CalTimeStep, timestep is too small."\
		 	<< " TimeStepIrr=" << TimeStepIrr	<< std::endl;
		TimeStepIrr  = std::max(dt_block*dt_min,             TimeStepIrr);
		TimeLevelIrr = std::max(dt_block_level+dt_level_min, TimeLevelIrr);
	}
	TimeStepIrr  = std::max(dt_min/2,             TimeStepIrr);
	TimeLevelIrr = std::max(dt_level_min-1, TimeLevelIrr);
}

// Update TimeStepReg // need to review
void Particle::calculateTimeStepReg() {
	//fprintf(stdout, "Number of AC=%d\n", NumberOfAC);
	//std::cout << NumberOfAC << std::flush;
	double TimeStepRegTmp, TimeStepRegTmp0;
	int TimeLevelTmp;

#ifdef time_trace
		_time.reg_dt1.markStart();
#endif
	getBlockTimeStep(getNewTimeStep(a_reg, a_reg), TimeLevelTmp, TimeStepRegTmp);

	//fprintf(stderr, "in CalReg, raw time step=%.2eMyr, ", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);
#ifdef time_trace
		_time.reg_dt1.markEnd();
		_time.reg_dt1.getDuration();

		_time.reg_dt2.markStart();
#endif

	//std::cout << "NBODY+: TimeStepRegTmp = " << TimeStepRegTmp << std::endl;

	if (TimeStepRegTmp > 2*TimeStepReg) {
		if (fmod(CurrentTimeReg, 2*TimeStepReg)==0 \
				&& CurrentTimeReg != 0) {
			TimeStepRegTmp = 2*TimeStepReg;
			TimeLevelTmp   = TimeLevelReg + 1;
			while ((TimeStepRegTmp0 > TimeStepRegTmp) \
					&& (fmod(CurrentTimeReg, 2*TimeStepRegTmp) == 0)) {
				TimeStepRegTmp = 2*TimeStepRegTmp;
				TimeLevelTmp   = TimeLevelTmp + 1;
			}
		}
		else {
			TimeStepRegTmp = TimeStepReg;
			TimeLevelTmp   = TimeLevelReg;
		}
	}
	else if (TimeStepRegTmp < TimeStepReg) {
		TimeStepRegTmp = TimeStepReg/2;
		TimeLevelTmp--;
		if (TimeStepRegTmp > TimeStepRegTmp0 ) {
			TimeStepRegTmp = TimeStepReg/4;
			TimeLevelTmp -= 2;
		}
	}
	else {
		TimeStepRegTmp  = TimeStepReg;
		TimeLevelTmp = TimeLevelReg;
	}


	// update needed. regcor_gpu.for:725 (Makino, ApJ, 369)
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

#ifdef time_trace
		_time.reg_dt2.markEnd();
		_time.reg_dt2.getDuration();

		_time.reg_dt3.markStart();
#endif

	//fprintf(stderr, " final time step=%.2eMyr\n", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);

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

	while (TimeStepReg < TimeStepIrr) {
		TimeStepIrr /= 2;
		TimeLevelIrr -= 1;
	}

#ifdef time_trace
		_time.reg_dt3.markEnd();
		_time.reg_dt3.getDuration();
#endif

	//std::cout << "NBODY+: TimeStepReg = " << TimeStepReg << std::endl;
}
