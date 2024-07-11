#include "global.h"
#include "defs.h"
#include <cmath>
#include <algorithm>
#include <iostream>




double getNewTimeStepReg(double v[3], double df[3][4]) {

	double v2, F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     = df[0][0]*df[0][0] + df[1][0]*df[1][0] + df[2][0]*df[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];
	v2     = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];


  TimeStep  = (std::sqrt(v2*Fdot2)+F2)/(std::sqrt(F2*F2dot2)+Fdot2);
	TimeStep  = std::sqrt(eta*TimeStep);

	return TimeStep;
}


double getNewTimeStepIrr(double f[3][4], double df[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     =   f[0][0]*f[0][0] +  f[1][0]*f[1][0]  + f[2][0]*f[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];



  TimeStep  = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
	TimeStep  = std::sqrt(eta*TimeStep);
	return TimeStep;
}

double getNewTimeStep(double f[3][4], double df[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     =   f[0][0]*f[0][0] +  f[1][0]*f[1][0]  + f[2][0]*f[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];


  TimeStep  = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
	TimeStep  = std::sqrt(eta*TimeStep);
	return TimeStep;
}

void getBlockTimeStep(double dt, int &TimeLevel, ULL &TimeBlock, double &TimeStep) {
	TimeLevel = static_cast<int>(floor(log(dt/EnzoTimeStep)/log(2.0)));
	//TimeLevel = static_cast<int>(ceil(log(dt/EnzoTimeStep)/log(2.0)));
	//std::cout << "NBODY+: TimeLevel = " << TimeLevel << std::endl;
	//std::cout << "NBODY+: TimeStep = " << TimeStep << std::endl;
	
	if (TimeLevel < time_block) {
		//std::cerr << "TimeLevel is less than time block!!" << std::endl;
		TimeLevel = time_block;
	}

	TimeStep = static_cast<double>(pow(2, TimeLevel));
	TimeBlock = static_cast<ULL>(pow(2, TimeLevel-time_block));
}


