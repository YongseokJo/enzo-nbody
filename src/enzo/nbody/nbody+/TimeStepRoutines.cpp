#include "global.h"
#include "defs.h"
#include <cmath>
#include <algorithm>
#include <iostream>


double getNewTimeStep(double f[3][4], double df[3][4], double dt) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     =   f[0][0]*f[0][0] +  f[1][0]*f[1][0]  + f[2][0]*f[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];


	/*
	else if  (F2 == 0 && F2dot2 != 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
	else if (F2 != 0 && F2dot2 == 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
													 */
	//fprintf(stdout, "NBODY+: f dots: %e, %e, %e, %e\n", F2, Fdot2, F2dot2, F3dot2);
	/*
	if (F2 != 0 && F2dot2 != 0 && dt != 0) {
		DivergentPrevent = std::sqrt(F2/Fdot2)/dt;
		fprintf(stdout, "F2 = %e, F2dot2 = %e, dt = %e\n", F2, F2dot2, dt);
		fprintf(stdout, "DivergentPrevent = %e\n", DivergentPrevent);
	}
	else {
		DivergentPrevent = 1;
	}
	*/

  TimeStep  = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
  //TimeStep  = F2/Fdot2;
	//TimeStep  = std::sqrt(DivergentPrevent*eta*TimeStep);
	TimeStep  = std::sqrt(eta*TimeStep);
	//std::cout<< TimeStep << " ";
	//exit(EXIT_FAILURE); 
	return TimeStep;
}

void getBlockTimeStep(double dt, int &TimeLevel, double &TimeStep) {
	TimeLevel = static_cast<int>(floor(log(dt/EnzoTimeStep)/log(2.0)));
	//std::cout << "NBODY+: TimeLevel = " << TimeLevel << std::endl;
	TimeStep = static_cast<double>(pow(2, TimeLevel));
	//std::cout << "NBODY+: TimeStep = " << TimeStep << std::endl;
}


