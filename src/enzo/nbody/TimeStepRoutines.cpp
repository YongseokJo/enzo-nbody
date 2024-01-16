#include "defs.h"
#include <cmath>
#include <algorithm>
#include <iostream>


double getNewTimeStep(double f[3][4], double df[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep;

	F2     =   f[0][0]*f[0][0] + f[1][0]*f[1][0]   + f[2][0]*f[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];

	TimeStep = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
	TimeStep = std::sqrt(eta*TimeStep);
	//std::cout<< TimeStep << " ";
	//exit(EXIT_FAILURE); 
	return TimeStep;
}

double getBlockTimeStep(double dt) {
	double normalized_dt = dt/EnzoTimeStep;
	double block_dt;


	block_dt = pow(2, floor(log(normalized_dt)/log(2.0)));

	if (block_dt < pow(2,-64)) {
		std::cerr << "The block time step is smaller than the smallest time step of the system. Program terminated.\n";
		exit(EXIT_FAILURE); 
	}

	return block_dt;
}



