#include <iostream>
#include "global.h"

void UpdateNextRegTime(std::vector<Particle*> &particle);

bool RegularAccelerationRoutine(std::vector<Particle*> &particle)
{
    std::cout << "Calculating regular force ...\n" << std::flush;
    
    // Calulating regular acceleration of the particles
    for (Particle *ptcl : particle)
        if (ptcl->isRegular)
        {
            fprintf(stdout, "Particle ID=%d, Time=%.4e, dtIrr=%.4e, dtReg=%.4e\n", ptcl->getPID(), ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->TimeStepReg);
            std::cerr << std::flush;
            // This only computes accelerations without updating particles.
            ptcl->calculateRegForce(particle);
        }
    // update the next regular time step
    UpdateNextRegTime(particle);
}



void UpdateNextRegTime(std::vector<Particle*> &particle) {

	int isRegular = 0;
	double time_tmp, time = 1e10;

    for (Particle *ptcl : particle)
    {
        // Next regular time step
        time_tmp = ptcl->CurrentTimeReg + ptcl->TimeStepReg;

        // Find the minum regular time step
        if (time > time_tmp)
            time = time_tmp;
    }
	NextRegTime = time;

    // Set isRegular of the particles that will be updated next to 1
	std::cerr << "Regular: ";
	for (Particle* ptcl: particle) {
		time_tmp = ptcl->CurrentTimeReg + ptcl->TimeStepReg;
		if (NextRegTime == time_tmp) {
			std::cerr << ptcl->getPID() << ' ';
			ptcl->isRegular = 1;
		}
	}
	std::cerr << '\n' << std::flush;
}

