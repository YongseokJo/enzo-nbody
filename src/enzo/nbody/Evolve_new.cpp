#include <cmath>
#include <iostream>
#include "defs.h"
#include "global.h"


int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int outputNum);
int ReceiveFromEzno(std::vector<Particle*> &particle);
int SendToEzno(std::vector<Particle*> &particle);
bool CreateComputationChain(std::vector<Particle*> &particle, std::vector<*Particle> &list);

bool IsOutput         = false;
double outputTime     = 0.;
double outputTimeStep = 0.;
double NextRegTime    = 0.;
std::vector<Particle*> &ComputationChain;


void Evolve(std::vector<Particle*> &particle) {

	std::cout << "nbody+: Evolve starting ..." << std::endl;
	int outNum = 0;
	int freq   = 0;

	if (NNB == 0) { 
		std::cout << "nbody+: no particle to be calculated ..." << std::endl;
		goto Communication;
	}

    CreateComputationChain(particle, ComputationChain);

	while (true)
	{
		while (global_time < 1)
		{
			// It's time to compute regular force.
			RegularAccelerationRoutine(particle); // do not update particles unless NNB=0
			IrregularAccelerationRoutine(particle);
		}

	Communication:
		do
		{
			SendToEzno(particle);
			ReceiveFromEzno(particle);
		} while (NNB == 0);
		global_time = 0.;
		NextRegTime = 0.;
	}
}



