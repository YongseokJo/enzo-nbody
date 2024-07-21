#include <cmath>
#include <iostream>
#include "defs.h"
#include "global.h"


int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int outputNum);
int ReceiveFromEnzo(std::vector<Particle*> &particle);
int SendToEnzo(std::vector<Particle*> &particle);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool RegularAccelerationRoutine(std::vector<Particle*> &particle);
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle);
void UpdateNextRegTime(std::vector<Particle*> &particle);



bool IsOutput           = false;
double binary_time      = 0;
double binary_time_prev = 0;
ULL binary_block        = 0;
double outputTime       = 0;
double outputTimeStep   = 0.;
int outNum              = 0;
double global_time_irr  = 0;
std::vector<Particle*> ComputationChain{};
#ifdef time_trace
TimeTracer _time;
#endif



void Evolve(std::vector<Particle*> &particle) {

	fprintf(nbpout, "Evolve Starts ...\n");

	if (NNB < 2) {
		fprintf(nbpout, "No particle to be calculated ...\n");
		fflush(nbpout);
		goto Communication;
	}

	//CreateComputationChain(particle);

	while (true)
	{
		//writeParticle(particle, EnzoCurrentTime, outNum++);
		while (global_time < 1)
		{
			// It's time to compute regular force.
#ifdef time_trace
		_time.irr.markStart();
#endif
			IrregularAccelerationRoutine(particle);
#ifdef time_trace
		_time.irr.markEnd();
		_time.reg.markStart();
#endif
			RegularAccelerationRoutine(particle); // do not update particles unless NNB=0
#ifdef time_trace
		_time.reg.markEnd();
		_time.reg.getDuration();
		_time.irr.getDuration();
		_time.output();
#endif
			/*
			std::cout << "CurrentTimeReg  =" << particle[0]->CurrentTimeReg << std::endl;
			std::cout << "CurrentTimeIrr  =" << particle[0]->CurrentTimeIrr << std::endl;
			std::cout << "TimeStepReg     =" << particle[0]->TimeStepReg << std::endl;
			std::cout << "NextRegTimeBlock=" << NextRegTimeBlock*time_step << std::endl;
			*/

		}

	Communication:
		do
		{
			// in case of nnb=1, only analytic solution be needed.
			fprintf(nbpout, "global time=%lf\n", global_time);
			SendToEnzo(particle);
			ReceiveFromEnzo(particle);
			UpdateNextRegTime(particle);
			FirstComputation = nullptr;
			ComputationList.clear();
			BinaryCandidateList.clear();
			BinaryList.clear();
			fprintf(nbpout, "RegularList size=%d\n", RegularList.size());
			fflush(nbpout);
		} while (NNB < 2);


		global_time      = 0.;
		global_time_irr  = 0.;
		//NextRegTimeBlock = 0.;
	}
}



