#include <cmath>
#include <iostream>
#include "defs.h"
#include "global.h"


void SortComputationChain(std::vector<Particle*> &particle);
void UpdateMinRegTime(std::vector<Particle*> &particle, double* MinRegTime);
void UpdateEvolveParticle(std::vector<Particle*> &particle, std::vector<Particle*> &list, double MinRegTime);
int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int outputNum);
int ReceiveFromEzno(std::vector<Particle*> &particle);
int SendToEzno(std::vector<Particle*> &particle);

bool IsEnzoCommunication = false;
bool IsOutput            = false;
double outputTime        = 0.;
double outputTimeStep    = 0.;


void Evolve(std::vector<Particle*> &particle) {

	std::cout << "Evolve starting ..." << std::endl;
	int outNum=0;
	int freq=0;
	double MinRegTime, EvolvedTime; // time where computing is actaully happening
	std::vector<Particle*> EvolveParticle{};
	std::vector<Particle*> EvolveParticleCopy{};


	std::cout << particle.size() << std::endl;
	UpdateMinRegTime(particle, &MinRegTime);

	// This part can be parallelized.
	while (true) {
		freq %= UpdateFrequency;
		if ((++freq == 0) || (EvolveParticle.size() == 0)) {
			// Find particles to compute
			UpdateEvolveParticle(particle, EvolveParticle, MinRegTime);

			// if no particles to compute, we have to do one of them: 
			// 1. regular force calculation; 2. enzo communication; 3. data dump
			// Even this can be parallelized. 
			if (EvolveParticle.size() == 0) {
				// It's time to compute regular force.
				if (global_time == MinRegTime) {
					std::cout <<  "Regular force calculating...\n" << std::flush;
					for (Particle* ptcl:particle) 
						if (ptcl->isRegular) {
							fprintf(stderr, "Particle ID=%d, Time=%.4e, dtIrr=%.4e, dtReg=%.4e\n",ptcl->getPID(), ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->TimeStepReg);
							std::cerr << std::flush;
							ptcl->calculateRegForce(particle, MinRegTime); // this only does acceleration computation without update.
						}
					UpdateMinRegTime(particle, &MinRegTime);
					if (MinRegTime >= 1) 
						exit(EXIT_FAILURE);
					// after regular force, find particles to evolve again
				}
				if (IsEnzoCommunication) {
					
					SendToEzno(particle);
					ReceiveFromEzno(particle);
				}
				if (IsOutput) {
					writeParticle(particle, MinRegTime, ++outNum);
					outputTime += outputTimeStep;
				}
				UpdateEvolveParticle(particle, EvolveParticle, MinRegTime); 
			}
		}


		// Irregular force calculation
		EvolveParticleCopy.clear();
		EvolveParticleCopy.assign(EvolveParticle.begin(), EvolveParticle.end());
		std::cout << "Irregular force calculating...\n" << std::flush;
		for (Particle* ptcl: EvolveParticleCopy) {
			EvolvedTime = ptcl->CurrentTimeIrr+ptcl->TimeStepIrr;
			fprintf(stderr, "PID=%d, MinRegTime= %e, EvloveTime = %e\n",
					ptcl->getPID(),MinRegTime, EvolvedTime);
			fprintf(stderr, "CurrentTimeIrr = %e, TimeStepIrr = %e, CurrentTimeReg=%e, TimeStepReg=%e\n",
					ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->CurrentTimeReg, ptcl->TimeStepReg);

			ptcl->calculateIrrForce(); // this includes particle position and time evolution.
			ptcl->updateEvolveParticle(EvolveParticle, MinRegTime);
		}
	}
}


void UpdateMinRegTime(std::vector<Particle*> &particle, double* MinRegTime) {

	int isRegular=0;
	double time_tmp, time=1e10;

	for (Particle* ptcl: particle) {
		isRegular += ptcl->isRegular;
	}

	// Find minium regular Time
	if (isRegular == 0) {
		for (Particle* ptcl: particle) {
			time_tmp = ptcl->CurrentTimeReg + ptcl->TimeStepReg;
			if (time > time_tmp)
				time = time_tmp;
		}
	}
	else
		return;

	// Set isRegular to 1
	std::cerr << "Regular: ";
	for (Particle* ptcl: particle) {
		time_tmp = ptcl->CurrentTimeReg + ptcl->TimeStepReg;
		if (time == time_tmp) {
			std::cerr << ptcl->getPID() << ' ';
			ptcl->isRegular = 1;
		}
	}
	std::cerr << '\n' << std::flush;
	*MinRegTime = time;
}


void UpdateEvolveParticle(std::vector<Particle*> &particle, std::vector<Particle*> &list, double MinRegTime) {
	double time, next_time;
	IsEnzoCommunication = true;
	list.clear();
	for (Particle* ptcl: particle) {
		time      = ptcl->CurrentTimeIrr;
		next_time = ptcl->TimeStepIrr;

		if ( (MinRegTime >= next_time) // Regular timestep
				 && ptcl->checkNeighborForEvolution()  // neighbor
				 && (EnzoTimeMark >= next_time) ) { // Enzo timestep
			ptcl->isEvolve = 1;
			list.push_back(ptcl);
		}
		//set global time as the time of a particle that has most advanced
		if (time > global_time)
			global_time = time;

		if (IsEnzoCommunication && EnzoTimeMark > next_time) 
			IsEnzoCommunication = false;
	}
	std::cerr << "EvolveParticle: ";
	std::cerr << list.size();
	std::cerr << '\n' << std::flush;
}

