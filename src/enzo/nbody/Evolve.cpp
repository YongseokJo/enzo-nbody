#include "Particle/Particle.h"
#include <cmath>


void SortComputationChain(std::vector<Particle*> &particle);
void UpdateMinRegTime(std::vector<Particle*> &particle, double* MinRegTime);
void UpdateEvolveParticle(std::vector<Particle*> &particle, std::vector<Particle*> &list, double MinRegTime);
int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int outputNum);



void Evolve(std::vector<Particle*> &particle) {

	std::cout << "Evolve starting ..." << std::endl;
	int outNum=0;
	Particle* ptcl;
	double MinRegTime, EvolvedTime; // time where computing is actaully happening
	std::vector<Particle*> EvolveParticle{};
	std::vector<Particle*> EvolveParticleCopy{};


	std::cout << particle.size() << std::endl;
	std::cout << LevelList[0]	<< std::endl;
	ptcl = particle[LevelList[0]];
	UpdateMinRegTime(particle, &MinRegTime);

	//int j=0; 
	// This part can be parallelized.
	while (true) {
		if (EvolveParticle.size() == 0)
			UpdateEvolveParticle(particle, EvolveParticle, MinRegTime);

		//Even this can be parallelized.
		if (EvolveParticle.size() == 0) {
			std::cout <<  "Regular force calculating...\n" << std::flush;
			for (Particle* ptcl:particle) {
				if (ptcl->isRegular) {//&& 
															//ptcl->CurrentTimeIrr==ptcl->CurrentTimeReg+ptcl->TimeStepReg) {
					fprintf(stderr, "Particle ID=%d, Time=%.4e, dtIrr=%.4e, dtReg=%.4e\n",ptcl->getPID(), ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->TimeStepReg);
					std::cerr << std::flush;
					ptcl->calculateRegForce(particle, MinRegTime); // this only does acceleration computation.
				}
			}
				writeParticle(particle, MinRegTime, ++outNum);
				UpdateMinRegTime(particle, &MinRegTime);
				UpdateEvolveParticle(particle, EvolveParticle, MinRegTime);
			}
			if (MinRegTime >= 1) 
				exit(EXIT_FAILURE);

			/*
			if (++j==4) {
				exit(EXIT_FAILURE);
			}
			*/

		// check if the particle is subject to calculation this level.
		EvolveParticleCopy.clear();
		EvolveParticleCopy.assign(EvolveParticle.begin(), EvolveParticle.end());

		fprintf(stderr, "Loop starting... size=%lu\n", EvolveParticleCopy.size());
		for (Particle* ptcl: EvolveParticleCopy) {
			EvolvedTime = ptcl->CurrentTimeIrr+ptcl->TimeStepIrr;
			fprintf(stderr, "PID=%d, MinRegTime= %e, EvloveTime = %e\n",
					ptcl->getPID(),MinRegTime, EvolvedTime);
			//if (ComputingTime >= EvolvedTime) { // Condition should be added more.
			//if (ptcl->checkNeighborForEvolution())
			//break;
			//if (i == NNB-1) {
			//std::cout << "Warning::Nothing to evolve!!!\n" <<std::flush;
			//std::cout << "Error::Nothing to evolve!!!\n" <<std::flush;
			//exit(EXIT_FAILURE);
			//}

			fprintf(stderr, "Particle ID=%d\n",ptcl->getPID());
			fprintf(stderr, "CurrentTimeIrr = %e, TimeStepIrr = %e, CurrentTimeReg=%e, TimeStepReg=%e\n",
					ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->CurrentTimeReg, ptcl->TimeStepReg);


			std::cout << "Irregular force calculating...\n" << std::flush;
			ptcl->calculateIrrForce(); // this includes particle and time advances.
			//ptcl = ptcl->NextParticle;
			std::cout << "Updating EvolveParticle...\n" << std::flush;
			ptcl->updateEvolveParticle(EvolveParticle, MinRegTime);

			for (Particle* elem: EvolveParticle) 
					std::cerr << elem->getPID() << ' ';
			std::cerr << '\n';
			std::cerr << ptcl->CurrentTimeIrr << '\n';
			std::cerr << ptcl->TimeStepIrr << '\n';
			std::cerr << '\n' << std::flush;
		}
		//SortComputationChain(particle);
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
	list.clear();
	for (Particle* ptcl: particle) {
		if (MinRegTime >= (ptcl->CurrentTimeIrr+ptcl->TimeStepIrr - mytolerance)
				&& (ptcl->checkNeighborForEvolution())) {
			ptcl->isEvolve = 1;
			list.push_back(ptcl);
			std::cerr << ptcl->getPID() << ' ';
		}
	}
	std::cerr << "EvolveParticle: ";
	std::cerr << list.size();
	std::cerr << '\n' << std::flush;
}

