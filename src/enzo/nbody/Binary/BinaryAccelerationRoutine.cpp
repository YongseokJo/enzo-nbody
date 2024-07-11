#include <iostream>
#include "global.h"
#include "defs.h"

void NewKSInitialization(Particle* ptclI, std::vector<Particle*> &particle, std::vector<Particle*> &ComputationList);
void KSTermination(Particle* ptclCM, std::vector<Particle*> &particle, double current_time, ULL current_block);

bool AddNewBinariesToList(std::vector<Particle*> &ComputationList, std::vector<Particle*> &particle) {

	fprintf(stdout, "Finding new binaries ...\n");
	// add new binaries
	for (Particle *ptcl : ComputationList) {
		// if the irregular time step is too short, check if it is binary
		if ((ptcl->TimeStepIrr*EnzoTimeStep*1e4<KSTime) && ( (ptcl->isBinary == false) && (ptcl->isCMptcl == false) )) {
			ptcl->isKSCandidate();
			if (ptcl->isBinary) {
				std::cout << "AddNewBinaries ... new binary pair found" << std::endl;
				fprintf(binout, "BinaryAccelerationRoutine.cpp: new binary particle found!\n");
				// the initialization of the binary counterpart will be done together in the following function.
				NewKSInitialization(ptcl,particle,ComputationList);
				fprintf(stdout, "New binary of (%d, %d) initialization finished ...\n",ptcl->PID, ptcl->BinaryPairParticle->PID);
				fprintf(binout,"\n After binary addition, the number of particles are... %d \n",int(particle.size()));
			}
		}
	}
	//fprintf(binout,"\n After binary addition, the number of particles are... %d \n",int(particle.size()));
	return true;
}

void BinaryAccelerationRoutine(double next_time, std::vector<Particle*> &particle) {

	int count;
	int bincount = 0;

	count = 0;

	if (next_time == 0) {
		return;
	}

	for (Binary* ptclBin: BinaryList) {

		ptclBin->KSIntegration(next_time, bincount);

		count += 1;

		fprintf(binout, "\nBinaryAccelerationRoutine.cpp: After KS Integration of %dth binary....\n", count);
		fprintf(binout, "The ID of ith particle is %d \n",ptclBin->ptclCM->BinaryParticleI->PID);
		fprintf(binout, "The ID of jth particle is %d \n",ptclBin->ptclCM->BinaryParticleJ->PID);
		//fflush(binout);	

		//if (bincount>0) {
		//	std::cout << "Integrating Binary ..." << std::endl;

		//	fprintf(binout, "KS coordinates - u1:%e, u2:%e, u3:%e, u4:%e\n", ptclBin->u[0], ptclBin->u[1], ptclBin->u[2], ptclBin->u[3]);
		//	fprintf(binout, "KS coordinates - udot1:%e, udot2:%e, udot3:%e, udot4:%e\n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
		//	fprintf(binout, "KS coordinates - udot1:%e, udot2:%e, udot3:%e, udot4:%e\n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
		//	fprintf(binout, "Other important KS variables - r:%e, h:%e, gamma: %e, tau:%e, step:%e, currentTime: %e \n", ptclBin->r, ptclBin->h, ptclBin->gamma, ptclBin->dTau, ptclBin->TimeStep, ptclBin->CurrentTime);
		//	fprintf(binout, "loop number = %d \n", bincount);
		//	fflush(binout);
		//}

	}
}
