#include "stdio.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"

bool CreateComputationList(Particle* ptcl);
bool CreateComputationChain(std::vector<Particle*> &particle);
void generate_Matrix(double a[3], double (&A)[3][4]);
void ReInitializeKSParticle(Particle* KSParticle, std::vector<Particle*> &particle);
void UpdateNextRegTime(std::vector<Particle*> &particle);

void KSTermination(Particle* ptclCM, std::vector<Particle*> &particle, double current_time, ULL current_block){

	double R[Dim], Rdot[Dim];
	double Rinv;
	double ratioM;
	double L[3][4];

	int ptclCMIndex;
	int ptclBinIndex;

	bool findPtclCM;

	Particle* ptclI;
	Particle* ptclJ;

	Binary* ptclBin;

	fprintf(stdout,"--------------------------------------\n");
	fprintf(stdout,"In KSRegularlizationTermination.cpp...\n\n");

	ptclI = ptclCM->BinaryParticleI;
	ptclJ = ptclCM->BinaryParticleJ;
	ptclBin = ptclCM->BinaryInfo;


	// convert partI and J's coordinates back to cartesian form KS coordinates
	ptclCM->convertBinaryCoordinatesToCartesian();


	// Initialize Neighbor list
	ptclI->ACList.clear();
	ptclI->NumberOfAC = 0;

	ptclJ->ACList.clear();
	ptclJ->NumberOfAC = 0;

	fprintf(stdout,"initialize particle I \n");
	ReInitializeKSParticle(ptclI, particle);
	fprintf(stdout,"initialize particle J \n");
	ReInitializeKSParticle(ptclJ, particle);



	//	InitializeTimeStep
	ptclI->calculateTimeStepReg();
	if (ptclI->TimeLevelReg <= ptclCM->TimeLevelReg-1 
			&& ptclI->TimeBlockReg/2+ptclI->CurrentBlockReg > ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
		ptclI->TimeLevelReg = ptclCM->TimeLevelReg-1;
	}
	else if  (ptclI->TimeLevelReg >= ptclCM->TimeLevelReg+1) {
		ptclI->TimeLevelReg = ptclCM->TimeLevelReg+1;
	}
	else 
		ptclI->TimeLevelReg = ptclCM->TimeLevelReg;

	ptclI->TimeStepReg  = static_cast<double>(pow(2, ptclI->TimeLevelReg));
	ptclI->TimeBlockReg = static_cast<ULL>(pow(2, ptclI->TimeLevelReg-time_block));
	ptclI->calculateTimeStepIrr(ptclI->a_tot, ptclI->a_irr);

	ptclJ->calculateTimeStepReg();
	if (ptclJ->TimeLevelReg <= ptclCM->TimeLevelReg-1 
			&& ptclJ->TimeBlockReg/2+ptclJ->CurrentBlockReg > ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
		ptclJ->TimeLevelReg = ptclCM->TimeLevelReg-1;
	}
	else if  (ptclJ->TimeLevelReg >= ptclCM->TimeLevelReg+1) {
		ptclJ->TimeLevelReg = ptclCM->TimeLevelReg+1;
	}
	else 
		ptclJ->TimeLevelReg = ptclCM->TimeLevelReg;
	ptclJ->TimeStepReg  = static_cast<double>(pow(2, ptclJ->TimeLevelReg));
	ptclJ->TimeBlockReg = static_cast<ULL>(pow(2, ptclJ->TimeLevelReg-time_block));
	ptclJ->calculateTimeStepIrr(ptclJ->a_tot, ptclJ->a_irr);



	// we also need to revert the neighbor list of Particles
	// assuming that all the neighbors are bidirectional
	// may need to update later if the radius for neighbor differs depending on the particle
	fprintf(stdout,"replacing CM particle in neighbor list to component particles \n");
	ptclCM->isErase = true;
	for (Particle* ptcl: particle) {
		auto it = ptcl->ACList.erase(
				std::remove_if(ptcl->ACList.begin(), ptcl->ACList.end(),
					[](Particle* p) {
					bool to_remove = p->isErase;
					return to_remove;
					}),
				ptcl->ACList.end());

		if (it != ptcl->ACList.end())  {
			ptcl->ACList.push_back(ptclI);
			ptcl->ACList.push_back(ptclJ);
		}
	}


	fprintf(stdout,"add the binary components to particle list\n");
	particle.push_back(ptclI);
	particle.push_back(ptclJ);


	// delete the original components from the list
	fprintf(stdout,"deleting CM particle from the particle list\n");


	fprintf(stderr,"PID of (CM, I, J) = (%d,%d,%d)\n",ptclCM->PID, ptclI->PID, ptclJ->PID);
	/*
	ComputationList.erase(
			std::remove_if(ComputationList.begin(), ComputationList.end(),
				[](Particle* p) {
				bool to_remove = p->isErase;
				//if (to_remove) delete p;
				return to_remove;
				}),
			ComputationList.end());
			*/
	// delete ptclCM from ComputationChain
	Particle* NextParticle=FirstComputation;
	while (NextParticle != nullptr) {
		if (NextParticle->NextParticleForComputation == ptclCM) {
			NextParticle->NextParticleForComputation = ptclCM->NextParticleForComputation;
		}
		NextParticle = NextParticle->NextParticleForComputation;
	}
	particle.erase(
			std::remove_if(particle.begin(), particle.end(),
				[](Particle* p) {
				bool to_remove = p->isErase;
				if (to_remove) delete p;
				return to_remove;
				}),
			particle.end());

	// we also need to delete it from the binary list
	fprintf(stdout,"deleting binary information from the BinaryList \n");
	ptclBin->isErase = true;
	BinaryList.erase(
			std::remove_if(BinaryList.begin(), BinaryList.end(),
				[](Binary* p) {
				bool to_remove = p->isErase;
				if (to_remove) delete p;
				return to_remove;
				}),
			BinaryList.end());


	//re-do UpdateNextRegTime
	UpdateNextRegTime(particle); //

	fprintf(stderr, "particle:");
	for (Particle* ptcl:particle) {
		fprintf(stderr, "%d, ", ptcl->PID);
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "ComputationList:");
	for (Particle* ptcl:ComputationList) {
		fprintf(stderr, "%d, ", ptcl->PID);
	}
	fprintf(stderr, "\n");

	// change the booleans and pointers for binary
	ptclI->isBinary = false;
	ptclI->BinaryPairParticle = nullptr;

	ptclJ->isBinary = false;
	ptclJ->BinaryPairParticle = nullptr;

	fprintf(stdout,"total number of particles = %lu, total number of binaries = %lu \n", particle.size(), BinaryList.size());
	fprintf(stdout,"total number of ComputationList = %lu\n", ComputationList.size());
	fprintf(stdout,"PID=%d\n", particle[0]->PID);


	fprintf(stdout,"end of KS Regularlization Termination \n ");

	fprintf(binout,"PID=%d\n", ptclI->PID);
	fprintf(binout, "\nPosition: ptclI - x:%e, y:%e, z:%e, \n", ptclI->Position[0], ptclI->Position[1], ptclI->Position[2]);
	fprintf(binout, "Velocity: ptclI - vx:%e, vy:%e, vz:%e, \n", ptclI->Velocity[0], ptclI->Velocity[1], ptclI->Velocity[2]);
	//fflush(binout);

	fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclI->a_tot[0][0], ptclI->a_tot[1][0], ptclI->a_tot[2][0]);
	fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclI->a_tot[0][1], ptclI->a_tot[1][1], ptclI->a_tot[2][1]);
	fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclI->a_tot[0][2], ptclI->a_tot[1][2], ptclI->a_tot[2][2]);
	fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclI->a_tot[0][3], ptclI->a_tot[1][3], ptclI->a_tot[2][3]);
	fprintf(binout, "Irr  Acceleration - ax:%e, ay:%e, az:%e, \n", ptclI->a_irr[0][0], ptclI->a_irr[1][0], ptclI->a_irr[2][0]);
	fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclI->a_irr[0][1], ptclI->a_irr[1][1], ptclI->a_irr[2][1]);
	fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclI->a_irr[0][2], ptclI->a_irr[1][2], ptclI->a_irr[2][2]);
	fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclI->a_irr[0][3], ptclI->a_irr[1][3], ptclI->a_irr[2][3]);
	fprintf(binout, "Time Steps - irregular:%e, regular:%e \n", ptclI->TimeStepIrr, ptclI->TimeStepReg);
	fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", ptclI->TimeBlockIrr, ptclI->TimeBlockReg);
	fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", ptclI->CurrentBlockIrr, ptclI->CurrentBlockReg);
	//fflush(binout);

	fprintf(binout,"PID=%d\n", ptclJ->PID);
	fprintf(binout, "\nPosition: ptclJ - x:%e, y:%e, z:%e, \n", ptclJ->Position[0], ptclJ->Position[1], ptclJ->Position[2]);
	fprintf(binout, "Velocity: ptclJ - vx:%e, vy:%e, vz:%e, \n", ptclJ->Velocity[0], ptclJ->Velocity[1], ptclJ->Velocity[2]);
	//fflush(binout);

	fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclI->a_tot[0][0], ptclI->a_tot[1][0], ptclJ->a_tot[2][0]);
	fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclJ->a_tot[0][1], ptclI->a_tot[1][1], ptclI->a_tot[2][1]);
	fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclJ->a_tot[0][2], ptclJ->a_tot[1][2], ptclJ->a_tot[2][2]);
	fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclJ->a_tot[0][3], ptclJ->a_tot[1][3], ptclJ->a_tot[2][3]);
	fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", ptclI->a_irr[0][0], ptclI->a_irr[1][0], ptclI->a_irr[2][0]);
	fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclI->a_irr[0][1], ptclI->a_irr[1][1], ptclI->a_irr[2][1]);
	fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclI->a_irr[0][2], ptclI->a_irr[1][2], ptclI->a_irr[2][2]);
	fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclI->a_irr[0][3], ptclI->a_irr[1][3], ptclI->a_irr[2][3]);
	fprintf(binout, "Time Steps - irregular:%e, regular:%e \n", ptclJ->TimeStepIrr, ptclJ->TimeStepReg);
	fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", ptclJ->TimeBlockIrr, ptclJ->TimeBlockReg);
	fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", ptclJ->CurrentBlockIrr, ptclJ->CurrentBlockReg);


	/*
	fprintf(binout, "\nPosition: ptclCM - x:%e, y:%e, z:%e, \n", ptclCM->Position[0], ptclCM->Position[1], ptclCM->Position[2]);
	fprintf(binout, "Velocity: ptclCM - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0], ptclCM->Velocity[1], ptclCM->Velocity[2]);
	fflush(binout);

	fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
	fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
	fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
	fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
	fprintf(binout, "Time Steps - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr, ptclCM->TimeStepReg);
	fprintf(binout, "Current time - irregular: %e, regular:%e \n", ptclCM->CurrentTimeIrr, ptclCM->CurrentTimeReg);
	*/
	//fflush(binout);

	fprintf(stdout,"--------------------------------------\n");
	//fflush(stdout); 

}
