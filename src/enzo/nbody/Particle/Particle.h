/*
 *  Purporse: Particle Class
 *
 *  Date    : 2024.01.02  by Yongseok Jo
 *
 */

#include <vector>
#include <iostream>
#include "../defs.h"
//#include "../global.h"

class Particle
{
	private:

		// Variables for KS Regularization (117p in GNS by Aarseth, will be added)	
		//
		//

	public:
		int PID;
		int ParticleType;
		double Mass;
		double InitialMass;
		double CreationTime;
		double DynamicalTime;
		double Velocity[Dim];
		double Position[Dim];
		double PredTime;
		double PredTimeIrr;
		double PredTimeReg;
		double CurrentTimeIrr;
		double CurrentTimeReg;
		double TimeStepIrr;
		double TimeStepReg;
		int TimeLevelIrr;
		int TimeLevelReg;
		double PredPosition[Dim];
		double PredVelocity[Dim];
		double a_tot[Dim][HERMITE_ORDER];
		double a_reg[Dim][HERMITE_ORDER];
		double a_irr[Dim][HERMITE_ORDER];
		double BackgroundAcceleration[Dim];
		Particle* NextParticleInEnzo;
		Particle* NextParticleForComputation;
		std::vector<Particle*> ACList;     // list of AC neighbor 
		int NumberOfAC; // number of neighbors
		double RadiusOfAC;
		int isEvolve;
		bool isRegular;
		bool isStarEvolution;

		// Constructor
		Particle(void) {
			//std::cout << "Constructor called" << std::endl;
			Mass            = 0;
			InitialMass     = 0;
			NumberOfAC      = 0; // number of neighbors
			RadiusOfAC      = InitialRadiusOfAC;
			ParticleType    = -9999;
			CurrentTimeIrr  = 0.; // consistent with actual current time
			CurrentTimeReg  = 0.;
			PredTimeIrr     = 0;
			PredTimeReg     = 0;
			TimeStepIrr     = 0;
			TimeStepReg     = 0;
			TimeLevelIrr    = 9999;
			TimeLevelReg    = 9999;
			isEvolve        = 0;
			isRegular       = false;
			isStarEvolution = true;
			for (int i=0; i<Dim; i++) {
				Velocity[i]     = 0;
				Position[i]     = 0;
				PredPosition[i] = 0;
				PredVelocity[i] = 0;
				BackgroundAcceleration[i] = 0;
				for (int j=0; j<HERMITE_ORDER; j++) {
					a_reg[i][j] = 0;
					a_irr[i][j] = 0;
				}
			}
			NextParticleInEnzo = nullptr;
			NextParticleForComputation = nullptr;
		}

		void updateParticle(double mass, double *vel, double pos[], int particletype) {

			Mass = mass;
			ParticleType = particletype;

			for (int i=0; i<Dim; i++) {
				Velocity[i] = vel[i];
				Position[i] = pos[i];
			}
		}
		double getDistanceTo(Particle *particle);
		void setParticleInfo(double *data, int PID);
		void setParticleInfo(double *data, int PID, Particle* NextParticleInEnzo);
		void setParticleInfo(int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
			 	double *Position[Dim], double *Velocity[Dim],
			 	double *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i);
		void setParticleInfo(int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
			 	double *Position[Dim], double *Velocity[Dim],
			 	double *BackgroundAcceleration[Dim],  int ParticleType, Particle* NextParticleInEnzo, int i);
		void setParticleInfo(int *PID, double *BackgroundAcceleration[Dim],
			 	Particle* NextParticleInEnzo, int i);
		void setParticleInfo(double *Mass, double *BackgroundAcceleration[Dim],
			 	Particle* NextParticleInEnzo, int i);
		void initializeTimeStep();
		int getPID() {return PID;};
		void calculateIrrForce();
		void calculateRegAccelerationSecondOrder(std::vector<Particle*> &particle);
		void calculateRegAccelerationFourthOrder(std::vector<Particle*> &particle);
		void predictParticleSecondOrder(double next_time);
		void correctParticleFourthOrder(double next_time, double a[3][4]);
		void normalizeParticle();
		void calculateTimeStepIrr(double f[3][4], double df[3][4]);
		void calculateTimeStepReg(double f[3][4], double df[3][4]);
		bool checkNeighborForEvolution();
		void updateEvolveParticle(std::vector<Particle*> &particle);
		void updateParticle(double next_time, double a[3][4]);
		double evolveStarMass(double t1, double t2);
};



