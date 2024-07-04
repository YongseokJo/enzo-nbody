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

class Binary;
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
		double Mdot;
		double InitialMass;
		double CreationTime;
		double DynamicalTime;
		double Velocity[Dim];
		double Position[Dim];
		double NewVelocity[Dim];
		double NewPosition[Dim];
		double PredTime;
		double PredTimeIrr;
		double PredTimeReg;
		double CurrentTimeIrr;
		double CurrentTimeReg;
		ULL CurrentBlockIrr;
		ULL CurrentBlockReg;
		double TimeStepIrr;
		double TimeStepReg;
		ULL TimeBlockIrr;
		ULL TimeBlockReg;
		int TimeLevelIrr;
		int TimeLevelReg;
		double PredMass;
		double PredPosition[Dim];
		double PredVelocity[Dim];
		double a_tot[Dim][HERMITE_ORDER];
		double a_reg[Dim][HERMITE_ORDER];
		double a_irr[Dim][HERMITE_ORDER];
		double BackgroundAcceleration[Dim];
		Particle* NextParticleInEnzo;
		Particle* NextParticleForComputation;
		Particle* BinaryPairParticle;
		Particle* BinaryParticleI;
		Particle* BinaryParticleJ;
		Binary* BinaryInfo;
		std::vector<Particle*> ACList;     // list of AC neighbor 
		int NumberOfAC; // number of neighbors
		double RadiusOfAC;
		bool isStarEvolution;
		bool isBinary; // check whether this is a member of the binary
		bool isCMptcl; // check if this particle is center-of-mass particle
		bool isErase;





		// Constructor
		Particle() {__initialize__();};
		Particle(int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
				double *Position[Dim], double *Velocity[Dim],
				double *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i);
		Particle(
				int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
				double *Position[Dim], double *Velocity[Dim],
				double *BackgroundAcceleration[Dim],  int ParticleType, Particle* NextParticleInEnzo, int i);

		void __initialize__(void) {
			//std::cout << "Constructor called" << std::endl;
			Mass            = 0;
			Mdot            = 0;
			InitialMass     = 0;
			NumberOfAC      = 0; // number of neighbors
			RadiusOfAC      = -1;
			ParticleType    = -9999;
			CurrentTimeIrr  = 0.; // consistent with actual current time
			CurrentTimeReg  = 0.;
			CurrentBlockIrr = 0; // consistent with actual current time
			CurrentBlockReg = 0;
			PredTimeIrr     = 0;
			PredTimeReg     = 0;
			TimeStepIrr     = 0;
			TimeStepReg     = 0;
			TimeLevelIrr    = 0;
			TimeLevelReg    = 0;
			TimeBlockIrr    = 0;
			TimeBlockReg    = 0;
			isStarEvolution = false;
			isBinary        = false;
			isCMptcl        = false;
			isErase         = false;
			PredMass        = 0.;
			for (int i=0; i<Dim; i++) {
				Velocity[i]     = 0.;
				Position[i]     = 0.;
				PredPosition[i] = 0.;
				PredVelocity[i] = 0.;
				BackgroundAcceleration[i] = 0.;
				for (int j=0; j<HERMITE_ORDER; j++) {
					a_tot[i][j] = 0.;
					a_reg[i][j] = 0.;
					a_irr[i][j] = 0.;
				}
			}
			NextParticleInEnzo = nullptr;
			NextParticleForComputation = nullptr;
			BinaryPairParticle = nullptr;
			BinaryParticleI = nullptr;
			BinaryParticleJ = nullptr;
			BinaryInfo      = nullptr;
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

		void predictParticleSecondOrder(double time);
		void predictParticleSecondOrderIrr(double time);
		void correctParticleFourthOrder(double current_time, double next_time, double a[3][4]);

		void normalizeParticle();
		void calculateTimeStepIrr(double f[3][4], double df[3][4]);
		void calculateTimeStepReg();
		bool checkNeighborForEvolution();
		void updateParticle();
		double evolveStarMass(double t1, double t2);
		void isKSCandidate();
		void convertBinaryCoordinatesToCartesian();
		void polynomialPrediction(double current_time);

		//destructor
    ~Particle() = default;
};



