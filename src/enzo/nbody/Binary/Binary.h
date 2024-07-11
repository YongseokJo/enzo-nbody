#include <vector>
#include <iostream>
#include "../defs.h"

class Particle;
class Binary
{
	private:


	public:

		// corresponding single particles and neighbor particles
		// overlapping variables with ptclCM are deleted. 
		// int NumberOfAC;     // number of neighbors - not currently in use
		// double RadiusOfAC;
		// std::vector<Particle*> ACList;     // list of AC neighbor - not currently in use

		Particle* ptclCM;
		bool isTerminate;
		bool isErase;

		// information of binary particles in cartesian coordinates
		double Mass;

		// double Position[Dim];
		// double Velocity[Dim];
		// double PredPosition[Dim];
		// double PredVelocity[Dim];
		double PredTime;
		double CurrentTime;  // this show how much the binary system has evolved
		double TimeStep;
		int TimeLevel;

		// double a_tot[Dim][HERMITE_ORDER];
		// double a_reg[Dim][HERMITE_ORDER];
		// double a_irr[Dim][HERMITE_ORDER];
		// double BackgroundAcceleration[Dim];


		// information of binary particles after Levi-civita transformation
		double u[4];
		double udot[4];
		double u_pred[4];
		double udot_pred[4];
		double u2dot[4];
		double u3dot[4];
		double u4dot[4];
		double u5dot[4];

		// Q = (L^T)F_perturb : Levi-civita transformation of perturbing force

		double Q[4];
		double Qdot[4];

		// r: distance between two binaries

		double r;
		double r0; // initial separation between binaries
		double rdot;
		double r2dot;
		double r3dot;
		double r4dot;
		double r5dot;

		// t: physical time. time derivatives with respect to tau
		// tdot = r, t2dot = rdot and so on

		double tdot;
		double t2dot;
		double t3dot;
		double t4dot;
		double t5dot;
		double t6dot;



		// h: energy of binary

		double h;
		double hdot;
		double h2dot;
		double h3dot;
		double h4dot;


		// tau: time unit for new coordinate system, tau = r*t

		double PredTau;
		double CurrentTau;
		double dTau;
		// double TauStep;

		// other variables

		double gamma;       // gamma: relative energy ratio of binary
		double a;       // semi-major axis of orbit

		// stumpff coefficient related variables

		double cn[stumpffN+1];
		double cn_4z[stumpffN+1];


		// Constructor
		Binary(void) {

			//std::cout << "Constructor called" << std::endl;

			isTerminate = false;
			isErase     = false;
			// NumberOfAC     = 0; // number of neighbors
			// RadiusOfAC     = InitialRadiusOfAC;
			// Mass           = 0;
			PredTime       = 0;
			CurrentTime    = 0;
			TimeStep       = 0;
			TimeLevel      = 0;

			PredTau       = 0;
			CurrentTau    = 0;
			dTau          = 9999;
			// TauStep       = 9999;

			gamma = 0;
			a     = 0;

			r0    = 0; // initial separation between binaries
			r     = 0;
			rdot  = 0;
			r2dot = 0;
			r3dot = 0;
			r4dot = 0;
			r5dot = 0;

			tdot = 0;
			t2dot = 0;
			t3dot = 0;
			t4dot = 0;
			t5dot = 0;
			t6dot = 0;


			h     = 0;
			hdot  = 0;
			h2dot = 0;
			h3dot = 0;
			h4dot = 0;

			// for (int i=0; i<Dim; i++) {

			// 	Position[i]     = 0;
			// 	Velocity[i]     = 0;
			// 	PredPosition[i] = 0;
			// 	PredVelocity[i] = 0;

			// 	BackgroundAcceleration[i] = 0;

			// 	for (int j=0; j<HERMITE_ORDER; j++) {
			// 		a_reg[i][j] = 0;
			// 		a_irr[i][j] = 0;
			// 		a_tot[i][j] = 0;
			// 	}
			// }  -> it is all saved in the CM particle information

			for (int i=0; i<4; i++){

				u[i]          = 0;
				udot[i]       = 0;
				u_pred[i]     = 0;
				udot_pred[i]  = 0;

				u2dot[i]      = 0;
				u3dot[i]      = 0;
				u4dot[i]      = 0;
				u5dot[i]      = 0;

				Q[i]           = 0;
				Qdot[i]        = 0;
			}

			for (int i=0; i<stumpffN+1; i++) {
				cn[i]    = 1.0;
				cn_4z[i] = 1.0;
			}
		}

		void InitializeBinary(double current_time);
		void getStumpffCoefficients(double z);
		void KSIntegration(double next_time, int &calnum);
		void predictBinary(double next_time);
		void IntegrateBinary(double next_time);

		~Binary() {};

};
