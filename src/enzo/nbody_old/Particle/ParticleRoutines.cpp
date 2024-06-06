#include "../global.h"
#include <cmath>


Particle::Particle(int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
			 	double *Position[Dim], double *Velocity[Dim],
			 	double *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i) {

			__initialize__();
	this->PID                        = PID[i];
	this->Mass                       = Mass[i]*EnzoMass;
	this->InitialMass                = this->Mass;
	this->CreationTime               = CreationTime[i]*EnzoTime;
	this->DynamicalTime              = DynamicalTime[i]*EnzoTime;
	this->Position[0]                = Position[0][i]*EnzoLength;
	this->Position[1]                = Position[1][i]*EnzoLength;
	this->Position[2]                = Position[2][i]*EnzoLength;
	this->Velocity[0]                = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]                = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]                = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->NextParticleInEnzo         = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	this->isRegular      = true;
	this->RadiusOfAC     = InitialRadiusOfAC;
}

Particle::Particle(
int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
			 	double *Position[Dim], double *Velocity[Dim],
			 	double *BackgroundAcceleration[Dim],  int ParticleType, Particle* NextParticleInEnzo, int i) {
			__initialize__();
	this->PID                        = PID[i];
	this->Mass                       = Mass[i]*EnzoMass;
	this->InitialMass                = this->Mass;
	this->CreationTime               = CreationTime[i]*EnzoTime;
	this->DynamicalTime              = DynamicalTime[i]*EnzoTime;
	this->Position[0]                = Position[0][i]*EnzoLength;
	this->Position[1]                = Position[1][i]*EnzoLength;
	this->Position[2]                = Position[2][i]*EnzoLength;
	this->Velocity[0]                = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]                = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]                = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->ParticleType               = ParticleType;
	this->NextParticleInEnzo         = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	this->isRegular      = true;
	this->RadiusOfAC     = InitialRadiusOfAC;
}


double Particle::getDistanceTo(Particle *particle) {
	double d=0;
	for (int i=0; i<Dim; i++)
		d += std::pow(this->Position[i]-particle->Position[i],2);
	return std::sqrt(d);
}

void Particle::setParticleInfo(double *data, int PID) {
	this->PID          = PID;
	this->Position[0]  = data[0];
	this->Position[1]  = data[1];
	this->Position[2]  = data[2];
	this->Velocity[0]  = data[3];
	this->Velocity[1]  = data[4];
	this->Velocity[2]  = data[5];
	this->Mass         = data[6];
	this->ParticleType = NormalStar+SingleParticle;
	this->RadiusOfAC     = InitialRadiusOfAC;
}

void Particle::setParticleInfo(double *data, int PID, Particle* NextParticleInEnzo) {
	this->PID          = PID;
	this->Position[0]  = data[0];
	this->Position[1]  = data[1];
	this->Position[2]  = data[2];
	this->Velocity[0]  = data[3];
	this->Velocity[1]  = data[4];
	this->Velocity[2]  = data[5];
	this->Mass         = data[6];
	this->ParticleType = NormalStar+SingleParticle;
	this->NextParticleInEnzo = NextParticleInEnzo;
	this->isRegular      = true;
	this->RadiusOfAC     = InitialRadiusOfAC;
}

void Particle::setParticleInfo(int *PID, double *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i) {
	this->PID                        = PID[i];
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->NextParticleInEnzo = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	this->isRegular      = true;
	this->RadiusOfAC     = InitialRadiusOfAC;
}

void Particle::setParticleInfo(double *Mass, double *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i) {
	this->Mass                       = Mass[i]*EnzoMass;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->NextParticleInEnzo = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	this->isRegular      = true;
	this->RadiusOfAC     = InitialRadiusOfAC;
}


void Particle::setParticleInfo(int *PID, double *Mass, double *CreationTime, double *DynamicalTime,
	 	double *Position[Dim], double *Velocity[Dim], double *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i) {
	this->PID                        = PID[i];
	this->Mass                       = Mass[i]*EnzoMass;
	this->InitialMass                = this->Mass;
	this->CreationTime               = CreationTime[i]*EnzoTime;
	this->DynamicalTime              = DynamicalTime[i]*EnzoTime;
	this->Position[0]                = Position[0][i]*EnzoLength;
	this->Position[1]                = Position[1][i]*EnzoLength;
	this->Position[2]                = Position[2][i]*EnzoLength;
	this->Velocity[0]                = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]                = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]                = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->ParticleType               = NormalStar+SingleParticle;
	this->NextParticleInEnzo         = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	this->isRegular      = true;
	this->RadiusOfAC     = InitialRadiusOfAC;
}


void Particle::setParticleInfo(int *PID, double *Mass, double *CreationTime, double *DynamicalTime, double *Position[Dim],
		double *Velocity[Dim], double *BackgroundAcceleration[Dim], int ParticleType, Particle* NextParticleInEnzo, int i) {
	this->PID          = PID[i];
	this->Mass         = Mass[i]*EnzoMass;
	this->InitialMass  = this->Mass;
	this->CreationTime               = CreationTime[i]*EnzoTime;
	this->DynamicalTime              = DynamicalTime[i]*EnzoTime;
	this->Position[0]  = Position[0][i]*EnzoLength;
	this->Position[1]  = Position[1][i]*EnzoLength;
	this->Position[2]  = Position[2][i]*EnzoLength;
	this->Velocity[0]  = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]  = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]  = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->ParticleType       = ParticleType;
	this->NextParticleInEnzo = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	this->isRegular      = true;
	this->RadiusOfAC     = InitialRadiusOfAC;
}


void Particle::normalizeParticle() {
	// pc to computing unit, km/s to computing unit
	Mass *= 1e9;
	Mass /= mass_unit;
	for (int dim=0; dim<Dim; dim++) {
		Position[dim] *= 1000; // kpc to pc
		Position[dim] /= position_unit;
		Velocity[dim] *= 1e5*yr/pc; // km/s to pc/yr
		Velocity[dim] /= velocity_unit;
	}
}




