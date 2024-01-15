#include "Particle.h"
#include <cmath>

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
	this->ParticleType = Star;
}


void Particle::setParticleInfo(int *PID, double *Mass, double *Position[Dim],
		double *Velocity[Dim], double *BackgroundAcceleration[Dim], int i) {
	this->PID          = PID[i];
	this->Mass         = Mass[i];
	this->Position[0]  = Position[0][i];
	this->Position[1]  = Position[1][i];
	this->Position[2]  = Position[2][i];
	this->Velocity[0]  = Velocity[0][i];
	this->Velocity[1]  = Velocity[1][i];
	this->Velocity[2]  = Velocity[2][i];
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i];
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i];
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i];
	this->ParticleType = Star;
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


