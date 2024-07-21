#include "../global.h"
#include <cmath>

void generate_Matrix(double a[3], double (&A)[3][4]);

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
	//this->RadiusOfAC     = InitialRadiusOfAC;
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
	//this->RadiusOfAC     = InitialRadiusOfAC;
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
	//this->RadiusOfAC     = InitialRadiusOfAC;
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
	//this->RadiusOfAC     = InitialRadiusOfAC;
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
	//this->RadiusOfAC     = InitialRadiusOfAC;
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
	//this->RadiusOfAC     = InitialRadiusOfAC;
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
	//this->RadiusOfAC     = InitialRadiusOfAC;
}


void Particle::setParticleInfo(int *PID, double *Mass, double *CreationTime, double *DynamicalTime, double *Position[Dim],
		double *Velocity[Dim], double *BackgroundAcceleration[Dim], int ParticleType, Particle* NextParticleInEnzo, int i) {
	this->PID												 = PID[i];
	this->Mass											 = Mass[i]*EnzoMass;
	this->InitialMass								 = this->Mass;
	this->CreationTime               = CreationTime[i]*EnzoTime;
	this->DynamicalTime              = DynamicalTime[i]*EnzoTime;
	this->Position[0]								 = Position[0][i]*EnzoLength;
	this->Position[1]                = Position[1][i]*EnzoLength;
	this->Position[2]								 = Position[2][i]*EnzoLength;
	this->Velocity[0]								 = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]								 = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]								 = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->ParticleType						   = ParticleType;
	this->NextParticleInEnzo				 = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
	//this->RadiusOfAC     = InitialRadiusOfAC;
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


void Particle::convertBinaryCoordinatesToCartesian() {
	if (!this->isCMptcl)  {
		fprintf(stderr,"This is NOT a CM particle!\n");
		return;
	}
	fprintf(stdout,"Converting the KS coordinates to physical coordinates of ptclI and ptclJ\n");
	fprintf(stderr,"Converting the KS coordinates to physical coordinates of ptclI and ptclJ\n");

	double R[Dim], Rdot[Dim];
	double Rinv;
	double ratioM;
	double L[3][4];

	// update the values of positions of BinaryParticleI and ptcl J
	R[0]   = BinaryInfo->u[0]*BinaryInfo->u[0] - BinaryInfo->u[1]*BinaryInfo->u[1] - BinaryInfo->u[2]*BinaryInfo->u[2] + BinaryInfo->u[3]*BinaryInfo->u[3];
	R[1]   = 2*(BinaryInfo->u[0]*BinaryInfo->u[1] - BinaryInfo->u[2]*BinaryInfo->u[3]);
	R[2]   = 2*(BinaryInfo->u[0]*BinaryInfo->u[2] + BinaryInfo->u[1]*BinaryInfo->u[3]);
	ratioM = BinaryParticleJ->Mass/this->Mass;


	for (int dim=0; dim<Dim; dim++) {
		BinaryParticleI->Position[dim] = this->Position[dim] + ratioM*R[dim];
		BinaryParticleJ->Position[dim] = this->Position[dim] - R[dim];
	}


	// do the same thing for velocity components
	generate_Matrix(BinaryInfo->u,L);

	Rinv = 1/(BinaryInfo->u[0]*BinaryInfo->u[0] + BinaryInfo->u[1]*BinaryInfo->u[1] + BinaryInfo->u[2]*BinaryInfo->u[2] + BinaryInfo->u[3]*BinaryInfo->u[3]) ;


	for (int dim=0; dim<Dim; dim++) {
		Rdot[dim] = 0.0;
		for (int dimu=0; dimu<4; dimu++) {
			Rdot[dim] += 2*L[dim][dimu]*BinaryInfo->udot[dim]*Rinv;
		}
	}


	for (int dim=0; dim<Dim; dim++) {
		BinaryParticleI->Velocity[dim] = this->Velocity[dim] + ratioM*Rdot[dim];
		BinaryParticleJ->Velocity[dim] = BinaryParticleI->Velocity[dim] - Rdot[dim];
	}

	BinaryParticleI->CurrentBlockIrr = this->CurrentBlockIrr;
	BinaryParticleI->CurrentBlockReg = this->CurrentBlockReg;
	BinaryParticleI->CurrentTimeIrr = this->CurrentBlockIrr*time_step;
	BinaryParticleI->CurrentTimeReg = this->CurrentBlockReg*time_step;

	BinaryParticleI->TimeStepIrr     = this->TimeStepIrr;
	BinaryParticleI->TimeBlockIrr    = this->TimeBlockIrr;
	BinaryParticleI->TimeLevelIrr    = this->TimeLevelIrr;

	BinaryParticleI->TimeStepReg     = this->TimeStepReg;
	BinaryParticleI->TimeBlockReg    = this->TimeBlockReg;
	BinaryParticleI->TimeLevelReg    = this->TimeLevelReg;

	BinaryParticleJ->CurrentBlockIrr = BinaryParticleI->CurrentBlockIrr;
	BinaryParticleJ->CurrentBlockReg = BinaryParticleI->CurrentBlockReg;
	BinaryParticleJ->CurrentTimeIrr = BinaryParticleI->CurrentTimeIrr;
	BinaryParticleJ->CurrentTimeReg = BinaryParticleI->CurrentTimeReg;

	BinaryParticleJ->TimeStepIrr     = this->TimeStepIrr;
	BinaryParticleJ->TimeBlockIrr    = this->TimeBlockIrr;
	BinaryParticleJ->TimeLevelIrr    = this->TimeLevelIrr;

	BinaryParticleJ->TimeStepReg     = this->TimeStepReg;
	BinaryParticleJ->TimeBlockReg    = this->TimeBlockReg;
	BinaryParticleJ->TimeLevelReg    = this->TimeLevelReg;

	fprintf(stdout,"END CONVERTING THE COORDINATES\n \n");
}


void Particle::ComputeAcceleration(Particle* ptcl, double (&a)[Dim], double (&adot)[Dim]) {
	double dx[Dim], dv[Dim];
	double dr2;
	double dxdv;
	double m_r3;

	if (this->PID == ptcl->PID) {
		return;
	}

	dr2  = 0.0;
	dxdv = 0.0;

	for (int dim=0; dim<Dim; dim++) {
		dx[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
		dv[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];
		dr2    += dx[dim]*dx[dim];
		dxdv   += dx[dim]*dv[dim];
	}


	if (dr2 < EPS2) {
		dr2 =  EPS2;
	}

	m_r3 = ptcl->Mass/dr2/sqrt(dr2);

	for (int dim=0; dim<Dim; dim++){
		a[dim]    = m_r3*dx[dim];
		adot[dim] = m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
	}

}



