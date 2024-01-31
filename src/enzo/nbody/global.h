#include <vector>
#include <mpi.h>
#include <iostream>
#include "Particle/Particle.h"

extern std::vector<int> LevelList;
extern int NNB;
extern int newNNB;
extern double global_time;
extern double NextRegTime;
extern const double dt_min;
extern const int dt_level_min;
extern double dt_block;  // this stores the minimum time step
extern int dt_block_level;
extern std::vector<Particle*> ComputationChain;
extern Particle* FirstComputation;


//extern bool debug;
extern char* fname;
extern bool restart;


// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTime;
extern double EnzoTimeStep;


// i/o
extern bool IsOutput;
extern double outputTime;
extern double outputTimeStep;
//
//


// MPI variables
extern MPI_Comm inter_comm;
extern MPI_Comm nbody_comm;


// ENZO variables
extern int StarParticleFeedback;
extern int HydroMethod;
extern double StarMassEjectionFraction;
extern double EnzoCurrentTime;


