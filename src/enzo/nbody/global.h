#include <vector>
#include <mpi.h>
#include <iostream>
#include "Particle/Particle.h"

extern std::vector<int> LevelList;
extern int NNB;
extern double global_time;
extern const double dt_min=0.03125;
extern const int dt_level_min=-5;
extern double dt_block;
extern int dt_block_level;


//extern bool debug;
extern char* fname;
extern bool restart;


// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTime;
extern bool IsEnzoCommunication;


// i/o
extern bool IsOutput;
extern double outputTime;
extern double outputTimeStep;
//
//


// MPI variables
extern MPI_Comm inter_comm;
extern MPI_Comm nbody_comm;

