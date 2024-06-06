#include <vector>
#include <mpi.h>
#include <iostream>
#include "Particle/Particle.h"
#include "defs.h"

#ifdef time_trace
#include "TimeTrace.h"
extern TimeTracer _time;
#endif


// General
extern int NNB;
extern int newNNB;
extern int NumNeighborMax;



//Time
extern double global_time;
extern ULL NextRegTimeBlock;
extern int time_block;
extern double time_step;
extern ULL block_max;

//Compuatation Chain
extern std::vector<Particle*> ComputationChain;
extern Particle* FirstComputation;
extern std::vector<Particle*> ComputationList;
extern std::vector<int> RegIndexList; // how about changing this to particle list 6/4/2024


//extern bool debug;
extern char* fname;
extern bool restart;


// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTime;
extern double EnzoTimeStep;
extern double ClusetrRadius2;
extern double ClusterAcceleration[Dim];
extern double ClusterPosition[Dim];
extern double ClusterVelocity[Dim];

extern double eta;
extern double EPS2;
extern double InitialRadiusOfAC;

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


