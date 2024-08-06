#include <vector>
#include <mpi.h>
#include <iostream>
#include "Particle/Particle.h"
#include "Binary/Binary.h"
#include "defs.h"
#include <stdio.h>
#include <stdexcept>

#ifdef time_trace
#include "TimeTrace.h"
extern TimeTracer _time;
#endif


/*********************************************************************
 *  General variables
 *********************************************************************/
extern int NNB;
extern int newNNB;
extern char* fname;
extern bool restart;
//extern bool debug;

extern double eta;
extern double EPS2;
extern double InitialRadiusOfAC;
extern int    FixNumNeighbor;
extern int    MaxNumNeighbor;
extern int    FixNumNeighbor0;

extern int BinaryRegularization;
extern double KSTime;
extern double KSDistance;

/*********************************************************************
 *  Time related variables
 *********************************************************************/
extern double global_time;
extern double global_time_irr;

extern int time_block;
extern double time_step;
extern ULL block_max;
extern ULL NextRegTimeBlock;

extern double binary_time;
extern double binary_time_prev;
extern ULL binary_block;

/*********************************************************************
 *  Computation chain variables
 *********************************************************************/
extern std::vector<Particle*> ComputationChain;
extern Particle* FirstComputation;
extern std::vector<Particle*> ComputationList;
extern int ComputationTimeMarker;
extern std::vector<Particle*> RegularList; // how about changing this to particle list
extern std::vector<Particle*> BinaryCandidateList;
extern std::vector<Binary*> BinaryList; // List of binaries to calculate


/*********************************************************************
 *  ENZO-NBODY+ Communication variables
 *********************************************************************/
extern Particle* FirstEnzoParticle;
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTime;
extern double EnzoTimeStep;
extern double ClusetrRadius2;
extern double ClusterAcceleration[Dim];
extern double ClusterPosition[Dim];
extern double ClusterVelocity[Dim];

extern MPI_Comm inter_comm;
extern MPI_Comm nbody_comm;



/*********************************************************************
 *  I/O variables
 *********************************************************************/
extern bool IsOutput;
extern double outputTime;
extern double outputTimeStep;
extern int outNum;

extern FILE* binout; // binary output
extern FILE* nbpout; // nbody+ output
extern FILE* gpuout; // cuda output
extern FILE* pfmout; // performance output



/*********************************************************************
 *  ENZO variables (mostly subgrid physics)
 *********************************************************************/
extern int StarParticleFeedback;
extern int HydroMethod;
extern double StarMassEjectionFraction;
extern double EnzoCurrentTime;
extern int IdentifyNbodyParticles;
extern int IdentifyOnTheFly;
extern double EnzoClusterPosition[Dim+1];


/*********************************************************************
 *  Cosmology
 *********************************************************************/
extern int ComovingCoordinates;

