#include <vector>
#include <mpi.h>

extern std::vector<int> LevelList;
extern int NNB;
extern double global_time;
extern double dt_min;
//extern bool debug;
extern char* fname;
extern bool restart;
// Enzo to Nbody
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTimeStep;
//extern double MinTimeStep;
//
//

// MPI variables
extern MPI_Comm inter_comm;
extern MPI_Comm nbody_comm;

