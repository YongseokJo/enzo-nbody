//#include <mpi.h>


using namespace std;

//int readData(std::vector<Particle*> &particle);
int createComputationChain(std::vector<Particle*> &particle);
int InitializeParticle(std::vector<Particle*> &particle);
void Evolve(std::vector<Particle*> &particle);
int InitialCommunication(std::vector<Particle*> &particle);
//int Parser(int argc, char* argv[]);

