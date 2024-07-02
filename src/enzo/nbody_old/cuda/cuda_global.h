#include "cuda_types.h"

#define THREAD 1024 // 64, 96, 128 or 192
#define BLOCK 64 // 8800GTS/512 has 16
//#define NIBLOCK 16 // 16 or 32 
//#define NIMAX (NTHREAD * NIBLOCK) // 1024

//#define NB_PER_BLOCK 256 // NNB per block
//#define NB_BUF_SIZE (1<<21)

/*
extern double time_send, time_grav, time_out, time_nb;
extern long long numInter;
extern int icall,ini,isend;

extern cudaPointer <BackgroundParticle> background;
//extern static cudaPointer <TargetParticle> target;
//static cudaPointer <uint16[NJBLOCK][NB_PER_BLOCK]>nbpart;
//extern static cudaPointer <int> nblist;
//extern static cudaPointer <int> nboff;

//static myvector<int> nblist;
extern int NNB, nbodymax;
extern int devid, numGPU;
extern bool is_open;
extern bool devinit;
// static int *nblist;
*/

extern int NNB;
extern double time_send, time_grav, time_out, time_nb;
extern long long numInter;
extern int icall,ini,isend;
extern int nbodymax;
extern int devid, numGPU;
extern bool is_open;
extern bool devinit;
extern cudaPointer <BackgroundParticle> background;

typedef unsigned short uint16;
