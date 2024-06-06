#define HERMITE_ORDER 4
#define Dim 3
#define NNB_MAX 1000000



// Particle Type Definition
// Particle type defined by combination of binary switches
// fictious (2^8), new (2^7), single (2^6), black hole (2^5), dead star (2^4)
// TypeIII (2^3), TypeII (2^2), TypeI (2^1), normal star or particle (2^0) 
// for example, normal signle particle is 001000001 = 2^6 + 2^0 = 65
// new black hole single particle is 011100000 = 2^7+2^6+2^5 = 
// fictious particles are usually COM particles for regularization
#define NormalStar 1
#define TypeIStar 2
#define TypeIIStar 4
#define TypeIIIStar 8
#define DeadStar 16
#define Blackhole 32
#define SingleParticle 64
#define NewParticle 128
#define FicticiousParticle 256


// Numerical 
#define MIN_LEVEL_BUFFER 5
#define FixNumNeighbor 0

#define DONE    1
#define FAIL   -1


// numerical values
#define mytolerance 5.4210109e-20


// Physical constants
#define G_cgs 6.67430e-8
#define G // pc, yr, Msun



// Physical units in cgs
#define pc 3.08567758149137e18
#define yr 3.1536e7
#define Msun 1.98847e33

// Code unit in pc, yr, Msun
#define time_unit 1e10 // in 1e10 yr
#define position_unit 4. // in 4 pc
#define velocity_unit 4e-10 // in 4e-10 pc/yr
//#define mass_unit 256e-20  // 256e-20 Msun in the unit that G = 1.
#define mass_unit 0.0001424198  // Msun in the unit that G = 1.

typedef unsigned long long ULL;

#define no_time_trace
