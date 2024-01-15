
#define HERMITE_ORDER 4
#define Dim 3

// Particle Type Definition
#define Star 0
#define TypeIStar 1
#define TypeIIStar 2
#define TypeIIIStar 3
#define DeadStar 3

#define Blackhole 9
// Fake Particles for KS regularization
#define BlackholeBiniary -2
#define Binary -1

#define DONE    1
#define FAIL   -1


// numerical values
#define mytolerance 5.4210109e-20

// Physics related parameters
#define eta 0.001 // by YS Jo subject to modifty
#define EPS2 0.001
#define InitialRadiusOfAC 0.01 // 0.04 pc

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
