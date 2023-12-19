#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
 

void  enzo_nb(int comm_enzo, int nbody_comm, int mpi_rank) {
	int EN, IS, IE, EID, I, J;	
	float * ENBODY, EX, EXODT, EF;
	float EMU, ELU, EVU, ETU;
	float LENGTHU,MASSU,VELU,TIMEU;
	float EDT;


	// MPI receive here?


	// unit conversion?
	
	// read input parameters  (start.F)

  // Adjust.for


	// Intgrt.for


	//Bunch of Iphases



}


