/***********************************************************************
	/
	/  FIND ALL STAR PARTICLES OVER ALL PROCESSORS
	/
	/  written by: Yongseok Jo
	/  date:       November, 2022
	/
 ************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"




int FindTotalNumberOfNbodyParticles(LevelHierarchyEntry *LevelArray[],int *LocalNumberOfNbodyParticles);
int FindTotalNumberOfNbodyParticles(LevelHierarchyEntry *LevelArray[],int *LocalNumberOfNbodyParticles, int *NewLocalNumberOfNbodyParticles);
int FindStartIndex(int* LocalNumberOfNbodyParticles);
void InitializeNbodyArrays(void);
void CopyNbodyArrayToOld(void);
void MatchAccelerationWithIndex(void);
void DeleteNbodyArrays(void);
void Scan(int *in, int *inout, int *len, MPI_Datatype *dptr);
void DeleteNbodyArrays(void);

