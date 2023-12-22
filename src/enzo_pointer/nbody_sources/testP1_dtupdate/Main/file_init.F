      SUBROUTINE FILE_INIT(ISTART)
*
*
*       Opening of Files with proper names.
*       ------------------------------------
*       ISTART=0: Open all other files after reading input data.
*       ISTART=1: Open only unit 1 for mydump.
*       ISTART=2: Open only unit 2 for mydump.
*
      INCLUDE 'common6.h'
*
*
      CHARACTER*10 FILE(100)
      CHARACTER*4 FRAG(100)
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif
*
*     unnecessary files deleted by sykim

      IF (KZ(22).GT.0)
     &OPEN (UNIT=10,STATUS='UNKNOWN',FORM='FORMATTED',FILE='dat.10')

*
      RETURN
      END
