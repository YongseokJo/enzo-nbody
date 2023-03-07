      SUBROUTINE MYDUMP(II,J)
*
*
*       COMMON save or read.
*       --------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      INCLUDE 'params.h'
      INCLUDE 'timing.h'
      INCLUDE 'mpi_base.h'
      INCLUDE 'tlist.h'
#ifdef TT
      INCLUDE 'tt.h'
#endif

*     deleted by sykim
*     we do not need mydump output
       
*
      RETURN
*
      END
