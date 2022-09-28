      PROGRAM ENZO_DEBUG

*
*
*       Energy check model of the enzo-nbody code
*       Written on 2022.09.02
*       ---------------------
*

      INCLUDE 'enzo_common6.h'
      INTEGER NXTLEN,N
      REAL*8 SCALE,DT,DTOUT
      INTEGER NUM

*     Get neccessary information for data

      open(1, FILE='input.dat', STATUS='OLD')
      read(1,*) NXTLEN,SCALE,DT,DTOUT,RVIR

      N = NXTLEN

*     set necessary parameters

      CALL ZERO

      DO NUM = 1,98
          CALL ENZO_NBODY6_DEBUG(NUM,N,NXTLEN,SCALE,DT,DTOUT)
      END DO

      close(1)
      close(9)
      close(11)
    
      END

