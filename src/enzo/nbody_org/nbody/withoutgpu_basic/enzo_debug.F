      PROGRAM ENZO_DEBUG

*
*
*       Energy check model of the enzo-nbody code
*       Written on 2022.09.02
*       ---------------------
*

      INCLUDE 'enzo_common6.h'

      INTEGER NXTLEN,N,J,L,SNAP
      REAL*8 SCALE,DT,DTOUT,RVIR
      INTEGER NUM

*     Get neccessary information for data

      open(1, FILE='input.dat', STATUS='OLD')
      read(1,*) NXTLEN,SCALE,DT,DTOUT,RVIR,SNAP

      N = NXTLEN
            
      CALL ZERO

*     set initial data to appropriate scale

*      open(12, FILE='init.dat', STATUS='OLD')
*      DO J = 1,NXTLEN
*          read(12,*) BODY(J),X0(1,J),X0(2,J),X0(3,J),
*     &               X0DOT(1,J),X0DOT(2,J),X0DOT(3,J)
*      END DO

*      CALL SCALEY(N,RVIR)
      
*      open(14, FILE='particle1.dat', STATUS='NEW')
      
*      do L = 1,NXTLEN
*        write(14,21) BODY(L),X0(1,L),X0(2,L),X0(3,L),
*     &               X0DOT(1,L),X0DOT(2,L),X0DOT(3,L)   
*      end do

*      close(14)

*   21 format (8(F20.10))  

*     start the repeating process

      DO NUM = 1,SNAP
          CALL ENZO_NBODY6_DEBUG(NUM,N,NXTLEN,SCALE,DT,DTOUT,RVIR)
      END DO

*      CALL UNSCALE(N)


      close(1)
      close(9)
*      close(12)
      close(11)
    
      END

