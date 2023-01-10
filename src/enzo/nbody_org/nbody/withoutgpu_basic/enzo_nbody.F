      PROGRAM ENZO_NBODY6

*
*             N B O D Y 6++
*             *************
*
*       Regularized AC N-body code with only Hermite Integration
*       Written on 2022.08.16


      INCLUDE 'enzo_common6.h'
      
      INTEGER NXTLEN
      REAL*8 SCALE,DT,DTOUT
      INTEGER I,J,K,L,M

*     Get neccessary information for data

      open(1, FILE='input.dat', STATUS='OLD')
      read(1,*) NXTLEN,SCALE,DT,DTOUT

*     Initialize the dataset

      DO K = 1, NXTLEN

          X0(1,K) = 0.0D0
          X0(2,K) = 0.0D0
          X0(3,K) = 0.0D0
          X0DOT(1,K) = 0.0D0
          X0DOT(2,K) = 0.0D0
          X0DOT(3,K) = 0.0D0

      END DO

*     Read Position and Velocity

      open(5, FILE='particle1.dat',STATUS='OLD')
      open(7, FILE='result.dat',STATUS='OLD')

      DO J = 1,NXTLEN
          read(5,*) BODY(J),X0(1,J),X0(2,J),X0(3,J),
     &              X0DOT(1,J),X0DOT(2,J),X0DOT(3,J)

      END DO


*     Carry out the NBODY integration
      
      CALL FORCEINT(SCALE,DT,NXTLEN)

      CALL FORCECOR(SCALE,DT,NXTLEN)

      CALL XVPRED(SCALE,DT,DTOUT,NXTLEN)
      
*     Write down the changes
      do L = 1,NXTLEN
        write(7,1000) BODY(L),XF(1,L),XF(2,L),XF(3,L),
     &                XFDOT(1,L),XFDOT(2,L),XFDOT(3,L)   
      end do
 1000 format (8(F10.5))



*     Close the files

      close(1)
      close(5)

      END