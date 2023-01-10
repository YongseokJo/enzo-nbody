      SUBROUTINE ENZO_NBODY6_DEBUG(NUM,N,NXTLEN,SCALE,DT,DTOUT)

*
*             N B O D Y 6++
*             *************
*
*       Regularized AC N-body code with only Hermite Integration with E check
*       Written on 2022.09.02


      INCLUDE 'enzo_common6.h'
      INTEGER NXTLEN,NUM,NEXTNUM,N
      REAL*8 SCALE,DT,DTOUT
      REAL*8 KE(NMAX),PE(NMAX),ME(NMAX)
      INTEGER I,J,K,L,M
      CHARACTER(len = 20) FILENAME,OLDFILE

*     Initialize the dataset

      NEXTNUM = NUM+1
      
      if (NEXTNUM.LT.10) then
          write (FILENAME, "(A8,I1,A4)") "particle", NEXTNUM, '.dat'
      else
          write (FILENAME, "(A8,I2,A4)") "particle", NEXTNUM, '.dat'
      end if


      if (NUM.LT.10) then
          write (OLDFILE, "(A8,I1,A4)") "particle", NUM, '.dat'
      else
          write (OLDFILE, "(A8,I2,A4)") "particle", NUM, '.dat'
      end if
      
      CALL ZERO

      DO K = 1, NXTLEN

          X0(1,K) = 0.0D0
          X0(2,K) = 0.0D0
          X0(3,K) = 0.0D0
          X0DOT(1,K) = 0.0D0
          X0DOT(2,K) = 0.0D0
          X0DOT(3,K) = 0.0D0

      END DO

*     Read Position and Velocity

      open(5, FILE=OLDFILE, STATUS='OLD')

      DO J = 1,NXTLEN
          read(5,*) BODY(J),X0(1,J),X0(2,J),X0(3,J),
     &              X0DOT(1,J),X0DOT(2,J),X0DOT(3,J)

      END DO


*     Carry out the NBODY integration
      
      CALL FORCEINT(SCALE,DT,NXTLEN)

      CALL FORCECOR(SCALE,DT,NXTLEN)

      CALL XVPRED(SCALE,DT,DTOUT,NXTLEN)

      CALL ENERGYCHECK(SCALE,NXTLEN,KE,PE,ME)
      
*     Write down the changes

      open(7, FILE=FILENAME,STATUS='NEW')
      open(9, FILE='energy.dat',STATUS='OLD')
      open(11, FILE='totalenergy.dat',STATUS='OLD')

      do L = 1,NXTLEN
        write(7,1000) BODY(L),XF(1,L),XF(2,L),XF(3,L),
     &                XFDOT(1,L),XFDOT(2,L),XFDOT(3,L)   
      end do
 1000 format (8(F12.8))
      
      write(9,*) 'step number',NUM
      write(9,2000) ALLKE,ALLPE,ALLME

      write(11,2000) ALLKE,ALLPE,ALLME

      do L = 1,NXTLEN
        write(9,2000) KE(L),PE(L),ME(L)  
      end do
 2000 format (3(F12.8))



*     Close the files

      close(5)
      close(7)

      END