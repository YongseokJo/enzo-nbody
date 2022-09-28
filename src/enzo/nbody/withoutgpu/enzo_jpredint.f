      SUBROUTINE JPRED_INT(SCALE,DT,NXTLEN)
*
*
*       Prediction of single particle, used in intgrt.F
*       Written on 2022.08.16
*       ----------------------------------------
*       Common variables (inputs and outputs)

      INCLUDE 'enzo_common6.h'
      INTEGER I,NXTLEN
      REAL*8  SCALE,DT

      S = DT
      S1 = 1.5*S
      S2 = 2.0*S

*     Note X may become corrected value X0 for zero interval.
      DO I=1,NXTLEN
      
         X(1,I) = ((F0DOT(1,I)*S + F0(1,I))*S + X0DOT(1,I))*S + X0(1,I)
         X(2,I) = ((F0DOT(2,I)*S + F0(2,I))*S + X0DOT(2,I))*S + X0(2,I)
         X(3,I) = ((F0DOT(3,I)*S + F0(3,I))*S + X0DOT(3,I))*S + X0(3,I)
         XDOT(1,I) = (F0DOT(1,I)*S1 + F0(1,I))*S2 + X0DOT(1,I)
         XDOT(2,I) = (F0DOT(2,I)*S1 + F0(2,I))*S2 + X0DOT(2,I)
         XDOT(3,I) = (F0DOT(3,I)*S1 + F0(3,I))*S2 + X0DOT(3,I)

*         PRINT *,'POSITION',X(1,I),XDOT(1,I)

      END DO
*
      RETURN
*
      END
