      SUBROUTINE JPRED_INT(I,X,XDOT,F,FD,SCALE,DT,NMAX,NXTLEN
     &     XP,XPDOT)
*
*
*       Prediction of single particle, used in intgrt.F
*       ----------------------------------------
*       Common variables (inputs and outputs)

      INTEGER I,NMAX,NXTLEN
      REAL*8  X(3,NMAX),XDOT(3,NMAX),F(3,NMAX),FD(3,NMAX),
      REAL*8  BODY(NMAX),SCALE,DT
      REAL*8  XP(3),XPDOT(3)

      S = DT
      S1 = 1.5*S
      S2 = 2.0*S

*     Note X may become corrected value X0 for zero interval.
      XP(1) = ((FDOT(1,I)*S + F(1,I))*S + XDOT(1,I))*S + X(1,I)
      XP(2) = ((FDOT(2,I)*S + F(2,I))*S + XDOT(2,I))*S + X(2,I)
      XP(3) = ((FDOT(3,I)*S + F(3,I))*S + XDOT(3,I))*S + X(3,I)
      XPDOT(1) = (FDOT(1,I)*S1 + F(1,I))*S2 + XDOT(1,I)
      XPDOT(2) = (FDOT(2,I)*S1 + F(2,I))*S2 + XDOT(2,I)
      XPDOT(3) = (FDOT(3,I)*S1 + F(3,I))*S2 + XDOT(3,I)
*
      RETURN
*
      END
