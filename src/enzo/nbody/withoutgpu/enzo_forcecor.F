      SUBROUTINE FORCECOR(I,X,XDOT,BODY,F,FD,SCALE,DT,NMAX,NXTLEN)

*
*
*     Force Correction of Particles w/o GPU
*     Modification complete on 2022.08.16
*     --------------------
*     Common variables (inputs and outputs)

      INTEGER I,NMAX,NXTLEN
      REAL*8  X(3,NMAX),XDOT(3,NMAX),F(3,NMAX),FD(3,NMAX),
      REAL*8  BODY(NMAX),SCALE,DT

*     Local variables

      REAL*8  XP(3),XPDOT(3),FP(3),FPDOT(3),

*     jpred_int -> predict position of single particle (future) this becomes X and Xdot

      call JPRED_INT(I,X,XDOT,F,FD,SCALE,DT,NMAX,NXTLEN,XP,XPDOT)

      call FPRED_INT(I,XP,XPDOT,BODY,SCALE,TIME,DT,NMAX,NXTLEN,FP,FPD)

      DTSQ = DT**2
      DT6 = 6.0/(DT*DTSQ)
      DT2 = 2.0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DT13 = ONE3*DT

      DO 80 K = 1,3
          DF = FP(K) - F(K,I)
          SUM = FPD(K) + FD(K,I)
          AT3 = 2.0*DF + DT*SUM
          BT2 = -3.0*DF - DT*(SUM + FID)

          D0(K,I) = F(K,I)
          D0R(K,I) = FP(K)
          D1R(K,I) = FPD(K)
          D3R(K,I) = AT3*DT6
          D2R(K,I) = (3.0*AT3 + BT2)*DT2

   80 continue

      RETURN
*
      END