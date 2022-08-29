      SUBROUTINE FORCECOR(SCALE,DT,NXTLEN)

*
*
*     Force Correction of Particles w/o GPU
*     Modification complete on 2022.08.16
*     --------------------
*     Common variables (inputs and outputs)

      INCLUDE 'enzo_common6.h'

      INTEGER I,K,NXTLEN
      REAL*8  SCALE,DT

*     Local variables

      REAL*8 DTSQ,DT6,DT2,DTSQ12,DT13
      REAL*8 DF,SUM,AT3,BT2

*     jpred_int -> predict position of single particle (future) this becomes X and Xdot

      call JPRED_INT(SCALE,DT,NXTLEN)

      call FPRED_INT(SCALE,DT,NXTLEN)

      DTSQ = DT**2
      DT6 = 6.0/(DT*DTSQ)
      DT2 = 2.0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DT13 = ONE3*DT
      
      DO I = 1,NXTLEN
          DO K = 1,3
             DF = F0(K,I) - F(K,I)
             SUM = F0DOT(K,I) + FDOT(K,I)
             AT3 = 2.0*DF + DT*SUM
             BT2 = -3.0*DF - DT*(SUM + FDOT(K,I))

             D0(K,I) = F0(K,I)
             D1(K,I) = F0DOT(K,I)
             D3(K,I) = AT3*DT6
             D2(K,I) = (3.0*AT3 + BT2)*DT2
          END DO
      END DO

      RETURN
*
      END