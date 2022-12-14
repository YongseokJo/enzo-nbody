      SUBROUTINE XVPRED(SCALE,DT,NXTLEN)
*
*
*     Prediction of coordinates & velocities w/o GPU
*     Modification complete on 2022.08.17   
*     ---------------------------------------

*     Common Variables

      INCLUDE 'enzo_common6.h'
      INTEGER NXTLEN
      REAL*8 SCALE,DT

*     Local Variables

      INTEGER I,K
      REAL*8 A1,A2,A3,A4
      REAL*8 FK,F1DOTK,F2DOTK,F3DOTK

      A1 = 0.05*DT
      A2 = 0.25*DT
      A3 = DT
      A4 = 0.5*A3
      
      DO I = 1, NXTLEN

         DO K = 1,3
             FK = ((ONE6*D3(K,I)*A3 + 0.5*D2(K,I))*A3
     &            + D1(K,I))*A3 + D0(K,I)
             F1DOTK = (D3(K,I)*A4 + D2(K,I))*A3 + D1(K,I)
             F2DOTK = 0.5*(D3(K,I)*A3 + 2*D2(K,I))
             F3DOTK = ONE6*(D3(K,I))
             XF(K,I) = ((((F3DOTK*A1 + ONE12*F2DOTK)*DT +
     &                 ONE6*F1DOTK)*DT + 0.5*FK)*DT + X0DOT(K,I))*DT +
     &                 X0(K,I)
             XFDOT(K,I) = (((F3DOTK*A2 + ONE3*F2DOTK)*DT +
     &                   0.5*F1DOTK)*DT + FK)*DT + X0DOT(K,I)
         END DO

      END DO

      RETURN

      END