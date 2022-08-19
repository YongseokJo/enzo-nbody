      SUBROUTINE XVPRED(NMAX,DT)
*
*
*     Prediction of coordinates & velocities w/o GPU
*     Modification complete on 2022.08.17   
*     ---------------------------------------

*     Common Variables

      INTEGER NMAX
      REAL*8 DT

*     Local Variables

      INTEGER I,K

      A1 = 0.05*DT
      A2 = 0.25*DT
      A3 = DT
      A4 = 0.5*A3
      
      DO I = 1, NXTLEN

         DO K = 1,3
             FK = ((ONE6*D3R(K,I)*A3 + 0.5*D2R(K,I))*A3
     &            + D1R(K,I))*A3 + D0R(K,I) + D0(K,I)
             F1DOTK = (D3R(K,I)*A4 + D2R(K,I))*A3 + D1R(K,I) + D1(K,I)
             F2DOTK = 0.5*(D3R(K,I)*A3 + D2R(K,I) + D2(K,I))
             F3DOTK = ONE6*(D3R(K,I) + D3(K,I))
             X(K,I) = ((((F3DOTK*A1 + ONE12*F2DOTK)*DT +
     &                 ONE6*F1DOTK)*DT + 0.5*FK)*DT + X0DOT(K,I))*DT +
     &                 X0(K,I)
             XDOT(K,I) = (((F3DOTK*A2 + ONE3*F2DOTK)*DT +
     &                   0.5*F1DOTK)*DT + FK)*DT + X0DOT(K,I)
         END DO

      END DO

      RETURN

      END