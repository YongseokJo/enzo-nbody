    SUBROUTINE FPERT(I,J,NP,PERT)


!       Perturbing force on dominant bodies.
!       ------------------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    REAL*8 ::  FP(3)


    DO 1 K = 1,3
        FP(K) = 0.0D0
    1 END DO
    CALL JPRED(I,TIME,TIME)
    CALL JPRED(J,TIME,TIME)


!       Obtain perturbation on body #I & J due to NP members of JLIST.
    DO 10 KK = 1,NP
        K = JLIST(KK)
        IF (K == I .OR. K == J) GO TO 10
        L = I
        CALL JPRED(K,TIME,TIME)
    !$$$          if (K.GE.IFIRST.and.TPRED(K).ne.TIME) then
    !$$$             S = TIME - T0(K)
    !$$$             X(1,K) = ((FDOT(1,K)*S + F(1,K))*S + X0DOT(1,K))*S +X0(1,K)
    !$$$             X(2,K) = ((FDOT(2,K)*S + F(2,K))*S + X0DOT(2,K))*S +X0(2,K)
    !$$$             X(3,K) = ((FDOT(3,K)*S + F(3,K))*S + X0DOT(3,K))*S +X0(3,K)
    !$$$          end if
    !$$$          if (L.GE.IFIRST.and.TPRED(L).ne.TIME) then
    !$$$             S = TIME - T0(L)
    !$$$             X(1,L) = ((FDOT(1,L)*S + F(1,L))*S + X0DOT(1,L))*S +X0(1,L)
    !$$$             X(2,L) = ((FDOT(2,L)*S + F(2,L))*S + X0DOT(2,L))*S +X0(2,L)
    !$$$             X(3,L) = ((FDOT(3,L)*S + F(3,L))*S + X0DOT(3,L))*S +X0(3,L)
    !$$$          end if
        5 A1 = X(1,K) - X(1,L)
        A2 = X(2,K) - X(2,L)
        A3 = X(3,K) - X(3,L)
        RIJ2 = A1**2 + A2**2 + A3**2
    !       Restrict evaluation to 100*RMIN (neighbour list may be used).
        IF (L == I .AND. RIJ2 > 1.0D+04*RMIN2) GO TO 10
        A4 = BODY(K)/(RIJ2*SQRT(RIJ2))
        IF (L == J) A4 = -A4
        FP(1) = FP(1) + A1*A4
        FP(2) = FP(2) + A2*A4
        FP(3) = FP(3) + A3*A4
        IF (L == I) THEN
            L = J
            GO TO 5
        END IF
    10 END DO

    PERT = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)

    RETURN

    END SUBROUTINE FPERT
