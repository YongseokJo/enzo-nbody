    SUBROUTINE KCPERT(I,I1,FIRR,FD)


!       Perturbation on KS components due to chain.
!       -------------------------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH, &
    LISTC(LMAX)
    REAL*8 ::  XI(3),XIDOT(3),FIRR(3),FD(3),FP(3),FPD(3),FPS(3),FDS(3)


    NNB1 = LIST(1,I1) + 1
!       See if perturber list contains chain c.m. body #ICH.
    DO 20 L = 2,NNB1
        J = LIST(L,I1)
    !       Finish search if perturber index exceeds chain c.m. index.
        IF (J > ICH) GO TO 30
    
        IF (J == ICH) THEN
        !       Obtain force & derivative on each KS component (saving the first).
            J1 = I1
            DO 10 KCOMP = 1,2
                DO 5 K = 1,3
                    FP(K) = 0.0D0
                    FPD(K) = 0.0D0
                    XI(K) = X(K,J1)
                    XIDOT(K) = XDOT(K,J1)
                5 END DO
                CALL FCHAIN(J1,1,XI,XIDOT,FP,FPD)
                IF (KCOMP == 1) THEN
                    DO 8 K = 1,3
                        FPS(K) = FP(K)
                        FDS(K) = FPD(K)
                    8 END DO
                END IF
                J1 = J1 + 1
            10 END DO
        
        !       Add mass-weighted contributions to irregular force & derivative.
            BODYIN = 1.0/BODY(I)
            DO 15 K = 1,3
                FIRR(K) = FIRR(K) + (BODY(I1)*FPS(K) + &
                BODY(J1)*FP(K))*BODYIN
                FD(K) = FD(K) + (BODY(I1)*FDS(K) + &
                BODY(J1)*FPD(K))*BODYIN
            15 END DO
        END IF
    20 END DO

    30 RETURN

    END SUBROUTINE KCPERT
