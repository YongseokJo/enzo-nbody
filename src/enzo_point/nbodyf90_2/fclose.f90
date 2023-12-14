    SUBROUTINE FCLOSE(I,NNB)


!       Force & first derivative from close bodies.
!       -------------------------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    REAL*8 ::  A(9)


!      CALL JPRED(I,TIME,TIME)
!       Initialize F & FDOT for body #I.
    DO 10 K = 1,3
        F(K,I) = 0.0D0
        FDOT(K,I) = 0.0D0
    10 END DO

!       Obtain F & FDOT due to NNB members of JLIST.
    DO 50 L = 1,NNB
        J = JLIST(L)
        IF (J == I) GO TO 50
    
        RIJ2 = 0.0D0
        RIJDOT = 0.0D0
    !          CALL JPRED(J,TIME,TIME)
        DO 20 K = 1,3
            A(K) = X(K,J) - X(K,I)
            A(K+3) = XDOT(K,J) - XDOT(K,I)
            RIJ2 = RIJ2 + A(K)**2
            RIJDOT = RIJDOT + A(K)*A(K+3)
        20 END DO
        A(8) = BODY(J)/(RIJ2*SQRT(RIJ2))
        A(9) = 3.0*RIJDOT/RIJ2
    
    !       Set approximate F & FDOT to be used by body #I in FPOLY2.
        DO 30 K = 1,3
            F(K,I) = F(K,I) + A(K)*A(8)
            FDOT(K,I) = FDOT(K,I) + (A(K+3) - A(K)*A(9))*A(8)
        30 END DO
    50 END DO

!       Initialize time of last force (prevents prediction in FPOLY2).
    T0(I) = TIME

    RETURN

    END SUBROUTINE FCLOSE
