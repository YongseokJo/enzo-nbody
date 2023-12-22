    SUBROUTINE FLYBY(I,ITERM)


!      Flyby check of KS orbit.
!       -----------------------

    USE POINTERS
    INCLUDE 'common6.h'
          


!       Set index of KS pair & first component of c.m. body #I.
    CALL JPRED(I,TIME,TIME)
    IPAIR = I - N
    I1 = 2*IPAIR - 1
    RJMIN2 = 1000.0
    ITERM = 0
    JCL = LIST(2,I1)

!       Find the closest body (JCL) and copy perturbers to JLIST.
    NNB = LIST(1,I1)
    DO 10 L = 2,NNB+1
        J = LIST(L,I1)
        JLIST(L-1) = J
        CALL JPRED(J,TIME,TIME)
        RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 + &
        (X(3,I) - X(3,J))**2
        IF (RIJ2 < RJMIN2) THEN
            RJMIN2 = RIJ2
            JCL = J
        END IF
    10 END DO

!       Skip rare case of merged binary or chain c.m. (denoted by NAME <= 0).
    IF (NAME(I) <= 0 .OR. NAME(JCL) <= 0) GO TO 20

    RDOT = (X(1,I) - X(1,JCL))*(XDOT(1,I) - XDOT(1,JCL)) + &
    (X(2,I) - X(2,JCL))*(XDOT(2,I) - XDOT(2,JCL)) + &
    (X(3,I) - X(3,JCL))*(XDOT(3,I) - XDOT(3,JCL))

    VREL2 = (XDOT(1,I) - XDOT(1,JCL))**2 + &
    (XDOT(2,I) - XDOT(2,JCL))**2 + &
    (XDOT(3,I) - XDOT(3,JCL))**2
     
!       Evaluate outer semi-major axis, eccentricity & impact parameter.
    RIJ = SQRT(RJMIN2)
    SEMI1 = 2.0/RIJ - VREL2/(BODY(I) + BODY(JCL))
    SEMI1 = 1.0/SEMI1
    ECC1 = SQRT((1.0D0 - RIJ/SEMI1)**2 + &
    RDOT**2/(SEMI1*(BODY(I) + BODY(JCL))))
    PMIN = SEMI1*(1.0D0 - ECC1)

!       Set semi-major axis & eccentricity of inner binary.
    SEMI = -0.5D0*BODY(I)/H(IPAIR)
    ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
    ECC = SQRT(ECC2)

    I2 = I1 + 1
    RJI1 = (X(1,JCL) - X(1,I1))**2 + (X(2,JCL) - X(2,I1))**2 + &
    (X(3,JCL) - X(3,I1))**2
    RJI2 = (X(1,JCL) - X(1,I2))**2 + (X(2,JCL) - X(2,I2))**2 + &
    (X(3,JCL) - X(3,I2))**2

!       Identify the dominant interaction (#I1 or #I2).
    FJI1 = (BODY(I1) + BODY(JCL))/RJI1
    FJI2 = (BODY(I2) + BODY(JCL))/RJI2
    NNB = NNB + 1
    IF (FJI1 > FJI2) THEN
        J = I1
        RJ2 = RJI1
    !       Add the other component to the perturber list.
        JLIST(NNB) = I2
    ELSE
        J = I2
        RJ2 = RJI2
        JLIST(NNB) = I1
    END IF

!       Evaluate radial velocity (routine SEARCH requires RD < 0).
    RD = 0.0
    RJJ = 0.0
    VJJ = 0.0
    DO 15 K = 1,3
        RD = RD + (X(K,JCL) - X(K,J))*(XDOT(K,JCL) - XDOT(K,J))
        RJJ = RJJ + (X(K,JCL) - X(K,J))**2
        VJJ = VJJ + (XDOT(K,JCL) - XDOT(K,J))**2
    15 END DO

!       Determine vectorial perturbation on intruder and closest component.
    CALL FPERT(JCL,J,NNB,PERT)
    GJ = PERT*RJ2/(BODY(JCL) + BODY(J))
    RJ = SQRT(RJJ)
!     if(rank.eq.0)
!    &WRITE (6,18) NAME(J),NAME(JCL),GAMMA(IPAIR),GJ,R(IPAIR),RJ,RD,
!    &             RDOT
!  18 FORMAT (' NMJ NMJC GI GJ R RJ RD RDI ',2I5,2F6.2,1P,4E9.1)

!       Terminate if perturber & I1/I2 provides dominant force.
    APO = ABS(SEMI)*(1.0 + ECC)
    RSUM = R(IPAIR) + SQRT(RJ2)
    XF = 2.0
!       Compare current radial velocity with #JCL & J.
    IF (TDOT2(IPAIR)/R(IPAIR) < RD/RJ) XF = 1.0

!       Increase the cross section for quadruples.
    IF (JCL > N) THEN
        SEMI2 = -0.5*BODY(JCL)/H(JCL-N)
        APO = APO + ABS(SEMI2)
    END IF

    SEMIJ = 2.0/RJ - VJJ/(BODY(JCL) + BODY(J))
    SEMIJ = 1.0/SEMIJ
!       Check for chain regularization test or standard termination.
    IF (PMIN < APO .AND. RSUM < XF*RMIN .AND. RDOT < 0.0) THEN
        IF (SEMI < 0.0 .AND. GAMMA(IPAIR) < 0.5) THEN
            ITERM = 0
            GO TO 20
        END IF
        ITERM = 1
    ELSE IF (GJ < 0.7*GAMMA(IPAIR) .AND. RD <= 0.0) THEN
        ITERM = 2
    END IF

    20 RETURN

    END SUBROUTINE FLYBY
