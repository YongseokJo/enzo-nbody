    SUBROUTINE FCORR(I,DM,KW)


!       Total force corrections due to masss loss.
!       ------------------------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    INCLUDE 'omp_lib.h'
    REAL*8 ::  A(6),VOLD(3)
    LOGICAL :: IKICK
!      REAL*8  XBACK(3,LMAX),XDBACK(3,LMAX)


!       Save the velocity components and square velocity.
    call xbpredall
!      call jpred(I,TIME,TIME)
    VI2 = 0.0
    DO 1 K = 1,3
        VOLD(K) = XDOT(K,I)
        VI2 = VI2 + XDOT(K,I)**2
    1 END DO

!       Include velocity kick in case of new WD, NS, BH or massless SN.
    IF (KW >= 10 .AND. KW <= 15) THEN
        IKICK = .TRUE.
    !       Distinguish between single star (first time only) and binary.
        IF (I <= N .AND. KW /= KSTAR(I)) THEN
            CALL KICK(I,1,KW,DM)
        ELSE IF (I > N) THEN
            IPAIR = I - N
            CALL KICK(IPAIR,0,KW,DM)
        END IF
    ELSE
        IKICK = .FALSE.
    END IF

!       Define consistent c.m. variables for KS mass loss (exclude Roche).
    VFAC = 0.0
    IF (I > N .AND. KSTAR(I) <= 10) THEN
        I2 = 2*(I - N)
        I1 = I2 - 1
        VF2 = 0.0
        DV2 = 0.0
        BODYI = BODY(I)
        IF (BODY(I) == 0.0D0) BODYI = BODY(I1) + BODY(I2)
        DO 8 K = 1,3
            X(K,I) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/BODYI
            XDOT(K,I) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/ &
            BODYI
            X0(K,I) = X(K,I)
            X0DOT(K,I) = XDOT(K,I)
            VF2 = VF2 + XDOT(K,I)**2
            DV2 = DV2 + (XDOT(K,I) - VOLD(K))**2
        8 END DO
        VFAC = SQRT(VF2/VI2)
    END IF

!       Correct potential energy, all forces & first derivatives.
    POTJ = 0.0D0
!$$$*       Backup all predicted x and xdot to xback and xdback
!$$$      NNB = LIST(1,I)
!$$$      DO II = 1, NNB
!$$$         IL = LIST(II+1,I)
!$$$         CALL JPRED(IL,TIME,TIME)
!$$$         XBACK(1,II) = X0(1,IL)
!$$$         XBACK(2,II) = X0(2,IL)
!$$$         XBACK(3,II) = X0(3,IL)
!$$$         X0(1,IL) = X(1,IL)
!$$$         X0(2,IL) = X(2,IL)
!$$$         X0(3,IL) = X(3,IL)
!$$$         XDBACK(1,II) = X0DOT(1,IL)
!$$$         XDBACK(2,II) = X0DOT(2,IL)
!$$$         XDBACK(3,II) = X0DOT(3,IL)
!$$$         X0DOT(1,IL) = XDOT(1,IL)
!$$$         X0DOT(2,IL) = XDOT(2,IL)
!$$$         X0DOT(3,IL) = XDOT(3,IL)
!$$$      END DO

! omp parallel do
! omp& private(J,NNB1,RIJ2,RIJDOT,DX,DXDOT,RIJ,
! omp& RDVDOT,A3,A4,A5,A6,A7,A,K)
! omp& reduction(+:POTJ)
    DO 40 J = IFIRST,NTOT
        #ifdef SIMD
        IF (J == I) GO TO 35
        #else
        IF (J == I) GO TO 40
        #endif
        RIJ2 = 0.0D0
        RIJDOT = 0.0D0
        RDVDOT = 0.0D0
    
        DO 10 K = 1,3
            A(K) = X(K,I) - X(K,J)
            A(K+3) = VOLD(K) - XDOT(K,J)
            RIJ2 = RIJ2 + A(K)**2
            RIJDOT = RIJDOT + A(K)*A(K+3)
            RDVDOT = RDVDOT + A(K)*(XDOT(K,I) - VOLD(K))
        10 END DO
    
        RIJ = SQRT(RIJ2)
        POTJ = POTJ + BODY(J)/RIJ
        A3 = 1.0/(RIJ2*RIJ)
        A4 = BODY(I)*A3
        A5 = DM*A3
        A6 = 3.0*RIJDOT/RIJ2
        A7 = 3.0*RDVDOT/RIJ2
    
        DO 15 K = 1,3
            A(K+3) = (A(K+3) - A6*A(K))*A5
            IF (IKICK) THEN
            !       Include FDOT corrections due to increased velocity.
                A(K+3) = A(K+3) + (XDOT(K,I) - VOLD(K))*A4
                A(K+3) = A(K+3) - A7*A(K)*A4
            END IF
        15 END DO
    
    !       Use neighbour list of #J to distinguish irregular & regular terms.
        NNB1 = LIST(1,J) + 1
        DO 25 L = 2,NNB1
            IF (LIST(L,J) == I) THEN
                DO 20 K = 1,3
                    F(K,J) = F(K,J) - 0.5*A(K)*A5
                    FI(K,J) = FI(K,J) - A(K)*A5
                    FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
                    D1(K,J) = D1(K,J) - A(K+3)
                    FIDOT(K,J) = FIDOT(K,J) - A(K+3)
                20 END DO
            !$$$                  IF(IKICK.AND.T0(J)+0.5D0*STEP(J).GE.TBLOCK) THEN
            !$$$                     STEP(J) = 0.5D0*STEP(J)
            !$$$                     IF(T0(J)+0.5D0*STEP(J).GE.TBLOCK) THEN
            !$$$                        STEP(J) = 0.5D0*STEP(J)
            !$$$                     END IF
            !$$$                     TIMENW(J) = T0(J) + STEP(J)
            !$$$                  END IF

            !     --03/19/14 22:48-lwang-debug--------------------------------------*
            !**** Note:------------------------------------------------------------**
            !$$$          if(name(J).eq.7) then
            !$$$             print*,rank,'J',J,'F',F(1,J),'FD',FDOT(1,J),
            !$$$     &            'FI',FI(1,J),'FID',D1(1,J),'FC',0.5*A(1)*A5,'FDC',
            !$$$     &            -ONE6*A(4),
            !$$$     &            't',time,'XDOT',XDOT(1,J)
            !$$$             call flush(6)
            !$$$          end if
            !     --03/19/14 22:48-lwang-end----------------------------------------*
                #ifdef SIMD
                GO TO 35
                #else
                GO TO 40
                #endif
            END IF
        25 END DO
        DO 30 K = 1,3
            F(K,J) = F(K,J) - 0.5*A(K)*A5
            FR(K,J) = FR(K,J) - A(K)*A5
            FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
            D1R(K,J) = D1R(K,J) - A(K+3)
            FRDOT(K,J) = FRDOT(K,J) - A(K+3)
        30 END DO
        #ifdef SIMD
    !     Reset AVX/SSE particle data
        35 call IRR_SIMD_SET_JP(J,X0(1,J),X0DOT(1,J),F(1,J),FDOT(1,J), &
        BODY(J),T0(J))
        #endif
    40 END DO
! omp end parallel do

!$$$*       Revert X0 and X0dot to original values
!$$$      DO II = 1, NNB
!$$$         IL = LIST(II+1,I)
!$$$         X0(1,IL) = XBACK(1,II)
!$$$         X0(2,IL) = XBACK(2,II)
!$$$         X0(3,IL) = XBACK(3,II)
!$$$         X0DOT(1,IL) = XDBACK(1,II)
!$$$         X0DOT(2,IL) = XDBACK(2,II)
!$$$         X0DOT(3,IL) = XDBACK(3,II)
!$$$      END DO

!       Update the potential and kinetic energy loss.
    EMDOT = EMDOT - DM*POTJ + 0.5*DM*VI2

!       Modify energy loss further for c.m. body (exclude Roche cases).
    IF (I > N .AND. KSTAR(I) <= 10) THEN
        ECDOT = ECDOT - 0.5*BODY(I)*VI2*(VFAC**2 - 1.0)
    END IF

!       See whether linearized tidal terms should be included.
    IF (KZ(14) > 0 .AND. KZ(14) < 3) THEN
        EMDOT = EMDOT - 0.5*DM*(TIDAL(1)*X(1,I)**2 + &
        TIDAL(3)*X(3,I)**2)
    END IF

!       Check optional Plummer potential.
    IF (KZ(14) == 4 .OR. KZ(14) == 3) THEN
        RI2 = AP2
        DO 50 K = 1,3
            RI2 = RI2 + X(K,I)**2
        50 END DO
        EMDOT = EMDOT - DM*MP/SQRT(RI2)
    END IF

!       Accumulate energy loss for conservation check (not used).
    E(12) = EMDOT

    RETURN

    END SUBROUTINE FCORR
