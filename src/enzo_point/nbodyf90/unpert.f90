    SUBROUTINE UNPERT(IPAIR)


!       Unperturbed two-body motion.
!       ----------------------------

    Include 'kspars.h'
    USE POINTERS
    INCLUDE 'common6.h'
          
    include 'timing.h'
    REAL*8 ::  UI(4),UIDOT(4)
!     Warning!: saved value is unsafe for parallel
    SAVE TCALL
    DATA TCALL /0.0D0/


!       Set first component & c.m. index and semi-major axis.
    I1 = 2*IPAIR - 1
    I = N + IPAIR
    SEMI = -0.5D0*BODY(I)/H(IPAIR)

!       Add elapsed unperturbed Kepler periods and update the time T0(I1).
    TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
    IF (TIME - T0(I1) < 2.0E+09*TK) THEN
        K = NINT((TIME - T0(I1))/TK)
    ELSE
        K = 0
        NPRECT = NPRECT + 1
    END IF
    T0(I1) = TIME

!       Reset unperturbed counter if > 2*10**9 and update the frequency.
    IF (NKSPER > 2000000000 .OR. NKSPER < 0) THEN
        NKSPER = 0
        NPRECT = NPRECT + 1
    END IF
    NKSPER = NKSPER + K

!       Include case of tidal dissipation or partial reflection (suppressed).
    IF (TDOT2(IPAIR) >= 0.0D0) THEN
    !         IF (KZ(25).GT.0.AND.ABS(TDOT2(IPAIR).GT.1.0E-10) THEN
    !             DT = 0.0
    !             GO TO 10
    !         END IF
    !       Ensure perturbation check at least once every c.m. step (KSTAR < 11).
        IF (KZ(27) > 0) THEN
            IF (KSTAR(I) < 11) THEN
                IF (STEP(I1) < TK) THEN
                    IF (STEP(I1) < 0.0001*STEP(I)) GO TO 9
                END IF
                IF (TIME - T0(I) > 2.0*STEP(I1)) THEN
                    DT = MIN(3.0D0*STEP(I1),STEP(I))
                    KPERT = 0
                    JCL = 0
                    GO TO 20
                END IF
            END IF
        END IF
    END IF

!       Evaluate interval for unperturbed motion (GAMMA < GMIN).
    9 KPERT = 1
    IP = IPAIR
!     Copy JCL from IPAIR and restore value
    call cputim(tt1)
    CALL TPERT(IPAIR,GMIN,DT)
    call cputim(tt2)
    tttp = tttp + (tt2-tt1)*60.
    JCL = IPAIR
    IPAIR = IP

!       Restore KS indicator and re-initialize if interval < period.
    IF (DT < TK) THEN
    !       Form perturber list and restart KS motion if required.
        CALL KSLIST(IPAIR)
        IF (LIST(1,I1) > 0) THEN
        !       Transform to apocentre variables in case of tidal dissipation.
            IF (R(IPAIR) < SEMI) THEN
                NP = LIST(1,I1)
                LIST(1,I1) = 0
            !       Do not allow backwards integration on switch from unperturbed motion.
                CALL KSPERI(IPAIR)
                LIST(1,I1) = NP
                CALL KSAPO(IPAIR)
            END IF
            KSLOW(IPAIR) = 1
            CALL RESOLV(IPAIR,1)
            CALL KSPOLY(IPAIR,1)
        !       Include differential correction (experimental 10/02).
            JLIST(1) = I
            NNB = LIST(1,I1)
            DO 12 L = 1,NNB
                JPERT(L) = LIST(L+1,I1)
            12 END DO
            CALL NBPOT(1,NNB,POT1)
            JLIST(1) = I1
            JLIST(2) = I1 + 1
            CALL NBPOT(2,NNB,POT2)
        !             if(rank.eq.0)
        !    &        WRITE (6,13)  NAME(I1), POT1-POT2, GAMMA(IPAIR)
        !  13         FORMAT ('  CORRECT    NAM DPHI G ',I6,1P,2E10.2)
            DETMP = POT2 - POT1
            EMERGE = EMERGE + DETMP
            BE(3) = BE(3) + DETMP
        !     --03/05/14 20:25-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$              print*,rank,'unpert',(POT2-POT1),'t',time
        !$$$              call flush(6)
        !     --03/05/14 20:25-lwang-end----------------------------------------*
        !     ks MPI communiction emerge
        !              call ksparmpi(K_store,K_real8,K_EMERGE,0,0,DETMP)
        !     ks MPI communiction BE(3)
        !              call ksparmpi(K_store,K_real8,K_BE3,0,0,DETMP)
            GO TO 30
        END IF
        IR = 1
        GO TO 28
    END IF

!       Check for tidal dissipation (peri & apocentre included; skip merger).
    20 IF (KZ(27) > 0 .AND. &
    (SEMI < RMIN .OR. TDOT2(IPAIR) > 0.0D0) .AND. &
    (KSTAR(I) == 0 .OR. KSTAR(I) == -1) .AND. &
    (NAME(I) > 0)) THEN
    
    !       Compare pericentre and effective capture distance.
        ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + &
        TDOT2(IPAIR)**2/(BODY(I)*SEMI)
        ECC = SQRT(ECC2)
        RP = SEMI*(1.0D0 - ECC)
        RT = 4.0*MAX(RADIUS(I1),RADIUS(I1+1))
    
    !       Specify circularization index (skip SPIRAL but include CHAOS).
        ICIRC = 0
        IF (KZ(27) == 1 .AND. TTOT > TCALL) THEN
            IF (RP < RT) ICIRC = 1
            TCALL = TTOT + 0.01
        ELSE IF (RP < 2.5*RT .AND. KSTAR(I) == 0 .AND. &
            TTOT > TCALL) THEN
        !       Perform occasional checks of TC (NB! small periods).
            CALL TCIRC(RP,ECC,I1,I1+1,ICIRC,TC)
            TCALL = TTOT + 0.01
        ELSE IF (KSTAR(I) == -1) THEN
            ICIRC = 1
        END IF
    
    !       See whether dissipation is active (A*(1 - E) < 4*MAX(R) & E > 0.002).
    !         IF (RP.LT.0.99*RT.AND.ECC.GT.0.002) THEN
        IF (ICIRC > 0 .AND. ECC > 0.002) THEN
        !       Transform to pericentre if R > A (KSAPO works both ways with PI/2).
            IF (R(IPAIR) > SEMI) THEN
                CALL KSAPO(IPAIR)
            END IF
        !       Implement energy loss at pericentre (exit on collision).
            CALL KSTIDE(IPAIR,RP)
            IF (IPHASE < 0) GO TO 30
            SEMI = -0.5D0*BODY(I)/H(IPAIR)
            TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
        !       Assign one Kepler period for active state and impose TDOT2 > 0.
            IF (R(IPAIR) < 0.99*RT .AND. KSTAR(I) /= 10) THEN
                STEP(I1) = TK
                TDOT2(IPAIR) = 1.0E-20
            ELSE
            !       Define apocentre and set 1/2 period for inactive or synchronous case.
                CALL KSAPO(IPAIR)
                STEP(I1) = 0.5*TK
            !       Note next unperturbed check at apocentre since T'' < 0 in KSAPO.
                IF (ABS(RT - RP)/RT > 0.1 .AND. KZ(27) <= 1) THEN
                    if(rank == 0) &
                    WRITE(6,25)  ECC, SEMI, R(IPAIR), RP, RADIUS(I1)
                    25 FORMAT (' INACTIVE PHASE    E A R RP R* ', &
                    F7.3,1P,4E10.2)
                END IF
                IF (ECC < 0.002 .AND. SEMI < 0.01*RMIN) THEN
                    KSTAR(I) = 10
                    TEV(I) = TIME
                !     ks MPI communicatin: KSTAR
                    call ksparmpi(K_store,K_int,K_KSTAR,I,0,KSTAR(I))
                !     ks MPI communicatin: TEV
                    call ksparmpi(K_store,K_real8,K_TEV,I,0,TEV(I))
                END IF
            END IF
            GO TO 30
        END IF
    END IF

!       Check for Roche overflow from synchronous orbit (but KSTAR < 13).
    IR = 0
    IF (KSTAR(I) == 13 .AND. MIN(TEV(I1),TEV(I1+1)) < TIME + DT) THEN
        IR = -1
    END IF
    IF (KSTAR(I) > 13) IR = 1
    IF (KZ(34) > 0 .AND. (KSTAR(I) >= 10 .AND. KSTAR(I) <= 12)) THEN
        TM = MIN(TEV(I1),TEV(I1+1),TEV(I))
        IF (TM < TIME) THEN
            CALL TRFLOW(IPAIR,DTR)
        !       Exit if ROCHE is indicated (change 16/08/2006).
        !             IF (DTR.LT.STEP(I1)) THEN
        !                 CALL ROCHE(IPAIR)
        !       Exit on coalescence (KPERT > 0 would cause trouble).
        !                 IF (IPHASE.LT.0) GO TO 30
        !             ELSE
        !                 IR = 1
        !             END IF
            IR = 1
            IF (DTR < STEP(I1)) GO TO 30
        
            IF (DTR > TK) THEN
                DT = MIN(DTR,DT)
                KPERT = 2
            END IF
        ELSE
            IR = 1
        END IF
    END IF

!       Perform general two-body collision test.
    IF (KZ(19) >= 3 .AND. NAME(I) > 0) THEN
        RI = 0.0
        DO 26 K = 1,4
            UI(K) = U0(K,IPAIR)
            UIDOT(K) = UDOT(K,IPAIR)
            RI = RI + UI(K)**2
        26 END DO
        CALL PERI(UI,UIDOT,RI,BODY(I1),BODY(I1+1),QPERI)
        IF (QPERI < 0.75*(RADIUS(I1) + RADIUS(I1+1))) THEN
        !$$$*       Obtain KS variables at pericentre before coalescence to one body.
        !$$$              CALL KSPERI(IPAIR)
        !$$$              KSPAIR = IPAIR
        !$$$*       Set indicator for skipping ECOLL updating in COAL.
        !$$$              IQCOLL = -2
        !$$$              CALL CMBODY(2)
        !$$$              IF (IPHASE.LT.0) GO TO 30
            IPHASE = -1
            GO TO 30
        END IF
    END IF

!       Specify a conservative number of unperturbed orbits from TPERT.
    IF (KPERT > 0 .AND. DT > 0.0) THEN
    !       Ensure that next look-up time is not exceeded (option #27).
        IF (KZ(27) > 0 .AND. KPERT > 1) THEN
            TM = MIN(TEV(I1),TEV(I1+1))
            IF (TM - T0(I1) > 0.0 .AND. TM - T0(I1) < 0.5*DT) THEN
                NWARN = NWARN + 1
            !                  IF (rank.eq.0.and.NWARN.LT.1000) THEN
            !                      WRITE (25,27)  IPAIR, KSTAR(I1), KSTAR(I1+1),
            !     &                               KSTAR(I), TM-T0(I1), STEP(I1)
            !   27                 FORMAT (' UNPERT WARNING!   IPAIR  K*(I1) K*(I2)',
            !     &                     ' K*(ICM) TEVNXT[NB] STEP(I1)[NB]',
            !     &                     I12,3I4,1P,2E26.17,0P)
            !                      CALL FLUSH(25)
            !                  END IF
                DT = MIN(2.0*(TM - T0(I1)),DT)
                DT = MAX(DT,TK,STEPX)
            END IF
        END IF
    !       Adopt c.m. step instead if integer argument exceeds 10**9.
        IF (DT < 2.0E+09*TK) THEN
            K = 1 + INT(0.5D0*DT/TK)
        !       Restrict Kepler period to c.m. step (case of very wide orbit).
            STEP(I1) = FLOAT(K)*MIN(TK,STEP(I))
        ELSE
            STEP(I1) = STEP(I)
        END IF
    !       Include optional treatment for spiralling of chaotic binary orbit.
        IF (KZ(27) > 1 .AND. KSTAR(I) == -2) THEN
        !       Ensure pericentre position after possible perturbed motion.
            IF (R(IPAIR) > SEMI) THEN
                CALL KSRECT(IPAIR)
            !       Reduce eccentric anomaly by pi for inward motion.
                IF (TDOT2(IPAIR) < 0.0D0) THEN
                    CALL KSAPO(IPAIR)
                END IF
            !       Transform from outward motion (anomaly < pi) to exact pericentre.
                CALL KSPERI(IPAIR)
            END IF
            CALL SPIRAL(IPAIR)
            IF (IPHASE < 0) GO TO 30
        END IF
    END IF

!       Check merger condition before continuing (skip small Roche steps).
    28 IF (KZ(15) > 0 .AND. IR >= 0 .AND. TIME+STEP(I1) > TBLOCK) THEN
        IF (STEP(I) < DTMIN) THEN
            CALL IMPACT(I)
        ELSE IF (JCL > 0 .AND. STEP(I) < 10.0*DTMIN) THEN
            CALL HISTAB(IPAIR,JCL,PMIN,RSTAB)
            IF (RSTAB < PMIN) THEN
                CALL IMPACT(I)
            END IF
        !       Include a more generous test for massive quadruples.
        ELSE IF (JCL > N) THEN
            FAC = 2.0*(BODY(I) + BODY(JCL))/BODYM
            IF (STEP(I) < FAC*DTMIN) THEN
                CALL HISTAB(IPAIR,JCL,PMIN,RSTAB)
                IF (RSTAB < PMIN) THEN
                    CALL IMPACT(I)
                END IF
            END IF
        END IF
    END IF

!       Check collision criterion for special case.
    IF (KZ(27) == -1) THEN
        ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + &
        TDOT2(IPAIR)**2/(BODY(I)*SEMI)
        ECC = SQRT(ECC2)
        QPERI = SEMI*(1.0 - ECC)
        RFAC = 2.0
        I2 = I1 + 1
        IF (QPERI < RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
            J1 = I1
            IF (RADIUS(I2) > RADIUS(I1)) J1 = I2
            FAC = 0.5*BODY(I)/BODY(J1)
        !       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
            RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
            IF (QPERI < RCOLL) THEN
                CALL TOUCH(IPAIR,I1,I2,RCOLL)
            END IF
        END IF
    END IF

    30 RETURN

    END SUBROUTINE UNPERT
