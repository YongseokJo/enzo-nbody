    SUBROUTINE BINDAT


!       Binary data bank.
!       -----------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
    COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX), &
    HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX), &
    NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
    REAL*8 ::  EB(KMAX),ECC(KMAX),RCM(KMAX),ECM(KMAX),PB(KMAX),AS(30)
    REAL*8 ::  XX(3,3),VV(3,3)
     
!      edited by sykim for outputs
    CHARACTER(27) :: OUTFILE1
    CHARACTER(29) :: OUTFILE2
    CHARACTER(20) :: OUTFORM1
    CHARACTER(20) :: OUTFORM2
    CHARACTER(20) :: TCHAR

    DATA BTCOUNT /1/ ! by sykim
    INTEGER :: BTCNUM ! by sykim
    SAVE BTCOUNT ! by sykim

!       added by sykim
!       write down appropriate time names for outputs

    write (6,*) 'file name output starts for binary - SY'

    BTCNUM = INT(LOG10(REAL(BTCOUNT))) + 1

    write(OUTFORM1,"(A5,I1,A1)") '(A7,I',BTCNUM,')'
    write(OUTFORM2,"(A5,I1,A1)") '(A9,I',BTCNUM,')'

    write(OUTFILE1,OUTFORM1) 'bdat.9_',INT(BTCOUNT)
    write(OUTFILE2,OUTFORM2) 'bwdat.19_',INT(BTCOUNT)

    write (6,*) 'file name output ends for binary - SY'

!       end added by sykim

!       Form binding energy and central distance for each KS pair.
    ZMBIN = 0.0
    DO 10 JPAIR = 1,NPAIRS
        J2 = 2*JPAIR
        J1 = J2 - 1
        ICM = N + JPAIR
        ZMBIN = ZMBIN + BODY(ICM)
        BODYCM = BODY(ICM)
    !       Determine merger & ghost index for negative c.m. name (skip ghost).
        IF (NAME(ICM) < 0 .AND. BODY(ICM) > 0.0) THEN
            CALL FINDJ(J1,J,IMERGE)
        !       Employ actual masses and two-body distance for energy & eccentricity.
            BODYCM = CM(1,IMERGE) + CM(2,IMERGE)
            EB(JPAIR) = CM(1,IMERGE)*CM(2,IMERGE)*HM(IMERGE)/BODYCM
            SEMI = -0.5*BODYCM/HM(IMERGE)
            RJ = SQRT(XREL(1,IMERGE)**2 + XREL(2,IMERGE)**2 + &
            XREL(3,IMERGE)**2)
        !       Assume that merged binary is near apo or peri (hence ignore TDOT2).
            ECC2 = (1.0 - RJ/SEMI)**2
        !       Include separate diagnostics for the hierarchy (inner comps J1 & J).
            SEMI1 = -0.5*BODY(ICM)/H(JPAIR)
            ECC1 = (1.0 - R(JPAIR)/SEMI1)**2 + &
            TDOT2(JPAIR)**2/(BODY(ICM)*SEMI1)
            E0 = SQRT(ECC2)
            E1 = SQRT(ECC1)
            IF (J < 0) J = J1
            RM = SEMI*(1.0 - E0)/MAX(RADIUS(J1),RADIUS(J),1.0D-20)
        !              RM = MIN(RM,99.9D0)
            P0 = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
            P1 = DAYS*SEMI1*SQRT(ABS(SEMI1)/BODY(ICM))
            DO 2 K = 1,3
                XX(K,1) = XREL(K,IMERGE)
                XX(K,2) = 0.0
                XX(K,3) = X(K,J2)
                VV(K,1) = VREL(K,IMERGE)
                VV(K,2) = 0.0
                VV(K,3) = XDOT(K,J2)
            2 END DO
            CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ALPH)
            PCR = stability(CM(1,IMERGE),CM(2,IMERGE),BODY(ICM),E0,E1, &
            ALPH)*SEMI
            PM = SEMI1*(1.0 - E1)/PCR
        !              if(rank.eq.0)
        !     &        WRITE (84,3) TTOT, NAME(J1), NAME(J), KSTAR(J1), KSTAR(J),
        !     &                     KSTARM(IMERGE), E0, E1, PM, RM, P0, P1, SEMI1
        !    3         FORMAT ('BINDAT: Time[NB] NAME[I1] NAME[I3] K*[I1]',
        !     &             'K*[I3] K*[IM] ECC0 ECC1 PERI(I3)/PCR ',
        !     &             'PERI(INCM)[RSM] P0[days] P1[days] SEMI1[NB] ',
        !     &             1P,E12.5,0P,2I12,3I4,2F7.3,2F12.5,1P,3E14.5,0P)
        !              CALL FLUSH(84)
        ELSE IF (BODY(J1) > 0.0D0) THEN
        !       Form binding energy and eccentricity for standard case.
            EB(JPAIR) = BODY(J1)*BODY(J2)*H(JPAIR)/ &
            (BODY(J1) + BODY(J2))
            SEMI = -0.5*BODY(ICM)/H(JPAIR)
            ECC2 = (1.0 - R(JPAIR)/SEMI)**2 + &
            TDOT2(JPAIR)**2/(BODY(ICM)*SEMI)
        ELSE
            IM = 0
        !       Search merger table to identify corresponding index of c.m. name.
            DO 5 K = 1,NMERGE
                IF (NAMEG(K) == NAME(ICM)) THEN
                    IM = K
                END IF
            5 END DO
            IF (IM == 0) GO TO 10
            BODYJ1 = CM(3,IM)
            BODYJ2 = CM(4,IM)
            BODYCM = BODYJ1 + BODYJ2
            BODYCM = MAX(BODYCM,1.0D-10)
            EB(JPAIR) = BODYJ1*BODYJ2*H(JPAIR)/BODYCM
            SEMI = -0.5*BODYCM/H(JPAIR)
            ECC2 = (1.0 - SEMI/R(JPAIR))**2
        END IF
        ECC(JPAIR) = SQRT(ECC2)
    !          EB(JPAIR) = MAX(EB(JPAIR),-9.99999)
        PB(JPAIR) = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
    !          PB(JPAIR) = MIN(PB(JPAIR),99999.9)
        IF (SEMI < 0.0) PB(JPAIR) = 0.0
    !       Obtain binding energy (per unit mass) of c.m. motion.
        VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2
        IF (BODY(ICM) == 0.0D0) VJ2 = 0.0
    !$$$          POTJ = 0.0
    !$$$          DO 9 J = IFIRST,NTOT
    !$$$              IF (J.EQ.ICM) GO TO 9
    !$$$              RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
    !$$$     &                                        (X(3,ICM) - X(3,J))**2
    !$$$              POTJ = POTJ + BODY(J)/SQRT(RIJ2)
    !$$$    9     CONTINUE
        ECM(JPAIR) = 0.5*VJ2 - phidbl(ICM)
    !       Check for external tidal field (note that HT includes mass).
        IF (KZ(14) > 0) THEN
            CALL XTRNLV(ICM,ICM)
            ECM(JPAIR) = ECM(JPAIR) + HT/(BODY(ICM) + 1.0E-20)
        END IF
        RCM(JPAIR) = SQRT((X(1,ICM) - RDENS(1))**2 + &
        (X(2,ICM) - RDENS(2))**2 + &
        (X(3,ICM) - RDENS(3))**2)
        RCM(JPAIR) = MIN(RCM(JPAIR),99.9)
    10 END DO

!       Copy relevant binary diagnostics to single precision.
    AS(1) = TIME + TOFF
    AS(2) = RSCALE
    AS(3) = RTIDE
    AS(4) = RC
    AS(5) = TPHYS
    AS(6) = -1.5*(TIDAL(1)*ZMASS**2)**0.3333
    AS(7) = 0.0
    DO 20 K = 1,10
        AS(K+7) = E(K)
    20 END DO
    AS(18) = SBCOLL
    AS(19) = BBCOLL
    AS(20) = ZKIN
    AS(21) = POT
    AS(22) = EBIN0
    AS(23) = EBIN
    AS(24) = ESUB
    AS(25) = EMERGE
    AS(26) = BE(3)
    AS(27) = ZMASS
    AS(28) = ZMBIN
    AS(29) = CHCOLL
    AS(30) = ECOLL

!       Write formatted data bank on unit 9.

    if(rank == 0)then
    !     Split the bdat.9 by time
    !         call string_left(TCHAR,TTOT,DELTAT)
    !         write(OUTFILE,118) TCHAR
    ! 118     format('bdat.9_',A20)
        write(6,*) 'bdat starts - sykim'
        OPEN (UNIT=9,STATUS='UNKNOWN',FORM='FORMATTED',FILE=OUTFILE1)
                 
        WRITE (9,30)  NPAIRS, MODEL, NRUN, N, NC, NMERGE, (AS(K),K=1,7)
        30 FORMAT (3I4,I6,2I4,2X,F7.1,2F7.2,F7.3,F8.1,2F9.4)
        WRITE (9,35)  (AS(K),K=8,17)
        35 FORMAT (10F11.6)
        WRITE (9,40)  (AS(K),K=18,30)
        40 FORMAT (13F10.5)
        WRITE (9,*) 'NAME(I1)    NAME(I2)    ', &
        'M1[M*]                    M2[M*]                    ', &
        'EB[NB]                    ECC                       ', &
        'P[Days]                   SEMI[AU]                  ', &
        'RI[PC]                    VI[km/s]                  ', &
        'K*(I1)      K*(I2)      ', &
        'ZN[NB]                    RP[NB]                    ', &
        'STEP(I1)[NB]              NAME(ICM)                 ', &
        'ECM[NB]                   K*(ICM)      '
    end if

    DO 50 JPAIR = 1,NPAIRS
        J1 = 2*JPAIR - 1
        J2 = 2*JPAIR
        KCM = KSTAR(N+JPAIR)
        IF (NAME(N+JPAIR) < 0) THEN
            KCM = -10
        END IF
        ICM = N + JPAIR
        SEMI = -0.5*(BODY(J1) + BODY(J2))/H(JPAIR)
        ZN = SQRT((BODY(J1) + BODY(J2))/SEMI**3)
        RP = U(1,JPAIR)**2 + U(2,JPAIR)**2 + U(3,JPAIR)**2 + &
        U(4,JPAIR)**2
        RI = SQRT((X(1,ICM) - RDENS(1))**2 + &
        (X(2,ICM) - RDENS(2))**2 + &
        (X(3,ICM) - RDENS(3))**2)
        VI = SQRT(XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2)
        if(rank == 0) &
        WRITE (9,*)  NAME(J1), NAME(J2), BODY(J1)*ZMBAR, &
        BODY(J2)*ZMBAR, EB(JPAIR), ECC(JPAIR), PB(JPAIR), &
        SEMI*RAU, RI*RBAR, VI*VSTAR, KSTAR(J1), KSTAR(J2), &
        ZN, RP, STEP(J1), NAME(N+JPAIR), ECM(JPAIR), KCM
        45 FORMAT (2I8,1P,8E13.5,0P,3I8,3I4)
    50 END DO
    CALL FLUSH(9)
    CLOSE(9)

!       Include optional table of wide binaries on fort.19.

    if(rank == 0) then
    !     Split the bwdat.9 by time
    !         write(OUTFILE2,119) TCHAR
    ! 119     format('bwdat.19_',A20)
        write(6,*) 'bwdat starts - sykim'
        OPEN (UNIT=19,STATUS='UNKNOWN',FORM='FORMATTED',FILE=OUTFILE2)


        WRITE (19,55)  TIME+TOFF, (TIME+TOFF)*TSTAR, N
        55 FORMAT(' WIDE PAIRS    T TPHYS N ',1P,E27.16,E27.16,0P,I12)
        WRITE (19,*) 'NAME(I1)    NAME(I2)    ', &
        'M(I1)[M*]                 M(I2)[M*]                 ', &
        'EB[NB]                    ECC                       ', &
        'P[Days]                   SEMI[AU]                  ', &
        'RI[PC]                    VI[km/s]                  ', &
        'K*(I1)      K*(I2)      '
    end if
!       Adopt a generous criterion for semi-major axis of wide binaries.
    RB1 = 0.1*RSCALE
    RB2 = RB1**2
    NEWI = 0
    DO 80 I = IFIRST,NTOT
        NNB = LIST(1,I)
        RCL2 = RB2
        JMIN = I
    !       Determine the closest particle.
        DO 65 L = 2,NNB+1
            J = LIST(L,I)
        !       Include fast skips on each dimension to reduce the effort.
            IF (ABS(X(1,I) - X(1,J)) > RB1) GO TO 65
            IF (ABS(X(2,I) - X(2,J)) > RB1) GO TO 65
            IF (ABS(X(3,I) - X(3,J)) > RB1) GO TO 65
            RIJ2 = 0.0
            DO 60 K = 1,3
                RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
            60 END DO
            IF (RIJ2 < RCL2) THEN
                JMIN = J
                RCL2 = RIJ2
            END IF
        65 END DO
        IF (RCL2 >= RB2) GO TO 80
        VREL2 = 0.0
        RDOT = 0.0
        DO 70 K = 1,3
            VREL2 = VREL2 + (XDOT(K,I) - XDOT(K,JMIN))**2
            RDOT = RDOT + (X(K,I)-X(K,JMIN))*(XDOT(K,I)-XDOT(K,JMIN))
        70 END DO
        RIJ = SQRT(RCL2)
        ZMB = BODY(I) + BODY(JMIN)
        SEMI = 2.0/RIJ - VREL2/ZMB
        SEMI = 1.0/SEMI
        IF (SEMI > 0.0 .AND. SEMI < RB1) THEN
        !     Exclude duplicates by examining current list of NEWI components.
            DO 72 L = 1,NEWI
                IF (I == JLIST(L) .OR. JMIN == JLIST(L)) GO TO 80
            72 END DO
            JLIST(NEWI+1) = I
            JLIST(NEWI+2) = JMIN
            NEWI = NEWI + 2
            ECC2 = (1.0 - RIJ/SEMI)**2 + RDOT**2/(ZMB*SEMI)
            ECC1 = SQRT(ECC2)
            HB = -0.5D0*BODY(I)*BODY(JMIN)/SEMI
                        
            TK = DAYS*SEMI*SQRT(SEMI/ZMB)
            RI = 0.0
            VI = 0.0
            DO K = 1,3
                XCM = (BODY(I)*X(K,I) + BODY(JMIN)*X(K,JMIN))/ZMB &
                - RDENS(K)
                VCM = (BODY(I)*XDOT(K,I) + BODY(JMIN)*XDOT(K,JMIN))/ZMB
                RI = RI + XCM**2
                VI = VI + VCM**2
            END DO
            RI = SQRT(RI)
            VI = SQRT(VI)
        !     Print basic binary parameters (SEMI in AU, period in years).
            if(rank == 0) then
                WRITE (19,*)  NAME(I), NAME(JMIN), BODY(I)*ZMBAR, &
                BODY(JMIN)*ZMBAR, HB, ECC1, TK, SEMI*RAU, RI*RBAR, &
                VI*VSTAR, KSTAR(I), KSTAR(JMIN)
            end if
            75 FORMAT (F8.3,F9.1,1P,E9.1,0P,2F6.1,2I7,2I4)
        END IF
    80 END DO
    CALL FLUSH(19)
    CLOSE(19)

    BTCOUNT = BTCOUNT + 1
     
    RETURN

    END SUBROUTINE BINDAT
