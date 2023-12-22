    SUBROUTINE START


!       Initialization of data & polynomials.
!       ------------------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    INCLUDE 'timing.h'
    EXTERNAL SCALE,MERGE
!      COMMON/NNBCOT/ nnbave
    PARAMETER  (NS=12)


!       Initialize global scalars, counters & useful constants.
    CALL ZERO

!       Read input parameters.
    CALL INPUT

!       Open all Files.
    CALL FILE_INIT(0)

!       Set initial conditions: BODY(I), X(K,I), XDOT(K,I); I=1,N & K=1,3.
    CALL DATA

!       Scale initial conditions to new units.
    CALL SCALE

!       Set total mass in case routines DATA & SCALE are not used.
    ZMASS = 0.0D0
    DO 10 I = 1,N
        ZMASS = ZMASS + BODY(I)
    10 END DO

!       Define mean mass in scaled units and solar mass conversion factor.
    BODYM = ZMASS/FLOAT(N)
!     Notice from now ZMBAR is mass scale factor instead of average mass of input data
    IF (KZ(22) /= 10) THEN
        ZMBAR = ZMBAR/BODYM
    ELSE
        ZMBAR = ZMBAR*FLOAT(N)
    END IF

!       Introduce scaling factors DAYS, YRS, SU, RAU, SMU, TSTAR & VSTAR.
    CALL UNITS

!       Check option for external force.
    IF (KZ(14) > 0) THEN
        #ifdef TT
    !** FlorentR - initialize the tensors
        IF (KZ(14) == 9) CALL TTINIT
    !** FRenaud
        #endif
        CALL XTRNL0
    END IF

!       Check optional scaling to hot system.
    IF (KZ(29) > 0) THEN
        CALL HOTSYS
    END IF

!       Check option for initial binaries.
    IF (KZ(8) == 1 .OR. KZ(8) >= 3) THEN
        CALL BINPOP
    END IF

!       Include stable primordial triples.
    IF (KZ(18) > 1 .AND. KZ(8) > 0) THEN
        CALL HIPOP
    END IF

!       Initial massive black hole input
    IF (KZ(24) == 1) CALL IMBHINIT

!       Set sequential name, maximum mass & primary velocity, initial single/binary total mass
    BODY1 = 0.0
    ZSMASS0 = 0.0
    ZBMASS0 = 0.0
    DO 20 I = 1,N
    !     exclude massive black hole mass in initial cluster mass
        IF(KZ(24) == 1 .AND. I == N) GO TO 19
        IF(I <= 2*NBIN0) THEN
            ZBMASS0 = ZBMASS0 + BODY(I)
        ELSE
            ZSMASS0 = ZSMASS0 + BODY(I)
        END IF
        19 NAME(I) = I
        BODY1 = MAX(BODY1,BODY(I))
        DO 15 K = 1,3
            X0DOT(K,I) = XDOT(K,I)
        15 END DO
    20 END DO

!       Initialize fixed block steps (64 levels).
    CALL IBLOCK

!       Create table of inverse Stumpff coefficients.
    DO 30 I = 1,NS
        SCOEFF(I) = 1.0/((I + 1)*(I + 2))
    30 END DO

!       Set optional stellar evolution parameters.
    IF (KZ(19) > 2) THEN
        CALL INSTAR
    ELSE IF (KZ(14) > 1) THEN
        DT = 1.0E-03/TSCALE
        CALL STEPK(DT,DTN)
        STEPX = DTN
    END IF

!       Initialize optional cloud parameters.
    IF (KZ(13) > 0) THEN
        CALL CLOUD0
    END IF

!       Regularize any hard primordial binaries (assume sequential ordering).
    IF (KZ(8) > 0 .AND. NBIN0 > 0) THEN
        SMMIN = 1.D30
        SMMAX = 0.D0
        XMMIN = 1.D30
        XMMAX = 0.D0
        TMMIN = 1.D30
        TMMAX = 0.D0
    
    !     temperately set kstart = -1 for KSREG
        KSDUM = KSTART
        KSTART = -1
        DO 50 IPAIR = 1,NBIN0
            ICOMP = 2*IPAIR - 1
            JCOMP = 2*IPAIR
            RIJ2 = 0.0
        !       Include standard distance criterion.
            DO 45 K = 1,3
                RIJ2 = RIJ2 + (X(K,ICOMP) - X(K,JCOMP))**2
            45 END DO
            IF(RIJ2 < SMMIN) SMMIN=RIJ2
            IF(RIJ2 > SMMAX) SMMAX=RIJ2
            XMBIN = BODY(ICOMP) + BODY(JCOMP)
            PERIOD = RIJ2**0.75/DSQRT(XMBIN)
            IF(XMBIN < XMMIN) XMMIN=XMBIN
            IF(XMBIN > XMMAX) XMMAX=XMBIN
            IF(PERIOD < TMMIN) TMMIN=PERIOD
            IF(PERIOD > TMMAX) TMMAX=PERIOD
            IF (RIJ2 < RMIN**2) THEN
                CALL KSPREG
            END IF
        50 END DO
    !     recover kstart
        KSTART = KSDUM
    END IF
!     Diagnostic
    if(rank == 0) write(6,51) SQRT(SMMIN)*RAU, SQRT(SMMAX)*RAU, &
    XMMIN*ZMBAR, XMMAX*ZMBAR, TMMIN*DAYS, TMMAX*DAYS
    51 format(/,' Binary parameters: R12(min)[AU]',F15.6, &
    ' R12(max)[AU]',F15.6,' MCM(min)[M*]',F15.6, &
    ' MCM(max)[M*]',F15.6,' P[min][days]',1P,E15.6, &
    ' P[max][days]',E15.6,0P)

!       Set initial neighbour list & corresponding radius.
    RS0 = RC
    call cputim(tt1)
! ifdef GPU
    CALL FPOLY0(RS0)
    call cputim(tt2)
    if(rank == 0) write(6,52) (tt2-tt1)*60.0
    52 format(/,' Fpoly0 + Neighbor List Timing: ',F15.6)
    call flush(6)
! else
!      NNB0 = 0
!*
!      DO 40 I = IFIRST,NTOT
!          CALL NBLIST(I,RS0)
!          NNB0 = NNB0 + LIST(1,I)
!   40 CONTINUE
!*
!*       Obtain force & first derivative.
! ifdef PARALLEL
!      CALL FPOLY1_MPI(IFIRST,NTOT,0)
! else
!      CALL FPOLY1(IFIRST,NTOT,0)
! endif
!      call cputim(tt2)
!      if(rank.eq.0) write(6,53) (tt2-tt1)*60.0
! 53   format(/,' Fpoly1 + Neighbor List Timing: ',F15.6)
!      call flush(6)
!*
! endif

    IF (KZ(40) > 0) THEN
    !     Obtain second & third force derivatives and set time-steps.
        call cputim(tt1)
        #ifdef PARALLEL
        CALL FPOLY2_MPI(IFIRST,NTOT,0)
        #else
        CALL FPOLY2(IFIRST,NTOT,0)
        #endif
        call cputim(tt2)
        if(rank == 0) write(6,54) (tt2-tt1)*60.0
        54 format(/,' Fpoly2 Timing: ',F15.6)
        call flush(6)
    ELSE
    !     Obtain time step only based on xdot, F and Fdot
        CALL STEPS2(IFIRST,NTOT)
    END IF

!     Generate perturber list for KS
    CALL KSPINIT
    call cputim(tt3)
    if(rank == 0) write(6,55) (tt3-tt2)*60.0
    55 format(/,' KS Perturbers List Initialization Timing: ',F15.6)
    call flush(6)

!     --11/08/14 10:12-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      DO II = IFIRST,NTOT
!$$$         write(113,*) 'II',II,'N',NAME(II),'X',X(1:3,II),
!$$$     &        'XD',XDOT(1:3,II),'F',F(1:3,II),'FD',FDOT(1:3,II),
!$$$     &        'D2',D2(1:3,II),'STEP',STEP(II),'NNB',LIST(1,II),
!$$$     &        'LIST',LIST(2:5,II)
!$$$      END DO
!$$$      call flush(113)
!$$$      DO II = 1, NPAIRS
!$$$         write(114,*) 'II',II,'H',H(II),'R',R(II),'STEP',STEP(2*II-1),
!$$$     &        'NP',LIST(1,2*II-1),'LIST',LIST(2:5,2*II-1)
!$$$      END DO
!$$$      call flush(114)
!$$$      call abort()
!     --11/08/14 10:12-lwang-end----------------------------------------*
!       Include optional regularization of primordial triples.
    IF (KZ(18) > 1 .AND. NHI0 > 0) THEN
        KSPAIR = 1
    !       Note that each KS pair will move to the top of the queue.
        60 ICOMP = 2*KSPAIR - 1
        ICM = KSPAIR + N
        RX2 = 1.0
    !       Find index of closest outer component without any assumption.
        DO 70 J = IFIRST,N
            RIJ2 = 0.0
            DO 65 K = 1,3
                RIJ2 = RIJ2 + (X(K,ICM) - X(K,J))**2
            65 END DO
            IF (RIJ2 < RX2) THEN
                RX2 = RIJ2
                JCOMP = J
            END IF
        70 END DO
    !       Increase KS pointer for rejected separation (include possible exit).
        IF (SQRT(RX2) > RMIN) THEN
            KSPAIR = KSPAIR + 1
            IF (KSPAIR > NPAIRS) GO TO 80
            GO TO 60
        END IF
    !       Evaluate PCRIT for R0(NPAIRS) in MERGE since IMPACT is bypassed.
        CALL HISTAB(KSPAIR,JCOMP,PMIN,RSTAB)
    !       Initialize the triple (constructed to be stable in HIPOP).
        IPHASE = 6
        CALL MERGE
        IF (NMERGE < NHI0) THEN
            GO TO 60
        END IF
    END IF

!       Adjust NNBMAX (R.Sp.)
    80 NNBMAX = MIN(N/2,LMAX - 3)
    ZNBMIN = MAX(0.01*FLOAT(NNBMAX),1.0)
    ZNBMAX = 0.9*FLOAT(NNBMAX)


!       Check the average neighbour number.
    NNB0 = 0
    DO 86 I = IFIRST,NTOT
        NNB0 = NNB0 + LIST(1,I)
    86 END DO
    ZNB = FLOAT(NNB0)/FLOAT(N)
!      nnbave = ZNB
    IF (ZNB < 0.25*ZNBMAX) THEN
        if(rank == 0)WRITE (6,90)  ZNB
        90 FORMAT (/,12X,'WARNING!   SMALL NEIGHBOUR NUMBERS   <NNB> =', &
        F5.1)
    END IF

!       Check option for writing the initial conditions on unit 10.
    IF (KZ(22) == 1 .AND. rank == 0) THEN
        DO 85 I = 1,N
            WRITE (10,*) BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3)
        85 END DO
    !  100 FORMAT(1X,1P,7(1X,E26.10))
    END IF
    CALL FLUSH(10)
    CLOSE(10)

    RETURN

    END SUBROUTINE START
