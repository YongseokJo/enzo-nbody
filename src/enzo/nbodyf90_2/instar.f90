    SUBROUTINE INSTAR


!       Initialization of stars.
!       ------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    REAL*8 :: TSCLS(20),LUMS(10),GB(10),TM,TN
    REAL*8 :: M0,M1,RM,LUM,AGE,MC,RCC
    REAL*8 :: MENV,RENV,K2,K3
    PARAMETER(K3=0.21D0)
    REAL*8 :: JSPIN,OSPIN
    REAL*8 :: VROTF
    EXTERNAL VROTF


!       Initialize mass loss variables & counters.
    TPHYS = 0.d0
    ZMRG = 0.d0
    ZMHE = 0.d0
    ZMRS = 0.d0
    ZMWD = 0.d0
    ZMSN = 0.d0
    ZMDOT = 0.d0
    NMDOT = 0
    NRG = 0
    NHE = 0
    NRS = 0
    NWD = 0
    NSN = 0
    NBH = 0
    NTZ = 0
    NBS = 0
    NKICK = 0
    NBKICK = 0

! B4 NEW
    ZMNH = 0.d0
    ZMBH = 0.d0
    ZMSY = 0.d0
    NROCHE = 0
    NCHAOS = 0
    IQCOLL = 0
    IBLUE = 0
    NCHA = 0
    NSP = 0
    NHG = 0
    NNH = 0
    NBR = 0
    NAS = 0
    NRO = 0
    NDD = 0
    NSPIR = 0
    NCIRC = 0
    NSLP = 0
    NCONT = 0
    NCOAL = 0
    NCE = 0
    INSTAB = 0
    NEINT = 0
    NEMOD = 0
    NHYP = 0
    NHYPC = 0
    NGB = 0
    NMS = N

    TMDOT = 1.0d+10
! B6  IF (KZ(27).EQ.0) THEN
!         TSTAR = TSCALE
!     END IF

!       Obtain optional stellar evolution parameters.
    IF (KZ(19) >= 3) THEN
        IF (ZMET > 0.03) ZMET = 0.03
        IF (ZMET < 0.0001) ZMET = 0.0001
        CALL zcnsts(ZMET,ZPARS)
        if(rank == 0) &
        WRITE (6,9)  ZPARS(11), ZPARS(12), ZMET
        9 FORMAT (//,12X,'ABUNDANCES:  X =',F6.3,'  Y =',F6.3, &
        '  Z =',F7.4)
    END IF
! B6  zpars(11)=xhyd
!     zpars(12)=yhel

!       Calculate scale factor for spin angular momentum.
    SPNFAC = ZMBAR*SU**2/(1.0D+06*TSTAR)

    EPOCH1 = EPOCH0
    DO 10 I = 1,N
    
    !       Obtain stellar parameters at current epoch.
        IF(KZ(19) == 2 .OR. KZ(19) == 4)THEN
        ! HECK 230109
        !             READ(21,*)M1,KW,M0,EPOCH1,OSPIN
            #ifdef PARALLEL
            CALL MPI_BCAST(M1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(KW,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(M0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(EPOCH1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(OSPIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            #endif
        ELSE
            M1 = BODY(I)*ZMBAR
            M0 = M1
            IF(KSTAR(I) > 1)THEN
                KW = KSTAR(I)
            ELSE
                KW = 1
            !               IF(M0.LE.0.01D0) KW = 10
            !               IF(M0.GE.100.D0) KW = 14
            ENDIF
        ENDIF
        MC = 0.D0
        AGE = TIME*TSTAR - EPOCH1
        CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
        CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS, &
        RM,LUM,KW,MC,RCC,MENV,RENV,K2)
    
    !       Assign initial spin angular momentum.
        IF(KZ(19) == 2 .OR. KZ(19) == 4)THEN
            JSPIN = (K2*(M1-MC)*RM**2 + K3*MC*RCC**2)*OSPIN
        ELSE
            OSPIN = 45.35d0*VROTF(M1)/RM
            JSPIN = K2*M1*RM**2*OSPIN
            IF(KW > 1) JSPIN = 0.D0
        ENDIF

    !       Convert from solar radii to scaled units (assume Sun = 1/215 AU).
        RADIUS(I) = RM/SU
        ZLMSTY(I) = LUM
        SPIN(I) = JSPIN/SPNFAC
    
        IF (KZ(27) == 3 .AND. KZ(28) == 4) THEN
            RADIUS(I) = 1.0D-12/(3.0*RBAR)
            KW = 13
            TEV(I) = 1000.0
            M0 = 25.0
            EPOCH(I) = -15.0/TSTAR
            IF (I == 1 .AND. rank == 0) WRITE (6,5) M1, KW, RADIUS(I)*SU
            5 FORMAT (/,12X,'FIRST STAR:    M K* R/SU ',F7.2,I4,1P,E9.1)
        END IF
    
    !       Initialize the stellar classification type (KW = 0 - 15).
        KSTAR(I) = KW
    
    !       Save the initial mass of all the stars in sequential order.
        BODY(I) = M1/ZMBAR
        BODY0(I) = M0/ZMBAR
    
    !       Set initial look-up time.
        EPOCH(I) = TIME*TSTAR - AGE
        TEV0(I) = TIME
        CALL TRDOT(I,DTM,M1)
        TEV(I) = DTM
        IF (KW >= 13) TEV(I) = 0.0
    
    !       Determine the time for next stellar evolution check.
        IF (TEV(I) < TMDOT) THEN
            TMDOT = TEV(I)
        END IF
    
    10 END DO

!       Define first quantized step < 1000 yrs (minimum interval for MDOT).
!       (changed to 100 yr to be consistent with trdot: 2/8/02)
    DT = 1.0d-04/TSCALE
    CALL STEPK(DT,DTN)
    IF(DTN*TSCALE < 100.0) DTN = 2.d0*DTN
    STEPX = DTN

!        Ensure binary components will be updated at the same time.
    DO 15 I = 1,NBIN0
        I1 = 2*I - 1
        I2 = I1 + 1
        TEV(I1) = MIN(TEV(I1),TEV(I2))
        TEV(I1) = MIN(TEV(I1),10.D0*STEPX)
        TMDOT = MIN(TMDOT,TEV(I1))
        TEV(I2) = TEV(I1)
    15 END DO

!       Initialize stellar collision matrix.

    ktype(0,0) = 1
    do 20 , j = 1,6
        ktype(0,j) = j
        ktype(1,j) = j
    20 END DO
    ktype(0,7) = 4
    ktype(1,7) = 4
    do 25 , j = 8,12
        if(j /= 10)then
            ktype(0,j) = 6
        else
            ktype(0,j) = 3
        endif
        ktype(1,j) = ktype(0,j)
    25 END DO
    ktype(2,2) = 3
    do 30 , i = 3,14
        ktype(i,i) = i
    30 END DO
    ktype(5,5) = 4
    ktype(7,7) = 1
    ktype(10,10) = 15
!       Change for CO+CO owing to Tout theory (21/11/08).
    ktype(11,11) = 12
    ktype(13,13) = 14
    do 35 , i = 2,5
        do 40 j = i+1,12
            ktype(i,j) = 4
        40 END DO
    35 END DO
    ktype(2,3) = 3
    ktype(2,6) = 5
    ktype(2,10) = 3
    ktype(2,11) = 5
    ktype(2,12) = 5
    ktype(3,6) = 5
    ktype(3,10) = 3
    ktype(3,11) = 5
    ktype(3,12) = 5
    ktype(6,7) = 4
    ktype(6,8) = 6
    ktype(6,9) = 6
    ktype(6,10) = 5
    ktype(6,11) = 6
    ktype(6,12) = 6
    ktype(7,8) = 8
    ktype(7,9) = 9
    ktype(7,10) = 7
    ktype(7,11) = 9
    ktype(7,12) = 9
    ktype(8,9) = 9
    ktype(8,10) = 7
    ktype(8,11) = 9
    ktype(8,12) = 9
    ktype(9,10) = 7
    ktype(9,11) = 9
    ktype(9,12) = 9
    ktype(10,11) = 9
    ktype(10,12) = 9
    ktype(11,12) = 12
    do 45 , i = 0,12
        ktype(i,13) = 13
        ktype(i,14) = 14
    45 END DO
    ktype(13,14) = 14

! Increase common-envelope cases by 100.
    do 50 , i = 0,9
        do 55 , j = i,14
            if(i <= 1 .OR. i == 7)then
                if(j >= 2 .AND. j <= 9 .AND. j /= 7)then
                    ktype(i,j) = ktype(i,j) + 100
                endif
            else
                ktype(i,j) = ktype(i,j) + 100
            endif
        55 END DO
    50 END DO

!       Assign the remaining values by symmetry.
    do 60 , i = 1,14
        do 65 , j = 0,i-1
            ktype(i,j) = ktype(j,i)
        65 END DO
    60 END DO

    if(rank == 0) &
    WRITE (6,75)  KTYPE
    75 FORMAT (/,11X,' KTYPE: ',15I4,14(/,19X,15I4))

    RETURN
    END SUBROUTINE INSTAR
