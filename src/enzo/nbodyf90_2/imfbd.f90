
! Note: this routine contains the correction from Carsten Weidner
!       (Nov.2004), that for LOM above a certain value the integral
!       was not computed correctly.


!===========================================================
    REAL*8 :: FUNCTION IMFBD(XX,LOM,UPM)
    IMPLICIT NONE
    INCLUDE 'mpi_base.h'


!        5-part-power-law IMF extending into BD regime.
!        ---------------------------------------------
!                 (P.Kroupa, Feb.99, Heidelberg)

! On input of a random deviate XX and the lower (LM) and upper (UM)
! mass limits (in Msun), IMFBD returns a mass in Msun.
! (For equal masses, LOM=UPM and IMFBD = UPM.)
! LOM and UPM must be kept unchanged.


    REAL*8 ::    ML,MH,M0,M1,M2,MU
    REAL*8 ::    ALPHA0,ALPHA1,ALPHA2,ALPHA3,ALPHA4
    REAL*8 ::    K,XH,X0,X1,X2,XU,MUP
    REAL*8 ::    MTOT,MASSH,MASS0,MASS1,MASS2,MASSU
    REAL*8 ::    XX,LOM,UPM,LM,UM,ZM
    REAL*8 ::    BETA,INTGR,MINTGR,MGEN,XIN
    EXTERNAL INTGR,MINTGR,MGEN
    LOGICAL ::  FIRSTBD

    SAVE  FIRSTBD,XH,X0,X1,X2,XU,K,LM,UM
    DATA  FIRSTBD /.TRUE./

! -------------------------------------------------------------------
!       MF (Kroupa, Tout & Gilmore, 1993, MNRAS, 262, 545)
!       extended into BD regime (masses in Msun).
!       (Note: ALPHA0 is between ML and MH, etc.)

!      DATA    ML,     MH,    M0,      M1,      M2,    MU/
!     &      0.01D0, 0.08D0, 0.5D0,  1.0D0,   8.0D0, 100.0D0/
!      DATA      ALPHA0, ALPHA1, ALPHA2, ALPHA3, ALPHA4/
!     &          0.3D0,  1.3D0,  2.2D0,   2.7D0,  2.7D0/

!..................................................................
!       Standard of average IMF: (Kroupa 2001)

    DATA    ML,     MH,    M0,      M1,      M2,    MU/ &
    0.01D0, 0.08D0, 0.5D0,  1.0D0,   8.0D0, 100.0D0/
    DATA      ALPHA0, ALPHA1, ALPHA2, ALPHA3, ALPHA4/ &
    0.3D0,  1.3D0,  2.3D0,   2.3D0,  2.3D0/



!      DATA    ML,     MH,    M0,      M1,      M2,    MU/
!     &      0.01D0, 0.08D0, 0.5D0,  1.0D0,   8.0D0, 150.0D0/
!      DATA      ALPHA0, ALPHA1, ALPHA2, ALPHA3, ALPHA4/
!     &          0.3D0,  1.3D0,  2.3D0,   2.35D0,  2.35D0/


!.....................................................................
!       Salpeter MF:

!      DATA   ML,       MH,    M0,      M1,      M2,    MU/
!     &     0.01D0,   0.08D0, 0.08D0,   1.0D0,  8.0D0, 100.0D0/
!      DATA    ALPHA0,  ALPHA1,   ALPHA2,  ALPHA3, ALPHA4/
!     &        0.3D0,   0.3D0,    1.2D0,   2.3D0,  2.3D0/

!.....................................................................
!       Flat MF:

!      DATA   ML,       MH,    M0,      M1,      M2,    MU/
!     &     0.01D0,   0.08D0, 0.5D0,  1.0D0,   5.0D0,  100.0D0/
!      DATA     ALPHA0,   ALPHA1, ALPHA2,  ALPHA3, ALPHA4/
!     &          1.D0,     1.D0,   1.D0,   1.D0,    1.D0/
! --------------------------------------------------------------------

!     Set initial ZM to negative value
    ZM = -1.D0

    IF ( .NOT. FIRSTBD .AND. UM == LM)THEN
        ZM = UM
        GOTO 500
    END IF

    IF (FIRSTBD) THEN
        LM = LOM
        UM = UPM
    
    !       Remind user of encoded MF:
    
        if(rank == 0)then
            WRITE(6,*)
            WRITE(6,'(a)')' The 5-part power-law MF:'
            WRITE(6,'(a)')'   ML      MH      M0      M1      M2      MU'
            WRITE(6,'(6F8.3)')ML,MH,M0,M1,M2,MU
            WRITE(6,'(a)')'      ALPHA0  ALPHA1  ALPHA2  ALPHA3  ALPHA4'
            WRITE(6,'(3x,6F8.3)')ALPHA0,ALPHA1,ALPHA2,ALPHA3,ALPHA4
            WRITE(6,*)
        end if
    
    !       Check MF:
    
        IF (rank == 0 .AND. ML > MH) THEN
            WRITE(6,*)
            WRITE(6,*)' FATAL: ML > MH'
            WRITE(6,*)' STOP'
            WRITE(6,*)
            STOP
        END IF
        IF (rank == 0 .AND. MH > M0) THEN
            WRITE(6,*)
            WRITE(6,*)' FATAL: MH > M0'
            WRITE(6,*)' STOP'
            WRITE(6,*)
            STOP
        END IF
        IF (rank == 0 .AND. M0 > M1) THEN
            WRITE(6,*)
            WRITE(6,*)' FATAL: M0 > M1'
            WRITE(6,*)' STOP'
            WRITE(6,*)
            STOP
        END IF
        IF (rank > 0 .AND. M1 > M2) THEN
            WRITE(6,*)
            WRITE(6,*)' FATAL: M1 > M2'
            WRITE(6,*)' STOP'
            WRITE(6,*)
            STOP
        END IF
        IF (rank > 0 .AND. M2 > MU) THEN
            WRITE(6,*)
            WRITE(6,*)' FATAL: M2 > MU'
            WRITE(6,*)' STOP'
            WRITE(6,*)
            STOP
        END IF
    
    !       Check for physical input mass range.
    
        IF (UM == LM) THEN
            if(rank == 0)then
                WRITE(6,*)
                WRITE(6,*)' WARNING: UM (= ',UM,' Msun) = LM= ', &
                LM,' Msun'
                WRITE(6,*)
            end if
            FIRSTBD = .FALSE.
            ZM = UM
            GOTO 500
        END IF
    
        IF (rank == 0 .AND. UM < LM) THEN
            WRITE(6,*)
            WRITE(6,*)' FATAL: UM (= ',UM,' Msun) < LM= ', &
            LM,' Msun'
            WRITE(6,*)' STOP'
            WRITE(6,*)
            STOP
        END IF
        IF (rank == 0 .AND. LM > MU) THEN
            WRITE(6,*)
            WRITE(6,*)' FATAL: BODYN (LM= ',LM,' Msun) > MU'
            WRITE(6,*)' STOP'
            WRITE(6,*)
            STOP
        END IF
        IF (LM < ML) THEN
            if(rank == 0)then
                WRITE(6,*)
                WRITE(6,*)' WARNING: BODYN (LM= ',LM,' Msun) < ML'
                WRITE(6,*)' BODYN = ',ML, ' Msun adopted.'
                WRITE(6,*)
            end if
            LM = ML
        END IF
        IF (UM > MU) THEN
            if(rank == 0)then
                WRITE(6,*)
                WRITE(6,*)' WARNING: BODY1 (UM= ',UM,' Msun) > MU= ', &
                MU,' Msun'
                WRITE(6,*)' BODY1 = ',MU, ' Msun adopted.'
                WRITE(6,*)
            end if
            UM = MU
        END IF
    
    !       Determine normalisation K for unit probablity (LM<=M<=UM).
    
        K = 0.0D0
        XH = 0.0D0
        X0 = 0.0D0
        X1 = 0.0D0
        X2 = 0.0D0
        XU = 0.0D0
    
        IF (LM >= ML .AND. LM <= MH) THEN
        
            IF (UM <= MH) THEN
                MUP = UM
            ELSE IF (UM > MH) THEN
                MUP = MH
            END IF
            BETA = MH**ALPHA0
            XH = INTGR(LM,MUP,ALPHA0,BETA)
        
            IF (UM > MH) THEN
                IF (UM <= M0) THEN
                    MUP = UM
                ELSE IF (UM > M0) THEN
                    MUP = M0
                END IF
                BETA = MH**ALPHA1
                X0 = INTGR(MH,MUP,ALPHA1,BETA)
            END IF
        
            IF (UM > M0) THEN
                IF (UM <= M1) THEN
                    MUP = UM
                ELSE IF (UM > M1) THEN
                    MUP = M1
                END IF
                BETA = M0**ALPHA2 * (MH/M0)**ALPHA1
                X1 = INTGR(M0,MUP,ALPHA2,BETA)
            END IF
        
            IF (UM > M1) THEN
                IF (UM <= M2) THEN
                    MUP = UM
                ELSE IF (UM > M2) THEN
                    MUP = M2
                END IF
                BETA = M1**ALPHA3 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2
                X2 = INTGR(M1,MUP,ALPHA3,BETA)
            END IF
        
            IF (UM > M2) THEN
                MUP = UM
                BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
                (M1/M2)**ALPHA3
                XU = INTGR(M2,MUP,ALPHA4,BETA)
            END IF
        
        ELSE IF (LM > MH .AND. LM <= M0) THEN
        
            IF (UM <= M0) THEN
                MUP = UM
            ELSE IF (UM > M0) THEN
                MUP = M0
            END IF
            BETA = MH**ALPHA1
            X0 = INTGR(LM,MUP,ALPHA1,BETA)
        
            IF (UM > M0) THEN
                IF (UM <= M1) THEN
                    MUP = UM
                ELSE IF (UM > M1) THEN
                    MUP = M1
                END IF
                BETA = M0**ALPHA2 * (MH/M0)**ALPHA1
                X1 = INTGR(M0,MUP,ALPHA2,BETA)
            END IF
        
            IF (UM > M1) THEN
                IF (UM <= M2) THEN
                    MUP = UM
                ELSE IF (UM > M2) THEN
                    MUP = M2
                END IF
                BETA = M1**ALPHA3 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2
                X2 = INTGR(M1,MUP,ALPHA3,BETA)
            END IF
        
            IF (UM > M2) THEN
                MUP = UM
                BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
                (M1/M2)**ALPHA3
                XU = INTGR(M2,MUP,ALPHA4,BETA)
            END IF
        
        ELSE IF (LM > M0 .AND. LM <= M1) THEN
        
            IF (UM <= M1) THEN
                MUP = UM
            ELSE IF (UM > M1) THEN
                MUP = M1
            END IF
            BETA = M0**ALPHA2 * (MH/M0)**ALPHA1
            X1 = INTGR(LM,MUP,ALPHA2,BETA)
        
            IF (UM > M1) THEN
                IF (UM <= M2) THEN
                    MUP = UM
                ELSE IF (UM > M2) THEN
                    MUP = M2
                END IF
                BETA = M1**ALPHA3 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2
                X2 = INTGR(M1,MUP,ALPHA3,BETA)
            END IF
        
            IF (UM > M2) THEN
                MUP = UM
                BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
                (M1/M2)**ALPHA3
                XU = INTGR(M2,MUP,ALPHA4,BETA)
            END IF
        
        ELSE IF (LM > M1 .AND. LM <= M2) THEN
        
            IF (UM <= M2) THEN
                MUP = UM
            ELSE IF (UM > M2) THEN
                MUP = M2
            END IF
            BETA = M1**ALPHA3 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2
            X2 = INTGR(LM,MUP,ALPHA3,BETA)
        
            IF (UM > M2) THEN
                MUP = UM
                BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
                (M1/M2)**ALPHA3
                XU = INTGR(M2,MUP,ALPHA4,BETA)
            END IF
        
        ELSE IF (LM > M2) THEN
        
            MUP = UM
            BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
            (M1/M2)**ALPHA3
            XU = INTGR(LM,MUP,ALPHA4,BETA)
        
        END IF
    
    
        K = XH + X0 + X1 + X2 + XU
        XH = XH/K
        X0 = X0/K
        X1 = X1/K
        X2 = X2/K
        XU = XU/K
    
    !       Write number fractions to log file.
    
        if(rank == 0)then
            WRITE(6,'(a,F7.3,a,F7.3,a)') &
            ' Number fractions for LM=',LM,', Msun, UM=',UM,' Msun:'
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') ML,' -',MH,' Msun:', XH
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') MH,' -',M0,' Msun:', X0
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') M0,' -',M1,' Msun:', X1
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') M1,' -',M2,' Msun:', X2
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') M2,' -',MU,' Msun:', XU
            WRITE(6,*)
        end if
    
    !       Determine average mass and mass fractions (for log-file).
    
        MTOT = 0.0D0
        MASSH = 0.0D0
        MASS0 = 0.0D0
        MASS1 = 0.0D0
        MASS2 = 0.0D0
        MASSU = 0.0D0
    
        IF (LM >= ML .AND. LM <= MH) THEN
        
            IF (UM <= MH) THEN
                MUP = UM
            ELSE IF (UM > MH) THEN
                MUP = MH
            END IF
            BETA = MH**ALPHA0
            MASSH = MINTGR(LM,MUP,ALPHA0,BETA)
        
            IF (UM > MH) THEN
                IF (UM <= M0) THEN
                    MUP = UM
                ELSE IF (UM > M0) THEN
                    MUP = M0
                END IF
                BETA = MH**ALPHA1
                MASS0 = MINTGR(MH,MUP,ALPHA1,BETA)
            END IF
        
            IF (UM > M0) THEN
                IF (UM <= M1) THEN
                    MUP = UM
                ELSE IF (UM > M1) THEN
                    MUP = M1
                END IF
                BETA = M0**ALPHA2 * (MH/M0)**ALPHA1
                MASS1 = MINTGR(M0,MUP,ALPHA2,BETA)
            END IF
        
            IF (UM > M1) THEN
                IF (UM <= M2) THEN
                    MUP = UM
                ELSE IF (UM > M2) THEN
                    MUP = M2
                END IF
                BETA = M1**ALPHA3 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2
                MASS2 = MINTGR(M1,MUP,ALPHA3,BETA)
            END IF
        
            IF (UM > M2) THEN
                MUP = UM
                BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
                (M1/M2)**ALPHA3
                MASSU = MINTGR(M2,MUP,ALPHA4,BETA)
            END IF
        
        ELSE IF (LM > MH .AND. LM <= M0) THEN
        
            IF (UM <= M0) THEN
                MUP = UM
            ELSE IF (UM > M0) THEN
                MUP = M0
            END IF
            BETA = MH**ALPHA1
            MASS0 = MINTGR(LM,MUP,ALPHA1,BETA)
        
            IF (UM > M0) THEN
                IF (UM <= M1) THEN
                    MUP = UM
                ELSE IF (UM > M1) THEN
                    MUP = M1
                END IF
                BETA = M0**ALPHA2 * (MH/M0)**ALPHA1
                MASS1 = MINTGR(M0,MUP,ALPHA2,BETA)
            END IF
        
            IF (UM > M1) THEN
                IF (UM <= M2) THEN
                    MUP = UM
                ELSE IF (UM > M2) THEN
                    MUP = M2
                END IF
                BETA = M1**ALPHA3 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2
                MASS2 = MINTGR(M1,MUP,ALPHA3,BETA)
            END IF
        
            IF (UM > M2) THEN
                MUP = UM
                BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
                (M1/M2)**ALPHA3
                MASSU = MINTGR(M2,MUP,ALPHA4,BETA)
            END IF
        
        ELSE IF (LM > M0 .AND. LM <= M1) THEN
        
            IF (UM <= M1) THEN
                MUP = UM
            ELSE IF (UM > M1) THEN
                MUP = M1
            END IF
            BETA = M0**ALPHA2 * (MH/M0)**ALPHA1
            MASS1 = MINTGR(LM,MUP,ALPHA2,BETA)
        
            IF (UM > M1) THEN
                IF (UM <= M2) THEN
                    MUP = UM
                ELSE IF (UM > M2) THEN
                    MUP = M2
                END IF
                BETA = M1**ALPHA3 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2
                MASS2 = MINTGR(M1,MUP,ALPHA3,BETA)
            END IF
        
            IF (UM > M2) THEN
                MUP = UM
                BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
                (M1/M2)**ALPHA3
                MASSU = MINTGR(M2,MUP,ALPHA4,BETA)
            END IF
        
        ELSE IF (LM > M1 .AND. LM <= M2) THEN
        
            IF (UM <= M2) THEN
                MUP = UM
            ELSE IF (UM > M2) THEN
                MUP = M2
            END IF
            BETA = M1**ALPHA3 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2
            MASS2 = MINTGR(LM,MUP,ALPHA3,BETA)
        
            IF (UM > M2) THEN
                MUP = UM
                BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
                (M1/M2)**ALPHA3
                MASSU = MINTGR(M2,MUP,ALPHA4,BETA)
            END IF
        
        ELSE IF (LM > M2) THEN
        
            MUP = UM
            BETA = M2**ALPHA4 * (MH/M0)**ALPHA1 * (M0/M1)**ALPHA2 * &
            (M1/M2)**ALPHA3
            MASSU = MINTGR(LM,MUP,ALPHA4,BETA)
        
        END IF
    
    
        MTOT = MASSH + MASS0 + MASS1 + MASS2 + MASSU
        MASSH = MASSH/MTOT
        MASS0 = MASS0/MTOT
        MASS1 = MASS1/MTOT
        MASS2 = MASS2/MTOT
        MASSU = MASSU/MTOT
    
    !       Write mass fractions to log file.
    
        if(rank == 0)then
            WRITE(6,'(a,F7.3,a,F7.3,a)') &
            ' Mass fractions for LM=',LM,', Msun, UM=',UM,' Msun:'
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') ML,' -',MH,' Msun:', MASSH
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') MH,' -',M0,' Msun:', MASS0
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') M0,' -',M1,' Msun:', MASS1
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') M1,' -',M2,' Msun:', MASS2
            WRITE(6,'(F7.3,a,F7.3,a,F7.4)') M2,' -',MU,' Msun:', MASSU
            WRITE(6,'(a,F7.3,a)')' Expected average stellar mass: ', &
            MTOT/K,' Msun'
            WRITE(6,*)
        end if
    
        FIRSTBD = .FALSE.
    END IF

!       Generate stellar mass from random deviate.

    IF (LM >= ML .AND. LM <= MH) THEN
    
        IF (XX <= XH) THEN
            XIN = XX
            BETA = K/MH**ALPHA0
            ZM = MGEN(XIN,LM,ALPHA0,BETA)
        ELSE IF (XX > XH .AND. XX <= (XH+X0)) THEN
            XIN = XX-XH
            BETA = K/MH**ALPHA1
            ZM = MGEN(XIN,MH,ALPHA1,BETA)
        ELSE IF (XX > (XH+X0) .AND. XX <= (XH+X0+X1)) THEN
            XIN = XX-XH-X0
            BETA = K/M0**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,M0,ALPHA2,BETA)
        ELSE IF (XX > (XH+X0+X1) .AND. XX <= (XH+X0+X1+X2)) THEN
            XIN = XX-XH-X0-X1
            BETA = K/M1**ALPHA3 * (M1/M0)**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,M1,ALPHA3,BETA)
        ELSE IF (XX > (XH+X0+X1+X2)) THEN
            XIN = XX-XH-X0-X1-X2
            BETA = K/M2**ALPHA4 * (M2/M1)**ALPHA3 * &
            (M1/M0)**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,M2,ALPHA4,BETA)
        END IF
    
    ELSE IF (LM > MH .AND. LM <= M0) THEN
    
        IF (XX > XH .AND. XX <= (XH+X0)) THEN
            XIN = XX-XH
            BETA = K/MH**ALPHA1
            ZM = MGEN(XIN,LM,ALPHA1,BETA)
        ELSE IF (XX > (XH+X0) .AND. XX <= (XH+X0+X1)) THEN
            XIN = XX-XH-X0
            BETA = K/M0**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,M0,ALPHA2,BETA)
        ELSE IF (XX > (XH+X0+X1) .AND. XX <= (XH+X0+X1+X2)) THEN
            XIN = XX-XH-X0-X1
            BETA = K/M1**ALPHA3 * (M1/M0)**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,M1,ALPHA3,BETA)
        ELSE IF (XX > (XH+X0+X1+X2)) THEN
            XIN = XX-XH-X0-X1-X2
            BETA = K/M2**ALPHA4 * (M2/M1)**ALPHA3 * &
            (M1/M0)**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,M2,ALPHA4,BETA)
        END IF
    
    ELSE IF (LM > M0 .AND. LM <= M1) THEN
    
        IF (XX > (XH+X0) .AND. XX <= (XH+X0+X1)) THEN
            XIN = XX-XH-X0
            BETA = K/M0**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,LM,ALPHA2,BETA)
        ELSE IF (XX > (XH+X0+X1) .AND. XX <= (XH+X0+X1+X2)) THEN
            XIN = XX-XH-X0-X1
            BETA = K/M1**ALPHA3 * (M1/M0)**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,M1,ALPHA3,BETA)
        ELSE IF (XX > (XH+X0+X1+X2)) THEN
            XIN = XX-XH-X0-X1-X2
            BETA = K/M2**ALPHA4 * (M2/M1)**ALPHA3 * &
            (M1/M0)**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,M2,ALPHA4,BETA)
        END IF
    
    ELSE IF (LM > M1 .AND. LM <= M2) THEN
    
        IF (XX > (XH+X0+X1) .AND. XX <= (XH+X0+X1+X2)) THEN
            XIN = XX-XH-X0-X1
            BETA = K/M1**ALPHA3 * (M1/M0)**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,LM,ALPHA3,BETA)
        ELSE IF (XX > (XH+X0+X1+X2)) THEN
            XIN = XX-XH-X0-X1-X2
            BETA = K/M2**ALPHA4 * (M2/M1)**ALPHA3 * &
            (M1/M0)**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,M2,ALPHA4,BETA)
        END IF
    
    ELSE IF (LM > M2) THEN
    
        IF (XX > (XH+X0+X1+X2)) THEN
            XIN = XX-XH-X0-X1-X2
            BETA = K/M2**ALPHA4 * (M2/M1)**ALPHA3 * &
            (M1/M0)**ALPHA2 * (M0/MH)**ALPHA1
            ZM = MGEN(XIN,LM,ALPHA4,BETA)
        END IF
    
    END IF

    500 IMFBD = ZM

!     Safe check
    IF(ZM < 0 .AND. rank == 0) THEN
        write(6,*) 'Error: mass not set in imfbd.f!'
        call flush(6)
        call abort()
    END IF

    RETURN
    END PROGRAM
!===========================================================
    real*8 :: FUNCTION INTGR(MA,MB,ALPHA,BETA)

!        Integral of mass function over mass interval MA to MB
!        -----------------------------------------------------

    IMPLICIT NONE
    REAL*8 ::   ALPHA,BETA,MA,MB


    IF (ALPHA /= 1.D0) THEN
        INTGR = (MB**(1.D0-ALPHA)-MA**(1.D0-ALPHA)) * BETA/(1.D0-ALPHA)
    ELSE IF (ALPHA == 1.D0) THEN
        INTGR = LOG(MB/MA) * BETA
    END IF

    RETURN
    END PROGRAM
!===========================================================
    real*8 :: FUNCTION MINTGR(MA,MB,ALPHA,BETA)

!        Total mass over mass interval MA to MB
!        --------------------------------------

    IMPLICIT NONE
    REAL*8 ::   ALPHA,BETA,MA,MB


    IF (ALPHA /= 2.D0) THEN
        MINTGR = (MB**(2.D0-ALPHA)-MA**(2.D0-ALPHA)) * BETA/(2.D0-ALPHA)
    ELSE IF (ALPHA == 2.D0) THEN
        MINTGR = LOG(MB/MA) * BETA
    END IF

    RETURN
    END PROGRAM
!===========================================================
    real*8 :: FUNCTION MGEN(XIN,MA,ALPHA,BETA)

!        Generate mass
!        -------------

    IMPLICIT NONE
    REAL*8 ::   ALPHA,BETA,XIN,MA


    IF (ALPHA /= 1.D0) THEN
        MGEN = (1.D0-ALPHA) * XIN * BETA
        MGEN = (MGEN + MA**(1.D0-ALPHA))**(1.D0/(1.D0-ALPHA))
    ELSE IF (ALPHA == 1.D0) THEN
        MGEN = MA * DEXP(XIN * BETA)
    END IF

    RETURN
    END PROGRAM
!===========================================================
