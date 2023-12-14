    SUBROUTINE sort_tlist(STEP,DTK,FLAG)


!     Sort block list array by steps
!     -----------------------------------

!     STEP: STEP array;  DTK: STEP of different level, FLAG: Full sorting flag

    include 'params.h'
    include 'tlist.h'
    REAL*8 :: STEP(NMAX),DTK(64)
    INTEGER :: J,I,K,L,IR,LTMP,NTC,LK
    LOGICAL :: FLAG

!     Heapsort method (Press p. 231).
    IF (NXTLIMIT <= 1) RETURN

!     If full sorting flag is false, jump to 30
    IF ( .NOT. FLAG) GO TO 30
          
    L = NXTLIMIT/2 + 1
    IR = NXTLIMIT
    10 CONTINUE
    IF (L > 1) THEN
        L = L - 1
        LTMP = NXTLST(L)
    ELSE
        LTMP = NXTLST(IR)
        NXTLST(IR) = NXTLST(1)
    !     Get criterion of steps
        IF (IR == NXTLIMIT) NTC = 1
        11 IF (STEP(NXTLST(IR)) < DTK(NTC)) THEN
            NTC = NTC + 1
            IF (NTC < 64) THEN
                NDTK(NTC) = IR
                GO TO 11
            ELSE
                write(6,*) 'Error!: Too small Step: I',NXTLST(IR), &
                'S',STEP(NXTLST(IR)),'DTK(64)',DTK(64)
                call flush(6)
                call abort()
            END IF
        END IF
        IF (IR == NXTLIMIT) THEN
            NDTMIN = NTC
            NDTK(NTC) = IR
            IF(NTC > 1) NDTK(1:NTC-1) = IR
        END IF
    !         NXTK(IR) = NTC
    !         NDTK(NTC) = IR
    
        IF (IR == 1) THEN
            NXTLST(1) = LTMP
            NDTMAX = NTC
        !     Reset unused level NDTK to zero
            IF(NTC < 64) NDTK(NTC+1:64) = 0
            RETURN
        END IF
        IR = IR - 1
    END IF
    I = L
    J = L + L
    20 IF (J <= IR) THEN
        IF (J < IR) THEN
            IF (STEP(NXTLST(J)) < STEP(NXTLST(J+1))) J = J + 1
        END IF
        IF (STEP(LTMP) < STEP(NXTLST(J))) THEN
            NXTLST(I) = NXTLST(J)
            I = J
            J = J + J
        ELSE
            J = IR + 1
        END IF
        GO TO 20
    END IF
    NXTLST(I) = LTMP
    GO TO 10

!     Only sorting before NXTLEN (NXTLEVEL), be careful not update NXTLEVEL earlier than sort
    30 DO 40 K = NDTMAX,NXTLEVEL,-1
        II = NDTK(K+1)+1
    !     Loop from NDTK(K+1)+1 to NDTK(K), since NDTK(K) will change, thus avoid use 'DO structure'
        39 IF(II > NDTK(K)) GO TO 40
        J = NXTLST(II)
    !     Detect step reduce
        NTC = K
        IF(STEP(J) < DTK(K)) THEN
            41 IF(STEP(J) < DTK(NTC)) THEN
                IF(NTC < 64) THEN
                    NTC = NTC + 1
                    GO TO 41
                ELSE
                    write(6,*) 'Error!: Too small Step: J',J,'K',K, &
                    'L',II,'S',STEP(J),'NDTK(K+1)',NDTK(K+1), &
                    'NDTK(K)',NDTK(K),'NDTMIN',NDTMIN,'NDTMAX', &
                    NDTMAX,'DTK(64)',DTK(64)
                    call flush(6)
                    call abort()
                END IF
            END IF
        !     Update smallest step level indicator
            IF(NDTMAX < NTC) NDTMAX = NTC
        !     Shift first index of current level to position J
            NXTLST(II) = NXTLST(NDTK(K+1)+1)
            DO LK = K+1, NTC-1
            !     Increase NDTK by one
                NDTK(LK) = NDTK(LK) + 1
            !     Shift first index to the end position of level LK
                NXTLST(NDTK(LK)) = NXTLST(NDTK(LK+1)+1)
            !     NXTK(NDTK(LK)) = LK
            END DO
        !     Increase target level NDTK by one
            NDTK(NTC) = NDTK(NTC) + 1
        !     Store index I in the last position of new level NL
            NXTLST(NDTK(NTC)) = J
        !     NXTK(NDTK(NTC)) = NTC
        ELSE IF(STEP(J) > DTK(K)) THEN
        !     Move ghost or special particle outside block list NDTK(1)
            IF(STEP(J) > DTK(1)) then
                NTC = 1
                NXTLIMIT = NXTLIMIT - 1
                NGHOSTS = NGHOSTS + 1
            ELSE
                42 IF(STEP(J) > DTK(NTC)) THEN
                    NTC = NTC - 1
                    GO TO 42
                END IF
            !     Update maximun step level indicator
                IF(NDTMIN > NTC) NDTMIN = NTC
            END IF
        !     Shift last index of current level to position J
            NXTLST(II) = NXTLST(NDTK(K))
            DO LK = K, NTC+1, -1
            !     Shift last index to the first position of level LK
                NXTLST(NDTK(LK)) = NXTLST(NDTK(LK-1))
            !     NXTK(NDTK(LK)) = LK-1
            !     Reduce NDTK(LK) by one
                NDTK(LK) = NDTK(LK) - 1
            !     Safe check
                IF(NDTK(LK) < 0) then
                    write(6,*) 'Error: NDTK LT ZERO! ', &
                    'K',K,'New K',NTC,'LOOP K',LK,'J',J, &
                    'L',II,'NDTMIN',NDTMIN,'NDTMAX',NDTMAX, &
                    'NDTK(1:NDTMAX)',NDTK(1:NDTMAX+1)
                    call flush(6)
                    call abort()
                END IF
            END DO
        !     Store index I in the last position of new level NL
            NXTLST(NDTK(NTC)) = J
        !     NXTK(NDTK(NTC)) = NTC
        !     Check same position since it's replaced by last particle index
            GO TO 39
        END IF
    !     Loop next
        II = II + 1
        GO TO 39

    40 END DO

    CALL shrink_tlist

    RETURN

    end SUBROUTINE sort_tlist

