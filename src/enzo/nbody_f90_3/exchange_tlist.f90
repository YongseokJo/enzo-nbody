    subroutine exchange_tlist(I,J,STEP,DTK)


!     Exchange I and J index in NXTLST
!     --------------------------------------------------------

    include 'params.h'
    include 'tlist.h'
    REAL*8 :: STEP(NMAX),DTK(64)
    INTEGER :: I,J,LI,LJ,L,LL

!     Get I and J step levels
    LI = k_step(STEP(I),DTK)
    LJ = k_step(STEP(J),DTK)

!     --07/08/14 16:49-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      print*,'exchange I',I,'J',J,'IK',LI,'JK',LJ,'NLSTDELAY',
!$$$     &     NLSTDELAY(1:NLSTDELAY(1)+1)
!     --07/08/14 16:49-lwang-end----------------------------------------*

!     IF I and J are same level, do nothing
    IF(LI == LJ) RETURN

!     Check ghost
    IF(LI == -1) THEN
        IF(NGHOSTS > 0) THEN
            L = NXTLIMIT + 1
            11 IF(NXTLST(L) == I) THEN
                NXTLST(L) =  J
                GO TO 30
            ELSE
                L = L + 1
                IF(L-NXTLIMIT <= NGHOSTS) GO TO 11
            END IF
        END IF
        IF(NLSTDELAY(1) > 0) GO TO 12
        write(6,*) 'Error: No ghost particle I',I,'STEP',STEP(I)
        call flush(6)
        call abort()
    ELSE
        LL = NDTK(LI+1)+1
        1 IF(NXTLST(LL) == I) THEN
            NXTLST(LL) = J
            GO TO 30
        ELSE
            IF(LL < NDTK(LI)) THEN
                LL = LL + 1
                GO TO 1
            ELSE
                IF(NLSTDELAY(1) > 0) GO TO 12
                write(6,*) 'Error: Index ',I,' not found in step level ', &
                LI,'!'
                call flush(6)
                call abort()
            END IF
        END IF
    END IF

!     Search delay list
    12 DO L = 2, NLSTDELAY(1)+1
        IF(NLSTDELAY(L) == I) THEN
            NLSTDELAY(L) = J
            GO TO 30
        END IF
    END DO
    write(6,*) 'Error: No delay particle I',I,'STEP',STEP(I)
    call flush(6)
    call abort()

!     Check ghost
    30 IF(LJ == -1) THEN
        IF(NGHOSTS > 0) THEN
            L = NXTLIMIT + 1
            21 IF(NXTLST(L) == J) THEN
                NXTLST(L) =  I
                RETURN
            ELSE
                L = L + 1
                IF(L-NXTLIMIT <= NGHOSTS) GO TO 21
            END IF
        END IF
        IF(NLSTDELAY(1) > 0) GO TO 22
        write(6,*) 'Error: No ghost particle J',J,'STEP',STEP(J)
        call flush(6)
        call abort()
    ELSE
        LL = NDTK(LJ+1)+1
        2 IF(NXTLST(LL) == J) THEN
            NXTLST(LL) = I
            RETURN
        ELSE
            IF(LL < NDTK(LJ)) THEN
                LL = LL + 1
                GO TO 2
            ELSE
                IF(NLSTDELAY(1) > 0) GO TO 22
                write(6,*) 'Error: Index ',J,' not found in step level ', &
                LJ,'!'
                call flush(6)
                call abort()
            END IF
        END IF
    END IF

    22 DO L = 2, NLSTDELAY(1)+1
        IF(NLSTDELAY(L) == J) THEN
            NLSTDELAY(L) = I
            RETURN
        END IF
    END DO

    write(6,*) 'Error: No delay particle J',J,'STEP',STEP(J)
    call flush(6)
    call abort()

    RETURN

    end subroutine exchange_tlist
