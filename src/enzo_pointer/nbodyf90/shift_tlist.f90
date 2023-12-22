    SUBROUTINE shift_tlist(I,DT,DL,DTK)


!     Shift particle in block list due to LIST_SF
!     -----------------------------------

!     I: particle index
!     DT: particle irregular step
!     DL: shifting step level

    include 'params.h'
    include 'tlist.h'
    REAL*8 :: DT,DTK(64)
    INTEGER :: I,DL
    INTEGER :: L,J,NL

!     Current step level
    L = K_STEP(DT,DTK)
!     escape for ghost case
    IF(L < 0) RETURN
!     Particle position in NXTLST
    J = NDTK(L+1)+1
    1 IF(NXTLST(J) /= I) THEN
        IF(J < NDTK(L)) THEN
            J = J + 1
            GO TO 1
        ELSE
            IF(NLSTDELAY(1) > 0) THEN
                J = 2
            !     If find in NLSTDELAY, just return
                2 IF (NLSTDELAY(J) == I) RETURN
                J = J + 1
                IF (J <= NLSTDELAY(1)+1) GO TO 2
            END IF
            write(6,*) 'Error: Index ',I,' not found in step level ', &
            L,'!'
            call flush(6)
            call abort()
        END IF
    END IF
!     New step level
    NL = L + DL
!     See whether the shift is reasonable
    IF(NL > 64) then
        write(6,*) 'Error!: Too small Step: I',I, &
        'Step level',NL
        call flush(6)
        call abort()
    END IF
!     See whether need to decrease NDTMIN or increase NDTMAX
    IF(NL < NDTMIN .AND. NL >= 1) NDTMIN = NL
    IF(NL > NDTMAX) NDTMAX = NL
!     If exceed maximum step, treat as special particle like ghost and move index outside nxtlimit
    IF(NL < 1) THEN
        NL = 1
        NXTLIMIT = NXTLIMIT - 1
    END IF
!     Shifting to larger step level (smaller time step)
    IF(DL > 0) THEN
    !     Shift first index of current level to position J
        NXTLST(J) = NXTLST(NDTK(L+1)+1)
        DO LK = L+1, NL-1
        !     Increase NDTK by one
            NDTK(LK) = NDTK(LK) + 1
        !     Shift first index to the end position of level LK
            NXTLST(NDTK(LK)) = NXTLST(NDTK(LK+1)+1)
        !            NXTK(NDTK(LK)) = LK
        END DO
    !     Increase target level NDTK by one
        NDTK(NL) = NDTK(NL) + 1
    !     Store index I in the last position of new level NL
        NXTLST(NDTK(NL)) = I
    !         NXTK(NDTK(NL)) = NL
    !     Shifting to smaller step level (larger time step)
    ELSE
    !     Shift last index of current level to position J
        NXTLST(J) = NXTLST(NDTK(L))
        DO LK = L, NL+1, -1
        !     Shift last index to the first position of level LK
            NXTLST(NDTK(LK)) = NXTLST(NDTK(LK-1))
        !            NXTK(NDTK(LK)) = LK-1
        !     Reduce NDTK(LK) by one
            NDTK(LK) = NDTK(LK) - 1
        END DO
    !     Store index I in the last position of new level NL
        NXTLST(NDTK(NL)) = I
    !         NXTK(NDTK(NL)) = NL
    END IF

!     Check whether need to modify NDTMAX and NDTMIN
    50 IF(NDTK(NDTMAX) == 0 .AND. NDTMAX > NDTMIN) THEN
        NDTMAX = NDTMAX - 1
        GO TO 50
    END IF

    51 IF(NDTK(NDTMIN) == NDTK(NDTMIN+1) &
     .AND. NDTK(NDTMIN) == NXTLIMIT) THEN
        NDTMIN = NDTMIN + 1
        GO TO 51
    END IF

    RETURN

    end SUBROUTINE shift_tlist
