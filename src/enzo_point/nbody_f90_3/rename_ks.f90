    SUBROUTINE RENAME_KS


!       Renaming of list arrays.
!       ------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    #ifdef SIMD
    LOGICAL :: FLAG_SEND
    FLAG_SEND = .false.
    #endif


!       Old & new names are shown explicitly for three different cases.

!       Standard sequence     ICOMP = 2*NPAIRS-1     ICOMP = 2*NPAIRS
!       -------------------------------------------------------------------
!       2*NPAIRS-1 ICOMP       2*NPAIRS-1 ICOMP       2*NPAIRS-1 JCOMP
!       2*NPAIRS   JCOMP       ICOMP      2*NPAIRS-1  ICOMP      2*NPAIRS-1
!       ICOMP      2*NPAIRS-1  2*NPAIRS   JCOMP       2*NPAIRS   2*NPAIRS-1
!       JCOMP      2*NPAIRS    JCOMP      2*NPAIRS    JCOMP      2*NPAIRS
!       -------------------------------------------------------------------

!       Form the sequential names of exchanged & regularized bodies.
    ILIST(1) = 2*NPAIRS - 1
    ILIST(2) = 2*NPAIRS
    ILIST(3) = ICOMP
    ILIST(4) = JCOMP
!       Save the corresponding new names for the renaming procedure.
    ILIST(5) = ILIST(3)
    ILIST(6) = ILIST(4)
    ILIST(7) = ILIST(1)
    ILIST(8) = ILIST(2)
!       Ensure consecutive search list and modify renaming list accordingly.
    IF (ICOMP <= ILIST(2)) THEN
        ILIST(3) = ILIST(2)
        ILIST(2) = ICOMP
    !       Note that all the new names are affected by the rearrangement.
        ILIST(5) = ILIST(2)
        ILIST(6) = ILIST(1)
        ILIST(7) = ILIST(4)
        ILIST(8) = ILIST(3)
        IF (ICOMP /= ILIST(1)) THEN
        !       First single particle is exchanged twice if ICOMP = 2*NPAIRS.
            ILIST(5) = JCOMP
            ILIST(7) = ILIST(1)
        !       New name of 2*NPAIRS is not really needed because of duplicity.
        END IF
    END IF

!       Rename neighbour lists with exchanged or regularized components.
    DO 40 J = 1,NTOT-1
        #ifdef SIMD
        FLAG_SEND = .false.
        #endif
    !       Skip modification if first neighbour comes after second component.
        IF (LIST(2,J) > JCOMP) GO TO 39
    !       No special precaution needed for case of no neighbours.
        NNB = LIST(1,J)
        NNB1 = 0
        NIJ = 0
    
    !       See whether any neighbours need to be renamed.
        K1 = 1
        DO 25 L = 2,NNB+1
            IF (LIST(L,J) > JCOMP) GO TO 30
        !       Start search loop above last identified value (bug fix 10/07).
            DO 20 K = K1,4
                IF (LIST(L,J) == ILIST(K)) THEN
                    K1 = K + 1
                    NNB1 = NNB1 + 1
                !       Save corresponding list location for the modification loop.
                    JLIST(NNB1) = K
                    JLIST(NNB1+4) = L
                !       Count number of identified KS components.
                    IF (ILIST(K) == ICOMP .OR. ILIST(K) == JCOMP) &
                    NIJ = NIJ + 1
                !       Advance neighbour list after each identification (note duplicity).
                    GO TO 25
                END IF
            20 END DO
        25 END DO
    
    !       Skip modification if no neighbours need to be renamed.
        30 IF (NNB1 == 0) GO TO 39
    
    !       Include c.m. renaming procedure if both components are neighbours.
        DO 38 K = 1,NNB1
            KTIME = JLIST(K)
        !       Indicator for identified members of search list.
            INEW = ILIST(KTIME+4)
        !       Replace identified components (or single perturber) by c.m.
            IF (INEW < IFIRST) THEN
                IF (NIJ == 2 .OR. J < IFIRST - 3) INEW = NTOT
            END IF
        !       Start with the saved list index (may differ by one either way).
            LS = JLIST(K+4)
        !       See if index needs adjusting after previous modification.
            IF (LIST(LS,J) < ILIST(KTIME)) LS = LS + 1
            IF (LIST(LS,J) > ILIST(KTIME)) LS = LS - 1
        
        !       Move list members up or down unless new name is identical to old.
            IF (LS <= 0) LS = 2
            LL = LS
            IF (INEW < ILIST(KTIME)) THEN
            !       Move list members down by one until location for new name is vacated.
                DO 32 L = LL,3,-1
                    IF (LIST(L-1,J) > INEW) THEN
                        LIST(L,J) = LIST(L-1,J)
                        LS = L - 1
                    ELSE
                        GO TO 36
                    END IF
                32 END DO
            ELSE IF (INEW > ILIST(KTIME)) THEN
            !       Move list members up by one until sequential location is reached.
                DO 34 L = LL,NNB
                    IF (LIST(L+1,J) < INEW) THEN
                        LIST(L,J) = LIST(L+1,J)
                        LS = L + 1
                    ELSE
                        GO TO 36
                    END IF
                34 END DO
            END IF
        !       Set renamed neighbour in sequential location.
            36 LIST(LS,J) = INEW
        38 END DO
    
    !       Reduce membership by one if c.m. set in last two locations.
        IF (NIJ == 2) LIST(1,J) = NNB - 1
        #ifdef SIMD
        IF (J >= IFIRST) FLAG_SEND = .TRUE. 
        #endif

        39 IF (J >= IFIRST-2) THEN
        !     Check replacing of single KS component by corresponding c.m.
            37 IF(LIST(1,J) > 0 .AND. LIST(2,J) < IFIRST) THEN
                JNEW = KVEC(LIST(2,J)) + N
                NNB = LIST(1,J)
                DO L = 2,NNB+1
                    IF (L <= NNB .AND. LIST(L+1,J) < JNEW) THEN
                        LIST(L,J) = LIST(L+1,J)
                    ELSE
                        LIST(L,J) = JNEW
                        #ifdef SIMD
                        IF(J >= IFIRST) FLAG_SEND = .TRUE. 
                        #endif
                    !     Check again until first neighbour > 2*NPAIRS.
                        GO TO 37
                    END IF
                END DO
            END IF
                         
        END IF
        #ifdef SIMD
    !     Update neighbor list in AVX/SSE library
        IF(TIME /= 0.0D0 .AND. FLAG_SEND) THEN
        !     --04/20/14 11:32-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            IF(J.EQ.7660.and.rank.eq.0) THEN
        !$$$               PRINT*,'J',J,'N',NAME(J),'LIST',LIST(1:LIST(1,J)+1,J),
        !$$$     &          'IFIRST',IFIRST
        !$$$               call flush(6)
        !$$$            end if
        !     --04/20/14 11:32-lwang-end----------------------------------------*
            CALL IRR_SIMD_SET_LIST(J,LIST(1,J))
        END IF
        #endif
                 
    40 END DO

!       Update the list of recently regularized particles.
    41 NNB = LISTR(1)
    DO 44 L = 2,NNB+1
    !       First see whether either component has been regularized before.
        IF (LISTR(L) == ICOMP .OR. LISTR(L) == JCOMP) THEN
        !       Remove corresponding pair even if only one component is present.
            J = 2*KVEC(L-1)
        !       Move up the subsequent pairs and reduce membership by two.
            DO 42 K = J,NNB-1
                LISTR(K) = LISTR(K+2)
            42 END DO
            LISTR(1) = LISTR(1) - 2
        !       Make a new search otherwise LISTR -> -2 if NNB = 2.
            GO TO 41
        END IF
    44 END DO

!       Rename exchanged components if present in the list.
    DO 48 KCOMP = 1,2
        DO 46 L = 2,NNB+1
            IF (LISTR(L) == 2*NPAIRS - 2 + KCOMP) THEN
                IF (KCOMP == 1) LISTR(L) = ICOMP
                IF (KCOMP == 2) LISTR(L) = JCOMP
            END IF
        46 END DO
    48 END DO

!       Update the list of high velocity particles.
    IF (LISTV(1) == 0) GO TO 70
    NNB = LISTV(1)
    L = 1
    50 L = L + 1
!       Check for removal of regularized component.
    IF (LISTV(L) == ICOMP .OR. LISTV(L) == JCOMP) THEN
        DO 52 K = L,NNB
            LISTV(K) = LISTV(K+1)
        52 END DO
        LISTV(1) = LISTV(1) - 1
        NNB = NNB - 1
    !       Consider the same location again.
        L = L - 1
    END IF
    IF (L <= NNB) GO TO 50

!       Rename exchanged components.
    DO 58 KCOMP = 1,2
        DO 56 L = 2,NNB+1
            IF (LISTV(L) == 2*NPAIRS - 2 + KCOMP) THEN
                IF (KCOMP == 1) LISTV(L) = ICOMP
                IF (KCOMP == 2) LISTV(L) = JCOMP
            END IF
        56 END DO
    58 END DO

!       Remove any fast particles which have slowed down or are outside 2<R>.
    L = 1
    60 L = L + 1
    IF (LISTV(1) == 0) GO TO 70
    J = LISTV(L)
    call jpred(J,TIME,TIME)
    A0 = XDOT(1,J)**2 + XDOT(2,J)**2 + XDOT(3,J)**2
    A2 = (X(1,J) - RDENS(1))**2 + (X(2,J) - RDENS(2))**2 + &
    (X(3,J) - RDENS(3))**2
    IF (A0 < 16.0*ECLOSE .OR. A2 > 4.0*RSCALE**2) THEN
        DO 65 K = L,NNB
            LISTV(K) = LISTV(K+1)
        65 END DO
        LISTV(1) = LISTV(1) - 1
        NNB = NNB - 1
    !       Consider the same location again.
        L = L - 1
    END IF
    IF (L <= NNB) GO TO 60

    70 RETURN

    END SUBROUTINE RENAME_KS
