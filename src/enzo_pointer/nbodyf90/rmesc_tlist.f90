    subroutine rmesc_tlist(I,B_FLAG)


!     Remove particle in NXTLST for escape

!     B_FLAG: binary remove flag

    include 'params.h'
    include 'tlist.h'
    INTEGER :: L,J,LK,K
    LOGICAL :: RM_FLAG,B_FLAG

!     Remove flag.
    RM_FLAG = .false.

!     Step level tracer
    K = NDTMAX
    L = 1
!     Here avoid use DO structure since NXTLIMIT is not constant
    1 IF(L <= NXTLIMIT) THEN
        J = NXTLST(L)
        10 IF(L > NDTK(K)) THEN
            K = K - 1
            GO TO 10
        END IF
        IF(J > I) THEN
            NXTLST(L) = NXTLST(L) - 1
        !     --07/10/14 13:03-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$         IF(J.GE.9951.and.j.LE.9960) print*,'J',J,'NDTK(K)'
        !$$$     &        ,NDTK(K),'K',K,'L',L,'NXTLST(L)',NXTLST(L),
        !$$$     &        'NDTK(K+1)',NDTK(K+1)
        !     --07/10/14 13:03-lwang-end----------------------------------------*
        ELSE IF(J == I) THEN
        !     decrease nxtlst ending point
        !     --07/08/14 16:49-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            print*,'REMOVE I',I,'NXTLIMIT',NXTLIMIT,'L',L,'K',K,
        !$$$     &           'NDTMIN',NDTMIN,'NDTK(K)',NDTK(K),'NDTMAX',NDTMAX,
        !$$$     &           'B_LFAG',B_FLAG
        !     --07/08/14 16:49-lwang-end----------------------------------------*

            IF(NXTLIMIT <= 1) THEN
                write(6,*) 'Error: No particle in NXTLST!'
                call flush(6)
                call abort()
            END IF
            NXTLIMIT = NXTLIMIT - 1
        !     Reduce all NDTK outside NDTMIN by one
            DO LK = 1,NDTMIN-1
                NDTK(LK) = NDTK(LK) - 1
            END DO
        !     Replace removing particle index by the end index of step level K
            IF(L /= NDTK(K)) THEN
                NXTLST(L) = NXTLST(NDTK(K))
                IF(NXTLST(L) > I) NXTLST(L) = NXTLST(L) - 1
            !     Avoid reduce index in position NXTLIMIT
            ELSE IF(K > NDTMIN .AND. NXTLST(NDTK(K-1)) > I) THEN
                NXTLST(NDTK(K-1)) = NXTLST(NDTK(K-1)) - 1
            END IF
            NDTK(K) = NDTK(K) - 1
            DO LK = K-1,NDTMIN,-1
            !     Shift last index position to beginning of level L
                NXTLST(NDTK(LK+1)+1) = NXTLST(NDTK(LK))
            !     Reduce step level L position by one
                NDTK(LK) = NDTK(LK) - 1
            END DO
            RM_FLAG = .true.
        END IF
        L = L + 1
        GO TO 1
    END IF

!     For ghost list
    IF(NGHOSTS > 0) THEN
        IF(RM_FLAG) NXTLST(NXTLIMIT+1) = NXTLST(NXTLIMIT+1+NGHOSTS)
        L = NXTLIMIT+1
        2 IF(NXTLST(L) == I) THEN
            IF(RM_FLAG) THEN
                write(6,*) 'Error: Particle I',I,'exist in both NXTLST ', &
                'and Ghost list! L',L
                call flush(6)
                call abort()
            END IF
            NXTLST(L) = NXTLST(NXTLIMIT+NGHOSTS)
            NGHOSTS = NGHOSTS - 1
            RM_FLAG = .true.
        END IF
        IF(NXTLST(L) > I) NXTLST(L) = NXTLST(L) - 1
        L = L + 1
        IF(L <= NXTLIMIT+NGHOSTS) GO TO 2
    END IF

!     For NLSTDELAY
    IF(NLSTDELAY(1) > 0) THEN
        L = 2
        3 IF(NLSTDELAY(L) == I) THEN
            IF(RM_FLAG) THEN
                write(6,*) 'Error: Particle I',I,'exist in both NXTLST ', &
                'and NLSTDELAY! L',L
                call flush(6)
                call abort()
            END IF
            NLSTDELAY(L) = NLSTDELAY(NLSTDELAY(1)+1)
            NLSTDELAY(1) = NLSTDELAY(1) - 1
            RM_FLAG = .true.
        END IF
        IF(NLSTDELAY(L) > I) NLSTDELAY(L) = NLSTDELAY(L) - 1
        L = L + 1
        IF(L <= NLSTDELAY(1)+1) GO TO 3
    END IF

!     Shift all index by 2 due to KS components removed
    IF(B_FLAG) THEN
        DO LK = 1, NXTLIMIT+NGHOSTS
            NXTLST(LK) = NXTLST(LK) - 2
        END DO
        IF(NLSTDELAY(1) > 0) THEN
            DO LK = 2, NLSTDELAY(1)+1
                NLSTDELAY(LK) = NLSTDELAY(LK) - 2
            END DO
        END IF
    END IF

    IF( .NOT. RM_FLAG) write(6,*) 'Warning!: particle',I,'not find in ', &
    'escape remove'
!     --07/10/14 13:13-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      print*,'-----------------'
!     --07/10/14 13:13-lwang-end----------------------------------------*
    call shrink_tlist

    RETURN

    end subroutine rmesc_tlist
