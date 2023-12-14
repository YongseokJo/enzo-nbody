    SUBROUTINE NBREST(ICM,NSYS,NP)


!       Restore ghosts in neighbour lists.
!       ----------------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          


!       Examine all members inside neighbour radius of body #ICM.
    DO 100 LL = 1,NP
        I = JPERT(LL)
        10 NNB1 = LIST(1,I) + 1
        IF (NNB1 == 1) GO TO 100
    
    !       First see whether body #ICM is a neighbour.
        DO 20 L = 2,NNB1
            IF (LIST(L,I) == ICM) THEN
                GO TO 30
            ELSE
                IF (LIST(L,I) > ICM) GO TO 100
            END IF
        20 END DO
    
    !       Add other members of the subsystem subject to limit.
        30 DO 60 K = 1,NSYS
            J = JLIST(K)
            IF (J == ICM) GO TO 60
        
        !       Skip addition if body #J is already a member.
            DO 40 L = 2,NNB1
                IF (LIST(L,I) == J) GO TO 60
            40 END DO
        
        !       Move members down until correct location identified.
            DO 50 L = NNB1,1,-1
                IF (LIST(L,I) > J .AND. L > 1) THEN
                    LIST(L+1,I) = LIST(L,I)
                ELSE
                !       Place body #J in vacated location and increase membership.
                    LIST(L+1,I) = J
                    LIST(1,I) = LIST(1,I) + 1
                    NNB1 = NNB1 + 1
                !       See whether neighbour list has reached maximum size.
                    IF (NNB1 > NNBMAX) GO TO 70
                    GO TO 60
                END IF
            50 END DO
        60 END DO
    
    !       Also check any KS perturber list.
        70 IF (I > N) THEN
            I = 2*(I - N) - 1
            IF (LIST(1,I) > 0) GO TO 10
        END IF
    100 END DO

    RETURN

    END SUBROUTINE NBREST
