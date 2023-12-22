    SUBROUTINE NBPOT(NB,NP,POTS)


!       Potential energy of subsystem.
!       ------------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          


!       Obtain potential energy of subsystem JLIST(NB) & JPERT(NP).
    POTS = 0.0D0

    DO K = 1,NP
        J = JPERT(K)
        call jpred(J,time,time)
    END DO
    DO 10 L = 1,NB
        I = JLIST(L)
        call jpred(I,time,time)
        DO 5 K = 1,NP
            J = JPERT(K)
            IF (J > N) THEN
            !       Skip outer binary c.m. or components during merger correction.
                IF (J == I) GO TO 5
                JP = J - N
                IF (I < IFIRST) THEN
                    IF (JP == KVEC(I)) GO TO 5
                END IF
            !       Resolve c.m. of other perturbed KS pairs.
                IF (LIST(1,2*JP-1) > 0) THEN
                    J = 2*JP - 1
                END IF
            END IF
            1 RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 + &
            (X(3,I) - X(3,J))**2
            POTS = POTS + BODY(I)*BODY(J)/SQRT(RIJ2)
        
        !       Check for second component of perturbed KS pair.
            IF (JPERT(K) > N) THEN
                IF (J == 2*JP - 1) THEN
                    J = J + 1
                    GO TO 1
                END IF
            END IF
        5 END DO
    10 END DO

!       Include any external potentials.
    IF (KZ(14) > 0) THEN
        DO 20 L = 1,NB
            I = JLIST(L)
            CALL XTRNLV(I,I)
            POTS = POTS - HT
        !       Note positive sign convention for potential energy.
        20 END DO
    END IF

    RETURN

    END SUBROUTINE NBPOT
