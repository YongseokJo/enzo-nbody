    SUBROUTINE VERIFY


!       Input validation.
!       -----------------

    USE POINTERS
    INCLUDE 'common6.h'
          


!       Check for unreasonable input parameters (initial & restart).
    IF (N >= NMAX - 2 .OR. NNBMAX > LMAX - 3 .OR. NNBOPT > NNBMAX) THEN
        if(rank == 0) &
        WRITE (6,10)  N, NNBMAX, NNBOPT
        10 FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I6, &
        ' NNBMAX =',I4' NNBOPT =',I4)
        STOP
    END IF

    IF (ETAI > 0.08 .OR. ETAR > 0.16) THEN
        if(rank == 0) &
        WRITE (6,20)  ETAI, ETAR
        20 FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   ETAI =',F6.2, &
        '  ETAR =',F6.2)
        STOP
    END IF

    IF (ETAU > 0.5 .OR. GMIN > 0.0001 .OR. GMAX > 0.10) THEN
        if(rank == 0) &
        WRITE (6,30)  ETAU, GMIN, GMAX
        30 FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   ETAU =',F6.2, &
        '  GMIN =',F11.7,'  GMAX =',F7.3)
        STOP
    END IF

!       Also check for zero or negative values.
    IF (N <= 0 .OR. NNBMAX <= 0 .OR. ETAI <= 0.0 .OR. ETAR <= 0.0) THEN
        if(rank == 0) &
        WRITE (6,40)  N, NNBMAX, ETAI, ETAR
        40 FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I5, &
        '  NNBMAX =',I4,'  ETAI =',F6.2,'  ETAR =',F6.2)
        STOP
    END IF

    IF (ETAU <= 0.0 .OR. GMIN <= 0.0 .OR. GMAX <= 0.0) THEN
        if(rank == 0) &
        WRITE (6,30)  ETAU, GMIN, GMAX
        STOP
    END IF

    IF (DTADJ <= 0.0 .OR. DELTAT <= 0.0 .OR. QE <= 0.0) THEN
        if(rank == 0) &
        WRITE (6,50)  DTADJ, DELTAT, QE
        50 FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   DTADJ =',F6.2, &
        '  DELTAT =',F6.2,'  QE =',1PE9.1)
        STOP
    END IF

    RETURN

    END SUBROUTINE VERIFY
