    SUBROUTINE SORT1(N,RA,RB)


!     Heapsort method (Press p. 231).
!     -------------------------------

    INTEGER ::  RB(N),RRB
    REAL*8 ::  RA(N),RRA


    IF (N <= 1) RETURN
!     !     bug fix Nov 2007.
    L = N/2+1
    IR=N
    10 CONTINUE
    IF (L > 1) THEN
        L=L-1
        RRA=RA(L)
        RRB=RB(L)
    ELSE
        RRA=RA(IR)
        RRB=RB(IR)
        RA(IR)=RA(1)
        RB(IR)=RB(1)
        IR=IR-1
        IF (IR == 1) THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
        END IF
    END IF
    I=L
    J=L+L
    20 IF (J <= IR) THEN
        IF (J < IR) THEN
            IF (RA(J) < RA(J+1)) J=J+1
        END IF
        IF (RRA < RA(J)) THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
        ELSE
            J=IR+1
        END IF
        GO TO 20
    END IF
    RA(I)=RRA
    RB(I)=RRB
    GO TO 10

    END SUBROUTINE SORT1

!*****-------------------------------------

    SUBROUTINE SORT1F(N,RA,RB)


!     Heapsort method (Press p. 231).
!     -------------------------------

    INTEGER ::  RB(N),RRB
    REAL ::  RA(N),RRA


    IF (N <= 1) RETURN
!     !     bug fix Nov 2007.
    L = N/2+1
    IR=N
    10 CONTINUE
    IF (L > 1) THEN
        L=L-1
        RRA=RA(L)
        RRB=RB(L)
    ELSE
        RRA=RA(IR)
        RRB=RB(IR)
        RA(IR)=RA(1)
        RB(IR)=RB(1)
        IR=IR-1
        IF (IR == 1) THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
        END IF
    END IF
    I=L
    J=L+L
    20 IF (J <= IR) THEN
        IF (J < IR) THEN
            IF (RA(J) < RA(J+1)) J=J+1
        END IF
        IF (RRA < RA(J)) THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
        ELSE
            J=IR+1
        END IF
        GO TO 20
    END IF
    RA(I)=RRA
    RB(I)=RRB
    GO TO 10

    END SUBROUTINE SORT1F
