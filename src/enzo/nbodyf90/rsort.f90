    SUBROUTINE RSORT(R,SWITCH,IND)


!       Sorting of mutual distances.
!       ----------------------------

    REAL*8 ::  R(6)
    INTEGER ::  IND(6)
    LOGICAL ::  SWITCH


!       Sort particle separations in increasing order.
    1 NSW = 0
    DO 2 K = 1,5
        IF (R(IND(K+1)) >= R(IND(K))) GO TO 2
        IND1 = IND(K)
        IND(K) = IND(K+1)
        IND(K+1) = IND1
        NSW = 1
    2 END DO

    IF (NSW > 0) GO TO 1

!       Check the switching conditions.
    SWITCH = .FALSE.
    IF (IND(1) > 3 .OR. IND(2) == 6) SWITCH = .TRUE. 

    RETURN

    END SUBROUTINE RSORT
