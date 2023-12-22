
!****************************************************************

    SUBROUTINE INDEXX(N,ARRIN,INDX)
!--
!--            INDEXX sorts ARRIN with increasing order
!--            ARRIN is saved during sorting
!--            ARRIN(INDX(I)) is in the i_th position
!--              of the sorted list
!--
!--                  subroutine from NUMERICAL RECIPES
!--
    DOUBLE PRECISION :: ARRIN(N), Q
    INTEGER :: INDX(N),INDXT,IR,L,I,J
    DO 11 J=1,N
        INDX(J)=J
    11 END DO
    L=N/2+1
    IR=N
    10 CONTINUE
    IF(L > 1)THEN
        L=L-1
        INDXT=INDX(L)
        Q=ARRIN(INDXT)
    ELSE
        INDXT=INDX(IR)
        Q=ARRIN(INDXT)
        INDX(IR)=INDX(1)
        IR=IR-1
        IF(IR == 1)THEN
            INDX(1)=INDXT
            RETURN
        ENDIF
    ENDIF
    I=L
    J=L+L
    20 IF(J <= IR)THEN
        IF(J < IR)THEN
            IF(ARRIN(INDX(J)) < ARRIN(INDX(J+1)))J=J+1
        ENDIF
        IF(Q < ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
        ELSE
            J=IR+1
        ENDIF
        GO TO 20
    ENDIF
    INDX(I)=INDXT
    GO TO 10
    END SUBROUTINE INDEXX
