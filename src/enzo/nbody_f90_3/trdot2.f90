!**
    SUBROUTINE TRDOT2(KW,AGE,TM,TN,TSCLS,DTM,DTR)


!       Time-scale for evolution changes.
!       ---------------------------------

    REAL*8 :: AGE,TM,TN,TSCLS(20),DTM,DTR
    REAL*8 :: PTS1,PTS2
    PARAMETER(PTS1=0.05D0,PTS2=0.02D0)

!       Base new time scale for changes in radius & mass on stellar type.
    if(kw <= 1)then
        dtm = pts1*tm
        dtr = tm - age
    elseif(kw >= 10)then
        dtm = 1.0d+02
        dtr = dtm
    elseif(kw == 2)then
        dtm = pts1*(tscls(1) - tm)
        dtr = tscls(1) - age
    elseif(kw == 3)then
        if(age < tscls(6))then
            dtm = pts2*(tscls(4) - age)
        else
            dtm = pts2*(tscls(5) - age)
        endif
        dtr = MIN(tscls(2),tn) - age
    elseif(kw == 4)then
        dtm = pts1*tscls(3)
        dtr = MIN(tn,tscls(2) + tscls(3)) - age
    elseif(kw == 5)then
        if(age < tscls(9))then
            dtm = pts2*(tscls(7) - age)
        else
            dtm = pts2*(tscls(8) - age)
        endif
        dtr = MIN(tn,tscls(13)) - age
    elseif(kw == 6)then
        if(age < tscls(12))then
            dtm = pts2*(tscls(10) - age)
        else
            dtm = pts2*(tscls(11) - age)
        endif
        dtm = MIN(dtm,0.005d0)
        dtr = tn - age
    elseif(kw == 7)then
        dtm = pts1*tm
        dtr = tm - age
    elseif(kw == 8 .OR. kw == 9)then
        if(age < tscls(6))then
            dtm = pts2*(tscls(4) - age)
        else
            dtm = pts2*(tscls(5) - age)
        endif
        dtr = tn - age
    endif
    dtm = MIN(dtm,dtr)
    dtm = MAX(dtm,1.0d-07)

    RETURN
    END SUBROUTINE TRDOT2
!**
