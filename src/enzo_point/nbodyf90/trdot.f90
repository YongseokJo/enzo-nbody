    SUBROUTINE TRDOT(I,DTM,M1)


!       Time-scale for expansion of radius.
!       -----------------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    REAL*8 :: TSCLS(20),LUMS(10),GB(10),TM,TN
    REAL*8 :: M0,M1,RM,LUM,AGE,MC,MC1,RCC,RM0,AGE0,M10
    REAL*8 :: menv,renv,k2
    REAL*8 :: pts1,pts2,eps,alpha2,tol
    PARAMETER(pts1=0.05d0,pts2=0.02d0)
    PARAMETER(eps=1.0d-06,alpha2=0.09d0,tol=1.0d-10)

!       Obtain stellar parameters at current epoch (body #I may be ghost).
    KW = KSTAR(I)
    M0 = BODY0(I)*ZMBAR
    IF(M1 <= 0.0) M1 = RADIUS(I)*SU
    M10 = M1
    MC = 0.D0
    AGE = TEV0(I)*TSTAR - EPOCH(I)
    CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
    CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS, &
    RM,LUM,KW,MC,RCC,MENV,RENV,K2)

!       Quit if there is a change of type at the current TEV.
    if((kstar(i) <= 6 .AND. kw > 6) .OR. &
    (kstar(i) <= 9 .AND. kw > 9))then
        m1 = m10
        dtm = 0.d0
        goto 10
    endif
!**

!       Base new time scale for changes in radius & mass on stellar type.
    if(kw <= 1)then
        dtm = pts1*tm
        dtr = tm - age
    elseif(kw >= 10)then
    !         dtm = 1.0d+02
        dtm = sqrt(0.1/lum)
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
    !     --09/27/13 23:37-lwang-improvements-------------------------------*
    !**** Note: dtm did not necessary become small-------------------------**
        dtm = MIN(dtm,0.005d0)
    !     --09/27/13 23:37-lwang-end----------------------------------------*
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

! Record radius.

    rm0 = rm
    if(kw >= 10) goto 30
    age0 = age
    kw0 = kw
    mc1 = mc

! Check for type change.

    it = 0
    if((dtr-dtm) <= tol)then
    
    ! Check final radius for too large a jump.
    
        age = MAX(age,age*(1.d0-eps) + dtr)
        CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS, &
        RM,LUM,KW,MC1,RCC,MENV,RENV,K2)
        dr = rm - rm0
        if(ABS(dr) > 0.1*rm0)then
            dtm = dtr - age0*eps
            dtdr = dtm/ABS(dr)
            dtm = alpha2*MAX(rm,rm0)*dtdr
            goto 20
        else
            dtm = dtr
            goto 30
        endif
    endif

! Limit to a 10% increase assuming no further mass loss
! and thus that the pertubation functions due to small envelope mass
! will not change the radius.

    20 age = age0 + dtm
    mc1 = mc
    CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS, &
    RM,LUM,KW,MC1,RCC,MENV,RENV,K2)
    dr = rm - rm0
    it = it + 1
    if(it == 20 .AND. kw == 4) goto 30
    IF(IT > 30)THEN
        if(rank == 0) WRITE (6,22) IT, KSTAR(I), M0, DR, RM0
        22 FORMAT (' DANGER!    TRDOT: IT K* M0 DR RM0 ',2I4,1P,3E10.2)
        goto 30
    ENDIF
    if(ABS(dr) > 0.1*rm0)then
        dtdr = dtm/ABS(dr)
        dtm = alpha2*MAX(rm0,rm)*dtdr
        if(it >= 20) dtm = 0.5d0*dtm
        goto 20
    endif

    30 continue

!       Impose a lower limit and convert time interval to scaled units.
    DTM = MAX(DTM,1.0D-04)/TSTAR

!       Use dummy for routine CHINIT (new chain or chain collision).
    10 IF(IPHASE == 8 .OR. IPHASE == 9)THEN
        M1 = RM
    ENDIF

    RETURN
    END SUBROUTINE TRDOT
