    subroutine ksparmpi(operator,par_type,par_index, &
    index2,index3,par_value)


!     Save and communicate variables inside parallel KS integration

    #ifdef PARALLEL
    Include 'kspars.h'
    Include 'common6.h'
    COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5), &
    NAMES(NCMAX,5),ISYS(5)
    COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX), &
    HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX), &
    NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
    COMMON/MODES/   EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX), &
    BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX), &
    RP(NTMAX),ES(NTMAX),CMM(2,NTMAX),IOSC(NTMAX), &
    NAMEC(NTMAX)

    INTEGER :: LISTILEN,LISTRLEN,par_type,par_index,index2,index3
    INTEGER :: LISTILEN_MPI(maxpe),LISTRLEN_MPI(maxpe)
    INTEGER :: LIDISPLACE(maxpe),LRDISPLACE(maxpe)
    INTEGER :: OPERATOR
    INTEGER :: Ibegin,Iend,Icounter
    INTEGER :: I_flag,I_index,I_index2,I_index3
    INTEGER :: I_LAST(3)
    REAL*8 :: par_value,par_last
    INTEGER :: IKSMPI(2*KMAX+10),IKSMERGE(2*KMAX+10)
    REAL*8 ::  RKSMPI(2*KMAX+10),RKSMERGE(2*KMAX+10)
    REAL*8 ::  EMERGE_MPI,BE3_MPI,ETIDE_MPI,ECOLL_MPI,EGRAV_MPI,E10_MPI
    REAL*8 ::  ECDOT_MPI,EKICK_MPI
    REAL*8 ::  EMERGEOLD,BEOLD,ETIDEOLD,ECOLLOLD,EGRAVOLD,EOLD10
    REAL*8 ::  ECDOTOLD,EKICKOLD
    SAVE IKSMPI,RKSMPI,LISTILEN,LISTRLEN
    SAVE EMERGEOLD,BEOLD,ETIDEOLD,ECOLLOLD,EGRAVOLD,EOLD10
    SAVE ECDOTOLD,EKICKOLD,I_LAST,par_last
    DATA I_LAST /3*-1/

!     --02/24/14 19:03-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      if(operator.eq.K_store) then
!$$$         print*,'Rank:',rank,'CALL KSPARMPI: TYPE',par_type,'Index ',
!$$$     &        par_index,index2,index3,'Value ',par_value,'t',time
!$$$         call flush(6)
!$$$      elseif(operator.eq.K_comm)then
!$$$         print*,'Rank:',rank,'CALL KSPARMPI COMM T:',time
!$$$         call flush(6)
!$$$      else
!$$$         print*,'Rank:',rank,'CALL KSPARMPI RESET T:',time
!$$$         call flush(6)
!$$$      end if
!     --02/24/14 19:03-lwang-end----------------------------------------*
!     save data
    IF (operator == K_store) then
    !     Check whether the data is dup.
        if(par_index == I_Last(1) .AND. index2 == I_last(2) &
         .AND. index3 == I_Last(3) .AND. par_last == par_value) &
        return
    !     Backup data
        I_LAST(1) = par_index
        I_LAST(2:3) = 0
        par_last = par_value
    !     integer variable
        IF (par_type == K_int) then
        !     save variable index
            LISTILEN = LISTILEN + 1
            IKSMPI(LISTILEN) = par_index
        !     check one dimensional array
            IF(par_index >= IISPLIT) then
                LISTILEN = LISTILEN + 1
                IKSMPI(LISTILEN) = index2
                I_LAST(2) = index2
            end if
        !     save value
            LISTILEN = LISTILEN + 1
            IKSMPI(LISTILEN) = int(par_value)
        !     real*8 variable
        ELSEIF (par_type == K_real8) then
        !     real*8 variable
            LISTRLEN = LISTRLEN + 1
        !     --04/08/14 11:39-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            if(LISTRLEN.GE.2*KMAX+10) then
        !$$$               write(130+rank,*) 'LEN',LISTRLEN,'RKS',RKSMPI(1:LISTRLEN)
        !$$$               call flush(130+rank)
        !$$$            end if
        !     --04/08/14 11:39-lwang-end----------------------------------------*
            RKSMPI(LISTRLEN) = DFLOAT(par_index)
        !     check one dimensional array
            IF(par_index >= IRSPLIT) then
                LISTRLEN = LISTRLEN + 1
                RKSMPI(LISTRLEN) = DFLOAT(index2)
                I_LAST(2) = index2
            !     check two dimensional array
                IF(par_index >= IRSPLIT2) then
                    LISTRLEN = LISTRLEN + 1
                    RKSMPI(LISTRLEN) = DFLOAT(index3)
                    I_LAST(3) = index3
                END IF
            END IF
        !     save value
            LISTRLEN = LISTRLEN + 1
            RKSMPI(LISTRLEN) = par_value
        END IF

    !     communication
    ELSEIF (operator == K_comm) then
    !     Gather all integer data
        call MPI_ALLGATHER(LISTILEN,1,MPI_INTEGER,LISTILEN_MPI(1), &
        1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(LISTRLEN,1,MPI_INTEGER,LISTRLEN_MPI(1), &
        1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        LIDISPLACE(1) = 0
        LRDISPLACE(1) = 0
        DO K=1,isize-1
            LIDISPLACE(K+1) = LIDISPLACE(K) + LISTILEN_MPI(K)
            LRDISPLACE(K+1) = LRDISPLACE(K) + LISTRLEN_MPI(K)
        END DO
    !     Integer part:
        IF(LIDISPLACE(isize) >= isize .OR. LISTILEN_MPI(isize) > 1) THEN
            IKSMPI(1) = LISTILEN
            CALL MPI_ALLGATHERV(IKSMPI(1),LISTILEN,MPI_INTEGER, &
            IKSMERGE(1),LISTILEN_MPI,LIDISPLACE,MPI_INTEGER, &
            MPI_COMM_WORLD,ierr)
        !     --02/25/14 16:39-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            LASTI = LIDISPLACE(isize)+LISTILEN_MPI(isize)
        !$$$            print*,'rank',rank,'T',time,'IKSMPI',IKSMERGE(1:LASTI)
        !$$$            call flush(6)
        !     --02/25/14 16:39-lwang-end----------------------------------------*
        !     Initial the flag and index
        !     I_flag: 0: read index, 1: read value, 2: read index2, 3: read index3
            I_FLAG = 0
            I_index = 0
            I_index2 = 0
            I_index3 = 0
        !     For MPI processors loop
            Icounter = 0
            IEND = 0
        !     LOOP
            1 IBEGIN = IEND + 2
            IEND = IKSMERGE(IEND+1) + IEND
            Icounter = Icounter + 1
            DO I=IBEGIN, IEND
                IF(I_flag == 0) then
                    I_index = IKSMERGE(I)
                    IF(I_index >= IISPLIT) then
                        I_FLAG = 2
                    ELSE
                        I_FLAG = 1
                    END IF
                ELSEIF(I_flag == 1) then
                    if(I_index == K_DELAY) then
                        call delay_bcast(icounter-1)
                    elseif(I_index == K_JCLOSE) then
                        JCLOSE = IKSMERGE(I)
                    elseif(I_index == K_JCMAX) then
                        JCMAX = IKSMERGE(I)
                    elseif(I_index == K_IQCOLL) then
                        IQCOLL = IKSMERGE(I)
                    elseif(I_index == K_IPHASE) then
                        IPHASE = IKSMERGE(I)
                    elseif(I_index == K_KSPAIR) then
                        KSPAIR = IKSMERGE(I)
                    elseif(I_index == K_KICK) then
                        call kick(icounter-1,-147,0,0.0)
                    elseif(I_index == K_KSTAR) then
                        KSTAR(I_index2) = IKSMERGE(I)
                    else
                        PRINT*,'Error! No Integer index ',I_index, &
                        ' For KS MPI communication'
                        call flush(6)
                        call abort()
                    end if
                    I_flag = 0
                ELSEIF(I_flag == 2) then
                    I_index2 = IKSMERGE(I)
                    I_FLAG = 1
                ELSE
                    PRINT*,'Error flag for reading data in ksparmpi (' &
                    ,I_FLAG,')!'
                    call flush(6)
                    call abort()
                END IF
            !     --02/25/14 17:00-lwang-debug--------------------------------------*
            !**** Note:------------------------------------------------------------**
            !$$$               print*,'rank',rank,'I',i,'FLAG',i_flag,'index',I_index,
            !$$$     &              I_index2,I_index3,'counter',icounter
            !$$$               call flush(6)
            !     --02/25/14 17:00-lwang-end----------------------------------------*
            END DO
            IF (Icounter < isize) GO TO 1
        END IF
    !     REAL*8 part
        IF(LRDISPLACE(isize) >= isize .OR. LISTRLEN_MPI(isize) > 1) THEN
            RKSMPI(1) = DFLOAT(LISTRLEN)
            CALL MPI_ALLGATHERV(RKSMPI(1),LISTRLEN,MPI_REAL8, &
            RKSMERGE(1),LISTRLEN_MPI,LRDISPLACE,MPI_REAL8, &
            MPI_COMM_WORLD,ierr)
        !     --02/25/14 16:39-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            LASTI = LRDISPLACE(isize)+LISTRLEN_MPI(isize)
        !$$$            print*,'rank',rank,'T',time,'RKSMPI',RKSMPI(1:LISTRLEN)
        !$$$            print*,'rank',rank,'T',time,'RKSMERGE',RKSMERGE(1:LASTI),
        !$$$     &           'IRSPLIT',IRSPLIT
        !$$$            call flush(6)
        !     --02/25/14 16:39-lwang-end----------------------------------------*
        !     Initial energy
            EMERGE_MPI = 0
            BE3_MPI = 0
            ETIDE_MPI = 0
            ECOLL_MPI = 0
            EGRAV_MPI = 0
            E10_MPI = 0
            ECDOT_MPI = 0
            EKICK_MPI = 0
                        
        !     Initial the flag and index
        !     I_flag: 0: read index, 1: read value, 2: read index2, 3: read index3
            I_FLAG = 0
            I_index = 0
            I_index2 = 0
            I_index3 = 0
        !     For MPI processors loop
            Icounter = 0
            IEND = 0
        !     LOOP
            2 IBEGIN = IEND + 2
            IEND = INT(RKSMERGE(IEND+1)) + IEND
            Icounter = Icounter + 1
            DO I=IBEGIN, IEND
                IF(I_flag == 0) then
                    I_index = Int(RKSMERGE(I))
                    IF(I_index >= IRSPLIT) then
                        I_FLAG = 2
                        IF(I_index >= IRSPLIT2) I_FLAG = 3
                    ELSE
                        I_FLAG = 1
                    END IF
                ELSEIF(I_flag == 1) then
                    if(I_index == K_EMERGE) then
                        EMERGE_MPI = EMERGE_MPI+RKSMERGE(I)
                    elseif(I_index == K_BE3) then
                        BE3_MPI = BE3_MPI+RKSMERGE(I)
                    elseif(I_index == K_ETIDE) then
                        ETIDE_MPI = ETIDE_MPI+RKSMERGE(I)
                    elseif(I_index == K_ECOLL) then
                        ECOLL_MPI = ECOLL_MPI+RKSMERGE(I)
                    elseif(I_index == K_EGRAV) then
                        EGRAV_MPI = EGRAV_MPI+RKSMERGE(I)
                    elseif(I_index == K_E10) then
                        E10_MPI = E10_MPI+RKSMERGE(I)
                    elseif(I_index == K_EBCH0) then
                        EBCH0 = RKSMERGE(I)
                    elseif(I_index == K_ECDOT) then
                        ECDOT_MPI = ECDOT_MPI+RKSMERGE(I)
                    elseif(I_index == K_EKICK) then
                        EKICK_MPI = EKICK_MPI+RKSMERGE(I)
                    elseif(I_index == K_RMAX) then
                        RMAX = MAX(RMAX,RKSMERGE(I))
                    elseif(I_index == K_TEV) then
                        TEV(I_index2) = RKSMERGE(I)
                    elseif(I_index == K_RADIUS) then
                        RADIUS(I_index2) = RKSMERGE(I)
                    elseif(I_index == K_SPIN) then
                        SPIN(I_index2) = RKSMERGE(I)
                    elseif(I_index == K_STEPS) then
                        STEPS(I_index2) = RKSMERGE(I)
                    elseif(I_index == K_TMDIS) then
                        TMDIS(I_index2) = RKSMERGE(I)
                    elseif(I_index == K_X0DOT) then
                        X0DOT(I_index3,I_index2) = RKSMERGE(I)
                    elseif(I_index == K_CM_B) then
                        CM(I_index3,I_index2) = RKSMERGE(I)
                    elseif(I_index == K_CM_M) then
                        CMM(I_index3,I_index2) = RKSMERGE(I)
                    else
                        PRINT*,'Error! No Integer index ',I_index, &
                        ' For KS MPI communication'
                        call flush(6)
                        CALL abort()
                    end if
                    I_flag = 0
                ELSEIF(I_flag == 2) then
                    I_index2 = int(RKSMERGE(I))
                    I_FLAG = 1
                ELSEIF(I_flag == 3) then
                    I_index3 = int(RKSMERGE(I))
                    I_FLAG = 2
                ELSE
                    PRINT*,'Error flag for reading data in ksparmpi (' &
                    ,I_FLAG,')!'
                    call flush(6)
                    CALL abort()
                END IF
            !     --02/25/14 17:00-lwang-debug--------------------------------------*
            !**** Note:------------------------------------------------------------**
            !$$$               print*,'rank',rank,'I',i,'FLAG',i_flag,'index',I_index,
            !$$$     &              I_index2,I_index3,'counter',icounter
            !$$$               call flush(6)
            !     --02/25/14 17:00-lwang-end----------------------------------------*
            END DO
            IF (Icounter < isize) GO TO 2
        !     Update energy
            IF(EMERGE_MPI /= 0) EMERGE = EMERGEOLD + EMERGE_MPI
        !     --03/05/14 21:26-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            if(EMERGE_MPI.GE.5) print*,'Warning! EMERGE_MPI Error'
        !     --03/05/14 21:26-lwang-end----------------------------------------*
            IF(BE3_MPI /= 0) BE(3) = BEOLD + BE3_MPI
            IF(ETIDE_MPI /= 0) ETIDE = ETIDEOLD + ETIDE_MPI
            IF(ECOLL_MPI /= 0) ECOLL = ECOLLOLD + ECOLL_MPI
            IF(EGRAV_MPI /= 0) EGRAV = EGRAVOLD + EGRAV_MPI
            IF(E10_MPI /= 0) E(10) = EOLD10 + E10_MPI
            IF(ECDOT_MPI /= 0) ECDOT = ECDOTOLD + ECDOT_MPI
            IF(EKICK_MPI /= 0) EKICK = EKICKOLD + EKICK_MPI
        END IF
    !     RESET
    ELSEIF (operator == K_reset) then
        LISTILEN = 1
        LISTRLEN = 1
        EMERGEOLD = EMERGE
        BEOLD = BE(3)
        ETIDEOLD = ETIDE
        ECOLLOLD = ECOLL
        EGRAVOLD = EGRAV
        EOLD10 = E(10)
        ECDOTOLD = ECDOT
        EKICKOLD = EKICK
        I_LAST(1:3) = -1
    END IF
    #endif

    RETURN
          
    end subroutine ksparmpi

!*************************************************************************
!     For delay
    subroutine delay_bcast(iproc)


!     bcast delay common variables

    #ifdef PARALLEL
    include 'mpi_base.h'
    COMMON/SAVEIT/ IPH0,KS0,KCH0,KS20,JC0,JCL0,PCR0
    integer :: iproc

    CALL MPI_BCAST(IPH0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(KS0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(KCH0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(KS20,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(JC0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(JCL0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(PCR0,1,MPI_REAL8,iproc,MPI_COMM_WORLD,ierr)
    #endif
    RETURN

    end subroutine delay_bcast
          
