    subroutine encorr_mdot(NPNUM,NDMLIST,DMLIST,GPUDM,GPUSW)


!     Correct the energy due to mass loss
!     -----------------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    INCLUDE 'timing.h'
    INCLUDE 'omp_lib.h'
    INTEGER :: NDMLIST(NMAX),NPNUM
!    &  ,NDCLIST(NMAX/10),NCNUM
    REAL*8 :: DMLIST(NMAX),DMPHI(NMAX)
    #ifdef PARALLEL
    REAL*8 :: DMPHI_L(NMAX)
    #endif
    LOGICAL :: GPUDM,GPUSW
!     &  ,DCLIST(NMAX/10)
          
!     --11/11/13 18:27-lwang-check--------------------------------------*
!**** Note: Energy correction for mass loss stars----------------------**
    if (NPNUM > 0) THEN
        #ifdef NGPU
        call cputim(tttpota)
        NNT = NTOT - IFIRST + 1
    !     Check whether need to send all particles to gpu
        if(GPUSW) then
            call gpupot_send(rank,NNT,BODY(IFIRST),X(1,IFIRST))
            GPUSW = .false.
        end if
    !     Calculation potential
        IF (NPNUM <= 2048) THEN
            IF(GPUDM) THEN
                CALL GPUPOT_DM(rank,NPNUM,NNT,1-IFIRST,NDMLIST(1), &
                DMLIST(1),DMPHI(1))
            ELSE
                CALL GPUPOT_FLOAT(rank,NPNUM,NNT,1-IFIRST,NDMLIST(1), &
                DMLIST(1),DMPHI(1))
            END IF
        !$$$               if(rank.eq.0) print*, 'GPU:',DMPHI(1:NPNUM)
        !$$$            CALL GPUPOT(rank,NDMLIST(1),NDMLIST(1),NNT,BODY(IFIRST),
        !$$$     *           X(1,IFIRST),DMPHI(1))
        !$$$               if(rank.eq.0) print*, 'GPU2:',DMPHI(1:NPNUM)
        else
            DO J=1,NPNUM,2048
                NN = 2048
                IF (NPNUM-J+1 < 2048) NN = NPNUM-J+1
                IF(GPUDM) THEN
                    CALL GPUPOT_DM(rank,NN,NNT,1-IFIRST,NDMLIST(J), &
                    DMLIST(J),DMPHI(J))
                ELSE
                    CALL GPUPOT_FLOAT(rank,NN,NNT,1-IFIRST,NDMLIST(J), &
                    DMLIST(J),DMPHI(J))
                END IF
            END DO
        end if
        call cputim(tttpotb)
        ttpotg = ttpotg + (tttpotb-tttpota)*60
        #else
        call cputim(tttpota)
        DMPHI(1:NPNUM) = 0.0
        #ifdef PARALLEL
        i_tot = NTOT - IFIRST + 1
        if(i_tot*NPNUM < 10000*icore*isize) then
            #endif
            DO III = 1, NPNUM
                DMPHI_T = 0.D0
            ! omp parallel do private(J,L,RIJ2)
            ! omp& reduction(+:DMPHI_T)
                DO J = IFIRST,NTOT
                    L = NDMLIST(III)
                    IF (J /= L) THEN
                        RIJ2 = 0.0D0
                        DO  K = 1,3
                            RIJ2 = RIJ2 + (X(K,L) - X(K,J))**2
                        END DO
                        DMPHI_T = DMPHI_T + BODY(J)/SQRT(RIJ2)
                    END IF
                END DO
            ! omp end parallel do
                DMPHI(III) = DMPHI_T
            END DO
            #ifdef PARALLEL
        else
            DMPHI_L(1:NPNUM) = 0.0
            i_block = i_tot/isize
            if(mod(i_tot,isize) /= 0) i_block = i_block + 1
            i_start = rank*i_block + ifirst
            i_end = (rank+1)*i_block + ifirst - 1
            if(rank == isize-1) i_end = NTOT
            DO III = 1, NPNUM
                DMPHI_T = 0.D0
            ! omp parallel do private(J,L,RIJ2) reduction(+:DMPHI_T)
                DO J = i_start,i_end
                    L = NDMLIST(III)
                    IF (J /= L) THEN
                        RIJ2 = 0.0D0
                        DO  K = 1,3
                            RIJ2 = RIJ2 + (X(K,L) - X(K,J))**2
                        END DO
                        DMPHI_T = DMPHI_T + BODY(J)/SQRT(RIJ2)
                    END IF
                END DO
            ! omp end parallel do
                DMPHI_L(III) = DMPHI_T
            END DO
            CALL MPI_ALLREDUCE(DMPHI_L,DMPHI,NPNUM,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        end if
        #endif
    !$$$            DO IIJ = 1, NPNUM
    !$$$               J = NDMLIST(IIJ)
    !$$$               DO III = 1, NPNUM
    !$$$                  IF (IIJ.NE.III) THEN
    !$$$                     L = NDMLIST(III)
    !$$$c$$$                     IF(J.EQ.L) print*,'Warning!: J.eq.L',J,L,IIJ,III
    !$$$                     RIJ2 = 0.0D0
    !$$$                     DO  K = 1,3
    !$$$                        RIJ2 = RIJ2 + (X(K,L) - X(K,J))**2
    !$$$                     END DO
    !$$$                    DMPHI(IIJ) = DMPHI(IIJ) + 0.5*DMLIST(III)/SQRT(RIJ2)
    !$$$                  END IF
    !$$$               END DO
    !$$$            END DO
    !$$$            IF (NCNUM.GT.0) THEN
    !$$$               DO IIJ = 1, NPNUM
    !$$$                  J = NDMLIST(IIJ)
    !$$$                  DO III = 1, NCNUM
    !$$$                     L = NDCLIST(III)
    !$$$                     IF(J.NE.L) THEN
    !$$$                        RIJ2 = 0.0D0
    !$$$                        DO  K = 1,3
    !$$$                           RIJ2 = RIJ2 + (X(K,L) - X(K,J))**2
    !$$$                        END DO
    !$$$                        DMPHI(IIJ) = DMPHI(IIJ) +
    !$$$     &                       0.5*DCLIST(III)/SQRT(RIJ2)
    !$$$                     END IF
    !$$$                  END DO
    !$$$               END DO
    !$$$            END IF
        call cputim(tttpotb)
        ttpot = ttpot + (tttpotb-tttpota)*60
    !     --05/20/14 22:10-lwang-debug--------------------------------------*
    !**** Note:------------------------------------------------------------**
    !         if(NPNUM.GE.2) then
    !$$$         print*,rank,'NPNUM',NPNUM,'DMPHI',DMPHI(1:NPNUM),'ttpot',ttpot
    !$$$         call flush(6)
    !$$$         stop
    !         end if
    !     --05/20/14 22:10-lwang-end----------------------------------------*
    !$$$            if(rank.eq.0)print*, 'HOST:',NPNUM,DMPHI(1:NPNUM)
    !$$$            call flush(6)
    !$$$            stop
        #endif
        call cputim(tttpota)
    !!$omp parallel do private(J,L,VI2) reduction(+:EMDOT)
        DO J = 1, NPNUM
            L = NDMLIST(J)
            VI2 = XDOT(1,L)**2 + XDOT(2,L)**2 + XDOT(3,L)**2
            EMDOT = EMDOT + DMLIST(J)*(0.5*VI2 - DMPHI(J))
        END DO
    !!$omp end parallel do
    !     --11/15/13 17:06-lwang-debug--------------------------------------*
    !**** Note:------------------------------------------------------------**
    !$$$         if(rank.eq.0) then
    !$$$            print*,'EMDOT',EMDOT,TIME,NDMLIST(1:NPNUM),DMPHI(1),VI2
    !$$$         end if
    !     --11/15/13 17:06-lwang-end----------------------------------------*
        NPNUM = 0
        call cputim(tttpotb)
        ttecor = ttecor + (tttpotb-tttpota)*60
    END IF
!     --11/11/13 18:27-lwang-end----------------------------------------*

    RETURN
          
    end subroutine encorr_mdot
