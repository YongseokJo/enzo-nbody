    SUBROUTINE SUBINT(IQ,I10)
!     label-ks/triple/quad

!       Decision-making for subsystems.
!       -------------------------------
!       R. Spurzem / S.J. Aarseth Experimental Version for parallel binaries
!       Sep 2001/July 2002

    Include 'kspars.h'
    USE POINTERS
    INCLUDE 'common6.h'
          
    INCLUDE 'timing.h'
    COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5), &
    NAMES(NCMAX,5),ISYS(5)
!$$$      COMMON/KSSTAT/ISTEPP,ISTEPU,IBINP,IBINU
    REAL*8 ::  TSLIST(10*KMAX),SBLOCK,DTBL,SPREV
    INTEGER ::  BSLIST(10*KMAX),IPAIR,KBLOCK
    SAVE  IRUN,LI,SPREV,KBLOCK,DTBL
    DATA  IRUN /0/
    DATA ICALL,KBSUM,ISBSUM,KSBSUM /0,0,0,0/
    #ifdef PARALLEL
    COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
    REAL*8 :: TPRED
    LOGICAL :: ipredall
    REAL*8 :: ZMPI(57,KMAX)
    Integer :: inum(maxpe),ista(maxpe),istart,iend
!      INTEGER CMBLIST(2,KMAX/maxpe),CMBLEN
    INTEGER :: JMPI(3,maxpe),JMPI_THR(3),IMPI_THR(LMAX*KMAX),NSTEPUOLD
    INTEGER :: NNBKTOT,NNBK_MPI(maxpe),IMPI_OFF(maxpe),IMPI(LMAX*KMAX)
!     JMPI: 1:NSTEPU; 2:Index of I particle; 3:Iphase of I particle
    REAL*8 :: EMPI(3,maxpe),EMPI_THR(3),EMERGEOLD,BEOLD
!     EMPI: 1:EMERGE; 2:BE(3); 3:RMAX
    LOGICAL :: CMFLAG
    #endif

!       Determine correct index after restart (NNTB = 0 initially).
    IF (IRUN == 0) THEN
        IRUN = 1
        TI = 1.0D+10
    !     Find smallest sum by looping backwards (avoids multiple entries).
        DO K = NNTB,1,-1
            J = KBLIST(K)
            TJ = T0(J) + STEP(J)
            IF (TJ < TI) THEN
                TI = TJ
                LI = K
            ELSE
            !     Adopt the previous index (increased by 1 below).
                LI = LI - 1
                GO TO 1
            END IF
        END DO
    END IF

    1 CONTINUE

!       See whether to advance any KS solutions at start of block-step.
    IF (NPAIRS > 0) THEN
        call cputim(ttks1)
    !       Obtain list of all KS pairs due in interval DTB.
        IF (TBLIST <= TBLOCK .OR. NPAIRS /= NBPREV .OR. &
        NNTB > KMAX) THEN
        !$$$            IF (DTB.EQ.0.0D0) THEN
        !$$$               DTB = MAX(DTMIN,TBLOCK - TPREV)
        !$$$            END IF
            DTBMIN = 1.D30
        !$$$ 2          TBLIST = TPREV + DTB
        !$$$            TBLIST = MAX(TBLOCK,TBLIST)
            TBLIST = TBLOCK
            NNTB = 0
            DO JPAIR = 1,NPAIRS
                J1 = 2*JPAIR - 1
            !$$$               TXLIST(J1) = T0(J1) + STEP(J1)
                IF (T0(J1) + STEP(J1) <= TBLIST) THEN
                    NNTB = NNTB + 1
                    KBLIST(NNTB) = J1
                    TSLIST(NNTB) = T0(J1) + STEP(J1)
                    IF(STEP(J1) < DTBMIN) DTBMIN = STEP(J1)
                END IF
            END DO
        !     Increase interval on zero membership.
            IF (NNTB == 0 .AND. rank == 0) THEN
            !$$$               DTB = 2.0*DTB
            !$$$               print*,'Warning! Non KS detected to be integrate',
            !$$$     &          ' but Npairs =',npairs,' TBlock =',tblock
            !$$$               call flush(6)
                GO TO 20
            END IF
        !     Determine minimum quantized step near binary steps.
            CALL STEPK(DTBMIN,DTNMIN)
        !$$$            ISTEPP=0
        !$$$            ISTEPU=0
        !$$$            IBINP=0
        !$$$            IBINU=0
            ICALL = ICALL + 1
            FAC = 1.0/LOG(1.9999999)
            JMM = INT(1 - LOG(DTNMIN)*FAC)
            DTBL = TBLOCK - TPREV
            JMN = INT(1 - LOG(DTBL)*FAC)
        
            IF(JMM < JMN)THEN
                KBLOCK = 1
            ELSE
                KBLOCK = 2**(JMM-JMN)
            END IF
        
            DTBL = DTBL/KBLOCK
            KBSUM = KBSUM + KBLOCK
        !$$$*     Stabilize interval on membership of 2*SQRT(NPAIRS).
        !$$$            NBTRY = 2*SQRT(FLOAT(NPAIRS))
        !$$$            IF (NNTB.GT.NBTRY)  DTB = 0.75*DTB
        !$$$            IF (NNTB.LT.NBTRY)  DTB = 1.25*DTB
        !     Sort the time-step list sequentially in KBLIST and reset pointer.
            IF (NNTB > 1) THEN
                CALL SORT1(NNTB,TSLIST,KBLIST)
            END IF
        
            LI = 0
        
            SPREV = TPREV
        
        !     --12/21/13 22:29-lwang-print--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$#ifdef PRINT
        !$$$            J1 = KBLIST(1)
        !$$$            JNNTB = KBLIST(NNTB)
        !$$$            Print*,' Kbl ',J1,Jnntb,T0(J1)+STEP(J1),
        !$$$     &           T0(JNNTB)+STEP(JNNTB)
        !$$$            PRINT*,' KBLOCK,SPREV,DTBL,NNTB,TBLIST ',KBLOCK,SPREV,
        !$$$     *           DTBL,NNTB,TBLIST
        !$$$            PRINT*,' 1-NNTB times ',J1,JNNTB,T0(J1)+STEP(J1),
        !$$$     *           T0(JNNTB)+STEP(JNNTB)
        !$$$            CALL FLUSH(6)
        !$$$            PRINT*,'--------------------'
        !$$$            PRINT*,' DTBL,DTBMIN,JMM,JMN=',DTBL,DTBMIN,JMM,JMN
        !$$$
        !$$$            IF(MOD(ICALL,1).EQ.0) PRINT*,' SUBINT: Av KBLOCK=',
        !$$$     *           REAL(KBSUM)/REAL(ICALL)
        !$$$#endif
        !     --12/21/13 22:28-lwang-end----------------------------------------*
        END IF
        call cputim(ttks2)
        ttksblist = ttksblist+(ttks2-ttks1)*60
    

    !     --12/21/13 22:29-lwang-end----------------------------------------*
        DO IX = 1,KBLOCK
        
            SBLOCK = SPREV + DTBL
        !     Form list of any KS pairs due in the current sub-block-step.
            KSB = 0
            DO  LJ = LI+1,NNTB
                J1 = KBLIST(LJ)
                TJ = T0(J1) + STEP(J1)
                IF (TJ <= SBLOCK) THEN
                !$$$                  IF(LIST(1,J1).GT.0)THEN
                !$$$                     IBINP = IBINP + 1
                !$$$                  ELSE
                !$$$                     IBINU = IBINU + 1
                !$$$                  END IF
                    KSB = KSB + 1
                    BSLIST(KSB) = J1
                ELSE
                !     Stop searching if sub-block level is reached.
                    GO TO 65
                END IF
            END DO
        
        !     Continue if no binaries to integrate
            65 IF (KSB == 0) GO TO 9
        
            ISBSUM = ISBSUM + 1
            KSBSUM = KSBSUM + KSB
        !     --03/14/14 12:27-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            if(rank.eq.0) write(6,141),TBlock,TBLIST,
        !$$$     ^           TIME,KBLOCK,IX,SBLOCK,KSB
        !$$$ 141        format('TBLOCK',E15.7,' TBLIST',E15.7,' T',E15.7,' KB',I6,
        !$$$     &           ' IX',I3,' SB',E15.7,' KSB',I6)
        !$$$            call flush(6)
        !$$$            do L= 1,ksb
        !$$$               i1 = bslist(l)
        !$$$               ipair = kvec(i1)
        !$$$               write(110+rank,*) 'i1',i1,'ip',ipair,'n',name(i1),
        !$$$     &              'u',u(1,ipair),'ud',udot(1,ipair),'t0',t0(i1),
        !$$$     &              'step',step(i1),'TBLOCK',tblock
        !$$$               call flush(110+rank)
        !$$$            end do
        !$$$            call mpi_barrier(MPI_COMM_WORLD,ierr)
        !     --03/14/14 12:27-lwang-end----------------------------------------*
        !     --12/21/13 22:29-lwang-print------------------------------------------*
        !**** Note:-----------------------------------------------------------------**
        !$$$  #ifdef PRINT
        !$$$  PRINT*,' Start Do Loop KSB, first min 5=',KSB,
        !$$$  *              ' BSLIST=',(BSLIST(L),L=1,MIN(5,KSB))
        !$$$  PRINT*,' Average KSB=',REAL(KSBSUM)/REAL(ISBSUM)
        !$$$  CALL FLUSH(6)
        !$$$  PRINT*,'--------------------'
        !$$$  #endif
        !     --12/21/13 22:30-lwang-end--------------------------------------------*
        !     Advance binaries due in the current sub-block-step.
            #ifdef PARALLEL
            IF(KSB <= iserks .OR. iserks == 0) THEN
            !     Reset ks par MPI
                call ksparmpi(K_reset,0,0,0,0,0.0)
                #endif
                call cputim(ttks5)
                DO L = 1,KSB
                
                    I1 = BSLIST(L)
                    IPAIR = KVEC(I1)
                    TIME = T0(I1) + STEP(I1)
                    IMULT = 0
                
                    10 continue
                !     --12/21/13 22:29-lwang-print------------------------------------------*
                !**** Note:-----------------------------------------------------------------**
                !$$$  #ifdef PRINT
                !$$$  PRINT*,' L, I1, IPAIR, TIME=',L,I1,IPAIR,
                !$$$  &                    TIME,IMULT,STEP(I1)
                !$$$  CALL FLUSH(6)
                !$$$  #endif
                !     --12/21/13 22:30-lwang-end----------------------------------------*
                
                !$$$*     --12/30/13 9:16-lwang-debug--------------------------------------*
                !$$$***** Note:------------------------------------------------------------**
                !$$$                  if(ipair.eq.1605.and.time.le.29.047283380920465) then
                !$$$                     print*,rank,'single ks t',time
                !$$$                  end if
                !$$$                  if(name(i1).eq.7173) then
                !$$$                  write(102+rank,*)'I1',I1,'IPAIR',ipair,'U',U(1,IPAIR),
                !$$$     &                 'U0',U0(1,IPAIR),'UDOT',UDOT(1,ipair),
                !$$$     &                 'FU',FU(1,ipair),'FT',FUDOT(1,ipair),
                !$$$     &                 'FT2',FUDOT2(1,ipair),'FT3',FUDOT3(1,ipair),
                !$$$     &                 'FP0',FP0(1,Ipair),'FD0',FD0(1,ipair),
                !$$$     &                 'H',H(ipair),'HT',HDOT(ipair),'HT2',HDOT2(ipair),
                !$$$     &                 'HT3',HDOT3(ipair),'HT4',HDOT4(ipair),
                !$$$     &                 'dt',DTAU(ipair),'TD2',TDOT2(ipair),
                !$$$     &                 'TD3',TDOT3(ipair),'R',R(ipair),'R0',r0(ipair),
                !$$$     &                 'Gamma',GAMMA(ipair),'H0',h0(ipair),
                !$$$     &                 'SF',sf(1,ipair),'kslow',kslow(ipair),
                !$$$     &                 'nnb',LIST(1,i1),'list',list(2,i1),
                !$$$     &                 'tblock',tblock,'ksb',ksb,'mass',body(i1),
                !$$$     *                 'r',radius(i1)
                !$$$                  call flush(102+rank)
                !$$$                  if(i1.eq.3327.and.tblock.ge.4.21E-3) then
                !$$$                     print*,rank,'SINGLE KSINT TBLOCK',tblock
                !$$$                     call flush(6)
                !$$$                  end if
                !$$$*     --12/30/13 9:16-lwang-end----------------------------------------*
                    CALL KSINT(I1)
                
                !     --03/17/14 20:22-lwang-debug--------------------------------------*
                !**** Note:------------------------------------------------------------**
                !$$$               if(name(i1).eq.7173) then
                !$$$                  write(102+rank,*)'I1',I1,'IPAIR',ipair,'U',U(1,IPAIR),
                !$$$     &                 'U0',U0(1,IPAIR),'UDOT',UDOT(1,ipair),
                !$$$     &                 'FU',FU(1,ipair),'FT',FUDOT(1,ipair),
                !$$$     &                 'FT2',FUDOT2(1,ipair),'FT3',FUDOT3(1,ipair),
                !$$$     &                 'FP0',FP0(1,Ipair),'FD0',FD0(1,ipair),
                !$$$     &                 'H',H(ipair),'HT',HDOT(ipair),'HT2',HDOT2(ipair),
                !$$$     &                 'HT3',HDOT3(ipair),'HT4',HDOT4(ipair),
                !$$$     &                 'dt',DTAU(ipair),'TD2',TDOT2(ipair),
                !$$$     &                 'TD3',TDOT3(ipair),'R',R(ipair),'R0',r0(ipair),
                !$$$     &                 'Gamma',GAMMA(ipair),'H0',h0(ipair),
                !$$$     &                 'SF',sf(1,ipair),'kslow',kslow(ipair),
                !$$$     &                 'nnb',LIST(1,i1),'list',list(2,i1),
                !$$$     &                 'tblock',tblock,'ksb',ksb,'mass',body(i1),
                !$$$     *                 'r',radius(i1),'T',TIME
                !$$$                  call flush(102+rank)
                !$$$                  if(i1.eq.3327.and.tblock.ge.4.21E-3) then
                !$$$                     print*,rank,'SINGLE KSINT TBLOCK',tblock
                !$$$                     call flush(6)
                !$$$               end if
                !$$$               print*,rank,'sig iq',iq,'t',time
                !$$$               call flush(6)
                !     --03/17/14 20:22-lwang-end----------------------------------------*
                !     --12/21/13 22:49-lwang-print--------------------------------------*
                !**** Note:------------------------------------------------------------**
                !$$$  #ifdef PRINT
                !$$$  PRINT*,' Ret KSINT T0+STEP STEP=',
                !$$$  &                    T0(I1)+STEP(I1),STEP(I1)
                !$$$  CALL FLUSH(6)
                !$$$  #endif
                !     --12/21/13 22:49-lwang-end----------------------------------------*

                !     Check for multiple calls of #I1 (saves using CALL INSERT).
                    IF (IPHASE == 0) THEN
                        TI = TIME + STEP(I1)
                        IF (TI <= SBLOCK) THEN
                            TIME = TI
                            IMULT = 1
                            GO TO 10
                        END IF
                    ELSE
                        I10 = I1
                        KSPAIR = KVEC(I1)
                        IQ = IPHASE
                        IF (IQ < 0) THEN
                            IQCOLL = -2
                            CALL CMBODY(2)
                            TIME = TBLOCK
                            RETURN
                        END IF
                    END IF
                END DO
                call cputim(ttks6)
                ttksints = ttksints+(ttks6-ttks5)*60.
                #ifdef PARALLEL
            !     Parallel KSINT
            ELSE
                call cputim(ttks7)
            !     Reset full prediction flag before enter into ksint.
                ipredall=.false.
            !     Reset ks par MPI
                call ksparmpi(K_reset,0,0,0,0,0.0)

                nl = KSB
            
                inl = nl/isize
                jsize = isize*inl
                idiff = nl - jsize
                irun = 0
            
                do ir = 1,isize
                    inum(ir)=inl
                    if(ir <= idiff)inum(ir) = inum(ir) + 1
                    ista(ir) = irun+1
                    if(ista(ir) > nl)inum(ir) = 0
                    irun = irun + inum(ir)
                end do
            
                istart = ista(rank+1)
                iend = ista(rank+1) + inum(rank+1) - 1
            !     reset cmbody counter to zero
            !$$$  CMBLEN = 1
            !$$$  CMBLIST(1,1) = 0
                CMFLAG = .false.
                NSTEPUOLD = NSTEPU
                EMERGEOLD = EMERGE
                BEOLD = BE(3)
                NNBKTOT = 0
                JMPI_THR(2:3) = 0

            !     --12/29/13 22:31-lwang-debug--------------------------------------*
            !**** Note:------------------------------------------------------------**
            !$$$  print*,'BSLIST ',rank,ksb,istart,iend,time
            !$$$  print*,rank,'B ',BSLIST(1:ksb)
            !$$$  call flush(6)
            !     --12/29/13 22:31-lwang-end----------------------------------------*
                DO L = istart,iend
                
                    I1 = BSLIST(L)
                    IPAIR = KVEC(I1)
                    TIME = T0(I1) + STEP(I1)
                    IMULT = 0
                
                    11 continue
                !     --12/21/13 22:29-lwang-print------------------------------------------*
                !**** Note:-----------------------------------------------------------------**
                !$$$  #ifdef PRINT
                !$$$  PRINT*,' L, I1, IPAIR, TIME=',L,I1,IPAIR,
                !$$$  &                    TIME,IMULT,STEP(I1)
                !$$$  CALL FLUSH(6)
                !$$$  #endif
                !     --12/21/13 22:30-lwang-end----------------------------------------*
                !     --03/14/14 11:43-lwang-debug--------------------------------------**
                !**** Note:------------------------------------------------------------**
                !$$$                  if(ipair.eq.1605.and.time.le.29.047283380920465) then
                !$$$                     print*,rank,'parallel ks t',time
                !$$$                  end if
                !     --03/14/14 11:43-lwang-end----------------------------------------*
                
                    CALL KSINT(I1)
                !     --12/21/13 22:49-lwang-print--------------------------------------*
                !**** Note:------------------------------------------------------------**
                !$$$  #ifdef PRINT
                !$$$  PRINT*,' Ret KSINT T0+STEP STEP=',
                !$$$  &                    T0(I1)+STEP(I1),STEP(I1)
                !$$$  CALL FLUSH(6)
                !$$$  #endif
                !     --12/21/13 22:49-lwang-end----------------------------------------*
                !     Check for multiple calls of #I1 (saves using CALL INSERT).
                    IF (IPHASE == 0) THEN
                        TI = TIME + STEP(I1)
                        IF (TI <= SBLOCK) THEN
                            TIME = TI
                            IMULT = 1
                            GO TO 11
                        END IF
                    ELSE
                        IF ( .NOT. cmflag) THEN
                            JMPI_THR(2) = I1
                            JMPI_THR(3) = IPHASE
                        !     If cmbody is needed, store this pair index and status till end of loop.
                            IF (IPHASE < 0) cmflag = .TRUE. 
                        !$$$  CMBLIST(1,1) = CMBLEN
                        !$$$  CMBLEN = CMBLEN + 1
                        !$$$  CMBLIST(1,CMBLEN) = IPAIR
                        !$$$  CMBSTAT(2,CMBLEN) = IPHASE
                        END IF
                    END IF

                    ZMPI(1:4,L) = U(1:4,IPAIR)
                    ZMPI(5:8,L) = U0(1:4,IPAIR)
                    ZMPI(9:12,L) = UDOT(1:4,IPAIR)
                    ZMPI(13:16,L) = FU(1:4,IPAIR)
                    ZMPI(17:20,L) = FUDOT(1:4,IPAIR)
                    ZMPI(21:24,L) = FUDOT2(1:4,IPAIR)
                    ZMPI(25:28,L) = FUDOT3(1:4,IPAIR)
                    ZMPI(29:32,L) = FP0(1:4,IPAIR)
                    ZMPI(33:36,L) = FD0(1:4,IPAIR)
                    ZMPI(37,L) = H(IPAIR)
                    ZMPI(38,L) = HDOT(IPAIR)
                    ZMPI(39,L) = HDOT2(IPAIR)
                    ZMPI(40,L) = HDOT3(IPAIR)
                    ZMPI(41,L) = HDOT4(IPAIR)
                    ZMPI(42,L) = DTAU(IPAIR)
                    ZMPI(43,L) = TDOT2(IPAIR)
                    ZMPI(44,L) = TDOT3(IPAIR)
                    ZMPI(45,L) = R(IPAIR)
                    ZMPI(46,L) = R0(IPAIR)
                    ZMPI(47,L) = GAMMA(IPAIR)
                    ZMPI(48,L) = H0(IPAIR)
                    ZMPI(49:55,L) = SF(1:7,IPAIR)
                    ZMPI(56,L) = STEP(I1)
                    ZMPI(57,L) = T0(I1)
                    NNBKS = LIST(1,I1) + 1
                    IMPI_THR(NNBKTOT+1:NNBKTOT+NNBKS)=LIST(1:NNBKS,I1)
                    IMPI_THR(NNBKTOT+NNBKS+1)=KSLOW(IPAIR)
                    NNBKTOT = NNBKTOT + NNBKS + 1
                !     --03/29/14 20:29-lwang-debug--------------------------------------*
                !**** Note:------------------------------------------------------------**
                !$$$                  if(tBLOCK.eq.6.05468750000000000E-002) then
                !$$$                     write(130+rank,*) 'I1',I1,NNBKS,
                !$$$     &                    'LIST',LIST(1:NNBKS,I1),'NTOT',NNBKTOT
                !$$$                     call flush(130+rank)
                !$$$                  end if
                !$$$                  write(120+rank,*),i1,list(1:NNBKS,i1)
                !$$$                  call flush(120+rank)
                !     --03/29/14 20:29-lwang-end----------------------------------------*
                !                  IMPI(NNBKS+1,L) = KSLOW(IPAIR)
                !     Do not forget to broadcast HT
                
                !     Finished ksint loop
                END DO

            !     Update Nstepu counter (only contain the total steps inside loop)
                JMPI_THR(1) = NSTEPU - NSTEPUOLD
            !     Update Emerge, BE(3) and RMAX
                EMPI_THR(1) = EMERGE - EMERGEOLD
                EMPI_THR(2) = BE(3) - BEOLD
                EMPI_THR(3) = RMAX

                call cputim(ttks8)
                ttksintp=ttksintp+(ttks8-ttks7)*60.

            !     Communication between nodes
                isend = rank + 1
                if(isend == isize)isend = 0
                irecv = rank - 1
                if(irecv == -1)irecv = isize - 1
            
                #ifdef PUREMPI
                do ir = 0,isize-2
                
                    irank = rank - ir
                    if(irank < 0)irank=irank+isize
                
                    istsen=ista(irank+1)
                    icnt = inum(irank+1)
                
                    if(irank == 0)irank=isize
                    istrec = ista(irank)
                    icnt2 = inum(irank)

                    call cputim(tta)
                !     --12/29/13 22:43-lwang-debug--------------------------------------*
                !**** Note:------------------------------------------------------------**
                !$$$  print*,'SERV',rank,istsen,istrec,icnt,icnt2
                !$$$  call flush(6)
                !     --12/29/13 22:43-lwang-end----------------------------------------*
                    CALL MPI_SENDRECV(ZMPI(1,istsen),57*icnt,MPI_REAL8, &
                    isend,rank,ZMPI(1,istrec),57*icnt2,MPI_REAL8, &
                    irecv,irecv,MPI_COMM_WORLD,status,ierr)
                !$$$                  CALL MPI_SENDRECV(IMPI(1,istsen),lmax*icnt,
                !$$$     *                 MPI_INTEGER,isend,rank,IMPI(1,istrec),
                !$$$     *                 lmax*icnt2,MPI_INTEGER,irecv,irecv,
                !$$$     *                 MPI_COMM_WORLD,status,ierr)
                    call cputim(ttb)
                    call mpi_barrier(MPI_COMM_WORLD,ierr)
                    call cputim(tt999)
                    ibarcount=ibarcount+1
                    ttbar = ttbar + (tt999-ttb)*60
                    ttksbar = ttksbar +(tt999-ttb)*60
                    xtsub3 = xtsub3 + dble((55*8+4)*(icnt+icnt2))
                    ttkssub = ttkssub + (ttb-tta)*60.
                !     --12/29/13 22:43-lwang-debug--------------------------------------*
                !**** Note:------------------------------------------------------------**
                !$$$  print*,'SERV DA',rank,ZMPI(1,1),ZMPI(1,ksb)
                !$$$  call flush(6)
                !     --12/29/13 22:43-lwang-end----------------------------------------*
                end do
                call cputim(tta)
                call MPI_ALLGATHER(NNBKTOT,1,MPI_INTEGER,NNBK_MPI(1), &
                1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
                IMPI_OFF(1) = 0
                DO K =1,isize-1
                    IMPI_OFF(K+1) = IMPI_OFF(K) + NNBK_MPI(K)
                END DO
                CALL MPI_ALLGATHERV(IMPI_THR(1),NNBKTOT,MPI_INTEGER, &
                IMPI(1),NNBK_MPI,IMPI_OFF,MPI_INTEGER, &
                MPI_COMM_WORLD,ierr)
                call cputim(ttb)
                ttkssub = ttkssub + (ttb-tta)*60.
                #endif
                call cputim(ttks9)
            !     --12/30/13 10:11-lwang-debug--------------------------------------*
            !**** Note:------------------------------------------------------------**
            !$$$               if(tBLOCK.eq.6.05468750000000000E-002) then
            !$$$                  write(140+rank,*)'NNBKTOT',NNBKTOT,
            !$$$     &             'NNBK_MPI',NNBK_MPI(1:isize)
            !$$$                  write(140+rank,*)'IMPI_THR',IMPI_THR(1:NNBKTOT)
            !$$$                  write(140+rank,*)'IMPI',
            !$$$     &                 IMPI(1:IMPI_OFF(isize)+NNBK_MPI(isize))
            !$$$                  call flush(140+rank)
            !$$$               end if
            !     --12/30/13 10:11-lwang-end----------------------------------------*
                NNBKTOT = 1
                DO L = 1,KSB
                    I1 = BSLIST(L)
                    IPAIR = KVEC(I1)
                    U(1:4,IPAIR) = ZMPI(1:4,L)
                    U0(1:4,IPAIR) = ZMPI(5:8,L)
                    UDOT(1:4,IPAIR) = ZMPI(9:12,L)
                    FU(1:4,IPAIR) = ZMPI(13:16,L)
                    FUDOT(1:4,IPAIR) = ZMPI(17:20,L)
                    FUDOT2(1:4,IPAIR) = ZMPI(21:24,L)
                    FUDOT3(1:4,IPAIR) = ZMPI(25:28,L)
                    FP0(1:4,IPAIR) = ZMPI(29:32,L)
                    FD0(1:4,IPAIR) = ZMPI(33:36,L)
                    H(IPAIR) = ZMPI(37,L)
                    HDOT(IPAIR) = ZMPI(38,L)
                    HDOT2(IPAIR) = ZMPI(39,L)
                    HDOT3(IPAIR) = ZMPI(40,L)
                    HDOT4(IPAIR) = ZMPI(41,L)
                    DTAU(IPAIR) = ZMPI(42,L)
                    TDOT2(IPAIR) = ZMPI(43,L)
                    TDOT3(IPAIR) = ZMPI(44,L)
                    R(IPAIR) = ZMPI(45,L)
                    R0(IPAIR) = ZMPI(46,L)
                    GAMMA(IPAIR) = ZMPI(47,L)
                    H0(IPAIR) = ZMPI(48,L)
                    SF(1:7,IPAIR) = ZMPI(49:55,L)
                    STEP(I1) = ZMPI(56,L)
                    T0(I1) = ZMPI(57,L)
                    NNBKS = IMPI(NNBKTOT) + 1
                    LIST(1:NNBKS,I1) = IMPI(NNBKTOT:NNBKTOT+NNBKS-1)
                    KSLOW(IPAIR) = IMPI(NNBKTOT+NNBKS)
                    NNBKTOT = NNBKTOT + NNBKS + 1
                !     --12/30/13 9:16-lwang-debug--------------------------------------*
                !**** Note:------------------------------------------------------------**
                !$$$                  if(tBLOCK.eq.6.05468750000000000E-002) then
                !$$$                     write(130+rank,*),'I1',I1,NNBKS,'NNBKTOT',NNBKTOT,
                !$$$     &                    LIST(1:NNBKS,I1),'T',TBLOCK
                !$$$                     call flush(130+rank)
                !$$$                  end if
                !$$$                 if(rank.eq.0) write(110,*),i1,list(1:NNBKS,i1)
                !$$$                  call flush(110)
                !$$$                  write(110+rank,*)'L',L,'I1',I1,'IPAIR',ipair,
                !$$$     &                 'U',U(1,IPAIR),'U0',U0(1,IPAIR),
                !$$$     &                 'UDOT',UDOT(1,ipair),'FU',FU(1,ipair),
                !$$$     &                 'FT',FUDOT(1,ipair),'FT2',FUDOT2(1,ipair),
                !$$$     &                 'FT3',FUDOT3(1,ipair),'FP0',FP0(1,Ipair),
                !$$$     &                 'FD0',FD0(1,ipair),'H',H(ipair),'HT',HDOT(ipair),
                !$$$     &                 'HT2',HDOT2(ipair),'HT3',HDOT3(ipair),
                !$$$     &                 'HT4',HDOT4(ipair),'dt',DTAU(ipair),
                !$$$     &                 'TD2',TDOT2(ipair),'TD3',TDOT3(ipair),
                !$$$     &                 'R',R(ipair),'R0',r0(ipair),
                !$$$     &                 'Gamma',GAMMA(ipair),'H0',h0(ipair),
                !$$$     &                 'SF',sf(1,ipair),'kslow',kslow(ipair),
                !$$$     &                 'nnb',LIST(1,i1),'list',list(2,i1),list(nnbks,i1)
                !$$$     &                 ,'step',step(i1),'t0',t0(i1),'tblock',tblock
                !$$$                  if(i1.eq.3327.and.tblock.ge.4.21E-3) then
                !$$$                     print*,rank,'PARALLEL KSINT TBLOCK',tblock
                !$$$                     call flush(6)
                !$$$                  end if
                !$$$                  call flush(110+rank)
                !$$$                  call mpi_barrier(MPI_COMM_WORLD,ierr)
                !$$$                  if(tblock.gt.6.30187988281250000E-002) stop
                !     --12/30/13 9:16-lwang-end----------------------------------------*
                END DO
                call cputim(ttks10)
                ttkscp=ttkscp+(ttks10-ttks9)*60.

                #ifdef PUREMPI
            !     Gather nstepu, specific iphase and ipair
                CALL MPI_ALLGATHER(JMPI_THR(1),3,MPI_INTEGER, &
                JMPI(1,1),3,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            !     Gather Emerge, BE(3) and RMAX
                CALL MPI_ALLGATHER(EMPI_THR(1),3,MPI_REAL8, &
                EMPI(1,1),3,MPI_REAL8,MPI_COMM_WORLD,ierr)
            !     --03/05/14 20:16-lwang-debug--------------------------------------*
            !$$$***** Note:------------------------------------------------------------**
            !$$$               print*,'r',rank,'THR',JMPI_THR(1:3),
            !$$$     &              'JMPI',JMPI(1:3,1:isize),'ib',ibarcount
            !$$$               call flush(6)
            !$$$*     --03/05/14 20:16-lwang-end----------------------------------------*
            !     Communicate specific value
                CALL ksparmpi(K_comm,0,0,0,0,0.0)
            !     --03/05/14 20:45-lwang-debug--------------------------------------*
            !**** Note:------------------------------------------------------------**
            !$$$               print*,rank,'rmax',rmax,'emerge',emerge,'BE3',BE(3),
            !$$$     *              'EGRAV',egrav,
            !$$$     &              't',tblock
            !$$$               call flush(6)
            !$$$               if(emerge.ge.100)  stop
            !     --03/05/14 20:45-lwang-end----------------------------------------*
                #endif
                NSTEPU = NSTEPUOLD
                EMERGE = EMERGEOLD
                BE(3) = BEOLD
                DO L = 1,isize
                !     Update step counter.
                    NSTEPU = NSTEPU + JMPI(1,L)
                !     Update Emerge.
                    EMERGE = EMERGE + EMPI(1,L)
                !     Update BE(3).
                    BE(3) = BE(3) + EMPI(2,L)
                !     Update RMAX
                    RMAX = MAX(RMAX,EMPI(3,L))
                
                    IF(JMPI(3,L) /= 0) THEN
                    !     Get ipair for specific treatment
                        I10 = JMPI(2,L)
                        KSPAIR = KVEC(I10)
                    !     Note IPHASE must be defined for chain (not needed for collision).
                        IPHASE = JMPI(3,L)
                        IQ = IPHASE
                    !     Treat collision explicitly before quitting (suppressed in KSINT).
                        IF (IQ < 0) THEN
                            IQCOLL = -2
                        !     Ensure block-step time for normal continuation.
                            TIME = TBLOCK
                            CALL CMBODY(2)
                            TIME = TBLOCK
                            RETURN
                        END IF
                    END IF
                END DO
                call cputim(ttks11)
                ttkscmb=ttkscmb+(ttks11-ttks10)*60.
            !     Reset full prediction flag at the end for safe
            !               ipredall=.false.
            !     End of parallel block
            END IF
            #endif
        
        
        !     Prepare pointer and time for INSERT routine
            LIOLD = LI
            LI = LI + KSB
        
        !     --12/21/13 22:31-lwang-print--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$  #ifdef PRINT
        !$$$  PRINT*,'--------------------'
        !$$$  PRINT*,' Before Insert LIOLD,KSB,LI,TBLIST ',
        !$$$  &              LIOLD,KSB,LI,TBLIST
        !$$$  CALL FLUSH(6)
        !$$$  J1 = KBLIST(LI)
        !$$$  JNNTB = KBLIST(NNTB)
        !$$$  PRINT*,' KBL ',J1,JNNTB,T0(J1)+STEP(J1),
        !$$$  &              T0(JNNTB)+STEP(JNNTB)
        !$$$  #endif
        !     --12/21/13 22:31-lwang-end----------------------------------------*
        !     See whether current pair is due before new KBLIST loop.
        !     Insert body #I1 in the correct sequential location.
            call cputim(ttks12)
            DO L = 1,KSB
            
                I1 = BSLIST(L)
            
                TIME = T0(I1) + STEP(I1)
                IF (TIME < TBLIST) THEN
                !$$$  TXLIST(I1) = T0(I1) + STEP(I1)
                                      
                !     --12/21/13 22:31-lwang-print--------------------------------------*
                !**** Note:------------------------------------------------------------**
                !$$$  #ifdef PRINT
                !$$$  PRINT*,' Call Insert I1,TXLIST=',I1,TXLIST(I1)
                !$$$  #endif
                !     --12/21/13 22:31-lwang-end----------------------------------------*
                    CALL INSERT(I1,LI)
                END IF
            
            END DO
        !     --12/30/13 22:55-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            DO L=LI+1,NNTB
        !$$$               K=KBLIST(L)
        !$$$               TTK1 = T0(K)+STEP(K)
        !$$$               TTK2 = T0(K+1)+STEP(K+1)
        !$$$               IF(TTK1.GT.TTK2) then
        !$$$                  PRINT*,'ERROR:TTK1>TTK2',ttk1,ttk2,time,LI,NNTB,K
        !$$$                  call flush(6)
        !$$$                  stop
        !$$$               end if
        !$$$            END DO
        !     --12/30/13 22:55-lwang-end----------------------------------------*
            call cputim(ttks13)
            ttksins=ttksins+(ttks13-ttks12)*60.
        
            9 SPREV = SBLOCK
        
        !     --03/14/14 12:27-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            do L= 1,ksb
        !$$$               i1 = bslist(l)
        !$$$               ipair = kvec(i1)
        !$$$               write(110+rank,*) 'aL',L,'i1',i1,'ip',ipair,'n',name(i1),
        !$$$     &              'u',u(1,ipair),'ud',udot(1,ipair),'t0',t0(i1),
        !$$$     &              'step',step(i1),'tb',tblock,'t',time
        !$$$               call flush(110+rank)
        !$$$            end do
        !$$$            call mpi_barrier(MPI_COMM_WORLD,ierr)
        !     --03/14/14 12:27-lwang-end----------------------------------------*
        END DO
    !     Copy original block time at end of KS treatment.
        TIME = TBLOCK
        NBPREV = NPAIRS
    END IF


!     Check time for advancing any triple, quad or chain regularization.
    20 call cputim(ttks14)
    IF (NSUB > 0) THEN
        30 TSUB = 1.0D+10
        DO L = 1,NSUB
            IF (TS(L) < TSUB) THEN
                ISUB = L
                TSUB = TS(L)
            END IF
        END DO
    
        IF (TSUB <= TBLOCK) THEN
            TIME = TSUB
        !       Decide between triple, quad or chain.
            IF (ISYS(ISUB) == 1) THEN
            !       Update unperturbed size of subsystem and copy c.m. step.
                CALL EXTEND(ISUB)
                CALL TRIPLE(ISUB)
            ELSE IF (ISYS(ISUB) == 2) THEN
                CALL EXTEND(ISUB)
                CALL QUAD(ISUB)
            ELSE
                IF (STEPS(ISUB) < 0.0D0) THEN
                    STEPS(ISUB) = 1.0D-10
                    GO TO 50
                END IF
                CALL CHAIN(ISUB)
                IF (ISUB > 0) THEN
                    IF(STEPS(ISUB) < 0.0D0) THEN
                        STEPS(ISUB) = 1.0D-10
                        GO TO 50
                    END IF
                END IF
            END IF
        
        !       Check for termination (set TPREV < TIME and set IQ < 0).
            IF (ISUB < 0 .OR. IPHASE < 0) THEN
            !     --03/07/14 23:10-lwang-note---------------------------------------*
            !**** Note: In parallel KS, when one particle is integrated in one processor,
            !****       U0 and U0dot is updated. In intgrt, if the time is same as
            !****       the previous time in intgrt, then if time not equal tprev,
            !****       subint will be called. If the particle predict to earlier time
            !****       in ksint in one processor, when go outside subint, the same particle
            !****       in other processors will not be predicted again, but this one will be
            !****       predicted with new U0, which will cause different X and Xdot for
            !****       its members, thus will bring MPI error. Here try to suppress this to avoid tprev lt. time
            !               TPREV = TIME - STEP(NTOT)
                TPREV = TIME
            !     --03/07/14 23:10-lwang-end----------------------------------------*
                IQ = -1
            END IF
            GO TO 30
        END IF
        50 TIME = TBLOCK
    END IF
    call cputim(ttks15)
    tttq = tttq+(ttks15-ttks14)*60.

    RETURN

    END SUBROUTINE SUBINT
