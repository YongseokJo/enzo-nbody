      SUBROUTINE INTGRT
*     label-intgrt
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      include 'tlist.h'
      include 'omp_lib.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/BLKLVL/  JHIST(0:NMAX),JHISTR(0:NMAX),JHISTU(0:NMAX)
      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
      COMMON/NEWXV/ XN(3,NMAX),XNDOT(3,NMAX)
      INTEGER JHIST,JHISTR,LSHORT(NMAX)
      INTEGER IREG(NMAX),IMINR_TMP(NMAX)
      REAL*8 FIRR_SIMD(3,NMAX),FD_SIMD(3,NMAX)

      LOGICAL LSTEPM,FNTLOOP,IPSTART,ipredall
      DATA LSHORT(1) /0/
      DATA IQ,ICALL,LSTEPM,STEPM,IPSTART /0,2,.FALSE.,0.03125,.true./
      SAVE IQ,ICALL,LSTEPM,STEPM,IPSTART,TMIN

      REAL*8 DTOUT
      SAVE DTOUT

*     included by YS Jo

      IF (TIME.EQ.0.0) THEN
          IPSTART = .true.
      END IF


*
*       Update quantized value of STEPM for large N (first time only).
      call cputim(tttbegin)
      IF (.NOT.LSTEPM.AND.NZERO.GT.1024) THEN
         K = INT((FLOAT(NZERO)/1024.0)**0.333333)
         STEPM = 0.03125D0/2**(K-1)
         LSTEPM = .TRUE.
      END IF
*
*       Open GPU for all single particles
      IF (IPSTART) THEN
         FNTLOOP = .false.

*     Initial output interval
         IF(KZ(47).GE.0) THEN
            DTOUT = DTK(KZ(47)+1)
         ELSE
*            DTOUT = DTK(64)
            DTOUT = DTK(1)
            DO I=KZ(47),-1
               DTOUT = DTOUT*2.0D0
            END DO
         END IF

*     Initial nxtlst
         IF(TIME+TOFF.EQ.0.D0) THEN
            NXTLIMIT = 0
            NGHOSTS = 0
            DO I = IFIRST,NTOT
               IF(BODY(I).NE.0.D0) THEN
                  NXTLIMIT = NXTLIMIT + 1
                  NXTLST(NXTLIMIT) = I
               ELSE
                  NXTLST(NMAX-NGHOSTS) = I
                  NGHOSTS = NGHOSTS + 1
               END IF
            END DO
            IF(NGHOSTS.GT.0) NXTLST(NXTLIMIT+1:NXTLIMIT+NGHOSTS) 
     &           = NXTLST(NMAX-NGHOSTS+1:NMAX)

            IF(KZ(46).GT.0) THEN
*     Initial NMINR
               IMINR(1:NTOT) = 1
               DO I = IFIRST,2*NBIN0
                  IF(NAME(I).LE.2*NBIN0) THEN
                     IF(MOD(I,2).EQ.0.AND.NAME(I-1).EQ.NAME(I)-1) THEN
                        IMINR(I) = I - 1
                     ELSE IF(MOD(I,2).EQ.1
     &                       .AND.NAME(I+1).EQ.NAME(I)+1) THEN
                        IMINR(I) = I + 1
                     END IF
                  END IF
               END DO

*     HDF5 format output initialization                  
*     At the first time, full dataset need to be stored,
               if(rank.eq.0) then
                 IF(KZ(46).GE.3.AND.KZ(46).LE.4.AND.SMAX.GT.
     &                 DTK(KZ(47)+1)) then
                     write(6,*) 'Warning!: Maximum time step ',SMAX,
     &                    ' > Output interval ',DTK(KZ(47)+1),
     &                    ' Wide bianries cannot be detected for ',
     &                    ' Non-active particles!'
                  END IF
               end if
               if(rank.eq.0) then
                  IF(KZ(46).EQ.1.OR.KZ(46).EQ.3) THEN
                     call custom_update_file(TTOT,DELTAT)
                  ELSE IF(KZ(46).EQ.2.OR.KZ(46).EQ.4) THEN
                     call custom_update_file(TTOT,DTOUT)
                  END IF
               end if
               CALL custom_output(NXTLEN,NXTLST,
     &              NXTLIMIT,NGHOSTS)
            END IF

*     Sorting NXTLST
            call sort_tlist(STEP,DTK,.true.)
            NXTLEVEL = NDTMAX

*     Initial tcheck block step level
         END IF

         IF(rank.eq.0.and.TOFF+TIME.NE.0.D0) THEN
            IF(KZ(46).EQ.1.OR.KZ(46).EQ.3) THEN
               call custom_update_file(TTOT,DELTAT)
            ELSE IF(KZ(46).EQ.2.OR.KZ(46).EQ.4) THEN
               call custom_update_file(TTOT,DTOUT)
            END IF
         END IF


#ifdef GPU
         NN = NTOT + 10
         CALL GPUNB_OPEN(NN,rank)
*         CALL GPUPOT_INIT_FLOAT(rank)
*         CALL GPUPOT_INIT(rank,NN)
#endif

         IPSTART = .false.
      END IF
      call cputim(tttiend)
      ttintb = ttintb +(tttiend-tttbegin)*60
*
*       Search for high velocities after escape or KS/chain termination.
  999 IF (KZ(37).GT.0.AND.(IPHASE.EQ.-1.OR.IPHASE.GE.2)) THEN
          CALL HIVEL(0)
      END IF
*
     call cputim(tttbegin)

*     Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
*
      IF (IQ.LT.0) ICALL = 0
      IQ = 0
      call cputim(tttiend)
      ttintb = ttintb +(tttiend-tttbegin)*60
*
*       Find all particles due at next block time.
    1 CONTINUE

*
*       Redetermine TMIN after main change to catch new small steps after chain.
      ICALL = ICALL + 1
      IF (ICALL.EQ.2) GO TO 999
*
      call cputim(tttbegin)

*     See whether need to add delay list into NXTLST
      IF(NLSTDELAY(1).GT.0) call delay_add_tlist(T0,STEP,DTK)

*     Determine NXTLEN, TMIN and NXTLEVEL
      call next_tlist(TMIN,DTK,T0)


*     added by sykim, to terminate at EDT(time from enzo)
*      IF (TMIN.GE.TCRIT) THEN

*          write(6,*) "TMIN",TMIN,"TCRIT",TCRIT
*          write (6,*) "TCRIT exceeded TMIN"
*     change variables so that every particle would be evaluated
*          NXTLEN = NXTLIMIT
*          NREG = NXTLIMIT
          
*          DO ICRIT = 1,NXTLEN
                     
*              TIMENW(ICRIT) = TCRIT : include this?
*              JCRIT = NXTLST(ICRIT)
*              IREG(ICRIT) = JCRIT

*              IF (JCRIT.GT.N) THEN
*                  IPAIR = JCRIT - N
*                  IF (LIST(1,2*IPAIR-1).GT.0) NSTEPB = NSTEPB + 1
*              END IF

*          END DO
          
*          LSKIP = NXTLEN
*          TIME = TCRIT

*          GO TO 31
*      END IF

*     end added by sykim
       

*     Determine MPI NXTLST and Regular block
      NREG = 0
      DO L = 1, NXTLEN
         J = NXTLST(L)
      
         TIMENW(J) = T0(J) + STEP(J)

         IF (TIMENW(J).GE.T0R(J) + STEPR(J)) THEN
            NREG = NREG + 1
            IREG(NREG) = J
         END IF
         IF (J.GT.N) THEN
            IPAIR = J - N
            IF (LIST(1,2*IPAIR-1).GT.0) NSTEPB = NSTEPB + 1
         END IF
      END DO
*
      LSKIP = NXTLEN 
      IF (LSKIP.LE.50) THEN
*       Update short timestep list for regularization search.
          CALL SHORT(NXTLEN,NXTLST,LSHORT)
      END IF
*
    
*       Set new time and save block time (for regularization terminations).
      TIME = TMIN

*       added by sykim, continue here when ending
*   31 CONTINUE
*       end added by sykim

      TBLOCK = TIME
      TTOT = TIME + TOFF
*
*     Check option for advancing interstellar clouds.
      IF (KZ(13).GT.0) THEN
          CALL CLINT
      END IF
*
       Check optional integration of cluster guiding centre.
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
          IF (KZ(14).EQ.3.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
             CALL GCINT
          END IF
*       Include mass loss by gas expulsion (Kroupa et al. MN 321, 699).
          IF (MPDOT.GT.0.0D0.AND.TIME + TOFF.GT.TDELAY) THEN
              MP = MP0/(1.0 + MPDOT*(TIME + TOFF - TDELAY))
          END IF
       END IF

      call cputim(tttiend)
      ttintb = ttintb +(tttiend-tttbegin)*60

*
*       Include commensurability test (may be suppressed if no problems).
*     IF (STEP(IMIN).LT.1.0E-15.OR.DMOD(TIME,STEP(IMIN)).NE.0.0D0) THEN
*         WRITE (6,1005) IMIN, NAME(IMIN), NSTEPI, TIME, STEP(IMIN), 
*    &                   TIME/STEP(IMIN)
*1005     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
*    &                        2I5,I11,F12.5,1P,E9.1,0P,F16.4)
*         CALL ABORT
*     END IF
*
*       Check for new regularization at end of block.
      IF (IKS.GT.0) THEN
*       added by sykim
          IF(TIME.GE.TCRIT) THEN
              write(6,*) "need to check IKS"
              IPHASE  = 3
              GO TO 100
          END IF
*       end added by sykim
          TIME = TPREV
          IPHASE = 1
          GO TO 100
      END IF
*
*       Check next adjust time before beginning a new block.
      IF (TIME.GT.TADJ).OR.(TIME.GE.TCRIT) THEN
          TIME = TADJ
          IPHASE = 3
          GO TO 100
      END IF
*
*       Also check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
          TIME = TNEXT
          CALL OUTPUT
          GO TO 1
      END IF
*
*       See whether to advance any close encounters at first new time.
      IF ((TIME.GT.TPREV).AND.(TIME.NE.TCRIT)) THEN
         call cputim(tt5)
         CALL SUBINT(IQ,I10)
         call cputim(tt6)
         ttksintgrt = ttksintgrt + (tt6-tt5)*60.
         ttks = ttks + (tt6-tt5)*60.
*
         IF (IQ.LT.0) GO TO 999
      END IF

*
*       Check regularization criterion for single particles.
      call cputim(tttbegin)
      IKS = 0
      ISMIN = 0
      DSMIN = DTMIN
      IF (LSKIP.LE.50) THEN
*       Search only in prepared list of short-step particles. (R.Sp.)
         ISHORT = LSHORT(1)
         DO L = 2,ISHORT+1
            I = LSHORT(L)
*       Search for minimum timestep candidate for not ordered steps (R.Sp.)
*       Beware that members of LSHORT may be members of KS pair (R.Sp.)
c$$$          If(step(i).LT.1.0D-8) THEN
c$$$             WRITE (6,24)  TIME, I, NAME(I), NXTLEN, NSTEPI,
c$$$     &            STEP(I), STEPR(I), NTOT, NPAIRS
c$$$ 24          FORMAT (' SMALL STEP!!  T I NAME LEN #I SI SR NT NP',F10.5,
c$$$     &            3I6,I11,1P,2E10.2,2I11)
c$$$             CALL FLUSH(6)
c$$$             STOP
c$$$  END IF
            IF (STEP(I).LT.DTMIN.AND.STEP(I).LT.DSMIN.AND.I.LE.N.AND.
     &           I.GE.IFIRST) THEN
               DSMIN = STEP(I)
               ISMIN = I
            END IF
         END DO
      END IF
*
*       See whether dominant body can be regularized.
      IF(ISMIN.GT.0) THEN
         CALL SEARCH(ISMIN,IKS)
*       Include close encounter search for low-eccentric massive binaries.
         IF (IKS.EQ.0.AND.NXTLEVEL.LE.4.AND.STEP(ISMIN).LT.4.0*DTMIN) 
     &        THEN
*       Consider massive single bodies in absence of subsystems.
            IF (ISMIN.LE.N.AND.BODY(I).GT.2.0*BODYM.AND.NSUB.EQ.0) THEN
*
*       Obtain two-body elements and relative perturbation.
               JMIN = 0
               CALL ORBIT(ISMIN,JMIN,SEMI,ECC,GI)
*     
               EB = -0.5*BODY(ISMIN)*BODY(JMIN)/SEMI
               IF (EB.LT.EBH.AND.GI.LT.0.25.AND.JMIN.GE.IFIRST) THEN
                  APO = SEMI*(1.0 + ECC)
*     Check eccentricity (cf. max perturbation) and neighbour radius.
                  IF (ECC.LT.0.25.AND.APO.LT.0.02*RS(ISMIN)) THEN
                     if(rank.eq.0) PRINT*, ' KS TRY: NAM E A EB ',
     *                    NAME(ISMIN), NAME(JMIN), ECC, SEMI, EB
                     CALL FLUSH(6)
                     IKS = IKS + 1
                     ICOMP = ISMIN
                     JCOMP = JMIN
                  END IF
               END IF
            END IF
         END IF
      END IF

      call cputim(tt1)
      ttintb=ttintb +(tt1-tttbegin)*60.
*
*       Initialize counters for irregular & regular integrations.
      TPREV = TIME
      NBLOCK = NBLOCK + 1

*     Initialize counters for irregular & regular integrations.
*      NREG = 0
*
*     write(6,*)' irr ',time,nxtlen,(nxtlst(k),step(k),k=1,5)
*     call flush(6)
*       Advance the irregular step for all particles in the current block.
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
      IF(KZ(33).GT.0)JHIST(NXTLEN) = JHIST(NXTLEN) + 1
*
*     Do full prediction when regular block step come
      ipredall = .false.
      IF(NREG.GT.0) THEN
         call cputim(tttpre1)
*       Predict all particles (except TPRED=TIME) in C++ on host. 
*         CALL CXVPRED(IFIRST,NTOT,TIME,T0,X0,X0DOT,F,FDOT,X,XDOT,TPRED)
         call xbpredall
         call cputim(tttpre2)
         ttpre = ttpre + (tttpre2-tttpre1)*60
         ttpreall = ttpreall + (tttpre2-tttpre1)*60
      END IF

*     Predict chain variables and perturbers at new block-time.
      IF (NCH.GT.0) THEN
          CALL JPRED_int(ICH,TIME)
          CALL XCPRED(2)
      END IF

         call cputim(ttnb1)
*#ifndef SIMD
*     Predict x and xdot
*         call xbpred(NREG,1,NXTLEN,NXTLST)
*#endif
*         call cputim(ttnb2)
*     OpenMP version of irregular force and fdot

         CALL NBINT(TIME,NXTLEN,NXTLST,X,XDOT,BODY,FIRR_SIMD,FD_SIMD,
     &        LIST,IMINR_TMP)

!$omp parallel do if(NXTLEN.GE.ITHREAD) private(L,I)
         DO L = 1,NXTLEN
            I = NXTLST(L)
            CALL NBINT_COR(I,FIRR_SIMD(1,L),FD_SIMD(1,L))
         END DO
!$omp end parallel do

!$omp parallel do if(NXTLEN.GE.ITHREAD) private(L,I,DTR)
         DO L = 1,NXTLEN
            I = NXTLST(L)
            X0(1:3,I) = XN(1:3,I)
c$$$            X(1:3,I) = XN(1:3,I)
            X0DOT(1:3,I) = XNDOT(1:3,I)
c$$$            XDOT(1:3,I) = XNDOT(1:3,I)
            D0(1:3,I) = FI(1:3,I)
            D1(1:3,I) = FIDOT(1:3,I)
            IMINR(I) = IMINR_TMP(L)
*     Save new block step and update T0 & next time
            T0(I) = TIME
            TIMENW(I) = T0(I) + STEP(I)
*     Set non-zero indicator for new regular force.
*     edited by sykim. may be wrong... need to check
            IF ((T0R(I) + STEPR(I).GT.TIME).OR.(TIME.EQ.TCRIT)) THEN
*     end edited by sykim
*     Extrapolate regular force & first derivatives to obtain F & FDOT.
              DTR = TIME - T0R(I)
              F(1,I) = 0.5*(FRDOT(1,I)*DTR + FR(1,I) + FI(1,I))
              F(2,I) = 0.5*(FRDOT(2,I)*DTR + FR(2,I) + FI(2,I))
              F(3,I) = 0.5*(FRDOT(3,I)*DTR + FR(3,I) + FI(3,I))
              FDOT(1,I) = ONE6*(FRDOT(1,I) + FIDOT(1,I))
              FDOT(2,I) = ONE6*(FRDOT(2,I) + FIDOT(2,I))
              FDOT(3,I) = ONE6*(FRDOT(3,I) + FIDOT(3,I))
           END IF
         END DO
!$omp end parallel do

         call cputim(tt3)
         ttirr = ttirr + (tt3-ttnb1)*60.
*


c$$$        call cputim(ttt32)
c$$$**!$omp parallel do if(NXTLEN.GE.ithread)
c$$$**!$omp& default(shared) private(L,I,DTR,K)
c$$$        DO L = 1,NXTLEN
c$$$           I = NXTLST(L)
c$$$*       Save new block step and update T0 & next time
c$$$           T0(I) = TIME
c$$$           TIMENW(I) = T0(I) + STEP(I)
c$$$           
c$$$*     Set non-zero indicator for new regular force.
c$$$           IF (T0R(I) + STEPR(I).GT.TIME) THEN
c$$$*     *!$omp critical 
c$$$*              NREG = NREG + 1
c$$$*              IREG(NREG) = I
c$$$**!$omp end critical
c$$$*          ELSE
c$$$*       Extrapolate regular force & first derivatives to obtain F & FDOT.
c$$$              DTR = TIME - T0R(I)
c$$$              F(1,I) = 0.5*(FRDOT(1,I)*DTR + FR(1,I) + FI(1,I))
c$$$              F(2,I) = 0.5*(FRDOT(2,I)*DTR + FR(2,I) + FI(2,I))
c$$$              F(3,I) = 0.5*(FRDOT(3,I)*DTR + FR(3,I) + FI(3,I))
c$$$              FDOT(1,I) = ONE6*(FRDOT(1,I) + FIDOT(1,I))
c$$$              FDOT(2,I) = ONE6*(FRDOT(2,I) + FIDOT(2,I))
c$$$              FDOT(3,I) = ONE6*(FRDOT(3,I) + FIDOT(3,I))
c$$$* Higher order extrapolation?
c$$$*                 F(K,I) = FI(K,I) + FR(K,I) + DTR*(FRDOT(K,I)
c$$$*    *                + DTR*(D2R(K,I)/2.D0 + DTR*D3R(K,I)/6.D0))
c$$$*                 FDOT(K,I) = FIDOT(K,I) + FRDOT(K,I)
c$$$*    *                + DTR*(D2R(K,I) + DTR*D3R(K,I)/2.D0)
c$$$*                 F(K,I) = F(K,I)/2.D0
c$$$*                 FDOT(K,I) = FDOT(K,I)/6.D0
c$$$           END IF
c$$$          
c$$$*     X0(K,I) = XN(K,I)
c$$$*     X0DOT(K,I) = XNDOT(K,I)
c$$$        END DO
c$$$**!$omp end parallel do
*
*
*       See whether any KS candidates are in the same block.
        IF (IKS.GT.0) THEN
*       Accept same time, otherwise reduce STEP(ICOMP) and/or delay.
           IF (T0(JCOMP).EQ.T0(ICOMP)) THEN
*             call cputim(ttnb1)
*             call xbpredall
*             call cputim(ttnb2)
*             ttpre = ttpre + (ttnb2-ttnb1)*60.
              I = ICOMP
              ICOMP = MIN(ICOMP,JCOMP)
              JCOMP = MAX(I,JCOMP)
              call jpred_int(icomp,time)
              call jpred_int(jcomp,time)
*     --09/25/13 21:44-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$              if(time.ge.0.3) then
c$$$                 print*,rank,'icomp',icomp,name(icomp),
c$$$     &            x(1,icomp),xdot(1,icomp),x0(1,icomp),x0dot(i,icomp)
c$$$              end if
*     --09/25/13 21:44-lwang-end----------------------------------------*
           ELSE IF (T0(JCOMP) + STEP(JCOMP).LT.T0(ICOMP)) THEN
              STEP(ICOMP) = 0.5D0*STEP(ICOMP)
              TIMENW(ICOMP) = T0(ICOMP) + STEP(ICOMP)
              IKS = 0
           ELSE
              IKS = 0
           END IF
        END IF
*
      NSTEPI = NSTEPI + NXTLEN

*
c$$$      call cputim(ttt33)
c$$$      ttnewt = ttnewt + (ttt33 - ttt32)*60.

***** --Regular force calculation--------------------------------------**
      IF(NREG.GT.0)THEN
      
#ifdef GPU
*       Send all single particles to GPU memory
         call cputim(tt51)
         NN = NTOT - IFIRST + 1
         CALL GPUNB_SEND(NN,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
         call cputim(tt52)
*     --09/26/13 16:58-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(rank.eq.0) print*,'NREG',NXTLEN,NREG,TIME,NGPUC,tt52-tt51
c$$$      NGPUC = NGPUC + 1
c$$$      if(rank.eq.0) then
c$$$         do i=1,nreg
c$$$            L=ireg(i)
c$$$            write(100,*) NREG, Name(L), STEPR(L),TIME
c$$$         end do
c$$$      end if
*     --09/26/13 16:58-lwang-end----------------------------------------*
         ttgrcomm = ttgrcomm + (tt52-tt51)*60.0
#endif

*
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
         IF(KZ(33).GT.1)JHISTR(NREG) = JHISTR(NREG) + 1
*

            call cputim(tt1)
            CALL CALC_REG_ON_GPU(IREG,1,NREG)
            call cputim(tt2)
            ttgrcalc = ttgrcalc + (tt2-tt1)*60.0
*
!$omp parallel do private(L,I)
            DO L = 1,NREG
               I = IREG(L)
               F(1,I) = 0.5D0*(FI(1,I) + FR(1,I))
               F(2,I) = 0.5D0*(FI(2,I) + FR(2,I))
               F(3,I) = 0.5D0*(FI(3,I) + FR(3,I))
               FDOT(1,I) = ONE6*(FIDOT(1,I) + FRDOT(1,I))
               FDOT(2,I) = ONE6*(FIDOT(2,I) + FRDOT(2,I))
               FDOT(3,I) = ONE6*(FIDOT(3,I) + FRDOT(3,I))               
               X0(1:3,I) = XN(1:3,I)
c$$$               X(1:3,I) = XN(1:3,I)
               X0DOT(1:3,I) = XNDOT(1:3,I)
c$$$               XDOT(1:3,I) = XNDOT(1:3,I)
*     Copy neighbor list
*               NNB = IMPI(1,L) + 1 
*               LIST(1:NNB,I) = IMPI(1:NNB,L)
               IF(LIST(1,I).GT.0) then
!$omp critical                   
                  RSMIN = MIN(RSMIN,RS(I))
!$omp end critical
               end if
           END DO
!$omp end parallel do
*     
            NBSUM = 0
            call cputim(tt2)
            ttreg = ttreg + (tt2-tt1)*60.


c$$$         call cputim(tt333) 
         NSTEPR = NSTEPR + NREG
         NBLCKR = NBLCKR + 1
c$$$*
c$$$c$$$**!$omp parallel do if(NREG.GE.ithread) private(L,I,K)      
c$$$c$$$         DO L = 1,NREG
c$$$c$$$            I = IREG(L)
c$$$c$$$*     
c$$$c$$$            LIST(1:LMAX,I) = IMPI(1:LMAX,L)
c$$$c$$$*
c$$$c$$$*       Check minimum neighbor sphere since last output
c$$$c$$$**!$omp critical          
c$$$c$$$            IF(LIST(1,I).GT.0)RSMIN = MIN(RSMIN,RS(I))
c$$$c$$$**!$omp end critical              
c$$$c$$$*
c$$$c$$$         END DO
c$$$c$$$**!$omp end parallel do      
c$$$         call cputim(tt334)
c$$$         ttmov = ttmov +(tt334-tt333)*60
*
*          OPEN(98,STATUS='OLD',ERR=123)
*         print*,' last reg block t=',time,' length=',nreg
*         print*,' first 10 =',(name(ireg(l)),l=1,min(nreg,10))
*         call flush(6)
*         CLOSE(98)
*123      CONTINUE
*
      END IF

     custom output
      call cputim(ttiout1)
      IF(KZ(46).GT.0.AND.DMOD(TIME,DTOUT).EQ.0) THEN
         call custom_output(NXTLEN,NXTLST,NXTLIMIT,NGHOSTS,IMINR)
      END IF
      call cputim(ttiout2)
      ttout = ttout + (ttiout2 - ttiout1)*60.

*     added by sykim, to end the calculation at TMIN
*     everything must be returned before new sorting of particles    
      IF (TIME.EQ.TCRIT) THEN
          IPHASE = 3
          GO TO 100
      END IF

*     Resorting the NXTLST
      call cputim(tt335)
      call sort_tlist(STEP,DTK,.false.)
      call cputim(tt336)
      ttintb =ttintb + (tt336-tt335)*60
      ttnewt = ttnewt + (tt336 - tt335)*60.      
     Output of Massive black hole data
      IF (KZ(24).EQ.1) then
         call BHPLOT
      END IF
c$$$*     --10/27/13 10:33-lwang-add----------------------------------------*
c$$$***** if ks need terminate, do prediction for icomp and jcomp-----------*
c$$$      if (iks.gt.0) then
c$$$         call xvpred(ICOMP,0)
c$$$         call xvpred(JCOMP,0)
c$$$      end if
c$$$*     --10/27/13 10:33-lwang-end----------------------------------------*
*
**!$omp parallel do if (nreg.ge.ithread) private(L,I)      
*       DO L =1,NREG
*          I = IREG(L)
*          DO 82 K = 1,3
*             X0(K,I) = XN(K,I)
*             X0DOT(K,I) = XNDOT(K,I)
* 82       CONTINUE
*       END DO
**!$omp end parallel do      
*     Needed to update only active particles, avoid sending all at reg block
*     Take care is this really efficient?
*         CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
*    &                                            BODY(I),T0(I))
      IF (NTAIL.GT.0) THEN
*     Allow large quantized interval with internal iteration.
         IF (DMOD(TIME,0.25D0).EQ.0.0D0) THEN
#ifndef TT
!$omp parallel do private(J)      
#endif
            DO 65 J = ITAIL0,NTTOT
               IF (TIMENW(J).LE.TIME) THEN
                  CALL NTINT(J)
               END IF
 65         CONTINUE
#ifndef TT
!$omp end parallel do
#endif
         END IF
      END IF
*
*       Exit on KS termination, new multiple regularization or merger.
      IF (IQ.GT.0) THEN
          NBPREV = 0
          IF (IQ.GE.4.AND.IQ.NE.7) THEN
              CALL DELAY(IQ,-1)
          ELSE
*       Ensure correct KS index (KSPAIR may denote second termination).
              KSPAIR = KVEC(I10)
              IPHASE = IQ
          END IF
          GO TO 100
      END IF
*
*       Perform optional check on high-velocity particles at major times.
      IF (KZ(37).GT.0.AND.LISTV(1).GT.0) THEN
          call cputim(tttsa)
          IF (DMOD(TIME,STEPM).EQ.0.0D0.AND.NREG.GT.0) THEN
              CALL SHRINK(TMIN)
              IF (LISTV(1).GT.0) THEN
                  CALL HIVEL(-1)
              END IF
          END IF
          call cputim(tttsb)
          ttshk = ttshk + (tttsb-tttsa)*60.0
      END IF
*
*       Check optional mass loss time.
      IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 1000-year step (new polynomials).
          IF (NREG.NE.0.AND.TIME.GT.TMDOT
     *        .AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
c$$$          IF (TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
             call cputim(tttmdota)
             IF (KZ(19).GE.3) THEN
                CALL MDOT
             ELSE
                CALL MLOSS
             END IF
             call cputim(tttmdotb)
             ttmdot = ttmdot +(tttmdotb-tttmdota)*60
             IF (IPHASE.LT.0.OR.IPHASE.EQ.1.OR.IPHASE.EQ.2) GO TO 999
          END IF
      END IF
*
*       Advance counters and check timer & optional COMMON save (NSUB = 0).
      NTIMER = NTIMER + NXTLEN
      IF (NTIMER.LT.NMAX) GO TO 1
      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.100*NMAX.AND.NSUB.EQ.0) THEN
         NSTEPS = 0
         IF (KZ(1).GT.1) CALL MYDUMP(1,1)
      END IF
*     
C*     Check option for general binary search.
C      IF (KZ(4).NE.0.AND.TIME - TLASTS.GT.DELTAS) THEN
C         CALL EVOLVE(0,0)
C      END IF
*     
*     Include facility for termination of run (create dummy file STOP).
      IF(rank.EQ.0)THEN
         OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
         IF (IO.EQ.0) THEN
            CLOSE (99)
            IF (NSUB.EQ.0.and.rank.eq.0) WRITE (6,90)
 90         FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
            CPU = 0.0
         END IF
      END IF
*     
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      TCOMP = (TCOMP-TTOTA)*60.
*     
      IF (TCOMP.LT.CPU) GO TO 1
*     
*     Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*     Specify zero step to enforce termination.
         STEPS(1:NSUB) = 0.0D0
         NTIMER = NMAX
         GO TO 1
      END IF
*     
*     Terminate run with optional COMMON save.
      IF (KZ(1).NE.0) THEN
         CPUTOT = CPUTOT + TCOMP - CPU0
         CALL MYDUMP(1,1)
         if(rank.eq.0)
     &        WRITE (6,98)  TOFF, TIME, TIME+TOFF, TCOMP, CPUTOT/60.0, 
     &        ERRTOT, DETOT
 98      FORMAT (//,9X,'COMMON SAVED AT TOFF/TIME/TTOT =',1P,E16.8,
     &        '  TCOMP =',E16.8,'  CPUTOT =',E16.8,
     &                  '  ERRTOT =',F10.6,'  DETOT =',F10.6)
      END IF
*
*     Determine time interval and step numbers per time unit
      TIMINT = TIME + TOFF - TINIT
*
#ifdef PARALLEL
      IF(rank.EQ.0)THEN
#endif
         WRITE (6,195)  rank,TIMINT,NSTEPI-NIR,NSTEPB-NIB,NSTEPR-NRGL,
     &        NSTEPU-NKS
 195     FORMAT (//,I9,' INTEGRATION INTERVAL =',F8.2,3X,' NIRR=',I11,
     &        ' NIRRB=',I11,' NREG=',I11,' NKS=',I11)
         WRITE (6,196)  (NSTEPI-NIR)/TIMINT,(NSTEPB-NIB)/TIMINT,
     &        (NSTEPR-NRGL)/TIMINT,(NSTEPU-NKS)/TIMINT
 196     FORMAT (//,9X,' PER TIME UNIT: NIRR=',1P,D12.5,' NIRRB=',
     &        D12.5,' NREG=',D12.5,' NKS=',D12.5)
#ifdef PARALLEL
      END IF
#endif
#ifdef GPU
      CALL GPUNB_CLOSE
*      CALL GPUPOT_CLOSE
*      CALL GPUPOT_CLOSE_FLOAT
#endif
#ifdef SIMD
      CALL IRR_SIMD_CLOSE(rank)
#endif
#ifdef PARALLEL
      call cputim(tt998)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call cputim(tt999)
      ibarcount=ibarcount+1
      ttbar = ttbar + (tt999-tt998)*60
      CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif
      STOP
*
 100  CONTINUE
*
*     Set current global time.
      TTOT = TIME + TOFF

      RETURN
*

      END
