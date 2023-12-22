    SUBROUTINE REGCOR_GPU(I,XI,XIDOT,FREG,FDR,NLIST)


!       Regular integration.
!       --------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    INCLUDE 'omp_lib.h'
!     Safe for parallel, used variable: Listc
    COMMON/NEWXV/ XN(3,NMAX),XNDOT(3,NMAX)
!      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
!     &                LISTC(LMAX)
    REAL*8 ::  XI(3),XIDOT(3),DFIRR(3),FREG(3),DV(3),DFD(3),FDR(3)
    REAL*8 ::  FRX(3),FDX(3)
    INTEGER :: NLIST(LMAX),JJLIST(2*LMAX)
!     REAL*8  DFI(3),DFRG(3)
!      LOGICAL LCHECK

!     Mass weighted neighbor list criterion flag
    IF(KZ(39) == 2 .OR. KZ(39) == 3) M_FLAG = 1
    IF(KZ(39) == 0 .OR. KZ(39) == 2) NB_FLAG = 1
          
!       Set neighbour number, time-step & choice of central distance.
    NNB = NLIST(1)
    NNB0 = LIST(1,I)
    DTR = TIME - T0R(I)
    IRSKIP = 0
!      NBSKIP = 0
    IF (NB_FLAG == 1) THEN
        RI2 = (XI(1) - RDENS(1))**2 + (XI(2) - RDENS(2))**2 + &
        (XI(3) - RDENS(3))**2
        RH2 = RSCALE**2
    ELSE
        RI2 = XI(1)**2 + XI(2)**2 + XI(3)**2
        RH2 = 9.0*RSCALE**2
    END IF

!       Obtain irregular & regular force and determine current neighbours.
    RS2 = RS(I)**2
!       Set radial velocity factor for the outer shell.
    VRFAC = -0.1*RS2/DTR

!     Choose appropriate force loop for single particle or c.m. body.
!$$$      IF (I.GT.N) THEN
!$$$*     Treat unperturbed KS in the single particle approximation.
!$$$         I1 = 2*(I - N) - 1
!$$$         IF (LIST(1,I1).GT.0) THEN
!$$$*     Correct irregular force on active c.m. body (same form as NBINT).
!$$$            CALL CMFIRR_COR(I,I1,FIRR,FD)
!$$$            GO TO 20
!$$$         END IF
!$$$      END IF

!     Add corrections from c.m. particles (singles or unperturbed).
!$$$      IF (NPAIRS.GT.0) THEN
!$$$          CALL CMFIRR_UCOR(I,NNB,NLIST,FIRR,FD)
!$$$      END IF
!     --10/17/13 19:07-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      if(rank.eq.0) print*,'CM FIRR',FIRR(1:3)
!$$$      call flush(6)
!$$$      if(time.ge.0.34741210937500000.and.name(i).eq.6193) then
!$$$         print*,rank,'ucm',i,name(i),firr,fd,time,x0dot(1,i),xdot(1,i)
!$$$         call flush(6)
!$$$      end if
!$$$      if(name(i).eq.1886) then
!$$$         write(105+rank,*),'reg c i',i,'n',name(i),'firr',firr(1),
!$$$     &        'fi',fi(1,I),'df',firr(1)-fi(1,i),
!$$$     &        't',time
!$$$         call flush(105+rank)
!$$$      end if
!     --10/27/13 18:21-lwang-end----------------------------------------*

!     Include treatment for regularized clump.
!$$$      IF (NCH.GT.0) THEN
!$$$*     Distinguish between chain c.m. and any other particle.
!$$$         IF (NAME(I).EQ.0) THEN
!$$$            CALL CHFIRR(I,1,XI,XIDOT,FIRR,FD)
!$$$         ELSE
!$$$*     Search the chain perturber list for #I.
!$$$            NP1 = LISTC(1) + 1
!$$$            DO 15 L = 2,NP1
!$$$               J = LISTC(L)
!$$$               IF (J.GT.I) GO TO 20
!$$$               IF (J.EQ.I) THEN
!$$$                  CALL FCHAIN(I,1,XI,XIDOT,FIRR,FD)
!$$$                  GO TO 20
!$$$               END IF
!$$$ 15         CONTINUE
!$$$         END IF
!$$$      END IF

!       Check whether cloud forces should be included.
    IF (KZ(13) > 0) THEN
        CALL FCLOUD(I,FREG,FDR,1)
    END IF

!       See whether an external force should be added.
    IF (KZ(14) > 0) THEN
        #ifdef TT
    !** FlorentR - include the case of tidal tensor
    !       Save current values for deriving work done by tides (#14 = 3 or 9).
        IF (KZ(14) == 3 .OR. KZ(14) == 9) THEN
        !** FRenaud
            #else
        !       Save current values for deriving work done by tides (#14 = 3).
            IF (KZ(14) == 3) THEN
                #endif
                DO 22 K = 1,3
                    FRX(K) = FREG(K)
                    FDX(K) = FDR(K)
                22 END DO
            END IF
        
        !       Obtain the tidal perturbation (force and first derivative).
            CALL XTRNLF(XI,XIDOT,DFIRR,FREG,DFD,FDR,1)
        
        END IF

    !       added on 23.04.24 by sykim
    !       add external force recieved from ENZO

        DO 23 K=1,3
            FREG(K) = FREG(K) + FENZO(K,I)
        23 END DO

    !       end added by sykim


    
    !       Check for zero neighbour number.
        IF (NNB == 0) THEN
        !       Assume small mass at centre for distant body or no neighbours.
            R2 = XI(1)**2 + XI(2)**2 + XI(3)**2
            FIJ = 0.01*BODYM/(R2*SQRT(R2))
            RDOT = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) + &
            XI(3)*XIDOT(3))/R2
            DO 25 K = 1,3
                DFIRR(K) = - FIJ*XI(K)
                DFD(K) = - (XIDOT(K) - RDOT*XI(K))*FIJ
            25 END DO
        !       Modify neighbour sphere gradually but allow encounter detection.
            IF (RDOT > 0.0) THEN
                RS(I) = MAX(0.75*RS(I),0.1*RSCALE)
            ELSE
                RS(I) = 1.1*RS(I)
            END IF
        ! omp critical
            NBVOID = NBVOID + 1
        ! omp end critical
            IRSKIP = 1
        !       Skip another full N loop (necessary for GPU/OpenMP case).
            GO TO 50
        !         GO TO 1
        END IF
    
    !       Stabilize NNB between ZNBMIN & ZNBMAX by square root of contrast.
    !       Include optional stabilization to increase neighbour number.
    !       Take input parameter NNBOPT as optimal neighbour number (R.Sp.)
    !       Note that it substitutes input parameter NNBMAX, which
    !       is now a parameter NNBMAX=LMAX-3
    
    !      FAC = 1.D0
    !      IF (KZ(40).GT.0) THEN
        FAC = 1.0 + 0.1*(FLOAT(NNBOPT) - FLOAT(NNB))/FLOAT(NNB)
    !      END IF
    
    !     A3 = ALPHA*FAC*SQRT(FLOAT(NNB)*RS(I))/RS2
        A3 = FLOAT(NNB)*FAC
        NBAVE = MIN(NNBOPT+150,MAX(NNBMAX-50,INT(ZNBMIN)))
        A3 = MIN(A3,0.9*NBAVE)
        NBP = INT(A3)
    !       Reduce predicted membership slowly outside half-mass radius.
        IF (RI2 > RH2 .AND. NB_FLAG == 1) THEN
            A3 = A3*MAX(RSCALE/SQRT(RI2),0.1)
        ELSE
            A3 = MAX(A3,ZNBMIN)
        END IF
        A4 = A3/FLOAT(NNB)
    
    !       Include inertial factor to prevent resonance oscillations in RS.
        IF ((A3 - FLOAT(NNB0))*(A3 - FLOAT(NNB)) < 0.0) A4 = SQRT(A4)
    
    !       Modify volume ratio by radial velocity factor outside the core.
        RIDOT = (XI(1) - RDENS(1))*XIDOT(1) + &
        (XI(2) - RDENS(2))*XIDOT(2) + &
        (XI(3) - RDENS(3))*XIDOT(3)
        IF (RI2 > RC2 .AND. NB_FLAG == 1 .AND. RI2 < 100.0*RH2) THEN
            A4 = A4*(1.0 + RIDOT*DTR/RI2)
        END IF
    
    !       See whether neighbour radius of c.m. body should be increased.
        IF (I > N .AND. RI2 <= RH2) THEN
        !       Set perturber range (soft binaries & H > 0 have short duration).
            A2 = 100.0*BODY(I)/ABS(H(I-N))
        !       Stabilize NNB on ZNBMAX if too few perturbers.
        !     --10/31/13 21:30-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$          print*,rank,'I',I,'NAME',NAME(I),'A2',A2,'RS',RS(I),'NNB',NNB
        !$$$          call flush(6)
        !     --10/31/13 21:30-lwang-end----------------------------------------*
            IF (A2 > RS(I)) THEN
                A3 = MAX(1.0 - FLOAT(NNB)/NBAVE,0.0D0)
            !       Modify volume ratio by approximate square root factor.
                A4 = 1.0 + 0.5*A3
            END IF
        END IF
    
    !       Limit change of RS for small steps (RSFAC = MAX(25/TCR,1.5*VC/RSMIN).
        A5 = MIN(RSFAC*DTR,0.25D0)
    !       Restrict volume ratio by inertial factor of 25 per cent either way.
        A4 = A4 - 1.0
        IF (A4 > A5) THEN
            A4 = A5
        ELSE
            A4 = MAX(A4,-A5)
        END IF
    
    !       Modify neighbour sphere radius by volume factor.
        IF (IRSKIP == 0) THEN
            A3 = ONE3*A4
        !       Use second-order cube root expansion (maximum volume error < 0.3 %).
            A1 = 1.0 + A3 - A3*A3
        !       Prevent reduction of small NNB if predicted value exceeds ZNBMIN.
            IF (NNB < ZNBMIN .AND. NBP > ZNBMIN) THEN
                IF (A1 < 1.0) A1 = 1.05
            END IF
        !       Skip modification for small changes (avoids oscillations in RS).
            IF (ABS(A1 - 1.0D0) > 0.003) THEN
            !     Restrict change in RS within 25 % of RI using approx (RS-RI)*(RS+RI).
                IF (RS(I)**2 > RI2 .AND. RS(I)**2 < 2.0*RI2 .AND. &
                NNB > 5) THEN
                    A1 = SQRT(A1)
                END IF
                RS(I) = A1*RS(I)
            END IF
        !     Employ extra reduction for inward motion and large membership.
            IF (NNB > NBAVE) THEN
                IF (RIDOT < 0.0) THEN
                    DRS = (FLOAT(NNB) - NBAVE)/(3.0*NNB)
                    RS(I) = RS(I)*(1.0 - DRS)
                END IF
            END IF
        END IF
    
        IF(M_FLAG == 0) THEN
        !       Calculate the radial velocity with respect to at most 3 neighbours.
            IF (NNB <= 3 .AND. RI2 < 100.0*RH2) THEN
                A1 = 2.0*RS(I)
                DO 45 L = 1,NNB
                    J = NLIST(L+1)
                    RIJ = SQRT((XI(1) - X(1,J))**2 + (XI(2) - X(2,J))**2 + &
                    (XI(3) - X(3,J))**2)
                    RSDOT = ((XI(1) - X(1,J))*(XIDOT(1) - XDOT(1,J)) + &
                    (XI(2) - X(2,J))*(XIDOT(2) - XDOT(2,J)) + &
                    (XI(3) - X(3,J))*(XIDOT(3) - XDOT(3,J)))/RIJ
                !     Find smallest neighbour distance assuming constant regular step.
                    A1 = MIN(A1,RIJ + RSDOT*DTR)
                45 END DO
            
            !     Increase neighbour sphere if all members are leaving inner region.
                RS(I) = MAX(A1,1.1*RS(I))
            END IF
        END IF
    
    !       Check optional procedures for adding neighbours.
        IF ((KZ(37) == 1 .AND. LISTV(1) > 0) .OR. KZ(37) > 1) THEN
        !         NBSKIP = NBSKIP - NNB
            CALL CHECKL(I,NNB,XI,XIDOT,RS2,DFIRR,FREG,DFD,FDR,NLIST)
        !         NBSKIP = NBSKIP + NNB
        END IF
    
    !     Initial irregular force change, Not use DFIRR and DFD from xtrnlf and checkl
        DFIRR(1:3) = 0.D0
        DFD(1:3) = 0.D0
                  
    !       Find loss or gain of neighbours at the same time.
        50 NBLOSS = 0
        NBGAIN = 0
    
    !       Check case of zero old or new membership (skip if both are zero).
        IF (NNB0 == 0) THEN
            NBGAIN = NNB
            DO 52 L=1,NNB
                JJLIST(L)=NLIST(L+1)
            52 END DO

            DF2 = 1.0
        !$$$          FI2 = 0.0
            FR2 = 0.0

            GO TO 70
        END IF
    
    !     Form expressions for deciding on corrections and time-step.
        DF2 = 0.0
        FR2 = 0.0
    !$$$      FI2 = 0.0
    !$$$      FDR2 = 0.0
        DO 55 K = 1,3
            DF2 = DF2 + (FREG(K) - FR(K,I))**2
            FR2 = FR2 + FREG(K)**2
        !$$$          FI2 = FI2 + FIRR(K)**2
        !$$$          FDR2 = FDR2 + FDR(K)**2
        55 END DO
    
    !       Skip neighbour list comparison without derivative corrections.
    !      IF (KZ(38).EQ.0.OR.(KZ(38).EQ.2.AND.DF2.LT.0.0001*FR2)) GO TO 70

        JMIN = 0
        L = 2
        LG = 2
    !       Set termination value in NLIST(NNB+2) and save last list member.
        NLIST(NNB+2) = NTOT + 1
        NLIST(1) = LIST(NNB0+1,I)
    
    !       Compare old and new list members in locations L & LG.
        56 IF (LIST(L,I) == NLIST(LG)) GO TO 58
    
    !       Now check whether inequality means gain or loss.
        IF (LIST(L,I) >= NLIST(LG)) THEN
            NBGAIN = NBGAIN + 1
            JJLIST(NNB0+NBGAIN) = NLIST(LG)
        !       Number of neighbour losses can at most be NNB0.
            L = L - 1
        !       The same location will be searched again after increasing L below.
        ELSE
            NBLOSS = NBLOSS + 1
            J = LIST(L,I)
            JJLIST(NBLOSS) = J
        !       Check SMIN step indicator (rare case permits fast skip below).
            IF (STEP(J) < SMIN) JMIN = J
            LG = LG - 1
        END IF
    
    !       See whether the last location has been checked.
        58 IF (L <= NNB0) THEN
            L = L + 1
            LG = LG + 1
        !       Last value of second search index is NNB + 2 which holds NTOT + 1.
            GO TO 56
        ELSE IF (LG <= NNB) THEN
            LG = LG + 1
            LIST(L,I) = NTOT + 1
        !       Last location of list holds termination value (saved in NLIST(1)).
            GO TO 56
        END IF
    
    !       See whether any old neighbour with small step should be retained.
        IF (JMIN == 0) GO TO 70
    
        K = 1
        60 IF (NNB > NNBMAX .OR. I > N) GO TO 70
        J = JJLIST(K)
    !       A single regularized component will be replaced by the c.m.
        IF (STEP(J) > SMIN .OR. J < IFIRST .OR. J > N) GO TO 68
    !       Retain old neighbour inside 2*RS to avoid large correction terms.
        RIJ2 = (XI(1) - X(1,J))**2 + (XI(2) - X(2,J))**2 + &
        (XI(3) - X(3,J))**2
        IF (RIJ2 > 4.0*RS2) GO TO 68
    
        L = NNB + 1
        62 IF (NLIST(L) < J) GO TO 64
        NLIST(L+1) = NLIST(L)
        L = L - 1
        IF (L > 1) GO TO 62
    
    !       Save index of body #J and update NNB & NBLOSS.
        64 NLIST(L+1) = J
    !      NBSKIP = 1
        NNB = NNB + 1
        NBLOSS = NBLOSS - 1
    !       Restore last old neighbour in case NBLOSS = 0 at end of search.
        LIST(NNB0+1,I) = NLIST(1)
    ! omp critical
        NBSMIN = NBSMIN + 1
    ! omp end critical
    
    !       Perform correction to irregular and regular force components.
        A1 = X(1,J) - XI(1)
        A2 = X(2,J) - XI(2)
        A3 = X(3,J) - XI(3)
        DV(1) = XDOT(1,J) - XIDOT(1)
        DV(2) = XDOT(2,J) - XIDOT(2)
        DV(3) = XDOT(3,J) - XIDOT(3)
    
        RIJ2 = A1*A1 + A2*A2 + A3*A3
        DR2I = 1.0/RIJ2
        DR3I = BODY(J)*DR2I*SQRT(DR2I)
        DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
        DRDP = 3.0*DRDV*DR2I
    
    !$$$      FIRR(1) = FIRR(1) + A1*DR3I
    !$$$      FIRR(2) = FIRR(2) + A2*DR3I
    !$$$      FIRR(3) = FIRR(3) + A3*DR3I
    !$$$      FD(1) = FD(1) + (DV(1) - A1*DRDP)*DR3I
    !$$$      FD(2) = FD(2) + (DV(2) - A2*DRDP)*DR3I
    !$$$      FD(3) = FD(3) + (DV(3) - A3*DRDP)*DR3I
        FREG(1) = FREG(1) - A1*DR3I
        FREG(2) = FREG(2) - A2*DR3I
        FREG(3) = FREG(3) - A3*DR3I
        FDR(1) = FDR(1) - (DV(1) - A1*DRDP)*DR3I
        FDR(2) = FDR(2) - (DV(2) - A2*DRDP)*DR3I
        FDR(3) = FDR(3) - (DV(3) - A3*DRDP)*DR3I
    
    !       Remove body #J from JJLIST unless it is the last or only member.
        IF (K > NBLOSS) GO TO 70
        DO 66 L = K,NBLOSS
            JJLIST(L) = JJLIST(L+1)
        66 END DO
    !       Index of last body to be moved up is L = NBLOSS + 1.
        K = K - 1
    !       Check the same location again since a new body has appeared.
        68 K = K + 1
    !       Last member to be considered is in location K = NBLOSS.
        IF (K <= NBLOSS) GO TO 60
    
    !       Form time-step factors and update regular force time.
        70 CONTINUE
    
    !     Irregular force correction due to neighbor changes
    !     Remove lost neighbor particle to irregular force differences
        DO L = 1,NBLOSS
            J = JJLIST(L)

            A1 = X(1,J) - XI(1)
            A2 = X(2,J) - XI(2)
            A3 = X(3,J) - XI(3)
            DV(1) = XDOT(1,J) - XIDOT(1)
            DV(2) = XDOT(2,J) - XIDOT(2)
            DV(3) = XDOT(3,J) - XIDOT(3)
        
            RIJ2 = A1*A1 + A2*A2 + A3*A3
            DR2I = 1.0/RIJ2
            DR3I = BODY(J)*DR2I*SQRT(DR2I)
            DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
            DRDP = 3.0*DRDV*DR2I
                  
            DFIRR(1) = DFIRR(1) - A1*DR3I
            DFIRR(2) = DFIRR(2) - A2*DR3I
            DFIRR(3) = DFIRR(3) - A3*DR3I
            DFD(1) = DFD(1) - (DV(1) - A1*DRDP)*DR3I
            DFD(2) = DFD(2) - (DV(2) - A2*DRDP)*DR3I
            DFD(3) = DFD(3) - (DV(3) - A3*DRDP)*DR3I
        END DO

    !     Add new neighbor particle to irregular force differences
        DO L = 1,NBGAIN
            J = JJLIST(NNB0+L)
                     
            A1 = X(1,J) - XI(1)
            A2 = X(2,J) - XI(2)
            A3 = X(3,J) - XI(3)
            DV(1) = XDOT(1,J) - XIDOT(1)
            DV(2) = XDOT(2,J) - XIDOT(2)
            DV(3) = XDOT(3,J) - XIDOT(3)
        
            RIJ2 = A1*A1 + A2*A2 + A3*A3
            DR2I = 1.0/RIJ2
            DR3I = BODY(J)*DR2I*SQRT(DR2I)
            DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
            DRDP = 3.0*DRDV*DR2I
                  
            DFIRR(1) = DFIRR(1) + A1*DR3I
            DFIRR(2) = DFIRR(2) + A2*DR3I
            DFIRR(3) = DFIRR(3) + A3*DR3I
            DFD(1) = DFD(1) + (DV(1) - A1*DRDP)*DR3I
            DFD(2) = DFD(2) + (DV(2) - A2*DRDP)*DR3I
            DFD(3) = DFD(3) + (DV(3) - A3*DRDP)*DR3I
        END DO
    
        DTSQ = DTR**2
        DT6 = 6.0/(DTR*DTSQ)
        DT2 = 2.0/DTSQ
        DTSQ12 = ONE12*DTSQ
        DTR13 = ONE3*DTR
        T0R(I) = TIME
    
    !       Suppress the corrector for large time-step ratios (experimental).
        IF (STEP(I) < 5.0D0*DTMIN .AND. DTR > 50.0*STEP(I)) THEN
            DTR13 = 0.0
            DTSQ12 = 0.0
        END IF
    
    !     LCHECK = NBLOSS+NBGAIN.EQ.0.AND.
    !    *   ((STEP(I).LT.1.D-09.AND.DFIRR.GT.1.0D-9).OR.
    !    *    (STEPR(I).LT.DTMIN/1.D1.AND.DFREG.GT.1.0D-9))
    !     IF (LCHECK) THEN
    !         DO K=1,3
    !             DFI(K) = (FIRR(K)-FI(K,I))**2
    !             DFRG(K) = (FREG(K)-FR(K,I))**2
    !         END DO
    !         DFIRR = DSQRT(DFI(1)+DFI(2)+DFI(3))**2
    !         DFREG = DSQRT(DFRG(1)+DFRG(2)+DFRG(3))**2
    !     IF (rank.eq.0)
    !    *WRITE (21,740)  TIME, I, NNB, DFIRR, DFREG, STEP(I), STEPR(I)
    ! 740 FORMAT ('-----------------------------------',/,
    !    *     ' CHECK! T I NB DFI DFR DTi,r',1P,E12.4,2I6,4E12.4)
    !         WRITE(21,*)NBLOSS,NBGAIN
    !         DO K=1,NNB
    !         J = NLIST(K)
    !         IF (J.GT.N) THEN
    !             JPAIR = J-N
    !             J2 = 2*JPAIR-1
    !         if(rank.eq.0)
    !    *    WRITE (21,742) J,J2,LIST(1,J2),GAMMA(JPAIR)
    ! 742 FORMAT (' CHECK  J J2 L Pert         ',3I6,1P,E12.4)
    !         END IF
    !         END DO
    !         IF (I.GT.N)THEN
    !             IPAIR = I-N
    !             I2 = 2*IPAIR-1
    !         if(rank.eq.0)
    !    *    WRITE (21,741) I,I2,LIST(1,I2),GAMMA(IPAIR)
    ! 741 FORMAT (' CHECK  I I2 L Pert          ',3I6,1P,E12.4)
    !     CALL FLUSH(21)
    !         END IF
    !     END IF
    !       Test for inconsistencies of nbint, regint and cmfreg Nov 2010
    !       Include the corrector and save F, FI, FR & polynomial derivatives.
        DO 75 K = 1,3
        !       Subtract change of neighbour force to form actual first derivative.
            DFR = FR(K,I) - FREG(K) - DFIRR(K)
            FDR0 = FDR(K) + DFD(K)
        !     --10/15/13 18:44-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$          if(name(i).eq.8689) then
        !$$$             write(102+rank,*)'K',K,'I',I,'N',NAME(I),'FI',FI(K,I),
        !$$$     &        'DFIRR',DFIRR(K),'FR',FR(K,I),'FREG',FREG(K),
        !$$$     &        'FR-FREG',FR(K,I)-FREG(K),'DFR',DFR,'FD',FIDOT(K,I),
        !$$$     &        'DFD',DFD(K),'FDR',FDR(K),'RS',RS(I),
        !$$$     &        'TIME',TIME,'NBGAIN',nbgain,'NBLOSS',nbloss
        !$$$             call flush(102+rank)
        !$$$          if(nbgain.eq.0.and.nbloss.eq.0.and.DFIRR(K).NE.0.D0)
        !$$$     &         call abort()
        !$$$          end if
        !     --10/15/13 18:44-lwang-end----------------------------------------*
        
            FRD = FRDOT(K,I)
            SUM = FRD + FDR0
            AT3 = 2.0*DFR + DTR*SUM
            BT2 = -3.0*DFR - DTR*(SUM + FRD)
        !       Use here new variables for consistency in parallel execution (R.Sp.)
            XN(K,I) = X0(K,I) + (0.6*AT3 + BT2)*DTSQ12
            XNDOT(K,I) = X0DOT(K,I) + (0.75*AT3 + BT2)*DTR13
        
            FI(K,I) = FI(K,I) + DFIRR(K)
            FR(K,I) = FREG(K)
            FIDOT(K,I) = FIDOT(K,I) + DFD(K)
            FRDOT(K,I) = FDR(K)

        !     Do tidal force correction for FIDOT
            IF (KZ(14) <= 2) THEN
                IF(K == 2) FIDOT(1,I) = FIDOT(1,I) - TIDAL(4)*DFR
                IF(K == 1) FIDOT(2,I) = FIDOT(2,I) + TIDAL(4)*DFR
            END IF
        
            D0(K,I) = FI(K,I)
            D0R(K,I) = FREG(K)
            D1R(K,I) = FDR(K)
            D3R(K,I) = AT3*DT6
            D2R(K,I) = (3.0*AT3 + BT2)*DT2
        75 END DO
    
    !     Form rate of tidal energy change during last regular step.
        #ifdef TT
    !** FlorentR - include the case of tidal tensor (really needed?)
        IF (KZ(14) == 3 .OR. KZ(14) == 9) THEN
        !** FRenaud
            #else
            IF (KZ(14) == 3) THEN
                #endif
                WDOT = 0.0
                W2DOT = 0.0
            !     W3DOT = 0.0
                DO 24 K = 1,3
                    PX = FREG(K) - FRX(K)
                    DPX = FDR(K) - FDX(K)
                    WDOT = WDOT + XIDOT(K)*PX
                    W2DOT = W2DOT + (FREG(K) + FI(K,I))*PX + XIDOT(K)*DPX
                !     W3DOT = W3DOT + 2.0*(FREG(K) + FIRR(K))*DPX +
                !     &                            (FDR(K) + FD(K))*PX
                24 END DO
            ! omp critical
            !       Note: second-order term derived by Douglas Heggie (Aug/03).
            !       Accumulate tidal energy change for general galactic potential.
            !       Note: Taylor series at end of interval with negative argument.
            !         ETIDE = ETIDE - BODY(I)*((ONE6*W3DOT*DTR - 0.5*W2DOT)*DTR +
            !    &                                                   WDOT)*DTR
                ETIDE = ETIDE + BODY(I)*(0.5*W2DOT*DTR - WDOT)*DTR
            ! omp end critical
            !       Note: integral of Taylor series for V*P using final values.
            END IF

        !     Check optional force polynomial corrections due to neighbour changes.
            IF (KZ(38) > 0) THEN
                IF (KZ(38) == 1) THEN
                    CALL FPCORR(I,NBLOSS,NBGAIN,XI,XIDOT,JJLIST)
                ELSE IF (KZ(38) == 2 .AND. DF2 > 0.0001*FR2) THEN
                    CALL FPCORR(I,NBLOSS,NBGAIN,XI,XIDOT,JJLIST)
                ELSE IF (NBLOSS > 0 .AND. KZ(38) == 3) THEN
                !       Ignore corrections for gains when treating the most critical losses.
                    NBGAIN = 0
                    CALL FPCORR(I,NBLOSS,NBGAIN,XI,XIDOT,JJLIST)
                END IF
            END IF
        !     --03/19/14 10:51-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$      if (name(i).eq.1886) then
        !$$$         j=i
        !$$$         write(105+rank,*),'reg f i',j,'n',name(j),'x0',x0(1,j),
        !$$$     *        'x0dot',x0dot(1,j),
        !$$$     *        'firr',firr(1),'fd',fd(1),'freg',freg(1),
        !$$$     &        'fdr',fdr(1),
        !$$$     *        't0',t0(j),'stepr',stepr(j),'dtr',dtr,
        !$$$     *        'body',body(j),'time',time,'list',list(1,j)
        !$$$         call flush(105+rank)
        !$$$      end if
        !     --03/19/14 10:51-lwang-end----------------------------------------*

        !     Copy new neighbour list if different from old list.
            IF (NBLOSS + NBGAIN > 0) THEN
                LIST(1,I) = NNB
                DO 80 L = 2,NNB+1
                    LIST(L,I) = NLIST(L)
                80 END DO
            ! omp critical
                NBFLUX = NBFLUX + NBLOSS + NBGAIN
            ! omp end critical
            END IF
        
        !       Obtain new regular integration step using composite expression.
        !       STEPR = (ETAR*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
            FR2 = FREG(1)**2 + FREG(2)**2 + FREG(3)**2
            W1 = FDR(1)**2 + FDR(2)**2 + FDR(3)**2
            W2 = D2R(1,I)**2 + D2R(2,I)**2 + D2R(3,I)**2
            W3 = D3R(1,I)**2 + D3R(2,I)**2 + D3R(3,I)**2
        
        !       Form new step by relative criterion.
            W0 = (SQRT(FR2*W2) + W1)/(SQRT(W1*W3) + W2)
            W0 = ETAR*W0
            TTMP = SQRT(W0)
            DT0 = TTMP
        
        !       Winston Sweatman's suggestion
        !     DVV = (XDOT(1,I)-X0DOT(1,I))**2 + (XDOT(2,I)-X0DOT(2,I))**2 +
        !    &     (XDOT(3,I)-X0DOT(3,I))**2
        !     FFD = FREG(1)**2 + FREG(2)**2 + FREG(3)**2
        !     ETARW = ETAR
        !     TTMPW = ETARW*DVV*BODY(I)/FFD
        
        !     PRINT*,' Reg I=',I,' TTMP,TTMPW,RATIO=',
        !    &  TTMP,TTMPW,TTMP/TTMPW
        
        !     IF(TTMP.GT.TTMPW)THEN
        !     IGT = IGT + 1
        !     ELSE
        !     ILE = ILE + 1
        !     END IF
        !     IF(MOD(IGT+ILE,100).EQ.0)PRINT*,' irr IGT,ILE=',IGT,ILE
        
        !     TTMP = MAX(TTMPW,TTMP)
        !     TTR = TSTEP(FREG,FDR,D2R(1,I),D3R(1,I),ETAR)
        
        !       Adopt FAC*MIN(FREG,FIRR) (or tidal force) for convergence test.
        !$$$      FAC = MIN(ETAR,0.04D0)
        !$$$      IF (KZ(14).EQ.0) THEN
        !$$$          FI2 = FIRR(1)**2 + FIRR(2)**2 + FIRR(3)**2
        !$$$          DF2 = FAC**2*MIN(FR2,FI2)
        !$$$*       Ignore irregular force criterion if no neighbours.
        !$$$          IF (NNB.EQ.0) DF2 = FAC**2*FR2
        !$$$      ELSE IF (KZ(14).EQ.1) THEN
        !$$$          W0 = (TIDAL(1)*XI(1))**2
        !$$$          DF2 = FAC**2*MAX(FR2,W0)
        !$$$      ELSE
        !$$$          DF2 = FAC**2*FR2
        !$$$      END IF
        !$$$*
        !$$$*       Obtain regular force change using twice the predicted step.
        !$$$      DTC = 2.0*TTMP
        !$$$      S2 = 0.5*DTC
        !$$$      S3 = ONE3*DTC
        !$$$      W0 = 0.0
        !$$$      DO 105 K = 1,3
        !$$$          W2 = ((D3R(K,I)*S3 + D2R(K,I))*S2 + FDR(K))*DTC
        !$$$          W0 = W0 + W2**2
        !$$$  105 CONTINUE
        !$$$*
        !$$$*       See whether regular step can be increased by factor 2.
        !$$$      IF (W0.LT.DF2) THEN
        !$$$          TTMP = DTC
        !$$$      END IF
        !$$$*
        !$$$*       Impose a smooth step reduction inside compact core.
        !$$$      IF (NC.LT.50.AND.RI2.LT.RC2) THEN
        !$$$          TTMP = TTMP*MIN(1.0D0,0.5D0*(1.0D0 + RI2*RC2IN))
        !$$$      END IF
        
        !       Select discrete value (increased by 2, decreased by 2 or unchanged).
        !     --09/26/13 21:14-lwang-change-------------------------------------*
        !**** Note: To increase the regular stepr if dt0 is much larger than ttmp
            IF (TTMP > 2.0*STEPR(I)) THEN
                IF (DMOD(TIME,2.0D0*DTR) == 0.0D0) THEN
                    TTMP = MIN(2.0*DTR,SMAX)
                !       Include factor 4 increase for rare cases of too small regular steps.
                    106 IF (DT0 > 2.0*TTMP .AND. 2.0*TTMP < SMAX .AND. &
                    DMOD(TIME,2.0D0*TTMP) == 0.0D0) THEN
                        TTMP = 2.0*TTMP
                        GO TO 106
                    END IF
                ELSE
                    TTMP = STEPR(I)
                END IF
            ELSE IF (TTMP < DTR) THEN
                TTMP = 0.5*DTR
            !       Allow a second reduction to prevent spurious contributions.
                IF (TTMP > DT0) THEN
                    TTMP = 0.5*TTMP
                END IF
            ELSE
                TTMP = STEPR(I)
            END IF

        !**** Note: For large N------------------------------------------------*
        !       Include convergence test for increasing step (Makino, ApJ, 369, 200).
            IF (TTMP > 0.1 .AND. TTMP > DTR) THEN
                DV2 = 0.0
                DO 90 K = 1,3
                    DV2 = DV2 + (XIDOT(K) - XNDOT(K,I))**2
                90 END DO
            !       Employ Jun Makino criterion to avoid over-shooting (cf. Book, 2.16).
                DTJ = DTR*(1.0D-04*DTR**2*FR2/DV2)**0.1
                IF (DTJ < TTMP) TTMP = DTR
            END IF
        !     --09/26/13 21:14-lwang-end----------------------------------------*
        !$$$      IF (TTMP .GT. 2.0*STEPR(I)) THEN
        !$$$         IF (DMOD(TIME,2.0*STEPR(I)) .EQ. 0.0D0) THEN
        !$$$              TTMP = MIN(2.0*STEPR(I),1.D0)
        !$$$          ELSE
        !$$$              TTMP = STEPR(I)
        !$$$          END IF
        !$$$      ELSE IF (TTMP .LT. STEPR(I)) THEN
        !$$$          TTMP = 0.5*STEPR(I)
        !$$$      ELSE
        !$$$          TTMP = STEPR(I)
        !$$$      END IF
        
        !     --11/04/13 23:10-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
            if (stepr(i) >= 1E2*dt0) then
                print*,rank,'Warning!: Regular step jump!: I',I,'N',name(i), &
                'FI',Fi(1,i),'FD',Fidot(1,i),'Freg',Freg(1),'FDR',FDR(1), &
                'ratio',dt0/stepr(i),'dt0',dt0,'stepr ',stepr(i), &
                't',TIME,'NNB',NNB
                call flush(6)
            !$$$         print*,'NNB',NNB,'LIST',(NAME(LIST(K,I)),K=2,NNB+1)
            !$$$         IF(NAME(i).EQ.1674) stop
            end if
        !$$$      if(name(i).eq.1674) then
        !$$$         print*,rank,'I',I,'N',name(i),'firr',fi(1,i),'fd',fidot(1,i),
        !$$$     &        'freg',freg(1),'fdr',fdr(1),
        !$$$     &        'stepr',stepr(i),'ttmp',ttmp,'t',time,'NNB',NNB
        !$$$         print*,'NNB',NNB,'LIST',(NAME(LIST(K,I)),K=2,NNB+1)
        !$$$         call flush(6)
        !$$$      end if
        !     --11/04/13 23:10-lwang-end----------------------------------------*
        
        !     ensure regular step not larger than smax
            TTMP = MIN(TTMP,SMAX)
        !       Set new regular step and reduce irregular step if STEPR < STEP.
            STEPR(I) = TTMP

        !     To avoid large Fdot causing large energy error for special case.(lwang)
            if(nnb0 == 0 .AND. nnb >= 0) then
                step(i) = 0.125*step(i)
                timenw(i) = t0(i) + step(i)
            end if

        !     --03/30/13 22:34-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$      if(rank.eq.0.and.name(i).eq.595) then
        !$$$         print*,'STEPR(595):',stepr(i),'time',time,'t0',t0(i),
        !$$$     &        'dm',dmod(time,stepr(i))
        !$$$         call flush(6)
        !$$$      end if
        !$$$      if (time.ge.60.37.and.name(i).eq.16107) then
        !$$$         j=i
        !$$$         write(105+rank,*),'reg',j,name(j),'x0',(x0(kk,j),kk=1,3),
        !$$$     *        'x0dot',(x0dot(kk,j),kk=1,3),'xn',(xn(kk,j),kk=1,3),
        !$$$     *        'xndot',(xndot(kk,j),kk=1,3),
        !$$$     *        't0',t0(j),'step',step(j),'dt0',dt0,
        !$$$     *        'stepr',stepr(j),'nbgain',nbgain,'nbloss',nbloss,
        !$$$     *        'f',(f(kk,j),kk=1,3),'fdot',(fdot(kk,j),kk=1,3),
        !$$$     *        'fi',(fi(kk,j),kk=1,3),'fidot',(fidot(kk,j),kk=1,3),
        !$$$     *        'fr',(fr(kk,j),kk=1,3),'frdot',(frdot(kk,j),kk=1,3),
        !$$$     *        'd0',(d0(kk,j),kk=1,3),'d1',(d1(kk,j),kk=1,3),
        !$$$     *        'd2',(d2(kk,j),kk=1,3),'d3',(d3(kk,j),kk=1,3),
        !$$$     *        'd0r',(d0r(kk,j),kk=1,3),'d1r',(d1r(kk,j),kk=1,3),
        !$$$     *        'd2r',(d2r(kk,j),kk=1,3),'d3r',(d3r(kk,j),kk=1,3),
        !$$$     *        'body',body(j),'time',time,'list',list(1,j)
        !$$$         call flush(105+rank)
        !$$$      end if
        !     --03/30/13 22:34-lwang-end----------------------------------------*

            110 IF (TTMP < STEP(I)) THEN
                STEP(I) = 0.5D0*STEP(I)
                TIMENW(I) = T0(I) + STEP(I)
            ! omp critical
                NICONV = NICONV + 1
            ! omp end critical
                GO TO 110
            END IF
        
            RETURN
        
            END SUBROUTINE REGCOR_GPU
