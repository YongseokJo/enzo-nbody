    SUBROUTINE KSINT(I1)


!       Regularized integration.
!       ------------------------

    Include 'kspars.h'
    USE POINTERS
    INCLUDE 'common6.h'
          
    include 'timing.h'
!     Safe for parallel, no value change, used value: tmdis,namem,nameg,kstarm
    COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX), &
    HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX), &
    NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
!     No value change, used value: TOSC,namec
    COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX), &
    BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX), &
    ZP(NTMAX),ES(NTMAX),CZ(2,NTMAX),IOSC(NTMAX), &
    NAMEC(NTMAX)
!     Safe for parallel, no value change
    COMMON/SLOW0/  RANGE,ISLOW(10)
!      COMMON/GAMDOT/  DGAM
    REAL*8 ::  UI(4),UIDOT(4),XI(6),VI(6),FP(6),FD(6)
    LOGICAL :: IQ
!     --01/02/14 20:43-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
!$$$      REAL*8 TPRED
!$$$      LOGICAL ipredall
!     --01/02/14 20:43-lwang-end----------------------------------------*


!       Set second component, pair index & c.m. index.
    I2 = I1 + 1
    IPAIR = KVEC(I1)
    I = N + IPAIR
    IPHASE = 0
    call jpred(i,time,time)
!     --09/25/13 14:22-lwang-debug--------------------------------------*
!**** Note: To test whether the value are same for all nodes at beginnging of ksint*
!$$$      if(i1.eq.15967.and.tblock.eq.0.73089027404785156) THEN
!$$$         write(120+rank,*),'KSINT I1',I1,'n',name(i1),'IPAIR',ipair,
!$$$     &        'U',U(1,IPAIR),'U0',U0(1,IPAIR),
!$$$     &        'UDOT',UDOT(1,ipair),'FU',FU(1,ipair),
!$$$     &        'FT',FUDOT(1,ipair),'FT2',FUDOT2(1,ipair),
!$$$     &        'FT3',FUDOT3(1,ipair),'FP0',FP0(1,Ipair),
!$$$     &        'FD0',FD0(1,ipair),'H',H(ipair),'HT',HDOT(ipair),
!$$$     &        'HT2',HDOT2(ipair),'HT3',HDOT3(ipair),
!$$$     &        'HT4',HDOT4(ipair),'dt',DTAU(ipair),
!$$$     &        'TD2',TDOT2(ipair),'TD3',TDOT3(ipair),
!$$$     &        'R',R(ipair),'R0',r0(ipair),
!$$$     &        'Gamma',GAMMA(ipair),'H0',h0(ipair),
!$$$     &        'SF',sf(1,ipair),'kslow',kslow(ipair),
!$$$     &        'nnb',LIST(1,i1),'list',list(2,i1)
!$$$     &        ,'step',step(i1),'t0',t0(i1),'tblock',tblock
!$$$         print*,rank,'ksint I1',i1,'X',x(1,i1),'xdot',xdot(1,i1),
!$$$     &        'UDOT',udot(1,ipair),'tblock',tblock
!$$$         call flush(120+rank)
!$$$         call mpi_barrier(MPI_COMM_WORLD,ierr)
!$$$      end if
!     --09/25/13 14:23-lwang-end----------------------------------------*

!       Define perturber membership & inverse c.m. mass.
    NNB0 = LIST(1,I1)
    BODYIN = 1.0/BODY(I)

!       Check for further unperturbed motion.
    IF (NNB0 == 0 .AND. H(IPAIR) < 0.0) THEN
    !$$$*       Include possible eccentricity modulation of hierarchical binary.
    !$$$          IF (NAME(I).LT.0) THEN
    !$$$              IM = 0
    !$$$              DO 1 K = 1,NMERGE
    !$$$                  IF (NAMEM(K).EQ.NAME(I)) IM = K
    !$$$    1         CONTINUE
    !$$$              IF (IM.GT.0.AND.TIME.GT.TMDIS(IM)) THEN
    !$$$                  SEMI = -0.5*BODY(I)/H(IPAIR)
    !$$$                  ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
    !$$$     &                                   TDOT2(IPAIR)**2/(BODY(I)*SEMI)
    !$$$                  RP = SEMI*(1.0 - SQRT(ECC2))*(1.0 - 2.0*GAMMA(IPAIR))
    !$$$                  IF (RP.LT.R0(IPAIR)) THEN
    !$$$                      if(rank.eq.0)
    !$$$     &                WRITE (77,4)  NAME(I1), SQRT(ECC2), RP, R0(IPAIR)
    !$$$    4                 FORMAT (' TMDIS TERM    NM E RP R0',I7,1P,3E10.2)
    !$$$                      CALL FLUSH(77)
    !$$$                      INSTAB = INSTAB + 1
    !$$$                      GO TO 90
    !$$$                  ELSE IF (KZ(27).EQ.2) THEN
    !$$$                      CALL ECCMOD(I,ITERM)
    !$$$*       Update time on termination to prevent looping.
    !$$$                      IF (ITERM.GT.0) THEN
    !$$$                          T0(I1) = TIME
    !$$$                          GO TO 90
    !$$$                      END IF
    !$$$                  END IF
    !$$$*       Check any inner and outer circularizing binary using NAMEM & NAMEG.
    !$$$                  DO 3 K = 1,NCHAOS
    !$$$                      IF (NAMEC(K).EQ.NZERO - NAMEM(IM)) THEN
    !$$$*       Update unperturbed binary if T - TOSC > 10 Myr (cf. IMPACT & DECIDE).
    !$$$                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) THEN
    !$$$                              T0(I1) = TIME
    !$$$                              GO TO 90
    !$$$                          END IF
    !$$$                      END IF
    !$$$                      IF (NAMEC(K).EQ.NAMEG(IM)) THEN
    !$$$                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) THEN
    !$$$                              T0(I1) = TIME
    !$$$                              GO TO 90
    !$$$                          END IF
    !$$$                      END IF
    !$$$    3             CONTINUE
    !$$$              END IF
    !$$$          END IF
        call cputim(tt1)
        CALL UNPERT(IPAIR)
        call cputim(tt2)
        ttup = ttup + (tt2-tt1)*60.
        GO TO 100
    END IF

!       Perform KS prediction of U, UDOT & H.
    CALL KSPRED(IPAIR,I1,I,BODYIN,UI,UIDOT,XI,VI)

!     --01/02/14 13:42-lwang-debug--------------------------------------*
!**** Note: check predicted value--------------------------------------**
!$$$      if(tblock.eq.29.053588867187500.and.name(i1).eq.1759) then
!$$$         write(125+rank,*)'KSPRED I1',I1,'IPAIR',ipair,'NNB',list(1,i1),
!$$$     &        'XI',XI(1),'VI',VI(1),'UI',UI(1),'UID',UIDOT(1),'B',BODYIN
!$$$     &        ,'t',time
!$$$         do L = 2, list(1,i1)+1
!$$$            K = list(l,i1)
!$$$            write(125+rank,*),k,'n',name(k),'x0',x0(1,k),'v0',x0dot(1,k)
!$$$     &           ,'f',f(1,k),'fdot',fdot(1,k),'x',x(1,k),'v',xdot(1,k),
!$$$     &           'nnb',list(1,k),'t0',t0(k),'tpred',tpred(k),
!$$$     &           'time',time,'flag',ipredall
!$$$            if(k.gt.n) then
!$$$               j = 2*(k - n)-1
!$$$               write(125+rank,*),j,'n',name(j),'x',x(1,j),'v',xdot(1,j),
!$$$     &              'nnb',list(1,k),'npert',list(1,j),'u',u(1,k-n),
!$$$     &              'ud',udot(1,k-n),
!$$$     &              'time',time
!$$$            end if
!$$$         end do
!$$$         call flush(125+rank)
!$$$      end if
!     --01/02/14 13:42-lwang-end----------------------------------------*
!       Obtain the perturbing force & derivative.
    CALL KSPERT(I1,NNB0,XI,VI,FP,FD)

!     --09/25/13 14:22-lwang-debug--------------------------------------*
!**** Note: To test after ksperts
!$$$      if(tblock.eq.29.053588867187500.and.name(i1).eq.1759) then
!$$$         write(125+rank,*)'KSPERT I1',I1,'IPAIR',ipair,
!$$$     &        'U',U(1,IPAIR),'U0',U0(1,IPAIR),
!$$$     &        'UDOT',UDOT(1,ipair),'FU',FU(1,ipair),
!$$$     &        'FT',FUDOT(1,ipair),'FT2',FUDOT2(1,ipair),
!$$$     &        'FT3',FUDOT3(1,ipair),'FP0',FP0(1,Ipair),
!$$$     &        'FD0',FD0(1,ipair),'H',H(ipair),'HT',HDOT(ipair),
!$$$     &        'HT2',HDOT2(ipair),'HT3',HDOT3(ipair),
!$$$     &        'HT4',HDOT4(ipair),'dt',DTAU(ipair),
!$$$     &        'TD2',TDOT2(ipair),'TD3',TDOT3(ipair),
!$$$     &        'R',R(ipair),'R0',r0(ipair),
!$$$     &        'Gamma',GAMMA(ipair),'H0',h0(ipair),
!$$$     &        'SF',sf(1,ipair),'kslow',kslow(ipair),
!$$$     &        'nnb',LIST(1,i1),'list',list(2,i1)
!$$$     &        ,'step',step(i1),'t0',t0(i1),'tblock',tblock
!$$$         write(125+rank,*),'KSPERT I1',i1,
!$$$     &        'UI',UI(1),'UIDOT',uidot(1),'fp',fp(1),
!$$$     &        'fd',fd(1),'tblock',tblock
!$$$         call flush(125+rank)
!$$$         call flush(6)
!$$$         call mpi_barrier(MPI_COMM_WORLD,ierr)
!$$$      end if
!     --09/25/13 14:23-lwang-end----------------------------------------*
!       Save old radial velocity & relative perturbation and set new GAMMA.
    RDOT = TDOT2(IPAIR)
!     G0 = GAMMA(IPAIR)
    GI = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)*R(IPAIR)**2*BODYIN
    GAMMA(IPAIR) = GI

!     --09/25/13 14:22-lwang-debug--------------------------------------*
!**** Note: To test after correct
!      if(i1.eq.1605.and.time.le.29.047285281042758) then
!$$$      if(tblock.eq.29.053588867187500.and.name(i1).eq.1759) then
!$$$         write(125+rank,*)'KSINT I1',I1,'IPAIR',ipair,
!$$$     &        'UI',UI(1),'U0',U0(1,IPAIR),
!$$$     &        'UIDOT',UIDOT(1),'FU',FU(1,ipair),
!$$$     &        'FT',FUDOT(1,ipair),'FT2',FUDOT2(1,ipair),
!$$$     &        'FT3',FUDOT3(1,ipair),'FP',FP(1),
!$$$     &        'FD',FD(1),'H',H(ipair),'HT',HDOT(ipair),
!$$$     &        'HT2',HDOT2(ipair),'HT3',HDOT3(ipair),
!$$$     &        'HT4',HDOT4(ipair),'dt',DTAU(ipair),
!$$$     &        'TD2',TDOT2(ipair),'TD3',TDOT3(ipair),
!$$$     &        'R',R(ipair),'R0',r0(ipair),
!$$$     &        'Gamma',GAMMA(ipair),'H0',h0(ipair),
!$$$     &        'SF',sf(1,ipair),'kslow',kslow(ipair),
!$$$     &        'nnb',LIST(1,i1),'list',list(2,i1)
!$$$     &        ,'step',step(i1),'t0',t0(i1),'tblock',tblock,'t',time
!$$$         print*,rank,'ksint I1',i1,'ipair',ipair,
!$$$     &        'H0',H0(ipair),'UDOT',udot(1,ipair),
!$$$     &        't',time,'t0',t0(i1)
!$$$         call flush(125+rank)
!$$$         call mpi_barrier(MPI_COMM_WORLD,ierr)
!$$$      end if
!       Apply the Hermite corrector.
    CALL KSCORR(IPAIR,UI,UIDOT,FP,FD,TD2,TDOT4,TDOT5,TDOT6)
!       Increase regularization time-step counter and update the time.
    NSTEPU = NSTEPU + 1
    T0(I1) = TIME
!     --09/25/13 14:22-lwang-debug--------------------------------------*
!**** Note: To test after correct
!$$$      if(tblock.eq.29.053588867187500.and.name(i1).eq.1759) then
!$$$         write(125+rank,*)'KSCORR I1',I1,'IPAIR',ipair,
!$$$     &        'U',U(1,IPAIR),'U0',U0(1,IPAIR),
!$$$     &        'UDOT',UDOT(1,ipair),'FU',FU(1,ipair),
!$$$     &        'FT',FUDOT(1,ipair),'FT2',FUDOT2(1,ipair),
!$$$     &        'FT3',FUDOT3(1,ipair),'FP0',FP0(1,Ipair),
!$$$     &        'FD0',FD0(1,ipair),'H',H(ipair),'HT',HDOT(ipair),
!$$$     &        'HT2',HDOT2(ipair),'HT3',HDOT3(ipair),
!$$$     &        'HT4',HDOT4(ipair),'dt',DTAU(ipair),
!$$$     &        'TD2',TDOT2(ipair),'TD3',TDOT3(ipair),
!$$$     &        'R',R(ipair),'R0',r0(ipair),
!$$$     &        'Gamma',GAMMA(ipair),'H0',h0(ipair),
!$$$     &        'SF',sf(1,ipair),'kslow',kslow(ipair),
!$$$     &        'nnb',LIST(1,i1),'list',list(2,i1)
!$$$     &        ,'step',step(i1),'t0',t0(i1),'tblock',tblock,'t',time
!$$$         print*,rank,'kscorr I1',i1,'ipair',ipair,
!$$$     &        'U',u(1,ipair),'UDOT',udot(1,ipair),
!$$$     &        't',time,'t0',t0(i1)
!$$$         call flush(125+rank)
!$$$         call mpi_barrier(MPI_COMM_WORLD,ierr)
!$$$         if(time.ge.29.047067570490370) stop
!$$$      end if
!     --09/25/13 14:23-lwang-end----------------------------------------*
!$$$*       Check for early return during termination (called from KSTERM).
!$$$      IF (IPHASE.NE.0) GO TO 100

!       Define useful scalars.
    RI = R(IPAIR)
    HI = H(IPAIR)

!       Initialize termination indicator and check for large perturbation.
    IQ = .FALSE.
!     Reduce the criterion from 0.25 to 0.1 to reduce the energy error
    IF (GI > 0.10) GO TO 2
    IF (GI < 0.03) THEN
        JCOMP = 0
        GO TO 20
    END IF
    CALL FLYBY(I,ITERM)
    IF (ITERM == 0 .AND. KSTAR(I) < 0) THEN
        IQ = .TRUE.
    END IF
    IF (ITERM == 1) THEN
    !       Delay chain regularization search until near end of block-step.
        IF (KZ(15) > 0 .AND. TIME + STEP(I1) > TBLOCK) THEN
            CALL IMPACT(I)
            IF (IPHASE > 0) GO TO 100
        END IF
    ELSE IF (ITERM == 2) THEN
        IQ = .TRUE.
        GO TO 20
    ELSE
        GO TO 20
    END IF

!       Find the dominant body for large perturbations.
    2 S = 4.0*STEP(I)
    FMAX = BODY(I)/RI**2
!       Initialize JCOMP for prediction and optional diagnostics in KSTERM.
    JCOMP = 0
    DO 10 L = 2,NNB0+1
        J = LIST(L,I1)
    !       Only search bodies within four times of the c.m. time-step.
        IF (STEP(J) > S) GO TO 10
    !       Compare strong perturber and either component with current pair.
        DO 5 K = I1,I2
            RIJ2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 + &
            (X(3,J) - X(3,K))**2
            IF (BODY(J) + BODY(K) > RIJ2*FMAX) JCOMP = J
        5 END DO
    10 END DO

!       Set termination if strong perturber <= N forms dominant pair.
    IF (JCOMP > 0 .OR. GI > 1.0) THEN
    !       Check optional binary search.
    !         IF (KZ(4).GT.0) THEN
    !             DGAM = GI - G0
    !             K = KZ(4)
    !             CALL EVOLVE(IPAIR,K)
    !         END IF
    !       Terminate the binary when close particle is KS pair
        IF (JCOMP > N) IQ = .TRUE. 
        IF (JCOMP <= N .OR. GI > 1.0) IQ = .TRUE. 
    END IF

!       Check termination of hyperbolic encounter (R > R0 or R > 2*RMIN).
    20 IF (HI > 0.0D0 .AND. NAME(I) > 0) THEN
        IF ((RI > R0(IPAIR) .AND. GI > GMAX) .OR. RI > 2.0*RMIN .OR. &
        (GI > 0.5 .AND. TD2 > 0.0)) THEN
        !       Skip termination delay in case of velocity kick (cf. routine KSTERM).
            IF (HI < 100.0 .OR. GI > 0.1 .OR. RI > 5.0*RMIN) THEN
                IQ = .TRUE.
            END IF
        END IF
    END IF

!       Choose basic regularized step using binding energy or lower limit.
    IF (ABS(HI) > ECLOSE) THEN
        W1 = 0.5/ABS(HI)
    ELSE
        W1 = R0(IPAIR)*BODYIN
        W2 = 0.5/ABS(HI)
        W1 = MIN(W1,W2)
        IF (RI > R0(IPAIR)) W1 = W1*R0(IPAIR)/RI
    !       Maximum square step for soft binaries & weak hyperbolic pairs.
        IF (HI < 0.0D0) THEN
        !       Include case of hard binary with massive components or merger.
            W2 = -0.5/HI
            IF (NAME(I) < 0) THEN
                W1 = RI*BODYIN
            END IF
            W1 = MIN(W2,W1)
        END IF
    END IF

!       Include perturbation factor in predicted step.
    IF (GI < 0.0005) THEN
    !       Use second-order expansion of cube root for small perturbations.
        W3 = 333.3*GI
        W2 = SQRT(W1)/(1.0 + W3*(1.0 - W3))
    ELSE
        W3 = 1.0 + 1000.0*GI
        W2 = SQRT(W1)/W3**0.3333
    END IF

!       Form new regularized step (include inertial factor).
    DTU = ETAU*W2
    DTU = MIN(1.2*DTAU(IPAIR),DTU)

!       Include convergence criterion DH = H'*DTU + H''*DTU**2/2 = 0.001*|H|.
    IF (GI > 1.0D-04) THEN
        DH = 1.0E-03*MAX(ABS(H(IPAIR)),ECLOSE)
        XF = 2.0*DH/ABS(HDOT2(IPAIR))
        YF = HDOT(IPAIR)/HDOT2(IPAIR)
        DTU1 = SQRT(XF + YF**2) - ABS(YF)
        DTU = MIN(DTU1,DTU)
    END IF

!       Check pericentre step reduction for perturbed spiral.
    IF (KSTAR(I) == -2 .AND. TD2 < 0.0) THEN
        SEMI = -0.5*BODY(I)/HI
        IF (RI < SEMI) THEN
        !       Predict radial velocity for step DTU (note scaled coefficients).
            RD = ((ONE3*TDOT5*DTU + TDOT4)*DTU + &
            TDOT3(IPAIR))*DTU + 2.0*TD2
        !       Adopt pericentre step with 1% safety factor (small TDOT4 is OK).
            IF (RD > 0.0) THEN
                A2 = 0.5*TDOT3(IPAIR)/TDOT4
                DTU1 = SQRT(A2**2 - 2.0*TD2/TDOT4) - A2
                DTU1 = 1.01*MAX(DTU1,1.0D-10)
                DTU = MIN(DTU1,DTU)
            END IF
        END IF
    END IF

!       Reset reference energy and generate new Stumpff coefficients.
    H0(IPAIR) = H(IPAIR)
    30 Z = -0.5D0*H(IPAIR)*DTU**2
    CALL STUMPF(IPAIR,Z)
    Z5 = SF(6,IPAIR)
    Z6 = SF(7,IPAIR)
    DT12 = ONE12*DTU*Z6

!       Convert to physical time units modified by Stumpff coefficients.
    STEP(I1) = (((((TDOT6*DT12 + TDOT5*Z5)*0.2*DTU + 0.5D0*TDOT4)*DTU &
    + TDOT3(IPAIR))*ONE6*DTU + TD2)*DTU + RI)*DTU

!     --09/25/13 14:26-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      IF(NAME(I1).EQ.37) THEN
!$$$         print*,'STEP(I1)',STEP(I1),'DTU',DTU,'H',H(IPAIR),'DT12',DT12,
!$$$     &        ''
!$$$         call flush(6)
!$$$      end if
!$$$      if(step(i1).lt.1/2.0**25.and.step(i1).gt.0) then
!$$$         print*,'I1',I1,'I',I,'N',NAME(I),'STEP',STEP(I1)
!$$$         call flush(6)
!$$$         call abort()
!$$$      end if
!$$$      if(tblock.eq.29.053588867187500.and.name(i1).eq.1759) then
!$$$         write(125+rank,*),'step i1',i1,'s',step(i1),'t0',t0(i1),
!$$$     &        'td6',tdot6,'dt12',dt12,'td5',tdot5,'td2',td2,'ri',ri,
!$$$     *        'z5',z5,'dtu',dtu,'td4',tdot4,'td3',tdot3(ipair),
!$$$     &        't',time
!$$$         call flush(125+rank)
!$$$      end if
!      if(tprev.ge.6.52122497558593750E-002.and.name(i1).eq.499) then
!         print*,rank,'step',i1,step(i1),t0(i1),time,tdot6,dt12,tdot5,
!     *        z5,dtu,tdot4,tdot3(ipair),one6,td2,ri
!         call flush(6)
!      end if
!     --09/25/13 14:26-lwang-end----------------------------------------*
!       Ensure that regularized step is smaller than the c.m. step.
    IF (STEP(I1) > STEP(I) .AND. HI < 0) THEN
        DTU = 0.5D0*DTU
        GO TO 30
    END IF
    DTAU(IPAIR) = DTU

!       See whether the KS slow-down procedure is activated.
    IMOD = KSLOW(IPAIR)
    ZMOD = 0
    IF (IMOD > 1) THEN
        ZMOD = FLOAT(ISLOW(IMOD))
        STEP(I1) = ZMOD*STEP(I1)
    END IF

!       Check diagnostic print option.
    IF (KZ(10) >= 3 .AND. rank == 0) THEN
        WRITE (6,40)  IPAIR, TIME+TOFF, H(IPAIR), RI, DTAU(IPAIR), &
        GI, STEP(I1), LIST(1,I1), IMOD
        40 FORMAT (3X,'KS MOTION',I6,2F10.4,1P,4E10.2,2I4)
    END IF

!       Employ special termination criterion in merger case.
    IF (NAME(I) < 0) THEN
    !       Terminate if apocentre perturbation > 0.25 (R > SEMI) or GI > 0.25.
        IF (HI < 0.0) THEN
            SEMI = -0.5*BODY(I)/HI
        !              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
        !              A0 = SEMI*(1.0 + SQRT(ECC2))/RI
        !       Replace eccentricity calculation with typical value.
            A0 = 1.5*SEMI/RI
            GA = GI*A0*A0*A0
            IF (GA > 0.25 .AND. RI > SEMI) IQ = .TRUE. 
        !       Terminate wide orbits if perturber number is relatively large.
            IF (RI > 20*RMIN .AND. NNB0 > 0.8*LIST(1,I)) IQ = .TRUE. 
        !              IF (GI.GT.0.1.AND.RI.GT.RMIN) IQ = .TRUE.
            IF (GI > 0.25) IQ = .TRUE. 
            IF (MIN(BODY(I1),BODY(I2)) < 0.05*BODYM) THEN
                IF (GI > 2.0D-04) IQ = .TRUE. 
            END IF
        ELSE
            IF (TD2 > 0.0 .AND. (GI > GMAX .OR. RI > RMIN)) IQ = .TRUE. 
        END IF
        IF ( .NOT. IQ) GO TO 60
    END IF

!     For massive black hole ks, try to avoid terminated with large seperation
    IF(KZ(24) == 1 .AND. BODY(I)/BODYM > 200 .AND. RI > 10.0*RMIN) THEN
        IQ = .true.
    END IF

!       Delay termination until end of block for large perturbations.
    IF (IQ) THEN
        DTR = TBLOCK - TIME
    !     --03/26/14 13:58-lwang-debug--------------------------------------*
    !**** Note:------------------------------------------------------------**
    !$$$         if(rank.eq.0)
    !$$$     &    WRITE (6,45)  IPAIR, TTOT, GI, RI, DTR, STEP(I1),TPREV,TBLOCK
    !$$$     *         ,TIME,T0(I1)
    !$$$ 45      FORMAT (' TERM TEST    KS'I4,' T',F10.4,' G',F7.3,' R',E10.2,
    !$$$     &        ' DTR',E15.7,' DT',E15.7,' TP',E15.7,' TB',E15.7,
    !$$$     &        ' T',E15.7,' T0',E15.7)
    !     --03/26/14 13:58-lwang-end----------------------------------------*
        IF (DTR < STEP(I1)) GO TO 90
    END IF

!       Check standard termination criterion (suppress on IQ = .true.).
!     IF (RI.GT.R0(IPAIR)) THEN
    IF (RI > R0(IPAIR) .AND. RI > 2.0*RMIN .AND. .NOT. IQ) THEN
    !       See whether termination can be delayed for sufficient perturbers.
        IF (NNB0 < 0.80*LIST(1,I) .AND. GI < 0.1) GO TO 60
    !     --11/20/13 16:26-lwang-debug--------------------------------------*
    !**** Note:------------------------------------------------------------**
    !$$$         if(rank.eq.0.and.name(i1).eq.230) then
    !$$$            print*,'230 ',RI,R0(IPAIR),RMIN,IQ
    !$$$         end if
    !     --11/20/13 16:26-lwang-end----------------------------------------*
    
    !       See whether termination can be delayed for intermediate energies.
    !$$$          IF (RI.GT.RMIN) THEN
    !$$$              A3 = RMIN*ABS(HI)*BODYIN
    !$$$              GIMAX = A3*A3*A3
    !$$$              IF (GI.LT.GIMAX) GO TO 60
    !       Check updating of R0 for newly hardened binary orbit.
        IF (HI < -ECLOSE) THEN
            SEMI = -0.5*BODY(I)/HI
            R0(IPAIR) = MAX(RMIN,2.0D0*SEMI)
            GO TO 70
        END IF
    !$$$          END IF
    
    !         IF (KZ(4).GT.0.AND.GI.GT.GPRINT(1)) THEN
    !             DGAM = GI - G0
    !             K = KZ(4)
    !             DO 55 L = 2,K
    !                 IF (GI.LT.GPRINT(L)) THEN
    !                     CALL EVOLVE(IPAIR,L-1)
    !                     GO TO 90
    !                 END IF
    !  55         CONTINUE
    !             CALL EVOLVE(IPAIR,K)
    !         END IF
        GO TO 90
    END IF

!       End integration cycle for hyperbolic motion.
    60 IF (HI >= 0.0D0) THEN
        IF (RDOT*TD2 < 0.0D0) THEN
        !       Determine pericentre for hyperbolic two-body motion.
            SEMI = -0.5D0*BODY(I)/HI
            ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
            QPERI = SEMI*(1.0D0 - SQRT(ECC2))
            DMIN2 = MIN(DMIN2,QPERI)
        
        !       Check optional tidal interaction or stellar collision.
            IF (KZ(19) >= 3 .AND. KSTAR(I) < 10) THEN
                VINF = SQRT(2.0*HI)*VSTAR
                KS1 = KSTAR(I1)
                KS2 = KSTAR(I2)
                RX = MAX(RADIUS(I1),RADIUS(I2))
            !       Determine maximum periastron factor for capture (VINF in km/sec).
                IF (KZ(27) <= 2) THEN
                    RFAC = RPMAX2(RADIUS(I1),RADIUS(I2),BODY(I1), &
                    BODY(I2),KS1,KS2,VINF)
                    RCAP = RFAC*RX
                ELSE
                    DV = SQRT(2.0*HI)
                !       Note that Quinlan & Shapiro function returns actual distance.
                    RCAP = RPMAX(BODY(I1),BODY(I2),VSTAR,DV,QPERI)
                END IF
            !                  IF (QPERI.LT.5.0*RX) THEN
            !                      if(rank.eq.0)
            !     &                WRITE (54,54)  TTOT, NAME(I1), NAME(I2), KS1,
            !     &                               KS2, VINF, RCAP*SU, RX*SU, QPERI*SU
            !   54                 FORMAT ('CLOSE   Time[NB] NAME(I1) NAME(I2) ',
            !     &                     'K*(I1) K*(I2) VINF[km/s] RCAP[R*] RX[R*] ',
            !     &                     'PERI[R*] ',1P,E26.17,0P,2I2,2I4,4F12.5)
            !                  END IF
                IF (QPERI < RCAP) THEN
                    J1 = I1
                    IF (RADIUS(I2) > RADIUS(I1)) J1 = I2
                    FAC = 0.5*BODY(I)/BODY(J1)
                !       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                    IF (KZ(27) <= 2) THEN
                        RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                    ELSE
                        RCOLL = 6.0*BODY(I)/CLIGHT**2
                    END IF
                    RI2 = 0.0
                    VI2 = 0.0
                    DO 61 K = 1,3
                        RI2 = RI2 + (X(K,I) - CMR(K))**2
                        VI2 = VI2 + XDOT(K,I)**2
                    61 END DO
                !                      if(rank.eq.0)then
                !                      WRITE (55,58)  TTOT, IPAIR, NAME(I1), NAME(I2),
                !     &                        KS1, KS2, KSTAR(I), VINF,
                !     &                        SQRT(ECC2), HI, R(IPAIR), SEMI,
                !     &                        QPERI, BODY(I1), BODY(I2),
                !     &                        BODY(I)*ZMBAR,
                !     &                        SQRT(RI2)/RC, SQRT(VI2)*VSTAR,
                !     &                        RHOD, RADIUS(I1)*SU, RADIUS(I2)*SU,
                !     &                        RCAP, RADIUS(J1)/QPERI, RCOLL/QPERI
                !   58                 FORMAT ('RPMAX:  Time[NB] IPAIR NAME(I1) ',
                !     &                     'NAME(I2) K*(I1) K*(I2) K*(ICM) VINF[km/s] ',
                !     &                     'ECC H[NB] R12[NB] SEMI[NB] PERI[NB] ',
                !     &                     'M(I1)[NB] M(I2)[NB] M(ICM)[M*] RI[RC] ',
                !     &                     'VI[km/s] RHOD RS(I1)[R*] RS(I2)[R*] ',
                !     &                     'RCAP[NB] RX/PERI RCOLL/PERI',1P,E25.17,0P,
                !     &                     3I12,3I4,2F12.5,1P,7E16.7,0P,8F12.6)
                !                      end if
                !                      CALL FLUSH(55)
                    IF (QPERI < RCOLL) THEN
                    !$$$*       Obtain KS variables at pericentre before merging into one body.
                    !$$$                          CALL KSPERI(IPAIR)
                    !$$$                          KSPAIR = IPAIR
                    !$$$                          IQCOLL = -2
                    !$$$                          CALL CMBODY(2)
                        IPHASE = -1
                    ELSE IF (KSTAR(I) >= 0 .AND. KZ(27) > 0) THEN
                        CALL KSTIDE(IPAIR,QPERI)
                    END IF
                END IF
            !       Check options for artificial collisions.
            ELSE IF (KZ(27) == -1) THEN
                RFAC = 2.0
                IF (QPERI < RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
                    J1 = I1
                    IF (RADIUS(I2) > RADIUS(I1)) J1 = I2
                    FAC = 0.5*BODY(I)/BODY(J1)
                !       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                    RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                    IF (QPERI < RCOLL) THEN
                        CALL TOUCH(IPAIR,I1,I2,RCOLL)
                    END IF
                END IF
            END IF
        END IF
        GO TO 100
    END IF

!       Check perturbation threshold (H < 0 & GAMMA > GMAX).
!     IF (KZ(4).EQ.0.OR.G0.LT.GMAX) GO TO 70

!     K = KZ(4)
!     DO 65 L = 1,K
!         IF ((G0 - GPRINT(L))*(GI - GPRINT(L)).LE.0.0) THEN
!             IF (L.EQ.1) THEN
!                 W1 = -0.5*BODY(I)/HI
!                 W2 = W1*BODYIN
!                 TK = TWOPI*W1*SQRT(W2)
!             END IF

!       Estimate smallest permitted output interval at new level.
!             DTCRIT = TK*ORBITS(L)
!             IF (TIME - TLASTB(L).LT.DTCRIT) GO TO 70
!             DGAM = GI - G0
!             CALL EVOLVE(IPAIR,L)
!             GO TO 70
!         END IF
!  65 CONTINUE

!       Check for partial reflection during approach (NB! only IMOD = 1).
!  70 IF (GI.LT.GMIN.AND.TD2.LT.0.0D0) THEN
!       Skip apocentre position itself.
!         IF (RDOT.LT.0.0D0.AND.IMOD.EQ.1) THEN
!             IF (KZ(25).GT.0) CALL FREEZE(IPAIR)
!             GO TO 100
!         END IF
!     END IF

!       Determine new perturbers for binary at apocentre turning point.
    70 IF (RDOT*TD2 >= 0.0D0) GO TO 100
    SEMI = -0.5D0*BODY(I)/HI

!       Check minimum two-body separation just after pericentre.
    IF (RDOT < 0.0D0) THEN
    !       Obtain pericentre by Mikkola's algorithm (GAMMA < 0.001).
        IF (GI < 0.001) THEN
            CALL PERI(UI,UIDOT,RI,BODY(I1),BODY(I2),QPERI)
        ELSE
            QPERI = RI
        END IF
        DMIN2 = MIN(DMIN2,QPERI)
    
    !       Check optional tidal interaction or stellar collision (skip merger).
        IF (KZ(19) >= 3 .AND. KSTAR(I) <= 10 .AND. NAME(I) > 0) THEN
            RFAC = 5.0
            IF (KZ(27) <= 2) THEN
                IF (KZ(27) == 1) RFAC = 4.0
                RX = RFAC*MAX(RADIUS(I1),RADIUS(I2))
            ELSE
                RX = RPMIN(BODY(I1),BODY(I2),VSTAR,HI,QPERI)
            END IF
            IF (QPERI < RX) THEN
                J1 = I1
                IF (RADIUS(I2) > RADIUS(I1)) J1 = I2
                FAC = 0.5*BODY(I)/BODY(J1)
            !       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                IF (KZ(27) <= 2) THEN
                    RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                ELSE
                    RCOLL = 6.0*BODY(I)/CLIGHT**2
                END IF
                IF (QPERI < RCOLL) THEN
                !$$$*       Obtain KS variables at pericentre before merging into one body.
                !$$$                      CALL KSPERI(IPAIR)
                !$$$                      KSPAIR = IPAIR
                !$$$                      IQCOLL = -2
                !$$$                      CALL CMBODY(2)
                    IPHASE = -1
                ELSE IF (KSTAR(I) >= 0) THEN
                !       Distinguish between sequential, standard and GR circularization.
                    IF (KZ(27) == 1) THEN
                        ICIRC = 1
                        TC = 0.0
                    ELSE IF (KZ(27) == 2) THEN
                        ECC2 = (1.0 - RI/SEMI)**2 + &
                        TDOT2(IPAIR)**2/(BODY(I)*SEMI)
                        ECC = SQRT(ECC2)
                        ICIRC = 0
                        CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
                    ELSE
                        ICIRC = 1
                        TC = 0.0
                    END IF
                    IF (KSTAR(I) >= 10) ICIRC = 0
                !       Skip tidal effects for circularization time above 100 Myr (07/08).
                    IF (ICIRC > 0 .AND. KZ(27) > 0 .AND. &
                    TC < 100.0) THEN
                        CALL KSTIDE(IPAIR,QPERI)
                    END IF
                END IF
            END IF
        !       Check for perturbed spiral or chaos case (skip collision).
            IF (KSTAR(I) == -2 .AND. IPHASE == 0) THEN
                CALL SPIRAL(IPAIR)
            ELSE IF (KSTAR(I) == -1 .AND. IPHASE == 0) THEN
                CALL KSTIDE(IPAIR,QPERI)
            END IF
        !       Check options for artificial collisions.
        ELSE IF (KZ(27) == -1) THEN
            RFAC = 2.0
            IF (QPERI < RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
                J1 = I1
                IF (RADIUS(I2) > RADIUS(I1)) J1 = I2
                FAC = 0.5*BODY(I)/BODY(J1)
            !       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                IF (QPERI < RCOLL) THEN
                    CALL TOUCH(IPAIR,I1,I2,RCOLL)
                END IF
            END IF
        END IF
        GO TO 100
    END IF

!       Save maximum separation of persistent binary.
    RMAX = MAX(RMAX,RI)
!     ks MPI communicate rmax
!      call ksparmpi(K_store,K_real8,K_RMAX,0,0,RMAX)
!     --03/05/14 20:34-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      print*,rank,'rmax',rmax
!     --03/05/14 20:34-lwang-end----------------------------------------*

!       Check binary reference radius or merger stability criterion.
    IF (NAME(I) > 0) THEN
    !       Update termination length scale in case of initial soft binary.
        EB = BODY(I1)*BODY(I2)*HI*BODYIN
        IF (EB < EBH) R0(IPAIR) = MAX(RMIN,2.0*SEMI)
    ELSE
        ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
        ECC = SQRT(ECC2)
        RP = SEMI*(1.0 - ECC)*(1.0 - 2.0*GI)
    !       Find merger index.
        IM = 0
        DO 72 K = 1,NMERGE
            IF (NAMEM(K) == NAME(I)) IM = K
        72 END DO
    !       Exclude inner planets from the general stability test.
        IF (MIN(BODY(I1),BODY(I2)) < 0.05*BODYM) THEN
            IF (RP < R0(IPAIR)) GO TO 90
        END IF
    !       Assess the stability inside critical pericentre (safety factor 1.04).
        IF (RP < 1.04*R0(IPAIR)) THEN
            EOUT = ECC
        !       Increase tolerance near sensitive stability boundary (RM 10/2008).
            IF (EOUT > 0.90) THEN
                DE = 0.5*(1.0 - EOUT)
                DE = MIN(DE,0.01D0)
            !       Add extra amount 0.011 to avoid switching.
                IF (ECC > 0.9) DE = DE + 0.011
                EOUT = EOUT - DE
            END IF
        !       Note: assessment needs to use same eccentricity as for acceptance.
            CALL ASSESS(IPAIR,IM,EOUT,SEMI,ITERM)
            IF (ITERM > 0) THEN
                INSTAB = INSTAB + 1
                GO TO 90
            END IF
        END IF
    !       Check possible eccentricity modulation or t_circ update.
        IF (IM > 0 .AND. (TIME > TMDIS(IM) .OR. &
        TMDIS(IM) > 1.0D+06)) THEN
            IF (KZ(27) == 2) THEN
                CALL ECCMOD(I,ITERM)
                IF (ITERM > 0) THEN
                !                     if(rank.eq.0)
                !    &                WRITE (6,76)  RP, R0(IPAIR)
                !  76                 FORMAT (' ECCMOD TERM    RP R0 ',1P,2E10.2)
                    GO TO 90
                END IF
            !       Consider both inner and possible outer circularizing binary.
                DO 78 K = 1,NCHAOS
                    IF (NAMEC(K) == NZERO - NAMEM(IM) .AND. &
                    KSTARM(IM) == -2) THEN
                    !       Update unperturbed binary if T - TOSC > 10 Myr (cf. IMPACT & DECIDE).
                        IF ((TIME - TOSC(K))*TSTAR > 10.0) GO TO 90
                    END IF
                    IF (NAMEC(K) == NAMEG(IM) .AND. &
                    KSTARM(IM) == -2) THEN
                        IF ((TIME - TOSC(K))*TSTAR > 10.0) GO TO 90
                    END IF
                !       Note: perturbed binary is treated if pericentre before next IMPACT.
                78 END DO
            END IF
        END IF
    END IF

!       Produce diagnostics for any circularizing perturbed binary.
    IF (KSTAR(I) == -2 .AND. GI > 0.01) THEN
        ECC = RI/SEMI - 1.0
        QPS = SEMI*(1.0 - ECC)/MAX(RADIUS(I1),RADIUS(I2))
        ZM = BODY(I1)/BODY(I2)
    !          if(rank.eq.0)
    !     &    WRITE (21,80)  TTOT, NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
    !     &                   LIST(1,I1), QPS, ZM, GI, ECC, SEMI
    !   80     FORMAT (' PERT SPIRAL    T NAM K* NP QP/S M1/M2 G E A ',
    !     &                             F11.4,2I6,3I3,2F5.1,2F8.4,1P,E10.2)
    !          CALL FLUSH(21)
    END IF

!       See whether KS slow-down procedure should be (re)-checked (no Chaos).
    IF (KZ(26) > 0 .AND. KSTAR(I) >= 0) THEN
        KMOD = INT(RANGE*GMIN/MAX(GI,1.0D-10))
        IF (KMOD > 1 .OR. IMOD > 1) THEN
            CALL KSMOD(IPAIR,KMOD)
            IF (KMOD < 0) GO TO 100
        END IF
    END IF

!       Set approximate value of next period with perturbation included.
    TK = TWOPI*SEMI*SQRT(SEMI*BODYIN)*(1.0 + GI)
    IF (IMOD > 1) THEN
        TK = ZMOD*TK
    END IF

!       Use old perturber list if next apocentre is before the c.m. step.
    IF (TIME + TK < T0(I) + STEP(I)) THEN
        GO TO 100
    END IF

!       Select new perturbers (J = N adopted for unperturbed Chaos).
    CALL KSLIST(IPAIR)

!       Check rectification of chaotic spiral at start of unperturbed motion.
    IF (KSTAR(I) == -2 .AND. LIST(1,I1) == 0) THEN
        DMR = 0.D0
        CALL CHRECT(IPAIR,DMR)
        IF (IPHASE < 0) GO TO 100
    ELSE
        CALL KSRECT(IPAIR)
    END IF

!       Check optional search criterion for multiple encounter or merger.
    IF (KZ(15) > 0 .AND. STEP(I) < DTMIN) THEN
        CALL IMPACT(I)
    END IF
    GO TO 100

!       Terminate regularization of current pair (IPAIR set in KSPAIR).
    90 KSPAIR = IPAIR
!     --03/10/14 22:02-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      print*,rank,'iq',iq,'t',time,'dtr',dtr,'step(i1)',step(i1),'iterm'
!$$$     *     ,iterm,'n',name(i),'GI',GI,'HI',hi,'JCOMP',jcomp,'N',N
!     --03/10/14 22:02-lwang-end----------------------------------------*
!       Set indicator for calling KSTERM in MAIN (permits phase overlay).
    IPHASE = 2
!       Check case of hierarchical binary.
    IF (NAME(I) < 0) IPHASE = 7

!     --09/25/13 14:22-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$ 100  if(tprev.ge.6.52122497558593750E-002.and.name(i1).eq.499) then
!$$$        print*,rank,'ksinte',i1,step(i1),t0(i1),time,tprev
!$$$        call flush(6)
!$$$      end if
!      print*,rank,'left subint',tblock
!      call flush(6)
!$$$      if (tprev.ge.6.80656433105468750E-002) stop
!     --09/25/13 14:23-lwang-end----------------------------------------*
          
    100 RETURN

    END SUBROUTINE KSINT
