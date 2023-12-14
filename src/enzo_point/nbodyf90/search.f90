    SUBROUTINE SEARCH(I,IKS)


!       Close encounter search.
!       -----------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5), &
    NAMES(NCMAX,5),ISYS(5)
    REAL*8 :: RJMIN2
!     Trace flag for debug
    INTEGER :: TRACE_FLAG
!     1. JCOMP<IFIRST || JCOMP>N
!     2. RJMIN2 > RMIN2
!     3. RDOT > 0.02*sqrt(BODY(I)+BODY(JCOMP)*RIJMIN)
!     4. GI>0.25
!     5. sub: NAMEI==NAME(I) || NAMEI==NAME(JCOMP)
!     6. ch: NAME(I)==0 || NAME(JCOMP)==0

    RJMIN2 = RMIN2 + RMIN2

!       Increase counter for regularization attempts.
    NKSTRY = NKSTRY + 1

!       Predict body #I to order FDOT for parallel prediction (not active RSp).
    call jpred(i,time,time)
!$$$      if(TPRED(I).ne.TIME) then
!$$$*         print*,rank,'Nopred in search',I
!$$$         S = TIME - T0(I)
!$$$         DO 1 K = 1,3
!$$$             X(K,I) = ((FDOT(K,I)*S + F(K,I))*S + X0DOT(K,I))*S +
!$$$     &                                                           X0(K,I)
!$$$             XDOT(K,I) = (3.0*FDOT(K,I)*S + 2.0*F(K,I))*S + X0DOT(K,I)
!$$$   1     CONTINUE
!$$$      end if

    FMAX = 0.0
    NCLOSE = 0
!       Find dominant neighbour by selecting all STEP(J) < 2*DTMIN.
    L = LIST(1,I) + 2
    2 L = L - 1
    IF (L < 2) GO TO 10

!       First see whether any c.m. with small step is within 2*RMIN.
    IF (LIST(L,I) <= N) GO TO 4

    J = LIST(L,I)
!      IF (STEP(J).GT.SMIN) GO TO 2
    IF (STEP(J) > 4*STEP(I) .AND. BODY(J) < 10*BODY(I)) GO TO 2
!       Prediction body #J to order FDOT for parallel prediction (not active RSp).
    call jpred(j,time,time)
!$$$      if(TPRED(J).ne.TIME) then
!$$$*         print*,rank,'Nopred in searchj',J
!$$$         S = TIME - T0(J)
!$$$         DO 225 K = 1,3
!$$$             X(K,J) = ((FDOT(K,J)*S + F(K,J))*S + X0DOT(K,J))*S +
!$$$     &                                                           X0(K,J)
!$$$             XDOT(K,J) = (3.0*FDOT(K,J)*S + 2.0*F(K,J))*S + X0DOT(K,J)
!$$$  225     CONTINUE
!$$$       end if
    A1 = X(1,J) - X(1,I)
    A2 = X(2,J) - X(2,I)
    A3 = X(3,J) - X(3,I)
    RIJ2 = A1*A1 + A2*A2 + A3*A3
    IF (RIJ2 > 2*RMIN22) GO TO 2

    FIJ = BODY(J)/RIJ2
    IF (FMAX < FIJ) FMAX = FIJ
!       Abandon further search if c.m. force exceeds half the total force.
    IF (FMAX**2 < F(1,I)**2 + F(2,I)**2 + F(3,I)**2) THEN
        NCLOSE = NCLOSE + 1
        JLIST(NCLOSE) = J
        GO TO 2
    ELSE
        GO TO 10
    END IF

!       Continue searching single particles with current value of FMAX.
    4 JCOMP = 0

    DO 6 K = L,2,-1
        J = LIST(K,I)
    !          IF (STEP(J).GT.SMIN) GO TO 6
        IF (STEP(J) > 8*STEP(I)) GO TO 6
    !       Prediction body #J to order FDOT for parallel prediction (not active RSp).
        call jpred(j,time,time)
    !$$$          if(TPRED(J).ne.TIME) then
    !$$$*             print*,rank,'Nopred in searchj2',j
    !$$$         S = TIME - T0(J)
    !$$$         DO 65 KK = 1,3
    !$$$         X(KK,J) = ((FDOT(KK,J)*S + F(KK,J))*S + X0DOT(KK,J))*S +
    !$$$     &                                                         X0(KK,J)
    !$$$         XDOT(KK,J) = (3.0*FDOT(KK,J)*S + 2.0*F(KK,J))*S + X0DOT(KK,J)
    !$$$  65     CONTINUE
    !$$$          end if
        A1 = X(1,J) - X(1,I)
        A2 = X(2,J) - X(2,I)
        A3 = X(3,J) - X(3,I)
        RIJ2 = A1*A1 + A2*A2 + A3*A3
        IF (RIJ2 < 2*RMIN22) THEN
            NCLOSE = NCLOSE + 1
            JLIST(NCLOSE) = J
        !       Remember index of every single body with small step inside 2*RMIN.
            FIJ = (BODY(I) + BODY(J))/RIJ2
            IF (FIJ > FMAX) THEN
                FMAX = FIJ
            !       Save square distance and global index of dominant body.
                RJMIN2 = RIJ2
                JCOMP = J
            END IF
        END IF
    6 END DO

!       See whether dominant component is a single particle inside RMIN.
    IF (JCOMP < IFIRST .OR. JCOMP > N) THEN
        TRACE_FLAG = 1
        GO TO 10
    END IF
    IF (RJMIN2 > RMIN2) THEN
        TRACE_FLAG = 2
        GO TO 10
    END IF

    RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) + &
    (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) + &
    (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))

!       Only select approaching particles (include nearly circular case).
    RIJMIN = SQRT(RJMIN2)
    IF (RDOT > 0.02*SQRT((BODY(I) + BODY(JCOMP))*RIJMIN)) THEN
        TRACE_FLAG = 3
        GO TO 10
    END IF

!       Ensure a massive neighbour is included in perturbation estimate.
    BCM = BODY(I) + BODY(JCOMP)
    IF (BODY1 > 10.0*BCM) THEN
        JBIG = 0
        BIG = BCM
        NNB1 = LIST(1,I) + 1
        DO 20 L = 2,NNB1
            J = LIST(L,I)
            IF (BODY(J) > BIG) THEN
                JBIG = J
                BIG = BODY(J)
            END IF
        20 END DO
    !       Check whether already present, otherwise add to JLIST.
        DO 25 L = 1,NCLOSE
            IF (JLIST(L) == JBIG) THEN
                JBIG = 0
            END IF
        25 END DO
        IF (JBIG > 0) THEN
            NCLOSE = NCLOSE + 1
            JLIST(NCLOSE) = JBIG
        END IF
    END IF

!       Evaluate vectorial perturbation due to the close bodies.
    CALL FPERT(I,JCOMP,NCLOSE,PERT)

!       Accept #I & JCOMP if the relative motion is dominant (GI < 0.25).
    GI = PERT*RJMIN2/BCM
    IF (GI > 0.25) THEN
    !         IF (KZ(4).GT.0.AND.TIME-TLASTT.GT.4.44*TCR/FLOAT(N))
    !    &                                             CALL EVOLVE(JCOMP,0)
        TRACE_FLAG=4
        GO TO 10
    END IF

!       Exclude any c.m. body of compact subsystem (TRIPLE, QUAD or CHAIN).
    DO 8 ISUB = 1,NSUB
        NAMEI = NAMES(1,ISUB)
        IF (NAMEI == NAME(I) .OR. NAMEI == NAME(JCOMP)) THEN
            TRACE_FLAG=5
            GO TO 10
        END IF
    8 END DO

!       Also check possible c.m. body of chain regularization (NAME = 0).
    IF (NCH > 0) THEN
        IF (NAME(I) == 0 .OR. NAME(JCOMP) == 0) THEN
            TRACE_FLAG=6
            GO TO 10
        END IF
    END IF

!       Save index and increase indicator to denote new regularization.
    ICOMP = I
    IKS = IKS + 1

    10 CONTINUE

!     --07/19/14 22:42-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!      print*,'SEARCH I',I,'NB',LIST(1,I),'L',LIST(2:LIST(1,I),I),
!     &     'TRACE_F',TRACE_FLAG,'IKS',IKS,'JCOMP',JCOMP
!     --07/19/14 22:42-lwang-end----------------------------------------*

          
    RETURN

    END SUBROUTINE SEARCH
