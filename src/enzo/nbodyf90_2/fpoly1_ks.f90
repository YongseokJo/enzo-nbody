    SUBROUTINE FPOLY1_KS(I1,I2,KCASE)


!       Force & first derivative.
!       -------------------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    REAL*8 ::  A(9),F1(3),F1DOT(3)
    REAL*8 ::  XBACK(3,LMAX),XDBACK(3,LMAX)
!     --03/05/14 15:20-lwang-debug--------------------------------------**
!**** Note:------------------------------------------------------------**
!$$$      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),predall
!$$$      REAL*8 TPRED
!$$$      LOGICAL predall
!     --03/05/14 15:20-lwang-end----------------------------------------*

!       Standard case, new c.m. or KS termination (KCASE = 0, 1, 2).
    JLAST = NTOT
!       Reduce loop size for new c.m. polynomial.
    IF (KCASE == 1) JLAST = NTOT - 1

!       Loop over all bodies, pair #ICOMP & JCOMP or one single body.
    DO 40 I = I1,I2
    
    !       Initialize forces & first differences for body #I.
        DO 10 K = 1,3
            FI(K,I) = 0.0D0
            FR(K,I) = 0.0D0
            D1(K,I) = 0.0D0
            D1R(K,I) = 0.0D0
        10 END DO
    
    !       Obtain force & first derivative by summing over all bodies.
        KDUM = 0
        NNB = LIST(1,I)
    !       Set index of first neighbour to be identified in force loop.
        L = 2
        NAMEJ = LIST(L,I)
    !       Backup all predicted x and xdot to xback and xdback
    !      CALL JPRED(I,TIME,TIME)
        DO II = 1, NNB
            IL = LIST(II+1,I)
        !         call jpred(IL,TIME,TIME)
            XBACK(1,II) = X0(1,IL)
            XBACK(2,II) = X0(2,IL)
            XBACK(3,II) = X0(3,IL)
            X0(1,IL) = X(1,IL)
            X0(2,IL) = X(2,IL)
            X0(3,IL) = X(3,IL)
            XDBACK(1,II) = X0DOT(1,IL)
            XDBACK(2,II) = X0DOT(2,IL)
            XDBACK(3,II) = X0DOT(3,IL)
            X0DOT(1,IL) = XDOT(1,IL)
            X0DOT(2,IL) = XDOT(2,IL)
            X0DOT(3,IL) = XDOT(3,IL)
        !     --09/25/13 19:20-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$      if(i.eq.14829) then
        !$$$         write(122+rank,*) 'I',I,'IL',IL,'X',XBACK(1,II),
        !$$$     &        'XD',XDBACK(1,II),'t',time,'NNB',nnb
        !$$$         call flush(122+rank)
        !$$$      end if
        !     --09/25/13 19:20-lwang-end----------------------------------------*
        END DO
    
    !       Sum over all other bodies.
        DO 30 JDUM = IFIRST,JLAST
            IF (JDUM == I) GO TO 30
            J = JDUM
            IF (J > N .AND. J == NAMEJ) THEN
                JPAIR = J - N
            !       Use c.m. approximation for unperturbed binary.
                IF (LIST(1,2*JPAIR-1) > 0) THEN
                    KDUM = 2*JPAIR - 1
                    J = KDUM
                END IF
            END IF
        
            12 If(J >= IFIRST) THEN
                DO 15 K = 1,3
                    A(K) = X0(K,J) - X(K,I)
                    A(K+3) = X0DOT(K,J) - XDOT(K,I)
                15 END DO
            ELSE
                DO K = 1,3
                    A(K) = X(K,J) - X(K,I)
                    A(K+3) = XDOT(K,J) - XDOT(K,I)
                END DO
            END IF
        !     --09/25/13 19:20-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$      if(i.eq.14829) then
        !$$$         if(j.ge.ifirst) then
        !$$$            write(120+rank,*) 'j',j,'n',name(j),'x0',x0(1,j),
        !$$$     &        'x0d',x0dot(1,j),'f',f(1,j),'fd',fdot(1,j),'t',time
        !$$$c$$$     &           'tpr',tpred(j)
        !$$$         call flush(120+rank)
        !$$$         else
        !$$$            write(120+rank,*) 'j',j,'n',name(j),'x',x(1,j),
        !$$$     &           'xd',xdot(1,j),'t',time
        !$$$         end if
        !$$$      end if
        !     --09/25/13 19:20-lwang-end----------------------------------------*
        
            A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
            A(8) = BODY(J)*A(7)*SQRT(A(7))
            A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
        
            DO 20 K = 1,3
                F1(K) = A(K)*A(8)
                F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
            20 END DO
        
        !       See whether summation index is equal to either component.
            IF (J == ICOMP .OR. J == JCOMP) THEN
                IF (KCASE == 1) GO TO 30
            !       Note that dominant terms cancel analytically in c.m. polynomial.
            END IF
        
            IF (JDUM /= NAMEJ) THEN
                DO 25 K = 1,3
                    FR(K,I) = FR(K,I) + F1(K)
                    D1R(K,I) = D1R(K,I) + F1DOT(K)
                25 END DO
            !     --11/24/13 15:40-lwang-debug--------------------------------------*
            !**** Note:------------------------------------------------------------**
            !$$$              if(rank.eq.0) then
            !$$$                 if(time.ge.3.11303138732910156E-003) then
            !$$$                    if(J.LE.300) then
            !$$$                       print*,'Body ',J,BODY(J)
            !$$$                    end if
            !$$$                 end if
            !$$$              end if
            !     --11/24/13 15:40-lwang-end----------------------------------------*
            ELSE
            !     --11/24/13 15:23-lwang-debug--------------------------------------*
            !**** Note:------------------------------------------------------------**
            !$$$             if(rank.eq.0) then
            !$$$                if(time.ge.3.11303138732910156E-003) then
            !$$$                   print*,'I NAMEJ D1',I,NAMEJ,D1(1,I),J,BODY(J),F1(1),
            !$$$     *                  F1DOT(1),X(1,j),XDOT(1,J)
            !$$$                end if
            !$$$             end if
            !$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)
            !     --11/24/13 15:23-lwang-end----------------------------------------*
                DO 28 K = 1,3
                    FI(K,I) = FI(K,I) + F1(K)
                    D1(K,I) = D1(K,I) + F1DOT(K)
                28 END DO
            
                IF (J == KDUM) THEN
                    J = J + 1
                    GO TO 12
                END IF
            
            !       Advance the neighbour list until last member has been considered.
                IF (L <= NNB) THEN
                    L = L + 1
                    NAMEJ = LIST(L,I)
                END IF
            END IF
        30 END DO
    
    !       Check option for interstellar clouds (force & first derivative).
        IF (KZ(13) > 0 .AND. TIME > 0.0D0) THEN
            CALL FCLOUD(I,F1,F1DOT,2)
        END IF

    !       Revert X0 and X0dot to original values
        DO II = 1, NNB
            IL = LIST(II+1,I)
            X0(1,IL) = XBACK(1,II)
            X0(2,IL) = XBACK(2,II)
            X0(3,IL) = XBACK(3,II)
            X0DOT(1,IL) = XDBACK(1,II)
            X0DOT(2,IL) = XDBACK(2,II)
            X0DOT(3,IL) = XDBACK(3,II)
        END DO
    40 END DO

!       Check option for external force.
    IF (KZ(14) > 0) THEN
        CALL XTRNLD(I1,I2,1)
    END IF

!       Set total force & first derivative.
    DO 50 I = I1,I2
    
        DO 45 K = 1,3
            F(K,I) = FI(K,I) + FR(K,I)
            FDOT(K,I) = D1(K,I) + D1R(K,I)
        45 END DO
    50 END DO

!       Check case of new c.m. (KCASE = 1 with I1 = ICOMP, I2 = JCOMP).
    IF (KCASE == 1) THEN
    !       Form total c.m. force & first derivative and also AC components.
        A1 = BODY(ICOMP)
        A2 = BODY(JCOMP)
        DO 60 K = 1,3
            F(K,NTOT) = (A1*F(K,ICOMP) + A2*F(K,JCOMP))/BODY(NTOT)
            FDOT(K,NTOT) = (A1*FDOT(K,ICOMP) + A2*FDOT(K,JCOMP))/ &
            BODY(NTOT)
            FI(K,NTOT) = (A1*FI(K,ICOMP) + A2*FI(K,JCOMP))/BODY(NTOT)
            D1(K,NTOT) = (A1*D1(K,ICOMP) + A2*D1(K,JCOMP))/BODY(NTOT)
            FR(K,NTOT) = (A1*FR(K,ICOMP) + A2*FR(K,JCOMP))/BODY(NTOT)
            D1R(K,NTOT) = (A1*D1R(K,ICOMP) + A2*D1R(K,JCOMP))/ &
            BODY(NTOT)
        60 END DO
        J1 = NTOT
        J2 = NTOT
    ELSE
        J1 = I1
        J2 = I2
    END IF

    DO 80 I = J1,J2
        DO 70 K = 1,3
            D0(K,I) = FI(K,I)
            D0R(K,I) = FR(K,I)
            FIDOT(K,I) = D1(K,I)
            FRDOT(K,I) = D1R(K,I)
        70 END DO
    80 END DO

    RETURN

    END SUBROUTINE FPOLY1_KS
