    subroutine xbpredall


!     Predict x and xdot. (L.WANG)

    USE POINTERS
    INCLUDE 'common6.h'
          
    INCLUDE 'omp_lib.h'
    COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
    REAL*8 :: TPRED
    LOGICAL :: iPREDALL

    IF (IPREDALL) RETURN

    NNPRED = NNPRED + 1
! omp parallel do private(J,S,S1,S2,JPAIR,J1,J2,ZZ)
    DO 40 J = IFIRST,NTOT
    !     IF(TPRED(J).NE.TIME) THEN
        S = TIME - T0(J)
        S1 = 1.5*S
        S2 = 2.0*S
        X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
        X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
        X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
        XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
        XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
        XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
        TPRED(J) = TIME
        IF (J > N) THEN
            JPAIR = J - N
            IF (LIST(1,2*JPAIR - 1) > 0) THEN
                ZZ = 1.0
                IF (GAMMA(JPAIR) > 1.0D-04) ZZ = 0.0
                CALL KSRES2(JPAIR,J1,J2,ZZ,TIME)
            END IF
        !     --03/07/14 21:22-lwang-debug--------------------------------------*
        !**** Note:------------------------------------------------------------**
        !$$$            IF(J.EQ.12195) THEN
        !$$$               print*,rank,'PA I',J,'N',NAME(J),'X',X(1,J),'X1',X(1,J1),
        !$$$     &              'X2',X(1,J2),'T',TIME,'TPR',TPRED(J),ipredall,
        !$$$     &              'T0',T0(J),'x0',x0(1,j),'F',F(1,j),'FD',FDOT(1,j)
        !$$$               call flush(6)
        !$$$            end if
        !     --03/07/14 21:22-lwang-end----------------------------------------*
        END IF
    40 END DO
! omp end parallel do
!     --05/12/14 17:26-lwang-debug--------------------------------------*
!**** Note:------------------------------------------------------------**
!$$$      do KK=ifirst,ntot
!$$$         write(102+rank,*) 'I',KK,'N',NAME(KK),'X',X(1:3,KK),
!$$$     &        'XD',XDOT(1:3,KK),'X0',X0(1:3,KK),'X0D',X0DOT(1:3,kk),
!$$$     &        'T0',T0(KK)
!$$$         if(KK.GT.N) then
!$$$            jpair=kk-n
!$$$            k1=2*(kk-n)-1
!$$$            k2=k1+1
!$$$            if(list(1,k1).gt.0) then
!$$$            write(102+rank,*) 'P',K1,K2,'N',NAME(K1),NAME(K2),'X1',
!$$$     &           X(1:3,K1),'X2',X(1:3,K2),'XD1',XDOT(1:3,K1),
!$$$     &              'XD2',XDOT(1:3,K2),'U0',U0(1:4,KK-n),'UDOT',
!$$$     &              UDOT(1:4,jpair),'FU',FU(1:4,jpair),'FUDOT',
!$$$     &              FUDOT(1:4,jpair),'FUDOT2',FUDOT2(1:4,jpair),
!$$$     &              'T',TIME
!$$$            end if
!$$$         end if
!$$$      end do
!     --05/12/14 17:26-lwang-end----------------------------------------*
    iPREDALL = .true.

    return

    end subroutine xbpredall
          
