    SUBROUTINE FINDM(I,ITERM,MG)


!       Find ghost mass.
!       ----------------

    USE POINTERS
    INCLUDE 'common6.h'
          
    COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX), &
    HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX), &
    NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
    REAL*8 :: MG


!       Distinguish between KS component and single particle or c.m.
    ITERM = 0
    IF (I < IFIRST) THEN
        IPAIR = KVEC(I)
        ICM = N + IPAIR
    !       Consider standard binary, simple hierarchy or high-order system.
        IF (NAME(ICM) > 0) THEN
            IM = 0
        !       Obtain merger index for later and identify component of CM.
            DO 3 K = 1,NMERGE
                IF (NAMEG(K) == NAME(ICM)) IM = K
            3 END DO
            IF (IM == 0) THEN
                ITERM = -1
                GO TO 50
            END IF
            J1 = 2*IPAIR - 1
            IF (NAME(J1) == NAME(I)) K1 = 3
            IF (NAME(J1+1) == NAME(I)) K1 = 4
        ELSE IF (NAME(ICM) >= -2*NZERO) THEN
            IM = 0
            DO 5 K = 1,NMERGE
                IF (NAMEM(K) == NAME(ICM)) IM = K
            5 END DO
            IF (IM == 0) THEN
                ITERM = -1
                GO TO 50
            END IF
            J1 = 2*IPAIR - 1
            K1 = 1
            IF (NAME(J1+1) == NAME(I)) THEN
                K1 = 2
            END IF
        ELSE IF (NAME(ICM) < -2*NZERO) THEN
            IM = 0
            DO 6 K = 1,NMERGE
                IF (NAMEM(K) == NAME(ICM)) IM = K
            6 END DO
            IF (IM == 0) THEN
                ITERM = -1
                GO TO 50
            END IF
            JH = 0
        !       Search c.m. ghost name to get KS pair index.
            DO 8 J = N+1,NTOT
                IF (NAME(J) == NAMEG(IM)) JH = J
            8 END DO
            IF (JH == 0) THEN
                ITERM = -1
                GO TO 50
            END IF
            JPAIR = JH - N
            J1 = 2*JPAIR - 1
        !       Compare component names in order to decide appropriate CM index.
            IF (NAME(J1) == NAME(I)) THEN
                K1 = 1
            ELSE IF (NAME(J1+1) == NAME(I)) THEN
                K1 = 2
            ELSE
                K1 = 2
            END IF
        ELSE
            ITERM = -1
            GO TO 50
        END IF
    ELSE
    !       Determine merger index of ghost particle.
        IM = 0
        DO 10 K = 1,NMERGE
            IF (NAMEG(K) == NAME(I)) IM = K
        10 END DO
        IF (IM == 0) THEN
            ITERM = -1
            GO TO 50
        END IF
        IF (I <= N) THEN
            J1 = 0
        !       Identify the location of the corresponding KS component.
            DO 15 J = 1,IFIRST
                IF (NAME(J) == -NAMEM(IM)) J1 = J
            15 END DO
            IF (J1 == 0) THEN
                ITERM = -1
                GO TO 50
            END IF
            IPAIR = KVEC(J1)
        !       Decide the relevant component by comparing look-up times.
            K1 = 1
            IF (TEV(2*IPAIR) < TEV(I)) K1 = 2
        ELSE
            ICM = 0
            DO 20 K = N+1,NTOT
                IF (NAMEM(IM) == NAME(K)) ICM = K
            20 END DO
            IF (ICM == 0) THEN
                ITERM = -1
                GO TO 50
            END IF
            IPAIR = ICM - N
            J = I
            K1 = 1
            IF (TEV(2*IPAIR) < TEV(J)) K1 = 2
        END IF
    END IF

!       Copy the desired mass from merger table.
    MG = CM(K1,IM)

    50 RETURN

    END SUBROUTINE FINDM
