      SUBROUTINE ADDINIT(INEW)
*
*
*       initialization of newly created particles
*       ----------------
*
      INCLUDE 'common6.h'

      INTEGER  I,INEW,NNBI
      INTEGER  J,J2,L,K

      REAL*8   FDUMI(3),JRKDUMI(3)

      REAL*8   XI,YI,ZI,VXI,VYI,VZI,MI
      REAL*8   DX,DY,DZ,DVX,DVY,DVZ

      REAL*8   R2,RINV1,RINV2,RI2,RSI
      REAL*8   MRINV1,MRINV3

      REAL*8   DV(3)
      REAL*8   A1,A2,A3

      I = INEW

      RI2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
      write(6,*) 'RI2 is ',RI2
      write(6,*) 'RS0 is',RS0
      RSI = 0.20D0*SQRT(1.0 + RI2)

      write(6,*) 'RSI is',RSI

    1 NNBI = 0

      FR(1:3,I) = 0.0D0
      D1R(1:3,I) = 0.0D0

      write(6,*) "NTOT is",NTOT
      write(6,*) "IFIRST is",IFIRST
      write(6,*) "NPAIRS is",NPAIRS


      DO J = (2*NPAIRS+1),N

        DX = X(1,J) - X(1,I)
        DY = X(2,J) - X(2,I)
        DZ = X(3,J) - X(3,I)

        DVX = XDOT(1,J) - XDOT(1,I)
        DVY = XDOT(2,J) - XDOT(2,I)
        DVZ = XDOT(3,J) - XDOT(3,I)

        R2 = DX*DX + DY*DY + DZ*DZ
        RINV1 = 1.0D0/SQRT(R2)
        RINV2 = 1.0D0/R2

        IF (R2.LE.RSI) THEN

            NNBI = NNBI + 1  
            LIST(1,I) = NNBI
            LIST((NNBI+1),I) = J
            
            RINV1 = 0.0D0
  
        END IF

        MRINV1 = BODY(J)*RINV1
        MRINV3 = MRINV1*RINV2
        RV = -3.0D0 * RINV2

        FDUMI(1) = MRINV3 * DX
        FDUMI(2) = MRINV3 * DY
        FDUMI(3) = MRINV3 * DZ

        JRKDUMI(1) = MRINV3 * (DVX + RV * DX)
        JRKDUMI(2) = MRINV3 * (DVY + RV * DY)
        JRKDUMI(3) = MRINV3 * (DVZ + RV * DZ)

        FR(1,I) = FR(1,I) + FDUMI(1)
        FR(1,I) = FR(2,I) + FDUMI(2)
        FR(1,I) = FR(3,I) + FDUMI(3)

        D1R(1,I) = D1R(1,I) + JRKDUMI(1)
        D1R(1,I) = D1R(2,I) + JRKDUMI(2)
        D1R(1,I) = D1R(3,I) + JRKDUMI(3)

      END DO

*     make sure that NNB is greater than 1

      IF (NNBI.LT.2) THEN
        write(6,*) 'resetting because of small NNB'
        write(6,*) 'current NNB is ',NNBI
        RSI = RSI*1.5
        GO TO 1
      END IF  

*     add particle to other interacting neighbor list as well

      DO J2 = 2,(NNBI+1)
        NNBJ2 = LIST(1,J2)
        LIST(1,J2) = NNBJ2 + 1
        LIST((NNBJ2+2),J2) = I
      END DO 


      FI(1:3,I) = 0.0D0
      D1(1:3,I) = 0.0D0

      DO L = 2,NNBI+1

        K = LIST(L,I)
        write(6,*) 'K is ',K
        A1 = X(1,K) - X(1,I)
        A2 = X(2,K) - X(2,I)
        A3 = X(3,K) - X(3,I)
        DV(1) = XDOT(1,K) - XDOT(1,I)
        DV(2) = XDOT(2,K) - XDOT(2,I)
        DV(3) = XDOT(3,K) - XDOT(3,I)
        RIJ2 = A1*A1 + A2*A2 + A3*A3

        write(6,*) 'RIJ2 = ',RIJ2    
 
        DR2I = 1.0D0/RIJ2
        DR3I = BODY(K)*DR2I*SQRT(DR2I)
        DRDV = 3.0D0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I

        write(6,*) 'DRs are ',DR2I,DR3I,DRDV

        FI(1,I) = FI(1,I) + A1*DR3I
        FI(2,I) = FI(2,I) + A2*DR3I
        FI(3,I) = FI(3,I) + A3*DR3I

        D1(1,I) = D1(1,I) + (DV(1) - A1*DRDV)*DR3I
        D1(2,I) = D1(2,I) + (DV(2) - A2*DRDV)*DR3I
        D1(3,I) = D1(3,I) + (DV(3) - A3*DRDV)*DR3I

      END DO

      F(K,I) = FI(K,I) + FR(K,I)
      FDOT(K,I) = D1(K,I) + D1R(K,I)
      D0(K,I) = FI(K,I)
      D0R(K,I) = FR(K,I)
      FIDOT(K,I) = D1(K,I)
      FRDOT(K,I) = D1R(K,I)


      write(6,*) 'initialization results'
      write(6,*) BODY(I),X(1:3,I),XDOT(1:3,I)
      write(6,*) FR(1:3,I),D1R(1:3,I)
      write(6,*) FI(1:3,I),D1(1:3,I)

      CALL STEPS2(I,I)

      write(6,*) 'calculated time step is'
      write(6,*) STEP(I),STEPR(I)


      DO 10 K = 1,3
        D2(K,I) = 0.0D0
        D3(K,I) = 0.0D0
        D2R(K,I) = 0.0D0
        D3R(K,I) = 0.0D0
   10 CONTINUE

*     need to add other things as well
      
      RS(I) = RSI

      RETURN

      END
