       SUBROUTINE ADD_STAR(INEW,KCASE,BODYNEW,XNEW,VNEW)
*
*
*       Particle removal.
*       -----------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tlist.h'

      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall

      INTEGER INEW,I,J,J1,J2,J3,K,K2,KCASE
      INTEGER IPAIR,I1,I2,ICM
      INTEGER J4,K1
      INTEGER L,LK,NNB

      REAL*8 BODYNEW,XNEW(3),VNEW(3)

      REAL*8 A(6)
      REAL*8 STEPI,STEPK

      LOGICAL RM_FLAG,B_FLAG

      write(6,*) 'ADD_STAR starts'

*       Move up all COMMON variables (escaper or old c.m. & KS comps) 
*       to make space for calculation of new star

      write(6,*) 'IFIRST = ',IFIRST
      write(6,*) 'NTOT = ',NTOT
      write(6,*) 'N = ',N
      write(6,*) 'NBIN = ',NBIN
      write(6,*) 'NGHOSTS = ',NGHOSTS
      write(6,*) 'NPAIRS = ',NPAIRS
      write(6,*) 'NXTLIMIT = ',NXTLIMIT 

      I = INEW

*    checking things related with KSPAIR

      IPAIR = KSPAIR
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ICM = N + IPAIR

      write(6,*) "KSPAIR = ",KSPAIR
      write(6,*) "I1 = ",I1,"I2 = ",I2
      write(6,*) "id of pair is ",NAME(I1),NAME(I2)
      write(6,*) "index of cm is",ICM
      write(6,*) "LIST of I1+1 is",LIST(2,I1+1)

*
      IF(NTOT.EQ.N) GO TO 6

      DO 5 J1 = NTOT,I,-1

          write(6,*) 'name of shifted particle is',NAME(J1)
          J = J1 + 1
          DO 2 K = 1,3
              X(K,J) = X(K,J1)
              X0(K,J) = X0(K,J1)
              X0DOT(K,J) = X0DOT(K,J1)
              XDOT(K,J) = XDOT(K,J1)
              F(K,J) = F(K,J1)
              FDOT(K,J) = FDOT(K,J1)
              FI(K,J) = FI(K,J1)
              FIDOT(K,J) = FIDOT(K,J1)
              D0(K,J) = D0(K,J1)
              D1(K,J) = D1(K,J1)
              D2(K,J) = D2(K,J1)
              D3(K,J) = D3(K,J1)
              FR(K,J) = FR(K,J1)
              FRDOT(K,J) = FRDOT(K,J1)
              D0R(K,J) = D0R(K,J1)
              D1R(K,J) = D1R(K,J1)
              D2R(K,J) = D2R(K,J1)
              D3R(K,J) = D3R(K,J1)
    2    CONTINUE
*
          TPRED(J) = TPRED(J1)
          BODY(J) = BODY(J1)
          RS(J) = RS(J1)
          RADIUS(J) = RADIUS(J1)
          TEV(J) = TEV(J1)
          TEV0(J) = TEV0(J1)
          BODY0(J) = BODY0(J1)
          EPOCH(J) = EPOCH(J1)
          SPIN(J) = SPIN(J1)
          ZLMSTY(J) = ZLMSTY(J1)
          KSTAR(J) = KSTAR(J1)
          NAME(J) = NAME(J1)+1
          STEP(J) = STEP(J1)
          STEPR(J) = STEPR(J1)
          T0(J) = T0(J1)
          T0R(J) = T0R(J1)

*
*       Transfer unmodified neighbour list.
          NNB = LIST(1,J1) + 1

          DO 3 L = 1,NNB
              LIST(L,J) = LIST(L,J1)
    3     CONTINUE
    5 CONTINUE

*      need to change the cm particle index in LIST

      DO J4 = 1,NTOT
         NNB = LIST(1,J1) + 1
         DO K1 = 2,NNB
           IF (LIST(K,J1).GE.I) THEN
              write(6,*) 'cm components are in LIST'
              write(6,*) 'need to change this!'
              write(6,*) 'LIST is ',LIST(K,J1)
              write(6,*) 'Name is ',NAME(LIST(K,J1))
              LIST(K,J1) = LIST(K,J1) + 1
           END IF
         END DO
      END DO

*      need to increase index of cmbody already existing in NXTLST

      DO J3 = 1,NXTLIMIT
         IF (NXTLST(J3).GT.N) THEN
            write(6,*) "cmbody already exists in NXTLST"
            write(6,*)  "index in NXTLST is ",J3
            write(6,*) "NAME is ",NAME(NXTLST(J3))
            NXTLST(J3) = NXTLST(J3) + 1
         END IF
      END DO

    6 CONTINUE

*       include the new particle information in the appropriate positions

       DO K = 1,3
         X(K,I) = XNEW(K)
         X0(K,I) = XNEW(K)
         X0DOT(K,I) = VNEW(K)
         XDOT(K,I) = VNEW(K)
       END DO
*
       TPRED(I) = TIME
       BODY(I) = BODYNEW
       RS(I) = RS0
       RADIUS(I) = 0.0D0
       TEV(I) = 0.0D0
       TEV0(I) = 0.0D0
       BODY0(I) = BODYNEW
       EPOCH(I) = 0.0D0
       SPIN(I) = 0.0D0
       KSTAR(I) = 0.0D0
       NAME(I) = I
       STEP(I) = 0.0D0
       STEPR(I) = 0.0D0
       T0(I) = TIME
       T0R(I) = TIME

*       calculate information related to intialization


      write(6,*) 'going into addinit'

      CALL ADDINIT(I)

      write(6,*) 'exiting addinit'


      write(6,*) 'adding new star into NXTLST'

      CALL add_tlist(I,STEP,DTK)

      write(6,*) 'end adding star to NXTLST'


*     For ghost list
      IF(NGHOSTS.GT.0) THEN

         write(6,*) "ghost not zero when adding - may be problematic"

         NXTLST(NXTLIMIT+1) = NXTLST(NXTLIMIT+1+NGHOSTS)
         L = NXTLIMIT + 1

    7    CONTINUE

         IF (NXTLST(L).GE.I) NXTLST(L) = NXTLST(L) + 1
         L = L + 1
         IF (L.LE.NXTLIMIT+NGHOSTS) GO TO 7
      END IF

*     For NLSTDELAY

      IF (NLSTDELAY(1).GT.0) THEN
         L = 2
    9    CONTINUE 
         IF (NLSTDELAY(L).GE.I) NLSTDELAY(L) = NLSTDELAY(L) + 1
         L = L + 1
         IF (L.LE.NLSTDELAY(1)+1) GO TO 9
      END IF

      

*       From now on it is code from remove.f

*       Correct force & first derivative of neighbours
      NNB1 = LIST(1,I) + 1
      DO 25 L = 2,NNB1
          J = LIST(L,I)
          RIJ2 = 0.0D0
          A7 = 0.0D0
          DO 21 K = 1,3
              A(K) = X(K,I) - X(K,J)
              A(K+3) = XDOT(K,I) - XDOT(K,J)
              RIJ2 = RIJ2 + A(K)**2
              A7 = A7 + A(K)*A(K+3)
   21     CONTINUE
          A8 = BODY(I)/(RIJ2*SQRT(RIJ2))
          DO 22 K = 1,3
              A(K+3) = (A(K+3) + 3.0*A7*A(K)/RIJ2)*A8
              F(K,J) = F(K,J) + 0.5*A(K)*A8
              FI(K,J) = FI(K,J) + A(K)*A8
              FDOT(K,J) = FDOT(K,J) + ONE6*A(K+3)
              D1(K,J) = D1(K,J) + A(K+3)
   22     CONTINUE
   25 CONTINUE  

      write(6,*) 'ADD_STAR finished'

      RETURN

      END
