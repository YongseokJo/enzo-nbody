      SUBROUTINE MYDUMP(II,J)
*
*
*       COMMON save or read.
*       --------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      INCLUDE 'params.h'
      INCLUDE 'timing.h'
      INCLUDE 'mpi_base.h'
      INCLUDE 'tlist.h'
#ifdef TT
      INCLUDE 'tt.h'
#endif
      PARAMETER  (NA=84,NB=168,NC=530,ND=392+MLR+MLD+MLV,NE=24,NM=40,
     &    NG=132,NL=99,NO=20*MCL+16,NP=32*NTMAX,NQ=31*MMAX,
     &    NS=44*MMAX)
      REAL*4  A,B,C,D,E,G,L,M,O,P,Q,S
      INTEGER K,I,NTSAVE,NPARTMP(10)
*
      COMMON/NAMES/  NTOT,NPAIRS,NTTOT,A(NA)
      COMMON/COUNTS/ B(NB)
      COMMON/PARAMS/ C(NC)
      COMMON/STARS/  D(ND)
      COMMON/PLPOT/  E(NE)
      COMMON/BLOCKS/ G(NG),KVEC(2*KMAX)
      COMMON/RAND2/  L(NL)
      COMMON/GALAXY/ M(NM)
      COMMON/CLOUDS/ O(NO)
      COMMON/MODES/  P(NP)
      COMMON/RCHE/   Q(NQ)
      COMMON/BINARY/ S(NS)

      COMMON/NBODY/  X(3,NMAX),XDOT(3,NMAX),X0(3,NMAX),X0DOT(3,NMAX),
     &               F(3,NMAX),FDOT(3,NMAX),BODY(NMAX),RS(NMAX),
     &               FI(3,NMAX),D1(3,NMAX),D2(3,NMAX),D3(3,NMAX),
     &               FR(3,NMAX),D1R(3,NMAX),D2R(3,NMAX),D3R(3,NMAX),
     &               STEP(NMAX),T0(NMAX),STEPR(NMAX),T0R(NMAX),
     &               TIMENW(NMAX),RADIUS(NMAX),TEV(NMAX),TEV0(NMAX),
     &               BODY0(NMAX),EPOCH(NMAX),SPIN(NMAX),XSTAR(NMAX),
     &               ZLMSTY(NMAX),FIDOT(3,NMAX),D0(3,NMAX),
     &               FRDOT(3,NMAX),D0R(3,NMAX),KSTAR(NMAX),NAME(NMAX),
     &               LIST(LMAX,NMAX)
*
      COMMON/PAIRS/  U(4,KMAX),U0(4,KMAX),UDOT(4,KMAX),FU(4,KMAX),
     &               FUDOT(4,KMAX),FUDOT2(4,KMAX),FUDOT3(4,KMAX),
     &               H(KMAX),HDOT(KMAX),HDOT2(KMAX),HDOT3(KMAX),
     &               HDOT4(KMAX),DTAU(KMAX),TDOT2(KMAX),TDOT3(KMAX),
     &               R(KMAX),R0(KMAX),GAMMA(KMAX),SF(7,KMAX),H0(KMAX),
     &               FP0(4,KMAX),FD0(4,KMAX),KBLIST(10*KMAX),
     &               KSLOW(KMAX),TBLIST
*

         
*
      RETURN
*
      END
