*       common6.
*       -------
*
      IMPLICIT REAL*8  (A-H,O-Z)

      INTEGER NMAX
      PARAMETER (NMAX=100)

      COMMON/NBODY/  X0(3,NMAX),X0DOT(3,NMAX),X(3,NMAX),XDOT(3,NMAX),
     &               F0(3,NMAX),F0DOT(3,NMAX),F(3,NMAX),FDOT(3,NMAX),
     &               BODY(NMAX),XF(3,NMAX),XFDOT(3,NMAX),
     &               D0(3,NMAX),D1(3,NMAX),D2(3,NMAX),D3(3,NMAX)

      COMMON/PARAMS/ TWOPI,ONE3,ONE6,ONE9,ONE12,TINY

