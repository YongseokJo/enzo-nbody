      SUBROUTINE UNSCALE(N)
*
*
*       Scaling to new units. 
*       Written on 2022.09.27
*       ---------------------
*
      INCLUDE 'enzo_common6.h'
      INTEGER N


*    We got the conversion factors, so let's convert it back to original scale
*    for conversion, I will assume input of mass = Msun, pos = pc, vel = km/s

      DO 10 I = 1,N
          BODY(I) = BODY(I)*ZMASS
 10   CONTINUE


      DO 20 I = 1,N
          DO 15 K = 1,3
              X(K,I) = X(K,I)*RVIR
              XDOT(K,I) = XDOT(K,I)*VELU
 15       CONTINUE
 20   CONTINUE

*       Adjust coordinates and velocities back to its original place (not cm frame)

      DO 30 I = 1,N
          DO 25 K = 1,3
              X(K,I) = X(K,I) + CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) + CMRDOT(K)/ZMASS
 25       CONTINUE
 30   CONTINUE

      RETURN

      END