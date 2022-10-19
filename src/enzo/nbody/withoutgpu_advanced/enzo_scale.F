      SUBROUTINE SCALE
*
*
*       Scaling to new units. 
*       Written on 2022.09.05
*       Modified on 2022.09.27
*       ---------------------
*
      INCLUDE 'common6.h'

      ZMASS = 0.0D0

      LENGTHU = 1.0D0
      MASSU = 1.0D0
      VELU = 1.0D0
      TIMEU = 1.0D0

      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
*
*     Form total mass and centre of mass displacements.
*     CMR, CMRDOT is in conventional units

      DO 30 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 25 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X0(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*X0DOT(K,I)
   25     CONTINUE
   30 CONTINUE

*       Adjust coordinates and velocities to c.m. rest frame.
*       (if start model data were read from other source do NOT scale)(R.Sp.)

      DO 40 I = 1,N
          DO 35 K = 1,3
              X0(K,I) = X0(K,I) - CMR(K)/ZMASS
              X0DOT(K,I) = X0DOT(K,I) - CMRDOT(K)/ZMASS
 35       CONTINUE
 40   CONTINUE

*     Scale masses to standard units of <M> = 1/N.

      DO 50 I = 1,N
          BODY(I) = BODY(I)/ZMASS
 50   CONTINUE


*     Originally, the scaling code calculates all virial energy, tidal energy et cetra
*     However, when I was working with NBODY6, the input was usually already scaled
*     Don't know what we should do yet, but since the energy calculation code is quite complex
*     I'll apply what I have done in my python code to convert the scales

*     for conversion, I will assume input of mass = Msun, pos = pc, vel = km/s


      LENTHU = RVIR
      MASSU = ZMASS
      VELU = (GRAV*MSUN*MASSU/(RVIR*PC))**(0.5)/(1e5)
      TIMEU = ((RVIR*PC)**3/(MASSU*MSUN*GRAV))**(0.5)/(YR)


*    We got the conversion factors, so let's convert!

      DO 60 I = 1,N
          DO 55 K = 1,3
              X0(K,I) = X0(K,I)/LENGTHU
              X0DOT(K,I) = X0DOT(K,I)/VELU
 55       CONTINUE
 60   CONTINUE

      RETURN

      END