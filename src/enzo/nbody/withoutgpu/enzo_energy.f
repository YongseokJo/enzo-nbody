      SUBROUTINE ENERGY

*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common6.h'
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)

*       Calculate the potential energy.
      ZKIN = 0.D00
      POT = 0.0
      VIR = 0.0

*       the origianal NBODY code uses GPU to calculate the KE and PE
*       later I would also have to include GPU
*       but just for a quick check I am using my knowledge to calcualte KE and PE

