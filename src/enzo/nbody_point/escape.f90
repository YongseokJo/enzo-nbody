      SUBROUTINE ESCAPE
*
*
*       Escaper removal.
*       ----------------
*
      USE POINTERS
      INCLUDE 'common6.h'
      
#ifdef TT
      include 'tt.h'
#endif
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
      CHARACTER*11  WHICH1
      LOGICAL FIRST
#ifdef TT
      LOGICAL RFLAG
#endif
      SAVE FIRST
      DATA FIRST/.true./
*
      END

