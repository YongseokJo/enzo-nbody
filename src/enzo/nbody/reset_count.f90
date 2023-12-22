      SUBROUTINE RESET_COUNT

*     routine for re-setting count variables after single run of NBODY
*     added by sykim


      USE POINTERS
      INCLUDE 'common6.h'
      

      INTEGER I,J


      NSTEPI = 0
      NSTEPR = 0
      NSTEPU = 0
      NNPRED = 0
      NBCORR = 0
      NBFULL = 0
      NBVOID = 0

      NNTB = 0
      NBSMIN = 0
      NLSMIN = 0
      NBDIS = 0
      NBDIS2 = 0
      NCMDER = 0
      NBDER = 0

      NFAST = 0
      NBFAST = 0
      NBLOCK = 0
      NBPRED = 0
      NICONV = 0
      NCHAIN = 0
      NSTEPC = 0

      NKSTRY = 0
      NKSREG = 0
      NKSHYP = 0
      NKSPER = 0
      NPRECT = 0
      NEWKS = 0
      NKSMOD = 0

      NTTRY = 0
      NTRIP = 0
      NQUAD = 0
      NMERG = 0
      NSTEPT = 0
      NSTEPQ = 0
      NDISS = 0
      NTIDE = 0

      NCOLL = 0
      NSYNC = 0
      NSESC = 0
      NBESC = 0
      NMESC = 0
      NTIMER = 0
      NSTEPS = 0
      NPRINT = 0

      NDUMP = 0
      NBPREV = 0
      NEWHI = 0
      NSTEPB = 0
      NBFLUX = 0
      NMTRY = 0
      NWARN = 0

      NIRECT = 0
      NURECT = 0
      NBRECT = 0
      NRRECT = 0
      KSMAG = 0

      NDUM(1) = 0
      NDUM(2) = 0

      DO I = 1,10
         NPOP(I) = 0
      END DO

      NBLCKR = 0
      NIRRF = 0
      KSTART = 0
      
      DO J = 1,97
         NDUMMY(J) = 0
      END DO

      RETURN

      END
