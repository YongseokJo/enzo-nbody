      SUBROUTINE INPUT
*
*
*       Parameter input.
*       ----------------
*
      USE POINTERS
      INCLUDE 'common6.h'
      
      EXTERNAL VERIFY
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif

*       input needed for input.F
!      N = 1000
      NFIX = 1
      NCRIT = 20
      NRAND = 10000
      NNBOPT = 31
      NRUN = 1 
*       1 is left out... why? 

      ETAI = 0.02D0
      ETAR = 0.02D0
      RS0 = 0.19587108D0
      DTADJ = 1.00D0
!      DELTAT = 1.00D0
      TCRIT = 1.0D10
      QE = 1.0D10
      RBAR = 2.58811228D0
!      ZMBAR = 100.00000000D0
      write(6,*) 'ZMBAR = ',ZMBAR

      KZ(1) = 2
      KZ(2) = 2
      KZ(3) = 1
      KZ(4) = 0
      KZ(5) = 1
      KZ(6) = 0
      KZ(7) = 4
      KZ(8) = 2
      KZ(9) = 0
      KZ(10) = 0

      KZ(11) = 0
      KZ(12) = 0
      KZ(13) = 0
      KZ(14) = 0 !edited
      KZ(15) = 2 !0 for starmake
      KZ(16) = 1
      KZ(17) = 1
      KZ(18) = 0
      KZ(19) = 0
      KZ(20) = 0

      KZ(21) = 1
      KZ(22) = 2
      KZ(23) = 0 !edited
      KZ(24) = 0
      KZ(25) = 0
      KZ(26) = 2
      KZ(27) = 0
      KZ(28) = 0
      KZ(29) = 0
      KZ(30) = 2 ! 0 for starmake

      KZ(31) = 0
      KZ(32) = 0
      KZ(33) = 0
      KZ(34) = 0
      KZ(35) = 1
      KZ(36) = 0
      KZ(37) = 1
      KZ(38) = 1
      KZ(39) = 0
      KZ(40) = 1

      KZ(41) = 0
      KZ(42) = 0
      KZ(43) = 0
      KZ(44) = 0
      KZ(45) = 0
      KZ(46) = 0
      KZ(47) = 0
      KZ(48) = 0
      KZ(49) = 0
      KZ(50) = 0

      DTMIN = 1.0D-05
      RMIN = 1.0D-04
      ETAU = 0.1D0
      ECLOSE = 1.0D0
      GMIN = 5.0D-05
      GMAX = 0.01D0
      SMAX = 1.0D0
*
*       Make a formal call to define input parameters & counters.
C      CALL DEFINE
*
*      IF(rank.eq.0)THEN
*       Read & print the main input parameters.
*         READ (5,*)  N, NFIX, NCRIT, NRAND, NNBOPT, NRUN
C Termination time in physical units, TCRITP, read in nbody6.F
*         READ (5,*)  ETAI, ETAR, RS0, DTADJ, DELTAT, TCRIT,
*     &               QE, RBAR, ZMBAR
*         READ (5,*)  (KZ(J),J=1,50)
*         READ (5,*)  DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX, SMAX
*      END IF

*     Check SMAX to make sure it have correct value
      if(rank.eq.0) THEN
         IF(SMAX.GT.1) THEN
            print*, 'Warning! SMAX > 1.0, reduce to 1.0.'
            SMAX = 1.0
         else
            DTNSMAX = 1.0
 1          IF(SMAX/DTNSMAX.LE.0.75) THEN
               DTNSMAX = 0.5D0*DTNSMAX
               IF(DTNSMAX.GT.1E-19) GO TO 1
            END IF
            SMAX = DTNSMAX
         END IF
      END IF

*
#if MPIINIT
      CALL MPI_BCAST(N,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NFIX,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NCRIT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NRAND,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NNBOPT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NRUN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
*
      CALL MPI_BCAST(KZ(1),50,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
*
      CALL MPI_BCAST(ETAI,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ETAR,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RS0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTADJ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DELTAT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TCRIT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(QE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RBAR,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ZMBAR,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTMIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RMIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      RMIN2 = RMIN**2 
      CALL MPI_BCAST(ETAU,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ECLOSE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(GMIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(GMAX,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(SMAX,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)      
#endif
*

#ifdef ENSEMBLE
      NRAND = NRAND + 17*rank
*      WRITE (6,11)rank,NRAND
*   11 FORMAT (//,' Run configured for ensemble average PE ',I5,
*     *   ' NRAND=',I9)
#endif

*      if(rank.eq.0)then
*         WRITE (6,10)
*   10    FORMAT (/////,15X,'N  NFIX  NCRIT  NRAND  NNBOPT  NRUN')
*         WRITE (6,12)  N, NFIX, NCRIT, NRAND, NNBOPT, NRUN
*   12    FORMAT (/,I16,I6,2I7,I8,I6)
*
C New: (Aug.1998, P.Kroupa)
*         WRITE(6,15)
*   15    FORMAT (//,12X,' ETAI      ETAR      RS0       DTADJ',
*     &                  '     DELTAT',
*     &                  '     TCRITP    TCRIT     QE', 
*     &                  '        RBAR       ZMBAR')
*         WRITE (6,20)  ETAI, ETAR, RS0, DTADJ, DELTAT, TCRITP, TCRIT, 
*     &              QE, RBAR,
*     &              ZMBAR
*   20    FORMAT (/,10X,1P10E10.1)
*
*         WRITE (6,22)
*   22    FORMAT (//,12X,'OPTIONS')
*         WRITE (6,24)  (J,J=1,50)
*   24    FORMAT (/,9X,50I3)
*         WRITE (6,26)  (KZ(J),J=1,50)
*   26    FORMAT (/,9X,50I3)
*         WRITE (6,28)
*   28    FORMAT (//,12X,'DTMIN     RMIN      ETAU      ECLOSE    GMIN',
*     &        '      GMAX     SMAX')
*         WRITE (6,30)  DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX, SMAX
*   30    FORMAT (/,9X,1P7E10.1)
*      end if
      call flush(6)
*
*       Define total particle number & neighbour membership range.
      NTOT = N
      NZERO = N
      NNBMAX = MIN(N/2,LMAX - 3)
      ZNBMIN = MAX(0.01*FLOAT(NNBMAX),1.0)
      ZNBMAX = 0.9*FLOAT(NNBMAX)
*       Save initial ETAI.
      ETA0 = ETAI
      RSMIN = RS0
      RC = RS0
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
      GPRINT(1) = 0.0
      DELTAS = 0.0
*     Suppress this KZ(4) since the output need special analysis tool
C      IF (KZ(4).GT.0) THEN
C*       Read parameters for binary evolution analysis.
C          K = KZ(4)
C          if(rank.eq.0)then
C          READ (5,*)  DELTAS, ORBITS(1), (GPRINT(J),J=1,K)
C          end if
C*
C#if MPIINIT
C      CALL MPI_BCAST(DELTAS,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C      CALL MPI_BCAST(ORBITS(1),9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C      CALL MPI_BCAST(GPRINT(1),9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C#endif
C*
C      if(rank.eq.0)WRITE (6,40)  DELTAS, ORBITS(1), (GPRINT(J),J=1,K)
C   40     FORMAT (//,12X,'DELTAS =',F6.2,'  ORBITS(1) =',F6.2,
C     &                                            '  GPRINT(J) =',9F7.3)
C*       Modify binary output factor by perturbation at different levels.
C          DO 50 L = 2,K
C              ORBITS(L) = ORBITS(1)*(GPRINT(1)/GPRINT(L))**0.3333
C   50     CONTINUE
C      END IF
*
C Old version:
*       Set random number skip for routine DATA.
c      IDUM1 = NRAND
C NEW version (14.08.98, P.Kroupa):
C*       Set random number SEED for routine DATA.
      IDUM1 = -1*NRAND
c+++ Notify others of this change on log file:
C      if(rank.eq.0)then
C      write(6,*)
C      write(6,*)' ****** NOTE: new random number seed initialisation!'
C      write(6,*)' ****** AND new ran2 from new ed. of Press et al.'
C      write(6,*)
C      end if
*
*
*       Save square of c.m. approximation parameter (coupled to GMIN).
      CMSEP2 = GMIN**(-0.666666667)
*
      RETURN
*
      END
