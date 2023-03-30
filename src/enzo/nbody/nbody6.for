      SUBROUTINE NBODY6(comm_enzo, nbody_comm, mpi_rank)
*
*             N B O D Y 6++
*             *************
*
*       Regularized AC N-body code with triple & binary collisions.
*       --------------------------------------------------------
*
*       Hermite integration scheme with block-steps (V 4.0.0 April/99).
*       ------------------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*       Message Passing Version NBODY6++ for Massively Parallel Systems
*       Developed by Rainer Spurzem, ARI, Heidelberg
*       
*       Hybrid parallelization (GPU, AVX/SSE, OpenMP + MPI) 
*       Developed by Long Wang, KIAA, Peking University
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      include 'omp_lib.h'
      
*     For MPI communication by YS Jo
      !USE ISO_C_BINDING
      integer :: nbody_comm
      INTEGER :: comm_enzo

*     added by sykim, parameters from enzo

      INTEGER, parameter:: EN_MAX = 2000
      INTEGER :: EN, IS,IE,EID,I,J
      integer :: istatus(MPI_STATUS_SIZE)

*      REAL*8 EBODY(EN_MAX),EX1(EN_MAX),EX2(EN_MAX),EX3(EN_MAX)
*      REAL*8 EXDOT1(EN_MAX),EXDOT2(EN_MAX),EXDOT3(EN_MAX)
*      REAL*8 EF1(EN_MAX),EF2(EN_MAX),EF3(EN_MAX)

*      REAL*8 EH11(EN_MAX),EH12(EN_MAX),EH13(EN_MAX)
*      REAL*8 EH21(EN_MAX),EH22(EN_MAX),EH23(EN_MAX)
*      REAL*8 EH31(EN_MAX),EH32(EN_MAX),EH33(EN_MAX)
*      REAL*8 EH41(EN_MAX),EH42(EN_MAX),EH43(EN_MAX)

      REAL*8, pointer :: EBODY(:),EX(:,:)
      REAL*8, pointer :: EXDOT(:,:), EF(:,:) !, EH(:,:,:)


*     conversion factors for enzo code unit -> cgs

      REAL*8 EMU,ELU,EVU,ETU


*     conversion factors for astronomical -> nbody

      REAL*8 LENGTHU,MASSU,VELU,TIMEU

      REAL*8 EDT

*     end added by sykim



      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
#ifdef DEBUG
*     --10/03/14 19:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      COMMON/adt/ adtime,dumptime,dprintt,dtprint,namep
*     --10/03/14 19:41-lwang-end----------------------------------------*
#endif
      EXTERNAL MERGE


#ifdef PARALLEL
      ! by YS should be 1 
#define MPIINIT 0
#else
#ifdef ENSEMBLE
      ! by YS should be 1
#define MPIINIT 0
#else
#define MPIINIT 0
#endif
#endif


*     added by sykim
*     recieve input from enzo and place it 
*     at necessary variables


*     conventional units are pc, Msun, km/s and Myr
*     the length unit is Rvir, mass :unit is total mass, et cetera


      write (0,*) 'In the Nbody6 '
*     Massive MPI Communication with Enzo code! (subroutine later) by YS
      write (0,*) 'Initialize Arrays'

*     Massive MPI Communication with Enzo code! (subroutine later) by YS
      write (0,*) 'Communication Starts'

      ! recieve the number of particles 
*----------------------------------------------------------------------------------*
*     MPI starts, refer to PrepareNbodyComputation.C for the counter part  by YS
      call MPI_RECV(EN, 1, MPI_INTEGER, 0, 100, comm_enzo, istatus,
     &           ierr)

      write (0,*) 'fortran: Number of Nbody Particles on Fortran', EN
      allocate(EBODY(EN))
      call MPI_RECV(EBODY, EN, MPI_DOUBLE_PRECISION, 0, 200, comm_enzo,
     &     istatus, ierr)
      allocate(EX(3,EN))
      allocate(EXDOT(3,EN))
      allocate(EF(3,EN))
      DO  I = 1,3
        call MPI_RECV(EX(I,:), EN, MPI_DOUBLE_PRECISION, 0, 300,
     &           comm_enzo, istatus,ierr)
        call MPI_RECV(EXDOT(I,:), EN, MPI_DOUBLE_PRECISION, 0, 400,
     &           comm_enzo, istatus,ierr)
        call MPI_RECV(EF(I,:), EN, MPI_DOUBLE_PRECISION, 0, 500,
     &            comm_enzo, istatus,ierr)
      END DO
      call MPI_RECV(EDT, 1, MPI_DOUBLE_PRECISION, 0, 600,
     &            comm_enzo, istatus,ierr)
      call MPI_RECV(ETU, 1, MPI_DOUBLE_PRECISION, 0, 700,
     &            comm_enzo, istatus,ierr)
      write (0,*) 'fortran: mass=', EBODY(1), 'X=', EX(1,1), 
     &             ', V=',EXDOT(1,1)
*     MPI done!
*----------------------------------------------------------------------------------*

      write (6,*) 'before conversion',EBODY(1),EX(1,1),EXDOT(1,1)

      MASSU = 0.0D0

      DO I = 1,EN
          MASSU = MASSU + EBODY(I)*EMU/(1.9891D33)
      END DO

*     need to calculate virial radius and put that into LENGTHU0
*     ** code for calculating virial radius** <- should be added

      LENGTHU = 2.58811228D0
      VELU = 6.557D0*((MASSU/LENGTHU)**(0.5D0))/(100.0D0)
      TIMEU = 14.94D0*(((LENGTHU)**3.0D0/MASSU)**(0.5D0))

      write (6,*) 'scaling',LENGTHU,MASSU,VELU,TIMEU

*     determine how much steps should be run depending on approximate
*     scaling
 
      N = EN
      TCRIT = EDT*ETU/(TIMEU*(3.1556952D13))

      write (6,*) 'timesteps',TCRIT


      DO 7 IS = 1,N
         BODY(IS) = EBODY(IS)*EMU/(MASSU*1.9891D33)
         DO J = 1,3
          X(J,IS) = EX(J,IS)*ELU/LENGTHU/(3.0857D18)
          XDOT(J,IS) = EXDOT(J,IS)*EVU/VELU/(1D5)
         END DO
    7 CONTINUE

      write (6,*) 'after conversion',BODY(1),X(1,1),XDOT(1,1)
*     end added by sykim


#if MPIINIT
*       Initialize MPI
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,group,ierr)
      CALL MPI_GROUP_SIZE(group,isize,ierr)
      CALL MPI_GROUP_RANK(group,rank,ierr)
      ibarcount=0
      write(6,11) rank,isize,group
 11   format('MPI-initial: This is rank=',I6,' size=',I6,' group=',I11)
#endif
*
*       Initialize the timer.
      CALL CPUTIM(ttota)

*       Get threads number
#ifdef OMP
!$omp parallel 
      icore=OMP_get_num_threads()
!$omp end parallel
      PRINT*,'RANK: ',rank,' OpenMP Number of Threads: ',icore
#else
      icore = 1
#endif

#ifdef PARALLEL
      call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif      
      call flush(6)
*

*      disabled by sykim, hard-code input parameters
*       Read start/restart indicator & CPU time.
*      IF(rank.eq.0)READ (5,*)  KSTART, TCOMP, TCRITP,
*     *    isernb,iserreg,iserks

*     added by sykim, code for defining input parameters

*     parameters needed for nbody6.F

      KSTART = 1
      TCOMP = 100000.0D0
      TCRITP = 1.0D8

      isernb = 40
      iserreg = 40
      iserks = 640

*    end added by sykim


#ifdef DEBUG
*     --10/03/14 19:41-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*      disabled by sykim, hard-code input parameters
*      if(rank.eq.0) read(5,*) adtime,dumptime,dprintt,dtprint,namep
#if MPIINIT
      CALL MPI_BCAST(adtime,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dumptime,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dprintt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dtprint,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(namep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif      
*     --10/03/14 19:41-lwang-end----------------------------------------*
#endif
*
#if MPIINIT
      CALL MPI_BCAST(isernb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iserreg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iserks,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KSTART,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TCOMP,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TCRITP,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
*
      isernb = max(isize,isernb*icore)
      iserreg = max(isize,iserreg*icore)
*      iserks = max(isize,iserks*icore)
      IF(rank.eq.0) then
        PRINT*,' iserreg,isernb,iserks,ithread=',iserreg,isernb,iserks,
     &        ithread
      end if
#endif
*
      IF (KSTART.EQ.1) THEN
*
*       Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          ttadj = ttadj + (tt8-tt7)*60.
      ELSE
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
*       Backup kstart value before call mydump
          KSTART0 = KSTART
          CALL MYDUMP(0,1)
*       Reset kstart to input value
          KSTART = KSTART0
*       
          IF (NDUMP.GE.3) STOP
*       Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0 
*       Set IPHASE = -1 for new NLIST in routine INTGRT (Hermite version).
          IPHASE = -1
*
*       Initialize evolution parameters which depend on metallicity.
          IF (KZ(19).GE.3) THEN
              CALL ZCNSTS(ZMET,ZPARS)
          END IF
*
*       Check reading modified restart parameters (KSTART = 3, 4 or 5).
          IF (KSTART.GT.2) THEN
              CALL MODIFY
          END IF
*
*       Open all other files.
          CALL FILE_INIT(0)
*
*       If no explicit new TCRIT given just go for another TCRIT of common block.
          TTOT = TIME + TOFF
          TCRIT = TTOT + TCRIT
          if(rank.eq.0)then
             WRITE (6,10) TTOT/TCR0, TIME/TCR0, TCRIT/TCR0, TTOT, TIME,
     &            TCRIT
             WRITE (6,20)  DTADJ, DELTAT, TADJ, TNEXT, TCRIT, QE
             WRITE (6,30)  ETAI, ETAR, ETAU, DTMIN, RMIN, NNBOPT
 10          FORMAT (' START AT TTOT/TIME ',2F16.8,' STOP INTENDED AT ',
     &            F16.8,' TCR0',/,' START AT TTOT/TIME ',2F16.8,
     &            ' STOP INTENDED AT ',F16.8,' NBODY-UNITS ',/)
 20          FORMAT (/,7X,'RESTART PARAMETERS:   DTADJ =',F7.3,
     &            '  DELTAT =',F7.3,'   TADJ =',F7.3,'   TNEXT =',
     &            F7.3,'  TCRIT =',F7.1,'  QE =',1PE9.1)
 30          FORMAT (/,7X,'                      ETAI =',F7.3,
     &            '  ETAR =',F7.3,'  ETAU =',F7.3,'  DTMIN =',1PE9.1,
     &            '  RMIN =',E9.1,' NNBOPT =',I5,/)
          end if
*
*       Find massive back hole after restart
          IF (KZ(24).EQ.1) call IMBHRESTART

      END IF
*
* (R.Sp.)Set time flag and step number flags for beginning of run
      TINIT = TTOT
      NIR = NSTEPI
      NIB = NSTEPB
      NRGL = NSTEPR
      NKS = NSTEPU
*
      call cputim(tt2)
      ttinitial = ttinitial + (tt2-ttota)*60.
*       Advance solutions until next output or change of procedure.
    1 CONTINUE

*     sykim: need to recieve force from enzo here!!
*     update the common variables to the recieved enzo variables. 

*      call MPI_RECV(EN, 1, MPI_INTEGER, 0, 100, comm_enzo, istatus,
*     &           ierr)
*      write (0,*) 'fortran: Number of Nbody Particles on Fortran', EN
*      allocate(EBODY(EN))
*      call MPI_RECV(EBODY, EN, MPI_DOUBLE_PRECISION, 0, 200, comm_enzo,
*     &     istatus, ierr)
*      allocate(EX(3,EN))
*      allocate(EXDOT(3,EN))

*----------------------------------------------------------------------------------*
*     MPI starts, refer to PrepareNbodyComputation.C:219 for the counter part  by YS
      allocate(EF(3,EN))
      DO I = 1,3 
        call MPI_RECV(EF(I,:), EN, MPI_DOUBLE_PRECISION, 0, 500,
     &            comm_enzo, istatus,ierr)
      END DO
      call MPI_RECV(EDT, 1, MPI_DOUBLE_PRECISION, 0, 600,
     &            comm_enzo, istatus,ierr)
      call MPI_RECV(ETU, 1, MPI_DOUBLE_PRECISION, 0, 700,
     &            comm_enzo, istatus,ierr)
*       TCRIT = TCRIT + dtENZO
*     for SY
      TCRIT = TCRIT + EDT*ETU/(TIMEU*(3.1556952D13)) ! should be fixed
      write (0,*) 'fortran: force=', EF(1,1)
*     MPI done!
*----------------------------------------------------------------------------------*


      call cputim(tt1)
*
*     --08/27/13 16:31-lwang-debug--------------------------------------*
***** Note: -----------------------------------------------------------**
*      if(time.ge.20.8) print*,rank,'aint ',time
*     --08/27/13 16:31-lwang-end-debug----------------------------------*
      CALL INTGRT
*     --08/27/13 16:32-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*      if(time.ge.20.8) print*,rank,'bint ',time
*     --08/27/13 16:32-lwang-end-debug----------------------------------*
*
      call cputim(tt2)
      ttint = ttint + (tt2-tt1)*60.
*



      IF (IPHASE.EQ.1) THEN
*       Prepare new KS regularization.
      call cputim(tt1)
          CALL KSREG
          CALL FLUSH(6)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
      ttksinit = ttksinit + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.2) THEN
*       Terminate KS regularization.
      call cputim(tt1)
          CALL KSTERM
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
      ttksterm = ttksterm + (tt2-tt1)*60.

      ELSE IF (IPHASE.EQ.3) THEN
*       Perform energy check & parameter adjustments and print
*       diagnostics.
        call cputim(tt7)
        CALL ADJUST

*       added by sykim. return datas to ENZO if TIME>TCRIT
        IF (IPHASE.EQ.13) THEN

*       need to make a new extrapolation scheme... in progress
*       after extrapolation,  update the EBODY, EX, EXDOT that would be passed onto ENZO
          DO EID = 1,N
            IE = NAME(EID)
              EBODY(IE) = BODY(IE)*MASSU*1.9891D33/EMU
              DO J = 1,3
                  EX(J,IE) = (X(J,IE)-RDENS(J))*LENGTHU*3.0857D18/ELU
                  EXDOT(J,IE) = XDOT(J,IE)*VELU*1D5/EVU
                END DO
                END DO 
*           17 CONTINUE ! this does not work for SY

*      sykim: need to add the MPI return scheme here!  
*      and do not return or end the program
*            RETURN
*            call MPI_RECV(EBODY, EN, MPI_DOUBLE_PRECISION, 0, 200, comm_enzo,
*     &     istatus, ierr)
*----------------------------------------------------------------------------------*
*     MPI starts, refer to FinalizeNbodyComputation.C for the counter part  by YS
          DO  I = 1,3
           call MPI_SSEND(EX(I,:), EN, MPI_DOUBLE_PRECISION, 0, 300,
     &           comm_enzo)
           call MPI_SSEND(EXDOT(I,:), EN, MPI_DOUBLE_PRECISION, 0, 400,
     &           comm_enzo)

          END DO
*          deallocate(EF)
          write (0,*) 'fortran: mass=', EBODY(1), 'X=', EX(1,1), 
     &             ', V=',EXDOT(1,1)
*            deallocate(EBODY(EN))
*            deallocate(EX(3,EN))
*            deallocate(EXDOT(3,EN))
*            deallocate(EH(3,4,EN))
*     MPI done!
*----------------------------------------------------------------------------------*
          IPHASE = 19 

          GO TO 1
        END IF
        call cputim(tt8)
        ttadj = ttadj + (tt8-tt7)*60.

*
      ELSE IF (IPHASE.EQ.4) THEN
*       Switch to unperturbed three-body regularization.
      call cputim(tt1)
          ISUB = 0 
          CALL TRIPLE(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.5) THEN
*       Switch to unperturbed four-body regularization.
      call cputim(tt1)
          ISUB = 0
          CALL QUAD(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
*       Adopt c.m. approximation for inner binary in hierarchical triple.
      ELSE IF (IPHASE.EQ.6) THEN
      call cputim(tt1)
          CALL MERGE
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.7) THEN
*       Restore old binary in hierarchical configuration.
      call cputim(tt1)
          CALL RESET
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
*
*       Begin chain regularization.
      ELSE IF (IPHASE.EQ.8) THEN
      call cputim(tt1)
          ISUB = 0
          CALL CHAIN(ISUB)
      call cputim(tt2)
      ttks = ttks + (tt2-tt1)*60.
      END IF
*
*       Continue integration.
      GO TO 1
*
      call MPI_FINALIZE(ierr)
      END
