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
      integer :: nbody_comm, NCOMM
      INTEGER :: comm_enzo, ECOMM
      LOGICAL :: FIRST=.true.
*     added by sykim, parameters from enzo

      INTEGER, parameter:: EN_MAX = 2000
      INTEGER :: EN, IS,IE,EID,I,J
      integer :: istatus(MPI_STATUS_SIZE)


      REAL*8, pointer :: EBODY(:),EX(:,:)
      REAL*8, pointer :: EXDOT(:,:), EF(:,:) !, EH(:,:,:)


*     conversion factors for enzo code unit -> cgs

      REAL*8 EMU,ELU,EVU,ETU


*     conversion factors for astronomical -> nbody

      REAL*8 LENGTHU,MASSU,VELU,TIMEU

      REAL*8 EDT

*     end added by sykim

*     common variables for MPI by YS
      COMMON/MPI/ ECOMM, NCOMM, FIRST

      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
#ifdef DEBUG
*     --10/03/14 19:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      COMMON/adt/ adtime,dumptime,dprintt,dtprint,namep
*     --10/03/14 19:41-lwang-end----------------------------------------*
#endif
      EXTERNAL MERGE


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
      ECOMM = comm_enzo
      NCOMM = nbody_comm
      call MPI_RECV(EN, 1, MPI_INTEGER, 0, 100, ECOMM, istatus,
     &           ierr)

      write (0,*) 'fortran: Number of Nbody Particles on Fortran', EN
      allocate(EBODY(EN))
      call MPI_RECV(EBODY, EN, MPI_DOUBLE_PRECISION, 0, 200, ECOMM,
     &     istatus, ierr)
      allocate(EX(3,EN))
      allocate(EXDOT(3,EN))
      allocate(EF(3,EN))
      DO  I = 1,3
        call MPI_RECV(EX(I,:), EN, MPI_DOUBLE_PRECISION, 0, 300,
     &           ECOMM, istatus,ierr)
        call MPI_RECV(EXDOT(I,:), EN, MPI_DOUBLE_PRECISION, 0, 400,
     &           ECOMM, istatus,ierr)
        call MPI_RECV(EF(I,:), EN, MPI_DOUBLE_PRECISION, 0, 500,
     &            ECOMM, istatus,ierr)
      END DO
      call MPI_RECV(EDT, 1, MPI_DOUBLE_PRECISION, 0, 600,
     &            ECOMM, istatus,ierr)
      call MPI_RECV(ETU, 1, MPI_DOUBLE_PRECISION, 0, 700,
     &            ECOMM, istatus,ierr)
      call MPI_RECV(ELU, 1, MPI_DOUBLE_PRECISION, 0, 800,
     &            ECOMM, istatus,ierr)
      call MPI_RECV(EMU, 1, MPI_DOUBLE_PRECISION, 0, 900,
     &            ECOMM, istatus,ierr)
      call MPI_RECV(EVU, 1, MPI_DOUBLE_PRECISION, 0, 1000,
     &            ECOMM, istatus,ierr)
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
*      CALL ENZO_TO_NB(EX, EXDOT)
      
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

      call flush(6)
*


*     added by sykim, code for defining input parameters

*     parameters needed for nbody6.F

      KSTART = 1
      TCOMP = 100000.0D0
      TCRITP = 1.0D8

      isernb = 40
      iserreg = 40
      iserks = 640

*    end added by sykim


*
*       Read input parameters, perform initial setup and obtain output.
      CPU = TCOMP
      CALL START
      call cputim(tt7)
      CALL ADJUST
      call cputim(tt8)
      ttadj = ttadj + (tt8-tt7)*60.

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

*      call MPI_RECV(EN, 1, MPI_INTEGER, 0, 100, ECOMM, istatus,
*     &           ierr)
*      write (0,*) 'fortran: Number of Nbody Particles on Fortran', EN
*      allocate(EBODY(EN))
*      call MPI_RECV(EBODY, EN, MPI_DOUBLE_PRECISION, 0, 200, ECOMM,
*     &     istatus, ierr)
*      allocate(EX(3,EN))
*      allocate(EXDOT(3,EN))



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

*           17 CONTINUE ! this does not work for SY ,why?

*      sykim: need to add the MPI return scheme here!  
*      and do not return or end the program
*            RETURN
*            call MPI_RECV(EBODY, EN, MPI_DOUBLE_PRECISION, 0, 200, ECOMM,
*     &     istatus, ierr)

          ! for SY enzo comm

          DO EID = 1,N
          IE = NAME(EID)
          EBODY(IE) = BODY(IE)*MASSU*1.9891D33/EMU
          DO J = 1,3
          EX(J,IE) = (X(J,IE)-RDENS(J))*LENGTHU*3.0857D18/ELU
          EXDOT(J,IE) = XDOT(J,IE)*VELU*1D5/EVU
          END DO
          END DO 
*----------------------------------------------------------------------------------*
      ! SEND
      allocate(EX(3,N))
      allocate(EDOT(3,N))
*     sorting
*     CALL NB_TO_ENZO(EX, EDOT)
*     MPI starts, refer to FinalizeNbodyComputation.C for the counter part  by YS
          write (0,*) 'fortran: write out' 
          write (6,*) 'fortran: write out' 

          DO  I = 1,3
          call MPI_SSEND(EX(I,:), N, MPI_DOUBLE_PRECISION, 0, 300,
     &           ECOMM, ierr)
          call MPI_SSEND(EXDOT(I,:), N, MPI_DOUBLE_PRECISION, 0, 400,
     &           ECOMM, ierr)
          END DO
*          deallocate(EF)
          write (0,*) 'fortran: X=', EX(1,1), ', V=',EXDOT(1,1)
*     MPI done!
*----------------------------------------------------------------------------------*

*----------------------------------------------------------------------------------*
          ! RECEIVE
*     MPI starts, refer to PrepareNbodyComputation.C:219 for the counter part  by YS

*        particle id
      call MPI_RECV(EN, 1, MPI_INTEGER, 0, 100, ECOMM, istatus,
     &           ierr)
      allocate(EF(3,EN))
        DO I = 1,3 
        call MPI_RECV(EF(I,:), EN, MPI_DOUBLE_PRECISION, 0, 500,
     &            ECOMM, istatus,ierr)
        END DO
        call MPI_RECV(EDT, 1, MPI_DOUBLE_PRECISION, 0, 600,
     &            ECOMM, istatus,ierr)
        call MPI_RECV(ETU, 1, MPI_DOUBLE_PRECISION, 0, 700,
     &            ECOMM, istatus,ierr)
*       TCRIT = TCRIT + dtENZO
*     for SY, here TIMEU changes adaptivelyi on the fly?
        write (0,*) 'fortran: force=', EF(1,1)
      END IF
*     MPI done!

*      CALL ENZO_TO_NB(EX, EXDOT)
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

      deallocate(EF(3,EN))



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