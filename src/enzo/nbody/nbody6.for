      SUBROUTINE NBODY6(world_comm, inter_comm, nbody_comm, mpi_rank)
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
      integer :: inter_comm, ICOMM
      INTEGER :: world_comm, WCOMM
      LOGICAL :: FIRST=.true.

*     added by sykim, parameters from enzo

      INTEGER, parameter:: EN_MAX = 200000
      INTEGER :: EN,IE
      integer :: istatus(MPI_STATUS_SIZE)


      REAL*8, pointer :: EBODY(:),EX(:,:)
      REAL*8, pointer :: EXDOT(:,:), EF(:,:) !, EH(:,:,:)
      INTEGER, pointer :: EID(:)
*     initial ID of particles recieved from ENZO

      INTEGER EIDINIT(NMAX)
      INTEGER I
      INTEGER OUTC, OUTCNUM
      INTEGER SHUFFLECNT
      CHARACTER*27 OUTFILE
      CHARACTER*27 OUTFORM

      INTEGER IS,JS  ! loop for starting
      INTEGER IP,JP,KP  ! loop for predingting/sending
      INTEGER IR,JR,KR  ! loop for recieving force

*     conversion factors for enzo code unit -> cgs

      REAL*8 EMU,ELU,EVU,ETU


*     conversion factors for astronomical <-> nbody

      REAL*8 LENGTHU,MASSU,VELU,TIMEU,FORCEU

*     conversion factors for enzo <-> nbody

      REAL*8 ELENGTHU,EMASSU,EVELU,ETIMEU,EFORCEU
      REAL*8 EDT

*     end added by sykim

*     common variables for MPI by YS
      COMMON/MPI/ WCOMM, ICOMM, NCOMM, FIRST

      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
#ifdef DEBUG
*     --10/03/14 19:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      COMMON/adt/ adtime,dumptime,dprintt,dtprint,namep
*     --10/03/14 19:41-lwang-end----------------------------------------*
#endif
      EXTERNAL MERGE


#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif


*-----MPI-scheme-of-nbody-added-by-YS-----------------------------------*
      OUTC = 1

      write (0,*) 'In the Nbody6 '
*     Massive MPI Communication with Enzo code! (subroutine later) by SY
      write (0,*) 'Initialize Arrays'

*     Massive MPI Communication with Enzo code! (subroutine later) by SY
      write (0,*) 'Communication Starts'

      ! recieve the number of particles
*     MPI starts, refer to PrepareNbodyComputation.C for the counter prt
*     by YS
      WCOMM = world_comm
      ICOMM = inter_comm
      NCOMM = nbody_comm


      call MPI_RECV(EN, 1, MPI_INTEGER, 0, 100, ICOMM, istatus,
     &           ierr)

      write (0,*) 'fortran: Number of Nbody Particles on Fortran', EN
      allocate(EBODY(EN))
      allocate(EID(EN))
      call MPI_RECV(EBODY, EN, MPI_DOUBLE_PRECISION, 0, 200, ICOMM,
     &     istatus, ierr)
      call MPI_RECV(EID, EN, MPI_INTEGER, 0, 250, ICOMM,istatus, ierr)
      allocate(EX(3,EN))
      allocate(EXDOT(3,EN))
      allocate(EF(3,EN))
      DO  I = 1,3
        call MPI_RECV(EX(I,:), EN, MPI_DOUBLE_PRECISION, 0, 300,
     &           ICOMM, istatus,ierr)
        call MPI_RECV(EXDOT(I,:), EN, MPI_DOUBLE_PRECISION, 0, 400,
     &           ICOMM, istatus,ierr)
        call MPI_RECV(EF(I,:), EN, MPI_DOUBLE_PRECISION, 0, 500,
     &            ICOMM, istatus,ierr)
      END DO
      call MPI_RECV(EDT, 1, MPI_DOUBLE_PRECISION, 0, 600,
     &            ICOMM, istatus,ierr)
      call MPI_RECV(ETU, 1, MPI_DOUBLE_PRECISION, 0, 700,
     &            ICOMM, istatus,ierr)
      call MPI_RECV(ELU, 1, MPI_DOUBLE_PRECISION, 0, 800,
     &            ICOMM, istatus,ierr)
      call MPI_RECV(EMU, 1, MPI_DOUBLE_PRECISION, 0, 900,
     &            ICOMM, istatus,ierr)
      call MPI_RECV(EVU, 1, MPI_DOUBLE_PRECISION, 0, 1000,
     &            ICOMM, istatus,ierr)
*     MPI done!

*-----end-added-by-YS---------------------------------------------------*



*-----unit-conversions-added-by-sykim-----------------------------------*

*     conventional units are pc, Msun, km/s and Myr
*     the length unit is Rvir, mass :unit is total mass, et cetera

      write (6,*) 'before conversion',EBODY(1),EX(1,1),EXDOT(1,1)

      MASSU = 0.0D0

      DO I = 1,EN
          MASSU = MASSU + EBODY(I)*EMU/(1.9891D33)
      END DO

*     need to calculate virial radius and put that into LENGTHU0
*     ** code for calculating virial radius** <- should be added

*     calculates nbody <-> astro units

      LENGTHU = 2.58811228D0
      VELU = 6.557D0*((MASSU/LENGTHU)**(0.5D0))/(100.0D0)
      TIMEU = 14.94D0*(((LENGTHU)**3.0D0/MASSU)**(0.5D0))
*     added on 0424 - sykim
*      EFORCEU = MASSU*VELU*1.9891D33*1D5/(TIMEU*3.1556952D13)*
*     &         EMU*EVU**2/ELU


*     calculates enzo <-> nbody units

      EMASSU = (MASSU*1.9891D33)/EMU
      ELENGTHU = (LENGTHU*3.0857D18)/ELU
      EVELU = (VELU*1d5)/EVU
      ETIMEU = (TIMEU*3.1556952D13)/ETU
*     what was force unit recieved from ENZO again? -sykim
      EFORCEU = EMASSU*EVELU*EVELU/ELENGTHU


      write (6,*) 'scaling',LENGTHU,MASSU,VELU,TIMEU

*-----end-added-by-sykim--------------------------------------------*



*----parameter-setting-and-initial-value-recieving-sykim------------*

*     determine how much steps should be run depending on approximate
*     scaling
 
      N = EN
      TCRIT = 1.0D5
      DELTAT = EDT/ETIMEU

      write (6,*) 'timesteps',DELTAT


*     move the variable recieved from ENZO to nbody

      DO IS = 1,N
         
         EIDINIT(IS) = EID(IS)
         BODY(IS) = EBODY(IS)/EMASSU
         
         DO JS = 1,3
           X(JS,IS) = EX(JS,IS)/ELENGTHU
           XDOT(JS,IS) = EXDOT(JS,IS)/EVELU
           FENZO(JS,IS) = EF(JS,IS)/EFORCEU
         END DO
         !FENZO(1,IS) = 0.05D0
         !FENZO(2,IS) = 0.0D0
         !FENZO(3,IS) = 0.0D0

      END DO

      write (6,*) 'after conversion',BODY(1),X(1,1),XDOT(1,1)
      write (6,*) 'initial enzo id',EIDINIT(1),EIDINIT(2),EIDINIT(3)
*     hard-code  parameters needed for nbody6.F

      KSTART = 1
      TCOMP = 100000.0D0
      TCRITP = 1.0D8
      isernb = 40
      iserreg = 40
      iserks = 640

*    end added by sykim


*-----end-added-by-sykim--------------------------------------------*



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




*       Read input parameters, perform initial setup and obtain output.

      CPU = TCOMP
      CALL START
      write(0,*) "start finished - SY"
      call cputim(tt7)
      CALL ADJUST
      write(0,*) "adjust finished - SY"
      call cputim(tt8)
      ttadj = ttadj + (tt8-tt7)*60.

      TNEXT = DELTAT ! added by  sykim


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

      write(0,*) "go into integration loop"

      call cputim(tt1)
*
      CALL INTGRT
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

*----communication-with-ENZO-YS------------------------------------------*
      IF ((EPHASE.EQ.2).OR.(EPHASE.EQ.3)) THEN

      SHUFFLECNT = 0

      DO 29 IP = 1,NTOT

         IF (NAME(IP).GT.N) GO TO 23 
         IE = EIDINIT(NAME(IP)) ! the id of the ipth particle 

         DO JP = 1,N
            IF (IE.EQ.EID(JP)) THEN
               EBODY(JP) = BODY(IP)*EMASSU
               DO KP = 1,3
                  EX(KP,JP) = X(KP,IP)*ELENGTHU
                  EXDOT(KP,JP) = XDOT(KP,IP)*EVELU
               END DO
            if (JP.LE.3) write(0,*) 'enzo id matches!'
            SHUFFLECNT = SHUFFLECNT + 1
            END IF
         END DO

   23   CONTINUE
   29   CONTINUE
         
          write (0,*) 'fortran: NTOT',NTOT
          write (0,*) 'fortran: X=', X(1,1), ', V=',XDOT(1,1)
          write (0,*) 'fortran: RDENS=',RDENS(1),RDENS(2),RDENS(3)
          write (0,*) '# of matched shuffles', SHUFFLECNT

*---------------------------------------------------------------------*
*    send particles to ENZO

*     MPI starts, refer to FinalizeNbodyComputation.C for the counter
*     part  by YS
          write (0,*) 'fortran: write out'
          write (6,*) 'fortran: write out'
          write (0,*) 'fortran: FN =', N

          call MPI_BARRIER(ICOMM, ierr)
          call MPI_SEND(EID, EN, MPI_INTEGER, 0, 250,
     &           ICOMM, ierr)
          DO  I = 1,3
          call MPI_SEND(EX(I,:), EN, MPI_DOUBLE_PRECISION, 0, 300,
     &           ICOMM, ierr)
          call MPI_SEND(EXDOT(I,:), EN, MPI_DOUBLE_PRECISION, 0, 400,
     &           ICOMM, ierr)
          END DO
          call MPI_BARRIER(ICOMM,ierr)
*          deallocate(EF)
          write (0,*) 'fortran: EX=', EX(1,1), ', EV=',EXDOT(1,1)
*     MPI done!

*----------------------------------------------------------------------*
          ! RECEIVE
*     MPI starts, refer to PrepareNbodyComputation.C:219 for the counter
*     part  by YS

*        particle id
*        force unit (likely cgs)
*          call MPI_RECV(EN, 1, MPI_INTEGER, 0, 100, ICOMM, istatus,
*     &         ierr)
          !allocate(EF(3,EN))
          write (0,*) 'fortran: read in'
          write (6,*) 'fortran: read in'

          call MPI_BARRIER(ICOMM,ierr)
          call MPI_RECV(EID, EN, MPI_INTEGER, 0, 250, ICOMM, istatus,
     &         ierr)
            DO I = 1,3
               call MPI_RECV(EF(I,:), EN, MPI_DOUBLE_PRECISION, 0, 500,
     &         ICOMM, istatus,ierr)
            END DO

          call MPI_RECV(EDT, 1, MPI_DOUBLE_PRECISION, 0, 600,
     &            ICOMM, istatus,ierr)
          call MPI_RECV(ETU, 1, MPI_DOUBLE_PRECISION, 0, 700,
     &            ICOMM, istatus,ierr)
          call MPI_BARRIER(ICOMM,ierr)
*       TCRIT = TCRIT + dtENZO
*     for SY, here TIMEU changes adaptivelyi on the fly?
          write (0,*) 'fortran: force=', EF(1,1)
*     MPI done!

*      Recieve forces from ENZO and update them
       
       DO IR = 1,NTOT

          IF (NAME(IR).GT.EN) GO TO 33
          IE = EIDINIT(NAME(IR)) 

          DO JR = 1,EN
             IF (EID(JR).EQ.IE) THEN
             DO KR = 1,3
                FENZO(KR,IR)= EF(KR,JR)/EFORCEU
             END DO
             END IF
          END DO

   33  CONTINUE

       END DO

          DELTAT = EDT/ETIMEU
          TNEXT = TNEXT + DELTAT  ! is it right? SY
          
          write (6,*) 'timesteps', DELTAT

          write (6,*) 'recieved and restarting'
         


        ! output
        OUTCNUM = INT(LOG10(REAL(OUTC)))+1
        write(OUTFORM,"(A5,I1,A4)") '(A7,I',OUTCNUM,',A4)'
        write(OUTFILE,OUTFORM) 'output_',OUTC
        OPEN (UNIT=3,STATUS='UNKNOWN',
     &         FILE=OUTFILE)

        OUTCNUM = NAME(J)
        DO J = 1, NTOT

          WRITE (3,
     &         '(f20.8,f20.8,f20.8,f20.8,f20.8,
     &          f20.8,f20.8,f20.8,f20.8,I5)')  
     &         (X(K,J),K=1,3), (XDOT(K,J),K=1,3),
     &         (FENZO(K,J),K=1,3), OUTCNUM
        END DO

        CLOSE(3)
      OUTC = OUTC + 1




*----end-added-by-YS-----------------------------------------------------*
      END IF


*       Continue integration.
      GO TO 1
*
      END
