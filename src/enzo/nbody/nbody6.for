      SUBROUTINE NBODY6(EN,EBODY,EX1,EX2,EX3,
     &            EXDOT1,EXDOT2,EXDOT3,
     &            EF1,EF2,EF3,
     &            EH11,EH12,EH13,
     &            EH21,EH22,EH23,
     &            EH31,EH32,EH33,
     &            EH41,EH42,EH43,EDT,
     &            EMU,ELU,EVU,ETU)
 
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

      INTEGER EN,IS,ENDSTEP,IE,EID,I

      REAL*8 EBODY(EN),EX1(EN),EX2(EN),EX3(EN)
      REAL*8 EXDOT1(EN),EXDOT2(EN),EXDOT3(EN)
      REAL*8 EF1(EN),EF2(EN),EF3(EN)

      REAL*8 EH11(EN),EH12(EN),EH13(EN)
      REAL*8 EH21(EN),EH22(EN),EH23(EN)
      REAL*8 EH31(EN),EH32(EN),EH33(EN)
      REAL*8 EH41(EN),EH42(EN),EH43(EN)

*     conversion factors for enzo code unit -> cgs

      REAL*8 EMU,ELU,EVU,ETU


*     conversion factors for astronomical -> nbody

      REAL*8 LENGTHU0,MASSU0,VELU0,TIMEU0
      REAL*8 LENGTHU,MASSU,VELU,TIMEU
      
      REAL*8 EDT

      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
#ifdef DEBUG
*     --10/03/14 19:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      COMMON/adt/ adtime,dumptime,dprintt,dtprint,namep
*     --10/03/14 19:41-lwang-end----------------------------------------*
#endif
      EXTERNAL MERGE
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

* recieve input from enzo and place it at necessary variables


*     conventional units are pc, Msun, km/s and Myr
*     the length unit is Rvir, mass :unit is total mass, et cetera
      
      write (6,*) 'before conversion',EBODY(1),EX1(1),EXDOT1(1)      

      MASSU0 = 0.0D0

      DO I = 1,EN
          MASSU0 = MASSU0 + EBODY(I)*EMU/(1.9891e33)
      END DO
*     need to calculate virial radius and put that into LENGTHU0
*     ** code for calculating virial radius** <- should be added
      LENGTHU0 = 2.58811
      VELU0 = 6.557*((MASSU0/LENGTHU0)**(0.5))/(100)
      TIMEU0 = 14.94*(((LENGTHU0)**3/MASSU0)**(0.5))
       
      write (6,*) 'initial scale',LENGTHU0,MASSU0,VELU0,TIMEU0

*     determine how much steps should be run depending on approximate
*     scaling

      ENDSTEP = INT(EDT*ETU/(TIMEU0*(3.1557E13)))+1
      TCRIT = REAL(ENDSTEP)
      
      write (6,*) 'timesteps',ENDSTEP,TCRIT

*     convert scaling according to new values

      TIMEU = EDT*ETU/(TCRIT*3.1557E13)
      MASSU = MASSU0
      LENGTHU = LENGTHU0*((TIMEU/TIMEU0)**(2.0/3.0))
      VELU = 6.557*((MASSU/LENGTHU)**(0.5))/(100)

      write (6,*) 'newscale',LENGTHU, MASSU, VELU, TIMEU       

      N = EN
      
      DO 7 IS = 1,N

         BODY(IS) = EBODY(IS)*EMU/(MASSU*1.9891e33)
         X(1,IS) = EX1(IS)*ELU/LENGTHU/(3.0857E18)
         X(2,IS) = EX2(IS)*ELU/LENGTHU/(3.0857E18)
         X(3,IS) = EX3(IS)*ELU/LENGTHU/(3.0857E18)
         XDOT(1,IS) = EXDOT1(IS)*EVU/VELU/(1e5)
         XDOT(2,IS) = EXDOT2(IS)*EVU/VELU/(1e5)
         XDOT(3,IS) = EXDOT3(IS)*EVU/VELU/(1e5)

    7 CONTINUE


      write (6,*) 'after conversion',BODY(1),X(1,1),XDOT(1,1)

*input needed for nbody6.F

      KSTART = 1
      TCOMP = 100000.0
      TCRITp = 1.E8
      isernb = 40
      iserreg = 40
      iserks = 640



#if MPIINIT
*       Initialize MPI
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,group,ierr)
      CALL MPI_GROUP_SIZE(group,isize,ierr)
      CALL MPI_GROUP_RANK(group,rank,ierr)
      ibarcount=0
*      write(6,11) rank,isize,group
* 11   format('MPI-initial: This is rank=',I6,' size=',I6,' group=',I11)
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
*      call flush(6)
*
*       Read start/restart indicator & CPU time.
*      IF(rank.eq.0)READ (5,*)  KSTART, TCOMP, TCRITP,
*     *    isernb,iserreg,iserks


#ifdef DEBUG
*     --10/03/14 19:41-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      if(rank.eq.0) read(5,*) adtime,dumptime,dprintt,dtprint,namep
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
*        PRINT*,' iserreg,isernb,iserks,ithread=',iserreg,isernb,iserks,
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
*          if(rank.eq.0)then
*             WRITE (6,10) TTOT/TCR0, TIME/TCR0, TCRIT/TCR0, TTOT, TIME,
*     &            TCRIT
*             WRITE (6,20)  DTADJ, DELTAT, TADJ, TNEXT, TCRIT, QE
*             WRITE (6,30)  ETAI, ETAR, ETAU, DTMIN, RMIN, NNBOPT
* 10          FORMAT (' START AT TTOT/TIME ',2F16.8,' STOP INTENDED AT ',
*     &            F16.8,' TCR0',/,' START AT TTOT/TIME ',2F16.8,
*     &            ' STOP INTENDED AT ',F16.8,' NBODY-UNITS ',/)
* 20          FORMAT (/,7X,'RESTART PARAMETERS:   DTADJ =',F7.3,
*     &            '  DELTAT =',F7.3,'   TADJ =',F7.3,'   TNEXT =',
*     &            F7.3,'  TCRIT =',F7.1,'  QE =',1PE9.1)
* 30          FORMAT (/,7X,'                      ETAI =',F7.3,
*     &            '  ETAR =',F7.3,'  ETAU =',F7.3,'  DTMIN =',1PE9.1,
*     &            '  RMIN =',E9.1,' NNBOPT =',I5,/)
*          end if
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
      DO 17 EID = 1,N
         
         IE = NAME(EID)
         EBODY(IE) = BODY(IE)*MASSU*1.9891e33/EMU
         EX1(IE) = X(1,IE)*LENGTHU*3.0857E18/ELU
         EX2(IE) = X(2,IE)*LENGTHU*3.0857E18/ELU
         EX3(IE) = X(3,IE)*LENGTHU*3.0857E18/ELU
         EXDOT1(IE) = XDOT(1,IE)*VELU*1e5/EVU
         EXDOT2(IE) = XDOT(2,IE)*VELU*1e5/EVU
         EXDOT3(IE) = XDOT(3,IE)*VELU*1e5/EVU

   17 CONTINUE

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
*
      ELSE IF (IPHASE.EQ.3) THEN
*       Perform energy check & parameter adjustments and print diagnostics.
          call cputim(tt7)
          CALL ADJUST
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
*       Continue integration
      GO TO 1
     
      END
