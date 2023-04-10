
      SUBROUTINE ENZO_COMMUNICATION
*
*
*       ENZO COMMUNICATION
*       ---------------------
*
      INCLUDE 'common6.h'



*----------------------------------------------------------------------------------*
      ! SEND
      allocate(EX(3,N))
      allocate(EDOT(3,N))
      CALL NB_TO_ENZO(EX, EDOT)
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

      CALL ENZO_TO_NB(EX, EXDOT)
*----------------------------------------------------------------------------------*





      deallocate(EF(3,EN))
*
      RETURN
*
      END
