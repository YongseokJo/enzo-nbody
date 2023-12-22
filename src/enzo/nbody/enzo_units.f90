      
*----------------------------------------------------------------------------------*
      SUBROUTINE  enzo_to_nb

        USE POINTERS
      INCLUDE 'common6.h'
      

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


        RETURN
        END

*----------------------------------------------------------------------------------*



*----------------------------------------------------------------------------------*
        SUBROUTINE  nb_to_enzo

          USE POINTERS
      INCLUDE 'common6.h'
      

*       need to make a new extrapolation scheme... in progress
*       after extrapolation,  update the EBODY, EX, EXDOT that would be passed onto ENZO


          RETURN
          END
*----------------------------------------------------------------------------------*
