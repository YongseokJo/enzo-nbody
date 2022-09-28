      SUBROUTINE ZERO
*
*
*     Initialization / Defining of global scalars.
*     Modification complete on 2022.08.17
*     ---------------------------------
*
      INCLUDE 'enzo_common6.h'

*
*       Set fractional constants & two PI.

      ONE3 = 1.0D0/3.0D0
      ONE6 = 1.0D0/6.0D0
      ONE9 = 1.0D0/9.0D0
      ONE12 = 1.0D0/12.0D0
      TWOPI = 8.0D0*ATAN(1.0D0)
      TINY = 1.0D-14

      RETURN
*
      END
