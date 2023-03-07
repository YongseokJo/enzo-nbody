      PROGRAM LINKSTART
      
      INTEGER EN,REPEATN,NUM

      EN = 1000

      READ (5,*) REPEATN

      DO NUM = 1, REPEATN
          CALL LINKTEST(EN,NUM)
      END DO

      END
