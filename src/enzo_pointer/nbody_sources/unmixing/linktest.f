      SUBROUTINE LINKTEST(EN)
     
      include 'common6.h'
 
      INTEGER I,J,EN
      
      REAL*8 EBODY(EN),EX1(EN),EX2(EN),EX3(EN)
      REAL*8 EXDOT1(EN),EXDOT2(EN),EXDOT3(EN)
      REAL*8 EF1(EN),EF2(EN),EF3(EN)

      REAL*8 EH11(EN),EH12(EN),EH13(EN)
      REAL*8 EH21(EN),EH22(EN),EH23(EN)
      REAL*8 EH31(EN),EH32(EN),EH33(EN)
      REAL*8 EH41(EN),EH42(EN),EH43(EN)

      REAL*8 EDT


*     read data from initial file

      OPEN (UNIT=10,FILE='dat.10')
   
      DO 5 I = 1,EN

          READ (10,*)  EBODY(I),EX1(I),EX2(I),EX3(I),
     &                 EXDOT1(I),EXDOT2(I),EXDOT3(I)

    5 CONTINUE


      DO 7 J = 1,EN

           EF1(J) = 0.0D0
           EF2(J) = 0.0D0
           EF3(J) = 0.0D0
           
           EH11(J) = 0.0D0
           EH12(J) = 0.0D0
           EH13(J) = 0.0D0
           EH21(J) = 0.0D0
           EH22(J) = 0.0D0
           EH23(J) = 0.0D0
           EH31(J) = 0.0D0
           EH32(J) = 0.0D0
           EH33(J) = 0.0D0
           EH41(J) = 0.0D0
           EH42(J) = 0.0D0
           EH43(J) = 0.0D0

    7 CONTINUE

*     the time in yr that the nbody simulation will run

      EDT = 311.0047

      CALL NBODY6(EN,EBODY,EX1,EX2,EX3,
     &            EXDOT1,EXDOT2,EXDOT3,
     &            EF1,EF2,EF3,
     &            EH11,EH12,EH13,
     &            EH21,EH22,EH23,
     &            EH31,EH32,EH33,
     &            EH41,EH42,EH43,EDT)     

      END
