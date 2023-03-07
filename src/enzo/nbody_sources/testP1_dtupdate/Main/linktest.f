      SUBROUTINE LINKTEST(EN,NUM)
     
      include 'common6.h'
 
      INTEGER I,J,J2,K,EN

      INTEGER NUM,NEXTNUM
      INTEGER INTNUM,INTNEXTNUM

      REAL*8 LOGNUM,LOGNEXTNUM
      CHARACTER (len=20) OLDFILE,NEWFILE
      CHARACTER (len=20) OLDFORM,NEWFORM
      
      REAL*8 EBODY(EN),EX1(EN),EX2(EN),EX3(EN)
      REAL*8 EXDOT1(EN),EXDOT2(EN),EXDOT3(EN)
      REAL*8 EF1(EN),EF2(EN),EF3(EN)

      REAL*8 EH11(EN),EH12(EN),EH13(EN)
      REAL*8 EH21(EN),EH22(EN),EH23(EN)
      REAL*8 EH31(EN),EH32(EN),EH33(EN)
      REAL*8 EH41(EN),EH42(EN),EH43(EN)

      REAL*8 EDT

      REAL*8 ELU,EVU,EMU,ETU

*     enzo code units. set to 1 for nbody-only runs

      ELU = 1.0D0
      EVU = 1.0D0
      EMU = 1.0D0
      ETU = 1.0D0


*     make appropriate file names

      NEXTNUM = NUM+1

      LOGNUM = LOG10(REAL(NUM))
      INTNUM = INT(LOGNUM)+1

      LOGNEXTNUM = LOG10(REAL(NEXTNUM))
      INTNEXTNUM = INT(LOGNEXTNUM)+1

      write (OLDFORM,"(A5,I1,A4)") '(A3,I',INTNUM,',A4)'
      write (NEWFORM,"(A5,I1,A4)") '(A3,I',INTNEXTNUM,',A4)'

      write (OLDFILE,OLDFORM) "out",NUM,'.dat'
      write (NEWFILE,NEWFORM) "out",NEXTNUM,'.dat'


*     open and read files

      OPEN (UNIT=10,STATUS ="UNKNOWN",FILE=OLDFILE)
      OPEN (UNIT=33,FILE='dt.dat')
   
      DO 5 I = 1,EN

          READ (10,*)  EBODY(I),EX1(I),EX2(I),EX3(I),
     &                 EXDOT1(I),EXDOT2(I),EXDOT3(I)

    5 CONTINUE

      READ (33,*) EDT

      EDT = EDT*3.1556952D13


*     input = in conventional astronomical units
*     Msun, pc, km/s, Myr



      DO 7 J = 1,EN

           EBODY(J) = EBODY(J)*1.9891D33

           EX1(J) = EX1(J)*3.0857D18
           EX2(J) = EX2(J)*3.0857D18
           EX3(J) = EX3(J)*3.0857D18
           EXDOT1(J) = EXDOT1(J)*1D5
           EXDOT2(J) = EXDOT2(J)*1D5
           EXDOT3(J) = EXDOT3(J)*1D5

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



*     the time in seconds that the nbody simulation will run


      CALL NBODY6(EN,EBODY,EX1,EX2,EX3,
     &            EXDOT1,EXDOT2,EXDOT3,
     &            EF1,EF2,EF3,
     &            EH11,EH12,EH13,
     &            EH21,EH22,EH23,
     &            EH31,EH32,EH33,
     &            EH41,EH42,EH43,EDT,
     &            EMU,ELU,EVU,ETU)     

       open (UNIT=17,STATUS ="UNKNOWN",FILE=NEWFILE)

      DO 9 J2 = 1,EN

           EBODY(J2) = EBODY(J2)/1.9891D33
           EX1(J2) = EX1(J2)/3.0857D18
           EX2(J2) = EX2(J2)/3.0857D18
           EX3(J2) = EX3(J2)/3.0857D18
           EXDOT1(J2) = EXDOT1(J2)/1D5
           EXDOT2(J2) = EXDOT2(J2)/1D5
           EXDOT3(J2) = EXDOT3(J2)/1D5

    9 CONTINUE


       DO K = 1, EN
            WRITE (17,'(7f20.8)')  EBODY(K),
     &      EX1(K),EX2(K),EX3(K),EXDOT1(K),
     &      EXDOT2(K),EXDOT3(K)
       END DO

      write(6,*) 'writing first output'

      close(10)
      close(17)
      close(33)

      RETURN

      END
