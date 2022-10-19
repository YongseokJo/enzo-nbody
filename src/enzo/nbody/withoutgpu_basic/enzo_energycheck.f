      SUBROUTINE ENERGYCHECK(SCALE,NXTLEN,
     &           KE,PE,ME)

*
*     Energy check of Particles w/o GPU
*     Modification complete on 2022.08.31
*     --------------------
*     Common variables (inputs and outputs)

      INCLUDE 'enzo_common6.h'

      INTEGER I,J,K,NXTLEN
      REAL*8 KE(NMAX),PE(NMAX),ME(NMAX)
      REAL*8 SCALE,RIJ2,DR2I,V2
      REAL*8 A1,A2,A3

      ALLKE = 0.0D0
      ALLPE = 0.0D0
      ALLME = 0.0D0

      DO I = 1, NXTLEN
      
         KE(I) = 0.0D0
         PE(I) = 0.0D0
         ME(I) = 0.0D0

*     Calculate the forces
*     I am going to calculate all the particles,
*     so I do not need LIST (neighbor list)
*     J covers all interacting particles

         DO J = 1, NXTLEN
            
            IF (J .LT. I) THEN

*     Relative position calculated (x2-x1)

               A1 = XF(1,J) - XF(1,I)
               A2 = XF(2,J) - XF(2,I)
               A3 = XF(3,J) - XF(3,I)

*     Calculate energy

               RIJ2 = A1*A1 + A2*A2 + A3*A3

               DR2I = 1.0/RIJ2

               PE(I) = PE(I) - G*BODY(I)*BODY(J)/SQRT(RIJ2)               

            END IF
         END DO
      END DO

      DO K = 1, NXTLEN

         V2 = XFDOT(1,K)*XFDOT(1,K)+XFDOT(2,K)*XFDOT(2,K)
     &   +XFDOT(3,K)*XFDOT(3,K)

         KE(K) = V2*0.5*BODY(K)
            
         ME(K) = KE(K) + PE(K)

         ALLKE = ALLKE + KE(K)
         ALLPE = ALLPE + PE(K)
         ALLME = ALLME + KE(K) + PE(K)

      END DO
      
      RETURN

      END
      