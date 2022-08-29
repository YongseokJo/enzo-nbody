      SUBROUTINE FORCEINT(SCALE,DT,NXTLEN)
*
*
*     Force Calculation of Particles w/o GPU
*     Modification complete on 2022.08.16
*     --------------------

*     Common variables (inputs and outputs)
*     
*     NMAX = max num of particles, NXTLEN = num of particles
*     X = position, XDOT = velocity
*     F = force, FD = force derivative, BODY = mass
*     SCALE = scale factor, TIME = current time, DT = time interval
      
      INCLUDE 'enzo_common6.h'
      INTEGER NXTLEN
      REAL*8  SCALE,DT

*     Local variables
*     A1-3 = relative position, DV(3) = relative velocity
*     RIJ2 = distance, DR2I = 1/r^2, DR3I = m/r^3, DRDV = r*v

      INTEGER I,J
      REAL*8  A1,A2,A3,DV(3),RIJ2,RIJMIN,DR2I,DR3I,DRDV,DRDP

*     Initialize arrays, I covers all particles

      DO I = 1, NXTLEN         
         F0(1,I) = 0.0D0
         F0(2,I) = 0.0D0
         F0(3,I) = 0.0D0
         F0DOT(1,I) = 0.0D0
         F0DOT(2,I) = 0.0D0
         F0DOT(3,I) = 0.0D0
         RIJMIN = 1.E20

*     Calculate the forces
*     I am going to calculate all the particles,
*     so I do not need LIST (neighbor list)
*     J covers all interacting particles

         DO J = 1, NXTLEN
            
            IF (J .NE. I) THEN

*     Relative position calculated (x2-x1)
               A1 = X0(1,J) - X0(1,I)
               A2 = X0(2,J) - X0(2,I)
               A3 = X0(3,J) - X0(3,I)

*     Predicted coordinates avoids spurious force differences.

               DV(1) = X0DOT(1,J) - X0DOT(1,I)
               DV(2) = X0DOT(2,J) - X0DOT(2,I)
               DV(3) = X0DOT(3,J) - X0DOT(3,I)

               RIJ2 = A1*A1 + A2*A2 + A3*A3

               IF(RIJ2.LT.RIJMIN) THEN
                  RIJMIN = RIJ2
               END IF

               DR2I = 1.0/RIJ2
               DR3I = BODY(J)*DR2I*SQRT(DR2I)
               DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
               DRDP = 3.0*DRDV*DR2I
    
               F0(1,I) = F0(1,I) + A1*DR3I
               F0(2,I) = F0(2,I) + A2*DR3I
               F0(3,I) = F0(3,I) + A3*DR3I
               F0DOT(1,I) = F0DOT(1,I) + (DV(1) - A1*DRDP)*DR3I
               F0DOT(2,I) = F0DOT(2,I) + (DV(2) - A2*DRDP)*DR3I
               F0DOT(3,I) = F0DOT(3,I) + (DV(3) - A3*DRDP)*DR3I

            END IF
         END DO
      END DO
      
      RETURN

      END
      