      SUBROUTINE FPRED_INT(SCALE,DT,NXTLEN)
*
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
         F(1,I) = 0.0D0
         F(2,I) = 0.0D0
         F(3,I) = 0.0D0
         FDOT(1,I) = 0.0D0
         FDOT(2,I) = 0.0D0
         FDOT(3,I) = 0.0D0
         RIJMIN = 1.E20

*     Calculate the forces
*     I am going to calculate all the particles,
*     so I do not need LIST (neighbor list)
*     J covers all interacting particles

         DO J = 1, NXTLEN
            
            IF (J .NE. I) THEN

*     Relative position calculated (x2-x1)
               A1 = X(1,J) - X(1,I)
               A2 = X(2,J) - X(2,I)
               A3 = X(3,J) - X(3,I)

*     Predicted coordinates avoids spurious force differences.

               DV(1) = XDOT(1,J) - XDOT(1,I)
               DV(2) = XDOT(2,J) - XDOT(2,I)
               DV(3) = XDOT(3,J) - XDOT(3,I)

               RIJ2 = A1*A1 + A2*A2 + A3*A3

               IF(RIJ2.LT.RIJMIN) THEN
                  RIJMIN = RIJ2
               END IF

               DR2I = 1.0/RIJ2
               DR3I = BODY(J)*DR2I*SQRT(DR2I)
               DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
               DRDP = 3.0*DRDV*DR2I
    
               F(1,I) = F(1,I) + A1*DR3I
               F(2,I) = F(2,I) + A2*DR3I
               F(3,I) = F(3,I) + A3*DR3I
               FDOT(1,I) = FDOT(1,I) + (DV(1) - A1*DRDP)*DR3I
               FDOT(2,I) = FDOT(2,I) + (DV(2) - A2*DRDP)*DR3I
               FDOT(3,I) = FDOT(3,I) + (DV(3) - A3*DRDP)*DR3I
            END IF
         END DO
*         PRINT *,'FORCE',F(1,I),FDOT(1,I)

      END DO
      
      RETURN

      END