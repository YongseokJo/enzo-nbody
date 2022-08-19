      SUBROUTINE FORCEINT(X,XDOT,BODY,SCALE,DT,NMAX,NXTLEN
     &     F,FD)
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

      INTEGER NMAX,NXTLEN
      REAL*8  F(3,NMAX),FD(3,NMAX),X(3,NMAX),XDOT(3,NMAX)
      REAL*8  BODY(NMAX),SCALE,DT

*     Local variables
*     A1-3 = relative position, DV(3) = relative velocity
*     RIJ2 = distance, DR2I = 1/r^2, DR3I = m/r^3, DRDV = r*v

      INTEGER I,J
      REAL*8  A1,A2,A3,DV(3),RIJ2,DR2I,DR3I,DRDV,DRDP

*     Initialize arrays, I covers all particles

      DO I = 1, NXTLEN         
         F(1:3,I) = 0.0D0
         FD(1:3,I) = 0.0D0
         RIJMIN = 1.E20

*     Calculate the forces
*     I am going to calculate all the particles,
*     so I do not need LIST (neighbor list)
*     J covers all interacting particles

         DO J = 1, NXTLEN

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
    
            F(1,II) = F(1,II) + A1*DR3I
            F(2,II) = F(2,II) + A2*DR3I
            F(3,II) = F(3,II) + A3*DR3I
            FD(1,II) = FD(1,II) + (DV(1) - A1*DRDP)*DR3I
            FD(2,II) = FD(2,II) + (DV(2) - A2*DRDP)*DR3I
            FD(3,II) = FD(3,II) + (DV(3) - A3*DRDP)*DR3I
            
         END DO
      END DO
      
      RETURN

      END
      