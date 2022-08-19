      SUBROUTINE FPRED_INT(I,XP,XPDOT,BODY,SCALE,TIME,DT,NMAX,NXTLEN
     &     FP,FPD)
*
*
*     Force Calculation of Particles in Integration w/o GPU
*     Modification complete on 2022.08.16
*     --------------------

*     Common variables (inputs and outputs)
*     
*     NMAX = max num of particles, NXTLEN = num of particles
*     X = position, XDOT = velocity
*     F = force, FD = force derivative, BODY = mass
*     SCALE = scale factor, TIME = current time, DT = time interval

      INTEGER I,NMAX,NXTLEN
      REAL*8  XP(3),XPDOT(3),FP(3),FPDOT(3),
      REAL*8  BODY(NMAX),SCALE,TIME,DT

*     Local variables
*     A1-3 = relative position, DV(3) = relative velocity
*     RIJ2 = distance, DR2I = 1/r^2, DR3I = m/r^3, DRDV = r*v

      INTEGER K
      REAL*8  A1,A2,A3,DV(3),RIJ2,DR2I,DR3I,DRDV,DRDP

*     Initialize arrays, I covers all particles

      DO K = 1, 3 

         FP(K) = 0.0D0
         FPD(K) = 0.0D0
      
      END DO

*     Relative position calculated (x2-x1)

      A1 = XP(1) - XP(1)
      A2 = XP(2) - XP(2)
      A3 = XP(3) - XP(3)

*     Predicted coordinates avoids spurious force differences.

      DV(1) = XPDOT(1) - XPDOT(1)
      DV(2) = XPDOT(2) - XPDOT(2)
      DV(3) = XPDOT(3) - XPDOT(3)

      RIJ2 = A1*A1 + A2*A2 + A3*A3

      DR2I = 1.0/RIJ2
      DR3I = BODY(I)*DR2I*SQRT(DR2I)
      DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
      DRDP = 3.0*DRDV*DR2I
    
      FP(1) = FP(1) + A1*DR3I
      FP(2) = FP(2) + A2*DR3I
      FP(3) = FP(3) + A3*DR3I
      FPD(1) = FPD(1) + (DV(1) - A1*DRDP)*DR3I
      FPD(2) = FPD(2) + (DV(2) - A2*DRDP)*DR3I
      FPD(3) = FPD(3) + (DV(3) - A3*DRDP)*DR3I
      
      RETURN

      END
      