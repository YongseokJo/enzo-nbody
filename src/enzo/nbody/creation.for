      SUBROUTINE CREATION(I,BODYNEW,XNEW,VNEW)
*
*
*       addition of newly created particles
*       ----------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tlist.h'
      INTEGER I,EN
      REAL*8 BODYNEW,XNEW(3),VNEW(3)
*
      write(6,*) 'creation starts'
      write(6,*) 'initial conditions are',BODYNEW,XNEW,VNEW

*       Update COMMON arrays to remove the escaper and correct F & FDOT.
      CALL ADD_STAR(I,1,BODYNEW,XNEW,VNEW)

      write(6,*) 'N = ',N
      write(6,*) 'NTOT = ',NTOT
      write(6,*) 'TIME = ',TIME
      write(6,*) 'NXTLIMIT = ',NXTLIMIT

*       increase particle number & total membership and check NNBMAX.
      N = N + 1
      NTOT = NTOT + 1
      NNBMAX = MIN(N/2,NNBMAX)
      ZNBMAX = 0.9*FLOAT(NNBMAX)
      ZNBMIN = MAX(0.01*FLOAT(NNBMAX),1.0)
*       Ensure NNBOPT is smaller than NNBMAX
      NNBOPT = MIN(NNBMAX-1,NNBOPT)
      NNBOPT = MAX(NNBOPT,1)

      write(6,*) 'creation routine ends'
*
      write(6,*) 'N = ',N
      write(6,*) 'NTOT = ',NTOT
      write(6,*) 'TIME = ',TIME
      write(6,*) 'NXTLIMIT = ',NXTLIMIT

      RETURN
*
      END
