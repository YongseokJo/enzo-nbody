
!***********************************************************************


    SUBROUTINE ellan


!***********************************************************************


!     Subroutine to analyze the ellipticity of the system


!=======================================================================

    USE POINTERS
    INCLUDE 'common6.h'
          
    COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
!   Declaration of local variables.
!   -------------------------------
    DOUBLE PRECISION :: etoti(1:nmax),mtot,poscm(3), &
    ti(3,3),tiwork(3,3),dwork(3),ework(3),lam(3)
    INTEGER :: i,j,k,ief,nbound,nstart,nnext,np,indexev(3)
    INTEGER :: index(1:nmax)

    INTEGER :: nef,nef1
    PARAMETER(nef=9,nef1=nef+1)
    DOUBLE PRECISION :: xf(nef1),ba(nef1),ca(nef1),taue(nef1), &
    evec(3,3,nef1)

    PARAMETER (tiny2=1.D-30)
    DATA (xf(i),i=1,nef1) /0.01D0,0.05D0,0.1D0,0.2D0,0.5D0, &
    0.8D0,0.9D0,0.95D0,1.0D0,99.D0/

!=======================================================================
            
!       calculate specific energy of particles
!       --------------------------------------

    DO 100 i=ifirst,ntot
        etoti(i-ifirst+1) = 0.5D0 * (xdot(1,i)**2 + xdot(2,i)**2 + &
        xdot(3,i)**2) - phidbl(i)
    100 END DO

!      calculate number of bound particles
!      -----------------------------------

    nbound = 0
    DO 150 i=1,ntot-ifirst+1
        IF(etoti(i) < 0.D0) nbound = nbound + 1
    150 END DO

!       sort for particle energy
!       ------------------------

    CALL indexx(ntot-ifirst+1,etoti,index)

!       initialize tensor of inertia
!       ----------------------------

    DO 210 i=1,3
        DO 200 k=1,3
            ti(i,k) = 0.D0
        200 END DO
    210 END DO

!       LOOP over fraction of most bound particles and all particles
!       ------------------------------------------------------------

    nstart   = 1
    mtot     = 0.D0
    poscm(1) = 0.D0
    poscm(2) = 0.D0
    poscm(3) = 0.D0
    DO 500 ief=1,nef1

        IF(ief <= nef) THEN
        !                                  only fraction of bound particles
        !                                  --------------------------------
            nnext = NINT(xf(ief) * nbound)
        ELSE
        !                                   all particles
        !                                   -------------
            nnext = ntot
        ENDIF

    !-----------------------------------------------------------------
    !--      at least two particles are required for ellipticity...
    !-----------------------------------------------------------------
        IF(nnext < 2) THEN
            ba(ief)  = 999.
            ca(ief)  = 999.
            taue(ief) = 999.
            DO 320 k=1,3
                DO 310 j=1,3
                    evec(k,j,ief) = 0.
                310 END DO
            320 END DO

        ELSE

        !       calculate tensor of inertia
        !       ---------------------------
        !             print*,' nstart,nnext=',nstart,nnext
            DO 400 i=nstart,nnext
                ipo = index(i) + ifirst - 1
                ti(1,1) = ti(1,1) + body(ipo) * &
                (x(2,ipo)*x(2,ipo) + x(3,ipo)*x(3,ipo))
                ti(2,2) = ti(2,2) + body(ipo) * &
                (x(1,ipo)*x(1,ipo) + x(3,ipo)*x(3,ipo))
                ti(3,3) = ti(3,3) + body(ipo) * &
                (x(1,ipo)*x(1,ipo) + x(2,ipo)*x(2,ipo))
                ti(1,2) = ti(1,2) - body(ipo) * x(1,ipo)*x(2,ipo)
                ti(1,3) = ti(1,3) - body(ipo) * x(1,ipo)*x(3,ipo)
                ti(2,3) = ti(2,3) - body(ipo) * x(2,ipo)*x(3,ipo)
            400 END DO
                  
        !       correct for center of mass
        !       --------------------------

        !   A) calculate center of mass data
        !   --------------------------------

        !--       remove previous correction for center of mass
            ti(1,1) = ti(1,1) + mtot * (poscm(2)**2+poscm(3)**2)
            ti(2,2) = ti(2,2) + mtot * (poscm(1)**2+poscm(3)**2)
            ti(3,3) = ti(3,3) + mtot * (poscm(1)**2+poscm(2)**2)
            ti(1,2) = ti(1,2) - mtot * poscm(1) * poscm(2)
            ti(1,3) = ti(1,3) - mtot * poscm(1) * poscm(3)
            ti(2,3) = ti(2,3) - mtot * poscm(2) * poscm(3)
            poscm(1) = poscm(1) * mtot
            poscm(2) = poscm(2) * mtot
            poscm(3) = poscm(3) * mtot
        !             xav = 0.d0
        !             yav = 0.d0
        !             zav = 0.d0

            DO 405 i=nstart,nnext
                ipo = index(i) + ifirst - 1
                poscm(1) = poscm(1) + body(ipo) * x(1,ipo)
                poscm(2) = poscm(2) + body(ipo) * x(2,ipo)
                poscm(3) = poscm(3) + body(ipo) * x(3,ipo)
                mtot     = mtot + body(ipo)
            !                xav = xav + abs(x(1,ipo))
            !                yav = yav + abs(x(2,ipo))
            !                zav = zav + abs(x(3,ipo))

            405 END DO
            poscm(1) = poscm(1) / mtot
            poscm(2) = poscm(2) / mtot
            poscm(3) = poscm(3) / mtot
        !             print*,' av=',xav,yav,zav

             
        !   B) apply correction
        !   -------------------
            ti(1,1) = ti(1,1) - mtot * (poscm(2)**2+poscm(3)**2)
            ti(2,2) = ti(2,2) - mtot * (poscm(1)**2+poscm(3)**2)
            ti(3,3) = ti(3,3) - mtot * (poscm(1)**2+poscm(2)**2)
            ti(1,2) = ti(1,2) + mtot * poscm(1) * poscm(2)
            ti(1,3) = ti(1,3) + mtot * poscm(1) * poscm(3)
            ti(2,3) = ti(2,3) + mtot * poscm(2) * poscm(3)

        !       set off-axis values by symmetry
        !       -------------------------------

            ti(2,1) = ti(1,2)
            ti(3,1) = ti(1,3)
            ti(3,2) = ti(2,3)
        
        !             print*,' mtot,poscm=',mtot,(poscm(k),k=1,3)
        !=======================================================
        !       determine eigenvalues and axis of inertia
        !=======================================================

        !------------------------------------------------------
        !--          copy tensor of inertia
        !------------------------------------------------------
            DO 420 i=1,3
                DO 410 k=1,3
                    tiwork(i,k) = ti(i,k)
                410 END DO
            420 END DO
            np = 3

        !------------------------------------------------------
        !--          calculate eigenvalues and eigenvectors
        !------------------------------------------------------
            CALL tred2(tiwork,np,np,dwork,ework)
            CALL tqli(dwork,ework,np,np,tiwork)

        !--               sort for increasing eigenvalues
            CALL indexx(np,dwork,indexev)
        !--               find eigenvectors
            DO 450 i=1,np
                lam(i) = dwork(indexev(i))
                DO 440 k=1,np
                    evec(k,i,ief) = tiwork(k,indexev(i))
                440 END DO
            450 END DO

            xhelp    = lam(3) + lam(2) - lam(1)
            xhelp1   = lam(2) - lam(3) + lam(1)
        !              IF(xhelp1.LT.0.D0) THEN
        !                 PRINT*,' ellan: xhelp1 < 0',xhelp1,tnow
        !                 xhelp1 = 0.D0
        !             ENDIF
            ba(ief)  = SQRT(MAX(tiny2,lam(3)-lam(2)+lam(1)) / xhelp)
            ca(ief)  = SQRT(MAX(tiny2,xhelp1) / xhelp)
            taue(ief) = (ba(ief)-ca(ief)) /MAX(tiny2,(1.D0 - ca(ief)))

            nstart = nnext + 1

        ENDIF

    500 END DO

!==================================================================
!==         OUTPUT of data
!===================================================================

!        DO 600 ief=1,nef1
!          if(rank.eq.0)
!     *     WRITE(60,910) ief, time,ba(ief),ca(ief),taue(ief),
!     *     mtot,(poscm(k),k=1,3)
! 00     CONTINUE

!     900     FORMAT(I4,1x,1p,e12.5,1x,0p,i9,1x,i9,1x,i2)
! 10     FORMAT(I4,1x,1P,8E12.5,0P)

    RETURN
            
    end SUBROUTINE ellan

