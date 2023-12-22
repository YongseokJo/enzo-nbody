!     general three-body stability algorithm

!     system is unstable if nstab=1 is returned
!     system is stable if nstab=0 is returned

!     Rosemary Mardling
!     School of Mathematical Sciences, Monash University

!     version as of 16-11-07
!     email mardling@sci.monash.edu.au to be added to updates email list
!     preprint on astro-ph by New Year :-)

!     sigma=period ratio (outer/inner) **should be > 1**
!     ai=inner semi-major axis
!     a0=outer semi-major axis
!     ei0=initial inner eccentricity
!     eo=outer eccentricity
!     relinc=relative inclination (radians)
!     m1, m2, m3=masses (any units; m3=outer body)

!     valid for all inclinations

!     MASS RATIO CONDITIONS
!     valid for systems with at least one of  m_2/m_1>0.05  OR  m_3/m_1>0.05
!     (so that one could have, for example,  m_2/m_1=0  and  m_3/m_1=0.1)
!     OR BOTH m2/m1>0.01  AND  m3/m1>0.01
!     **future version will include other resonances to cover smaller mass ratios

!     assumes resonance angle phi=0 because resonance overlap criterion doesn't recognize
!     instability outside separatrix.

!     system is unstable if nstab=1 is returned
!     system is stable if nstab=0 is returned

    integer function nstab(ai,a0,ei0,eo,relinc,m1,m2,m3)

    implicit real*8 (a-h,m,o-z)
    common/params2/mm1,mm2,mm3
    common/savepi/pi
    save itime
    data itime/0/

!     reject outer pericentre inside inner apocentre
    if(a0*(1.-eo) < ai*(1.+ei0))then
        nstab=1
        return
    endif

    if(itime == 0)then
        pi=4.d0*datan(1.0d0)
        itime=1
    endif

    mm1=m1
    mm2=m2
    mm3=m3
          
    m12=m1+m2
    m123=m12+m3

!     set period ratio (outer/inner)
    a0ai=a0/ai
    sigma=sqrt(a0ai**3*m12/m123)
!     do not allow period ratio < 1
    if(sigma < 1.0)then
        nstab=1
        return
    endif
          
    Mi2=m3/m123
    Mo2=(m1*m2/m12**2)*(m12/m123)**(2./3.)
    Mi3=(m3/m12)*(m12/m123)**(4./3.)*(m1-m2)/m12
    Mo3=(m1*m2/m12**2)*(m12/m123)*(m1-m2)/m12
          
    c22=3./8.
    c20=0.25
    c31=sqrt(3.)/4.
    c33=-sqrt(5.)/4.
          
    e=eo
          
!     inclination coefficients

    win=0
          
    A=sqrt(1-ei0**2)*cos(relinc)
    Z=(1-ei0**2)*(1+sin(relinc)**2)+5*ei0**2* &
    (sin(win)*sin(relinc))**2
    Del=z**2+25+16*A**4-10*Z-20*A**2-8*A**2*Z
          
    eK=sqrt(abs((Z+1-4*A**2+sqrt(Del))/6.))
    cosIK=A/sqrt(1-eK**2)
    sinIK=sqrt(1-cosIK**2)
          
    gam222=0.25*(1+cosIK)**2
    gam22m2=0.25*(1-cosIK)**2
    gam220=0.5*sqrt(1.5)*sinIK**2
    gam200=0.5*(3*cosIK**2-1)
               
!     induced inner eccentricity
    ei=ein_induced(sigma,ei0,e,relinc)
               
!     octopole emax
    if(m1 /= m2)then
        eoctmax=eoct(sigma,ei0,e)
        ei=max(eoctmax,ei)
    endif
               
    ei=max(eK,ei)
    ei=min(ei,1.0d0)
               
    n=sigma
    nstab=0
               
!     [n:1](222) resonance
    s221=-3*ei+(13./8.)*ei**3+(5./192.)*ei**5

    f22n=flmn(2,2,n,e)/(1-e)**3

    An=abs(6*c22*s221*f22n*(Mi2+Mo2*sigma**0.666)*gam222)
    phi=0
    En=0.5*(sigma-n)**2-An*(1+cos(phi))

!     [n+1:1](222) resonance
    f22n=flmn(2,2,n+1,e)/(1-e)**3

    An=abs(6*c22*s221*f22n*(Mi2+Mo2*sigma**0.666)*gam222)
               
    Enp1=0.5*(sigma-(n+1))**2-An*(1+cos(phi))
    if(En < 0 .AND. Enp1 < 0)nstab=1
               
!     [n:1](22-2) resonance
    s22m1=-(ei**3*(4480 + 1880*ei**2 + 1091*ei**4))/15360.
    f22n=flmn(2,2,n,e)/(1-e)**3

    An=abs(6*c22*s22m1*f22n*(Mi2+Mo2*sigma**0.666)*gam22m2)
    phi=0
    En=0.5*(sigma-n)**2-An*(1+cos(phi))

!     [n+1:1](22-2) resonance
    f22n=flmn(2,2,n+1,e)/(1-e)**3

    An=abs(6*c22*s22m1*f22n*(Mi2+Mo2*sigma**0.666)*gam22m2)
               
    Enp1=0.5*(sigma-(n+1))**2-An*(1+cos(phi))
    if(En < 0 .AND. Enp1 < 0)nstab=1

!     [n:1](202) resonance
    s201=(ei*(-9216 + 1152*ei**2 - 48*ei**4 + ei**6))/9216.
    f22n=flmn(2,2,n,e)/(1-e)**3

    An=abs(6*sqrt(c20*c22)*s201*f22n* &
    (Mi2+Mo2*sigma**0.666)*gam220)
               
    phi=0
    En=0.5*(sigma-n)**2-An*(1+cos(phi))

!     [n+1:1](202) resonance
    f22n=flmn(2,2,n+1,e)/(1-e)**3

    An=abs(6*sqrt(c20*c22)*s201*f22n &
    *(Mi2+Mo2*sigma**0.666)*gam220)
    Enp1=0.5*(sigma-(n+1))**2-An*(1+cos(phi))
    if(En < 0 .AND. Enp1 < 0)nstab=1

!     [n:1](002) resonance
    s201=(ei*(-9216 + 1152*ei**2 - 48*ei**4 + ei**6))/9216.
    f20n=flmn(2,0,n,e)/(1-e)**3

    An=abs(3*c20*s201*f20n*(Mi2+Mo2*sigma**0.666)*gam200)
               
    phi=0
    En=0.5*(sigma-n)**2-An*(1+cos(phi))

!     [n+1:1](002) resonance
    f20n=flmn(2,0,n+1,e)/(1-e)**3

    An=abs(3*c20*s201*f20n*(Mi2+Mo2*sigma**0.666)*gam200)
               
    Enp1=0.5*(sigma-(n+1))**2-An*(1+cos(phi))
    if(En < 0 .AND. Enp1 < 0)nstab=1
               
    end function nstab

!     -----------------------------------------
!     Asymptotic expression for f^(lm)_n(e) for all e<1 and n.

    real*8 :: function flmn(l,m,n,e)

    implicit real*8 (a-h,o-z)
    common/savepi/pi

    if(e < 5.e-3)then
        if(m == n)then
            flmn=1
        else
            flmn=0
        endif
        return
    endif
          
    rho=n*(1-e)**1.5
          
    xi=(a_cosh(1/e)-sqrt(1-e**2))/(1-e)**1.5
          
    flmn=(1/(2*pi*n))*2.0**m*(sqrt(2*pi)/facfac(l,m))* &
    ((1+e)**(real(3*m-l-1)/4.)/e**m)* &
    (rho**(real(l+m+1)/2.))* &
    exp(-rho*xi)

    END PROGRAM


!     -----------------------------------------
    real*8 :: function ein_induced(sigma,ei0,e,relinc)

    implicit real*8 (a-h,m,o-z)
    common/params2/m1,m2,m3
    common/savepi/pi

    m123=m1+m2+m3
    n=sigma

    gam222=0.25*(1+cos(relinc))**2
    gam220=0.5*sqrt(1.5)*sin(relinc)**2
    gam200=0.5*(3*cos(relinc)**2-1)
                    
    f22n=flmn(2,2,n,e)/(1-e)**3
    f20n=flmn(2,0,n,e)/(1-e)**3
                    
    prod222=f22n*gam222
    prod220=f22n*gam220
    prod200=f20n*gam200
                              
    prod=max(prod222,prod220,prod200)
                                   
    a=4.5*(m3/m123)*(2*pi*n)*prod/sigma**2
                                   
    ein_induced=sqrt(ei0**2+a**2)
                                   
    END PROGRAM
!     -----------------------------------------
!     eoct.f

!     calculates maximum eccentricity for arbitrary coplanar system
!     using Mardling (2007) MNRAS in press

    real*8 :: function eoct(sigma,ei0,eo)
    implicit real*8 (a-h,m,o-z)
    common/params2/m1,m2,m3
    common/savepi/pi

    m12=m1+m2
    m123=m12+m3
    aoai=((m123/m12)*sigma**2)**0.3333
    al=1/aoai
          
    epso=sqrt(1-eo**2)
          
    eeq=1.25*al*eo/epso**2/abs(1-sqrt(al)*(m2/m3)/epso)

    AA=abs(1-ei0/eeq)
          
    if(AA < 1)then
        eoct=(1+AA)*eeq
    else
        eoct=ei0+2*eeq
    endif
          
    END PROGRAM

!     -----------------------------------------
    real*8 :: function a_cosh(x)
    real*8 :: x

    a_cosh=dlog(x+dsqrt(x**2-1.d0))

    END PROGRAM
!     -----------------------------------------
    real*8 :: function sgn(x)
    real*8 :: x

    if(x < 0)then
        sgn=-1
    else
        sgn=1
    endif

    END PROGRAM
!     -----------------------------------------
    real*8 :: function facfac(l,m)
    implicit real*8 (a-h,o-z)

    prod=1

    n=l+m-1

    do i=1,n,2
        prod=prod*i
    enddo

    facfac=prod

    END PROGRAM
