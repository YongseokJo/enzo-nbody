!***************************************************************

    SUBROUTINE tred2(a,n,np,d,e)
    INTEGER :: n,np
    DOUBLE PRECISION :: a(np,np),d(np),e(np)
    INTEGER :: i,j,k,l
    DOUBLE PRECISION :: f,g,h,hh,scale
    do 18 i=n,2,-1
        l=i-1
        h=0.d0
        scale=0.d0
        if(l > 1)then
            do 11 k=1,l
                scale=scale+abs(a(i,k))
            11 END DO
            if(scale == 0.d0)then
                e(i)=a(i,l)
            else
                do 12 k=1,l
                    a(i,k)=a(i,k)/scale
                    h=h+a(i,k)**2
                12 END DO
                f=a(i,l)
                g=-sign(sqrt(h),f)
                e(i)=scale*g
                h=h-f*g
                a(i,l)=f-g
                f=0.d0
                do 15 j=1,l
                !     Omit following line if finding only eigenvalues
                    a(j,i)=a(i,j)/h
                    g=0.d0
                    do 13 k=1,j
                        g=g+a(j,k)*a(i,k)
                    13 END DO
                    do 14 k=j+1,l
                        g=g+a(k,j)*a(i,k)
                    14 END DO
                    e(j)=g/h
                    f=f+e(j)*a(i,j)
                15 END DO
                hh=f/(h+h)
                do 17 j=1,l
                    f=a(i,j)
                    g=e(j)-hh*f
                    e(j)=g
                    do 16 k=1,j
                        a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                    16 END DO
                17 END DO
            endif
        else
            e(i)=a(i,l)
        endif
        d(i)=h
    18 END DO
!     Omit following line if finding only eigenvalues.
    d(1)=0.d0
    e(1)=0.d0
    do 24 i=1,n
    !     Delete lines from here ...
        l=i-1
        if(d(i) /= 0.d0)then
            do 22 j=1,l
                g=0.d0
                do 19 k=1,l
                    g=g+a(i,k)*a(k,j)
                19 END DO
                do 21 k=1,l
                    a(k,j)=a(k,j)-g*a(k,i)
                21 END DO
            22 END DO
        endif
    !     ... to here when finding only eigenvalues.
        d(i)=a(i,i)
    !     Also delete lines from here ...
        a(i,i)=1.d0
        do 23 j=1,l
            a(i,j)=0.d0
            a(j,i)=0.d0
        23 END DO
    !     ... to here when finding only eigenvalues.
    24 END DO
    return
    end SUBROUTINE tred2
!  (C) Copr. 1986-92 Numerical Recipes Software 5.d0).

!***************************************************************

    SUBROUTINE tqli(d,e,n,np,z)
    INTEGER :: n,np
    DOUBLE PRECISION :: d(np),e(np),z(np,np)
!U    USES pythag
    INTEGER :: i,iter,k,l,m
    DOUBLE PRECISION :: b,c,dd,f,g,p,r,s,pythag
    do 11 i=2,n
        e(i-1)=e(i)
    11 END DO
    e(n)=0.d0
    do 15 l=1,n
        iter=0
        1 do 12 m=l,n-1
            dd=abs(d(m))+abs(d(m+1))
            if (abs(e(m))+dd == dd) goto 2
        12 END DO
        m=n
        2 if(m /= l)then
        !         if(iter.eq.30)pause 'too many iterations in tqli'
            iter=iter+1
            g=(d(l+1)-d(l))/(2.d0*e(l))
            r=pythag(g,1.d0)
            g=d(m)-d(l)+e(l)/(g+sign(r,g))
            s=1.d0
            c=1.d0
            p=0.d0
            do 14 i=m-1,l,-1
                f=s*e(i)
                b=c*e(i)
                r=pythag(f,g)
                e(i+1)=r
                if(r == 0.d0)then
                    d(i+1)=d(i+1)-p
                    e(m)=0.d0
                    goto 1
                endif
                s=f/r
                c=g/r
                g=d(i+1)-p
                r=(d(i)-g)*s+2.d0*c*b
                p=s*r
                d(i+1)=g+p
                g=c*r-b
            !     Omit lines from here ...
                do 13 k=1,n
                    f=z(k,i+1)
                    z(k,i+1)=s*z(k,i)+c*f
                    z(k,i)=c*z(k,i)-s*f
                13 END DO
            !     ... to here when finding only eigenvalues.
            14 END DO
            d(l)=d(l)-p
            e(l)=g
            e(m)=0.d0
            goto 1
        endif
    15 END DO
    return
    end SUBROUTINE tqli
!  (C) Copr. 1986-92 Numerical Recipes Software 5.d0).

!***************************************************************

    FUNCTION pythag(a,b)
    DOUBLE PRECISION :: a,b,pythag
    DOUBLE PRECISION :: absa,absb
    absa=abs(a)
    absb=abs(b)
    if(absa > absb)then
        pythag=absa*sqrt(1.d0+(absb/absa)**2)
    else
        if(absb == 0.d0)then
            pythag=0.d0
        else
            pythag=absb*sqrt(1.d0+(absa/absb)**2)
        endif
    endif
    return
    end FUNCTION pythag
!  (C) Copr. 1986-92 Numerical Recipes Software 5.d0).

