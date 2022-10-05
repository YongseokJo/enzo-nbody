      PROGRAM PARTICLE
      implicit none   

      real, dimension(5) :: m, x, y, z  
      real, dimension(5) :: vx, vy, vz
      integer :: i  
    
      do i = 1,5  
        m(i) = 1
        x(i) = i * 2 / 10 
        y(i) = sin(x(i)) * (1-cos(x(i)/3.0))
        z(i) = sin(x(i))
        vx(i) = i * 2 / 10
        vy(i) = sin(x(i))
        vz(i) = sin(x(i)) * (1-cos(x(i)/3.0))   
      end do  
   
      open(1, file = 'particle1.dat', status='new')  
        do i = 1,5
        write(1,1000) m(i),x(i), y(i), z(i), vx(i), vy(i), vz(i)   
 1000 format (8(F10.5))
      end do  

      close(1) 
   
      END