!
function diffuse(x,t,D,p,x0)
!
   implicit none
   double precision :: diffuse
   double precision :: x,t,D
   integer          :: p
!
  ! diffuse = 1.d0/(2.0**p * (3.14159265d0*D*t)**(0.5d0*p)) * exp(-)  
!
   return diffuse
!
end function diffuse
!
program diffusion
!
   implicit none
!
   integer, parameter          :: nx = 100
   integer, parameter          :: nt = 10
   double precision, parameter :: x0 = 40.0
   double precision, parameter :: x1 = 60.0
   double precision, parameter :: dx = 1.0
   double precision, parameter :: dt = 1.e-1
   double precision, parameter :: q0 = 1.d-1
   double precision, parameter :: q1 = 1.d1
   double precision, parameter :: D  = 1.0
   !
   double precision   :: q(nx),q_old(nx)
   double precision   :: ix,it,x,t
   !
   do ix=1,nx
      x = dx*ix
      if(x<x0.or.x>x1) then
        q(ix) = q0
      else
        q(ix) = q1
      endif 
   enddo
!
   do it=1,nt
      q_old = q
      q(1)  = q_old(2) 
      q(nx) = q_old(nx-1) 
      do ix=2,nx-1
!        q(ix) = q_old(ix) + D * dt/(dx*dx) * (q_old(ix+1)-2.0*q_old(ix)+q_old(ix-1))    
        q(ix) = q_old(ix) + D * dt/(dx*dx) * (q_old(ix+1)-2.0*q_old(ix)+q_old(ix-1))    
      enddo
   enddo
!
   open(unit=1,file='q_final.out')
   do ix=1,nx
      write(1,*) q(ix)
   enddo 
   close(1)
!
!
end program diffusion
