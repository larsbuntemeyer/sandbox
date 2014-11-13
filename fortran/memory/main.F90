program main
!
   implicit none
!
   real :: x,y,z
   integer :: i
   integer, parameter :: nmax = 200000000
   real, save :: a(nmax)
!
   do i=1,nmax
!
      a(i) = '1.0'
!
   enddo
!
   write(*,*),'done'
!
!  ----------------------------------------------
   contains
!  ----------------------------------------------
!
!  ============================================== 
   subroutine test(a,b,c)
!  ============================================== 

   implicit none

   real, intent(in)      :: a,b
   real, intent(out)     :: c
!
   c = 0.0
!
   end subroutine test
!  ============================================== 
!
end program main
