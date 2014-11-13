program main
!
   implicit none
!
   real :: x,y,z,a,b
   integer :: i,j,k
   integer :: ib=1
   integer :: ie=4
   integer :: jb=1
   integer :: je=2
   integer :: kb=1
   integer :: ke 
!
   write (*,*) 'FORTRAN test program'
!
!
   write(*,*) 'Modulo Test'
!   read(*,*) i
!   read(*,*) j
!   write(*,*) 'mod(',i,',',j,'):',mod(i,j) 
   do i=ib,ie
      do j=jb,je
         write(*,*) 'i',i         
         write(*,*) 'j',j         
         write(*,*) '--------------------------'
      enddo
   enddo
!
   ke = ie*je
!
   write(*,*) '=========================='
   do k=1,ke
      i = 1+int((k-1)/je)
      j = k-(i-1)*je
      write(*,*) 'k',k
      write(*,*) 'i',i
      write(*,*) 'j',j
      write(*,*) '--------------------------'
   enddo
   
!   call test(x,y,z)
!
!  ----------------------------------------------
   contains
!  ----------------------------------------------
!
!  ============================================== 
   subroutine test(a,b,c)
!  ============================================== 

   implicit none

   real, intent(in)    :: a,b
   real, intent(out)   :: c

   write(*,*) 'writing from subroutine test:',c
   c = a*b*c
   write(*,*) 'writing from subroutine test:',c

   end subroutine test
!  ============================================== 
!
end program main
