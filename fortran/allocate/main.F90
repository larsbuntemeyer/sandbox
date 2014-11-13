program main
!
   implicit none
!
   real :: x,y,z,a,b
   integer :: i,j
   integer, parameter       :: n=100
   integer, allocatable, dimension(:,:)   :: integerList
!
   write (*,*) 'FORTRAN test program'
!
!
   allocate(integerList(n,n))
   integerList = 0
   allocate(integerList(n,n))
   deallocate(integerList)



end program main
