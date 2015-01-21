program main
!
   use test_module_interface, only: testA
!
   implicit none
!
   real :: x,y,z,a,b
   real, dimension(4,4) :: arrayA
   real, dimension(4)   :: arrayB
   integer :: i,j
!
   call testA()
   call one_d_motion 
!
end program main
