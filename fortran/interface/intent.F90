program main
!
implicit none
!
real :: w,x,y,z
!
w = 1.0
x = 2.0
!
call sub1(w,x,y,z)
!
write(*,*) 'w,x:  ',x,y
write(*,*) 'y,z:  ',y,z
!
contains
 !
 subroutine sub1(a,b,c,d)
 !
 implicit none
 !
 real, intent(in)    :: a,b
 real, intent(out)   :: c
 real, intent(inout) :: d
 !
 write(*,*) 'sub1: a,b:  ',a,b
 write(*,*) 'sub1: c,d:  ',c,d
 !
 end subroutine sub1
 !
end program main
