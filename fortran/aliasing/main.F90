program xxxx
  
  implicit none

  real :: x,y,z

  x = 1.0
  y = 2.0
  z = 3.0

  call sub1(x,y,z)

end program xxxx


subroutine sub1(a,b,c)

  implicit none

  real, intent(in) :: a,b,c,d

  write(*,*) 'a:',a
  write(*,*) 'b:',b
  write(*,*) 'c:',c

end subroutine sub1
