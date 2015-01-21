program xxx
  implicit none
  real :: x,y,z
  x = 1.0
  y = 2.0
  z = 3.0
  call sub(x,y,z)
  write(*,*) x,y,z
end program xxx

subroutine sub(a,b,c)
  implicit none
  real, intent(in)  ::  a,b
  real, intent(out) ::  c
  c = a+b
  call sub1(a,b,c)
end subroutine sub

subroutine sub1(a,b,c)
  implicit none
  real, intent(in) :: a,b
  real, intent(out):: c
  c = c+a+b
end subroutine sub1
