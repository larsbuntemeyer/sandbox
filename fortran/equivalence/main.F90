program main
  !
  implicit none
  !
  real :: x
  character(10) :: y = 'Hallo'
  !
  equivalence(x,y)
  !
  x = 1.0
  !
  write(*,*) x,y
  !
end program main
