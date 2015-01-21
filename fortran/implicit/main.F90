program main

  implicit none

  ! a=10
  ! b=20

  write(*,*) square(2.)

  contains

  real function square(x)
     real x
     square = x*x
  end function

end program main
