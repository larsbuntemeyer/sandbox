program bsp
  implicit none
  
  integer :: zahl = 5
 
  call unterprogramm(zahl)
end program bsp


subroutine unterprogramm(a)
  implicit none
  
  real, intent(in) :: a
  
  write (*,*) a  
end subroutine unterprogramm
