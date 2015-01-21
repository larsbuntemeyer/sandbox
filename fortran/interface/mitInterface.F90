program bsp
  implicit none
  
  integer :: zahl = 5
  
  interface
    subroutine unterprogramm(a)
      real, intent(in) :: a
    end subroutine unterprogramm    
  end interface
  
  call unterprogramm(zahl)
end program bsp


subroutine unterprogramm(a)
  implicit none
  
  real, intent(in) :: a
  
  write (*,*) a  
end subroutine unterprogramm
