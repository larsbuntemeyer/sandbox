program bsp
  implicit none
  
  interface gensub
    subroutine writeReal(val)
      real, intent(in) ::val
    end subroutine writeReal 

    subroutine writeInteger(val)
      integer, intent(in) ::val
    end subroutine writeInteger

    subroutine writeCharacter(val)
      character, intent(in) :: val
    end subroutine writeCharacter   
  end interface gensub
 
  call gensub(5.5)
  call gensub(3)
  call gensub("H")
  call writeCharacter("X") 

! Ausgabe:
!   Real-Wert =    5.500000
!   Integer-Wert =            3
!   Zeichen = H
!   Zeichen = X
end program bsp


subroutine writeReal(val)
  real, intent(in) ::val
  write (*,*) "Real-Wert = ", val
end subroutine writeReal 

subroutine writeInteger(val)
  integer, intent(in) ::val
  write (*,*) "Integer-Wert = ", val
end subroutine writeInteger

subroutine writeCharacter(val)
  character, intent(in) ::val
  write (*,*) "Zeichen = ", val
end subroutine writeCharacter
