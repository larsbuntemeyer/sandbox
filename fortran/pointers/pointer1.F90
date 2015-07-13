program bsp
 implicit none
 
  integer, pointer       :: ptr1 => null()
  character(20), pointer :: ptr2
  character(20), target  :: str
 
  str = "Hallo, Welt!"
  ptr2 => str
 
  write(*,*) associated(ptr1)
  ! Ausgabe: F
 
  write(*,*) associated(ptr2)
  ! Ausgabe: T
 
  write(*,*) associated(ptr2, str)
  ! Ausgabe: F
 
  write(*,*) associated(ptr2, str)
  ! Ausgabe: T
end program bsp
