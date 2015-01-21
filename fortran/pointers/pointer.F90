program bsp
  implicit none
 
  real, pointer :: ptr
  real, target  :: trg
 
  trg = 5.5
  ptr => trg
  write(*,*) ptr,trg
  ! Ausgabe: 5.50000
  ! Zeiger werden bei Bedarf automatisch dereferenziert
end program bsp

