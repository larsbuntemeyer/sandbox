program bsp
  implicit none
 
  real, pointer :: ptr(:)
  real, target  :: trg(1000),trg1(1000)
 
  trg = 5.5
!  write(*,*) sizeof(trg)
!  write(*,*) sizeof(ptr)
  ptr => trg
!  write(*,*) ptr,trg
!  write(*,*) sizeof(trg)
!  write(*,*) sizeof(ptr)
  ! Ausgabe: 5.50000
  ! Zeiger werden bei Bedarf automatisch dereferenziert
end program bsp

