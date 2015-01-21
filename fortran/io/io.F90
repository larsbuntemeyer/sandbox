program bsp
  implicit none
 
  integer :: stat
  character(80) :: str
 
  open(20, file='./testdatei.txt', iostat=stat)
 
  if(stat == 0) then
    write(*,*) 'Das Öffnen der Datei war ein voller Erfolg'
 
    do 
      read(20, '(A)', iostat=stat) str
      ! Bei EOF wird stat /= 0 
      if (stat /= 0) exit
      write(*,*) str
    end do  
  else 
    write(*,*) 'Datei konnte nicht geöffnet werden'
  end if
   
  close(20)
end program bsp
