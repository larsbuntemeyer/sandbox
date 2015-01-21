program bsp
  implicit none
 
  character (len = 80) :: a
  integer              :: st = 0            
 
  open (20, file='./testdatei.txt', status='OLD', iostat=st)
  
  if (st /= 0) then    
    stop "open-Fehler!"
  end if  
 
  ! Aus Datei lesen
  do
    read (20, 888, iostat=st) a
    ! Auf Standardausgabe schreiben   
    
    if (st == 0) then
      write (*, 888) a
    else if (st > 0) then  
      write (*,*) "read-Fehler!"
      exit
    else if (st < 0) then
      exit
    end if  
  end do
 
  close(20)
 
 888  format(A) 
end program bsp
