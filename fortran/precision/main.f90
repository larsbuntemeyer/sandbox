program bsp
  implicit none
 
  real(kind=8) :: variable
  ! variable ist nun vom Datentyp "double precision"
 
  variable = 1.55555555555_8
 
  write(*,*) variable
  ! Ausgabe: 1.55555555555000
  
  ! Hier wird nur eine gewöhnliche real-Zahl mit 7 Nachkommastellen 
  ! zugewiesen
  variable = 1.55555555555
 
  write(*,*) variable
  ! Ausgabe: 1.55555558204651
end program bsp
