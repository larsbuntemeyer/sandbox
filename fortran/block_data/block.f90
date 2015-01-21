module common

implicit none

INTEGER A,B,C
REAL    I,J,K,L
COMMON  /AREA1/ A,B,C
COMMON  /AREA2/ I,J,K,L

end module common

program xxxx

use common

implicit none

!EXTERNAL SETUP

write(*,*) a,b,c,i,j,k,l

end program xxxx

BLOCK DATA SETUP
IMPLICIT NONE
INTEGER A,B,C
REAL    I,J,K,L
COMMON  /AREA1/ A,B,C
COMMON  /AREA2/ I,J,K,L
DATA    A,B,C,I,J,K,L/0,1,2,10.0,-20.0,30.0,-40.0/
END
