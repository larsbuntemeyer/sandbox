program test

implicit none

real :: a,b,c,d
real :: x,y 

a=1.13423
b=2.32232
c=3.34504
d=4.0

x = a*(b+c)
y = a*b+a*c

write(*,*) 'x:',x
write(*,*) 'y:',y
write(*,'(a10,e20.10)') '|x-y|/x:',abs(x-y)/x

end program test
