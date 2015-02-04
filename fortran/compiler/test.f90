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

write(*,100) 'x:',x
write(*,100) 'y:',y
write(*,100) '|x-y|/x:',abs(x-y)/x

100 format(a10,e20.10)
101 format(e20.10)

end program test
