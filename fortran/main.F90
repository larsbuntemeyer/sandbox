program main
!
   implicit none
!
   real :: x,y,z,a,b
   integer :: i,j
   integer, parameter       :: n=100
   real, dimension(n,n,n)   :: grid
   integer, allocatable, dimension(:)   :: integerList
   integer, pointer, dimension(:,:,:)   :: pointerGrid
!
   character :: text*10
!
   write (*,*) 'FORTRAN test program'
!
   x = 1.0
   y = -1.0
   z = -0.5
   a = 0.5
   b = 2.5
!   
   write (*,*) 'int',x, int(x)
   write (*,*) 'int',y, int(y)
   write (*,*) 'int',z, int(z)
   write (*,*) 'int',b, int(b)
   write (*,*) 'int(2.5)', int(2.5)
   write (*,*) 'int(2.8)', int(2.8)
   write (*,*) 'int(2.3)', int(2.3)
   write (*,*) 'sign(0.5,-0.5)', sign(0.5,-0.5)
   write (*,*) 'ceiling(sign(0.5,-0.5))', ceiling(sign(0.5,-0.5))
   write (*,*) 'sign(0.5,0.5)', sign(0.5,0.5)
   write (*,*) 'ceiling(sign(0.5,0.5))', ceiling(sign(0.5,0.5))
   write (*,*) 'anint',z, anint(z)
   write (*,*) 'anint',a, anint(a)
   write (*,*) 'nint',z, nint(z)
   write (*,*) 'nint',a, nint(a)
   write (*,*) 'ceiling',z, ceiling(z)
   write (*,*) 'ceiling',a, ceiling(a)
   write (*,*) 'ceiling',x, ceiling(x)
   write (*,*) 'ceiling',y, ceiling(y)
   write (*,*) '-1.0 .lt. 0.0',(-1.d0.lt.0.0)
   write (*,*) '-1.0  <   0.0',(-1.d0<0.0)
   write (*,*) '8.0.eq.8?', 8.0.eq.8
   write (*,*) '8.0.eq.1.0*8?', 8.0.eq.1.0*8
   write (*,*) '18.0.ne.18.0?', 18.0.ne.18.0
   write (*,*) 'sin(40.0/180.0*3.141)', sin(40.0/180.0*3.141)
   write (*,*) 'mod(4,4)', mod(4,4)
   write (*,*) 'mod(4.0,4.0)', mod(4.0,4.0)
   write (*,*) 'mod(3.0+1,4.0)', mod(3.0+1,4.0)
   write (*,*) 'mod(3.0,4.0)', mod(3.0,4.0)
   write (*,*) 'equal(0.0,0.0)', equal(0.0,0.0)
   write (*,*) 'equal(2.0,2.0)', equal(2.0,2.0)
   write (*,*) 'equal(2.0,2.00001)', equal(2.0,2.00001)
   write (*,*) 'equal(2.0,2.0000000001)', equal(2.0,2.0000000001)
   write (*,*) 'equal(1.12345678e18,1.12345677e18)', equal(1.12345678e18,1.12345677e18)
   write (*,*) 'equal(1.1234567e18,1.1234566e18)', equal(1.1234567e18,1.1234566e18)
   write (*,*) 'equal(0.e0,1.e-17)', equal(0.e0,1.e-17)
   write (*,*) '1.e-8*1.e-17', 1.e-8*1.e-17
   write (*,*) '1.e-8*1.e-17.lt.1.e-8', 1.e-8*1.e-17.lt.1.e-8
   write (*,*) '1.e-17-0.e0', 1.e-17-0.e0
   write (*,*) '1.e-17+0.e0', 1.e-17+0.e0
   write (*,*) 'abs(1.e-17-0.e0)', abs(1.e-17-0.e0)
   write (*,*) '1.e-17', 1.e-17
   write (*,*) '1.e-8*max(abs(1.e-17),abs(0.e0)', 1.e-8*max(abs(1.e-17),abs(0.e0))
   write (*,*) 'abs(1.e-17-0.e0) .lt. 1.e-8*max(abs(1.e-17),abs(0.e0)', abs(1.e-17-0.e0) .lt. 1.e-8*max(abs(1.e-17),abs(0.e0))
   write (*,*) 'atan2(0.5,0.5)', atan2(0.5,0.5)
   write (*,*) 'atan2(-0.5,0.5)', atan2(-0.5,0.5)
   write (*,*) 'atan2(0.5,-0.5)', atan2(0.5,-0.5)
   write (*,*) 'atan2(-0.5,-0.5)', atan2(-0.5,-0.5)
   write(*,*) 'sizeOf(a):',sizeOf(a)
   write(*,*) 'sizeOf(1.d0):',sizeOf(1.d0)
   write(*,*) 'mod(2.0,2.0)',mod(2.0,2.0)
   write(*,*) 'mod(3.0,2.0)',mod(3.0,2.0)
   write(*,*) 'mod(2,2)',mod(2,2)
   write(*,*) 'mod(3,2)',mod(3,2)
   write(*,*) 'mod(1,2)',mod(1,2)
   write(*,*) '2/2',2/2
   write(*,*) '3/2',3/2
   write(*,*) '5/2',5/2
   write(*,*) '1/2',1/2
   write(*,*) 'sizeof(1.e0)',sizeof(1.e0)
   write(*,*) 'mod(0.5,0.5)',mod(1.5,0.5)
!   call workOnLargeArray(grid) 
!
!  ----------------------------------------------
   contains
!  ----------------------------------------------
!
!  ============================================== 
   subroutine test(a,b,c)
!  ============================================== 

   implicit none

   real, intent(in)    :: a,b
   real, intent(out)   :: c

   write(*,*) 'writing from subroutine test:',c
   c = a*b*c
   write(*,*) 'writing from subroutine test:',c

   end subroutine test
!  ============================================== 
   subroutine workOnLargeArray(a)
!  ============================================== 

   implicit none

   real, dimension(:,:,:), intent(inout)  :: a

   write(*,*) 'writing from subroutine workOnLargeArray:'
   write(*,*) 'a:',a

   end subroutine workOnLargeArray
!  ============================================== 
!
logical function equal(a,b)
 
   implicit none
   
   real :: a,b
   real :: x 
   real :: precision = 1.e-8
 
   equal =  abs(a-b) .le. (precision * max(abs(a),abs(b)))  
 
end function equal





end program main
