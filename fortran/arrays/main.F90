program main
!
   implicit none
!
   real :: x,y,z,a,b
   real, dimension(4,4) :: arrayA
   real, dimension(4)   :: arrayB
   integer :: i,j
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
!
   do i=0,3
      do j=0,3
         arrayA(i,j) = i*j
      enddo
      arrayB(i) = i
   enddo
!
   do i=0,3
      do j=0,3
         write(*,*)'A',i,j,arrayA(i,j)
      enddo
   enddo
!
   write(*,*)
!
   arrayB = arrayA(0,:)
!
   do i=0,3
      write(*,*) 'B',i,arrayB(i)
   enddo
!
   write(*,*) 'arrayA',arrayA
   write(*,*) 'arrayB',arrayB
!
!   call test(x,y,z)
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
!
end program main
