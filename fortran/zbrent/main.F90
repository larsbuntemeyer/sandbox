program main
!
   implicit none
!
   real :: a,b,c,x
   real :: x1,x2,tol
   real :: test_result, check_result 
   integer :: i
!
   write(*,*), '------------------------------'
   write (*,*) 'FORTRAN test program'
   write(*,*), '------------------------------'
!
   a = 1.5
   b = 5.0
   c = -10.0
   x1  = 0.0
   x2  = 100.0
   tol = 0.000001
!
!   test_result = test_function(a,b,c,x)
!
   write(*,*), 'a: ', a
   write(*,*), 'b: ', b
   write(*,*), 'c: ', c
   write(*,*), 'x: ', x
   write(*,*), 'calling zbrent...'
!
   test_result = zbrent(test_function,x1,x2,tol)
   check_result = test_function(test_result)
!
   write(*,*), '------------------------------'
   write(*,*), 'test_result', test_result
   write(*,*), 'check_result', check_result
   write(*,*), '------------------------------'
!
!   call test(x,y,z)
!
!  ----------------------------------------------
   contains
!  ----------------------------------------------
!
!  ============================================== 
   subroutine test_routine(a,b,c)
!  ============================================== 

   implicit none

   real, intent(in)    :: a,b
   real, intent(out)   :: c

   write(*,*) 'writing from subroutine test:',c
   c = a*b*c
   write(*,*) 'writing from subroutine test:',c

   end subroutine test_routine
!
!  ============================================== 
   real function test_function(x)
!  ============================================== 

   implicit none

   real :: h,i,j,x
!
   h = 1.5
   i = 5.0
   j = -10.0
!
   write(*,*), '------------------------------'
   write(*,*), 'test function called with x=',x
   write(*,*), '------------------------------'
!
   test_function = h*x*x + i*x + j
!
   end function test_function
!  ============================================== 
!
      FUNCTION zbrent(func,x1,x2,tol)
      INTEGER ITMAX
      REAL zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause &
      'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      end function zbrent
end program main
