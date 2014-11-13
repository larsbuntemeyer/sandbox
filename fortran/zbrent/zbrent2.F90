!********************************!
!                                !
! File contains:                 !
!               function zbrent. !
!                                !
!********************************!

!*******************************************************************
!
      FUNCTION ZBRENT(func,x1,x2,tol)
!
! This is a bisection routine. When ZBRENT is called, we provide a
! reference to a particular function and also two values which bound
! the arguments for the function of interest. ZBRENT finds a root of
! the function (i.e. the point where the function equals zero), that
! lies between the two bounds.  For a full description see Press et
! al. (1986).
!
!*******************************************************************
      real,intent(in)   :: tol,x1,x2
      real :: func !OLI: added to satisfy the implicit none compiler option.
      external          :: func
      real              :: zbrent

      ! internal variables...
      integer           :: iter
      integer,parameter :: ITMAX=30
      real              :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      real,parameter    :: EPS=3.e-8

      ! calculations...
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if ((fa.gt.0..and.fb.gt.0.) .or. (fa.lt.0..and.fb.lt.0.)) then
                fa=func(a)
                fb=func(b)
            write(*,*)'             fa            fb              x1                x2'
                write(*,*)fa,fb,x1,x2
                write(*,*)"ZBRENT.F90: root must be bracketed"
                fa=func(a)
                fb=func(b)
          endif
      c=b
      fc=fb
      do iter=1,ITMAX
        if ((fb.gt.0..and.fc.gt.0.) .or. (fb.lt.0..and.fc.lt.0.)) then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if (abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if (abs(xm).le.tol1 .or. fb.eq.0.) then
          zbrent=b
          return
        endif
        if (abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if (a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if (p.gt.0.) q=-q
          p=abs(p)
          if (2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
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
        if (abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
      enddo
      write(*,*)"ZBRENT: exceeding maximum iterations"
      zbrent=b

END FUNCTION
!
!********************************************* 
