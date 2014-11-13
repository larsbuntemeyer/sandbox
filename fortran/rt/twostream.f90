
!-------------------------------------------------------------------------
!             SUBROUTINE FOR OLSON & KUNASZ (1987) QUADRATURE
!-------------------------------------------------------------------------
function qdr_olsonkunasz(int,src_prev,src_curr,src_next, &
                         dtau_prev,dtau_next,order,      &
                         u,v,w)
  implicit none
  double precision :: qdr_olsonkunasz
  double precision :: int,src_prev,src_curr,src_next
  double precision :: dtau_prev,dtau_next,u,v,w
  integer :: order
  double precision :: xp,e0,e1,e2
  !
  ! Precalculate e^{-\Delta\tau}
  !
  xp = exp(-dtau_prev)
  !
  ! Calculate
  !
  select case(order) 
  case(1)
     !
     ! First order integration
     !
     u  = 0.5d0*(1.d0-xp)
     v  = 0.5d0*(1.d0-xp)
     w  = 0.d0
  case(2)
     !
     ! Second order integration
     !
     u  = (1.d0-(1.d0+dtau_prev)*xp)/dtau_prev
     v  = (dtau_prev-1+xp)/dtau_prev
     w  = 0.d0
  case(3)
     !
     ! Third order integration
     !
     e0 = 1.d0-xp
     e1 = dtau_prev-e0
     e2 = dtau_prev*dtau_prev-2.d0*e1
     u  = e0+(e2-(2.d0*dtau_prev+dtau_next)*e1)/(dtau_prev*(dtau_prev+dtau_next))
     v  = ((dtau_prev+dtau_next)*e1-e2)/(dtau_prev*dtau_next)
     w  = (e2-dtau_prev*e1)/(dtau_next*(dtau_prev+dtau_next))
  case default
     write(*,*) "ERROR: Do not know order = ",order
     stop
  end select
  !
  ! Now return the new intensity
  !
  qdr_olsonkunasz = int*xp + u*src_prev + v*src_curr + w*src_next
  return
end function qdr_olsonkunasz


!-------------------------------------------------------------------------
!                           NG ACCELERATOR
!-------------------------------------------------------------------------
subroutine ng_accel(sppprev,spprev,sprev,scur,n)
  implicit none
  double precision :: sppprev(n),spprev(n),sprev(n),scur(n)
  double precision :: q1(n),q2(n),q3(n)
  double precision :: a1,b1,c1,a2,b2,c2,a,b
  integer :: n,i
  do i=1,n
     q1(i) = scur(i)-2*sprev(i)+spprev(i)
     q2(i) = scur(i)-sprev(i)-spprev(i)+sppprev(i)
     q3(i) = scur(i)-sprev(i)
  enddo
  a1 = 0.d0
  b1 = 0.d0
  c1 = 0.d0
  a2 = 0.d0
  b2 = 0.d0
  c2 = 0.d0
  do i=1,n
     a1 = a1 + q1(i)*q1(i)
     b1 = b1 + q1(i)*q2(i)
     c1 = c1 + q1(i)*q3(i)
     a2 = a2 + q2(i)*q1(i)
     b2 = b2 + q2(i)*q2(i)
     c2 = c2 + q2(i)*q3(i)
  enddo
  a  = (c1*b2-c2*b1)/(a1*b2-a2*b1)
  b  = (c2*a1-c1*a2)/(a1*b2-a2*b1)
  do i=1,n
     scur(i) = (1.d0-a-b)*scur(i)+a*sprev(i)+b*spprev(i)
  enddo
end subroutine ng_accel


!-------------------------------------------------------------------------
!                TRIDIAGONAL SOLVER OF "NUMERICAL RECIPES"
!-------------------------------------------------------------------------
subroutine tridag(a,b,c,r,u,n)
  implicit none
  integer :: n
  double precision :: a(n),b(n),c(n),r(n),u(n)
  integer :: j
  double precision :: bet,gam(n)
  if(b(1).eq.0.d0) stop 'tridag: rewrite equations'
  bet=b(1)
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j)*gam(j)
     if(bet.eq.0.) stop 'tridag failed'
     u(j)=(r(j)-a(j)*u(j-1))/bet
  enddo
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  enddo
  return
end subroutine tridag


!-------------------------------------------------------------------------
!                             MAIN PROGRAM
!
! The equation to solve is the following radiative transfer equation:
!
!   S = eps * B + (1-eps)*Lambda[S]
!
! where Lambda is the Lambda Operator. 
!-------------------------------------------------------------------------
program twostream
  implicit none
  integer, parameter :: nzc=64
  integer, parameter :: nzi=nzc+1
  double precision :: mu=1.0
  double precision, parameter :: alpha0=1d7,alpha1=1d-1
  double precision :: bnu(nzi),dtau(nzc),dtaui(nzi)
  double precision :: zi(nzi),zc(nzc),dz(nzc),alpha(nzc)
  double precision :: xp(nzc),taui(nzi)
  double precision :: eps,int,u,v,w
  double precision :: dmt_a(nzi),dmt_b(nzi),dmt_c(nzi),drs(nzi)
  double precision :: mat_a(nzi),mat_b(nzi),mat_c(nzi),rhs(nzi)
  double precision :: lambda_a(nzi),lambda_b(nzi),lambda_c(nzi)
  double precision :: j_diff(nzi),s_diff(nzi)
  double precision :: j_iter(nzi),s_iter(nzi),s_bk(nzi,4)
  double precision :: j_tri(nzi)
  double precision :: sp,sc,sn,dtp,dtn
  double precision :: qdr_olsonkunasz
  integer :: niter,order,ali,ng,iter
  integer :: i,ing
  !
  ! The input parameters
  !
  write(*,*) 'Photon destruction probability (epsilon)?'
  read(*,*) eps
  write(*,*) 'Nr of iterations?'
  read(*,*) niter
  write(*,*) 'Order of integration (1, 2 or 3)?'
  read(*,*) order
  write(*,*) 'Iteration method (0=LI, 1=ALI local, 2=ALI tridiag)?'
  read(*,*) ali
  write(*,*) 'Ng acceleration (0=no, 1=yes)?'
  read(*,*) ng
  !
  ! Make the grid from 0 to 1
  !
  ! ...First the cell wall (interfaces)
  !
  do i=1,nzi
     zi(i) = (i-1.d0)/(nzi-1.d0)
  enddo
  !
  ! ...Then the cell centers
  !
  do i=1,nzc
     zc(i) = 0.5d0*(zi(i)+zi(i+1))
  enddo
  !
  ! ...Then the cell widths Delta z
  !
  do i=1,nzc
     dz(i) = zi(i+1)-zi(i)
  enddo
  !
  ! Make an exponential atmosphere
  !
  do i=1,nzc
     alpha(i) = exp(log(alpha0)+(log(alpha1)-log(alpha0))*zc(i))
  enddo
  !
  ! Compute the dtau over each cell, including the mu factor
  !
  do i=1,nzc
     dtau(i)  = alpha(i)*dz(i)/mu
  enddo
  !
  ! For the diffusion equation we need also the 
  ! Delta tau from cell-center to cell-center
  !
  do i=2,nzi-1
     dtaui(i) = 0.5d0*(dtau(i)+dtau(i-1))
  enddo
  dtaui(1)  = 2*dtaui(2)-dtaui(3)
  dtaui(nzi) = 2*dtaui(nzi-1)-dtaui(nzi-2)
  !
  ! Put the Planck function to 1 everywhere
  !
  do i=1,nzi
     bnu(i) = 1.d0
  enddo
  !
  ! Precalculate exp(-Deltatau)
  !
  do i=1,nzc
     xp(i) = exp(-dtau(i))
  enddo
  !
  ! For your convenience: The integrated tau from top to bottom 
  ! along the direction cos(theta)=mu
  !
  taui(nzi) = dtaui(nzi)
  do i=nzi-1,1,-1 
     taui(i) = taui(i+1) + dtaui(i)
  enddo
  !
  ! Write the tau grid to file
  !
  open(unit=1,file='taui.out')
  do i=1,nzi
     write(1,*) taui(i)
  enddo
  close(1)
  open(unit=1,file='dtau.out')
  do i=1,nzc
     write(1,*) dtau(i)
  enddo
  close(1)
  !
  !==========================================================
  !
  ! First solve the diffusion equation, so get the "true" solution
  ! to the two-stream equations
  !
  !   d^2J/dz^2 = 3*alpha^2*eps* (J-B)
  !
  ! The Matrix elements
  ! 
  do i=2,nzi-1
     dmt_a(i) = -1.d0/(dtaui(i)*dtau(i-1))
     dmt_b(i) = eps+1.d0/(dtaui(i)*dtau(i-1))+1.d0/(dtaui(i)*dtau(i))
     dmt_c(i) = -1.d0/(dtaui(i)*dtau(i))
  enddo
  !
  ! Right hand side
  !
  do i=2,nzi-1
     drs(i)   = eps*bnu(i)
  enddo
  !
  ! Boundary condition left
  !
  dmt_a(1) = 0.d0
  dmt_b(1) = 1.d0
  dmt_c(1) = 0.d0
  drs(1)   = bnu(1)
  !
  ! Boundary condition right
  !
  dmt_a(nzi) = -1.d0/dtau(nzi-1)
  dmt_b(nzi) = 1.d0/dtau(nzi-1)+1.d0    ! Note: No 1/sqrt(3.d0) because tau is along mu=1/sqrt(3.d0) already
  dmt_c(nzi) = 0.d0
  drs(nzi)   = 0.d0
  !
  ! Solve for J
  !
  call tridag(dmt_a,dmt_b,dmt_c,drs,j_diff,nzi)
  !
  ! Compute S
  !
  do i=1,nzi
     s_diff(i) = eps*bnu(i) + (1.d0-eps)*j_diff(i)
  enddo
  !
  ! Write the solution to file
  !
  open(unit=1,file='s_diff.out')
  do i=1,nzi
     write(1,*) s_diff(i)
  enddo
  close(1)
  !
  !==========================================================
  !
  ! Now (A)LI 
  !
  ! Open file
  !
  open(unit=2,file='s_iter.out')
  !
  ! Initial guess
  !
  do i=1,nzi
     s_iter(i) = eps*bnu(i)
  enddo
     do i=1,nzi
        write(2,*) s_iter(i)
     enddo
     write(2,*)
  !
  ! Reset Ng counter
  !
  ing = -1
  !
  ! Start the (A)LI iteration
  !
  do iter=1,niter 
     !
     ! Reset some stuff
     !
     do i=1,nzi
        lambda_a(i) = 0.d0
        lambda_b(i) = 0.d0
        lambda_c(i) = 0.d0
        mat_a(i)    = 0.d0
        mat_b(i)    = 0.d0
        mat_c(i)    = 0.d0
        j_iter(i)   = 0.d0
     enddo
     !
     ! Starting condition
     !
     int       = bnu(1)
     j_iter(1) = j_iter(1) + 0.5d0*int
     !
     ! Integrate from bottom to top
     !
     do i=2,nzi
        !
        ! Prepare the arguments for the qdr_olsonkunasz subroutine
        !
        sp  = s_iter(i-1)
        sc  = s_iter(i)
        dtp = dtau(i-1)
        if(i.ne.nzi) then 
           sn  = s_iter(i+1)
           dtn = dtau(i)
        else
           sn  = sc
           dtn = dtaui(i)
        end if
        !
        ! Call the qdr_olsonkunasz subroutine to make one quadrature step
        !
        int = qdr_olsonkunasz(int,sp,sc,sn,dtp,dtn,order,u,v,w)
        !
        ! Update the mean intensity
        !
        j_iter(i) = j_iter(i) + 0.5d0*int
        !
        ! Update the Lambda^* matrix elements
        !
        lambda_a(i) = lambda_a(i) + 0.5d0*u
        lambda_b(i) = lambda_b(i) + 0.5d0*v
        lambda_c(i) = lambda_c(i) + 0.5d0*w
     enddo
     ! 
     ! Starting condition at the top
     !
     int         = 0.d0
     j_iter(nzi) = j_iter(nzi) + 0.5d0*int
     !
     ! Integrate from top to bottom
     !
     do i=nzi-1,1,-1
        !
        ! Prepare the arguments for the qdr_olsonkunasz subroutine
        !
        sp  = s_iter(i+1)
        sc  = s_iter(i)
        dtp = dtau(i)
        if(i.ne.1) then
           sn  = s_iter(i-1)
           dtn = dtau(i-1)
        else
           sn  = sc
           dtn = dtaui(1)
        endif
        !
        ! Call the qdr_olsonkunasz subroutine to make one quadrature step
        !
        int = qdr_olsonkunasz(int,sp,sc,sn,dtp,dtn,order,u,v,w)
        !
        ! Update the mean intensity
        !
        j_iter(i) = j_iter(i) + 0.5d0*int
        !
        ! Update the Lambda^* matrix elements
        !
        lambda_a(i) = lambda_a(i) + 0.5d0*w
        lambda_b(i) = lambda_b(i) + 0.5d0*v
        lambda_c(i) = lambda_c(i) + 0.5d0*u
     enddo
     !
     ! Now create the M^* matrix
     !
     do i=1,nzi 
        mat_a(i) = -(1.d0-eps)*lambda_a(i)
        mat_b(i) = 1.d0-(1.d0-eps)*lambda_b(i)
        mat_c(i) = -(1.d0-eps)*lambda_c(i)
     enddo
     !
     ! Now do LI, ALI-local or ALI-tridiag
     !
     select case(ali)
     case(0)
        !
        ! Lambda Iteration
        !
        do i=1,nzi
           s_iter(i) = eps*bnu(i) + (1.d0-eps)*j_iter(i)
        enddo
     case(1)
        !
        ! ALI with local operator
        !
        do i=1,nzi
           s_iter(i) = ( eps*bnu(i) + (1.d0-eps)*                     &
                         (j_iter(i)-lambda_b(i)*s_iter(i)) ) / mat_b(i)
        enddo
     case(2)
        !
        ! ALI with tridiagonal operator
        !
        ! First multiply Lambda^* with S
        !
        do i=2,nzi-1
           j_tri(i) = lambda_a(i)*s_iter(i-1) + &
                      lambda_b(i)*s_iter(i) +   &
                      lambda_c(i)*s_iter(i+1)
        enddo
        j_tri(1)   = lambda_b(1)*s_iter(1)+lambda_c(1)*s_iter(2)
        j_tri(nzi) = lambda_a(nzi)*s_iter(nzi-1)+lambda_b(nzi)*s_iter(nzi)
        !
        ! Make the right-hand side
        !
        do i=1,nzi
           rhs(i) = eps*bnu(i) + (1.d0-eps)*(j_iter(i)-j_tri(i))
        enddo
        !
        ! Solve the tridiagonal matrix equation
        !
        call tridag(mat_a,mat_b,mat_c,rhs,s_iter,nzi)
     end select
     !
     ! The Ng acceleration
     !
     if(ng.eq.1) then
        ing = ing+1
     endif
     if((ing.ge.0).and.(ing.lt.4)) then 
        do i=1,nzi
           s_bk(i,ing+1) = s_iter(i)
        enddo
     endif
     if(ing.eq.3) then
        call ng_accel(s_bk(:,1),s_bk(:,2),s_bk(:,3),s_iter(:),nzi)
        ing = -1
     endif
     !
     ! End of iteration loop
     !
     ! Write this snapshot
     !
     do i=1,nzi
        write(2,*) s_iter(i)
     enddo
     write(2,*)
  enddo
  !
  ! Close the file
  !
  close(2)
  !
  ! For convenience also write the final answer
  !
  open(unit=1,file='s_final.out')
  do i=1,nzi
     write(1,*) j_iter(i)
  enddo
  close(1)
end program twostream
