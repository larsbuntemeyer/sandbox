!
!
!  Analytical Solver for a Radiation Hydrodynamics Shock
!
!
program RHD_Shock
!
   implicit none
!
   integer, parameter   :: n     = 100 
   real                 :: M0    = 50.0
   real                 :: x0    = -0.4
   real                 :: x1    = 0.0
   real, dimension(n)   :: x,temp,rhoT
   integer              :: i 
   real                 :: dx
!
   real                 :: gasConst    = 8.31447E+07
   real                 :: gamma       = 5.0/3.0
   real                 :: dens_ref    = 1.0
   real                 :: clight      = 2.99792e+10
   real                 :: L_ref       = 1.0
   real                 :: radConst    = 7.56577e-15
  
   real                 :: trans_opac
   real                 :: kappa       = 1.e-4
   real                 :: cs_ref,cs_star 
!
   real                 :: P0          = 1.e-4
   real                 :: T0          = 1.1604505e6 ! 100 eV
   real                 :: T1          = 1.0 
!
!   T1                = 2.07756999533*T0
   T0 = 1.0
   cs_ref            = sqrt(gamma*gasConst*T0)
   !cs_star           = (1.0+4.0/9.0*P0*(gamma-1.0)*(4.0*P0*gamma+3.0*(5.0-3.0*gamma)/(1.0+4.0*P0*gamma*(gamma-1.0))))**0.25 
   cs_star           = 4.0/9.0 * P0 
   dens_ref          = radConst*(T0**4)/(P0*cs_ref**2)
   trans_opac        = (1.0/kappa)*(clight*radConst*(T0**4)/(3.0*(cs_ref**3)*L_ref*dens_ref))
   T1                = ((1.0-gamma+2.0*gamma*M0*M0)*(2.0+(gamma-1.0)*M0*M0))/(((gamma+1.0)**2)*M0*M0)
!
   write(*,*) ''
   write(*,*) ''
   write(*,*) '====================================================='
   write(*,*) 'Analytical Solver for a Radiation Hydrodynamics Shock'
   write(*,*) '====================================================='
   write(*,*) ''
   write(*,*) ' Shock Parameter'
   write(*,*) ''
   write(*,*) ' gamma  :', gamma
   write(*,*) ' P0     :', P0
   write(*,*) ' M0     :', M0
   write(*,*) ' T0     :', T0
   write(*,*) ' T1     :', T1
   write(*,*) ' kappa  :',kappa
   write(*,*) ''
   write(*,*) ' Grid Parameters'
   write(*,*) ' n      :', n
   write(*,*) ' x0     :', x0
   write(*,*) ' x1     :', x1
   write(*,*) ' dx     :', abs(x1-x0)/n
   write(*,*) ''
   write(*,*) '-----------------------------------------------------'
   write(*,*) ''
   write(*,*) 'Starting Integrations...'
   write(*,*) ''
!
!
   temp(n) = T1
   rhoT(n) = rho(T1)
   dx      = (x1-x0)/n
!
!  prepare x-coordinates
!
   x(1) = x0
   do i=2,n
      x(i) = x(i-1) + dx
   enddo
!
!  integrate temperature backwards from x1 to x0
!
   call integrate_precursor(rhoT,temp,x0,x1)
!
!  update density
!
   do i=1,n
      rhoT(i) = rho(temp(i))
   enddo
!
!
!  output
!
   open(unit=1,file='data.out')
   do i=1,n
      write(1,*) x(i), rhoT(i), temp(i)
   enddo
   close(1)
!
   write(*,*) 'wrote output to data.out'
   write(*,*) 'SUCCESS'
   write(*,*) ''
   write(*,*) ''
   write(*,*) ''
!
!
!
!  ----------------------------------------------
   contains
!  ----------------------------------------------
!
!  ============================================== 
   real function m(T)
!  ============================================== 
!
   implicit none
!
   real :: T
!
   m = 0.5*(gamma*M0**2+1.0) + &
       gamma*P0/6.0*(1.0-T**4)
!
   return
!
   end function m
!  ============================================== 
!
!  ============================================== 
   real function rho(T)
!  ============================================== 
!
   implicit none
!
   real :: T
   real :: mT,discr
!
   mT = m(T)
!
   if(mT**2>gamma*T*M0**2) then
      discr = mT**2-gamma*T*M0**2
   else
      discr = 0.0
   endif
!
   if(T>0.0) then
      rho = (mT-sqrt(discr))/T
   else
      rho = 0.0
      write(*,*) 'WARNING, T=0 and rho=0'
   endif
!
   return
!
   end function rho
!  ============================================== 
!
!  ============================================== 
   real function f1(T)
!  ============================================== 
!
   implicit none
!
   real :: T
!
   f1 = 3.0*(gamma+1)*(T-1.0)- &
        P0*gamma*(gamma-1.0)*(7.0+T**4)
!
   return
!
   end function f1
!  ============================================== 
!
!  ============================================== 
   real function f2(T)
!  ============================================== 
!
   implicit none
!
   real :: T
!
   f2 = 12.0*(gamma-1)**2 &
        *T*(3.0+gamma*P0*(1.0+7.0*T**4))
!
   return
!
   end function f2
!  ============================================== 
!
!  ============================================== 
   real function f3(rho,T)
!  ============================================== 
!
   implicit none
!
   real :: T,rho
!
   f3 = 6.0*rho*rho*(T-1.0)/(gamma-1.0) + &
        3.0*(1.0-rho*rho)*M0*M0       + &
        8.0*P0*(T**4-rho)*rho         
!
   return
!
   end function f3
!  ============================================== 
!
!  ============================================== 
   subroutine integrate_precursor(rhoT,temp,x0,x1)
!  ============================================== 
!
!  integrates the temperature backwards from
!  x1 to x0
!
   implicit none
!
   real, dimension(n), intent(inout) :: temp
   real, dimension(n), intent(inout) :: rhoT
   real, intent(in)                  :: x0,x1
   real                              :: dx
   integer                           :: i
!
   dx = (x0-x1)/n
!
!  integration loop
!
   do i=n-1,1,-1
!
      rhoT(i) = rho(temp(i+1))
! 
      temp(i) = temp(i+1) + dx * M0*f3(rhoT(i+1),temp(i+1))     / &
                (24.0*kappa*rhoT(i+1)*temp(i+1)*temp(i+1)*temp(i+1))
!       
   enddo
!
   end subroutine integrate_precursor
!  ============================================== 
!
end program RHD_Shock
