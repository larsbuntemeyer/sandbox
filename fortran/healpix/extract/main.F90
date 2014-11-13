module healpix

  implicit none


  contains

  subroutine pix2ang_ring_external  (nside, ipix, theta, phi)

    implicit none

    INTEGER, INTENT(IN)  :: nside
    INTEGER, INTENT(IN)  :: ipix
    double precision,     INTENT(OUT) :: theta, phi

    INTEGER ::  nl2, nl4, iring, iphi
    INTEGER ::  npix, ncap, ip
    double precision ::  fodd, dnside
    INTEGER, PARAMETER :: ns_max4=8192     ! 2^13
    INTEGER, PARAMETER :: ns_max=ns_max4 ! largest nside available
    double precision, parameter :: half = 0.5d0
    double precision, parameter :: PI = 3.14159265359d0
    double precision, parameter :: HALFPI = 1.57079632679d0 
    !real(kind=dp), parameter :: one  = 1.000000000000000_dp
    !real(kind=dp), parameter :: three = 3.00000000000000_dp
    double precision, parameter :: threehalf = 1.5d0
    character(len=*), parameter :: code = "pix2ang_ring"
    !-----------------------------------------------------------------------
    if (nside > ns_max4) then
       print*,code,"> nside out of range"
       stop
    endif
    npix = nside2npix(nside)       ! total number of points
    if (ipix <0 .or. ipix>npix-1) then 
       print*,code,"> ipix out of range"
       stop
    endif

    nl2  = 2*nside
    ncap = nl2*(nside-1) ! points in each polar cap, =0 for nside =1
    dnside = real(nside)

    if (ipix < ncap) then ! North Polar cap -------------

!       iring = nint( sqrt( (ipix+1) * half ), kind=MKD) ! counted from North pole
       iring = (cheap_isqrt(2*ipix+2) + 1)/2
       iphi  = ipix - 2*iring*(iring - 1)

!       theta = ACOS( one - (iring/dnside)**2 / three )
       theta = 2.d0 * asin(iring / (sqrt(6.d0)*dnside))
       phi   = (real(iphi) + half) * HALFPI/iring

    elseif (ipix < npix-ncap) then ! Equatorial region ------

       ip    = ipix - ncap
       nl4   = 4*nside
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = iand(ip, nl4-1)

       fodd  = half * ( iand(iring+nside+1,1) )  ! 0 if iring+nside is odd, 1/2 otherwise
       theta = ACOS( (nl2 - iring) / (threehalf*dnside) )
       phi   = (real(iphi) + fodd) * HALFPI / dnside

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix
!       iring = nint( sqrt( ip * half ), kind=MKD)     ! counted from South pole
       iring = (cheap_isqrt(2*ip) + 1) / 2
       iphi  = 2*iring*(iring + 1) - ip

!       theta = ACOS( (iring/dnside)**2 / three  - one)
       theta = PI - 2.d0 * asin(iring / (sqrt(6.d0)*dnside))
       phi   = (real(iphi) + half) * HALFPI/iring

    endif

    return

    contains


  function nside2npix(nside) result(npix_result)
    !=======================================================================
    ! given nside, returns npix such that npix = 12*nside^2
    !  nside should be a power of 2 smaller than ns_max
    !  if not, -1 is returned
    ! EH, Feb-2000
    ! 2009-03-04: returns i8b result, faster
    !=======================================================================
    INTEGER             :: npix_result
    INTEGER, INTENT(IN) :: nside

    INTEGER :: npix
    CHARACTER(LEN=*), PARAMETER :: code = "nside2npix"
    !=======================================================================

    npix = (12*nside)*nside
    if (nside < 1 .or. nside > ns_max .or. iand(nside-1,nside) /= 0) then
       print*,code,": Nside=",nside," is not a power of 2."
       npix = -1
    endif
    npix_result = npix

    return
  end function nside2npix

  function cheap_isqrt(lin) result (lout)
    integer, intent(in) :: lin
    integer :: lout, diff
    real :: dout, din
    lout = floor(sqrt(dble(lin))) ! round-off error may offset result
    diff = lin - lout*lout ! test Eq (1)
    if (diff <0)      lout = lout - 1
    if (diff >2*lout) lout = lout + 1
    return
  end function cheap_isqrt


  end subroutine pix2ang_ring_external


end module healpix




program main
!
   use pix_tools
   use healpix
!
   implicit none
!
   real :: x,y,z,a,b
   double precision :: theta,phi
   double precision :: thetad,phid
   integer :: nside,ipix,npix
!
   write (*,*) 'Simple Healpix test'
!
!
!
   write(*,*) ''
   write(*,*) 'nside?'
   read(*,*) nside
   npix = 12*nside**2
!
   do ipix=0,npix-1
      call pix2ang_ring(nside,ipix,thetad,phid)
      call pix2ang_ring_external(nside,ipix,theta,phi)
      write(*,*) '----------Original------------------'
      write(*,*) 'ipix',ipix
      write(*,*) 'theta',thetad/3.14159
      write(*,*) 'phi',phid/3.14159
      write(*,*) 'theta',thetad/3.14159*180.0
      write(*,*) 'phi',phid/3.14159*180.0
      write(*,*) 'cos(theta)',cos(thetad)
      write(*,*) '----------External------------------'
      write(*,*) 'ipix',ipix
      write(*,*) 'theta',theta/3.14159
      write(*,*) 'phi',phi/3.14159
      write(*,*) 'theta',theta/3.14159*180.0
      write(*,*) 'phi',phi/3.14159*180.0
      write(*,*) 'cos(theta)',cos(theta)
      if(abs(thetad-theta)/thetad>1.e-8) stop
      if(phi>0.d0.and.abs(phid-phi)/phid>1.e-8) stop
   enddo
   write(*,*) ''
   write(*,*) '------------------------------------'
   write(*,*) 'Total number of pixels:',npix
   write(*,*) ''
   write(*,*) 'DONE'
   write(*,*) ''
!

end program main
