program main
!
   use pix_tools
!
   implicit none
!
   real :: x,y,z,a,b
   real*8 :: theta,phi
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
      call pix2ang_ring(nside,ipix,theta,phi)
      write(*,*) '------------------------------------'
      write(*,*) 'ipix',ipix
      write(*,*) 'theta',theta/3.14159
      write(*,*) 'phi',phi/3.14159
      write(*,*) 'theta',theta/3.14159*180.0
      write(*,*) 'phi',phi/3.14159*180.0
      write(*,*) 'cos(theta)',cos(theta)
   enddo
   write(*,*) ''
   write(*,*) '------------------------------------'
   write(*,*) 'Total number of pixels:',npix
   write(*,*) ''
   write(*,*) 'DONE'
   write(*,*) ''
!

end program main


