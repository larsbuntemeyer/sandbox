program main
   !
   implicit none
   !
   integer,parameter :: n=2
   real              :: arrayA(n,n,n,3)
   integer           :: i,j,k,l
   integer,parameter :: ijk=n*n*n
   !
   write(*,*),'---------------'
   do k=1,n
      do j=1,n
         do i=1,n
            arrayA(i,j,k,1) = i
            arrayA(i,j,k,2) = j
            arrayA(i,j,k,3) = k
            write(*,*),arrayA(i,j,k,:)
         enddo
      enddo
   enddo 
   write(*,*),'---------------'
   !
   call sub(arrayA)
   !
   contains
   !
   subroutine sub(arrayB)
      !
      implicit none
      !
      real, intent(inout) :: arrayB(ijk,3)
      integer             :: i
      !
      write(*,*),'---------------'
      do i=1,ijk
         write(*,*),arrayB(i,:)
      enddo
      write(*,*),'---------------'
      write(*,*),'---------------'
      write(*,*),arrayA
      write(*,*),'---------------'
      write(*,*),arrayB
      write(*,*),'---------------'
      write(*,*),'---------------'
      !
      !
   end subroutine
   !
end program main

