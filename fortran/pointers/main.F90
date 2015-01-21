program main
!
   implicit none
!
   real :: x,y,z,a,b
   integer :: i,j,k
   integer, parameter       :: n=100
   integer, parameter       :: l=10
   real, dimension(n,n,n)   :: grid
   integer, allocatable, dimension(:)        :: integerList
   integer, allocatable, dimension(:,:,:,:)  :: listOfIntegerLists
!
   write (*,*) 'FORTRAN test program'
!
   allocate(integerList(l))
   allocate(listOfIntegerLists(n,n,n,l))
!
   do i=1,l
      integerList(i) = i 
   enddo
! 
   do k=1,n
   do j=1,n
   do i=1,n
      listOfIntegerLists(i,j,k,:) = integerList
   enddo
   enddo
   enddo
!
   integerList = listOfIntegerLists(1,1,1,:)
!
   write(*,*), integerList
   call writeIntegerList(listOfIntegerLists(1,1,1,:))
!
!  ----------------------------------------------
   contains
!  ----------------------------------------------
!  ============================================== 
   subroutine workOnLargeArray(a)
!  ============================================== 

   implicit none

   real, dimension(:,:,:), intent(inout)  :: a

   write(*,*) 'writing from subroutine workOnLargeArray:'

   end subroutine workOnLargeArray
!  ============================================== 
   subroutine writeIntegerList(list)
!  ============================================== 

   implicit none

   integer, dimension(:), intent(in)  :: list

   write(*,*) 'integerList:',list

   end subroutine writeIntegerList
!  ============================================== 

end program main
