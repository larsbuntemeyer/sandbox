module test_module
!
   implicit none
!
   contains
!
   subroutine testA()
      implicit none
      write(*,*) "testA is called..."
   end subroutine testA
!
end module test_module
