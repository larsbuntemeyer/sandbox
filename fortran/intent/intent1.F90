program xxxx

    implicit none  
    integer :: i,j  
    real    :: x,y 

  !  interface
  !     subroutine sub(i)
  !     implicit none
  !     integer :: i
  !     end subroutine sub
  !  end interface

    i = 9
    call sub(i)  
    call sub3(i,i)  
    print*,i ! will print 7 on all compilers I checked  
end  
subroutine sub(i)  
    implicit none
    integer, intent(in) :: i  
    call sub2(i)  
end  
subroutine sub2(i)  
    implicit none  
    integer, intent(inout) :: i  
    i = 7  ! This works since the "intent" information was lost.  
end
subroutine sub3(i,j)  
    implicit none  
    integer, intent(in)  :: i  
    integer, intent(out) :: j  
end
