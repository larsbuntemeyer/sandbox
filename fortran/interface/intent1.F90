program xxxx

    implicit none  
    integer :: i,j  
    real    :: x,y 

    interface
       subroutine sub(i)
       implicit none
       integer :: i
       end subroutine sub
    end interface

    i = 9
    call sub(i)  
    print*,i ! will print 7 on all compilers I checked  
end  
subroutine sub(i)  
    implicit none
    integer, intent(in) :: i  
    call sub2(i)  
end  
subroutine sub2(i)  
    implicit none  
    !integer :: i  
    integer :: i  
    i = 7  ! This works since the "intent" information was lost.  
end
