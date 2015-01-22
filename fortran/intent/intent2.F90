         program intent_gotcha
         type mytype
           integer :: x
           real :: y
         end type mytype

         type (mytype) :: a
         a%x = 1 ; a%y = 2.
         print *, a
         call assign(a)
! a%y COULD BE UNDEFINED HERE
         print *, a

         contains

         subroutine assign(this)
         type (mytype), intent (out) :: this
! THIS IS THE WRONG WAY
         this%x = 2
         end subroutine assign

!         subroutine assign(this)
!         type (mytype), intent (out) :: this
!! THIS IS THE RIGHT WAY
!         this%x = 2 ; this%y = 2.
!         end subroutine assign
         end program intent_gotcha

