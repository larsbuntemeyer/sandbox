module commondat
	implicit none
	real, public :: a,b,c
	data a,b,c / -1.0, -1.0, -1.0/
contains
	subroutine outdat()
	write(6,'(1x,"a,b,c=",3f7.2)') a,b,c
	end subroutine outdat
end module commondat

program tcommon
	use commondat
	implicit none
	call suba()
	call subb()
	call subc()
	call outdat()
end program tcommon

subroutine suba()
! _X(USE the module ... it handles all)
! _X(the necessary declarations)
	use commondat, only : x=>a
	implicit none
	real :: a = 1.0
	x = a
	return
end subroutine suba

subroutine subb()
! _X(ONLY use the variable we need ...)
! _X(even rename it, all other are 'invisible')
	use commondat, only : x=>b
	implicit none
	real :: a = 2.0
	x = a
	return
end subroutine subb

subroutine subc()
	use commondat, only : x=>c
! _X(note the rename syntax is same as)
! _X(pointer assignment)
	implicit none
	real :: a = 3.0
	x = a
	return
end subroutine subc
