PROGRAM test_transfer
integer :: x = 2143289344
print *, transfer(x, 1.0)    ! prints "NaN" on i686
END PROGRAM
