program main

  !implicit none

  ! a=10
  ! b=20

  call implicit_decl(a,b,c,i,j,k,l,m,n,x,y,z)

  contains

     subroutine implicit_decl(a,b,c,i,j,k,l,m,n,x,y,z)

     write(*,*) 'a,b,c',a,b,c
     write(*,*) 'i,j,k',i,j,k
     write(*,*) 'l,m,n',l,m,n
     write(*,*) 'x,y,z',x,y,z
 
     end subroutine implicit_decl

end program main
