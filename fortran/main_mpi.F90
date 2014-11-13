program main
!
   implicit none
!
   integer, parameter                 :: sizeOfMessage = 4
   integer, parameter                 :: io = 0
   integer, allocatable, dimension(:) :: req
   real, allocatable, dimension(:,:)  :: receiveArray
   real, allocatable, dimension(:)    :: sendArray
   character (len=40)                 :: rankChar
   integer        :: i,j
   integer        :: ierr,rank,nrOfProcs,tag
   logical        :: nonBlocking = .false. 
   !MPI_REQUEST    :: request
!
#include "mpif.h"
!
   write (*,*) 'FORTRAN for MPI Routines'
!
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr ) 
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nrOfProcs, ierr ) 
!
   allocate(req(nrOfProcs-1))
   allocate(sendArray(sizeOfMessage))
   allocate(receiveArray(0:nrOfProcs,sizeOfMessage))
!
   write (*,*)
   if(rank.eq.0) write (*,*) 'nr of processes:', nrOfProcs
   write (*,*) 'my rank:', rank
!
   write(rankChar,'(I4.4)'),rank
   open(io,file='data_'//trim(rankChar)//'.out')
!
   write(io,*) '------------------------------'
   write(io,*) 'sendArray before Communication'
   write(io,*) '------------------------------'
!
   do i=1,sizeOfMessage
      sendArray(i) = real((rank+1)*i)
      receiveArray(rank,i) = sendArray(i)
      write(io,*) i,sendArray(i)
   enddo
!
!  do it the usual way
!
   if(not(nonBlocking)) then
       do j=0,nrOfProcs-1
          if(j.ne.rank) then
             write(*,*),'Rank',rank,' sending to', j
             call MPI_SEND(sendArray,sizeOfMessage,MPI_DOUBLE_PRECISION,j,1,MPI_COMM_WORLD,ierr) 
          endif
       enddo
!
!      receive message 
!
       do j=0,nrOfProcs-1
          if(j.ne.rank) then
             write(*,*),'Rank',rank,' receiving from', j
             call MPI_RECV(receiveArray(j,:),sizeOfMessage,MPI_DOUBLE_PRECISION,j,1,MPI_COMM_WORLD,ierr) 
          endif
       enddo
!
!  do it the non-blocking way
!
   else 
       do j=0,nrOfProcs-1
          if(j.ne.rank) then
             write(*,*),'Rank',rank,' sending to', j
             call MPI_SEND(sendArray,sizeOfMessage,MPI_DOUBLE_PRECISION,j,1,MPI_COMM_WORLD,ierr) 
          endif
       enddo
!
!      receive message 
!
       do j=0,nrOfProcs-1
          if(j.ne.rank) then
             write(*,*),'Rank',rank,' receiving from', j
             call MPI_RECV(receiveArray(j,:),sizeOfMessage,MPI_DOUBLE_PRECISION,j,1,MPI_COMM_WORLD,ierr) 
          endif
       enddo
   endif
!
   write(io,*) '--------------------------------'
   write(io,*) 'receiveArray after Communication'
   write(io,*) '--------------------------------'
!
   do j=0,nrOfProcs-1
      write(io,*),'rank',j
      do i=1,sizeOfMessage
         write(io,*) i,receiveArray(j,i)
      enddo
   enddo
!
   close(io)
!
   call MPI_FINALIZE(ierr)
   if(ierr==0) write(*,*) 'rank',rank,'finalized'
!  ----------------------------------------------
   contains
!  ----------------------------------------------
!
!  ============================================== 
   subroutine test(a,b,c)
!  ============================================== 

   implicit none

   real, intent(in)    :: a,b
   real, intent(out)   :: c

   write(*,*) 'writing from subroutine test:',c
   c = a*b*c
   write(*,*) 'writing from subroutine test:',c

   end subroutine test
!  ============================================== 
!
end program main
