program main
!
!   use mpi
   implicit none
!
#include "mpif.h"
!
   integer(kind=8), parameter                            :: sizeOfMessage = 10
   integer, parameter                                    :: nrOfMessage = 2
   integer, parameter                                    :: io = 0
   logical, parameter                                    :: nonBlocking = .true. 
   logical, parameter                                    :: check = .true. 
   logical, parameter                                    :: allgather = .false. 
   logical, parameter                                    :: wf = .true. 
   integer, save                                         :: nrOfProcs, myRank, request, status
   integer, allocatable, save, dimension(:)              :: req
   integer, allocatable, save, dimension(:)              :: sendReq
   integer, allocatable, save, dimension(:)              :: receiveReq
   integer, allocatable, save, dimension(:,:)            :: combinedArray
   real,    allocatable, save, dimension(:,:)            :: sendArray,receiveArray
   real,    allocatable, save, dimension(:,:,:)          :: hugeGrid
   real                                                  :: elapsed(2),eTime,totalTime,tOld,tStart,dt
   character (len=40)                                    :: rankChar
   integer, allocatable, dimension(:,:)                  :: sendStats,receiveStats,stats
   integer                                               :: i,j,k,ni,nj,nk,n,l
   integer                                               :: sender,receiver
   integer                                               :: ierr,tag
   real                                                  :: x,y,z
   real                                                  :: wt1,wt2,wtOld
!
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr ) 
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nrOfProcs, ierr ) 
!
   dt = 0.0
   tOld = 0.0
   totalTime = 0.0
!
   ni = 1000
   nj = 1000
   nk = 100
!
   allocate(req(2*(nrOfProcs-1)))
   allocate(sendReq(0:nrOfProcs-1))
   allocate(receiveReq(0:nrOfProcs-1))
   allocate(sendArray(sizeOfMessage,nrOfMessage))
   allocate(receiveArray(sizeOfMessage,nrOfMessage*nrOfProcs))
   allocate(combinedArray(sizeOfMessage,0:nrOfProcs-1))
   allocate(sendStats(MPI_STATUS_SIZE,(nrOfProcs-1)))
   allocate(receiveStats(MPI_STATUS_SIZE,(nrOfProcs-1)))
   allocate(stats(MPI_STATUS_SIZE,2*(nrOfProcs-1)))
   allocate(hugeGrid(ni,nj,nk))
!
   totalTime = 0.d0
   tOld = 0.d0
!
   if(myRank.eq.0) then
!
      write(*,*) 'FORTRAN for MPI Routines'
      write(*,*) 
      write(*,*) 'nr of processes:', nrOfProcs
      write(*,*) 'size of message:', sizeOfMessage
      write(*,*) 'size of message [MB]:', 8.d0*real(sizeOfMessage)/1.d6
!
   endif
!
!  Open individual file for myRank
!
   write(rankChar,'(I4.4)'),myRank
   open(io,file='data_'//trim(rankChar)//'.out')
   write(io,*) 'size of message [MB]:', 8.d0*real(sizeOfMessage)/1.d6
   write(io,*) '---------------------------------------------------------------'
!
!  Fill the messages
!
   receiveArray = 0.d0
   do j=1,nrOfMessage
      do i=1,sizeOfMessage
         sendArray(i,j) = dble(myRank*sizeOfMessage+i)
         receiveArray(i,myRank*nrOfMessage+j) = sendArray(i,j)
         if(wf) write(io,*) i,j,sendArray(i,j)
      enddo
   enddo
!
   if(wf) then
   do k=0,nrOfProcs-1
      write(io,*) 'rank',k 
      do j=1,nrOfMessage
        do i=1,sizeOfMessage
         write(io,*) i,j,receiveArray(i,k*nrOfMessage+j)
        enddo
      enddo
   enddo
   endif
!
   req = MPI_REQUEST_NULL
!
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
   call cpu_time(totalTime)
   tOld = totalTime
   tStart = totalTime
   wt1 = MPI_WTIME()
   wtOld = wt1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  let's do it manually non-blocking
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   do n = 0, nrOfProcs-1
   
    receiveReq(n)=MPI_REQUEST_NULL
    sender = n

    if(myRank .eq. sender) then
      do l = 0, nrOfProcs-1
         receiver = l
         sendReq(l)=MPI_REQUEST_NULL
         if(receiver .ne. sender) then
           !write(*,*)'rank',myRank,'starting communication'
           call MPI_ISEND(sendArray,sizeOfMessage*nrOfMessage,MPI_DOUBLE_PRECISION,receiver,1,MPI_COMM_WORLD,sendReq(l),ierr)
           !call MPI_SSEND(sendArray(1),sizeOfMessage,MPI_DOUBLE_PRECISION,receiver,1,MPI_COMM_WORLD,ierr)
           !call MPI_SEND(sendArray(1),sizeOfMessage,MPI_DOUBLE_PRECISION,receiver,1,MPI_COMM_WORLD,ierr)
           !req(n) = sendReq(l)
         end if
       end do
    else
            !!THIS IS RECIEVER
       call MPI_IRECV(receiveArray(:,1+n*nrOfMessage),sizeOfMessage*nrOfMessage,&
            MPI_DOUBLE_PRECISION,sender,1,MPI_COMM_WORLD,receiveReq(n),ierr)
       !req(nrOfProcs-1+n) = receiveReq(n)
    end if
   end do
!
   call cpu_time(totalTime)
   write(io,*) 'time for starting communication:',totalTime-tOld
   tOld = totalTime
!
!  walk throug a huge grid which is independent of communicated data
!
!   call grid_walk()
!
   call cpu_time(totalTime)
   write(io,*) 'time for grid grid loop:',totalTime-tOld
   tOld = totalTime
!
!
   write(*,*)'rank',myRank,'starting MPI_Wait'
   do n = 0, nrOfProcs-1
      if (n .ne. myRank) then
   !====================
   ! WAIT for message
   !====================
   
          call MPI_WAIT(receiveReq(n),status,ierr)
            !! do something with recieved message
   
      end if
   end do
!
   call cpu_time(totalTime)
   write(io,*) 'time for MPI_WAIT(Receive):',totalTime-tOld
   tOld = totalTime
!
   do n = 0, nrOfProcs-1
      if (n .ne. myRank) then
   !====================
   ! WAIT for message
   !====================
   
          call MPI_WAIT(sendReq(n),status,ierr)
            !! do something with recieved message
   
      end if
   end do
!
   call cpu_time(totalTime)
   write(io,*) 'time for MPI_WAIT(Send):',totalTime-tOld
   tOld = totalTime
!
   write(*,*)'rank',myRank,'finished MPI_WAIT'
!
!  Check received messages
!
   do k=0,nrOfProcs-1
      if(wf) write(io,*),'rank',k
      do j=1,nrOfMessage
      do i=1,sizeOfMessage
         l = k*sizeOfMessage+i
         if(check.and.receiveArray(i,nrOfMessage*k+j).ne.dble(l)) then
           write(io,*) 'ERROR in Communication with MPI_ALLGATHER'
           write(*,*) 'ERROR in Communication with MPI_ALLGATHER'
           write(*,*) 'i',i
           write(*,*) 'j',j
           write(*,*) 'l',l
           write(*,*) 'receiveArray(i,j+nrOfMessage*k)',receiveArray(i,j+nrOfMessage*k)
           stop
         endif
         if(wf) write(io,*)i,j,receiveArray(i,nrOfMessage*k+j)
      enddo
      enddo
   enddo
!
   if(check) write(*,*) 'rank',myRank, 'check successfull'
!
   call cpu_time(totalTime)
   write(io,*) 'time for check:',totalTime-tOld
   tOld = totalTime
!
   wt1 = MPI_WTIME()
!
   write(io,*) 'total time for manually non-blocking:',totalTime-tStart
   write(io,*) 'total time for manually non-blocking (wtime):',wt1-wtOld
   write(io,*) '---------------------------------------------------------------'
!
!  Reset messages
!
   receiveArray = 0.d0
   do i=1,sizeOfMessage
      do j=1,nrOfMessage
         sendArray(i,j) = dble(myRank*sizeOfMessage+i)
         receiveArray(i,myRank*nrOfMessage+j) = sendArray(i,j)
      enddo
   enddo
!
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
   call cpu_time(totalTime)
   tOld = totalTime
   tStart = totalTime
   wt1 = MPI_WTIME()
   wtOld = wt1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  now do it blocking with allgather
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
   call MPI_ALLGATHER(sendArray,sizeOfMessage*nrOfMessage,MPI_DOUBLE_PRECISION,      &
                      receiveArray,sizeOfMessage*nrOfMessage,MPI_DOUBLE_PRECISION,   &
                      MPI_COMM_WORLD,ierr)
!
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
   call cpu_time(totalTime)
   write(io,*) 'time for MPI_ALLGATHER:',totalTime-tOld
   tOld = totalTime
!
!  walk throug a huge grid which is independent of communicated data
!
!   call grid_walk()
!
   call cpu_time(totalTime)
   write(io,*) 'time for grid loop:',totalTime-tOld
   tOld = totalTime
!
!  Check received messages
!
   do k=0,nrOfProcs-1
      if(wf) write(io,*),'rank',k
      do j=1,nrOfMessage
      do i=1,sizeOfMessage
         l = k*sizeOfMessage+i
         if(check.and.receiveArray(i,nrOfMessage*k+j).ne.dble(l)) then
           write(io,*) 'ERROR in Communication with MPI_ALLGATHER'
           write(*,*) 'ERROR in Communication with MPI_ALLGATHER'
           write(*,*) 'i',i
           write(*,*) 'j',j
           write(*,*) 'l',l
           write(*,*) 'receiveArray(i,j+nrOfMessage*k)',receiveArray(i,j+nrOfMessage*k)
           stop
         endif
         if(wf) write(io,*)i,j,receiveArray(i,nrOfMessage*k+j)
      enddo
      enddo
   enddo
!
   if(check) write(*,*) 'rank',myRank, 'check successfull'
!
   call cpu_time(totalTime)
   write(io,*) 'time for check:',totalTime-tOld
   tOld = totalTime
!
   wt1 = MPI_WTIME()
!
   call cpu_time(totalTime)
   write(io,*) 'total time using MPI_ALLGATHER (blocking):',totalTime-tStart
   write(io,*) 'total time using MPI_ALLGATHER (blocking) (wtime):',wt1-wtOld
   tOld = totalTime
   write(io,*) '---------------------------------------------------------------'
!
   write(*,*) 'rank',myRank, 'waiting'
!
!  Reset messages
!
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   write(*,*) 'rank',myRank, 'passed Barrier...'
!
   call cpu_time(totalTime)
   write(io,*) '---------------------------------------------------------------'
   write(io,*) 'total runtime:',totalTime
!
   write(*,*) 'rank',myRank, 'closing file...'
   close(io)
!
   write(*,*) 'rank',myRank, 'finalize'
!
   call MPI_FINALIZE(ierr)
!
   contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  large grid wald
!
   subroutine grid_walk()
!
      implicit none
!
      do i=1,ni
        do j=1,nj
          do k=1,nk
            hugeGrid(i,j,k) = 1.0
          enddo 
        enddo 
      enddo 
!
      return
!
   end subroutine grid_walk
!
end program main
