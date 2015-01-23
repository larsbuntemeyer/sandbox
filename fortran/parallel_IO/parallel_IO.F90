! example of parallel MPI write into a single file, in Fortran 
PROGRAM main 
    ! Fortran 90 users can (and should) use 
    !     use mpi 
    ! instead of include 'mpif.h' if their MPI implementation provides a 
    ! mpi module. 
    !
    implicit none
    !
    include 'mpif.h' 
     
    integer :: ierr, i, myrank, BUFSIZE, thefile 
    integer, allocatable :: buf(:) 
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer :: count
    real    :: t1,t2,time,total_time

    total_time = 0.0
    
    call cpu_time(t1) 
    call MPI_INIT(ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, count)  
    call cpu_time(t2)
    if(myrank==0) write(*,*) 'time for MPI_INIT:',t2-t1
     
    BUFSIZE = int(960000000./count) 

    if(myrank==0) then
       write(*,*) 'count:',count
       write(*,*) 'BUFSIZE:',BUFSIZE
    endif

    call cpu_time(t1)
    allocate(buf(BUFSIZE))
    call cpu_time(t2)
    if(myrank==0) write(*,*) 'time for allocatio:',t2-t1
 
    call cpu_time(t1)
    do i = 0, BUFSIZE 
        buf(i) = myrank * BUFSIZE + i
    enddo 
    call cpu_time(t2)
    if(myrank==0) write(*,*) 'time for buffer filling:',t2-t1

    call cpu_time(t1)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'testfile', & 
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                       MPI_INFO_NULL, thefile, ierr) 
    call cpu_time(t2)
    if(myrank==0) write(*,*) 'time for MPI_FILE_OPEN:',t2-t1

    ! assume 4-byte integers 
    disp = myrank * BUFSIZE * 4 
    call cpu_time(t1)
    call MPI_FILE_SET_VIEW(thefile, disp, MPI_INTEGER, & 
                           MPI_INTEGER, 'native', & 
                           MPI_INFO_NULL, ierr) 
    call cpu_time(t2)
    if(myrank==0) write(*,*) 'time for MPI_FILE_SET_VIEW:',t2-t1

    call cpu_time(t1)
    call MPI_FILE_WRITE(thefile, buf, BUFSIZE, MPI_INTEGER, & 
                        MPI_STATUS_IGNORE, ierr) 
    call cpu_time(t2)
    if(myrank==0) write(*,*) 'time for MPI_FILE_WRITE:',t2-t1

    call cpu_time(t1)
    call MPI_FILE_CLOSE(thefile, ierr) 
    call cpu_time(t2)
    if(myrank==0) write(*,*) 'time for MPI_FILE_CLOSE:',t2-t1
    call MPI_FINALIZE(ierr) 

    call cpu_time(t1)
    if(myrank==0) write(*,*) 'total time:',t1
 
END PROGRAM main 
