   program dfft
   use omp_lib
   use mpi
   implicit none
   integer::sizet,sectionsize,index1,index2
   double precision,allocatable::x(:),f(:)
   double precision::PI,W,sum1,sum2,t0,t1,t2,t3,tempsum1,tempsum2
   integer::i,j,k,n,iterations
   integer, parameter :: dp = mpi_double_precision, comm=mpi_comm_world
   integer :: p, rank, status(mpi_status_size),ierror,tag=1
  call mpi_init(ierror)
  call mpi_comm_size(comm, p, ierror)
  call mpi_comm_rank(comm, rank, ierror)
   call omp_set_nested(0)
 !print*,'number of processes=',p
   sizet=1
do iterations=1,25
   sizet=sizet*2
   sectionsize=sizet/p
   call mpi_Bcast(sectionsize,1,mpi_integer,0,comm,ierror)
  allocate(x(sizet)) 
  allocate(f(sizet))     
   PI=4.d0*atan(1.d0)
   W=2*PI/(sizet+1)
 
   if (rank .eq.0)then
   call random_number(f)
   x(1:sizet)=f(1:sizet)
   endif
   !if(rank .eq. 0)then
   !     t0=omp_get_wtime()
   !endif
   call mpi_bcast(x,sizet,dp,0,comm,ierror)  
   call mpi_barrier(comm,ierror)    
   if(rank .eq. 0)then
        t0=omp_get_wtime()
   endif   
    
   if(sectionsize*rank .le. sizet)then
        index1=sectionsize*rank
        index2=index1+sectionsize

    call omp_set_num_threads(16)
    !print*,rank,index1,index2
    !$omp parallel do schedule(auto) 
    do i=index1+1,index2
       ! print*,rank,i,omp_get_thread_num()
       sum1=0.0
       sum2=0.0

       do j=1,sizet
          sum1=x(j)*cos(W*j*i)+sum1
          sum2=x(j)*sin(W*j*i)+sum2
       end do
   !print*,rank,i
   end do

   !$omp end parallel do
   !print*,rank,index1,index2,'done'
   call mpi_gather(x(index1),sectionsize,dp,x,sectionsize,dp,0,comm,ierror)
   !call mpi_barrier(comm,ierror)
  if(rank .eq. 0)then
        t1=omp_get_wtime()
        print*,sizet,(t1-t0)
   endif
end if
   !t3=omp_get_wtime()
 ! print*,sizet,(t1-t0)
  !print*,'time differe=',(t1-t0)-(t3-t2)
  deallocate(x)
  deallocate(f)
    
end do
end program
