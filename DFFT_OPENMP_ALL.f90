   program dfft
   use omp_lib
   implicit none
   integer::sizet
   double precision,allocatable::x(:),f(:),flush(:)
   double precision::PI,W,sum1,sum2,t0,t1,t2,t3,tempsum1,tempsum2
   integer::i,j,k,n,iterations
   !call omp_set_nested(0)
   call omp_set_num_threads(32)
   sizet=1
do iterations=1,20
   sizet=sizet*2

  allocate(x(sizet)) 
  allocate(f(sizet))     
  allocate(flush(4*sizet))
   PI=4.d0*atan(1.d0)
   W=2*PI/(sizet+1)
 

   call random_number(f)
   x(1:sizet)=f(1:sizet)

   t0=omp_get_wtime()
    !$omp parallel do schedule(static) 
    do i=1,(sizet+1)
       sum1=0.0
       sum2=0.0
       do j=1,sizet
          sum1=x(j)*cos(W*j*i)+sum1
          sum2=x(j)*sin(W*j*i)+sum2
       end do
   end do
   !$omp end parallel do
  t1=omp_get_wtime()
  call random_number(flush)
  call omp_set_num_threads(4)
  !call omp_set_nested(1)
    t2=omp_get_wtime()
    !$omp parallel do schedule(static) 
    do i=1,(sizet+1)
       sum1=0.0
       sum2=0.0
       do j=1,sizet
          sum1=f(j)*cos(W*j*i)+sum1
          sum2=f(j)*sin(W*j*i)+sum2
       end do
   end do
   !$omp end parallel do
   t3=omp_get_wtime()
  print*,sizet,(t3-t2),(t1-t0)
  !print*,'time differe=',(t1-t0)-(t3-t2)
  deallocate(x)
  deallocate(f)
  deallocate(flush)  
end do
end program
