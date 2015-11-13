module help
   contains


    subroutine initializeAngles(sizet,angles)
         use omp_lib
         implicit none
         integer(kind=4)::sizeT,i
         real(kind =8) angles(sizeT),T
         T=2.d0*3.141592653589793D+00/(real(sizeT))
        !$omp parallel shared (sizeT,angles)private (i)
        !$omp do
           do i = 1, sizeT/2
           angles(2*i-1) = cos ( T * real ( i - 1 ) )
           angles(2*i)   = sin (T * real ( i - 1 ))
         end do
       !$omp end do
       !$omp end parallel  
    end subroutine         






    subroutine calculateFFTS(sizet,input,output,angles,sgn)
    use omp_lib
    implicit none
    integer(kind=4) sizet,sectionSize,steps,j
    logical change
    real(kind=8) input(2*sizet), output(2*sizet),angles(sizet),sgn
    steps=int(log(real(sizet))/log(real(2)))
    sectionSize=1
    change=.true.

    call oneheadS(sizet,sectionSize,input(1),input((sizet/2)*2+1),output(1),output(sectionSize*2+1),angles,sgn)


    if (sizet == 2 ) then
    return
    end if

    do j=1,steps-2
        sectionSize=sectionSize*2
        if(change)then
           call oneheadS(sizet,sectionSize,output(1),output((sizet/2)*2+1),input(1),input(sectionSize*2+1),angles,sgn)
           change=.false.
        else
            call oneheadS(sizet,sectionSize,input(1),input((sizet/2)*2+1),output(1),output(sectionSize*2+1),angles,sgn)
            change=.true.
        endif
    end do



   if(change)then
        input(1:2*sizet)=output(1:2*sizet)
   end if


  sectionSize=sizet/2
  call oneheadS(sizet,sectionSize,input(1),input((sizet/2)*2+1),output(1),output(sectionSize*2+1),angles,sgn)
  return
    end subroutine








    

    subroutine calculateFFT(sizet,input,output,angles,sgn)
    use omp_lib
    implicit none
    integer(kind=4) sizet,sectionSize,steps,j
    logical change
    real(kind=8) input(2*sizet), output(2*sizet),angles(sizet),sgn
    steps=int(log(real(sizet))/log(real(2)))
    sectionSize=1
    change=.true.
    
    call onehead(sizet,sectionSize,input(1),input((sizet/2)*2+1),output(1),output(sectionSize*2+1),angles,sgn)


    if (sizet == 2 ) then
    return
    end if

    do j=1,steps-2
        sectionSize=sectionSize*2
        if(change)then
           call onehead(sizet,sectionSize,output(1),output((sizet/2)*2+1),input(1),input(sectionSize*2+1),angles,sgn)
           change=.false.
        else
            call onehead(sizet,sectionSize,input(1),input((sizet/2)*2+1),output(1),output(sectionSize*2+1),angles,sgn)
            change=.true.
        endif
    end do



   if(change)then
        input(1:2*sizet)=output(1:2*sizet)
   end if


  sectionSize=sizet/2
  call onehead(sizet,sectionSize,input(1),input((sizet/2)*2+1),output(1),output(sectionSize*2+1),angles,sgn)
  return
    end subroutine





   subroutine oneheadS(sizet,sectionSize,input1,input2,output1,output2,angles,sgn)
      use omp_lib
      implicit none
      integer(kind=4)sizet,sectionSize,j,a,b,c,d,w,k,sizeofsection,sections
      real(kind=8)input1(sizet),input2(sizet),output1(sizet),output2(sizet),angles(sizet),sgn,tempR,tempI,angle1,angle2
      sections=sectionSize*2
      sizeofsection=sizet/sections
    
       do j=0,sizeofsection-1
                w=j*sectionSize
                a=w
                b=a
                c=j*sections
                d=c
                angle1=angles(w*2+1)
                angle2=angles(w*2+2)

               if ( sgn < 0.0D+00 ) then
                 angle2 = - angle2
               end if


                do k=0,sectionSize-1
                   output1((c+k)*2+1)=input1((a+k)*2+1)+input2((b+k)*2+1)
                   output1((c+k)*2+2)=input1((a+k)*2+2)+input2((b+k)*2+2)

                   tempr=input1((a+k)*2+1)-input2((b+k)*2+1)
                   tempi=input1((a+k)*2+2)-input2((b+k)*2+2)

                   output2((d+k)*2+1)=angle1*tempr-angle2*tempi
                   output2((d+k)*2+2)=angle2*tempr-angle1*tempi
                end do
        end do
      return
   end subroutine



 

   subroutine onehead(sizet,sectionSize,input1,input2,output1,output2,angles,sgn)
      use omp_lib
      implicit none
      integer(kind=4)sizet,sectionSize,j,a,b,c,d,w,k,sizeofsection,sections
      real(kind=8)input1(sizet),input2(sizet),output1(sizet),output2(sizet),angles(sizet),sgn,tempR,tempI,angle1,angle2
      sections=sectionSize*2
      sizeofsection=sizet/sections
     !$omp parallel shared ( input1, input2, output1, output2, sizeofsection, sectionSize, sections, sgn, angles ) private (tempr,tempi,j,a,b,c,d,w,k,angle1,angle2)
     !$omp do 
       do j=0,sizeofsection-1
                
                w=j*sectionSize
                a=w
                b=a
                c=j*sections
                d=c
                angle1=angles(w*2+1)
                angle2=angles(w*2+2)

               if ( sgn < 0.0D+00 ) then
                 angle2 = - angle2
               end if


                do k=0,sectionSize-1
                   output1((c+k)*2+1)=input1((a+k)*2+1)+input2((b+k)*2+1)
                   output1((c+k)*2+2)=input1((a+k)*2+2)+input2((b+k)*2+2)
        
                   tempr=input1((a+k)*2+1)-input2((b+k)*2+1)
                   tempi=input1((a+k)*2+2)-input2((b+k)*2+2)

                   output2((d+k)*2+1)=angle1*tempr-angle2*tempi
                   output2((d+k)*2+2)=angle2*tempr-angle1*tempi
                end do
!        print*,'done=',j
        end do

     !$omp end do
     !$omp end parallel
      return
   end subroutine





   subroutine onehead2(sizet,sectionSize,input1,input2,output1,output2,angles,sgn)
      use omp_lib
      implicit none
      integer(kind=4)sizet,sectionSize,j,a,b,c,d,w,k,sizeofsection,sections,index1,index2
      real(kind=8)input1(sizet),input2(sizet),output1(sizet),output2(sizet),angles(sizet),sgn,tempR,tempI,angle1,angle2
      sections=sectionSize*2
      sizeofsection=sizet/sections
     !$omp parallel shared ( input1, input2, output1, output2, sizeofsection,sectionSize, sections, sgn, angles ) private(tempr,tempi,j,a,b,c,d,w,k,angle1,angle2)
     !$omp do 
       do k=0,sectionSize-1

                w=j*sectionSize
                a=w
                b=a
                c=j*sections
                d=c
                angle1=angles(w*2+1)
                angle2=angles(w*2+2)
                index1=(a+k)*2
                index2=(c+k)*2
               if ( sgn < 0.0D+00 ) then
                 angle2 = - angle2
               end if


                do j=0,sizeofsection-1
                   output1(index2+1)=input1(index1+1)+input2(index1+1)
                   output1(index2+2)=input1(index1+2)+input2(index1+2)

                   tempr=input1(index1+1)-input2(index1+1)
                   tempi=input1(index1+2)-input2(index1+2)

                   output2(index2+1)=angle1*tempr-angle2*tempi
                   output2(index2+2)=angle2*tempr-angle1*tempi
                end do
!        print*,'done=',j
        end do

     !$omp end do
     !$omp end parallel
      return
   end subroutine





end module





program FFT
use omp_lib
use help


implicit none

real(kind=8),allocatable,dimension(:) :: input,output,angles,copy,flush
real(kind =8) t0,t1,sgn,temp0,temp1,seed,fnml,ggl,flopsP,flops,t2,t3
integer(kind=4) icase,nits,sizet,i,j,max_size,nt
logical first,switch
integer(kind=4), parameter::iterations=32
nits=10000
nt=32
first=.true.
seed  = 331.0D+00
sizet=1
call omp_set_num_threads(nt)
print*,'Number of threads=',omp_get_max_threads()

do i=1,iterations
   sizet=sizet*2
   max_size=sizet*2
   

   !!!!!!!!!!!!!!!Allocation!!!!!!!!!!!!!!!!
    allocate ( input(1:max_size) )
    allocate ( output(1:max_size) )
    allocate ( copy(1:max_size) )
    allocate ( angles(1:sizet) )
    allocate ( flush(1: 6*sizet))
   ! print*,'Allocation completed'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          do j=1,max_size-1,2
               ! temp0=randomNum(seed)
               ! temp1=randomNum(seed)
                input(j)=temp0
                copy(j)=temp0
                input(j+1)=temp1
                copy(j+1)=temp1
         end do    
        
      call initializeAngles(sizet,angles)
   t0=omp_get_wtime()
        sgn=+1.d0
        call calculateFFT(sizet,input,output,angles,sgn)
        sgn=-1.d0
        call calculateFFT(sizet,output,input,angles,sgn)
   t1=omp_get_wtime()

call random_number(flush)  
!call omp_set_num_threads(4*nt)
   t2=omp_get_wtime()
        sgn=+1.d0
        call calculateFFT(sizet,copy,output,angles,sgn)
        sgn=-1.d0
        call calculateFFT(sizet,output,copy,angles,sgn)
   t3=omp_get_wtime()

   !!!!!!!!!!!!!!!!DeAllocation!!!!!!!!!!!

    deallocate(input)
    deallocate(output)
    deallocate(angles)
    deallocate(copy)
    deallocate(flush)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*,sizet,t1-t0,t3-t2
end do !iterations end here


end program






