! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module fft_routines
  ! FFT_ROUTINES contains variables and routines related to taking
  ! Fourier transforms of fields using the FFTW 2.x MPI library
  ! 
  ! Modified: 12 September 2007

  use constants, only : pr,pin,int8

  implicit none

  save

  ! FFTW constants
  
  integer FFTW_FORWARD,FFTW_BACKWARD
  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
  
  integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
  parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)
  
  integer FFTW_ESTIMATE,FFTW_MEASURE
  parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
  
  integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
  parameter (FFTW_OUT_OF_PLACE=0)
  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
  
  integer FFTW_THREADSAFE
  parameter (FFTW_THREADSAFE=128)
  
  integer FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
  integer FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
  parameter(FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
  parameter(FFTW_SCRAMBLED_INPUT=8192)
  parameter(FFTW_SCRAMBLED_OUTPUT=16384)

  ! Global Parameters:
  ! USE_WORK = use work array to speed up calculation
  ! LEVEL = plan level
  ! PLANF = used by FFTW for forward transforms
  ! PLANI = used by FFTW for inverse transforms
  ! MX = lower index of x on this process
  ! PX = upper index of x on this process
  ! LNX = length of x on this process
  ! MK = lower index of k on this process (for complex-packed data)
  ! PK = upper index of k on this process (for complex-packed data)
  ! LNK = length of k on this process (for complex-packed data)
  ! LN = length of complex data vector on this process
  ! ...0=...
  ! MXG = array holding lower index of (complex) x for all processes
  ! PXG = array holding upper index of (complex) x for all processes
  ! MKG = array holding lower index of k for all processes
  ! PKG = array holding upper index of k for all processes
  ! MX0G = array holding lower index of (real) x for all processes
  ! PX0G = array holding upper index of (real) x for all processes
  ! KID = rank of process holding particular mode
  ! FACI = normalization factor for IFFT(FFT(.))=a (usually 1/N or similar)
  ! EF = exponential for forward transform (single Fourier mode)
  ! EI = exponential for inverse transform (single Fourier mode)
  ! R = real data array
  ! C = complex data array
  ! D = data array used by FFTW
  ! DX = data array with more intuitive indices
  ! DK = data array with more intuitive indices
  ! W = work array used by FFTW
  ! WK = work array with more intuitive indices
  ! COUNTS_X = MPI message sizes (x distributed)
  ! DISPLS_X = MPI message displacements (x distributed)

  integer(pin),parameter :: use_work=1
  integer(int8) :: planf,plani
  integer(pin) :: level,nx,mx,px,lnx,nk,mk,pk,lnk,ln,mx0,px0,lnx0
  integer(pin),dimension(:),allocatable :: mxg,pxg,mkg,pkg,mx0g,px0g,kid
  real(pr) :: faci
  complex(pr) :: Ef,Ei
  real(pr),dimension(:),allocatable :: R
  complex(pr),dimension(:),allocatable :: C,D,Dx,Dk,W,Wk
  integer(pin),dimension(:),allocatable :: counts_x,displs_x


contains


  subroutine init_fft(mdl,cnv)
    ! INIT_FFT initializes fft_routines variables
    ! 
    ! Modified: 31 July 2007

    use constants, only : one,two,img
    use model, only : model_type
    use convolution, only : convolution_type
    use mpi_routines, only : nprocs,my_rank
    use mpi

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables

    type(model_type),intent(inout) :: mdl
    type(convolution_type),intent(inout) :: cnv

    ! Internal Parameters:
    ! IK = index in k
    ! IP = index for process
    ! IERR = MPI error flag
    ! CD = counter for cumulative displacements

    integer(pin) :: ik,ip,ierr,cd

    ! set data sizes

    nx = mdl%nx
    nk = nx/2+1

    ! initialize FFT variables

    if (mdl%dim==1) then

       Ef = exp(-img*mdl%k*mdl%x(1))
       Ei = exp( img*mdl%k*mdl%x(1))

    else

       ! assign the factor needed to get the proper inverse transform
       
       faci = two/(real(mdl%nx,pr)*real(mdl%ny,pr))
       
       ! plan fourier transforms (required by FFTW)
       
       call plan_fft(mdl)

    end if

    ! set packed array bounds

    call fftw_f77_mpi_local_sizes(planf,lnx,mx,lnk,mk,ln)

    if (lnx==0) then
       ! if not holding x data, set mx=px=-1
       mx = -1
       px = -1
    else
       px = mx+lnx-1
    end if

    if (lnk==0) then
       ! if not holding k data, set mk=pk=-1
       mk = -1
       pk = -1
    else
       pk = mk+lnk-1
    end if

    ! set unpacked array bounds

    lnx0 = 2*lnx
    if (lnx0==0) then
       mx0 = -1
       px0 = -1
    else
       mx0 = 2*mx
       px0 = mx0+lnx0-1
    end if

    ! store necessary data in derived types (shifting indices by 1 to start at 1 instead of 0)

    mdl%mx = mx0+1
    mdl%px = px0+1
    mdl%holds_x = (lnx/=0)

    cnv%nkx = nk
    cnv%mkx = mk+1
    if (pk==nk-2) then
       ! add extra value for Nyquist
       cnv%pkx = pk+2
    else
       cnv%pkx = pk+1
    end if
    cnv%holds_k = (lnk/=0)

    ! MPI information

    ! counts/displacements for scatterv/gatherv

    allocate(counts_x(0:nprocs-1),displs_x(0:nprocs-1))

    call MPI_AllGather(lnx0,1,MPI_INTEGER,counts_x,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

    cd = 0
    do ip = 0,nprocs-1
       displs_x(ip) = cd
       cd = cd+counts_x(ip)
    end do

    ! write data limits and sizes

    !print *
    !print *, my_rank,'total: ',ln
    !print *, my_rank,'x:     ',mx,px,lnx
    !print *, my_rank,'x0:    ',mx0,px0,lnx0
    !print *, my_rank,'k:     ',mk,pk,lnk
    !print *, my_rank,'mdl%x: ',mdl%mx,mdl%px

    ! share data bounds with all processes

    allocate(mxg(0:nprocs-1),pxg(0:nprocs-1))
    allocate(mkg(0:nprocs-1),pkg(0:nprocs-1))
    allocate(mx0g(0:nprocs-1),px0g(0:nprocs-1))

    call MPI_AllGather(mx,1,MPI_INTEGER,mxg,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_AllGather(px,1,MPI_INTEGER,pxg,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_AllGather(mk,1,MPI_INTEGER,mkg,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_AllGather(pk,1,MPI_INTEGER,pkg,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_AllGather(mx0,1,MPI_INTEGER,mx0g,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_AllGather(px0,1,MPI_INTEGER,px0g,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

    ! establish which process holds which Fourier mode (of complex-packed data)

    allocate(kid(0:nk-1))

    do ik = 0,nk-1
       do ip = 0,nprocs-1
          if (mkg(ip)<=ik.and.ik<=pkg(ip)) then
             kid(ik) = ip
             exit
          end if
       end do
    end do

    ! allocate memory to arrays

    if (lnx0/=0) allocate(R(mx0:px0))

    if (lnk/=0) then
       if (pk==nk-2) then
          ! include Nyquist at end
          allocate(C(mk:pk+1))
       else
          allocate(C(mk:pk))
       end if
    end if

    allocate(D(0:ln-1))
    allocate(W(0:ln-1))

    if (lnx/=0) allocate(Dx(mx:px))
    if (lnk/=0) allocate(Dk(mk:pk),Wk(mk:pk))

  end subroutine init_fft


  subroutine destroy_fft
    ! DESTROY_FFT deallocates temporary arrays
    ! 
    ! Modified: 2 August 2007

    implicit none

    ! deallocate temporary arrays

    if (allocated(mxg)) deallocate(mxg) 
    if (allocated(pxg)) deallocate(pxg) 
    if (allocated(mx0g)) deallocate(mx0g) 
    if (allocated(px0g)) deallocate(px0g) 
    if (allocated(mkg)) deallocate(mkg) 
    if (allocated(pkg)) deallocate(pkg) 
    if (allocated(kid)) deallocate(kid) 

    if (allocated(counts_x)) deallocate(counts_x)  
    if (allocated(displs_x)) deallocate(displs_x)  

    if (allocated(R)) deallocate(R)
    if (allocated(C)) deallocate(C) 
    if (allocated(D)) deallocate(D) 
    if (allocated(Dx)) deallocate(Dx) 
    if (allocated(Dk)) deallocate(Dk) 
    if (allocated(W)) deallocate(W) 
    if (allocated(Wk)) deallocate(Wk) 

  end subroutine destroy_fft


  subroutine plan_fft(mdl)
    ! PLAN_FFT plans FFTs, sometimes by running several test problems
    ! 
    ! Modified: 12 September 2007

    use constants, only : real4,real8
    use model, only : model_type
    use io, only : error
    use mpi

    implicit none    

    ! I/O Parameters:
    ! MDL = model variables

    type(model_type),intent(in) :: mdl

    ! FFTW requires that you 'plan' the FFTs before executing them.  It often runs a few tests
    ! to pick the fastest routine given the array dimensions.  There is obviously a trade-off
    ! between the time spent planning and the time saved later, which is controlled by the following 
    ! variables.  From least to most planning are: fftw_estimate, fftw_measure.

    ! assign FFTW plan level

    level = fftw_estimate
    !level = fftw_measure

    ! create forward and backward plans for complex data of half the size

    select case(pr)

    case(real4) ! single-precision transforms

       call error('Single precision FFTs not implemented','plan_fft')

    case(real8) ! double-precision transforms

       select case(mdl%dim)
       case(2)
          call fftw_f77_mpi_create_plan(planf,MPI_COMM_WORLD,nx/2,fftw_forward ,level)
          call fftw_f77_mpi_create_plan(plani,MPI_COMM_WORLD,nx/2,fftw_backward,level)
       case(3)
          call error('2D MPI FFTs not implemented','plan_fft')
       end select

    end select

  end subroutine plan_fft


  subroutine forward_fft(mdl,cnv,S,F) 
    ! FORWARD_FFT forward FFTs a spatial field to Fourier domain
    ! 
    ! Modified: 12 September 2007

    use model, only : model_type
    use convolution, only : convolution_type
    use constants, only : img,zero,half,twopi,zeroc
    use mpi_routines, only : MPI_CPLX_PR
    use mpi

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! S = field in space
    ! F = field in Fourier domain

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(in) :: cnv
    real(pr),dimension(mdl%mx:mdl%px,mdl%ny),intent(in) :: S
    complex(pr),dimension(cnv%mkx:cnv%pkx,cnv%nky),intent(out) :: F

    ! Internal Parameters:
    ! IERR = MPI error flag
    ! I = index
    ! C0 = zero mode
    ! CN = Nyquist mode
    ! STATUS = MPI status

    integer(pin) :: ierr,i,status(MPI_STATUS_SIZE)
    complex(pr) :: C0,CN

    ! special treatment of low dimensional cases

    select case(mdl%dim)
    case(0)
       F = cmplx(S,zero,pr)
       return
    case(1)
       F = Ef*S
       return
    end select

    ! pack data

    if (lnx/=0) then

       R = S(:,1)

       do i = mx,px
          Dx(i) = cmplx(R(2*i),R(2*i+1),pr)
       end do
     
       D(0:lnx-1) = Dx

    end if

    ! transform

    call fftw_f77_mpi(planf,1,D,W,use_work)

    ! move D to Dk, which has more intuitive indices
    
    if (lnk/=0) Dk = D(0:lnk-1)
    
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    ! pass information between different processes to set Wk(i)=Dk(nx/2-i)
    
    if (lnk/=0) call reverse_vector(Dk,Wk)
    
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    if (lnk/=0) then
    
    ! unpack transformed data
       
       ! set Wk(i)=Dk(nx/2-i)*

       Wk = conjg(Wk)

       ! unpack

       do i = mk,pk

          if (i==0) then

             ! zero and Nyquist extracted together (since both are real)

             C(0) = cmplx( &
                  real(half*(Dk(0)+Wk(0))-half*img*(Dk(0)-Wk(0))) , &
                  real(half*(Dk(0)+Wk(0))+half*img*(Dk(0)-Wk(0))) ,pr)

             ! zero

             C0 = cmplx(real(half*(Dk(i)+Wk(i))-half*img*(Dk(i)-Wk(i))),zero,pr)

             ! Nyquist

             CN = cmplx(real(half*(Dk(i)+Wk(i))+half*img*(Dk(i)-Wk(i))),zero,pr)

          else

             C(i) = half*(Dk(i)+Wk(i))-half*img*(Dk(i)-Wk(i))*exp(-twopi*img*real(i,pr)/real(nx,pr))

          end if

       end do

       ! store zero and Nyquist on appropriate processes:
       ! C(0) = C0, C(nk-1) = CN
       
       if (mk==0) then
          
          ! store zero mode
          
          C(0) = C0
          
          ! send Nyquist mode
          
          call MPI_Send(CN,1,MPI_CPLX_PR,kid(nk-2),7777,MPI_COMM_WORLD,ierr)
 
       end if
       
       ! receive Nyquist mode and store it

       if (pk==nk-2) &
            call MPI_Recv(C(nk-1),1,MPI_CPLX_PR,kid(0),7777,MPI_COMM_WORLD,status,ierr)
       
    end if

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! place data in output array

    if (lnk==0) then
       F = zeroc ! need to return something even if not used
    else
       F(:,1) = C
    end if

  end subroutine forward_fft


  subroutine inverse_fft(mdl,cnv,F,S)
    ! INVERSE_FFT inverse FFTs a Fourier-domain field to space
    ! 
    ! Modified: 12 September 2007

    use model, only : model_type
    use convolution, only : convolution_type
    use constants, only : img,zero,half,twopi
    use mpi_routines, only : MPI_CPLX_PR
    use mpi

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! F = field in Fourier domain
    ! S = field in space

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(in) :: cnv
    complex(pr),dimension(cnv%mkx:cnv%pkx,cnv%nky),intent(in) :: F
    real(pr),dimension(mdl%mx:mdl%px,mdl%ny),intent(out) :: S

    ! Internal Parameters:
    ! IERR = MPI error flag
    ! I = index
    ! C0 = zero mode
    ! CN = Nyquist mode
    ! STATUS = MPI status

    integer(pin) :: ierr,i,status(MPI_STATUS_SIZE)
    complex(pr) :: C0,CN

    ! special treatment of low dimensional cases

    select case(mdl%dim)
    case(0)
       S = real(F,pr)
       return
    case(1)
       S = real(Ei*F,pr)
       return
    end select

    if (lnk/=0) then

       ! store F in C

       C = F(:,1)
              
       ! send Nyquist mode
       if (pk==nk-2) &
            call MPI_Send(C(nk-1),1,MPI_CPLX_PR,kid(0),8888,MPI_COMM_WORLD,ierr)
          
       if (mk==0) then
          
          ! receive Nyquist mode
          
          call MPI_Recv(CN,1,MPI_CPLX_PR,kid(nk-2),8888,MPI_COMM_WORLD,status,ierr)
          
          ! store zero and Nyquist modes together
          
          C0 = C(0)
          
          C(0) = cmplx(real(C0),real(CN),pr)
          
       end if
       
    end if

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    ! move C to Dk (note that Nyquist index is never used)
    
    if (lnk/=0) Dk = C(mk:pk)
    
    ! pass information between different processes to set Wk(i)=Dk(nx/2-i)
    
    if (lnk/=0) call reverse_vector(Dk,Wk)
    
    if (lnk/=0) then

       ! pack transformed data
    
       do i = mk,pk
          if (i==0) then
             ! zero mode
             ! Dk(i) holds zero mode (real)
             ! Wk(i) holds Nyquist mode (real)
             Dk(i) = cmplx(real( Dk(i)),zero,pr)
             Wk(i) = cmplx(aimag(Wk(i)),zero,pr)
             Dk(i) = cmplx( &
                  real(half*(Dk(i)+Wk(i))), &
                  real(half*(Dk(i)-Wk(i))), pr)
          else
             ! Wk(i) holds Dk(nx/2-i)*
             Wk(i) = conjg(Wk(i))
             Dk(i) = half*(Dk(i)+Wk(i))+half*img*(Dk(i)-Wk(i))*exp(twopi*img*real(i,pr)/real(nx,pr))
          end if
       end do
       
       D(0:lnk-1) = Dk
       
    end if
       
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! inverse transform

    call fftw_f77_mpi(plani,1,D,W,use_work)
    
    ! unpack transform

    if (lnx/=0) then

       Dx = D(0:lnx-1)

       do i = mx,px
          R(2*i) = real(Dx(i))
          R(2*i+1) = aimag(Dx(i))
       end do
     
       ! multiply by necessary factor
       
       R = R*faci
       
       ! place data in output array
       
    end if

    if (lnx==0) then
       S = zero ! need to return something even if not used
    else
       S(:,1) = R
    end if

  end subroutine inverse_fft


  subroutine reverse_vector(A,B)
    
    use mpi_routines, only : MPI_CPLX_PR
    use mpi

    implicit none

    complex(pr),dimension(mk:pk),intent(in)  :: A
    complex(pr),dimension(mk:pk),intent(out) :: B

    integer(pin) :: i,ierr,status(MPI_STATUS_SIZE)

    ! first half, ascending 

    do i = mk,min(pk,nx/4) ! make sure process holds lower half modes
       
       ! zero,Nyquist treated as special case
       
       if (i==0) then
          B(0) = A(0)
          cycle
       end if
       
       ! don't send what process already has, but do store
       
       if (kid(i)==kid(nx/2-i)) then
          B(i) = A(nx/2-i)
          cycle
       end if
       
       ! send and receive to/from same location:
       ! send A(i) to kid(nx/2-i) and store in B(nx/2-i)
       ! receive A(nx/2-i) from kid(nx/2-i) and store in B(i)

       call MPI_SendRecv(A(i),1,MPI_CPLX_PR, &
            kid(nx/2-i),nx/2-i,B(i),1,MPI_CPLX_PR, &
            kid(nx/2-i),i,MPI_COMM_WORLD,status,ierr)

    end do

    ! upper half, descending

    if (nx/4<=pk) then ! make sure process holds upper half modes

       do i = pk,max(nx/4,mk),-1
          
          ! zero,Nyquist treated as special case
          
          if (i==0) then
             B(0) = A(0)
             cycle
          end if
          
          ! don't send what process already has, but do store
          
          if (kid(i)==kid(nx/2-i)) then
             B(i) = A(nx/2-i)
             cycle
          end if
          
          ! send and receive to/from same location:
          ! send A(i) to kid(nx/2-i) and store in B(nx/2-i)
          ! receive A(nx/2-i) from kid(nx/2-i) and store in B(i)
          
          call MPI_SendRecv(A(i),1,MPI_CPLX_PR, &
               kid(nx/2-i),nx/2-i,B(i),1,MPI_CPLX_PR, &
               kid(nx/2-i),i,MPI_COMM_WORLD,status,ierr)
          
       end do

    end if

  end subroutine reverse_vector


  subroutine test_reverse_vector
    
    use mpi_routines, only : my_rank
    use constants, only : zero,one
    
    implicit none
        
    ! test subroutine reverse_vector: B(0)=A(0), B(1:nk-1)=A(nk-1:1:-1)
    
    complex(pr),dimension(mk:pk) :: A,B
    
    integer(pin) :: i
    
    ! initial data
    
    do i = mk,pk
       A(i) = cmplx(real(i,pr),zero,pr)
    end do
    
    B = -one ! make sure B gets set to proper value

    ! reverse vector

    call reverse_vector(A,B)

    ! output data

    write(6,'(a)') ''
    do i = mk,pk
       write(6,'(a,i6,i6,f6.0,i6,f6.0)') 'reverse_vector:  ',my_rank,i,real(A(i)),nx/2-i,real(B(i))
    end do

  end subroutine test_reverse_vector


  function proc_for_task(t,n,np) result(ip)

    use io, only : error

    implicit none

    ! T = specific task
    ! N = number of tasks
    ! NP = number of processes

    integer(pin),intent(in) :: t,n,np
    integer(pin) :: ip

    integer(pin),dimension(0:np-1) :: l,u

    if (n==nx) then
       ! add 1 since location in x lies between 1 and nx
       l = mx0g+1
       u = px0g+1
    elseif (n==nk) then
       ! add 1 since location in k lies between 1 and nk
       l = mkg+1
       u = pkg+1
    else
       call error('Only x and k are distributed','proc_for_task')
    end if

    do ip = 0,np-1
       if (l(ip)<=t.and.t<=u(ip)) return
    end do

    call error('Task not assigned to any process','proc_for_task')

  end function proc_for_task


  function proc_holds_x()

    implicit none

    logical :: proc_holds_x

    if (lnx==0) then
       proc_holds_x = .false.
    else
       proc_holds_x = .true.
    end if

  end function proc_holds_x


end module fft_routines
