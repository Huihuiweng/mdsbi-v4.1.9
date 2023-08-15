! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module fft_routines
  ! FFT_ROUTINES contains variables and routines related to taking
  ! Fourier transforms of fields using the FFTW 2.x MPI library with
  ! complex transforms (sacrificing a factor of two in wavenumber space)
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
  ! MK = lower index of k on this process
  ! PK = upper index of k on this process
  ! LNK = length of k on this process
  ! LN = length of complex data vector on this process
  ! MXG = array holding lower index of x for all processes
  ! PXG = array holding upper index of x for all processes
  ! MKG = array holding lower index of k for all processes
  ! PKG = array holding upper index of k for all processes
  ! FACI = normalization factor for IFFT(FFT(.))=a (usually 1/N or similar)
  ! EF = exponential for forward transform (single Fourier mode)
  ! EI = exponential for inverse transform (single Fourier mode)
  ! D = data array used by FFTW
  ! W = work array used by FFTW
  ! COUNTS_X = MPI message sizes (x distributed)
  ! DISPLS_X = MPI message displacements (x distributed)

  integer(pin),parameter :: use_work=1
  integer(int8) :: planf,plani
  integer(pin) :: level,nx,mx,px,lnx,nk,mk,pk,lnk,ln
  integer(pin),dimension(:),allocatable :: mxg,pxg,mkg,pkg
  real(pr) :: faci
  complex(pr) :: Ef,Ei
  complex(pr),dimension(:),allocatable :: D,W
  integer(pin),dimension(:),allocatable :: counts_x,displs_x


contains


  subroutine init_fft(mdl,cnv)
    ! INIT_FFT initializes fft_routines variables
    ! 
    ! Modified: 1 August 2007

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
    ! IP = index for process
    ! NP = number of processes
    ! MYID = rank of current process
    ! IERR = MPI error flag
    ! CD = counter for cumulative displacements

    integer(pin) :: ip,np,myid,ierr,cd

    ! set data sizes

    nx = mdl%nx
    nk = mdl%nx

    ! initialize FFT variables

    if (mdl%dim==1) then

       Ef = exp(-img*mdl%k*mdl%x(1))
       Ei = exp( img*mdl%k*mdl%x(1))

    else

       ! assign the factor needed to get the proper inverse transform
       
       faci = one/(real(mdl%nx,pr)*real(mdl%ny,pr))
       
       ! plan fourier transforms (required by FFTW)
       
       call plan_fft(mdl)

    end if

    ! set packed array bounds

    call fftw_f77_mpi_local_sizes(planf,lnx,mx,lnk,mk,ln)

    if (lnx==0) then
       ! if not holding x data, set mx=px=-1
       mx = -1
       px = -1
       print *, myid,'holds no x data (use fewer processors)'
    else
       px = mx+lnx-1
    end if

    if (lnk==0) then
       ! if not holding k data, set mk=pk=-1
       mk = -1
       pk = -1
       print *, myid,'holds no k data (use fewer processors)'
    else
       pk = mk+lnk-1
    end if

    ! store necessary data in derived types (shifting indices by 1 to start at 1 instead of 0)

    mdl%mx = mx+1
    mdl%px = px+1
    mdl%holds_x = (lnx/=0)

    cnv%nkx = nk
    cnv%mkx = mk+1
    cnv%pkx = pk+1
    cnv%holds_k = (lnk/=0)

    ! MPI information

    np = nprocs
    myid = my_rank

    ! counts/displacements for scatterv/gatherv

    allocate(counts_x(0:np-1),displs_x(0:np-1))

    call MPI_AllGather(lnx,1,MPI_INTEGER,counts_x,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

    cd = 0
    do ip = 0,np-1
       displs_x(ip) = cd
       cd = cd+counts_x(ip)
    end do

    ! write data limits and sizes

    !print *
    !print *, myid,'total: ',ln
    !print *, myid,'x:     ',mx,px,lnx
    !print *, myid,'k:     ',mk,pk,lnk
    !print *, myid,'mdl%x: ',mdl%mx,mdl%px
    !print *, myid,'cnv%k: ',cnv%mkx,cnv%pkx

    ! share data bounds with all processes

    allocate(mxg(0:np-1),pxg(0:np-1))
    allocate(mkg(0:np-1),pkg(0:np-1))

    call MPI_AllGather(mx,1,MPI_INTEGER,mxg,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_AllGather(px,1,MPI_INTEGER,pxg,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_AllGather(mk,1,MPI_INTEGER,mkg,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_AllGather(pk,1,MPI_INTEGER,pkg,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

    ! allocate memory to arrays

    allocate(D(0:ln-1))
    allocate(W(0:ln-1))

  end subroutine init_fft


  subroutine destroy_fft
    ! DESTROY_FFT deallocates temporary arrays
    ! 
    ! Modified: 2 August 2007

    implicit none

    ! deallocate temporary arrays

    if (allocated(mxg)) deallocate(mxg) 
    if (allocated(pxg)) deallocate(pxg) 
    if (allocated(mkg)) deallocate(mkg) 
    if (allocated(pkg)) deallocate(pkg) 

    if (allocated(counts_x)) deallocate(counts_x)  
    if (allocated(displs_x)) deallocate(displs_x)  

    if (allocated(D)) deallocate(D) 
    if (allocated(W)) deallocate(W) 

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

    !level = fftw_estimate
    level = fftw_measure

    ! create forward and backward plans for complex data of half the size

    select case(pr)

    case(real4) ! single-precision transforms

       call error('Single precision FFTs not implemented','plan_fft')

    case(real8) ! double-precision transforms

       select case(mdl%dim)
       case(2)
          call fftw_f77_mpi_create_plan(planf,MPI_COMM_WORLD,nx,fftw_forward ,level)
          call fftw_f77_mpi_create_plan(plani,MPI_COMM_WORLD,nx,fftw_backward,level)
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
    use constants, only : zero,zeroc

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

    ! special treatment of low dimensional cases

    select case(mdl%dim)
    case(0)
       F = cmplx(S,zero,pr)
       return
    case(1)
       F = Ef*S
       return
    end select

    ! store data from input array

    if (lnx/=0) D(0:lnx-1) = cmplx(S(:,1),zero,pr)

    ! transform

    call fftw_f77_mpi(planf,1,D,W,use_work)

    ! unpack data in output array

    if (lnk==0) then
       F = zeroc ! need to return something even if not used
    else
       F(:,1) = D(0:lnk-1)
    end if

  end subroutine forward_fft


  subroutine inverse_fft(mdl,cnv,F,S)
    ! INVERSE_FFT inverse FFTs a Fourier-domain field to space
    ! 
    ! Modified: 12 September 2007

    use model, only : model_type
    use convolution, only : convolution_type
    use constants, only : img,zero,half,twopi

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

    ! special treatment of low dimensional cases

    select case(mdl%dim)
    case(0)
       S = real(F,pr)
       return
    case(1)
       S = real(Ei*F,pr)
       return
    end select

    ! store data from input array

    if (lnk/=0) D(0:lnk-1) = F(:,1)
      
    ! inverse transform

    call fftw_f77_mpi(plani,1,D,W,use_work)
    
    ! unpack data in output array, multiplying by appropriate factor

    if (lnx==0) then
       S = zero ! need to return something even if not used
    else
       S(:,1) = real(D(0:lnx-1),pr)*faci
    end if

  end subroutine inverse_fft


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
       l = mxg+1
       u = pxg+1
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
