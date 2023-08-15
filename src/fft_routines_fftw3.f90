! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module fft_routines
  ! FFT_ROUTINES contains variables and routines related to taking
  ! Fourier transforms of fields using the FFTW 3.x library
  ! 
  ! Modified: 19 July 2010

  use constants, only : pr,pin,int8

  implicit none

  save

  ! FFTW constants

  INTEGER FFTW_R2HC
  PARAMETER (FFTW_R2HC=0)
  INTEGER FFTW_HC2R
  PARAMETER (FFTW_HC2R=1)
  INTEGER FFTW_DHT
  PARAMETER (FFTW_DHT=2)
  INTEGER FFTW_REDFT00
  PARAMETER (FFTW_REDFT00=3)
  INTEGER FFTW_REDFT01
  PARAMETER (FFTW_REDFT01=4)
  INTEGER FFTW_REDFT10
  PARAMETER (FFTW_REDFT10=5)
  INTEGER FFTW_REDFT11
  PARAMETER (FFTW_REDFT11=6)
  INTEGER FFTW_RODFT00
  PARAMETER (FFTW_RODFT00=7)
  INTEGER FFTW_RODFT01
  PARAMETER (FFTW_RODFT01=8)
  INTEGER FFTW_RODFT10
  PARAMETER (FFTW_RODFT10=9)
  INTEGER FFTW_RODFT11
  PARAMETER (FFTW_RODFT11=10)
  INTEGER FFTW_FORWARD
  PARAMETER (FFTW_FORWARD=-1)
  INTEGER FFTW_BACKWARD
  PARAMETER (FFTW_BACKWARD=+1)
  INTEGER FFTW_MEASURE
  PARAMETER (FFTW_MEASURE=0)
  INTEGER FFTW_DESTROY_INPUT
  PARAMETER (FFTW_DESTROY_INPUT=1)
  INTEGER FFTW_UNALIGNED
  PARAMETER (FFTW_UNALIGNED=2)
  INTEGER FFTW_CONSERVE_MEMORY
  PARAMETER (FFTW_CONSERVE_MEMORY=4)
  INTEGER FFTW_EXHAUSTIVE
  PARAMETER (FFTW_EXHAUSTIVE=8)
  INTEGER FFTW_PRESERVE_INPUT
  PARAMETER (FFTW_PRESERVE_INPUT=16)
  INTEGER FFTW_PATIENT
  PARAMETER (FFTW_PATIENT=32)
  INTEGER FFTW_ESTIMATE
  PARAMETER (FFTW_ESTIMATE=64)
  INTEGER FFTW_ESTIMATE_PATIENT
  PARAMETER (FFTW_ESTIMATE_PATIENT=128)
  INTEGER FFTW_BELIEVE_PCOST
  PARAMETER (FFTW_BELIEVE_PCOST=256)
  INTEGER FFTW_NO_DFT_R2HC
  PARAMETER (FFTW_NO_DFT_R2HC=512)
  INTEGER FFTW_NO_NONTHREADED
  PARAMETER (FFTW_NO_NONTHREADED=1024)
  INTEGER FFTW_NO_BUFFERING
  PARAMETER (FFTW_NO_BUFFERING=2048)
  INTEGER FFTW_NO_INDIRECT_OP
  PARAMETER (FFTW_NO_INDIRECT_OP=4096)
  INTEGER FFTW_ALLOW_LARGE_GENERIC
  PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
  INTEGER FFTW_NO_RANK_SPLITS
  PARAMETER (FFTW_NO_RANK_SPLITS=16384)
  INTEGER FFTW_NO_VRANK_SPLITS
  PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
  INTEGER FFTW_NO_VRECURSE
  PARAMETER (FFTW_NO_VRECURSE=65536)
  INTEGER FFTW_NO_SIMD
  PARAMETER (FFTW_NO_SIMD=131072)
  INTEGER FFTW_NO_SLOW
  PARAMETER (FFTW_NO_SLOW=262144)
  INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
  PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
  INTEGER FFTW_ALLOW_PRUNING
  PARAMETER (FFTW_ALLOW_PRUNING=1048576)
  
  ! Global Parameters:
  ! PLAN_LEVEL = how much overhead spent on optimizing
  ! FFT routines before problem is run
  ! PF_XY = used by FFTW for forward transforms
  ! PI_XY = used by FFTW for inverse transforms
  ! PF_SY_X = used by FFTW for forward transform along x, symmetry in y direction
  ! PI_SY_X = used by FFTW for inverse transform along x, symmetry in y direction
  ! PF_SY_YE = used by FFTW for forward transform along y, even symmetry in y direction
  ! PF_SY_YO = used by FFTW for forward transform along y, odd symmetry in y direction
  ! PI_SY_YE = used by FFTW for inverse transform along y, even symmetry in y direction
  ! PI_SY_YO = used by FFTW for inverse transform along y, odd symmetry in y direction
  ! INV_FFT_FAC = normalization factor for inverse FFT
  ! EXP_FWRD = exponential for forward transform (single Fourier mode)
  ! EXP_BWRD = exponential for backward transform (single Fourier mode)

  integer(int8) :: plan_level
  integer(int8) :: pf_xy,pi_xy,pf_sy_x,pf_sy_ye,pf_sy_yo,pi_sy_x,pi_sy_ye,pi_sy_yo
  real(pr) :: inv_fft_fac
  complex(pr) :: exp_fwrd,exp_bwrd

  ! COUNTS_X = MPI message sizes (x distributed)
  ! DISPLS_X = MPI message displacements (x distributed)
  ! COUNTS_K = MPI message sizes (k distributed)
  ! DISPLS_K = MPI message displacements (k distributed)

  integer(pin),dimension(:),allocatable :: counts_x,displs_x,counts_k,displs_k

  ! S = temporary array for spatial-domain field
  ! F = temporary array for Fourier-domain field
  ! SX = temporary array for 1D spatial-domain field in x direction
  ! SY = temporary array for 1D spatial-domain field in y direction
  ! FX = temporary array for 1D Fourier-domain field in x direction
  ! FY = temporary array for 1D Fourier-domain field in y direction

  real(pr),dimension(:,:),allocatable :: S
  complex(pr),dimension(:,:),allocatable :: F
  real(pr),dimension(:),allocatable :: Sx,Sy,Fy
  complex(pr),dimension(:),allocatable :: Fx

contains


  subroutine init_fft(mdl,cnv)
    ! INIT_FFT initializes fft_routines variables
    ! 
    ! Modified: 2 August 2007

    use constants, only : one,two,img
    use model, only : model_type
    use convolution, only : convolution_type
    use mpi_routines, only : nprocs,my_rank
    use io, only : message

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables

    type(model_type),intent(inout) :: mdl
    type(convolution_type),intent(inout) :: cnv

    character(64) :: str1,str2

    ! allocate memory for temporary arrays

    select case(mdl%dim)
    case(1,2)
       allocate(S(mdl%nx,mdl%ny))
       allocate(F(mdl%nx/2+1,mdl%ny))
    case(3)
       allocate(S(mdl%nx,mdl%ny))
       allocate(F(mdl%nx/2+1,mdl%ny))
       if (mdl%sym_y) then
          allocate(Sx(mdl%nx))
          allocate(Fx(mdl%nx/2+1))
          allocate(Sy(mdl%ny))
          allocate(Fy(mdl%ny))
       end if
    end select

    ! initialize FFT variables

    if (mdl%dim==1) then

       exp_fwrd = exp(-img*mdl%k*mdl%x(1))
       exp_bwrd = exp( img*mdl%k*mdl%x(1))

    else

       ! assign the factor needed to get the proper inverse transform
       
       if (mdl%sym_y) then 
          inv_fft_fac = one/(real(mdl%nx,pr)*two*real(mdl%ny,pr))
       else
          inv_fft_fac = one/(real(mdl%nx,pr)*real(mdl%ny,pr))
       end if
       
       ! plan fourier transforms (required by FFTW)
       
       call plan_fft(mdl,mdl%sym_y)

    end if

    ! set array bounds (for particular process)

    allocate(counts_x(0:nprocs-1),displs_x(0:nprocs-1))
    call decompose1d_mpi(mdl%nx,nprocs,my_rank,mdl%mx,mdl%px,counts_x,displs_x)

    cnv%nkx = mdl%nx/2+1
    allocate(counts_k(0:nprocs-1),displs_k(0:nprocs-1))
    call decompose1d_mpi(cnv%nkx,nprocs,my_rank,cnv%mkx,cnv%pkx,counts_k,displs_k)

    mdl%holds_x = .true.
    cnv%holds_k = .true.

    !write(str1,'(a,i6,a,i6,a,i6)') &
    !     'x :',mdl%mx ,'<=i<=',mdl%px ,' cnt=',counts_x(my_rank)
    !write(str2,'(a,i6,a,i6,a,i6)') &
    !     'kx:',cnv%mkx,'<=i<=',cnv%pkx,' cnt=',counts_k(my_rank)
    !write(str1,'(a,i6,a,i6,a,i6,a,i6)') &
    !     'x :',mdl%mx ,'<=i<=',mdl%px ,' cnt=',counts_x(my_rank),' displ=',displs_x(my_rank)
    !write(str2,'(a,i6,a,i6,a,i6,a,i6)') &
    !     'kx:',cnv%mkx,'<=i<=',cnv%pkx,' cnt=',counts_k(my_rank),' displ=',displs_k(my_rank)
    !call message(str1)
    !call message(str2)

  end subroutine init_fft


  subroutine destroy_fft
    ! DESTROY_FFT deallocates temporary arrays
    ! 
    ! Modified: 21 July 2007

    implicit none

    ! deallocate temporary arrays

    if (allocated(S)) deallocate(S)
    if (allocated(F)) deallocate(F) 

    if (allocated(counts_x)) deallocate(counts_x)  
    if (allocated(displs_x)) deallocate(displs_x)  
    if (allocated(counts_k)) deallocate(counts_k)  
    if (allocated(displs_k)) deallocate(displs_k)  

    if (allocated(Sx)) deallocate(Sx)
    if (allocated(Sy)) deallocate(Sy)
    if (allocated(Fx)) deallocate(Fx)  
    if (allocated(Fy)) deallocate(Fy)  

  end subroutine destroy_fft

  
  subroutine plan_fft(mdl,sym_y)
    ! PLAN_FFT runs several test problems to choose the optimum FFT method
    ! 
    ! Modified: 12 September 2007

    use constants, only : real4,real8
    use model, only : model_type
    use io, only : error

    implicit none    

    ! I/O Parameters:
    ! MDL = model variables
    ! SYM_Y = symmetry in y direction

    type(model_type),intent(in) :: mdl
    logical,intent(in) :: sym_y

    ! FFTW requires that you 'plan' the FFTs before executing them.
    ! It often runs a few tests to pick the fastest routine given
    ! the array dimensions.  There is obviously a trade-off between
    ! the time spent planning and the time saved later, which is controlled
    ! by the following variables.  From least to most planning are:
    ! fftw_estimate, fftw_measure, fftw_patient, fftw_exhaustive

    ! assign FFTW plan level

    plan_level = fftw_measure

    ! single precision options - this requires installing the single-precision
    ! version of FFTW in addition to the standard double-precision version.
    ! This section is included if you wish to do this.  It should be as simple
    ! as copying what was done below for double precision, but replacing
    ! any instances of dfftw_ with dfftwf_.  This also requires linking to
    ! the single precision library in the Makefile.

    select case(pr)

    case(real4) ! single-precision transforms

       call error('Single precision FFTs not implemented','plan_fft')

    case(real8) ! double-precision transforms

       select case(mdl%dim)
       case(2)
          call dfftw_plan_dft_r2c_1d(pf_xy,mdl%nx,S(:,1),F(:,1),plan_level)
          call dfftw_plan_dft_c2r_1d(pi_xy,mdl%nx,F(:,1),S(:,1),plan_level)
       case(3)
          if (sym_y) then
             call dfftw_plan_r2r_1d(pf_sy_ye,mdl%ny,Sy,Fy,FFTW_REDFT10,plan_level)
             call dfftw_plan_r2r_1d(pi_sy_ye,mdl%ny,Fy,Sy,FFTW_REDFT01,plan_level)
             call dfftw_plan_r2r_1d(pf_sy_yo,mdl%ny,Sy,Fy,FFTW_RODFT10,plan_level)
             call dfftw_plan_r2r_1d(pi_sy_yo,mdl%ny,Fy,Sy,FFTW_RODFT01,plan_level)
             call dfftw_plan_dft_r2c_1d(pf_sy_x,mdl%nx,Sx,Fx,plan_level)
             call dfftw_plan_dft_c2r_1d(pi_sy_x,mdl%nx,Fx,Sx,plan_level)
          else
             call dfftw_plan_dft_r2c_2d(pf_xy,mdl%nx,mdl%ny,S,F,plan_level)
             call dfftw_plan_dft_c2r_2d(pi_xy,mdl%nx,mdl%ny,F,S,plan_level) 
          end if
       end select

    end select

  end subroutine plan_fft


  subroutine forward_fft(mdl,cnv,R,C,sym) 
    ! FORWARD_FFT forward FFTs a spatial field to Fourier domain
    ! 
    ! Modified: 12 September 2007

    use model, only : model_type
    use convolution, only : convolution_type
    use constants, only : zero,img,onec
    use mpi_routines, only : is_master,master,MPI_REAL_PR,MPI_CPLX_PR,my_rank
    use mpi

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! R = field in space
    ! C = field in Fourier domain
    ! SYM = even or odd symmetry

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(in) :: cnv
    real(pr),dimension(mdl%mx:mdl%px,mdl%ny),intent(in) :: R
    complex(pr),dimension(cnv%mkx:cnv%pkx,cnv%nky),intent(out) :: C
    character(1),intent(in),optional :: sym

    ! Internal Parameters:
    ! PY = plan for transform in y direction
    ! I = index in x and kx direction
    ! J = index in y and ky direction
    ! IERR = MPI error flag
    ! T = temporary array
    ! FAC = factor used for DCT and DST
    ! RM = array in physical domain
    ! CM = array in Fourier domain

    integer(int8) :: py
    integer(pin) :: i,j,ierr
    complex(pr) :: fac
    real(pr),dimension(:,:),allocatable :: Rm,T
    complex(pr),dimension(:,:),allocatable :: Cm

    ! special treatment of low dimensional cases

    select case(mdl%dim)
    case(0)
       C = cmplx(R,zero,pr)
       return
    case(1)
       C = exp_fwrd*R
       return
    end select

    ! allocate memory

    allocate(Rm(mdl%nx,mdl%ny),T(mdl%nx,mdl%ny))
    allocate(Cm(cnv%nkx,cnv%nky))

    ! gather R from slaves to Rm on master

    do j = 1,mdl%ny
       call MPI_Gatherv(R(:,j),counts_x(my_rank),MPI_REAL_PR,Rm(:,j),counts_x,displs_x,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
    end do

    ! transform (on master)

    if (is_master) then

       if (mdl%sym_y) then
          select case(sym)
          case('e')
             py = pf_sy_ye
             fac = onec ! since exp(-ikx)=cos(kx)-i*sin(kx) (and taking coefficient of cosine)
          case('o')
             py = pf_sy_yo
             fac = -img ! since exp(-ikx)=cos(kx)-i*sin(kx) (and taking coefficient of sine)
          end select

          do i = 1,mdl%nx
             Sy = Rm(i,:)
             call dfftw_execute(py) ! Sy -> Fy
             T(i,:) = Fy
          end do

          do j = 1,mdl%ny
             Sx = T(:,j)
             call dfftw_execute(pf_sy_x) ! Sx -> Fx
             Cm(:,j) = Fx*fac
          end do

       else

          S = Rm
          call dfftw_execute(pf_xy)
          Cm = F

       end if
       
    end if

    ! scatter Cm from master to C on slaves

    do j = 1,cnv%nky
       call MPI_Scatterv(Cm(:,j),counts_k,displs_k,MPI_CPLX_PR,C(:,j),counts_k(my_rank),MPI_CPLX_PR,master,MPI_COMM_WORLD,ierr)
    end do

    ! deallocate temporary arrays

    deallocate(Rm,T,Cm)

  end subroutine forward_fft


  subroutine inverse_fft(mdl,cnv,C,R,sym)
    ! INVERSE_FFT inverse FFTs a Fourier-domain field to space
    ! 
    ! Modified: 12 September 2007

    use model, only : model_type
    use convolution, only : convolution_type
    use constants, only : img,onec,zero
    use mpi_routines, only : is_master,master,MPI_REAL_PR,MPI_CPLX_PR,my_rank
    use mpi

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! C = field in Fourier domain
    ! R = field in space
    ! SYM = even or odd symmetry

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(in) :: cnv
    complex(pr),dimension(cnv%mkx:cnv%pkx,cnv%nky),intent(in) :: C
    real(pr),dimension(mdl%mx:mdl%px,mdl%ny),intent(inout) :: R
    character(1),intent(in),optional :: sym

    ! Internal Parameters:
    ! PY = plan for transform in y direction
    ! I = index in x and kx direction
    ! J = index in y and ky direction
    ! IERR = MPI error flag
    ! T = temporary array
    ! FAC = factor used for DCT and DST
    ! RM = array in physical domain
    ! CM = array in Fourier domain

    integer(int8) :: py
    integer(pin) :: i,j,ierr
    complex(pr) :: fac
    real(pr),dimension(:,:),allocatable :: Rm,T
    complex(pr),dimension(:,:),allocatable :: Cm

    ! special treatment of low dimensional cases

    select case(mdl%dim)
    case(0)
       R = real(C,pr)
       return
    case(1)
       R = real(exp_bwrd*C,pr)
       return
    end select

    ! allocate memory

    allocate(Rm(mdl%nx,mdl%ny),T(mdl%nx,mdl%ny))
    allocate(Cm(cnv%nkx,cnv%nky))

    ! gather C from slaves to Cm on master

    do j = 1,cnv%nky
       call MPI_Gatherv(C(:,j),counts_k(my_rank),MPI_CPLX_PR,Cm(:,j),counts_k,displs_k,MPI_CPLX_PR,master,MPI_COMM_WORLD,ierr)
    end do

    ! transform (on master)
    
    if (is_master) then

       if (mdl%sym_y) then
          
          select case(sym)
          case('e')
             py = pi_sy_ye
             fac = onec
          case('o')
             py = pi_sy_yo
             fac = img
          end select

          do j = 1,mdl%ny
             Fx = Cm(:,j)*fac
             call dfftw_execute(pi_sy_x) ! Fx -> Sx
             T(:,j) = Sx
          end do

          do i = 1,mdl%nx
             Fy = T(i,:)
             call dfftw_execute(py) ! Fy -> Sy
             S(i,:) = Sy*inv_fft_fac
          end do
          
       else

          ! DEBUG: CHECK IF NYQUIST AND ZERO MODES HAVE ZERO IMAG PART
          
          if (aimag(Cm(1,1))/=zero.or.aimag(Cm(cnv%nkx,1))/=zero) then
             print *, 'zero or Nyquist mode has nonzero imaginary part'
             print *, 'zero:',Cm(1,1)
             print *, 'Nyquist:',Cm(cnv%nkx,1)
             stop
          end if

          F = Cm
          call dfftw_execute(pi_xy)
          S = S*inv_fft_fac
          
       end if
    
       ! S now holds new value inverse transformed field and R 
       ! holds old value, so overwrite old value with new one
       
       Rm = S

    end if

    ! scatter Rm from master to R on slaves

    do j = 1,mdl%ny
       call MPI_Scatterv(Rm(:,j),counts_x,displs_x,MPI_REAL_PR,R(:,j),counts_x(my_rank),MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
    end do

    ! deallocate temporary arrays

    deallocate(Rm,T,Cm)

  end subroutine inverse_fft


  subroutine decompose1d_mpi(n,p,i,l,u,counts,displs)

    use mpi_routines, only : decompose1d
    use mpi

    implicit none
    
    ! N = number of tasks
    ! P = number of processors
    ! I = rank of process (0<=i<=p-1)
    ! L = starting index
    ! U = ending index

    integer(pin),intent(in) :: n,p,i
    integer(pin),intent(out) :: l,u
    integer(pin),dimension(0:p-1),intent(out) :: counts,displs

    integer :: c,ierr

    call decompose1d(n,p,i,l,u,c)

    call MPI_Allgather(c,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_Allgather(l,1,MPI_INTEGER,displs,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    displs = displs-1

  end subroutine decompose1d_mpi


  subroutine decompose1d_old(n,p,i,l,u,counts,displs)

    use mpi

    implicit none
    
    ! N = number of tasks
    ! P = number of processors
    ! I = rank of process (0<=i<=p-1)
    ! L = starting index
    ! U = ending index

    integer(pin),intent(in) :: n,p,i
    integer(pin),intent(out) :: l,u
    integer(pin),dimension(0:p-1),intent(out) :: counts,displs

    integer(pin) :: c,j,datatype,ierr

    ! have MPI perform decomposition (for consistency with data output routines)

    call MPI_Type_create_darray(p,i,1,(/n/),(/MPI_DISTRIBUTE_BLOCK/), &
         (/MPI_DISTRIBUTE_DFLT_DARG/),(/p/),MPI_ORDER_FORTRAN,MPI_BYTE,datatype,ierr)

    ! determine size of local datatype, which is number of tasks (count, c) for process i

    call MPI_Type_size(datatype,c,ierr)

    ! get counts from all processes

    call MPI_Allgather(c,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

    ! cumulative sum of counts gives displacements and lower/upper indices

    displs(0) = 0
    l = 1
    do j = 0,p-2
       displs(j+1) = displs(j)+counts(j)
       if (i>j) l = l+counts(j)
    end do

    u = l+c-1

  end subroutine decompose1d_old


  function proc_holds_x()

    implicit none

    logical :: proc_holds_x

    proc_holds_x = .true.

  end function proc_holds_x


end module fft_routines
