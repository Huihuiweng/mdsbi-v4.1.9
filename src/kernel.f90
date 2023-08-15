! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module kernel

  ! KERNEL contains variables and routines pertaining to the convolution kernel
  ! 
  ! Modified: 20 July 2010

  use constants, only : pr,pin
  implicit none

  ! KERNEL_TYPE is a derived type containing kernel variables
  !
  ! Parameters:
  ! LOADED = flag indicating if kernel table is loaded (default=.false.)
  ! METHOD = method for calculating kernel at arbitrary T
  !      linear = linear interpolation (fast)
  !      spline = cubic spline interpolation (slow)
  ! PRECALC = precalculate kernel or not
  ! NT = number of points in kernel file 
  ! T = values of T=k*cs*t at which kernel has been evaluated in file
  ! C1  = mode I   kernel
  ! C2  = mode II  kernel
  ! C3  = mode III kernel
  ! C11 = in-plane horizontal-horizontal kernel
  ! C12 = in-plane horizontal-vertical kernel
  ! C22 = in-plane vertical-vertical kernel
  ! C33 = anti-plane kernel
  ! C1S  = mode I   kernel (static)
  ! C2S  = mode II  kernel (static)
  ! C3S  = mode III kernel (static)
  ! C11S = in-plane horizontal-horizontal kernel (static)
  ! C12S = in-plane horizontal-vertical kernel (static)
  ! C22S = in-plane vertical-vertical kernel (static)
  ! C33S = anti-plane kernel (static)
  ! S1  = mode I   kernel (spline coefficient*)
  ! S2  = mode II  kernel (spline coefficient*)
  ! S3  = mode III kernel (spline coefficient*)
  ! S11 = in-plane horizontal-horizontal kernel (spline coefficient*)
  ! S12 = in-plane horizontal-vertical kernel (spline coefficient*)
  ! S22 = in-plane vertical-vertical kernel (spline coefficient*)
  ! S33 = anti-plane kernel (spline coefficient*)
  ! *one-half of second derivative of convolution kernel with respect to T

  type kernel_type  
     logical :: loaded,precalc
     integer(pin) :: nT
     character(64) :: method
     real(pr) :: C1S,C2S,C3S,C11S,C12S,C22S,C33S
     real(pr),dimension(:),allocatable ::  T, &
          C1,C2,C3,C11,C12,C22,C33, &
          S1,S2,S3,S11,S12,S22,S33
  end type kernel_type


contains


  subroutine read_kernel(ninput,krn)
    ! READ_KERNEL reads in kernel variables from file
    ! 
    ! Modified: 20 July 2010

    use io, only : error

    implicit none

    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! KRN = kernel variables

    integer(pin),intent(in) :: ninput
    type(kernel_type),intent(out) :: krn

    ! Internal Parameters:
    ! PRECALC = precalculate kernel or not
    ! METHOD  = method for evaluating kernel at arbitrary T
    ! STAT = I/O error flag

    logical :: precalc
    character(64) :: method
    integer(pin) :: stat

    ! make namelist of user input variables

    namelist /kernel_list/ precalc,method

    ! defaults
    
    precalc = .true.
    method = 'linear'
    
    ! read namelist from input file, call error routine if needed
    
    rewind(ninput)
    read(ninput,nml=kernel_list,iostat=stat)
    if (stat>0) call error("Error reading namelist 'kernel_list' in .in file",'read_kernel')
    
    ! check if input variable values are acceptable
    
    if (.not.(method=='linear'.or.method=='spline')) &
         call error("Invalid parameter 'method' in namelist 'kernel_list' in .in file",'read_kernel')
    
    ! assign input variables to components of derived type
    
    krn%precalc = precalc
    krn%method = method

  end subroutine read_kernel


  subroutine init_kernel(necho,mdl,krn,formulation)
    ! INIT_KERNEL initializes kernel variables
    ! 
    ! Modified: 30 May 2007

    use constants, only : one,two
    use io, only : write_matlab,error
    use model, only : model_type
    use mpi_routines, only : is_master

    implicit none

    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! MDL = model variables
    ! KRN = kernel variables
    ! FORMULATION = convolution scheme

    integer(pin),intent(in) :: necho
    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(inout) :: krn
    character(*),intent(in) :: formulation
       
    ! static convolution kernels
    
    if (mdl%bm) then
       ! check that etap=etam (not set up to handle etap/=etam case, which would
       ! require separate storage of convolution kernels on each side)
       if (mdl%etap/=mdl%etam) &
            call error("Error: Poisson's ratio must be identical on both sides",'init_kernel')
       krn%C11S = two/(one+one/mdl%etap**2) ! H11
       krn%C12S = mdl%etap*(mdl%etap-one)**2/(mdl%etap**2+one) ! H12
       krn%C22S = two/(one+one/mdl%etap**2) ! H22
       krn%C33S = one ! H33
    else
       krn%C1S = two*(one-one/mdl%eta**2) !C1 (not yet for combined bimaterial)
       krn%C2S = two*(one-one/mdl%eta**2) !C2 (not yet for combined bimaterial)
       krn%C3S = two/(one+mdl%rrho*mdl%rcs**2) ! C3 (combined bimaterial)
    end if
    
    ! load convolution kernels using one of two options (comment one)

    ! option 1: master reads table, distributes kernels

    !krn%loaded = .false.
    !if (is_master) call load_kernel(mdl,krn,formulation)
    !call distribute_kernel(mdl,krn,formulation)

    ! option 2: each slave reads table, no need to distribute kernels

    call load_kernel(mdl,krn,formulation)

    ! output variable values into matlab file

    if (is_master) then

       call write_matlab(necho,'precalc',krn%precalc,'krn')
       call write_matlab(necho,'method',krn%method,'krn')

    end if

  end subroutine init_kernel


  subroutine load_kernel(mdl,krn,formulation)
    ! LOAD_KERNEL loads convolution kernels
    ! 
    ! Modified: 28 May 2007

    use constants, only : zero
    use utilities, only : init_spline,cumquad_trap,cumquad_spline
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! KRN = kernel variables
    ! FORMULATION = convolution scheme

    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(inout) :: krn
    character(*),intent(in) :: formulation

    ! return if no kernels are needed

    if (formulation=='static'.or.formulation=='quasidynamic') return

    ! read kernel table from file

    if (mdl%bm) then
       call get_kernel(mdl,krn,'kernel/kernel_db.dat')
    else
       call get_kernel(mdl,krn,'kernel/kernel_di.dat')
    end if

    ! initialize the interpolation method of choice
    ! (NOTE: H12 has non-zero slope at T=0; ideally this would be calculated and used
    ! as a boundary condition on the spline fit, but instead use a natural spline;
    ! all other kernels have zero slope at T=0, so use that as boundary condition)

    if (krn%method=='spline') then

       if (mdl%bm) then

          allocate(krn%S11(krn%nT))
          allocate(krn%S12(krn%nT))
          allocate(krn%S22(krn%nT))
          allocate(krn%S33(krn%nT))
          
          call init_spline(krn%nT,krn%T,krn%C11,krn%S11,zero)
          call init_spline(krn%nT,krn%T,krn%C12,krn%S12) ! natural spline
          call init_spline(krn%nT,krn%T,krn%C22,krn%S22,zero)
          call init_spline(krn%nT,krn%T,krn%C33,krn%S33,zero)

       else

          allocate(krn%S1(krn%nT))
          allocate(krn%S2(krn%nT))
          allocate(krn%S3(krn%nT))

          call init_spline(krn%nT,krn%T,krn%C1,krn%S1,zero)
          call init_spline(krn%nT,krn%T,krn%C2,krn%S2,zero)
          call init_spline(krn%nT,krn%T,krn%C3,krn%S3,zero)

       end if

    end if

    ! calculate cumulative integral of displacement formulation kernel to obtain 
    ! velocity formulation kernel (if needed)

    if (formulation=='velocity') then

       select case(krn%method)

       case('spline')

          if (mdl%bm) then
             krn%C11 = krn%C11S-cumquad_spline(krn%T,krn%C11,krn%S11)
             krn%C12 = krn%C12S-cumquad_spline(krn%T,krn%C12,krn%S12)
             krn%C22 = krn%C22S-cumquad_spline(krn%T,krn%C22,krn%S22)
             krn%C33 = krn%C33S-cumquad_spline(krn%T,krn%C33,krn%S33)
          else
             krn%C1 = krn%C1S-cumquad_spline(krn%T,krn%C1,krn%S1)
             krn%C2 = krn%C2S-cumquad_spline(krn%T,krn%C2,krn%S2)
             krn%C3 = krn%C3S-cumquad_spline(krn%T,krn%C3,krn%S3)
          end if

       case default

          if (mdl%bm) then
             krn%C11 = krn%C11S-cumquad_trap(krn%T,krn%C11)
             krn%C12 = krn%C12S-cumquad_trap(krn%T,krn%C12)
             krn%C22 = krn%C22S-cumquad_trap(krn%T,krn%C22)
             krn%C33 = krn%C33S-cumquad_trap(krn%T,krn%C33)
          else
             krn%C1 = krn%C1S-cumquad_trap(krn%T,krn%C1)
             krn%C2 = krn%C2S-cumquad_trap(krn%T,krn%C2)
             krn%C3 = krn%C3S-cumquad_trap(krn%T,krn%C3)
          end if

       end select

    end if

    ! set flag to indicate that kernel is loaded

    krn%loaded = .true.

  end subroutine load_kernel


  subroutine get_kernel(mdl,krn,filename)
    ! GET_KERNEL loads table of kernel values from file
    ! 
    ! Modified: 30 May 2007

    use io, only : new_io_unit,error
    use model, only : model_type

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! KRN = kernel variables
    ! FILENAME = name of kernel file

    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(inout) :: krn
    character(*),intent(in) :: filename

    ! Internal Parameters:
    ! HEADER = file header (description of kernel file)
    ! IT =  index of T vector (and kernels)
    ! NT = number of points in kernel file
    ! STAT = I/O error flag

    character(512) :: header
    integer(pin) :: iT,nT,nkernel,stat

    ! open unit to kernel file, call error routine if needed

    nkernel = new_io_unit()
    open(nkernel,file=filename,status='old',iostat=stat)
    if (stat/=0) call error("Error opening file: '" // trim(adjustl(filename)) // "'",'get_kernel')

    ! read header of kernel file

    read(nkernel,'(a)',iostat=stat) header
    if (stat/=0) call error("Error reading kernel file: '" // trim(adjustl(filename)) // "'",'get_kernel')
    read(nkernel,*,iostat=stat) nT
    if (stat/=0) call error("Error reading kernel file: '" // trim(adjustl(filename)) // "'",'get_kernel')

    ! assign length of T to derived type

    krn%nT = nT

    ! allocate memory to arrays that store kernel data

    allocate(krn%T(nT))

    if (mdl%bm) then
       allocate(krn%C11(krn%nT))
       allocate(krn%C12(krn%nT))
       allocate(krn%C22(krn%nT))
       allocate(krn%C33(krn%nT))
    else
       allocate(krn%C1(krn%nT))
       allocate(krn%C2(krn%nT))
       allocate(krn%C3(krn%nT))
    end if

    ! read kernel data from kernel file

    if (mdl%bm) then
       do iT = 1,nT
          read(nkernel,*,iostat=stat) krn%T(iT),krn%C11(iT),krn%C12(iT),krn%C22(iT),krn%C33(iT)
          if (stat/=0) call error("Error reading kernel file: '" // trim(adjustl(filename)) // "'",'get_kernel')
       end do
    else
       do iT = 1,nT
          read(nkernel,*,iostat=stat) krn%T(iT),krn%C1(iT),krn%C2(iT),krn%C3(iT)
          if (stat/=0) call error("Error reading kernel file: '" // trim(adjustl(filename)) // "'",'get_kernel')
       end do
    end if

    ! close kernel file

    close(nkernel,iostat=stat)
    if (stat/=0) call error("Error closing file: '" // trim(adjustl(filename)) // "'",'get_kernel')

  end subroutine get_kernel


  subroutine distribute_kernel(mdl,krn,formulation)
    ! DISTRIBUTE_KERNEL distributes convolution kernel from master to slaves
    ! 
    ! Modified: 30 May 2007

    use model, only : model_type
    use mpi_routines, only : master,MPI_REAL_PR
    use mpi

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! KRN = kernel variables
    ! FORMULATION = convolution scheme

    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(inout) :: krn
    character(*),intent(in) :: formulation

    ! Internal Parameters:
    ! IERR = MPI error flag

    integer(pin) :: ierr

    ! return if no kernels are needed

    if (formulation=='static'.or.formulation=='quasidynamic') return

    ! distribute basic data

    call MPI_Bcast(krn%nT,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

    ! allocate arrays using shared data

    if (.not.krn%loaded) then

       allocate(krn%T(krn%nT))

       if (mdl%bm) then
          allocate(krn%C11(krn%nT))
          allocate(krn%C12(krn%nT))
          allocate(krn%C22(krn%nT))
          allocate(krn%C33(krn%nT))
       else
          allocate(krn%C1(krn%nT))
          allocate(krn%C2(krn%nT))
          allocate(krn%C3(krn%nT))
       end if
       
       if (krn%method=='spline') then
          if (mdl%bm) then
             allocate(krn%S11(krn%nT))
             allocate(krn%S12(krn%nT))
             allocate(krn%S22(krn%nT))
             allocate(krn%S33(krn%nT))
          else
             allocate(krn%S1(krn%nT))
             allocate(krn%S2(krn%nT))
             allocate(krn%S3(krn%nT))
          end if
       end if

    end if

    ! distribute kernels

    call MPI_Bcast(krn%T,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)

    if (mdl%bm) then
       call MPI_Bcast(krn%C11,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(krn%C12,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(krn%C22,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(krn%C33,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
    else
       call MPI_Bcast(krn%C1,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(krn%C2,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(krn%C3,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
    end if

    if (krn%method=='spline') then
       if (mdl%bm) then
          call MPI_Bcast(krn%S11,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
          call MPI_Bcast(krn%S12,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
          call MPI_Bcast(krn%S22,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
          call MPI_Bcast(krn%S33,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
       else
          call MPI_Bcast(krn%S1,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
          call MPI_Bcast(krn%S2,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
          call MPI_Bcast(krn%S3,krn%nT,MPI_REAL_PR,master,MPI_COMM_WORLD,ierr)
       end if
    end if

    ! set flag to indicate that kernel is loaded

    krn%loaded = .true.

  end subroutine distribute_kernel


  subroutine precalc_kernel(C,n,k,mdl,krn)
    ! PRECALC_KERNEL precalculates convolution kernel for a particular mode
    ! 
    ! Modified: 7 August 2007

    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! C = kernel array
    ! N = number of time steps in history array (at which kernel will be calculated)
    ! MDL = model variables
    ! KRN = kernel variables
    ! K = wavenumber amplitude

    real(pr),dimension(:,:,:),intent(out) :: C
    integer(pin),intent(in) :: n
    real(pr),intent(in) :: k
    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(in) :: krn

    ! Internal Parameters:
    ! IT = counter for time steps
    ! T = non-dimensional time vector
    ! TP = non-dimensional time vector, plus side
    ! TM = non-dimensional time vector, minus side

    integer(pin) :: iT
    real(pr),dimension(0:n) :: T,Tp,Tm

    if (mdl%bm) then

       Tp = k*mdl%csp*mdl%dt*real( (/ (iT,iT=0,n) /) ,pr)
       Tm = k*mdl%csm*mdl%dt*real( (/ (iT,iT=0,n) /) ,pr)
       
       C(:,1,1) = calc_kernel(Tp,11,krn)
       C(:,1,2) = calc_kernel(Tm,11,krn)
       C(:,2,1) = calc_kernel(Tp,12,krn)
       C(:,2,2) = calc_kernel(Tm,12,krn)
       C(:,3,1) = calc_kernel(Tp,22,krn)
       C(:,3,2) = calc_kernel(Tm,22,krn)
       C(:,4,1) = calc_kernel(Tp,33,krn)
       C(:,4,2) = calc_kernel(Tm,33,krn)
       
    else

       T = k*mdl%cs*mdl%dt*real( (/ (iT,iT=0,n) /) ,pr)
       
       select case(mdl%dim)
       case default
          select case(mdl%mode)
          case(1)
             C(:,1,1) = calc_kernel(T,1,krn)
          case(2)
             C(:,1,1) = calc_kernel(T,2,krn)
          case(3)
             C(:,1,1) = calc_kernel(T,3,krn)
          case default
             C(:,1,1) = calc_kernel(T,2,krn)
             C(:,2,1) = calc_kernel(T,3,krn)
          end select
       case(3)
          C(:,1,1) = calc_kernel(T,2,krn)
          C(:,2,1) = calc_kernel(T,3,krn)
       end select

    end if

  end subroutine precalc_kernel


  subroutine destroy_kernel(krn)
    ! DESTROY_KERNEL destroys derived type krn
    ! 
    ! Modified: 25 January 2007
    
    implicit none

    ! I/O Parameters:
    ! KRN = kernel variables

    type(kernel_type),intent(inout) :: krn

    ! deallocate memory

    if (allocated(krn%T)) deallocate(krn%T)

    if (allocated(krn%C1)) deallocate(krn%C1)
    if (allocated(krn%C2)) deallocate(krn%C2)
    if (allocated(krn%C3)) deallocate(krn%C3)
    if (allocated(krn%C11)) deallocate(krn%C11)
    if (allocated(krn%C12)) deallocate(krn%C12)
    if (allocated(krn%C22)) deallocate(krn%C22)
    if (allocated(krn%C33)) deallocate(krn%C33)

    if (allocated(krn%S1)) deallocate(krn%S1)
    if (allocated(krn%S2)) deallocate(krn%S2)
    if (allocated(krn%S3)) deallocate(krn%S3)
    if (allocated(krn%S11)) deallocate(krn%S11)
    if (allocated(krn%S12)) deallocate(krn%S12)
    if (allocated(krn%S22)) deallocate(krn%S22)
    if (allocated(krn%S33)) deallocate(krn%S33)

    ! set flag to indicate that kernel is not loaded

    krn%loaded = .false.

  end subroutine destroy_kernel


  function calc_kernel(T,mode,krn) result(C)
    ! CALC_KERNEL returns interpolated kernel C, given vector T
    ! 
    ! Modified: 23 January 2007

    use utilities, only : interp_linear,interp_spline

    implicit none
    
    ! I/O Parameters:
    ! T = values of T at which kernel is to be interpolated
    ! MODE = crack geometry mode
    ! KRN = kernel variables
    ! C = interpolated values of kernel corresponding to input T

    real(pr),dimension(:),intent(in) :: T
    integer(pin),intent(in) :: mode
    type(kernel_type),intent(in) :: krn
    real(pr) :: C(size(T))
    
    ! call interpolation routines

    select case(krn%method)

    case('linear') ! linear interpolation

       select case(mode)
       case(1)
          C = interp_linear(krn%T,krn%C1,T)
       case(2)
          C = interp_linear(krn%T,krn%C2,T)
       case(3)
          C = interp_linear(krn%T,krn%C3,T)
       case(11)
          C = interp_linear(krn%T,krn%C11,T)
       case(12,21)
          C = interp_linear(krn%T,krn%C12,T)
       case(22)
          C = interp_linear(krn%T,krn%C22,T)
       case(33)
          C = interp_linear(krn%T,krn%C33,T)
       end select

    case('spline') ! cubic spline

       select case(mode)
       case(1)
          C = interp_spline(krn%T,krn%C1,krn%S1,T)
       case(2)
          C = interp_spline(krn%T,krn%C2,krn%S2,T)
       case(3)
          C = interp_spline(krn%T,krn%C3,krn%S3,T)
       case(11)
          C = interp_spline(krn%T,krn%C11,krn%S11,T)
       case(12,21)
          C = interp_spline(krn%T,krn%C12,krn%S12,T)
       case(22)
          C = interp_spline(krn%T,krn%C22,krn%S22,T)
       case(33)
          C = interp_spline(krn%T,krn%C33,krn%S33,T)
       end select

    end select
    
  end function calc_kernel


end module kernel
