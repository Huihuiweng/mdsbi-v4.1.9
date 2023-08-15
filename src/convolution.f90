! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module convolution

  ! CONVOLUTION contains variables and routines to calculate stress transfer, including
  ! Fourier coefficients of slip and/or slip velocity history used within convolutions
  ! 
  ! Modified: 20 July 2010
  
  use constants, only : pr,pin
  use history, only : history_type

  implicit none

  ! CONVOLUTION_TYPE is a derived type containing convolution variables
  !
  ! Parameters:
  ! FORMULATION = convolution scheme
  !      static = static
  !      quasidynamic = static+radiation damping
  !      displacement = convolution using slip history
  !      velocity = convolution using slip velocity history
  ! TRUNCATE = convolution truncation
  !      T = truncate convolution at Tmax (see below)
  !      F = do not truncate convolution
  ! TMAX = maximum value of T=k*cs*t beyond which convolution is truncated, where
  !      k = wavenumber of Fourier mode
  !      cs = S-wave speed
  !      t = time
  ! HOLDS_K = flag indicating if this process holds k data
  ! MKX = lower bound on i (for particular process)
  ! PKX = upper bound on i (for particular process)
  ! NKX = length of kx
  ! NKY = length of ky
  ! KX = wavenumber vector in x direction
  ! KY = wavenumber vector in y direction
  ! KN = Nyquist wavenumber
  ! DX = Fourier coefficients of x component of slip ux
  ! DY = Fourier coefficients of y component of slip uy
  ! DXP = Fourier coefficients of x component of displacement ux, plus side
  ! DXM = Fourier coefficients of x component of displacement ux, minus side
  ! DYP = Fourier coefficients of y component of displacement uy, plus side
  ! DYM = Fourier coefficients of y component of displacement uy, minus side
  ! DZP = Fourier coefficients of z component of displacement uz, plus side
  ! DZM = Fourier coefficients of z component of displacement uz, minus side
  ! FX = Fourier coefficients of x component of stress transfer functional
  ! FY = Fourier coefficients of y component of stress transfer functional
  ! FXH = contribution of past time steps (excluding the current time step) to Fx
  ! FYH = contribution of past time steps (excluding the current time step) to Fy
  ! FXP = Fourier coefficients of x component of stress transfer functional, plus side
  ! FXM = Fourier coefficients of x component of stress transfer functional, minus side
  ! FYP = Fourier coefficients of y component of stress transfer functional, plus side
  ! FYM = Fourier coefficients of y component of stress transfer functional, minus side
  ! FZP = Fourier coefficients of z component of stress transfer functional, plus side
  ! FZM = Fourier coefficients of z component of stress transfer functional, minus side
  ! FXHP = contribution of past time steps (excluding the current time step) to Fxp
  ! FXHM = contribution of past time steps (excluding the current time step) to Fxm
  ! FYHP = contribution of past time steps (excluding the current time step) to Fyp
  ! FYHM = contribution of past time steps (excluding the current time step) to Fym
  ! FZHP = contribution of past time steps (excluding the current time step) to Fzp
  ! FZHM = contribution of past time steps (excluding the current time step) to Fzm
  ! HX = Fourier coefficients of x component of history (current step)
  ! HY = Fourier coefficients of y component of history (current step)
  ! HXP = Fourier coefficients of x component of history, plus side (current step)
  ! HXM = Fourier coefficients of x component of history, minus side (current step)
  ! HYP = Fourier coefficients of y component of history, plus side (current step)
  ! HYM = Fourier coefficients of y component of history, minus side (current step)
  ! HZP = Fourier coefficients of z component of history, plus side (current step)
  ! HZM = Fourier coefficients of z component of history, minus side (current step)
  ! H = history of Fourier coefficients for each mode, with indices i,j of H(i,j) labelling the mode

  type convolution_type
     character(64) :: formulation
     logical :: truncate,holds_k
     real(pr) :: kN,Tmax
     integer(pin) :: nkx,nky,mkx,pkx
     real(pr),dimension(:),allocatable :: kx,ky
     complex(pr),dimension(:,:),allocatable :: Dx,Dy, &
          Dxp,Dxm,Dyp,Dym,Dzp,Dzm
     complex(pr),dimension(:,:),allocatable :: Fx,Fy,FxH,FyH, &
          Fxp ,Fxm ,Fyp ,Fym ,Fzp ,Fzm, &
          FxHp,FxHm,FyHp,FyHm,FzHp,FzHm
     complex(pr),dimension(:,:),allocatable :: Hx,Hy, &
          Hxp,Hxm,Hyp,Hym,Hzp,Hzm
     type(history_type),dimension(:,:),pointer :: H=>null()
  end type convolution_type


contains
  

  subroutine read_convolution(ninput,cnv)
    ! READ_CONVOLUTION reads in convolution variables from file
    ! 
    ! Modified: 20 July 2010
    
    use constants, only : zero
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! CNV = convolution variables

    integer(pin),intent(in) :: ninput
    type(convolution_type),intent(out) :: cnv

    ! Internal Parameters:
    ! FORMULATION = convolution scheme
    !      static = static
    !      quasidynamic = static+radiation damping
    !      displacement = convolution using slip history
    !      velocity = convolution using slip velocity history
    ! TRUNCATE = convolution truncation
    !      T = truncate convolution at Tmax
    !      F = do not truncate convolution
    ! TMAX = maximum value of T=k*cs*t beyond which convolution is truncated
    ! STAT = I/O error flag

    character(64) :: formulation
    logical :: truncate
    real(pr) :: Tmax
    integer(pin) :: stat

    ! make namelist of user input variables

    namelist /convolution_list/ formulation,truncate,Tmax

    ! defaults
    
    formulation = 'displacement'
    truncate = .false.
    Tmax = zero
    
    ! read namelist from input file, call error routine if needed
    
    rewind(ninput)
    read(ninput,nml=convolution_list,iostat=stat)
    if (stat/=0) call error("Error reading namelist 'convolution_list' in .in file",'read_convolution')
    
    ! assign input variables to components of derived type
    
    cnv%formulation = formulation
    cnv%truncate = truncate
    cnv%Tmax = Tmax
    
  end subroutine read_convolution

 
  subroutine init_convolution(necho,mdl,krn,cnv)
    ! INIT_CONVOLUTION initializes convolution variables
    ! 
    ! Modified: 21 July 2007
    
    use constants, only : pi
    use model, only : model_type
    use io, only : write_matlab
    use kernel, only : kernel_type
    use mpi_routines, only : is_master

    implicit none

    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! KRN = kernel variables
    ! MDL = model variables
    ! CNV = convolution variables

    integer(pin),intent(in) :: necho
    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(inout) :: krn
    type(convolution_type),intent(inout) :: cnv

    ! set size of wavenumber and Fourier-domain arrays;
    ! only the non-negative portion of kx is stored due to symmetry of transforms of real data
    ! nkx is set in fft_routines (because it relies on how FFT procedure distributes arrays)

    if (mdl%sym_y) then
       cnv%nky = mdl%ny+1
    else
       cnv%nky = mdl%ny
    end if

    ! assign Nyquist wavenumber using
    ! fN=1/(2*h) => kN=2*pi*fN=pi/h

    cnv%kN = pi/mdl%h

    ! initialize wavenumber

    call init_wavenumber(mdl,cnv)

    ! initialize variables specific to identical materials or bimaterial case

    if (mdl%bm) then
       call init_convolution_bm(cnv)
    else
       call init_convolution_im(cnv)
    end if

    ! initialize history vectors (if k data stored for this process)

    if (cnv%holds_k) call init_history(mdl,krn,cnv)

    ! output variable values into matlab file
    
    if (is_master) then

       call write_matlab(necho,'formulation',cnv%formulation,'cnv')
       call write_matlab(necho,'truncate',cnv%truncate,'cnv')
       if (cnv%truncate) call write_matlab(necho,'Tmax',cnv%Tmax,'cnv')

    end if

  end subroutine init_convolution


  subroutine init_wavenumber(mdl,cnv)
    ! INIT_WAVENUMBER initializes wavenumber
    ! 
    ! Modified: 6 July 2008

    use constants, only : zero,two
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(inout) :: cnv

    ! Internal Parameters:
    ! I = index along kx direction
    ! J = index along ky direction
    ! K = temporary wavenumber vector

    integer(pin) :: i,j
    real(pr),dimension(:),allocatable :: k

    ! allocate memory for wavenumber arrays and assign values
    ! In general, kx designates the wavenumber along the spatial x axis,
    ! and ky designates the wavenumber along the spatial y axis.  However,
    ! this is not true for the 1D/2D mixed-mode case.  For this case, the kx-ky
    ! coordinate system is rotated at mdl%angle with respect to the x-y
    ! spatial coordinate system.

    select case(mdl%dim)

    case(1) ! 1D (single Fourier mode)

       ! set bounds on j

       allocate(cnv%kx(1))
       allocate(cnv%ky(1))

       if (mdl%mode==0) then
          ! mixed mode
          cnv%kx = cnv%kN*cos(mdl%angle)
          cnv%ky = cnv%kN*sin(mdl%angle)
       else
          ! pure mode
          cnv%kx = cnv%kN
          cnv%ky = zero
       end if

    case(2)

       ! note that ky is same size as kx for this case

       allocate(     k(cnv%nkx))
       allocate(cnv%kx(cnv%nkx))
       allocate(cnv%ky(cnv%nkx))

       if (cnv%nkx==mdl%nx) then
          ! use positive and negative wavenumbers
          k(1:cnv%nkx/2+1) = two*real( (/ (i, i=0,cnv%nkx/2) /) ,pr)/real(cnv%nkx,pr)
          k(cnv%nkx:cnv%nkx/2+2:-1) = -k(2:cnv%nkx/2)
          !k(cnv%nkx/2+1) = zero ! set Nyquist to zero
       else
          ! only use non-negative half of wavenumbers
          k = real( (/ (i, i=0,cnv%nkx-1) /) ,pr)/real(cnv%nkx-1,pr)
          !k(cnv%nkx) = zero ! set Nyquist to zero
       end if

       k = k*cnv%kN

       if (mdl%mode==0) then
          ! mixed mode
          cnv%kx = k*cos(mdl%angle)
          cnv%ky = k*sin(mdl%angle)
       else
          ! pure mode
          cnv%kx = k
          cnv%ky = zero
       end if
       deallocate(k)

    case(3)

       allocate(cnv%kx(cnv%nkx))
       allocate(cnv%ky(cnv%nky))

       if (cnv%nkx==mdl%nx) then
          ! use positive and negative wavenumbers
          cnv%kx(1:cnv%nkx/2+1) = two*real( (/ (i, i=0,cnv%nkx/2) /) ,pr)/real(cnv%nkx,pr)
          cnv%kx(cnv%nkx:cnv%nkx/2+2:-1) = -cnv%kx(2:cnv%nkx/2)
          !cnv%kx(cnv%nkx/2+1) = zero ! set Nyquist to zero
       else
          ! only use non-negative half of wavenumbers
          cnv%kx = real( (/ (i, i=0,cnv%nkx-1) /) ,pr)/real(cnv%nkx-1,pr)
          !cnv%kx(cnv%nkx) = zero ! set Nyquist to zero
       end if
       cnv%kx = cnv%kx*cnv%kN

       if (mdl%sym_y) then
          cnv%ky = real( (/ (j, j=0,cnv%nky-1) /) ,pr)/real(cnv%nky-1,pr)
          !cnv%ky(cnv%nky) = zero ! set Nyquist to zero
       else
          cnv%ky(1:cnv%nky/2+1) = two*real( (/ (j, j=0,cnv%nky/2) /) ,pr)/real(cnv%nky,pr)
          cnv%ky(cnv%nky:cnv%nky/2+2:-1) = -cnv%ky(2:cnv%nky/2)
          !cnv%ky(cnv%nky/2+1) = zero ! set Nyquist to zero
       end if
       cnv%ky = cnv%ky*cnv%kN

    end select

  end subroutine init_wavenumber


  subroutine init_convolution_im(cnv)
    ! INIT_CONVOLUTION_IM initializes convolution variables, identical materials case
    ! 
    ! Modified: 12 June 2007
    
    use constants, only : zeroc

    implicit none

    ! I/O Parameters:
    ! CNV = convolution variables

    type(convolution_type),intent(inout) :: cnv

    ! allocate memory to arrays of Fourier coefficients of slip

    allocate(cnv%Dx(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Dy(cnv%mkx:cnv%pkx,cnv%nky))    

    cnv%Dx = zeroc
    cnv%Dy = zeroc

    ! allocate memory to arrays of Fourier coefficients of stress transfer 

    allocate(cnv%Fx(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Fy(cnv%mkx:cnv%pkx,cnv%nky))
    
    cnv%Fx = zeroc
    cnv%Fy = zeroc

    ! allocate memory to arrays of Fourier coefficients of current step in history 

    if (.not.(cnv%formulation=='static'.or.cnv%formulation=='quasidynamic')) then

       allocate(cnv%Hx(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%Hx = zeroc
       allocate(cnv%Hy(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%Hy = zeroc

    end if
    
    ! allocate memory to arrays of Fourier coefficients of stress transfer
    ! (past history excluding current time step)

    if (.not.(cnv%formulation=='static'.or.cnv%formulation=='quasidynamic')) then

       allocate(cnv%FxH(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%FxH = zeroc
       allocate(cnv%FyH(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%FyH = zeroc

    end if
    
  end subroutine init_convolution_im


  subroutine init_convolution_bm(cnv)
    ! INIT_CONVOLUTION_BM initializes convolution variables, bimaterial case
    ! 
    ! Modified: 12 June 2007
    
    use constants, only : zeroc

    implicit none

    ! I/O Parameters:
    ! CNV = convolution variables

    type(convolution_type),intent(inout) :: cnv

    ! allocate Fourier transform of displacement arrays

    allocate(cnv%Dxp(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Dxm(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Dyp(cnv%mkx:cnv%pkx,cnv%nky))       
    allocate(cnv%Dym(cnv%mkx:cnv%pkx,cnv%nky)) 
    allocate(cnv%Dzp(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Dzm(cnv%mkx:cnv%pkx,cnv%nky))
    
    cnv%Dxp = zeroc
    cnv%Dxm = zeroc
    cnv%Dyp = zeroc
    cnv%Dym = zeroc
    cnv%Dzp = zeroc
    cnv%Dzm = zeroc

    ! allocate memory to arrays of Fourier coefficents of stress transfer

    allocate(cnv%Fxp(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Fxm(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Fyp(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Fym(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Fzp(cnv%mkx:cnv%pkx,cnv%nky))
    allocate(cnv%Fzm(cnv%mkx:cnv%pkx,cnv%nky))

    cnv%Fxp = zeroc
    cnv%Fxm = zeroc
    cnv%Fyp = zeroc
    cnv%Fym = zeroc
    cnv%Fzp = zeroc
    cnv%Fzm = zeroc

    ! allocate memory to arrays of Fourier coefficients of current step in history 

    if (.not.(cnv%formulation=='quasidynamic'.or.cnv%formulation=='static')) then

       allocate(cnv%Hxp(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%Hxp = zeroc
       allocate(cnv%Hxm(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%Hxm = zeroc
       allocate(cnv%Hyp(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%Hyp = zeroc
       allocate(cnv%Hym(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%Hym = zeroc
       allocate(cnv%Hzp(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%Hzp = zeroc
       allocate(cnv%Hzm(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%Hzm = zeroc

    end if
    
    ! allocate memory to arrays of Fourier coefficients of stress transfer
    ! (past history excluding current time step)

    if (.not.(cnv%formulation=='quasidynamic'.or.cnv%formulation=='static')) then

       allocate(cnv%FxHp(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%FxHp = zeroc
       allocate(cnv%FxHm(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%FxHm = zeroc
       allocate(cnv%FyHp(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%FyHp = zeroc
       allocate(cnv%FyHm(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%FyHm = zeroc
       allocate(cnv%FzHp(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%FzHp = zeroc
       allocate(cnv%FzHm(cnv%mkx:cnv%pkx,cnv%nky))
       cnv%FzHm = zeroc

    end if
    
  end subroutine init_convolution_bm


  subroutine init_history(mdl,krn,cnv)
    ! INIT_HISTORY initializes history arrays
    ! 
    ! Modified: 24 January 2007

    use constants, only : zero
    use model, only : model_type
    use kernel, only : kernel_type,precalc_kernel,destroy_kernel
    use history, only : create_history
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! KRN = kernel variables
    ! CNV = convolution variables

    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(inout) :: krn
    type(convolution_type),intent(inout) :: cnv

    ! Internal Parameters:
    ! I = index along kx direction
    ! J = index along ky direction
    ! NH = number of time steps saved for a given Fourier mode
    ! K = wavenumber
    ! KM = wavenumber product array (not used)
    ! KX = wavenumber in x direction (not used)
    ! KY = wavenumber in y direction (not used)
    ! TMAX = maximum value of T=k*cs*t at which history is saved
    ! KMAX = maximum wavenumber
    ! CSMAX = maximum S-wave speed
    ! STR = character string
    ! STR1 = character string
    ! STR2 = character string

    integer(pin) :: i,j,nH
    real(pr) :: k,km(2,2),kx,ky,Tmax,kmax,csmax
    character(64) :: str,str1,str2

    ! return if histories will not be used

    if (cnv%formulation=='quasidynamic'.or.cnv%formulation=='static') return

    ! check if maximum T=k*cs*t exceeds maximum T in kernel data file
    
    if (mdl%bm) then
       csmax = max(mdl%csp,mdl%csm)
    else
       csmax = mdl%cs
    end if
    
    kmax = sqrt(maxval(cnv%kx)**2+maxval(cnv%ky)**2)
    Tmax = kmax*csmax*real(mdl%nt,pr)*mdl%dt
    if (mdl%tmax/=zero) Tmax = min(Tmax,kmax*csmax*mdl%tmax)
    if (cnv%truncate) Tmax = min(Tmax,cnv%Tmax)
    
    if (Tmax>krn%T(krn%nT)) then
       write(str1,'(f12.2)') krn%T(krn%nT)
       write(str2,'(f12.2)') Tmax
       str = 'Tabulated kernel T=k*cs*t=' // trim(adjustl(str1)) // ', but max(T)=' // trim(adjustl(str2))
       call error(str,'init_history')
    end if
    
    ! allocate memory to arrays of Fourier mode history by looping over all modes and
    ! (1) calculating appropriate wavenumber for the mode
    ! (2) calculate the number of time steps to be saved, which is either the number of time steps
    !     to be run or a value derived from placing a limit on T=k*cs*t (the non-dimensional
    !     argument of the convolution kernel), which is done when truncating the convolution
    ! (3) create a history derived type using these values
    
    allocate(cnv%H(cnv%mkx:cnv%pkx,cnv%nky))
    
    do j = 1,cnv%nky
       do i = cnv%mkx,cnv%pkx
          
          ! calculate wavenumber
          
          call wavenumber(mdl,cnv,i,j,k,km,kx,ky)
          
          if (k==zero) then
             
             nH = 0
             cnv%H(i,j)%n = nH
             
          else
             
             if (cnv%truncate) then
                ! truncate at constant T (nondimensional kernel argument):
                ! since T=k*cs*t and t=n*dt, then n=Tmax/(k*cs*dt)
                nH = min(floor(cnv%Tmax/(k*csmax*mdl%dt)),mdl%nt)
             else
                nH = mdl%nt
             end if
             
             ! limit nH using tmax (if this is set to stop calculation)
             
             if (mdl%tmax/=zero) nH = min(nH,ceiling(mdl%tmax/mdl%dt))
             
             ! initialize history vector
             
             call create_history(mdl,cnv%H(i,j),nH,krn%precalc)
             
             ! precalculate kernel (if needed)
             
             if (krn%precalc) call precalc_kernel(cnv%H(i,j)%C,cnv%H(i,j)%n,k,mdl,krn)
             
          end if
          
       end do
    end do

    ! destroy tabulated kernel vectors in derived type krn if precalculating kernels
    
    if (krn%precalc) call destroy_kernel(krn)
    
  end subroutine init_history


  subroutine destroy_convolution(cnv)
    ! DESTROY_CONVOLUTION destroys derived type cnv
    ! 
    ! Modified: 23 January 2007

    use history, only : destroy_history

    implicit none

    ! I/O Parameters:
    ! CNV = convolution variables

    type(convolution_type),intent(inout) :: cnv

    ! Internal Parameters:
    ! I = index in kx direction
    ! J = index in ky direction

    integer(pin) :: i,j
    
    ! deallocate memory assigned to allocatable arrays

    if (allocated(cnv%kx)) deallocate(cnv%kx)
    if (allocated(cnv%ky)) deallocate(cnv%ky)

    if (allocated(cnv%Fx)) deallocate(cnv%Fx)
    if (allocated(cnv%Fy)) deallocate(cnv%Fy)
    if (allocated(cnv%FxH)) deallocate(cnv%FxH)
    if (allocated(cnv%FyH)) deallocate(cnv%FyH)

    if (allocated(cnv%Dx)) deallocate(cnv%Dx)
    if (allocated(cnv%Dy)) deallocate(cnv%Dy)       

    if (allocated(cnv%Fxp)) deallocate(cnv%Fxp)       
    if (allocated(cnv%Fxm)) deallocate(cnv%Fxm)       
    if (allocated(cnv%Fyp)) deallocate(cnv%Fyp)       
    if (allocated(cnv%Fym)) deallocate(cnv%Fym)       
    if (allocated(cnv%Fzp)) deallocate(cnv%Fzp)       
    if (allocated(cnv%Fzm)) deallocate(cnv%Fzm)       
    if (allocated(cnv%FxHp)) deallocate(cnv%FxHp)       
    if (allocated(cnv%FxHm)) deallocate(cnv%FxHm)       
    if (allocated(cnv%FyHp)) deallocate(cnv%FyHp)       
    if (allocated(cnv%FyHm)) deallocate(cnv%FyHm)       
    if (allocated(cnv%FzHp)) deallocate(cnv%FzHp)       
    if (allocated(cnv%FzHm)) deallocate(cnv%FzHm)       

    if (allocated(cnv%Dxp)) deallocate(cnv%Dxp)       
    if (allocated(cnv%Dxm)) deallocate(cnv%Dxm)       
    if (allocated(cnv%Dyp)) deallocate(cnv%Dyp)       
    if (allocated(cnv%Dym)) deallocate(cnv%Dym)       
    if (allocated(cnv%Dzp)) deallocate(cnv%Dzp)       
    if (allocated(cnv%Dzm)) deallocate(cnv%Dzm)       

    ! loop over modes, and destroy the history type associated with
    ! each mode, then deallocate/nullify pointer holding them

    if (associated(cnv%H)) then
       do j = 1,cnv%nky
          do i = cnv%mkx,cnv%pkx
             call destroy_history(cnv%H(i,j))
          end do
       end do       
       deallocate(cnv%H)       
       cnv%H => null()
    end if

  end subroutine destroy_convolution


  subroutine wavenumber(mdl,cnv,i,j,k,km,kx,ky)
    ! WAVENUMBER evaluates wavenumbers associated with a given Fourier mode
    ! 
    ! Modified: 23 January 2007

    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! I = index in the kx direction
    ! J = index in the ky direction
    ! K = wavenumber amplitude (k**2=kx**2+ky**2)
    ! KM = wavenumber product array 
    !      km(1,1) = kx**2
    !      km(1,2) = km(2,1) = kx*ky
    !      km(2,2) = ky**2
    ! KX = wavenumber in x direction
    ! KY = wavenumber in y direction

    ! Note that wavenumber storage is different for 2D mixed mode case.
    ! In this case, the x direction is the only spatial direction.  However, kx
    ! and ky do not measure wavenumbers parallel and perpendicular to
    ! this spatial direction, but instead measure wavenumbers in a coordinate
    ! system that is inclined at an angle mdl%angle with respect to it.  In this
    ! case, the index i measures distance along the (spatial) x axis, and is the
    ! appropriate index for both kx and ky.

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(in) :: cnv
    integer(pin),intent(in) :: i,j
    real(pr),intent(out) :: k,km(2,2),kx,ky

    kx = cnv%kx(i)

    if(mdl%dim/=3) then
       ky = cnv%ky(i)
    else
       ky = cnv%ky(j)
    end if

    k = sqrt(kx**2+ky**2)

    km(1,1) = kx**2
    km(1,2) = kx*ky
    km(2,1) = km(1,2)
    km(2,2) = ky**2

  end subroutine wavenumber


end module convolution
