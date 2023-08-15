! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module model

  ! MODEL contains variables and routines associated with model:
  ! grid spacing, number of grid points, elastic constants, spatial 
  ! coordinates, etc.
  ! 
  ! Modified: 6 August 2010

  use constants, only : pr,pin
  use runge_kutta, only : rk_type

  implicit none
  
  ! MODEL_TYPE is a derived type containing model variables
  !
  ! Parameters:
  ! DIM = dimension
  !      0 = spring-block model
  !      1 = single Fourier mode
  !      2 = 2D
  !      3 = 3D
  ! TRANSVERSE_SLIP = allow vector slip
  !      T = vector slip (in both x and y directions)
  !      F = scalar slip (in only x direction)
  ! OPENING = allow opening (only for bimaterial problem)
  !      T = allow fault opening
  !      F = do not allow fault opening
  ! MODE = crack geometry mode, valid only for dim=1 or 2
  !      2 = mode II
  !      3 = mode III
  !      0 = mixed mode
  ! ANGLE = mixed mode angle (0<=angle<=pi/2)
  !      0 = mode II
  !      pi/2 = mode III
  ! BM = bimaterial problem
  !      T = bimaterial
  !      F = identical materals
  ! SYM_Y = symmetry in y-direction
  !      T = symmetry about y=y0
  !      F = no symmetry about y=y0
  ! MU = shear modulus
  ! CS = S-wave speed
  ! CP = P-wave speed
  ! ETA = ratio of P-wave speed to S-wave speed
  ! NU = Poisson's ratio
  ! MUP = shear modulus, plus side
  ! MUM = shear modulus, minus side
  ! CSP = S-wave speed, plus side
  ! CSM = S-wave speed, minus side
  ! CPP = S-wave speed, plus side
  ! CPM = S-wave speed, minus side
  ! NUP = Poisson's ratio, plus side
  ! NUM = Poisson's ratio, minus side
  ! ETAP = ratio of P-wave speed to S-wave speed, plus side
  ! ETAM = ratio of P-wave speed to S-wave speed, minus side
  ! RRHO = ratio of densities, rhom/rhop
  ! RCS = ratio of S-wave speeds, csm/csp
  ! K = wavenumber of single Fourier mode (for dim=1 only)
  ! M = spring constant (for dim=0 spring-block model only)
  ! X0 = center of coordinate system in x direction
  ! Y0 = center of coordinate system in y direction
  ! NX = number of grid points in x direction (before refinement)
  ! NY = number of grid points in y direction (before refinement)
  ! NZ = number of grid points in z direction (before refinement)
  ! MX = lower bound on i (for particular process)
  ! PX = upper bound on i (for particular process)
  ! H = distance between spatial grid points in x and y directions
  ! DX = distance between spatial grid points in x direction
  ! DY = distance between spatial grid points in y direction
  ! DZ = distance between spatial grid points in z direction
  ! REFINE = refinement factor
  ! CASE =  labels some cases to simplify convolution routines and history
  ! vector storage (see init_model for details)
  ! MIXED2D = flag for 2D mixed mode, since wavenumber storage is different for this case
  !      T = 2D mixed mode
  !      F = not 2D mixed mode
  ! CV = coefficient of radiation damping term
  ! CVSP = coefficient of shear radiation damping term, plus side
  ! CVSM = coefficient of shear radiation damping term, minus side
  ! CVNP = coefficient of normal radiation damping term, plus side
  ! CVNM = coefficient of normal radiation damping term, minus side
  ! X = spatial coordinate in x direction
  ! Y = spatial coordinate in y direction
  ! Z = spatial coordinate in z direction
  ! POINT_AT_ZERO = place grid point exactly at zero (rather than centered around it)
  ! ZMAX = maximm value of z
  ! N = current time step
  ! NT = maximum number of time steps
  ! DT = time step
  ! DT_OLD = time step, previous time step (only used for substepping)
  ! DTMIN = minimum time step (for adaptive integration and substepping schemes)
  ! DTMAX = maximum time step (for adaptive integration and substepping schemes)
  ! T0 = initial time
  ! TMAX = maximum time
  ! T = current time
  ! CFL = CFL parameter
  ! METHOD = time integration method
  !      static = static elasticity solver
  !      iterative = iterative method
  !      substep = substepping within elastodynamic time step
  !      adaptive = embedded Runge-Kutta (for quasidynamic only)
  ! FRICTION_METHOD = friction law implementation
  !      strength = solve elasticity coupled with strength equation (using state variable)
  !      strength_rate = solve elasticity coupled with d(strength)/dt equation (no state variable)
  ! SUBSTEP_METHOD = substepping integration method
  !      iterative = iterative method
  !      adaptive = embedded Runge-Kutta
  !      rk = fixed steplength Runge-Kutta
  ! INTERP_METHOD = stress transfer interpolation method
  !      forward = forward predictor-corrector
  !      backward = backward predictor-corrector
  ! PREDICT_F_METHOD = method to predict stress transfer
  !      extrapolate = extrapolate from past values
  !      old = use old value
  !      integrate = explicitly integrate displacement/slip and convolve
  ! NS = number of stages/subdivisions within time step
  ! MS = lower bound on stress transfer stage index
  ! PS = upper bound on stress transfer stage index
  ! SUBSTEP_NT = number of substeps within elastodynamic time step for method='substep'
  ! NITER = number of iterations of elastodynamic time step
  ! ITER = iteration of current elastodynamic time step
  ! SUBSTEP_NITER = number of iterations of substep
  ! RTOL = relative error tolerance
  ! ATOL = absolute error tolerance
  ! SCHEME = Runge-Kutta scheme (see rk.f90 for schemes)
  ! RK = Runge-Kutta variables
  ! HOLDS_X = flag indicating if process holds x data
  ! SUBSTEP_INFO = print information about substepping

  type model_type
     integer(pin) :: dim,mode,nx,ny,nz,mx,px,case,n,nt,ns,ms,ps,substep_nt,niter,iter,substep_niter
     character(64) :: method,substep_method,interp_method,friction_method,predict_f_method,scheme
     real(pr) :: refine,k,M,h,dx,dy,dz,angle,x0,y0,dt,dt_old,dtmin,dtmax,cfl,t,t0,tmax, &
          mu,cs,cp,nu,eta,rrho,rcs,cV,mup,mum,csp,csm,cpp,cpm,nup,num,etap,etam, &
          cvsp,cvsm,cvnp,cvnm,rtol,atol,zmax
     logical :: transverse_slip,opening,mixed2d,bm,sym_y,holds_x,point_at_zero,substep_info
     real(pr),dimension(:),allocatable :: x,y,z
     type(rk_type),pointer :: rk=>null()
  end type model_type

contains
  

  subroutine read_model(ninput,mdl)
    ! READ_MODEL reads in model variables from file
    ! 
    ! Modified: 19 July 2010

    use constants, only : zero,half,one,three
    use io, only : error

    implicit none
    
    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! MDL = model variables

    type(model_type),intent(out) :: mdl
    integer(pin),intent(in) :: ninput
    
    ! Internal Parameters:
    ! DIM = dimension
    ! TRANSVERSE_SLIP = allow vector slip
    ! OPENING = allow opening (only for bimaterial problem)
    ! MODE = crack geometry mode (for dim=1 or 2)
    ! ANGLE = mixed mode angle
    ! BM = bimaterial problem
    ! SYM_Y = symmetry in y-direction
    ! MU = shear modulus
    ! CS = S-wave speed
    ! ETA = ratio of P-wave speed to S-wave speed
    ! MUP = shear modulus, plus side
    ! MUM = shear modulus, minus side
    ! CSP = S-wave speed, plus side
    ! CSM = S-wave speed, minus side
    ! ETAP = ratio of P-wave speed to S-wave speed, plus side
    ! ETAM = ratio of P-wave speed to S-wave speed, minus side
    ! RRHO = ratio of densities, rhom/rhop
    ! RCS = ratio of S-wave speeds, csm/csp
    ! K = wavenumber of single Fourier mode (for dim=1 only)
    ! M = spring constant (for dim=0 spring-block model only)
    ! CV = coefficient of radiation damping term
    ! X0 = center of coordinate system in x direction
    ! Y0 = center of coordinate system in y direction
    ! NX = number of grid points in x direction (before refinement)
    ! NY = number of grid points in y direction (before refinement)
    ! NZ = number of grid points in z direction (before refinement)
    ! H = distance between spatial grid points in x and y direction
    ! POINT_AT_ZERO = place grid point exactly at zero (rather than centered around it)
    ! DZ = distance between spatial grid points in z direction at the center
    ! REFINE = refinement factor
    ! METHOD = time integration method
    ! FRICTION_METHOD = friction law implementation
    ! SUBSTEP_METHOD = substepping integration method
    ! INTERP_METHOD = stress transfer interpolation method
    ! PREDICT_F_METHOD = method to predict stress transfer
    ! DT = time step
    ! DTMIN = minimum time step
    ! DTMAX = maximum time step
    ! NT = maximum number of time steps
    ! CFL = CFL parameter
    ! T0 = initial time
    ! TMAX = maximum time
    ! STAT = I/O error flag
    ! SUBSTEP_NT = number of substeps within elastodynamic time step for method='substep'
    ! NITER = maximum number of iterations of elastodynamic time step
    ! SUBSTEP_NITER = number of iterations of substep
    ! RTOL = relative error tolerance
    ! ATOL = absolute error tolerance
    ! SCHEME = Runge-Kutta scheme
    ! SUBSTEP_INFO = print information about substepping

    integer(pin) :: dim,mode,nx,ny,nz,nt,substep_nt,niter,substep_niter,stat
    character(64) :: method,substep_method,interp_method,friction_method,predict_f_method,scheme
    real(pr) :: refine,k,M,cV,h,dz,angle,x0,y0,dt,dtmin,dtmax,cfl,t0,tmax, &
         mu,cs,eta,mup,mum,csp,csm,etap,etam,rrho,rcs,rtol,atol,zmax
    logical :: transverse_slip,opening,bm,sym_y,point_at_zero,substep_info
    
    ! make namelist of user input variables

    namelist /model_list/ dim,bm,transverse_slip,opening,sym_y, &
         mode,angle,mu,cs,eta,mup,mum,csp,csm,etap,etam,rrho,rcs, &
         k,M,cV,x0,y0,nx,ny,nz,h,dz,refine,method,substep_method,friction_method,predict_f_method, &
         dt,dtmin,dtmax,nt,cfl,t0,tmax, &
         interp_method,substep_nt,niter,substep_niter,rtol,atol,scheme,point_at_zero,substep_info
    
    ! defaults
    
    dim = 2
    bm = .false.
    transverse_slip = .false.
    opening = .false.
    sym_y = .false.
    mode = 2
    angle = zero
    nx = 1
    ny = 1
    nz = 1
    mu = one
    cs = one
    eta = sqrt(three)
    mup = one
    mum = one
    csp = zero
    csm = zero
    etap = sqrt(three)
    etam = sqrt(three)
    rrho = one
    rcs = one
    k = one
    M = one
    cV = zero
    x0 = zero
    y0 = zero
    nx = 1
    ny = 1
    nz = 1
    h = zero
    dz = zero
    refine = one
    method = 'substep'
    substep_method = 'adaptive'
    interp_method = 'forward'
    friction_method = 'strength'
    predict_f_method = 'integrate'
    dt = zero
    dtmin = zero
    dtmax = huge(zero)
    nt = 0
    cfl = zero
    t0 = zero
    tmax = zero
    substep_nt = 1
    niter = 1
    substep_niter = 1
    rtol = 1.d-3
    atol = zero
    scheme = 'RK3(2)3'
    zmax = zero
    point_at_zero = .true.
    substep_info = .false.
    
    ! read namelist from input file, call error routine if needed
    
    rewind(ninput)
    read(ninput,nml=model_list,iostat=stat)
    if (stat/=0) call error("Error reading namelist 'model_list' in .in file",'read_model')
    
    ! assign input variables to components of derived type
    
    mdl%dim = dim
    mdl%bm = bm
    mdl%transverse_slip = transverse_slip
    mdl%opening = opening
    mdl%sym_y = sym_y
    mdl%mode = mode
    mdl%angle = angle
    mdl%mu = mu
    mdl%cs = cs
    mdl%eta = eta
    mdl%mup = mup
    mdl%mum = mum
    mdl%csp = csp
    mdl%csm = csm
    mdl%etap = etap
    mdl%etam = etam
    mdl%rrho = rrho
    mdl%rcs = rcs
    mdl%k = k
    mdl%M = M
    mdl%cV = cV
    mdl%x0 = x0
    mdl%y0 = y0
    mdl%nx = nx
    mdl%ny = ny
    mdl%nz = nz
    mdl%h = h
    mdl%dz = dz
    mdl%refine = refine
    mdl%method = method
    mdl%substep_method = substep_method
    mdl%interp_method = interp_method
    mdl%friction_method = friction_method
    mdl%predict_f_method = predict_f_method
    mdl%dt = dt
    mdl%dtmin = dtmin
    mdl%dtmax = dtmax
    mdl%nt = nt
    mdl%cfl = cfl
    mdl%t0 = t0
    mdl%tmax = tmax
    mdl%substep_nt = substep_nt
    mdl%niter = niter
    mdl%substep_niter = substep_niter
    mdl%rtol = rtol
    mdl%atol = atol
    mdl%scheme = scheme
    mdl%zmax = zmax
    mdl%point_at_zero = point_at_zero
    mdl%substep_info = substep_info

  end subroutine read_model


  subroutine init_model(necho,mdl)
    ! INIT_MODEL initializes model variables
    ! 
    ! Modified: 6 August 2010

    use constants, only : zero,half,one,two,pi
    use io, only : write_matlab,error
    use runge_kutta, only : init_rk
    use mpi_routines, only : is_master

    implicit none
    
    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! MDL = model variables

    integer(pin),intent(in) :: necho
    type(model_type),intent(inout) :: mdl
    
    ! Internal Parameters:
    ! I = index in x direction
    ! J = index in y direction
    ! K = index in z direction

    integer(pin) :: i,j,k

    ! set wave speeds for bimaterial problem based on density and wave-speed
    ! ratios, but only if these values were specified within the input file
    
    if (mdl%bm) then
       if (mdl%csm==zero.and.mdl%csp==zero) then
          mdl%csm = mdl%cs
          mdl%csp = mdl%cs/mdl%rcs
          mdl%mum = mdl%mu
          mdl%mup = mdl%mu/(mdl%rrho*mdl%rcs**2)
       else
          mdl%rcs = mdl%csm/mdl%csp
          mdl%rrho = mdl%mum/mdl%mup/mdl%rcs**2
       end if
    end if
    
    ! calculate wave speeds and Poisson's ratio
    
    if (mdl%bm) then
       mdl%cpp = mdl%etap*mdl%csp
       mdl%cpm = mdl%etam*mdl%csm
       mdl%nup = (mdl%etap**2-two)/(two*(mdl%etap**2-one))
       mdl%num = (mdl%etam**2-two)/(two*(mdl%etam**2-one))
    else
       mdl%cp = mdl%eta*mdl%cs
       mdl%nu = (mdl%eta**2-two)/(two*(mdl%eta**2-one))
    end if

    ! set plus and minus side material constants to identical materials values

    if (.not.mdl%bm) then
       mdl%csp = mdl%cs
       mdl%csm = mdl%cs
       mdl%cpp = mdl%cp
       mdl%cpm = mdl%cp
       mdl%mup = mdl%mu
       mdl%mum = mdl%mu
       mdl%nup = mdl%nu
       mdl%num = mdl%nu
       mdl%etap = mdl%eta
       mdl%etam = mdl%eta
    end if

    ! assign radiation damping coefficient
    
    mdl%cvsp = mdl%mup/mdl%csp
    mdl%cvsm = mdl%mum/mdl%csm
    mdl%cvnp = mdl%etap*mdl%mup/mdl%csp
    mdl%cvnm = mdl%etam*mdl%mum/mdl%csm
    if (mdl%dim/=0.and.(.not.mdl%bm)) mdl%cV = mdl%mu/mdl%cs/(one+mdl%rrho*mdl%rcs)
    
    ! if not 3D, set sym_y=F
    
    if (mdl%dim/=3) mdl%sym_y = .false.
    
    ! set mode=0 if 3D
    
    if (mdl%dim==3) mdl%mode = 0
    
    ! for the 1D and 2D pure mode cases,  set transverse_slip to 'F'
    
    if (mdl%dim/=3.and.mdl%mode/=0) mdl%transverse_slip = .false.
    
    ! assign convolution cases
    
    if (mdl%bm) then
       
       ! to simplify some cases for convolution, set mdl%case depending
       ! on the existence of transverse slip (uym/=uyp) and opening (uzm/=uzp)
       ! according to the following:
       ! transverse_slip=T,opening=T => case=1
       ! transverse_slip=T,opening=F => case=2
       ! transverse_slip=F,opening=T => case=3
       ! transverse_slip=F,opening=F => case=4
       
       if (mdl%transverse_slip) then
          if (mdl%opening) then
             mdl%case = 1
          else
             mdl%case = 2
          end if
       else
          if (mdl%opening) then
             mdl%case = 3
          else
             mdl%case = 4
          end if
       end if
       
    else
       
       ! to simplify some cases for convolution, set mdl%case depending
       ! on the existence of transverse slip (uy,Dy) and transverse stress
       ! transfer (fy,Fy) according to the following:
       ! Dy=T,Fy=T => case=1
       ! Dy=F,Fy=T => case=2
       ! Dy=F,Fy=F => case=3
       ! see the logic that follows for when these conditions arise
       
       select case(mdl%dim)
       case(1,2) ! 1D (single Fourier mode) and 2D
          select case(mdl%mode)
          case(0) ! mixed mode
             if (mdl%transverse_slip) then
                ! with transverse slip
                mdl%case = 1
             else
                ! without transverse slip
                mdl%case = 2
             end if
          case(2,3) ! pure mode
             mdl%case = 3
          end select
       case(3) ! 3D
          if (mdl%transverse_slip) then
             ! with transverse slip
             mdl%case = 1
          else
             ! without transverse slip
             mdl%case = 2
          end if
       end select
       
    end if
    
    ! set special flag for 2D mixed mode since wavenumber storage 
    ! is different in this case
    
    if (mdl%dim==2.and.mdl%mode==0) then
       mdl%mixed2d = .true.
    else
       mdl%mixed2d = .false.
    end if
    
    ! for single Fourier mode (dim=1), ensure consistency of h and k
    ! (such that k is Nyquist wavenumber for grid spacing h)
    
    if (mdl%dim==1) then
       if (mdl%h==zero) then
          ! set h from k
          mdl%h = pi/mdl%k
       else
          ! set k from h
          mdl%k = pi/mdl%h
       end if
    end if
    
    ! modify fault grid spacing according to refine and sym_y
    
    if (mdl%nx/=1) mdl%nx = ceiling(real(mdl%nx,pr)*mdl%refine)
    if (mdl%ny/=1) mdl%ny = ceiling(real(mdl%ny,pr)*mdl%refine)
    if (mdl%sym_y) then
       if (mod(mdl%ny,2)==0) then
          mdl%ny = mdl%ny/2
       else
          call error('ny must be even for sym_y to be used','init_model')
       end if
    end if
    if (mdl%dim/=1) mdl%h = mdl%h/mdl%refine
    mdl%dx = mdl%h
    mdl%dy = mdl%h
    if (mdl%nz/=1) mdl%nz = ceiling(real(mdl%nz-1,pr)*mdl%refine)+1
    mdl%dz = mdl%dz/mdl%refine
    
    ! set appropriate dimensions to unity for 1D and 2D
    
    select case(mdl%dim)
    case(0,1)
       mdl%nx = 1
       mdl%ny = 1
    case(2)
       mdl%ny = 1
    end select
    
    ! allocate memory to coordinate vectors
    
    allocate(mdl%x(mdl%nx))
    allocate(mdl%y(mdl%ny))
    allocate(mdl%z(mdl%nz))
    
    ! initialize spatial coordinate vectors with coordinates centered around x0 and y0
    
    ! x direction
    
    select case(mdl%dim)
       
    case(0,1)
       
       mdl%x(1) = mdl%x0

    case(2,3)

       if (mdl%nx==1) then
          mdl%x = mdl%x0
       else
          if (mod(mdl%nx,2)==0) then
             ! nx even
             mdl%x = real( (/ (i, i=-mdl%nx/2+1,mdl%nx/2) /) ,pr)
             if (.not.mdl%point_at_zero) mdl%x = mdl%x-half
          else
             ! nx odd
             mdl%x = real( (/ (i, i=-(mdl%nx-1)/2,(mdl%nx-1)/2) /) ,pr)
          end if
          mdl%x = mdl%x*mdl%dx+mdl%x0
       end if

    end select

    ! y direction

    select case(mdl%dim)
       
    case(0,1,2)
       
       mdl%y(1) = mdl%y0
       
    case(3)

       if (mdl%sym_y) then
          mdl%y = real( (/ (j, j=1,mdl%ny) /) ,pr)-half
       else
          if (mod(mdl%ny,2)==0) then
             ! ny even
             mdl%y = real( (/ (j, j=-mdl%ny/2+1,mdl%ny/2) /) ,pr)
             if (.not.mdl%point_at_zero) mdl%y = mdl%y-half
          else
             ! ny odd
             mdl%y = real( (/ (j, j=-(mdl%ny-1)/2,(mdl%ny-1)/2) /) ,pr)
          end if
       end if
       mdl%y = mdl%y*mdl%dy+mdl%y0
       
    end select

    ! initialize z direction coordinate vector
    
    mdl%z = real( (/ (k, k=0,mdl%nz-1) /) ,pr)*mdl%dz
    
    ! ensure that time step will be nonzero
    
    if (mdl%cfl==zero.and.mdl%dt==zero) mdl%cfl = 0.25_pr
    
    ! set time step or CFL ratio
    
    if (mdl%cfl/=zero) then
       ! set time step from CFL ratio
       if (mdl%bm) then
          mdl%dt = mdl%cfl*mdl%h/max(mdl%csp,mdl%csm)
       else
          mdl%dt = mdl%cfl*mdl%h/mdl%cs
       end if
    else
       ! set CFL from time step
       if (mdl%h==zero) mdl%h = one ! default grid spacing
       mdl%dt = mdl%dt/mdl%refine
       if (mdl%bm) then
          mdl%cfl = mdl%dt*max(mdl%csp,mdl%csm)/mdl%h
       else
          mdl%cfl = mdl%dt*mdl%cs/mdl%h
       end if
    end if
    
    if (mdl%dim==0.or.mdl%dim==1) then
       mdl%dt = mdl%dt/mdl%refine
       mdl%cfl = mdl%cfl/mdl%refine
    end if
    
    ! modify time steps according to refine
    
    if (mdl%method/='static') &
         mdl%nt = ceiling(real(mdl%nt,pr)*mdl%refine)
    
    ! make sure number of time steps and tmax are consistent
    
    if (mdl%nt==0) then
       if (mdl%tmax-mdl%t0<zero) then
          call error('Number of time steps or simulation time cannot be negative','init_model')
       else
          ! tmax given, set nt (but prevent numerical problems)
          if ((mdl%tmax-mdl%t0)/mdl%dt>real(huge(pin),pr)) then
             mdl%nt = huge(pin)
          else
             mdl%nt = ceiling((mdl%tmax-mdl%t0)/mdl%dt)
          end if
       end if
    else
       if (mdl%tmax==zero) then
          ! nt given, set tmax
          mdl%tmax = real(mdl%nt,pr)*mdl%dt+mdl%t0
       else
          ! both given, do nothing
       end if
    end if
    
    ! initialize time to initial time
    
    mdl%n = 0
    mdl%t = mdl%t0
    
    ! set number of subdivisions for time step and indices of stress transfer based on integration method
    
    select case(mdl%method)
    case('static')
       mdl%ns = 1
       mdl%ms = 0
       mdl%ps = 0
    case('iterative')
       mdl%ns = 2
       mdl%ms = 0
       mdl%ps = 0
    case('substep')
       select case(mdl%substep_method)
       case('iterative')
          mdl%ns = 3
       case('adaptive','rk')
          allocate(mdl%rk)
          call init_rk(mdl%scheme,mdl%rk)
          mdl%ns = mdl%rk%s+1
       end select
       mdl%ms = -1
       mdl%ps = 1
    case('adaptive')
       allocate(mdl%rk)
       call init_rk(mdl%scheme,mdl%rk)
       mdl%ns = mdl%rk%s
       mdl%ms = 0
       mdl%ps = 0
    end select

    if (is_master) then
       
       ! output variable values into matlab file
       
       call write_matlab(necho,'dim',mdl%dim,'mdl')
       call write_matlab(necho,'bm',mdl%bm,'mdl')
       call write_matlab(necho,'sym_y',mdl%sym_y,'mdl')
       call write_matlab(necho,'transverse_slip',mdl%transverse_slip,'mdl')
       call write_matlab(necho,'mode',mdl%mode,'mdl')
       if (mdl%mode==0.and.mdl%dim/=3) &
            call write_matlab(necho,'angle',mdl%angle,'mdl')
       if (mdl%bm) then
          call write_matlab(necho,'mup',mdl%mup,'mdl')
          call write_matlab(necho,'mum',mdl%mum,'mdl')
          call write_matlab(necho,'csp',mdl%csp,'mdl')
          call write_matlab(necho,'csm',mdl%csm,'mdl')
          call write_matlab(necho,'cpp',mdl%cpp,'mdl')
          call write_matlab(necho,'cpm',mdl%cpm,'mdl')
          call write_matlab(necho,'nup',mdl%nup,'mdl')
          call write_matlab(necho,'num',mdl%num,'mdl')
          call write_matlab(necho,'etap',mdl%etap,'mdl')
          call write_matlab(necho,'etam',mdl%etam,'mdl')
          call write_matlab(necho,'opening',mdl%opening,'mdl')
          call write_matlab(necho,'cvsp',mdl%cvsp,'mdl')
          call write_matlab(necho,'cvsm',mdl%cvsm,'mdl')
          call write_matlab(necho,'cvnp',mdl%cvnp,'mdl')
          call write_matlab(necho,'cvnm',mdl%cvnm,'mdl')
       else
          call write_matlab(necho,'mu',mdl%mu,'mdl')
          call write_matlab(necho,'cs',mdl%cs,'mdl')
          call write_matlab(necho,'cp',mdl%cp,'mdl')
          call write_matlab(necho,'eta',mdl%eta,'mdl')
          call write_matlab(necho,'nu',mdl%nu,'mdl')
          call write_matlab(necho,'rrho',mdl%rrho,'mdl')
          call write_matlab(necho,'rcs',mdl%rcs,'mdl')
          call write_matlab(necho,'cV',mdl%cV,'mdl')
       end if
       call write_matlab(necho,'x0',mdl%x0,'mdl')
       call write_matlab(necho,'y0',mdl%y0,'mdl')
       call write_matlab(necho,'refine',mdl%refine,'mdl')
       call write_matlab(necho,'nx',mdl%nx,'mdl')
       call write_matlab(necho,'ny',mdl%ny,'mdl')
       call write_matlab(necho,'nz',mdl%nz,'mdl')
       call write_matlab(necho,'nt',mdl%nt,'mdl')
       if (mdl%dim==0) call write_matlab(necho,'M',mdl%M,'mdl')
       if (mdl%dim==1) call write_matlab(necho,'k',mdl%k,'mdl')
       if (mdl%dim>1) then
          call write_matlab(necho,'cfl',mdl%cfl,'mdl')
          call write_matlab(necho,'h',mdl%h,'mdl')
          call write_matlab(necho,'dx',mdl%dx,'mdl')
       end if
       if (mdl%dim==3) then
          call write_matlab(necho,'dy',mdl%dy,'mdl')
       end if
       call write_matlab(necho,'dz',mdl%dz,'mdl')
       call write_matlab(necho,'zmax',mdl%zmax,'mdl')
       call write_matlab(necho,'dt',mdl%dt,'mdl')
       if (mdl%dtmin/=zero) call write_matlab(necho,'dtmin',mdl%dtmin,'mdl')
       if (mdl%dtmax/=huge(zero)) call write_matlab(necho,'dtmax',mdl%dtmax,'mdl')
       call write_matlab(necho,'nt',mdl%nt,'mdl')
       call write_matlab(necho,'t0',mdl%t0,'mdl')
       call write_matlab(necho,'tmax',mdl%tmax,'mdl')
       call write_matlab(necho,'method',mdl%method,'mdl')
       select case(mdl%method)
       case('iterative')
          call write_matlab(necho,'interp_method',mdl%interp_method,'mdl')
          call write_matlab(necho,'niter',mdl%niter,'mdl')
       case('substep')
          call write_matlab(necho,'interp_method',mdl%interp_method,'mdl')
          call write_matlab(necho,'substep_method',mdl%substep_method,'mdl')
          call write_matlab(necho,'predict_f_method',mdl%predict_f_method,'mdl')
          call write_matlab(necho,'substep_nt',mdl%substep_nt,'mdl')
          call write_matlab(necho,'niter',mdl%niter,'mdl')
          call write_matlab(necho,'dtmin',mdl%dtmin,'mdl')
          select case(mdl%substep_method)
          case('iterative')
             call write_matlab(necho,'substep_niter',mdl%niter,'mdl')
          case('adaptive')
             call write_matlab(necho,'scheme',mdl%scheme,'mdl')
             call write_matlab(necho,'rtol',mdl%rtol,'mdl')
             call write_matlab(necho,'atol',mdl%atol,'mdl')
          end select
       case('adaptive')
          call write_matlab(necho,'scheme',mdl%scheme,'mdl')
          call write_matlab(necho,'rtol',mdl%rtol,'mdl')
          call write_matlab(necho,'atol',mdl%atol,'mdl')
          call write_matlab(necho,'dtmin',mdl%dtmin,'mdl')
       case('rk')
          call write_matlab(necho,'scheme',mdl%scheme,'mdl')
       end select
       call write_matlab(necho,'friction_method',mdl%friction_method,'mdl')

    end if

  end subroutine init_model


  subroutine destroy_model(mdl)
    ! DESTROY_MODEL destroys derived type mdl
    ! 
    ! Modified: 21 August 2007

    use runge_kutta, only : destroy_rk

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables

    type(model_type),intent(inout) :: mdl
    
    ! deallocate memory assigned to allocatable arrays and pointers

    if (allocated(mdl%x)) deallocate(mdl%x)
    if (allocated(mdl%y)) deallocate(mdl%y)
    if (allocated(mdl%z)) deallocate(mdl%z)

    if (associated(mdl%rk)) then
       call destroy_rk(mdl%rk)
       deallocate(mdl%rk)
    end if
    mdl%rk => null()

  end subroutine destroy_model


  subroutine get_bound_point(field,dir,val,ots_min,ots_max,mdl,imin,imax,frac)
    ! GET_BOUND_POINT returns indices imin and imax of an array x(low:high) 
    ! as well as frac and val such that one of the following is true:
    ! 1. val lies within range - imin and imax bound location of val,
    !    0<=frac<=1
    ! 2. val lies below range - imin=low and imax=low+1 and either
    !    a. ots_min=T => extrapolation: frac<0
    !    b. ots_min=F => set val=x(low): frac=0
    ! 3. val lies above range - imin=high-1 and imax=high and either
    !    a. ots_max=T => extrapolation: frac>1
    !    b. ots_max=F => set val=x(high): frac=1
    ! 
    ! Modified: 19 July 2010

    use constants, only : zero
    use utilities, only : search_binary

    implicit none

    ! I/O Parameters:
    ! FIELD = field
    ! DIR = coordinate direction
    ! VAL = desired value of coordinate
    ! OTS_MIN = allow extrapolation below coordinate range
    ! OTS_MAX = allow extrapolation above coordinate range
    ! FRAC = fraction of grid spacing at which value lies
    ! MDL = model parameters
    ! IMIN = index of coordinate providing lower bound
    ! IMAX = index of coordinate providing upper bound

    character(*),intent(in) :: field,dir
    logical,intent(in) :: ots_min,ots_max
    real(pr),intent(inout) :: val
    real(pr),intent(out) :: frac
    type(model_type),intent(in) :: mdl
    integer(pin),intent(out) :: imin,imax

    ! Internal Parameters:
    ! OTS_FLR = flag indicating whether or not val lies outside of range,
    ! for floor mode
    ! OTS_CLG = flag indicating whether or not val lies outside of range,
    ! for ceiling mode
    ! VAL_FLR = closest value of coordinate to val in floor mode
    ! VAL_CLG = closest value of coordinate to val in ceiling mode
    ! LOW = lowest index of coordinate vector
    ! HIGH = highest index of coordinate vector
    ! DATA = temporary vector holding values of coordinate vector
    
    logical :: ots_flr,ots_clg
    real(pr) :: val_flr,val_clg
    integer(pin) :: low,high
    real(pr),dimension(:),allocatable :: data

    if (field=='rerr'.or.field=='aerr') then
       imin = 0
       imax = 0
       val = zero
       frac = zero
       return
    end if

    select case(dir)
    case('x')
       low = 1
       high = mdl%nx
       allocate(data(low:high))
       data = mdl%x
    case('y')
       low = 1
       high = mdl%ny
       allocate(data(low:high))
       data = mdl%y
    case('z')
       select case(field)
       case('p','T')
          low = 1
          high = mdl%nz
          allocate(data(low:high))
          data = mdl%z
       case default
          imin = 1
          imax = 1
          val = zero
          frac = zero
          return
       end select
    end select

    imin = search_binary(data,val,'floor',ots_flr,z=val_flr)+low-1
    imax = search_binary(data,val,'ceiling',ots_clg,z=val_clg)+low-1
    if (ots_flr.and.(.not.ots_min)) val = val_flr ! case 2b
    if (ots_clg.and.(.not.ots_max)) val = val_clg ! case 3b

    if (val_clg/=val_flr) then
       frac = (val-val_flr)/(val_clg-val_flr)
    else
       frac = zero
    end if

    if (allocated(data)) deallocate(data)

  end subroutine get_bound_point


  subroutine get_bound_range(field,dir,val_min,val_max,ots_min,ots_max,mdl,imin,imax)
    ! GET_BOUND_RANGE returns indices imin and imax of an array x(low:high) such that
    ! x(imin:imax) has as its lower and upper bounds the values of x(:) closest to
    ! the desired bounds val_min and val_max.  Rounding up or down in the location of
    ! imin and imax is controlled by ots_min and ots_max (if F, then enforce strict bounds).
    ! If the desired bounds lie outside of the range of x(:), then use low or high
    ! for imin or imax.
    !
    ! Modified: 16 June 2006

    use constants, only : zero
    use utilities, only : search_binary

    implicit none

    ! I/O Parameters:
    ! FIELD = field
    ! DIR = coordinate direction
    ! VAL_MIN = minimum desired value of coordinate
    ! VAL_MAX = maximum desired value of coordinate
    ! OTS_MIN = allow extrapolation below coordinate range
    ! OTS_MAX = allow extrapolation above coordinate range
    ! MDL = model parameters
    ! IMIN = index of coordinate providing lower bound
    ! IMAX = index of coordinate providing upper bound

    character(*),intent(in) :: field,dir
    logical,intent(in) :: ots_min,ots_max
    real(pr),intent(inout) :: val_min,val_max
    type(model_type),intent(in) :: mdl
    integer(pin),intent(out) :: imin,imax

    ! Internal Parameters:
    ! OTS_FLR = flag indicating whether or not val lies outside of range,
    ! for floor mode
    ! OTS_CLG = flag indicating whether or not val lies outside of range,
    ! for ceiling mode
    ! VAL_FLR = closest value of coordinate to val in floor mode
    ! VAL_CLG = closest value of coordinate to val in ceiling mode
    ! LOW = lowest index of coordinate vector
    ! HIGH = highest index of coordinate vector
    ! DATA = temporary vector holding values of coordinate vector

    logical :: ots_flr,ots_clg
    real(pr) :: val_flr,val_clg
    integer(pin) :: low,high
    real(pr),dimension(:),allocatable :: data

    if (field=='rerr'.or.field=='aerr') then
       imin = 0
       imax = 0
       val_min = zero
       val_max = zero
       return
    end if

    select case(dir)
    case('x')
       low = 1
       high = mdl%nx
       allocate(data(low:high))
       data = mdl%x
    case('y')
       low = 1
       high = mdl%ny
       allocate(data(low:high))
       data = mdl%y
    case('z')
       select case(field)
       case('p','T')
          low = 1
          high = mdl%nz
          allocate(data(low:high))
          data = mdl%z
       case default
          imin = 0
          imax = 0
          val_min = zero
          val_max = zero
          return
       end select
    end select

    ! lower bound
    if (ots_min) then
       imin = search_binary(data,val_min,'floor',ots_flr,z=val_flr)+low-1
       val_min = val_flr
    else
       imin = search_binary(data,val_min,'ceiling',ots_clg,z=val_clg)+low-1
       if (ots_clg) then
          imin = search_binary(data,val_min,'floor',ots_flr,z=val_flr)+low-1
          val_min = val_flr
       else
          val_min = val_clg
       end if
    end if

    ! upper bound
    if (ots_max) then
       imax = search_binary(data,val_max,'ceiling',ots_clg,z=val_clg)+low-1
       val_max = val_clg
    else
       imax = search_binary(data,val_max,'floor',ots_flr,z=val_flr)+low-1
       if (ots_flr) then
          imax = search_binary(data,val_max,'ceiling',ots_clg,z=val_clg)+low-1
          val_max = val_clg
       else
          val_max = val_flr
       end if
    end if

    if (allocated(data)) deallocate(data)

  end subroutine get_bound_range


end module model
