! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module fault

  ! FAULT contains arrays of all of the fields (slip, stress, etc.) that 
  ! are defined on the fault.  It also holds the initial values of the fields 
  ! and calls routines associated with asperities in these fields.
  ! 
  ! Modified: 20 July 2010

  use constants, only : pin,pr
  use asperity, only : asperity_list
  use preslip, only : preslip_type

  implicit none

  ! FAULT_TYPE is a derived type containing fault variables
  !
  ! Parameters:
  ! PRESLIPPED = initiate with nonzero initial slip or not
  ! THERMPRES = use thermal pressurization model or not
  ! DYNLOAD = include dynamic (time-dependent) load or not
  ! POROELASTIC = calculate poroelastic response or not
  ! PLASTIC = permit plastic deformation of fault zone
  ! TEMPERATURE = whether or not temperature field will be used by friction law
  ! STATE = whether or not state variable will be used by friction law
  ! RERR = relative error; has three indices, from 1 (step n+1) to -1 (step n-1) 
  ! AERR = absolute error; has three indices, from 1 (step n+1) to -1 (step n-1)
  ! RERR_OLD = relative error, previous time step
  ! AERR_OLD = absolute error, previous time step
  ! QM = combination of poroelastic properties, minus side 
  ! QP = combination of poroelastic properties, plus side
  ! Q0 = combination of poroelastic properties
  ! QM0 = combination of poroelastic properties, minus side (initial value)
  ! QP0 = combination of poroelastic properties, plus side (initial value)
  ! Q00 = combination of poroelastic properties (initial value)
  ! MUFP = shear modulus within fault, plus side
  ! MUFM = shear modulus within fault, minus side
  ! NUFP = Poisson ratio within fault, plus side
  ! NUFM = Poisson ratio within fault, minus side
  ! BP = Skempton's coefficient, plus side
  ! BM = Skempton's coefficient, minus side
  ! ZP = sqrt(permeability*storage coefficient), plus side
  ! ZM = sqrt(permeability*storage coefficient), minus side
  ! MUFP0 = shear modulus within fault, plus side (initial value)
  ! MUFM0 = shear modulus within fault, minus side (initial value)
  ! NUFP0 = Poisson ratio within fault, plus side (initial value)
  ! NUFM0 = Poisson ratio within fault, minus side (initial value)
  ! BP0 = Skempton's coefficient, plus side (initial value)
  ! BM0 = Skempton's coefficient, minus side (initial value)
  ! ZP0 = sqrt(permeability*storage coefficient), plus side (initial value)
  ! ZM0 = sqrt(permeability*storage coefficient), minus side (initial value)
  ! YP = yielded or not, plus side
  ! YM = yielded or not, minus side
  ! MUIP = internal friction coefficient, plus side
  ! MUIM = internal friction coefficient, minus side
  ! COHP = cohesion, plus side
  ! COHM = cohesion, minus side
  ! BETAP = dilatancy factor, plus side
  ! BETAM = dilatancy factor, minus side
  ! WP = half-width of fault zone, plus side
  ! WM = half-width of fault zone, minus side
  ! SXZ0 = initial stress, xz component
  ! SYZ0 = initial stress, yz component
  ! SZZ0 = initial stress, zz component
  ! SXX0 = initial stress, xx component
  ! SYY0 = initial stress, yy component
  ! SXY0 = initial stress, xy component
  ! UX = slip in x direction
  ! UY = slip in y direction
  ! VX = slip velocity in x direction
  ! VY = slip velocity in y direction
  ! U = cumulative slip (integrated as line integral)
  ! V = magnitude of slip velocity
  ! O = opening velocity
  ! N = effective normal stress, positive in compression
  ! DV = magnitude of slip acceleration
  ! UXP = displacement in x direction, plus side
  ! UXM = displacement in x direction, minus side
  ! UYP = displacement in y direction, plus side
  ! UYM = displacement in y direction, minus side
  ! UZP = displacement in z direction, plus side
  ! UZM = displacement in z direction, minus side
  ! VXP = particle velocity in x direction, plus side
  ! VXM = particle velocity in x direction, minus side
  ! VYP = particle velocity in y direction, plus side
  ! VYM = particle velocity in y direction, minus side
  ! VZP = particle velocity in z direction, plus side
  ! VZM = particle velocity in z direction, minus side
  ! SX = stress in x direction
  ! SY = stress in y direction
  ! SZ = stress in z direction
  ! SX0 = load in x direction
  ! SY0 = load in y direction
  ! SZ0 = load in z direction, including pore pressure
  ! FXI = stress transfer in x direction (interpolated in time)
  ! FYI = stress transfer in y direction (interpolated in time)
  ! FXPI = stress transfer in x direction, plus side (interpolated in time)
  ! FXMI = stress transfer in x direction, minus side (interpolated in time)
  ! FYPI = stress transfer in y direction, plus side (interpolated in time)
  ! FYMI = stress transfer in y direction, minus side (interpolated in time)
  ! FZPI = stress transfer in z direction, plus side (interpolated in time)
  ! FZMI = stress transfer in z direction, minus side (interpolated in time)
  ! FX = stress transfer in x direction
  ! FY = stress transfer in y direction
  ! FXP = stress transfer in x direction, plus side
  ! FXM = stress transfer in x direction, minus side
  ! FYP = stress transfer in y direction, plus side
  ! FYM = stress transfer in y direction, minus side
  ! FZP = stress transfer in z direction, plus side
  ! FZM = stress transfer in z direction, minus side
  ! EXXP = fault-parallel strain, xx component, plus side
  ! EXXM = fault-parallel strain, xx component, minus side
  ! EXYP = fault-parallel strain, xy component, plus side
  ! EXYM = fault-parallel strain, xy component, minus side
  ! EYYP = fault-parallel strain, yy component, plus side
  ! EYYM = fault-parallel strain, yy component, minus side
  ! TFRONT = time at which the rupture front passes each point
  ! Q = state variable
  ! DQ = rate of change of state variable
  ! S = strength
  ! DS = rate of change of strength
  ! T0 = temperature on fault
  ! P0 = pore pressure on fault
  ! T = temperature
  ! DT = rate of change of temperature
  ! P = pore pressure
  ! DP = rate of change of pore pressure
  ! The first two indices of each array correspond to grid points in the
  ! x and y directions (except for temperature and pressure, wheter the first index runs in z direction
  ! and the next two correspond to x and y). Some fields have a final index, which is used to store
  ! values of the field at time locations between time steps (used in higher-order time stepping methods).
  ! SCALE_U = scale for displacement (used for relative error estimates)
  ! SCALE_Q = scale for state (used for relative error estimates)
  ! SCALE_T = scale for temperature (used for relative error estimates)
  ! SCALE_P = scale for pore pressure (used for relative error estimates)
  ! SCALE_S = scale for stress (used for relative error estimates)
  ! SCALE_VMIN = scale for velocity, minimum  (used for relative error estimates)
  ! SCALE_VMAX = scale for velocity, maximum  (used for relative error estimates)
  ! Set above SCALE_XXX to zero if current value of XXX is to be used for typical scale,
  ! but note that this causes problems when value of XXX approaches zero.
  ! V0 = initial slip velocity (used as guess in some methods, or as actual initial value in others)
  ! VXL = load point slip velocity in x direction
  ! VYL = load point slip velocity in y direction
  ! PS = preslip variables
  ! ASPERITY_FILE = add heterogeneity to fault fields
  !      T = use asperity array data file
  !      F = do not use asperity array data file
  ! FILENAME = asperity file name
  ! LIST = singly-linked list of asperities

  type fault_type
     logical :: preslipped,thermpres,dynload,poroelastic,plastic,temperature,state,asperity_file
     real(pr) :: qm0,qp0,q00,mufp0,mufm0,nufp0,nufm0,Bp0,Bm0,Zp0,Zm0, &
          muip,muim,cohp,cohm,betap,betam,wp,wm, &
          sxz0,syz0,szz0,sxx0,syy0,sxy0, &
          rerr(-1:1),aerr(-1:1),rerr_old(-1:1),aerr_old(-1:1),V0,Vxl,Vyl, &
          scale_U,scale_Q,scale_T,scale_p,scale_s,scale_Vmin,scale_Vmax
     real(pr),dimension(:,:),allocatable :: fxi,fyi,fxpi,fxmi,fypi,fymi,fzpi,fzmi, &
          exxp,exxm,exyp,exym,eyyp,eyym
     real(pr),dimension(:,:,:),allocatable :: U,V,dV, &
          Ux,Uy,Vx,Vy, & 
          uxp,vxp,uyp,vyp,uzp,vzp, &
          uxm,vxm,uym,vym,uzm,vzm, &
          sx,sy,sz, &
          fx,fy,fxp,fxm,fyp,fym,fzp,fzm, &
          Q,dQ,T0,p0,S,dS,O,N
     real(pr),dimension(:,:,:,:),allocatable :: T,dT,p,dp
     real(pr),dimension(:,:),allocatable :: sx0,sy0,sz0,tfront,qm,qp,q0,mufp,mufm,nufp,nufm,Bp,Bm,Zp,Zm,Yp,Ym
     character(64) :: filename
     type(asperity_list) :: list
     type(preslip_type),pointer :: ps=>null()
  end type fault_type


contains


  subroutine read_fault(ninput,mdl,flt)
    ! READ_FAULT reads in fault variables from file
    ! 
    ! Modified: 20 July 2010

    use constants, only : zero
    use asperity, only : read_asperity_list
    use io, only : error
    use model, only : model_type
    use preslip, only : read_preslip

    implicit none

    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! MDL = model variables
    ! FLT = fault variables

    integer(pin),intent(in) :: ninput
    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt

    ! Internal Parameters:
    ! PRESLIPPED = initiate with nonzero initial slip or not
    ! THERMPRES = use thermal pressurization model or not
    ! DYNLOAD = including dynamic (time-dependent) load or not
    ! POROELASTIC = calculate poroelastic response or not
    ! PLASTIC = permit plastic deformation of fault zone
    ! QM = combination of poroelastic properties, minus side
    ! QP = combination of poroelastic properties, plus side
    ! Q0 = combination of poroelastic properties
    ! MUFP = shear modulus within fault, plus side
    ! MUFM = shear modulus within fault, minus side
    ! NUFP = Poisson ratio within fault, plus side
    ! NUFM = Poisson ratio within fault, minus side
    ! MUIP = internal friction coefficient, plus side
    ! MUIM = internal friction coefficient, minus side
    ! COHP = cohesion, plus side
    ! COHM = cohesion, minus side
    ! BETAP = dilatancy factor, plus side
    ! BETAM = dilatancy factor, minus side
    ! WP = half-width of fault zone, plus side
    ! WM = half-width of fault zone, minus side
    ! SXZ0 = initial stress, xz component
    ! SYZ0 = initial stress, yz component
    ! SZZ0 = initial stress, zz component
    ! SXX0 = initial stress, xx component
    ! SYY0 = initial stress, yy component
    ! SXY0 = initial stress, xy component
    ! SCALE_U = scale for displacement (used for relative error estimates)
    ! SCALE_Q = scale for state (used for relative error estimates)
    ! SCALE_T = scale for temperature (used for relative error estimates)
    ! SCALE_P = scale for pore pressure (used for relative error estimates)
    ! SCALE_S = scale for stress (used for relative error estimates)
    ! SCALE_VMIN = scale for velocity, minimum  (used for relative error estimates)
    ! SCALE_VMAX = scale for velocity, maximum  (used for relative error estimates)
    ! V0 = initial slip velocity
    ! VXL = load point slip velocity in x direction
    ! VYL = load point slip velocity in y direction
    ! STAT = I/O error flag
    ! NDATA = number of fault perturbations
    ! ASPERITY_LIST = add heterogeneity with list
    !      T = use asperity list
    !      F = do not use asperity list
    ! ASPERITY_FILE = add heterogeneity with file
    ! FILENAME = asperity file name

    logical :: preslipped,thermpres,dynload,poroelastic,plastic,asperity_list,asperity_file
    real(pr) :: qm,qp,q0,mufp,mufm,nufp,nufm,Bp,Bm,Zp,Zm, &
         muip,muim,cohp,cohm,betap,betam,wp,wm, &
         sxz0,syz0,szz0,sxx0,syy0,sxy0, &
         scale_U,scale_Q,scale_T,scale_p,scale_s,scale_Vmin,scale_Vmax,V0,Vxl,Vyl
    integer(pin) :: stat,ndata
    character(64) :: filename

    ! make namelist of user input variables

    namelist /fault_list/ preslipped,thermpres,dynload,poroelastic,plastic, &
         qm,qp,q0,mufp,mufm,nufp,nufm,Bp,Bm,Zp,Zm, &
         muip,muim,cohp,cohm,betap,betam,wp,wm, &
         sxz0,syz0,szz0,sxx0,syy0,sxy0, &
         scale_U,scale_Q,scale_T,scale_p,scale_s,scale_Vmin,scale_Vmax,V0,Vxl,Vyl, &
         asperity_list,asperity_file,filename

    ! set default values
    
    preslipped = .false.
    dynload = .false.
    thermpres = .false.
    poroelastic = .false.
    plastic = .false.
    qm = zero
    qp = zero
    q0 = zero
    mufp = zero
    mufm = zero
    nufp = zero
    nufm = zero
    Bp = zero
    Bm = zero
    Zp = zero
    Zm = zero
    muip = zero
    muim = zero
    cohp = zero
    cohm = zero
    betap = zero
    betam = zero
    wp = zero
    wm = zero
    sxz0 = zero
    syz0 = zero
    szz0 = zero
    sxx0 = zero
    syy0 = zero
    sxy0 = zero
    scale_U = zero
    scale_Q = zero
    scale_T = zero
    scale_p = zero
    scale_s = zero
    scale_Vmin = zero
    scale_Vmax = zero
    V0 = zero
    Vxl = zero
    Vyl = zero
    asperity_list = .false.
    asperity_file = .false.
    
    ! read namelist from input file, call error routine if needed
    
    rewind(ninput)
    read(ninput,nml=fault_list,iostat=stat)
    if (stat/=0) call error("Error reading namelist 'fault_list' in .in file",'read_fault')
    
    ! assign input variables to components of derived type
    
    flt%preslipped = preslipped
    flt%dynload = dynload
    flt%thermpres = thermpres
    flt%poroelastic = poroelastic
    flt%plastic = plastic
    flt%qm0 = qm
    flt%qp0 = qp
    flt%q00 = q0
    flt%mufp0 = mufp
    flt%mufm0 = mufm
    flt%nufp0 = nufp
    flt%nufm0 = nufm
    flt%Bp0 = Bp
    flt%Bm0 = Bm
    flt%Zp0 = Zp
    flt%Zm0 = Zm
    flt%muip = muip
    flt%muim = muim
    flt%cohp = cohp
    flt%cohm = cohm
    flt%betap = betap
    flt%betam = betam
    flt%wp = wp
    flt%wm = wm
    flt%sxz0 = sxz0
    flt%syz0 = syz0
    flt%szz0 = szz0
    flt%sxx0 = sxx0
    flt%syy0 = syy0
    flt%sxy0 = sxy0
    flt%scale_U = scale_U
    flt%scale_Q = scale_Q
    flt%scale_T = scale_T
    flt%scale_p = scale_p
    flt%scale_s = scale_s
    flt%scale_Vmin = scale_Vmin
    flt%scale_Vmax = scale_Vmax
    flt%V0 = V0
    flt%Vxl = Vxl
    flt%Vyl = Vyl
    flt%asperity_file = asperity_file
    flt%filename = filename
    
    ! asperity list
    
    ! set ndata (sx0,sy0,sz0, and optionally Q,T)
    
    ndata = 3
    if (flt%state) ndata = ndata+1
    if (flt%temperature.or.flt%thermpres) ndata = ndata+1
    
    if (asperity_list) &
         flt%list = read_asperity_list(ninput,'fault',ndata)

    ! read static crack variables if needed
    
    if (flt%preslipped) then
       allocate(flt%ps)
       call read_preslip(ninput,mdl,flt%ps)
    end if
    
  end subroutine read_fault


  subroutine init_fault(necho,mdl,flt)
    ! INIT_FAULT initializes fault variables
    ! 
    ! Modified: 20 July 2010

    use constants, only : zero,one
    use model, only : model_type
    use asperity, only : assign_list_data,destroy_asperity_list
    use io, only : write_matlab,file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use mpi_routines, only : is_master,MPI_REAL_PR,subarray
    use mpi

    implicit none

    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! MDL = model variables
    ! FLT = fault variables

    integer(pin),intent(in) :: necho
    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt

    ! Internal Parameters:
    ! IFIELD = index of field in asperity list input
    ! DARRAY = distributed array type
    ! FH = file handle

    integer(pin) :: ifield,darray
    type(file_distributed) :: fh

    ! allocate memory to fault arrays

    allocate(flt%sx0(mdl%mx:mdl%px,mdl%ny))
    allocate(flt%sy0(mdl%mx:mdl%px,mdl%ny))
    allocate(flt%sz0(mdl%mx:mdl%px,mdl%ny))

    allocate(flt%sx(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%sy(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%sz(mdl%mx:mdl%px,mdl%ny,mdl%ns))

    allocate(flt%U(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%V(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%O(mdl%mx:mdl%px,mdl%ny,mdl%ns))

    allocate(flt%tfront(mdl%mx:mdl%px,mdl%ny)) ! MAKE TFRONT OPTIONAL

    allocate(flt%Q( mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%dQ(mdl%mx:mdl%px,mdl%ny,mdl%ns))

    allocate(flt%dV(mdl%mx:mdl%px,mdl%ny,mdl%ns))

    allocate(flt%S( mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%dS(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%N( mdl%mx:mdl%px,mdl%ny,mdl%ns))

    if (flt%thermpres.or.flt%temperature) allocate(flt%T0(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    if (flt%thermpres.or.flt%poroelastic) allocate(flt%p0(mdl%mx:mdl%px,mdl%ny,mdl%ns))

    if (flt%poroelastic) then
       allocate(flt%exxp(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%exxm(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%exyp(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%exym(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%eyyp(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%eyym(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%qp  (mdl%mx:mdl%px,mdl%ny))
       allocate(flt%qm  (mdl%mx:mdl%px,mdl%ny))
       allocate(flt%q0  (mdl%mx:mdl%px,mdl%ny))
       allocate(flt%mufp(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%mufm(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%nufp(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%nufm(mdl%mx:mdl%px,mdl%ny))
       allocate(flt%Bp  (mdl%mx:mdl%px,mdl%ny))
       allocate(flt%Bm  (mdl%mx:mdl%px,mdl%ny))
       allocate(flt%Zp  (mdl%mx:mdl%px,mdl%ny))
       allocate(flt%Zm  (mdl%mx:mdl%px,mdl%ny))
       if (flt%plastic) then
          allocate(flt%Yp(mdl%mx:mdl%px,mdl%ny))
          allocate(flt%Ym(mdl%mx:mdl%px,mdl%ny))
       end if
    end if

    ! initialize fields

    flt%rerr = zero
    flt%aerr = zero

    flt%sx0 = zero
    flt%sy0 = zero
    flt%sz0 = zero

    flt%sx = zero
    flt%sy = zero
    flt%sz = zero

    if (allocated(flt%tfront)) flt%tfront = zero

    if (allocated(flt%T0)) flt%T0 = zero
    if (allocated(flt%p0)) flt%p0 = zero

    if (allocated(flt%exxp)) flt%exxp = zero
    if (allocated(flt%exxm)) flt%exxm = zero
    if (allocated(flt%exyp)) flt%exyp = zero
    if (allocated(flt%exym)) flt%exym = zero
    if (allocated(flt%eyyp)) flt%eyyp = zero
    if (allocated(flt%eyym)) flt%eyym = zero

    if (allocated(flt%qp)) flt%qp = flt%qp0
    if (allocated(flt%qm)) flt%qm = flt%qm0
    if (allocated(flt%q0)) flt%q0 = flt%q00

    ! set fault-zone properties to far-field values (unless already set)

    if (flt%poroelastic) then

       if (flt%mufp0==zero) flt%mufp0 = mdl%mup
       flt%mufp = flt%mufp0

       if (flt%mufm0==zero) flt%mufm0 = mdl%mum
       flt%mufm = flt%mufm0

       if (flt%nufp0==zero) flt%nufp0 = mdl%nup
       flt%nufp = flt%nufp0

       if (flt%nufm0==zero) flt%nufm0 = mdl%num
       flt%nufm = flt%nufm0

       flt%Bp = flt%Bp0
       flt%Bm = flt%Bm0
       flt%Zp = flt%Zp0
       flt%Zm = flt%Zm0

    end if

    if (allocated(flt%Yp)) flt%Yp = -huge(zero)
    if (allocated(flt%Ym)) flt%Ym = -huge(zero)

    flt%Q  = zero
    flt%dQ = zero

    flt%dV = zero

    flt%S  = zero
    flt%dS = zero

    flt%O = zero

    ! read asperity arrays from file
       
    if (flt%asperity_file) then
       call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
       call open_file_distributed(fh,flt%filename,'read',MPI_COMM_WORLD,darray,pr)
       call read_file_distributed(fh,flt%sx0)
       call read_file_distributed(fh,flt%sy0)
       call read_file_distributed(fh,flt%sz0)
       if (flt%state) call read_file_distributed(fh,flt%Q(:,:,1))
       if (allocated(flt%T0)) call read_file_distributed(fh,flt%T0(:,:,1))
       call close_file_distributed(fh)
    end if
       
    ! use list data to perturb fields
    
    call assign_list_data(flt%sx0,flt%list,1,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(flt%sy0,flt%list,2,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(flt%sz0,flt%list,3,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    ifield = 4
    if (flt%state) then
       call assign_list_data(flt%Q(:,:,1),flt%list,ifield,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       ifield = ifield+1
    end if
    if (allocated(flt%T0)) then
       call assign_list_data(flt%T0(:,:,1),flt%list,ifield,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       ifield = ifield+1
    end if
       
    ! deallocate memory for asperity list
    
    call destroy_asperity_list(flt%list)

    ! spread values across all stages (if needed)
    
    if (flt%state) flt%Q = spread(flt%Q(:,:,1),3,mdl%ns)
    if (allocated(flt%T0)) flt%T0  = spread(flt%T0(:,:,1),3,mdl%ns)

    ! add uniform contributions to initial stress

    flt%sx0 = flt%sx0+flt%sxz0
    flt%sy0 = flt%sy0+flt%syz0
    flt%sz0 = flt%sz0+flt%szz0

    ! allocate and initialize fields specific to formulation

    if (mdl%bm) then
       call init_fault_bm(mdl,flt)
    else
       call init_fault_im(mdl,flt)
    end if

    ! initialize strength and effective normal stress

    flt%S = sqrt(flt%sx**2+flt%sy**2)
    flt%N = -flt%sz
    if (flt%thermpres) flt%N = flt%N-flt%p0

    ! deallocate memory to preslip type
    
    if (associated(flt%ps)) deallocate(flt%ps)
    flt%ps => null()

    ! output variable values into matlab file

    if (is_master) then

       call write_matlab(necho,'preslipped',flt%preslipped,'flt')
       call write_matlab(necho,'dynload',flt%dynload,'flt')
       call write_matlab(necho,'thermpres',flt%thermpres,'flt')
       call write_matlab(necho,'poroelastic',flt%poroelastic,'flt')
       call write_matlab(necho,'plastic',flt%plastic,'flt')
       call write_matlab(necho,'scale_s',flt%scale_s,'flt')
       call write_matlab(necho,'scale_Vmin',flt%scale_Vmin,'flt')
       call write_matlab(necho,'scale_Vmax',flt%scale_Vmax,'flt')
       call write_matlab(necho,'V0',flt%V0,'flt')
       call write_matlab(necho,'sxz0',flt%sxz0,'flt')
       call write_matlab(necho,'syz0',flt%syz0,'flt')
       call write_matlab(necho,'szz0',flt%szz0,'flt')
       call write_matlab(necho,'sxx0',flt%sxx0,'flt')
       call write_matlab(necho,'syy0',flt%syy0,'flt')
       call write_matlab(necho,'szz0',flt%szz0,'flt')
       call write_matlab(necho,'muip',flt%muip,'flt')
       call write_matlab(necho,'muim',flt%muim,'flt')
       call write_matlab(necho,'cohp',flt%cohp,'flt')
       call write_matlab(necho,'cohm',flt%cohm,'flt')
       call write_matlab(necho,'betap',flt%betap,'flt')
       call write_matlab(necho,'betam',flt%betam,'flt')
       call write_matlab(necho,'wp',flt%wp,'flt')
       call write_matlab(necho,'wm',flt%wm,'flt')       
       if (mdl%dim==0) then
          call write_matlab(necho,'Vxl',flt%Vxl,'flt')
          call write_matlab(necho,'Vyl',flt%Vyl,'flt')
       end if
       if (flt%poroelastic.or.flt%plastic) then
          call write_matlab(necho,'qm',flt%qm0,'flt')
          call write_matlab(necho,'qp',flt%qp0,'flt')
          call write_matlab(necho,'q0',flt%q00,'flt')
          call write_matlab(necho,'mufp',flt%mufp0,'flt')
          call write_matlab(necho,'mufm',flt%mufm0,'flt')
          call write_matlab(necho,'nufp',flt%nufp0,'flt')
          call write_matlab(necho,'nufm',flt%nufm0,'flt')
          call write_matlab(necho,'Bp',flt%Bp0,'flt')
          call write_matlab(necho,'Bm',flt%Bm0,'flt')
          call write_matlab(necho,'Zp',flt%Zp0,'flt')
          call write_matlab(necho,'Zm',flt%Zm0,'flt')
      end if
       if (mdl%method=='adaptive'.or. &
          (mdl%method=='substep'.and. &
          (mdl%substep_method=='adaptive'.or.mdl%substep_method=='rk'))) then
          call write_matlab(necho,'scale_U',flt%scale_U,'flt')
          call write_matlab(necho,'scale_Q',flt%scale_Q,'flt')
          call write_matlab(necho,'scale_T',flt%scale_T,'flt')
          call write_matlab(necho,'scale_p',flt%scale_p,'flt')
       end if

    end if

  end subroutine init_fault


  subroutine init_fault_im(mdl,flt)
    ! INIT_FAULT_IM initializes fault variables, identical materials case
    ! 
    ! Modified: 3 December 2007

    use constants, only : zero
    use model, only : model_type
    use preslip, only : init_preslip_im

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt

    ! allocate memory to fault arrays

    allocate(flt%Ux(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%Uy(mdl%mx:mdl%px,mdl%ny,mdl%ns))

    allocate(flt%Vx(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%Vy(mdl%mx:mdl%px,mdl%ny,mdl%ns))

    allocate(flt%fx(mdl%mx:mdl%px,mdl%ny,mdl%ms:mdl%ps))
    allocate(flt%fy(mdl%mx:mdl%px,mdl%ny,mdl%ms:mdl%ps))

    allocate(flt%fxi(mdl%mx:mdl%px,mdl%ny))
    allocate(flt%fyi(mdl%mx:mdl%px,mdl%ny))

    ! initialize fields

    flt%Ux = zero
    flt%Uy = zero

    flt%Vx = flt%V0
    flt%Vy = zero

    flt%fx = zero
    flt%fy = zero

    flt%fxi = zero
    flt%fyi = zero

    ! calculate preslip fields

    if (flt%preslipped) &
         call init_preslip_im(mdl,flt%ps,flt%Ux(:,:,1),flt%Uy(:,:,1))

    ! spread values of initial slip and slip velocity

    flt%Ux = spread(flt%Ux(:,:,1),3,mdl%ns)
    flt%Uy = spread(flt%Uy(:,:,1),3,mdl%ns)

    flt%Vx = spread(flt%Vx(:,:,1),3,mdl%ns)
    flt%Vy = spread(flt%Vy(:,:,1),3,mdl%ns)

    flt%U = sqrt(flt%Ux**2+flt%Uy**2)
    flt%V = sqrt(flt%Vx**2+flt%Vy**2)

    ! set stress to ensure that fx=fy=0 at t=0-

    flt%sx = spread(flt%sx0,3,mdl%ns)-mdl%cV*flt%Vx
    flt%sy = spread(flt%sy0,3,mdl%ns)-mdl%cV*flt%Vy
    flt%sz = spread(flt%sz0,3,mdl%ns)

  end subroutine init_fault_im


  subroutine init_fault_bm(mdl,flt)
    ! INIT_FAULT_BM initializes fault variables, bimaterial case
    ! 
    ! Modified: 3 December 2007

    use constants, only : zero,one
    use model, only : model_type
    use preslip, only : init_preslip_bm

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt

    ! allocate memory to fault arrays

    allocate(flt%uxp(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%uxm(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%uyp(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%uym(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%uzp(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%uzm(mdl%mx:mdl%px,mdl%ny,mdl%ns))

    allocate(flt%vxp(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%vxm(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%vyp(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%vym(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%vzp(mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%vzm(mdl%mx:mdl%px,mdl%ny,mdl%ns))

    allocate(flt%fxp(mdl%mx:mdl%px,mdl%ny,mdl%ms:mdl%ps))
    allocate(flt%fxm(mdl%mx:mdl%px,mdl%ny,mdl%ms:mdl%ps))
    allocate(flt%fyp(mdl%mx:mdl%px,mdl%ny,mdl%ms:mdl%ps))
    allocate(flt%fym(mdl%mx:mdl%px,mdl%ny,mdl%ms:mdl%ps))
    allocate(flt%fzp(mdl%mx:mdl%px,mdl%ny,mdl%ms:mdl%ps))
    allocate(flt%fzm(mdl%mx:mdl%px,mdl%ny,mdl%ms:mdl%ps))

    allocate(flt%fxpi(mdl%mx:mdl%px,mdl%ny))
    allocate(flt%fxmi(mdl%mx:mdl%px,mdl%ny))
    allocate(flt%fypi(mdl%mx:mdl%px,mdl%ny))
    allocate(flt%fymi(mdl%mx:mdl%px,mdl%ny))
    allocate(flt%fzpi(mdl%mx:mdl%px,mdl%ny))
    allocate(flt%fzmi(mdl%mx:mdl%px,mdl%ny))

    ! initialize fields

    flt%uxp = zero
    flt%uxm = zero
    flt%uyp = zero
    flt%uym = zero
    flt%uzp = zero
    flt%uzm = zero

    flt%vxp =  flt%V0/(one+mdl%cvsp/mdl%cvsm)
    flt%vxm = -flt%V0/(one+mdl%cvsm/mdl%cvsp)
    flt%vyp = zero
    flt%vym = zero
    flt%vzp = zero
    flt%vzm = zero

    flt%fxp = zero
    flt%fxm = zero
    flt%fyp = zero
    flt%fym = zero
    flt%fzp = zero
    flt%fzm = zero

    flt%fxpi = zero
    flt%fxmi = zero
    flt%fypi = zero
    flt%fymi = zero
    flt%fzpi = zero
    flt%fzmi = zero

    ! calculate preslip fields

    if (flt%preslipped) &
         call init_preslip_bm(mdl,flt%ps,flt%uxp(:,:,1),flt%uxm(:,:,1),flt%uyp(:,:,1),flt%uym(:,:,1),flt%uzp(:,:,1),flt%uzm(:,:,1))

    ! spread values of initial slip and slip velocity

    flt%uxp = spread(flt%uxp(:,:,1),3,mdl%ns)
    flt%uxm = spread(flt%uxm(:,:,1),3,mdl%ns)
    flt%uyp = spread(flt%uyp(:,:,1),3,mdl%ns)
    flt%uym = spread(flt%uym(:,:,1),3,mdl%ns)
    flt%uzp = spread(flt%uzp(:,:,1),3,mdl%ns)
    flt%uzm = spread(flt%uzm(:,:,1),3,mdl%ns)

    flt%vxp = spread(flt%vxp(:,:,1),3,mdl%ns)
    flt%vxm = spread(flt%vxm(:,:,1),3,mdl%ns)
    flt%vyp = spread(flt%vyp(:,:,1),3,mdl%ns)
    flt%vym = spread(flt%vym(:,:,1),3,mdl%ns)
    flt%vzp = spread(flt%vzp(:,:,1),3,mdl%ns)
    flt%vzm = spread(flt%vzm(:,:,1),3,mdl%ns)

    flt%U = sqrt((flt%uxp-flt%uxm)**2+(flt%uyp-flt%uym)**2+(flt%uzp-flt%uzm)**2)
    flt%V = sqrt((flt%vxp-flt%vxm)**2+(flt%vyp-flt%vym)**2+(flt%vzp-flt%vzm)**2)

    ! set stress to ensure that fx=fy=fz=0 at t=0-

    flt%sx = spread(flt%sx0,3,mdl%ns)-mdl%cvsp*flt%vxp
    flt%sy = spread(flt%sy0,3,mdl%ns)-mdl%cvsp*flt%vyp
    flt%sz = spread(flt%sz0,3,mdl%ns)-mdl%cvnp*flt%vzp

  end subroutine init_fault_bm
  

  subroutine destroy_fault(flt)
    ! DESTROY_FAULT destroys derived type flt
    ! 
    ! Modified: 5 December 2007

    implicit none

    ! I/O Parameters:
    ! FLT = fault variables

    type(fault_type),intent(inout) :: flt

    ! deallocate memory assigned to allocatable arrays

    if (allocated(flt%sx)) deallocate(flt%sx)
    if (allocated(flt%sy)) deallocate(flt%sy)
    if (allocated(flt%sz)) deallocate(flt%sz)

    if (allocated(flt%sx0)) deallocate(flt%sx0)
    if (allocated(flt%sy0)) deallocate(flt%sy0)
    if (allocated(flt%sz0)) deallocate(flt%sz0)

    if (allocated(flt%U)) deallocate(flt%U)
       
    if (allocated(flt%Ux)) deallocate(flt%Ux)
    if (allocated(flt%Uy)) deallocate(flt%Uy)

    if (allocated(flt%Vx)) deallocate(flt%Vx)
    if (allocated(flt%Vy)) deallocate(flt%Vy)

    if (allocated(flt%fx)) deallocate(flt%fx)
    if (allocated(flt%fy)) deallocate(flt%fy)

    if (allocated(flt%fxi)) deallocate(flt%fxi)
    if (allocated(flt%fyi)) deallocate(flt%fyi)

    if (allocated(flt%uxp)) deallocate(flt%uxp)
    if (allocated(flt%uxm)) deallocate(flt%uxm)
    if (allocated(flt%uyp)) deallocate(flt%uyp)
    if (allocated(flt%uym)) deallocate(flt%uym)
    if (allocated(flt%uzp)) deallocate(flt%uzp)
    if (allocated(flt%uzm)) deallocate(flt%uzm)

    if (allocated(flt%vxp)) deallocate(flt%vxp)
    if (allocated(flt%vxm)) deallocate(flt%vxm)
    if (allocated(flt%vyp)) deallocate(flt%vyp)
    if (allocated(flt%vym)) deallocate(flt%vym)
    if (allocated(flt%vzp)) deallocate(flt%vzp)
    if (allocated(flt%vzm)) deallocate(flt%vzm)

    if (allocated(flt%fxp)) deallocate(flt%fxp)
    if (allocated(flt%fxm)) deallocate(flt%fxm)
    if (allocated(flt%fyp)) deallocate(flt%fyp)
    if (allocated(flt%fym)) deallocate(flt%fym)
    if (allocated(flt%fzp)) deallocate(flt%fzp)
    if (allocated(flt%fzm)) deallocate(flt%fzm)

    if (allocated(flt%fxpi)) deallocate(flt%fxpi)
    if (allocated(flt%fxmi)) deallocate(flt%fxmi)
    if (allocated(flt%fypi)) deallocate(flt%fypi)
    if (allocated(flt%fymi)) deallocate(flt%fymi)
    if (allocated(flt%fzpi)) deallocate(flt%fzpi)
    if (allocated(flt%fzmi)) deallocate(flt%fzmi)

    if (allocated(flt%exxp)) deallocate(flt%exxp)
    if (allocated(flt%exxm)) deallocate(flt%exxm)
    if (allocated(flt%exyp)) deallocate(flt%exyp)
    if (allocated(flt%exym)) deallocate(flt%exym)
    if (allocated(flt%eyyp)) deallocate(flt%eyyp)
    if (allocated(flt%eyym)) deallocate(flt%eyym)

    if (allocated(flt%qp)) deallocate(flt%qp)
    if (allocated(flt%qm)) deallocate(flt%qm)
    if (allocated(flt%q0)) deallocate(flt%q0)

    if (allocated(flt%mufp)) deallocate(flt%mufp)
    if (allocated(flt%mufm)) deallocate(flt%mufm)
    if (allocated(flt%nufp)) deallocate(flt%nufp)
    if (allocated(flt%nufm)) deallocate(flt%nufm)

    if (allocated(flt%Bp)) deallocate(flt%Bp)
    if (allocated(flt%Bm)) deallocate(flt%Bm)
    if (allocated(flt%Zp)) deallocate(flt%Zp)
    if (allocated(flt%Zm)) deallocate(flt%Zm)

    if (allocated(flt%Yp)) deallocate(flt%Yp)
    if (allocated(flt%Ym)) deallocate(flt%Ym)

    if (allocated(flt%tfront)) deallocate(flt%tfront)

    if (allocated(flt%V)) deallocate(flt%V)
    if (allocated(flt%dV)) deallocate(flt%dV)

    if (allocated(flt%S)) deallocate(flt%S)
    if (allocated(flt%dS)) deallocate(flt%dS)

    if (allocated(flt%Q)) deallocate(flt%Q)
    if (allocated(flt%dQ)) deallocate(flt%dQ)

    if (allocated(flt%T0)) deallocate(flt%T0)
    if (allocated(flt%p0)) deallocate(flt%p0)

  end subroutine destroy_fault


  subroutine set_stress_transfer(mdl,flt,stage)
    ! SET_STRESS_TRANSFER sets current value of stress transfer
    ! 
    ! Modified: 7 February 2007

    use model, only : model_type

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! STAGE = stage at which f is returned

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    integer(pin),intent(in) :: stage

    if (mdl%bm) then

       flt%fxpi = flt%fxp(:,:,stage)
       flt%fxmi = flt%fxm(:,:,stage)
       flt%fypi = flt%fyp(:,:,stage)
       flt%fymi = flt%fym(:,:,stage)
       flt%fzpi = flt%fzp(:,:,stage)
       flt%fzmi = flt%fzm(:,:,stage)

    else

       flt%fxi = flt%fx(:,:,stage)
       flt%fyi = flt%fy(:,:,stage)

    end if

  end subroutine set_stress_transfer


  subroutine interp_stress_transfer(flt,mdl,t1,t2,t,method)
    ! INTERP_STRESS_TRANSFER interpolates stress transfer over times
    ! within a single elastodynamic time step
    ! 
    ! Modified: 14 February 2007

    use constants, only : half,two
    use io, only : error
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! FLT = fault variables
    ! MDL = model variables
    ! T1 = time at beginning of elastodynamic time step
    ! T2 = time at end of elastodynamic time step
    ! T = current time (between t1 and t2)
    ! METHOD = interpolation method

    type(fault_type),intent(inout) :: flt
    type(model_type),intent(in) :: mdl
    real(pr),intent(in) :: t1,t2,t
    character(*),intent(in) :: method

    ! Internal Parameters:
    ! Q = fraction of elastodynamic time step completed

    real(pr) :: q

    q = (t-t1)/(t2-t1)

    ! indices of f are:
    ! ... t-2dt t-dt t t+dt
    ! ...  -2    -1  0   1

    select case(method)
       
    case default

       call error('invalid interpolation method','interp_stress_transfer')

    case('constant') ! constant at f(t)
          
       if (mdl%bm) then
          flt%fxpi = flt%fxp(:,:,0)
          flt%fxmi = flt%fxm(:,:,0)
          flt%fypi = flt%fyp(:,:,0)
          flt%fymi = flt%fym(:,:,0)
          flt%fzpi = flt%fzp(:,:,0)
          flt%fzmi = flt%fzm(:,:,0)
       else
          flt%fxi  = flt%fx( :,:,0)
          flt%fyi  = flt%fy( :,:,0)
       end if
       
    case('linear_backward') ! linear using f(t-dt) and f(t)
          
       if (mdl%bm) then
          flt%fxpi = flt%fxp(:,:,0)+q*(flt%fxp(:,:,0)-flt%fxp(:,:,-1))
          flt%fxmi = flt%fxm(:,:,0)+q*(flt%fxm(:,:,0)-flt%fxm(:,:,-1))
          flt%fypi = flt%fyp(:,:,0)+q*(flt%fyp(:,:,0)-flt%fyp(:,:,-1))
          flt%fymi = flt%fym(:,:,0)+q*(flt%fym(:,:,0)-flt%fym(:,:,-1))
          flt%fzpi = flt%fzp(:,:,0)+q*(flt%fzp(:,:,0)-flt%fzp(:,:,-1))
          flt%fzmi = flt%fzm(:,:,0)+q*(flt%fzm(:,:,0)-flt%fzm(:,:,-1))
       else
          flt%fxi  = flt%fx( :,:,0)+q*(flt%fx( :,:,0)-flt%fx( :,:,-1))
          flt%fyi  = flt%fy( :,:,0)+q*(flt%fy( :,:,0)-flt%fy( :,:,-1))
       end if
       
    case('linear_forward') ! linear using f(t) and f(t+dt)
       
       if (mdl%bm) then
          flt%fxpi = flt%fxp(:,:,0)+q*(flt%fxp(:,:,1)-flt%fxp(:,:,0))
          flt%fxmi = flt%fxm(:,:,0)+q*(flt%fxm(:,:,1)-flt%fxm(:,:,0))
          flt%fypi = flt%fyp(:,:,0)+q*(flt%fyp(:,:,1)-flt%fyp(:,:,0))
          flt%fymi = flt%fym(:,:,0)+q*(flt%fym(:,:,1)-flt%fym(:,:,0))
          flt%fzpi = flt%fzp(:,:,0)+q*(flt%fzp(:,:,1)-flt%fzp(:,:,0))
          flt%fzmi = flt%fzm(:,:,0)+q*(flt%fzm(:,:,1)-flt%fzm(:,:,0))
       else
          flt%fxi  = flt%fx( :,:,0)+q*(flt%fx( :,:,1)-flt%fx( :,:,0))
          flt%fyi  = flt%fy( :,:,0)+q*(flt%fy( :,:,1)-flt%fy( :,:,0))
       end if
       
    case('quadratic') ! quadratic using f(t-dt), f(t), f(t+dt)
       
       if (mdl%bm) then
          flt%fxpi = flt%fxp(:,:,0)+half*q*(flt%fxp(:,:,1)-flt%fxp(:,:,-1))+ &
               half*q**2*(flt%fxp(:,:,-1)-two*flt%fxp(:,:,0)+flt%fxp(:,:,1))
          flt%fxmi = flt%fxm(:,:,0)+half*q*(flt%fxm(:,:,1)-flt%fxm(:,:,-1))+ &
               half*q**2*(flt%fxm(:,:,-1)-two*flt%fxm(:,:,0)+flt%fxm(:,:,1))
          flt%fypi = flt%fyp(:,:,0)+half*q*(flt%fyp(:,:,1)-flt%fyp(:,:,-1))+ &
               half*q**2*(flt%fyp(:,:,-1)-two*flt%fyp(:,:,0)+flt%fyp(:,:,1))
          flt%fymi = flt%fym(:,:,0)+half*q*(flt%fym(:,:,1)-flt%fym(:,:,-1))+ &
               half*q**2*(flt%fym(:,:,-1)-two*flt%fym(:,:,0)+flt%fym(:,:,1))
          flt%fzpi = flt%fzp(:,:,0)+half*q*(flt%fzp(:,:,1)-flt%fzp(:,:,-1))+ &
               half*q**2*(flt%fzp(:,:,-1)-two*flt%fzp(:,:,0)+flt%fzp(:,:,1))
          flt%fzmi = flt%fzm(:,:,0)+half*q*(flt%fzm(:,:,1)-flt%fzm(:,:,-1))+ &
               half*q**2*(flt%fzm(:,:,-1)-two*flt%fzm(:,:,0)+flt%fzm(:,:,1))
       else
          flt%fxi  = flt%fx( :,:,0)+half*q*(flt%fx( :,:,1)-flt%fx( :,:,-1))+ &
               half*q**2*(flt%fx( :,:,-1)-two*flt%fx( :,:,0)+flt%fx( :,:,1))
          flt%fyi  = flt%fy( :,:,0)+half*q*(flt%fy( :,:,1)-flt%fy( :,:,-1))+ &
               half*q**2*(flt%fy( :,:,-1)-two*flt%fy( :,:,0)+flt%fy( :,:,1))
       end if
       
    end select
    
  end subroutine interp_stress_transfer
 

  subroutine stress_invariants(mu,nu,B,sxz,syz,szz,szz0,exx,exy,eyy,sxx0,sxy0,syy0,N,T)
    ! STRESS_INVARIANTS calculates mean normal stress and deviatoric stress
    ! from set of stresses and strains and elastic properties
    ! 
    ! Modified: 6 December 2007

    use constants, only : one,two,three

    implicit none

    ! I/O Parameters:
    ! MU = shear modulus
    ! NU = Poisson's ratio
    ! B = Skempton's coefficient
    ! SXZ = stress, xz component
    ! SYZ = stress, yz component
    ! SZZ = stress, zz component
    ! SZZ0 = initial stress, zz component
    ! EXX = strain, xx component
    ! EXY = strain, xy component
    ! EYY = strain, yy component
    ! SXX0 = initial stress, xx component
    ! SXY0 = initial stress, xy component
    ! SYY0 = initial stress, yy component
    ! N = mean normal stress (from first invariant)
    ! T = deviatoric stress (from second invariant)

    real(pr),dimension(:,:),intent(in) :: mu,nu,B,sxz,syz,szz,szz0,exx,exy,eyy
    real(pr),intent(in) :: sxx0,sxy0,syy0
    real(pr),dimension(:,:),intent(out) :: N,T

    ! Internal Parameters:
    ! NX = number of grid points in x direction
    ! NY = number of grid points in y direction
    ! SXX = stress, xx component
    ! SYY = stress, yy component
    ! SXY = stress, xy component

    integer(pin) :: nx,ny
    real(pr),dimension(:,:),allocatable :: sxx,sxy,syy

    ! allocate arrays
    
    nx = size(sxz,1)
    ny = size(sxz,2)

    allocate(sxx(nx,ny),sxy(nx,ny),syy(nx,ny))

    ! calculate unknown stress components (initial+change)
    ! note use of szz-szz0 for change in szz on right-hand side

    sxx = sxx0+mu*two/(one-nu)*exx+nu/(one-nu)*((szz-szz0)+two*mu*eyy)
    syy = syy0+mu*two/(one-nu)*eyy+nu/(one-nu)*((szz-szz0)+two*mu*exx)
    sxy = sxy0+two*mu*exy

    ! mean normal stress, positive in compression
    ! N = -(1/3)*(skk0+(1-B)*(skk-skk0))
    ! where sij0 is initial effective stress,
    ! simplify using skk0+(1-B)*(skk-skk0) = (1-B)*skk+B*skk0

    N = -((one-B)*(sxx+syy+szz)+B*(sxx0+syy0+szz0))/three

    ! deviatoric stress

    ! 2*T^2 = s_ij*s_ij, sij=sigma_ij+N*delta_ij (deviatoric stress tensor)

    T = (sxx+N)**2+(syy+N)**2+(szz+N)**2+sxy**2+sxz**2+syz**2
    T = sqrt(T/two)

    ! deallocate arrays

    deallocate(sxx,sxy,syy)

  end subroutine stress_invariants


end module fault
    
