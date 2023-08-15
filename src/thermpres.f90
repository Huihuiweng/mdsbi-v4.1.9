! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module thermpres

  ! THERMPRES contains variables and routines related to the
  ! thermal pressurization friction law
  ! 
  ! Modified: 12 October 2011

  use constants, only : pr,pin

  implicit none


  ! THERMPRES_TYPE is a derived type containing variables related to 
  ! thermal pressurization friction law
  !
  ! Parameters:
  ! RHOC = volumetric heat capacity
  ! K = thermal conductivity
  ! ALPHATH = thermal diffusivity
  ! ALPHAHY = hydraulic diffusivity
  ! LAMBDA = pore pressure rise per unit temprature rise
  ! SHAPE = shape of strain distribution
  ! W = width of deformation zone
  ! PHI = Maximum inelastic porosity change at large slip
  ! DELTAD = Characteristic slip distance associated with the change in porosity
  ! BETA =  volumetric pore fluid storage coefficient beta = n*(beta_n+beta_f)
  ! ASPERITY_FILE = add heterogeneity to initial fields
  !      T = use asperity array data file
  !      F = do not use asperity array data file
  ! FILENAME = asperity file name

  type thermpres_type
     real(pr) :: rhoc,K,alphath,alphahy,Lambda,W
     real(pr) :: beta,deltaD,Phi
     character(64) :: shape
     logical :: asperity_file
     character(64) :: filename
  end type thermpres_type


contains


  subroutine read_thermpres(ninput,tp)
    ! READ_THERMPRES reads in thermpres variables from file
    ! 
    ! Modified: 12 October 2011

    use constants, only : zero,half,one
    use io, only : error

    implicit none
    
    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! TP = thermal pressurization variables

    integer(pin),intent(in) :: ninput
    type(thermpres_type),intent(out) :: tp

    ! Internal Parameters:
    ! STAT = I/O error flag
    ! RHOC = volumetric heat capacity
    ! K = thermal conductivity
    ! ALPHATH = theamal diffusivity
    ! ALPHAHY = hydraulic diffusivity
    ! LAMBDA = pore pressure rise per unit temerature rise
    ! SHAPE = shape of strain distribution
    ! W = width of deformation zone
    ! T_BND = temperature at remote boundary
    ! P_BND = pressure at remote boundary
    ! ASPERITY_FILE = add heterogeneity to initial fields
    !      T = use asperity array data file
    !      F = do not use asperity array data file
    ! FILENAME = asperity file name
    
    integer(pin) :: stat
    real(pr) :: rhoc,K,alphath,alphahy,Lambda,W
    real(pr) :: beta,deltaD,Phi
    character(64) :: shape,filename
    logical :: asperity_file

    ! make namelist of user input variables

    namelist /thermpres_list/ rhoc,K,alphath,alphahy,Lambda,shape,W,beta,deltaD,Phi,asperity_file,filename

    ! assign default values
    
    rhoc = zero
    K = zero
    alphath = one
    alphahy = one
    Lambda = one
    shape = 'plane'
    W = half
    beta   = one
    Phi    = zero
    deltaD = one
    asperity_file = .false.

    ! read namelist from input file, call error routine if needed
    !write(*,*) beta,Phi,deltaD
    rewind(ninput)
    read(ninput,nml=thermpres_list,iostat=stat)
    if (stat>0) call error("Error reading namelist 'thermpres_list' in .in file",'read_thermpres')
    
    ! assign input variables to components of derived type

    ! user only specifies one of rhoc and K (rhoc takes precedence if both given)
    if (rhoc/=zero) then
       tp%rhoc = rhoc
       tp%K = rhoc*alphath
    else
       tp%K = K
       tp%rhoc = K/alphath
    end if

    tp%alphath = alphath
    tp%alphahy = alphahy
    tp%Lambda = Lambda
    tp%shape = shape
    tp%W = W
    tp%beta   = beta
    tp%Phi    = Phi
    tp%deltaD = deltaD
    tp%asperity_file = asperity_file
    tp%filename = filename
    !write(*,*) beta,Phi,deltaD

  end subroutine read_thermpres


  subroutine init_thermpres(necho,mdl,flt,tp)
    ! INIT_THERMPRES initializes thermpres variables
    ! 
    ! Modified: 12 October 2011

    use constants, only : zero
    use io, only : write_matlab,file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use model, only : model_type
    use fault, only : fault_type
    use mpi_routines, only : is_master,MPI_REAL_PR,subarray
    use mpi

    implicit none
    
    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! MDL = model variables
    ! FLT = fault variables
    ! TP = thermal pressurization variables

    integer(pin),intent(in) :: necho
    type(model_type),intent(in) :: mdl
    type(thermpres_type),intent(in) :: tp
    type(fault_type),intent(inout) :: flt

    ! Internal Parameters:
    ! K = index in z direction
    ! DARRAY = distributed array type
    ! FH = file handle

    integer(pin) :: k,darray
    type(file_distributed) :: fh

    ! allocate and initialize diffusion grids

    allocate(flt%T( mdl%nz,mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%dT(mdl%nz,mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%p( mdl%nz,mdl%mx:mdl%px,mdl%ny,mdl%ns))
    allocate(flt%dp(mdl%nz,mdl%mx:mdl%px,mdl%ny,mdl%ns))
    
    ! initialize temperature with linear profile and pressure as zero

    flt%T = zero
    flt%dT = zero
    flt%p = zero
    flt%dp = zero

    ! read asperity arrays from file
    
    if (tp%asperity_file) then
       call subarray(mdl%nz,mdl%nx,mdl%ny,1,mdl%nz,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
       call open_file_distributed(fh,tp%filename,'read',MPI_COMM_WORLD,darray,pr)
       call read_file_distributed(fh,flt%T(:,:,:,1))
       call read_file_distributed(fh,flt%p(:,:,:,1))
       call close_file_distributed(fh)
    end if
       
    ! spread values across all stages
    flt%T = spread(flt%T(:,:,:,1),4,mdl%ns)
    flt%p = spread(flt%p(:,:,:,1),4,mdl%ns)

    ! add constant background

    do k = 1,mdl%nz
       flt%T(k,:,:,:) = flt%T(k,:,:,:)+flt%T0
    end do

    ! output variable values into matlab file

    if (is_master) then

       call write_matlab(necho,'K',tp%K,'tp')
       call write_matlab(necho,'alphath',tp%alphath,'tp')
       call write_matlab(necho,'alphahy',tp%alphahy,'tp')
       call write_matlab(necho,'Lambda',tp%Lambda,'tp')
       call write_matlab(necho,'shape',tp%shape,'tp')
       call write_matlab(necho,'W',tp%W,'tp')
       call write_matlab(necho,'beta',tp%beta,'tp')
       call write_matlab(necho,'deltaD',tp%deltaD,'tp')
       call write_matlab(necho,'Phi',tp%Phi,'tp')

    end if

  end subroutine init_thermpres


  subroutine destroy_thermpres(flt)
    ! DESTROY_THERMPRES destroys derived type tp
    ! 
    ! Modified: 19 September 2007

    use fault, only : fault_type

    implicit none

    ! I/O Parameters:
    ! FLT = fault variables

    type(fault_type),intent(inout) :: flt

    if (allocated(flt%T)) deallocate(flt%T)
    if (allocated(flt%dT)) deallocate(flt%dT)
    if (allocated(flt%p)) deallocate(flt%p)
    if (allocated(flt%dp)) deallocate(flt%dp)

  end subroutine destroy_thermpres


  subroutine thermpres_rate(mdl,flt,tp,stage)
    ! THERMPRES_RATE solves the 1D diffusion equations for thermal pressurization
    ! using explicit finite differences
    ! 
    ! Modified: 21 July 2007

    use model, only : model_type
    use fault, only : fault_type

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! TP = thermal pressurization variables
    ! STAGE = integration step

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(thermpres_type),intent(in) :: tp
    integer(pin),intent(in) :: stage

    ! Internal Parameters:
    ! I = index in x direction
    ! J = index in y direction
    ! Q = shear heating rate

    integer(pin) :: i,j
    real(pr),dimension(:,:),allocatable :: Q

    ! return if not including thermal pressurization model
    
    if (.not.flt%thermpres) return
       
    ! return if process holds no x data 

    if (.not.mdl%holds_x) return

    ! calculate shear heat

    allocate(Q(mdl%mx:mdl%px,mdl%ny))
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)

    !$OMP DO

    do j = 1,mdl%ny
       do i = mdl%mx,mdl%px
          
          Q(i,j) = shear_heat(flt,i,j,stage)
          
       end do
    end do
    
    !$OMP END DO

    !$OMP END PARALLEL
    
    ! loop over the fault plane and solve diffusion equations
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)

    !$OMP DO

    do j = 1,mdl%ny
       do i = mdl%mx,mdl%px
                 
          call thermpres_rate_point(mdl,flt%dT(:,i,j,stage),flt%dp(:,i,j,stage), &
               flt%T(:,i,j,stage),flt%p(:,i,j,stage),Q(i,j),flt%V(i,j,stage),flt%U(i,j,stage),tp)
          
       end do
    end do

    !$OMP END DO

    !$OMP END PARALLEL
    
    ! deallocate memory

    deallocate(Q)

  end subroutine thermpres_rate


  subroutine thermpres_rate_point(mdl,dT,dp,T,p,Q,V,U,tp)
    ! THERMPRES_RATE_POINT solves the 1D diffusion equations for thermal pressurization
    ! using explicit finite differences at one point on fault
    ! 
    ! Modified: 12 October 2011

    use model, only : model_type
    use constants, only : zero,half,one,two,twopi
    use asperity, only : asperity_type,ampl
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! DT = rate of change of temperature
    ! DP = rate of change of pore pressure
    ! T = temperature
    ! P = pore pressure
    ! Q = shear heating rate (slip velocity times shear stress)
    ! TP = thermal pressurization variables

    type(model_type),intent(in) :: mdl
    real(pr),dimension(:),intent(out) :: dT,dp
    real(pr),dimension(:),intent(in) :: T,p
    real(pr),intent(in) :: Q
    real(pr),intent(in) :: V,U
    type(thermpres_type),intent(in) :: tp

    ! Internal Parameters:
    ! M = length of T,p
    ! K = index in z direction
    ! SN = effective normal stress, positive in compression
    ! RT = alpha_th/dz**2
    ! RP = alpht_hy/dz**2
    ! DTDZ = derivative of T with respect to z at left boundary
    ! G = nondimensional shear strain rate
    ! A = normalization for strain rate distribution
    ! Z = distance normal to fault
    ! ASP = asperity, used to evaluate strain rate distribution
    ! C2 = coefficients of second-order approximation of Laplacian

    integer(pin) :: M,k
    real(pr) :: dTdz,A,rt,rp,rd
    real(pr),dimension(:),allocatable :: G,G_temp
    type(asperity_type),pointer :: asp=>null()

    M = size(T)
    ! shear heating, either as boundary condition for slip-on-plane model or as
    ! volumetric heat source term for finite width shear zone model

    select case(tp%shape)

    case('plane')

       ! heat flux at boundary

       dTdz = -half*Q/tp%K

       ! no volumetric source term

       dT = zero
   
    case default

       ! symmetry BC at z=0

       dTdz = zero

       ! volumetric source term

       allocate(G(M))
       allocate(G_temp(M))
       allocate(asp)

       asp%type = tp%shape
       asp%inside = .true. 
       asp%x0 = zero

       ! calculate nondimensionalized strain rate, G=W*(strain rate)/V,
       ! note that integral over z of G(z) must equal W, in order that
       ! integral over z of strain rate equals V

       ! calculate coefficient to properly normalize strain rate distribution

       select case(tp%shape)
       case default
          call error('Invalid shear zone shape','thermpres_rate')
       case('box')
          ! W is half-width for consistency with 'box' in asperity routines
          A = half
          asp%length_x = tp%W
          do k = 1,M
             G(k) = A*ampl(mdl%z(k),zero,mdl%dz,one,asp,dim=2)
          end do
       case('gaussian')
          ! W is standard deviation of gaussian (see asperity.f90 for exact formula)
          A = one/sqrt(twopi)
          asp%length_x = tp%W !length of asperity in x direction
          do k = 1,M
             G(k) = A*ampl(mdl%z(k),zero,mdl%dz,one,asp,dim=2)
          end do
       end select

       dT = G*Q/(tp%W*tp%rhoc)
       G_temp = G
       deallocate(G)
       deallocate(asp)
       asp => null()

    end select

    ! special case of adiabatic, undrained response (only works for finite width shear zone)
    ! case for nz = 1    
    if (M==1) then
       dp = tp%Lambda*dT
       return
    else
       dp = zero
    end if

    ! add contribution from diffusion of heat and fluid mass

    ! parameter combinations


    rt = tp%alphath/mdl%dz**2
    rp = tp%alphahy/mdl%dz**2
    rd = (1/tp%beta)*tp%Phi*(V/tp%deltaD)*exp(-U/tp%deltaD)
        
    ! diffusion terms in interior

    do k = 2,M-1
       dT(k) = dT(k)+rt*(T(k+1)-two*T(k)+T(k-1))
       dp(k) = dp(k)+rp*(p(k+1)-two*p(k)+p(k-1)) - rd*G_temp(k)
    end do

    ! and at left boundary

    dT(1) = dT(1)+rt*two*(T(2)-T(1)-mdl%dz*dTdz)
    dp(1) = dp(1)+rp*two*(p(2)-p(1)) - rd*G_temp(1)  ! dp/dz = 0 on fault (no fluid flow across it by symmetry)

    ! thermal pressurization term in pressure rate

    dp = dp+tp%Lambda*dT

    ! fix T and p at right boundary

    dT(M) = zero
    dp(M) = zero

  end subroutine thermpres_rate_point


  function shear_heat(flt,i,j,stage) result(Q)
    ! SHEAR_HEAT calculates shear heating rate
    ! 
    ! Modified: 6 November 2006

    use fault, only : fault_type

    implicit none

    ! I/O Parameters:
    ! FLT = fault variables   
    ! I = index in x direction
    ! J = index in y direction
    ! STAGE = integration sub-step

    type(fault_type),intent(inout) :: flt
    integer(pin),intent(in) :: i,j,stage
    real(pr) :: Q

    ! Internal Parameters:
    ! S = amplitude of shear stress

    real(pr) :: S

    S = sqrt(flt%sx(i,j,stage)**2+flt%sy(i,j,stage)**2)

    Q = s*flt%V(i,j,stage)

  end function shear_heat


end module thermpres

