! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module ratestate

  ! RATESTATE contains variables and routines related to the
  ! rate-and-state friction laws
  ! 
  ! Modified: 20 July 2010

  use constants, only : pr,pin
  use asperity, only : asperity_list

  implicit none

  ! RATESTATE_TYPE is a derived type containing variables related to 
  ! rate-and-state friction law
  !
  ! Parameters:
  ! ASPERITY_FILE = add heterogeneity to initial rate-and-state fields
  !      T = use asperity array data file
  !      F = do not use asperity array data file
  ! FILENAME = asperity file name
  ! LIST = singly-linked list of asperities
  ! FORM = which form of rate-and-state to use, options are:
  !      LN = linearized Dieterich-Ruina
  !      AL = Dieterich-Ruina
  !      RAL = regularized Dieterich-Ruina
  !      SL = slip law with low-speed weakening
  !      SF = slip law with flash heating (Vweak as parameter)
  !      RSF = regularized slip law with flash heating (Vweak as parameter)
  !      EW = Dieterich-Ruina with enhanced weakening
  !      FL = slip law with flash heating (Vweak from actual T)
  !      RFL = regularized slip law with flash heating (Vweak from actual T)
  ! A_OF_T = temperature-dependent direct effect parameter a
  !      T = a increases linearly with absolute T, a(i,j) is value at room temperature
  !      F = constant a, independent of temperature
  ! F0 = reference (or unweakened) friction coefficient
  ! V0 = reference sliding velocity
  ! A = direct effect parameter
  ! B = healing parameter b
  ! L = state evolution distance
  ! VW = weakening velocity
  ! FW = weakened friction coefficient
  ! TC = temperature scale, (shear contact strength)/(specific heat)
  ! TW = weakening temperature
  ! VTH = thermal velocity (thermal diffusivity)/(contact dimension)

  type ratestate_type
     logical :: asperity_file,a_of_T
     character(64) :: filename,form
     type(asperity_list) :: list
     real(pr),dimension(:,:),allocatable :: f0,V0,a,b,L,Vw,fw,Tc,Tw,Vth
  end type ratestate_type

contains


  subroutine read_ratestate(ninput,rs,temperature)
    ! READ_RATESTATE reads in rate-and-state variables from file
    ! 
    ! Modified: 20 July 2010

    use io, only : error
    use asperity, only : read_asperity_list

    implicit none
    
    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! RS = rate-and-state variables
    ! TEMPERATURE = whether or not friction law requires temperature field

    integer(pin),intent(in) :: ninput
    logical,intent(out) :: temperature
    type(ratestate_type),intent(out) :: rs

    ! Internal Parameters:
    ! STAT = I/O error flag
    ! ASPERITY_LIST = add heterogeneity to initial rs fields with list
    ! ASPERITY_FILE = add heterogeneity to initial rs fields from file
    ! FILENAME = asperity file name
    ! NDATA = number of field perturbations
    ! FORM = which form of rate-and-state to use
    ! A_OF_T = temperature-dependent direct effect parameter a

    integer(pin) :: stat,ndata
    logical :: asperity_list,asperity_file,a_of_T
    character(64) :: filename,form

    ! make namelist of user input variables

    namelist /ratestate_list/ asperity_list,asperity_file,filename,form,a_of_T

    ! default values
    
    form = 'SL'
    a_of_T = .false.
    asperity_list = .false.
    asperity_file = .false.
    
    ! read namelist from input file, call error routine if needed
    
    rewind(ninput)
    read(ninput,nml=ratestate_list,iostat=stat)
    if (stat/=0) call error("Error reading namelist 'ratestate_list' in .in file",'read_ratestate')
    
    ! assign input variables to components of derived type
    
    rs%asperity_file = asperity_file
    rs%filename = filename
    rs%form = adjustl(form)
    rs%a_of_T = a_of_T
    
    ! asperity list, set ndata based on form of friction law
    
    select case(rs%form)
    case('LN','AL','RAL','SL','RSL')
       ndata = 5 ! for f0,V0,a,b,L
    case('EW')
       ndata = 6 ! for f0,V0,a,b,L,Vw
    case('SF','RSF')
       ndata = 7 ! for f0,V0,a,b,L,fw,Vw
    case('FL','RFL')
       ndata = 9 ! for f0,V0,a,b,L,fw,Tc,Tw,Vth
    case default
       call error('Invalid form of rate-and-state friction','read_ratestate')
    end select
    
    if (asperity_list) &
         rs%list = read_asperity_list(ninput,'ratestate',ndata)

    ! set flag indicating if temperature will be required for particular friction law
    
    select case(rs%form)
    case('FL','RFL')
       temperature = .true.
    case default
       temperature = .false.
    end select
       
  end subroutine read_ratestate


  subroutine init_ratestate(necho,mdl,rs)
    ! INIT_RATESTATE initializes rate-and-state variables
    ! 
    ! Modified: 20 July 2010

    use constants, only : zero
    use io, only : write_matlab,file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use model, only : model_type
    use asperity, only : assign_list_data,destroy_asperity_list
    use mpi_routines, only : is_master,MPI_REAL_PR,subarray
    use mpi

    implicit none
    
    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! MDL = model variables
    ! RS = rate-and-state variables

    integer(pin),intent(in) :: necho
    type(model_type),intent(in) :: mdl
    type(ratestate_type),intent(inout) :: rs

    ! Internal Parameters:
    ! DARRAY = distributed array type
    ! FH = file handle

    integer(pin) :: darray
    type(file_distributed) :: fh

    ! allocate arrays

    select case(rs%form)

    case('LN','AL','RAL','SL','RSL')

       ! allocate memory for rate-and-state parameter arrays
       
       allocate(rs%f0(mdl%mx:mdl%px,mdl%ny))
       allocate(rs%V0(mdl%mx:mdl%px,mdl%ny))
       allocate(rs%a( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%b( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%L( mdl%mx:mdl%px,mdl%ny))
       
       ! initialize to zero

       rs%f0 = zero
       rs%V0 = zero
       rs%a  = zero
       rs%b  = zero
       rs%L  = zero
       
       ! read asperity arrays from file
       
       if (rs%asperity_file) then
          call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
          call open_file_distributed(fh,rs%filename,'read',MPI_COMM_WORLD,darray,pr)
          call read_file_distributed(fh,rs%f0)
          call read_file_distributed(fh,rs%V0)
          call read_file_distributed(fh,rs%a)
          call read_file_distributed(fh,rs%b)
          call read_file_distributed(fh,rs%L)
          call close_file_distributed(fh)
       end if
       
       ! use list data to perturb fields
       
       call assign_list_data(rs%f0,rs%list,1,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%V0,rs%list,2,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%a ,rs%list,3,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%b ,rs%list,4,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%L ,rs%list,5,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)

    case('EW')

       ! allocate memory for rate-and-state parameter arrays
       
       allocate(rs%f0(mdl%mx:mdl%px,mdl%ny))
       allocate(rs%V0(mdl%mx:mdl%px,mdl%ny))
       allocate(rs%a( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%b( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%L( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%Vw(mdl%mx:mdl%px,mdl%ny))
       
       ! initialize to zero

       rs%f0 = zero
       rs%V0 = zero
       rs%a  = zero
       rs%b  = zero
       rs%L  = zero
       rs%Vw = zero
       
       ! read asperity arrays from file
          
       if (rs%asperity_file) then
          call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
          call open_file_distributed(fh,rs%filename,'read',MPI_COMM_WORLD,darray,pr)
          call read_file_distributed(fh,rs%f0)
          call read_file_distributed(fh,rs%V0)
          call read_file_distributed(fh,rs%a)
          call read_file_distributed(fh,rs%b)
          call read_file_distributed(fh,rs%L)
          call read_file_distributed(fh,rs%Vw)
          call close_file_distributed(fh)
       end if
          
       ! use list data to perturb fields
          
       call assign_list_data(rs%f0,rs%list,1,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%V0,rs%list,2,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%a ,rs%list,3,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%b ,rs%list,4,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%L ,rs%list,5,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%Vw,rs%list,6,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)

    case('SF','RSF')

       ! allocate memory for rate-and-state parameter arrays
       
       allocate(rs%f0(mdl%mx:mdl%px,mdl%ny))
       allocate(rs%V0(mdl%mx:mdl%px,mdl%ny))
       allocate(rs%a( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%b( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%L( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%fw(mdl%mx:mdl%px,mdl%ny))
       allocate(rs%Vw(mdl%mx:mdl%px,mdl%ny))
       
       ! initialize to zero

       rs%f0 = zero
       rs%V0 = zero
       rs%a  = zero
       rs%b  = zero
       rs%L  = zero
       rs%fw = zero
       rs%Vw = zero

       ! read asperity arrays from file
       
       if (rs%asperity_file) then
         call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
          call open_file_distributed(fh,rs%filename,'read',MPI_COMM_WORLD,darray,pr)
          call read_file_distributed(fh,rs%f0)
          call read_file_distributed(fh,rs%V0)
          call read_file_distributed(fh,rs%a)
          call read_file_distributed(fh,rs%b)
          call read_file_distributed(fh,rs%L)
          call read_file_distributed(fh,rs%fw)
          call read_file_distributed(fh,rs%Vw)
          call close_file_distributed(fh)
       end if
          
       ! use list data to perturb fields
       
       call assign_list_data(rs%f0,rs%list,1,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%V0,rs%list,2,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%a ,rs%list,3,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%b ,rs%list,4,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%L ,rs%list,5,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%fw,rs%list,6,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%Vw,rs%list,7,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)

    case('FL','RFL')

       ! allocate memory for rate-and-state parameter arrays
       
       allocate(rs%f0( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%V0( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%a(  mdl%mx:mdl%px,mdl%ny))
       allocate(rs%b(  mdl%mx:mdl%px,mdl%ny))
       allocate(rs%L(  mdl%mx:mdl%px,mdl%ny))
       allocate(rs%fw( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%Tc( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%Tw( mdl%mx:mdl%px,mdl%ny))
       allocate(rs%Vth(mdl%mx:mdl%px,mdl%ny))

       ! initialize to zero

       rs%f0  = zero
       rs%V0  = zero
       rs%a   = zero
       rs%b   = zero
       rs%L   = zero
       rs%fw  = zero
       rs%Tc  = zero
       rs%Tw  = zero
       rs%Vth = zero
       
       ! read asperity arrays from file
          
       if (rs%asperity_file) then
          call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
          call open_file_distributed(fh,rs%filename,'read',MPI_COMM_WORLD,darray,pr)
          call read_file_distributed(fh,rs%f0)
          call read_file_distributed(fh,rs%V0)
          call read_file_distributed(fh,rs%a)
          call read_file_distributed(fh,rs%b)
          call read_file_distributed(fh,rs%L)
          call read_file_distributed(fh,rs%fw)
          call read_file_distributed(fh,rs%Tc)
          call read_file_distributed(fh,rs%Tw)
          call read_file_distributed(fh,rs%Vth)
          call close_file_distributed(fh)
       end if
          
       ! use list data to perturb fields
       
       call assign_list_data(rs%f0 ,rs%list,1,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%V0 ,rs%list,2,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%a  ,rs%list,3,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%b  ,rs%list,4,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%L  ,rs%list,5,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%fw ,rs%list,6,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%Tc ,rs%list,7,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%Tw ,rs%list,8,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
       call assign_list_data(rs%Vth,rs%list,9,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)

    end select

    ! deallocate memory for asperity list
       
    call destroy_asperity_list(rs%list)

    ! output variable values into matlab file

    if (is_master) then
       
       call write_matlab(necho,'form',rs%form,'rs')
    
    end if

  end subroutine init_ratestate


  subroutine destroy_ratestate(rs)
    ! DESTROY_RATESTATE destroys derived type rs
    ! 
    ! Modified: 10 June 2007

    implicit none

    ! I/O Parameters:
    ! RS = rate-and-state variables

    type(ratestate_type),intent(inout) :: rs

    ! deallocate memory assigned to allocatable arrays

    if (allocated(rs%f0)) deallocate(rs%f0)
    if (allocated(rs%V0)) deallocate(rs%V0)
    if (allocated(rs%a)) deallocate(rs%a)
    if (allocated(rs%b)) deallocate(rs%b)
    if (allocated(rs%L)) deallocate(rs%L)
    if (allocated(rs%Vw)) deallocate(rs%Vw)
    if (allocated(rs%fw)) deallocate(rs%fw)
    if (allocated(rs%Tc)) deallocate(rs%Tc)
    if (allocated(rs%Tw)) deallocate(rs%Tw)
    if (allocated(rs%Vth)) deallocate(rs%Vth)

  end subroutine destroy_ratestate


  subroutine fric_rs(rs,V,Q,T,i,j,f,dfdV,err)
    ! FRIC_RS returns the coefficient of friction for a rate-and-state friction law
    ! 
    ! Modified: 20 July 2010

    use constants, only : zero,half,one,two,four
    use utilities, only : arcsinh

    implicit none

    ! I/O Parameters:
    ! RS = rate-and-state variables
    ! V = slip velocity
    ! Q = state variable
    ! T = temperature
    ! I = index in x direction
    ! J = index in y direction
    ! F = coefficient of friction
    ! DFDV = derivative of friction with respect to velocity
    ! ERR = error flag (true if friction isn't defined at particular V)

    type(ratestate_type),intent(in) :: rs
    real(pr),intent(in) :: V,Q,T
    integer(pin),intent(in) :: i,j
    real(pr),intent(out) :: f
    real(pr),intent(out),optional :: dfdV
    logical,intent(out),optional :: err

    ! Internal  Parameters:
    ! WK = parameter combination in EW law
    ! O = parameter combination in regularized law
    ! A = direct effect parameter

    real(pr) :: wk,O,a

    ! check if V is acceptable (some laws require positive V)

    if (present(err)) then
       select case(rs%form)
       case default
          err = .false.
       case('AL','SL','SF','EW','FL')
          if (V<=zero) then
             err = .true.
             return
          else
             err = .false.
          end if
       end select
    end if

    ! set magnitude of direct effect

    if (rs%a_of_T) then
       a = rs%a(i,j)*(T/293._pr)
    else
       a = rs%a(i,j)
    end if

    ! evaluate friction coefficient

    select case(rs%form)
    case('LN')
       f = rs%f0(i,j)+a*V/rs%V0(i,j)+rs%b(i,j)*rs%V0(i,j)*Q/rs%L(i,j)
       if (present(dfdV)) dfdV = a/rs%V0(i,j)
    case('EW')
       wk = one+rs%L(i,j)/(Q*rs%Vw(i,j))
       f = (rs%f0(i,j)+a*log(V/rs%V0(i,j))+rs%b(i,j)*log(rs%V0(i,j)*Q/rs%L(i,j)))/wk
       if (present(dfdV)) dfdV = a/(V*wk)
    case('SL','SF','FL','AL')
       f = a*log(V/rs%V0(i,j))+Q
       if (present(dfdV)) dfdV = a/V
    case('RSL','RSF','RFL','RAL')
       O = half/rs%V0(i,j)*exp(Q/a)
       f = a*arcsinh(O*V)
       if (present(dfdV)) dfdV = a*O/sqrt(one+(V*O)**2)
    end select

  end subroutine fric_rs


  subroutine strength_rs(rs,V,Q,T,sz,i,j,SV,dSdV,err)
    ! STRENGTH_RS returns the strength for rate-and-state law
    ! 
    ! Modified: 19 September 2007

    use constants, only : zero

    implicit none

    ! I/O Parameters:
    ! RS = rate-and-state variables
    ! V = slip velocity
    ! Q = state variable
    ! T = temperature
    ! SZ = stress in z direction
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! SV = velocity-dependent strength
    ! DSDV = derivative of strength with respect to slip velocity
    ! ERR = error flag (true if strength isn't defined at particular V)

    type(ratestate_type),intent(in) :: rs
    real(pr),intent(in) :: V,Q,T,sz
    integer(pin),intent(in) :: i,j
    real(pr),intent(out) :: SV,dSdV
    logical,intent(out) :: err

    ! Internal  Parameters:
    ! F = coefficient of friction
    ! DFDV = derivative of f with respect to V
    ! N = normal stress, positive in compression

    real(pr) :: f,dfdV,N

    ! calculate normal stress

    N = max(-sz,zero)

    ! calculate coefficient of friction and derivative with respect to V

    call fric_rs(rs,V,Q,T,i,j,f,dfdV,err)

    ! set strength to product of coefficient of friction and normal stress

    SV = f*N
    dSdV = dfdV*N
    
  end subroutine strength_rs


  function evolution_rs(V,Q,T,i,j,rs) result(dQ)
    ! EVOLUTION_RS returns the state-rate
    ! 
    ! Modified: 20 July 2010

    use constants, only : one,two,pi

    implicit none

    ! I/O Parameters:
    ! V = slip velocity
    ! Q = state variable
    ! T = temperature
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! RS = rate-and-state variables
    ! DQ = state-rate

    real(pr),intent(in) :: V,Q,T
    integer(pin),intent(in) :: i,j
    type(ratestate_type),intent(in) :: rs
    real(pr) :: dQ

    ! Internal  Parameters:
    ! F = friction
    ! FSS = steady state friction
    ! VW = weakening velocity
    ! QSS = steady state state
    ! A = direct effect parameter

    real(pr) :: f,fss,Vw,Qss,a

    ! set magnitude of direct effect

    if (rs%a_of_T) then
       a = rs%a(i,j)*(T/293._pr)
    else
       a = rs%a(i,j)
    end if

    ! evaluate state evolution rate

    select case(rs%form)
    case('LN')
       dQ = -V/rs%V0(i,j)-rs%V0(i,j)*Q/rs%L(i,j)
    case('EW')
       dQ = one-V*Q/rs%L(i,j)
    case('AL','RAL') ! uses Q = f0+b*log(V0*theta/L)
       dQ = rs%b(i,j)*rs%V0(i,j)/rs%L(i,j)*(exp((rs%f0(i,j)-Q)/rs%b(i,j))-V/rs%V0(i,j))
    case('SL','RSL','SF','RSF','FL','RFL')
       call fric_rs(rs,V,Q,T,i,j,f)
       call fss_rs(rs,V,T,i,j,fss)
       dQ = V/rs%L(i,j)*(fss-f)
    !case('RSF') ! used in SCEC code validation problems
    !   call fss_rs(rs,V,T,i,j,fss)
    !   Qss = rs%a(i,j)*log(two*rs%V0(i,j)/V*sinh(fss/rs%a(i,j)))
    !   dQ = V/rs%L(i,j)*(Qss-Q)
    end select

  end function evolution_rs


  function check_state_rs(Q,rs) result(accept)
    ! CHECK_STATE_RS checks if state variable value is acceptable
    ! for the particular rate-and-state law being used (some do not
    ! permit negative state values)
    ! 
    ! Modified: 20 July 2010

    use constants, only : zero

    implicit none

    ! I/O Parameters:
    ! Q = state variable
    ! RS = rate-and-state variables
    ! ACCEPT = accept value of state variable or not

    real(pr),intent(in) :: Q
    type(ratestate_type),intent(in) :: rs
    logical :: accept

    accept = .true.
    
    select case(rs%form)
    case('EW') ! only laws that use "theta" as state variable
       if (Q<=zero) then
          accept = .false.
       end if
    end select

  end function check_state_rs
  

  subroutine directeffect_rs(rs,V,S,T,sz,i,j,D,dDdV,dDdS)
    ! DIRECTEFFECT_RS returns direct effect, D, of a rate-and-state law written in the form:
    ! dS/dt=D(V,S)*dV/dt+R(V,S)
    ! 
    ! Modified: 19 September 2007

    use constants, only : zero,one
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! RS = rate-and-state variables
    ! V = slip velocity
    ! S = strength
    ! T = temperature
    ! SZ = stress in z direction
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! D = direct effect (coefficient of dV/dt)
    ! DDDV = derivative of D(V,S) with respect to V
    ! DDDS = derivative of D(V,S) with respect to S

    type(ratestate_type),intent(in) :: rs
    real(pr),intent(in) :: V,S,T,sz
    integer(pin),intent(in) :: i,j
    real(pr),intent(out) :: D
    real(pr),intent(out),optional :: dDdV,dDdS

    ! Internal  Parameters:
    ! N = normal stress, positive in compression
    ! A = direct effect parameter
    ! AN = a*N

    real(pr) :: N,a,aN

    ! set magnitude of direct effect

    if (rs%a_of_T) then
       a = rs%a(i,j)*(T/293._pr)
    else
       a = rs%a(i,j)
    end if

    ! calculate normal stress

    N = max(-sz,zero)

    aN = a*N

    ! evaluate direct effect term

    select case(rs%form)
    case('SL','FL','SF')
       D = aN/V
       if (present(dDdV)) dDdV = -D/V
       if (present(dDdS)) dDdS = zero
    case('RFL','RSF')
       D = aN*tanh(S/aN)/V
       if (present(dDdV)) dDdV = -D/V
       if (present(dDdS)) dDdS = (aN/V)*(one-tanh(S/aN)**2)
    case default
       call error('Direct effect not yet defined for friction law: ' // &
            trim(adjustl(rs%form)),'directeffect_rs')
    end select
    
  end subroutine directeffect_rs


  subroutine directeffect_log_rs(rs,V,S,T,sz,i,j,D,dDdV,dDdS)
    ! DIRECTEFFECT_LOG_RS returns direct effect, D, of a rate-and-state law written in the form:
    ! dS/dt=D(V,S)*dlogV/dt+R(V,S)
    ! 
    ! Modified: 19 September 2007

    use constants, only : zero,one
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! RS = rate-and-state variables
    ! V = slip velocity
    ! S = strength
    ! T = temperature
    ! SZ = stress in z direction
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! D = direct effect (coefficient of dlogV/dt)
    ! DDDV = derivative of D(V,S) with respect to V
    ! DDDS = derivative of D(V,S) with respect to S

    type(ratestate_type),intent(in) :: rs
    real(pr),intent(in) :: V,S,T,sz
    integer(pin),intent(in) :: i,j
    real(pr),intent(out) :: D
    real(pr),intent(out),optional :: dDdV,dDdS

    ! Internal  Parameters:
    ! N = normal stress, positive in compression
    ! A = direct effect parameter
    ! AN = a*N

    real(pr) :: N,a,aN

    ! set magnitude of direct effect

    if (rs%a_of_T) then
       a = rs%a(i,j)*(T/293._pr)
    else
       a = rs%a(i,j)
    end if

    ! calculate normal stress

    N = max(-sz,zero)

    aN = a*N

    ! evaluate direct effect term

    select case(rs%form)
    case('SL','FL','SF')
       D = aN
       if (present(dDdV)) dDdV = zero
       if (present(dDdS)) dDdS = zero
    case('RFL','RSF')
       D = A*tanh(S/A)
       if (present(dDdV)) dDdV = zero
       if (present(dDdS)) dDdS = aN*(one-tanh(S/aN)**2)
    case default
       call error('Direct effect not yet defined for friction law: ' // &
            trim(adjustl(rs%form)),'directeffect_log_rs')
    end select
    
  end subroutine directeffect_log_rs


  subroutine rate_rs(rs,V,S,T,sz,i,j,R,dRdV,dRdS,err)
    ! RATE_RS returns rate, R, of a rate-and-state law written in the form: 
    ! dS/dt=D(V,S)*dV/dt+R(V,S)
    ! 
    ! Modified: 19 September 2007

    use constants, only : zero
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! RS = rate-and-state variables
    ! V = slip velocity
    ! S = strength
    ! T = temperature
    ! SZ = stress in z direction
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! R = rate (how strength evolves, excluding direct effect)
    ! DRDV = derivative of R(V,S) with respect to V
    ! DRDS = derivative of R(V,S) with respect to S
    ! ERR = error flag (true if friction isn't defined at particular V)

    type(ratestate_type),intent(in) :: rs
    real(pr),intent(in) :: V,S,T,sz
    integer(pin),intent(in) :: i,j
    real(pr),intent(out) :: R
    real(pr),intent(out),optional :: dRdV,dRdS
    logical,intent(out),optional :: err

    ! Internal  Parameters:
    ! F = steady state coefficient of friction
    ! DFDV = derivative of f with respect to V
    ! N = normal stress, positive in compression

    real(pr) :: f,dfdV,N

    ! calculate normal stress

    N = max(-sz,zero)

    ! calculate steady state coefficient of friction and its derivative with respect to V

    call fss_rs(rs,V,T,i,j,f,dfdV,err)

    if (err) return

    select case(rs%form)
    case('SL','FL','RFL','SF','RSF')
       R = V/rs%L(i,j)*(f*N-S)
       if (present(dRdV)) dRdV = R/V+V/rs%L(i,j)*dfdV*N
       if (present(dRdS)) dRdS = -V/rs%L(i,j)
    case default
       call error('Rate not yet defined for friction law: ' // &
            trim(adjustl(rs%form)),'rate_rs')
    end select
    
  end subroutine rate_rs


  subroutine fss_rs(rs,V,T,i,j,f,dfdV,err)
    ! FSS_RS returns the steady state friction coefficient
    ! 
    ! Modified: 19 September 2007
    
    use constants, only : zero,one,pi
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! RS = rate-and-state variables
    ! V = slip velocity
    ! T = temperature
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! F = steady state friction coefficient
    ! DFDV = derivative of f with respect to V
    ! ERR = error flag (true if strength isn't defined at particular V)

    real(pr),intent(in) :: V,T
    integer(pin),intent(in) :: i,j
    type(ratestate_type),intent(in) :: rs
    real(pr),intent(out) :: f
    real(pr),intent(out),optional :: dfdV
    logical,intent(out),optional :: err

    ! Internal  Parameters:
    ! FL = low-speed steady-state friction
    ! VW = weakening velocity
    ! WK = parameter combination in EW law
    ! DFLDV = derivative of fl with respect to V
    ! A = direct effect parameter

    real(pr) :: fl,Vw,wk,dfldV,a
    real(pr),parameter :: n=8._pr

    ! default is no error

    if (present(err)) err = .false.

    ! set magnitude of direct effect

    if (rs%a_of_T) then
       a = rs%a(i,j)*(T/293._pr)
    else
       a = rs%a(i,j)
    end if

    ! linearized law (easier to treat separately)

    if (rs%form=='LN') then
       f = rs%f0(i,j)+(a-rs%b(i,j))*V/rs%V0(i,j)
       if (present(dfdV)) &
            dfdV = (a-rs%b(i,j))/rs%V0(i,j)
       return
    end if

    ! check for possible problems

    if (V<=zero) then
       if (present(err)) err = .true.
       return
    end if

    ! evaluate friction coefficient
    
    fl = rs%f0(i,j)+(a-rs%b(i,j))*log(V/rs%V0(i,j))
    if (present(dfdV)) dfldV = (a-rs%b(i,j))/V

    select case(rs%form)
    case('AL','RAL','SL','RSL')
       f = fl
       if (present(dfdV)) dfdV = dfldV
    case('EW')
       wk = one+V/rs%Vw(i,j)
       f = fl/wk
       if (present(dfdV)) dfdV = dfldV/wk-fl/rs%Vw(i,j)/wk**2
    case('FL','RFL','SF','RSF')
       if (rs%form=='FL'.or.rs%form=='RFL') then
          Vw = pi*rs%Vth(i,j)*((rs%Tw(i,j)-T)/rs%Tc(i,j))**2
       else
          Vw = rs%Vw(i,j)
       end if
       f = rs%fw(i,j)+(fl-rs%fw(i,j))/(one+(V/Vw)**n)**(one/n)
       if (present(dfdV)) dfdV = (dfldV-(fl-rs%fw(i,j))/(V*(one+(Vw/V)**n)))/ &
            (one+(V/Vw)**n)**(one/n)
       !if (V<=Vw) then
       !   f = fl
       !   if (present(dfdV)) dfdV = dfldV
       !else
       !   f = rs%fw(i,j)+(fl-rs%fw(i,j))*Vw/V
       !   if (present(dfdV)) dfdV = (dfldV-(fl-rs%fw(i,j))/V)*Vw/V
       !end if
    case default
       call error('Steady state friction not defined for friction law: ' // &
            trim(adjustl(rs%form)),'fss_rs')
    end select

  end subroutine fss_rs


end module ratestate
