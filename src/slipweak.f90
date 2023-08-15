! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module slipweak

  ! SLIPWEAK contains variables and routines related to the
  ! slip-weakening friction laws, in which shear strength depends
  ! only on slip
  ! 
  ! Modified: 20 July 2010

  use constants, only : pr,pin
  use asperity, only : asperity_list

  implicit none

  ! SLIPWEAK_TYPE is a derived type containing variables related to 
  ! slip-weakening friction law
  !
  ! Parameters:
  ! FORM = form of slip-weakening
  !      linear = linear slip-weakening
  !      exponential = exponential slip-weakening
  !      dugdale = Dugdale slip-weakening
  !      power = power-law slip-weakening
  !      U-shape = nonmonotonic slip-weakening
  !      kinematic = prescribed rupture speed
  ! ASPERITY_FILE = add heterogeneity to initial slip-weakening fields
  !      T = use asperity array data file
  !      F = do not use asperity array data file
  ! FILENAME = asperity file name
  ! LIST = singly-linked list of asperities
  ! P = power-law exponent
  ! VM = prescribed rupture velocity, left tip
  ! VP = prescribed rupture velocity, right tip
  ! DFDX = increase of coefficient of friction with distance
  ! XM = prescribed position of left tip
  ! XP = prescribed position of right tip
  ! DFDXM = increase of coefficient of friction with distance to left
  ! DFDXP = increase of coefficient of friction with distance to right
  ! RGLR = type of regularization
  !     N = none
  !     LT = length (slip) and time
  !     T = time only
  !     L = length (slip) only
  ! L = regularizing slip scale
  ! V = regularizing slip-rate scale
  ! T = regularizing time scale
  ! DC = slip-weakening distance at each point
  ! FD = dynamic coefficient of friction at each point
  ! FS = static coefficient of friction at each point

  type slipweak_type
     logical :: asperity_file
     character(64) :: filename,form,rglr
     type(asperity_list) :: list
     real(pr) :: p,L,V,T,vm,vp,dfdx,xm,xp,dfdxm,dfdxp
     real(pr),dimension(:,:),allocatable :: Dc,fd,fs
  end type slipweak_type

contains


  subroutine read_slipweak(ninput,sw)
    ! READ_SLIPWEAK reads in slipweak variables from file
    ! 
    ! Modified: 20 July 2010

    use constants, only : zero,one
    use io, only : error
    use asperity, only : read_asperity_list

    implicit none
    
    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! SW = slip-weakening variables

    integer(pin),intent(in) :: ninput
    type(slipweak_type),intent(out) :: sw

    ! Internal Parameters:
    ! STAT = I/O error flag
    ! FORM = form of slip-weakening
    ! ASPERITY_LIST = add heterogeneity to initial sw fields with list
    ! ASPERITY_FILE = add heterogeneity to initial sw fields from file
    ! FILENAME = asperity file name
    ! VM = prescribed rupture velocity, left tip
    ! VP = prescribed rupture velocity, right tip
    ! DFDX = increase of coefficient of friction with distance
    ! P = power-law exponent
    ! XM = prescribed position of left tip
    ! XP = prescribed position of right tip
    ! DFDXM = increase of coefficient of friction with distance to left
    ! DFDXP = increase of coefficient of friction with distance to right
    ! RGLR = type of regularization
    ! L = regularizing slip scale
    ! V = regularizing slip-rate scale
    ! T = regularizing time scale
    ! NDATA = number of field perturbations

    integer(pin) :: stat,ndata
    logical :: asperity_list,asperity_file
    character(64) :: filename,form,rglr
    real(pr) :: p,L,V,T,vm,vp,dfdx,xm,xp,dfdxm,dfdxp

    ! make namelist of user input variables

    namelist /slipweak_list/ form,asperity_list,asperity_file,filename,p,rglr,L,V,T, &
         vm,vp,dfdx,xm,xp,dfdxm,dfdxp

    ! default values
    
    form = 'linear'
    vm = zero
    vp = zero
    dfdx = one
    xm  = -one
    xp = one
    dfdxm = one
    dfdxp = one
    p = one
    rglr = 'N'
    L = zero
    V = zero
    T = zero
    asperity_list = .false.
    asperity_file = .false.
    
    ! read namelist from input file, call error routine if needed
    rewind(ninput)
    read(ninput,nml=slipweak_list,iostat=stat)
    if (stat/=0) call error("Error reading namelist 'slipweak_list' in .in file",'read_slipweak')
    
    ! assign input variables to components of derived type
    
    sw%form = form
    sw%asperity_file = asperity_file
    sw%filename = filename
    sw%vm = vm
    sw%vp = vp
    sw%dfdx = dfdx
    sw%rglr = rglr
    sw%p = p
    sw%L = L
    sw%V = V
    sw%T = T
    sw%xm = xm
    sw%xp = xp
    sw%dfdxm = dfdxm
    sw%dfdxp = dfdxp
    
    ! asperity list, set ndata = 3 (for Dc,fd,fs)
    
    ndata = 3
    if (asperity_list) &
         sw%list = read_asperity_list(ninput,'slipweak',ndata)
    
    ! set rglr='N' if appropriate
    
    select case(sw%rglr)
    case('LT','L')
       if (sw%L==zero) sw%rglr = 'N'
    case('T')
       if (sw%T==zero) sw%rglr = 'N' 
    end select

  end subroutine read_slipweak


  subroutine init_slipweak(necho,mdl,sw)
    ! INIT_SLIPWEAK initializes slipweak variables
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
    ! SW = slip-weakening variables

    integer(pin),intent(in) :: necho
    type(model_type),intent(in) :: mdl
    type(slipweak_type),intent(inout) :: sw

    ! Internal Parameters:
    ! DARRAY = distributed array type
    ! FH = file handle

    integer(pin) :: darray
    type(file_distributed) :: fh

    ! allocate arrays

    allocate(sw%Dc(mdl%mx:mdl%px,mdl%ny))
    allocate(sw%fd(mdl%mx:mdl%px,mdl%ny))
    allocate(sw%fs(mdl%mx:mdl%px,mdl%ny))

    sw%Dc = zero
    sw%fd = zero
    sw%fs = zero

    ! read asperity arrays from file
    
    if (sw%asperity_file) then
       call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
       call open_file_distributed(fh,sw%filename,'read',MPI_COMM_WORLD,darray,pr)
       call read_file_distributed(fh,sw%Dc)
       call read_file_distributed(fh,sw%fd)
       call read_file_distributed(fh,sw%fs)
       call close_file_distributed(fh)
    end if
       
    ! use list data to perturb fields
    
    call assign_list_data(sw%Dc,sw%list,1,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(sw%fd,sw%list,2,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(sw%fs,sw%list,3,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    
    ! deallocate memory for asperity list
    
    call destroy_asperity_list(sw%list)

    ! output variable values into matlab file

    if (is_master) then
       
       call write_matlab(necho,'form',sw%form,'sw')
       
       select case(sw%form)
       case('power','powerlaw')
          call write_matlab(necho,'p',sw%p,'sw')
       case('linear')
       case('time-weakening')
       case('exponential')
       case('dugdale','Dugdale')
       case('kinematic')
          call write_matlab(necho,'vm',sw%vm,'sw')
          call write_matlab(necho,'vp',sw%vp,'sw')
          call write_matlab(necho,'dfdx',sw%dfdx,'sw')
       case('static')
          call write_matlab(necho,'xm',sw%xm,'sw')
          call write_matlab(necho,'xp',sw%xp,'sw')
          call write_matlab(necho,'dfdxm',sw%dfdxm,'sw')
          call write_matlab(necho,'dfdxp',sw%dfdxp,'sw')
       end select
       
       call write_matlab(necho,'rglr',sw%rglr,'sw')
       
       select case(sw%rglr)
       case('LT')
          call write_matlab(necho,'L',sw%L,'sw')
          call write_matlab(necho,'V',sw%V,'sw')
       case('L')
          call write_matlab(necho,'L',sw%L,'sw')
       case('T')
          call write_matlab(necho,'T',sw%T,'sw')
       end select

    end if

  end subroutine init_slipweak


  subroutine destroy_slipweak(sw)
    ! DESTROY_SLIPWEAK destroys derived type sw
    ! 
    ! Modified: 26 November 2005

    implicit none

    ! I/O Parameters:
    ! SW = slip-weakening variables

    type(slipweak_type),intent(inout) :: sw

    ! deallocate memory assigned to allocatable arrays

    if (allocated(sw%Dc)) deallocate(sw%Dc)
    if (allocated(sw%fd)) deallocate(sw%fd)
    if (allocated(sw%fs)) deallocate(sw%fs)

  end subroutine destroy_slipweak


  function fric_sw(U,i,j,r,t,sw) result(f)
    ! FRIC_SW returns the coefficient of friction for a slip-dependent friction law, given the slip
    ! 
    ! Modified: 23 July 2008

    use constants, only : zero,half,one,two,four
        
    implicit none

    ! I/O Parameters:
    ! U = cumulative slip (line integral)
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! R = ridus 
    ! T = time
    ! SW = slip-weakening variables
    ! F = coefficient of friction

    real(pr),intent(in) :: U,r,t
    integer(pin),intent(in) :: i,j
    type(slipweak_type),intent(in) :: sw
    real(pr) :: f

    ! Internal  Parameters:
    ! XM = position of left crack-tip
    ! XP = position of right crack-tip
    ! FK = friction coefficient from kinematic (forced) rupture
    ! FKM = friction coefficient from kinematic (forced) rupture, minus branch
    ! FKP = friction coefficient from kinematic (forced) rupture, plus branch
    ! FSW = friction coefficient from slip-weakening

    real(pr) :: xm,xp,fk,fkm,fkp,fsw

    ! evaluate coefficient of friction

    if (sw%Dc(i,j)==zero) then
       f = sw%fd(i,j)
       return
    end if

    select case(sw%form)

    case('power','powerlaw') ! power-law weakening

       f = max(sw%fs(i,j)-(sw%fs(i,j)-sw%fd(i,j))*(sw%p/(sw%p+one)*U/sw%Dc(i,j))**sw%p,sw%fd(i,j))

    case('linear') ! linear weakening

       f = max(sw%fs(i,j)-(sw%fs(i,j)-sw%fd(i,j))*U/sw%Dc(i,j),sw%fd(i,j))

    case('time-weakening') ! time weakening

       fk = max(sw%fs(i,j)-(sw%fs(i,j)-sw%fd(i,j))*U/sw%Dc(i,j),sw%fd(i,j))
        
       ! Prescribed rupture tip location
       xm = sw%V*t
       if(abs(r)<sw%V*sw%T) then
             if (abs(r)<=xm-sw%L) then
                fkm = sw%fd(i,j)
             elseif((abs(r)<xm).and.(abs(r)>xm-sw%L)) then
                fkm = sw%fs(i,j) - (sw%fs(i,j)-sw%fd(i,j))*(xm-abs(r))/sw%L
             else
                fkm = fk
             end if
       else
             fkm=fk
       end if
       ! combine by taking smaller of two friction coefficients
       f = min(fk,fkm)

    case('exponential') ! exponential weakening

       f = sw%fd(i,j)+(sw%fs(i,j)-sw%fd(i,j))*exp(-U/sw%Dc(i,j))

    case('dugdale','Dugdale') ! Dugdale model

       if (U>=sw%Dc(i,j)) then
          f = sw%fd(i,j)
       elseif (U<sw%Dc(i,j)) then
          f = sw%fs(i,j)
       else ! U=Dc
          f = half*(sw%fd(i,j)+sw%fs(i,j))
       end if

    case('U-shape') ! U-shape (weakening and strengthing)

       f = sw%fd(i,j)+(sw%fs(i,j)-sw%fd(i,j))* &
            exp(two*((two*U/sw%Dc(i,j)-one)**four-(one+U/sw%Dc(i,j))))

    case('kinematic') ! prescribed rupture speed

       ! left tip

       xm = sw%vm*t

       if (r<xm) then
          fkm = sw%fd(i,j)+sw%dfdx*(xm-r)
       else
          fkm = sw%fd(i,j)
       end if

       ! right tip

       xp = sw%vp*t

       if (r>xp) then
          fkp = sw%fd(i,j)+sw%dfdx*(r-xp)
       else
          fkp = sw%fd(i,j)
       end if

       ! combine

       fk = max(fkm,fkp)

       ! slip-weakening friction coefficient

       fsw = max(sw%fs(i,j)-(sw%fs(i,j)-sw%fd(i,j))*U/sw%Dc(i,j),sw%fd(i,j))       

       ! combine by taking smaller of two friction coefficients

       f = min(fk,fsw)

    case('static') ! static crack

       xm = sw%xm
       xp = sw%xp

       if (r<xm) then ! left of crack
          f = sw%fd(i,j)+sw%dfdxm*(xm-r)
       elseif (r>xp) then ! right of crack
          f = sw%fd(i,j)+sw%dfdxp*(r-xp)
       else
          f = sw%fd(i,j)
       end if

    end select

  end function fric_sw


  function strength_sw(U,Q,sz,i,j,x,t,sw) result(st)
    ! STRENGTH_SW returns the fault strength
    ! 
    ! Modified: 22 August 2006

    use constants, only : zero

    implicit none

    ! I/O Parameters:
    ! U = slip
    ! Q = state variable (strength)
    ! SZ = stress in z direction
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! X = distance in x direction
    ! T = time
    ! SW = slip-weakening variables
    ! ST = strength

    real(pr),intent(in) :: U,Q,sz,x,t
    integer(pin),intent(in) :: i,j
    type(slipweak_type),intent(in) :: sw
    real(pr) :: st

    ! Internal  Parameters:
    ! F = coefficient of friction
    ! SN = normal stress, positive in compression

    real(pr) :: f,sn

    if (sw%rglr=='N') then
       sn = max(-sz,zero)
       f = fric_sw(U,i,j,x,t,sw)
       st = f*sn
    else
       st = Q
    end if

  end function strength_sw


  function set_strength_sw(U,sz,i,j,x,t,sw) result(st)
    ! SET_STRENGTH_SW returns the fault strength for use as initial value of 'state' variable
    ! 
    ! Modified: 14 August 2007

    use constants, only : zero

    implicit none

    ! I/O Parameters:
    ! U = slip
    ! SZ = stress in z direction
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! X = distance in x direction
    ! T = time
    ! SW = slip-weakening variables
    ! ST = strength

    real(pr),intent(in) :: U,sz,x,t
    integer(pin),intent(in) :: i,j
    type(slipweak_type),intent(in) :: sw
    real(pr) :: st

    ! Internal  Parameters:
    ! SN = normal stress, positive in compression
    ! F = coefficient of friction

    real(pr) :: sn,f

    ! calculate normal stress

    sn = max(-sz,zero)

    ! calculate coefficient of friction

    f = fric_sw(U,i,j,x,t,sw)

    ! set strength to product of coefficient of friction and normal stress

    st = f*sn

  end function set_strength_sw


  function evolution_sw(U,V,Q,sz,i,j,x,t,sw) result(dQ)
    ! EVOLUTION_SW returns the state-rate, where 'state' is shear strength
    ! 
    ! Modified: 22 August 2006

    use constants, only : zero

    implicit none

    ! I/O Parameters:
    ! U = slip
    ! V = slip velocity
    ! Q = state variable
    ! SZ = stress in z direction
    ! I = index in x direction at which friction is applied
    ! J = index in y direction at which friction is applied  
    ! X = distance in x direction
    ! T = time
    ! SW = slip-weakening variables
    ! DQ = state-rate

    real(pr),intent(in) :: U,V,Q,sz,x,t
    integer(pin),intent(in) :: i,j
    type(slipweak_type),intent(in) :: sw
    real(pr) :: dQ

    ! Internal  Parameters:
    ! SN = normal stress, positive in compression
    ! F = coefficient of friction
    ! ST = shear strength

    real(pr) :: sn,f,st

    ! calculate strength

    sn = max(-sz,zero)
    f = fric_sw(U,i,j,x,t,sw)
    st = f*sn

    ! calculate state-rate

    select case(sw%rglr)
    case('N')
       dQ = zero
    case('LT')
       dQ = -(Q-st)*(abs(V)+sw%V)/sw%L
    case('L')
       dQ = -(Q-st)*abs(V)/sw%L
    case('T')
       dQ = -(Q-st)/sw%T
    end select

  end function evolution_sw


end module slipweak
