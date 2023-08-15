! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module load

  ! LOAD contains routines to manipulate time-dependent loads
  ! 
  ! Modified: 13 July 2012

  use constants, only : pin,pr
  use asperity, only : asperity_list

  implicit none

  ! LOAD_TYPE is a derived type containing load variables
  !
  ! Parameters:
  ! ASPERITY_FILE = add heterogeneity to initial slip-weakening fields
  !      T = use asperity array data file
  !      F = do not use asperity array data file
  ! FILENAME = asperity file name
  ! FORM = form of time-dependent load
  ! T = loading time parameter (time beyond which load is held fixed, period of loading, etc.)
  ! SX0 = load in x direction
  ! SY0 = load in y direction
  ! SZ0 = load in z direction
  ! AX = loading parameter (loading rate, perturbation amplitude, etc.) in x direction 
  ! AY = loading parameter in y direction
  ! AZ = loading parameter in z direction
  ! WAVE_LOAD = dynamic load from passage of seismic wave
  ! KX = wavenumber in x direction (wave load)
  ! KY = wavenumber in y direction (wave load)
  ! OMEGA = angular frequency (wave load)
  ! WX = amplitude in x direction (wave load)
  ! WY = amplitude in y direction (wave load)
  ! WZ = amplitude in z direction (wave load)

  type load_type
     logical :: asperity_file,wave_load
     character(64) :: filename,form
     type(asperity_list) :: list
     real(pr) :: T,kx,ky,omega,Wx,Wy,Wz
     real(pr),dimension(:,:),allocatable :: sx0,sy0,sz0,Ax,Ay,Az
  end type load_type

contains


  subroutine read_load(ninput,ld)
    ! READ_LOAD reads in load variables from file
    ! 
    ! Modified: 13 July 2012

    use constants, only : zero,one
    use io, only : error
    use asperity, only : read_asperity_list

    implicit none

    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! LD = load variables

    integer(pin),intent(in) :: ninput
    type(load_type),intent(out) :: ld

    ! Internal Parameters:
    ! STAT = I/O error flag
    ! ASPERITY_LIST = add heterogeneity to initial ld fields with list
    ! ASPERITY_FILE = add heterogeneity to initial ld fields from file
    ! FILENAME = asperity file name
    ! NDATA = number of field perturbations
    ! T = loading time parameter
    ! FORM = form of time-dependent load
    ! WAVE_LOAD = dynamic load from passage of seismic wave
    ! KX = wavenumber in x direction (wave load)
    ! KY = wavenumber in y direction (wave load)
    ! OMEGA = angular frequency (wave load)
    ! WX = amplitude in x direction (wave load)
    ! WY = amplitude in y direction (wave load)
    ! WZ = amplitude in z direction (wave load)
    
    integer(pin) :: stat,ndata
    logical :: asperity_list,asperity_file,wave_load
    character(64) :: filename,form
    real(pr) :: T,kx,ky,omega,Wx,Wy,Wz

    ! make namelist of user input variables

    namelist /load_list/ asperity_list,asperity_file,filename,T,form, &
         wave_load,kx,ky,omega,Wx,Wy,Wz

    ! default values
    
    asperity_list = .false.
    asperity_file = .false.
    T = zero
    form = 'step'
    wave_load = .false.
    kx = zero
    ky = zero
    omega = zero
    Wx = zero
    Wy = zero
    Wz = zero

    ! read namelist from input file, call error routine if needed
    
    rewind(ninput)
    read(ninput,nml=load_list,iostat=stat)
    if (stat/=0) call error("Error reading namelist 'load_list' in .in file",'read_load')
    
    ! assign input variables to components of derived type
    
    ld%asperity_file = asperity_file
    ld%filename = filename
    ld%T = T
    ld%form = form
    ld%wave_load = wave_load
    ld%kx = kx
    ld%ky = ky
    ld%omega = omega
    ld%Wx = Wx
    ld%Wy = Wy
    ld%Wz = Wz

    ! asperity list, set ndata = 3 (for Ax,Ay,Az)
    
    ndata = 3
    if (asperity_list) &
         ld%list = read_asperity_list(ninput,'load',ndata)

  end subroutine read_load


  subroutine init_load(necho,mdl,ld)
    ! INIT_LOAD initializes load variables
    ! 
    ! Modified: 13 July 2012

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
    ! LD = load variables

    integer(pin),intent(in) :: necho
    type(model_type),intent(in) :: mdl
    type(load_type),intent(inout) :: ld

    ! Internal Parameters:
    ! DARRAY = distributed array type
    ! FH = file handle
    
    integer(pin) :: darray
    type(file_distributed) :: fh

    ! allocate arrays

    allocate(ld%sx0(mdl%mx:mdl%px,mdl%ny))
    allocate(ld%sy0(mdl%mx:mdl%px,mdl%ny))
    allocate(ld%sz0(mdl%mx:mdl%px,mdl%ny))
    allocate(ld%Ax(mdl%mx:mdl%px,mdl%ny))
    allocate(ld%Ay(mdl%mx:mdl%px,mdl%ny))
    allocate(ld%Az(mdl%mx:mdl%px,mdl%ny))

    ld%sx0 = zero
    ld%sy0 = zero
    ld%sz0 = zero
    ld%Ax = zero
    ld%Ay = zero
    ld%Az = zero

    ! read asperity arrays from file
    
    if (ld%asperity_file) then
       call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
       call open_file_distributed(fh,ld%filename,'read',MPI_COMM_WORLD,darray,pr)
       call read_file_distributed(fh,ld%Ax)
       call read_file_distributed(fh,ld%Ay)
       call read_file_distributed(fh,ld%Az)
       call close_file_distributed(fh)
    end if
       
    ! use list data to perturb fields
    
    call assign_list_data(ld%Ax,ld%list,1,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(ld%Ay,ld%list,2,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(ld%Az,ld%list,3,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    
    ! deallocate memory for asperity list
    
    call destroy_asperity_list(ld%list)

    ! output variable values into matlab file

    if (is_master) then
       call write_matlab(necho,'form',ld%form,'ld')
       call write_matlab(necho,'T',ld%T,'ld')
       if (ld%wave_load) then
          call write_matlab(necho,'kx',ld%kx,'ld')
          call write_matlab(necho,'ky',ld%ky,'ld')
          call write_matlab(necho,'omega',ld%omega,'ld')
          call write_matlab(necho,'Wx',ld%Wx,'ld')
          call write_matlab(necho,'Wy',ld%Wy,'ld')
          call write_matlab(necho,'Wz',ld%Wz,'ld')
       end if
    end if

  end subroutine init_load


  subroutine destroy_load(ld)
    ! DESTROY_LOAD destroys derived type ld
    ! 
    ! Modified: 20 April 2007

    implicit none

    ! I/O Parameters:
    ! LD = load variables

    type(load_type),intent(inout) :: ld

    ! deallocate memory assigned to allocatable arrays

    if (allocated(ld%sx0)) deallocate(ld%sx0)
    if (allocated(ld%sy0)) deallocate(ld%sy0)
    if (allocated(ld%sz0)) deallocate(ld%sz0)
    if (allocated(ld%Ax)) deallocate(ld%Ax)
    if (allocated(ld%Ay)) deallocate(ld%Ay)
    if (allocated(ld%Az)) deallocate(ld%Az)

  end subroutine destroy_load


  subroutine set_load(mdl,ld,flt)
    ! SET_LOAD sets values of the load
    ! 
    ! Modified: 13 July 2012

    use constants, only : zero,half,one,two,pi,twopi
    use model, only : model_type
    use fault, only : fault_type
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! LD = load variables
    ! FLT = fault variables

    type(model_type),intent(in) :: mdl
    type(load_type),intent(inout) :: ld
    type(fault_type),intent(inout) :: flt

    ! Internal Parameters:
    ! F = fraction of load to be applied

    real(pr) :: f

    ! return if not including time-dependent load

    if (.not.flt%dynload) return

    ! subtract old time-dependent load from load
    
    flt%sx0 = flt%sx0-ld%sx0
    flt%sy0 = flt%sy0-ld%sy0
    flt%sz0 = flt%sz0-ld%sz0

    ! calculate current time-dependent load

    select case(ld%form)
    case default
       call error('Invalid form of time-dependent load','set_load')
    case('step') ! discontinuous
       if (mdl%t<ld%T) then
          f = zero
       else
          f = one
       end if
    case('ramp') ! discontinuous first derivative (C_1)
       f = max(min(mdl%t,ld%T),zero)/ld%T
    case('smooth1')
       if (mdl%t<=zero) then
          f = zero
       elseif (mdl%t<ld%T) then
          f = half*(one-cos(pi*mdl%t/ld%T))
       else
          f = one
       end if
    case('smooth2') ! discontinuous second derivative (C_2)
       if (mdl%t<=zero) then
          f = zero
       elseif (mdl%t<half*ld%T) then
          f = mdl%t/ld%T-sin(twopi*mdl%t/ld%T)/twopi
       elseif (mdl%t<ld%T) then
          f = mdl%t/ld%T+sin(twopi*(ld%T-mdl%t)/ld%T)/twopi
       else
          f = one
       end if
    case('smoothI') ! all derivative continuous (C_infinity)
       if (mdl%t<=zero) then
          f = zero
       elseif (mdl%t<ld%T) then
          f = half*(one+tanh(-ld%T/mdl%t+ld%T/(ld%T-mdl%t)))
       else
          f = one
       end if
    case('smooth') ! all derivative continuous (C_infinity)
       if (mdl%t<=zero) then
          f = zero
       elseif (mdl%t<ld%T) then
          f = exp((mdl%t-ld%T)**2/(mdl%t*(mdl%t-two*ld%T)))
       else
          f = one
       end if
    case('sine') ! sine with period T
       f = sin(twopi*mdl%t/ld%T)
    end select

    ld%sx0 = ld%Ax*f
    ld%sy0 = ld%Ay*f
    ld%sz0 = ld%Az*f

    ! wave load

    if (ld%wave_load) call set_wave_load(mdl,ld)

    ! add time-dependent load to load

    flt%sx0 = flt%sx0+ld%sx0
    flt%sy0 = flt%sy0+ld%sy0
    flt%sz0 = flt%sz0+ld%sz0

  end subroutine set_load


  subroutine set_wave_load(mdl,ld)
    ! SET_WAVE_LOAD adds perturbation to load from wave
    ! 
    ! Modified: 13 July 2012

    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! LD = load variables

    type(model_type),intent(in) :: mdl
    type(load_type),intent(inout) :: ld

    ! Internal Parameters:
    ! I = index in x direction
    ! J = index in y direction
    ! PHASE_FACTOR = phase factor of sinusoidal wave

    integer(pin) :: i,j
    real :: phase_factor

    do j = 1,mdl%ny
       do i = mdl%mx,mdl%px
          phase_factor = sin(ld%kx*mdl%x(i)+ld%ky*mdl%y(j)-ld%omega*mdl%t)
          ld%sx0(i,j) = ld%sx0(i,j)+ld%Wx*phase_factor
          ld%sy0(i,j) = ld%sy0(i,j)+ld%Wy*phase_factor
          ld%sz0(i,j) = ld%sz0(i,j)+ld%Wz*phase_factor
       end do
    end do

  end subroutine set_wave_load


end module load
