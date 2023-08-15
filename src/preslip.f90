! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module preslip

  ! PRESLIP contains parameters and routines related to pre-existing crack solutions
  ! that can be used to nucleate ruptures from pre-existing cracks.
  ! 
  ! Modified: 20 July 2010

  use constants, only : pin,pr
  use asperity, only : asperity_list

  implicit none

  ! PRESLIP_TYPE is a derived type containing pre-existing crack variables
  !
  ! Parameters:
  ! FORM = friction law model
  !      dugdale = Dugdale cohesive zone model
  !      linear = linear slip-weakening cohesive zone model
  ! C = crack length, including cohesive zone
  ! A = crack length, excluding cohesive zone (a<c)
  ! SP = peak strength
  ! SR = residual strength
  ! ASPERITY_FILE = add heterogeneity to initial slip-weakening fields
  !      T = use asperity array data file
  !      F = do not use asperity array data file
  ! FILENAME = asperity file name
  ! LIST = singly-linked list of asperities

  type preslip_type
     logical :: asperity_file
     character(64) :: filename,form
     real(pr) :: c,a,sp,sr
     type(asperity_list) :: list
  end type preslip_type


contains


  subroutine read_preslip(ninput,mdl,ps)
    ! READ_PRESLIP reads in pre-existing crack variables from file
    ! 
    ! Modified: 20 July 2010

    use constants, only : zero,one,two
    use io, only : error
    use model, only : model_type
    use asperity, only : read_asperity_list

    implicit none

    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! MDL = model variables
    ! PS = pre-existing crack variables

    integer(pin),intent(in) :: ninput
    type(model_type),intent(in) :: mdl
    type(preslip_type),intent(out) :: ps

    ! Internal Parameters:
    ! FORM = friction law model
    ! ASPERITY_LIST = add heterogeneity to initial slip with list
    ! ASPERITY_FILE = add heterogeneity to initial slip from file
    ! FILENAME = asperity file name
    ! C = crack length, including cohesive zone
    ! A = crack length, excluding cohesive zone
    ! TP = peak strength
    ! TR = residual strength
    ! STAT = I/O error flag
    ! NDATA = number of field perturbations

    integer(pin) :: stat,ndata
    logical :: asperity_list,asperity_file
    character(64) :: filename,form
    real(pr) :: c,a,sp,sr

    ! make namelist of user input variables

    namelist /preslip_list/ form,asperity_list,asperity_file,filename,c,a,sp,sr

    ! default values
    
    form = 'none'
    c = one
    a = two
    sp = one
    sr = zero
    asperity_list = .false.
    asperity_file = .false.
    
    ! read namelist from input file, call error routine if needed
    
    rewind(ninput)
    read(ninput,nml=preslip_list,iostat=stat)
    if (stat/=0) call error("Error reading namelist 'preslip_list' in .in file",'read_preslip')
    
    ! assign input variables to components of derived type
    
    ps%form = form
    ps%c = c
    ps%a = a
    ps%sp = sp
    ps%sr = sr
    ps%filename = filename
    ps%asperity_file = asperity_file
    
    ! asperity list, set ndata
    
    if (mdl%bm) then
       ndata = 6 ! for uxp,uxm,uyp,uym,uzp,uzm
    else
       ndata = 2 ! for Ux,Uy 
    end if
    if (asperity_list) &
         ps%list = read_asperity_list(ninput,'preslip',ndata)

  end subroutine read_preslip


  subroutine init_preslip_im(mdl,ps,Ux,Uy)
    ! INIT_PRESLIP_IM initializes preslip variables, identical materials case
    ! 
    ! Modified: 20 July 2010

    use model, only : model_type
    use asperity, only : assign_list_data,destroy_asperity_list
    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use mpi_routines, only : MPI_REAL_PR,subarray
    use mpi

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! PS = preslip variables
    ! UX = slip in x direction
    ! UY = slip in y direction

    type(model_type),intent(in) :: mdl
    type(preslip_type),intent(inout) :: ps
    real(pr),dimension(:,:),intent(out) :: Ux,Uy

    ! Internal Parameters
    ! I = index in x direction
    ! J = index in y direction
    ! DARRAY = distributed array type
    ! FH = file handle

    integer(pin) :: i,j,darray
    type(file_distributed) :: fh

    ! read asperity arrays from file
    
    if (ps%asperity_file) then
       call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
       call open_file_distributed(fh,ps%filename,'read',MPI_COMM_WORLD,darray,pr)
       call read_file_distributed(fh,Ux)
       call read_file_distributed(fh,Uy)
       call close_file_distributed(fh)
    end if
       
    ! use list data to perturb fields
    
    call assign_list_data(Ux,ps%list,1,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(Uy,ps%list,2,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    
    ! deallocate memory for asperity list
       
    call destroy_asperity_list(ps%list)
    
    ! add (2D) analytical solution if required
    
    if (ps%form/='none') then
       call param_preslip(mdl,ps)
       do j = 1,mdl%ny
          do i = mdl%mx,mdl%px
             Ux(i,j) = Ux(i,j)+eval_preslip(mdl%x(i),mdl,ps)
          end do
       end do
    end if

  end subroutine init_preslip_im


  subroutine init_preslip_bm(mdl,ps,uxp,uxm,uyp,uym,uzp,uzm)
    ! INIT_PRESLIP_BM initializes preslip variables, bimaterial case
    ! 
    ! Modified: 20 July 2010

    use model, only : model_type
    use asperity, only : assign_list_data,destroy_asperity_list
    use io, only :file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use mpi_routines, only : MPI_REAL_PR,subarray
    use mpi 

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! PS = preslip variables
    ! UXP = displacement in x direction, plus side
    ! UXM = displacement in x direction, minus side
    ! UYP = displacement in y direction, plus side
    ! UYM = displacement in y direction, minus side
    ! UZP = displacement in z direction, plus side
    ! UZM = displacement in z direction, minus side

    type(model_type),intent(in) :: mdl
    type(preslip_type),intent(inout) :: ps
    real(pr),dimension(:,:),intent(out) :: uxp,uxm,uyp,uym,uzp,uzm

    ! Internal Parameters:
    ! DARRAY = distributed array type
    ! FH = file handle

    integer(pin) :: darray
    type(file_distributed) :: fh

    ! read asperity arrays from file
    
    if (ps%asperity_file) then
       call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PR,darray)
       call open_file_distributed(fh,ps%filename,'read',MPI_COMM_WORLD,darray,pr)
       call read_file_distributed(fh,uxp)
       call read_file_distributed(fh,uxm)
       call read_file_distributed(fh,uyp)
       call read_file_distributed(fh,uym)
       call read_file_distributed(fh,uzp)
       call read_file_distributed(fh,uzm)
       call close_file_distributed(fh)
    end if
       
    ! use list data to perturb fields
    
    call assign_list_data(uxp,ps%list,1,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(uxm,ps%list,2,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(uyp,ps%list,3,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(uym,ps%list,4,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(uzp,ps%list,5,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    call assign_list_data(uzm,ps%list,6,mdl%x(mdl%mx:mdl%px),mdl%y,mdl%dx,mdl%dy,mdl%dim)
    
    ! deallocate memory for asperity list
    
    call destroy_asperity_list(ps%list)

  end subroutine init_preslip_bm


  function eval_preslip(x,mdl,ps) result(U)
    ! EVAL_PRESLIP calculates slip for a pre-existing crack solution
    ! 
    ! Modified: 24 October 2006

    use constants, only : zero,half,one,two,pi
    use utilities, only : arccosh
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! X = distance along fault
    ! U = slip
    ! MDL = model variables
    ! PS = pre-existing crack variables

    real(pr),intent(in) :: x
    type(model_type),intent(in) :: mdl
    type(preslip_type),intent(in) :: ps
    real(pr) :: U

    ! Internal Parameters:
    ! MUE = effective shear modulus (depends on loading mode)

    real(pr) :: mue

    select case(mdl%mode)
    case default
       mue = mdl%mu/(one-mdl%nu)
    case(3)
       mue = mdl%mu
    end select

    if (abs(x)<ps%c) then
       select case(ps%form)
       case('linear')
          U = two*(ps%sp-ps%sr)/mue/(pi*(ps%c-ps%a))*( &
               sqrt(ps%c**2-x**2)*sqrt(ps%c**2-ps%a**2)-half*(x**2+ps%a**2)*log(abs((sqrt(ps%c**2-x**2)+ &
               sqrt(ps%c**2-ps%a**2))/(sqrt(ps%c**2-x**2)-sqrt(ps%c**2-ps%a**2))))+ &
               x*ps%a*log(abs((ps%a*sqrt(ps%c**2-x**2)+x*sqrt(ps%c**2-ps%a**2))/(ps%a*sqrt(ps%c**2-x**2)- &
               x*sqrt(ps%c**2-ps%a**2)))))
       case('dugdale','Dugdale')
          U = two*(ps%sp-ps%sr)/mue/pi*( &
               (x+ps%a)*arccosh(abs((ps%c**2+ps%a*x)/(ps%c*(ps%a+x))))- &
               (x-ps%a)*arccosh(abs((ps%c**2-ps%a*x)/(ps%c*(ps%a-x)))) )
       end select
    else
       U = zero
    end if

  end function eval_preslip


  subroutine param_preslip(mdl,ps)
    ! PARAM_PRESLIP calculates parameters for a pre-existing crack solution
    ! 
    ! Modified: 27 August 2006

    use constants, only : one,two,pi
    use io, only : message
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! PS = pre-existing crack variables

    type(model_type),intent(in) :: mdl
    type(preslip_type),intent(in) :: ps

    ! Internal Parameters:
    ! MUE = effective shear modulus (depends on loading mode)
    ! E = ratio of a to c
    ! G = stress ratio
    ! S = another stress ratio
    ! ASE = arcsin(e)
    ! S0 = initial stress for equilibrium
    ! D0 = weakening distance (slip at edge of cohesive zone)
    ! STR = character string

    real(pr) :: mue,e,g,S,ase,s0,d0
    character(64) :: str

    select case(mdl%mode)
    case default
       mue = mdl%mu/(one-mdl%nu)
    case(3)
       mue = mdl%mu
    end select

    e = ps%a/ps%c

    call message('pre-existing crack parameters:')

    select case(ps%form)
    case('linear')
       g = two*(sqrt(one-e**2)-e*acos(e))/(pi*(one-e))
       s0 = ps%sr+(ps%sp-ps%sr)*g
       d0 = two*(ps%sp-ps%sr)/mue/pi*ps%c*(one-e**2+two*e**2*log(e))/(one-e)
    case('dugdale','Dugdale')
       ase = asin(e)
       S = ase/(pi/two-ase)
       write(str,'(a,f14.7)') 'S = ',S
       call message(str)
       s0 = (ps%sp+S*ps%sr)/(one+S)
       d0 = 4._pr*(ps%sp-ps%sr)/mue/pi*ps%a*log(one/e)
    end select

    write(str,'(a,f14.7)') 's0 = ',s0
    call message(str)
    write(str,'(a,f14.7)') 'd0 = ',d0
    call message(str)

  end subroutine param_preslip


end module preslip
