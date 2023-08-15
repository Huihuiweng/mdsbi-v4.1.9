! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module init

  ! INIT contains routines to initialize a problem
  ! 
  ! Modified: 19 July 2010

  use constants, only : pin

  implicit none

contains


  subroutine initialize(pb,ninput)
    ! INITIALIZE initializes a problem
    ! 
    ! Modified: 19 July 2010
    
    use io, only : new_io_unit,error
    use problem, only : problem_type
    use model, only : read_model,init_model
    use convolution, only : read_convolution,init_convolution
    use kernel, only : read_kernel,init_kernel
    use fields, only : read_fields,init_fields
    use fft_routines, only : init_fft
    use friction, only : read_friction,init_friction
    use mesh, only : read_mesh,init_mesh
    use front, only : read_front,init_front
    use mpi_routines, only : is_master
    use mpi

    implicit none
    
    ! I/O Parameters:
    ! PB = problem
    ! NINPUT = unit for input file

    type(problem_type),intent(inout) :: pb
    integer(pin),intent(inout) :: ninput

    ! Internal Parameters:
    ! NAME_INPUT = name of input file
    ! NAME_ECHO = name of file in which to echo all input variables
    ! (this is set up for use with matlab)
    ! NECHO = unit for echo (matlab) file
    ! STAT = I/O error flag
    ! IERR = MPI error flag

    character(64) :: name_echo
    integer(pin) :: necho,stat,ierr
    
    if (is_master) then

       ! open file unit for echoing input variables
       
       necho = new_io_unit()
       name_echo = 'data/' // trim(adjustl(pb%name)) // '.m'
       open(necho,file=name_echo,iostat=stat)
       if (stat/=0) call error("Error opening file: '" // &
            trim(adjustl(name_echo)) // "'",'initialize')
       
    else

       ! set default value indicating that unit is not connected to file

       necho = -1

    end if

    ! read in configuration variables set in *.in input file
       
    call read_model(ninput,pb%mdl)
    call read_kernel(ninput,pb%krn)
    call read_convolution(ninput,pb%cnv)
    call read_friction(ninput,pb%mdl,pb%fld%flt,pb%fri)
    call read_fields(ninput,pb%mdl,pb%fld)

    call read_mesh(ninput,pb%msh)
    call read_front(ninput,pb%frt)

    ! close input file

    close(ninput,iostat=stat)
    if (stat/=0) call error('Error closing input file','initialize')

    ! allocate and initialize main data types
    ! echo input variables to a single matlab file stored in data directory

    call init_model(necho,pb%mdl)
    call init_fft(pb%mdl,pb%cnv)
    call init_fields(necho,pb%mdl,pb%fld)
    call init_friction(necho,pb%mdl,pb%fri,pb%fld%flt)
    call init_kernel(necho,pb%mdl,pb%krn,pb%cnv%formulation)
    call init_convolution(necho,pb%mdl,pb%krn,pb%cnv)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! initialize output routines
       
    call init_mesh(necho,pb%name,pb%msh,pb%mdl,pb%fld)
    call init_front(necho,pb%name,pb%frt,pb%mdl,pb%fld)
       
    ! close echo file

    if (is_master) then
       
       close(necho,iostat=stat)
       if (stat/=0) call error("Error closing file: '" // &
            trim(adjustl(name_echo)) // "'",'initialize')

    end if

  end subroutine initialize
  

end module init
