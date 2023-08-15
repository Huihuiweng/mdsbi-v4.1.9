! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module front

  ! FRONT contains variables and routines for outputting rupture front history
  ! 
  ! Modified: 6 August 2010

  use io, only : file_distributed
  use constants, only : pr,pin

  implicit none

  ! FRONT_TYPE is a derived type containing rupture front variables
  !
  ! Parameters:
  ! FRONT = output rupture front history
  !      T = output rupture front history
  !      F = do not output rupture front history
  ! NAME = name of data file
  ! FIELD = field component to be contoured
  ! VAL = minimum value field must exceed (defines rupture front)
  ! XMIN = minimum value of x
  ! XMAX = maximum value of x
  ! YMIN = minimum value of y
  ! YMAX = maximum value of y
  ! OTS_XMIN = whether or not to include the grid point just outside of bound set by xmin
  ! OTS_XMAX = same as above for xmax
  ! OTS_YMIN = same as above for ymin
  ! OTS_YMAX = same as above for ymax
  ! FH = file handle
  ! MX = index corresponding to xmin
  ! PX = index corresponding to xmax
  ! MY = index corresponding to ymin
  ! PY = index corresponding to ymax
  ! NX = number of grid points in x direction
  ! NY = number of grid points in y direction
  ! DX = spacing in x direction
  ! DY = spacing in y direction
  ! STRIDE_X = stride in x direction
  ! STRIDE_Y = stride in y direction
  ! T = array containing the time at which the rupture front passes each point
  ! PROC_HOLDS_DATA = flag indicating if process holds data to be output
  ! DARRAY = MPI distributed array type

  type front_type
     logical :: front,proc_holds_data,ots_xmin,ots_xmax,ots_ymin,ots_ymax
     character(64) :: name
     character(4) :: field
     integer(pin) :: mx,px,my,py,nx,ny,stride_x,stride_y, &
          darray
     real(pr) :: val,xmin,xmax,ymin,ymax,dx,dy
     real(pr),dimension(:,:),allocatable :: t
     type(file_distributed) :: fh
  end type front_type

contains
  

  subroutine read_front(ninput,frt)
    ! READ_FRONT reads in front variables from file
    ! 
    ! Modified: 20 July 2010

    use io, only : error

    implicit none

    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! FRT = rupture front variables

    integer(pin),intent(in) :: ninput
    type(front_type),intent(out) :: frt
    
    ! Internal Parameters:
    ! FRONT = output rupture front history
    ! FIELD = field component to be contoured
    ! XMIN = minimum value of x
    ! XMAX = maximum value of x
    ! YMIN = minimum value of y
    ! YMAX = maximum value of y
    ! OTS_XMIN = whether or not to include the grid point 
    ! just outside of bound set by xmin
    ! OTS_XMAX = same as above for xmax
    ! OTS_YMIN = same as above for ymin
    ! OTS_YMAX = same as above for ymax
    ! STRIDE_X = stride in x direction
    ! STRIDE_Y = stride in y direction
    ! VAL = minimum value a point must exceed to denote 
    ! the passage of the rupture front
    ! STAT = I/O error flag

    logical :: front,ots_xmin,ots_xmax,ots_ymin,ots_ymax
    character(4) :: field
    integer(pin) :: stride_x,stride_y
    real(pr) :: val,xmin,xmax,ymin,ymax
    integer(pin) :: stat

    ! make namelist of user input variables

    namelist /front_list/ front,field,val,xmin,xmax,ymin,ymax, &
         ots_xmin,ots_xmax,ots_ymin,ots_ymax,stride_x,stride_y

    ! defaults
    
    front = .false.
    field = 'V'
    val = 1.d-3
    xmin = -15._pr
    xmax = 15._pr
    ymin = -15._pr
    ymax = 0._pr
    ots_xmin = .false.
    ots_xmax = .false.
    ots_ymin = .false.
    ots_ymax = .false.
    stride_x = 1
    stride_y = 1
    
    ! read namelist from input file, call error routine if needed
    
    rewind(ninput)
    read(ninput,nml=front_list,iostat=stat)
    if (stat>0) call error("Error reading namelist 'front_list' in .in file",'read_front')
    
    ! assign input variables to components of derived type
    
    frt%front = front
    frt%field = field
    frt%val = val
    frt%xmin = xmin
    frt%xmax = xmax
    frt%ymin = ymin
    frt%ymax = ymax
    frt%ots_xmin = ots_xmin
    frt%ots_xmax = ots_xmax
    frt%ots_ymin = ots_ymin
    frt%ots_ymax = ots_ymax
    frt%stride_x = stride_x
    frt%stride_y = stride_y

  end subroutine read_front

  
  subroutine init_front(necho,pb_name,frt,mdl,fld)
    ! INIT_FRONT initializes front variables
    ! 
    ! Modified: 20 July 2010

    use model, only : model_type
    use io, only : write_matlab
    use fields, only : fields_type,check_field
    use fft_routines, only : proc_holds_x
    use mpi_routines, only : is_master,MPI_REAL_PSAV,psav,subarray

    implicit none

    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! PB_NAME = problem name
    ! FRT = rupture front variables
    ! MDL = model variables
    ! FLD = fields variables

    integer(pin),intent(in) :: necho
    character(*),intent(in) :: pb_name
    type(front_type),intent(inout) :: frt
    type(model_type),intent(in) :: mdl
    type(fields_type),intent(in) :: fld

    ! return if not using this output method

    if (.not.frt%front) return

     ! initialize distributed output 

    call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PSAV,frt%darray)

    ! check if memory is allocated to field

    call check_field(frt%field,fld)

    ! get integer bounds and update min,max

    ! check if process holds data
    frt%proc_holds_data = proc_holds_x()
    if (.not.frt%proc_holds_data) return

    ! for now, include all points in x
    frt%mx = mdl%mx
    frt%px = mdl%px
    frt%xmin = mdl%x(frt%mx)
    frt%xmax = mdl%x(frt%px)
    frt%nx = 1+frt%px-frt%mx
    frt%dx = mdl%dx
    frt%stride_x = 1

    ! for now, include all points in y
    frt%my = 1
    frt%py = mdl%ny
    frt%ymin = mdl%y(frt%my)
    frt%ymax = mdl%y(frt%py)
    frt%ny = 1+frt%py-frt%my
    frt%dy = mdl%dy
    frt%stride_y = 1

    ! allocate memory for rupture front time array,
    ! initialize array with t=0

    allocate(frt%t(frt%mx:frt%px,frt%ny))
    frt%t = 0._pr

    ! name of data file

    frt%name = 'data/' // trim(adjustl(pb_name)) // '_front.dat'

    ! output variable values into matlab file

    if (is_master) then
       
       call write_matlab(necho,'name','front','frt')
       call write_matlab(necho,'field','t','frt')
       call write_matlab(necho,'val',frt%val,'frt')
       call write_matlab(necho,'xmin',mdl%x(1),'frt')
       call write_matlab(necho,'xmax',mdl%x(mdl%nx),'frt')
       call write_matlab(necho,'ymin',mdl%y(1),'frt')
       call write_matlab(necho,'ymax',mdl%y(mdl%ny),'frt')
       call write_matlab(necho,'zmin',0._pr,'frt')
       call write_matlab(necho,'zmax',0._pr,'frt')
       call write_matlab(necho,'nx',mdl%nx,'frt')
       call write_matlab(necho,'ny',mdl%ny,'frt')
       call write_matlab(necho,'nz',1,'frt')
       call write_matlab(necho,'dx',frt%dx,'frt')
       call write_matlab(necho,'dy',frt%dy,'frt')
       call write_matlab(necho,'dz',0._pr,'frt')
       call write_matlab(necho,'tmin',0._pr,'frt')
       call write_matlab(necho,'tmax',0._pr,'frt')

    end if

  end subroutine init_front

 
  subroutine destroy_front(frt)
    ! DESTROY_FRONT destroys derived type frt
    ! 
    ! Modified: 6 August 2010

    use io, only : close_file_distributed

    implicit none

    ! I/O Parameters:
    ! FRT = rupture front variables

    type(front_type),intent(inout) :: frt

    ! return if not using this output method

    if (.not.frt%front) return

    ! close output file (NOT REQUIRED, AS FILE IS ALWAYS CLOSED AFTER WRITE)

    !if(frt%proc_holds_data) call close_file_distributed(frt%fh)

    ! deallocate memory assigned to allocatable arrays

    if (allocated(frt%t)) deallocate(frt%t)

  end subroutine destroy_front


  subroutine output_front(frt)
    ! OUTPUT_FRONT writes rupture front data to output file
    ! 
    ! Modified: 6 August 2010

    use io, only : open_file_distributed,write_file_distributed,close_file_distributed
    use constants, only : psav
    use mpi

    implicit none

    ! I/O Parameters:
    ! FRT = rupture front variables

    type(front_type),intent(inout) :: frt

    ! return if not using this output method

    if (.not.frt%front) return

    ! proceed only if this process holds data

    if (.not.frt%proc_holds_data) return

    ! open file, clear contents

    call open_file_distributed(frt%fh,frt%name,'write',MPI_COMM_WORLD,frt%darray,psav)

    ! write data

    call write_file_distributed(frt%fh,frt%t)

    ! close file

    call close_file_distributed(frt%fh)

  end subroutine output_front


  subroutine set_front(frt,mdl,fld,t,dt)
    ! SET_FRONT marks progress of rupture front, assigns new
    ! value of time if front has passed a particular point
    ! 
    ! Modified: 3 March 2008

    use constants, only : zero
    use model, only : model_type
    use fields, only : fields_type,get_field,get_output_time

    implicit none

    ! I/O Parameters:
    ! FRT = rupture front variables
    ! MDL = model variables
    ! FLD = fields variables
    ! T = time of current time step 
    ! DT = time step

    type(model_type),intent(in) :: mdl
    type(fields_type),intent(in) :: fld
    type(front_type),intent(inout) :: frt
    real(pr),intent(in) :: t,dt

    ! Internal Parameters:
    ! TIM = time corresponding to particular field
    ! DATA = strided section of field array
    
    real(pr) :: tim
    real(pr),dimension(:,:),allocatable :: data

    ! return if not using this output method

    if (.not.frt%front) return

    ! search everywhere that has not already been set
    ! if rupture front condition is met, then store current time
    ! in the appropriate spot in the array

    allocate(data(frt%nx,frt%ny))

    call get_field(data,frt%field,mdl,fld,frt%nx,frt%ny,1, &
         frt%mx,frt%px,frt%my,frt%py,0,0,frt%stride_x,frt%stride_y,1)

    tim = get_output_time(frt%field,t,dt)

    where (data>=frt%val.and.frt%t==zero) frt%t = tim

    deallocate(data)

  end subroutine set_front


end module front

