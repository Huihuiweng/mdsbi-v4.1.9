! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module mesh

  ! MESH contains variables and routines pertaining to outputting data on a mesh
  ! 
  ! Modified: 19 July 2010

  use constants, only : pr,pin
  use io, only : file_distributed

  implicit none

  ! MESH_ITEM is a derived type containing mesh variables
  !
  ! Parameters:
  ! NAME = name of item
  ! FIELD = field components to output
  ! FH = file handle
  ! FHT = file handle for time
  ! SX = stride in x direction 
  ! (e.g., 2 means output every other point)
  ! SY = stride in y direction
  ! SZ = stride in z direction
  ! ST = stride in time
  ! N = number of time steps within time output window
  ! MX = minimum index in x direction
  ! PX = maximum index in x direction
  ! MY = minimum index in y direction
  ! PY = maximum index in y direction
  ! MZ = minimum index in z direction
  ! PZ = maximum index in z direction
  ! NX = number of strided grid points to be output in x direction
  ! NY = number of strided grid points to be output in y direction
  ! NZ = number of strided grid points to be output in z direction
  ! IX = interpolate in x direction
  ! IY = interpolate in y direction
  ! IZ = interpolate in z direction
  ! XMIN = minimum value of x
  ! XMAX = maximum value of x
  ! YMIN = minimum value of y
  ! YMAX = maximum value of y
  ! ZMIN = minimum value of z
  ! ZMAX = maximum value of z
  ! TMIN = minimum value of t
  ! TMAX = maximum value of t
  ! DX = spatial step in x direction
  ! DY = spatial step in y direction
  ! DZ = spatial step in z direction
  ! FX = normalized position of output in x direction (used in interpolation)
  ! FY = normalized position of output in y direction
  ! FZ = normalized position of output in z direction
  ! OTS_XMIN = whether or not to include the grid point just outside of bound set by xmin
  ! OTS_XMAX = same as above for xmax
  ! OTS_YMIN = same as above for ymin
  ! OTS_YMAX = same as above for ymax
  ! OTS_ZMIN = same as above for zmin
  ! OTS_ZMAX = same as above for zmax
  ! DISTRIBUTED = flag indicating if data is distributed across multiples processes or not
  ! PROC_HOLDS_DATA = flag indicating if process holds data to be output
  ! DARRAY = MPI distributed array type

  type :: mesh_item
     character(64) :: name
     character(5) :: field
     logical :: distributed,proc_holds_data,ix,iy,iz, &
          ots_xmin,ots_xmax,ots_ymin,ots_ymax,ots_zmin,ots_zmax
     integer(pin) :: n,sx,sy,sz,st, &
          mx,px,my,py,mz,pz,nx,ny,nz, &
          darray
     real(pr) :: xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz,tmin,tmax, &
          fx,fy,fz
     type(file_distributed) :: fh,fht
  end type mesh_item

  ! MESH_NODE is a derived type for basic node in singly-linked list 
  ! (single mesh item variables and pointer to next node)
  !
  ! Parameters:
  ! MS = mesh variables
  ! NEXT = pointer to next node

  type :: mesh_node
     private
     type(mesh_item),pointer :: ms=>null()
     type(mesh_node),pointer :: next=>null()
  end type mesh_node

  ! MESH_TYPE is a derived type containing all mesh variables
  !
  ! Parameters:
  ! ROOT = start of singly-linked list of mesh items

  type :: mesh_type
     type(mesh_node),pointer :: root=>null()
  end type mesh_type

contains


  subroutine read_mesh(ninput,msh)
    ! READ_MESH reads in mesh variables from file
    ! 
    ! Modified: 19 July 2010

    use utilities, only : word_count,extract_words
    use io, only : new_io_unit,error

    implicit none

    ! I/O Parameters:
    ! NINPUT = unit for input file
    ! MSH = mesh variables

    integer(pin),intent(in) :: ninput
    type(mesh_type),intent(out) :: msh

    ! Internal Parameters:
    ! NAME = name of mesh item
    ! TEMP = temporary character string
    ! IFLD = index of field
    ! NFLD = number of fields
    ! FIELDS = field names
    ! MS = mesh item

    integer(pin) :: ifld,nfld
    character(64) :: name,temp
    character(5),dimension(:),allocatable :: fields
    type(mesh_item) :: ms

    ! move to start of mesh list

    rewind(ninput)
    do
       read(ninput,*,end=100) temp
       if (temp=='!---begin:mesh_list---') exit
    end do
       
    ! read mesh data and add nodes to list

    do

       ! read item data

       read(ninput,*,end=100) name

       if (name=='!---end:mesh_list---') return
       
       read(ninput,'(a)') temp
       nfld = word_count(temp)

       allocate(fields(nfld))
       call extract_words(temp,fields)

       read(ninput,*) ms%xmin,ms%xmax,ms%sx,ms%ots_xmin,ms%ots_xmax
       read(ninput,*) ms%ymin,ms%ymax,ms%sy,ms%ots_ymin,ms%ots_ymax
       read(ninput,*) ms%zmin,ms%zmax,ms%sz,ms%ots_zmin,ms%ots_zmax
       read(ninput,*) ms%tmin,ms%tmax,ms%st

       ! store as nodes in list

       do ifld = 1,nfld
          ms%name = trim(name) // '_' // fields(ifld)
          ms%field = fields(ifld)
          call create_mesh_node(msh%root,ms)
       end do

       deallocate(fields)

    end do
    
100 call error("Error reading list 'mesh' in .in file",'read_mesh')

  end subroutine read_mesh


  recursive subroutine create_mesh_node(node,ms)
    ! CREATE_MESH_NODE recursively traverses list and adds node at end
    ! 
    ! Modified: 4 August 2007

    implicit none

    ! I/O Parameters
    ! NODE = mesh node
    ! MS = mesh item

    type(mesh_node),pointer :: node
    type(mesh_item),intent(in) :: ms
    
    if (.not.associated(node)) then
       allocate(node)
       nullify(node%next)
       allocate(node%ms)
       node%ms = ms
    else
        call create_mesh_node(node%next,ms)
    endif

  end subroutine create_mesh_node


  subroutine init_mesh(necho,pb_name,msh,mdl,fld)
    ! INIT_MESH initializes mesh variables and creates mesh items
    ! 
    ! Modified: 13 September 2007

    use model, only : model_type
    use fields, only : fields_type

    implicit none
    
    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! PB_NAME = name of problem
    ! MSH = mesh variables
    ! MDL = model variables
    ! FLD = fields variables

    integer(pin),intent(in) :: necho
    character(*),intent(in) :: pb_name
    type(mesh_type),intent(inout) :: msh
    type(model_type),intent(in) :: mdl
    type(fields_type),intent(in) :: fld

    ! traverse list and initialize all nodes

    call init_mesh_node(msh%root,necho,pb_name,mdl,fld)

  end subroutine init_mesh

  
  recursive subroutine init_mesh_node(node,necho,pb_name,mdl,fld)
    ! INIT_MESH_NODE recursively initializes all mesh nodes
    ! 
    ! Modified: 4 August 2007

    use model, only : model_type
    use fields, only : fields_type

    implicit none
    
    ! I/O Parameters:
    ! NODE = mesh node
    ! NECHO = unit number for output file
    ! PB_NAME = name of problem
    ! MDL = model variables
    ! FLD = fields variables

    type(mesh_node),pointer :: node
    integer(pin),intent(in) :: necho
    character(*),intent(in) :: pb_name
    type(model_type),intent(in) :: mdl
    type(fields_type),intent(in) :: fld

    ! exit recursion if at end of list

    if (.not.associated(node)) return

    ! allocate and initialize mesh item

    call init_mesh_item(necho,node%ms,mdl,fld)

    ! open file for output
    
    call open_mesh_file(node%ms,pb_name)

    ! repeat for next node

    call init_mesh_node(node%next,necho,pb_name,mdl,fld)    

  end subroutine init_mesh_node


  subroutine init_mesh_item(necho,ms,mdl,fld)
    ! INIT_MESH_ITEM creates and initializes a mesh item
    ! 
    ! Modified: 19 July 2010

    use model, only : model_type,get_bound_point
    use io, only : write_matlab,error
    use fields, only : fields_type,check_field
    use constants, only : psav,zero,half
    use fft_routines, only : proc_holds_x
    use mpi_routines, only : is_master,MPI_REAL_PSAV,subarray
    use mpi

    implicit none
    
    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! MS = mesh item
    ! MDL = model variables
    ! FLD = fields variables

    integer(pin),intent(in) :: necho
    type(mesh_item),intent(inout) :: ms
    type(model_type),intent(in) :: mdl
    type(fields_type),intent(in) :: fld

    integer(pin) :: ierr

    ! check if memory is allocated to field

    call check_field(ms%field,fld)

    ! no interpolation

    ms%ix = .false.
    ms%iy = .false.
    ms%iz = .false.

    ! check if data is from single point, otherwise initialize distributed output 

    if (ms%xmin==ms%xmax.and.ms%ymin==ms%ymax) then
       ms%distributed = .false.
       select case(ms%field)
       case('p','T')
          call MPI_Type_contiguous(mdl%nz,MPI_REAL_PSAV,ms%darray,ierr)
          call MPI_Type_commit(ms%darray,ierr)
       case default
          ms%darray = MPI_REAL_PSAV
       end select
    else
       ms%distributed = .true.
       select case(ms%field)
       case('p','T')
          call subarray(mdl%nx,mdl%ny,mdl%nz,mdl%mx,mdl%px,1,mdl%ny,1,mdl%nz,MPI_REAL_PSAV,ms%darray)
       case default
          call subarray(mdl%nx,mdl%ny,mdl%mx,mdl%px,1,mdl%ny,MPI_REAL_PSAV,ms%darray)
       end select
    end if

    ! get integer bounds and update min,max when not interpolating
    ! NO INTERPOLATION NOW -- FEATURE WAS DISABLED WHEN MPI SUPPORT ADDED

    if (ms%xmin==ms%xmax) then
       ! for now, find closest point (no interpolation)
       call get_bound_point(ms%field,'x',ms%xmin,ms%ots_xmin,ms%ots_xmax,mdl,ms%mx,ms%px,ms%fx)
       ! ensure mx==px (single point)
       if (ms%mx/=ms%px) then
          if (ms%fx<half) then
             ms%px = ms%mx ! round down
          else
             ms%mx = ms%px ! round up
          end if
       end if
       ms%xmin = mdl%x(ms%mx)
       ms%xmax = mdl%x(ms%px)
       ms%nx = 1
       ms%dx = zero
       ms%proc_holds_data = (mdl%mx<=ms%mx.and.ms%mx<=mdl%px)
    else
       ! check if process holds data
       ms%proc_holds_data = proc_holds_x()
       if (.not.ms%proc_holds_data) return
       ! for now, include all points in x
       ms%mx = mdl%mx
       ms%px = mdl%px
       ms%xmin = mdl%x(ms%mx)
       ms%xmax = mdl%x(ms%px)
       ms%nx = 1+ms%px-ms%mx
       ms%dx = mdl%dx
       ms%fx = zero
    end if

    if (ms%ymin==ms%ymax) then
       ! for now, find closest point (no interpolation)
       call get_bound_point(ms%field,'y',ms%ymin,ms%ots_ymin,ms%ots_ymax,mdl,ms%my,ms%py,ms%fy)
       ! ensure my==py (single point)
       if (ms%my/=ms%py) then
          if (ms%fy<half) then
             ms%py = ms%my ! round down
          else
             ms%my = ms%py ! round up
          end if
       end if
       ms%ymin = mdl%y(ms%my)
       ms%ymax = mdl%y(ms%py)
       ms%ny = 1
       ms%dy = zero
   else
       ! for now, include all points in y
       ms%my = 1
       ms%py = mdl%ny
       ms%ymin = mdl%y(ms%my)
       ms%ymax = mdl%y(ms%py)
       ms%ny = 1+ms%py-ms%my
       ms%dy = mdl%dy
       ms%fy = zero
    end if

    if (ms%zmin==ms%zmax) then
       call get_bound_point(ms%field,'z',ms%zmin,ms%ots_zmin,ms%ots_zmax,mdl,ms%mz,ms%pz,ms%fz)
       ! ensure mz==pz (single point)
       if (ms%mz/=ms%pz) then
          if (ms%fz<half) then
             ms%pz = ms%mz ! round down
          else
             ms%mz = ms%pz ! round up
          end if
       end if
       ms%zmin = mdl%z(ms%mz)
       ms%zmax = mdl%z(ms%pz)
       ms%nz = 1
       ms%dz = zero
    else
       ! for now, include all points in z
       ms%mz = 1
       ms%pz = mdl%nz
       ms%zmin = mdl%z(ms%mz)
       ms%zmax = mdl%z(ms%pz)
       ms%nz = 1+ms%pz-ms%mz
       ms%dz = mdl%dz
       ms%fz = zero
    end if

    ! zero number of elapsed time steps within output window

    ms%n = 0

    if (is_master) then

       call write_matlab(necho,'name',ms%name,ms%name)

       call write_matlab(necho,'field',ms%field,ms%name)

       if (ms%distributed) then
          call write_matlab(necho,'nx',mdl%nx,ms%name)
          call write_matlab(necho,'xmin',mdl%x(1),ms%name)
          call write_matlab(necho,'xmax',mdl%x(mdl%nx),ms%name)
       else
          call write_matlab(necho,'nx',ms%nx,ms%name)
          call write_matlab(necho,'xmin',ms%xmin,ms%name)
          call write_matlab(necho,'xmax',ms%xmax,ms%name)
       end if
       call write_matlab(necho,'dx',ms%dx,ms%name)
       call write_matlab(necho,'stride_x',ms%sx,ms%name)

       if (ms%distributed) then
          call write_matlab(necho,'ny',mdl%ny,ms%name)
          call write_matlab(necho,'ymin',mdl%y(1),ms%name)
          call write_matlab(necho,'ymax',mdl%y(mdl%ny),ms%name)
       else
          call write_matlab(necho,'ny',ms%ny,ms%name)
          call write_matlab(necho,'ymin',ms%ymin,ms%name)
          call write_matlab(necho,'ymax',ms%ymax,ms%name)
       end if
       call write_matlab(necho,'dy',ms%dy,ms%name)
       call write_matlab(necho,'stride_y',ms%sy,ms%name)

       if (ms%distributed) then
          call write_matlab(necho,'nz',ms%nz,ms%name)
          call write_matlab(necho,'zmin',mdl%z(1),ms%name)
          call write_matlab(necho,'zmax',mdl%z(mdl%nz),ms%name)
       else
          call write_matlab(necho,'nz',ms%nz,ms%name)
          call write_matlab(necho,'zmin',ms%zmin,ms%name)
          call write_matlab(necho,'zmax',ms%zmax,ms%name) 
       end if
       call write_matlab(necho,'stride_z',ms%sz,ms%name)
       call write_matlab(necho,'dz',ms%dz,ms%name)

       call write_matlab(necho,'tmin',ms%tmin,ms%name)
       call write_matlab(necho,'tmax',ms%tmax,ms%name)
       call write_matlab(necho,'stride_t',ms%st,ms%name)

    end if

  end subroutine init_mesh_item


  subroutine open_mesh_file(ms,pb_name)
    ! OPEN_MESH_FILE opens a mesh file
    ! 
    ! Modified: 19 July 2010

    use io, only : open_file_distributed
    use mpi_routines, only : is_master,MPI_REAL_PSAV,psav
    use mpi

    implicit none

    ! I/O Parameters:
    ! MS = mesh item
    ! PB_NAME = problem name
    
    type(mesh_item),intent(inout) :: ms
    character(*),intent(in) :: pb_name
    
    ! Internal Parameters:
    ! NAME = file name
    ! COMM = MPI communicator

    character(64) :: name
    integer(pin) :: comm

    ! open output time file and set file size to zero (to remove any pre-existing files of same name)

    if (is_master) then
       name = 'data/' // trim(adjustl(pb_name)) // '_' // trim(adjustl(ms%name)) // '.t'
       call open_file_distributed(ms%fht,name,'write',MPI_COMM_SELF,MPI_REAL_PSAV,psav)
    end if

    ! open data file and set file size to zero (to remove any pre-existing files of same name)

    if (ms%distributed) then
       comm = MPI_COMM_WORLD
    else
       comm = MPI_COMM_SELF
       if (.not.ms%proc_holds_data) return
    end if

    name = 'data/' // trim(adjustl(pb_name)) // '_' // trim(adjustl(ms%name)) // '.dat'
    call open_file_distributed(ms%fh,name,'write',comm,ms%darray,psav)

  end subroutine open_mesh_file


  subroutine output_mesh(mdl,fld,msh,t)
    ! OUTPUT_MESH obtains array of field values, appropriately strided, and calls
    ! a routine to write them to a data file
    !
    ! Modified: 13 September 2007

    use constants, only : psav
    use model, only : model_type
    use fields, only : fields_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLD = fields variables
    ! MSH = mesh variables
    ! T = current time

    type(model_type),intent(in) :: mdl
    type(fields_type),intent(in) :: fld
    type(mesh_type),intent(inout) :: msh
    real(pr),intent(in) :: t

    ! recursively traverse list to output all mesh items

    call output_mesh_node(msh%root,mdl,fld,t)

  end subroutine output_mesh


  recursive subroutine output_mesh_node(node,mdl,fld,t)
    ! OUTPUT_MESH_NODE recursively outputs data for nodes in list
    !
    ! Modified: 19 July 2010

    use model, only : model_type
    use fields, only : fields_type,write_field
    use io, only : write_file_distributed
    use mpi_routines, only : is_master

    implicit none

    ! I/O Parameters:
    ! NODE = mesh node
    ! MDL = model variables
    ! FLD = fields variables
    ! T = current time

    type(mesh_node),pointer :: node
    type(model_type),intent(in) :: mdl
    type(fields_type),intent(in) :: fld
    real(pr),intent(in) :: t

    ! proceed only if node exists

    if (.not.associated(node)) return

    ! proceed only if current time is within output time window
    
    if (node%ms%tmin<=t.and.t<=node%ms%tmax) then
       
       ! tmin<=t<=tmax, so increment potential time output index
       
       node%ms%n = node%ms%n+1
       
       ! proceed only if stride of potential time output index is acceptable
       
       if (mod(node%ms%n-1,node%ms%st)==0) then
          
          ! write output time
          
          if (is_master) call write_file_distributed(node%ms%fht,(/t/))
          
          ! write data
          
          if (node%ms%proc_holds_data) &
               call write_field(node%ms%fh,node%ms%field,mdl,fld,node%ms%nx,node%ms%ny,node%ms%nz, &
               node%ms%mx,node%ms%px,node%ms%my,node%ms%py,node%ms%mz,node%ms%pz, &
               node%ms%sx,node%ms%sy,node%ms%sz)

       end if

    end if

    ! move to next mesh node
       
    call output_mesh_node(node%next,mdl,fld,t)

  end subroutine output_mesh_node


  subroutine destroy_mesh(msh)
    ! DESTROY_MESH destroys derived type msh
    ! 
    ! Modified: 13 September 2007

    implicit none

    ! I/O Parameters:
    ! MSH = mesh variables

    type(mesh_type),intent(inout) :: msh

    ! recursively traverse list, destroying each node (slow implementation)

    do while(associated(msh%root)) 
       call destroy_mesh_node(msh%root)
    end do

  end subroutine destroy_mesh


  recursive subroutine destroy_mesh_node(node)
    ! DESTROY_MESH_NODE recursively destroys all nodes in list
    ! 
    ! Modified: 21 September 2007

    use io, only : close_file_distributed
    use mpi_routines, only : is_master

    implicit none

    ! I/O Parameters:
    ! NODE = mesh node

    type(mesh_node),pointer :: node

    ! move to end of list

    if (associated(node%next)) then

       ! not at end, move to next node

       call destroy_mesh_node(node%next)

    else

       ! end of list, destroy node

       if (associated(node%ms)) then
          
          ! close data file
          
          if(node%ms%proc_holds_data) call close_file_distributed(node%ms%fh)
          
          ! close time file
          
          if (is_master) call close_file_distributed(node%ms%fht)
          
          ! deallocate and nullify mesh item
          
          deallocate(node%ms)
          node%ms => null()
          
       end if
       
       ! deallocate and nullify mesh node

       deallocate(node)
       node => null()
       
       return

    end if

  end subroutine destroy_mesh_node


end module mesh

