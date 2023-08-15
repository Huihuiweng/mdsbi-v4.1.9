! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module fields

  ! FIELDS contains fields defined on fault and on hydraulic
  ! grid used for thermal pressurization calculation, as well as
  ! routines for manipulating these fields
  ! 
  ! Modified: 19 July 2010

  use constants, only : pin,pr
  use fault, only : fault_type
  use thermpres, only : thermpres_type
  use load, only : load_type

  implicit none

  ! FIELDS_TYPE is a derived type containing fields variables
  !
  ! Parameters:
  ! FLT = fault variables
  ! TP = thermal pressurization variables
  ! LD = load variables
  ! BM = bimaterial or not

  type fields_type
     type(fault_type) :: flt
     type(thermpres_type) :: tp
     type(load_type) :: ld
     logical :: bm
  end type fields_type

contains


  subroutine read_fields(ninput,mdl,fld)
    ! READ_FIELDS reads in fields variables from file
    ! 
    ! Modified: 19 October 2007

    use model, only : model_type
    use fault, only : read_fault
    use thermpres, only : read_thermpres
    use load, only : read_load

    implicit none

    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! MDL = model variables
    ! FLD = fields variables

    integer(pin),intent(in) :: ninput
    type(model_type),intent(in) :: mdl
    type(fields_type),intent(inout) :: fld

    fld%bm = mdl%bm
    
    call read_fault(ninput,mdl,fld%flt)
    if (fld%flt%thermpres) call read_thermpres(ninput,fld%tp)
    if (fld%flt%dynload) call read_load(ninput,fld%ld)

  end subroutine read_fields


  subroutine init_fields(necho,mdl,fld)
    ! INIT_FIELDS initializes derived type fld
    ! 
    ! Modified: 19 October 2007

    use model, only : model_type
    use fault, only : init_fault
    use thermpres, only : init_thermpres
    use load, only : init_load

    implicit none

    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! MDL = model variables
    ! FLD = fields variables

    integer(pin),intent(in) :: necho
    type(model_type),intent(in) :: mdl
    type(fields_type),intent(inout) :: fld

    call init_fault(necho,mdl,fld%flt)
    if (fld%flt%dynload) call init_load(necho,mdl,fld%ld)
    if (fld%flt%thermpres) call init_thermpres(necho,mdl,fld%flt,fld%tp)

  end subroutine init_fields


  subroutine destroy_fields(fld)
    ! DESTROY_FIELDS destroys derived type fld
    ! 
    ! Modified: 19 September 2007

    use fault, only : destroy_fault
    use thermpres, only : destroy_thermpres
    use load, only : destroy_load

    implicit none

    ! I/O Parameters:
    ! FLD = fields variables

    type(fields_type),intent(inout) :: fld

    call destroy_load(fld%ld)
    call destroy_thermpres(fld%flt)
    call destroy_fault(fld%flt)

  end subroutine destroy_fields


  subroutine check_field(field,fld)
    ! CHECK_FIELD checks whether or not field is allocated
    ! 
    ! Modified: 19 July 2010

    use io, only : error

    implicit none

    ! I/O Parameters:
    ! FIELD = field
    ! FLD = fields variables

    character(*),intent(in) :: field
    type(fields_type),intent(in) :: fld

    ! Internal Parameters:
    ! ALO = flag to test whether field is allocated or not

    logical :: alo

    select case(field)

    case('U')
       alo = allocated(fld%flt%U)
    case('V')
       alo = allocated(fld%flt%V)
    case('O')
       alo = allocated(fld%flt%O)

    case('S')
       alo = allocated(fld%flt%S)
    case('N')
       alo = allocated(fld%flt%N)

    case('Ux')
       if (fld%bm) then
          alo = allocated(fld%flt%uxp).and.allocated(fld%flt%uxm)
       else
          alo = allocated(fld%flt%Ux)
       end if
    case('Uy')
       if (fld%bm) then
          alo = allocated(fld%flt%uyp).and.allocated(fld%flt%uym)
       else
          alo = allocated(fld%flt%Uy)
       end if
    case('Uz')
       if (fld%bm) then
          alo = allocated(fld%flt%uzp).and.allocated(fld%flt%uzm)
       else
          alo = .false.
       end if

    case('Vx')
       if (fld%bm) then
          alo = allocated(fld%flt%vxp).and.allocated(fld%flt%vxm)
       else
          alo = allocated(fld%flt%Vx)
       end if
    case('Vy')
       if (fld%bm) then
          alo = allocated(fld%flt%vyp).and.allocated(fld%flt%vym)
       else
          alo = allocated(fld%flt%Vy)
       end if
    case('Vz')
       if (fld%bm) then
          alo = allocated(fld%flt%vzp).and.allocated(fld%flt%vzm)
       else
          alo = .false.
       end if

    case('uxp')
       alo = allocated(fld%flt%uxp)
    case('uxm')
       alo = allocated(fld%flt%uxm)
    case('uyp')
       alo = allocated(fld%flt%uyp)
    case('uym')
       alo = allocated(fld%flt%uym)
    case('uzp')
       alo = allocated(fld%flt%uzp)
    case('uzm')
       alo = allocated(fld%flt%uzm)

    case('vxp')
       alo = allocated(fld%flt%vxp)
    case('vxm')
       alo = allocated(fld%flt%vxm)
    case('vyp')
       alo = allocated(fld%flt%vyp)
    case('vym')
       alo = allocated(fld%flt%vym)
    case('vzp')
       alo = allocated(fld%flt%vzp)
    case('vzm')
       alo = allocated(fld%flt%vzm)

    case('sx')
       alo = allocated(fld%flt%sx)
    case('sy')
       alo = allocated(fld%flt%sy)
    case('sz')
       alo = allocated(fld%flt%sz)

    case('sx0')
       alo = allocated(fld%flt%sx0)
    case('sy0')
       alo = allocated(fld%flt%sy0)
    case('sz0')
       alo = allocated(fld%flt%sz0)

    case('fx')
       alo = allocated(fld%flt%fx)
    case('fy')
       alo = allocated(fld%flt%fy)

    case('fxp')
       alo = allocated(fld%flt%fxp)
    case('fxm')
       alo = allocated(fld%flt%fxm)
    case('fyp')
       alo = allocated(fld%flt%fyp)
    case('fym')
       alo = allocated(fld%flt%fym)
    case('fzp')
       alo = allocated(fld%flt%fzp)
    case('fzm')
       alo = allocated(fld%flt%fzm)

    case('Q','logQ')
       alo = allocated(fld%flt%Q)
    case('dQ')
       alo = allocated(fld%flt%dQ)

    case('T0')
       alo = allocated(fld%flt%T0)
    case('p0')
       alo = allocated(fld%flt%p0)

    case('exxp')
       alo = allocated(fld%flt%exxp)
    case('exxm')
       alo = allocated(fld%flt%exxm)
    case('exyp')
       alo = allocated(fld%flt%exyp)
    case('exym')
       alo = allocated(fld%flt%exym)
    case('eyyp')
       alo = allocated(fld%flt%eyyp)
    case('eyym')
       alo = allocated(fld%flt%eyym)

    case('Yp')
       alo = allocated(fld%flt%Yp)
    case('Ym')
       alo = allocated(fld%flt%Ym)

    case('p')
       alo = allocated(fld%flt%p)
    case('dp')
       alo = allocated(fld%flt%dp)
    case('T')
       alo = allocated(fld%flt%T)
    case('dT')
       alo = allocated(fld%flt%dT)

    case default
       call error("Invalid field '" // field // "' to be output",'check_field')

    end select

    if (.not.alo) call error("Field '" // field // "' is not allocated and cannot be output",'check_field')

  end subroutine check_field


  subroutine get_field(data,field,mdl,fld,nx,ny,nz,mx,px,my,py,mz,pz,sx,sy,sz)
    ! GET_FIELD returns a section of a field array, primarily for data output routines
    ! 
    ! Modified: 19 July 2010

    use constants, only : zero
    use io, only : error
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! DATA = field array to be returned
    ! FIELD = field to be returned
    ! MDL = model variables
    ! FLD = fields variables
    ! NX = number of strided grid points in x direction
    ! NY = number of strided grid points in y direction
    ! NZ = number of strided grid points in z direction
    ! MX = minimum index of field array in x direction
    ! PX = maximum index of field array in x direction
    ! MY = minimum index of field array in y direction
    ! PY = maximum index of field array in y direction
    ! MZ = minimum index of field array in z direction
    ! PZ = maximum index of field array in z direction
    ! SX = stride in x direction
    ! SY = stride in y direction
    ! SZ = stride in z direction

    integer(pin),intent(in) :: nx,ny,nz,mx,px,my,py,mz,pz,sx,sy,sz
    real(pr),dimension(nx,ny,nz),intent(out) :: data
    character(*),intent(in) :: field
    type(model_type),intent(in) :: mdl
    type(fields_type),intent(in) :: fld

    select case(field)

    case('U')
       data(:,:,1) = fld%flt%U(mx:px:sx,my:py:sy,1)
    case('V')
       data(:,:,1) = fld%flt%V(mx:px:sx,my:py:sy,1)
    case('O')
       data(:,:,1) = fld%flt%O(mx:px:sx,my:py:sy,1)

    case('S')
       data(:,:,1) = fld%flt%S(mx:px:sx,my:py:sy,1)
    case('N')
       data(:,:,1) = fld%flt%N(mx:px:sx,my:py:sy,1)

    case('Ux')
       if (fld%bm) then
          data(:,:,1) = fld%flt%uxp(mx:px:sx,my:py:sy,1)-fld%flt%uxm(mx:px:sx,my:py:sy,1)
       else
          data(:,:,1) = fld%flt%Ux(mx:px:sx,my:py:sy,1)
       end if
    case('Uy')
       if (fld%bm) then
          data(:,:,1) = fld%flt%uyp(mx:px:sx,my:py:sy,1)-fld%flt%uym(mx:px:sx,my:py:sy,1)
       else
          data(:,:,1) = fld%flt%Uy(mx:px:sx,my:py:sy,1)
       end if
    case('Uz')
       if (fld%bm) then
          data(:,:,1) = fld%flt%uzp(mx:px:sx,my:py:sy,1)-fld%flt%uzm(mx:px:sx,my:py:sy,1)
       else
          data(:,:,1) = zero
       end if

    case('Vx')
       if (fld%bm) then
          data(:,:,1) = fld%flt%vxp(mx:px:sx,my:py:sy,1)-fld%flt%vxm(mx:px:sx,my:py:sy,1)          
       else
          data(:,:,1) = fld%flt%Vx(mx:px:sx,my:py:sy,1)
       end if
    case('Vy')
       if (fld%bm) then
          data(:,:,1) = fld%flt%vyp(mx:px:sx,my:py:sy,1)-fld%flt%vym(mx:px:sx,my:py:sy,1)
       else
          data(:,:,1) = fld%flt%Vy(mx:px:sx,my:py:sy,1)       
       end if
    case('Vz')
       if (fld%bm) then
          data(:,:,1) = fld%flt%vzp(mx:px:sx,my:py:sy,1)-fld%flt%vzm(mx:px:sx,my:py:sy,1)
       else
          data(:,:,1) = zero
       end if

    case('uxp')
       data(:,:,1) = fld%flt%uxp(mx:px:sx,my:py:sy,1)
    case('uxm')
       data(:,:,1) = fld%flt%uxm(mx:px:sx,my:py:sy,1)
    case('uyp')
       data(:,:,1) = fld%flt%uyp(mx:px:sx,my:py:sy,1)
    case('uym')
       data(:,:,1) = fld%flt%uym(mx:px:sx,my:py:sy,1)
    case('uzp')
       data(:,:,1) = fld%flt%uzp(mx:px:sx,my:py:sy,1)
    case('uzm')
       data(:,:,1) = fld%flt%uzm(mx:px:sx,my:py:sy,1)

    case('vxp')
       data(:,:,1) = fld%flt%vxp(mx:px:sx,my:py:sy,1)
    case('vxm')
       data(:,:,1) = fld%flt%vxm(mx:px:sx,my:py:sy,1)
    case('vyp')
       data(:,:,1) = fld%flt%vyp(mx:px:sx,my:py:sy,1)
    case('vym')
       data(:,:,1) = fld%flt%vym(mx:px:sx,my:py:sy,1)
    case('vzp')
       data(:,:,1) = fld%flt%vzp(mx:px:sx,my:py:sy,1)
    case('vzm')
       data(:,:,1) = fld%flt%vzm(mx:px:sx,my:py:sy,1)

    case('sx')
       data(:,:,1) = fld%flt%sx(mx:px:sx,my:py:sy,1)
    case('sy')
       data(:,:,1) = fld%flt%sy(mx:px:sx,my:py:sy,1)
    case('sz')
       data(:,:,1) = fld%flt%sz(mx:px:sx,my:py:sy,1)

    case('sx0')
       data(:,:,1) = fld%flt%sx0(mx:px:sx,my:py:sy)
    case('sy0')
       data(:,:,1) = fld%flt%sy0(mx:px:sx,my:py:sy)
    case('sz0')
       data(:,:,1) = fld%flt%sz0(mx:px:sx,my:py:sy)

    case('fx')
       data(:,:,1) = fld%flt%fx(mx:px:sx,my:py:sy,0)
    case('fy')
       data(:,:,1) = fld%flt%fy(mx:px:sx,my:py:sy,0)

    case('fxp')
       data(:,:,1) = fld%flt%fxp(mx:px:sx,my:py:sy,0)
    case('fxm')
       data(:,:,1) = fld%flt%fxm(mx:px:sx,my:py:sy,0)
    case('fyp')
       data(:,:,1) = fld%flt%fyp(mx:px:sx,my:py:sy,0)
    case('fym')
       data(:,:,1) = fld%flt%fym(mx:px:sx,my:py:sy,0)
    case('fzp')
       data(:,:,1) = fld%flt%fzp(mx:px:sx,my:py:sy,0)
    case('fzm')
       data(:,:,1) = fld%flt%fzm(mx:px:sx,my:py:sy,0)

    case('Q')
       data(:,:,1) = fld%flt%Q(mx:px:sx,my:py:sy,1)
    case('dQ')
       data(:,:,1) = fld%flt%dQ(mx:px:sx,my:py:sy,1)
    case('logQ')
       data(:,:,1) = log10(fld%flt%Q(mx:px:sx,my:py:sy,1))

    case('T0')
       data(:,:,1) = fld%flt%T0(mx:px:sx,my:py:sy,1)
    case('p0')
       data(:,:,1) = fld%flt%p0(mx:px:sx,my:py:sy,1)

    case('exxp')
       data(:,:,1) = fld%flt%exxp(mx:px:sx,my:py:sy)
    case('exxm')
       data(:,:,1) = fld%flt%exxm(mx:px:sx,my:py:sy)
    case('exyp')
       data(:,:,1) = fld%flt%exyp(mx:px:sx,my:py:sy)
    case('exym')
       data(:,:,1) = fld%flt%exym(mx:px:sx,my:py:sy)
    case('eyyp')
       data(:,:,1) = fld%flt%eyyp(mx:px:sx,my:py:sy)
    case('eyym')
       data(:,:,1) = fld%flt%eyym(mx:px:sx,my:py:sy)

    case('Yp')
       data(:,:,1) = fld%flt%Yp(mx:px:sx,my:py:sy)
    case('Ym')
       data(:,:,1) = fld%flt%Ym(mx:px:sx,my:py:sy)

    ! note for arrays below, z-dimension is switched from first (most efficient in calculations)
    ! to last (for compatibility with output routines)

    case('p')
       data = reshape(fld%flt%p( mz:pz:sz,mx:px:sx,my:py:sy,1),shape(data),order=(/2,3,1/))
    case('T')
       data = reshape(fld%flt%T( mz:pz:sz,mx:px:sx,my:py:sy,1),shape(data),order=(/2,3,1/))
    case('dp')
       data = reshape(fld%flt%dp(mz:pz:sz,mx:px:sx,my:py:sy,1),shape(data),order=(/2,3,1/))
    case('dT')
       data = reshape(fld%flt%dT(mz:pz:sz,mx:px:sx,my:py:sy,1),shape(data),order=(/2,3,1/))

    case('rerr')
       data(:,:,1) = fld%flt%rerr(1)
    case('aerr')
       data(:,:,1) = fld%flt%aerr(1)

    case default
       call error("Invalid field '" // field // "' to be output",'get_field')

    end select

  end subroutine get_field


  subroutine write_field(fh,field,mdl,fld,nx,ny,nz,mx,px,my,py,mz,pz,sx,sy,sz)
    ! WRITE_FIELD writes a section of a field array
    ! 
    ! Modified: 19 July 2010

    use constants, only : zero
    use io, only : error,file_distributed,write_file_distributed
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! FH = file handle
    ! FIELD = field to be returned
    ! MDL = model variables
    ! FLD = fields variables
    ! NX = number of strided grid points in x direction
    ! NY = number of strided grid points in y direction
    ! NZ = number of strided grid points in z direction
    ! MX = minimum index of field array in x direction
    ! PX = maximum index of field array in x direction
    ! MY = minimum index of field array in y direction
    ! PY = maximum index of field array in y direction
    ! MZ = minimum index of field array in z direction
    ! PZ = maximum index of field array in z direction
    ! SX = stride in x direction
    ! SY = stride in y direction
    ! SZ = stride in z direction

    type(file_distributed),intent(in) :: fh
    integer(pin),intent(in) :: nx,ny,nz,mx,px,my,py,mz,pz,sx,sy,sz
    character(*),intent(in) :: field
    type(model_type),intent(in) :: mdl
    type(fields_type),intent(in) :: fld

    select case(field)

    case('U')
       call write_file_distributed(fh,fld%flt%U(mx:px:sx,my:py:sy,1))
    case('V')
       call write_file_distributed(fh,fld%flt%V(mx:px:sx,my:py:sy,1))
    case('O')
       call write_file_distributed(fh,fld%flt%O(mx:px:sx,my:py:sy,1))

    case('S')
       call write_file_distributed(fh,fld%flt%S(mx:px:sx,my:py:sy,1))
    case('N')
       call write_file_distributed(fh,fld%flt%N(mx:px:sx,my:py:sy,1))

    case('Ux')
       if (fld%bm) then
          call write_file_distributed(fh, &
               fld%flt%uxp(mx:px:sx,my:py:sy,1)- &
               fld%flt%uxm(mx:px:sx,my:py:sy,1))
       else
          call write_file_distributed(fh,fld%flt%Ux(mx:px:sx,my:py:sy,1))
       end if
    case('Uy')
       if (fld%bm) then
          call write_file_distributed(fh, &
               fld%flt%uyp(mx:px:sx,my:py:sy,1)- &
               fld%flt%uym(mx:px:sx,my:py:sy,1))
       else
          call write_file_distributed(fh,fld%flt%Uy(mx:px:sx,my:py:sy,1))
       end if
    case('Uz')
       if (fld%bm) then
          call write_file_distributed(fh, &
               fld%flt%uzp(mx:px:sx,my:py:sy,1)- &
               fld%flt%uzm(mx:px:sx,my:py:sy,1))
       else
          call write_file_distributed(fh,spread(spread(zero,1,nx),2,ny))
       end if

    case('Vx')
       if (fld%bm) then
          call write_file_distributed(fh, &
               fld%flt%vxp(mx:px:sx,my:py:sy,1)- &
               fld%flt%vxm(mx:px:sx,my:py:sy,1))          
       else
          call write_file_distributed(fh,fld%flt%Vx(mx:px:sx,my:py:sy,1))
       end if
    case('Vy')
       if (fld%bm) then
          call write_file_distributed(fh, &
               fld%flt%vyp(mx:px:sx,my:py:sy,1)- &
               fld%flt%vym(mx:px:sx,my:py:sy,1))
       else
          call write_file_distributed(fh,fld%flt%Vy(mx:px:sx,my:py:sy,1))      
       end if
    case('Vz')
       if (fld%bm) then
          call write_file_distributed(fh, &
               fld%flt%vzp(mx:px:sx,my:py:sy,1)- &
               fld%flt%vzm(mx:px:sx,my:py:sy,1))
       else
          call write_file_distributed(fh,spread(spread(zero,1,nx),2,ny))
       end if

    case('uxp')
       call write_file_distributed(fh,fld%flt%uxp(mx:px:sx,my:py:sy,1))
    case('uxm')
       call write_file_distributed(fh,fld%flt%uxm(mx:px:sx,my:py:sy,1))
    case('uyp')
       call write_file_distributed(fh,fld%flt%uyp(mx:px:sx,my:py:sy,1))
    case('uym')
       call write_file_distributed(fh,fld%flt%uym(mx:px:sx,my:py:sy,1))
    case('uzp')
       call write_file_distributed(fh,fld%flt%uzp(mx:px:sx,my:py:sy,1))
    case('uzm')
       call write_file_distributed(fh,fld%flt%uzm(mx:px:sx,my:py:sy,1))

    case('vxp')
       call write_file_distributed(fh,fld%flt%vxp(mx:px:sx,my:py:sy,1))
    case('vxm')
       call write_file_distributed(fh,fld%flt%vxm(mx:px:sx,my:py:sy,1))
    case('vyp')
       call write_file_distributed(fh,fld%flt%vyp(mx:px:sx,my:py:sy,1))
    case('vym')
       call write_file_distributed(fh,fld%flt%vym(mx:px:sx,my:py:sy,1))
    case('vzp')
       call write_file_distributed(fh,fld%flt%vzp(mx:px:sx,my:py:sy,1))
    case('vzm')
       call write_file_distributed(fh,fld%flt%vzm(mx:px:sx,my:py:sy,1))

    case('sx')
       call write_file_distributed(fh,fld%flt%sx(mx:px:sx,my:py:sy,1))
    case('sy')
       call write_file_distributed(fh,fld%flt%sy(mx:px:sx,my:py:sy,1))
    case('sz')
       call write_file_distributed(fh,fld%flt%sz(mx:px:sx,my:py:sy,1))

    case('sx0')
       call write_file_distributed(fh,fld%flt%sx0(mx:px:sx,my:py:sy))
    case('sy0')
       call write_file_distributed(fh,fld%flt%sy0(mx:px:sx,my:py:sy))
    case('sz0')
       call write_file_distributed(fh,fld%flt%sz0(mx:px:sx,my:py:sy))

    case('fx')
       call write_file_distributed(fh,fld%flt%fx(mx:px:sx,my:py:sy,0))
    case('fy')
       call write_file_distributed(fh,fld%flt%fy(mx:px:sx,my:py:sy,0))

    case('fxp')
       call write_file_distributed(fh,fld%flt%fxp(mx:px:sx,my:py:sy,0))
    case('fxm')
       call write_file_distributed(fh,fld%flt%fxm(mx:px:sx,my:py:sy,0))
    case('fyp')
       call write_file_distributed(fh,fld%flt%fyp(mx:px:sx,my:py:sy,0))
    case('fym')
       call write_file_distributed(fh,fld%flt%fym(mx:px:sx,my:py:sy,0))
    case('fzp')
       call write_file_distributed(fh,fld%flt%fzp(mx:px:sx,my:py:sy,0))
    case('fzm')
       call write_file_distributed(fh,fld%flt%fzm(mx:px:sx,my:py:sy,0))

    case('Q')
       call write_file_distributed(fh,fld%flt%Q(mx:px:sx,my:py:sy,1))
    case('dQ')
       call write_file_distributed(fh,fld%flt%dQ(mx:px:sx,my:py:sy,1))
    case('logQ')
       call write_file_distributed(fh,log10(fld%flt%Q(mx:px:sx,my:py:sy,1)))

    case('T0')
       call write_file_distributed(fh,fld%flt%T0(mx:px:sx,my:py:sy,1))
    case('p0')
       call write_file_distributed(fh,fld%flt%p0(mx:px:sx,my:py:sy,1))

    case('exxp')
       call write_file_distributed(fh,fld%flt%exxp(mx:px:sx,my:py:sy))
    case('exxm')
       call write_file_distributed(fh,fld%flt%exxm(mx:px:sx,my:py:sy))
    case('exyp')
       call write_file_distributed(fh,fld%flt%exyp(mx:px:sx,my:py:sy))
    case('exym')
       call write_file_distributed(fh,fld%flt%exym(mx:px:sx,my:py:sy))
    case('eyyp')
       call write_file_distributed(fh,fld%flt%eyyp(mx:px:sx,my:py:sy))
    case('eyym')
       call write_file_distributed(fh,fld%flt%eyym(mx:px:sx,my:py:sy))

    case('Yp')
       call write_file_distributed(fh,fld%flt%Yp(mx:px:sx,my:py:sy))
    case('Ym')
       call write_file_distributed(fh,fld%flt%Ym(mx:px:sx,my:py:sy))

    ! note for arrays below, z-dimension is switched from first (most efficient in calculations)
    ! to last (for compatibility with output routines)

    case('p')
       call write_file_distributed(fh, &
            reshape(fld%flt%p( mz:pz:sz,mx:px:sx,my:py:sy,1),(/nx,ny,nz/),order=(/2,3,1/)))
    case('T')
       call write_file_distributed(fh, &
            reshape(fld%flt%T( mz:pz:sz,mx:px:sx,my:py:sy,1),(/nx,ny,nz/),order=(/2,3,1/)))
    case('dp')
       call write_file_distributed(fh, &
            reshape(fld%flt%dp(mz:pz:sz,mx:px:sx,my:py:sy,1),(/nx,ny,nz/),order=(/2,3,1/)))
    case('dT')
       call write_file_distributed(fh, &
            reshape(fld%flt%dT(mz:pz:sz,mx:px:sx,my:py:sy,1),(/nx,ny,nz/),order=(/2,3,1/)))

    case('rerr')
       call write_file_distributed(fh,(/fld%flt%rerr(1)/))
    case('aerr')
       call write_file_distributed(fh,(/fld%flt%aerr(1)/))

    case default
       call error("Invalid field '" // field // "' to be output",'write_field')

    end select

  end subroutine write_field


  function get_output_time(field,t,dt) result(time)
    ! GET_OUTPUT_TIME returns a current value of time for a particular field,
    ! trivial for this code but included for compatility with staggered
    ! grid codes (for which fields are staggered in time)
    ! 
    ! Modified: 24 October 2006

    use constants, only : zero

    implicit none

    ! I/O Parameters:
    ! FIELD = field
    ! T = current time
    ! DT = time step
    ! TIME = time returned by function

    character(*),intent(in) :: field
    real(pr),intent(in) :: t,dt
    real(pr) :: time

    select case(field)
    case default
       time = t+zero*dt
    end select

  end function get_output_time


end module fields
    
