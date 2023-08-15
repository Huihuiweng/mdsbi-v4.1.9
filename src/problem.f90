! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module problem

  ! PROBLEM contains routines for creating and destroying a problem and its
  ! constituent objects
  ! 
  ! Modified: 9 July 2007

  use model, only : model_type
  use convolution, only : convolution_type
  use fields, only : fields_type
  use kernel, only : kernel_type
  use friction, only : friction_type
  use mesh, only : mesh_type
  use front, only : front_type

  implicit none

  ! PROBLEM_TYPE is a derived type containing problem structure
  !
  ! Parameters:
  ! NAME = name of problem
  ! MDL = model variables
  ! CNV = convolution variables
  ! FLD = fields variables
  ! KRN = kernel variables
  ! FRI = friction variables
  ! MSH = mesh variables
  ! FRT = front variables

  type problem_type
     character(64) :: name
     type(model_type),pointer :: mdl=>null()
     type(convolution_type),pointer :: cnv=>null()
     type(fields_type),pointer :: fld=>null()
     type(kernel_type),pointer :: krn=>null()
     type(friction_type),pointer :: fri=>null()
     type(mesh_type),pointer :: msh=>null()
     type(front_type),pointer :: frt=>null()
  end type problem_type

contains


  subroutine init_problem(pb,name)
    ! INIT_PROBLEM initializes a problem
    ! 
    ! Modified: 9 July 2007

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! NAME = name of problem

    type(problem_type),intent(inout) :: pb
    character(*),intent(in) :: name

    pb%name = name

    ! allocate derived types used by slaves
    
    allocate(pb%mdl)
    allocate(pb%cnv)
    allocate(pb%krn)
    allocate(pb%fld)
    allocate(pb%fri)
    allocate(pb%msh)
    allocate(pb%frt)
    
  end subroutine init_problem


  subroutine destroy_problem(pb)
    ! DESTROY_PROBLEM destroys all problem objects
    ! 
    ! Modified: 9 July 2007

    use model, only : destroy_model
    use convolution, only : destroy_convolution
    use fields, only : destroy_fields
    use fft_routines, only : destroy_fft
    use kernel, only : destroy_kernel
    use friction, only : destroy_friction
    use mesh, only : destroy_mesh
    use front, only : destroy_front

    implicit none

    ! I/O Parameters:
    ! PB = problem structure

    type(problem_type),intent(inout) :: pb

    ! destroy constituent objects

    if (associated(pb%mdl)) then
       call destroy_model(pb%mdl)
       deallocate(pb%mdl)
    end if
    pb%mdl => null()

    if (associated(pb%cnv)) then
       call destroy_convolution(pb%cnv)
       deallocate(pb%cnv)
    end if
    pb%cnv => null()

    if (associated(pb%fld)) then
       call destroy_fields(pb%fld)
       deallocate(pb%fld)
    end if
    pb%fld => null()

    call destroy_fft

    if (associated(pb%krn)) then
       call destroy_kernel(pb%krn)
       deallocate(pb%krn)
    end if
    pb%krn => null()

    if (associated(pb%fri)) then
       call destroy_friction(pb%fri)
       deallocate(pb%fri)
    end if
    pb%fri => null()

    if (associated(pb%msh)) then
       call destroy_mesh(pb%msh)
       deallocate(pb%msh)
    end if
    pb%msh => null()

    if (associated(pb%frt)) then
       call destroy_front(pb%frt)
       deallocate(pb%frt)
    end if
    pb%frt => null()

  end subroutine destroy_problem


end module problem
