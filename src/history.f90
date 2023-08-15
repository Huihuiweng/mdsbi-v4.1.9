! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module history

  ! HISTORY constains variables for storing history of Fourier coefficients
  ! of slip, slip velocity, displacement, or particle velocity along with
  ! precalculated convolution kernels
  ! 
  ! Modified: 6 August 2007

  use constants, only : pr,pin
  implicit none

  ! HISTORY_TYPE is a derived type containing history of Fourier coefficients for one mode
  !
  ! Parameters:
  ! N = number of time steps saved for each Fourier mode
  ! HX = history of Fourier coefficients of x component of slip (Ux) or slip velocity (Vx)
  ! HY = history of Fourier coefficients of y component of slip (Uy) or slip velocity (Vy)
  ! HXP = history of Fourier coefficients of x component of displacement (uxp) or particle velocity (vxp), plus side
  ! HXM = history of Fourier coefficients of x component of displacement (uxm) or particle velocity (vxm), minus side
  ! HYP = history of Fourier coefficients of y component of displacement (uyp) or particle velocity (vyp), plus side
  ! HYM = history of Fourier coefficients of y component of displacement (uym) or particle velocity (vym), minus side
  ! HZP = history of Fourier coefficients of z component of displacement (uzp) or particle velocity (vzp), plus side
  ! HZM = history of Fourier coefficients of z component of displacement (uzm) or particle velocity (vzm), minus side
  ! Note: This implementation allows different modes to have different numbers of time steps saved in their history. 
  ! The current time step is stored as the first element; ie, Hx(0) and Hy(0) and a shift is applied each 
  ! time step so that the current time step always holds the first spot.
  ! C = precalculated convolution kernel (first index is time step, second index is mode,
  ! third indicate plus or minus side for bimaterial problem)

  type history_type
     integer(pin) :: n
     complex(pr),dimension(:),pointer :: Hx=>null(),Hy=>null(), &
          Hxp=>null(),Hxm=>null(),Hyp=>null(),Hym=>null(),Hzp=>null(),Hzm=>null()
     real(pr),dimension(:,:,:),allocatable :: C
  end type history_type


contains


  subroutine create_history(mdl,his,n,precalc)
    ! CREATE_HISTORY creates derived type his
    ! 
    ! Modified: 6 August 2007

    use constants, only : zeroc
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! HIS = history variables
    ! N = number of time steps at which history is saved
    ! PRECALC = precalculate convolution kernels or not

    type(model_type),intent(in) :: mdl
    type(history_type),intent(inout) :: his
    integer(pin),intent(in) :: n
    logical,intent(in) :: precalc

    ! assign input variables to components of derived type

    his%n = n

    ! allocate memory for history arrays and initialize

    if (mdl%bm) then

       allocate(his%Hxp(0:n))
       his%Hxp = zeroc
       allocate(his%Hxm(0:n))
       his%Hxm = zeroc

       allocate(his%Hyp(0:n))
       his%Hyp = zeroc

       if (mdl%transverse_slip) then
          allocate(his%Hym(0:n))
          his%Hym = zeroc
       else
          his%Hym => his%Hyp
       end if

       allocate(his%Hzp(0:n))
       his%Hzp = zeroc

       if (mdl%opening) then
          allocate(his%Hzm(0:n))
          his%Hzm = zeroc
       else
          his%Hzm => his%Hzp
       end if

       if (precalc) allocate(his%C(0:n,4,2))

    else

       allocate(his%Hx(0:n))
       his%Hx = zeroc
       if (mdl%transverse_slip) then
          allocate(his%Hy(0:n))
          his%Hy = zeroc
       end if

       if (precalc) then
          select case(mdl%dim)
          case default
             select case(mdl%mode)
             case(1) ! mode I only
                allocate(his%C(0:n,1:1,1))
             case(2) ! mode II only
                allocate(his%C(0:n,2:2,1))
             case(3) ! mode III only
                allocate(his%C(0:n,3:3,1))
             case default ! mixed mode (II and III)
                allocate(his%C(0:n,2:3,1))
             end select
          case(3)
             allocate(his%C(0:n,2:3,1))
          end select
       end if

    end if

  end subroutine create_history


  subroutine destroy_history(his)
    ! DESTROY_HISTORY destroys derived type his
    ! 
    ! Modified: 20 January 2007

    implicit none

    ! I/O Parameters:
    ! HIS = history variables

    type(history_type),intent(inout) :: his

    ! deallocate memory assigned to allocatable arrays

    if (associated(his%Hx)) deallocate(his%Hx)
    his%Hx => null()
    if (associated(his%Hy)) deallocate(his%Hy)
    his%Hy => null()

    if (associated(his%Hxp)) deallocate(his%Hxp)
    his%Hxp => null()
    if (associated(his%Hxm)) deallocate(his%Hxm)
    his%Hxm => null()

    if (associated(his%Hym).and.(.not.associated(his%Hym,his%Hyp))) deallocate(his%Hym)
    his%Hym => null()
    if (associated(his%Hyp)) deallocate(his%Hyp)
    his%Hyp => null()

    if (associated(his%Hzm).and.(.not.associated(his%Hzm,his%Hzp))) deallocate(his%Hzm)
    his%Hzm => null()
    if (associated(his%Hzp)) deallocate(his%Hzp)
    his%Hzp => null()

    if (allocated(his%C)) deallocate(his%C)

  end subroutine destroy_history


end module history
