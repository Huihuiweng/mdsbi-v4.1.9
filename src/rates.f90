! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module rates

  ! RATES contains routines to set the rates (velocity, etc.) in the system of
  ! functional differential equations to be solved
  ! 
  ! Modified: 10 August 2010

  use constants, only : pr,pin

  implicit none

contains


  subroutine update_stress_transfer(pb,stage,update_history)
    ! UPDATE_STRESS_TRANSFER performs convolution and calculates stress transfer f
    ! 
    ! Modified: 23 August 2007

    use problem, only : problem_type
    use fourier, only : fft_current_U,fft_history,inverse_fft_F
    use fault, only : set_stress_transfer
    use convolution_routines, only : do_convolution

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! STAGE = integration stage
    ! UPDATE_HISTORY = update history for convolution or overwrite

    type(problem_type),intent(inout) :: pb
    integer(pin),intent(in) :: stage
    logical,intent(in) :: update_history

    ! FFT displacement/slip U(t+dt)
    
    call fft_current_U(pb%mdl,pb%fld%flt,pb%cnv,stage)
    
    ! forward FFT displacement/slip U(t+dt) or velocity V(t+dt) for use in convolution routines
    
    call fft_history(pb%mdl,pb%fld%flt,pb%cnv,stage,update_history)
    
    ! perform time convolution to get F(t+dt)
    
    call do_convolution(pb%mdl,pb%krn,pb%cnv,update_history)
    
    ! inverse FFT F(t+dt) to f(t+dt)
    
    select case(pb%mdl%method)
    case default
       call inverse_fft_F(pb%mdl,pb%fld%flt,pb%cnv,stage=0)
       call set_stress_transfer(pb%mdl,pb%fld%flt,stage=0)
    case('substep')
       call inverse_fft_F(pb%mdl,pb%fld%flt,pb%cnv,stage=1)
       call set_stress_transfer(pb%mdl,pb%fld%flt,stage=1)
    end select

  end subroutine update_stress_transfer


  subroutine init_rates(pb,accept_step)
    ! INIT_RATES calculates initial rates (velocity, state-rate, rates of temperature and pressure change, etc.)
    ! and stores initial values of convolution history
    ! 
    ! Modified: 10 August 2010

    use problem, only : problem_type
    use friction_routines, only : solve_friction,state_rate,strength_rate
    use fourier, only : fft_current_U,fft_history
    use load, only : set_load
    use thermpres, only : thermpres_rate
    use io, only : error
    use mpi

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! ACCEPT_STEP = accept integration step

    type(problem_type),intent(inout) :: pb
    logical,intent(out) :: accept_step

    ! Internal Parameters:
    ! MY_ACCEPT_STEP = accept integration step (local to each slave)
    ! IERR = MPI error flag  

    logical :: my_accept_step
    integer(pin) :: ierr

    ! set time-dependent loads s0(t)
    
    call set_load(pb%mdl,pb%fld%ld,pb%fld%flt)
    
    ! FFT displacement/slip U(t)
    
    call fft_current_U(pb%mdl,pb%fld%flt,pb%cnv,stage=1)

    ! set stress transfer f(t)

    call update_stress_transfer(pb,stage=1,update_history=.false.)

    ! initialize T and p on fault

    if (pb%fld%flt%thermpres) then
       pb%fld%flt%T0(:,:,1) = pb%fld%flt%T(1,:,:,1)
       pb%fld%flt%p0(:,:,1) = pb%fld%flt%p(1,:,:,1)
    end if

    ! if solving static elasticity problem, set all rates to zero and return
       
    if (pb%mdl%method=='static') then
          
       call zero_rates(pb%mdl,pb%fld%flt)
          
       accept_step = .true.
       
       return
       
    end if

    select case(pb%mdl%friction_method)

    case default

       call error('Invalid friction_method','init_rates')

    case('strength')

       ! solve elasticity equation coupled with friction for V(t) and s(t) at each point using
       ! U(t), f(t), Q(t), p(t), T(t)
       
       call solve_friction(pb%mdl,pb%cnv,pb%fld%flt,pb%fri,s=1,accept=my_accept_step)

       ! if unacceptable, exit routine
       
       call MPI_AllReduce(my_accept_step,accept_step,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
       
       if (.not.accept_step) return
       
       ! calculate state-rate dQ(t) (if needed)
    
       call state_rate(pb%mdl,pb%fld%flt,pb%fri,s=1)
       
    case('strength_rate')

       ! solve elasticity equation coupled with friction for V(t) and s(t) at each point using
       ! U(t), f(t), Q(t), p(t), T(t)

       call solve_friction(pb%mdl,pb%cnv,pb%fld%flt,pb%fri,s=1,accept=my_accept_step,dt=pb%mdl%dt,s0=1,direct_only=.true.)

       ! if unacceptable, exit routine
       
       call MPI_AllReduce(my_accept_step,accept_step,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
       
       if (.not.accept_step) return
       
       call strength_rate(pb%mdl,pb%fld%flt,pb%fri,stage=1)

       accept_step = .true.

    end select

    ! forward FFT displacement/slip U(t) or velocity V(t) for use in convolution routines
    
    call fft_history(pb%mdl,pb%fld%flt,pb%cnv,stage=1,update_history=.true.)
    
    ! calculate rates dT(t) and dp(t) in thermal pressurization equations
    
    call thermpres_rate(pb%mdl,pb%fld%flt,pb%fld%tp,stage=1)

  end subroutine init_rates


  subroutine set_rates(pb,stage,stress_transfer,accept_step,update_history,t1,t2)
    ! SET_RATES calculates rates (velocity, state-rate, rates of temperature and pressure change, etc.)
    ! at some stage.  Stress transfer is calculated by either 1. FFT, convolution, inverse FFT, or 
    ! 2. interpolation.
    ! 
    ! Modified: 5 December 2007

    use problem, only : problem_type
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! STAGE = integration stage
    ! STRESS_TRANSFER = method used to calculate stress transfer
    !      convolution = perform convolution
    !      interpolation = interpolate between known values
    ! ACCEPT_STEP = accept integration step
    ! UPDATE_HISTORY = update history for convolution or overwrite
    ! T1 = time at beginning of elastodynamic step
    ! T2 = time at end of elastodynamic step

    type(problem_type),intent(inout) :: pb
    integer(pin),intent(in) :: stage
    character(*),intent(in) :: stress_transfer
    logical,intent(out) :: accept_step
    logical,intent(in),optional :: update_history
    real(pr),intent(in),optional :: t1,t2

    ! special case for static elasticity
    
    if (pb%mdl%method=='static') then
       call set_rates_static(pb,stage)
       accept_step = .true.
       return
    end if

    ! check inputs and call set_rates_ex with appropriate arguments

    select case(stress_transfer)
    case('convolution')
       if (.not.(present(update_history))) &
            call error('update_history must be present for convolution','set_rates')
       call set_rates_ex(pb,stage,stress_transfer,accept_step,update_history)
    case ('interpolate')
       if (.not.(present(t1).and.present(t2))) &
            call error('t1 and t2 must be present for interpolation','set_rates')
       call set_rates_ex(pb,stage,stress_transfer,accept_step,t1=t1,t2=t2)
    end select

  end subroutine set_rates


  subroutine set_rates_static(pb,stage)
    ! SET_RATES_STATIC sets stress transfer and zeros rates.  This subroutine is for static elasticity.
    ! 
    ! Modified: 23 August 2007

    use problem, only : problem_type
    use convolution_routines, only : do_convolution
    use fourier, only : fft_current_U,inverse_fft_F,set_normal_displacement
    use fault, only : set_stress_transfer
    use load, only : set_load

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! STAGE = integration stage

    type(problem_type),intent(inout) :: pb
    integer(pin),intent(in) :: stage

    ! set time-dependent loads s0(t)
    
    call set_load(pb%mdl,pb%fld%ld,pb%fld%flt)
    
    ! FFT displacement/slip U(t)
    
    call fft_current_U(pb%mdl,pb%fld%flt,pb%cnv,stage)
    
    ! set fault-normal displacements (when opening is not permitted)
    
    call set_normal_displacement(pb%mdl,pb%krn,pb%fld%flt,pb%cnv,stage)
    
    ! perform time convolution to get F(t)
    
    call do_convolution(pb%mdl,pb%krn,pb%cnv,update_history=.true.)
    
    ! inverse FFT F(t) to f(t)
    
    call inverse_fft_F(pb%mdl,pb%fld%flt,pb%cnv,stage=0)
    
    call set_stress_transfer(pb%mdl,pb%fld%flt,stage=0)
    
    ! set all rates to zero
    
    call zero_rates(pb%mdl,pb%fld%flt)
    
  end subroutine set_rates_static


  subroutine set_rates_ex(pb,stage,stress_transfer,accept_step,update_history,t1,t2)
    ! SET_RATES_EX calculates rates (velocity, state-rate, rates of temperature and pressure change, etc.)
    ! at some stage.  Stress transfer is calculated by either 1. FFT, convolution, inverse FFT, or 
    ! 2. interpolation.  This subroutine is for explicit integration.
    ! 
    ! Modified: 3 December 2007

    use problem, only : problem_type
    use fault, only : interp_stress_transfer
    use friction_routines, only : solve_friction,state_rate
    use fourier, only : fft_current_U
    use load, only : set_load
    use thermpres, only : thermpres_rate
    use io, only : error
    use mpi

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! STAGE = integration stage
    ! STRESS_TRANSFER = method used to calculate stress transfer
    !      convolution = perform convolution
    !      interpolation = interpolate between known values
    ! ACCEPT_STEP = accept integration step
    ! UPDATE_HISTORY = update history for convolution or overwrite
    ! T1 = time at beginning of elastodynamic step
    ! T2 = time at end of elastodynamic step

    type(problem_type),intent(inout) :: pb
    integer(pin),intent(in) :: stage
    character(*),intent(in) :: stress_transfer
    logical,intent(out) :: accept_step
    logical,intent(in),optional :: update_history
    real(pr),intent(in),optional :: t1,t2

    ! Internal Parameters:
    ! METHOD = interpolation method
    ! MY_ACCEPT_STEP = accept integration step (local to each process)
    ! IERR = MPI error flag  

    character(64) :: method
    integer(pin) :: ierr
    logical :: my_accept_step

    ! set time-dependent loads s0(t)
    
    call set_load(pb%mdl,pb%fld%ld,pb%fld%flt)
    
    ! set stress transfer f(t)

    select case(stress_transfer)
       
    case('convolution')

       ! FFT necessary fields and inverse FFT to get f(t)

      call update_stress_transfer(pb,stage,update_history)

    case ('interpolate')

       ! FFT displacement/slip U(t), but only if needed (for poroelastic effect)
          
       if (pb%fld%flt%poroelastic) call fft_current_U(pb%mdl,pb%fld%flt,pb%cnv,stage)
    
       ! interpolate stress transfer f(t) after setting interpolation method
              
       call set_interp_method(pb,method)

       call interp_stress_transfer(pb%fld%flt,pb%mdl,t1,t2,pb%mdl%t,method)
       
    end select

    ! solve elasticity equation coupled with friction for V(t) and s(t) at each point using
    ! U(t), f(t), Q(t), p(t), T(t), as well as fault-normal velocity vz(t) and effective normal stress sz(t)
    
    call solve_friction(pb%mdl,pb%cnv,pb%fld%flt,pb%fri,stage,my_accept_step)  

    call MPI_AllReduce(my_accept_step,accept_step,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)

    ! if unacceptable, reject this step and exit routine
    ! (usually to try again with smaller step size)
    
    if (.not.accept_step) return
    
    ! calculate state-rate dQ(t)
       
    call state_rate(pb%mdl,pb%fld%flt,pb%fri,stage)

    ! calculate rates dT(t) and dp(t) in thermal pressurization equations
    
    call thermpres_rate(pb%mdl,pb%fld%flt,pb%fld%tp,stage)
    
  end subroutine set_rates_ex


  subroutine zero_rates(mdl,flt)
    ! ZERO_RATES sets all rates to zero
    ! 
    ! Modified: 22 August 2007

    use constants, only : zero
    use model, only : model_type
    use fault, only : fault_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt

    if (mdl%bm) then
       flt%vxp = zero
       flt%vxm = zero
       flt%vyp = zero
       flt%vym = zero
       flt%vzp = zero
       flt%vzm = zero
    else
       flt%Vx = zero
       flt%Vy = zero
    end if
       
    flt%dQ = zero
    flt%dS = zero
    flt%dV = zero
       
    if (allocated(flt%dT)) flt%dT = zero
    if (allocated(flt%dp)) flt%dp = zero

  end subroutine zero_rates


  subroutine set_interp_method(pb,method)
    ! SET_INTERP_METHOD returns the interpolation method to be used based on current time step
    ! and other parameters
    ! 
    ! Modified: 22 August 2007
    
    use problem, only : problem_type

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! METHOD = interpolation method

    type(problem_type),intent(in) :: pb
    character(*),intent(out) :: method

    if (pb%mdl%iter==0) then
       
       ! explicit or predictor
       
       select case(pb%mdl%interp_method)
       case('forward')
          if(pb%mdl%n==1) then
             method = 'linear_forward'
          else
             method = 'quadratic'
          end if
       case('backward')
          if(pb%mdl%n==1) then
             method = 'constant'
          else
             method = 'linear_backward'
          end if
       case('safe1','safe2')
          method = 'constant'
       end select
       
    else
       
       ! corrector
       
       select case(pb%mdl%interp_method)
       case('forward')
          if(pb%mdl%n==1) then
             method = 'linear_forward'
          else
             method = 'quadratic'
          end if
       case('backward','safe2')
          if(pb%mdl%n==1) then
             method = 'linear_forward'
          else
             method = 'quadratic'
          end if
       case('safe1')
          method = 'linear_forward'
       end select
       
    end if
    
  end subroutine set_interp_method


end module rates
