! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module time_step

  ! TIME_STEP contains routines to advance problem through one time step
  ! 
  ! Modified: 3 December 2007

  use constants, only : pr,pin

  implicit none


contains


  subroutine predict_stress_transfer(pb,dt)
    ! PREDICT_STRESS_TRANSFER explicitly integrates V to predict U, then calculates stress transfer f
    ! 
    ! Modified: 25 August 2007

    use problem, only : problem_type
    use integration, only : integrate_rate_eqns,move_rates
    use rates, only : update_stress_transfer
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! DT = time step

    type(problem_type),intent(inout) :: pb
    real(pr),intent(in) :: dt

    select case(pb%mdl%predict_f_method)

    case default

       call error('Invalid method to predict stress transfer','predict_stress_transfer')

    case('integrate')
       
       ! predict rates at t+dt (stage=2) are those at t (stage=1)
       
       call move_rates(pb%fld%flt,pb%mdl,n1=2,n2=1)
       
       ! integrate to get U(t+dt)
       
       call integrate_rate_eqns(pb%mdl,pb%fld%flt,'U',dt,'explicit',s0=1)
       
       ! FFT necessary fields and inverse FFT to get f(t+dt)
       
       call update_stress_transfer(pb,stage=2,update_history=.true.)

    case('old')

       ! set f(t+dt) = f(t)

    case('extrapolate')

       ! interpolate stress transfer f(t+dt) from f(t-dt) and f(t)

    end select

  end subroutine predict_stress_transfer


  subroutine step_iterative(pb,dt,s0,stress_transfer,first_step,accept_step,update_history,t1,t2)
    ! STEP_ITERATIVE advances problem through one time step by basic integration (explicit or trapezoid)
    ! 
    ! Modified: 23 August 2007

    use problem, only : problem_type
    use integration, only : integrate_rate_eqns,move_rates
    use rates, only : set_rates
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! DT = time step
    ! S0 = stage at beginning of time step
    ! STRESS_TRANSFER = method used to calculate stress transfer
    !      convolution = perform convolution
    !      interpolation = interpolate between known values
    ! FIRST_STEP = flag indicating if this is first time step, in which case additional tasks are performed
    ! ACCEPT_STEP = accept integration step
    ! UPDATE_HISTORY = update history for convolution or overwrite
    ! T1 = time at beginning of elastodynamic step
    ! T2 = time at end of elastodynamic step

    type(problem_type),intent(inout) :: pb
    real(pr),intent(in) :: dt
    integer(pin),intent(in) :: s0
    character(*),intent(in) :: stress_transfer
    logical,intent(in) :: first_step
    logical,intent(out) :: accept_step
    logical,intent(in),optional :: update_history
    real(pr),intent(in),optional :: t1,t2

    ! input checks

    select case(stress_transfer)
    case('convolution')
       if (.not.(present(update_history))) &
            call error('update_history must be present for convolution','step_iterative')
    case('interpolate')
       if (.not.(present(t1).and.present(t2))) &
            call error('t1 and t2 must be present for interpolation','step_iterative')
    end select

    ! check if start of new time step

    if (first_step) then

       ! advance current time from t to t+dt
       
       pb%mdl%t = pb%mdl%t+dt
       
       ! predict rates at t+dt are those at t
       
       call move_rates(pb%fld%flt,pb%mdl,s0+1,s0)
    
    end if

    ! integrate rate equations

    if (first_step) then
       call integrate_rate_eqns(pb%mdl,pb%fld%flt,'all',dt,'explicit',s0)
    else
       call integrate_rate_eqns(pb%mdl,pb%fld%flt,'all',dt,'trapezoid',s0)
    end if

    ! set rates at t+dt

    select case(stress_transfer)
    case('convolution')
       call set_rates(pb,s0+1,stress_transfer,accept_step,update_history)
    case('interpolate')
       call set_rates(pb,s0+1,stress_transfer,accept_step,t1=t1,t2=t2)
    end select

  end subroutine step_iterative


  subroutine step_rk(pb,dt,s0,stress_transfer,accept_step,attempts,t0,t1,t2)
    ! STEP_RK advances problem through one time step using explicit Runge-Kutta
    ! time stepping with embedded error estimation formulas
    ! 
    ! Modified: 21 September 2007
    
    use problem, only : problem_type
    use integration, only : integrate_rate_eqns
    use rates, only : set_rates
    use io, only : error,warning
    use mpi_routines, only : is_master
    use mpi

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! DT = time step
    ! S0 = stage at beginning of time step
    ! STRESS_TRANSFER = method used to calculate stress transfer
    !      convolution = perform convolution
    !      interpolation = interpolate between known values
    ! ACCEPT_STEP = accept integration step
    ! ATTEMPTS = number of attempts to take acceptable time step
    ! T0 = initial time
    ! T1 = time at beginning of elastodynamic time step
    ! T2 = time at end of elastodynamic time step

    type(problem_type),intent(inout) :: pb
    real(pr),intent(inout) :: dt
    integer(pin),intent(in) :: s0,attempts
    character(*),intent(in) :: stress_transfer
    logical,intent(out) :: accept_step
    real(pr),intent(in) :: t0
    real(pr),intent(in),optional :: t1,t2

    ! Internal Parameters:
    ! STR = string used in output
    ! M = Runge-Kutta stage
    ! STAGE = current integration stage (for storage)
    ! FAILED = flag indicating if error exceeds tolerance so time step failed

    character(64) :: str
    integer(pin) :: m,stage
    logical :: failed,accept_step2

    ! input checks

    select case(stress_transfer)
    case('interpolate')
       if (.not.(present(t1).and.present(t2))) &
            call error('t1 and t2 must be present for interpolation','step_rk')
    end select

    ! set number of stages

    ! ---stage 1---

    ! keep values from end of previous time step

    ! ---stages m=2,...,s---

    do m = 2,pb%mdl%rk%s

       ! set storage stage

       stage = s0-1+m

       ! advance current time
          
       pb%mdl%t = t0+pb%mdl%rk%c(m)*dt

       ! integrate rate equations
       
       call integrate_rate_eqns(pb%mdl,pb%fld%flt,'all',dt,'rk_stage',s0,m)
       
       ! set rates at new time

       select case(stress_transfer)

       case('convolution')
          call set_rates(pb,stage,stress_transfer,accept_step,update_history=.true.)
       case('interpolate')
          call set_rates(pb,stage,stress_transfer,accept_step,t1=t1,t2=t2)
       end select

       ! if unacceptable, reject this step and exit routine
       
       if (.not.accept_step) then
          call set_time_step(pb,controller='min',dt=dt)
          return
       end if

    end do

    ! ---error estimate---
    
    call integrate_rate_eqns(pb%mdl,pb%fld%flt,'all',dt,'rk_error',s0)

    ! reject step if error greater than tolerance

    failed = (pb%fld%flt%rerr(1)>pb%mdl%rtol)

    if (failed) then
          
       if (attempts==10) then
             
          ! accept this step, even though error exceeds tolerance

          if (is_master) then
             write(str,'(a,i14)') 'Accepting unacceptable time step: ',pb%mdl%n
             call warning(str,'step_rk')
             write(str,'(a,e14.6,a,e14.6)') 'relative error = ',pb%fld%flt%rerr(1),', tolerance = ',pb%mdl%rtol
             call warning(str)
          end if

          accept_step = .true.

       else

          ! reject unless using fixed steplength RK method

          if (pb%mdl%method=='substep'.and.pb%mdl%substep_method=='rk') then
             accept_step = .true.
          else
             accept_step = .false.
          end if

       end if

    end if

    ! ---use higher-order formula for final values---

    ! set storage stage

    stage = pb%mdl%ns

    ! advance current time
       
    pb%mdl%t = t0+dt

    call integrate_rate_eqns(pb%mdl,pb%fld%flt,'all',dt,'rk_final',s0)

    ! set final rates

    select case(stress_transfer)
    case('convolution')
       call set_rates(pb,stage,stress_transfer,accept_step2,update_history=.true.)
    case('interpolate')
       call set_rates(pb,stage,stress_transfer,accept_step2,t1=t1,t2=t2)
    end select
    if (.not.accept_step2) accept_step = .false.

    ! if unacceptable, reject this step and exit routine

    if (.not.accept_step) then
       call set_time_step(pb,controller='min',dt=dt)
       return
    end if
    
    ! predict next time step
              
    if (failed) then

       if (attempts>=8) then
          call set_time_step(pb,controller='min',dt=dt)
       else
          call set_time_step(pb,controller='PID',dt=dt)
       end if

    else

       call set_time_step(pb,controller='PID',dt=dt)

    end if
    
  end subroutine step_rk


  subroutine step_rk_imex(pb,dt,s0,stress_transfer,accept_step,attempts,t0,t1,t2)
    ! STEP_RK_IMEX advances problem through one time step using implicit-explicit 
    ! Runge-Kutta time stepping with embedded error estimation formulas
    ! 
    ! Modified: 3 December 2007
    
    use problem, only : problem_type
    use integration, only : integrate_rate_eqns,predict_rk
    use load, only : set_load
    use rates, only : update_stress_transfer,set_interp_method
    use friction_routines, only : solve_friction
    use fourier, only : fft_current_U
    use fault, only : interp_stress_transfer
    use thermpres, only : thermpres_rate
    use io, only : error,warning
    use mpi_routines, only : is_master
    use mpi

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! DT = time step
    ! S0 = stage at beginning of time step
    ! STRESS_TRANSFER = method used to calculate stress transfer
    !      convolution = perform convolution
    !      interpolation = interpolate between known values
    ! ACCEPT_STEP = accept integration step
    ! ATTEMPTS = number of attempts to take acceptable time step
    ! T0 = initial time
    ! T1 = time at beginning of elastodynamic time step
    ! T2 = time at end of elastodynamic time step

    type(problem_type),intent(inout) :: pb
    real(pr),intent(inout) :: dt
    integer(pin),intent(in) :: s0,attempts
    character(*),intent(in) :: stress_transfer
    logical,intent(out) :: accept_step
    real(pr),intent(in) :: t0
    real(pr),intent(in),optional :: t1,t2

    ! Internal Parameters:
    ! METHOD = interpolation method
    ! STR = string used in output
    ! M = Runge-Kutta stage
    ! STAGE = current integration stage (for storage)
    ! FAILED = flag indicating if error exceeds tolerance so time step failed
    ! MY_ACCEPT_STEP = accept integration step (local to each process)
    ! IERR = MPI error flag  

    character(64) :: method,str
    integer(pin) :: m,stage,ierr
    logical :: failed,my_accept_step

    ! input checks

    select case(stress_transfer)
    case('interpolate')
       if (.not.(present(t1).and.present(t2))) &
            call error('t1 and t2 must be present for interpolation','step_rk')
    end select

    ! set number of stages

    ! ---stage 1---

    ! keep values from end of previous time step

    ! ---stages m=2,...,s---

    do m = 2,pb%mdl%rk%s

       ! set storage stage

       stage = s0-1+m

       ! advance current time
          
       pb%mdl%t = t0+pb%mdl%rk%c(m)*dt

       ! integrate explicit rate equations
       
       call integrate_rate_eqns(pb%mdl,pb%fld%flt,'U',dt,'rk_stage',s0,m)
       call integrate_rate_eqns(pb%mdl,pb%fld%flt,'Tp',dt,'rk_stage',s0,m)
       
       ! set time-dependent loads s0(t)
    
       call set_load(pb%mdl,pb%fld%ld,pb%fld%flt)
    
       ! set stress transfer f(t)
       
       select case(stress_transfer)
          
       case('convolution')
          
          ! FFT necessary fields and inverse FFT to get f(t)
          
          call update_stress_transfer(pb,stage,update_history=.true.)
          
       case ('interpolate')
          
          ! FFT displacement/slip U(t), but only if needed (for poroelastic effect)
          
          if (pb%fld%flt%poroelastic) call fft_current_U(pb%mdl,pb%fld%flt,pb%cnv,stage)
          
          ! interpolate stress transfer f(t) after setting interpolation method
          
          call set_interp_method(pb,method)
          
          call interp_stress_transfer(pb%fld%flt,pb%mdl,t1,t2,pb%mdl%t,method)
          
       end select
    
       ! predict stage values to speed up Newton's method
       
       call predict_rk(pb%mdl,pb%fld%flt,dt,s0,stage)
       
       ! solve elasticity equation coupled with friction for V(t) and s(t) at each point 
       ! using U(t), f(t), Q(t), p(t), T(t)
       
       call solve_friction(pb%mdl,pb%cnv,pb%fld%flt,pb%fri,stage,my_accept_step,dt,s0)

       call MPI_AllReduce(my_accept_step,accept_step,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)

       ! if unacceptable, reject this step and exit routine
       
       if (.not.accept_step) then
          call set_time_step(pb,controller='min',dt=dt)
          return
       end if

       ! calculate rates dT(t) and dp(t) in thermal pressurization equations
       
       call thermpres_rate(pb%mdl,pb%fld%flt,pb%fld%tp,stage)
    
    end do

    ! ---error estimate---
    
    call integrate_rate_eqns(pb%mdl,pb%fld%flt,'all',dt,'rk_error',s0)

    ! reject step if error greater than tolerance

    failed = (pb%fld%flt%rerr(1)>pb%mdl%rtol)

    if (failed) then
       
       if (attempts==10) then
          
          ! accept this step, even though error exceeds tolerance

          if (is_master) then
             write(str,'(a,i14)') 'Accepting unacceptable time step: ',pb%mdl%n
             call warning(str,'step_rk_imex')
             write(str,'(a,e14.6,a,e14.6)') 'relative error = ',pb%fld%flt%rerr(1),', tolerance = ',pb%mdl%rtol
             call warning(str)
          end if
          
          accept_step = .true.
          
       else

          ! reject unless using fixed steplength RK method

          if (pb%mdl%method=='substep'.and.pb%mdl%substep_method=='rk') then
             accept_step = .true.
          else
             accept_step = .false.
          end if

       end if

    end if

    ! ---use higher-order formula for final values---

    ! set storage stage

    stage = pb%mdl%ns

    ! advance current time
       
    pb%mdl%t = t0+dt

    ! set final values

    call integrate_rate_eqns(pb%mdl,pb%fld%flt,'U',dt,'rk_final',s0)
    call integrate_rate_eqns(pb%mdl,pb%fld%flt,'Tp',dt,'rk_final',s0)

    ! predict next time step
       
    if (failed) then

       if (attempts>=8) then
          call set_time_step(pb,controller='min',dt=dt)
       else
          call set_time_step(pb,controller='PID',dt=dt)
       end if

    else

       call set_time_step(pb,controller='PID',dt=dt)

    end if
    
  end subroutine step_rk_imex


  subroutine set_time_step(pb,controller,dt)
    ! SET_TIME_STEP sets adaptive time step
    ! 
    ! Modified: 28 May 2007

    use constants, only : zero,one
    use problem, only : problem_type
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! CONTROLLER = step-size controller
    ! DT = time step

    type(problem_type),intent(in) :: pb
    character(*),intent(in) :: controller
    real(pr),intent(inout) :: dt

    ! Internal Parameters:
    ! CONTROL = actual step-size controller used
    ! ERR = error
    ! TOL = tolerance
    ! KAPPA = safety factor
    ! ALPHA = exponent in step-size formula
    ! BETA = exponent in step-size formula
    ! GAMMA = exponent in step-size formula
    ! STEP_RATIO = step-size ratio (new step-size/current step-size)
    ! STEP_RATIO_MIN = minimum step-size ratio
    ! STEP_RATIO_MAX = maximum step-size ratio
    ! ERROR_TYPE = which error (relative or absolute) to control

    character(64) :: control
    real(pr) :: err(-1:1),tol,step_ratio,alpha,beta,gamma
    real(pr),parameter :: step_ratio_min=0.01_pr, step_ratio_max=10._pr, kappa=0.9_pr
    character(64),parameter :: error_type = 'relative'

    ! index of error corresponds to following times:
    ! time:  t-2*dt t-dt t
    ! index:   -1     0  1

    ! set error and tolerance
    
    select case(error_type)
    case('relative')
       err = pb%fld%flt%rerr
       tol = pb%mdl%rtol
    case('absolute')
       err = pb%fld%flt%aerr
       tol = pb%mdl%atol
    end select
    
    ! set default controller to input value
    
    control = controller
    
    ! switch controller for special cases that might cause problems
    
    if (controller/='min') then
       if (err(1)==zero) control = 'max'
       if (err(0)==zero.and.(control=='PID'.or.control=='PI')) control = 'I'
       if (err(-1)==zero.and.(control=='PID')) control = 'PI'
    end if
    
    ! set time step with controller
    
    select case(control)
       
    case default
       
       call error('invalid step-size controller','set_time_step')
       
    case ('min') ! minimum step-size ratio
       
       step_ratio = step_ratio_min
       
    case ('max') ! maximum step-size ratio
       
       step_ratio = step_ratio_max
       
    case('I') ! integral feedback (I) controller
       
       alpha = one/(one+pb%mdl%rk%p)
       step_ratio = kappa*(tol/err(1))**alpha
       
    case('PI') ! proportional+integral feedback (PI) controller
       
       alpha = 0.7_pr/pb%mdl%rk%p
       beta = 0.4_pr/pb%mdl%rk%p
       step_ratio = kappa*(tol/err(1))**alpha*(err(0)/tol)**beta
       
    case('PID') ! proportional+integral+derivative feedback (PID) controller
       
       alpha = 0.49_pr/pb%mdl%rk%p
       beta = 0.34_pr/pb%mdl%rk%p
       gamma = 0.1_pr/pb%mdl%rk%p
       step_ratio = kappa*(tol/err(1))**alpha*(err(0)/tol)**beta*(tol/err(-1))**gamma
       
    end select
    
    ! bound step-size ratio (prevent overly small or large changes)
    
    dt = dt*min(step_ratio_max,max(step_ratio_min,step_ratio))
    
    ! bound step size
    
    dt = min(pb%mdl%dtmax,max(pb%mdl%dtmin,dt))
    
  end subroutine set_time_step


end module time_step
