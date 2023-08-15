! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module problem_routines

  ! PROBLEM_ROUTINES contains routines to solve a problem
  ! 
  ! Modified: 6 August 2010

  use constants, only : pr,pin

  implicit none

contains


  subroutine solve_problem(pb,ninput,start_time,ninfo)
    ! SOLVE_PROBLEM solves a problem
    ! 
    ! Modified: 6 August 2010

    use constants, only : zero
    use problem, only : problem_type
    use init, only : initialize
    use mesh, only : output_mesh
    use front, only : output_front,set_front
    use io, only : message,error
    use rates, only : init_rates
    use static, only : iterate_static
    use mpi_routines, only : is_master,clock_split

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! NINPUT = unit for input file

    type(problem_type),intent(inout) :: pb
    integer(pin),intent(inout) :: ninput
    real(pr),intent(in) :: start_time
    integer(pin),intent(in) :: ninfo

    ! Internal Parameters:
    ! DT = time step
    ! DTS = substep length
    ! STR = character string used for writing message to status file
    ! ACCEPT_STEP = accept integration step
    ! DONE = flag indicating when all time steps are complete

    real(pr) :: dt,dts,time_per_step,end_time
    character(128) :: str
    logical :: accept_step,done

    ! initialize everything at t=0
    
    call initialize(pb,ninput)

    ! solve static or dynamic problem

    select case(pb%mdl%method)
     
    case default
       
       ! set initial rates (at t=0)
       
       call init_rates(pb,accept_step)
       
       if (.not.accept_step) call error('Initial conditions rejected','solve_problem')
       
       ! print initialized statement
       
       if (is_master) call message('initialized')
          
       ! output fields on mesh
          
       call output_mesh(pb%mdl,pb%fld,pb%msh,pb%mdl%t)
          
       ! update position of rupture front
       
       call set_front(pb%frt,pb%mdl,pb%fld,pb%mdl%t,dt)
       
       ! output rupture front position
       
       call output_front(pb%frt)
              
       ! set initial time step (to elastodynamic time step)
       
       dt = pb%mdl%dt
       dts = pb%mdl%dt
                
       ! loop over time steps

       done = .false.
       
       do
          
          ! advance one time step (all processes participate)
          
          call advance_time_step(pb,dt,dts)

          ! output fields on mesh
          
          call output_mesh(pb%mdl,pb%fld,pb%msh,pb%mdl%t)          
          
          ! update position of rupture front
          
          call set_front(pb%frt,pb%mdl,pb%fld,pb%mdl%t,dt)
          
          ! output rupture front position
          
          call output_front(pb%frt)
             
          ! status update
          
          if (mod(pb%mdl%n-1,ninfo)==0) then
             call clock_split(time_per_step,end_time)
             time_per_step = time_per_step/dble(min(ninfo,pb%mdl%n))
             if (is_master) then
                write(str,'(a,i0,a,i0,a,f10.5,a,f10.5,a,f0.6,a,f0.6,a)') &
                     'time step ',pb%mdl%n,' of ',pb%mdl%nt,': t = ',pb%mdl%t,', dt = ',dt, &
                     ' (wall clock: ',end_time-start_time,' s, ',time_per_step,' s/step)'
                !if (pb%mdl%method=='adaptive') then
                !   write(str,'("simulation time: ",e20.10," out of ",e20.10," at ",i12, &
                !        " time steps")') pb%mdl%t,pb%mdl%tmax,pb%mdl%n
                !else
                !   write(str,'("time steps: ",i12," out of ",i12)') pb%mdl%n,pb%mdl%nt
                !end if
                call message(str)
             end if
          end if
  
          ! done if time reaches/exceeds maximum time
          
          if (pb%mdl%t>=pb%mdl%tmax.and.pb%mdl%tmax/=zero) done = .true.
          
          ! done if maximum number of time steps is exceeded
          
          if (pb%mdl%n>=pb%mdl%nt) done = .true.
          
          if (done) exit

       end do

    case('static')
       
       if (is_master) then

          ! print initialized statement
          
          call message('initialized')
          
          ! solve static problem
          
          call iterate_static(pb)
          
          ! output fields on mesh
          
          call output_mesh(pb%mdl,pb%fld,pb%msh,pb%mdl%t)
          
          ! update position of rupture front
          
          call set_front(pb%frt,pb%mdl,pb%fld,pb%mdl%t,dt)
          
          ! output rupture front position
          
          call output_front(pb%frt)
     
       end if

    end select

  end subroutine solve_problem


  subroutine advance_time_step(pb,dt,dts)
    ! ADVANCE_TIME_STEP advances problem through one time step
    ! 
    ! Modified: 18 September 2007

    use constants, only : zero
    use problem, only : problem_type
    use integration, only : finalize_fields,move_error
    use rates, only : set_rates,update_stress_transfer
    use thermpres, only : thermpres_rate
    use substep, only : integrate_substep
    use time_step, only : predict_stress_transfer,step_iterative,step_rk,step_rk_imex
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! DT = time step
    ! DTS = substep length

    type(problem_type),intent(inout) :: pb
    real(pr),intent(inout) :: dt,dts

    ! Internal Parameters:
    ! ITER = iteration of elastodynamic time step
    ! T0 = initial time
    ! ATTEMPTS = number of attempts to take acceptable time step
    ! ACCEPT_STEP = accept integration step

    integer(pin) :: iter,attempts
    real(pr) :: t0
    logical :: accept_step

    ! update time step index
    
    pb%mdl%n = pb%mdl%n+1
    
    ! store initial time
       
    t0 = pb%mdl%t
    
    ! iterate over step sizes until acceptable

    attempts = 0

    do

       ! update step size attempt counter

       attempts = attempts+1

       ! start time from beginning of time step
       
       pb%mdl%t = t0

       ! prevent overshooting tmax by reducing time step to stop at exactly at tmax if dt too large
       
       if (pb%mdl%tmax/=zero.and.pb%mdl%t+dt>pb%mdl%tmax) dt = pb%mdl%tmax-pb%mdl%t
       
       ! take time step

       select case(pb%mdl%method)

       case ('iterative')

          ! take explicit step

          call step_iterative(pb,dt,s0=1,stress_transfer='convolution', &
               first_step=.true.,accept_step=accept_step,update_history=.true.)

          ! repeat step using trapezoid rule some number of times

          do iter = 1,pb%mdl%niter
             
             pb%mdl%iter = iter

             call step_iterative(pb,dt,s0=1,stress_transfer='convolution', &
                  first_step=.false.,accept_step=accept_step,update_history=.false.)

             if (.not.accept_step) exit

          end do
             
          ! exit loop searching for time step if step is acceptable
          
          if (accept_step) exit
          
          ! time step failed

          call error('Time step failed, try substepping','advance_time_step')

       case ('substep')

          ! predict stress transfer f(t+dt)
          
          call predict_stress_transfer(pb,dt)

          ! integrate over substeps, iterating as desired to improve stress transfer estimate
          ! (note that first time, iter=0, must always be done; iter>=1 are corrections)

          do iter = 0,pb%mdl%niter
             
             pb%mdl%iter = iter
             
             ! revert current time back to t0
             
             pb%mdl%t = t0
             
             ! integrate elasticity+friction law over substeps (length dts) within elastodynamic 
             ! step (length dt) with interpolation of stress transfer between t and t+dt

             call integrate_substep(pb,dt,dts,accept_step)

             ! update stress transfer f(t+dt) and (if appropriate) set rates at t+dt
             
             select case(pb%mdl%friction_method)
             case default
                call set_rates(pb,stage=2,stress_transfer='convolution', &
                     accept_step=accept_step,update_history=.false.)
             case('strength_rate')
                if (iter/=pb%mdl%niter) &
                     call update_stress_transfer(pb,stage=2,update_history=.false.)
                call thermpres_rate(pb%mdl,pb%fld%flt,pb%fld%tp,stage=2)
             end select

             if (.not.accept_step) exit
             
          end do

          ! exit loop searching for time step if step is acceptable
             
          if (accept_step) exit

          call error('Time step failed','advance_time_step')
          
       case('adaptive')

          ! advance one time step, estimating error to predict new time step
       
          select case(pb%mdl%friction_method)
          case('strength')
             call step_rk(pb,dt,s0=1,stress_transfer='convolution', &
                  accept_step=accept_step,attempts=attempts,t0=t0)
          case('strength_rate')
             call step_rk_imex(pb,dt,s0=1,stress_transfer='convolution', &
                  accept_step=accept_step,attempts=attempts,t0=t0)
          end select

          ! exit loop searching for time step if step is acceptable; move error estimate
          
          if (accept_step) then
             call move_error(pb%fld%flt)
             exit
          end if

       end select

    end do

    ! accept step; overwrite fields at time t with fields at time t+dt
       
    call finalize_fields(pb%fld%flt,pb%mdl)
    
  end subroutine advance_time_step


end module problem_routines

