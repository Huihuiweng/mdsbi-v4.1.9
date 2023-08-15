! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module substep
  ! SUBSTEP contains routines for integrating over several small substeps within
  ! an elastodynamic time step using an interpolated value of stress transfer
  ! 
  ! Modified: 22 September 2007

  use constants, only : pr,pin

  implicit none


contains


  subroutine integrate_substep(pb,dte,dt,accept_step)
    ! INTEGRATE_SUBSTEP integrates diffusion and state-evolution equations over
    ! multiple substeps within one elastodynamic time step
    ! 
    ! Modified: 21 September 2007

    use problem, only : problem_type
    use integration, only : move_fields

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! DTE = elastodynamic time step
    ! DT = substep length
    ! ACCEPT_STEP = accept integration step

    type(problem_type),intent(inout) :: pb
    real(pr),intent(in) :: dte
    real(pr),intent(inout) :: dt
    logical,intent(out) :: accept_step

    ! transfer values from stage=1 to stage=2 (to start integration)
       
    call move_fields(pb%fld%flt,pb%mdl,n1=2,n2=1)
        
    ! store or reset error estimate and time step
    
    if (pb%mdl%iter==0) then

       ! store values from previous time step

       pb%fld%flt%aerr_old = pb%fld%flt%aerr
       pb%fld%flt%rerr_old = pb%fld%flt%rerr
       pb%mdl%dt_old = dt

    else

       ! reset values to those from previous time step
       ! (ensures that each iteration has same initial conditions)

       pb%fld%flt%aerr = pb%fld%flt%aerr_old
       pb%fld%flt%rerr = pb%fld%flt%rerr_old
       dt = pb%mdl%dt_old

    end if

    ! integrate over substeps

    select case(pb%mdl%substep_method)
    case('iterative')
       call substep_iterative(pb,dte,dt,accept_step)
    case('rk')
       call substep_rk(pb,dte,dt,accept_step)
    case('adaptive')
       call substep_adaptive(pb,dte,dt,accept_step)
    end select

  end subroutine integrate_substep


  subroutine substep_iterative(pb,dte,dt,accept_step)
    ! SUBSTEP_ITERATIVE integrates diffusion and state-evolution equations over constant substeps
    ! within one elastodynamic time step using iterative predictor-corrector method
    ! 
    ! Modified: 20 September 2007

    use problem, only : problem_type
    use integration, only : move_fields
    use time_step, only : step_iterative

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! DTE = elastodynamic time step
    ! DT = substep time step
    ! ACCEPT_STEP = accept integration step

    type(problem_type),intent(inout) :: pb
    real(pr),intent(in) :: dte
    real(pr),intent(inout) :: dt
    logical,intent(out) :: accept_step

    ! Internal variables:
    ! T1 = time at beginning of elastodynamic time step
    ! T2 = time at end of elastodynamic time step
    ! ITER = iteration index
    ! N = index of substep
    ! NT = number of substeps in elastodynamic time step

    real(pr) :: t1,t2
    integer(pin) :: iter,n,nt

    ! fix number of substeps

    nt = pb%mdl%substep_nt

    ! calculate substep length
    
    dt = dte/real(nt,pr)
    
    ! store initial and final time

    t1 = pb%mdl%t
    t2 = pb%mdl%t+dte
       
    ! loop over substeps
    
    do n = 1,nt
       
       ! take explicit step

       call step_iterative(pb,dt,s0=2,stress_transfer='interpolate', &
            first_step=.true.,accept_step=accept_step,t1=t1,t2=t2)

       ! if unacceptable, reject this step and exit routine
    
       if (.not.accept_step) return

       ! iterate to improve solution

       do iter = 1,pb%mdl%substep_niter

          ! repeat step using trapezoid rule
          
          call step_iterative(pb,dt,s0=2,stress_transfer='interpolate', &
               first_step=.false.,accept_step=accept_step,t1=t1,t2=t2)
          
          ! if unacceptable, reject this step and exit routine
          
          if (.not.accept_step) return

       end do

       ! move field values from end of current time step (stage=3) to start of next (stage=2)
       
       call move_fields(pb%fld%flt,pb%mdl,n1=2,n2=3)
          
    end do

  end subroutine substep_iterative


  subroutine substep_rk(pb,dte,dt,accept_step)
    ! SUBSTEP_RK integrates diffusion and state-evolution equations over
    ! constant substeps within one elastodynamic time step (but using Runge-Kutta methods)
    ! 
    ! Modified: 21 September 2007

    use problem, only : problem_type
    use integration, only : move_fields
    use time_step, only : step_rk,step_rk_imex
    use mpi_routines, only : is_master
    use constants, only : zero

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! DTE = elastodynamic time step
    ! DT = substep length
    ! ACCEPT_STEP = accept integration step

    type(problem_type),intent(inout) :: pb
    real(pr),intent(in) :: dte
    real(pr),intent(inout) :: dt
    logical,intent(out) :: accept_step

    ! Internal variables:
    ! T0 = initial time (at beginning of substep)
    ! T1 = time at beginning of elastodynamic time step
    ! T2 = time at end of elastodynamic time step
    ! DTO = time step (original)
    ! N = index of substep
    ! NT = number of substeps in elastodynamic time step
    ! ATTEMPTS = number of attempts to take acceptable time step

    real(pr) :: t0,t1,t2,dto
    integer(pin) :: n,nt,attempts

    ! fix number of substeps

    nt = pb%mdl%substep_nt

    ! calculate substep length
    
    dto = dte/real(nt,pr)
    
    ! store initial and final time

    t1 = pb%mdl%t
    t2 = pb%mdl%t+dte
       
    !if (is_master.and.pb%mdl%substep_info) &
    !     write(6,'(a,i6,a,i1,a,e20.10,a)') '----Starting substepping, n = ',pb%mdl%n,', iter = ',pb%mdl%iter,', dt = ',dto,'----'

    ! loop over substeps
    
    do n = 1,nt

       ! only one attempt

       attempts = 1

       ! store initial time
       
       t0 = pb%mdl%t
       
       ! set time step to specified value (rather than predicted value from error estimate)

       dt = dto
    
       ! Runge-Kutta step
       
       select case(pb%mdl%friction_method)
       case('strength')
          call step_rk(pb,dt,s0=2,stress_transfer='interpolate', &
               accept_step=accept_step,attempts=attempts,t0=t0,t1=t1,t2=t2)
       case('strength_rate')
          call step_rk_imex(pb,dt,s0=2,stress_transfer='interpolate', &
               accept_step=accept_step,attempts=attempts,t0=t0,t1=t1,t2=t2)
       end select
       
       if (is_master.and.pb%mdl%substep_info) &
            write(6,'(a,i3,a,e20.10,a,e20.10)') 'step ',n,', err = ',pb%fld%flt%rerr(1),', dt = ',dt

       ! if unacceptable, reject this step and exit routine
    
       if (.not.accept_step) return

       ! move field values from end of current time step [stage=rk%s+1=ns] to start of next [stage=2]
       
       call move_fields(pb%fld%flt,pb%mdl,n1=2,n2=pb%mdl%ns)
          
    end do

    if (is_master.and.pb%mdl%substep_info) &
         write(6,'(a,i3,a)') '----total number of substeps = ',nt,'----'

  end subroutine substep_rk


  subroutine substep_adaptive(pb,dte,dt,accept_step)
    ! SUBSTEP_ADAPTIVE integrates diffusion and state-evolution equations over
    ! adaptively chosen substeps within one elastodynamic time step
    ! 
    ! Modified: 22 September 2007

    use constants, only : zero
    use problem, only : problem_type
    use integration, only : move_fields,move_error
    use time_step, only : step_rk,step_rk_imex
    use mpi_routines, only : is_master

    implicit none

    ! I/O Parameters:
    ! PB = problem variables
    ! DTE = elastodynamic time step
    ! DT = substep length
    ! ACCEPT_STEP = accept integration step

    type(problem_type),intent(inout) :: pb
    real(pr),intent(in) :: dte
    real(pr),intent(inout) :: dt
    logical,intent(out) :: accept_step

    ! Internal variables:
    ! T0 = initial time (at beginning of substep)
    ! T1 = time at beginning of elastodynamic time step
    ! T2 = time at end of elastodynamic time step
    ! NT = number of substeps
    ! ATTEMPTS = number of attempts to take acceptable time step
    ! REJECTIONS = number of step rejections

    real(pr) :: t0,t1,t2
    integer(pin) :: nt,attempts,rejections
    character(16) :: str1,str2,str3,str4

    ! store initial and final time

    t1 = pb%mdl%t
    t2 = pb%mdl%t+dte

    ! set substep counter and substep rejection counter to zero

    nt = 0
    rejections = 0

    !if (is_master.and.pb%mdl%substep_info) &
    !     write(6,'(a,i6,a,i1,a)') '----Starting substepping, n = ',pb%mdl%n,', iter = ',pb%mdl%iter,'----'

    ! take multiple substeps
    
    do
       
       ! update substep counter

       nt = nt+1

       ! store initial time
       
       t0 = pb%mdl%t
    
       ! iterate over step sizes until acceptable

       attempts = 0

       do 

          ! update step size attempt counter
          
          attempts = attempts+1
                 
          ! start time from beginning of substep
          
          pb%mdl%t = t0

          ! prevent overshooting t2 by reducing time step to stop at exactly at t2 if dt too large
          
          if (pb%mdl%t+dt>t2) dt = t2-pb%mdl%t

          ! advance one time step, estimating error to predict new time step
       
          select case(pb%mdl%friction_method)
          case('strength')
             call step_rk(pb,dt,s0=2,stress_transfer='interpolate', &
                  accept_step=accept_step,attempts=attempts,t0=t0,t1=t1,t2=t2)
          case('strength_rate')
             call step_rk_imex(pb,dt,s0=2,stress_transfer='interpolate', &
                  accept_step=accept_step,attempts=attempts,t0=t0,t1=t1,t2=t2)
          end select

          ! exit loop searching for time step if step is acceptable
          
          if (accept_step) exit

          rejections = rejections+1

          !if (is_master) &
          !     write(6,'(a,i3,a,i2,a,e20.10,a,e20.10)') 'reject step ',nt,'(',attempts,'), err = ',pb%fld%flt%rerr(1),', dt = ',dt
          
       end do

       ! accept step

       !if (is_master) &
       !     write(6,'(a,i3,a,i2,a,e20.10,a,e20.10)') 'accept step ',nt,'(',attempts,'), err = ',pb%fld%flt%rerr(1),', dt = ',dt

       ! move error estimate
       
       call move_error(pb%fld%flt)
       
       ! move field values from end of current time step [stage=rk%s+1=ns] to start of next [stage=2]
       
       call move_fields(pb%fld%flt,pb%mdl,n1=2,n2=pb%mdl%ns)
          
       ! exit when end of elastodynamic time step is reached
       
       if (pb%mdl%t>=t2) exit

    end do

    if (is_master.and.pb%mdl%substep_info) then
       write(str1,'(i0)') pb%mdl%n
       write(str2,'(i0)') pb%mdl%iter
       write(str3,'(i0)') nt
       write(str4,'(i0)') rejections
       write(6,'(a,e20.10)') '---Substep info: step(iter) = ' // trim(str1) // &
            '(' // trim(str2) // ')   substeps(rejections) = ' // &
            trim(str3) // '(' // trim(str4) // ')   err = ',pb%fld%flt%rerr(0)
       !write(6,'(a,i6,a,i1,a)') '----Starting substepping, n = ',pb%mdl%n,', iter = ',pb%mdl%iter,'----'
       !write(6,'(a,i3,a,i3,a,e20.10,a)') '----total number of substeps(rejections) = ',nt, &
       ! '(',rejections,'), err = ',pb%fld%flt%rerr(0),'----'
    end if

  end subroutine substep_adaptive


end module substep

