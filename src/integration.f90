! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module integration

  ! INTEGRATION contains routines to integrate elasticity and state-evolution equations
  ! 
  ! Modified: 25 August 2007

  use constants, only : pr,pin

  implicit none

contains


  subroutine integrate_rate_eqns(mdl,flt,field_name,dt,step,s0,s)
    ! INTEGRATE_RATE_EQNS integrates rate equations: dU/dt=V, dQ/dt=dQ,
    ! dT/dt=dT, dp/dt=dp to update displacement/slip, state, temperature, pressure
    ! 
    ! Modified: 25 August 2007

    use constants, only : zero,half,one
    use model, only : model_type
    use fault, only : fault_type
    use mpi_routines, only : MPI_REAL_PR
    use io, only : error
    use mpi

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! FIELD_NAME = names of fields to update
    ! DT = time step
    ! STEP = integration step
    ! S0 = location of beginning of time step (in Runge-Kutta method)
    ! S = current stage (in Runge-Kutta method)

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    character(*),intent(in) :: field_name
    real(pr),intent(in) :: dt
    character(*),intent(in) :: step
    integer(pin),intent(in) :: s0
    integer(pin),intent(in),optional :: s

    ! Internal Parameters:
    ! K = index in z direction
    ! N = index of new field values
    ! O = index of old field values
    ! R1 = lower index of rates
    ! R2 = upper index of rates
    ! NR = number of rates used in update
    ! R = coefficients of rate vector
    ! NX = number of points in x direction (for particular process)
    ! MY_ABS_ERR = absolute error (for particular process)
    ! MY_REL_ERR = relative error (for particular process)
    ! ABS_ERR = absolute error
    ! REL_ERR = relative error
    ! IERR = MPI error flag

    integer(pin) :: k,n,o,r1,r2,nr,nx,ierr
    real(pr),dimension(:),allocatable :: R
    real(pr) :: my_abs_err,my_rel_err,abs_err,rel_err

    ! return if process holds no x data 
    
    if (.not.mdl%holds_x) return

    select case(step)
       
    case default

       call error('invalid integration step','integrate_rate_eqns')

    case('explicit')
       
       o = s0   ! beginning of time step
       n = s0+1 ! end of time step
       r1 = s0
       r2 = s0
       allocate(R(1:1))
       R = one
  
    case('trapezoid')
       
       o = s0   ! beginning of time step
       n = s0+1 ! end of time step
       r1 = s0
       r2 = s0+1
       allocate(R(1:2))
       R = half

    case('rk_stage')
       
       if (.not.present(s)) &
            call error('current stage location must be given','integrate_rate_eqns')

       o = s0       ! beginning of time step
       n = (s0-1)+s ! current stage
       r1 = s0
       r2 = r1+s-2
       allocate(R(1:s-1))
       R = mdl%rk%a(s,1:s-1)
       
    case('rki_stage')
       
       if (.not.present(s)) &
            call error('current stage location must be given','integrate_rate_eqns')

       o = s0       ! beginning of time step
       n = (s0-1)+s ! current stage
       r1 = s0
       r2 = r1+s-1
       allocate(R(1:s))
       R = mdl%rk%ai(s,1:s)
       
    case('rk_final')

       o = s0              ! beginning of time step
       n = (s0-1)+mdl%rk%s ! end of time step
       r1 = s0
       r2 = (s0-1)+mdl%rk%s
       allocate(R(1:mdl%rk%s))
       R = mdl%rk%b

    case('rk_error')
       
       o = s0 ! beginning of time step
       r1 = s0
       r2 = (s0-1)+mdl%rk%s
       allocate(R(1:mdl%rk%s))
       R = mdl%rk%bh-mdl%rk%b
       
    end select
    
    nr = r2-r1+1
    nx = mdl%px-mdl%mx+1

    select case(step)

    case default

       !---update---
       
       if (field_name=='U'.or.field_name=='all') then
          
          if (mdl%bm) then
             call update_field(nx,mdl%ny,nr,flt%uxm(:,:,n),flt%uxm(:,:,o),flt%vxm(:,:,r1:r2),dt,R)
             call update_field(nx,mdl%ny,nr,flt%uxp(:,:,n),flt%uxp(:,:,o),flt%vxp(:,:,r1:r2),dt,R)
             call update_field(nx,mdl%ny,nr,flt%uym(:,:,n),flt%uym(:,:,o),flt%vym(:,:,r1:r2),dt,R)
             call update_field(nx,mdl%ny,nr,flt%uyp(:,:,n),flt%uyp(:,:,o),flt%vyp(:,:,r1:r2),dt,R)
             call update_field(nx,mdl%ny,nr,flt%uzm(:,:,n),flt%uzm(:,:,o),flt%vzm(:,:,r1:r2),dt,R)
             call update_field(nx,mdl%ny,nr,flt%uzp(:,:,n),flt%uzp(:,:,o),flt%vzp(:,:,r1:r2),dt,R)
          else
             call update_field(nx,mdl%ny,nr,flt%Ux(:,:,n),flt%Ux(:,:,o),flt%Vx(:,:,r1:r2),dt,R)
             call update_field(nx,mdl%ny,nr,flt%Uy(:,:,n),flt%Uy(:,:,o),flt%Vy(:,:,r1:r2),dt,R)
          end if
          
          call update_field(nx,mdl%ny,nr,flt%U(:,:,n),flt%U(:,:,o),flt%V(:,:,r1:r2),dt,R)
          
       end if
       
       if ((field_name=='V'.or.field_name=='all').and.mdl%friction_method=='strength_rate') &
            call update_field(nx,mdl%ny,nr,flt%V(:,:,n),flt%V(:,:,o),flt%dV(:,:,r1:r2),dt,R)
       
       if ((field_name=='S'.or.field_name=='all').and.mdl%friction_method=='strength_rate') &
            call update_field(nx,mdl%ny,nr,flt%S(:,:,n),flt%S(:,:,o),flt%dS(:,:,r1:r2),dt,R)
       
       if ((field_name=='Q'.or.field_name=='all').and.mdl%friction_method=='strength') &
            call update_field(nx,mdl%ny,nr,flt%Q(:,:,n),flt%Q(:,:,o),flt%dQ(:,:,r1:r2),dt,R)
       
       if ((field_name=='Tp'.or.field_name=='all').and.flt%thermpres) then

          do k = 1,mdl%nz
             call update_field(nx,mdl%ny,nr,flt%T(k,:,:,n),flt%T(k,:,:,o),flt%dT(k,:,:,r1:r2),dt,R)
             call update_field(nx,mdl%ny,nr,flt%p(k,:,:,n),flt%p(k,:,:,o),flt%dp(k,:,:,r1:r2),dt,R)
          end do
             
          flt%T0(:,:,n) = flt%T(1,:,:,n)
          flt%p0(:,:,n) = flt%p(1,:,:,n)
       
       end if

    case ('rk_error')

       !---error estimate---
    
       my_abs_err = zero
       my_rel_err = zero
       
       if (field_name=='U'.or.field_name=='all') &
            call error_field(nx,mdl%ny,nr,flt%U(:,:,o),flt%V(:,:,r1:r2),dt,R,my_abs_err,my_rel_err,flt%scale_U)

       if ((field_name=='S'.or.field_name=='all').and.mdl%friction_method=='strength_rate') &
            call error_field(nx,mdl%ny,nr,flt%S(:,:,o),flt%dS(:,:,r1:r2),dt,R,my_abs_err,my_rel_err,flt%scale_s)
       
       if ((field_name=='V'.or.field_name=='all').and.mdl%friction_method=='strength_rate') &
            call error_field(nx,mdl%ny,nr,flt%V(:,:,o),flt%dV(:,:,r1:r2),dt,R,my_abs_err,my_rel_err,flt%scale_Vmin,flt%scale_Vmax)

       if ((field_name=='Q'.or.field_name=='all').and.mdl%friction_method=='strength') &
            call error_field(nx,mdl%ny,nr,flt%Q(:,:,o),flt%dQ(:,:,r1:r2),dt,R,my_abs_err,my_rel_err,flt%scale_Q)

       if ((field_name=='Tp'.or.field_name=='all').and.flt%thermpres) then

          do k = 1,mdl%nz
             call error_field(nx,mdl%ny,nr,flt%T(k,:,:,o),flt%dT(k,:,:,r1:r2),dt,R,my_abs_err,my_rel_err,flt%scale_T)
             call error_field(nx,mdl%ny,nr,flt%p(k,:,:,o),flt%dp(k,:,:,r1:r2),dt,R,my_abs_err,my_rel_err,flt%scale_p)
          end do

       end if

       call MPI_AllReduce(my_abs_err,abs_err,1,MPI_REAL_PR,MPI_MAX,MPI_COMM_WORLD,ierr)
       call MPI_AllReduce(my_rel_err,rel_err,1,MPI_REAL_PR,MPI_MAX,MPI_COMM_WORLD,ierr)
       
       ! store error estimates
       
       flt%rerr(1) = rel_err
       flt%aerr(1) = abs_err

    end select
    
    if (allocated(R)) deallocate(R)
    
  end subroutine integrate_rate_eqns


  subroutine update_field(nx,ny,nr,f,f0,dfdt,dt,R)
    ! UPDATE_FIELD updates a field from f0 to f as f = f0+dt*R*dfdt
    ! 
    ! Modified: 7 August 2007

    implicit none

    ! I/O Parameters:
    ! NX = length in x direction
    ! NY = length in y direction
    ! NR = number of rates used in update
    ! F = final value of field
    ! F0 = initial value of field
    ! DFDT = rate of change of field
    ! DT = time step
    ! R = weighting coefficients for rates

    integer(pin),intent(in) :: nx,ny,nr
    real(pr),dimension(nx,ny),intent(out) :: f
    real(pr),dimension(nx,ny),intent(in) :: f0
    real(pr),dimension(nx,ny,nr),intent(in) :: dfdt
    real(pr),intent(in) :: dt
    real(pr),dimension(nr),intent(in) :: R

    ! Internal Parameters:
    ! IR = index of rates
    ! I = index of x
    ! J = index of y

    integer(pin) :: ir!,i,j

    ! update field

    ! option 1 (vectorized)

    f = f0+dt*dfdt(:,:,1)*R(1)
    do ir = 2,nr
       f = f+dt*dfdt(:,:,ir)*R(ir)
    end do

    ! option 2 (OpenMP)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
    !$OMP DO
    !do j = 1,ny
    !   do i = 1,nx
    !      f(i,j) = f0(i,j)+dt*dot_product(dfdt(i,j,:),R)
    !   end do
    !end do
    !$OMP END DO
    !$OMP END PARALLEL
    
  end subroutine update_field


  subroutine error_field(nx,ny,nr,f0,dfdt,dt,R,abs_err,rel_err,scale_f,scale_fmax)
    ! ERROR_FIELD estimates error from embedded Runge-Kutta formula
    ! 
    ! Modified: 27 August 2007

    use constants, only : zero

    implicit none

    ! I/O Parameters:
    ! NX = length in x direction
    ! NY = length in y direction
    ! NR = number of rates used in update
    ! F0 = initial value of field
    ! DFDT = rate of change of field
    ! DT = time step
    ! R = weighting coefficients for rates
    ! ABS_ERR = absolute error
    ! REL_ERR = relative error
    ! SCALE_F = scale of field (for relative error estimate)
    ! SCALE_FMAX = maximum scale of field (for relative error estimate)

    integer(pin),intent(in) :: nx,ny,nr
    real(pr),dimension(nx,ny),intent(in) :: f0
    real(pr),dimension(nx,ny,nr),intent(in) :: dfdt
    real(pr),intent(in) :: dt,scale_f
    real(pr),dimension(nr),intent(in) :: R
    real(pr),intent(inout) :: abs_err,rel_err
    real(pr),intent(in),optional :: scale_fmax

    ! Internal Parameters:
    ! IR = index of rates
    ! I = index of x
    ! J = index of y
    ! DF = update to field
    ! ERR = error in field (absolute value of update df)
    ! SCALE = scale of field (for relative error)
    
    integer(pin) :: ir !,i,j
    real(pr),dimension(nx,ny) :: df,err,scale

    ! calculate error with embedded Runge-Kutta formula

    ! option 1 (vectorized)

    df = dt*dfdt(:,:,1)*R(1)
    do ir = 2,nr
       df = df+dt*dfdt(:,:,ir)*R(ir)
    end do

    ! option 2 (OpenMP)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
    !$OMP DO
    !do j = 1,ny
    !   do i = 1,nx
    !      df(i,j) = dt*dot_product(dfdt(i,j,:),R)
    !   end do
    !end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! calculate absolute error
  
    err = abs(df)
    abs_err = max(abs_err,maxval(err))

    ! set scale for relative error estimates

    if (scale_f/=zero) then
       ! if scale is given (as some non-zero value), use that scale
       if (present(scale_fmax)) then
          ! use 'average' field value as scale, but only if between scale_f and scale_fmax
          scale = min(max(abs(f0)+dt*maxval(abs(dfdt),3),scale_f),scale_fmax)
       else
          scale = scale_f
       end if
    else
       ! otherwise determine an 'average' value of f over step
       scale = abs(f0)+dt*maxval(abs(dfdt),3)
    end if

    ! calculate relative error

    where(scale==zero) 
       err = zero ! no relative error estimate possible
    elsewhere
       err = err/scale ! relative error
    end where

    rel_err = max(rel_err,maxval(err))

  end subroutine error_field


  subroutine finalize_fields(flt,mdl)
    ! FINALIZE_FIELDS overwrites fields at time t with values from time t+dt
    ! 
    ! Modified: 30 May 2007

    use fault, only : fault_type
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! FLT = fault variables
    ! MDL = model variables

    type(fault_type),intent(inout) :: flt
    type(model_type),intent(in) :: mdl

    ! Internal Parameters:
    ! N1 = index of step containing fields at t
    ! N2 = index of step containing fields at t+dt

    integer(pin) :: n1,n2

    ! transfer field values from end of current time step (n2) to beginning of next (n1)

    n1 = 1

    select case(mdl%method)
    case('iterative','substep')
       n2 = 2
    case('adaptive')
       n2 = mdl%ns
    end select
    
    call move_fields(flt,mdl,n1,n2)
    call shift_stress_transfer(flt,mdl)

  end subroutine finalize_fields


  subroutine move_fields(flt,mdl,n1,n2)
    ! MOVE_FIELDS overwrites fields at index n1 with those at n2
    ! 
    ! Modified: 21 August 2007

    use fault, only : fault_type
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! FLT = fault variables
    ! MDL = model variables
    ! N1 = index of field to be overwritten
    ! N2 = index of field to replace old value

    type(fault_type),intent(inout) :: flt
    type(model_type),intent(in) :: mdl
    integer(pin),intent(in) :: n1,n2

    if (mdl%bm) then
       
       flt%uxp(:,:,n1) = flt%uxp(:,:,n2)
       flt%uxm(:,:,n1) = flt%uxm(:,:,n2)
       flt%uyp(:,:,n1) = flt%uyp(:,:,n2)
       flt%uym(:,:,n1) = flt%uym(:,:,n2)
       flt%uzp(:,:,n1) = flt%uzp(:,:,n2)
       flt%uzm(:,:,n1) = flt%uzm(:,:,n2)
       
       flt%vxp(:,:,n1) = flt%vxp(:,:,n2)
       flt%vxm(:,:,n1) = flt%vxm(:,:,n2)
       flt%vyp(:,:,n1) = flt%vyp(:,:,n2)
       flt%vym(:,:,n1) = flt%vym(:,:,n2)
       flt%vzp(:,:,n1) = flt%vzp(:,:,n2)
       flt%vzm(:,:,n1) = flt%vzm(:,:,n2)
       
    else
       
       flt%Ux(:,:,n1) = flt%Ux(:,:,n2)
       flt%Uy(:,:,n1) = flt%Uy(:,:,n2)
       
       flt%Vx(:,:,n1) = flt%Vx(:,:,n2)
       flt%Vy(:,:,n1) = flt%Vy(:,:,n2)
       
    end if
    
    flt%U(:,:,n1) = flt%U(:,:,n2)
    
    flt%V( :,:,n1) = flt%V( :,:,n2)
    flt%dV(:,:,n1) = flt%dV(:,:,n2)
    
    flt%sx(:,:,n1) = flt%sx(:,:,n2)
    flt%sy(:,:,n1) = flt%sy(:,:,n2)
    flt%sz(:,:,n1) = flt%sz(:,:,n2)
    
    flt%Q( :,:,n1) = flt%Q( :,:,n2)
    flt%dQ(:,:,n1) = flt%dQ(:,:,n2)
    
    flt%S( :,:,n1) = flt%S( :,:,n2)
    flt%dS(:,:,n1) = flt%dS(:,:,n2)
    
    if (allocated(flt%T0)) flt%T0(:,:,n1) = flt%T0(:,:,n2)
    if (allocated(flt%p0)) flt%p0(:,:,n1) = flt%p0(:,:,n2)
    
    if (allocated(flt%T )) flt%T( :,:,:,n1) = flt%T( :,:,:,n2)
    if (allocated(flt%dT)) flt%dT(:,:,:,n1) = flt%dT(:,:,:,n2)

    if (allocated(flt%p )) flt%p( :,:,:,n1) = flt%p(:,:,:,n2)
    if (allocated(flt%dp)) flt%dp(:,:,:,n1) = flt%dp(:,:,:,n2)

  end subroutine move_fields


  subroutine move_rates(flt,mdl,n1,n2)
    ! MOVE_RATES overwrites rates at index n1 with those at n2
    ! 
    ! Modified: 21 August 2007

    use fault, only : fault_type
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! FLT = fault variables
    ! MDL = model variables
    ! N1 = index of field to be overwritten
    ! N2 = index of field to replace old value

    type(fault_type),intent(inout) :: flt
    type(model_type),intent(in) :: mdl
    integer(pin),intent(in) :: n1,n2

    if (mdl%bm) then
       
       flt%vxp(:,:,n1) = flt%vxp(:,:,n2)
       flt%vxm(:,:,n1) = flt%vxm(:,:,n2)
       flt%vyp(:,:,n1) = flt%vyp(:,:,n2)
       flt%vym(:,:,n1) = flt%vym(:,:,n2)
       flt%vzp(:,:,n1) = flt%vzp(:,:,n2)
       flt%vzm(:,:,n1) = flt%vzm(:,:,n2)
       
    else
       
       flt%Vx(:,:,n1) = flt%Vx(:,:,n2)
       flt%Vy(:,:,n1) = flt%Vy(:,:,n2)
       
    end if
    
    flt%V( :,:,n1) = flt%V( :,:,n2)
    flt%dV(:,:,n1) = flt%dV(:,:,n2)
    
    flt%sx(:,:,n1) = flt%sx(:,:,n2)
    flt%sy(:,:,n1) = flt%sy(:,:,n2)
    flt%sz(:,:,n1) = flt%sz(:,:,n2)
    
    flt%dQ(:,:,n1) = flt%dQ(:,:,n2)
    flt%dS(:,:,n1) = flt%dS(:,:,n2)

    if (allocated(flt%dT)) flt%dT(:,:,:,n1) = flt%dT(:,:,:,n2)
    
    if (allocated(flt%dp)) flt%dp(:,:,:,n1) = flt%dp(:,:,:,n2)

  end subroutine move_rates


  subroutine shift_stress_transfer(flt,mdl)
    ! SHIFT_STRESS_TRANSFER shifts stress transfer to make room for new time step
    ! 
    ! Modified: 14 February 2007

    use fault, only : fault_type
    use model, only : model_type

    implicit none

    ! I/O Parameters:
    ! FLT = fault variables
    ! MDL = model variables

    type(fault_type),intent(inout) :: flt
    type(model_type),intent(in) :: mdl

    ! return if no shift required

    if (mdl%method/='substep') return

    ! shift stress transfer

    if (mdl%bm) then

       flt%fxp = eoshift(flt%fxp,1,dim=3)
       flt%fxm = eoshift(flt%fxm,1,dim=3)
       flt%fyp = eoshift(flt%fyp,1,dim=3)
       flt%fym = eoshift(flt%fym,1,dim=3)
       flt%fzp = eoshift(flt%fzp,1,dim=3)
       flt%fzm = eoshift(flt%fzm,1,dim=3)
       
    else
       
       flt%fx  = eoshift(flt%fx, 1,dim=3)
       flt%fy  = eoshift(flt%fy, 1,dim=3)

    end if

  end subroutine shift_stress_transfer


  subroutine move_error(flt)
    ! MOVE_ERROR shifts error estimate array
    ! 
    ! Modified: 21 September 2007

    use constants, only : zero
    use fault, only : fault_type

    implicit none

    ! I/O Parameters:
    ! FLT = fault variables

    type(fault_type),intent(inout) :: flt

    flt%rerr = eoshift(flt%rerr,1,zero)
    flt%aerr = eoshift(flt%aerr,1,zero)

  end subroutine move_error


  subroutine predict_rk(mdl,flt,dt,s0,s)
    ! PREDICT_RK predicting stage values of V and S (to speed convergence
    ! of Newton's method for implicit Runge-Kutta solver)
    ! 
    ! Modified: 27 August 2007

    use constants, only : half
    use model, only : model_type
    use fault, only : fault_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! DT = time step
    ! S0 = first integration stage
    ! S = stage at which values are to be predicted

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    real(pr),intent(in) :: dt
    integer(pin),intent(in) :: s0,s

    ! Internal Parameters:
    ! SRK = Runge-Kutta stage
    ! DTS = time step over which interpolation occurs

    integer(pin) :: srk
    real(pr) :: dts

    ! set Runge-Kutta stage index

    srk = s-s0+1

    ! set time step between stages s and s-1

    dts = (mdl%rk%c(srk)-mdl%rk%c(srk-1))*dt

    ! predict stage values during next step
    
    if (mdl%bm) then
       flt%vxp(:,:,s) = flt%vxp(:,:,s-1)+half*flt%dV(:,:,s-1)*dts
       flt%vxm(:,:,s) = flt%vxm(:,:,s-1)-half*flt%dV(:,:,s-1)*dts
       flt%vyp(:,:,s) = flt%vyp(:,:,s-1)
       flt%vym(:,:,s) = flt%vym(:,:,s-1)
    else
       flt%Vx(:,:,s) = flt%Vx(:,:,s-1)+flt%dV(:,:,s-1)*dts
       flt%Vy(:,:,s) = flt%Vy(:,:,s-1)
    end if
    flt%V(:,:,s) = flt%V(:,:,s-1)+flt%dV(:,:,s-1)*dts
    
    flt%sx(:,:,s) = flt%Sx(:,:,s-1)+flt%dS(:,:,s-1)*dts
    flt%sy(:,:,s) = flt%Sy(:,:,s-1)
    flt%S( :,:,s) = flt%S( :,:,s-1)+flt%dS(:,:,s-1)*dts
       
  end subroutine predict_rk


end module integration
