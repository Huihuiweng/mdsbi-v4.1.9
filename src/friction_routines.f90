! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module friction_routines

  ! FRICTION_ROUTINES contains subroutines to call friction law
  ! solvers for identical material or bimaterial problems
  ! 
  ! Modified: 6 August 2010

  use constants, only : pr,pin

  implicit none

contains


  subroutine solve_friction(mdl,cnv,flt,fri,s,accept,dt,s0,direct_only)
    ! SOLVE_FRICTION solves the friction law coupled with elasticity
    ! for unknown velocity and stress
    ! 
    ! Modified: 5 December 2007

    use model, only : model_type
    use fault, only : fault_type
    use friction, only : friction_type
    use convolution, only : convolution_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! FLT = fault variables
    ! FRI = friction variables
    ! S = integration stage for storage
    ! ACCEPT = accept solution
    ! DT = time step
    ! S0 = initial integration stage (corresponding to start of step)
    ! DIRECT_ONLY = solve for response by direct effect only (input value)

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(in) :: cnv
    type(fault_type),intent(inout) :: flt
    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: s
    logical,intent(out) :: accept
    real(pr),intent(in),optional :: dt
    integer(pin),intent(in),optional :: s0
    logical,intent(in),optional :: direct_only

    ! Separate normal and shear solvers:
       
    ! set fault-normal velocity vz(t) and effective normal stress sz(t)
    
    call solve_friction_normal(mdl,flt,cnv,s)
    
    ! solve elasticity equation coupled with friction for V(t) and s(t) at each point using
    ! U(t), f(t), Q(t), p(t), T(t)
    
    if (present(dt).and.present(s0)) then
       if (present(direct_only)) then
          call solve_friction_shear(mdl,flt,fri,s,accept,dt,s0,direct_only)
       else
          call solve_friction_shear(mdl,flt,fri,s,accept,dt,s0)
       end if
    else
       call solve_friction_shear(mdl,flt,fri,s,accept)
    end if

    ! Combined normal and shear solvers:
    ! solve elasticity equation coupled with friction for V(t) and s(t) at each point using
    ! U(t), f(t), Q(t), p(t), T(t), as well as fault-normal velocity vz(t) and effective normal stress sz(t)
       
    ! solve elasticity equation coupled with friction for V(t) and s(t) at each point using
    ! U(t), f(t), Q(t), p(t), T(t), as well as fault-normal velocity vz(t) and effective normal stress sz(t)
    
    !call solve_friction_combined(mdl,cnv,flt,fri,s,accept)
    
  end subroutine solve_friction


  subroutine solve_friction_normal(mdl,flt,cnv,s)
    ! SOLVE_FRICTION_NORMAL calculates effective normal stress by adding contributions from
    ! 1. normal stress loads (in absence of slip)
    ! 2. normal stress changes due to elastic mismatch, assuming fault walls remain in contact
    ! 3. pore pressure changes due to poroelastic response
    ! 4. pore pressure changes due to thermal pressurization.  
    ! and also fault-normal velocities
    ! 
    ! Modified: 6 December 2007

    use constants, only : zero,half,one,two,three
    use model, only : model_type
    use fault, only : fault_type,stress_invariants
    use convolution, only : convolution_type
    use fourier, only : strain
    use io, only : error
    
    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! CNV = convolution variables
    ! S = integration substep

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(convolution_type),intent(in) :: cnv
    integer(pin),intent(in) :: s

    ! Internal Parameters:
    ! NP = mean normal stress, plus side
    ! NM = mean normal stress, minus side
    ! TP = deviatoric stress, plus side
    ! TM = deviatoric stress, minus side

    real(pr),dimension(:,:),allocatable :: qp,qm,q0,QQp,QQm,zp,zm,Np,Nm,Tp,Tm

    ! return if process holds no x data 

    if (.not.mdl%holds_x) return

    ! fault-normal velocity, assuming fault walls remain in contact

    if (mdl%bm) then
       flt%vzp(:,:,s) = (flt%fzpi-flt%fzmi)/(mdl%cvnp+mdl%cvnm)
       flt%vzm(:,:,s) = flt%vzp(:,:,s)
    end if
       
    ! (total) normal stress changes from slip, assuming fault walls remain in contact

    if (mdl%bm) then
       flt%sz(:,:,s) = flt%fzpi-mdl%cvnp*flt%vzp(:,:,s)
    else
       flt%sz(:,:,s) = zero
    end if

    ! poroelastic fault zone response

    if (flt%poroelastic) then

       ! allocate arrays

       allocate( &
            qp(mdl%mx:mdl%px,mdl%ny),qm(mdl%mx:mdl%px,mdl%ny),q0(mdl%mx:mdl%px,mdl%ny), &
            QQp(mdl%mx:mdl%px,mdl%ny),QQm(mdl%mx:mdl%px,mdl%ny), &
            zp(mdl%mx:mdl%px,mdl%ny),zm(mdl%mx:mdl%px,mdl%ny) )

       ! calculate fault-parallel strains 
       
       call strain(mdl,cnv,flt%exxp,flt%exxm,flt%exyp,flt%exym,flt%eyyp,flt%eyym)
       
       ! check for plastic yielding

       if (flt%plastic) then
          
          ! allocate arrays
          
          allocate( &
               Np(mdl%mx:mdl%px,mdl%ny),Nm(mdl%mx:mdl%px,mdl%ny), &
               Tp(mdl%mx:mdl%px,mdl%ny),Tm(mdl%mx:mdl%px,mdl%ny) )
          
          ! calculate stress invariants
          
          call stress_invariants(flt%mufp,flt%nufp,flt%Bp, &
               flt%sx(:,:,s),flt%sy(:,:,s),flt%sz(:,:,s)+flt%sz0,flt%sz0, &
               flt%exxp,flt%exyp,flt%eyyp,flt%sxx0,flt%sxy0,flt%syy0,Np,Tp)
          
          call stress_invariants(flt%mufm,flt%nufm,flt%Bm, &
               flt%sx(:,:,s),flt%sy(:,:,s),flt%sz(:,:,s)+flt%sz0,flt%sz0, &
               flt%exxm,flt%exym,flt%eyym,flt%sxx0,flt%sxy0,flt%syy0,Nm,Tm)

          ! check for yielding and update maximum 

          flt%Yp = max(flt%Yp,Tp-flt%muip*Np+flt%cohp)
          flt%Ym = max(flt%Ym,Tm-flt%muim*Nm+flt%cohm)

          ! degrade material instantaneously upon reaching yield

          where (flt%Yp>zero)
             flt%mufp = half*flt%mufp0
             flt%nufp = 0.3_pr
             flt%Bp = 0.8_pr
             flt%Zp = sqrt(20._pr)*flt%Zp0
          end where

          where (flt%Ym>zero)
             flt%mufm = half*flt%mufm0
             flt%nufm = 0.3_pr
             flt%Bm = 0.8_pr
             flt%Zm = sqrt(20._pr)*flt%Zm0
          end where

          !if (any(flt%Yp>zero.or.flt%Ym>zero)) then
          !   print *
          !   print *, 'n=',mdl%n
          !end if

          !j = 1
          !do i = mdl%mx,mdl%px
          !   if (flt%Yp(i,j)>zero.or.flt%Ym(i,j)>zero) then
          !      print *, i,flt%Yp(i,j),flt%mufp(i,j),flt%nufp(i,j),flt%Bp(i,j),flt%Zp(i,j)
          !      print *, i,flt%Ym(i,j),flt%mufm(i,j),flt%nufm(i,j),flt%Bm(i,j),flt%Zm(i,j)
          !   end if
          !end do

          ! deallocate arrays

          deallocate(Np,Nm,Tp,Tm)
       
       end if

       ! calculate poroelastic properties

       qp = flt%qp
       qm = flt%qm
       q0 = flt%q0

       !QQp = two*(one+flt%nufp)/(three*(one-flt%nufp))*flt%Bp
       !QQm = two*(one+flt%nufm)/(three*(one-flt%nufm))*flt%Bm
       !zp = flt%Zp/(flt%Zp+flt%Zm)
       !zm = flt%Zm/(flt%Zp+flt%Zm)
       !qp = (flt%mufp/mdl%mup)*zp*QQp
       !qm = (flt%mufm/mdl%mum)*zm*QQm
       !q0 = half*(zp*QQp+zm*QQm)

       ! calculate pore pressure change from fault-parallel strains and normal stress

       flt%p0(:,:,s) = -q0*flt%sz(:,:,s) &
            -qp*mdl%mup*(flt%exxp+flt%eyyp) &
            -qm*mdl%mum*(flt%exxm+flt%eyym)

       ! add this to effective normal stress
    
       flt%sz(:,:,s) = flt%sz(:,:,s)+flt%p0(:,:,s)

       ! deallocate arrays

       deallocate(qp,qm,q0,QQp,QQm,zp,zm)

    end if

    ! add pore pressure changes from thermal pressurization

    if (flt%thermpres) &
         flt%sz(:,:,s) = flt%sz(:,:,s)+flt%p0(:,:,s)

    ! add normal load

    flt%sz(:,:,s) = flt%sz(:,:,s)+flt%sz0

    ! switch boundary conditions to free surface to prevent
    ! effective normal stress from becoming tensile

    if (mdl%bm) then
       if (any(flt%sz(:,:,s)>zero).and.(.not.mdl%opening)) &
            call error('Tensile normal stress but opening not permitted','solve_friction_normal')
       where (flt%sz(:,:,s)>zero)
          flt%sz(:,:,s) = zero
          flt%vzp(:,:,s) =  (flt%sz0+flt%fzpi)/mdl%cvnp
          flt%vzm(:,:,s) = -(flt%sz0+flt%fzmi)/mdl%cvnm
       end where
    end if
  
  end subroutine solve_friction_normal
 

  subroutine solve_friction_shear(mdl,flt,fri,s,accept,dt,s0,direct_only)
    ! SOLVE_FRICTION_SHEAR solves the friction law coupled with elasticity
    ! for unknown velocity and stress
    ! 
    ! Modified: 6 December 2007

    use model, only : model_type
    use fault, only : fault_type
    use friction, only : friction_type
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! FRI = friction variables
    ! S = integration stage for storage
    ! ACCEPT = accept solution
    ! DT = time step
    ! S0 = initial integration stage (corresponding to start of step)
    ! DIRECT_ONLY = solve for response by direct effect only (input value)

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: s
    logical,intent(out) :: accept
    real(pr),intent(in),optional :: dt
    integer(pin),intent(in),optional :: s0
    logical,intent(in),optional :: direct_only

    ! Internal Parameters:
    ! I = index in x direction
    ! J = index in y direction
    ! DIRECT = solve for response by direct effect only

    integer(pin) :: i,j
    logical :: direct

    ! return if process holds no x data 

    if (.not.mdl%holds_x) then
       accept = .true.
       return
    end if

    ! loop over fault

    do j = 1,mdl%ny
       do i = mdl%mx,mdl%px

          if (allocated(flt%T0)) fri%T = flt%T0(i,j,s)
          
          select case(mdl%friction_method)
             
          case default

             call error('Invalid friction_method:' // trim(adjustl(mdl%friction_method)) ,'')

          case('strength')

             call friction_strength(mdl,flt,fri,i,j,s,accept)

          case('strength_rate')

             if (.not.(present(dt))) &
                  call error('Time step, dt, needed for strength-rate form of friction law','solve_friction_shear')
             if (.not.(present(s0))) &
                  call error('Initial integration stage, s0, needed for strength-rate form of friction law','solve_friction_shear')

             if (present(direct_only)) then
                direct = direct_only
             else
                direct = .false.
             end if

             call friction_strength_rate(mdl,flt,fri,i,j,s,accept,dt,s0,direct)

          end select

          ! if unacceptable, reject this step and exit routine
          ! (usually to try again with smaller step size)
          
          if (.not.accept) return
          
       end do
    end do

  end subroutine solve_friction_shear


  subroutine friction_strength(mdl,flt,fri,i,j,s,accept)
    ! FRICTION_STRENGTH solves the friction law coupled with elasticity
    ! for unknown velocity and stress, strength formulation
    ! 
    ! Modified: 6 December 2007

    use constants, only : zero,half,one
    use model, only : model_type
    use fault, only : fault_type
    use friction, only : friction_type
    use friction_im, only : couple_elastic_strength_im
    use friction_bm, only : couple_elastic_strength_bm
    use ratestate, only : check_state_rs
    use slipweak, only : strength_sw,evolution_sw
    use io, only : warning

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! FRI = friction variables
    ! I = index in x direction
    ! J = index in y direction
    ! S = integration stage for storage
    ! ACCEPT = accept solution

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: i,j,s
    logical,intent(out) :: accept

    ! Internal Parameters:
    ! SO = velocity-independent portion of strength
    ! ERR = flag that indicates error in solving friction law
    ! INFO = print information from friction solver

    real(pr) :: So
    logical :: err,info

    ! no information by default

    info = .false.

    ! check if value of state is acceptable, otherwise reject step
    
    select case(fri%friction)
    case default
       accept = .true.
    case('ratestate')
       accept = check_state_rs(flt%Q(i,j,s),fri%rs)
    end select
    
    ! if unacceptable, reject this step and exit routine
    ! (usually to try again with smaller step size)
    
    if (.not.accept) return
    
    ! calculate strength, but only if velocity-independent
    
    select case(fri%friction)
    case('slipweak')
       So = strength_sw(flt%U(i,j,s),flt%Q(i,j,s),flt%sz(i,j,s), &
            i,j,(mdl%x(i)**2+mdl%y(j)**2)**0.5,mdl%t,fri%sw)
    case('constfV')
       ! constant f,V: set V=1 and calculate new shear resistance
       ! (assuming f=1 and no resistance in y direction)
       if (mdl%bm) then
          flt%vxp(i,j,s) =  half
          flt%vxm(i,j,s) = -half
          flt%vyp(i,j,s) = zero
          flt%vym(i,j,s) = zero
          flt%V(i,j,s) = sqrt( &
               (flt%vxp(i,j,s)-flt%vxm(i,j,s))**2+ &
               (flt%vyp(i,j,s)-flt%vym(i,j,s))**2)
       else
          flt%Vx(i,j,s) = one
          flt%Vy(i,j,s) = zero
          flt%V(i,j,s) = sqrt(flt%Vx(i,j,s)**2+flt%Vy(i,j,s)**2)
       end if
       flt%sx(i,j,s) = -flt%sz(i,j,s)
       flt%sy(i,j,s) = zero
       return
    case default
       So = zero
    end select
             
    ! solve elasticity coupled with possibly velocity-dependent friction law;
    ! perturb initial guess at Vx if originally zero and using rate-and-state law
             
    if (mdl%bm) then
       
       if (fri%friction=='ratestate'.and.flt%vxp(i,j,s)-flt%vxm(i,j,s)==zero) then
          if (flt%scale_Vmin<=zero) &
               call warning('Friction solver may have problems, set scale_Vmin>0','friction_strength')
          flt%vxp(i,j,s) =  half*flt%scale_Vmin
          flt%vxm(i,j,s) = -half*flt%scale_Vmin
       end if
       
       call couple_elastic_strength_bm(fri,i,j, &
            So,flt%Q(i,j,s),flt%sz(i,j,s), &
            flt%vxp(i,j,s),flt%vxm(i,j,s), &
            flt%vyp(i,j,s),flt%vym(i,j,s), &
            flt%scale_Vmin,flt%scale_Vmax, &
            flt%sx(i,j,s),flt%sy(i,j,s), &
            flt%scale_s, &
            flt%sx0(i,j),flt%sy0(i,j), &
            flt%fxpi(i,j),flt%fxmi(i,j), &
            flt%fypi(i,j),flt%fymi(i,j), &
            mdl%cvsp,mdl%cvsm,mdl%transverse_slip,err)
       
       flt%V(i,j,s) = sqrt( &
            (flt%vxp(i,j,s)-flt%vxm(i,j,s))**2+ &
            (flt%vyp(i,j,s)-flt%vym(i,j,s))**2)
       
    else
                
       if (fri%friction=='ratestate'.and.flt%Vx(i,j,s)==zero) then
          if (flt%scale_Vmin<=zero) &
               call warning('Friction solver may have problems, set scale_Vmin>0','friction_strength')
          flt%Vx(i,j,s) = flt%scale_Vmin
       end if

!!$       if (abs(mdl%x(i))>=24d0) then
!!$
!!$          ! force sliding below some depth
!!$          flt%Vx(i,j,s) = 1.109106761539448d-9
!!$          flt%sx(i,j,s) = flt%sx0(i,j)+flt%fxi(i,j)-mdl%cV*flt%Vx(i,j,s)
!!$          err = .false.
!!$          flt%Vy(i,j,s) = 0d0
!!$          flt%sy(i,j,s) = flt%sy0(i,j)+flt%fyi(i,j)-mdl%cV*flt%Vy(i,j,s)
!!$
!!$       else

          call couple_elastic_strength_im(fri,i,j, &
               So,flt%Q(i,j,s),flt%sz(i,j,s), &
               flt%Vx(i,j,s),flt%Vy(i,j,s), &
               flt%scale_Vmin,flt%scale_Vmax, &
               flt%sx(i,j,s),flt%sy(i,j,s), &
               flt%scale_s, &
               flt%sx0(i,j),flt%sy0(i,j), &
               flt%fxi(i,j),flt%fyi(i,j), &
               mdl%cV,mdl%transverse_slip,err,info)
          
!!$       end if

       flt%V(i,j,s) = sqrt(flt%Vx(i,j,s)**2+flt%Vy(i,j,s)**2)
       
    end if
    
    ! if no solution can be found, reject this step and exit routine
    ! (usually to try again with smaller step size)
    
    if (err) then
       accept = .false.
       return
    end if
             
    ! set rate of change of strength for regularized slip-weakening law

    if (fri%friction=='slipweak') &
         flt%dQ(i,j,s) = evolution_sw(flt%U(i,j,s),flt%V(i,j,s),flt%Q(i,j,s), &
         flt%sz(i,j,s),i,j,(mdl%x(i)**2+mdl%y(j)**2)**0.5,mdl%t,fri%sw)

  end subroutine friction_strength


  subroutine friction_strength_rate(mdl,flt,fri,i,j,s,accept,dt,s0,direct_only)
    ! FRICTION_STRENGTH_RATE solves the friction law coupled with elasticity
    ! for unknown velocity and stress, strength rate formulation
    ! 
    ! Modified: 1 March 2007

    use constants, only : zero,half,one
    use model, only : model_type
    use fault, only : fault_type
    use friction, only : friction_type
    use friction_im, only : couple_elastic_strength_rate_im
    use friction_bm, only : couple_elastic_strength_rate_bm
    use io, only : warning

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! FRI = friction variables
    ! I = index in x direction
    ! J = index in y direction
    ! S = integration stage for storage
    ! ACCEPT = accept solution
    ! DT = time step
    ! S0 = initial integration stage (corresponding to start of step)
    ! DIRECT_ONLY = solve for response by direct effect only (input value)

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: i,j,s
    logical,intent(out) :: accept
    real(pr),intent(in) :: dt
    integer(pin),intent(in) :: s0
    logical,intent(in) :: direct_only

    ! Internal Parameters:
    ! A = coefficient of current rate
    ! DSO = contribution to rate of change of strength from previous integration stages
    ! DVO = contribution to rate of change of strength from previous integration stages
    ! SRK = current Runge-Kutta stage
    ! ERR = flag that indicates error in solving friction law

    integer(pin) :: srk
    real(pr) :: a,dVo,dSo
    logical :: err

    accept = .true.
    
    ! RK stages are 1:ns and corresponding storage stages are s0:ns+s0-1;
    ! current storage stage is s and current RK stage is srk=s-s0+1
    
    srk = s-s0+1
    
    ! implicit RK formula is
    ! Y(i) = y+dt*[a(i,1)*dY(1)+...+a(i,i-1)*dY(i-1)]+dt*a(i,i)*dY(i)
    !      = y+dt*dYo                                +dt*a     *dY(i)
    ! => dY(i) = [(Y(i)-y)/dt-dYo]/a (used in Newton's method solver)
    
    if (s==s0) then
       dVo = zero
       dSo = zero
    else
       dVo = dot_product(mdl%rk%ai(srk,1:srk-1),flt%dV(i,j,s0:s-1))
       dSo = dot_product(mdl%rk%ai(srk,1:srk-1),flt%dS(i,j,s0:s-1))
    end if
    a = mdl%rk%ai(srk,srk)
    
    if (mdl%bm) then
       
       ! perturb initial guess at Vx if originally zero and using rate-and-state law
       
       if (fri%friction=='ratestate'.and.flt%vxp(i,j,s)-flt%vxm(i,j,s)==zero) then
          if (flt%scale_Vmin<=zero) &
               call warning('Friction solver may have problems, set scale_Vmin>0','friction_strength_rate')
          flt%vxp(i,j,s) =  half*flt%scale_Vmin
          flt%vxm(i,j,s) = -half*flt%scale_Vmin
       end if
       
       ! solve elasticity coupled with velocity-dependent friction law
       
       call couple_elastic_strength_rate_bm(fri,i,j, &
            flt%sz(i,j,s),dt,a, &
            flt%V(i,j,s0),flt%S(i,j,s0),dVo,dSo, &
            flt%vxp(i,j,s),flt%vxm(i,j,s), &
            flt%vyp(i,j,s),flt%vym(i,j,s), &
            flt%scale_Vmin,flt%scale_Vmax, &
            flt%sx(i,j,s),flt%sy(i,j,s), &
            flt%scale_s, &
            flt%sx0(i,j),flt%sy0(i,j), &
            flt%fxpi(i,j),flt%fxmi(i,j), &
            flt%fypi(i,j),flt%fymi(i,j), &
            mdl%cvsp,mdl%cvsm,mdl%transverse_slip,direct_only,err)
       
       flt%V(i,j,s) = sqrt( &
            (flt%vxp(i,j,s)-flt%vxm(i,j,s))**2+ &
            (flt%vyp(i,j,s)-flt%vym(i,j,s))**2)
       
    else
       
       ! perturb initial guess at Vx if originally zero and using rate-and-state law
       
       if (fri%friction=='ratestate'.and.flt%Vx(i,j,s)==zero) then
          if (flt%scale_Vmin<=zero) &
               call warning('Friction solver may have problems, set scale_Vmin>0','friction_strength_rate')
          flt%Vx(i,j,s) = flt%scale_Vmin
       end if
       
       ! solve elasticity coupled with velocity-dependent friction law
       
       call couple_elastic_strength_rate_im(fri,i,j, &
            flt%sz(i,j,s),dt,a,&
            flt%V(i,j,s0),flt%S(i,j,s0),dVo,dSo, &
            flt%Vx(i,j,s),flt%Vy(i,j,s), &
            flt%scale_Vmin,flt%scale_Vmax, &
            flt%sx(i,j,s),flt%sy(i,j,s), &
            flt%scale_s, &
            flt%sx0(i,j),flt%sy0(i,j), &
            flt%fxi(i,j),flt%fyi(i,j), &
            mdl%cV,mdl%transverse_slip,direct_only,err)
       
       flt%V(i,j,s) = sqrt(flt%Vx(i,j,s)**2+flt%Vy(i,j,s)**2)
       
    end if
    
    ! if no solution can be found, reject this step and exit routine
    ! (usually to try again with smaller step size)
             
    if (err) then
       accept = .false.
       return
    end if
             
    flt%S(i,j,s) = sqrt(flt%sx(i,j,s)**2+flt%sy(i,j,s)**2)

    if (direct_only) return

    flt%dV(i,j,s) = ((flt%V(i,j,s)-flt%V(i,j,s0))/dt-dVo)/a
    flt%dS(i,j,s) = ((flt%S(i,j,s)-flt%S(i,j,s0))/dt-dSo)/a

  end subroutine friction_strength_rate


  subroutine state_rate(mdl,flt,fri,s)
    ! STATE_RATE solves for the state-rate
    ! 
    ! Modified: 19 September 2007

    use model, only : model_type
    use fault, only : fault_type
    use friction, only : friction_type
    use ratestate, only : evolution_rs

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! FRI = friction variables
    ! S = integration sub-step

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: s

    ! Internal Parameters:
    ! I = index in x direction
    ! J = index in y direction

    integer(pin) :: i,j

    ! return if process holds no x data 

    if (.not.mdl%holds_x) return

    ! return if not using rate-and-state friction

    if (fri%friction/='ratestate') return
    
    ! return if not using state variable formulation

    if (mdl%friction_method/='strength') return

    ! loop over fault

    do j = 1,mdl%ny
       do i = mdl%mx,mdl%px

          if (allocated(flt%T0)) fri%T = flt%T0(i,j,s)

          flt%dQ(i,j,s) = evolution_rs(flt%V(i,j,s),flt%Q(i,j,s),fri%T,i,j,fri%rs)
          
       end do
    end do

  end subroutine state_rate


  subroutine strength_rate(mdl,flt,fri,stage)
    ! STRENGTH_RATE solves for the strength-rate
    ! 
    ! Modified: 19 September 2007

    use model, only : model_type
    use fault, only : fault_type
    use friction, only : friction_type
    use ratestate, only : directeffect_rs,rate_rs

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! FRI = friction variables
    ! STAGE = integration stage

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: stage

    ! Internal Parameters:
    ! I = index in x direction
    ! J = index in y direction
    ! D = direct effect (coefficient of dV/dt)
    ! R = rate (how strength evolves, excluding direct effect)
    ! ERR = error flag (true if error)

    integer(pin) :: i,j
    real(pr) :: D,R
    logical :: err

    ! return if process holds no x data 

    if (.not.mdl%holds_x) return

    ! return if not using rate-and-state friction

    if (fri%friction/='ratestate') return

    ! loop over fault

    do j = 1,mdl%ny
       do i = mdl%mx,mdl%px

          if (allocated(flt%T0)) fri%T = flt%T0(i,j,stage)

          call directeffect_rs(fri%rs,flt%V(i,j,stage),flt%S(i,j,stage),fri%T,flt%sz(i,j,stage),i,j,D)
          call rate_rs(fri%rs,flt%V(i,j,stage),flt%S(i,j,stage),fri%T,flt%sz(i,j,stage),i,j,R,err=err)

          if (err) return

          flt%dS(i,j,stage) = D*flt%dV(i,j,stage)+R
          
       end do
    end do

  end subroutine strength_rate


  subroutine static_friction(flt,mdl,fri,s,g,m,p)
    ! STATIC_FRICTION evaluates the misfit between friction law and elasticity
    ! 
    ! Modified: 2 January 2007

    use constants, only : zero
    use model, only : model_type
    use fault, only : fault_type
    use friction, only : friction_type
    use slipweak, only : strength_sw
    use friction_static, only : static_elasticity_im,static_elasticity_bm

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! FRI = friction variables
    ! S = integration sub-step
    ! G = misfit
    ! M = minimum index of slipping zone in x direction
    ! P = maximum index of slipping zone in x direction

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: s
    real(pr),dimension(:),intent(out) :: g
    integer(pin),intent(in) :: m,p

    ! Internal Parameters:
    ! I = index in x direction
    ! J = index in y direction
    ! ST = strength
    ! MISFIT = misfit function (in array format)
    ! N = size of misfit

    integer(pin) :: i,j,n
    real(pr) :: st
    real(pr),dimension(:,:,:),allocatable :: misfit

    ! allocate memory to misfit array

    if (mdl%bm) then
       allocate(misfit(m:p,mdl%ny,4))
    else
       allocate(misfit(m:p,mdl%ny,2))
    end if

    ! set stresses

    if (mdl%bm) then
       flt%sx(:,:,s) = flt%sx0+flt%fxpi
       flt%sy(:,:,s) = flt%sy0+flt%fypi
    else
       flt%sx(:,:,s) = flt%sx0+flt%fxi
       flt%sy(:,:,s) = flt%sy0+flt%fyi
    end if

    ! set vector slip

    if (mdl%bm) then
       flt%U(:,:,s) = sqrt( &
            (flt%uxp(:,:,s)-flt%uxm(:,:,s))**2+ &
            (flt%uyp(:,:,s)-flt%uym(:,:,s))**2+ &
            (flt%uzp(:,:,s)-flt%uzm(:,:,s))**2 )
    else
       flt%U(:,:,s) = sqrt(flt%Ux(:,:,s)**2+flt%Uy(:,:,s)**2)
    end if

    ! loop over slipping portion of fault and calculate misfit

    do j = 1,mdl%ny
       do i = m,p

          ! calculate strength

          select case(fri%friction)
          case('slipweak')
             st = strength_sw(flt%U(i,j,s),flt%Q(i,j,s),flt%sz(i,j,s), &
                 i,j,(mdl%x(i)**2+mdl%y(j)**2)**0.5,mdl%t,fri%sw)
          case default
             st = zero
          end select
          
          ! evaluate misfit

          if (mdl%bm) then
             call static_elasticity_bm(st, &
                  flt%sx(i,j,s),flt%sy(i,j,s), &
                  flt%fxpi(i,j),flt%fxmi(i,j), &
                  flt%fypi(i,j),flt%fymi(i,j), &
                  flt%uxp(i,j,s),flt%uxm(i,j,s), &
                  flt%uyp(i,j,s),flt%uym(i,j,s), &
                  mdl%transverse_slip,misfit(i,j,:))
          else
             call static_elasticity_im(st, &
                  flt%sx(i,j,s),flt%sy(i,j,s), &
                  flt%Ux(i,j,s),flt%Uy(i,j,s), &
                  mdl%transverse_slip,misfit(i,j,:))
          end if
          
       end do
    end do

    ! pack misfit array into vector

    n = (p-m+1)*mdl%ny

    if (mdl%bm) then
       if (mdl%transverse_slip) then
          g(    1:  n) = reshape(misfit(:,:,1),shape(g(1:n)))
          g(  n+1:2*n) = reshape(misfit(:,:,2),shape(g(1:n)))
          g(2*n+1:3*n) = reshape(misfit(:,:,3),shape(g(1:n)))
          g(3*n+1:4*n) = reshape(misfit(:,:,4),shape(g(1:n)))
       else
          g(    1:  n) = reshape(misfit(:,:,1),shape(g(1:n)))
          g(  n+1:2*n) = reshape(misfit(:,:,2),shape(g(1:n)))
       end if
    else
       if (mdl%transverse_slip) then
          g(  1:  n) = reshape(misfit(:,:,1),shape(g(1:n)))
          g(n+1:2*n) = reshape(misfit(:,:,2),shape(g(1:n)))
       else
          g(1:n) = reshape(misfit(:,:,1),shape(g(1:n)))
       end if
    end if
    
    ! deallocate memory for misfit array
    
    deallocate(misfit)

  end subroutine static_friction


  subroutine solve_friction_combined(mdl,cnv,flt,fri,s,accept)
    ! SOLVE_FRICTION_COMBINED solves the friction law coupled with elasticity
    ! for unknown velocity and stress, treats both tangential and normal components
    ! 
    ! Modified: 5 December 2007

    use constants, only : zero,half,one,two,four
    use model, only : model_type
    use convolution, only : convolution_type
    use fault, only : fault_type
    use friction, only : friction_type
    use fourier, only : strain
    use ratestate, only : check_state_rs!,rate_rs
    use slipweak, only : fric_sw!,evolution_sw
    use io, only : error,warning

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! FLT = fault variables
    ! FRI = friction variables
    ! S = integration stage for storage
    ! ACCEPT = accept solution

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(in) :: cnv
    type(fault_type),intent(inout) :: flt
    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: s
    logical,intent(out) :: accept

    ! Internal Parameters:
    ! P = pore pressure
    ! MUP = (far-field) shear modulus, plus side
    ! MUM = (far-field) shear modulus, minus side
    ! V = slip velocity
    ! O = opening velocity
    ! N = effective normal stress, positive in compression
    ! T = shear stress
    ! I = index in x direction
    ! J = index in y direction
    ! ERR = flag that indicates error in solving friction law

    real(pr),dimension(mdl%mx:mdl%px,mdl%ny) :: p
    real(pr) :: mum,mup,V,O,N,T
    integer(pin) :: i,j
    logical :: err

    ! calculate fault-parallel strains if needed

    if (flt%poroelastic) &
         call strain(mdl,cnv,flt%exxp,flt%exxm,flt%exyp,flt%exym,flt%eyyp,flt%eyym)

    ! return if process holds no x data 

    if (.not.mdl%holds_x) then
       accept = .true.
       return
    end if

    ! calculate contributions to pore pressure change

    p = zero

    ! poroelastic fault zone response

    if (flt%poroelastic) then
    
       ! (far-field) shear moduli of both sides
    
       if (mdl%bm) then
          mum = mdl%mum
          mup = mdl%mup
       else
          mum = mdl%mu
          mup = mdl%mu
       end if

       ! calculate poro pressure change from fault-parallel strains and normal stress

       print *, 'flt%sz not known--fix this routine'

       flt%p0(:,:,s) = -flt%q0*flt%sz(:,:,s)-flt%qp*mup*(flt%exxp+flt%eyyp)-flt%qm*mum*(flt%exxm+flt%eyym)

       ! add this contribution
    
       p = p+flt%p0(:,:,s)

    end if

    ! thermal pressurization

    if (flt%thermpres) p = p+flt%p0(:,:,s)

    ! loop over fault

    do j = 1,mdl%ny
       do i = mdl%mx,mdl%px
    
          ! store necessary information in derived type fri

          if (allocated(flt%T0)) fri%T = flt%T0(i,j,s)
          fri%p = p(i,j)

          select case(mdl%friction_method)
             
          case default

             call error('Invalid friction_method:' // trim(adjustl(mdl%friction_method)) ,'')

          case('strength')
             
             ! check if value of state is acceptable, otherwise reject step
             
             select case(fri%friction)
             case default
                accept = .true.
             case('ratestate')
                accept = check_state_rs(flt%Q(i,j,s),fri%rs)
             end select
             
             ! if unacceptable, reject this step and exit routine
             ! (usually to try again with smaller step size)
             
             if (.not.accept) return
             
             ! calculate friction coefficient, but only if velocity-independent
             
             select case(fri%friction)
             case('slipweak')
                !fri%f0 = fric_sw(flt%U(i,j,s),i,j,mdl%x(i),mdl%t,fri%sw)     
               fri%f0 = fric_sw(flt%U(i,j,s),i,j,(mdl%x(i)**2+mdl%y(j)**2)**0.5,mdl%t,fri%sw)
             case('constfV')
                ! constant f,V: set V=1 and calculate new shear resistance
                ! (assuming f=1 and no resistance in y direction)
                if (mdl%bm) then
                   flt%vxp(i,j,s) =  half
                   flt%vxm(i,j,s) = -half
                   flt%vyp(i,j,s) = zero
                   flt%vym(i,j,s) = zero
                else
                   flt%Vx(i,j,s) = one
                   flt%Vy(i,j,s) = zero
                end if
                flt%sx(i,j,s) = -flt%sz(i,j,s)
                flt%sy(i,j,s) = zero
                cycle
             case default
                fri%f0 = zero
             end select
             
             ! solve elasticity coupled with possibly velocity-dependent friction law;
             
             ! first assume no opening
                
             fri%opening = .false.

             do

                ! lock fault in y direction (for now)

                flt%vyp(i,j,s) = (flt%fypi(i,j)-flt%fymi(i,j))/(mdl%cvsp+mdl%cvsm)
                flt%vym(i,j,s) = flt%vyp(i,j,s)
                flt%sy(i,j,s) = flt%sy0(i,j)+flt%fypi(i,j)-mdl%cvsp*flt%vyp(i,j,s)

                ! set parameters required by friction solver

                fri%i = i
                fri%j = j

                fri%sx0 = flt%sx0(i,j)
                fri%fx = fri%cV*(flt%fxpi(i,j)/mdl%cvsp+flt%fxmi(i,j)/mdl%cvsm)
                fri%sz0 = flt%sz0(i,j)
                fri%fz = fri%cO*(flt%fzpi(i,j)/mdl%cvnp+flt%fzmi(i,j)/mdl%cvnm)

                fri%Q = flt%Q(i,j,s)

                V = flt%V(i,j,s)
                O = flt%O(i,j,s)
                T = flt%S(i,j,s)
                N = flt%N(i,j,s)

                ! perturb initial guess at V if originally zero and using rate-and-state law

                if (fri%friction=='ratestate'.and.V==zero) then
                   if (flt%scale_Vmin<=zero) &
                        call warning('Friction solver may have problems, set scale_Vmin>0','solve_friction_combined')
                   V = flt%scale_Vmin
                end if
                
                ! solve friction

                call couple_elastic_strength_combined(fri, &
                     V,O,T,N,flt%scale_Vmin,flt%scale_Vmax,flt%scale_s,err)
                
                ! if no solution can be found, reject this step and exit routine
                ! (usually to try again with smaller step size)
                
                if (err) then
                   accept = .false.
                   return
                end if
                
                if (.not.fri%opening) then
                   ! no opening permitted, check sign of effective normal stress
                   if (N>=zero) then
                      ! compressive effective normal stress, solution is acceptable
                      exit
                   else
                      ! tensile normal stress, solve again permitting fault opening
                      if (.not.mdl%opening) call error('Tensile normal stress but opening not permitted','solve_friction_combined')
                      fri%opening = .true.
                      cycle
                   end if
                end if

             end do

             ! particle velocities

             flt%vxp(i,j,s) =  (-flt%sx(i,j,s)+flt%sx0(i,j)+flt%fxpi(i,j))/mdl%cvsp
             flt%vxm(i,j,s) = -(-flt%sx(i,j,s)+flt%sx0(i,j)+flt%fxmi(i,j))/mdl%cvsm
             flt%vzp(i,j,s) =  (-flt%sz(i,j,s)+flt%sz0(i,j)+flt%fzpi(i,j))/mdl%cvnp
             flt%vzm(i,j,s) = -(-flt%sz(i,j,s)+flt%sz0(i,j)+flt%fzmi(i,j))/mdl%cvnm

             ! set magnitude of slip and opening velocity, shear and effective normal stress
             
             flt%V(i,j,s) = V
             flt%O(i,j,s) = O
             flt%S(i,j,s) = T
             flt%N(i,j,s) = N

          end select

       end do
    end do

  end subroutine solve_friction_combined


  subroutine couple_elastic_strength_combined(fri,V,O,T,N,scale_Vmin,scale_Vmax,scale_s,err)
    ! COUPLE_ELASTIC_STRENGTH_COMBINED simultaneously solves elasticity and friction
    ! for shear and normal components of motion
    ! 
    ! Modified: 6 August 2010

    use constants, only : zero,one,two,three
    use friction, only : friction_type,elastic_strength_combined
    use utilities, only : newton

    implicit none

    ! I/O Parameters:
    ! FRI = friction variables
    ! V = slip velocity
    ! O = opening velocity
    ! T = shear stress
    ! N = effective normal stress, positive in compression
    ! SCALE_VMIN = scale of velocity, minimum
    ! SCALE_VMAX = scale of velocity, maximum
    ! SCALE_S = scale of stress
    ! ERR = flag that indicates error in solving friction law

    type(friction_type),intent(inout) :: fri
    real(pr),intent(in) :: scale_s,scale_Vmin,scale_Vmax
    real(pr),intent(inout) :: V,O,T,N
    logical,intent(out) :: err

    ! X = independent variables for Newton's method
    ! XSCALE = scale of independent variables for convergence of Newton's method
    ! FSCALE = scale of function for convergence of Newton's method
    ! PARAM = parameters required for function that is solved using Newton's method
    ! SIZE_PARAM = size of param
    ! NMAX = maximum number of iterations allowed to solve equations with Newton's method
    ! TOLX = maximum tolerated error in independent variables
    ! TOLF = maximum tolerated error in dependent variables

    real(pr),dimension(:),allocatable :: x,xscale,fscale
    integer(pin),dimension(:),allocatable :: param
    integer(pin) :: nmax,size_param
    real(pr) :: tolx,tolf
    
    ! set iteration limit and tolerances

    nmax = 1000
    tolx = epsilon(one)
    tolf = epsilon(one)**(two/three)

    allocate(x(2))
    
    x(1) = V
    x(2) = O
    
    allocate(xscale(2),fscale(2))
    
    xscale = min(max(abs(V),scale_Vmin),scale_Vmax)
    
    if (scale_s==zero) then
       fscale = abs(T)
    else
       fscale = scale_s
    end if
    
    size_param = size(transfer(fri,param))
    allocate(param(size_param))
    param = transfer(fri,param)
    
    call newton(x,nmax,tolx,tolf,param, &
         elastic_strength_combined,xscale,fscale,err)
    
    V = x(1)
    O = x(2)

    T =   fri%sx0+fri%fx-fri%cV*V
    N = -(fri%sz0+fri%fz-fri%cO*O)-fri%p
    
    deallocate(x)
    deallocate(xscale)
    deallocate(fscale)
    deallocate(param)

  end subroutine couple_elastic_strength_combined


end module friction_routines
