! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module friction

  ! FRICTION contains variables and routines related to applying
  ! friction laws on the fault surface, coupling slip, slip velocity,
  ! state variable, and stress fields
  ! 
  ! Modified: 6 August 2010

  use constants, only : pr,pin
  use slipweak, only : slipweak_type
  use ratestate, only : ratestate_type

  implicit none

  ! FRICTION_TYPE is a derived type containing variables related to friction law
  !
  ! Parameters:
  ! FRICTION = type of friction to be applied
  ! SW = slip-weakening variables
  ! RS = rate-and-state variables
  ! SLIP_X = extent of slipping region in x direction
  ! SLIP_X0 = center of slipping region in x direction
  ! (for these derived types, memory is only allocated if needed)
  ! The following variables are used for temporary storage when using Newton's method:
  ! S0 = velocity-independent strength
  ! F0 = velocity-independent friction coefficient
  ! Q = state variable
  ! SZ = stress in z-direction
  ! CV = shear radiation damping coefficient
  ! CO = normal radiation damping coefficient
  ! SFX = sum of load and stress transfer in x direction (sfx = sx0+fx)
  ! SFY = sum of load and stress transfer in y direction (sfy = sy0+fy)
  ! SFZ = sum of load and stress transfer in y direction (sfy = sy0+fy)
  ! SX0 = xz component of load
  ! SZ0 = zz component of load
  ! FX = xz component of stress transfer 
  ! FZ = zz component of stress transfer
  ! VO = velocity at previous time step
  ! SO = strength at previous time step
  ! DVO = contribution to rate of change of velocity from previous integration stages
  ! DSO = contribution to rate of change of strength from previous integration stages
  ! A = Runge-Kutta coefficient of current stage
  ! DT = time step
  ! T = temperature
  ! P = pore pressure change
  ! I = index in x direction
  ! J = index in y direction
  ! DIRECT_ONLY = solve for direct effect only
  ! OPENING = permit opening
  
  type friction_type
     character(64) :: friction
     type(slipweak_type),pointer :: sw=>null()
     type(ratestate_type),pointer :: rs=>null()
     logical :: direct_only,opening
     real(pr) :: slip_x,slip_x0,S0,f0,Q,sz,cV,cO,sfx,sfy,sx0,sz0,fx,fz,Vo,So,dVo,dSo,a,dt,T,p
     integer(pin) :: i,j
  end type friction_type

contains


  subroutine read_friction(ninput,mdl,flt,fri)
    ! READ_FRICTION reads in friction variables from file
    ! 
    ! Modified: 20 July 2010

    use constants, only : zero,one
    use io, only : error
    use model, only : model_type
    use fault, only : fault_type
    use slipweak, only : read_slipweak
    use ratestate, only : read_ratestate

    implicit none

    ! I/O Parameters:
    ! NINPUT = unit number for *.in input file
    ! MDL = model variables
    ! FLT = fault variables
    ! FRI = friction variables

    integer(pin),intent(in) :: ninput
    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(friction_type),intent(out) :: fri

    ! Internal Parameters:
    ! FRICTION = type of friction to be applied
    ! SLIP_X = extent of slipping region in x direction
    ! SLIP_X0 = center of slipping region in x direction
    ! STAT = I/O error flag

    character(64) :: friction
    real(pr) :: slip_x,slip_x0
    integer(pin) :: stat

    ! make namelist of user input variables

    namelist /friction_list/ friction

    ! set default values
    
    friction = 'slipweak'
    slip_x = huge(one)
    slip_x0 = zero
    
    ! read namelist from input file, call error routine if needed
    
    rewind(ninput)
    read(ninput,nml=friction_list,iostat=stat)
    if (stat/=0) call error("Error reading namelist 'friction_list' in .in file",'read_friction')
    
    ! assign input variables to components of derived type
    
    fri%friction = friction
    fri%slip_x = slip_x
    fri%slip_x0 = slip_x0

    ! read specific friction law, after checking if this is a valid law;
    ! also set whether or not state and temperature fields are required by friction law
    
    select case(fri%friction)
    case default
       call error("Invalid parameter 'friction' in namelist 'friction_list' in .in file",'read_friction')
    case('none','constfV')
       flt%state = .false.
       flt%temperature = .false.
    case('slipweak')
       allocate(fri%sw)
       call read_slipweak(ninput,fri%sw)
       flt%state = .false.
       flt%temperature = .false.
    case('ratestate')
       allocate(fri%rs)
       call read_ratestate(ninput,fri%rs,flt%temperature)
       select case(mdl%friction_method)
       case('strength')
          flt%state = .true.
       case('strength_rate')
          flt%state = .false.
       end select
    end select

  end subroutine read_friction


  subroutine init_friction(necho,mdl,fri,flt)
    ! INIT_FRICTION initializes friction variables
    ! 
    ! Modified: 3 December 2007

    use constants, only : one,two
    use model, only : model_type
    use fault, only : fault_type
    use io, only : write_matlab
    use slipweak, only : init_slipweak,set_strength_sw
    use ratestate, only : init_ratestate
    use mpi_routines, only : is_master

    implicit none

    ! I/O Parameters:
    ! NECHO = unit number for output file
    ! MDL = model variables
    ! FRI = friction variables

    integer(pin),intent(in) :: necho
    type(model_type),intent(inout) :: mdl
    type(friction_type),intent(inout) :: fri
    type(fault_type),intent(inout) :: flt

    ! Internal Parameters:
    ! I = index in x direction
    ! J = index in y direction
    ! STAGE = integration stage

    integer(pin) :: i,j,stage

    ! set default absolute temperature (used by some rate-and-state laws) to room temperature

    fri%T = 293._pr

    ! set parameters required by friction solver (using material properties, etc.)

    fri%cV = mdl%cvsp*mdl%cvsm/(mdl%cvsp+mdl%cvsm)
    fri%cO = mdl%cvnp*mdl%cvnm/(mdl%cvnp+mdl%cvnm)

    ! initialize specific friction law

    select case(fri%friction)

    case('slipweak')

       call init_slipweak(necho,mdl,fri%sw)

       ! set fault strength in regularized friction law

       if (.not.mdl%holds_x) return

       do j = 1,mdl%ny
          do i = mdl%mx,mdl%px
             do stage = 1,mdl%ns
                flt%Q(i,j,stage) = set_strength_sw(flt%U(i,j,stage),flt%sz(i,j,stage),i,j,(mdl%x(i)**2+mdl%y(j)**2)**0.5,mdl%t,fri%sw)
             end do
          end do
       end do

    case('ratestate')

       call init_ratestate(necho,mdl,fri%rs)

    end select

    ! output variable values into matlab file

    if (is_master) then
       
       call write_matlab(necho,'friction',fri%friction,'fri')
       call write_matlab(necho,'slip_x',fri%slip_x,'fri')
       call write_matlab(necho,'slip_x0',fri%slip_x0,'fri')
       
    end if

  end subroutine init_friction


  subroutine destroy_friction(fri)
    ! DESTROY_FRICTION destroys derived type fri
    ! 
    ! Modified: 28 May 2007

    use slipweak, only : destroy_slipweak
    use ratestate, only : destroy_ratestate

    implicit none
    
    ! I/O Parameters:
    ! FRI = friction variables    
    
    type(friction_type),intent(inout) :: fri
    
    ! destroy derived types for specific friction law and
    ! deallocate memory associated with them

    select case(fri%friction)
    case('slipweak')
       call destroy_slipweak(fri%sw)
       if (associated(fri%sw)) deallocate(fri%sw)
       fri%sw => null()
    case('ratestate')
       call destroy_ratestate(fri%rs)
       if (associated(fri%rs)) deallocate(fri%rs)
       fri%rs => null()
    end select
    
  end subroutine destroy_friction


  subroutine strengthV(fri,V,SV,dSdV,err)
    ! STRENGTHV returns the velocity-dependent strength and derivative
    ! of strength with respect to slip velocity at a point
    ! 
    ! Modified: 19 September 2007

    use constants, only : zero
    use ratestate, only : strength_rs

    implicit none

    ! I/O Parameters:
    ! FRI = friction variables
    ! V = slip velocity
    ! SV = velocity-dependent strength
    ! DSDV = derivative of strength with respect to slip velocity
    ! ERR = error flag (true if strength isn't defined at particular V)

    type(friction_type),intent(in) :: fri
    real(pr),intent(in) :: V
    real(pr),intent(out) :: SV,dSdV
    logical,intent(out) :: err

    select case(fri%friction)
    case default
       SV = zero
       dSdV = zero
       err = .false.
    case('ratestate')
       call strength_rs(fri%rs,V,fri%Q,fri%T,fri%sz,fri%i,fri%j,SV,dSdV,err)
    end select

  end subroutine strengthV


  subroutine fricV(fri,V,f,dfdV,err)
    ! FRICV returns the velocity-dependent friction coefficient and derivative
    ! of friction coefficient with respect to slip velocity at a point
    ! 
    ! Modified: 19 October 2007

    use constants, only : zero
    use ratestate, only : fric_rs

    implicit none

    ! I/O Parameters:
    ! FRI = friction variables
    ! V = slip velocity
    ! F = velocity-dependent friction coefficient
    ! DFDV = derivative of friction coefficient with respect to slip velocity
    ! ERR = error flag (true if friction isn't defined at particular V)

    type(friction_type),intent(in) :: fri
    real(pr),intent(in) :: V
    real(pr),intent(out) :: f,dfdV
    logical,intent(out) :: err

    select case(fri%friction)
    case default
       f = zero
       dfdV = zero
       err = .false.
    case('ratestate')
       call fric_rs(fri%rs,V,fri%Q,fri%T,fri%i,fri%j,f,dfdV,err)
    end select
    
  end subroutine fricV
  

  subroutine elastic_strength(x,f,param,err,jac)
    ! ELASTIC_STRENGTH solves elasticity equations for shear stress and slip velocity given
    ! shear strength.  This version applies to both scalar and vector slip (in which 
    ! shear stress opposes the slip velocity in direction).  The equations are:
    ! (1) S = St(V)
    ! (2) sx = sx0+fx-cV*Vx
    ! (3) sy = sy0+fy-cV*Vy
    ! (4) angle_V = angle_S
    ! where (2) and (3) are the elastodynamic equations, (1) says that the magnitude of 
    ! the shear stress is equal to the strength, and (4) says that sliding occurs in the same
    ! direction as the shear stress.  There are either two (V,S) or four unknowns (V,S,angle_V,angle_S), 
    ! which are solved for using Newton's method.
    ! 
    ! Modified: 6 August 2010

    use constants, only : zero,one

    implicit none
    
    ! I/O Parameters:
    ! X = independent variable for Newton's method
    ! F = dependent variable for Newton's method
    ! PARAM = variables required for evaluating elastodynamic equation and strength
    ! ERR = error flag (true if error)
    ! JAC = Jacobian df/dx

    real(pr),dimension(:),intent(in) :: x
    real(pr),intent(out) :: f(size(x))
    integer(pin),dimension(:),intent(in) :: param
    logical,intent(out) :: err
    real(pr),intent(out),optional :: jac(size(x),size(x))

    ! Internal Parameters:
    ! VX = slip velocity in x direction
    ! VY = slip velocity in y direction
    ! V = magnitude of slip velocity
    ! SX = stress in x direction
    ! SY = stress in y direction
    ! S = magnitude of shear stress
    ! ANGLE_S = angle of shear stress vector (w.r.t. x-direction)
    ! ANGLE_V = angle of slip velocity vector (w.r.t. x-direction)
    ! ST = shear strength
    ! SV = velocity-dependent strength
    ! DSTDV = derivative of shear strength with respect to slip velocity
    ! FRI = friction variables

    real(pr) :: Vx,Vy,V,sx,sy,S,angle_S,angle_V,St,SV,dStdV
    type(friction_type) :: fri

    fri = transfer(param,fri)

    ! determine if scalar or vector slip, solve appropriate equations

    select case(size(x))

    case(2) ! scalar

       V = x(1)
       S = x(2)

       call strengthV(fri,V,SV,dStdV,err)

       if (err) return

       St = fri%S0+SV

       f(1) = -S+St
       f(2) = -S+fri%sfx-fri%cV*V

       if (present(jac)) then
          
          jac(1,1) = dStdV
          jac(1,2) = -one
          
          jac(2,1) = -fri%cV
          jac(2,2) = -one
          
       end if
    
    case(4) ! vector
       
       V = x(1)
       S = x(2)
       angle_V = x(3)
       angle_S = x(4)

       Vx = V*cos(angle_V)
       Vy = V*sin(angle_V)
       sx = S*cos(angle_S)
       sy = S*sin(angle_S)

       call strengthV(fri,V,SV,dStdV,err)

       if (err) return

       St = fri%S0+SV

       f(1) = -S+St
       f(2) = -sx+fri%sfx-fri%cV*Vx
       f(3) = -sy+fri%sfy-fri%cV*Vy
       f(4) = -angle_V+angle_S
       
       if (present(jac)) then
          
          jac(1,1) = dStdV
          jac(1,2) = -one
          jac(1,3) = zero
          jac(1,4) = zero
          
          jac(2,1) = -fri%cV*Vx/V
          jac(2,2) = -sx/S
          jac(2,3) = fri%cV*Vy
          jac(2,4) = sy
          
          jac(3,1) = -fri%cV*Vy/V
          jac(3,2) = -sy/S
          jac(3,3) = -fri%cV*Vx
          jac(3,4) = -sx
          
          jac(4,1) = zero
          jac(4,2) = zero
          jac(4,3) = -one
          jac(4,4) = one
          
       end if

    end select

  end subroutine elastic_strength


  subroutine elastic_strength_rate(x,f,param,err,jac)
    ! ELASTIC_RATE solves elasticity and friction equations for shear stress and velocity.
    ! This version applies to both scalar and vector slip (in which shear stress opposes
    ! slip velocity in direction).  The equations are:
    ! (1) dS/dt = D(V,S)*dV/dt+R(V,S) with
    !     S = abs(sx) OR S = sqrt(sx**2+sy**2)
    !     V = abs(Vx) OR V = sqrt(vx**2+vy**2)
    !     dS/dt = [(S-So)/dt-dSo]/a from implicit RK formula S = So+dt*(dSo+a*dS/dt)
    !     dV/dt = [(V-Vo)/dt-dVo]/a (from implicit RK formula)
    ! (2) sx = sx0+fx-cV*Vx
    ! (3) sy = sy0+fy-cV*Vy
    ! (4) Vx*sy = Vy*sx
    ! where (2) and (3) are the elastodynamic equations, (1) is the time derivative of the 
    ! friction law and (4) says that shear stress opposes slip velocity.  There are either 
    ! two (Vx,sx) or four unknowns (Vx,Vy,sx,sy), which are solved for using Newton's method.
    ! 
    ! Modified: 6 August 2010

    use constants, only : zero,one,two
    use ratestate, only : directeffect_rs,directeffect_log_rs,rate_rs

    implicit none
    
    ! I/O Parameters:
    ! X = independent variable for Newton's method
    ! F = dependent variable for Newton's method
    ! PARAM = variables required for evaluating elastodynamic equation and strength
    ! ERR = error flag (true if error)
    ! JAC = Jacobian df/dx

    real(pr),dimension(:),intent(in) :: x
    real(pr),intent(out) :: f(size(x))
    integer(pin),dimension(:),intent(in) :: param
    logical,intent(out) :: err
    real(pr),intent(out),optional :: jac(size(x),size(x))

    ! Internal Parameters:
    ! VX = slip velocity in x direction
    ! VY = slip velocity in y direction
    ! V = slip velocity (current time step)
    ! VO = slip velocity (previous time step)
    ! DV = rate of change of slip velocity (times dt)
    ! DVO = rate of change of slip velocity (contribution from previous stages) (times dt)
    ! SX = stress in x direction
    ! SY = stress in y direction
    ! S = strength (current time step)
    ! SO = strength (previous time step)
    ! DS = rate of change of strength (times dt)
    ! DSO = rate of change of strength (contribution from previous stages) (times dt)
    ! DT = time step
    ! D = direct effect (coefficient of dV/dt)
    ! R = rate (how strength evolves, excluding direct effect)
    ! A = coefficient of current rate
    ! DDDV = derivative of D(V,S) with respect to V
    ! DDDS = derivative of D(V,S) with respect to S
    ! DRDV = derivative of R(V,S) with respect to V
    ! DRDS = derivative of R(V,S) with respect to S
    ! DVDS = derivative of V with respect to t
    ! DIRECT_ONLY = solve for direct effect only
    ! DLOGV = rate of change of log(slip velocity) (times dt)
    ! FRI = friction variables

    logical :: direct_only
    real(pr) :: Vx,Vy,V,dV,Vo,dVo,sx,sy,S,dS,So,dSo,dt,D,R,a,dDdV,dDdS,dRdV,dRdS,dlogV
    type(friction_type) :: fri

    print *, 'STRENGTH-RATE SOLVER NEEDS TO BE UPDATED FOR VECTOR FRICTION'

    fri = transfer(param,fri)

    dt = fri%dt
    Vo = fri%Vo
    So = fri%So
    dVo = fri%dVo
    dSo = fri%dSo
    a = fri%a
    direct_only = fri%direct_only

    ! solve direct effect only: dS = D(V,S)*dlogV, where D(V,S) is coefficient of dlogV/dt
    ! instead of dV/dt in this case

    if (direct_only) then

       Vx = x(1)
       sx = x(2)
       
       V = Vx ! assumes positive Vx
       S = abs(sx)
       
       call directeffect_log_rs(fri%rs,V,S,fri%T,fri%sz,fri%i,fri%j,D,dDdV,dDdS)
       
       dlogV = log(V/Vo) ! = log(V)-log(Vo)
       dS = S-So
       
       f(1) = -dS+D*dlogV
       f(2) = -sx+fri%sfx-fri%cV*Vx
       
       if (present(jac)) then
          
          jac(1,1) =    D/V+dDdV*dlogV
          jac(1,2) = -one+dDdS*dlogV
          
          jac(2,1) = -fri%cV
          jac(2,2) = -one
          
       end if
       
       err = .false.

       return

    end if

    ! determine if scalar or vector slip, solve appropriate equations

    select case(size(x))

    case(2) ! scalar

       Vx = x(1)
       sx = x(2)

       V = Vx ! assumes positive Vx
       S = abs(sx)

       call directeffect_rs(fri%rs,V,S,fri%T,fri%sz,fri%i,fri%j,D,dDdV,dDdS)
       call rate_rs(fri%rs,V,S,fri%T,fri%sz,fri%i,fri%j,R,dRdV,dRdS,err)

       if (err) return

       dV = (V-Vo-dt*dVo)/a
       dS = (S-So-dt*dSo)/a

       f(1) = -dS+D*dV+dt*R
       f(2) = -sx+fri%sfx-fri%cV*Vx

       if (present(jac)) then
          
          jac(1,1) =    D/a+dDdV*dV+dt*dRdV
          jac(1,2) = -one/a+dDdS*dV+dt*dRdS
          
          jac(2,1) = -fri%cV
          jac(2,2) = -one
          
       end if
    
    case(4) ! vector
       
       Vx = x(1)
       sx = x(2)
       Vy = x(3)
       sy = x(4)

       V = sqrt(Vx**2+Vy**2)
       S = sqrt(sx**2+sy**2)

       call directeffect_rs(fri%rs,V,S,fri%T,fri%sz,fri%i,fri%j,D,dDdV,dDdS)
       call rate_rs(fri%rs,V,S,fri%T,fri%sz,fri%i,fri%j,R,dRdV,dRdS,err)

       if (err) return

       dV = (V-Vo-dt*dVo)/a
       dS = (S-So-dt*dSo)/a

       f(1) = -dS+D*dV+dt*R
       f(2) = -sx+fri%sfx-fri%cV*Vx
       f(3) = -sy+fri%sfy-fri%cV*Vy
       f(4) = -Vy*sx+Vx*sy
       
       if (present(jac)) then
          
          jac(1,1) =    (D/a+dDdV*dV+dt*dRdV)*Vx/V
          jac(1,3) =    (D/a+dDdV*dV+dt*dRdV)*Vy/V
          jac(1,2) = (-one/a+dDdS*dV+dt*dRdS)*sx/S
          jac(1,4) = (-one/a+dDdS*dV+dt*dRdS)*sy/S
          
          jac(2,1) = -fri%cV
          jac(2,2) = -one
          jac(2,3) = zero
          jac(2,4) = zero
          
          jac(3,1) = zero
          jac(3,2) = zero
          jac(3,3) = -fri%cV
          jac(3,4) = -one
          
          jac(4,1) = sy
          jac(4,2) = -Vy
          jac(4,3) = -sx
          jac(4,4) = Vx
          
       end if

    end select

  end subroutine elastic_strength_rate


  subroutine elastic_strength_combined(x,R,param,err,jac)
    ! ELASTIC_STRENGTH_COMBINED solves elasticity equations for shear stress, T, normal stress, N,
    ! slip velocity, V, and opening velocity, O, given
    ! shear strength. The equations are:
    ! (1) T =  (sx0+fx-cV*V)+m*[2*b*(sz0+fz-cO*O)+cV*O-A]
    ! (2) N = -(sz0+fz-cO*O)+m*[2*  (sx0+fx-cV*V)+cO*V  ]
    ! (3) T = f(V)*N
    ! (4) N*O = 0, N>0, O>0 => O=0 if N>=0 or N=0 if O>=0
    ! where (1) and (2) are the elastodynamic equations, (3) is the friction law, and
    ! (4) prevents opening unless normal stress would otherwise become tensile. V and O
    ! are solved for using Newton's method.
    ! 
    ! Modified: 6 August 2010

    use constants, only : zero,one,two

    implicit none
    
    ! I/O Parameters:
    ! X = independent variable for Newton's method
    ! R = residual (function to be set equal to zero) for Newton's method
    ! PARAM = variables required for evaluating elastodynamic equation and strength
    ! ERR = error flag (true if error)
    ! JAC = Jacobian dR/dx

    real(pr),dimension(:),intent(in) :: x
    real(pr),intent(out) :: R(size(x))
    integer(pin),dimension(:),intent(in) :: param
    logical,intent(out) :: err
    real(pr),intent(out),optional :: jac(size(x),size(x))

    ! Internal Parameters:
    ! V = slip velocity
    ! O = opening velocity
    ! F = friction coefficient
    ! FV = velocity-dependent friction coefficient
    ! DFDV = derivative of friction coefficient with respect to slip velocity
    ! T = stress stress
    ! N = effective normal stress, positive in compression
    ! DTDV = derivative of T with respect to V
    ! DTDO = derivative of T with respect to O
    ! DNDV = derivative of N with respect to V
    ! DNDO = derivative of N with respect to O
    ! FRI = friction variables

    real(pr) :: V,O,f,fV,dfdV,T,N,dTdV,dTdO,dNdV,dNdO
    type(friction_type) :: fri

    fri = transfer(param,fri)

    V = x(1)
    O = x(2)
    
    call fricV(fri,V,fV,dfdV,err)
    
    if (err) return
    
    f = fri%f0+fV
    
    ! elastic relations and partial derivatives

    T =   fri%sx0+fri%fx-fri%cV*V
    N = -(fri%sz0+fri%fz-fri%cO*O)-fri%p

    dTdV = -fri%cV
    dNdO =  fri%cO

    dTdO = zero
    dNdV = zero

    ! equation 1: R(V,O) = -T(V,O)+f(V)*N(V,O)

    R(1) = -T+f*N

    ! equation 2: R(V,O) = N(V,O) or O

    if (fri%opening) then
       R(2) = N
    else
       R(2) = O
    end if
    
    if (present(jac)) then
       
       jac(1,1) = -dTdV+dfdV*N+f*dNdV
       jac(1,2) = -dTdO       +f*dNdO
       
       if (fri%opening) then
          jac(2,1) = dNdV
          jac(2,2) = dNdO
       else       
          jac(2,1) = zero
          jac(2,2) = one
       end if

    end if
    
  end subroutine elastic_strength_combined


end module friction
