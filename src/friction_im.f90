! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module friction_im

  ! FRICTION_IM contains routines to couple friction laws and elasticity for identical materials
  ! 
  ! Modified: 6 August 2010

  use constants, only : pr,pin

  implicit none


contains


  subroutine couple_elastic_strength_im(fri,i,j,S0,Q,sz,Vx,Vy,scale_Vmin,scale_Vmax, &
       sx,sy,scale_s,sx0,sy0,fx,fy,cV,transverse_slip,err,info)
    ! COUPLE_ELASTIC_STRENGTH_IM solves elasticity equation with stress<=strength
    ! 
    ! Modified: 6 August 2010

    use constants, only : zero,one,two,three
    use friction, only : friction_type,elastic_strength
    use utilities, only : newton

    implicit none

    ! I/O Parameters:
    ! FRI = friction variables
    ! I = index in x direction
    ! J = index in y direction
    ! S0 = portion of shear strength independent of slip velocity
    ! Q = state variable
    ! SZ = stress in z direction
    ! VX = slip velocity in x direction
    ! VY = slip velocity in y direction
    ! SCALE_VMIN = scale of velocity, minimum
    ! SCALE_VMAX = scale of velocity, maximum
    ! SX = stress in x direction
    ! SY = stress in y direction
    ! SCALE_S = scale of stress
    ! SX0 = load in x direction
    ! SY0 = load in y direction
    ! FX = stress transfer in x direction
    ! FY = stress transfer in y direction
    ! CV = coefficient of radiation damping
    ! TRANSVERSE_SLIP = scalar or vector slip
    ! ERR = flag that indicates error in solving friction law
    ! INFO = print details of Newton's method solution

    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: i,j
    real(pr),intent(in) :: Q,S0,sz,sx0,sy0,fx,fy,cV,scale_s,scale_Vmin,scale_Vmax
    real(pr),intent(inout) :: Vx,Vy,sx,sy
    logical,intent(in) :: transverse_slip,info
    logical,intent(out) :: err

    ! Internal Parameters:
    ! SL = stress if fault is locked
    ! NEWTON_X = independent variables for Newton's method
    ! NEWTON_XSCALE = scale of independent variables for convergence of Newton's method
    ! NEWTON_FSCALE = scale of function for convergence of Newton's method
    ! PARAM = parameters required for function that is solved using Newton's method
    ! SIZE_PARAM = size of param
    ! NEWTON_NMAX = maximum number of iterations allowed to solve equations with Newton's method
    ! NEWTON_TOLX = maximum tolerated error in independent variables
    ! NEWTON_TOLF = maximum tolerated error in dependent variables

    real(pr) :: sl
    real(pr),dimension(:),allocatable :: newton_x,newton_xscale,newton_fscale
    integer(pin),dimension(:),allocatable :: param
    integer(pin) :: newton_nmax,size_param
    real(pr) :: newton_tolx,newton_tolf
   
    ! no error

    err = .false.

    ! set iteration limit and tolerances

    newton_nmax = 1000
    newton_tolx = epsilon(one)
    newton_tolf = epsilon(one)**(two/three)
 
    ! store variables needed elsewhere

    fri%i = i
    fri%j = j
    fri%S0 = S0
    fri%Q = Q
    fri%sz = sz
    fri%cV = cV
    fri%sfx = sx0+fx
    fri%sfy = sy0+fy

    if (transverse_slip) then

       ! assume fault is locked, calculate stress
       
       sl = sqrt((sx0+fx)**2+(sy0+fy)**2)
    
       ! compare strength to elastic shear stress if fault were locked,
       ! allow failure if stress exceeds strength
          
       if (fri%friction/='ratestate'.and.sl<S0) then
          
          ! stress less than strength, fault locked
          
          Vx = zero
          sx = sx0+fx
          Vy = zero
          sy = sy0+fy
          
       else
          
          ! stress set to strength, fault slips;
          ! solution via Newton's method--see the subroutine elastic_strength for further description
          
          allocate(newton_x(4))

          newton_x(1) = max(sqrt(Vx**2+Vy**2),epsilon(one))
          newton_x(2) = sqrt(sx**2+sy**2)
          newton_x(3) = atan2(Vy,Vx)
          newton_x(4) = atan2(sy,sx)
          
          allocate(newton_xscale(4))
          
          newton_xscale(1) = min(max(abs(newton_x(1)),scale_Vmin),scale_Vmax)
          
          if (scale_s==zero) then
             newton_xscale(2) = newton_x(2)
          else
             newton_xscale(2) = scale_s
          end if
          
          newton_xscale(3) = one
          newton_xscale(4) = one
          
          allocate(newton_fscale(4))
          
          newton_fscale(1:3) = newton_xscale(2)
          newton_fscale(4) = one
       
          size_param = size(transfer(fri,param))
          allocate(param(size_param))
          param = transfer(fri,param)

          call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
               elastic_strength,newton_xscale,newton_fscale,err)
          
          Vx = newton_x(1)*cos(newton_x(3))
          Vy = newton_x(1)*sin(newton_x(3))
          sx = newton_x(2)*cos(newton_x(4))
          sy = newton_x(2)*sin(newton_x(4))
       
          deallocate(newton_x)
          deallocate(newton_xscale)
          deallocate(newton_fscale)
          deallocate(param)

       end if

    else

       ! lock fault in y-direction
       
       Vy = zero
       sy = sy0+fy
       
       ! assume fault is locked, calculate stress
       
       sl = abs(sx0+fx)
       
       ! compare strength to elastic shear stress if fault were locked,
       ! allow failure if stress exceeds strength
       
       if (fri%friction/='ratestate'.and.sl<S0) then
          
          ! stress less than strength, fault locked
          
          Vx = zero
          sx = sx0+fx
          
       else
          
          ! stress set to strength, fault slips;
          ! solution via Newton's method--see the subroutine elastic_strength for further description
          
          allocate(newton_x(2))

          newton_x(1) = Vx
          newton_x(2) = sx
          
          allocate(newton_xscale(2))

          newton_xscale(1) = min(max(abs(newton_x(1)),scale_Vmin),scale_Vmax)

          if (scale_s==zero) then
             newton_xscale(2) = abs(sx)
          else
             newton_xscale(2) = scale_s
          end if

          allocate(newton_fscale(2))

          newton_fscale = newton_xscale(2)

          size_param = size(transfer(fri,param))
          allocate(param(size_param))
          param = transfer(fri,param)
       
          if (info) then
             print *, 'before'
             print *, newton_x
             print *, newton_xscale
             print *, newton_fscale
          end if

          call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
               elastic_strength,newton_xscale,newton_fscale,err)

          if (info) then
             print *, 'after, err=',err
             print *, newton_x
          end if

          Vx = newton_x(1)
          sx = newton_x(2)

          deallocate(newton_x)
          deallocate(newton_xscale)
          deallocate(newton_fscale)
          deallocate(param)

       end if
       
    end if

  end subroutine couple_elastic_strength_im


  subroutine couple_elastic_strength_rate_im(fri,i,j,sz,dt,a,Vo,So,dVo,dSo,Vx,Vy,scale_Vmin,scale_Vmax, &
       sx,sy,scale_s,sx0,sy0,fx,fy,cV,transverse_slip,direct_only,err)
    ! COUPLE_ELASTIC_STRENGTH_RATE_IM solves elasticity equation with implicitly discretized friction law
    ! 
    ! Modified: 6 August 2010

    use constants, only : zero,one,two,three
    use friction, only : friction_type,elastic_strength_rate
    use utilities, only : newton

    implicit none

    ! I/O Parameters:
    ! FRI = friction variables
    ! I = index in x direction
    ! J = index in y direction
    ! SZ = stress in z direction
    ! DT = time step
    ! A = coefficient of current rate
    ! VO = slip velocity at previous time step
    ! SO = strength at previous time step
    ! DVO = contribution to rate of change of velocity from previous integration stages
    ! DSO = contribution to rate of change of strength from previous integration stages
    ! VX = slip velocity in x direction
    ! VY = slip velocity in y direction
    ! SCALE_VMIN = scale of velocity, minimum
    ! SCALE_VMAX = scale of velocity, maximum
    ! SX = stress in x direction
    ! SY = stress in y direction
    ! SCALE_S = scale of stress
    ! SX0 = load in x direction
    ! SY0 = load in y direction
    ! FX = stress transfer in x direction
    ! FY = stress transfer in y direction
    ! CV = coefficient of radiation damping
    ! TRANSVERSE_SLIP = scalar or vector slip
    ! DIRECT_ONLY = solve for response by direct effect only
    ! ERR = flag that indicates error in solving friction law

    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: i,j
    real(pr),intent(in) :: sz,dt,a,Vo,So,dVo,dSo,sx0,sy0,fx,fy,cV,scale_s,scale_Vmin,scale_Vmax
    real(pr),intent(inout) :: Vx,Vy,sx,sy
    logical,intent(in) :: transverse_slip,direct_only
    logical,intent(out) :: err

    ! Internal Parameters:
    ! NEWTON_X = independent variables for Newton's method
    ! NEWTON_XSCALE = scale of independent variables for convergence of Newton's method
    ! NEWTON_FSCALE = scale of function for convergence of Newton's method
    ! PARAM = parameters required for function that is solved using Newton's method
    ! SIZE_PARAM = size of param
    ! NEWTON_NMAX = maximum number of iterations allowed to solve equations with Newton's method
    ! NEWTON_TOLX = maximum tolerated error in independent variables
    ! NEWTON_TOLF = maximum tolerated error in dependent variables
    ! MSG = error message
    ! N = iteration on initial guess

    real(pr),dimension(:),allocatable :: newton_x,newton_xscale,newton_fscale
    integer(pin),dimension(:),allocatable :: param
    integer(pin) :: newton_nmax,n,size_param
    real(pr) :: newton_tolx,newton_tolf
    character(64) :: msg

    ! no error

    err = .false.

    ! set iteration limit and tolerances

    newton_nmax = 1000
    newton_tolx = epsilon(one)
    newton_tolf = epsilon(one)**(two/three)
 
    ! store variables needed elsewhere

    fri%i = i
    fri%j = j
    fri%sz = sz
    fri%cV = cV
    fri%sfx = sx0+fx
    fri%sfy = sy0+fy
    fri%Vo = Vo
    fri%So = So
    fri%dVo = dVo
    fri%dSo = dSo
    fri%dt = dt
    fri%a = a
    fri%direct_only = direct_only

    if (transverse_slip) then

       ! stress set to strength, fault slips;
       ! solution via Newton's method--see the subroutine elastic_strength_rate for further description
       
       allocate(newton_x(4))
       
       newton_x(1) = max(sqrt(Vx**2+Vy**2),epsilon(one))
       newton_x(2) = sqrt(sx**2+sy**2)
       newton_x(3) = atan2(Vy,Vx)
       newton_x(4) = atan2(sy,sx)
       
       allocate(newton_xscale(4))
       
       newton_xscale(1) = min(max(abs(newton_x(1)),scale_Vmin),scale_Vmax)

       if (scale_s==zero) then
          newton_xscale(2) = newton_x(2)
       else
          newton_xscale(2) = scale_s
       end if
       
       newton_xscale(3) = one
       newton_xscale(4) = one

       allocate(newton_fscale(4))
       
       newton_fscale(1:3) = newton_xscale(2)
       newton_fscale(4) = one
       
       size_param = size(transfer(fri,param))
       allocate(param(size_param))
       param = transfer(fri,param)
       
       call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
            elastic_strength_rate,newton_xscale,newton_fscale,err)
       
       Vx = newton_x(1)*cos(newton_x(3))
       Vy = newton_x(1)*sin(newton_x(3))
       sx = newton_x(2)*cos(newton_x(4))
       sy = newton_x(2)*sin(newton_x(4))
       
       deallocate(newton_x)
       deallocate(newton_xscale)
       deallocate(newton_fscale)
       deallocate(param)
       
    else
       
       ! lock fault in y-direction
       
       Vy = zero
       sy = sy0+fy
       
       ! stress set to strength, fault slips;
       ! solution via Newton's method--see the subroutine elastic for further description
       
       allocate(newton_x(2))
       
       newton_x(1) = Vx
       newton_x(2) = sx
       
       allocate(newton_xscale(2))
       
       newton_xscale(1) = min(max(abs(newton_x(1)),scale_Vmin),scale_Vmax)
       
       if (scale_s==zero) then
          newton_xscale(2) = abs(sx)
       else
          newton_xscale(2) = scale_s
       end if
       
       allocate(newton_fscale(2))
       
       if (scale_s==zero) then
          newton_fscale(1) = max(abs(fx),abs(cV*Vx))
          newton_fscale(2) = newton_fscale(1)
       else
          newton_fscale = scale_s
       end if
       
       size_param = size(transfer(fri,param))
       allocate(param(size_param))
       param = transfer(fri,param)
       
       n = 0
       do
          
          n = n+1
          
          call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
               elastic_strength_rate,newton_xscale,newton_fscale,err,msg)
          
          if (err) then
             ! repeat Newton's method with different initial guess
             if (n>=100) exit
             print *, msg
             print *, i,n
             print *, Vx,sx
             print *, newton_x,'(failed)'
             if (.not.direct_only) fri%direct_only = .true.
             call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
                  elastic_strength_rate,newton_xscale,newton_fscale,err,msg)
             if (.not.direct_only) fri%direct_only = .false.
             print *, newton_x,'(direct only)'
          else
             if (n/=1) then
                print *, newton_x,'accepted'
                print *
             end if
             exit
          end if

       end do

       Vx = newton_x(1)
       sx = newton_x(2)
       
       deallocate(newton_x)
       deallocate(newton_xscale)
       deallocate(newton_fscale)
       deallocate(param)
       
    end if

  end subroutine couple_elastic_strength_rate_im


end module friction_im
