! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module friction_solver

  ! FRICTION_SOLVER contains routines to couple friction laws and elasticity
  ! 
  ! Modified: 11 October 2007

  use constants, only : pr,pin

  implicit none


contains


  subroutine couple_elastic_strength(mdl,fri,i,j,S0,Q,Vx,Vy,Vz,scale_Vmin,scale_Vmax, &
       Tx,Ty,Tz,scale_s,sfx,sfy,sfz,cVs,cVn,mx,my,err)
    ! COUPLE_ELASTIC_STRENGTH solves elasticity equation with stress<=strength
    ! 
    ! Modified: 11 October 2007

    use constants, only : zero,one,two,three
    use model, only : model_type
    use friction, only : friction_type,elastic_strength
    use utilities, only : newton

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! FRI = friction variables
    ! I = index in x direction
    ! J = index in y direction
    ! S0 = portion of shear strength independent of slip velocity
    ! Q = state variable
    ! VX = slip velocity in x direction
    ! VY = slip velocity in y direction
    ! VZ = opening velocity in z direction
    ! SCALE_VMIN = scale of velocity, minimum
    ! SCALE_VMAX = scale of velocity, maximum
    ! TX = shear traction in x direction
    ! TY = shear traction in y direction
    ! N = normal stress, positive in compression
    ! SCALE_S = scale of stress
    ! SFX = load+stress transfer in x direction
    ! SFY = load+stress transfer in y direction
    ! SFZ = load+stress transfer in z direction
    ! CVS = coefficient of shear radiation damping
    ! CVN = coefficient of normal radiation damping
    ! TRANSVERSE_SLIP = scalar or vector slip
    ! OPENING = permit opening or not
    ! ERR = flag that indicates error in solving friction law

    type(model_type),intent(in) :: mdl
    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: i,j
    real(pr),intent(in) :: Q,S0,sfx,sfy,sfz,cVs,cVn,scale_s,scale_Vmin,scale_Vmax
    real(pr),intent(inout) :: Vx,Vy,Vz,sx,sy,sz
    logical,intent(in) :: transverse_slip,opening
    logical,intent(out) :: err
    
    ! Internal Parameters:
    ! TL = shear stress if fault is locked
    ! NEWTON_X = independent variables for Newton's method
    ! NEWTON_XSCALE = scale of independent variables for convergence of Newton's method
    ! NEWTON_FSCALE = scale of function for convergence of Newton's method
    ! PARAM = parameters required for function that is solved using Newton's method
    ! NEWTON_NMAX = maximum number of iterations allowed to solve equations with Newton's method
    ! NEWTON_TOLX = maximum tolerated error in independent variables
    ! NEWTON_TOLF = maximum tolerated error in dependent variables

    real(pr) :: Tl
    real(pr),dimension(:),allocatable :: newton_x,newton_xscale,newton_fscale,param
    integer(pin) :: newton_nmax
    real(pr) :: newton_tolx,newton_tolf
   
    ! set iteration limit and tolerances

    newton_nmax = 1000
    newton_tolx = epsilon(one)
    newton_tolf = epsilon(one)**(two/three)
 
    ! store variables needed elsewhere

    fri%i = i
    fri%j = j
    fri%S0 = S0
    fri%Q = Q
    fri%cVs = cVs
    fri%cVn = cVn
    fri%sfx = sfx
    fri%sfy = sfy
    fri%sfz = sfz
    fri%mx = mx
    fri%my = my

    if (mdl%transverse_slip) then

       ! assume fault is locked, calculate shear stress
       
       Tl = sqrt(sfx**2+sfy**2)
    
       ! compare strength to elastic shear stress if fault were locked,
       ! allow failure if stress exceeds strength
          
       if (fri%friction/='ratestate'.and.Tl<S0) then
          
          ! stress less than strength, fault locked
          
          Vx = zero
          Tx = sfx
          Vy = zero
          Ty = sfy
          
       else
          
          ! stress set to strength, fault slips;
          ! solution via Newton's method--see the subroutine elastic_strength for further description
          
          allocate(newton_x(4))

          newton_x(1) = Vx
          newton_x(2) = Tx
          newton_x(3) = Vy
          newton_x(4) = Ty
          
          allocate(newton_xscale(4))

          newton_xscale(1) = min(max(abs(newton_x(1)),scale_Vmin),scale_Vmax)
          newton_xscale(3) = newton_xscale(1)

          if (scale_s==zero) then
             newton_xscale(2) = max(abs(Tx),abs(Ty))
             newton_xscale(4) = max(abs(Tx),abs(Ty))
          else
             newton_xscale(2) = scale_s
             newton_xscale(4) = scale_s
          end if

          allocate(newton_fscale(4))

          if (scale_s==zero) then
             newton_fscale(2) = max(abs(fx),abs(cV*Vx))
             newton_fscale(3) = max(abs(fy),abs(cV*Vy))
             newton_fscale(1) = newton_fscale(2)**2+newton_fscale(3)**2
             newton_fscale(4) = abs(newton_fscale(2)*Vy)+abs(newton_fscale(3)*Vx)
          else
             newton_fscale = scale_s
          end if

          allocate(param(size(transfer(fri,param))))
          param = transfer(fri,param)

          call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
               elastic_strength,newton_xscale,newton_fscale,err)
          
          Vx = newton_x(1)
          Tx = newton_x(2)
          Vy = newton_x(3)
          Ty = newton_x(4)
          
          deallocate(newton_x)
          deallocate(newton_xscale)
          deallocate(newton_fscale)
          deallocate(param)

       end if

    else

       ! lock fault in y-direction
       
       Vy = zero
       Ty = sfy
       
       ! first attempt solution without opening

       Vz = zero

       ! assume fault is locked, calculate stress
       
       Tl = abs(sfx)
       
       ! compare strength to elastic shear stress if fault were locked,
       ! allow failure if stress exceeds strength
       
       if (fri%friction/='ratestate'.and.Tl<S0) then
          
          ! stress less than strength, fault locked
          
          Vx = zero
          Tx = sfx
          
       else
          
          ! stress set to strength, fault slips;
          ! solution via Newton's method--see the subroutine elastic_strength for further description
          
          allocate(newton_x(2))

          newton_x(1) = Vx
          newton_x(2) = Tx
          
          allocate(newton_xscale(2))

          newton_xscale(1) = min(max(abs(newton_x(1)),scale_Vmin),scale_Vmax)

          if (scale_s==zero) then
             newton_xscale(2) = abs(Tx)
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

          allocate(param(size(transfer(fri,param))))
          param = transfer(fri,param)

          call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
               elastic_strength,newton_xscale,newton_fscale,err)
          
          Vx = newton_x(1)
          Tx = newton_x(2)
          
          deallocate(newton_x)
          deallocate(newton_xscale)
          deallocate(newton_fscale)
          deallocate(param)

       end if
       
    end if

  end subroutine couple_elastic_strength


  subroutine couple_elastic_strength_rate_im(fri,i,j,sz,dt,a,Vo,So,dVo,dSo,Vx,Vy,scale_Vmin,scale_Vmax, &
       sx,sy,scale_s,sx0,sy0,fx,fy,cV,transverse_slip,direct_only,err)
    ! COUPLE_ELASTIC_STRENGTH_RATE_IM solves elasticity equation with implicitly discretized friction law
    ! 
    ! Modified: 18 September 2007

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
    ! NEWTON_NMAX = maximum number of iterations allowed to solve equations with Newton's method
    ! NEWTON_TOLX = maximum tolerated error in independent variables
    ! NEWTON_TOLF = maximum tolerated error in dependent variables
    ! MSG = error message
    ! N = iteration on initial guess

    real(pr),dimension(:),allocatable :: newton_x,newton_xscale,newton_fscale,param
    integer(pin) :: newton_nmax,n
    real(pr) :: newton_tolx,newton_tolf
    character(64) :: msg

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
       
       newton_x(1) = Vx
       newton_x(2) = sx
       newton_x(3) = Vy
       newton_x(4) = sy
       
       allocate(newton_xscale(4))
       
       newton_xscale(1) = min(max(abs(newton_x(1)),scale_Vmin),scale_Vmax)
       newton_xscale(3) = newton_xscale(1)
       
       if (scale_s==zero) then
          newton_xscale(2) = max(abs(sx),abs(sy))
          newton_xscale(4) = max(abs(sx),abs(sy))
       else
          newton_xscale(2) = scale_s
          newton_xscale(4) = scale_s
       end if
       
       allocate(newton_fscale(4))
       
       if (scale_s==zero) then
          newton_fscale(2) = max(abs(fx),abs(cV*Vx))
          newton_fscale(3) = max(abs(fy),abs(cV*Vy))
          newton_fscale(1) = newton_fscale(2)**2+newton_fscale(3)**2
          newton_fscale(4) = abs(newton_fscale(2)*Vy)+abs(newton_fscale(3)*Vx)
       else
          newton_fscale = scale_s
       end if
       
       allocate(param(size(transfer(fri,param))))
       param = transfer(fri,param)
       
       call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
            elastic_strength_rate,newton_xscale,newton_fscale,err)
       
       Vx = newton_x(1)
       sx = newton_x(2)
       Vy = newton_x(3)
       sy = newton_x(4)
       
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
       
       allocate(param(size(transfer(fri,param))))
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


end module friction_solver
