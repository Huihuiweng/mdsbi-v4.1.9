! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module friction_bm

  ! FRICTION_BM contains routines to couple friction laws and elasticity for bimaterials
  ! 
  ! Modified: 6 August 2010

  use constants, only : pr,pin

  implicit none


contains


  subroutine couple_elastic_strength_bm(fri,i,j,S0,Q,sz,vxp,vxm,vyp,vym,scale_Vmin,scale_Vmax, &
       sx,sy,scale_s,sx0,sy0,fxp,fxm,fyp,fym,cvsp,cvsm,transverse_slip,err)
    ! COUPLE_ELASTIC_STRENGTH_BM solves elasticity equation with stress<=strength
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
    ! SZ = stress in z direction
    ! Q = state variable
    ! VXP = particle velocity in x direction, plus side
    ! VXM = particle velocity in x direction, minus side
    ! VYP = particle velocity in y direction, plus side
    ! VYM = particle velocity in y direction, minus side
    ! SCALE_VMIN = scale of velocity, minimum
    ! SCALE_VMAX = scale of velocity, maximum
    ! SX = stress in x direction
    ! SY = stress in y direction
    ! SCALE_S = scale of stress
    ! SX0 = load in x direction
    ! SY0 = load in y direction
    ! FXP = stress transfer in x direction, plus side
    ! FXM = stress transfer in x direction, minus side
    ! FYP = stress transfer in y direction, plus side
    ! FYM = stress transfer in y direction, minus side
    ! CVSP = coefficient of shear radiation damping, plus side
    ! CVSM = coefficient of shear radiation damping, minus side
    ! TRANSVERSE_SLIP = scalar or vector slip
    ! ERR = flag that indicates error in solving friction law

    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: i,j
    real(pr),intent(in) :: S0,Q,sz,sx0,sy0,fxp,fxm,fyp,fym,cvsp,cvsm,scale_s,scale_Vmin,scale_Vmax
    real(pr),intent(inout) :: vxp,vxm,vyp,vym,sx,sy
    logical,intent(in) :: transverse_slip
    logical,intent(out) :: err

    ! Internal Parameters:
    ! VXPL = particle velocity in x direction if fault is locked, plus side
    ! VXML = particle velocity in x direction if fault is locked, minus side
    ! VYPL = particle velocity in y direction if fault is locked, plus side
    ! VYML = particle velocity in y direction if fault is locked, minus side
    ! SXL = stress in x direction if fault is locked
    ! SYL = stress in y direction if fault is locked
    ! SL = magnitude of shear stress vector if fault is locked
    ! NEWTON_X = independent variables for Newton's method
    ! NEWTON_XSCALE = scale of independent variables for convergence of Newton's method
    ! NEWTON_FSCALE = scale of function for convergence of Newton's method
    ! SIZE_PARAM = size of param
    ! PARAM = parameters required for function that is solved using Newton's method
    ! NEWTON_NMAX = maximum number of iterations allowed to solve equations with Newton's method
    ! NEWTON_TOLX = maximum tolerated error in independent variables
    ! NEWTON_TOLF = maximum tolerated error in dependent variables

    real(pr) :: vxpl,vxml,vypl,vyml,sxl,syl,sl
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
    fri%cV = cvsp*cvsm/(cvsp+cvsm)
    fri%sfx = (sx0+fxp)/(one+cvsp/cvsm)+(sx0+fxm)/(one+cvsm/cvsp)
    fri%sfy = (sy0+fyp)/(one+cvsp/cvsm)+(sy0+fym)/(one+cvsm/cvsp)

    if (transverse_slip) then

       ! assume fault is locked, calculate velocities and stresses
       
       vxpl = (fxp-fxm)/(cvsp+cvsm)
       vxml = vxpl
       sxl = sx0+fxp-cvsp*vxpl
       vypl = (fyp-fym)/(cvsp+cvsm)
       vyml = vypl
       syl = sy0+fyp-cvsp*vypl
       sl = sqrt(sxl**2+syl**2)
       
       ! compare strength to elastic shear stress if fault were locked,
       ! allow failure if stress exceeds strength
       
       if (fri%friction/='ratestate'.and.sl<S0) then
          
          ! stress less than strength, fault locked
          
          vxp = vxpl
          vxm = vxml
          vyp = vypl
          vym = vyml
          sx = sxl
          sy = syl
          
       else
          
          ! stress set to strength, fault slips;
          ! solution via Newton's method--see the subroutine elastic for further description
          
          allocate(newton_x(4))

          newton_x(1) = max(vxp-vxm,epsilon(one))
          newton_x(2) = sx
          newton_x(3) = vyp-vym
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
             newton_fscale(2) = max(abs(fxp),abs(fxm),abs(cvsp*vxp),abs(cvsm*vxm))
             newton_fscale(3) = max(abs(fyp),abs(fym),abs(cvsp*vyp),abs(cvsm*vym))
             newton_fscale(1) = newton_fscale(2)**2+newton_fscale(3)**2
             newton_fscale(4) = abs(newton_fscale(2)*(vyp-vym))+abs(newton_fscale(3)*(vxp-vxm))
          else
             newton_fscale = scale_s
          end if

          size_param = size(transfer(fri,param))
          allocate(param(size_param))
          param = transfer(fri,param)
       
          call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
               elastic_strength,newton_xscale,newton_fscale,err)
          
          vxp =  (-newton_x(2)+sx0+fxp)/cvsp
          vxm = -(-newton_x(2)+sx0+fxm)/cvsm
          sx = newton_x(2)
          vyp =  (-newton_x(4)+sy0+fyp)/cvsp
          vym = -(-newton_x(4)+sy0+fym)/cvsm
          sy = newton_x(4)
          
          deallocate(newton_x)
          deallocate(newton_xscale)
          deallocate(newton_fscale)
          deallocate(param)

       end if

    else

       ! lock fault in y-direction
       
       vyp = (fyp-fym)/(cvsp+cvsm)
       vym = vyp
       sy = sy0+fyp-cvsp*vyp
       
       ! assume fault is locked, calculate velocities and stresses
       
       vxpl = (fxp-fxm)/(cvsp+cvsm)
       vxml = vxpl
       sxl = sx0+fxp-cvsp*vxpl
       
       sl = abs(sxl)

       ! compare strength to elastic shear stress if fault were locked,
       ! allow failure if stress exceeds strength
       
       if (fri%friction/='ratestate'.and.sl<S0) then
          
          ! stress less than strength, fault locked
          
          vxp = vxpl
          vxm = vxml
          sx = sxl
          
       else
          
          ! stress set to strength, fault slips;
          ! solution via Newton's method--see the subroutine elastic for further description
          
          allocate(newton_x(2))

          newton_x(1) = vxp-vxm
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
             newton_fscale(1) = max(abs(fxp),abs(fxm),abs(cvsp*vxp),abs(cvsm*vxm))
             newton_fscale(2) = newton_fscale(1)
          else
             newton_fscale = scale_s
          end if

          size_param = size(transfer(fri,param))
          allocate(param(size_param))
          param = transfer(fri,param)
       
          call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
               elastic_strength,newton_xscale,newton_fscale,err)
          
          vxp =  (-newton_x(2)+sx0+fxp)/cvsp
          vxm = -(-newton_x(2)+sx0+fxm)/cvsm
          sx = newton_x(2)
          
          deallocate(newton_x)
          deallocate(newton_xscale)
          deallocate(newton_fscale)
          deallocate(param)

       end if

    end if

  end subroutine couple_elastic_strength_bm


  subroutine couple_elastic_strength_rate_bm(fri,i,j,sz,dt,a,Vo,So,dVo,dSo,vxp,vxm,vyp,vym,scale_Vmin,scale_Vmax, &
       sx,sy,scale_s,sx0,sy0,fxp,fxm,fyp,fym,cvsp,cvsm,transverse_slip,direct_only,err)
    ! COUPLE_ELASTIC_STRENGTH_RATE_BM solves elasticity equation with implicitly discretized friction law
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
    ! VXP = particle velocity in x direction, plus side
    ! VXM = particle velocity in x direction, minus side
    ! VYP = particle velocity in y direction, plus side
    ! VYM = particle velocity in y direction, minus side
    ! SCALE_VMIN = scale of velocity, minimum
    ! SCALE_VMAX = scale of velocity, maximum
    ! SX = stress in x direction
    ! SY = stress in y direction
    ! SCALE_S = scale of stress
    ! SX0 = load in x direction
    ! SY0 = load in y direction
    ! FXP = stress transfer in x direction, plus side
    ! FXM = stress transfer in x direction, minus side
    ! FYP = stress transfer in y direction, plus side
    ! FYM = stress transfer in y direction, minus side
    ! CVSP = coefficient of shear radiation damping, plus side
    ! CVSM = coefficient of shear radiation damping, minus side
    ! TRANSVERSE_SLIP = scalar or vector slip
    ! DIRECT_ONLY = solve for response by direct effect only
    ! ERR = flag that indicates error in solving friction law

    type(friction_type),intent(inout) :: fri
    integer(pin),intent(in) :: i,j
    real(pr),intent(in) :: sz,dt,a,Vo,So,dVo,dSo,sx0,sy0,fxp,fxm,fyp,fym,cvsp,cvsm,scale_s,scale_Vmin,scale_Vmax
    real(pr),intent(inout) :: vxp,vxm,vyp,vym,sx,sy
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
    fri%cV = cvsp*cvsm/(cvsp+cvsm)
    fri%sfx = (sx0+fxp)/(one+cvsp/cvsm)+(sx0+fxm)/(one+cvsm/cvsp)
    fri%sfy = (sy0+fyp)/(one+cvsp/cvsm)+(sy0+fym)/(one+cvsm/cvsp)
    fri%Vo = Vo
    fri%So = So
    fri%dVo = dVo
    fri%dSo = dSo
    fri%dt = dt
    fri%a = a
    fri%direct_only = direct_only

    if (transverse_slip) then

       ! stress set to strength, fault slips;
       ! solution via Newton's method--see the subroutine elastic for further description
       
       allocate(newton_x(4))
       
       newton_x(1) = max(vxp-vxm,epsilon(one))
       newton_x(2) = sx
       newton_x(3) = vyp-vym
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
          newton_fscale(2) = max(abs(fxp),abs(fxm),abs(cvsp*vxp),abs(cvsm*vxm))
          newton_fscale(3) = max(abs(fyp),abs(fym),abs(cvsp*vyp),abs(cvsm*vym))
          newton_fscale(1) = newton_fscale(2)**2+newton_fscale(3)**2
          newton_fscale(4) = abs(newton_fscale(2)*(vyp-vym))+abs(newton_fscale(3)*(vxp-vxm))
       else
          newton_fscale = scale_s
       end if
       
       size_param = size(transfer(fri,param))
       allocate(param(size_param))
       param = transfer(fri,param)
       
       call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param, &
            elastic_strength_rate,newton_xscale,newton_fscale,err)
       
       ! set values to ensure that elasticity is solved exactly
       ! (either accept Vx from Newton's method and set sx, or vice versa)

       !sx = newton_x(2)
       !sy = newton_x(4)
       !vxp =  (-sx+sx0+fxp)/cvsp
       !vxm = -(-sx+sx0+fxm)/cvsm
       !vyp =  (-sy+sy0+fyp)/cvsp
       !vym = -(-sy+sy0+fym)/cvsm

       vxp = (fxp-fxm+cvsm*newton_x(1))/(cvsp+cvsm)
       vxm = (fxp-fxm-cvsp*newton_x(1))/(cvsp+cvsm)
       vyp = (fyp-fym+cvsm*newton_x(3))/(cvsp+cvsm)
       vym = (fyp-fym-cvsp*newton_x(3))/(cvsp+cvsm)
       sx = sx0+fxp-cvsp*vxp
       sy = sy0+fyp-cvsp*vyp
       
       deallocate(newton_x)
       deallocate(newton_xscale)
       deallocate(newton_fscale)
       deallocate(param)

    else
       
       ! lock fault in y-direction
       
       vyp = (fyp-fym)/(cvsp+cvsm)
       vym = vyp
       sy = sy0+fyp-cvsp*vyp
       
       ! stress set to strength, fault slips;
       ! solution via Newton's method--see the subroutine elastic for further description
       
       allocate(newton_x(2))
              
       newton_x(1) = vxp-vxm
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
          newton_fscale(1) = max(abs(fxp),abs(fxm),abs(cvsp*vxp),abs(cvsm*vxm))
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
             print *, vxp-vxm,sx
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

       ! set values to ensure that elasticity is solved exactly
       ! (either accept Vx from Newton's method and set sx, or vice versa)

       !sx = newton_x(2)
       !vxp =  (-sx+sx0+fxp)/cvsp
       !vxm = -(-sx+sx0+fxm)/cvsm

       vxp = (fxp-fxm+cvsm*newton_x(1))/(cvsp+cvsm)
       vxm = (fxp-fxm-cvsp*newton_x(1))/(cvsp+cvsm)
       sx = sx0+fxp-cvsp*vxp
       
       deallocate(newton_x)
       deallocate(newton_xscale)
       deallocate(newton_fscale)
       deallocate(param)
       
    end if

  end subroutine couple_elastic_strength_rate_bm


end module friction_bm

