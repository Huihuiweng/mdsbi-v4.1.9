! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module static

  ! STATIC contains routines to solve a static elasticity problem
  ! 
  ! Modified: 3 December 2007

  use constants, only : pr,pin

  implicit none


contains


  subroutine iterate_static(pb)
    ! ITERATE_STATIC solves static elasticity problems, iteratively varying initial
    ! stress levels to remove singularity at crack tips.
    ! 
    ! Modified: 26 January 2007

    use constants, only : zero,half,one
    use problem, only : problem_type
    use model, only : get_bound_range

    implicit none

    ! I/O Parameters:
    ! PB = problem variables

    type(problem_type),intent(inout) :: pb

    ! Internal Parameters:
    ! M = minimum index of slipping zone in x direction
    ! P = maximum index of slipping zone in x direction
    ! SLIP_XMIN = minimum location of slipping zone in x direction
    ! SLIP_XMAX = maximum location of slipping zone in x direction
    ! NM = index of left-tip iteration
    ! NP = index of right-tip iteration
    ! EPS = offset used to calculate finite-difference Jacobian
    ! XMN = left-tip location, new iteration
    ! XMO = left-tip location, old iteration
    ! XPN = right-tip location, new iteration
    ! XPO = right-tip location, old iteration
    ! E = error at right crack tip
    ! E = error at left crack tip
    ! XM = left side of fully weakened region
    ! XP = right side of fully weakened region

    integer(pin) :: m,p,nm,np
    real(pr) :: slip_xmin,slip_xmax,eps,xmn,xmo,xpn,xpo,emo,emn,epo,epn,dedxm,dedxp

    ! set eps to square-root of machine precision

    eps = sqrt(epsilon(one))

    ! calculate size of x (assuming slip zone is restricted in x-direction only)

    slip_xmin = pb%fri%slip_x0-pb%fri%slip_x
    slip_xmax = pb%fri%slip_x0+pb%fri%slip_x

    if (pb%mdl%bm) then
       call get_bound_range(field='uxp',dir='x',val_min=slip_xmin,val_max=slip_xmax, &
            ots_min=.false.,ots_max=.false.,mdl=pb%mdl,imin=m,imax=p)
    else
       call get_bound_range(field='Ux',dir='x',val_min=slip_xmin,val_max=slip_xmax, &
            ots_min=.false.,ots_max=.false.,mdl=pb%mdl,imin=m,imax=p)
    end if

    ! set initial slip
    
    if (pb%mdl%nz/=0) then
       if (pb%mdl%bm) then
          pb%fld%flt%uxp(m:p,:,1) =  half*half*pb%fri%sw%Dc(m:p,:)
          pb%fld%flt%uxm(m:p,:,1) = -half*half*pb%fri%sw%Dc(m:p,:)
       else
          pb%fld%flt%Ux(m:p,:,1) = half*pb%fri%sw%Dc(m:p,:)
       end if
    end if

    ! set initial trial values

    xmn = pb%fri%sw%xm
    xmo = pb%fri%sw%xm
    xpn = pb%fri%sw%xp
    xpo = pb%fri%sw%xp

    ! determine location of tips

    do np = 0,pb%mdl%niter

       ! set initial position of right tip
       
       pb%fri%sw%xp = xpn
       pb%fri%sw%dfdxp = 0.1_pr/(pb%fri%slip_x-xpn)

       do nm = 0,pb%mdl%niter

          ! set initial position of left tip
          
          pb%fri%sw%xm = xmn
          pb%fri%sw%dfdxm = 0.1_pr/(pb%fri%slip_x+xmn)

          ! solve static problem
          
          call solve_static(pb)
          
          ! calculate error in stress at left tip

          ! stress at node beyond tips equals peak strength
          epn = pb%fld%flt%sx(p+1,1,1)+pb%fri%sw%fs(p+1,1)*pb%fld%flt%sz(p+1,1,1)
          emn = pb%fld%flt%sx(m-1,1,1)+pb%fri%sw%fs(m-1,1)*pb%fld%flt%sz(m-1,1,1)

          ! continuity of stress
          !epn = pb%fld%flt%sx(p+1,1,1)-pb%fld%flt%sx(p,1,1)
          !emn = pb%fld%flt%sx(m-1,1,1)-pb%fld%flt%sx(m,1,1)

          if (nm==0) emo = emn
          if (np==0) epo = epn

          print *
          print *, 'np=',np,' nm=',nm
          print *, 'xm=',xmn,' xp=',xpn
          print *, 'dfdxm=',pb%fri%sw%dfdxm,' dfdxp=',pb%fri%sw%dfdxp
          print *, 'em=',emn,' ep=',epn
          
          if (pb%mdl%niter==0) exit ! exit without iterating

          ! exit when error is sufficiently small

          if (abs(emn)<pb%mdl%rtol) exit

          ! secant approximation to derivative
       
          if (emn/=emo) then
             dedxm = (emn-emo)/(xmn-xmo) 
          else
             dedxm = -emn/(eps*xmn)
          end if
       
          ! update
          
          xmo = xmn
          emo = emn
          xmn = xmo-emo/dedxm

       end do

       if (pb%mdl%niter==0) exit ! exit without iterating
       
       ! exit when error is sufficiently small
       
       if (abs(epn)<pb%mdl%rtol) exit
       
       ! secant approximation to derivative
       
       if (epn/=epo) then
          dedxp = (epn-epo)/(xpn-xpo) 
       else
          dedxp = -epn/(eps*xpn)
       end if
       
       ! update
       
       xpo = xpn
       epo = epn
       xpn = xpo-epo/dedxp

    end do

  end subroutine iterate_static


  subroutine solve_static(pb)
    ! SOLVE_STATIC solves static elasticity problem, using Newton's method, 
    ! combined with a forward-difference approximation to calculate the 
    ! Jacobian (the routines are similar to "newton" and "jacfd" in utilities.f90,
    ! but are coded explicitly here)
    ! 
    ! Modified: 4 May 2007

    use constants, only : zero,one
    use problem, only : problem_type
    use model, only : get_bound_range
    use io, only : warning,message

    implicit none

    ! I/O Parameters:
    ! PB = problem variables

    type(problem_type),intent(inout) :: pb

    ! uncomment to use LAPACK's linear solver
    interface
       subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
         ! DGESV solves ax=b for x, which overwrites b on output
         ! (routine from LAPACK library)
         use constants, only : pr,pin
         implicit none
         integer(pin),intent(in) :: n,nrhs,lda,ldb
         integer(pin),intent(out) :: info
         integer(pin),dimension(n),intent(out) :: ipiv
         real(pr),dimension(ldb,nrhs),intent(inout) :: b
         real(pr),dimension(lda,n),intent(inout) :: a
       end subroutine dgesv
    end interface

    ! Internal Parameters:
    ! M = minimum index of slipping zone in x direction
    ! P = maximum index of slipping zone in x direction
    ! SLIP_XMIN = minimum location of slipping zone in x direction
    ! SLIP_XMAX = maximum location of slipping zone in x direction
    ! N = index of Newton iteration
    ! I = index of x
    ! NX = length of x
    ! EPS = offset used to calculate finite-difference Jacobian
    ! G0 = vector containing misfit function (at current iteration)
    ! G = vector containing misfit function (offset from current iteration values)
    ! X = vector of slip or displacement (at current iteration)
    ! C = correction to x from Newton's method
    ! TEMP = temporary variable holding component of slip/displacement vector
    ! DX = offset of x
    ! JAC = Jacobian
    ! IPIV = index of pivot (used to diagnose errors)
    ! INFO = error flag (non-zero if error)
    ! STR = string used in output

    integer(pin) :: m,p,n,i,nx,info
    integer(pin),dimension(:),allocatable :: ipiv
    character(64) :: str
    real(pr) :: slip_xmin,slip_xmax,eps,temp,dx
    real(pr),dimension(:),allocatable :: x,g0,g,c
    real(pr),dimension(:,:),allocatable :: jac

    ! set eps to square-root of machine precision

    eps = sqrt(epsilon(one))

    ! calculate size of x (assuming slip zone is restricted in x-direction only)

    slip_xmin = pb%fri%slip_x0-pb%fri%slip_x
    slip_xmax = pb%fri%slip_x0+pb%fri%slip_x

    if (pb%mdl%bm) then
       call get_bound_range(field='uxp',dir='x',val_min=slip_xmin,val_max=slip_xmax, &
            ots_min=.false.,ots_max=.false.,mdl=pb%mdl,imin=m,imax=p)
    else
       call get_bound_range(field='Ux',dir='x',val_min=slip_xmin,val_max=slip_xmax, &
            ots_min=.false.,ots_max=.false.,mdl=pb%mdl,imin=m,imax=p)
    end if

    if (pb%mdl%bm) then
       if (pb%mdl%transverse_slip) then
          nx = 4*(p-m+1)*pb%mdl%ny
       else
          nx = 2*(p-m+1)*pb%mdl%ny
       end if
    else
       if (pb%mdl%transverse_slip) then
          nx = 2*(p-m+1)*pb%mdl%ny
       else
          nx = (p-m+1)*pb%mdl%ny
       end if
    end if

    ! allocate memory to arrays

    allocate(x(nx),g(nx),g0(nx),c(nx),jac(nx,nx),ipiv(nx))

    ! convert slip/displacement into x vector

    call pack_U(x,pb%mdl,pb%fld%flt,m,p)

    ! solve static problem with Newton's method

    do n = 1,pb%mdl%nt
       
       ! write iteration step to status file
       
       write(str,'("iteration step: ",i6," out of ",i6)') n,pb%mdl%nt
       call message(str)
       
       ! evaluate misfit at initial guess
       
       call evaluate_static(x,g0,pb,m,p)
       
       ! calculate error in g
       
       pb%fld%flt%aerr = maxval(abs(g0))
       
       ! write error in g to status file
       
       write(str,'(".   error in g: ",e12.7," tolerance: ",e12.7)') pb%fld%flt%aerr(1),pb%mdl%atol
       call message(str)
       
       ! exit if error in g is sufficiently small
       
       if (pb%fld%flt%aerr(1)<=pb%mdl%atol) exit
       
       ! calculate Jacobian with forward finite-differences
       
       do i = 1,nx
          temp = x(i)
          dx = eps*abs(temp)
          if (dx==zero) dx = eps
          x(i) = temp+dx
          dx = x(i)-temp
          call evaluate_static(x,g,pb,m,p)
          x(i) = temp
          jac(:,i) = (g-g0)/dx
       end do
       
       ! solve linear system for correction
       
       ! uncomment to use LAPACK's linear solver
       call dgesv(nx,1,jac,nx,ipiv,g0,nx,info)
       if (info/=0) then
          call warning('Warning: Problem solving linear system','solve_static')
          write(str,'("info= ",i10)') info
          call message(str)
          exit
       end if
       c = -g0
        
       ! uncomment to use Gaussian elimination
       !call gauss(jac,-g0,c)
  
       ! calculate error in x
       
       pb%fld%flt%aerr = maxval(abs(c))
       
       ! apply correction
       
       x = x+c
       
       ! write error in x to status file
       
       write(str,'(".   error in U: ",e12.7," tolerance: ",e12.7)') pb%fld%flt%aerr(1),pb%mdl%atol
       call message(str)
       
       ! exit if error in x is sufficiently small
       
       if (pb%fld%flt%aerr(1)<=pb%mdl%atol) exit
       
    end do
    
    ! print warning message if convergence not achieved in Newton's method
    
    if (n>=pb%mdl%nt) then
       str = "Warning: Newton's method did not converge"
       call warning(trim(adjustl(str)),'solve_static')
       write(str,'(a,i5)') 'Iterations: nmax=',pb%mdl%nt
       call message(str)
    end if

    ! convert final x into slip/displacement
    
    call unpack_U(x,pb%mdl,pb%fld%flt,m,p)
    
    ! deallocate memory to arrays
    
    deallocate(x,g,g0,c,jac,ipiv)
    
  end subroutine solve_static


  subroutine evaluate_static(x,g,pb,m,p)
    ! EVALUATE_STATIC computes the misfit function (how well a particular slip or displacement
    ! solution satisfies frictional boundary conditions and static elasticity) for an input 
    ! slip/displacement solution
    !
    ! Modified: 3 December 2007

    use problem, only : problem_type
    use rates, only : set_rates
    use friction_routines, only : solve_friction_normal,static_friction

    ! I/O Parameters:
    ! X = vector of slip or displacement
    ! G = vector containing misfit function
    ! PB = problem variables
    ! M = minimum index of slipping zone in x direction
    ! P = maximum index of slipping zone in x direction

    real(pr),dimension(:),intent(in) :: x
    real(pr),dimension(size(x)),intent(out) :: g
    type(problem_type),intent(inout) :: pb
    integer(pin),intent(in) :: m,p

    ! Internal Parameters:
    ! ACCEPT_STEP = accept integration step
    ! STAGE = integration stage (stage=1 here)
    ! UPDATE_HISTORY = update or overwrite history of Fourier coefficients
    ! STRESS_TRANSFER = method used to calculate stress transfer
    !      convolution = perform convolution
    !      interpolation = interpolate between known values

    logical :: accept_step
    integer(pin),parameter :: stage = 1
    logical,parameter :: update_history = .false.
    character(64),parameter :: stress_transfer = 'convolution'

    ! convert input x into slip/displacement

    call unpack_U(x,pb%mdl,pb%fld%flt,m,p)

    ! set load, calculate fault-normal displacements, solve for stress transfer
    
    call set_rates(pb,stage,stress_transfer,accept_step,update_history)

    ! set normal stress, assuming no opening
    
    call solve_friction_normal(pb%mdl,pb%fld%flt,pb%cnv,stage)

    ! set shear stress, evaluate frictional boundary conditions to calculate misfit

    call static_friction(pb%fld%flt,pb%mdl,pb%fri,stage,g,m,p)

  end subroutine evaluate_static


  subroutine pack_U(x,mdl,flt,m,p)
    ! PACK_U packs slip/displacement arrays into vector x
    ! 
    ! Modified: 27 November 2006

    use model, only : model_type
    use fault, only : fault_type

    implicit none

    ! I/O Parameters:
    ! X = vector of slip/displacement
    ! MDL = model variables
    ! FLT = fault variables
    ! M = minimum index of slipping zone in x direction
    ! P = maximum index of slipping zone in x direction

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(in) :: flt
    real(pr),dimension(:),intent(out) :: x
    integer(pin),intent(in) :: m,p

    ! Internal Parameters:
    ! N = product of x and y dimensions
    ! S = stage (=1 here)

    integer(pin) :: n
    integer(pin),parameter :: s=1

    n = (p-m+1)*mdl%ny

    if (mdl%bm) then

       if (mdl%transverse_slip) then

          x(    1:  n) = reshape(flt%uxp(m:p,:,s),shape(x(1:n)))
          x(  n+1:2*n) = reshape(flt%uxm(m:p,:,s),shape(x(1:n)))
          x(2*n+1:3*n) = reshape(flt%uyp(m:p,:,s),shape(x(1:n)))
          x(3*n+1:4*n) = reshape(flt%uym(m:p,:,s),shape(x(1:n)))

       else

          x(    1:  n) = reshape(flt%uxp(m:p,:,s),shape(x(1:n)))
          x(  n+1:2*n) = reshape(flt%uxm(m:p,:,s),shape(x(1:n)))

       end if

    else
       
       if (mdl%transverse_slip) then
          
          x(  1:  n) = reshape(flt%Ux(m:p,:,s),shape(x(1:n)))
          x(n+1:2*n) = reshape(flt%Uy(m:p,:,s),shape(x(1:n)))
          
       else
          
          x(1:n) = reshape(flt%Ux(m:p,:,s),shape(x(1:n)))
          
       end if

    end if
    
  end subroutine pack_U


  subroutine unpack_U(x,mdl,flt,m,p)
    ! UNPACK_U unpacks vector x into slip/displacement arrays
    ! 
    ! Modified: 27 November 2006

    use model, only : model_type
    use fault, only : fault_type

    implicit none

    ! I/O Parameters:
    ! X = vector of slip/displacement
    ! MDL = model variables
    ! FLT = fault variables
    ! M = minimum index of slipping zone in x direction
    ! P = maximum index of slipping zone in x direction

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    real(pr),dimension(:),intent(in) :: x
    integer(pin),intent(in) :: m,p

    ! Internal Parameters:
    ! N = product of x and y dimensions
    ! S = stage (=1 here)

    integer(pin) :: n
    integer(pin),parameter :: s=1

    n = (p-m+1)*mdl%ny

    if (mdl%bm) then

       if (mdl%transverse_slip) then

          flt%uxp(m:p,:,s) = reshape(x(    1:  n),shape(flt%uxp(m:p,:,s)))
          flt%uxm(m:p,:,s) = reshape(x(  n+1:2*n),shape(flt%uxm(m:p,:,s)))
          flt%uyp(m:p,:,s) = reshape(x(2*n+1:3*n),shape(flt%uyp(m:p,:,s)))
          flt%uym(m:p,:,s) = reshape(x(3*n+1:4*n),shape(flt%uym(m:p,:,s)))

       else

          flt%uxp(m:p,:,s) = reshape(x(    1:  n),shape(flt%uxp(m:p,:,s)))
          flt%uxm(m:p,:,s) = reshape(x(  n+1:2*n),shape(flt%uxm(m:p,:,s)))

       end if

    else
    
       if (mdl%transverse_slip) then

          flt%Ux(m:p,:,s) = reshape(x(  1:  n),shape(flt%Ux(m:p,:,s)))
          flt%Uy(m:p,:,s) = reshape(x(n+1:2*n),shape(flt%Uy(m:p,:,s)))

       else

          flt%Ux(m:p,:,s) = reshape(x(  1:  n),shape(flt%Ux(m:p,:,s)))

       end if

    end if

  end subroutine unpack_U


end module static
