! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module utilities

  ! UTILITIES contains variables and routines for miscellaneous tasks
  ! 
  ! Modified: 6 August 2010

  use constants, only : pin,pr

  implicit none

  interface interp_linear
     ! INTERP_LINEAR overloads a linear interpolation routines
     module procedure interp_linear_1d,interp_linear_2d,interp_linear_3d,interp_linear_vector
  end interface

  interface interp_spline
     ! INTERP_SPLINE overloads a cubic spline interpolation routines
     module procedure interp_spline_scalar,interp_spline_vector
  end interface

contains


  function search_binary(x,y,mode,out_range,ig,z) result(mid)
    ! SEARCH_BINARY performs a binary search of ordered list x(1:n), 
    ! returning index i of x whose value is closest to y and satisfies either
    ! x(i)<=y<x(i+1) (floor mode) or x(i-1)<y<=x(i) (ceiling mode)
    ! as well as z = x(i).  The values of i must satisfy 1<=i<=n.  If
    ! the value y lies below (or equal to, for ceiling mode) x(1) then 
    ! i=1 is returned.  If the value y lies above (or equal to, for floor
    ! mode) x(n) then i=n is returned.  In these cases, out_range is
    ! set to T to indicate that y lies outside range of x.
    ! 
    ! Modified: 17 February 2006

    use io, only : error

    implicit none

    ! I/O Parameters:
    ! X = ordered list
    ! Y = value to be found in list
    ! MID = mean of upper and low bounds
    ! MODE = whether x(i) is greater or lesser than y
    !      'floor' = x(i)<=y<x(i+1)
    !      'ceiling' = x(i-1)<y<=x(i)
    ! OUT_RANGE = flag used to signify that y lies outside range of x
    !      T = outside range of x
    !      F = inside range of x
    ! IG = guess at index i
    ! Z = value of x(i) that satisfies inequality

    real(pr),dimension(:),intent(in) :: x
    real(pr),intent(in) :: y
    character(len=*),intent(in) :: mode
    logical,intent(out) :: out_range
    integer(pin),intent(in),optional :: ig
    real(pr),intent(out),optional :: z
    integer(pin) :: mid

    ! Internal Parameters:
    ! N = length of x
    ! LOW = lower bound on index
    ! HIGH = upper bound on index
    ! I = used in expanding binary search
    ! J = used in expanding binary search

    integer(pin) :: n,low,high,i,j

    n = size(x)

    ! special case of size(x)==1

    if (n==1) then
       mid = 1
       if (y==x(1)) then
          out_range = .false.
       else
          out_range = .true.
       end if
       if (present(z)) z = x(1)
       return
    end if

    ! use initial guess (if given) to isolate range containing index,
    ! otherwise start from full range

    if (present(ig)) then

       ! make sure initial guess is within bounds

       i = min(max(ig,1),n)
       
       ! determine if y is above or below guess x(i), set low and high

       if (y==x(i)) then
          mid = i
          if (present(z)) z = x(i)
          out_range = .false.
          return
       elseif (y<x(i)) then
          ! y less than x(i) => i is upper bound
          high = i
          ! find lower bound, moving out in powers of 2
          j = 0
          do
             low = max(1,high-2**j)
             if (x(low)<y.or.low==1) exit
             j = j+1
          end do
       else ! y>x(i)
          ! y greater than x(i) => i is lower bound
          low = i
          ! find upper bound, moving out in powers of 2
          j = 0
          do
             high = min(low+2**j,n)
             if (y<x(high).or.high==n) exit
             j = j+1
          end do
       end if

    else

       low = 1
       high = n

    end if

    ! perform binary search to locate index

    select case(mode)

    case('floor')

       if (y<x(low)) then
          ! out of bounds below
          mid = low
          out_range = .true.
       elseif (y>=x(high)) then
          ! out of bounds above
          mid = high-1
          out_range = .true.
       else
          mid = search_binary_floor(x,y,low,high)
          out_range = .false.
       end if

    case('ceiling')

       if (y<=x(low)) then
          ! out of bounds below
          mid = low+1
          out_range = .true.
       elseif (y>x(high)) then
          ! out of bounds above
          mid = high
          out_range = .true.
       else
          mid = search_binary_ceiling(x,y,low,high)
          out_range = .false.
       end if

    case default

       call error('Error: Invalid mode in binary search','search_binary')

    end select

    if (present(z)) z = x(mid)

  end function search_binary


  recursive function search_binary_floor(x,y,low,high) result(mid)
    ! SEARCH_BINARY_FLOOR performs a binary search of ordered list x, 
    ! returning index i of x whose value is closest to y and satisfies x(i)<=y<x(i+1)
    ! 
    ! Modified: 22 July 2005

    implicit none

    ! I/O Parameters:
    ! X = ordered list
    ! Y = value to be found in list
    ! LOW = lower bound on index
    ! HIGH = upper bound on index
    ! MID = mean of upper and low bounds

    real(pr),dimension(:),intent(in) :: x
    real(pr),intent(in) :: y
    integer(pin),intent(in) :: low,high
    integer(pin) :: mid

    ! find mean of low and high, rounding down
    mid = (low+high)/2

    if (x(mid)<=y) then
       ! value greater than or equal to mid
       if (y<x(mid+1)) then
          ! value less than mid+1, properly bracketed
          return
       else
          ! value greater than mid+1, make this new lower bound
          mid = search_binary_floor(x,y,mid+1,high)
       end if
    else
       ! value less than mid, make this new upper bound
       mid = search_binary_floor(x,y,low,mid)
    end if

  end function search_binary_floor


  recursive function search_binary_ceiling(x,y,low,high) result(mid)
    ! SEARCH_BINARY_CEILING performs a binary search of ordered list x, 
    ! returning index i of x whose value is closest to y and satisfies x(i-1)<y<=x(i)
    ! 
    ! Modified: 22 July 2005

    implicit none

    ! I/O Parameters:
    ! X = ordered list
    ! Y = value to be found in list
    ! LOW = lower bound on index
    ! HIGH = upper bound on index
    ! MID = mean of upper and low bounds

    real(pr),dimension(:),intent(in) :: x
    real(pr),intent(in) :: y
    integer(pin),intent(in) :: low,high
    integer(pin) :: mid

    ! find mean of low and high, rounding up
    mid = ceiling(real(low+high,pr)/2._pr)

    if (x(mid)>=y) then
       ! value less than mid
       if (y>x(mid-1)) then
          ! value less than or equal to mid-1, properly bracketed
          return
       else
          ! value less than mid-1, make this new upper bound
          mid = search_binary_ceiling(x,y,low,mid-1)
       end if
    else
       ! value greater than or equal to mid, make this new lower bound
       mid = search_binary_ceiling(x,y,mid,high)
    end if

  end function search_binary_ceiling


  function interp_linear_vector(x,f,xi) result(fi)
    ! INTERP_LINEAR_VECTOR performs a linear interpolation of a vector
    ! 
    ! Modified: 17 February 2006

    use io, only : error

    implicit none

    ! I/O Parameters:
    ! X = independent variable (data)
    ! F = dependent variable (data)
    ! XI = independent variable (points at which function is to be interpolated)
    ! FI = dependent variable (linearly interpolated values of f(xi))

    real(pr),dimension(:),intent(in) :: x,f
    real(pr),dimension(:),intent(in) :: xi
    real(pr),dimension(size(xi)) :: fi

    ! Internal Parameters:
    ! I = index of xi vector
    ! J = holds index of x(i) such that x(j)<=xi<x(j+1)
    ! OTS = flag denoting if xi lies outside range of x

    integer(pin) :: i
    integer(pin),dimension(size(xi)) :: j
    logical :: ots

    ! find index j such that x(j)<=xi(i)<x(j+1)

    do i = 1,size(xi)

       if (i==1) then
          ! use binary search to find starting value
          j(i) = search_binary(x,xi(i),'floor',ots)
       else
          ! use previous index as starting value for binary search
          j(i) = search_binary(x,xi(i),'floor',ots,j(i-1))
       end if

       if (ots) call error('Error: Attempting to extrapolating data outside range','linear_interp_vector')

    end do

    ! linear interpolation
    
    fi = f(j)+(f(j+1)-f(j))/(x(j+1)-x(j))*(xi-x(j))

  end function interp_linear_vector


  function interp_linear_1d(x,f,xi) result(fi)
    ! INTERP_LINEAR_1D performs linear interpolation of a single point
    ! 
    ! Modified: 20 February 2005

    implicit none

    ! I/O Parameters:
    ! X = independent variable (data)
    ! F = dependent variable (data)
    ! XI = independent variable (point at which function is to be interpolated)
    ! FI = dependent variable (linearly interpolated value of f(xi))

    real(pr),dimension(2),intent(in) :: x,f
    real(pr),intent(in) :: xi
    real(pr) :: fi

    ! note that this version assumes that xi has been bracketed already:
    ! x(1)<=xi<=x(2) so only uses two data points (although formula can
    ! be used to extrapolate outside of range in x)

    fi = f(1)+(f(2)-f(1))*(xi-x(1))/(x(2)-x(1))

  end function interp_linear_1d


  function interp_linear_2d(x,y,f,xi,yi) result(fi)
    ! INTERP_LINEAR_2D performs bilinear interpolation in 2D
    ! by consecutive linear interpolation in each direction
    ! 
    ! Modified: 23 December 2005

    implicit none

    ! I/O Parameters:
    ! X = data in x direction
    ! Y = data in y direction
    ! F = field values at points corresponding to x and y
    ! XI = location in x direction at which f is to be interpolated
    ! YI = location in y direction at which f is to be interpolated
    ! FI = interpolated value of f(xi,yi)

    real(pr),intent(in) :: x(2),y(2),f(2,2),xi,yi
    real(pr) :: fi
    
    ! Internal Parameters:
    ! FIY = field values after interpolation in y direction

    real(pr),dimension(2) :: fiy

    fiy(1) = interp_linear_1d(y,f(1,1:2),yi)
    fiy(2) = interp_linear_1d(y,f(2,1:2),yi)

    fi = interp_linear_1d(x,fiy,xi)

  end function interp_linear_2d


  function interp_linear_3d(x,y,z,f,xi,yi,zi) result(fi)
    ! INTERP_LINEAR_3D performs trilinear interpolation in 3D
    ! by consecutive linear interpolation in each direction
    ! 
    ! Modified: 23 December 2005

    implicit none

    ! I/O Parameters:
    ! X = data in x direction
    ! Y = data in y direction
    ! Z = data in z direction
    ! F = field values at points corresponding to x,y,z
    ! XI = location in x direction at which f is to be interpolated
    ! YI = location in y direction at which f is to be interpolated
    ! YI = location in z direction at which f is to be interpolated
    ! FI = interpolated value of f(xi,yi,zi)

    real(pr),intent(in) :: x(2),y(2),z(2),f(2,2,2),xi,yi,zi
    real(pr) :: fi
    
    ! Internal Parameters:
    ! FIZ = field values after interpolation in z direction
    ! FIZY = field values after interpolation in z and y directions

    real(pr) :: fiz(2,2),fizy(2)

    fiz(1,1) = interp_linear_1d(z,f(1,1,1:2),zi)
    fiz(1,2) = interp_linear_1d(z,f(1,2,1:2),zi)
    fiz(2,1) = interp_linear_1d(z,f(2,1,1:2),zi)
    fiz(2,2) = interp_linear_1d(z,f(2,2,1:2),zi)

    fizy(1) = interp_linear_1d(y,fiz(1,1:2),yi)
    fizy(2) = interp_linear_1d(y,fiz(2,1:2),yi)

    fi = interp_linear_1d(x,fizy,xi)

  end function interp_linear_3d


  subroutine newton1d(x,nmax,tolx,tolf,param,func)
    ! NEWTON1D solves equation f(x)=0 by Newton's method, for scalar x
    !
    ! Modified: 6 August 2010

    use io, only : warning

    implicit none

    ! I/O Parameters:
    ! NMAX = maximum number of iterations
    ! TOLX = tolerance for convergence of x
    ! TOLF = tolerance for convergence of f
    ! X = independent variable
    ! PARAM = parameters passed to function f(x)

    integer(pin),intent(in) :: nmax
    real(pr),intent(in) :: tolx,tolf
    real(pr),intent(inout) :: x
    integer(pin),dimension(:),intent(in) :: param
    
    interface
       subroutine func(x,f,dfdx,param)
         ! FUNC returns f(x) and the derivative df/dx
         use constants, only : pr,pin
         implicit none
         real(pr),intent(in) :: x
         real(pr),intent(out) :: f,dfdx
         integer(pin),dimension(:),intent(in) :: param
       end subroutine func
    end interface

    ! Internal Parameters:
    ! N = number of iterations
    ! F = function value f(x)
    ! DFDX = derivative df/dx
    ! DX = x(new)-x(old)
    ! STR = string used in output

    integer(pin) :: n
    real(pr) :: f,dx,dfdx
    character(64) :: str

    do n=1,nmax
       call func(x,f,dfdx,param)
       if (abs(f)<=tolf) return
       dx = -f/dfdx
       x = x+dx
       if (abs(dx)<=tolx) return
    end do

    str = "Warning: Newton's method did not converge"
    call warning(str,'newton')
    write(str,'(a,i14)') 'Iterations: nmax=',nmax
    call warning(str)
    write(str,'(a,e14.6,a,e14.6)') 'Function: errf=',abs(f),' tolf=',tolf
    call warning(str)
    write(str,'(a,e14.6,a,e14.6)') 'Independent variables: errx=',abs(dx),' tolx=',tolx
    call warning(str)

  end subroutine newton1d


  subroutine newton(x,nmax,tolx,tolf,param,func,xscale,fscale,errn,msg,method)
    ! NEWTON solves equation f(x)=0 by Newton's method, optionally 
    ! (although this is default) performing a backtracking line search
    ! along Newton direction to ensure global convergence.  
    ! See Section 6.3 in Dennis and Schnabel, Numerical Methods for 
    ! Unconstrained Optimization and Nonlinear Equations, 1996.
    !
    ! Modified: 6 August 2010

    use constants, only : zero,half,one
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! NMAX = maximum number of iterations
    ! TOLX = tolerance for convergence of x
    ! TOLF = tolerance for convergence of f
    ! X = independent variables
    ! PARAM = parameters passed to function f(x)
    ! XSCALE = scale for independent variables x
    ! FSCALE = scale for function f
    ! ERRN = error flag (true if error)
    ! MSG = message describing error
    ! METHOD = solution method, input (default is linesearch)
    !      standard = standard Newton's method
    !      linesearch = Newton's method with backtracking line search

    integer(pin),intent(in) :: nmax
    real(pr),intent(in) :: tolx,tolf
    real(pr),dimension(:),intent(inout) :: x
    integer(pin),dimension(:),intent(in) :: param
    real(pr),dimension(size(x)),intent(in),optional :: xscale,fscale
    logical,intent(out),optional :: errn
    character(*),intent(out),optional :: msg
    character(*),intent(in),optional :: method

    interface

       subroutine func(x,f,param,err,jac)
         ! FUNC returns f(x) and (optionally) the Jacobian df/dx
         use constants, only : pr,pin
         implicit none
         real(pr),dimension(:),intent(in) :: x
         real(pr),intent(out) :: f(size(x))
         integer(pin),dimension(:),intent(in) :: param
         logical,intent(out) :: err
         real(pr),intent(out),optional :: jac(size(x),size(x))
       end subroutine func

       ! uncomment to use LAPACK's linear solver
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
    ! N = number of iterations
    ! NX = length of x
    ! X0 = original x value
    ! A = fraction of Newton correction
    ! CMAX = maximum permissible (scaled) correction to x
    ! H = function to be minimized, h(x)=0.5*dot_prod(f(x),f(x))
    ! H0 = h(x0)
    ! F = function f(x)
    ! P = Newton correction
    ! G = gradient of h(x)
    ! XS = scale for x
    ! FS = scale for f
    ! JAC = Jacobian df/dx
    ! ERRX = error in x
    ! ERRF = error in f
    ! B = backtrack iteration
    ! BMAX = maximum number of backtrack iterations
    ! IPIV = index of pivot (used to diagnose errors)
    ! INFO = error flag (nonzero if error)
    ! ERR = error flag (true if error)
    ! STR = string used in output
    ! CHECK = flag indicating that solution needs to be verified
    ! MTHD = solution method
    !      standard = standard Newton's method
    !      linesearch = Newton's method with backtracking line search

    integer(pin) :: n,nx,b,info
    integer(pin),dimension(size(x)) :: ipiv
    real(pr) :: a,cmax,h,h0,errx,errf
    real(pr),dimension(size(x)) :: x0,f,p,g,xs,fs
    real(pr),dimension(size(x),size(x)) :: jac
    character(64) :: str,mthd
    logical :: err,check
    integer(pin),parameter :: bmax=100

    ! assign solution method

    if (present(method)) then
       mthd = method
    else
       mthd = 'linesearch'
    end if

    ! default error and message
    
    if (present(errn)) errn = .false.
    if (present(msg)) msg = ''

    ! size of system to be solved

    nx = size(x)

    ! set scale of x

    if (present(xscale)) then
       xs = abs(xscale)
       where (xs==zero) xs = one
    else
       xs = one
    end if

    ! set scale of f

    if (present(fscale)) then
       fs = abs(fscale)
       where (fs==zero) fs = one
    else
       fs = one
    end if

    ! evaluate function and Jacobian at initial guess

    call func(x,f,param,err,jac)
    
    if (err) then
       if (present(errn)) then
          errn = .true.
          if (present(msg)) msg = 'Function not defined at initial guess'
          return
       else
          call error('Function not defined at initial guess','newton')
       end if
    end if
     
    ! check for convergence on initial guess

    errf = maxval(abs(f/fs))
    if (errf<=0.01_pr*tolf) return
    
    ! bound search area

    cmax = 100._pr*max(norm2(x/xs),norm2(one/xs))

    ! iterate to improve guess at x

    do n = 1,nmax

       ! evaluate function to be minimized

       h = fmin(f,fs)

       ! save old values

       x0 = x
       h0 = h

       ! evaluate (scaled) gradient in h(x)

       g = matmul(f/fs**2,jac)

       ! solve Newton system for correction p

       ! uncomment to use LAPACK's linear solver
       p = -f
       call dgesv(nx,1,jac,nx,ipiv,p,nx,info)
       if (info/=0) then
          write(str,'(a,i5)') 'Error in DGESV, info=',info
          if (present(errn)) then
             errn = .true.
             if (present(msg)) msg = str
             return
          else
             call error(str,'newton')
          end if
       end if

       ! uncomment to use Gaussian elimination (not recommended)
       !call gauss(jac,-f,p)

       ! update x, calculate new f(x) and Jacobian

       select case(mthd)

       case default 

          call error('Invalid method','newton')

       case('standard')
          
          ! if |p|<tolx, then converged
          
          errx = maxval(abs(p)/xs)
          if (errx<=tolx) return

          ! correct x
          
          x = x0+p
          
          ! evaluate function and Jacobian at x
          
          call func(x,f,param,err,jac)
                    
          ! if function is not defined at x, then backtrack along Newton direction
          ! by factors of two until it is defined
          
          a = one
          do b = 1,bmax
             if (.not.err) exit
             a = half*a
             x = x0+a*p
             call func(x,f,param,err,jac)
          end do
          if (b==bmax) then
             if (present(errn)) then
                errn = .true.
                if (present(msg)) msg = 'Backtracked many times, but function never defined'
                return
             else
                call error('Backtracked many times, but function never defined','newton')
             end if
          end if

       case ('linesearch')
          
          ! perform line search to correct x, also evaluate f(x) and Jacobian at x

          call line_search(x0,h0,g,p,x,h,f,jac,cmax,xs,fs,tolx,check,func,param)

          ! check against false convergence by verifying that gradient of h is zero

          if (check) then
             check = (maxval(abs(g)*xs/max(h,half*size(x)))<tolx)
             if (check) then
                if (present(errn)) then
                   errn = .true.
                   if (present(msg)) msg = 'False convergence, gradient nonzero'
                   return
                else
                   call error('False convergence, gradient nonzero','newton')
                end if
             end if
          end if

          ! if |x-x0|<tolx, then converged

          errx = maxval(abs(x-x0)/xs)
          if (errx<=tolx) then
             return
          end if

       end select

       ! if |f|<tolf, then converged

       errf = maxval(abs(f/fs))
       if (errf<=tolf) then
          return
       end if

    end do

    if (present(errn)) then
       errn = .true.
       if (present(msg)) msg = "Newton's method did not converge"
    else
       call error("Newton's method did not converge",'newton')
    end if

  end subroutine newton


  function fmin(f,fs) result(h)
    ! FMIN returns scalar function that, when minimized, solves
    ! the nonlinear system of equations f=0.  Scaling is included.

    use constants, only : half

    implicit none

    ! I/O Parameters:
    ! F = vector of function values
    ! FS = scale of f
    ! H = function to be minimized

    real(pr),dimension(:),intent(in) :: f
    real(pr),dimension(size(f)),intent(in) :: fs
    real(pr) :: h

    ! Internal Parameters:
    ! FSCL = scaled f

    real(pr),dimension(size(f)) :: fscl

    fscl = f/fs
    h = half*dot_product(fscl,fscl)

  end function fmin


  subroutine line_search(x0,h0,g,p,x,h,f,jac,cmax,xs,fs,tolx,check,func,param)
    ! LINE_SEARCH is used with Newton's method.  It performs a backtracking 
    ! line search along the Newton direction until certain criteria are met.
    ! See Section 6.3 in Dennis and Schnabel, Numerical Methods for 
    ! Unconstrained Optimization and Nonlinear Equations, 1996.

    ! Modified: 6 August 2010

    use constants, only : zero,half,one,two,three
    use io, only : error,warning

    implicit none

    ! I/O Parameters:
    ! X0 = initial guess at solution
    ! H0 = initial value of scalar function to be minimized, h(x0)
    ! G = (scaled) gradient of h at x0
    ! P = search direction from x0 (Newton correction)
    ! X = improved solution
    ! H = updated value of scalar function to be minimized, h(x)
    ! F = vector function, f(x)
    ! JAC = Jacobian, df/dx
    ! CMAX = maximum permissible (scaled) correction to x
    ! XS = scale of x
    ! FS = scale for f
    ! TOLX = tolerance in x
    ! CHECK = flag used to indicate potential false convergence
    ! PARAM = parameters for function
    ! FUNC = function subroutine

    real(pr),dimension(:),intent(in) :: x0
    integer(pin),dimension(:),intent(in) :: param
    real(pr),dimension(size(x0)),intent(in) :: g,xs,fs
    real(pr),intent(in) :: cmax,h0,tolx
    real(pr),dimension(size(x0)),intent(inout) :: p
    real(pr),intent(out) :: h
    real(pr),dimension(size(x0)),intent(out) :: x,f
    real(pr),dimension(size(x0),size(x0)),intent(out) :: jac
    logical,intent(out) :: check
    
    interface

       subroutine func(x,f,param,err,jac)
         ! FUNC returns f(x) and (optionally) the Jacobian df/dx
         use constants, only : pr,pin
         implicit none
         real(pr),dimension(:),intent(in) :: x
         real(pr),intent(out) :: f(size(x))
         integer(pin),dimension(:),intent(in) :: param
         logical,intent(out) :: err
         real(pr),intent(out),optional :: jac(size(x),size(x))
       end subroutine func

    end interface

    ! Internal Parameters:
    ! N = iteration number
    ! NMAX = maximum number of iterations
    ! AL = constant used in step-acceptance test
    ! C = (scaled) length of corrective step
    ! ERRX = error in x
    ! A = fraction of corrective step
    ! In determining optimum value of a, we seek to minimize the function
    ! H(a) = h(x0+a*p), and will approximate H(a) by a cubic:
    ! H(a) = h0+h1*a+h2*a**2+h3*a**3
    ! H1 = coefficient of linear term
    ! H2 = coefficient of quadratic term
    ! H3 = coefficient of cubic term
    ! D = discriminant of quadratic equation to be solved
    ! AO = old/previous value of a
    ! AN = new/trial value of a
    ! HO = function to be minimized, evaluated at ao, H(ao)
    ! B = backtrack iteration
    ! BMAX = maximum number of backtrack iterations
    ! M = matrix used to solve linear system
    ! R = vector used to solve linear system
    ! ERR = error flag (nonzero if function not defined for chosen x)

    real(pr),parameter :: al = 1.d-4
    real(pr) :: c,errx,a,ao,an,ho,h1,h2,h3,D,m(2,2),r(2)
    integer(pin) :: n,b
    integer(pin),parameter :: bmax=100, nmax=100
    logical :: err

    ! enforce maximum limit on step length

    c = norm2(p/xs)
    if (c>cmax) p = p*cmax/c

    ! calculate linear coefficient [slope, H'(0)]

    h1 = dot_product(g,p)
    if (h1>zero) call error('Slope is positive--any corrections will increase error','line_search')

    ! start with full correction

    a = one

    ! iterate until permissible step is found

    do n = 1,nmax
       
       ! update x
       
       x = x0+a*p
       
       ! update f(x), Jacobian, and h(x)
       
       call func(x,f,param,err,jac)
       
       ! if function is not defined at x, then backtrack by factors of two until it is defined
       
       do b = 1,bmax
          if (.not.err) exit
          a = half*a
          x = x0+a*p
          call func(x,f,param,err,jac)
       end do
       if (b==bmax) call error('Backtracked many times, but function never defined','linesearch')
       
       h = fmin(f,fs)
              
       ! if |x-x0|<tolx, then converged (false convergence possible, so check)
       
       errx = maxval(abs(x-x0)/xs)
       if (errx<=tolx) then
          check = .true.
          return
       end if
       
       ! if h has sufficiently decreased, then converged
       
       if (h<=h0+al*a*h1) then
          check = .false.
          return
       end if
       
       ! backtrack by finding a that minimizes H(a)
       
       if (n==1) then ! first time, quadratic fit to H(a)
          
          ! H(a)=h0+h1*a+h2*a**2 with following information:
          ! H(0)=h0, H'(0)=h1, H(1)=h => h2=h-h0-h1, 
          ! minimizer given by:
          ! H'(an)=0=h1+2*h2*an => an=-h1/(2*h2)
          
          h2 = h-h0-h1
          an = -h1/(two*h2)
                    
       else ! subsequent backtracks, cubic fit to H(a)
          
          ! H(a)=h0+h1*a+h2*a**2+h3*a**3 with following information:
          ! H(0)=h0, H'(0)=h1, H(ao)=ho, H(a)=h
          ! => know h0 and h1, need to solve following for h2,h3:
          ! h0+h1*ao+h2*ao**2+h3*ao**3=ho
          ! h0+h1*a +h2*a**2 +h3*a**3 =h
          ! or in matrix form:
          ! [1 ao][h2]   [(ho-h0-h1*ao)/ao**2]
          ! [1 a ][h3] = [(h -h0-h1*a )/a**2 ]
          ! inverse of matrix on LHS times determinant is m,
          ! vector on RHS divided by determinant is r

          m(1,1) = a
          m(1,2) = -ao
          m(2,1) = -one
          m(2,2) = one
          r(1) = (ho-h0-ao*h1)/ao**2
          r(2) = (h -h0-a *h1)/a**2
          r = r/(a-ao) ! divide by determinant
          h2 = dot_product(m(1,:),r)
          h3 = dot_product(m(2,:),r)
          
          if (h3==zero) then

             ! cubic is really quadratic, minimizer given by
             ! H'(an)=0=h1+2*h2*an => an=-h1/(2*h2)

             an = -h1/(two*h2)

          else

             ! cubic with several possible minimizers given by
             ! H'(an)=0=h1+2*h2*an+3*h3*an**2 (quadratic)
             ! use discriminant 4*D of quadratic to distinguish types of solutions
             ! (note that D is really discriminant/4 to simplify formulas)

             D = h2**2-three*h3*h1

             if (D<zero) then
                ! two complex conjugate solutions (not useful here),
                ! so just decrease a by factor of two
                an = half*a
             elseif (h2<zero) then
                ! minimizing solution to cubic (making sure an>0) 
                an = (-h2+sqrt(D))/(three*h3)
             else
                ! use this formula to ensure an>0
                an = -h1/(h2+sqrt(D))
             end if

          end if
                    
          ! must decrease a by at least a factor of two

          if (an>half*a) an = half*a
                    
       end if
       
       ! save these values, since they are used in cubic fit

       ao = a
       ho = h
       
       ! cannot decrease a by more than a factor of ten

       a = max(an,0.1_pr*a)

    end do

    call warning('Maximum number of iterations reached','linesearch')

  end subroutine line_search


  subroutine jacfd(x,f0,param,func,jac)
    ! JACFD estimates the Jacobian of a function by forward finite differences,
    ! based on SLATEC library's dfdjc1.f and Numerical Recipes' fdjac
    !
    ! Modified: 6 August 2010

    use constants, only : zero,one

    implicit none

    ! I/O Parameters:
    ! X = independent variables
    ! F0 = f(x) evaluated at input value of x
    ! PARAM = parameters passed to function f(x)
    ! JAC = estimate of Jacobian df/dx

    real(pr),dimension(:),intent(inout) :: x
    real(pr),dimension(:),intent(in) :: f0
    integer(pin),dimension(:),intent(in) :: param
    real(pr),dimension(size(x),size(x)),intent(out) :: jac

    interface

       subroutine func(x,f,param)
         ! FUNC returns f(x)
         use constants, only : pr,pin
         implicit none
         real(pr),dimension(:),intent(in) :: x
         real(pr),intent(out) :: f(size(x))
         integer(pin),dimension(:),intent(in) :: param
       end subroutine func

    end interface

    ! Internal Parameters:
    ! I = index of x
    ! N = size of x
    ! F = function values f(x)
    ! TEMP = temporary value of x
    ! DX = step size used in finite-difference approximation to derivative
    
    integer(pin) :: i,n
    real(pr) :: temp,dx,eps
    real(pr),dimension(size(x)) :: f

    ! set eps to square-root of machine precision

    eps = sqrt(epsilon(one))

    n = size(x)

    ! calculate Jacobian with forward finite-differences

    do i = 1,n
       temp = x(i)
       dx = eps*abs(temp)
       if (dx==zero) dx = eps
       x(i) = temp+dx
       dx = x(i)-temp
       call func(x,f,param)
       x(i) = temp
       jac(:,i) = (f-f0)/dx
    end do

  end subroutine jacfd


  function cumquad_trap(x,y) result(f)
    ! CUMQUAD_TRAP returns an approximation to the cumulative integral of y with respect to x
    ! from min(x) to each value of x using trapezoidal rule
    !
    ! Modified: 23 September 2006

    use constants, only : zero,half

    implicit none

    ! I/O Parameters:
    ! X = integration variable
    ! Y = function values
    ! F = cumulative integral

    real(pr),dimension(:),intent(in) :: x
    real(pr),dimension(size(x)),intent(in) :: y
    real(pr),dimension(size(x)) :: f
    
    ! Internal Parameters:
    ! I = counter for x
    ! N = number of data points
    ! H = grid spacing

    integer(pin) :: i,n
    real(pr),dimension(size(x)-1) :: h

    n = size(x)
    h = x(2:n)-x(1:n-1)

    f(1) = zero

    do i = 2,n
       f(i) = f(i-1)+half*h(i-1)*(y(i-1)+y(i))
    end do

  end function cumquad_trap


  subroutine gauss(A,b,x)
    ! GAUSS solves linear system Ax=b for x using Gaussian elimination with pivoting
    ! 
    ! Modified: 3 May 2007

    use constants, only : one
    use io, only : error

    implicit none

    real(pr),dimension(:,:),intent(in) :: A
    real(pr),dimension(:),intent(in) ::  b
    real(pr),dimension(:),intent(out) :: x
    
    integer(pin) :: i,j   ! Local index variables.
    integer(pin) :: N      ! Order of the linear system.
    real(pr) :: m         ! Multiplier.
    real(pr) :: smallestPivot
    
    ! Pointers to the appropriate rows of the matrix during the elmination.
    real(pr),dimension(:), pointer :: pivotRow
    real(pr),dimension(:), pointer :: currentRow
    
    ! Copies of the input arguments.  These copies are modified during
    ! the computation.
    ! The target attribute is used to indicate that the specified 
    ! variable may be the target of a pointer.  Rows of ACopy are targets
    ! of pivotRow and currentRow, defined above.
    
    real(pr),dimension(size(A,1),size(A,2)), target :: ACopy
    real(pr),dimension(size(b)) :: bCopy
    
    ! Status of the computation.  The return value of the function.
    logical :: successful
    
    smallestPivot = epsilon(one)
    N = size(b)   
    ACopy = A
    bCopy = b
    successful = .true.
    
    ! put system into upper triangular form

    i = 1
    do while (successful.and.(i<=N-1))
       
       pivotRow => ACopy(i,:)
       
       successful = (abs(pivotRow(i))>=smallestPivot)
       
       if (successful) then
          
          ! eliminate entries in pivot column below pivot row
          
          do j = i+1,N
             currentRow => ACopy(j,:)
             ! calculate multiplier
             m = currentRow(i)/pivotRow(i) 
             ! perform elimination step on currentRow and right hand side, bCopy
             currentRow = m*pivotRow-currentRow
             bCopy(j) = m*bCopy(i)-bCopy(j)
          end do
          
       end if
       
       ! next row
       i = i+1
       
    end do
    
    ! check last pivot
    pivotRow => ACopy(N,:)
    if (successful) successful = (abs(pivotRow(N))>=smallestPivot)
    
    if (successful) then
       do i = N,2,-1   ! backward substitution
          ! determine ith unknown, x(i).
          x(i) = bCopy(i)/ACopy(i,i)
          ! substitute known value of x(i), reducing order of system by 1
          bCopy = bCopy-x(i)*ACopy(:,i)
       end do
    end if
    
    ! determine the value of x(1) as special case
    if (successful) x(1) = bCopy(1)/ACopy(1,1 )
    
    ! return value of function
    if (.not.successful) then

       print *, abs(pivotRow(N)),' < ',smallestPivot
       print *
       do i = 1,N
          write(6,'(2f20.10,a,f20.10)') A(i,:),'        ',b(i)
       end do

       call error('Solving linear system with Gaussian elimination failed','gauss')
    
    end if

  end subroutine gauss


  subroutine tridiag(n,l,d,u,b)
    ! TRIDIAG solves tridiagonal system Ax=b for x, where A is tridiagonal matrix
    ! with diagonal d, superdiagonal u, and subdiagonal l.
    ! 
    ! Modified: 27 August 2006

    implicit none

    ! I/O Parameters:
    ! N = length of x
    ! L = subdiagonal elements of tridiagonal matrix
    ! D = diagonal elements of tridiagonal matrix
    ! U = superdiagonal elements of tridiagonal matrix
    ! Note that for the linear system, the index labels the row of the system

    integer(pin),intent(in) :: n
    real(pr),dimension(2:n),intent(inout) :: l
    real(pr),dimension(1:n),intent(inout) :: d,b
    real(pr),dimension(1:n-1),intent(inout) :: u

    ! Internal Parameters:
    ! I = loop index
    ! M = parameter combination

    integer(pin) :: i
    real(pr) :: m

    ! forward

    do i = 2,n
       m = l(i)/d(i-1)
       d(i) = d(i)-m*u(i-1)
       b(i) = b(i)-m*b(i-1)
    end do

    ! backward

    b(n) = b(n)/d(n)

    do i = n-1,1,-1
       b(i) = (b(i)-u(i)*b(i+1))/d(i)
    end do

  end subroutine tridiag


  subroutine init_spline(n,x,y,c,dydxm,dydxp)
    ! INIT_SPLINE calculates second derivatives of cubic spline.  If first derivatives
    ! are given for the endpoints, these values are used to set the boundary conditions.
    ! Otherwise, natural splines (with vanishing second derivatives) are used.
    !
    ! Modified: 1 February 2007

    use constants, only : zero,one,two,three

    implicit none

    ! I/O Parameters:
    ! N = length of x
    ! X = values at which function y(x) is tabulated
    ! Y = values of function y(x)
    ! C = one-half the second derivative of y with respect to x, d^2y/dx^2,
    ! also used for right-hand side of tridiagonal system
    ! DYDXM = dy/dx at first endpoint
    ! DYDXP = dy/dx at last endpoint

    integer(pin),intent(in) :: n
    real(pr),dimension(n),intent(in) :: x,y
    real(pr),dimension(n),intent(out) :: c
    real(pr),intent(in),optional :: dydxm,dydxp

    ! Internal Parameters:
    ! H = grid spacings, h(i)=x(i+1)-x(i)
    ! L = subdiagonal elements of tridiagonal matrix
    ! D = diagonal elements of tridiagonal matrix
    ! U = superdiagonal elements of tridiagonal matrix
    ! Note that for the linear system, the index labels the row of the system

    real(pr),dimension(1:n-1) :: h
    real(pr),dimension(3:n-1) :: l
    real(pr),dimension(2:n-1) :: d
    real(pr),dimension(2:n-2) :: u

    ! calculate grid spacing

    h(1:n-1) = x(2:n)-x(1:n-1)
    
    ! generate linear system to be solved

    u(2:n-2) = h(2:n-2)
    d(2:n-1) = two*(h(1:n-2)+h(2:n-1))
    l(3:n-1) = h(2:n-2)

    c(2:n-1) = three*(y(3:n)-y(2:n-1))/h(2:n-1)-three*(y(2:n-1)-y(1:n-2))/h(1:n-2)

    ! modify i=2 equation if first derivative is given

    if (present(dydxm)) then
       d(2) = 1.5_pr*h(1)+two*h(2)
       c(2) = three*(y(3)-y(2))/h(2)-1.5_pr*(three*(y(2)-y(1))/h(1)-dydxm)
    end if

    ! modify i=n-1 equation if first derivative is given

    if (present(dydxp)) then
       d(n-1) = two*h(n-2)+1.5_pr*h(n-1)
       c(n-1) = 1.5_pr*(three*(y(n)-y(n-1))/h(n-1)-dydxp)-three*(y(n-1)-y(n-2))/h(n-2)
    end if

    ! solve tridiagonal system

    call tridiag(n-2,l,d,u,c(2:n-1))

    ! complete solution by specifying second derivatives at endpoints

    if (present(dydxm)) then
       c(1) = (three*(y(2)-y(1))/h(1)-three*dydxm-c(2)*h(1))/(two*h(1))
    else
       c(1) = zero
    end if

    if (present(dydxp)) then
       c(n) = -(three*(y(n)-y(n-1))/h(n-1)-three*dydxp+c(n-1)*h(n-1))/(two*h(n-1))
    else
       c(n) = zero
    end if

  end subroutine init_spline


  function interp_spline_scalar(x,y,c,xi,ix) result(yi)
    ! INTERP_SPLINE_SCALAR evaluates a cubic spline y=a+b*x+c*x**2+d*x**3
    !
    ! Modified: 23 September 2006
    
    use constants, only : two,three
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! X = values at which function y(x) is tabulated
    ! Y = values of function y(x)
    ! C = coefficient of quadratic term
    ! XI = value of x at which spline is to be evaluated
    ! IX = guess at index i satisfying x(i)<=xi<x(i+1)
    ! YI = approximation of y(xi)

    real(pr),dimension(:),intent(in) :: x
    real(pr),dimension(size(x)),intent(in) :: y
    real(pr),dimension(size(x)),intent(in) :: c
    real(pr),intent(in) :: xi
    integer(pin),intent(in),optional :: ix
    real(pr) :: yi

    ! Internal Parameters:
    ! I = index satisfying x(i)<=xi<x(i+1)
    ! H = grid spacing
    ! B = coefficient of linear term
    ! D = coefficient of cubic term
    ! OTS = flag denoting if xi lies outside range of x

    integer(pin) :: i
    real(pr) :: h,b,d
    logical :: ots

    ! obtain index of interval containing xi

    if (present(ix)) then
       i = search_binary(x,xi,'floor',ots,ix)
    else
       i = search_binary(x,xi,'floor',ots)
    end if

    if (ots) call error('Error: Attempting to extrapolate data outside range','eval_spline')

    ! compute coefficients and evaluate cubic function

    h = x(i+1)-x(i)
    b = (y(i+1)-y(i))/h-h/three*(c(i+1)+two*c(i))
    d = (c(i+1)-c(i))/(three*h)

    yi = y(i)+b*(xi-x(i))+c(i)*(xi-x(i))**2+d*(xi-x(i))**3

  end function interp_spline_scalar


  function interp_spline_vector(x,y,c,xi) result(yi)
    ! INTERP_SPLINE_VECTOR evaluates a cubic spline y=a+b*x+c*x**2+d*x**3
    ! for a vector xi 
    !
    ! Modified: 23 September 2006
    
    use constants, only : two,three
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! X = values at which function y(x) is tabulated
    ! Y = values of function y(x)
    ! C = coefficient of quadratic term
    ! XI = values of x at which spline is to be evaluated
    ! YI = approximation of y(xi)

    real(pr),dimension(:),intent(in) :: x
    real(pr),dimension(size(x)),intent(in) :: y
    real(pr),dimension(size(x)),intent(in) :: c
    real(pr),dimension(:),intent(in) :: xi
    real(pr),dimension(size(xi)) :: yi

    ! Internal Parameters:
    ! I = index of xi vector
    ! J = index satisfying x(j)<=xi<x(j+1)
    ! H = grid spacing
    ! B = coefficient of linear term
    ! D = coefficient of cubic term
    ! OTS = flag denoting if xi lies outside range of x

    integer(pin) :: i
    integer(pin),dimension(size(xi)) :: j
    real(pr),dimension(size(xi)) :: h,b,d
    logical :: ots

    ! find index j such that x(j)<=xi(i)<x(j+1)

    do i = 1,size(xi)

       if (i==1) then
          ! use binary search to find starting value
          j(i) = search_binary(x,xi(i),'floor',ots)
       else
          ! use previous index as starting value for binary search
          j(i) = search_binary(x,xi(i),'floor',ots,j(i-1))
       end if

       if (ots) call error('Error: Attempting to extrapolate data outside range','eval_spline')

    end do

    ! compute coefficients and evaluate cubic function

    h = x(j+1)-x(j)
    b = (y(j+1)-y(j))/h-h/three*(c(j+1)+two*c(j))
    d = (c(j+1)-c(j))/(three*h)

    yi = y(j)+b*(xi-x(j))+c(j)*(xi-x(j))**2+d*(xi-x(j))**3

  end function interp_spline_vector


  function cumquad_spline(x,y,c) result(f)
    ! CUMQUAD_SPLINE returns an approximation to the cumulative integral of y with respect to x
    ! from min(x) to each value of x using cubic splines
    !
    ! Modified: 23 September 2006

    use constants, only : zero,one,two,three,four

    implicit none

    ! I/O Parameters:
    ! X = integration variable
    ! Y = function values
    ! C = coefficient of quadratic term
    ! F = cumulative integral

    real(pr),dimension(:),intent(in) :: x
    real(pr),dimension(size(x)),intent(in) :: y,c
    real(pr),dimension(size(x)) :: f
    
    ! Internal Parameters:
    ! I = counter for x
    ! N = number of data points
    ! H = grid spacing
    ! B = coefficient of linear term
    ! D = coefficient of cubic term

    integer(pin) :: i,n
    real(pr),dimension(size(x)-1) :: h,b,d

    n = size(x)
    h = x(2:n)-x(1:n-1)
    b = (y(2:n)-y(1:n-1))/h-h/three*(c(2:n)+two*c(1:n-1))
    d = (c(2:n)-c(1:n-1))/(three*h)

    f(1) = zero

    do i = 2,n
       f(i) = f(i-1)+ &
            h(i-1)         *y(i-1)+ &
            h(i-1)**2/two  *b(i-1)+ &
            h(i-1)**3/three*c(i-1)+ &
            h(i-1)**4/four *d(i-1)
    end do

  end function cumquad_spline


  function log2(x)
    ! LOG2 returns base two logarithm
    !
    ! Modified: 17 February 2006

    use constants, only : logtwo

    implicit none

    real(pr),intent(in) :: x
    real(pr) :: log2

    log2 = log(x)/logtwo

  end function log2


  function arcsinh(x) result(f)
    ! ARCSINH returns the inverse hyperbolic sine of x
    !
    ! Modified: 23 September 2006

    use constants,only : one

    implicit none

    ! I/O Parameters:
    ! X = argument of arcsinh
    ! F = arcsinh(x)

    real(pr),intent(in) :: x
    real(pr) :: f

    ! Internal Parameters:
    ! Y = absolute value of x

    real(pr) :: y

    y = abs(x)
    f = log(y+sqrt(y**2+one))
    f = sign(f,x)

  end function arcsinh


  function arccosh(x) result(f)
    ! ARCCOSH returns the inverse hyperbolic cosine of x
    !
    ! Modified: 23 September 2006

    use constants, only : one
    use io, only : error

    implicit none

    ! I/O Parameters:
    ! X = argument of arccosh
    ! F = arccosh(x)

    real(pr),intent(in) :: x
    real(pr) :: f

    if (x<one) call error('Error: argument of arccosh cannot be less than 1','arccosh')

    f = log(x+sqrt(x**2-one))

  end function arccosh

  
  function norm2(x) result (n)
    ! NORM2 returns the two-norm (length) of a vector x
    !
    ! Modified: 25 September 2006

    ! I/O Parameters:
    ! X = vector
    ! N = norm (length) of vector

    real(pr),dimension(:),intent(in) :: x
    real(pr) :: n

    n = sqrt(dot_product(x,x))

  end function norm2


  pure function convolve(C,H,term) result(CH)
    ! CONVOLVE evaluates convolution of C and H
    !
    ! Modified: 28 December 2006

    use constants, only : half

    implicit none

    ! I/O Parameters:
    ! C = convolution kernel
    ! H = history
    ! TERM = term in convolution to compute (history or current step)
    ! CH = convolution of C and H
    
    real(pr),dimension(:),intent(in) :: C
    complex(pr),dimension(:),intent(in) :: H
    character(*),intent(in) :: term
    complex(pr) :: CH

    ! Internal Parameters:
    ! N = length of C and H

    integer(pin) :: n

    select case(term)
    case('history')
       n = size(C)
       CH = dot_product(C(1:n-1),H(1:n-1))+half*C(n)*H(n)
    case('current')
       CH = half*C(1)*H(1)
    end select

  end function convolve


  function word_count(str) result(nw)
    ! WORD_COUNT counts the number of words (delimited by spaces) in a string
    !
    ! Modified: 16 June 2007

    implicit none
    
    ! I/O Parameters:
    ! STR = string to be processed
    ! NW = number of words in string

    character(*),intent(in) :: str
    integer(pin) :: nw
    
    ! Internal Parameters:
    ! I = location of space in string
    ! C = current position in string

    integer(pin) :: i,c

    nw = 0
    c = 1
    do
       i = index(trim(str(c:)),' ')
       nw = nw+1
       if (i==0) return
       c = c+i
    end do

  end function word_count


  subroutine extract_words(str,words)
    ! EXTRACT_WORDS extracts space-delimited words from string
    !
    ! Modified: 16 June 2007
    
    implicit none

    ! I/O Parameters:
    ! STR = string to be processed
    ! WORDS = extracted words

    character(*),intent(in) :: str
    character(*),dimension(:),intent(out) :: words
    
    ! Internal Parameters:
    ! I = location of space in string
    ! N = length of trimmed string
    ! C = current position in string
    ! NW = number of words

    integer(pin) :: i,n,c,nw

    n = len_trim(str)

    nw = 0
    c = 1
    do
       i = index(trim(str(c:)),' ')
       nw = nw+1
       if (i==0) then
          words(nw) = str(c:n)
          return
       else
          words(nw) = str(c:c+i-2)
          c = c+i
       end if
    end do

  end subroutine extract_words
  

  elemental subroutine convert_time(total_sc,hr,mn,sc)

    real(pr),intent(in) :: total_sc
    integer(pin),intent(out) :: hr,mn
    real(pr),intent(out) :: sc

    hr = total_sc/3600
    mn = (total_sc-3600*hr)/60
    sc = total_sc-3600d0*dble(hr)-60d0*dble(mn)

  end subroutine convert_time


end module utilities
