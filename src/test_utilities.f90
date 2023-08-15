! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module test_utilities
  
  use constants, only : pr,pin

  implicit none


contains


  subroutine test_gauss

    use constants, only : zero
    use utilities, only : gauss

    integer(pin),parameter :: n=6
    integer(pin) :: i
    real(pr) :: A(n,n),b(n),x(n)

    call random_number(A)
    call random_number(b)

    print *
    write(6,'(a)') 'A='
    do i = 1,n
       write(6,'(6f20.10)') A(i,:)
    end do

    print *
    write(6,'(a)') 'b='
    do i = 1,n
       write(6,'(f20.10)') b(i)
    end do

    call gauss(A,b,x)

    print *
    write(6,'(a)') 'x='
    do i = 1,n
       write(6,'(f20.10)') x(i)
    end do

  end subroutine test_gauss


  subroutine test_update

    use constants, only : zero,half,one
    use integration, only : update_field

    implicit none

    integer(pin),parameter :: nx=2048, ny=1, nr=2, nt=10000000
    real(pr),parameter :: dt=one, R(nr)=half
    real(pr),dimension(nx,ny) :: f0,f
    real(pr),dimension(nx,ny,nr) :: dfdt

    integer(pin) :: n

    f = zero
    dfdt = one

    do n = 1,nt
       f0 = f
       call update_field(nx,ny,nr,f,f0,dfdt,dt,R)
    end do

  end subroutine test_update


  subroutine test_newton
    ! TEST_NEWTON tests line search Newton's method with the example given by
    ! Dennis and Schnabel, Numerical Methods for Uncontrained Optimization and Nonlinear
    ! Equations, 1996 on p. 149-151 (example 6.5.1).  Also test regular Newton's method
    ! using the same equation but with different starting value (example 5.4.2, p. 98).
    ! Perform second test where function has compact support (and is undefined elsewhere).
    ! 
    ! Modified: 10 August 2010
    
    use constants, only : half,one,two,three
    use utilities, only : newton
    
    implicit none
    
    ! Internal Parameters:
    ! NEWTON_X = independent variables for Newton's method
    ! PARAM = parameters required for function that is solved using Newton's method
    ! NEWTON_NMAX = maximum number of interations allowed to solve equations with Newton's method
    ! NEWTON_TOLX = maximum tolerated error in independent variables
    ! NEWTON_TOLF = maximum tolerated error in dependent variables
    
    real(pr) :: newton_x(2),newton_tolx,newton_tolf
    integer(pin) :: newton_nmax,param(0)

    newton_nmax = 1000
    newton_tolx = epsilon(one)
    newton_tolf = epsilon(one)**(two/three)
    
    write(6,'(a)') ''
    write(6,'(a)') "Newton's method, test 1:"
    write(6,'(a,f20.10,a,f20.10)') 'exact:         ',one,'  ',one
    
    newton_x(1) = two
    newton_x(2) = half
    
    call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param,fcn_test1)
    
    write(6,'(a,f20.10,a,f20.10)') 'numerical1:    ',newton_x(1),'  ',newton_x(2)
    
    newton_x(1) = two
    newton_x(2) = three
    
    call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param,fcn_test1)
    
    write(6,'(a,f20.10,a,f20.10)') 'numerical2:    ',newton_x(1),'  ',newton_x(2)
    
    write(6,'(a)') ''
    write(6,'(a)') "Newton's method, test 2:"
    write(6,'(a,f20.10,a,f20.10)') 'exact:         ',0.999876605423074_pr,'  ',0.000123394576925610_pr
    
    newton_x(1) = one
    newton_x(2) = one
    
    call newton(newton_x,newton_nmax,newton_tolx,newton_tolf,param,fcn_test2)
    
    write(6,'(a,f20.10,a,f20.10)') 'numerical:     ',newton_x(1),'  ',newton_x(2)
    write(6,'(a)') ''
    
  end subroutine test_newton
  
  
  subroutine fcn_test1(x,f,param,err,jac)
    ! FCN_TEST1 is a function used to test line search Newton's method
    !
    ! Modified: 6 August 2010
    
    use constants, only : one,two,three
  
    implicit none
    
    ! I/O Parameters:
    ! X = independent variable for Newton's method
    ! F = dependent variable for Newton's method
    ! PARAM = variables required for evaluating elastodynamic equation and strength
    ! ERR = error flag (non-zero if function not defined)
    ! JAC = Jacobian df/dx
    
    real(pr),dimension(:),intent(in) :: x
    real(pr),intent(out) :: f(size(x))
    integer(pin),dimension(:),intent(in) :: param
    logical,intent(out) :: err
    real(pr),intent(out),optional :: jac(size(x),size(x))
    
    err = .false.
    
    f(1) = x(1)**2+x(2)**2-two
    f(2) = exp(x(1)-one)+x(2)**3-two
    
    if (present(jac)) then
       
       jac(1,1) = two*x(1)
       jac(1,2) = two*x(2)
       jac(2,1) = exp(x(1)-one)
       jac(2,2) = three*x(2)**2
       
    end if
    
  end subroutine fcn_test1
  
  
  subroutine fcn_test2(x,f,param,err,jac)
    ! FCN_TEST2 is a function used to test bounded domain Newton's method
    !
    ! Modified: 6 August 2010
    
    use constants, only : zero,one
    
    implicit none
    
    ! I/O Parameters:
    ! X = independent variable for Newton's method
    ! F = dependent variable for Newton's method
    ! PARAM = variables required for evaluating elastodynamic equation and strength
    ! ERR = error flag (non-zero if function not defined)
    ! JAC = Jacobian df/dx
    
    real(pr),dimension(:),intent(in) :: x
    real(pr),intent(out) :: f(size(x))
    integer(pin),dimension(:),intent(in) :: param
    logical,intent(out) :: err
    real(pr),intent(out),optional :: jac(size(x),size(x))
    
    if (x(2)<=zero) then
       err = .true.
       return
    else
       err = .false.
    end if
    
    f(1) = -x(1)+one-x(2)
    f(2) = -x(1)+10_pr+log(x(2))
    
    if (present(jac)) then
       
       jac(1,1) = -one
       jac(1,2) = -one
       jac(2,1) = -one
       jac(2,2) = one/x(2)
       
    end if
    
  end subroutine fcn_test2
  

  subroutine test_interp

    use constants, only : one
    use utilities, only : init_spline,interp_linear,interp_spline, &
         cumquad_trap,cumquad_spline

    implicit none

    integer(pin) :: n,i
    real(pr),dimension(:),allocatable :: x,y,c,yct,ycs,yc
    real(pr) :: xi(3),yl(3),ys(3),ylm(3),ysm(3)
    real(pr),parameter :: dydxm=0.00125664_pr,dydxp=0.0001_pr

    ! test 1: natural spline

    n = 15
    allocate(x(n),y(n),c(n))

    x = real((/(i,i=1,n)/),pr)
    y = real((/ 7,6,4,4,5,4,2,3,5,7,6,4,4,5,7 /),pr)

    call init_spline(n,x,y,c)
    
    write(6,*)
    write(6,'(a)') 'Interpolation routines: Test 1 (linear and natural spline):'

    !write(6,*)
    !write(6,'(a)') 'x y c'

    !do i = 1,n
    !   write(6,'(3f20.10)') x(i),y(i),2._pr*c(i)
    !end do
    
    xi = (/ 1.2,7.31,12.2 /)
    ylm = (/ 6.8,2.31,4. /)

    yl = interp_linear(x,y,xi)
    ys = interp_spline(x,y,c,xi)

    write(6,*)
    write(6,'(a)') '        x                   y(x)-linear         y(x)-matlab         y(x)-spline         y(x)-matlab'
    do i = 1,3
       write(6,'(5f20.10)') xi(i),yl(i),ylm(i),ys(i),ysm(i)
    end do

    deallocate(x,y,c)

    ! test 2: fixed derivative spline

    write(6,*)
    write(6,'(a)') 'Interpolation Routines: Test 2 (linear and fixed derivative spline)'
    write(6,*)

    n = 9
    allocate(x(n),y(n),c(n))

    x = (/ 0.,8.2,14.7,17.0,21.1,35.0,54.1,104.,357. /)
    y = (/ 0.,0.5,1.,1.1,1.2,1.4,1.5,1.6,1.7 /)

    call init_spline(n,x,y,c,dydxm,dydxp)

    !write(6,*)
    !write(6,'(a)') 'x y c'

    !do i = 1,n
    !   write(6,'(3f20.10)') x(i),y(i),2._pr*c(i)
    !end do

    xi = (/ 1.2,79.1,121. /)
    ylm = (/ 0.07317073170732,1.55010020040080,1.60671936758893 /)
    ysm = (/ 0.01667585025790,1.55847653649851,1.62246538476981 /)

    yl = interp_linear(x,y,xi)
    ys = interp_spline(x,y,c,xi)

    write(6,'(a)') '        x                   y(x)-linear         y(x)-matlab         y(x)-spline         y(x)-matlab'
    do i = 1,3
       write(6,'(5f20.10)') xi(i),yl(i),ylm(i),ys(i),ysm(i)
    end do

    deallocate(x,y,c)

    ! test 3: natural spline, cumulative quadrature

    n = 10
    allocate(x(n),y(n),c(n),yct(n),ycs(n),yc(n))

    x = real((/(i,i=0,n-1)/),pr)
    y = exp(-x)

    call init_spline(n,x,y,c)

    yct = cumquad_trap(x,y)
    ycs = cumquad_spline(x,y,c)
    yc = one-y

    write(6,*)
    write(6,'(a)') 'Interpolation Routines: Test 3 (cumulative quadrature with trapezoidal rule and natural spline)'
    write(6,*)

    write(6,'(a)') '        x                   Y(x)-trap           Y(x)-spline         Y(x)-exact'
    do i = 1,n
       write(6,'(4f20.10)') x(i),yct(i),ycs(i),yc(i)
    end do

    deallocate(x,y,c,yct,ycs,yc)

  end subroutine test_interp



end module test_utilities
