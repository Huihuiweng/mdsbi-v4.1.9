! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module convolution_im

  ! CONVOLUTION_IM contains subroutines to perform the convolution over 
  ! past time of the slip or slip velocity history with the appropriate convolution
  ! kernel.  This version is for identical materials.
  ! 
  ! Modified: 14 September 2007

  use constants, only : pr,pin

  implicit none


contains


  subroutine do_convolution_im(mdl,krn,cnv,update_history)
    ! DO_CONVOLUTION_IM performs convolution, identical materials case
    ! 
    ! Modified: 12 June 2007

    use model, only : model_type
    use kernel, only : kernel_type
    use convolution, only : convolution_type

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! KRN = kernel variables
    ! CNV = convolution variables
    ! UPDATE_HISTORY = perform convolution over past history instead of just updating the portion 
    ! of the convolution pertaining to the current time step

    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(in) :: krn
    type(convolution_type),intent(inout) :: cnv
    logical,intent(in) :: update_history

    ! Internal Parameters:
    ! I = index in the kx direction
    ! J = index in the ky direction
    ! Note that for 2D mixed mode case, both kx and ky are indexed by i.  In this case the x in kx 
    ! does not coincide with the x direction in space (likewise with y).  See the user guide for 
    ! further explanation.

    integer(pin) :: i,j

    ! loop over Fourier modes and perform appropriate convolution
    ! for each mode

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)

    !$OMP DO

    do j = 1,cnv%nky

       do i = cnv%mkx,cnv%pkx

          ! perform convolution on this mode

          call do_single_convolution(i,j,mdl,krn,cnv,update_history)

       end do

    end do

    !$OMP END DO

    !$OMP END PARALLEL
    
  end subroutine do_convolution_im


  subroutine do_single_convolution(i,j,mdl,krn,cnv,update_history)
    ! DO_SINGLE_CONVOLUTION performs convolution, identical materials case
    ! 
    ! Modified: 13 September 2007

    use constants, only : zeroc,zero,one
    use model, only : model_type
    use kernel, only : kernel_type
    use convolution, only : convolution_type,wavenumber
    
    implicit none

    ! I/O Parameters:
    ! I = index in the kx direction
    ! J = index in the ky direction
    ! Note that for 2D mixed mode case, both kx and ky are indexed by i.  In this case the x in kx 
    ! does not coincide with the x direction in space (likewise with y).  See the user guide for 
    ! further explanation.
    ! MDL = model variables
    ! KRN = kernel variables
    ! CNV = convolution variables
    ! UPDATE_HISTORY = perform convolution over past history instead of just updating the portion of
    ! the convolution pertaining to the current time step

    integer(pin),intent(in) :: i,j
    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(in) :: krn
    type(convolution_type),intent(inout) :: cnv
    logical,intent(in) :: update_history

    ! Internal Parameters:
    ! LB = lower index/bound of history vectors Hx and Hy
    ! UB = upper index/bound of history vectors Hx and Hy
    ! FXS = static contribution to stress transfer functional Fx
    ! FYS = static contribution to stress transfer functional Fy
    ! FXC = contribution of current time step to Fx
    ! FYC = contribution of current time step to Fy
    ! K = wavenumber amplitude (k**2=kx**2+ky**2)
    ! KM = wavenumber product array 
    !      km(1,1) = kx**2
    !      km(1,2) = km(2,1) = kx*ky
    !      km(2,2) = ky**2
    ! KX = wavenumber in x direction
    ! KY = wavenumber in y direction
    ! TERM = term in convolution to compute (history or current step)

    integer(pin) :: lb,ub
    complex(pr) :: FxS,FyS,FxC,FyC
    real(pr) :: k,km(2,2),kx,ky
    character(64) :: term

    ! special case for spring-block model

    if (mdl%dim==0) then

       cnv%Fx = -mdl%M*cnv%Dx
       cnv%Fy = -mdl%M*cnv%Dy

       return

    end if

    ! calculate wavenumber array
    
    call wavenumber(mdl,cnv,i,j,k,km,kx,ky)

    ! treat k=0 case separately to avoid useless work and possible division by zero
    
    if (k==zero) then

       cnv%Fx(i,j) = zeroc
       cnv%Fy(i,j) = zeroc

       return

    end if
    
    ! calculate static contribution to stress transfer
    
    if (cnv%formulation=='displacement') then
       FxS = zeroc
       FyS = zeroc
    else
       call conv_static(mdl,krn,k,km,cnv%Dx(i,j),cnv%Dy(i,j),FxS,FyS)
    end if

    ! sum contributions to stress transfer and return if static/quasidynamic (no history contribution)

    if (cnv%formulation=='static'.or.cnv%formulation=='quasidynamic') then

       cnv%Fx(i,j) = FxS
       cnv%Fy(i,j) = FyS

       return

    end if

    ! calculate contribution to stress transfer from current time step
    
    term = 'current'
    
    lb = 0
    ub = 0

    call conv_history(mdl,krn,cnv,term,k,kx,ky,lb,ub,cnv%H(i,j),FxC,FyC)
    
    ! update contribution to stress transfer from past history, excluding current time step
    
    if (update_history) then
       
       term = 'history'

       lb = 1
       ub = min(cnv%H(i,j)%n,mdl%n)

       if (ub==0) then
          cnv%FxH(i,j) = zeroc
          cnv%FyH(i,j) = zeroc
       else
          call conv_history(mdl,krn,cnv,term,k,kx,ky,lb,ub,cnv%H(i,j),cnv%FxH(i,j),cnv%FyH(i,j))
       end if

    end if

    ! sum contributions to stress transfer
    
    cnv%Fx(i,j) = cnv%FxH(i,j)+FxC+FxS
    cnv%Fy(i,j) = cnv%FyH(i,j)+FyC+FyS
    
  end subroutine do_single_convolution


  subroutine conv_static(mdl,krn,k,km,Dx,Dy,Fx,Fy)
    ! CONV_STATIC calculates static contribution to stress transfer,
    ! identical materials case
    ! 
    ! Modified: 23 January 2007

    use constants, only : half,zeroc
    use model, only : model_type 
    use kernel, only : kernel_type

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! KRN = kernel variables
    ! K = wavenumber amplitude
    ! KM = wavenumber product array
    ! DX = Fourier transform of slip in x direction
    ! DY = Fourier transform of slip in y direction
    ! FX = Fourier transform of stress transfer in x direction
    ! FY = Fourier transform of stress transfer in y direction

    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(in) :: krn
    real(pr),intent(in) :: k,km(2,2)
    complex(pr),intent(in) :: Dx,Dy
    complex(pr),intent(out) :: Fx,Fy
    
    ! form the appropriate products of Fourier coefficients of slip
    ! with convolution kernel and wavenumbers to get static 
    ! stress transfer functional

    select case(mdl%case)
    case(1)
       Fx = -half*mdl%mu/k* &
            ((km(1,1)*krn%C2S+km(2,2)*krn%C3S)*Dx+ &
             (km(1,2)*krn%C2S-km(2,1)*krn%C3S)*Dy)
       Fy = -half*mdl%mu/k* &
            ((km(2,1)*krn%C2S-km(1,2)*krn%C3S)*Dx+ &
             (km(2,2)*krn%C2S+km(1,1)*krn%C3S)*Dy)    
    case(2)
       Fx = -half*mdl%mu/k* &
            (km(1,1)*krn%C2S+km(2,2)*krn%C3S)*Dx
       Fy = -half*mdl%mu/k* &
            (km(2,1)*krn%C2S-km(1,2)*krn%C3S)*Dx
    case(3)
       select case(mdl%mode)
       case(2)
          Fx = -half*mdl%mu*k*krn%C2S*Dx
       case(3)
          Fx = -half*mdl%mu*k*krn%C3S*Dx
       end select
       Fy = zeroc
    end select

  end subroutine conv_static


  subroutine conv_history(mdl,krn,cnv,term,k,kx,ky,lb,ub,H,Fx,Fy)
    ! CONV_HISTORY calculates contribution to stress transfer 
    ! from past history, identical materials case
    ! 
    ! Modified: 19 January 2007

    use constants, only : zeroc,half
    use model, only : model_type 
    use kernel, only : kernel_type,calc_kernel
    use convolution, only : convolution_type
    use history, only : history_type
    use utilities, only : convolve

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! KRN = kernel variables
    ! CNV = convolution variables
    ! TERM = term in convolution to compute (history or current step)
    ! K = wavenumber amplitude
    ! KX = wavenumber  in x direction
    ! KY = wavenumber  in y direction
    ! LB = lower index/bound of history vectors Hx and Hy
    ! UB = upper index/bound of history vectors Hx and Hy
    ! H = history variables
    ! FX = Fourier transform of stress transfer in x direction
    ! FY = Fourier transform of stress transfer in y direction    

    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(in) :: krn
    type(convolution_type),intent(in) :: cnv
    character(*),intent(in) :: term
    real(pr),intent(in) :: k,kx,ky
    integer(pin),intent(in) :: lb,ub
    type(history_type),intent(in) :: H
    complex(pr),intent(out) :: Fx,Fy

    ! Internal Parameters:
    ! N = counter for time steps
    ! CF = coefficient of stress transfer
    ! T = non-dimensional time vector
    ! C2HX = dot product of C2 and Hx
    ! C2HY = dot product of C2 and Hy
    ! C3HX = dot product of C3 and Hx
    ! C3HY = dot product of C3 and Hy
    ! C2 = values of mode II convolution kernel evaluated at T
    ! C3 = values of mode III convolution kernel evaluated at T

    integer(pin) :: n
    real(pr) :: cF
    complex(pr) :: C2Hx,C2Hy,C3Hx,C3Hy
    real(pr),dimension(:),allocatable :: T,C2,C3

    ! evaluate products of kernels and history

    if (krn%precalc) then

       ! evaluate stress transfer by forming appropriate products
       
       select case(mdl%case)

       case(1)

          C2Hx = convolve(H%C(lb:ub,2,1),H%Hx(lb:ub),term)
          C3Hx = convolve(H%C(lb:ub,3,1),H%Hx(lb:ub),term)
          C2Hy = convolve(H%C(lb:ub,2,1),H%Hy(lb:ub),term)
          C3Hy = convolve(H%C(lb:ub,3,1),H%Hy(lb:ub),term)
          
       case(2)
          
          C2Hx = convolve(H%C(lb:ub,2,1),H%Hx(lb:ub),term)
          C3Hx = convolve(H%C(lb:ub,3,1),H%Hx(lb:ub),term)
          
       case(3)
          
          select case(mdl%mode)
          case(2)
             C2Hx = convolve(H%C(lb:ub,2,1),H%Hx(lb:ub),term)
          case(3)
             C3Hx = convolve(H%C(lb:ub,3,1),H%Hx(lb:ub),term)
          end select
          
       end select

    else

       ! initialize non-dimensional time vector T=k*cs*t
       
       allocate(T(lb:ub))

       T = k*mdl%cs*mdl%dt*real( (/ (n,n=lb,ub) /) ,pr)
       
       select case(mdl%case)
          
       case(1)

          allocate(C2(lb:ub),C3(lb:ub))

          C2 = calc_kernel(T,2,krn)
          C3 = calc_kernel(T,3,krn)
          
          C2Hx = convolve(C2,H%Hx(lb:ub),term)
          C3Hx = convolve(C3,H%Hx(lb:ub),term)
          C2Hy = convolve(C2,H%Hy(lb:ub),term)
          C3Hy = convolve(C3,H%Hy(lb:ub),term)
          
       case(2)
          
          allocate(C2(lb:ub),C3(lb:ub))

          C2 = calc_kernel(T,2,krn)
          C3 = calc_kernel(T,3,krn)
          
          C2Hx = convolve(C2,H%Hx(lb:ub),term)
          C3Hx = convolve(C3,H%Hx(lb:ub),term)
          
       case(3)
          
          select case(mdl%mode)
          case(2)
             allocate(C2(lb:ub))
             C2 = calc_kernel(T,2,krn)
             C2Hx = convolve(C2,H%Hx(lb:ub),term)
          case(3)
             allocate(C2(lb:ub),C3(lb:ub))
             C3 = calc_kernel(T,3,krn)
             C3Hx = convolve(C3,H%Hx(lb:ub),term)
          end select
          
       end select

       deallocate(T)
       if (allocated(C2)) deallocate(C2)
       if (allocated(C3)) deallocate(C3)

    end if

    ! evaluate stress transfer by forming appropriate products
    
    select case(mdl%case)
    case(1)
       Fx = kx**2*C2Hx+kx*ky*C2Hy+ky**2*C3Hx-kx*ky*C3Hy
       Fy = ky**2*C2Hy+kx*ky*C2Hx+kx**2*C3Hy-kx*ky*C3Hx
    case(2)
       Fx = kx**2*C2Hx+ky**2*C3Hx
       Fy = kx*ky*C2Hx-kx*ky*C3Hx
    case(3)
       select case(mdl%mode)
       case(2)
          Fx = k**2*C2Hx
       case(3)
          Fx = k**2*C3Hx
       end select
       Fy = zeroc
    end select

    ! multiply stress transfer by appropriate factor

    select case(cnv%formulation)
    case('displacement')
       cF = -half*mdl%mu*mdl%cs*mdl%dt
    case('velocity')
       cF = half*mdl%mu*mdl%dt/k
    end select

    Fx = cF*Fx
    Fy = cF*Fy

  end subroutine conv_history


end module convolution_im
