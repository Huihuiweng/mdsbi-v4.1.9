! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module convolution_bm

  ! CONVOLUTION_BM contains subroutines to perform the convolution over 
  ! past time of the slip or slip velocity history with the appropriate convolution 
  ! kernel.  This version is for the bimaterial problem.
  ! 
  ! Modified: 12 June 2007

  use constants, only : pr,pin

  implicit none


contains


  subroutine do_convolution_bm(mdl,krn,cnv,update_history)
    ! DO_CONVOLUTION_BM performs convolution, bimaterial case
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

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)

    do j = 1,cnv%nky
       do i = cnv%mkx,cnv%pkx

          ! perform convolution on this mode

          call do_single_convolution(i,j,mdl,krn,cnv,update_history)

       end do
    end do

    !$OMP END PARALLEL DO
    

  end subroutine do_convolution_bm


  subroutine do_single_convolution(i,j,mdl,krn,cnv,update_history)
    ! DO_CONVOLUTION performs convolution, bimaterial case
    ! 
    ! Modified: 23 January 2007

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
    ! FXSP = static contribution to stress transfer functional Fxp
    ! FXSM = static contribution to stress transfer functional Fxm
    ! FYSP = static contribution to stress transfer functional Fyp
    ! FYSM = static contribution to stress transfer functional Fym
    ! FZSP = static contribution to stress transfer functional Fzp
    ! FZSM = static contribution to stress transfer functional Fzm
    ! FXPP = Poisson contribution to stress transfer functional Fxp
    ! FXPM = Poisson contribution to stress transfer functional Fxm
    ! FYPP = Poisson contribution to stress transfer functional Fyp
    ! FYPM = Poisson contribution to stress transfer functional Fym
    ! FZPP = Poisson contribution to stress transfer functional Fzp
    ! FZPM = Poisson contribution to stress transfer functional Fzm
    ! FXCP = contribution of current time step to Fxp
    ! FXCM = contribution of current time step to Fxm
    ! FYCP = contribution of current time step to Fyp
    ! FYCM = contribution of current time step to Fym
    ! FZCP = contribution of current time step to Fzp
    ! FZCM = contribution of current time step to Fzm
    ! K = wavenumber amplitude (k**2=kx**2+ky**2)
    ! KM = wavenumber product array 
    !      km(1,1) = kx**2
    !      km(1,2) = km(2,1) = kx*ky
    !      km(2,2) = ky**2
    ! KX = wavenumber in x direction
    ! KY = wavenumber in y direction
    ! TERM = term in convolution to compute (history or current step)

    integer(pin) :: lb,ub
    complex(pr) :: &
         FxSp,FxSm,FySp,FySm,FzSp,FzSm, &
         FxPp,FxPm,FyPp,FyPm,FzPp,FzPm, &
         FxCp,FxCm,FyCp,FyCm,FzCp,FzCm
    real(pr) :: k,km(2,2),kx,ky
    character(64) :: term

    ! calculate wavenumber array
    
    call wavenumber(mdl,cnv,i,j,k,km,kx,ky)

    ! treat k=0 case separately to avoid useless work and possible division by zero
    
    if (k==zero) then

       cnv%Fxp(i,j) = zeroc
       cnv%Fxm(i,j) = zeroc
       cnv%Fyp(i,j) = zeroc
       cnv%Fym(i,j) = zeroc
       cnv%Fzp(i,j) = zeroc
       cnv%Fzm(i,j) = zeroc

       return

    end if
    
    ! calculate static contribution to stress transfer

    if (cnv%formulation=='displacement') then
       FxSp = zeroc
       FxSm = zeroc
       FySp = zeroc
       FySm = zeroc
       FzSp = zeroc
       FzSm = zeroc
    else
       call conv_static(k,km,kx,ky,krn,mdl, &
            cnv%Dxp(i,j),cnv%Dxm(i,j),cnv%Dyp(i,j),cnv%Dym(i,j),cnv%Dzp(i,j),cnv%Dzm(i,j), &
            FxSp,FxSm,FySp,FySm,FzSp,FzSm)
    end if

    ! calculate Poisson contribution to stress transfer
    
    call conv_poisson(k,kx,ky,mdl, &
         cnv%Dxp(i,j),cnv%Dxm(i,j),cnv%Dyp(i,j),cnv%Dym(i,j),cnv%Dzp(i,j),cnv%Dzm(i,j), &
         FxPp,FxPm,FyPp,FyPm,FzPp,FzPm)

    ! sum contributions to stress transfer and return if static/quasidynamic (no history contribution)

    if (cnv%formulation=='static'.or.cnv%formulation=='quasidynamic') then

       cnv%Fxp(i,j) = FxSp+FxPp
       cnv%Fxm(i,j) = FxSm+FxPm
       cnv%Fyp(i,j) = FySp+FyPp
       cnv%Fym(i,j) = FySm+FyPm
       cnv%Fzp(i,j) = FzSp+FzPp
       cnv%Fzm(i,j) = FzSm+FzPm

       return

    end if

    ! calculate contribution to stress transfer from current time step
    
    term = 'current'

    lb = 0
    ub = 0

    call conv_history(mdl,krn,cnv,term,k,kx,ky,lb,ub,cnv%H(i,j), &
         FxCp,FxCm,FyCp,FyCm,FzCp,FzCm)


    ! update contribution to stress transfer from past history, excluding current time step

    if (update_history) then
       
       term = 'history'

       
       lb = 1
       ub = min(cnv%H(i,j)%n,mdl%n)
       
       call conv_history(mdl,krn,cnv,term,k,kx,ky,lb,ub,cnv%H(i,j), &
            cnv%FxHp(i,j),cnv%FxHm(i,j), &
            cnv%FyHp(i,j),cnv%FyHm(i,j), &
            cnv%FzHp(i,j),cnv%FzHm(i,j))

    end if

    ! sum contributions to stress transfer
    
    cnv%Fxp(i,j) = cnv%FxHp(i,j)+FxCp+FxSp+FxPp
    cnv%Fxm(i,j) = cnv%FxHm(i,j)+FxCm+FxSm+FxPm
    cnv%Fyp(i,j) = cnv%FyHp(i,j)+FyCp+FySp+FyPp
    cnv%Fym(i,j) = cnv%FyHm(i,j)+FyCm+FySm+FyPm
    cnv%Fzp(i,j) = cnv%FzHp(i,j)+FzCp+FzSp+FzPp
    cnv%Fzm(i,j) = cnv%FzHm(i,j)+FzCm+FzSm+FzPm

  end subroutine do_single_convolution


  subroutine conv_static(k,km,kx,ky,krn,mdl,Dxp,Dxm,Dyp,Dym,Dzp,Dzm,Fxp,Fxm,Fyp,Fym,Fzp,Fzm)
    ! CONV_STATIC calculates static contribution to stress transfer, bimaterial case
    ! 
    ! Modified: 23 January 2007

    use constants, only : img,zeroc
    use model, only : model_type 
    use kernel, only : kernel_type

    implicit none
    
    ! I/O Parameters:
    ! KRN = kernel variables
    ! MDL = model variables
    ! K = wavenumber amplitude
    ! KM = wavenumber product array
    ! KX = wavenumber in x direction
    ! KY = wavenumber in y direction
    ! DXP = Fourier transform of displacement in x direction, plus side
    ! DXM = Fourier transform of displacement in x direction, minus side
    ! DYP = Fourier transform of displacement in y direction, plus side
    ! DYM = Fourier transform of displacement in y direction, minus side
    ! DZP = Fourier transform of displacement in z direction, plus side
    ! DZM = Fourier transform of displacement in z direction, minus side
    ! FXP = Fourier transform of stress transfer in x direction, plus side
    ! FXM = Fourier transform of stress transfer in x direction, minus side
    ! FYP = Fourier transform of stress transfer in y direction, plus side
    ! FYM = Fourier transform of stress transfer in y direction, minus side
    ! FZP = Fourier transform of stress transfer in z direction, plus side
    ! FZM = Fourier transform of stress transfer in z direction, minus side

    type(kernel_type),intent(in) :: krn
    type(model_type),intent(in) :: mdl
    real(pr),intent(in) :: k,km(2,2),kx,ky
    complex(pr),intent(in)  :: Dxp,Dxm,Dyp,Dym,Dzp,Dzm
    complex(pr),intent(out) :: Fxp,Fxm,Fyp,Fym,Fzp,Fzm
    
    ! form the appropriate products of Fourier coefficients of slip with convolution kernel 
    ! and wavenumbers to get static stress transfer functional

    select case(mdl%mode)
    case(0)
       Fxp = -mdl%mup/k* &
            ((km(1,1)*krn%C11S+km(2,2)*krn%C33S)*Dxp+ &
             (km(1,2)*krn%C11S-km(2,1)*krn%C33S)*Dyp- &
             img*k*kx*krn%C12S*Dzp)
       Fxm =  mdl%mup/k* &
            ((km(1,1)*krn%C11S+km(2,2)*krn%C33S)*Dxm+ &
             (km(1,2)*krn%C11S-km(2,1)*krn%C33S)*Dym+ &
             img*k*kx*krn%C12S*Dzm)
       Fyp = -mdl%mup/k* &
            ((km(2,1)*krn%C11S-km(1,2)*krn%C33S)*Dxp+ &
             (km(2,2)*krn%C11S+km(1,1)*krn%C33S)*Dyp- &
             img*k*ky*krn%C12S*Dzp)
       Fym =  mdl%mup/k* &
            ((km(2,1)*krn%C11S-km(1,2)*krn%C33S)*Dxm+ &
             (km(2,2)*krn%C11S+km(1,1)*krn%C33S)*Dym+ &
             img*k*ky*krn%C12S*Dzm)
       Fzp = -mdl%mup* &
            (k*krn%C22S*Dzp+img*krn%C12S*(kx*Dxp+ky*Dyp))
       Fzm =  mdl%mum* &
            (k*krn%C22S*Dzm-img*krn%C12S*(kx*Dxm+ky*Dym))
    case(1,2)
       Fxp = -mdl%mup*k*(krn%C11S*Dxp-img*krn%C12S*Dzp)
       Fxm =  mdl%mum*k*(krn%C11S*Dxm+img*krn%C12S*Dzm)
       Fyp = zeroc
       Fym = zeroc
       Fzp = -mdl%mup*k*(krn%C22S*Dzp+img*krn%C12S*Dxp)
       Fzm =  mdl%mum*k*(krn%C22S*Dzm-img*krn%C12S*Dxm)
    case(3)
       Fxp = -mdl%mup*k*krn%C33S*Dxp
       Fxm =  mdl%mum*k*krn%C33S*Dxm
       Fyp = zeroc
       Fym = zeroc
       Fzp = zeroc
       Fzm = zeroc
    end select

  end subroutine conv_static


  subroutine conv_poisson(k,kx,ky,mdl,Dxp,Dxm,Dyp,Dym,Dzp,Dzm,Fxp,Fxm,Fyp,Fym,Fzp,Fzm)
    ! CONV_POISSON calculates instantaneous Poisson contribution to stress transfer
    ! 
    ! Modified: 24 October 2006

    use constants, only : img,zeroc,two
    use model, only : model_type 

    implicit none
    
    ! I/O Parameters:
    ! KRN = kernel variables
    ! MDL = model variables
    ! K = wavenumber amplitude
    ! KX = wavenumber in x direction
    ! KY = wavenumber in y direction
    ! DXP = Fourier transform of displacement in x direction, plus side
    ! DXM = Fourier transform of displacement in x direction, minus side
    ! DYP = Fourier transform of displacement in y direction, plus side
    ! DYM = Fourier transform of displacement in y direction, minus side
    ! DZP = Fourier transform of displacement in z direction, plus side
    ! DZM = Fourier transform of displacement in z direction, minus side
    ! FXP = Fourier transform of stress transfer in x direction, plus side
    ! FXM = Fourier transform of stress transfer in x direction, minus side
    ! FYP = Fourier transform of stress transfer in y direction, plus side
    ! FYM = Fourier transform of stress transfer in y direction, minus side
    ! FZP = Fourier transform of stress transfer in z direction, plus side
    ! FZM = Fourier transform of stress transfer in z direction, minus side

    type(model_type),intent(in) :: mdl
    real(pr),intent(in) :: k,kx,ky
    complex(pr),intent(in)  :: Dxp,Dxm,Dyp,Dym,Dzp,Dzm
    complex(pr),intent(out) :: Fxp,Fxm,Fyp,Fym,Fzp,Fzm
    
    ! form the appropriate products of Fourier coefficients of slip with convolution kernel 
    ! and wavenumbers to get instantaneous Poisson contribution to stress transfer functional

    select case(mdl%mode)
    case(0)
       Fxp =  img*(two-mdl%etap)*mdl%mup*kx*Dzp
       Fxm =  img*(two-mdl%etam)*mdl%mum*kx*Dzm
       Fyp =  img*(two-mdl%etap)*mdl%mup*ky*Dzp
       Fym =  img*(two-mdl%etam)*mdl%mum*ky*Dzm
       Fzp = -img*(two-mdl%etap)*mdl%mup*(kx*Dxp+ky*Dyp)
       Fzm = -img*(two-mdl%etam)*mdl%mum*(kx*Dxm+ky*Dym)
    case(1,2)
       Fxp =  img*(two-mdl%etap)*mdl%mup*k*Dzp
       Fxm =  img*(two-mdl%etam)*mdl%mum*k*Dzm
       Fyp = zeroc
       Fym = zeroc
       Fzp = -img*(two-mdl%etap)*mdl%mup*k*Dxp
       Fzm = -img*(two-mdl%etam)*mdl%mum*k*Dxm
    case(3)
       Fxp = zeroc
       Fxm = zeroc
       Fyp = zeroc
       Fym = zeroc
       Fzp = zeroc
       Fzm = zeroc
    end select

  end subroutine conv_poisson


  subroutine conv_history(mdl,krn,cnv,term,k,kx,ky,lb,ub,H,Fxp,Fxm,Fyp,Fym,Fzp,Fzm)
    ! CONV_HISTORY calculates contribution to stress transfer from past history, bimaterial case
    ! 
    ! Modified: 22 January 2007

    use constants, only : img,zeroc
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
    ! KX = wavenumber in x direction
    ! KY = wavenumber in y direction
    ! LB = lower index/bound of history vectors Hx and Hy
    ! UB = upper index/bound of history vectors Hx and Hy
    ! H = history variables
    ! FXP = Fourier transform of stress transfer in x direction, plus side
    ! FXM = Fourier transform of stress transfer in x direction, minus side
    ! FYP = Fourier transform of stress transfer in y direction, plus side
    ! FYM = Fourier transform of stress transfer in y direction, minus side
    ! FZP = Fourier transform of stress transfer in z direction, plus side
    ! FZM = Fourier transform of stress transfer in z direction, minus side

    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(in) :: krn
    type(convolution_type),intent(in) :: cnv
    character(*),intent(in) :: term
    real(pr),intent(in) :: k,kx,ky
    integer(pin),intent(in) :: lb,ub
    type(history_type),intent(in) :: H
    complex(pr),intent(out) :: Fxp,Fxm,Fyp,Fym,Fzp,Fzm

    ! Internal Parameters:
    ! N = counter for time steps
    ! CFP = coefficient of stress transfer, plus side
    ! CFM = coefficient of stress transfer, minus side
    ! TP = non-dimensional time vector, plus side
    ! TM = non-dimensional time vector, minus side
    ! C11HXP = dot product of C11 and Hx, plus side
    ! C11HXM = dot product of C11 and Hx, minus side
    ! C..H.. = ...
    ! C11P = values of convolution kernel C11 evaluated at Tp
    ! C11M = values of convolution kernel C11 evaluated at Tm
    ! C... = ...
    ! note that C11,C12,C22,C33 are called H11,H12,H22,H33 in
    ! papers by Breitenfeld and Geubelle as well as Cochard and Rice

    integer(pin) :: n
    real(pr) :: cFp,cFm
    real(pr),dimension(:),allocatable :: Tp,Tm, &
         C11p,C11m,C12p,C12m,C22p,C22m,C33p,C33m
    complex(pr) :: &
         C11Hxp,C11Hxm,C11Hyp,C11Hym, &
         C12Hxp,C12Hxm,C12Hyp,C12Hym,C12Hzp,C12Hzm, &
         C22Hzp,C22Hzm, &
         C33Hxp,C33Hxm,C33Hyp,C33Hym

    ! evaluate products of kernels and history

    if (krn%precalc) then
       
       select case(mdl%mode)
       case(0)
          C11Hxp = convolve(H%C(lb:ub,1,1),H%Hxp(lb:ub),term)
          C12Hxp = convolve(H%C(lb:ub,2,1),H%Hxp(lb:ub),term)
          C33Hxp = convolve(H%C(lb:ub,4,1),H%Hxp(lb:ub),term)
          C11Hxm = convolve(H%C(lb:ub,1,2),H%Hxm(lb:ub),term)
          C12Hxm = convolve(H%C(lb:ub,2,2),H%Hxm(lb:ub),term)
          C33Hxm = convolve(H%C(lb:ub,4,2),H%Hxm(lb:ub),term)
          C11Hyp = convolve(H%C(lb:ub,1,1),H%Hyp(lb:ub),term)
          C12Hyp = convolve(H%C(lb:ub,2,1),H%Hyp(lb:ub),term)
          C33Hyp = convolve(H%C(lb:ub,4,1),H%Hyp(lb:ub),term)
          C11Hym = convolve(H%C(lb:ub,1,2),H%Hym(lb:ub),term)
          C12Hym = convolve(H%C(lb:ub,2,2),H%Hym(lb:ub),term)
          C33Hym = convolve(H%C(lb:ub,4,2),H%Hym(lb:ub),term)
          C12Hzp = convolve(H%C(lb:ub,2,1),H%Hzp(lb:ub),term)
          C22Hzp = convolve(H%C(lb:ub,3,1),H%Hzp(lb:ub),term)
          C12Hzm = convolve(H%C(lb:ub,2,2),H%Hzm(lb:ub),term)
          C22Hzm = convolve(H%C(lb:ub,3,2),H%Hzm(lb:ub),term)
       case(1,2)
          C11Hxp = convolve(H%C(lb:ub,1,1),H%Hxp(lb:ub),term)
          C12Hxp = convolve(H%C(lb:ub,2,1),H%Hxp(lb:ub),term)
          C11Hxm = convolve(H%C(lb:ub,1,2),H%Hxm(lb:ub),term)
          C12Hxm = convolve(H%C(lb:ub,2,2),H%Hxm(lb:ub),term)
          C12Hzp = convolve(H%C(lb:ub,2,1),H%Hzp(lb:ub),term)
          C22Hzp = convolve(H%C(lb:ub,3,1),H%Hzp(lb:ub),term)
          C12Hzm = convolve(H%C(lb:ub,2,2),H%Hzm(lb:ub),term)
          C22Hzm = convolve(H%C(lb:ub,3,2),H%Hzm(lb:ub),term)
       case(3)
          C33Hxp = convolve(H%C(lb:ub,4,1),H%Hxp(lb:ub),term)
          C33Hxm = convolve(H%C(lb:ub,4,2),H%Hxm(lb:ub),term)
       end select

    else

       ! initialize non-dimensional time vector T=k*cs*t
       
       allocate(Tp(lb:ub),Tm(lb:ub))

       Tp = k*mdl%csp*mdl%dt*real( (/ (n,n=lb,ub) /) ,pr)
       Tm = k*mdl%csm*mdl%dt*real( (/ (n,n=lb,ub) /) ,pr)
       
       select case(mdl%mode)

       case(0)

          allocate(C11p(lb:ub),C11m(lb:ub))
          allocate(C12p(lb:ub),C12m(lb:ub))
          allocate(C22p(lb:ub),C22m(lb:ub))
          allocate(C33p(lb:ub),C33m(lb:ub))

          C11p = calc_kernel(Tp,11,krn)
          C11m = calc_kernel(Tm,11,krn)
          C12p = calc_kernel(Tp,12,krn)
          C12m = calc_kernel(Tm,12,krn)
          C22p = calc_kernel(Tp,22,krn)
          C22m = calc_kernel(Tm,22,krn)
          C33p = calc_kernel(Tp,33,krn)
          C33m = calc_kernel(Tm,33,krn)

          C11Hxp = convolve(C11p,H%Hxp(lb:ub),term)
          C12Hxp = convolve(C12p,H%Hxp(lb:ub),term)
          C33Hxp = convolve(C33p,H%Hxp(lb:ub),term)
          C11Hxm = convolve(C11m,H%Hxm(lb:ub),term)
          C12Hxm = convolve(C12m,H%Hxm(lb:ub),term)
          C33Hxm = convolve(C33m,H%Hxm(lb:ub),term)
          C11Hyp = convolve(C11p,H%Hyp(lb:ub),term)
          C12Hyp = convolve(C12p,H%Hyp(lb:ub),term)
          C33Hyp = convolve(C33p,H%Hyp(lb:ub),term)
          C11Hym = convolve(C11m,H%Hym(lb:ub),term)
          C12Hym = convolve(C12m,H%Hym(lb:ub),term)
          C33Hym = convolve(C33m,H%Hym(lb:ub),term)
          C12Hzp = convolve(C12p,H%Hzp(lb:ub),term)
          C22Hzp = convolve(C22p,H%Hzp(lb:ub),term)
          C12Hzm = convolve(C12m,H%Hzm(lb:ub),term)
          C22Hzm = convolve(C22m,H%Hzm(lb:ub),term)

       case(1,2)


          allocate(C11p(lb:ub),C11m(lb:ub))
          allocate(C12p(lb:ub),C12m(lb:ub))
          allocate(C22p(lb:ub),C22m(lb:ub))

          C11p = calc_kernel(Tp,11,krn)
          C11m = calc_kernel(Tm,11,krn)
          C12p = calc_kernel(Tp,12,krn)
          C12m = calc_kernel(Tm,12,krn)
          C22p = calc_kernel(Tp,22,krn)
          C22m = calc_kernel(Tm,22,krn)
          
          C11Hxp = convolve(C11p,H%Hxp(lb:ub),term)
          C12Hxp = convolve(C12p,H%Hxp(lb:ub),term)
          C11Hxm = convolve(C11m,H%Hxm(lb:ub),term)
          C12Hxm = convolve(C12m,H%Hxm(lb:ub),term)
          C12Hzp = convolve(C12p,H%Hzp(lb:ub),term)
          C22Hzp = convolve(C22p,H%Hzp(lb:ub),term)
          C12Hzm = convolve(C12m,H%Hzm(lb:ub),term)
          C22Hzm = convolve(C22m,H%Hzm(lb:ub),term)

       case(3)

          allocate(C33p(lb:ub),C33m(lb:ub))

          C33p = calc_kernel(Tp,33,krn)
          C33m = calc_kernel(Tm,33,krn)
             
          C33Hxp = convolve(C33p,H%Hxp(lb:ub),term)
          C33Hxm = convolve(C33m,H%Hxm(lb:ub),term)

       end select

    end if

    ! evaluate stress transfer by forming appropriate products

    select case(mdl%mode)
    case(0)
       Fxp = kx**2*C11Hxp+kx*ky*C11Hyp+ky**2*C33Hxp-kx*ky*C33Hyp-img*k*kx*C12Hzp
       Fxm = kx**2*C11Hxm+kx*ky*C11Hym+ky**2*C33Hxm-kx*ky*C33Hym+img*k*kx*C12Hzm
       Fyp = ky**2*C11Hyp+kx*ky*C11Hxp+kx**2*C33Hyp-kx*ky*C33Hxp-img*k*ky*C12Hzp
       Fym = ky**2*C11Hym+kx*ky*C11Hxm+kx**2*C33Hym-kx*ky*C33Hxm+img*k*ky*C12Hzm
       Fzp = k**2*C22Hzp+img*k*(kx*C12Hxp+ky*C12Hyp)
       Fzm = k**2*C22Hzm-img*k*(kx*C12Hxm+ky*C12Hym)
    case(1,2)
       Fxp = k**2*(C11Hxp-img*C12Hzp)
       Fxm = k**2*(C11Hxm+img*C12Hzm)
       Fyp = zeroc
       Fym = zeroc
       Fzp = k**2*(C22Hzp+img*C12Hxp)
       Fzm = k**2*(C22Hzm-img*C12Hxm)
    case(3)
       Fxp = k**2*C33Hxp
       Fxm = k**2*C33Hxm
       Fyp = zeroc
       Fym = zeroc
       Fzp = zeroc
       Fzm = zeroc
    end select

    ! multiply stress transfer by appropriate factor

    select case(cnv%formulation)
    case('displacement')
       cFp = -mdl%mup*mdl%csp*mdl%dt
       cFm =  mdl%mum*mdl%csm*mdl%dt
    case('velocity')
       cFp =  mdl%mup*mdl%dt/k
       cFm = -mdl%mum*mdl%dt/k
    end select

    Fxp = cFp*Fxp
    Fxm = cFm*Fxm
    Fyp = cFp*Fyp
    Fym = cFm*Fym
    Fzp = cFp*Fzp
    Fzm = cFm*Fzm

  end subroutine conv_history


end module convolution_bm
