! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module fourier

  ! FOURIER contains routines to manipulate fields in Fourier domain (to store history
  ! for use in convolutions, to calculate strains, etc.)
  ! 
  ! Modified: 4 December 2007

  use constants, only : pr,pin

  implicit none

contains


  subroutine test_fft(mdl,cnv)

    use model, only : model_type
    use convolution, only : convolution_type
    use fft_routines, only : forward_fft,inverse_fft
    use constants, only : zero

    implicit none

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(in) :: cnv

    ! test subroutine forward_fft by transforming Ux to Dx

    real(pr) :: err
    real(pr),dimension(mdl%mx:mdl%px,mdl%ny) :: A
    complex(pr),dimension(cnv%mkx:cnv%pkx,cnv%nky) :: B

    integer(pin) :: i

    ! initial data

    do i = mdl%mx,mdl%px
       A(i,:) = real(i,pr)
    end do

    ! output data

    !write(6,'(a)') ''
    !do i = mdl%mx,mdl%px
    !   write(6,'(a,i6,i6,f6.0)') 'A',my_rank,i,A(i,1)
    !end do

    ! transform

    call forward_fft(mdl,cnv,A,B)

    ! output data

    !write(6,'(a)') ''
    !do i = cnv%mkx,cnv%pkx
    !   write(6,'(a,i6,i6,2f20.10)') 'B',my_rank,i,B(i,1)
    !end do

    ! inverse transform

    call inverse_fft(mdl,cnv,B,A)

    ! output data

    !write(6,'(a)') ''
    !do i = mdl%mx,mdl%px
    !   write(6,'(a,i6,i6,f6.0)') 'FFT test',my_rank,i,A(i,1)
    !end do

    err = zero
    do i = mdl%mx,mdl%px
       err = err+abs(A(i,1)-real(i,pr))
    end do

    !print *
    !print *, 'rank ',my_rank,'FFT error=',err

  end subroutine test_fft


  subroutine fft_history(mdl,flt,cnv,stage,update_history)
    ! FFT_HISTORY Fourier transforms displacement/slip or velocity
    ! for use in performing convolution, only history portion here
    ! 
    ! Modified: 21 July 2007

    use model, only : model_type
    use fault, only : fault_type
    use convolution, only : convolution_type
    use fft_routines, only : forward_fft

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! CNV = convolution variables
    ! STAGE = integration stage
    ! UPDATE_HISTORY = overwrite history or shift it

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(in) :: flt
    type(convolution_type),intent(inout) :: cnv
    integer(pin),intent(in) :: stage
    logical,intent(in) :: update_history

    ! return if history isn't stored

    if (cnv%formulation=='static'.or.cnv%formulation=='quasidynamic') return

    ! forward FFT displacement/slip U or velocity V
    
    if (mdl%bm) then
       
       select case(cnv%formulation)
       case('displacement')
          call forward_fft(mdl,cnv,flt%uxp(:,:,stage),cnv%Hxp)
          call forward_fft(mdl,cnv,flt%uxm(:,:,stage),cnv%Hxm)
          call forward_fft(mdl,cnv,flt%uyp(:,:,stage),cnv%Hyp)
          if (mdl%transverse_slip) call forward_fft(mdl,cnv,flt%uym(:,:,stage),cnv%Hym)
          call forward_fft(mdl,cnv,flt%uzp(:,:,stage),cnv%Hzp)
          if (mdl%opening) call forward_fft(mdl,cnv,flt%uzm(:,:,stage),cnv%Hzm)
       case('velocity')
          call forward_fft(mdl,cnv,flt%vxp(:,:,stage),cnv%Hxp)
          call forward_fft(mdl,cnv,flt%vxm(:,:,stage),cnv%Hxm)
          call forward_fft(mdl,cnv,flt%vyp(:,:,stage),cnv%Hyp)
          if (mdl%transverse_slip) call forward_fft(mdl,cnv,flt%vym(:,:,stage),cnv%Hym)
          call forward_fft(mdl,cnv,flt%vzp(:,:,stage),cnv%Hzp)
          if (mdl%opening) call forward_fft(mdl,cnv,flt%vzm(:,:,stage),cnv%Hzm)
       end select
       
    else
       
       select case(cnv%formulation)
       case('displacement')
          call forward_fft(mdl,cnv,flt%Ux(:,:,stage),cnv%Hx)
          if (mdl%transverse_slip) call forward_fft(mdl,cnv,flt%Uy(:,:,stage),cnv%Hy)
       case('velocity')
          call forward_fft(mdl,cnv,flt%Vx(:,:,stage),cnv%Hx)
          if (mdl%transverse_slip) call forward_fft(mdl,cnv,flt%Vy(:,:,stage),cnv%Hy)
       end select
       
    end if

    ! return if process holds no k data

    if (.not.cnv%holds_k) return

    ! store result in history vector, either shifting or overwriting with new value

    if (mdl%bm) then
       call store_history_bm(mdl,cnv,update_history)
    else
       call store_history_im(mdl,cnv,update_history)
    end if

  end subroutine fft_history


  subroutine store_history_im(mdl,cnv,update)
    ! STORE_HISTORY_IM stores slip or slip velocity history in history vector
    !
    ! Modified: 12 June 2007

    use model, only : model_type
    use convolution, only : convolution_type

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! UPDATE = whether or not to update/shift the history vector
    !      T = shift contents of history vector, adding new Fourier
    !          coefficients as first element and moving all the other
    !          elements backward in the vector
    !      F = do not shift elements, instead just overwrite the first
    !          element corresponding to the current time step

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(inout) :: cnv
    logical,intent(in) :: update

    ! Internal Parameters:
    ! I = index in kx direction
    ! J = index in ky direction

    integer(pin) :: i,j

    ! For each mode, store new Fourier coefficients in history vector.  If the Fourier coefficient for this
    ! time step is new, then shift the contents of the history vector by one element and place the new value
    ! as the first element.  This is accomplished using the intrinsic eoshift (end-off shift). Otherwise,
    ! overwrite the first element of the history vector with the new value.

    do j = 1,cnv%nky
       do i = cnv%mkx,cnv%pkx

          if (cnv%H(i,j)%n==0) cycle
          
          if (update) then

             cnv%H(i,j)%Hx = eoshift(cnv%H(i,j)%Hx,-1,cnv%Hx(i,j))
             if (mdl%transverse_slip) &
                  cnv%H(i,j)%Hy = eoshift(cnv%H(i,j)%Hy,-1,cnv%Hy(i,j))

          else

             cnv%H(i,j)%Hx(0) = cnv%Hx(i,j)
             if (mdl%transverse_slip) &
                  cnv%H(i,j)%Hy(0) = cnv%Hy(i,j)

          end if

       end do
    end do

  end subroutine store_history_im


  subroutine store_history_bm(mdl,cnv,update)
    ! STORE_HISTORY_BM stores slip or slip velocity history in history vector
    ! 
    ! Modified: 12 June 2007

    use model, only : model_type
    use convolution, only : convolution_type

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! UPDATE = whether or not to update/shift the history vector
    !      T = shift contents of history vector, adding new Fourier coefficients as first element 
    !          and moving all the other elements backward in the vector
    !      F = do not shift elements, instead just overwrite the first element corresponding to the 
    !          current time step

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(inout) :: cnv
    logical,intent(in) :: update

    ! Internal Parameters:
    ! I = index in kx direction
    ! J = index in ky direction

    integer(pin) :: i,j

    ! For each mode, store new Fourier coefficients in history vector.  If the Fourier coefficient for this
    ! time step is new, then shift the contents of the history vector by one element and place the new value
    ! as the first element.  This is accomplished using the intrinsic eoshift (end-off shift). Otherwise,
    ! overwrite the first element of the history vector with the new value.

    do j = 1,cnv%nky
       do i = cnv%mkx,cnv%pkx
          
          if (cnv%H(i,j)%n==0) cycle
          
          if (update) then

             cnv%H(i,j)%Hxp = eoshift(cnv%H(i,j)%Hxp,-1,cnv%Hxp(i,j))
             cnv%H(i,j)%Hxm = eoshift(cnv%H(i,j)%Hxm,-1,cnv%Hxm(i,j))

             if (mdl%transverse_slip) then
                cnv%H(i,j)%Hyp = eoshift(cnv%H(i,j)%Hyp,-1,cnv%Hyp(i,j))
                cnv%H(i,j)%Hym = eoshift(cnv%H(i,j)%Hym,-1,cnv%Hym(i,j))
             else
                cnv%H(i,j)%Hyp = eoshift(cnv%H(i,j)%Hyp,-1,cnv%Hyp(i,j))
             end if

             if (mdl%opening) then
                cnv%H(i,j)%Hzp = eoshift(cnv%H(i,j)%Hzp,-1,cnv%Hzp(i,j))
                cnv%H(i,j)%Hzm = eoshift(cnv%H(i,j)%Hzm,-1,cnv%Hzm(i,j))
             else
                cnv%H(i,j)%Hzp = eoshift(cnv%H(i,j)%Hzp,-1,cnv%Hzp(i,j))
             end if

          else

             cnv%H(i,j)%Hxp(0) = cnv%Hxp(i,j)
             cnv%H(i,j)%Hxm(0) = cnv%Hxm(i,j)

             if (mdl%transverse_slip) then
                cnv%H(i,j)%Hyp(0) = cnv%Hyp(i,j)
                cnv%H(i,j)%Hym(0) = cnv%Hym(i,j)
             else
                cnv%H(i,j)%Hyp(0) = cnv%Hyp(i,j)
             end if

             if (mdl%opening) then
                cnv%H(i,j)%Hzp(0) = cnv%Hzp(i,j)
                cnv%H(i,j)%Hzm(0) = cnv%Hzm(i,j)
             else
                cnv%H(i,j)%Hzp(0) = cnv%Hzp(i,j)
             end if

          end if
          
       end do
    end do
  
  end subroutine store_history_bm
        

  subroutine fft_current_U(mdl,flt,cnv,stage)
    ! FFT_CURRENT_U Fourier transforms displacement/slip
    ! 
    ! Modified: 22 September 2007

    use constants, only : zero,zeroc
    use model, only : model_type
    use fault, only : fault_type
    use convolution, only : convolution_type
    use fft_routines, only : forward_fft
    use io, only : error

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! CNV = convolution variables
    ! STAGE = integration stage

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(in) :: flt
    type(convolution_type),intent(inout) :: cnv
    integer(pin),intent(in) :: stage

    if (mdl%dim==0) then

       ! spring-block model, so subtract load point displacement

       if (mdl%bm) then
          call error('No spring-block model for bimaterials','fft_current_U')
       else
          call forward_fft(mdl,cnv,flt%Ux(:,:,stage)-flt%Vxl*mdl%t,cnv%Dx)
          call forward_fft(mdl,cnv,flt%Uy(:,:,stage)-flt%Vyl*mdl%t,cnv%Dy)
       end if

    else

       if (mdl%bm) then
          
          ! bimaterial version
          
          call forward_fft(mdl,cnv,flt%uxp(:,:,stage),cnv%Dxp)
          call forward_fft(mdl,cnv,flt%uxm(:,:,stage),cnv%Dxm)
          call forward_fft(mdl,cnv,flt%uyp(:,:,stage),cnv%Dyp)
          call forward_fft(mdl,cnv,flt%uym(:,:,stage),cnv%Dym)
          call forward_fft(mdl,cnv,flt%uzp(:,:,stage),cnv%Dzp)
          call forward_fft(mdl,cnv,flt%uzm(:,:,stage),cnv%Dzm)
          
       else
          
          ! identical materials version
          
          call forward_fft(mdl,cnv,flt%Ux(:,:,stage),cnv%Dx)
          if (mdl%transverse_slip) call forward_fft(mdl,cnv,flt%Uy(:,:,stage),cnv%Dy)
          
       end if

    end if

  end subroutine fft_current_U


  subroutine inverse_fft_F(mdl,flt,cnv,stage)
    ! INVERSE_FFT_F inverse Fourier transforms stress transfer functional
    ! 
    ! Modified: 22 September 2007

    use model, only : model_type
    use fault, only : fault_type
    use convolution, only : convolution_type
    use fft_routines, only : inverse_fft

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! FLT = fault variables
    ! CNV = convolution variables   
    ! STAGE = stage at which f is returned

    type(model_type),intent(in) :: mdl
    type(fault_type),intent(inout) :: flt
    type(convolution_type),intent(in) :: cnv
    integer(pin),intent(in) :: stage

    if (mdl%bm) then
       
       call inverse_fft(mdl,cnv,cnv%Fxp,flt%fxp(:,:,stage))
       call inverse_fft(mdl,cnv,cnv%Fxm,flt%fxm(:,:,stage))
       call inverse_fft(mdl,cnv,cnv%Fyp,flt%fyp(:,:,stage))
       call inverse_fft(mdl,cnv,cnv%Fym,flt%fym(:,:,stage))
       call inverse_fft(mdl,cnv,cnv%Fzp,flt%fzp(:,:,stage))
       call inverse_fft(mdl,cnv,cnv%Fzm,flt%fzm(:,:,stage))
       
    else
       
       call inverse_fft(mdl,cnv,cnv%Fx,flt%fx(:,:,stage))
       call inverse_fft(mdl,cnv,cnv%Fy,flt%fy(:,:,stage))
       
    end if
    
  end subroutine inverse_fft_F


  subroutine strain(mdl,cnv,exxp,exxm,exyp,exym,eyyp,eyym)
    ! STRAIN calculates fault-parallel strains using
    ! spectral methods to take spatial derivatives
    ! 
    ! Modified: 4 December 2007

    use constants, only : zeroc,img,zero,half,one,two
    use model, only : model_type
    use convolution, only : convolution_type,wavenumber
    use fft_routines, only : inverse_fft

    implicit none
    
    ! I/O Parameters:
    ! MDL = model variables
    ! CNV = convolution variables
    ! EXXP = fault-parallel strain, xx component, plus side
    ! EXXM = fault-parallel strain, xx component, minus side
    ! EXYP = fault-parallel strain, xy component, plus side
    ! EXYM = fault-parallel strain, xy component, minus side
    ! EYYP = fault-parallel strain, yy component, plus side
    ! EYYM = fault-parallel strain, yy component, minus side

    type(model_type),intent(in) :: mdl
    type(convolution_type),intent(in) :: cnv
    real(pr),dimension(mdl%mx:mdl%px,mdl%ny),intent(out) :: exxp,exxm,exyp,exym,eyyp,eyym

    ! Internal Parameters:
    ! I = index in the kx direction
    ! J = index in the ky direction
    ! K = wavenumber amplitude (k**2=kx**2+ky**2)
    ! KM = wavenumber product array 
    !      km(1,1) = kx**2
    !      km(1,2) = km(2,1) = kx*ky
    !      km(2,2) = ky**2
    ! KX = wavenumber in x direction
    ! KY = wavenumber in y direction
    ! CXXP = transform of exxp
    ! CXXM = transform of exxm
    ! CXYP = transform of exyp
    ! CXYM = transform of exym
    ! CYYP = transform of eyyp
    ! CYYM = transform of eyym

    integer(pin) :: i,j
    real(pr) :: k,km(2,2),kx,ky
    complex(pr),dimension(cnv%mkx:cnv%pkx,cnv%nky) :: Cxxp,Cxxm,Cxyp,Cxym,Cyyp,Cyym

    if (mdl%bm) then
       
       ! bimaterial version
       
       do j = 1,cnv%nky
          do i = cnv%mkx,cnv%pkx
             
             ! calculate wavenumber array
             
             call wavenumber(mdl,cnv,i,j,k,km,kx,ky)

             ! calculate Fourier transforms of strains
             
             Cxxp(i,j) = img*kx*cnv%Dxp(i,j)
             Cxxm(i,j) = img*kx*cnv%Dxm(i,j)

             if (mdl%transverse_slip) then
                Cxyp(i,j) = half*img*(kx*cnv%Dyp(i,j)+ky*cnv%Dxp(i,j))
                Cxym(i,j) = half*img*(kx*cnv%Dym(i,j)+ky*cnv%Dxm(i,j))
                Cyyp(i,j) = img*ky*cnv%Dyp(i,j)
                Cyym(i,j) = img*ky*cnv%Dym(i,j)
             else
                Cxyp(i,j) = half*img*ky*cnv%Dxp(i,j)
                Cxym(i,j) = half*img*ky*cnv%Dxm(i,j)
                Cyyp(i,j) = zeroc
                Cyym(i,j) = zeroc
             end if

          end do
       end do

    else

       ! identical materials version

       do j = 1,cnv%nky
          do i = cnv%mkx,cnv%pkx
             
             ! calculate wavenumber array
             
             call wavenumber(mdl,cnv,i,j,k,km,kx,ky)

             ! calculate Fourier transforms of strains
             
             Cxxp(i,j) =  half*img*kx*cnv%Dx(i,j)
             Cxxm(i,j) = -half*img*kx*cnv%Dx(i,j)
             Cxyp(i,j) =  half**2*img*(kx*cnv%Dy(i,j)+ky*cnv%Dx(i,j))
             Cxym(i,j) = -half**2*img*(kx*cnv%Dy(i,j)+ky*cnv%Dx(i,j))
             Cyyp(i,j) =  half*img*ky*cnv%Dy(i,j)
             Cyym(i,j) = -half*img*ky*cnv%Dy(i,j)

          end do
       end do

    end if

    ! invert Fourier transforms to obtain strains
    
    call inverse_fft(mdl,cnv,Cxxp,exxp)
    call inverse_fft(mdl,cnv,Cxxm,exxm)
    call inverse_fft(mdl,cnv,Cxyp,exyp)
    call inverse_fft(mdl,cnv,Cxym,exym)
    call inverse_fft(mdl,cnv,Cyyp,eyyp)
    call inverse_fft(mdl,cnv,Cyym,eyym)

  end subroutine strain


  subroutine set_normal_displacement(mdl,krn,flt,cnv,stage)
    ! SET_NORMAL_DISPLACEMENT sets normal displacement in bimaterial problem
    ! using no opening constraint (in Fourier domain)
    ! 
    ! Modified: 6 July 2007

    use constants, only : zeroc,img,two
    use io, only : error
    use model, only : model_type
    use kernel, only : kernel_type
    use fault, only : fault_type
    use convolution, only : convolution_type
    use fft_routines, only : inverse_fft

    implicit none

    ! I/O Parameters:
    ! MDL = model variables
    ! KRN = kernel variables
    ! FLT = fault variables
    ! CNV = convolution variables
    ! STAGE = integration stage

    type(model_type),intent(in) :: mdl
    type(kernel_type),intent(in) :: krn
    type(fault_type),intent(inout) :: flt
    type(convolution_type),intent(inout) :: cnv
    integer(pin),intent(in) :: stage

    ! return if identical materials

    if (.not.mdl%bm) return

    ! return if opening is permitted

    if (mdl%opening) return

    ! solve for fault-normal displacement in Fourier domain

    select case(mdl%mode)
    case(0)
       call error('setting normal displacement not yet implemented in mixed mode or 3D','normal_displacement')
    case(1,2)
       cnv%Dzp = -img* &
            (mdl%mup*(krn%C12S+two-mdl%etap)*cnv%Dxp-mdl%mum*(krn%C12S+two-mdl%etam)*cnv%Dxm)/ &
            ((mdl%mup+mdl%mum)*krn%C22S)
       cnv%Dzm = cnv%Dzp
    case(3)
       cnv%Dzp = zeroc
       cnv%Dzm = zeroc
    end select

    ! invert transforms

    select case(mdl%dim)
    case(1,2)
       call inverse_fft(mdl,cnv,cnv%Dzp,flt%uzp(:,:,stage))
       call inverse_fft(mdl,cnv,cnv%Dzm,flt%uzm(:,:,stage))
    case(3)
       call error('setting normal displacement not yet implemented in 3D','normal_displacement')
    end select
    
  end subroutine set_normal_displacement


end module fourier
