! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module convolution_routines

  ! CONVOLUTION_ROUTINES contains subroutines to call convolution routines
  ! for identical material or bimaterial problems
  ! 
  ! Modified: 21 July 2007

  implicit none

  
contains

  
  subroutine do_convolution(mdl,krn,cnv,update_history)
    ! DO_CONVOLUTION performs convolution
    ! 
    ! Modified: 21 July 2007
    
    use model, only : model_type
    use kernel, only : kernel_type
    use convolution, only : convolution_type
    use convolution_im, only : do_convolution_im
    use convolution_bm, only : do_convolution_bm

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

    ! return if no k data stored for this process

    if (.not.cnv%holds_k) return

    ! perform convolution

    if (mdl%bm) then
       call do_convolution_bm(mdl,krn,cnv,update_history)
    else
       call do_convolution_im(mdl,krn,cnv,update_history)
    end if

  end subroutine do_convolution


end module convolution_routines
