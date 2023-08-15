! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module friction_static

  ! FRICTION_STATIC contains routines to solve friction laws with static elasticity
  ! 
  ! Modified: 29 November 2006

  use constants, only : pr,pin

  implicit none


contains


  subroutine static_elasticity_im(st,sx,sy,Ux,Uy,transverse_slip,g)
    ! STATIC_ELASTICITY_IM returns misfit between elasticity equation and frictional strength.
    ! This version is for identical materials.
    ! 
    ! Modified: 29 November 2006

    implicit none

    ! I/O Parameters:
    ! ST = shear strength
    ! SX = stress in x direction
    ! SY = stress in y direction
    ! UX = slip in x direction
    ! UY = slip in y direction
    ! TRANSVERSE_SLIP = scalar or vector slip
    ! G = misfit
    
    real(pr),intent(in) :: st,Ux,Uy,sx,sy
    logical,intent(in) :: transverse_slip
    real(pr),intent(out) :: g(:)

    if (transverse_slip) then

       ! 1. shear stress equals strength
       
       g(1) = sqrt(sx**2+sy**2)-st

       ! 2. stress and slip are antiparallel

       g(2) = sx*Uy-sy*Ux

    else

       ! 1. shear stress equals strength
       
       g(1) = sx-st

       ! 2. no slip in y direction

       g(2) = Uy

    end if


  end subroutine static_elasticity_im


  subroutine static_elasticity_bm(st,sx,sy,fxp,fxm,fyp,fym,uxp,uxm,uyp,uym,transverse_slip,g)
    ! STATIC_ELASTICITY_BM returns misfit between elasticity equation and frictional strength.
    ! This version is for bimaterial problems.
    ! 
    ! Modified: 29 November 2006

    implicit none

    ! I/O Parameters:
    ! ST = shear strength
    ! SX = stress in x direction
    ! SY = stress in y direction
    ! FXP = stress transfer in x direction, plus side
    ! FXM = stress transfer in x direction, minus side
    ! FYP = stress transfer in y direction, plus side
    ! FYM = stress transfer in y direction, minus side
    ! UXP = displacement in x direction, plus side
    ! UXM = displacement in x direction, minus side
    ! UYP = displacement in y direction, plus side
    ! UYM = displacement in y direction, minus side
    ! TRANSVERSE_SLIP = scalar or vector slip
    ! G = misfit
    
    real(pr),intent(in) :: st,sx,sy,fxp,fxm,fyp,fym,uxp,uxm,uyp,uym
    logical,intent(in) :: transverse_slip
    real(pr),intent(out) :: g(:)

    if (transverse_slip) then

       ! 1. shear stress equals strength
       
       g(1) = sqrt(sx**2+sy**2)-st

       ! 2. stress and slip are antiparallel

       g(2) = sx*(uyp-uym)-sy*(uxp-uxm)

       ! 3. continuity of stress in x direction

       g(3) = fxp-fxm

       ! 4. continuity of stress in y direction

       g(4) = fyp-fym

    else

       ! 1. shear stress equals strength
       
       g(1) = sx-st

       ! 2. continuity of stress in x direction

       g(2) = fxp-fxm

       ! 3. no slip in y direction

       g(3) = uyp-uym

       ! 4. continuity of stress in y direction

       g(4) = fyp-fym

    end if


  end subroutine static_elasticity_bm


end module friction_static
