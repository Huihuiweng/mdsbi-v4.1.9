! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module runge_kutta

  ! RUNGE_KUTTA contains variables and routines related to Runge-Kutta time stepping
  !
  ! Schemes:
  ! 2(1)2           Huen's method + Euler
  ! RK2(1)3S        Dormand (1996), stable
  ! RK3(2)3         Dormand (1996)
  ! RK3(2)4BS       Bogacki and Shampine (1989)
  ! RK3(2)4[2R+]C   Kennedy et al. (2000)
  ! RK4(3)5M        Dormand (1996)
  ! RK4(3)5[2R+]C   Kennedy et al. (2000)
  ! RK4(3)5[4R+]M   Kennedy et al. (2000)
  ! RK5(4)6         Cash and Karp (1990)
  ! RK5(4)7S        Dormand and Prince (1980)
  ! RK5(4)7M        (DOPRI5) Dormand and Prince (1980)
  ! ARK3(2)4        Kennedy and Carpenter (2003)
  ! ARK4(3)6        Kennedy and Carpenter (2003)
  ! LIRK3           Calvo et al. (2001)
  ! RK3TVD          Shu (TVD method used in WENO)
  !
  ! Modified: 24 October 2008

  use constants, only : pr,pin

  implicit none

  ! Runge-Kutta coefficients
  ! SCHEME = name of scheme
  ! Q = order of accuracy of updating method
  ! P = order of accuracy of embedded method
  ! S = number of stages
  ! A = coefficient matrix a(i,j)
  ! AI = coefficient matrix a(i,j) (implicit method)
  ! B = coefficient array for updating method (order q)
  ! BH = coefficient array for embedded method (order p)
  ! C = coefficient array for time steps

  type rk_type
     character(64) :: scheme
     integer(pin) :: q,p,s
     real(pr),dimension(:,:),allocatable :: a,ai
     real(pr),dimension(:),allocatable :: b,bh,c
  end type rk_type

contains


  subroutine init_rk(scheme,rk)
    ! INIT_RK initializes Runge-Kutta variables.  Schemes come from:
    !
    ! Dormand, J. R. and P. J. Prince (1980) A family of embedded Runge-Kutta formulae,
    ! Journal of Computational and Applied Mathematics, 6(1), 19-26.
    !
    ! Bogaki, P. and L. F. Shampine (1989) A 3(2) pair of Runge-Kutta formulas, Applied
    ! Mathematics Letters, 2(4), 321-325
    !
    ! Cash, J. R. and A. H. Karp (1990) A variable order Runge-Kutta method for initial value
    ! problems with rapidly varying right-hand sides, ACM Transactions on Mathematical Software.
    ! 16(3), 201-222.
    !
    ! Dormand, J. R. (1996) Numerical Methods for Differential Equations: A Computational Approach. 
    ! Boca Raton, FL: CRC Press, 1996
    !
    ! Kennedy et al. (2000) Low-storage, explicit Runge-Kutta schemes for the compressible
    ! Navier-Stokes equations, Applied Numerical Mathematics, 35, 177-219.
    !
    ! Kennedy, C. A. and M. H. Carpenter (2003) Additive Runge-Kutta schemes for 
    ! convection-diffusion-reaction equations, Applied Numerical Mathematics, 44,  139-181.
    !
    ! Calvo, M. P. et al. (2001) Linearly implicit Runge-Kutta methods for advection-reaction-diffusion 
    ! equations, Applied Numerical Mathematics, 32, 535-549.
    ! 
    ! Modified: 23 October 2008

    use constants, only : zero
    use io, only : error,warning

    implicit none

    ! I/O Parameters:
    ! SCHEME = Runge-Kutta scheme
    ! RK = Runge-Kutta variables

    character(*),intent(in) :: scheme
    type(rk_type),intent(out) :: rk

    ! Internal Parameters:
    ! M = stage

    integer(pin) :: m

    rk%scheme = scheme

    select case(scheme)

    case default

       call error('invalid Runge-Kutta scheme','init_rk')

    case('RK2(1)2') ! Huen's method + Euler

       rk%q = 2
       rk%p = 1
       rk%s = 2

       allocate(rk%b(2))
       rk%b(1) = 1._pr/2._pr
       rk%b(2) = 1._pr/2._pr

       allocate(rk%a(2,1))
       rk%a = zero
       rk%a(2,1) = 1._pr

       allocate(rk%bh(2))
       rk%bh(1) = 1._pr
       rk%bh(2) = 0._pr

    case('RK2(1)3S') ! Dormand (1996), p. 131

       rk%q = 2
       rk%p = 1
       rk%s = 3

       allocate(rk%b(3))
       rk%b(1) = 0._pr
       rk%b(2) = 0._pr
       rk%b(3) = 1._pr

       allocate(rk%a(3,2))
       rk%a = zero
       rk%a(2,1) = 1._pr/8._pr
       rk%a(3,2) = 1._pr/2._pr

       allocate(rk%bh(3))
       rk%bh(1) = 19._pr/243._pr
       rk%bh(2) = 608._pr/729._pr
       rk%bh(3) = 64._pr/729._pr

    case('RK3(2)3') ! Dormand (1996), p. 79

       rk%q = 3
       rk%p = 2
       rk%s = 3

       allocate(rk%b(3))
       rk%b(1) = 1._pr/6._pr
       rk%b(2) = 2._pr/3._pr
       rk%b(3) = 1._pr/6._pr

       allocate(rk%a(3,2))
       rk%a = zero
       rk%a(2,1) = 1._pr/2._pr
       rk%a(3,1) = -1._pr
       rk%a(3,2) = 2._pr

       allocate(rk%bh(3))
       rk%bh(1) = 1._pr/2._pr
       rk%bh(2) = zero
       rk%bh(3) = 1._pr/2._pr

    case('RK3(2)4BS') ! Bogacki and Shampine (1989)

       rk%q = 3
       rk%p = 2
       rk%s = 4

       allocate(rk%b(4))
       rk%b(1) = 2._pr/9._pr
       rk%b(2) = 1._pr/3._pr
       rk%b(3) = 4._pr/9._pr
       rk%b(4) = 0._pr

       allocate(rk%a(4,3))
       rk%a = zero
       rk%a(2,1) = 1._pr/2._pr
       rk%a(3,1) = 0._pr
       rk%a(3,2) = 3._pr/4._pr
       rk%a(4,1) = rk%b(1)
       rk%a(4,2) = rk%b(2)
       rk%a(4,3) = rk%b(3)

       allocate(rk%bh(4))
       rk%bh(1) = 7._pr/24._pr
       rk%bh(2) = 1._pr/4._pr
       rk%bh(3) = 1._pr/3._pr
       rk%bh(4) = 1._pr/8._pr

    case('RK3(2)4[2R+]C') ! Kennedy et al. (2000)

       rk%q = 3
       rk%p = 2
       rk%s = 4

       allocate(rk%b(4))
       rk%b(1) = 1017324711453._pr/9774461848756._pr
       rk%b(2) = 8237718856693._pr/13685301971492._pr
       rk%b(3) = 57731312506979._pr/19404895981398._pr
       rk%b(4) = -101169746363290._pr/37734290219643._pr

       allocate(rk%a(4,3))
       rk%a = zero
       rk%a(2,1) = 11847461282814._pr/36547543011857._pr
       rk%a(3,1) = rk%b(1)
       rk%a(3,2) = 3943225443063._pr/7078155732230._pr
       rk%a(4,1) = rk%b(1)
       rk%a(4,2) = rk%b(2)
       rk%a(4,3) = -346793006927._pr/4029903576067._pr

       allocate(rk%bh(4))
       rk%bh(1) = 15763415370699._pr/46270243929542._pr
       rk%bh(2) = 514528521746._pr/5659431552419._pr
       rk%bh(3) = 27030193851939._pr/9429696342944._pr
       rk%bh(4) = -69544964788955._pr/30262026368149._pr

    case('RK4(3)5M') ! Dormand (1996) p. 83

       rk%q = 4
       rk%p = 3
       rk%s = 5

       allocate(rk%b(5))
       rk%b(1) = 13._pr/96._pr
       rk%b(2) = 0._pr
       rk%b(3) = 25._pr/48._pr
       rk%b(4) = 25._pr/96._pr
       rk%b(5) = 1._pr/12._pr

       allocate(rk%a(5,4))
       rk%a = zero
       rk%a(2,1) = 1._pr/5._pr
       rk%a(3,1) = 0._pr
       rk%a(3,2) = 2._pr/5._pr
       rk%a(4,1) = 6._pr/5._pr
       rk%a(4,2) = -12._pr/5._pr
       rk%a(4,3) = 2._pr
       rk%a(5,1) = -17._pr/8._pr
       rk%a(5,2) = 5._pr
       rk%a(5,3) = -5._pr/2._pr
       rk%a(5,4) = 5._pr/8._pr

       allocate(rk%bh(5))
       rk%bh(1) = 23._pr/192._pr
       rk%bh(2) = 0._pr
       rk%bh(3) = 55._pr/96._pr
       rk%bh(4) = 35._pr/192._pr
       rk%bh(5) = 1._pr/8._pr

    case('RK4(3)5[2R+]C') ! Kennedy et al. (2000)

       rk%q = 4
       rk%p = 3
       rk%s = 5

       allocate(rk%b(5))
       rk%b(1) = 1153189308089._pr/22510343858157._pr
       rk%b(2) = 1772645290293._pr/4653164025191._pr
       rk%b(3) = -1672844663538._pr/4480602732383._pr
       rk%b(4) = 2114624349019._pr/3568978502595._pr
       rk%b(5) = 5198255086312._pr/14908931495163._pr

       allocate(rk%a(5,4))
       rk%a = zero
       rk%a(2,1) = 970286171893._pr/4311952581923._pr
       rk%a(3,1) = rk%b(1)
       rk%a(3,2) = 6584761158862._pr/12103376702013._pr
       rk%a(4,1) = rk%b(1)
       rk%a(4,2) = rk%b(2)
       rk%a(4,3) = 2251764453980._pr/15575788980749._pr
       rk%a(5,1) = rk%b(1)
       rk%a(5,2) = rk%b(2)
       rk%a(5,3) = rk%b(3)
       rk%a(5,4) = 26877169314380._pr/34165994151039._pr

       allocate(rk%bh(5))
       rk%bh(1) = 1016888040809._pr/7410784769900._pr
       rk%bh(2) = 11231460423587._pr/58533540763752._pr
       rk%bh(3) = -1563879915014._pr/6823010717585._pr
       rk%bh(4) = 606302364029._pr/971179775848._pr
       rk%bh(5) = 1097981568119._pr/3980877426909._pr

    case('RK4(3)5[4R+]M') ! Kennedy et al. (2000)

       rk%q = 4
       rk%p = 3
       rk%s = 5

       allocate(rk%b(5))
       rk%b(1) = 297809._pr/2384418._pr
       rk%b(2) = 0._pr
       rk%b(3) = 156250000._pr/270591503._pr
       rk%b(4) = 5030000._pr/888933._pr
       rk%b(5) = -2927._pr/546._pr

       allocate(rk%a(5,4))
       rk%a = zero
       rk%a(2,1) = 7142524119._pr/20567653057._pr
       rk%a(3,1) = 15198616943._pr/89550000000._pr
       rk%a(3,2) = 20567653057._pr/89550000000._pr
       rk%a(4,1) = 9890667227._pr/80359280000._pr
       rk%a(4,2) = -226244183627._pr/80359280000._pr
       rk%a(4,3) = 7407775._pr/2008982._pr
       rk%a(5,1) = rk%b(1)
       rk%a(5,2) = -20567653057._pr/6979191486._pr
       rk%a(5,3) = 33311687500._pr/8703531091._pr
       rk%a(5,4) = -4577300._pr/867302297._pr

       allocate(rk%bh(5))
       rk%bh(1) = 121286694859._pr/931793198518._pr
       rk%bh(2) = 0._pr
       rk%bh(3) = 9680751416357._pr/17201392077364._pr
       rk%bh(4) = 6633076090000._pr/1042143269349._pr
       rk%bh(5) = -127961558623._pr/21123456354._pr

    case('RK5(4)6') ! Cash and Karp (1990)

       rk%q = 5
       rk%p = 4
       rk%s = 6

       allocate(rk%b(6))
       rk%b(1) = 37._pr/378._pr
       rk%b(2) = 0._pr
       rk%b(3) = 250._pr/621._pr
       rk%b(4) = 125._pr/594._pr
       rk%b(5) = 0._pr
       rk%b(6) = 512._pr/1771._pr

       allocate(rk%a(6,5))
       rk%a = zero
       rk%a(2,1) = 1._pr/5._pr
       rk%a(3,1) = 3._pr/40._pr
       rk%a(3,2) = 9._pr/40._pr
       rk%a(4,1) = 3._pr/10._pr
       rk%a(4,2) = -9._pr/10._pr
       rk%a(4,3) = 6._pr/5._pr
       rk%a(5,1) = -11._pr/54._pr
       rk%a(5,2) = 5._pr/2._pr
       rk%a(5,3) = -70._pr/27._pr
       rk%a(5,4) = 35._pr/27._pr
       rk%a(6,1) = 1631._pr/55296._pr
       rk%a(6,2) = 175._pr/512._pr
       rk%a(6,3) = 575._pr/13824._pr
       rk%a(6,4) = 44275._pr/110592._pr
       rk%a(6,5) = 253._pr/4096._pr

       allocate(rk%bh(6))
       rk%bh(1) = 2825._pr/27648._pr
       rk%bh(2) = 0._pr
       rk%bh(3) = 18575._pr/48384._pr
       rk%bh(4) = 13525._pr/55296._pr
       rk%bh(5) = 277._pr/14336._pr
       rk%bh(6) = 1._pr/4._pr

    case('RK5(4)7S') ! Dormand and Prince (1980)

       rk%q = 5
       rk%p = 4
       rk%s = 7

       allocate(rk%b(7))
       rk%b(1) = 19._pr/200._pr
       rk%b(2) = 0._pr
       rk%b(3) = 3._pr/5._pr
       rk%b(4) = -243._pr/400._pr
       rk%b(5) = 33._pr/40._pr
       rk%b(6) = 7._pr/80._pr
       rk%b(7) = 0._pr

       allocate(rk%a(7,6))
       rk%a = zero
       rk%a(2,1) = 2._pr/9._pr
       rk%a(3,1) = 1._pr/12._pr
       rk%a(3,2) = 1._pr/4._pr
       rk%a(4,1) = 55._pr/324._pr
       rk%a(4,2) = -25._pr/108._pr
       rk%a(4,3) = 50._pr/81._pr
       rk%a(5,1) = 83._pr/330._pr
       rk%a(5,2) = -13._pr/22._pr
       rk%a(5,3) = 61._pr/66._pr
       rk%a(5,4) = 9._pr/110._pr
       rk%a(6,1) = -19._pr/28._pr
       rk%a(6,2) = 9._pr/4._pr
       rk%a(6,3) = 1._pr/7._pr
       rk%a(6,4) = -27._pr/7._pr
       rk%a(6,5) = 22._pr/7._pr
       rk%a(7,1) = 19._pr/200._pr
       rk%a(7,2) = 0._pr
       rk%a(7,3) = 3._pr/5._pr
       rk%a(7,4) = -243._pr/400._pr
       rk%a(7,5) = 33._pr/40._pr
       rk%a(7,6) = 7._pr/80._pr

       allocate(rk%bh(7))
       rk%bh(1) = 431._pr/5000._pr
       rk%bh(2) = 0._pr
       rk%bh(3) = 333._pr/500._pr
       rk%bh(4) = -7857._pr/10000._pr
       rk%bh(5) = 957._pr/1000._pr
       rk%bh(6) = 193._pr/2000._pr
       rk%bh(7) = -1._pr/50._pr

    case('RK5(4)7M') ! (DOPRI5) Dormand and Prince (1980)

       call warning('preliminary tests of this RK scheme show errors','init_rk')

       rk%q = 5
       rk%p = 4
       rk%s = 7

       allocate(rk%b(7))
       rk%b(1) = 35._pr/384._pr
       rk%b(2) = 0._pr
       rk%b(3) = 500._pr/1113._pr
       rk%b(4) = 125._pr/192._pr
       rk%b(5) = -2187._pr/6784._pr
       rk%b(6) = 11._pr/84._pr
       rk%b(7) = 0._pr

       allocate(rk%a(7,6))
       rk%a = zero
       rk%a(2,1) = 0.2_pr
       rk%a(3,1) = 3._pr/40._pr
       rk%a(3,2) = 9._pr/40._pr
       rk%a(4,1) = 44._pr/45._pr
       rk%a(4,2) = -56._pr/15._pr
       rk%a(4,3) = 32._pr/9._pr
       rk%a(5,1) = 19372._pr/6561._pr
       rk%a(5,2) = -25360._pr/2187._pr
       rk%a(5,3) = 64448._pr/6561._pr
       rk%a(5,4) = -212._pr/729._pr
       rk%a(6,1) = 9017._pr/3168._pr
       rk%a(6,2) = -355._pr/33._pr
       rk%a(6,3) = 46732._pr/5247._pr
       rk%a(6,4) = 49._pr/176._pr
       rk%a(6,5) = -5103._pr/18656._pr
       rk%a(7,1) = 35._pr/384._pr
       rk%a(7,2) = 0._pr
       rk%a(7,3) = 500._pr/1113._pr
       rk%a(7,4) = 125._pr/192._pr
       rk%a(7,5) = -2187._pr/6784._pr
       rk%a(7,6) = 11._pr/84._pr

       allocate(rk%bh(7))
       rk%bh(1) = 5197._pr/57600._pr
       rk%bh(2) = 0._pr
       rk%bh(3) = 7571._pr/16695._pr
       rk%bh(4) = 393._pr/640._pr
       rk%bh(5) = -92097._pr/339200._pr
       rk%bh(6) = 187._pr/2100._pr
       rk%bh(7) = 1._pr/40._pr

    case('ARK3(2)4') ! Kennedy and Carpenter (2003)
       
       rk%q = 3
       rk%p = 2
       rk%s = 4

       allocate(rk%b(4))
       rk%b(1) = 1471266399579._pr/7840856788654._pr
       rk%b(2) = -4482444167858._pr/7529755066697._pr
       rk%b(3) = 11266239266428._pr/11593286722821._pr
       rk%b(4) = 1767732205903._pr/4055673282236._pr

       allocate(rk%bh(4))
       rk%bh(1) = 2756255671327._pr/12835298489170._pr
       rk%bh(2) = -10771552573575._pr/22201958757719._pr
       rk%bh(3) = 9247589265047._pr/10645013368117._pr
       rk%bh(4) = 2193209047091._pr/5459859503100._pr

       allocate(rk%a(4,4))
       rk%a = zero
       rk%a(2,1) = 1767732205903._pr/2027836641118._pr
       rk%a(3,1) = 5535828885825._pr/10492691773637._pr
       rk%a(3,2) = 788022342437._pr/10882634858940._pr
       rk%a(4,1) = 6485989280629._pr/16251701735622._pr
       rk%a(4,2) = -4246266847089._pr/9704473918619._pr
       rk%a(4,3) = 10755448449292._pr/10357097424841._pr

       allocate(rk%ai(4,4))
       rk%ai = zero
       rk%ai(2,1) = 1767732205903._pr/4055673282236._pr
       rk%ai(2,2) = 1767732205903._pr/4055673282236._pr
       rk%ai(3,1) = 2746238789719._pr/10658868560708._pr
       rk%ai(3,2) = -640167445237._pr/6845629431997._pr
       rk%ai(3,3) = 1767732205903._pr/4055673282236._pr
       rk%ai(4,1) = 1471266399579._pr/7840856788654._pr
       rk%ai(4,2) = -4482444167858._pr/7529755066697._pr
       rk%ai(4,3) = 11266239266428._pr/11593286722821._pr
       rk%ai(4,4) = 1767732205903._pr/4055673282236._pr

    case('ARK4(3)6') ! Kennedy and Carpenter (2003)
       
       rk%q = 4
       rk%p = 3
       rk%s = 6

       allocate(rk%b(6))
       rk%b(1) = 82889._pr/524892._pr
       rk%b(2) = 0._pr
       rk%b(3) = 15625._pr/83664._pr
       rk%b(4) = 69875._pr/102672._pr
       rk%b(5) = -2260._pr/8211._pr
       rk%b(6) = 1._pr/4._pr

       allocate(rk%bh(6))
       rk%bh(1) = 4586570599._pr/29645900160._pr
       rk%bh(2) = 0._pr
       rk%bh(3) = 178811875._pr/945068544._pr
       rk%bh(4) = 814220225._pr/1159782912._pr
       rk%bh(5) = -3700637._pr/11593932._pr
       rk%bh(6) = 61727._pr/225920._pr

       allocate(rk%a(6,6))
       rk%a = zero
       rk%a(2,1) = 1._pr/2._pr
       rk%a(3,1) = 13861._pr/62500._pr
       rk%a(3,2) = 6889._pr/62500._pr
       rk%a(4,1) = -116923316275._pr/2393684061468._pr
       rk%a(4,2) = -2731218467317._pr/15368042101831._pr
       rk%a(4,3) = 9408046702089._pr/11113171139209._pr
       rk%a(5,1) = -451086348788._pr/2902428689909._pr
       rk%a(5,2) = -2682348792572._pr/7519795681897._pr
       rk%a(5,3) = 12662868775082._pr/11960479115383._pr
       rk%a(5,4) = 3355817975965._pr/11060851509271._pr
       rk%a(6,1) = 647845179188._pr/3216320057751._pr
       rk%a(6,2) = 73281519250._pr/8382639484533._pr
       rk%a(6,3) = 552539513391._pr/3454668386233._pr
       rk%a(6,4) = 3354512671639._pr/8306763924573._pr
       rk%a(6,5) = 4040._pr/17871._pr

       allocate(rk%ai(6,6))
       rk%ai = zero
       rk%ai(2,1) = 1._pr/4._pr
       rk%ai(2,2) = 1._pr/4._pr
       rk%ai(3,1) = 8611._pr/62500._pr
       rk%ai(3,2) = -1743._pr/31250._pr
       rk%ai(3,3) = 1._pr/4._pr
       rk%ai(4,1) = 5012029._pr/34652500._pr
       rk%ai(4,2) = -654441._pr/2922500._pr
       rk%ai(4,3) = 174375._pr/388108._pr
       rk%ai(4,4) = 1._pr/4._pr
       rk%ai(5,1) = 15267082809._pr/155376265600._pr
       rk%ai(5,2) = -71443401._pr/120774400._pr
       rk%ai(5,3) = 730878875._pr/902184768._pr
       rk%ai(5,4) = 2285395._pr/8070912._pr
       rk%ai(5,5) = 1._pr/4._pr
       rk%ai(6,1) = 82889._pr/524892._pr
       rk%ai(6,2) = 0._pr
       rk%ai(6,3) = 15625._pr/83664._pr
       rk%ai(6,4) = 69875._pr/102672._pr
       rk%ai(6,5) = -2260._pr/8211._pr
       rk%ai(6,6) = 1._pr/4._pr

    case('LIRK3') ! Calvo et al. (2001)

       rk%q = 3
       rk%p = 2
       rk%s = 5

       allocate(rk%b(5))
       rk%b(1) = zero
       rk%b(2) = 1.208496649176010070_pr
       rk%b(3) = -0.644363170684469070_pr
       rk%b(4) = 0.4358665215084589994_pr
       rk%b(5) = zero

       allocate(rk%bh(5))
       rk%bh = rk%b
       print *, 'no error estimate for LIRK3 yet'

       allocate(rk%a(5,5))
       rk%a = zero
       rk%a(2,1) = 0.4358665215084589994_pr
       rk%a(3,1) = 1.0679332607542294997_pr
       rk%a(3,2) = -0.35_pr
       rk%a(4,2) = 1.98917572467984610_pr
       rk%a(4,3) = -0.98917572467984610_pr
       rk%a(5,2) = 1.208496649176010070_pr
       rk%a(5,3) = -0.644363170684469070_pr
       rk%a(5,4) = 0.4358665215084589994_pr

       allocate(rk%ai(5,5))
       rk%ai = zero
       rk%ai(2,2) = 0.4358665215084589994_pr
       rk%ai(3,2) = 0.2820667392457705003_pr
       rk%ai(3,3) = 0.4358665215084589994_pr
       rk%ai(4,2) = 1.208496649176010070_pr
       rk%ai(4,3) = -0.644363170684469070_pr
       rk%ai(4,4) = 0.4358665215084589994_pr
       rk%ai(5,2) = 1.208496649176010070_pr
       rk%ai(5,3) = -0.644363170684469070_pr
       rk%ai(5,5) = 0.4358665215084589994_pr
       
    case('Ascher(1,2,2)') ! implicit-explicit midpoint

       rk%q = 2
       rk%p = 1
       rk%s = 2

       allocate(rk%b(2))
       rk%b(1) = 0._pr
       rk%b(2) = 1._pr

       allocate(rk%a(2,2))
       rk%a = zero
       rk%a(2,1) = 0.5_pr

       allocate(rk%ai(2,2))
       rk%ai = zero
       rk%ai(2,2) = 0.5_pr

       allocate(rk%bh(2))
       rk%bh = rk%b

    case('Ascher(2,3,2)')

       rk%q = 2
       rk%p = 1
       rk%s = 3

       allocate(rk%b(3))
       rk%b(1) = 0._pr
       rk%b(2) = 0.29289321881345_pr
       rk%b(3) = 1._pr

       allocate(rk%a(3,3))
       rk%a = zero
       rk%a(2,1) = 0.29289321881345_pr
       rk%a(3,1) = -0.94280904158206_pr
       rk%a(3,2) = 1.94280904158206_pr

       allocate(rk%ai(3,3))
       rk%ai = zero
       rk%ai(2,2) = 0.29289321881345_pr
       rk%ai(3,2) = 0.70710678118655_pr
       rk%ai(3,3) = 0.29289321881345_pr

       allocate(rk%bh(3))
       rk%bh = rk%b ! no error estimate

    case('IE') ! implicit Euler

       rk%q = 2
       rk%p = 1
       rk%s = 2

       allocate(rk%b(2))
       rk%b(1) = 0.5_pr
       rk%b(2) = 0.5_pr

       allocate(rk%a(2,2))
       rk%a = zero
       rk%a(2,1) = 1._pr

       allocate(rk%ai(2,2))
       rk%ai = zero
       rk%ai(2,2) = 1._pr

       allocate(rk%bh(2))
       rk%bh = rk%b ! no error estimate

    case('RK3TVD') ! Shu (TVD for WENO)

       rk%q = 3
       rk%p = 2
       rk%s = 3

       allocate(rk%b(3))
       rk%b(1) = 1._pr/6._pr
       rk%b(2) = 1._pr/6._pr
       rk%b(3) = 2._pr/3._pr

       allocate(rk%a(3,2))
       rk%a = zero
       rk%a(2,1) = 1._pr
       rk%a(3,1) = 1._pr/4._pr
       rk%a(3,2) = 1._pr/4._pr

       allocate(rk%bh(3))
       rk%bh = rk%b ! no error estimate

    end select

    ! use consistency condition to set coefficients c(m)

    allocate(rk%c(rk%s))
    rk%c(1) = zero
    do m = 2,rk%s
       rk%c(m) = sum(rk%a(m,:))
    end do
    
  end subroutine init_rk


  subroutine destroy_rk(rk)
    ! DESTROY_RK destroys derived type rk
    ! 
    ! Modified: 21 August 2007

    implicit none

    ! I/O Parameters:
    ! RK = Runge-Kutta variables

    type(rk_type),intent(inout) :: rk

    ! deallocate memory assigned to allocatable arrays

    if (allocated(rk%a)) deallocate(rk%a)
    if (allocated(rk%ai)) deallocate(rk%ai)
    if (allocated(rk%b)) deallocate(rk%b)
    if (allocated(rk%bh)) deallocate(rk%bh)
    if (allocated(rk%c)) deallocate(rk%c)

  end subroutine destroy_rk


end module runge_kutta
