! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module constants

  ! CONSTANTS contains compile-time constants, including the precision 
  ! of real variables used for computations and for writing binary data to files
  ! 
  ! Modified: 30 August 2006

  implicit none

  ! Global Parameters:
  ! INT4 = 4 byte integer
  ! INT8 = 8 byte integer
  ! PIN = precision of integer variables for calculations
  ! REAL4 = single precision
  ! REAL8 = double precision
  ! REAL16 = quad precision
  ! PSAV = precision in which to output binary data
  ! PR = precision of real variables for calculations
  ! PI = the constant pi
  ! TWOPI = 2*pi
  ! PIOVERTWO = pi/2
  ! ZERO = zero in real form
  ! ONE = one in real form
  ! TWO = two in real form
  ! THREE = three in real form
  ! FOUR = four in real form
  ! FIVE = five in real form
  ! SIX = six in real form
  ! HALF = 1/2 in real form
  ! IMG = square root of -1
  ! ZEROC = zero in complex form
  ! ONEC = one in complex form
  ! LOGTWO = log(2) in real form
  ! MINUS = -1 in real form

  integer,parameter :: int4 = selected_int_kind(9)
  integer,parameter :: int8 = selected_int_kind(18)
  integer,parameter :: real4  = selected_real_kind(6,37)
  integer,parameter :: real8  = selected_real_kind(15,307)

  integer,parameter :: pin = int4
  integer,parameter :: pr = real8
  integer,parameter :: psav = real4

  real(pr),parameter :: zero = 0._pr
  real(pr),parameter :: one = 1._pr
  real(pr),parameter :: two = 2._pr
  real(pr),parameter :: three = 3._pr
  real(pr),parameter :: four = 4._pr
  real(pr),parameter :: five = 5._pr
  real(pr),parameter :: six = 6._pr

  real(pr),parameter :: half = 0.5_pr
  real(pr),parameter :: minus = -1._pr

  real(pr),parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pr
  real(pr),parameter :: twopi = 2._pr*pi
  real(pr),parameter :: piovertwo = pi/2._pr
  real(pr),parameter :: logtwo = 0.69314718055994530941723212145817656807550013436025_pr
  complex(pr),parameter :: img = (0._pr,1._pr)
  complex(pr),parameter :: zeroc = (0._pr,0._pr)
  complex(pr),parameter :: onec = (1._pr,0._pr)

end module constants
