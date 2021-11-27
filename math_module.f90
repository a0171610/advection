module math_module
  implicit none

! Mathematical constants
  real(8), parameter, public :: &
    math_pi = acos(-1.0d0), &
    math_pir = 1.0d0/math_pi, &
    math_pi2 = 2.0d0*math_pi, math_pih = 0.5d0*math_pi, &
    math_deg2rad = math_pi/180.0d0

end module math_module