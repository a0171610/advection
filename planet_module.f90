module planet_module
  use math_module, only: pi=>math_pi
  implicit none

! Geophysical constants
  real(8), public :: &
    planet_radius = 6.371d6, &
    day_in_sec = 86400.0d0, angular_velocity

  public :: planet_init

contains

  subroutine planet_init()
    implicit none

    angular_velocity= 2.0d0*pi/day_in_sec

  end subroutine planet_init

end module planet_module