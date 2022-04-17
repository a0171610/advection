module time_module
  implicit none
  private

  integer(8), public ::  nstep = 80, hstep = 5
  real(8), public :: deltat = 21600.00
  ! クーラン数を1にするときはdeltatを18515.908736445574d0に、5.2にする時は96282.72542951698d0にする

  character(len=10), public :: &
    model = "direction", imethod = "sph", velocity = "nodiv", field = "ghill"
  logical, public :: conserve = .false., local_conserve = .true.
end module time_module