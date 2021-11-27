module time_module
  implicit none
  private

  integer(8), public ::  nstep = 80, hstep = 10
  real(8), public :: deltat = 21600.0d0, kappa = 0.0d0

  character(len=7), public :: &
    model = "slag", imethod = "sph", imethoduv = "bilin"
end module time_module
