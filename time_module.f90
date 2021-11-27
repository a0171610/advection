module time_module
  implicit none
  private

  integer(8), public ::  nstep = 320, hstep = 10
  real(8), public :: deltat = 5400.0d0, kappa = 0.0d0

  character(len=7), public :: &
    model = "euler", imethod = "polint2", imethoduv = "bilin"
end module time_module
