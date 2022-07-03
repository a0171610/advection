module time_module
  implicit none
  private

  integer(8), public ::  nstep = 120, hstep = 5
  !integer(8), public ::  nstep = 80, hstep = 5

  real(8), public :: deltat = 0.0833333333333333/2.0d0
  !real(8), public :: deltat = 21600.0d0

  character(len=10), public :: model = "nisl_2step", imethod = "fd", velocity = "nodiv", field = "cbell2"
  logical, public :: conserve = .false., local_conserve = .false.
end module time_module