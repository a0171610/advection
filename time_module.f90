module time_module
  implicit none
  private

  integer(8), public ::  nstep = 80, hstep = 10
  real(8), public :: deltat = 21600.0d0, etf = 0.0d0, kappa = 0.0d0

  character(len=7), public :: &
    model = "nisl", imethod = "polint2", imethoduv = "bilin"


  public :: time_init

contains

  subroutine time_init()
    use grid_module, only: ntrunc
    use planet_module, only: d=>day_in_sec
    implicit none

    real(8) :: tau = 0.0d0

    if (tau<=0.0d0) then
      kappa = 0.0d0
    else
      kappa = 1.0d0/(tau*d*(ntrunc*(ntrunc+1.0d0))**2)
    end if
    print *, "tau=", tau, " kappa=", kappa

  end subroutine time_init

end module time_module
