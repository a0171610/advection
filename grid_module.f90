module grid_module
  implicit none
  private

  integer(8), parameter, public ::  ntrunc = 42, nlon = 128, nlat = 64

  complex(8), dimension(:,:), allocatable, public :: sphi, sphi_old
  real(8), dimension(:,:), allocatable, public :: gphi, gphi_initial, gu, gv
  real(8), dimension(:), allocatable, public :: lon, lat, coslat, coslatr, wgt

  public :: grid_init, grid_clean

contains

  subroutine grid_init()
    use math_module, only: pi2=>math_pi2
    use legendre_transform_module, only: &
      legendre_init, legendre_analysis
    use init_module, only: &
      init_ghill, init_ghill2, init_cbell2, init_scyli2, init_ccbel2
    use uv_module, only: uv_sbody, uv_nodiv, uv_div
    implicit none


    integer(8) :: i, m, n
    real(8) :: dlon

    allocate(lon(nlon), lat(nlat), coslat(nlat), coslatr(nlat), wgt(nlat))
    allocate(gphi(nlon,nlat), gphi_initial(nlon, nlat), gu(nlon,nlat), gv(nlon,nlat))
    allocate(sphi(0:ntrunc,0:ntrunc), sphi_old(0:ntrunc,0:ntrunc))
 
    dlon = pi2/nlon
    do i=1, nlon
      lon(i) = dlon * dble(i-1)
    end do
    call legendre_init(nlon,nlat,ntrunc,lat,wgt)
    coslat(:) = cos(lat(:))
    coslatr(:) = 1.0d0 / coslat(:)

    call init_ghill(lon,lat,gphi)
    gphi_initial(:, :) = gphi(:, :)
    call legendre_analysis(gphi, sphi)
    do m = 0, ntrunc
      do n = m, ntrunc
        sphi_old(n, m) = sphi(n, m)
      enddo
    enddo

    call uv_sbody(lon, lat, gu, gv)
      

  end subroutine grid_init

  subroutine grid_clean
    implicit none

    deallocate(lon, lat, coslat, coslatr, gphi, gu, gv, sphi, sphi_old)

  end subroutine grid_clean

end module grid_module