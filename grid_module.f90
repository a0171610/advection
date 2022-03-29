module grid_module
  implicit none
  private

  integer(8), parameter, public ::  ntrunc = 42, nlon = 128, nlat = 64

  complex(8), dimension(:,:), allocatable, public :: sphi, sphi_old
  real(8), dimension(:,:), allocatable, public :: gphi, gphi_initial, gu, gv
  real(8), dimension(:), allocatable, public :: lon, lat, coslat, coslatr, wgt

  public :: grid_init, grid_clean, pole_regrid, grid_id

contains

  subroutine grid_init()
    use math_module, only: pi2=>math_pi2
    use legendre_transform_module, only: &
      legendre_init, legendre_analysis
    use init_module, only: &
      init_ghill, init_ghill2, init_cbell2, init_scyli2, init_ccbel2
    use uv_module, only: uv_sbody, uv_nodiv, uv_div
    use time_module, only: velocity, field
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

    select case(field)
      case("ghill")
        call init_ghill(lon,lat,gphi)
      case("ghill2")
        call init_ghill2(lon,lat,gphi)
      case("cbell2")
        call init_cbell2(lon,lat,gphi)
      case("scyli2")
        call init_scyli2(lon,lat,gphi)
      case("ccbel2")
        call init_ccbel2(lon,lat,gphi)
      case default
        print *, "No matching initial field"
      stop
    end select

    gphi_initial(:, :) = gphi(:, :)
    call legendre_analysis(gphi, sphi)
    do m = 0, ntrunc
      do n = m, ntrunc
        sphi_old(n, m) = sphi(n, m)
      enddo
    enddo

    select case(velocity)
    case("sbody")
      call uv_sbody(lon, lat, gu, gv)
    case("nodiv")
      call uv_nodiv(0.0d0, lon, lat, gu, gv)
    case("div")
      call uv_div(0.0d0, lon, lat, gu, gv)
    case default
      print *, "No matching initial wind"
      stop
  end select
      

  end subroutine grid_init

  subroutine grid_clean
    implicit none

    deallocate(lon, lat, coslat, coslatr, gphi, gu, gv, sphi, sphi_old)

  end subroutine grid_clean

  subroutine pole_regrid(x, y)
    implicit none
    integer(8), intent(inout) :: x, y
    if (x < 1) then
        x = x + nlon
    endif
    if (x > nlon) then
        x = x - nlon
    endif
    if ( y < 1 ) then
        y = 1 - y
        if ( x > nlon/2) then
            x = x - nlon/2
        else
            x = x + nlon/2
        endif
    endif
    if ( y > nlat ) then
        y = 2 * nlat + 1 - y
        if ( x > nlon/2) then
            x = x - nlon/2
        else
            x = x + nlon/2
        endif
    endif
end subroutine pole_regrid

  ! 経度lon, 緯度latが与えられた時、 lon(lo1) <= lon < lon(lo2), lat(la1) <= lat < lat(la2) を満たすlo1, lo2, la1, la2を返す
  subroutine grid_id(long, lati, lo, la)
    implicit none
    real(8), intent(in) :: long, lati
    integer(8), intent(out) :: lo, la
    integer(8) :: i

    lo = -1
    do i = 1, nlon - 1
      if (lon(i) <= long .and. long < lon(i + 1)) then
        lo = i
        exit
      endif
    end do

    if (lo == -1) then
      lo = nlon
    endif

    la = -1
    do i = 1, nlat - 1
      if (lat(i + 1) <= lati .and. lati <lat(i)) then
        la = i
        exit
      endif
    end do

    if (la == -1 .and. lati < lat(nlat)) then
      la = 0
    endif

    if (la == -1 .and. lat(nlat) < lati) then
      la = nlat
    endif

  end subroutine grid_id

end module grid_module