module direction16_module

  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat, wgt
  use field_module, only : X, Y
  use mass_module, only: mass_correct
  use time_module, only: conserve
  private
  
  integer(8), dimension(:, :, :), allocatable, private :: p, q
  real(8), allocatable, private :: weight(:, :, :)
  real(8), dimension(:,:), allocatable, private :: &
    gphi_old, dgphi, gphim, gphix, gphiy, gphixy, deplon, deplat
  real(8), dimension(:, :, :), allocatable, private :: midlon, midlat
  real(8), dimension(:, :, :), allocatable, private :: gum, gvm
  real(8), dimension(:, :, :), allocatable, private :: dgphim
  real(8), dimension(:, :), allocatable, private :: gmin, gmax, w
  complex(8), dimension(:,:), allocatable, private :: sphi1

  private :: update, bicubic_interpolation_set
  public :: direction16_init, direction16_timeint, direction16_clean

contains

  subroutine direction16_init()
    use time_module, only: deltat
    use interpolate16_module, only: interpolate16_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i,j

    allocate(sphi1(0:ntrunc, 0:ntrunc),gphi_old(nlon, nlat))
    allocate(gphim(nlon, nlat),dgphi(nlon, nlat))
    allocate(midlon(16, nlon, nlat), midlat(16, nlon, nlat))
    allocate(deplon(nlon, nlat), deplat(nlon, nlat), p(16, nlon, nlat), q(16, nlon, nlat))
    allocate(weight(16, nlon, nlat))
    allocate(gum(16, nlon, nlat), gvm(16, nlon, nlat), dgphim(16, nlon, nlat))

    allocate(gphix(nlon, nlat), gphiy(nlon, nlat), gphixy(nlon, nlat))
    if (conserve) then
      allocate(gmax(nlon,nlat),gmin(nlon,nlat),w(nlon,nlat))
      do j=1, nlat
        w(:,j) = wgt(j)
      end do
      gmin(:,:) = minval(gphi)
      gmax(:,:) = maxval(gphi)
    endif

    call interpolate16_init(gphi)

    call legendre_synthesis(sphi_old,gphi_old)
    gphi(:, :) = gphi_old(:, :)

    open(11, file="animation.txt")
    do i = 1, nlon
      do j = 1, nlat
          write(11,*) X(i, j), Y(i, j), gphi(i, j)
      end do        
    end do
    call update(deltat)
    write(*, *) 'step = 1 ', "maxval = ", maxval(gphi), 'minval = ', minval(gphi)

  end subroutine direction16_init

  subroutine direction16_clean()
    use interpolate16_module, only: interpolate16_clean
    implicit none

    deallocate(sphi1, gphi_old, gphim, dgphi, deplon, deplat)
    deallocate(midlon, midlat, gum, gvm)
    deallocate(weight, p, q)
    call interpolate16_clean()

  end subroutine direction16_clean

  subroutine direction16_timeint()
    use time_module, only: nstep, deltat, hstep
    implicit none

    integer(8) :: i, j, k

    do i = 2, nstep
      call update(2.0d0*deltat)
      write(*, *) 'step = ', i, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
      if ( mod(i, hstep) == 0 ) then
        do j = 1, nlon
            do k = 1, nlat
              write(11,*) X(j, k), Y(j, k), gphi(j, k)
            end do
        end do
      endif
    end do
    close(11)
    
  end subroutine direction16_timeint

  subroutine update(dt)
    use upstream_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat
    use interpolate16_module, only: interpolate16_set, interpolate16_setd, find_stencil_16
    use interpolate16_module, only: interpolate16_bicubic, interpolate16_dist, interpolate16_dist_ratio
    use interpolate16_module, only: interpolate16_dist, interpolate16_dist_ratio
    implicit none

    integer(8) :: i, j, m, k
    real(8), intent(in) :: dt

    call find_points(gu, gv, 0.5*dt, midlon(1,:,:), midlat(1,:,:), deplon, deplat)
    ! dtに0.5をかけているのは引数のdtが最初のステップ以外は2.0*deltatを渡しているから

    do i = 1, nlon
        do j = 1, nlat
          call interpolate16_dist_ratio(deplon(i, j), deplat(i, j), weight(:, i, j))
        end do
    end do

    call set_niuv(dt)

    call legendre_synthesis(sphi_old, gphi_old)

    if (conserve) then
      gmin(:, :) = min(0.0d0, minval(gphi))
      gmax(:, :) = max(0.0d0, maxval(gphi))
    end if

    ! まずはgphiにbilinear法で求めた値を詰めていく
    call interpolate16_set(gphi_old)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate16_dist(deplon(i, j), deplat(i, j), gphi(i, j))
      end do
    end do

    ! dF/dlon
    call legendre_synthesis_dlon(sphi, dgphi)
    call bicubic_interpolation_set(dgphi)
    call interpolate16_set(dgphi)
    call interpolate16_setd(gphix, gphiy, gphixy)
    gphim(:, :) = 0.0d0
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, 16
          call interpolate16_bicubic(midlon(k, i, j), midlat(k, i, j), dgphim(k, i, j))
          gphim(i, j) = gphim(i, j) + weight(k, i, j) * gum(k, i, j) * dgphim(k, i, j) / cos(latitudes(j))
        end do
      enddo
    enddo

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi) 
    call bicubic_interpolation_set(dgphi)
    call interpolate16_set(dgphi)
    call interpolate16_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        do k = 1, 16
          call interpolate16_bicubic(midlon(k, i, j), midlat(k, i, j), dgphim(k, i, j))
          gphim(i, j) = gphim(i, j) + weight(k, i, j) * gvm(k, i, j) * dgphim(k, i, j) / cos(latitudes(j))
        end do
      enddo
    enddo

    gphi(:, :) = gphi(:, :) + dt * gphim(:, :)

    if (conserve) then
      call mass_correct(gphi,gphi_old,gmax,gmin,w)
    end if
    
! time step
    call legendre_analysis(gphi, sphi1)
    do m = 0, ntrunc
      sphi_old(m : ntrunc, m) = sphi(m : ntrunc, m)       
      sphi(m : ntrunc, m) = sphi1(m : ntrunc, m)
    enddo

  end subroutine update

  subroutine set_niuv(dt)
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    use uv_module, only: uv_sbody_calc
    use interpolate16_module, only: find_stencil_16
    implicit none

    real(8), intent(in) :: dt

    integer(8) :: i, j, k
    real(8) :: dlonr

    dlonr = 0.5d0 * nlon / math_pi
    call set_velocity_zero()
    do j = 1, nlat
      do i = 1, nlon
        ! find grid points near departure points
        call find_stencil_16(deplon(i, j), deplat(i, j), p(:, i, j), q(:, i, j))
        do k = 1, 16
          call calc_niuv(dt, p(k, i, j), q(k, i, j), longitudes(i), latitudes(j), midlon(k, i, j), midlat(k, i, j), &
            gum(k, i, j), gvm(k, i, j))
        end do
      end do
    end do

  end subroutine  set_niuv

  subroutine bicubic_interpolation_set(f)
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis_dlat, legendre_synthesis_dlon, &
      legendre_synthesis_dlonlat
    implicit none
    integer(8) :: j
    real(8), intent(in) :: f(nlon, nlat)

    call legendre_analysis(f, sphi1)
    call legendre_synthesis_dlon(sphi1, gphix)
    call legendre_synthesis_dlat(sphi1, gphiy)
    call legendre_synthesis_dlonlat(sphi1, gphixy)
    do j = 1, nlat
      gphiy(:,j) = gphiy(:,j) / cos(latitudes(j))
      gphixy(:,j) = gphixy(:,j) / cos(latitudes(j))
    end do

  end subroutine bicubic_interpolation_set

  subroutine set_velocity_zero
    implicit none
    gum(:, :, :) = 0.0d0
    gvm(:, :, :) = 0.0d0
  end subroutine set_velocity_zero

  ! Ritchie1987 式(45)のu^* - Uをgumに詰める(gvmも)
  subroutine calc_niuv(dt, p1, q1, lon, lat, midlon1, midlat1, gum1, gvm1)
    use grid_module, only: latitudes => lat, longitudes => lon
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    use uv_module, only: uv_sbody_calc
    implicit none
    real(8), intent(in) :: dt
    integer(8), intent(in) :: p1, q1
    real(8), intent(in) :: lon, lat
    real(8), intent(out) :: gum1, gvm1
    real(8), intent(out) :: midlon1, midlat1

    real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xdot, ydot, zdot, lon_grid, lat_grid, u, v, b


    lon_grid = longitudes(p1)
    lat_grid = latitudes(q1)
    call lonlat2xyz(lon_grid, lat_grid, xr, yr, zr)
    ! arrival points
    call lonlat2xyz(lon, lat, xg, yg, zg)

    b = 1.0d0 / sqrt( 2.0d0 * (1.0d0 + (xg*xr + yg*yr + zg*zr)) ) ! Ritchie1987 式(44)
    xm = b * (xg + xr)
    ym = b * (yg + yr)
    zm = b * (zg + zr)
    midlon1 = modulo(atan2(ym, xm) + pi2, pi2)
    midlat1 = asin(zm)

    xdot = (xg - xr) / dt
    ydot = (yg - yr) / dt
    zdot = (zg - zr) / dt
    call xyz2uv(xdot, ydot, zdot, midlon1, midlat1, u, v)  !Richie1987式(49)
    gum1 = gum1 + u
    gvm1 = gvm1 + v

    call uv_sbody_calc(midlon1, midlat1, u, v)
    gum1 = gum1 - u
    gvm1 = gvm1 - v

  end subroutine calc_niuv

end module direction16_module