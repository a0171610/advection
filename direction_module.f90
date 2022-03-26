module direction_module

  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat, wgt
  use mass_module, only: mass_correct
  use time_module, only: conserve, local_conserve, velocity
  private
  
  integer(8), dimension(:, :), allocatable, private :: pa, qa, pb, qb, pc, qc, pd, qd
  real(8), allocatable, private :: A(:, :), B(:, :), C(:, :), D(:, :) 
  real(8), dimension(:,:), allocatable, private :: &
    gphi_old, dgphi, gphim, gphix, gphiy, gphixy, deplon, deplat
  real(8), dimension(:, :), allocatable, private :: midlonA, midlatA, midlonB, midlatB, midlonC, midlatC, midlonD, midlatD
  real(8), dimension(:, :), allocatable, private :: gumA, gvmA, gumB, gvmB, gumC, gvmC, gumD, gvmD
  real(8), dimension(:, :), allocatable, private :: dgphimA, dgphimB, dgphimC, dgphimD
  real(8), dimension(:, :), allocatable, private :: gmin, gmax, w
  complex(8), dimension(:,:), allocatable, private :: sphi1

  private :: update, bicubic_interpolation_set
  public :: direction_init, direction_timeint, direction_clean

contains

  subroutine direction_init()
    use time_module, only: deltat
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i,j

    allocate(sphi1(0:ntrunc, 0:ntrunc),gphi_old(nlon, nlat))
    allocate(gphim(nlon, nlat),dgphi(nlon, nlat))
    allocate(midlonA(nlon, nlat), midlatA(nlon, nlat), midlonB(nlon, nlat), midlatB(nlon, nlat))
    allocate(midlonC(nlon, nlat), midlatC(nlon, nlat), midlonD(nlon, nlat), midlatD(nlon, nlat))
    allocate(deplon(nlon, nlat), deplat(nlon, nlat), pa(nlon, nlat), qa(nlon, nlat), pb(nlon, nlat), qb(nlon, nlat))
    allocate(pc(nlon, nlat), qc(nlon, nlat), pd(nlon, nlat), qd(nlon, nlat))
    allocate(A(nlon, nlat), B(nlon, nlat), C(nlon, nlat), D(nlon, nlat))
    allocate(guma(nlon, nlat), gvma(nlon, nlat), gumb(nlon, nlat), gvmb(nlon, nlat))
    allocate(gumc(nlon, nlat), gvmc(nlon, nlat), gumd(nlon, nlat), gvmd(nlon, nlat))
    allocate(dgphimA(nlon, nlat), dgphimB(nlon, nlat), dgphimC(nlon, nlat), dgphimD(nlon, nlat))
    allocate(gphix(nlon, nlat), gphiy(nlon, nlat), gphixy(nlon, nlat))
    if (conserve) then
      allocate(gmax(nlon,nlat),gmin(nlon,nlat),w(nlon,nlat))
      do j=1, nlat
        w(:,j) = wgt(j)
      end do
      gmin(:,:) = minval(gphi)
      gmax(:,:) = maxval(gphi)
    endif

    call interpolate_init(gphi)

    call legendre_synthesis(sphi_old,gphi_old)
    gphi(:, :) = gphi_old(:, :)

    open(11, file="animation.txt")
    do i = 1, nlon
      do j = 1, nlat
          write(11,*) longitudes(i), latitudes(j), gphi(i, j)
      end do        
    end do
    call update(0.5d0*deltat, deltat)
    write(*, *) 'step = 1 ', "maxval = ", maxval(gphi), 'minval = ', minval(gphi)

  end subroutine direction_init

  subroutine direction_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(sphi1, gphi_old, gphim, dgphi, deplon, deplat)
    deallocate(midlonA, midlatA, midlonB, midlatB, midlonC, midlatC, midlonD, midlatD)
    deallocate(guma, gvma, gumb, gvmb, gumc, gvmc, gumd, gvmd)
    deallocate(A, B, C, D, pa, qa, pb, qb, pc, qc, pd, qd)
    deallocate(dgphimA, dgphimB, dgphimC, dgphimD)
    call interpolate_clean()

  end subroutine direction_clean

  subroutine direction_timeint()
    use time_module, only: nstep, deltat, hstep
    implicit none

    integer(8) :: i, j, k

    do i = 2, nstep
      call update((i-1)*deltat, 2.0d0*deltat)
      write(*, *) 'step = ', i, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
      if ( mod(i, hstep) == 0 ) then
        do j = 1, nlon
            do k = 1, nlat
              write(11,*) longitudes(j), latitudes(k), gphi(j, k)
            end do
        end do
      endif
    end do
    close(11)
    
  end subroutine direction_timeint

  subroutine update(t, dt)
    use uv_module, only: uv_nodiv, uv_div, uv_sbody
    use grid_module, only: lat_search
    use upstream_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat
    use interpolate_module, only: interpolate_set, interpolate_setd, find_stencil_
    use interpolate_module, only: interpolate_bicubic, interpolate_bilinear, interpolate_bilinear_ratio
    use interpolate_module, only: interpolate_dist, interpolate_dist_ratio, interpolate_setuv
    use interpolate_module, only: record_departure_point, record_i, record_j, record_d
    implicit none

    integer(8) :: i, j, k, dir, m, ilat
    real(8), intent(in) :: t, dt
    real(8) :: sum_tmp = 0.0d0

    select case(velocity)
    case("nodiv")
      call uv_nodiv(t,longitudes,latitudes,gu,gv)
    case("div")
      call uv_div(t,longitudes,latitudes,gu,gv)
    case("sbody")
      call uv_sbody(longitudes, latitudes, gu, gv)
    case default
      print *, "No matching model for", velocity
    end select
    call find_points(gu, gv, 0.5d0*dt, midlonA, midlatA, deplon, deplat)
    ! dtに0.5をかけているのは引数のdtが最初のステップ以外は2.0*deltatを渡しているから

    call interpolate_setuv(gu, gv)

    do i = 1, nlon
      do j = 1, nlat
        call interpolate_bilinear_ratio(deplon(i, j), deplat(i, j), A(i, j), B(i, j), C(i, j), D(i, j))
      end do
    end do

    if (local_conserve) then
      call record_departure_point(deplon, deplat)
      do i = 2, nlon
        do j = 2, nlat
          ilat = lat_search(deplat(i, j))
          !sum_tmp = A(i, j) + B(i, j - 1) + C(i - 1, j - 1) + D(i - 1, j)
          do k = 1, 5
            dir = record_d(i, j, k)
            if (dir == 1) then
              sum_tmp = sum_tmp + A(record_i(i, j, k), record_j(i, j, k))
            else if (dir == 2) then
              sum_tmp = sum_tmp + B(record_i(i, j, k), record_j(i, j, k))
            else if (dir == 3) then
              sum_tmp = sum_tmp + C(record_i(i, j, k), record_j(i, j, k))
            else if (dir == 4) then
              sum_tmp = sum_tmp + D(record_i(i, j, k), record_j(i, j, k))
            endif
          end do
          if (dir == 1) then
            A(record_i(i, j, k), record_j(i, j, k)) = A(record_i(i, j, k), record_j(i, j, k)) * wgt(ilat) / (wgt(i) * sum_tmp)
          else if (dir == 2) then
            B(record_i(i, j, k), record_j(i, j, k)) = B(record_i(i, j, k), record_j(i, j, k)) * wgt(ilat) / (wgt(i) * sum_tmp)
          else if (dir == 3) then
            C(record_i(i, j, k), record_j(i, j, k)) = C(record_i(i, j, k), record_j(i, j, k)) * wgt(ilat) / (wgt(i) * sum_tmp)
          else if (dir == 4) then
            D(record_i(i, j, k), record_j(i, j, k)) = D(record_i(i, j, k), record_j(i, j, k)) * wgt(ilat) / (wgt(i) * sum_tmp)
          endif
          !A(i, j) = A(i , j) * wgt(ilat) / wgt(j)
          !!B(i, j - 1) = B(i, j - 1) * wgt(ilat) / wgt(j)
          !C(i - 1, j - 1) = C(i - 1, j - 1) * wgt(ilat) / wgt(j)
          !D(i - 1, j) = D(i - 1, j) * wgt(ilat) / wgt(j)
        end do
      end do
    endif

    call set_niuv(dt)

    call legendre_synthesis(sphi_old, gphi_old)

    if (conserve) then
      gmin(:, :) = min(0.0d0, minval(gphi))
      gmax(:, :) = max(0.0d0, maxval(gphi))
    end if

    ! まずはgphiにbilinear法で求めた値を詰めていく
    call interpolate_set(gphi_old)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bilinear(deplon(i, j), deplat(i, j), gphi(i, j))
      end do
    end do

    ! dF/dlon
    call legendre_synthesis_dlon(sphi, dgphi)
    call bicubic_interpolation_set(dgphi) 
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    gphim(:, :) = 0.0d0
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(midlonA(i, j), midlatA(i, j), dgphimA(i, j))
        call interpolate_bicubic(midlonB(i, j), midlatB(i, j), dgphimB(i, j))
        call interpolate_bicubic(midlonC(i, j), midlatC(i, j), dgphimC(i, j))
        call interpolate_bicubic(midlonD(i, j), midlatD(i, j), dgphimD(i, j))

        gphim(i, j) = gphim(i, j) + A(i, j) * gumA(i, j) * dgphimA(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + B(i, j) * gumB(i, j) * dgphimB(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + C(i, j) * gumC(i, j) * dgphimC(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + D(i, j) * gumD(i, j) * dgphimD(i, j) / cos(latitudes(j))
      enddo
    enddo

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi) 
    call bicubic_interpolation_set(dgphi)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(midlonA(i, j), midlatA(i, j), dgphimA(i, j))
        call interpolate_bicubic(midlonB(i, j), midlatB(i, j), dgphimB(i, j))
        call interpolate_bicubic(midlonC(i, j), midlatC(i, j), dgphimC(i, j))
        call interpolate_bicubic(midlonD(i, j), midlatD(i, j), dgphimD(i, j))

        gphim(i, j) = gphim(i, j) + A(i, j) * gvmA(i, j) * dgphimA(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + B(i, j) * gvmB(i, j) * dgphimB(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + C(i, j) * gvmC(i, j) * dgphimC(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + D(i, j) * gvmD(i, j) * dgphimD(i, j) / cos(latitudes(j))
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
    use interpolate_module, only: find_stencil_
    implicit none

    real(8), intent(in) :: dt

    integer(8) :: i, j
    real(8) :: dlonr
    integer(8), dimension(4) :: tmp1, tmp2

    dlonr = 0.5d0 * nlon / math_pi
    call set_velocity_zero()
    do j = 1, nlat
      do i = 1, nlon
        ! find grid points near departure points
        call find_stencil_(deplon(i, j), deplat(i, j), tmp1, tmp2)
        pa(i, j) = tmp1(1); qa(i, j) = tmp2(1)
        pb(i, j) = tmp1(2); qb(i, j) = tmp2(2)
        pc(i, j) = tmp1(3); qc(i, j) = tmp2(3)
        pd(i, j) = tmp1(4); qd(i, j) = tmp2(4)


        call calc_niuv(dt, pa(i, j), qa(i, j), longitudes(i), latitudes(j), midlonA(i, j), midlatA(i, j), gumA(i, j), gvmA(i, j))
        call calc_niuv(dt, pb(i, j), qb(i, j), longitudes(i), latitudes(j), midlonB(i, j), midlatB(i, j), gumB(i, j), gvmB(i, j))
        call calc_niuv(dt, pc(i, j), qc(i, j), longitudes(i), latitudes(j), midlonC(i, j), midlatC(i, j), gumC(i, j), gvmC(i, j))
        call calc_niuv(dt, pd(i, j), qd(i, j), longitudes(i), latitudes(j), midlonD(i, j), midlatD(i, j), gumD(i, j), gvmD(i, j))
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
    gumA(:, :) = 0.0d0
    gumB(:, :) = 0.0d0
    gumC(:, :) = 0.0d0
    gumD(:, :) = 0.0d0
    gvmA(:, :) = 0.0d0
    gvmB(:, :) = 0.0d0
    gvmC(:, :) = 0.0d0
    gvmD(:, :) = 0.0d0
  end subroutine set_velocity_zero

  ! Ritchie1987 式(45)のu^* - Uをgumに詰める(gvmも)
  subroutine calc_niuv(dt, p, q, lon, lat, midlon, midlat, gum, gvm)
    use grid_module, only: latitudes => lat, longitudes => lon
    use interpolate_module, only: interpolate_bilinearuv
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    implicit none
    real(8), intent(in) :: dt
    integer(8), intent(in) :: p, q
    real(8), intent(in) :: lon, lat
    real(8), intent(out) :: gum, gvm
    real(8), intent(out) :: midlon, midlat

    real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xdot, ydot, zdot, lon_grid, lat_grid, u, v, b1


    lon_grid = longitudes(p)
    lat_grid = latitudes(q)
    call lonlat2xyz(lon_grid, lat_grid, xr, yr, zr)
    ! arrival points
    call lonlat2xyz(lon, lat, xg, yg, zg)

    b1 = 1.0d0 / sqrt( 2.0d0 * (1.0d0 + (xg*xr + yg*yr + zg*zr)) ) ! Ritchie1987 式(44)
    xm = b1 * (xg + xr)
    ym = b1 * (yg + yr)
    zm = b1 * (zg + zr)
    midlon = modulo(atan2(ym, xm) + pi2, pi2)
    midlat = asin(zm)

    xdot = (xg - xr) / dt
    ydot = (yg - yr) / dt
    zdot = (zg - zr) / dt
    call xyz2uv(xdot, ydot, zdot, midlon, midlat, u, v)  !Richie1987式(49)
    gum = gum + u
    gvm = gvm + v

    call interpolate_bilinearuv(midlon, midlat, u, v)

    gum = gum - u
    gvm = gvm - v

  end subroutine calc_niuv

end module direction_module