module nisl_2step_module

  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat, wgt
  private
  
  integer(8), allocatable, private :: p(:,:), q(:,:)
  real(8), dimension(:,:), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphix, gphiy, gphixy, &
    midlon, midlat, deplon, deplat, gum, gvm
  complex(8), dimension(:,:), allocatable, private :: sphi1

  private :: update
  public :: nisl_2step_init, nisl_2step_timeint, nisl_2step_clean

contains

  subroutine nisl_2step_init()
    use time_module, only: deltat
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i,j

    allocate(sphi1(0:ntrunc, 0:ntrunc),gphi_old(nlon, nlat))
    allocate(gphim(nlon, nlat),dgphi(nlon, nlat),dgphim(nlon, nlat))
    allocate(midlon(nlon, nlat), midlat(nlon, nlat))
    allocate(deplon(nlon, nlat), deplat(nlon, nlat), p(nlon, nlat), q(nlon, nlat))
    allocate(gum(nlon, nlat), gvm(nlon, nlat))
    allocate(gphix(nlon, nlat), gphiy(nlon, nlat), gphixy(nlon, nlat))

    call interpolate_init(gphi)

    call legendre_synthesis(sphi_old,gphi_old)
    gphi(:, :) = gphi_old(:, :)

    open(11, file="animation.txt")
    do i = 1, nlon
      do j = 1, nlat
          write(11,*) longitudes(i), latitudes(j), gphi(i, j)
      end do        
    end do
    call update(deltat)
    write(*, *) 'step = ', 1, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)

  end subroutine nisl_2step_init

  subroutine nisl_2step_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(sphi1,gphi_old,gphim,dgphi,dgphim,gum,gvm, &
      midlon,midlat,deplon,deplat,p,q)
    call interpolate_clean()

  end subroutine nisl_2step_clean

  subroutine nisl_2step_timeint()
    use time_module, only: nstep, deltat, hstep, field
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    do i = 2, nstep
      call update(2.0d0*deltat)
      write(*, *) 'step = ', i, "maxval = ", maxval(gphi), 'minval = ', minval(gphi)
      if ( mod(i, hstep) == 0 ) then
        do j = 1, nlon
            do k = 1, nlat
              write(11,*) longitudes(j), latitudes(k), gphi(j, k)
            end do
        end do
      endif
      if (i == nstep / 2 .and. field == "cbell2") then
        open(10, file="log_cbell.txt")
        do j = 1, nlon
          do k = 1, nlat
            write(10,*) gphi(j, k)
          enddo
        enddo
        close(10)
      endif
      if (i == nstep / 2 .and. field == "ccbell2") then
        open(10, file="log_ccbell.txt")
        do j = 1, nlon
          do k = 1, nlat
            write(10,*) wgt(k), gphi(j, k)
          enddo
        enddo
        close(10)
      endif
    end do
    close(11)
    
  end subroutine nisl_2step_timeint

  subroutine update(dt)
    use uv_module, only: uv_nodiv, uv_div
    use upstream_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate_module, only: interpolate_set, interpolate_bilinear, interpolate_setd, interpolate_bicubic
    use interpolate_module, only: interpolate_setuv
    implicit none

    integer(8) :: i, j, m
    real(8), intent(in) :: dt
    real(8) :: b(nlon * nlat), x(nlon * nlat)

    call interpolate_setuv(gu, gv)

    call find_points(gu, gv, 0.5d0*dt, midlon, midlat, deplon, deplat)   ! dtに0.5をかけているのは引数のdtが最初のステップ以外は2.0*deltatを渡しているから
    call set_niuv(dt)

    call legendre_synthesis(sphi_old, gphi_old)

    do j = 1, nlat
      do i = 1, nlon
        gphi(i,j) = gphi_old(p(i, j), q(i, j))
      end do
    end do

    ! dF/dlon
    call legendre_synthesis_dlon(sphi_old, dgphi)
    do j = 1, nlat
      do i = 1, nlon
        gphim(i, j) = gum(i, j) * dgphi(p(i, j), q(i, j)) / cos(latitudes(j)) ! gum: -u'
      enddo
    enddo

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi_old, dgphi)
    do j = 1, nlat
      do i = 1, nlon
        gphim(i, j) = gphim(i, j) + gvm(i, j) * dgphi(p(i, j), q(i, j)) / cos(latitudes(j)) ! gvm: -v'
      enddo
    enddo

    gphi = gphi + dt * gphim * 0.5d0

    do i = 1, nlon
      do j = 1, nlat
        m = i + (j - 1) * nlon
        b(m) = gphi(i, j)
      end do
    end do

    call solve_sparse_matrix(dt, b, x)
    do m = 1, nlon*nlat
      i = mod(m, nlon)
      if (i == 0) then
        i = nlon
      endif
      j = int((m-i)/nlon) + 1
      gphi(i, j) = x(m)
    end do

! time filter
    call legendre_analysis(gphi, sphi1)
    do m = 0, ntrunc
      sphi_old(m : ntrunc, m) = sphi(m : ntrunc, m)       
      sphi(m : ntrunc, m) = sphi1(m : ntrunc, m)
    enddo

  end subroutine update

  subroutine solve_sparse_matrix(dt, b, x)
    use lsqr_module, only: lsqr_solver_ez
    use sphere_module, only: orthodrome
    use grid_module, only: pole_regrid
    use math_module, only : pih=>math_pih, pir=>math_pir
    implicit none

    integer(8), parameter :: sz = nlat * nlon
    real(8), intent(in) :: dt
    integer(8) :: x1, y1, x2, y2, x3, y3, x4, y4
    real(8), intent(in) :: b(sz)
    real(8), intent(out) :: x(sz)
    type(lsqr_solver_ez) :: solver
    integer :: istop
    integer :: i, j, id, row, col
    integer, allocatable :: icol(:), irow(:)
    real(8), allocatable :: a(:)
    real(8) :: val, dlonr, dlat, eps

    allocate( icol(sz * 5), irow(sz * 5), a(sz * 5) )

    dlonr = longitudes(3) - longitudes(1)

    id = 1
    do i = 1, nlon
      do j = 1, nlat
        row = i + (j-1) * int(nlon)
        x1 = i + 1; y1 = j
        x2 = i - 1; y2 = j
        x3 = i; y3 = j + 1
        x4 = i; y4 = j - 1
        call pole_regrid(x1, y1)
        call pole_regrid(x2, y2)
        call pole_regrid(x3, y3)
        call pole_regrid(x4, y4)

        val = gum(i, j) * dt / (2.0d0 * dlonr * cos(latitudes(j)))
        col = int(x1 + (y1 - 1) * nlon)
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1

        val = -gum(i, j) * dt / (2.0d0 * dlonr * cos(latitudes(j)))
        col = int(x2 + (y2 - 1) * nlon)
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1 

        eps = pih - latitudes(1)
        if (j == 1) then
          dlat = pih + eps - latitudes(2)
        elseif (j == nlat) then
          dlat = -pih-eps-latitudes(nlat-1)
        else
          dlat = latitudes(j + 1) - latitudes(j - 1)
        endif
        val = gvm(i, j) * dt / (2.0d0 * dlat)
        col = int(x3 + (y3 - 1) * nlon)
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1

        val = -gvm(i, j) * dt / (2.0d0 * dlat)
        col = int(x4 + (y4 - 1) * nlon)
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1
      end do
    end do

    do i = 4*sz+1, 5*sz
      irow(i) = int(i - 4*sz)
      icol(i) = int(i - 4*sz)
      a(i) = 1.0d0
    end do

    call solver%initialize(int(sz), int(sz), a, irow, icol) ! use defaults for other optional inputs
    call solver%solve(b, 0.0d0, x, istop)       ! solve the linear system
    write(*,*) 'istop = ', istop
  end subroutine solve_sparse_matrix

  subroutine set_niuv(dt)
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    use grid_module, only: pole_regrid
    implicit none

    real(8), intent(in) :: dt

    integer(8) :: i, j
    real(8) :: dlonr

    dlonr = 0.5d0 * nlon / math_pi
    gum(:, :) = 0.0d0
    gvm(:, :) = 0.0d0
    do j = 1, nlat
      do i = 1, nlon
        ! find grid points near departure points

        p(i, j) = int(anint( deplon(i, j) * dlonr + 1.0d0 ))
        if ( p(i,j) > nlon ) then
          p(i, j) = p(i, j) - nlon
        end if
        ! lat = (J+1-2j)pi/(2J+1)
        q(i, j) = int(anint( 0.5d0 * (nlat + 1.0d0 - (2.0d0*dble(nlat)+1.0d0)*deplat(i, j) / math_pi) ))  !latitudesは大きい順で詰められているので注意
        call pole_regrid(p(i, j), q(i, j))
        call calc_niuv(dt, p(i, j), q(i, j), longitudes(i), latitudes(j), midlon(i, j), midlat(i, j), gum(i, j), gvm(i, j))
      end do
    end do
        
  end subroutine  set_niuv

  ! Ritchie1987 式(45)のu^* - Uをgumに詰める(gvmも)
  subroutine calc_niuv(dt, p1, q1, lon, lat, midlon1, midlat1, gum1, gvm1)
    use grid_module, only: latitudes => lat, longitudes => lon
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    use interpolate_module, only: interpolate_bilinearuv
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

    call interpolate_bilinearuv(midlon1, midlat1, u, v)
    gum1 = gum1 - u
    gvm1 = gvm1 - v

  end subroutine calc_niuv

end module nisl_2step_module