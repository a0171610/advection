module nisl_2step_module

  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat, wgt
  use time_module, only: deltat, velocity
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
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i,j

    allocate(sphi1(0:ntrunc, 0:ntrunc), gphi_old(nlon, nlat))
    allocate(gphim(nlon, nlat),dgphi(nlon, nlat),dgphim(nlon, nlat))
    allocate(midlon(nlon, nlat), midlat(nlon, nlat))
    allocate(deplon(nlon, nlat), deplat(nlon, nlat), p(nlon, nlat), q(nlon, nlat))
    allocate(gum(nlon, nlat), gvm(nlon, nlat))
    allocate(gphix(nlon, nlat), gphiy(nlon, nlat), gphixy(nlon, nlat))

    call interpolate_init(gphi)

    call legendre_synthesis(sphi_old,gphi_old)
    gphi(:, :) = gphi_old(:, :)

    do i = 1, nlon
      do j = 1, nlat
        midlat(i, j) = latitudes(j)
        midlon(i, j) = longitudes(i)
      end do
    end do

    open(11, file="animation.txt")
    do i = 1, nlon
      do j = 1, nlat
          write(11,*) longitudes(i), latitudes(j), gphi(i, j)
      end do        
    end do
  end subroutine nisl_2step_init

  subroutine nisl_2step_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(sphi1, gphi_old,gphim,dgphi,dgphim,gum,gvm, &
      midlon,midlat,deplon,deplat,p,q)
    call interpolate_clean()

  end subroutine nisl_2step_clean

  subroutine nisl_2step_timeint()
    use time_module, only: nstep, hstep, field
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i, j, k

    do i = 1, nstep
      call update((i- 1.0d0)*deltat)
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

  ! dt = leapfrog法の+と-の時刻差, t=中央の時刻
  subroutine update(t)
    use grid_module, only: pole_regrid
    use uv_module, only: uv_nodiv, uv_div
    use sphere_module, only: orthodrome
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate_module, only: interpolate_set, interpolate_bilinear, interpolate_setd, interpolate_bicubic
    use interpolate_module, only: interpolate_setuv
    use math_module, only : pir=>math_pir, pih=>math_pih
    use upstream_module, only: find_points
    implicit none

    integer(8) :: i, j, m
    real(8) :: b(nlon * nlat), x(nlon * nlat)
    real(8), intent(in) :: t
    real(8) :: eps
    real(8), dimension(nlon) :: gphitmp

    call legendre_synthesis(sphi_old, gphi_old)

    call calculate_residual_velocity(t)

    do j = 1, nlat
      do i = 1, nlon
        gphi(i,j) = gphi_old(p(i, j), q(i, j))
      end do
    end do

    !gphix(1,:) = (gphi_old(2,:) - gphi_old(nlon,:)) / (longitudes(3) - longitudes(1))
    !gphix(nlon,:) = (gphi_old(1,:) - gphi_old(nlon-1,:)) / (longitudes(3) - longitudes(1))
    !do i=2, nlon-1
    !  gphix(i,:) = (gphi_old(i+1,:) - gphi_old(i-1,:)) / (longitudes(3) - longitudes(1))
    !end do
    ! d/dphi
    !eps = pih-latitudes(1)
    !gphitmp = cshift(gphi_old(:,1),nlon/2)
    !gphiy(:,1) = (gphitmp-gphi_old(:,2))/(pih+eps-latitudes(2))
    !gphitmp = cshift(gphi_old(:,nlat),nlon/2)
    !gphiy(:,nlat) = (gphitmp-gphi_old(:,nlat-1))/(-pih-eps-latitudes(nlat-1))
    !do j=2, nlat-1
    !  gphiy(:,j) = (gphi_old(:,j+1)-gphi_old(:,j-1))/(latitudes(j+1)-latitudes(j-1))
    !end do

    !do i = 1, nlon
    !  do j = 1, nlat
    !    gphim(i, j) = gum(i, j) * gphix(i, j) / cos(deplat(i, j)) +  gvm(i, j) * gphiy(i, j)
    !    gphi(i, j) = gphi(i, j) + deltat * gphim(i, j) * 0.5d0
    !  end do
    !end do

        ! dF/dlon
    call legendre_synthesis_dlon(sphi_old, dgphi)
    call bicubic_interpolation_set(dgphi)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(longitudes(p(i, j)), latitudes(q(i, j)), dgphim(i, j))
      enddo
      gphim(:, j) = gum(:, j) * dgphim(:,j) / cos(latitudes(j)) ! gum: -u'
    enddo

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi_old, dgphi) 
    call bicubic_interpolation_set(dgphi)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(longitudes(p(i, j)), latitudes(q(i, j)), dgphim(i, j))
      enddo
      gphim(:, j) = gphim(:, j) + gvm(:, j) * dgphim(:,j) / cos(latitudes(j)) * 0.5d0 ! gvm: -v'
    enddo

    do i = 1, nlon
      do j = 1, nlat
        gphi(i, j) = gphi(i, j) + deltat * gphim(i, j)
      end do
    end do

    do i = 1, nlon
      do j = 1, nlat
        m = i + (j - 1) * nlon
        b(m) = gphi(i, j)
      end do
    end do

    call solve_sparse_matrix(t, b, x)
    do m = 1, nlon*nlat
      i = mod(m, nlon)
      if (i == 0) then
        i = nlon
      endif
      j = int((m-i)/nlon) + 1
      gphi(i, j) = x(m)
    end do

    call legendre_analysis(gphi, sphi)
    do m = 0, ntrunc
      sphi_old(m : ntrunc, m) = sphi(m : ntrunc, m)
    enddo

  end subroutine update

  subroutine solve_sparse_matrix(t, b, x)
    use lsqr_module, only: lsqr_solver_ez
    use sphere_module, only: orthodrome
    use grid_module, only: pole_regrid
    use math_module, only : pih=>math_pih, pir=>math_pir
    implicit none

    integer(8), parameter :: sz = nlat * nlon
    real(8), intent(in) :: t
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

    call calculate_residual_velocity(t + deltat)

    id = 1
    do i = 1, nlon
      do j = 1, nlat
        col = i + (j-1) * int(nlon)
        x1 = i + 1; y1 = j
        x2 = i - 1; y2 = j
        x3 = i; y3 = j + 1
        x4 = i; y4 = j - 1
        call pole_regrid(x1, y1)
        call pole_regrid(x2, y2)
        call pole_regrid(x3, y3)
        call pole_regrid(x4, y4)

        val = gum(i, j) * deltat / (2.0d0 * dlonr * cos(latitudes(j)))
        row = int(x1 + (y1 - 1) * nlon)
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1

        val = -gum(i, j) * deltat / (2.0d0 * dlonr * cos(latitudes(j)))
        row = int(x2 + (y2 - 1) * nlon)
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
        val = gvm(i, j) * deltat / (2.0d0 * dlat)
        row = int(x3 + (y3 - 1) * nlon)
        irow(id) = row
        icol(id) = col
        a(id) = val
        id = id + 1

        val = -gvm(i, j) * deltat / (2.0d0 * dlat)
        row = int(x4 + (y4 - 1) * nlon)
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
  end subroutine solve_sparse_matrix

  subroutine calculate_residual_velocity(t)
    use upstream_module, only: find_points
    use math_module, only: pi2=>math_pi2
    use interpolate_module, only: interpolate_bilinearuv, interpolate_setuv
    use sphere_module, only: xyz2uv, lonlat2xyz
    use uv_module, only: uv_div, uv_nodiv
    implicit none
    real(8), intent(in) :: t
    integer(8) :: i, j
    real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xdot, ydot, zdot, u, v, b, midlon1, midlat1

    select case(velocity)
    case("nodiv")
      call uv_nodiv(t,longitudes,latitudes,gu,gv)
    case("div")
      call uv_div(t,longitudes,latitudes,gu,gv)
    end select
    call interpolate_setuv(gu, gv)

    call find_points(gu, gv, deltat, midlon, midlat, deplon, deplat)
    call find_nearest_grid

    do i = 1, nlon
      do j = 1, nlat
        call lonlat2xyz(longitudes(p(i, j)), latitudes(q(i , j)), xr, yr, zr)
        ! arrival points
        call lonlat2xyz(longitudes(i), latitudes(j), xg, yg, zg)

        b = 1.0d0 / sqrt( 2.0d0 * (1.0d0 + (xg*xr + yg*yr + zg*zr)) ) ! Ritchie1987 式(44)
        xm = b * (xg + xr)
        ym = b * (yg + yr)
        zm = b * (zg + zr)
        midlon1 = modulo(atan2(ym, xm) + pi2, pi2)
        midlat1 = asin(zm)

        xdot = (xg - xr) / (2.0d0 * deltat)
        ydot = (yg - yr) / (2.0d0 * deltat)
        zdot = (zg - zr) / (2.0d0 * deltat)
        call xyz2uv(xdot, ydot, zdot, midlon1, midlat1, u, v)  !Richie1987式(49)
        gum(i, j) = gum(i, j) + u
        gvm(i, j) = gvm(i, j) + v

        call interpolate_bilinearuv(midlon1, midlat1, u, v)
        gum(i, j) = gum(i, j) - u
        gvm(i, j) = gvm(i, j) - v
      enddo
    enddo

  end subroutine calculate_residual_velocity

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

  subroutine find_nearest_grid
    use math_module, only: math_pi, pi2=>math_pi2
    use sphere_module, only: xyz2uv, lonlat2xyz
    use grid_module, only: pole_regrid
    implicit none

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
      end do
    end do
  end subroutine  find_nearest_grid

end module nisl_2step_module