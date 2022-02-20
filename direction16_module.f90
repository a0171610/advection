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
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i,j

    allocate(sphi1(0:ntrunc, 0:ntrunc),gphi_old(nlon, nlat))
    allocate(gphim(nlon, nlat),dgphi(nlon, nlat))
    allocate(midlon(4, nlon, nlat), midlat(4, nlon, nlat))
    allocate(deplon(nlon, nlat), deplat(nlon, nlat), p(4, nlon, nlat), q(4, nlon, nlat))
    allocate(weight(4, nlon, nlat))
    allocate(gum(4, nlon, nlat), gvm(4, nlon, nlat), dgphim(4, nlon, nlat))

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
          write(11,*) X(i, j), Y(i, j), gphi(i, j)
      end do        
    end do
    call update(deltat)
    write(*, *) 'step = 1 ', "maxval = ", maxval(gphi), 'minval = ', minval(gphi)

  end subroutine direction16_init

  subroutine direction16_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(sphi1, gphi_old, gphim, dgphi, deplon, deplat)
    deallocate(midlon, midlat, gum, gvm)
    deallocate(weight, p, q)
    call interpolate_clean()

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
    use interpolate_module, only: interpolate_set, interpolate_setd, find_stencil_
    use interpolate_module, only: interpolate_bicubic, interpolate_bilinear, interpolate_bilinear_ratio
    use interpolate_module, only: interpolate_dist, interpolate_dist_ratio
    implicit none

    integer(8) :: i, j, m
    real(8), intent(in) :: dt

    call find_points(gu, gv, 0.5*dt, midlon(1,:,:), midlat(1,:,:), deplon, deplat)
    ! dtに0.5をかけているのは引数のdtが最初のステップ以外は2.0*deltatを渡しているから

    do i = 1, nlon
        do j = 1, nlat
          call interpolate_bilinear_ratio(deplon(i, j), deplat(i, j),&
           weight(1, i, j), weight(2, i, j), weight(3, i, j), weight(4, i, j))
        end do
    end do

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
        call interpolate_bicubic(midlon(1, i, j), midlat(1, i, j), dgphim(1, i, j))
        call interpolate_bicubic(midlon(2, i, j), midlat(2, i, j), dgphim(2, i, j))
        call interpolate_bicubic(midlon(3, i, j), midlat(3, i, j), dgphim(3, i, j))
        call interpolate_bicubic(midlon(4, i, j), midlat(4, i, j), dgphim(4, i, j))

        gphim(i, j) = gphim(i, j) + weight(1, i, j) * gum(1, i, j) * dgphim(1, i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + weight(2, i, j) * gum(2, i, j) * dgphim(2, i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + weight(3, i, j) * gum(3, i, j) * dgphim(3, i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + weight(4, i, j) * gum(4, i, j) * dgphim(4, i, j) / cos(latitudes(j))
      enddo
    enddo

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi) 
    call bicubic_interpolation_set(dgphi)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(midlon(1, i, j), midlat(1, i, j), dgphim(1, i, j))
        call interpolate_bicubic(midlon(2, i, j), midlat(2, i, j), dgphim(2, i, j))
        call interpolate_bicubic(midlon(3, i, j), midlat(3, i, j), dgphim(3, i, j))
        call interpolate_bicubic(midlon(4, i, j), midlat(4, i, j), dgphim(4, i, j))

        gphim(i, j) = gphim(i, j) + weight(1, i, j) * gvm(1, i, j) * dgphim(1, i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + weight(2, i, j) * gvm(2, i, j) * dgphim(2, i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + weight(3, i, j) * gvm(3, i, j) * dgphim(3, i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + weight(4, i, j) * gvm(4, i, j) * dgphim(4, i, j) / cos(latitudes(j))
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
    use interpolate_module, only: find_stencil_
    use upstream_module, only: calc_niuv
    implicit none

    real(8), intent(in) :: dt

    integer(8) :: i, j, k
    real(8) :: dlonr
    integer(8), dimension(4) :: tmp1, tmp2

    dlonr = 0.5d0 * nlon / math_pi
    call set_zero()
    do j = 1, nlat
      do i = 1, nlon
        ! find grid points near departure points
        call find_stencil_(deplon(i, j), deplat(i, j), tmp1, tmp2)
        do k = 1, 4
            p(k, i, j) = tmp1(k)
            q(k, i, j) = tmp2(k)
        end do

        call calc_niuv(dt, p(1,i,j), q(1,i,j), longitudes(i), latitudes(j), midlon(1, i,j), midlat(1, i, j), gum(1,i,j), gvm(1,i,j))
        call calc_niuv(dt, p(2,i,j), q(2,i,j), longitudes(i), latitudes(j), midlon(2, i,j), midlat(2, i, j), gum(2,i,j), gvm(2,i,j))
        call calc_niuv(dt, p(3,i,j), q(3,i,j), longitudes(i), latitudes(j), midlon(3, i,j), midlat(3, i, j), gum(3,i,j), gvm(3,i,j))
        call calc_niuv(dt, p(4,i,j), q(4,i,j), longitudes(i), latitudes(j), midlon(4, i,j), midlat(4, i, j), gum(4,i,j), gvm(4,i,j))
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
    call legendre_synthesis_dlat(sphi1, gphix)
    call legendre_synthesis_dlat(sphi1, gphiy)
    call legendre_synthesis_dlonlat(sphi1, gphixy)
    do j = 1, nlat
      gphiy(:,j) = gphiy(:,j) / cos(latitudes(j))
      gphixy(:,j) = gphixy(:,j) / cos(latitudes(j))
    end do

  end subroutine bicubic_interpolation_set

  subroutine set_zero
    implicit none
    gum(:, :, :) = 0.0d0
    gvm(:, :, :) = 0.0d0
  end subroutine set_zero

end module direction16_module