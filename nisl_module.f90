module nisl_module

  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat
  use field_module, only : X, Y
  private
  
  integer(8), allocatable, private :: p(:,:), q(:,:)
  real(8), dimension(:,:), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphix, gphiy, gphixy, &
    midlon, midlat, deplon, deplat, gum, gvm
  complex(8), dimension(:,:), allocatable, private :: sphi1

  private :: update, bicubic_interpolation_set
  public :: nisl_init, nisl_timeint, nisl_clean

contains

  subroutine nisl_init()
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
          write(11,*) X(i, j), Y(i, j), gphi(i, j)
      end do        
    end do
    call update(deltat)

  end subroutine nisl_init

  subroutine nisl_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(sphi1,gphi_old,gphim,dgphi,dgphim,gum,gvm, &
      midlon,midlat,deplon,deplat,p,q)
    call interpolate_clean()

  end subroutine nisl_clean

  subroutine nisl_timeint()
    use time_module, only: nstep, deltat, hstep
    use legendre_transform_module, only: legendre_synthesis
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
    
  end subroutine nisl_timeint

  subroutine update(dt)
    use uv_module, only: uv_nodiv, uv_div
    use upstream_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate_module, only: interpolate_set, interpolate_bilinear, interpolate_setd, interpolate_bicubic
    implicit none

    integer(8) :: i, j, m
    real(8), intent(in) :: dt

    call find_points(gu, gv, 0.5d0*dt, midlon, midlat, deplon, deplat)   ! dtに0.5をかけているのは引数のdtが最初のステップ以外は2.0*deltatを渡しているから
    call set_niuv(dt)

    call legendre_synthesis(sphi_old, gphi_old)

    do j = 1, nlat
      do i = 1, nlon
        gphi(i,j) = gphi_old(p(i, j), q(i, j))
      end do
    end do

    ! dF/dlon
    call legendre_synthesis_dlon(sphi, dgphi)
    call bicubic_interpolation_set(dgphi) 
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(midlon(i, j), midlat(i, j), dgphim(i, j))
      enddo
      gphim(:, j) = gum(:, j) * dgphim(:,j) / cos(latitudes(j)) ! gum: -u'
    enddo

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi) 
    call bicubic_interpolation_set(dgphi)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(midlon(i, j), midlat(i, j), dgphim(i, j))
      enddo
      gphim(:, j) = gphim(:, j) + gvm(:, j) * dgphim(:,j) / cos(latitudes(j)) ! gvm: -v'
    enddo

    gphi = gphi + dt * gphim

! time filter
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
    use upstream_module, only: calc_niuv
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

        p(i, j) = anint( deplon(i, j) * dlonr + 1.0d0 )
        if ( p(i,j) > nlon ) then
          p(i, j) = p(i, j) - nlon
        end if
        ! lat = (J+1-2j)pi/(2J+1)
        q(i, j) = anint( 0.5d0 * (nlat + 1.0d0 - (2.0d0*dble(nlat)+1.0d0)*deplat(i, j) / math_pi) )  !latitudesは大きい順で詰められているので注意
        call calc_niuv(dt, p(i, j), q(i, j), longitudes(i), latitudes(j), midlon(i, j), midlat(i, j), gum(i, j), gvm(i, j))
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

end module nisl_module