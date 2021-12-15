module nisl_module

  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, sphi_old, sphi, longitudes=>lon, latitudes=>lat, coslatr
  use field_module, only : X, Y
  private
  
  integer(8), allocatable, private :: p(:,:), q(:,:)
  real(8), dimension(:,:), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphi_initial, &
    midlon, midlat, deplon, deplat, gum, gvm
  complex(8), dimension(:,:), allocatable, private :: sphi1

  private :: update
  public :: nisl_init, nisl_timeint, nisl_clean

contains

  subroutine nisl_init()
    use time_module, only: deltat
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(8) :: i,j

    allocate(sphi1(0:ntrunc,0:ntrunc),gphi_old(nlon,nlat), &
             gphim(nlon,nlat),dgphi(nlon,nlat),dgphim(nlon,nlat), &
             midlon(nlon,nlat),midlat(nlon,nlat), gphi_initial(nlon, nlat), &
             deplon(nlon,nlat),deplat(nlon,nlat), p(nlon,nlat), q(nlon,nlat), &
             gum(nlon,nlat),gvm(nlon,nlat))
    call interpolate_init(gphi)

    call legendre_synthesis(sphi_old,gphi_old)
    gphi = gphi_old
    gphi_initial(:, :) = gphi(:, :)

    do i=1, nlon
      midlon(i,:) = longitudes(i)
    end do
    do j=1, nlat
      midlat(:,j) = latitudes(j)
    end do

    open(11, file="animation.txt")
    do i = 1, nlon
      do j = 1, nlat
          write(11,*) X(i, j), Y(i, j), gphi(i, j)
      end do        
    end do
    call update(0.25d0*deltat)
    print *, "step=1", " t=", real(deltat)
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

    do i=2, nstep
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
    
    open(10, file="log.txt")
    do i = 1, nlon
      do j = 1, nlat
        write(10,*) X(i, j), Y(i, j), gphi(i, j)
      enddo
    enddo
    close(10)

    open(12, file="error.txt")
    do j = 1, nlon
        do k = 1, nlat
          write(12,*) X(j, k), Y(j, k), gphi_initial(j, k) - gphi(j, k)
        end do
    end do
    close(12)
  end subroutine nisl_timeint

  subroutine update(dt)
    use time_module, only: imethod
    use uv_module, only: uv_nodiv, uv_div
    use upstream_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate_module, only: interpolate_set, interpolate_bilinear, interpolate_polin2
    implicit none

    integer(8) :: i, j, m
    real(8), intent(in) :: dt

    call find_points(gu, gv, 0.5d0*dt, midlon, midlat, deplon, deplat)
    call calc_niuv(dt)

    call legendre_synthesis(sphi_old, gphi_old)
    do j = 1, nlat
      do i = 1, nlon
        gphi(i,j) = gphi_old(p(i,j), q(i,j))
      end do
    end do

    ! dF/dlon
    call legendre_synthesis_dlon(sphi, dgphi)
    call interpolate_set(dgphi)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bilinear(midlon(i,j), midlat(i,j), dgphim(i,j))
      enddo
      gphim(:,j) = gum(:,j)*coslatr(j)*dgphim(:,j) ! gum: -u'
    enddo

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi) 
    call interpolate_set(dgphi)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bilinear(midlon(i,j), midlat(i,j), dgphim(i,j))
      enddo
      gphim(:,j) = gphim(:,j) + gvm(:,j)*coslatr(j)*dgphim(:,j) ! gvm: -v'
    enddo

    gphi = gphi + dt*gphim

! time filter
    call legendre_analysis(gphi, sphi1)
    do m = 0, ntrunc
      sphi_old(m : ntrunc, m) = sphi(m : ntrunc, m)       
      sphi(m : ntrunc, m) = sphi1(m : ntrunc, m)
    enddo

  end subroutine update

  subroutine calc_niuv(dt)
    use math_module, only: math_pi, pi2=>math_pi2
    use time_module, only: imethoduv
    use sphere_module, only: xyz2uv, lonlat2xyz, lonlat2uv
    use uv_module, only: uv_sbody_calc
    implicit none

    real(8), intent(in) :: dt

    integer(8) :: i,j
    real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xdot, ydot, zdot, lon, lat, lon_grid, lat_grid, u, v, dlonr, b

    dlonr = 0.5d0 * nlon / math_pi
    do j=1, nlat
      lat = latitudes(j)
      do i=1, nlon
! find grid points near departure points
        p(i,j) = anint(deplon(i,j)*dlonr+1.0d0)
        if (p(i,j)>nlon) then
          p(i,j) = p(i,j)-nlon
        end if
! lat = (J+1-2j)pi/(2J+1)
        q(i,j) = anint(0.5d0*(nlat+1-(2.0d0*nlat+1.0d0)*deplat(i,j)/math_pi))
        lon_grid = longitudes(p(i,j))
        lat_grid = latitudes(q(i,j))  
        call lonlat2xyz(lon_grid, lat_grid, xr, yr, zr)
! arrival points
        lon = longitudes(i)
        call lonlat2xyz(lon,lat,xg,yg,zg)
! calculate midpoints between integer departure points and arrival points
        b = 1.0d0/sqrt(2.0d0*(1.0d0+(xg*xr+yg*yr+zg*zr))) ! Ritchie1987 式(44)
        xm = b*(xg + xr)
        ym = b*(yg + yr)
        zm = b*(zg + zr)
        midlon(i,j) = modulo(atan2(ym,xm)+pi2,pi2)
        midlat(i,j) = asin(zm)
!       print *, real(lon*180/pi), real(midlon(i,j)*180/pi), real(lon_grid*180/pi), &
!         real(lat*180/pi), real(midlat(i,j)*180/pi), real(lat_grid*180/pi)
! calculate integer velocities at midpoints
        xdot = (xg - xr) / dt
        ydot = (yg - yr) / dt
        zdot = (zg - zr) / dt
        call xyz2uv(xdot, ydot, zdot, midlon(i,j), midlat(i,j), u, v)  !Richie1987式(49)
!       print *, real(lon_grid*180/pi), real(lat_grid*180/pi), real(u*a), real(v*a)
!       print *, real(zg), real(zr), real(zd), real(v*a)
        gum(i,j) = u
        gvm(i,j) = v
! calculate velocity at midpoints and -residual velocities
        call uv_sbody_calc(midlon(i, j), midlat(i, j), u, v)
        gum(i,j) = gum(i,j) - u
        gvm(i,j) = gvm(i,j) - v
      end do
    end do
        
  end subroutine  calc_niuv

end module nisl_module