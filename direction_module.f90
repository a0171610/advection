module direction_module

  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, gphi_initial, sphi_old, sphi, longitudes=>lon, latitudes=>lat
  use field_module, only : X, Y
  private
  
  integer(8), dimension(:, :), allocatable, private :: pa, qa, pb, qb, pc, qc, pd, qd
  real(8), allocatable, private :: A(:, :), B(:, :), C(:, :), D(:, :) 
  real(8), dimension(:,:), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim, gphix, gphiy, gphixy, deplon, deplat
  real(8), dimension(:, :), allocatable, private :: midlonA, midlatA, midlonB, midlatB, midlonC, midlatC, midlonD, midlatD
  real(8), dimension(:, :), allocatable, private :: guma, gvma, gumb, gvmb, gumc, gvmc, gumd, gvmd
  complex(8), dimension(:,:), allocatable, private :: sphi1
  integer(8), dimension(:, :, :), allocatable,  private :: is, js

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
    allocate(gphim(nlon, nlat),dgphi(nlon, nlat),dgphim(nlon, nlat))
    allocate(midlonA(nlon, nlat), midlatA(nlon, nlat), midlonB(nlon, nlat), midlatB(nlon, nlat))
    allocate(midlonC(nlon, nlat), midlatC(nlon, nlat), midlonD(nlon, nlat), midlatD(nlon, nlat))
    allocate(deplon(nlon, nlat), deplat(nlon, nlat), pa(nlon, nlat), qa(nlon, nlat), pb(nlon, nlat), qb(nlon, nlat))
    allocate(pc(nlon, nlat), qc(nlon, nlat), pd(nlon, nlat), qd(nlon, nlat))
    allocate(A(nlon, nlat), B(nlon, nlat), C(nlon, nlat), D(nlon, nlat))
    allocate(guma(nlon, nlat), gvma(nlon, nlat), gumb(nlon, nlat), gvmb(nlon, nlat))
    allocate(gumc(nlon, nlat), gvmc(nlon, nlat), gumd(nlon, nlat), gvmd(nlon, nlat))
    allocate(gphix(nlon, nlat), gphiy(nlon, nlat), gphixy(nlon, nlat))
    allocate(is(nlon, nlat, 4), js(nlon, nlat, 4))

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

  end subroutine direction_init

  subroutine direction_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(sphi1, gphi_old, gphim, dgphi, dgphim, is, js, deplon, deplat)
    deallocate(midlonA, midlatA, midlonB, midlatB, midlonC, midlatC, midlonD, midlatD)
    deallocate(guma, gvma, gumb, gvmb, gumc, gvmc, gumd, gvmd)
    deallocate(A, B, C, D, pa, qa, pb, qb, pc, qc, pd, qd)
    call interpolate_clean()

  end subroutine direction_clean

  subroutine direction_timeint()
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
    
  end subroutine direction_timeint

  subroutine update(dt)
    use upstream_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat
    use interpolate_module, only: interpolate_set, interpolate_setd, find_stencil_
    use interpolate_module, only: interpolate_bicubic, interpolate_bilinear_ratio, interpolate_bilinear
    implicit none

    integer(8) :: i, j, m
    real(8), intent(in) :: dt

    call find_points(gu, gv, 0.5d0*dt, deplon=deplon, deplat=deplat)
    ! dtに0.5をかけているのは引数のdtが最初のステップ以外は2.0*deltatを渡しているから
    
    do i = 1, nlon
        do j = 1, nlat
            call find_stencil_(deplon(i, j), deplat(i, j), is(i, j, :), js(i, j, :))
            call interpolate_bilinear_ratio(deplon(i, j), deplat(i, j), A(i, j), B(i, j), C(i, j), D(i, j))
        end do
    end do
    call set_niuv(dt)

    call legendre_synthesis(sphi_old, gphi_old)

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
        call interpolate_bicubic(midlon(i, j), midlat(i, j), dgphim(i, j))
        gphim(i, j) = gphim(i, j) + A(i, j) * gumm(is(i,j,1), js(i,j,1)) * dgphim(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + B(i, j) * gumm(is(i,j,2), js(i,j,2)) * dgphim(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + C(i, j) * gumm(is(i,j,3), js(i,j,3)) * dgphim(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + D(i, j) * gumm(is(i,j,4), js(i,j,4)) * dgphim(i, j) / cos(latitudes(j))
      enddo
    enddo

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi) 
    call bicubic_interpolation_set(dgphi)
    call interpolate_set(dgphi)
    call interpolate_setd(gphix, gphiy, gphixy)
    do j = 1, nlat
      do i = 1, nlon
        call interpolate_bicubic(midlon(i, j), midlat(i, j), dgphim(i, j))
        gphim(i, j) = gphim(i, j) + A(i, j) * gvmm(is(i,j,1), js(i,j,1)) * dgphim(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + B(i, j) * gvmm(is(i,j,2), js(i,j,2)) * dgphim(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + C(i, j) * gvmm(is(i,j,3), js(i,j,3)) * dgphim(i, j) / cos(latitudes(j))
        gphim(i, j) = gphim(i, j) + D(i, j) * gvmm(is(i,j,4), js(i,j,4)) * dgphim(i, j) / cos(latitudes(j))
      enddo
    enddo

    gphi(:, :) = gphi(:, :) + dt * gphim(:, :)

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

  subroutine regrid(f, ff)
    implicit none

    real(8), dimension(:, :), intent(in) :: f
    real(8), dimension(0:, -1:), intent(out) :: ff

    integer(8) :: i, j
    integer(8) :: nx, ny

    nx = size(f, 1)
    ny = size(f, 2)

    ff(1:nx,1:ny) = f
    do j = 1, 2
      ff(1:nx,1-j) = cshift(ff(1:nx,j),nx/2)
      ff(1:nx,ny+j) = cshift(ff(1:nx,ny-(j-1)),nx/2)
    end do
    do i = 1, 1
      ff(1-i,:) = ff(nx-(i-1),:)
    end do
    do i = 1, 2
      ff(nx+i,:) = ff(1+(i-1),:)
    end do

  end subroutine regrid

end module direction_module