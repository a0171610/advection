module analysis_module
  use grid_module, only: gphi, gphi_initial, lat, wgt, lon, ntrunc, Umax
  use field_module, only: X, Y
  use legendre_transform_module, only: legendre_analysis
  use math_module, only: math_pi
  use planet_module, only: planet_radius
  use time_module, only: deltat
  implicit none

contains

  subroutine error_log()
    implicit none
    integer(8) :: i, j, nlon, nlat
    real(8), allocatable :: w(:, :)
    real(8), allocatable :: l1(:, :), l2(:, :)
    complex(8), allocatable :: l1_t(:, :), l2_t(:, :)
    real(8) :: dq, dqp, rmse
    real(8) :: sum_g1, sum_g2
    real(8) :: d_lamda

    nlat = size(gphi, 2)
    nlon = size(gphi, 1)
    allocate(w(nlon, nlat))

    do j=1, nlat
      w(:, j) = wgt(j)
    end do

    open(10, file="log.txt")
    do i = 1, nlon
      do j = 1, nlat
        write(10,*) lon(i), lat(j), gphi(i, j)
      enddo
    enddo
    close(10)

    open(12, file="error.txt")
    do i = 1, nlon
        do j = 1, nlat
          write(12,*) lon(i), lat(j), gphi_initial(i, j) - gphi(i, j)
        end do
    end do
    close(12)

    open(14, file = "error_equator.txt")
    do i = 1, nlat
      write(14,*) lat(i), gphi(1, i) - gphi_initial(1, i)
    end do

    dq = sum(gphi(:, :) - gphi_initial(:, :)) / dble(nlat * nlon)

    dqp = 0.0d0
    do i = 1, nlon
        do j = 1, nlat
            if(gphi(i, j) > gphi_initial(i, j)) then
              dqp = dqp + gphi(i, j) - gphi_initial(i, j)
              dqp = dqp / dble(nlon * nlat)
            endif
        end do
    end do

    rmse = 0.0d0
    do i = 1, nlon
      do j = 1, nlat
        rmse = rmse + (gphi(i, j) - gphi_initial(i, j)) ** 2
        rmse = rmse / dble(nlon * nlat)
      end do
    end do
    rmse = sqrt(rmse)

    write(*,*) 'error = ', dq, 'positive error = ', dqp, 'RMSE = ', rmse

    sum_g1 = 0.0d0
    sum_g2 = 0.0d0

    do i = 1, nlon
        do j = 1, nlat
            sum_g1 = sum_g1 + gphi_initial(i, j) * w(i, j)
            sum_g2 = sum_g2 + gphi(i, j) * w(i, j)
        end do
    end do
    write(*,*) "initial global mass sum", sum_g1, "final global mass sum", sum_g2

    ! l2ノルムを求める
    allocate(l1(nlon, nlat), l2(nlon, nlat))
    allocate(l1_t(0:ntrunc,0:ntrunc), l2_t(0:ntrunc,0:ntrunc))
    do i = 1, nlon
      do j = 1, nlat
        l1(i, j) = ((gphi(i, j) - gphi_initial(i, j)) ** 2)
        l2(i, j) = (gphi_initial(i, j) ** 2)
      end do
    end do

    sum_g1 = 0.0d0
    sum_g2 = 0.0d0
    do i = 1, nlon
      do j = 1, nlat
        sum_g1 = sum_g1 + l1(i, j) * wgt(j)
        sum_g2 = sum_g2 + l2(i, j) * wgt(j)
      end do
    end do

    write(*,*) "l2 norm = ", sqrt(sum_g1 / sum_g2)
    write(*,*) 'l1', sum(l1), 'l2', sum(l2)
    write(*,*) sum_g1, sum_g2

    d_lamda = 360.0d0 / dble(nlon)

    write(*,*) 'Courant Number = ', Umax * deltat / (d_lamda * math_pi * planet_radius / 180.0d0)

  end subroutine error_log
end module analysis_module