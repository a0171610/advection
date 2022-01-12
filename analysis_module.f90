module analysis_module
  use grid_module, only: gphi, gphi_initial, lat
  use field_module, only: X, Y
  implicit none

contains

  subroutine error_log()
    implicit none
    integer(8) :: i, j, nlon, nlat
    real(8) :: dq, dqp

    nlat = size(gphi, 2)
    nlon = size(gphi, 1)

    open(10, file="log.txt")
    do i = 1, nlon
      do j = 1, nlat
        write(10,*) X(i, j), Y(i, j), gphi(i, j)
      enddo
    enddo
    close(10)

    open(12, file="error.txt")
    do i = 1, nlon
        do j = 1, nlat
          write(12,*) X(i, j), Y(i, j), gphi_initial(i, j) - gphi(i, j)
        end do
    end do
    close(12)

    open(14, file = "error_equator.txt")
    do i = 1, nlat
      write(14,*) lat(i), gphi(1, i) - gphi_initial(1, i)
    end do

    dq = sum(gphi(:, :) - gphi_initial(:, :))

    dqp = 0.0d0
    do i = 1, nlon
        do j = 1, nlat
            if(gphi(i, j) > gphi_initial(i, j)) then
              dqp = dqp + gphi(i, j) - gphi_initial(i, j)
            endif
        end do
    end do
    write(*,*) '△q = ', dq, '△q+', dqp
  end subroutine error_log
end module analysis_module