module interpolate16_module
  ! interpolate in a stencil
    use grid_module, only: latitudes=>lat, wgt, longitudes=>lon
    use sphere_module, only: lon2i, lat2j
    implicit none
    private
  
    integer(8), private :: nx, ny, n = 3, nh, nhalo, nx1, nx2, ny1, ny2 
    real(8), private :: u, t, dlon                            ! tとuは線分比
    real(8), dimension(:), allocatable, public :: lon_extend, lat_extend
    real(8), dimension(:,:), allocatable, private :: ff, ffx, ffy, ffxy, fu, fv, ffxl, ffyl
    integer(8), dimension(4), private :: is_, js_

    public :: interpolate16_init, interpolate16_clean, interpolate16_set, find_stencil_16, &
              interpolate16_dist, interpolate16_dist_ratio, interpolate16_setuv, &
              interpolate16_bicubic, interpolate16_setd, interpolate16_bilinearuv, &
              interpolate12_dist_ratio
  
  contains
  
    subroutine interpolate16_init(f)
      use math_module, only: pi2=>math_pi2, pih=>math_pih
      implicit none
  
      real(8), dimension(:,:) :: f
  
      integer(8) :: i, j
  
      nh = n/2
      nhalo = nh
    
      nx = size(f,1)
      ny = size(f,2)
      nx1 = 1 - nhalo
      nx2 = nx + nhalo + 1
      ny1 = 1 - nhalo - 1
      ny2 = ny + nhalo + 1
  
      allocate(lon_extend(nx1:nx2), lat_extend(ny1:ny2), ff(nx1:nx2,ny1:ny2), &
               ffx(nx1:nx2,ny1:ny2), ffy(nx1:nx2,ny1:ny2), ffxy(nx1:nx2,ny1:ny2), &
               ffxl(nx1:nx2,ny1:ny2), ffyl(nx1:nx2,ny1:ny2), &
               fu(nx1:nx2,ny1:ny2),fv(nx1:nx2,ny1:ny2))
  
      dlon = pi2/nx
      do i = nx1, nx2
        lon_extend(i) = dlon * (i-1)
      end do
      lat_extend(1 : ny) = latitudes
      do j = 1, nhalo + 1
        lat_extend(1 - j)   = pih + (pih - latitudes(j))
        lat_extend(ny + j)   = -pih + (-pih - latitudes(ny-j+1))
      end do
  
    end subroutine interpolate16_init
  
    subroutine interpolate16_clean()
      implicit none
  
      deallocate(lon_extend, lat_extend, ff, ffx, ffy, ffxy, fu, fv, ffxl, ffyl)
  
    end subroutine  interpolate16_clean

    subroutine interpolate16_dist_ratio(lon, lat, weight)
      use sphere_module, only: orthodrome
      use grid_module, only: pole_regrid
      implicit none
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: weight(16)
      integer(8) lon_g(16), lat_g(16)
      real(8) :: dist(16), dist_sum
      integer(8) :: i

      call find_stencil_16(lon, lat, lon_g, lat_g)

      do i = 1, 16
          dist(i) = 1.0d0 / (orthodrome(lon, lat, longitudes(lon_g(i)), latitudes(lat_g(i))) ** 5)
      end do
      dist_sum = sum(dist(:))

      do i = 1, 16
          weight(i) = dist(i) / dist_sum
      end do

    end subroutine interpolate16_dist_ratio

    subroutine interpolate12_dist_ratio(lon, lat, weight)
      use sphere_module, only: orthodrome
      use grid_module, only: pole_regrid
      implicit none
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: weight(16)
      integer(8) lon_g(16), lat_g(16)
      real(8) :: dist(16), dist_sum
      integer(8) :: i

      call find_stencil_16(lon, lat, lon_g, lat_g)

      do i = 1, 16
          dist(i) = 1.0d0 / (orthodrome(lon, lat, longitudes(lon_g(i)), latitudes(lat_g(i))) ** 5)
      end do
      !dist(1) = 0.0d0; dist(4) = 0.0d0
      !dist(13) = 0.0d0; dist(16) = 0.0d0
      call change_array_min_N_to_zero(12, dist)
      dist_sum = sum(dist(:))

      do i = 1, 16
          weight(i) = dist(i) / dist_sum
      end do

    end subroutine interpolate12_dist_ratio

    subroutine interpolate16_dist(lon, lat, ans)
      implicit none

      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: ans
      real(8) :: weight(16)
      real(8), dimension(16) :: fs
      integer(8), dimension(16) :: is, js
      integer(8) :: i

      call find_stencil_16(lon, lat, is, js)
      do i = 1, 16
        fs(i) = ff(is(i), js(i))
      end do

      call interpolate16_dist_ratio(lon, lat, weight)

      ans = 0.0d0
      do i = 1, 16
          ans = ans + fs(i) * weight(i)
      end do
    end subroutine interpolate16_dist
  
    subroutine interpolate16_set(f)
      implicit none
  
      real(8), dimension(:, :), intent(in) :: f
  
      integer(8) :: i, j
  
      ff(1:nx, 1:ny) = f
      do j = 1, nh+1
        ff(1:nx, 1-j) = cshift(ff(1:nx,j), nx/2)
        ff(1:nx, ny+j) = cshift(ff(1:nx,ny-(j-1)), nx/2)
      end do
      do i = 1, nh
        ff(1-i, :) = ff(nx-(i-1), :)
      end do
      do i = 1, nh+1
        ff(nx+i, :) = ff(1+(i-1), :)
      end do
  
    end subroutine interpolate16_set

    subroutine interpolate16_bilinearuv(lon, lat, fiu, fiv)
      implicit none

      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: fiu, fiv
      real(8), dimension(4) :: fsu, fsv
      integer(8) :: k

      call find_stencil(lon, lat)
      do k=1, 4
        fsu(k) = fu(is_(k),js_(k))
        fsv(k) = fv(is_(k),js_(k))
      end do
      fiu = (1.0d0-u)*((1.0d0-t)*fsu(1)+t*fsu(2)) + u*(t*fsu(3)+(1.0d0-t)*fsu(4))
      fiv = (1.0d0-u)*((1.0d0-t)*fsv(1)+t*fsv(2)) + u*(t*fsv(3)+(1.0d0-t)*fsv(4))

    end subroutine interpolate16_bilinearuv

    subroutine find_stencil(lon, lat)
      implicit none

      real(8), intent(in) :: lon, lat

      integer(8) :: j

      is_(1) = lon2i(lon,nx)
      is_(2) = is_(1) + 1
      t = lon/dlon - is_(1) + 1.0d0 ! t = (lon - dlon*(i-1))/dlon
      is_(3:4) = is_(2:1:-1)

      j = lat2j(lat, ny)
      if (lat > lat_extend(j)) then
        j = j - 1
      end if
      js_(1:2) = j
      js_(3:4) = j + 1
      u = (lat-lat_extend(j))/(lat_extend(j+1)-lat_extend(j))

    end subroutine find_stencil
  
    subroutine find_stencil_4(lon, lat, is, js)
      use grid_module, only: pole_regrid
      implicit none
      real(8), intent(in) :: lon, lat
      integer(8), dimension(:), intent(out) :: is, js
   
      integer(8) :: j
  
      is(1) = lon2i(lon, nx)
      is(2) = is(1) + 1
      t = lon/dlon - is(1) + 1.0d0 ! t = (lon - dlon*(i-1))/dlon
      is(3:4) = is(2:1:-1)
  
      j = lat2j(lat, ny)
      if (lat > lat_extend(j)) then
        j = j - 1
      end if
      js(1 : 2) = j
      js(3 : 4) = j + 1
      u = (lat - lat_extend(j)) / (lat_extend(j+1) - lat_extend(j))

    end subroutine find_stencil_4

    subroutine find_stencil_16(lon, lat, is, js)
      use grid_module, only: pole_regrid
      implicit none
      real(8), intent(in) :: lon, lat
      integer(8), intent(out) :: is(16), js(16)
      integer(8) :: tmpi(4), tmpj(4)
      integer(8) :: i

      call  find_stencil_4(lon, lat, tmpi, tmpj)

      is(7) = tmpi(1); js(7) = tmpj(1)
      is(6) = tmpi(2); js(6) = tmpj(2)
      is(11) = tmpi(3); js(11) = tmpj(3)
      is(10) = tmpi(4); js(10) = tmpj(4)

      is(1) = is(7) - 1; js(1) = js(7) - 1
      is(2) = is(7); js(2) = js(7) - 1
      is(3) = is(6); js(3) = js(7) - 1
      is(4) = is(6) + 1; js(4) = js(7) - 1
      is(5) = is(6) + 1; js(5) = js(6)
      is(8) = is(7) - 1; js(8) = js(6)
      is(9) = is(7) - 1; js(9) = js(10)
      is(12) = is(6) + 1; js(12) = js(10)
      is(13) = is(6) + 1; js(13) = js(10) + 1
      is(14) = is(6); js(14) = js(10) + 1
      is(15) = is(7); js(15) = js(10) + 1
      is(16) = is(7) - 1; js(16) = js(10) + 1

      do i = 1, 16
        call pole_regrid(is(i), js(i))
      end do

    end subroutine find_stencil_16

    subroutine interpolate16_setuv(gu,gv)
      implicit none
  
      real(8), dimension(:,:), intent(in) :: gu, gv
  
      integer(8) :: i, j
  
      do j=1, ny
        fu(1:nx,j) = gu(:,j)
        fv(1:nx,j) = gv(:,j)
      end do
  ! direction of u, v is reversed beyond poles
      do j=1, nh+1
        fu(1:nx,1-j) = -cshift(fu(1:nx,j),nx/2)
        fu(1:nx,ny+j) = -cshift(fu(1:nx,ny-(j-1)),nx/2)
        fv(1:nx,1-j) = -cshift(fv(1:nx,j),nx/2)
        fv(1:nx,ny+j) = -cshift(fv(1:nx,ny-(j-1)),nx/2)
      end do
      do i=1, nh
        fu(1-i,:) = fu(nx-(i-1),:)
        fv(1-i,:) = fv(nx-(i-1),:)
      end do
      do i=1, nh+1
        fu(nx+i,:) = fu(1+(i-1),:)
        fv(nx+i,:) = fv(1+(i-1),:)
      end do
  
    end subroutine interpolate16_setuv

    subroutine interpolate16_bicubic(lon, lat, fi)
      use bicubic_module, only: bcucof, bcuint, bcuintp
      implicit none
  
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: fi
  
      real(8), dimension(4) :: z, zx, zy, zxy
      integer(8), dimension(4) :: xs, ys
      integer(8) :: k
      real(8) :: dlat
  
      call find_stencil_4(lon, lat, xs, ys)

      dlat = lat_extend(ys(4)) - lat_extend(ys(1))
      do k = 1, 4
        z(k) = ff(xs(k), ys(k))
        zx(k) = ffx(xs(k), ys(k))
        zy(k) = ffy(xs(k), ys(k))
        zxy(k) = ffxy(xs(k), ys(k))
      end do
  
      call bcucof(z, zx, zy, zxy, dlon, dlat)
      fi = bcuint(t, u)
  
    end subroutine interpolate16_bicubic

    subroutine interpolate16_setd(fx,fy,fxy)
      implicit none
  
      real(8), dimension(:, :), intent(in) :: fx, fy, fxy
  
      integer(8) :: i, j
  
      ffx(1:nx, 1:ny) = fx
      ffy(1:nx, 1:ny) = fy
      ffxy(1:nx, 1:ny) = fxy
  ! directions of d/dx and d/dy are reversed beyond poles
      do j=1, nh+1
        ffx(1:nx, 1-j) = -cshift(ffx(1:nx,j), nx/2)
        ffx(1:nx, ny+j) = -cshift(ffx(1:nx,ny-(j-1)), nx/2)
        ffy(1:nx, 1-j) = -cshift(ffy(1:nx,j), nx/2)
        ffy(1:nx, ny+j) = -cshift(ffy(1:nx,ny-(j-1)), nx/2)
        ffxy(1:nx, 1-j) = cshift(ffxy(1:nx, j), nx/2)
        ffxy(1:nx, ny+j) = cshift(ffxy(1:nx, ny-(j-1)), nx/2)
      end do
      do i=1, nh
        ffx(1-i, :) = ffx(nx-(i-1), :)
        ffy(1-i, :) = ffy(nx-(i-1), :)
        ffxy(1-i, :) = ffxy(nx-(i-1), :)
      end do
      do i=1, nh+1
        ffx(nx+i, :) = ffx(1+(i-1), :)
        ffy(nx+i, :) = ffy(1+(i-1), :)
        ffxy(nx+i, :) = ffxy(1+(i-1), :)
      end do
  
    end subroutine interpolate16_setd

    subroutine change_array_min_N_to_zero(id, arr)
      use sort_module, only: bubblesort
      implicit none
      integer, intent(in) :: id
      real(8), intent(inout) :: arr(16)
      real(8) :: arr1(16)
      real(8) :: limit
      integer(8) :: i
      real(8), parameter :: eps = 0.0000001d0

      arr1(:) = arr(:)
      call bubblesort(16, arr1)
      limit = arr1(id)

      arr1(:) = arr(:)
      do i = 1, 16
        if (arr1(i) < limit .or. abs(arr1(i) - limit) < eps) then
          arr1(i) = 0.0d0
        endif
      enddo
      arr(:) = arr1(:)

    end subroutine change_array_min_N_to_zero

  end module interpolate16_module