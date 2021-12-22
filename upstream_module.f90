module upstream_module
  ! finds departure and mid-points
    use grid_module, only: latitudes=>lat
    use time_module, only: imethoduv
    use interpolate_module, only: &
      interpolate_setuv, interpolate_bilinearuv, interpolate_polin2uv
    use sphere_module, only: lonlat2xyz, uv2xyz
    implicit none
    private
  
    integer(8), public :: itermax = 5
    real(8), public :: small = 1.0d-10
  
    public :: find_points, calc_niuv
  
  contains
  !! midlon, midlatはdt前の出発点(middle), deplon, deplatは2*dt前
    subroutine find_points(u, v, dt, midlon, midlat, deplon, deplat)
      use math_module, only: pi2=>math_pi2
      use uv_module, only: uv_sbody_calc
      implicit none
  
      real(8), dimension(:,:), intent(in) :: u, v
      real(8), intent(in) :: dt
      real(8), dimension(:,:), intent(inout) :: midlon, midlat, deplon, deplat
  
      integer(8) :: nx, ny, i, j, step
      real(8) :: un, vn, &     ! normalised velocity
                       bk, &          ! correction factor
                       xd, yd, zd, & ! Cartesian velocity
                       xg, yg, zg, & ! arival point in Cartesian coordinates
                       x0, y0, z0, & ! present point in Cartesian coordinates
                       x1, y1, z1, & ! updated point in Cartesian coordinates
                       lon, lat, err
  
      nx = size(u,1)
      ny = size(u,2)
  
      call interpolate_setuv(u,v)
      do j = 1, ny
        do i = 1, nx
          ! calculate initial values
          un = u(i,j)
          vn = v(i,j)
          lon = pi2*dble(i-1)/dble(nx) ! calculate (lon,lat) from (i,j)
          lat = latitudes(j)
          call lonlat2xyz(lon, lat, xg, yg, zg) ! transform into Cartesian coordinates
          ! r = g as an initial point for the 1st time step, 最初midlat, midlonには格子点上の値が入っている
          call lonlat2xyz(midlon(i,j), midlat(i,j), x0, y0, z0) 
          step = 1
          do 
            call uv2xyz(un,vn,lon,lat,xd,yd,zd) ! normalized Cartesian velocity
            ! correction factor
            bk = 1.0d0/sqrt(1.0d0+dt*dt*(xd*xd+yd*yd+zd*zd)-2.0d0*dt*(xd*xg+yd*yg+zd*zg))
            x1 =  bk*(xg - dt*xd) ! calculate new points
            y1 =  bk*(yg - dt*yd)
            z1 =  bk*(zg - dt*zd)
            ! calculate (lon,lat) from (x,y,z)
            lat = asin(z1)
            lon = modulo(atan2(y1,x1)+pi2,pi2)
            call uv_sbody_calc(lon, lat, un, vn)
            err = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0)) ! calculate error
            x0 = x1 ! save as the current point
            y0 = y1
            z0 = z1
            step = step + 1
            if ((err<small).or.(step>itermax)) then
              exit
            end if
          end do
          midlon(i,j) = lon ! store as the mid-point
          midlat(i,j) = lat
          bk = 2.0d0*(x0*xg+y0*yg+z0*zg) ! calculate the departure point
          x1 = bk*x0 - xg
          y1 = bk*y0 - yg
          z1 = bk*z0 - zg
          deplon(i,j) = modulo(atan2(y1,x1)+pi2,pi2)
          deplat(i,j) = asin(z1)
        end do
      end do
      
    end subroutine find_points

    subroutine calc_niuv(dt, p, q, lon, lat, midlon, midlat, gum, gvm)
      use grid_module, only: latitudes => lat, longitudes => lon
      use math_module, only: math_pi, pi2=>math_pi2
      use sphere_module, only: xyz2uv, lonlat2xyz
      use uv_module, only: uv_sbody_calc
      implicit none
      real(8), intent(in) :: dt
      integer(8), intent(in) :: p, q
      real(8), intent(in) :: lon, lat
      real(8), intent(out) :: gum, gvm
      real(8), intent(out) :: midlon, midlat

      real(8) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xdot, ydot, zdot, lon_grid, lat_grid, u, v, b


      lon_grid = longitudes(p)
      lat_grid = latitudes(q)
      call lonlat2xyz(lon_grid, lat_grid, xr, yr, zr)
      ! arrival points
      call lonlat2xyz(lon, lat, xg, yg, zg)

      b = 1.0d0 / sqrt( 2.0d0 * (1.0d0 + (xg*xr + yg*yr + zg*zr)) ) ! Ritchie1987 式(44)
      xm = b * (xg + xr)
      ym = b * (yg + yr)
      zm = b * (zg + zr)
      midlon = modulo(atan2(ym, xm) + pi2, pi2)
      midlat = asin(zm)

      xdot = (xg - xr) / dt
      ydot = (yg - yr) / dt
      zdot = (zg - zr) / dt
      call xyz2uv(xdot, ydot, zdot, midlon, midlat, u, v)  !Richie1987式(49)
      gum = gum + u
      gvm = gvm + v

      call uv_sbody_calc(midlon, midlat, u, v)
      gum = gum - u
      gvm = gvm - v


    end subroutine calc_niuv
  
  end module upstream_module