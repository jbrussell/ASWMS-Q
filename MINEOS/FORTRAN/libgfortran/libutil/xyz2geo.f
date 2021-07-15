      subroutine xyz2geo(iflag,x,lat,lon,hght)

*     Convert Cartesian to geodetic coordinates assuming a given
*     reference ellipsoid (finv,semi) and geocenter offset t(3)
*       iflag = 1  xyz to geodetic
*             = 2  geodetic to xyz

*       lat and lon in degrees
*       hght in meters

      integer*4 i,iflag
      real*8    lat,lon,hght,x(3),semi,finv,t(3)
      real*8    f,e2,sinlat,coslat,sinlon,coslon,curvn
      real*8    sqr,lat0,cutoff
      real*8    pi,rad_to_deg

*     Numerical constants
      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( rad_to_deg    = 180.d0   /pi         )

*     GRS80 ellipsoid used by NAD83 and WGS84
      data semi,finv,(t(i),i=1,3)/6378137.,298.257222101,0.,0.,0./

      f= 1.d0/finv
      e2= 2.d0*f - f*f

      if( iflag.eq.1) then
*        xyz to geodetic: requires iterations on latitude
         do i = 1, 3
            x(i) = x(i) - t(i)
         enddo

         if (x(1)+1.d0.eq.1.d0 .and. x(2)+1.d0.eq.1.d0) then
            lon = 0.d0
         else
            lon= datan2(x(2),x(1))
         endif

         if( lon.lt.0d0 ) lon=lon + 2.d0*pi

*        starting value for latitude iteration
         sqr= dsqrt(x(1)**2 + x(2)**2)
         lat0= datan2(x(3)/sqr,1.d0-e2)
         lat= lat0

   40    sinlat= dsin(lat)
         curvn= semi/(dsqrt(1.d0-e2*sinlat*sinlat))
         lat= datan2((x(3)+e2*curvn*sinlat),sqr)
*        iterate to millimeter level
         if( dabs(lat-lat0).lt.1.d-10) goto 30
         lat0= lat
         goto 40

   30    continue
         cutoff= 80.d0/rad_to_deg
         if(lat.le.cutoff) then
            hght= (sqr/dcos(lat))-curvn
         else
            hght= z/dsin(lat)-curvn+e2*curvn
         endif
         lat=lat*rad_to_deg
         lon=lon*rad_to_deg
         if(lon.gt.180.d0) lon=lon-360.d0

      else
*        geodetic to xyz
         lat = lat/rad_to_deg
         lon = lon/rad_to_deg
         sinlat= dsin(lat)
         coslat= dcos(lat)
         sinlon= dsin(lon)
         coslon= dcos(lon)
         curvn= semi/(dsqrt(1.d0-e2*sinlat*sinlat))

         x(1) = (curvn+hght)*coslat*coslon + t(1)
         x(2) = (curvn+hght)*coslat*sinlon + t(2)
         x(3) = (curvn*(1.d0-e2)+hght)*sinlat + t(3)

      endif

      return
      end


