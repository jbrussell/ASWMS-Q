c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine convert_geoc(lat,lon,glat,glon)
c-----------------------------------------------------------------------
c converts geographic to geocentric coordinates
c of course lon is not changed
c taken from subroutine azimth (Mark Riedesel)
c-----------------------------------------------------------------------
      real*8 lat,lon,glat,glon,rlat
      real*8 pi,eps,radcon,s,w

      pi = 4.d0*datan(1.d0)
      radcon = pi/180.d0
      eps = 1.d0/298.247d0

      rlat = lat*radcon
      w = dsin(rlat)
      s = ((2.d0-eps)*w +4.d0*eps*(w**3))*eps*cos(rlat)

      glat = (rlat-s)/radcon
      glon = lon

      return
      end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine azi(slat1,slon1,slat2,slon2,azim)
c-----------------------------------------------------------------------
c azimuth clockwise from north from location 1 to location 2
c input in geocentric coordinates (degree)
c output in degrees
c based on program azimth due to Mark Riedesel
c-----------------------------------------------------------------------
      real*8 slat1,slon1,slat2,slon2
      real*8 rcolat1,rlon1,rcolat2,rlon2
      real*8 pi,radcon,eps,azim,delta
      real*8 c1,s1,c2,s2,s3,x0,y0,z0,x1,y1,z1,arg
      real*8 x2,y2,z2

      eps = 1.d-5
      pi = 4.d0*datan(1.d0)
      radcon = pi/180.d0

      rcolat1 = pi/2.d0 - slat1*radcon
      rlon1   = slon1*radcon
      rcolat2 = pi/2.d0 - slat2*radcon
      rlon2   = slon2*radcon

      c2 = dcos(rcolat1)
      s2 = dsin(rcolat1)
      c1 = dcos(rlon1)
      s1 = dsin(rlon1)
      s3 = dsin(rcolat2)

c  find the azimuth and distance by rotating the source to the
c  North pole

      x0= s3*cos(rlon2)
      y0= s3*sin(rlon2)
      z0= dcos(rcolat2)
      x1 =  c1*x0 + s1*y0
      y1 = -s1*x0 + c1*y0
      z1 =  z0
      x2=  c2*x1 - s2*z1
      y2=  y1
      z2=  c2*z1 + s2*x1

      arg = dsqrt(x2*x2+y2*y2)
      delta= datan2(arg,z2)
      if(dabs(x2).le.eps.and.dabs(y2).le.eps) then
         azim=0.d0
      else
         azim=datan2(y2,x2)
      end if
      azim  = azim/radcon
      delta = delta/radcon

      azim=180.d0-azim

      return
      end

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      double precision function dist_a(slat1,slong1,slat2,slong2)
c----------------------------------------------------------------------
c     computes the ANGLE in degrees between two points on a sphere 
c     given by their lat,long
c----------------------------------------------------------------------
      real*8 slat1,slong1,slat2,slong2
      real*8 slat11,slong11,slat22,slong22
      real*8 pi,conv,temp

      pi     = 4.d0*datan(1.d0)
      conv   = pi/180.d0

      if (slat1.eq.slat2.and.slong1.eq.slong2) then
        dist_a = 0.d0
        return
      endif

      slat11  = (90.d0-dble(slat1))*conv
      slong11 = dble(slong1)*conv
      slat22  = (90.d0-dble(slat2))*conv
      slong22 = dble(slong2)*conv

      temp = dcos(slat11)*dcos(slat22)+
     &           dsin(slat11)*dsin(slat22)*dcos(slong11-slong22)

      dist_a = dacos(temp)/conv
      return
      end            
