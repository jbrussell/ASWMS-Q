      subroutine azimth(slat,slon,rlat,rlon,delta,azim,bazim)
c
c   This routine uses Euler angles to find the geocentric 
c   distance, azimuth, and back azimuth for a source-reciever
c   pair.
c
c                      Input
c   slat  - source geographic latitude in decimal degrees
c   slon  - source longitude in decimal degrees
c   rlat  - receiver geographic latitude in decimal degrees
c   rlon  - receiver longitude in decimal degrees
c
c                     Output
c   delta - source-reciever distance in decimal degrees of arc
c   azim  - azimuth from the source to the reciever in degrees
c   bazim - back azimuth from the receiver to the source in degrees
c
c   Mark Riedesel, January 30, 1986
c
      include 'numerical.h'
      dtor = drad
      e = 1.0/flt
c
c    convert to geocentric coordinates and from latitude to 
c    colatitude
c
      slatra = dtor*slat
      w = sin(slatra)
      s = ((2.-e)*w +4.*e*(w**3))*e*cos(slatra)
      scolat = pi2 - slatra + s
c
      rlatra = dtor*rlat
      w = sin(rlatra)
      s = ((2.-e)*w +4.*e*(w**3))*e*cos(rlatra)
      rcolat = pi2 - rlatra + s
c
      slonra=slon*dtor
      rlonra=rlon*dtor
      c2=cos(scolat)
      s2=sin(scolat)
      c1=cos(slonra)
      s1=sin(slonra)
c
c  find the azimuth and distance by rotating the source to the
c  North pole
c
      slatrc=sin(rcolat)
      x0=slatrc*cos(rlonra)
      y0=slatrc*sin(rlonra)
      z0=cos(rcolat)
c
      x1=c1*x0+s1*y0
      y1=-s1*x0+c1*y0
      z1=z0
      x2=c2*x1-s2*z1
      y2=y1
      z2=c2*z1+s2*x1
c
      call angles(x2,y2,z2,delta,azim)
      azim = 180.-azim
c
c  find the back azimuth by rotating the reciever to the 
c  North pole
c
      c2=cos(rcolat)
      s2=sin(rcolat)
      c1=cos(rlonra)
      s1=sin(rlonra)
c
      slatrc=sin(scolat)
      x0=slatrc*cos(slonra)
      y0=slatrc*sin(slonra)
      z0=cos(scolat)
c
      x1=c1*x0+s1*y0
      y1=-s1*x0+c1*y0
      z1=z0
      x2=c2*x1-s2*z1
      y2=y1
      z2=c2*z1+s2*x1
c
      call angles(x2,y2,z2,delta,bazim)
      bazim = 180.-bazim
      return
      end
