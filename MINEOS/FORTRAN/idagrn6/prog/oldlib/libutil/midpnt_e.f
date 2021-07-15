      subroutine midpnt_e(ealat,ealong,slat,slong,theta,phi)
c
c     Subroutine to determine the mid-point between two points on the
c     earth's surface.
c
c     This version uses geocentric coordinates
c
c     inputs
c       ealat is earthquake latitude
c       ealong is earthquake longitude
c       slat is station latitiude
c       slong is station longitude
c     outputs
c       theta is mid-point latitude
c       phi is mid-point longitude
c
c     all latitudes and longitudes are in degrees
c       latitude is positive north, negative south
c       longitude is positive east, negative west
c
      real*4 n, l1, l2, l3, mm1, mm2, mm3
c
      common/b2/l1,l2,l3,mm1,mm2,mm3,factor
c
      data n/1hn/,s/1hs/,e/1he/,w/1hw/
c
      include 'numerical.h'
c
      factor = drad
      e = 1.0/flt
      a = (1.0 - e)**2
      b = 1.0/a
c
c     Aa1,aa2,aa3 are cartesian coordinates of projected earthquake sour
c     c1,c2,c3 are cartesian coordinates of station
c
      eathet=ealat*factor
ccc
      beta = atan( a * tan(eathet))
      eathet = beta
ccc
      eaphi=ealong*factor
      aa1= cos(eathet)* cos(eaphi)
      aa2= cos(eathet)* sin(eaphi)
      aa3= sin(eathet)
      stheta=slat*factor
ccc
      beta = atan( a * tan(stheta))
      stheta = beta
ccc
      sphi=slong*factor
      c1= cos(stheta)* cos(sphi)
      c2= cos(stheta)* sin(sphi)
      c3= sin(stheta)
      mm1=aa3*c2-aa2*c3
      mm2=aa1*c3-aa3*c1
      mm3=aa2*c1-aa1*c2
      l1=aa1-c1
      l2=aa2-c2
      l3=aa3-c3
      r1=l3*mm2-l2*mm3
      r2=l1*mm3-l3*mm1
c     Test to find right phi
      z=r2/r1
      ph1= atan(z)
      phi1=ph1/factor
      if(phi1.le.0.) phi2=phi1+180.
      if(phi1.gt.0.) phi2=phi1-180.
      clat1=fcn(phi1)
      clat2=fcn(phi2)
ccc
      beta = clat1*factor
      clat1 = atan(b*tan(beta))/factor
      beta = clat2*factor
      clat2 = atan(b*tan(beta))/factor
ccc
c
c     Pick right solution
c
      q1=clat1*factor
      q2=phi1*factor
      q3=clat2*factor
      q4=phi2*factor
      b11= cos(q1)* cos(q2)
      b21= cos(q1)* sin(q2)
      b31= sin(q1)
      b12= cos(q3)* cos(q4)
      b22= cos(q3)* sin(q4)
      b32= sin(q3)
      dot1=aa1*b11+aa2*b21+aa3*b31
      dot2=aa1*b12+aa2*b22+aa3*b32
      if(dot1.gt.dot2) phi=phi1
      if(dot2.gt.dot1) phi=phi2
      if(phi.eq.phi1) theta=clat1
      if(phi.eq.phi2) theta=clat2
 1000 continue
      return
      end
