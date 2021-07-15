c
      subroutine moment_a(xm,azm,a)
c
c     compute "a" functions for adding "fundamental faults"
c     ala Jost and Herrmann, SRL, 60, 37-57, 1989
c
      real*4 xm(*), a(*), azm, az
      real*4 mxx, myy, mzz, mxy, mxz, myz
      include 'numerical.h'
c
c     convert from spherical coordinate system to Aki and Richards
c
      mzz = xm(1)
      mxx = xm(2)
      myy = xm(3)
      mxz = xm(4)
      myz = -xm(5)
      mxy = -xm(6)
c
      az = azm * drad
      c1 = cos(az)
      c2 = cos(2.*az)
      s1 = sin(az)
      s2 = sin(2.*az)
c
      a(1) = 0.5*(mxx - myy)*c2 + mxy*s2
      a(2) = mxz*c1 + myz*s1
      a(3) = -0.5*(mxx + myy)
      a(4) = 0.5*(mxx - myy)*s2 - mxy*c2
      a(5) = mxz*s1 - myz*c1
      a(6) = third*(mxx + myy + mzz)
c
      return                                                        
      end 
