c
      subroutine moment_rtf(strike,dip,rake,mo,xm)
c
c     determine components of moment tensor from strike, dip, and rake
c     ala Aki and Richards, p. 117
c
      real*4 xm(*), mo
      real*4 mxx, myy, mzz, mxy, mxz, myz
      include 'numerical.h'
c
      st = strike*drad
      di = dip*drad
      ra = rake*drad
c
      sid = sin(di)
      si2d = sin(2.*di)
      cod = cos(di)
      co2d = cos(2.*di)
      sir = sin(ra)
      cor = cos(ra)
      sins = sin(st)
      sin2s = sin(2.*st)
      coss = cos(st)
      cos2s = cos(2.*st)
c
      mxx = -mo*(sid*cor*sin2s + si2d*sir*sins*sins)
      mxy = mo*(sid*cor*cos2s + 0.5*si2d*sir*sin2s)
      mxz = -mo*(cod*cor*coss + co2d*sir*sins)
      myy = mo*(sid*cor*sin2s - si2d*sir*coss*coss)
      myz = -mo*(cod*cor*sins - co2d*sir*coss)
      mzz = mo*si2d*sir
c
c     convert to spherical coordinate system
c
c     mzz = mrr
c     mxx = mtt
c     myy = mff
c     mxy = -mtf
c     mxz = mrt
c     myz = -mrf
c
      xm(1) = mzz
      xm(2) = mxx
      xm(3) = myy
      xm(4) = mxz
      xm(5) = -myz
      xm(6) = -mxy
c
      return                                                        
      end 
