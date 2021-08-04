      subroutine zfcns(c,s)
c
c     zfcns computes z(m,l,theta) and dz(m,l,theta)/dtheta,
c     denoted z and p resp. all functions for 0 le m le
c     max(2,l) and 0 le l le lmax are computed. c is cos(theta)
c     and s is sin(theta).
c
      implicit real*8(a-h,o-z)
      real*8 c, s
c
      include 'parameter5.h'
c
      common/zf/z(3,maxl), p(3,maxl)
      common/mode_hdr/llmax,llmin
c
      data pi/3.14159265358979d0/,tol/1.d-13/
c
      lmax=llmax
      sos=s
      if(s.lt.tol) s=tol
      pi4=4d0*pi
      pi8=8d0*pi
      z(1,1)=1d0/pi4
      z(2,1)=0.d0
      z(3,1)=0.d0
      z(1,2)=3d0*c/pi4
      z(2,2)=3d0*s/pi8
      z(3,2)=0.d0
      z(1,3)=5d0*(3d0*c*c-1d0)/pi8
      z(2,3)=15d0*c*s/pi8
      z(3,3)=15d0*s*s/(4d0*pi8)
      p(1,1)=0.d0
      p(2,1)=0.d0
      p(3,1)=0.d0
      p(1,2)=-3d0*s/pi4
      p(2,2)=3d0*c/pi8
      p(3,2)=0.d0
      p(1,3)=-15d0*c*s/pi4
      p(2,3)=15d0*(c*c-s*s)/pi8
      p(3,3)=15d0*c*s/(2d0*pi8)
      lmaxp1=lmax+1
      tlp1=5d0
      do 10 lp1=4,lmaxp1
         l=lp1-1
         lm1=l-1
         elmm=l
         elpmm1=lm1
         tlp1=tlp1+2d0
         tlm3=tlp1-4d0
         do 10 mp1=1,3
            r=tlp1/elmm
            q=elpmm1/tlm3
            z(mp1,lp1)=r*(c*z(mp1,l)-q*z(mp1,lm1))
            p(mp1,lp1)=r*(c*p(mp1,l)-s*z(mp1,l)-q*p(mp1,lm1))
            elmm=elmm-1d0
            elpmm1=elpmm1+1d0
   10    continue
      s=sos
      return
      end
