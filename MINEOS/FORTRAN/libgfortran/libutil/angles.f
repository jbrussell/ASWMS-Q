


      subroutine angles(x,y,z,theta,phi)
c
c  finds the angles theta and phi of a spherical polar coordinate
c  system from the cartesion coordinates x, y, and z.
c
c  Mark Riedesel, 1983
c
      include 'numerical.h'
c
      eps = 1.e-4
      rtod = radd
c
      arg1=sqrt(x*x+y*y)
      theta=atan2(arg1,z)
      if(abs(x).le.eps.and.abs(y).le.eps) then
         phi=0.
      else
         phi=atan2(y,x)
      end if
      phi=phi*rtod
      theta=theta*rtod
      return
      end
