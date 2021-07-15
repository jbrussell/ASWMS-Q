      subroutine spread(p,dpdd,rs,rr,s_slow,r_slow,delta,b)
c
c     subroutine to calculate geometrical spreading coefficient b
c
c     p - ray parameter in s/deg
c     dpdd - dp/ddelta in s/rad/rad
c     rs - source radius in km
c     rr - receiver radius in km
c     s_slow - source slowness in s/km (scaled)
c     r_slow - receiver slowness in s/km (scaled)
c     delta - distance in deg
c     b - geometrical spreading coefficient
c
c     spreading checked 6/9/92 against figure from Kovach and Anderson
c     1/r**2 scaling is removed for stability!
c
c     input assumed from ttimes
c
      real*4 p, dpdd, rs, rr, s_slow, r_slow, delta
      real*4 s_vel, r_vel
c
      include 'numerical.h'
c
c     convert ray parameter and its derivative to s
c
      pk = p * radd
c
c     unscale the slowness and convert to velocity
c
      s_vel = (1.0/s_slow)*(rs/rnk)
      r_vel = (1.0/r_slow)*(rr/rnk)
c
c     5/29/92 - my thinking on units may be incorrect
c     p is in sec/deg but dpdd may be in sec/rad/rad already
c
c     dpk = dpdd * radd * radd
c
      dpk = dpdd
c
      dist = delta * drad
c
c     calculate the take off angles
c
      xs = asin(pk*s_vel/rs)
      xr = asin(pk*r_vel/rr)
c
      xa = tan(xs)/(sin(dist)*cos(xr))
      b = sqrt(xa * (s_vel/rs) * (abs(dpk)/(rr**2)))
c
c     1/r**2 scaling removed for stability
c
ccccccc
      b = sqrt(xa * (s_vel/rs) * (abs(dpk)))
ccccccc
c
c     print*,' geometrical spreading info'
c     print*, pk, dpk, xs*radd, xr*radd, delta, s_vel, r_vel
c
      return
      end
