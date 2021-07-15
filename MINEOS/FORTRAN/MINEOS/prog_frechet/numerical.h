c
c     include file for useful constants
c
      pi = 3.14159265350
      tpi = 2.0*pi
      pi2 = 0.5*pi
      pi4 = 0.25*pi
c
      drad = pi/180.
c
      rmhz = 1000./tpi
      rad = 1.0/rmhz
c
      rn = 6371000.0
      rnk = 6371.
c
      dkm = rnk*drad
c
      bigg = 6.6732e-11
      rhobar = 5515.0
c
      third = 1.0/3.0
      tthird = 2.0/3.0
      fthird = 4.0/3.0
c
c     for Barbara's mode scalings
c
      gn = pi*bigg*rhobar*rn
      vn2 = gn*rn
      vn = sqrt(vn2)
      wn = vn/rn
c
