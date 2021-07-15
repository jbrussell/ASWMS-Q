c-----------------------------------------------------------
c
c parameter for lateral heterogeneous models 
c M84A, M84C, SH8....
c
c nord = total angular oder for complex spherical harmonics
c        routines
c kmax = order of Legendre polynomials (radial dependency)
c        (M84..)
c kmax1 = order of Chebychev polynomials (radial dependency) SH8-Whole mantle
c kmax2 = SH8-U lower mantle
c kmax3 = SH8-U upper mantle      
c------------------------------------------------------------

      integer*4 nord,kmax,kmax1,kmax2,kmax3
      
      parameter (nord  = 12)
      parameter (kmax  =  3)
      parameter (kmax1 = 13)
      parameter (kmax2 =  8)
      parameter (kmax3 =  4)

c...........................................................
c  mode handling for Spheroidal modes up to 25/50 mHz
c...........................................................

      integer*4 maxn,maxl,maxmodes
    
C      parameter (maxn =  75)
C      parameter (maxl = 275)
C      parameter (maxmodes = 4300)

      parameter (maxn = 145)
      parameter (maxl = 580)
      parameter (maxmodes = 16800)

c...........................................................
c  stations,nknots
c...........................................................

      integer*4 maxstat,nknot,nknot6,maxbyte,maxtime

c      parameter (maxstat =   40)
c      parameter (maxtime = 3610)

      parameter (maxstat =  100)
      parameter (maxtime = 7210)

      parameter (nknot   = 1000)
      parameter (nknot6  = 6000)
      parameter (maxbyte = 6005) 
      parameter (nbranch = 145)
      parameter (maxcomp = 6)


