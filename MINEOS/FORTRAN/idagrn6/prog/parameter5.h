c pamameter file for idagrn5 -- merging of iidagrn4 param file
c (hetero.param.h) with additional params need for the frechet
c files (parameter.h)
c-----------------------------------------------------------
c
c parameter for lateral heterogeneous models 
c M84A, M84C, SH8....
c
c nord = total angular oder for complex spherical harmonics
c        routines (also defined in /usr/local/src/Libs/Harm/harm_param.h;
c        make sure that nord in harm_param.h .ge. nord here)
c kmax = order of Legendre polynomials (radial dependency)
c        (M84..)
c kmax1 = order of Chebychev polynomials (radial dependency) SH8-Whole mantle
c kmax2 = SH8-U lower mantle
c kmax3 = SH8-U upper mantle      
c------------------------------------------------------------

      integer*4 nord,kmax,kmax1,kmax2,kmax3
      
      parameter (nord  = 36)
      parameter (kmax  =  3)
      parameter (kmax1 = 13)
      parameter (kmax2 =  8)
      parameter (kmax3 =  4)

c...........................................................
c  mode handling for Spheroidal modes up to 25/50 mHz
c...........................................................

      integer*4 maxn,maxl,maxmodes,nbranch
    
C      parameter (maxn =  75)
C      parameter (maxl = 275)
C      parameter (maxmodes = 4300)

      parameter (maxn = 4000)
C      parameter (maxl = 700)
      parameter (maxl = 10000)
      parameter (maxmodes = 400000)
      parameter (nbranch = 4000)

c...........................................................
c  stations,nknots
c...........................................................

      integer*4 maxstat,maxtime,maxcomp

c      parameter (maxstat =  100)
c      parameter (maxtime = 3610)

c      parameter (maxstat =  100)
      parameter (maxstat =  2000)

c      parameter (maxtime = 7210)
      parameter (maxtime = 260000)

      parameter (maxcomp = 6)

c additional mode parameters for frechet files, etc.

       integer*4 nknot_t,nknot_s,nknot,nknot3,nknot4,nknot5,nknot6
       integer*4 nknot10,nknot14,maxbyte3,maxbyte4,maxbyte5,maxbyte
       integer*4 maxdisc,maxdh,mbuf,mfrechet
       real*4    rfrctn
c
c     more realistic estimates of the number of knots
c
      parameter (nknot_t = 600)
      parameter (nknot_s = 800)
c
c     knot definitions for all raw mineos programs
c
      parameter (nknot = 1000)
      parameter (nknot3 = 3*nknot)
      parameter (nknot4 = 4*nknot)
      parameter (nknot5 = 5*nknot)
      parameter (nknot6 = 6*nknot)
      parameter (nknot9 = 9*nknot)
      parameter (nknot10 = 10*nknot)
      parameter (nknot14 = 14*nknot)
      parameter (nknot18 = 18*nknot)
      parameter (maxbyte3 = nknot3+5)
      parameter (maxbyte4 = nknot4+5)
      parameter (maxbyte5 = nknot5+5)
      parameter (maxbyte = nknot6+5) 

      parameter (maxdisc = 20)
      parameter (maxdh = 50000)
      parameter (rfrctn = 10.)

      parameter (mbuf = 6+6*nknot_s)
      parameter (mfrechet = 6+6*nknot_s + maxdisc)


c     parameter for crust5.1 correcitions
      parameter (nlith=30)



