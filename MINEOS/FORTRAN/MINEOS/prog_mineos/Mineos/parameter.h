       integer*4 nknot_t,nknot_s,nknot,nknot3,nknot4,nknot5,nknot6
       integer*4 nknot10,nknot14,maxbyte3,maxbyte4,maxbyte5,maxbyte
       integer*4 maxmodes,maxmodes1,maxl,maxll,maxtime,maxcomp,maxstat
       integer*4 maxn,maxold,maxold2,maxold3,nbranch,nprop,maxdisc,maxdh
       integer*4 mbuf,mfrechet,lhdr
       real*4    rfrctn
c
c     more realistic estimates of the number of knots
c     need more t knots for iasp91 model (5 km mantle)
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
c
c     other parameters for mineos programs & idagrn
c

c      parameter (nbranch = 150)
      parameter (nbranch = 2000)
      parameter (nbranch2 = 2*nbranch)
      parameter (maxn = 2000)
      parameter (maxl = 50000)
      parameter (maxll = 50000)
      parameter (maxold = maxn*maxll)
      parameter (maxmodes = 300000)
      parameter (maxmodes1 = 300001)

c      parameter (nbranch = 150)
c      parameter (nbranch2 = 2*nbranch)
c      parameter (maxn = 150)
c      parameter (maxl = 700)
c      parameter (maxll = 1100)
c      parameter (maxold = maxn*maxl)
c      parameter (maxmodes = 25800)
c      parameter (maxmodes1 = 25801)
c
c     parameters for idagrn 
c
      parameter (maxtime = 7210)
      parameter (maxstat = 100)
      parameter (maxcomp = 6)
c
c     parameter for plot_wk
c
      parameter (nprop = 29)
c
c     parameters of mineos_frechet_new
c
      parameter (maxdisc = 30)
      parameter (maxdh = 50000)
      parameter (rfrctn = 10.)
c
c     6 + 1000*3 + 20 
c    yielding a desire to make the arrays smaller:
c     6 + 800*3 + 20
c    but need to consider anisotropy
c     6 + 800*6 + 20 
c
      parameter (mbuf = 6+6*nknot_s)
      parameter (mfrechet = mbuf + maxdisc)
c
c     mineos_partial
c
      parameter (lhdr = 100)
c
c  nfreq = 2*number of modes - maximum
c  nknot = number of knots in model
c  nknot3 = 3*number of knots
c  nknot6 = 6*number of knots
c  and so forth
c 
c  maxbyte = 6*knot + 5
c
c  if you change nknot, be sure to change the other parameters as well
c

