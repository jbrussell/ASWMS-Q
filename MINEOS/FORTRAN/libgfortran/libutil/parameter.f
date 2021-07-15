      integer*4 npmax,npmax2,npmax4,nsmax,ncmax,nlmax,nphase
      integer*4 maxcal,mkine,mrayp
c
c     npmax is the standard dimension for arrays -- it is used in
c     some programs that DOUBLE the trace (i.e. xcor) so it needs to
c     be twice what you need -- kind of a space waste, should be 
c     changed eventually...
c
      parameter (npmax  = 400000)
c
c     npmax2 is (npmax+2)/2 for the fft
c
      parameter (npmax2 = 200001)
c
c     what is this used for -- can't find it?
c
      parameter (npmax4 = 640004)
      parameter (nsmax  = 120)
      parameter (ncmax  = 3)
      parameter (nlmax  = 100)
      parameter (nphase = 100)
      parameter (maxcal = 5000)
       parameter (mkine = 10)
       parameter (mrayp = 100)
