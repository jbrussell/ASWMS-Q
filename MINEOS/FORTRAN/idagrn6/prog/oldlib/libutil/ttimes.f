      subroutine ttimes(zs,delta,plist,n,tt,dtdd,dtdh,dddp,phcd)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c     zs - source depth (km)
c     delta - distance in degrees
c     plist - phase request
c    
c     n - number of phases
c     tt - array of travel times
c     dtdd - array of d time /d distance values
c     dtdh - array of d time /d depth values
c     dddp - array of d distane /d ray parameter values
c     phcd - array of phase names
c
      save
c
      include 'parameter.f'
c
      logical prnt(3), lprem
      character*8 phcd(nphase),phlst(1), plist, plast
      character*80 modnam
      character*60 seismdir
      real*4 tt(nphase),dtdd(nphase),dtdh(nphase),dddp(nphase)
      real*4 usrc(2)
c
      data plast /' '/
c
      common /t_slow/ usrc
      common /t_model/ lprem
c
      in = 1
c
      call getenv('SEISM',seismdir)
      call klen(seismdir,kl)
c
      if (lprem) then
        modnam = seismdir(1:kl)//'/prem'
      else
        modnam = seismdir(1:kl)//'/iasp91'
      end if
c
      do ii = 1, 3
        prnt(ii) = .false.
      end do
      do ii = 1, nphase
        tt(ii) = 0.0
        dtdd(ii) = 0.0
        dtdh(ii) = 0.0
        dddp(ii) = 0.0
        phcd(ii) = ' '
      end do
      usrc(1) = 0.0
      usrc(2) = 0.0
      phlst(1) = plist
      n = 0
      nc = 1
c
      if (plast .ne. plist) then
        close(in)
        call tabin(in,modnam)
        call brnset(nc,phlst,prnt)
        plast = plist
      endif
c 
      call depset(zs,usrc)
c 
      call trtm(delta,nphase,n,tt,dtdd,dtdh,dddp,phcd)
c
      end
