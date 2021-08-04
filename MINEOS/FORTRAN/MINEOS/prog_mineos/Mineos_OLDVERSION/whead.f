      subroutine whead(ioeig, jcoma, wmina, wmaxa, lmin, lmax, wgrava)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c   writes two records of header to the formatted output file
c   on the first entry
c
c    calls no other routines
c
      implicit real*8 (a-h,o-z)
c
      include 'parameter.h'
c
      real*4 wminl, wmaxl, wgravl, trefl, tdum(nknot,9)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/head/ tdum
c
      data ifirst /0/
c
      save ifirst
c
c   convert double precision to single precision
c
      wminl = wmina
      wmaxl = wmaxa
      wgravl = wgrava
      trefl = tref
c
      if (ifirst .eq. 0) then
        write(ioeig) jcom, wminl, wmaxl, lmin, lmax, wgravl
        write(ioeig) n, nic, noc, ifanis, trefl, (tdum(i,1),i=1,n), 
     &       (tdum(i,2),i=1,n), (tdum(i,3),i=1,n), (tdum(i,4),i=1,n),
     &       (tdum(i,5),i=1,n), (tdum(i,6),i=1,n), (tdum(i,7),i=1,n),
     &       (tdum(i,8),i=1,n), (tdum(i,9),i=1,n)
      endif
      ifirst = 1
c
      return
      end
