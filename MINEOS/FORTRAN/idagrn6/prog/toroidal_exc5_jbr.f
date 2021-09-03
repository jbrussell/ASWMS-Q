c***********************************************************************
c               TOROIDAL
c
c contains the modified subroutine toroidal and the version of it
c that implements a given laterally heterogenous model,i.e
c (source-station sensitive)
c
c      subroutine toroidal_exc_h2(lcomp,num_modes,ista,
c     &                          del,azm,moment,lhetero)
c
c common/vecnl/ changed to common/vecnl_h
c
c PI 1/21/92
c
c jbr 09/01/21 - modified slightly to output ascii file rather than binary.
c              - ak: excitation amplitude
c              - phi2: excitation phase
c jbr - I think the excitation equations follow 10.42-10.48 of Dahlen and Tromp
c
c***********************************************************************

      subroutine toroidal_exc_h2_jbr(lcomp,num_modes,ista,
     &               del,azm,moment,lhetero,sta)
c
c     calculates toroidal excitation kernels - for source phase
c     more or less identical to subroutine toroidal
c
c     z functions - these are the product of bcoeff and x functions
c       note to those who care:
c         the arrays in zfcns are from 1:maxl, which means l = l + 1
c         in the common blocks in idagrn, the arrays are 0:maxl, which
c         means that the zfcns can be correctly retrieved by l
c
      implicit real*8 (a-h,o-z)
c
      include 'parameter5.h'
c
      real*8 co,si,c1,c2,s1,s2,z,zp,azimdp
      real*8 wt1, wt2, wt3
      real*8 c1sp, s1sp, c2sp, s2sp
      real*8 e1, e2, e3, e4, au, av
      real*8 vecnl,vecnl_h,pi
      real*8 sigma1, sigma2, lam, aj
      real*8 mm1, mm2
      real*8 a, phi2, phip, ak, ae
      real*8 sign1

c     common block params needed to make w available for dd_c51
      real*4 w, ql, alpha, phi
      real*4 buff(maxmodes,6)
      real*4 abuf(6)
c
      real*4 del(maxstat), azm(maxstat)
      real*4 moment(6)
      real*4 wref, dep

      logical lhetero
      
c jbr - define station variable
      character*10 sta
c
      character*1 lcomp
      character*4 stn(maxstat), comps
      character*8 pre
      character*12 modeid
      character*80 filep
      character*256 exf
c
      integer*4 scnl,n, l, nn, ll
      integer*4 num_modes
      integer*4 evt(5), jcom
c
c common block for crust5.1/dbdb options

      logical lcrust51, ldbdb5
      character c51opt*2, cdbopt*2 
      common /coption/ lcrust51,ldbdb5,c51opt,cdbopt

c other commons

      common/vecnl_h/scnl(maxmodes,2),vecnl(maxmodes,4),
     &               vecnl_h(maxmodes,5:6)
      common/zf/ z(3,maxl), zp(3,maxl)
      common/ahdr/ jcom, knots
      common/exc/ evt, filep, modeid, stn, comps, dep
      common/noname/buff
      common/buf/nn, ll, w, ql, alpha, phi
      
c jbr - read in frequencies     
      real*8 wsave, wdb
      common/wsave_jbr/wsave(maxmodes)
c
      data rad/57.29578d0/
      data pi/3.14159265359d0/
c
      equivalence (nn,abuf)
      pi4 = pi/4.d0
c
c     set up a few things for the file name  -- if standard 
c     p(8 char date).dir format, then use date.  If not, set
c     a temp name
c
      call kblnk(filep,kk)
      do i = 1,kk
       if (filep(i:i+3) .eq. '.dir') then
        if (filep(i-9:i-9) .eq. 'p') then
         pre = filep(i-8:i-1)
        else
         pre = 'tempname'
        endif
       else
        pre = 'tempname'
       endif
      enddo
c      if (filep(1:4) .eq. 'data') then
c        ks = 10
c      else
c        ks = 1
c      endif
c      call kblnk(modeid,km)
      nup = 0
      ndown = 99999
      lup = 0
      ldown = 99999
      do j = 1, num_modes
        nup = max0(scnl(j,1),nup)
        ndown = min0(scnl(j,1),ndown)
        lup = max0(scnl(j,2),lup)
        ldown = min0(scnl(j,2),ldown)
      end do
c
      ddelta=dble(del(ista))/rad
      azimdp=(180.d0-dble(azm(ista)))/rad
      c1=dcos(azimdp)
      s1=dsin(azimdp)
      c2=dcos(azimdp*2.d0)
      s2=dsin(azimdp*2.d0)
      c1sp=2.d0*c1
      s1sp=2.d0*s1
      c2sp=2.d0*c2
      s2sp=4.d0*s2
      if (stn(ista)(4:4) .eq. ' ') then
        kl = 3
      else
        kl = 4
      endif
c
      mm1 = moment(4)*s1 - moment(5)*c1
      mm2 = 0.5d0*s2*(moment(2)-moment(3)) - c2*moment(6)
c
c     test here for component
c 
ccc   transverse or phi component
c
      if (lcomp .eq. 'T') then
c          exf = filep(ks:kk)//'.'//stn(ista)(1:kl)//'.'//modeid(1:km)/
c     &          /'.'//lcomp//'.exc'
C          exf = pre(1:8)//'.'//stn(ista)(1:kl)//'.'//lcomp//'.exc'
          exf = trim(sta)//'.'//lcomp//'.excite.asc'
          print*,'opening output file ',exf
C          open(unit=4,file=exf,form='unformatted',access='sequential')
          open(unit=4,file=exf,form='formatted',access='sequential')
C          write(4) jcom, num_modes, nup, ndown, lup, ldown
C          write(4) (evt(jj),jj=1,5),stn(ista),comps,del(ista),dep
c
          do 27 j = 1, num_modes
cc...computes apparent shift in epicentral distance
            if (lhetero) then
              nn = scnl(j,1) 
              ll = scnl(j,2)
              call calc_dd(nn,ll,ista,dd)
              co=dcos(ddelta+dd)
              si=dsin(ddelta+dd)
            elseif (lcrust51 .and. c51opt(1:2).eq.'gc') then
              do ii = 1, 6
                abuf(ii) = buff(j,ii)
              end do
              call dd_c51(nn,ll,w,ddelta,ista,dd)
              co=dcos(ddelta+dd)
              si=dsin(ddelta+dd)
            elseif (ldbdb5 .and. cdbopt(1:2).eq.'gc') then
              do ii = 1, 6
                abuf(ii) = buff(j,ii)
              end do
              call dd_dbd(nn,ll,w,ddelta,ista,dd)
              co=dcos(ddelta+dd)
              si=dsin(ddelta+dd)
            else
              nn = scnl(j,1) 
              ll = scnl(j,2)
              co=dcos(ddelta)
              si=dsin(ddelta)
            endif
            wdb = wsave(j)
c
            call zfcns(co,si)
c
            e1 = vecnl(j,1)
            e2 = vecnl(j,2)
            e3 = vecnl(j,3)
            e4 = vecnl(j,4)
            au = vecnl_h(j,5)
            av = vecnl_h(j,6)
            lp1 = ll + 1
c
c           weight the strains by the appropriate Legendre polynomial
c
            wt1 = 0.d0
            wt2 = e1*zp(2,lp1)
            wt3 = e2*zp(3,lp1)
c
c           form the excitation kernel - au is the acceleration scaler
c           (aw in this case)
c
            a  = (wt1*dble(moment(1))
     &         + wt3*s2*dble(moment(2)-moment(3))
     &         + wt2*(dble(moment(4))*s1sp-dble(moment(5))*c1sp)
     &         - wt3*c2sp*dble(moment(6)))*au
c
            e2 = e2*0.25d0
            lam = dfloat(ll) + 0.5d0
            aj = (dsqrt((lam)**5)/(2.d0*pi))*dsqrt(2.d0/(pi*si))
            sigma1 = e1*mm1
            sigma2 = lam*e2*mm2
            if ((sigma1 .eq. 0.d0) .and. (sigma2 .eq. 0.d0)) then
              ak = -9999.d0
              phi2 = 0.d0
            else
              ak   = dsqrt(sigma1**2 + sigma2**2)*(au*aj*1.d06)
              phi2  = datan2(sigma2,sigma1)
              sign1 = dsign(1.d0,ak)
              if (sign1 .lt. 0.d0) then
                ak = dabs(ak)                     
                phi2  = phi2 + pi
              endif
            endif
            phip = 0.d0
            ae = a*1.d06
c
c           note that these are now real numbers (not double precision)
c
C            write(4) nn, ll, sngl(ak), sngl(phi2), sngl(phip), sngl(ae)
            write(4,10) nn, ll, wdb, sngl(ak), sngl(phi2), sngl(phip), sngl(ae)
  10        format(I3,I5,F15.8,F25.8,F25.8,F25.8,F25.8)
ccc
ccc         optional check to amplitude
ccc
c           if (ae .ne. 0.0) then
c             aj = ak*dcos(lam*ddelta - pi4 - phi2)
c             write(6,'(2i4,5f15.8,2g16.8)') nn, ll, ak, phi2, aj, ae,
c    &          (aj - ae)*100./ae, sigma1, sigma2
c           endif
ccc
ccc         end optional check
ccc
c
c           end of mode loop
c
  27      continue
c
          close(4)
c
c     end of test on component
c
      endif
c
c     head back
c
      return
      end
