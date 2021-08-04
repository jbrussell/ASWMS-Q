c***********************************************************************
c               SPHEROIDAL_EXC
c
c contains the modified subroutine spheroidal_exc and the version of it
c that implements a given laterally heterogenous model,i.e
c (source-station sensitive)
c
c      subroutine spheroidal_exc_h2(lcomp,num_modes,ista,
c     &                            del,azm,moment,lhetero)
c
c common/vecnl/ changed to common/vecnl_h
c
c PI 1/21/92
c
c PP 05/08/95 bug for case dd.ne.0.d0 fixed
c
c***********************************************************************
      subroutine spheroidal_exc_h2(lcomp,num_modes,ista,
     &             del,azm,moment,lhetero)
c
c     calculates spheroidal excitation kernels - for source phase
c     more or less identical to subroutine spheroidal
c
      implicit real*8 (a-h,o-z)
c
      include 'parameter5.h'
c
      real*8 co,si,c1,c2,s1,s2,z,zp,cz,sz,caz,saz,azimdp
      real*8 wt1, wt2, wt3, wt4 
      real*8 c1sp, s1sp, c2sp, s2sp, dsi
      real*8 e1, e2, e3, e4, au, av
      real*8 vecnl,vecnl_h, pi
      real*8 sigma1, sigma2, lam, aj
      real*8 mm1, mm2, mm3, mm4
      real*8 a, phi2, phip, ak, ae, sign1
c
c     common block params needed to make w available for dd_c51
      integer*4 scnl, nn, ll
      real*4 w, ql, alpha, phi
      real*4 buff(maxmodes,6)
      real*4 abuf(6)

      real*4 del(maxstat), azm(maxstat)
      real*4 moment(6)
      real*4 wref,dep

      logical lhetero
c
      character*1 lcomp
      character*4 stn(maxstat), comps
      character*8 pre
      character*12 modeid
      character*80 filep
      character*256 exf
c
c      integer*4 scnl,n, l, nn, ll
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
c
      data rad/57.29578d0/
      data pi/3.14159265359d0/
c
      equivalence (nn,abuf)
c
      pi4 = pi/4.d0
c
c     set up a few things for the file name
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
c      call kblnk(filep,kk)
c      if (filep(kk-3:kk) .eq. '.dir') then
c        kk = kk - 4           
c      endif
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
      mm1 = moment(1)
      mm2 = moment(2) + moment(3)
      mm3 = c1*moment(4) + s1*moment(5)
      mm4 = c2*(moment(2) - moment(3)) + 2.d0*s2*moment(6)
c
c     test here for component
c 
ccc   vertical component
c
      if (lcomp .eq. 'Z') then
c
c          exf = filep(ks:kk)//'.'//stn(ista)(1:kl)//'.'//modeid(1:km)/
c     &          /'.'//lcomp//'.exc'
          exf = pre(1:8)//'.'//stn(ista)(1:kl)//'.'//lcomp//'.exc'
          print*,'opening output file ',exf
          open(unit=4,file=exf,form='unformatted',access='sequential')
          write(4) jcom, num_modes, nup, ndown, lup, ldown
          write(4) (evt(jj),jj=1,5),stn(ista),comps,del(ista),dep
c
          do 25 j = 1, num_modes
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
c           weight the strains the appropriate Legendre polynomials
c
            wt1 = e1*z(1,lp1)
            wt2 = e2*z(1,lp1)
            wt3 = e3*z(2,lp1)
            wt4 = e4*z(3,lp1)
c
c           form the excitation kernel - au is the acceleration scaler
c
            a  = (wt1*dble(moment(1))
     &         + (wt2 + wt4*c2sp)*dble(moment(2))
     &         + (wt2 - wt4*c2sp)*dble(moment(3))           
     &         + wt3*(c1sp*dble(moment(4)) + s1sp*dble(moment(5)))
     &         + wt4*s2sp*dble(moment(6)))*au
c
            lam = float(ll) + 0.5d0
            aj = dsqrt(lam)/(2.d0*pi)*dsqrt(2.d0/(pi*si))
            sigma1 = (e1*mm1 + e2*mm2) - 0.25d0*e4*mm4*lam**2
            sigma2 = e3*mm3*lam
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
            write(4) nn, ll, sngl(ak), sngl(phi2), sngl(phip), sngl(ae)
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
  25      continue
c
          close(4)
c
c
ccc    radial or theta component
c
      else if (lcomp .eq. 'R') then
c
c          exf = filep(ks:kk)//'.'//stn(ista)(1:kl)//'.'//modeid(1:km)/
c     &          /'.'//lcomp//'.exc'
c          open(unit=4,file=exf,form='unformatted',access='sequential')
          write(4) jcom, num_modes, nup, ndown, lup, ldown
          write(4) (evt(jj),jj=1,5),stn(ista),comps,del(ista)
c
          do 26 j = 1, num_modes
            nn = scnl(j,1) 
            ll = scnl(j,2)
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
c           weight the strains the appropriate Legendre polynomials
c
            wt1 = e1*zp(1,lp1)
            wt2 = e2*zp(1,lp1)
            wt3 = e3*zp(2,lp1)
            wt4 = e4*zp(3,lp1)
c
c           form the excitation kernel - av is the acceleration scaler
c
            a  = (wt1*dble(moment(1))
     &            + (wt2 + wt4*c2sp)*dble(moment(2))
     &            + (wt2 - wt4*c2sp)*dble(moment(3))
     &            + wt3*(c1sp*dble(moment(4)) + s1sp*dble(moment(5)))
     &            + wt4*s2sp*dble(moment(6)))*av
c

            lam = float(ll) + 0.5d0
            aj = dsqrt(lam**3)/(2.d0*pi)*dsqrt(2.d0/(pi*si))
            sigma1 = -((e1*mm1 + e2*mm2) - 0.25d0*e4*mm4*lam**2)
            sigma2 = e3*mm3*lam
            if ((sigma1 .eq. 0.d0) .and. (sigma2 .eq. 0.d0)) then
              ak = -9999.d0
              phi2 = 0.d0
            else
              ak   = dsqrt(sigma1**2 + sigma2**2)*(av*aj*1.d06)
              phi2  = datan2(sigma1,sigma2)
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
            write(4) nn, ll, sngl(ak), sngl(phi2), sngl(phip), sngl(ae)
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
  26      continue
c
          close(4)
c
c
c     end of test on component
c
      endif
c
c     head back
c
      return
      end
