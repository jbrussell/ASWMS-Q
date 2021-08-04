c***********************************************************************
c               SPHEROIDAL
c
c contains the modified subroutine spheroidal and the version of it
c that implements a given laterally heterogenous model,i.e
c (source-station sensitive)
c
c      subroutine spheroidal_h2(lcomp,num_modes,ista,
c     &                        del,azm,moment,index,lhetero)
c
c common/vecnl/ changed to common/vecnl_h
c
c PI 1/21/92
c
c***********************************************************************

      subroutine spheroidal_h2(lcomp,num_modes,ista,
     &                        del,azm,moment,index,lhetero)
c
c     calculates spheroidal excitation kernels - for summation
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
      real*8 a
      real*8 wt1, wt2, wt3, wt4 
      real*8 c1sp, s1sp, c2sp, s2sp, dsi
      real*8 e1, e2, e3, e4, au, av
      real*8 vecnl,vecnl_h
      real*8 m1, m2, m3, m4, m5, m6
      real*8 rad,ddelta,dd

c     common block params needed to make w available for dd_c51
      real*4 w, ql, alpha, phi
      real*4 buff(maxmodes,6)
      real*4 abuf(6)
c
      real*4 del(maxstat), azm(maxstat)
      real*4 moment(6)
c
      character*1 lcomp
c
      integer*4 scnl,n, l,nn,ll
      integer*4 num_modes, nentry, index

      logical lhetero

c common block for crust5.1/dbdb options

      logical lcrust51, ldbdb5
      character c51opt*2, cdbopt*2 
      common /coption/ lcrust51,ldbdb5,c51opt,cdbopt

c other commons

      common/vecnl_h/scnl(maxmodes,2),vecnl(maxmodes,4),
     &               vecnl_h(maxmodes,5:6)
      common/zf/z(3,maxl), zp(3,maxl)
      common/aa/ a(maxmodes,maxcomp)
      common/noname/buff
      common/buf/nn, ll, w, ql, alpha, phi
c
      data rad/57.2957795d0/,icomp/0/
c
      equivalence (nn,abuf)
c
      m1 = dble(moment(1))
      m2 = dble(moment(2))
      m3 = dble(moment(3))
      m4 = dble(moment(4))
      m5 = dble(moment(5))
      m6 = dble(moment(6))
c
c     test here for component
c 
ccc   vertical component
c
      if (lcomp .eq. 'Z') then

          azimdp=(180.d0-dble(azm(ista)))/rad
          c1=dcos(azimdp)
          s1=dsin(azimdp)
          c2=dcos(azimdp*2.d0)
          s2=dsin(azimdp*2.d0)
          c1sp=2.d0*c1
          s1sp=2.d0*s1
          c2sp=2.d0*c2
          s2sp=4.d0*s2

          do 25 j = 1, num_modes

            ddelta=dble(del(ista))/rad
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

            dsi = 1.d0/si
            call zfcns(co,si)

            e1 = vecnl(j,1)
            e2 = vecnl(j,2)
            e3 = vecnl(j,3)
            e4 = vecnl(j,4)
            au = vecnl_h(j,5)
            lp1 = scnl(j,2) + 1
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
            a(j,1) = wt1*au
            a(j,2) = (wt2 + wt4*c2sp)*au
            a(j,3) = (wt2 - wt4*c2sp)*au
            a(j,4) = wt3*c1sp*au       
            a(j,5) = wt3*s1sp*au       
            a(j,6) = wt4*s2sp*au
c
c           end of mode loop
c
  25      continue
c
ccc    radial or theta component
c
      else if (lcomp .eq. 'R') then

          azimdp=(180.d0-dble(azm(ista)))/rad
          c1=dcos(azimdp)
          s1=dsin(azimdp)
          c2=dcos(azimdp*2.d0)
          s2=dsin(azimdp*2.d0)
          c1sp=2.d0*c1
          s1sp=2.d0*s1
          c2sp=2.d0*c2
          s2sp=4.d0*s2

          do 26 j = 1, num_modes

            ddelta=dble(del(ista))/rad
cc...computes apparent shif in epicentral distance
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

            dsi = 1.d0/si
            call zfcns(co,si)

            e1 = vecnl(j,1)
            e2 = vecnl(j,2)
            e3 = vecnl(j,3)
            e4 = vecnl(j,4)
            av = vecnl_h(j,6)
            lp1 = scnl(j,2) + 1
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
            a(j,1) = wt1*av
            a(j,2) = (wt2 + wt4*c2sp)*av
            a(j,3) = (wt2 - wt4*c2sp)*av
            a(j,4) = wt3*c1sp*av
            a(j,5) = wt3*s1sp*av       
            a(j,6) = wt4*s2sp*av

c
  26      continue
c
ccc   transverse or phi component
c
      else if (lcomp .eq. 'T') then

          azimdp=(180.d0-dble(azm(ista)))/rad
          c1=dcos(azimdp)
          s1=dsin(azimdp)
          c2=dcos(azimdp*2.d0)
          s2=dsin(azimdp*2.d0)
          c1sp=2.d0*c1
          s1sp=2.d0*s1
          c2sp=2.d0*c2
          s2sp=4.d0*s2

          do 27 j = 1, num_modes

            ddelta=dble(del(ista))/rad
cc...computes apparent shif in epicentral distance
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

            dsi = 1.d0/si
            call zfcns(co,si)

            e3 = vecnl(j,3)
            e4 = vecnl(j,4)
            av = vecnl_h(j,6)
            lp1 = scnl(j,2) + 1
c
c           weight the strains the appropriate Legendre polynomials
c
            wt1 = 0.d0
            wt2 = e4*z(3,lp1)
            wt3 = e3*z(2,lp1)
c
c           form the excitation kernel - av is the acceleration scaler
c           (aw in this case)
c
            a(j,1) = wt1*av*dsi
            a(j,2) = -wt2*s2sp*av*dsi
            a(j,3) = wt2*s2sp*av*dsi
            a(j,4) = -wt3*s1sp*av*dsi
            a(j,5) = wt3*c1sp*av*dsi
            a(j,6) = wt2*4.d0*c2sp*av*dsi
c
c           end of mode loop
c
  27      continue
c
c     end of test on component
c
      endif
c
c     dot in moment tensor if this is a structure problem
c
      if (index .eq. 1) then
          do j = 1, num_modes
           a(j,1) = a(j,1)*m1 + a(j,2)*m2 + a(j,3)*m3
     &               + a(j,4)*m4 + a(j,5)*m5 + a(j,6)*m6
          end do
      end if
c
c     head back
c
      return
      end
c
