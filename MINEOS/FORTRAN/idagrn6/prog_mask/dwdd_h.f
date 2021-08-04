c***********************************************************************
c                 PERTURBATIONS 2
c
c contains the  subroutines needed to compute dw, the frequency shift
c due to aphericity and dd, the apparent shift in epicentral distance
c
c      subroutine calc_dw(nn,ll,nsta,dw)
c      subroutine calc_dd(nn,ll,nsta,dw)
c
c      subroutine ang_param(elat,elon,slat,slon,
c     &                    delta,azmp,theta,avg_great,avg_diff,ismax)
c
c***********************************************************************

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine calc_dw(nn,ll,nsta,dw)
c-----------------------------------------------------------------------
c computes the shift in eigenfrenquency due to asphericity
c (lateral heterogeneity and ellipticity of figure)
c
c nn = branch
c ll = total angular oder of the mode
c nsta = id # of the station
c dw = eigenfrequency shift 
c-----------------------------------------------------------------------
      include 'parameter5.h'
c-----------------------------------------------------------------------
c declaration for radial integrals and ellipticity correction factor
c spheroidal as well as toroidal

      integer*4 jcom_h,ifanis_h,imode
      integer*4 nnll(0:maxn,0:maxl)
      real*4 tref_h,rn_h
      real*8 int_ker(maxmodes,0:kmax1)
      real*8 int_cru(maxmodes)
      real*8 wa_t(maxmodes),gvelo(maxmodes)
      character*80 model_h

      common /integ_const_h/model_h,tref_h,rn_h,jcom_h,ifanis_h,imode
      common /integ_h/int_ker,int_cru,wa_t,gvelo,nnll
c-----------------------------------------------------------------------
c declaration for M84 coefficients

      integer*4 icrust,ismax
      complex*16 CA(-nord:nord,0:nord,0:kmax)
      complex*16 CC(-nord:nord,0:nord,0:kmax)
      complex*16 crustCC(-nord:nord,0:nord)

      common /coeff_m84/CA,CC,crustCC,icrust,ismax
c-----------------------------------------------------------------------
c declaration for SH8-12 coefficients

      complex*16 manW(-nord:nord,0:nord,0:kmax1)
      complex*16 manL(-nord:nord,0:nord,0:kmax2)
      complex*16 manU(-nord:nord,0:nord,0:kmax3)
      complex*16 cruC(-nord:nord,0:nord)

      common /coeff_sh8/manW,manL,manU,cruC
c-----------------------------------------------------------------------
c declaration for angular parameters and averages over complex
c spherical harmonics

      real*8 dist(maxstat),azmp(maxstat),theta(maxstat)
      complex*16 avg_great(-nord:nord,0:nord,maxstat)
      complex*16 avg_diff(-nord:nord,0:nord,maxstat)

      common /ang_h/avg_great,avg_diff,dist,azmp,theta
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      integer*4 nn,ll,nsta
      real*8 dw,pi,radcon
      complex*16 sum  

      pi = 4.d0*datan(1.d0)
      radcon = pi/180.d0

cc...pointer for integrated kernels
      ip = nnll(nn,ll)

cc^^^^^^^^^^^^^^^^^^^    KMAX   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      if (model_h(1:3) .eq. 'M84') then

cc...model M84A
        if (icrust.eq.0) then
         dw = 0.d0
         do kk = 0,kmax
           sum = (0.d0,0.d0)
           do iss = 0,ismax
            do itt = -iss,iss 
             sum = sum + avg_great(itt,iss,nsta)*CA(itt,iss,kk)
            enddo
           enddo
           dw = dw + dble(sum)*int_ker(ip,kk)
         enddo    
cc...model M84C
        elseif (icrust.eq.1) then
         dw = 0.d0
         do kk = 0,kmax
           sum = (0.d0,0.d0)
           do iss = 0,ismax
             do itt = -iss,iss 
               sum = sum + avg_great(itt,iss,nsta)*CC(itt,iss,kk)
             enddo
           enddo
           dw = dw + dble(sum)*int_ker(ip,kk)
         enddo  
         sum = (0.d0,0.d0)
         do iss = 0,ismax
           do itt = -iss,iss 
             sum = sum + avg_great(itt,iss,nsta)*crustCC(itt,iss)
           enddo
         enddo 
         dw = dw + dble(sum)*int_cru(ip)
        endif
 
cc^^^^^^^^^^^^^^^^^^^    KMAX1  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      elseif (model_h(1:4).eq.'SH8W'.or.model_h(1:4).eq.'SH12') then

        dw = 0.d0
        do kk = 0,kmax1
          sum = (0.d0,0.d0)
          do iss = 0,ismax
            do itt = -iss,iss 
              sum = sum + avg_great(itt,iss,nsta)*manW(itt,iss,kk)
            enddo
          enddo
          dw = dw + dble(sum)*int_ker(ip,kk)
        enddo   
 
        sum = (0.d0,0.d0)
        do iss = 0,ismax
          do itt = -iss,iss 
            sum = sum + avg_great(itt,iss,nsta)*cruC(itt,iss)
          enddo
        enddo 
        dw = dw + dble(sum)*int_cru(ip)

cc^^^^^^^^^^^^^^^^^^^    KMAX2 KMAX3 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      elseif (model_h(1:4) .eq. 'SH8U') then

        dw = 0.d0
        do kk = 0,kmax2
          sum = (0.d0,0.d0)
          do iss = 0,ismax
            do itt = -iss,iss 
              sum = sum + avg_great(itt,iss,nsta)*manL(itt,iss,kk)
            enddo
          enddo
          dw = dw + dble(sum)*int_ker(ip,kk)
        enddo

        do kk = 0,kmax3
          sum = (0.d0,0.d0)
          do iss = 0,ismax
            do itt = -iss,iss 
              sum = sum + avg_great(itt,iss,nsta)*manU(itt,iss,kk)
            enddo
          enddo
          dw = dw + dble(sum)*int_ker(ip,kk+9)
        enddo
    
        sum = (0.d0,0.d0)
        do iss = 0,ismax
          do itt = -iss,iss 
            sum = sum + avg_great(itt,iss,nsta)*cruC(itt,iss)
          enddo
        enddo 
        dw = dw + dble(sum)*int_cru(ip)
      endif

cc...ellipticity correction
      dw = dw + wa_t(ip)*(1.d0-3.d0*dcos(theta(nsta)*radcon)**2)

      return
      end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine calc_dd(nn,ll,nsta,dd)
c-----------------------------------------------------------------------
c computes the apparent shift in epicentral distance due to asphericity
c (lateral heterogeneity and ellipticity of figure) 
c-----------------------------------------------------------------------
      include 'parameter5.h'
c-----------------------------------------------------------------------
c declaration for radial integrals and ellipticity correction factor
c spheroidal as well as toroidal

      integer*4 jcom_h,ifanis_h,imode
      integer*4 nnll(0:maxn,0:maxl)
      real*4 tref_h,rn_h
      real*8 int_ker(maxmodes,0:kmax1)
      real*8 int_cru(maxmodes)
      real*8 wa_t(maxmodes),gvelo(maxmodes)
      character*80 model_h

      common /integ_const_h/model_h,tref_h,rn_h,jcom_h,ifanis_h,imode
      common /integ_h/int_ker,int_cru,wa_t,gvelo,nnll
c-----------------------------------------------------------------------
c declaration for M84 coefficients

      integer*4 icrust,ismax
      complex*16 CA(-nord:nord,0:nord,0:kmax)
      complex*16 CC(-nord:nord,0:nord,0:kmax)
      complex*16 crustCC(-nord:nord,0:nord)

      common /coeff_m84/CA,CC,crustCC,icrust,ismax
c-----------------------------------------------------------------------
c declaration for SH8-12 coefficients

      complex*16 manW(-nord:nord,0:nord,0:kmax1)
      complex*16 manL(-nord:nord,0:nord,0:kmax2)
      complex*16 manU(-nord:nord,0:nord,0:kmax3)
      complex*16 cruC(-nord:nord,0:nord)

      common /coeff_sh8/manW,manL,manU,cruC
c-----------------------------------------------------------------------
c declaration for angular parameters and averages over complex
c spherical harmonics

      real*8 dist(maxstat),azmp(maxstat),theta(maxstat)
      complex*16 avg_great(-nord:nord,0:nord,maxstat)
      complex*16 avg_diff(-nord:nord,0:nord,maxstat)

      common /ang_h/avg_great,avg_diff,dist,azmp,theta
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      integer*4 nn,ll,nsta
      real*8 dd,pi,radcon,fact,fact2
      real*8 delt,azim,thet,rnd,gv
      complex*16 sum  

      pi = 4.d0*datan(1.d0)
      radcon = pi/180.d0

cc...pointer for integrated radial kernels
      ip = nnll(nn,ll)

      delt = dist(nsta)*radcon
      azim = azmp(nsta)*radcon
      thet = theta(nsta)*radcon
      rnd = dble(rn_h)
      gv = gvelo(ip)*1000.d0

      if (gv.eq.0.d0) then
        dd = 0.d0
        return
      endif

cc^^^^^^^^^^^^^^^^^^^    KMAX   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      if (model_h(1:3) .eq. 'M84') then

        if (icrust.eq.0) then
cc...model M84A mantle only
          dd = 0.d0
          do kk = 0,kmax
            sum = (0.d0,0.d0)
            do iss = 0,ismax
              do itt = -iss,iss 
                sum = sum + avg_diff(itt,iss,nsta)*CA(itt,iss,kk)
              enddo
            enddo
            dd = dd + dble(sum)*int_ker(ip,kk)
          enddo    

        elseif (icrust.eq.1) then
cc...model M84C mantle at first
          dd = 0.d0
          do kk = 0,kmax
            sum = (0.d0,0.d0)
            do iss = 0,ismax
              do itt = -iss,iss 
                sum = sum + avg_diff(itt,iss,nsta)*CC(itt,iss,kk)
              enddo
            enddo
            dd = dd + dble(sum)*int_ker(ip,kk)
          enddo
cc...model M84C then crust 
          sum = (0.d0,0.d0)
          do iss = 0,ismax
            do itt = -iss,iss 
              sum = sum + avg_diff(itt,iss,nsta)*crustCC(itt,iss)
            enddo
          enddo 
          dd = dd + dble(sum)*int_cru(ip)
        endif

cc^^^^^^^^^^^^^^^^^^^    KMAX1  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      elseif (model_h(1:4).eq.'SH8W'.or.model_h(1:4).eq.'SH12') then

        dd = 0.d0
        do kk = 0,kmax1
          sum = (0.d0,0.d0)
          do iss = 0,ismax
            do itt = -iss,iss 
              sum = sum + avg_diff(itt,iss,nsta)*manW(itt,iss,kk)
            enddo
          enddo
          dd = dd + dble(sum)*int_ker(ip,kk)
        enddo   
 
        sum = (0.d0,0.d0)
        do iss = 0,ismax
          do itt = -iss,iss 
            sum = sum + avg_diff(itt,iss,nsta)*cruC(itt,iss)
          enddo
        enddo 
        dd = dd + dble(sum)*int_cru(ip)

cc^^^^^^^^^^^^^^^^^^^    KMAX2 KMAX3 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      elseif (model_h(1:4) .eq. 'SH8U') then

        dd = 0.d0
        do kk = 0,kmax2
          sum = (0.d0,0.d0)
          do iss = 0,ismax
            do itt = -iss,iss 
              sum = sum + avg_diff(itt,iss,nsta)*manL(itt,iss,kk)
            enddo
          enddo
          dd = dd + dble(sum)*int_ker(ip,kk)
        enddo

        do kk = 0,kmax3
          sum = (0.d0,0.d0)
          do iss = 0,ismax
            do itt = -iss,iss 
              sum = sum + avg_diff(itt,iss,nsta)*manU(itt,iss,kk)
            enddo
          enddo
          dd = dd + dble(sum)*int_ker(ip,kk+9)
        enddo
    
        sum = (0.d0,0.d0)
        do iss = 0,ismax
          do itt = -iss,iss 
            sum = sum + avg_diff(itt,iss,nsta)*cruC(itt,iss)
          enddo
        enddo 
        dd = dd + dble(sum)*int_cru(ip)

      endif

cc...value of dd
      fact = rnd/((dble(ll)+0.5d0)*gv)
      dd = dd*fact*delt
cc...ellipticity correction
      fact2 = -3.d0*fact*dsin(delt)*dcos(2.d0*azim)*dsin(thet)**2
      dd = dd + wa_t(ip)*fact2

      return
      end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine ang_param(elat,elon,slat,slon,
     &                    delta,azmp,theta,avg_great,avg_diff,ismax)
c-----------------------------------------------------------------------
c computes all angular quantities and averages over great-circle 
c and difference between great-circle and minor-arc averages 
c needed to compute delta-omega and delta-delta
c
c   elat,elon = event lat, lon (geographic coord)
c   slat,slon = station lat,lon (geographic coord)
c
c   delta = epicentral distance
c   azmp  = azimuth from north of midpoint of minor-arc between event 
c           and station from pole to great circle
c   theta = colat of great circle pole
c   avg_great = average of complex spherical harmonics up to degree 
c               ismax over great circle
c   avg_diff  = difference between great circle and minor arc average
c----------------------------------------------------------------------                  
      include 'parameter5.h'

cc...global variables    
      real*8 elat,elon,slat,slon
      real*8 delta,azmp,theta,dist_a
      complex*16 avg_great(-nord:nord,0:nord)
      complex*16 avg_diff(-nord:nord,0:nord)

cc...local variable
      real*8 eglat,eglon,sglat,sglon
      real*8 polelat,polelon,dlat,dlon
      complex*16 avg_g(0:nord,-nord:nord)
      complex*16 avg_m(0:nord,-nord:nord)

c-----------------------------------------------------------------------
cc...converts to geocentric coordinates
      call convert_geoc(elat,elon,eglat,eglon)
      call convert_geoc(slat,slon,sglat,sglon)

cc...epicentral distance
      delta = dist_a(eglat,eglon,sglat,sglon)

cc...pole to great circle
      call pole(eglat,eglon,sglat,sglon,polelat,polelon)
      theta = 90.d0-polelat

cc...azimuth pole to midpoint
      call midpt(eglat,eglon,sglat,sglon,dlat,dlon)
      call azi(polelat,polelon,dlat,dlon,azmp)

cc...great circle and minor arc averages
      call avg_shm('great ',eglat,eglon,sglat,sglon,ismax,avg_g)
      call avg_shm('minor ',eglat,eglon,sglat,sglon,ismax,avg_m)

cc...compute avg_diff and swap storage space for use in main program
      do iss = 0, ismax
        do itt=-iss,iss
            avg_great(itt,iss)= avg_g(iss,itt)
            avg_diff(itt,iss) = avg_g(iss,itt)-avg_m(iss,itt)
        enddo
      enddo

      return
      end



