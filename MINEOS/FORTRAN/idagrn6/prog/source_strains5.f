c***********************************************************************
c                   SOURCE STRAINS 4.0
c
c contains the modified subroutine source_strains and the version of it
c that implements a given laterally heterogenous model and another
c version that is station sensitive
c 2 mode tables : 1) source strains  2) propagation and receiver
c
c      subroutine source_strains_h2(source_file,d0,num_modes,ista,lread)
c
c***********************************************************************

       subroutine source_strains_h2(source_file,d0,num_modes,ista,
     &                              lread,lhetero)
c
c     this subroutine handles the heavy io of the source mode file
c     it assumes a reformatted mineos output file
c     d0 is the depth of the event in kilometers
c
c     02/13/90 - ### introduced sign change in synthetics ####
c     04/27/90 - added branch option
c     07/02/90 - new mode table header format
c     01/20/92 - added aspherical representation,modified common/modes/ 
c                to common/modes_h/ and common/vecnl/ to 
c                common/vecnl_h/

      implicit real*8 (a-h,o-z)
c
      include 'parameter5.h'
c
      character*256 source_file(3), m_file, b_file
      character*1 lcomp
c
      integer*4 ip_s(0:maxl,0:maxn),ip_r(0:maxl,0:maxn)
c
      real*4 d0, dt
      real*4 abuf(6), r1
      real*4 w, ql, alpha, phi
      real*4 f0(maxmodes),f1(maxmodes),f(4,3,maxmodes)
      real*4 f0_temp(maxmodes),f1_temp(maxmodes),f_temp(4,3,maxmodes)
      real*4 x1, x2
      real*4 wmhz
      real*4 buff(maxmodes,6)
      real*4 sort(0:maxl,0:maxn)
      real*4 ak(maxmodes)
      real*4 dum1, dum2, dum3, dum4, dum5, tref
      real*4 digar(maxstat)
      real*4 qalpha(10000), qbeta(10000), rq(10000)
c
      real*8 fl, fl1(0:maxl), fl3(0:maxl), fl4(0:maxl)
      real*8 ddr, ddi, qinv
      real*8 vecnl,vecnl_h, pi, wref
      real*8 rs, r0
      real*8 wd, wdb, wdt, decay, qdb
      real*8 u, v, up, vp
      real*8 e14, e01, e02, e22, au, av
      real*8 wr, wp, e15, e26, aw
      real*8 u0, v0, w0
      real*8 au1, au2, au3, av1, av2, av3
      real*8 scale1, scale2, scale3
      real*8 radian, mhz
      real*8 dw,wddw,wddwsq    

c
      integer*4 scnl, nn, ll
      integer*4 jcom, modes, llmin, llmax, knots
      integer*4 jcom2,jcom3,jcom4
      integer*4 nbr
      integer*4 nptsar(maxstat)
c
      logical ocean, lmask, ldisp, lida, lbranch, lpart
      logical lread,lhetero

c common block for crust5.1/dbdb options

      logical lcrust51, ldbdb5
      character c51opt*2, cdbopt*2

      common /coption/ lcrust51,ldbdb5,c51opt,cdbopt

c other commons
      common/ahdr/jcom,knots
      common/scl/x1,r1,x2,f
      common/vecnl_h/scnl(maxmodes,2),vecnl(maxmodes,4),
     &               vecnl_h(maxmodes,5:6)
      common/mode_h/ddr(maxmodes),ddi(maxmodes),qinv(maxmodes)
      common/mode_hdr/llmax,llmin
      common/buf/nn, ll, w, ql, alpha, phi
      common/limits/dt,wmhz,lcomp,ldisp,lida
      common/mask/ m_file, lmask
      common/branch/ b_file, lbranch, nbr
      common /stasen/digar,nptsar

      common/ipoint/ip_s,ip_r
      common/tempor/f_temp,f0_temp,f1_temp
      common/noname/buff
c
      equivalence (nn,abuf)
c
      data pi/3.14159265359d0/
      data rn /6371000.d0/
      data bigg /6.6723d-11/
      data g0 /9.780317d0/
      data rhobar/5515.d0/
      data ocean /.false./
 
c
      mhz = 500.d0/pi
      radian = 1.0/mhz
      lpart = .false.

cc...some constants
cc
      scale1 = 1.d0/(rn*dsqrt(rn*pi*bigg)*rhobar)
      scale2 = dsqrt(pi*bigg/rn)
      scale3 = 1.d20
      wref = dble(wmhz) * radian
      rs = rn - dble(d0)*1000.d0
      r0 = rs/rn
      r1 = sngl(r0)
      io1 = 1
      io2 = 2
      io3 = 3

cc...loads source and receiver eigenfunctions and pointers file
cc   only first time through

      if (.not.lread) then

        call load_source(io1,source_file(1),jcom,modes,llmin,llmax)
        call load_receiver(io2,source_file(3),jcom2,modes2,
     &                     llmin2,llmax2)

cc...open propagation file,get reference frequency
cc   the propagation file controls modes,llmin and llmax

        call kblnk(source_file(2),k)
        open(io3,form='unformatted',file=source_file(2)(1:k)//'_hdr',
     &       access = 'sequential')
        read(io3) jcom4, modes,dum4,dum5, llmin, llmax
        read(io3) dum1, dum2, dum3, ifanis, tref
        read(io3)
        read(io3)
        read(io3) nqk, (rq(i),i=1,nqk), (qalpha(i),i=1,nqk), (qbeta(i),i
     &       =1,nqk)
c        do i=1,nqk
c           write(*,*) i, rq(i), qalpha(i), qbeta(i)
c        enddo
              
cc...scale1 contains the normalization for u,v,w, and their derivatives
cc   scale2 contains the normalization for phi and phip

        do lc = llmin, llmax
          fl = dble(lc)
          fl1(lc) = fl + 1.d0
          fl3(lc) = fl*fl1(lc)
          fl4(lc) = dsqrt(fl3(lc))
        enddo

cc...load mask or branch file
        jcom3 = jcom
        if (lmask.or.lbranch) then
          call load_bm(modes,jcom3,sort,lpart)
        endif

cc...test jcom
        if (jcom.ne.jcom2.and.jcom.ne.jcom4) then
          print*,' Mode tables for source and receiver incompatible'
          print*,' jcom ',jcom,jcom2,jcom4
          stop
        endif
        if (jcom.ne.jcom3) then
          print*,' WARNING'
          print*,' Mask or branch file incompatible with mode table'
        endif

cc...read propagation part
cc...check availability of modes from source and receiver mode tables
cc...if this is a partial sum, form a reduced mode set
c
        num_modes = 0
        if (lpart) then
          print*,' forming partial sums'
          do ii = 1, modes
            read(io3) (abuf(i),i=1,6)

            if (sort(ll,nn) .ne. 0.) then
            if (ip_s(ll,nn) .ne. 0 ) then
            if (ip_r(ll,nn) .ne. 0 ) then
              num_modes = num_modes + 1
              if (num_modes .gt. ii) then
                print*,' problem with partial sum: ', num_modes, ii
              endif
cc   strain pointer
              ip = ip_s(ll,nn)
              do jj = 1, 3
                do kk = 1, 4
                  f(kk,jj,num_modes) = f_temp(kk,jj,ip)
                end do
              end do
cc   receiver pointer
              ip = ip_r(ll,nn)
              f0(num_modes) = f0_temp(ip)
              f1(num_modes) = f1_temp(ip)
              ak(num_modes) = sort(ll,nn)
              do mm = 1, 6
                buff(num_modes, mm) = abuf(mm)
              end do
            endif
            endif
            endif
          end do
        else
          print*,' forming complete synthetic'
          do ii = 1, modes
            read(io3) (abuf(i),i=1,6)
c            if (ll.eq.20) write(*,*) nn, ll, w
            if (ip_s(ll,nn) .ne. 0 ) then
            if (ip_r(ll,nn) .ne. 0 ) then
              num_modes = num_modes + 1
              if (num_modes .gt. ii) then
                print*,' problem with partial sum: ', num_modes, ii
             endif
cc   strain pointer
             ip = ip_s(ll,nn)
             do jj = 1, 3
               do kk = 1, 4
                 f(kk,jj,num_modes) = f_temp(kk,jj,ip)
               end do
             end do
cc   receiver pointer
             ip = ip_r(ll,nn)
             f0(num_modes) = f0_temp(ip)
             f1(num_modes) = f1_temp(ip)
             ak(num_modes) = 1.0
             do mm = 1, 6
               buff(num_modes, mm) = abuf(mm)
             end do
           endif
           endif
         end do

        end if
        close(io3)

c
c   cubic interpolation to source depth
c
        call cubic(jcom,num_modes)

      endif

      lread = .true.
c
c   here begins a vast if statement
c   test for mode type and then loop over the modes
c   check to see if Q perturbation to eigenfrequency 
c   was calculated in mineos
cc..................................................................
cc
cc    SPHEROIDAL STRAINS 
cc..................................................................

      print*,' '
      if (jcom .eq. 3) then 
        print*,' spheroidal mode strains'
        do nm = 1, num_modes
            do ii = 1, 6
              abuf(ii) = buff(nm,ii)
            end do
            wdb = dble(w)
            qdb = dble(ql)

cc{{ pfi compute for each station the delta-w
cc   jbg add in alternative call for crust/dbdb5 corrections

            if (lhetero) then
               call calc_dw(nn,ll,ista,dw)
            elseif (lcrust51) then
               call dw_c51(nn,ll,wdb,ista,dw)
            elseif (ldbdb5) then
               call dw_dbd(nn,ll,wdb,ista,dw)
c               if(nn.eq.0 .and. mod(ll,100).eq.0) then
c                 print*,nn,ll,wdb,dw
c               endif
            else
               dw = 0d0
            endif

            if (tref .le. 0.) then
              wd = wdb*(1.d0 + (1.d0/(qdb*pi))*dlog(wdb/wref))
            else
              wd = wdb
            endif
          
            qinv(nm) = 0.5d0/qdb
            wddw   = wd + dw
            wddwsq = wddw*wddw
            dt     = digar(ista)
            wdt    = wddw*dble(dt)
            decay  = dexp(-wdt*qinv(nm))
            ddr(nm)= decay*dcos(wdt)
            ddi(nm)= decay*dsin(wdt)
cc }}
c
c           up and vp are normalized by 1/rn in addition to scale1
c
            up = dble(f(2,2,nm))/rn
            vp = dble(f(4,2,nm))/rn
c
c           dividing u and v here by rs simplifies the algebra for the
c             strains
c
            v = dble(f(3,2,nm))/rs
            u = dble(f(1,2,nm))/rs
c
c           strains a la Gilbert and Dziewonski, 1978
c
            e01 = up
            e02 = u - .5d0 * fl3(ll) * v
            e14 = vp + u - v
            e22 = 2.d0 * v
c
c           au is the vertical acceleration scalar
c           av is the horizontal acceleration scalar
c           refer to Gilbert, 1980
c
            u0 = dble(f0(nm))*scale1
            v0 = dble(f1(nm))*scale1
            phi0 = dble(phi)*scale2
c
c           note that there is a different normalization
c             between G&D (1975) and Gilbert (1980)
c             one could be hopeless confused by this
c             hopefully, I am not
c
            au2 = 2.d0*g0*u0/rn
            au3 = fl1(ll)*phi0/rn

            av2 = -g0*u0/rn
            av3 = -phi0/rn

            au1 = wddwsq*u0
            au  = au1 + au2 + au3
            av1 = wddwsq*v0
            av  = av1 + av2 + av3
            if (ldisp) then
              au = -u0
              av = -v0
            endif
            vecnl_h(nm,5) = au*scale3*ak(nm)
            vecnl_h(nm,6) = av*scale3*ak(nm)

            scnl(nm,1)  = nn
            scnl(nm,2)  = ll
            vecnl(nm,1) = e01*scale1
            vecnl(nm,2) = e02*scale1
            vecnl(nm,3) = e14*scale1
            vecnl(nm,4) = e22*scale1
c
        end do

cc..................................................................
cc
cc    TOROIDAL STRAINS 
cc..................................................................

      elseif (jcom .eq. 2) then
        print*,' toroidal mode strains'
        do nm = 1, num_modes
            do ii = 1, 6
              abuf(ii) = buff(nm,ii)
            end do
            wdb = dble(w)
            qdb = dble(ql)

cc  pfi compute for each station the delta-w
c            print*,'here x....'
            if (lhetero) then
               call calc_dw(nn,ll,ista,dw)
            elseif (lcrust51) then
c               print*, 'here lcrust51...a', ista
               call dw_c51(nn,ll,wdb,ista,dw)
c               print*, 'here lcrust51...b'
c               if(nn.eq.0 .and. mod(ll,100).eq.0) then
c     print*,nn,ll,wdb,dw
c               endif
            elseif (ldbdb5) then
               call dw_dbd(nn,ll,wdb,ista,dw)
            else
               dw = 0d0
            endif
c            print*,'here y....'
            if (tref .le. 0.) then
               wd = wdb*(1.d0 + (1.d0/(qdb*pi))*dlog(wdb/wref))
            else
               wd = wdb
            endif

            qinv(nm) = 0.5d0/qdb
            wddw   = wd + dw
            wddwsq = wddw*wddw
            dt     = digar(ista)
            wdt    = wddw*dble(dt)
            decay  = dexp(-wdt*qinv(nm))
            ddr(nm)= decay*dcos(wdt)
            ddi(nm)= decay*dsin(wdt)
c
c           wp is normalized by 1/rn in addition to scale1
c
            wp = dble(f(2,2,nm))/rn
c
c           dividing wr by rs simplifies the algebra of the strains
c
            wr = dble(f(1,2,nm))/rs
c
c           aw is the horizontal acceleration scalar
c           refer to Gilbert, 1980
c
            w0 = dble(f0(nm))*scale1

            e15 = -(wp - wr)
            e26 = -4.d0*wr

            aw = wddwsq*w0
            if (ldisp) then
              aw = -w0
            endif
            vecnl_h(nm,5) = aw*scale3*ak(nm)
            vecnl_h(nm,6) = 0.d0

            scnl(nm,1)  = nn
            scnl(nm,2)  = ll
            vecnl(nm,1) = e15*scale1
            vecnl(nm,2) = e26*scale1
            vecnl(nm,3) = 0.d0
            vecnl(nm,4) = 0.d0

        end do
        print*,'End toroidal mode strains'
cc..................................................................
cc
cc    RADIAL STRAINS 
cc..................................................................

      elseif (jcom .eq. 1) then
        print*,' radial mode strains'
        do nm = 1, num_modes
            do ii = 1, 6
              abuf(ii) = buff(nm,ii)
            end do
            wdb = dble(w)
            qdb = dble(ql)

cc{{ pfi compute for each station the delta-w
cc does not make sense, since first order perturbation theory
cc is not valid for radial modes, implemented for reasons of
cc symmetry

            if (lhetero) then
               call calc_dw(nn,ll,ista,dw)
            elseif (lcrust51) then
               call dw_c51(nn,ll,wdb,ista,dw)
            elseif (ldbdb5) then
               call dw_dbd(nn,ll,wdb,ista,dw)
            else
               dw = 0d0
            endif

            if (tref .le. 0.) then
              wd = wdb*(1.d0 + (1.d0/(qdb*pi))*dlog(wdb/wref))
            else
              wd = wdb
            endif

            qinv(nm) = 0.5d0/qdb
            wddw   = wd + dw
            wddwsq = wddw*wddw
            dt     = digar(ista)
            wdt    = wddw*dble(dt)
            decay  = dexp(-wdt*qinv(nm))
            ddr(nm)= decay*dcos(wdt)
            ddi(nm)= decay*dsin(wdt)
cc }}
c
c           up is normalized by 1/rn in addition to scale1
c
            up = dble(f(2,2,nm))/rn
c
c           dividing u by rs here simplifies the algebra of the strains
c
            u = dble(f(1,2,nm))/rs
c
c           au is the vertical acceleration scalar
c           refer to Gilbert, 1980
c
            u0 = dble(f0(nm))*scale1
            e01 = up
            e02 = u 
            au2 = 2.d0*g0*u0/rn

            au1 = wddwsq*u0
            au  = au1 + au2
            if (ldisp) then
              au = -u0
            endif
            vecnl_h(nm,5) = au*scale3*ak(nm)
            vecnl_h(nm,6) = 0.d0

            scnl(nm,1)  = nn
            scnl(nm,2)  = ll
            vecnl(nm,1) = e01*scale1
            vecnl(nm,2) = e02*scale1
            vecnl(nm,3) = 0.d0
            vecnl(nm,4) = 0.d0

        end do
      endif
      if (lida) then
        do nm = 1, num_modes
          vecnl_h(nm,5) = -vecnl_h(nm,5)
          vecnl_h(nm,6) = -vecnl_h(nm,6)
        end do
      endif
c
      return
      end

