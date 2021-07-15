      program mineos_qcorrectphv
c
c     mineos_q performs the modal q integration on a mineos output file
c     it leaves the eigenfunctions untouched, outputting an ascii file
c     which contains the mode header, including the calculated Q
c
c     will integrate both spheroidal and toroidal modes
c
c     mineos_q expects a raw mineos output file containing:
c       two record header
c       then 1 record for every mode, depending on type
c        spheroidal modes: nn, ll, w, q, gv, u(1:knots), up(1:knots),
c                          v(1:knots), vp(1:knots), 
c                          phi(1:knots), phip(1:knots)
c        radial modes:     nn, ll, w, q, gv, u(1:knots), up(1:knots)
c        toroidal modes:   nn, ll, w, q, gv, w(1:knots), wp(1:knots)
c    
c     03/24/87 spheroidal, radial, and toroidal calculation checked 
c              against q values from Sailor and Dziewonski, 1978       
c     11/04/88 rechecked calculations - error ~1.e-3 for certain modes
c     07/02/90 test added for Q calculation in mineos
c     01/09/91 check for overlapping eigenfrequency added (thanks, Jim)
c     04/22/91 minor bugs cleaned up (concerning anisotropic prem, 
c              overlapping eigenfrequency test,gravity scalar on ocean 
c              floor for oceanic model) pp
c     07/20/91 updated to handle anisotropic models.  Calculated q is 
c              still equivelent isotropic Q,as in PREM
c              Tested by comparing Q results from an isotropic model 
c              run here as both isotropic and anisotropic.
c              Occasionally Q values differ by 2 part in 10**8 
c              (i.e. 272.43210 vs 272.43212), which seems OK
c        01/92 w,h1,gv from mineos are now real*8 and nvec 
c              is increased by 3
c        04/92 option to overwrite modal Q's if already calculated
c              in mineos
c
c     12/10/14 adding 3 output colloums for physicail dispersion correction
c              for phase velocity. 
c              original: write(3,120) nn,ll,w,qq,phi,cv,gv
c              updated:  write(3,120) nn,ll,w,qq,phi,cv,gv,cvq,Tq,T
c              the equation is got from 
c              pylin.patty 

      include 'parameter.h'
c
      real*8 intg(nknot),x(nknot),dr(nknot),t
      real*8 u, up, v , vp, wl, wp, f
      real*8 ff, mm, kk, scale, pi, rn, rj, third
      real*8 fl1(0:maxll), fl2(0:maxll),fl3(0:maxll)
      real*8 mu(nknot),muq(nknot)
      real*8 kappa(nknot),kapq(nknot)
      
      real*8 wd,h1d,gvd
c
      real*4 w,h1,gv
      real*4 bulkq(nknot),shearq(nknot),rq(nknot)
      real*4 mbq(nknot),bbq(nknot),msq(nknot),bsq(nknot)
      real*4 rad(nknot)
      real*4 interple, wlast
      real*4 abuf(maxbyte5+3),radius(nknot)
      real*4 dn(nknot),alpha(nknot),beta(nknot),qa(nknot),qb(nknot)
      real*4 alphah,betah,eta
      real*4 temp(nknot,5)
      real*4 old(maxn, maxll)
      real*8 radian, mhz,wwd,wref
    
c
      integer*4 nvec
      integer*4 nn, ll, knot, ifirst, obs
c
      character*1 comp(3),ans,dtest,set
      character*256 q_file,m_file,o_file
c
      logical lq
c
      common/ablk/nn,ll,wd,h1d,gvd,buf(nknot5)
      equivalence(nn,abuf)
c
      data comp /'S','T','S'/
      data pi/3.1415926535d0/
      data rn/6371000.d0/
      data old /maxold*0.0/
      data lq /.false./
      data ifirst /0/
      data qflag /999999.0/
c
      third = 1.d0/3.d0
      fthird = 4.d0/3.d0
c
c     open and read q file
c
      print*,'enter name of q model'
      read(*,'(a)') q_file
      open(unit=2,file=q_file)
      read(2,*) nq
      if (nq .gt. nknot) then
        print*,' array limits of q exceeded'
        print*,' reset nqknot'
        stop
      endif
      read(2,*)(rq(i),shearq(i),bulkq(i),i=1,nq)
      print*,' Q Model: ', q_file
      print*,'    Radius      Shear       Bulk' 
      write(*,130)(rq(i),shearq(i),bulkq(i),i=1,nq)
      close(2)
c
c     convert q radii to meters
c
      do 10 j=1,nq
        rq(j)=rq(j)*1000.
 10   continue
c
c     linear fit the kappa q's
c 
      call interpol(1,nq,rq,bulkq,mbq,bbq)
      call interpol(1,nq,rq,shearq,msq,bsq)
c
c     open ascii output file
c
      print*,' Enter name of ascii output file'
      read(*,'(a)')o_file
      open(unit=3, file=o_file)
c
c     print out q model to ascii output file
c
      write(3,129) nq
      write(3,130)(rq(i),shearq(i),bulkq(i),i=1,nq)
c
c     open mode file and read header records - 
c     beginning of loop over mode files
c
  5   print*,' Enter name of input mode file'
      read(*,'(a)')m_file 
      if (m_file .eq. ' ') then
        close(3)
        stop
      endif
      open(unit=2,file=m_file,form='unformatted',access='sequential')
csun     +    recl=36000)
      read(2) jcom, wmin, wmax, lmin, lmax, wgrav
      read(2) knot, nic, noc, ifanis, tref, (rad(i), i=1,knot),
     &    (dn(i),i=1,knot), (alpha(i),i=1,knot), (beta(i),i=1,knot),
     &    ((temp(i,j),i=1, knot), j = 1, 5)

      ifirst = ifirst + 1
c
c TESTS  1) array dimensions
c
      if (lmax.gt.maxl) then
        print*,' ****************************************************'
        print*,' WARNING lmax > maxl  !   '
        print*,' mode file header(your guess in mineos) lmax = ',lmax
        print*,' from parameter statement maxl = ',maxl
        print*
        print*,' Do you want to set lmax = maxl < y/n > ?'
        read(*,'(a)') set
        if (set.eq.'y') then
           lmax=maxl
           print*,' Overwrite lmax in mode file with maxl < y/n > ??'
           read(*,'(a)') dtest
           if (dtest.eq.'y') then
             close(2)
             open(unit=2,file=m_file,form='unformatted',
     &            access='sequential')
             write(2) jcom, wmin, wmax, lmin, maxl, wgrav
             read(2)
           endif
        else
          print*,' Sorry, the program has to abort'
          print*,' Please modify the parameter statement'
          close(2)
          close(3)
          stop
        endif
        print*,' ***************************************************'
      endif
c
c TESTS  2) q calculation in mineos
c        3) oceanic model
c
      if (ifirst .eq. 1) then
        if (tref .le. 0.) then
c....     check if a Q-model was input in mineos -> 
c         calc. of modal Q's done there
          lq = .true.
          do i = 1,knot
            qa(i) = temp(i,1)
            qb(i) = temp(i,2)
            if (qa(i) .gt. 0.e0 .and. qa(i) .lt. qflag) lq = .false.
            if (qb(i) .gt. 0.e0 .and. qb(i) .lt. qflag) lq = .false.
          enddo
        endif
c
c....   give user a chance to override and calculate Q's anyway 
c       iff no dispersion correction had been made in mineos!
        if (.not. lq .and. tref .le. 0) then
          print*,' '
          print*,'modal Q already calculated in mineos'
          print*,'calculate Q again here < Y/N > ?'
          read(*,'(a1)') ans
          print*,' '
         if (ans .eq. 'y' .or. ans .eq. 'Y') lq = .true.
        endif
c
c....   if this is an oceanic model we need to know 
c       the number of knots under water (obs)
        obs = 0
        if (beta(knot) .eq. 0.e0) then
          do i = knot,1,-1
            if (beta(i) .ne. 0.e0) then
              obs = knot - i
              go to 15
            endif
          enddo
   15     continue
        endif
      endif
cpp <<
c
      print*, ' mode type = ',jcom,' lmin = ',lmin,' lmax = ', lmax
      print*, ' wmin = ', wmin,' wmax = ', wmax
      print*, ' knots = ',knot

      if (lq) then
        print*, ' '
        print*,' Calculating Q'
        print*, ' '
      else
        print*, ' '
        print*,' Saving Q from mineos'
        print*, ' '
      endif
      nocor = knot - noc
c
c     zero knot+1 radius BEFORE resetting knots for toroidal case
c        jbg 6/92
c
      rad(knot+1) = 0.0
      rad(1)=1.0
c
      if (jcom .eq. 2) then
        knot = nocor
      endif
ccc
ccc   big test to calculate Q or save mineos values
ccc
      if (lq) then
c
c        rad(knot+1)=0.0
c        rad(1)=1.0
c
c       DEBUGGING!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c       ifanis=0
c
c       set up intergration parameters and constants for 
c       spheroidal and radial modes
c
        if (jcom .ne. 2) then
          rl = rad(1)
          do j = 1, knot
             radius(j) = rad(j)
             muq(j) = dble(interple(1,nq,rq,rad(j),rl,shearq,msq,bsq))
             kapq(j) = dble(interple(1,nq,rq,rad(j),rl,bulkq,mbq,bbq))
             if (ifanis.eq.1) then
               alphah=temp(j,3)
               betah=temp(j,4)
               eta=temp(j,5)
               mu(j) = dble(1./15.*dn(j) * ((1.-2.*eta)*alphah**2. + 
     &                 alpha(j)**2. + 5.*betah**2 +
     &                 (6.+4.*eta)*beta(j)**2.))
               kappa(j) = dble(1./9.*dn(j) * ((4.+4.*eta)*alphah**2.+
     &                    alpha(j)**2. - 8.*eta*beta(j)**2. -
     &                    4.*betah**2.)) 
             else
               mu(j) = dble(beta(j)**2*dn(j))
               kappa(j) = dble((alpha(j)**2 - fthird*beta(j)**2)*dn(j))
             end if
             dr(j) = rad(j+1) - rad(j)
             rl = rad(j)
          end do
c
c       set up integration parameters and constants for toroidal modes
c
        else
          rl = rad(noc + 1)
          do j = 1, knot
             radius(j) = rad(noc + j)
             muq(j)=dble(interple(1,nq,rq,rad(noc+j),rl,shearq,msq,bsq))
             if (ifanis.eq.1) then
               alphah=temp(noc+j,3)
               betah=temp(noc+j,4)
               eta=temp(noc+j,5)
c               print*,rad(noc+j),dn(noc+j),alpha(noc+j),beta(noc+j),
c     &                alphah,betah,eta
               mu(j) = dble(1./15.*dn(noc+j)*((1.-2.*eta)*alphah**2. + 
     &                 alpha(noc+j)**2. + 5.*betah**2 +
     &                 (6.+4.*eta)*beta(noc+j)**2.))
             else
               mu(j) = dble(dn(noc +j) * beta(noc +j)**2)
             end if
             dr(j) = rad(noc + j + 1) - rad(noc + j)
c             if (dr(j).lt.0.) then
c                print*,j,dr(j),rad(noc + j + 1),rad(noc + j)
c             end if
             rl = rad(noc + j)
          end do
        endif
        print*,' knots in integration: ', knot
ccc
      endif
ccc
c
c     set up integration parameters
c
      do 20 j = lmin, lmax
        f = dble(j)
        fl1(j) = f * (f + 1.d0)
        fl2(j) = dsqrt(fl1(j))
        fl3(j) = (f - 1.d0) * (f + 2.d0)
 20   continue
c
c     scale is the normalization for the eigenfunctions - 
c     checked and rechecked 11/5/88
c
      scale = 1.d0/(rn*dsqrt(rn*pi*6.6723d-11)*5515.d0)
c
      if (jcom .eq. 3) then
        nvec=5*knot+8
      else
        nvec=2*knot+8
      endif
      k4 = 4*knot
      k3 = 3*knot
      k2 = 2*knot
      k1 = knot-1
c
c     begin reading mode files
c
 30   continue
      read(2,end=999)(abuf(i),i=1,nvec)
c      print*,'ll= ',ll,' nn= ',nn,' w= ',wd,' h1= ',h1d,' gv= ',gvd
c
c     look for spurious mode
c
c.... put w,q,gv in real*4 array
      w  = sngl(wd)
      h1 = sngl(h1d)
      gv = sngl(gvd)
c      print*,'ll= ',ll,' nn= ',nn,' w= ',w,' h1= ',h1,' gv= ',gv
      if (nn .lt. 0) then
        print *, ' apparent error in mode calculation'
        print *, nn,' ',comp(jcom),' ',ll
        if (nn.lt.-1) then
          print*,'n=-10: skipping next mode as well'
          read(2,end=999)
        end if
        go to 30 
      else if (nn .gt. 0) then
        wlast = old(nn,ll+1)
        if (abs(w - wlast) .lt. 1.e-7*w) then
          print*,' skipping overlapping eigenfrequency'
          print *, nn,' ',comp(jcom),' ',ll
          go to 30 
        endif 
      endif
c
c     check to see if this mode has been read already
c
      if (old(nn+1,ll+1) .ne. w) then
        old(nn+1,ll+1) = w
ccc
        if (lq) then        
c
c         equations for Q calculation taken from p. 266-267 
c         of Woodhouse, 1980
c
c         calculate compressional and shear energy densities
c         NORMALIZATIONS:
c         up, vp, and wp are normalized by 1/rn through (rj/rn) 
c                 in the equations
c         v, vp, w, and wp all carry  an implicit sqrt(ll(ll+1)) 
c                 normalization
c
          if (jcom .eq. 3) then
            do j = 1, knot
              rj = dble(radius(j))
              u = dble(buf(j))*scale
              up = dble(buf(j + knot))*scale
              v = dble(buf(j + k2))*scale
              vp = dble(buf(j + k3))*scale
              ff = 2.d0 * u - fl2(ll) * v
              kk = (up*(rj/rn) + ff)**2
              mm = third*(2.0*up*(rj/rn) - ff)**2
     &           + (vp*(rj/rn) - v + fl2(ll)*u)**2
     &           + fl3(ll) * v**2
              intg(j) = (kappa(j)/kapq(j))*kk + (mu(j)/muq(j))*mm
            enddo
          elseif (jcom .eq. 1) then
            do j= 1, knot
              rj = dble(radius(j))
              u = dble(buf(j))*scale
              up = dble(buf(j + knot))*scale
              v = 0.d0
              vp = 0.d0
              ff = 2.d0 * u - fl2(ll) * v
              kk = (up*(rj/rn) + ff)**2
              mm = third*(2.0*up*(rj/rn) - ff)**2
     &           + (vp*(rj/rn) - v + fl2(ll)*u)**2
     &           + fl3(ll) * v**2
              intg(j) = (kappa(j)/kapq(j))*kk + (mu(j)/muq(j))*mm
            enddo
          elseif (jcom .eq. 2) then
            do j = 1, knot
              rj = dble(radius(j))
              wl = dble(buf(j))*scale
              wp = dble(buf(j + knot))*scale
              mm = ((fl3(ll) * wl**2) + ((rj/rn)*wp - wl)**2)
c              if (nn.eq.34 .and. ll.eq.1) then
c                 print*,nn,ll,mm,mu(j),muq(j)
c              end if
              intg(j) = (mu(j)/muq(j))*mm
            enddo
          endif
 35       continue
c
c         form the integrand - two-point trapezoidal integration
c
          t = 0.d0
          do i = 1, k1
c            if (nn.eq.34 .and. ll.eq.1) then
c              print*,nn,ll,intg(i),intg(i+1),dr(i),t
c            end if
            t = t + 0.5d0*(intg(i)+intg(i+1))*dr(i)
          enddo
          qq = 1.d0/t
ccc
ccc     save values from mineos
ccc
        else
          qq = h1
ccc
        endif
ccc
c
c       save the value of the gravity eigenfunction at the 
c       earth's surface
c       or - if this an oceanic model - on the ocean floor
c
        if (jcom .eq. 3) then
          phi = buf(5*knot - obs)
        else
          phi = 0.
        endif
c
c       calculate phase velocity
c
        if (ll .ne. 0) then
          mhz = 35.d0
          wref = dble(mhz) * 2.d0 *pi / 1000
c         wref = 0.3500140
		  print *, w,wref,mhz,1.d0 + (1.d0/(qq*pi))*dlog(w/wref)
          wwd = w*(1.d0 + (1.d0/(qq*pi))*dlog(w/wref))
          cv = w * (rn/1000.) / fl2(ll)
          cvq = wwd * (rn/1000.) / fl2(ll)
		  T = 1. / (w/(2.*pi)) 
          Tq =  1. / (wwd/(2.*pi))
        else
          cv = 99999999.
        endif
c
c       write out q and rest of the mode header 
c   
        write(3,120) nn,ll,w,qq,phi,cv,gv,cvq,Tq,T
      else
        print*, 
     &   ' mode ',nn,comp(jcom),ll,' overlaps with previously read mode'
      endif
      goto 30
c
c     exit this program
c
999   continue
      close(2)
c
c     loop to continue q integration on another mode file
c
      go to 5
c
c     formats and other such garbage
c
110   format(f8.0,3f9.0)
120   format(2i6,f15.5,f20.5,6f15.5)
129   format(i4)
130   format(3f11.0)
      end
c
c
c
      subroutine interpol(n1, n2, x, y, m, b)
c
c     computes the coefficients for linear interpolation
c     y = mx + b
c
c     inputs:
c       n1:      lower bound for interpolation
c       n2:      upper bound for interpolation
c       x(n):    points at which the function is evaluated
c       y(n):    function to be interpolated
c     outputs:
c       m(n):    slopes of lines
c       b(n):    intercepts
c
      save
      parameter (n=1000)
      real x(n), y(n)
      real b(n), m(n)
c
      if ((n2-n1) .gt. n) then
        print*,' array limits exceeded in interpl'
        stop
      endif
      do i = n1, n2-1
        dx = x(i+1) - x(i)
        dy = y(i+1) - y(i)
        if (dx .eq. 0.) then
          m(i) = 999.0
        else
          m(i) = dy/dx
        endif
        b(i) = y(i) - m(i)*x(i)
      end do
      return
      end
c
c
c
      real function interple(n1, n2, x, dx, xlast, y, m, b)
c
c     given the coefficients for linear interpolation
c     this routine calculates y for an input x
c
c     inputs:
c       n1:      lower bound
c       n2:      upper bound
c       x(n):    array of x-values
c       dx:      point a which the function is to be evaluated
c       y(n):    function to be interpolated
c       m(n-1):  slopes
c       b(n-1):  intercepts
c     returned
c       y:       interpolated value
c
      parameter (n=1000)
      real x(n), dx, y(n)
      real b(n), m(n), xlast
c
      if ((n2-n1) .gt. n) then
        print*,' array limits exceeded in interpl'
        stop
      endif
c
      do i = n1, n2
        if (dx .eq. x(i)) then
          if (dx .eq. x(i+1)) then
            if (xlast .eq. 0.) then
              interple = y(i+1)
              return
            elseif (xlast .lt. x(i)) then
              interple = y(i)
              return
            else
              interple = y(i+1)
              return
            endif
          else
            interple = y(i)
            return
          endif
        elseif ((dx .gt. x(i)) .and. (dx .lt. x(i+1))) then
          if (m(i) .ge. 999.0) then
            if (xlast .lt. dx) then
              interple = y(i)
            else
              interple = y(i+1)
            endif
          else
            interple = m(i)*dx + b(i)
          endif
          return
        endif
      end do
20    continue
c
c     outside array bounds - extrapolate
c
      if (dx .lt. x(n1)) then
        interple = m(n1)*dx + b(n1)
      elseif (dx .gt. x(n2)) then
        interple = m(n2)*dx + b(n2)
      else
        print*,' error in interpolation'
      endif
      return
      end

