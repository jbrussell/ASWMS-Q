      program eigen_asc
      
c Written/Modified by Josh Russell: 9/13/15
c
c       Input: (1) name of input branch file (input)
c              (2) name of output ascii file (output)
c			   (3) name of raw eigenfunction file (input)
c              (4) number of discontinuities
c
c		Program to convert Mineos unformatted eigenfunction file to 
c		ascii file to be plotted
c
c       This code was adapted from frechet.f which reads in the unformatted
c       eigenfunction file and calculates frechet derivatives
c

c      Raw mineos output eigenfunction file contains:
c       two record header
c       then 1 record for every mode, depending on type
c        spheroidal modes: nn, ll, w, q, gv, u(1:knots), up(1:knots),
c                    v(1:knots), vp(1:knots), phi(1:knots), phip(1:knots)
c        radial modes:     nn, ll, w, q, gv, u(1:knots), up(1:knots)
c        toroidal modes:   nn, ll, w, q, gv, w(1:knots), wp(1:knots)
c

      implicit real*4 (a-h,o-z)
c
      include 'parameter.h'
c                
      integer*4 nqknot
      parameter (nqknot = 20)
c                 
      real*8 temp (nknot), t,dr(nknot), wd,qd,gvd
c
      real*4 u, up, v,vp, phi, phip, wl,wp, w, wmin, wmax
c
      real*4 ff, mm, mmbar, kk, kkbar, rr, rr1, rr2, rr3,
     &      rr4, qq, mvs, xaa, xcc, xff, xll, xnn
c
      real*4 abuf(maxbyte+3),stuff(4), intg(nknot,7),
     &       disc (maxdisc), discsv
c                                                  
      real*4 rad(nknot), radius(nknot), rn, rj, rm, rl
c
      real*4 bigg, rhobar, pi, gg (nknot)
c
      real*4 bulkq(nqknot), shearq(nqknot), rq(nqknot),
     &       mbq(nknot), bbq(nknot), msq(nknot), bsq(nknot)
c
      real*4 kappa(nknot),kapq(nknot),kap(nknot),mu(nknot),muq(nknot)
c
      real*4 dn(nknot), alphav(nknot),alphah(nknot),
     &       betav(nknot),betah(nknot), eta(nknot), dnn(nknot)
c                                       
      real*4 wgrav, tref, q, gv, cv   
c
      real*4 scale, scale2, scale3, third, tthird, fthird
c
      real*4 f,fl1(0:maxll),fl2(0:maxll),fl3(0:maxll),fl4(0:maxll)
c
      real*4 interple
c             
      integer*4 nn,ni,ll,li,im,nnmin,nnmax,llmin,llmax,lmin,lmax,
     &     nnb(0:maxl,nbranch),llb(0:maxl,nbranch), numb(nbranch)
c
      integer*4 old(0:maxl,0:maxn),kntdsc(maxdisc),knewd(maxdisc)
c                                           
      integer*4 ifirst, knot, nq, nb, ii, jj, il, jcom, nmodes,
     &    ksave, nic, noc, nocor, ifanis
c                                                                           
      integer*4 ndisc, dscnum

      integer*4 nvec, newvec,nrec,irec,jrec,k1,k2,k3,k4,k5
c
      integer*4 index, ic, i, j

      logical lexist, isdisc(0:nknot), force
c
      character*40 fmt
      character*256 q_file, m_file, o_file, b_file
c
      common/ablk/nn,ll,wd,qd,gvd,buf(nknot6)
      equivalence(nn,abuf)
c equivalence statement means abuf(i) contains entire 'ablk'
c This includes nn,ll,wd,qd,gvd,buf(nknot6)
c

c
c     open branch file
c
      print*,' enter name of branch file'
      read(*,'(256a)') b_file
cad      open(unit=1,file=b_file,form='unformatted',access='sequential',
cad     +    recl=24000)

      open(unit=1,file=b_file,form='unformatted',access='sequential')
cc

      read(1) jcom, nmodes, nnmin, nnmax, llmin, llmax
      print*,'llmax = ',llmax
      read(1) nb, (numb(ii), ii = 1, nb)
      do ii = 1, nb
	read(1) (nnb(jj,ii), llb(jj,ii), 
     &          (stuff(il),il = 1,4),jj = 1,numb(ii))
      end do
      close(1)

c	  open output ascii file
      print*,' Enter name of output ascii file'
      read(*,'(256a)') o_file
      open(unit=3,file=o_file,form='formatted',access='sequential')

c
c     open eigenfunction file and read header records 
c     - beginning of loop over eigenfunction files
c
      ifirst = 0
  5   print*,' Enter name of input eigenfunction file'
      read(*,'(a)') m_file 
      if (m_file .eq. ' ') then
        close(2)
c        go to 1000
      elseif (ifirst .ne. 0) then
        close(2)
      endif
      open(unit=2,file=m_file,form='unformatted',access='sequential',
     +    recl=36000)
      ifirst = ifirst + 1
c
c     note that mineos writes out 9 model parameters even for isotropic
c     models, so we can read in without checking for anisotropy first.
c

c read in first line of header
      read(2) jcom, wmin, wmax, lmin, lmax, wgrav
      print*, 'wmax = ',wmax
      print*, 'lmax = ',lmax
      
c      print*, ' mode type = ',jcom,' lmin = ',lmin,' lmax = ', lmax
c      print*, ' wmin = ', wmin,' wmax = ', wmax
c      write(1,10)
c  10  format(I,3X,F,F,F,F
c
c     if not the first file, skip the model information.  Otherwise,
c     do some work the first time through - assuming all files are 
c     for the same model
c
      if (ifirst .ne. 1) then
        read(2)
        print*, ' knots in output= ',ksave
      else
c
c read in second (final) line of header
        read(2) ksave,nic,noc,ifanis,tref,(rad(i), i=1,ksave),
     &      (dn(i), i=1,ksave),(alphav(i), i=1,ksave),
     &      (betav(i), i=1,ksave), (muq(i), i=1,ksave),
     &      (kapq(i), i=1,ksave),(alphah(i), i=1,ksave),
     &      (betah(i), i=1,ksave),(eta(i), i=1,ksave)
c        print*, ' knots nic, and noc in input= ',ksave,nic,noc
c	print*, ' ifanis, force= ',ifanis,force
c	do i=1,ksave
c	 print*, 'rad(i) = ',rad(i),'dn(i) = ',dn(i)
c	 print*, 'alphah(i) = ',alphah(i),'betah(i) = ',betah(i)
c	 print*, 'alphav(i) = ',alphav(i),'betav(i) = ',betav(i)
c	 print*, 'muq(i) = ',muq(i),'eta(i) = ',eta(i)
c	end do
	
	
c       store original knot structure for input
c
        knotos=ksave
        knotot=ksave-noc
c
        if((ifanis.ne.1) .and. (force)) then
          print*,'CALCULATE ANISOTROPIC PARTIALS WITH ISOTROPIC MODEL'
          ifanis=1
        end if
c
c       get discontinuities to add to model.  Directly adds them to
c       input model immediately, and later duplicates eigenfunction
c       for new knot when it's read in from the buffer.  Allows several
c       new discos, dependant only on # in starting model and value of
c       maxdisc.  For my models, maxdisc=20, number of discs ~ 10.
c
        print*,'Enter # of new discontinuities to add to model'
        read(5,*) nnew
        do jj=1,nnew
          print*,'  '
42        print*,'Enter DEPTH (in KM) of new discontinuity'
          print*,'If adding more than one, add deepest first'
          print*,'maxdisc = ',maxdisc,' so BEWARE'
          read(5,*) dnew
          if(dnew.lt.1. .or. dnew.ge.rn) then
             print*,'bad depth: try again'
             go to 42
          end if
          rnew=rn - dnew*1000.
c
          do ii=ksave,1,-1
             rad(ii+1)=rad(ii)
             dn(ii+1)=dn(ii)
             alphav(ii+1)=alphav(ii)
             betav(ii+1)=betav(ii)
             muq(ii+1)=muq(ii)
             kapq(ii+1)=kapq(ii)
             alphah(ii+1)=alphah(ii)
             betah(ii+1)=betah(ii)
             eta(ii+1)=eta(ii)
             if (rnew.ge.rad(ii)) then
               print*,'Added discontinuity: ',(rn-rad(ii))/1000.,' km'
               ksave=ksave+1
               if(rnew.lt.ric) then
                  nic=nic+1
                  noc=noc+1
               elseif(rnew.lt.roc) then
                  noc=noc+1
               endif
               knewd(jj)=ii
               if (jj.gt.1) then
                 if (knewd(jj).lt.knewd(jj-1)) then
                    print*,'ERROR: order new discos deep to shallow'
                    print*,'Stopped.  Try again.'
                    stop 
                 end if
               end if
               if(rnew.eq.rad(ii-1)) then
                 print*,'WARNING -- DISCONTINUITY ALREADY EXISTED'
               end if
               go to 45
             end if
          end do
45        continue
        end do
        knewd(nnew+1)=-1
c
c       debugging
c 
c        do i=1,ksave
c          write(6,'(F8.0,3F9.2,2F9.1,2F9.2,F9.5)')
c     &     rad(i),dn(i),alphav(i),betav(i),muq(i),kapq(i),alphah(i),
c     &     betah(i),eta(i)
c        end do
c
        nocor = ksave - noc
        rad(ksave+1)=0.0
        rad(1)=1.0
        print*,'  '
        print*,'knots, nic, and noc in output',ksave,nic,noc
c
c       count number of discontinuities
c
        ndisc = 0
        isdisc (0) = .false.
c
c
c       set up integration parameters and constants for spheroidal 
c       and radial modes
c                                                           
        print*, 'The discontinuity radii (km) are:'
c       pattylin for spheroidal and radial modes        
		if (jcom .ne. 2) then
          knot = ksave
          rl = rad(1)
          do j = 1, knot
             radius(j) = rad(j)
c
c            equivilent isotropic mu, kappa from PREM, if model is anisotropic
c
             if (ifanis.eq.1) then 
                mu(j) = 1./15.*dn(j) * ((1.-2.*eta(j))*alphah(j)**2. + 
     &                 alphav(j)**2. + 5.*betah(j)**2 +
     &                 (6.+4.*eta(j))*betav(j)**2.)
                kappa(j) = 1./9.*dn(j) * ((4.+4.*eta(j))*alphah(j)**2.+
     &                    alphav(j)**2. - 8.*eta(j)*betav(j)**2. -
     &                    4.*betah(j)**2.)
             else
                mu(j) = betav(j)**2*dn(j)
                kappa(j) = alphav(j)**2*dn(j) - fthird*mu(j)
             end if
             kap(j) = kappa(j)/dn(j)
c
c            if not done in mineos, Q subscript mu interpolated at the 
c            eigenfunction file's knot j
c
             if(tref.lt.0) 
     &          muq(j) = interple(1,nq,rq,rad(j),rl,shearq,msq,bsq)
c
             kap(j) = alphav(j)**2 - fthird*betav(j)**2
c
c            if not done in mineos, Q subscript kappa interpolated at 
c            the eigenfunction file's knot j
c
             if(tref.lt.0)
     &          kapq(j) = interple(1,nq,rq,rad(j),rl,bulkq,mbq,bbq)
c
             dnn(j) = dn(j)
             dr(j) = rad(j+1) - rad(j)
c
c            dr ~ 0 => a discontinuity             
c
             if (abs (dr (j)) .le. rfrctn) then 
                ndisc = ndisc + 1
                isdisc (j) = .true.
                kntdsc (ndisc) = j
                print*, ndisc, rad (j) / 1000.0, rad (j+1) / 1000.0
             else
                isdisc (j) = .false.
             end if
c
             temp(j) = dnn(j)*radius(j)**2
             rl = rad(j)                  
          end do
          temp(knot+1) = 0.d0
c
c         integrate density structure for gravitational acceleration 
c         as a function of depth
c    
          do ii = 1, knot
            t = 0.0d0
            do j = 1, ii - 1
              t = t + 0.5d0*(temp(j) + temp(j+1)) * dr(j)
            end do
            gg(ii) = 4.0*pi*bigg*(t/(radius(ii)**2))
          end do
c
c         set up integration parameters and constants for toroidal modes
c
        else
          knot = nocor                    
          rl = rad(noc + 1)
          do j = 1, knot
             radius(j) = rad(noc + j)
c
c            if necessary, mu for anisotropic structure from PREM
c
             if (ifanis.eq.1) then 
                mu(j) = 1./15.*dn(noc+j) * ((1.-2.*eta(noc+j))*
     &                 alphah(noc+j)**2. + 
     &                 alphav(noc+j)**2. + 5.*betah(noc+j)**2 +
     &                 (6.+4.*eta(noc+j))*betav(noc+j)**2.)
             else
                mu(j) = dn(noc +j) * betav(noc +j)**2
             end if
c
c            if not done in mineos, Q subscript mu interpolated at 
c            the eigenfunction file's knot j
c
             if(tref.lt.0) 
     &          muq(j) = interple(1,nq,rq,rad(noc+j),rl,shearq,msq,bsq)
             dnn(j) = dn(noc + j)
             dr(j) = rad(noc + j + 1) - rad(noc + j)
c
c            dr ~ 0 => a discontinuity             
c
             if (abs (dr (j)) .le. rfrctn) then
                ndisc = ndisc + 1
                isdisc (j) = .true.
                kntdsc (ndisc) = noc + j
                print*, ndisc, rad (noc + j) / 1000.0, 
     &                  rad (noc + j + 1) / 1000.0
             else
                isdisc (j) = .false.
             end if
c
             rl = rad(noc + j)
          end do
        endif
c
c       determine the number of bytes for reading and writing buffers.
c       note that number of knots in and out could differ due to 
c       added discontinuities.
c       pattylin for radial and toroidal modes 
        if (jcom .ne. 3) then
          if (ifanis.eq.1) then
             index = 3
          else
             index = 2
          end if
          knoto=knotot
          newvec = index*nknot_t + 6 + maxdisc
          nvec=2*knoto + 8
        else
          if (ifanis.eq.1) then
             index = 6
          else
             index = 3
          end if
          knoto=knotos
          newvec = index*nknot_s + 6 + maxdisc
          nvec=6*knoto + 8
        endif
        k5 = 5*knoto
        k4 = 4*knoto
        k3 = 3*knoto
        k2 = 2*knoto
        k1 = knoto - 1
c
        nrec = newvec * 4
c


c================================================
c
c     begin reading eigenfunction file
c
c

c		Spheroidal Mode
      if (jcom .eq. 3) then
        print*,'Spheroidal Mode'
100     continue
        read(2,end=5) (abuf(i), i = 1, nvec)
c
c       transfer double precision values to single precision
c
        w = sngl(wd)
        q = sngl(qd)
        gv = sngl(gvd)
c
c       look for spurious mode or previously read mode
c
        if (nn .lt. 0) then
          print*, ' apparent error in mode calculation: ',nn, ll
          if (nn.lt.-1) then
            print*,'n=-10: skipping next mode as well'
            read(2,end=5)
          end if
          go to 100
        elseif (old(ll,nn) .gt. 0) then
          print*,'old(ll,nn) = ', old(ll,nn)
          print*,' skipping previously stored mode: ', nn, ll
          go to 100
        endif
        old(ll,nn) = 1
c
c       keep track of records
c
        ic = 0
        do ii = 1, nb
          do im = 1, numb(ii)
            ic = ic + 1
            ni = nnb(im,ii)
            li = llb(im,ii)
            if (ni .eq. nn) then
              if (li .eq. ll) then
                go to 200
              endif
            endif
          end do
        end do
        print*,' mode not found in branch file: ', nn, ll
        go to 100
200     continue
        irec = nextrec + nb + ic
c
c       need to reuse buf for knots with new discontinuities
c
        iold=1
        kkk=1
        do ii = 1, knot
          u = buf(iold)
          v = buf(iold + k2)
          ff = 2.0 * u - fl2(ll) * v
          temp(ii) = u*ff*(dnn(ii)/radius(ii))
c 
c         debugging
c
c          if(rad(ii).eq.rnew) then
c             print*,'radius= ',rnew,' should be there'
c             print*,'knewd= ',knewd(kkk),' knot= ',ii
c          end if              
c
          if(knewd(kkk).eq.ii) then
            kkk=kkk+1
          else 
            iold=iold+1
          end if
        end do
c                                   
        iold=1
        kkk=1
        dscnum = 0
        do j = 1, knot
          u = buf(iold)
          up = buf(iold + knoto)
          v = buf(iold + k2)
          vp = buf(iold + k3)
          phi = buf(iold + k4)*scale3
          phip = buf(iold + k5)*scale3
          if(knewd(kkk).eq.j) then
            kkk=kkk+1
          else 
            iold=iold+1
          end if
c
          rj = radius(j)
          rm = rj/rn
          ff = 2.0 * u - fl2(ll) * v
          
          r = rj/1000
          
          
          write(3,11) ll,nn,w,radius(j),kappa(j),mu(j),dnn(j),u,up,v,vp,phi,phip
11        format(I5,I5,F15.8,F15.3,F20.3,F20.3,F20.3,F20.8,F20.8,F20.8,F20.8,F25.3,F25.3)
        end do
        
        go to 100
        
        
c=================================================        
        
c       Toroidal Mode
      elseif (jcom .eq. 2) then
        print*,'Toroidal Mode'
101     continue
        read(2,end=5) (abuf(i), i = 1, nvec)
c
c       transfer double precision values to single precision
c
        w = sngl(wd)
        q = sngl(qd)
        gv = sngl(gvd)
c
c       look for spurious mode or previously read mode
c
        if (nn .lt. 0) then
          print*, ' apparent error in mode calculation: ',nn, ll
          if (nn.lt.-1) then
            print*,'n=-10: skipping next mode as well'
            read(2,end=5)
          end if
          go to 101
        elseif (old(ll,nn) .gt. 0) then
          print*,' skipping previously stored mode: ', nn, ll
          go to 101
        endif
        old(ll,nn) = 1
c
c       keep track of records
c
        ic = 0
        do ii = 1, nb
          do im = 1, numb(ii)
            ic = ic + 1
            ni = nnb(im,ii)
            li = llb(im,ii)
            if (ni .eq. nn) then
              if (li .eq. ll) then
                go to 201
              endif
            endif
          end do
        end do
        print*,' mode not found in branch file: ', nn, ll
        go to 101
201     continue
        irec = nextrec + nb + ic
c
c
        iold=1
        kkk=1
        dscnum = 0
        do j = 1, knot
          wl = buf(iold)
          wp = buf(iold + knoto)
          rj = radius(j)                        
          rm = rj / rn 
          
          r = rj/1000
          
          
          write(3,10) ll,nn,w,radius(j),kappa(j),mu(j),dnn(j),wl,wp
10        format(I5,I5,F15.8,F15.3,F20.3,F20.3,F20.3,F20.8,F20.8)
c          print*,'nn= ',nn,'ll= ',ll,'wl= ',wl,'wp= ',wp,
c     &     'rj= ',rj,'rm= ',rm
c 
c         debugging
c
c          if(rad(j+noc).eq.rnew) then
c             print*,'radius= ',rnew,' should be there'
c             print*,'knewd= ',knewd(kkk),' knot= ',j+noc
c          end if 
c              
          if(knewd(kkk).eq.(j+noc)) then
            kkk=kkk+1
          else 
            iold=iold+1
          end if
        end do

        go to 101


c===========================================

c		Radial Mode       
      elseif (jcom .eq. 1) then
        print*,'Radial Mode'
102     continue
        read(2,end=5) (abuf(i), i = 1, nvec)
c
c       transfer double precision values to single precision
c
        w = sngl(wd)
        q = sngl(qd)
        gv = sngl(gvd)
c
c       look for spurious mode or previously read mode
c
        if (nn .lt. 0) then
          print*, ' apparent error in mode calculation: ',nn, ll
          if (nn.lt.-1) then
            print*,'n=-10: skipping next mode as well'
            read(2,end=5)
          end if
          go to 102
        elseif (old(ll,nn) .gt. 0) then
          print*,' skipping previously stored mode: ', nn, ll
          go to 102
        endif
        old(ll,nn) = 1
c
c       keep track of records
c
        ic = 0
        do ii = 1, nb
          do im = 1, numb(ii)
            ic = ic + 1
            ni = nnb(im,ii)
            li = llb(im,ii)
            if (ni .eq. nn) then
              if (li .eq. ll) then
                go to 202
              endif
            endif
          end do
        end do
        print*,' mode not found in branch file: ', nn, ll
        go to 102
202     continue
        irec = nextrec + nb + ic
c
c
        iold=1
        kkk=1
        do ii = 1, knot
          u = buf(iold)
          v = buf(iold + k2)
          ff = 2.0 * u - fl2(ll) * v
          temp(ii) = u*ff*(dnn(ii)/radius(ii))
          if(knewd(kkk).eq.ii) then
            kkk=kkk+1
          else 
            iold=iold+1
          end if
        end do
c
        discnum=0
        iold=1
        kkk=1
        do j= 1, knot
          u = buf(iold)
          up = buf(iold + knoto)
          v = 0.0
          vp = 0.0
          phi = 0.0
          phip = 0.0
          if(knewd(kkk).eq.j) then
            kkk=kkk+1
          else 
            iold=iold+1
          end if
c
          rj = radius(j)
          rm = rj/rn
          ff = 2.0 * u - fl2(ll) * v
        end do
        
        go to 102
      
      endif  
      endif    
      end

c************************** FUNCTION ******************************
      real function interple(n1, n2, x, dx, xlast, y, m, b)
c
c     given the coefficients for linear interpolation
c     this routine calculates y for an input x
c
c     inputs:
c       n1:      lower bound
c       n2:      upper bound
c       x(n):    array of x-values
c       dx:      point at which the function is to be evaluated
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