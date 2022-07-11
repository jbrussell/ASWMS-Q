      program frechet_G
c
c     WARNING currently hardwired to compute anisotropic partials
c     irregardless of input model.  Reset variable "force" to
c     ".false." to turn off this option.
c
c     program to produce complete Frechet kernels in velocity
c
c     will integrate spheroidal, toroidal, and radial modes
c
c     frechet expects a raw mineos output file containing:
c       two record header
c       then 1 record for every mode, depending on type
c        spheroidal modes: nn, ll, w, q, gv, u(1:knots), up(1:knots),
c                    v(1:knots), vp(1:knots), phi(1:knots), phip(1:knots)
c        radial modes:     nn, ll, w, q, gv, u(1:knots), up(1:knots)
c        toroidal modes:   nn, ll, w, q, gv, w(1:knots), wp(1:knots)
c
c     expressions from Backus & Gilbert (1967)(B&G) & Woodhouse & Dahlen (1978) 
c
c     06/14/89 - more or less tested out
c                error in dw of the order 0.015%
c
c     08/02/89 - this version does not rescale the eigenfunctions in the 
c                interest of numerical stability 
c
c     04/30/90 - stores velocity kernels - rather than energy densities  
c                this version combines mineos_frechet and mineos_frechet_sort
c
c     10/08/90 - Modified by Steven S. Shapiro 8 October 1990      
C                to handle perturbations in discontinuity locations. Relevant 
C                equations can be found in, "The effect of a general aspherical 
C                perturbation onthe free oscillations of the earth" by J.H. 
C                Woodhouse and F.A. Dahlen,GJRAS, vol. 53, pp. 335-354, 1978
C                (W&D). The spheroidal and toroidal modes have been tested 
C                using "mineos_frechet_new_sss.ftn". The radial mode section
c                has not been explicitely checked but should be O.K. since 
C                the spheroidal section has been tested and these sections 
c                are essentially identical.
c
c     06/13/91 - now can handle anisotropic models.  Calculates the anisotropic
c                kernels for SV, SH, PV, PH, and eta, in addition to (isotropic) 
c                Q, density.  Primary reference for correct kernels is Anderson 
c                and Dziewonski (1982), alternatively Dziewonski and Anderson 
c                (1981) (PREM), although note that the "k twiddle" term in PREM 
c                is incorrect, it is correct in A&D. Much of the kernel 
c                calculation section of the code was rewritten at this time. 
c                Renamed from mineos_frechet to frechet-- Jim Gaherty
c
c     07/22/92 - minor modification to handle the real*8 variables in the
c                eigenfunction files -- increased size of abuf is the 
c                important change  -- JBG
c
c     08/12/92 - Modified to allow addition of "discontinuities" to model.
c                Prompts user for depth to add, puts it in model header,
c                and calculates a kernel for it.  Allows a "new" disco
c                to develop via inversion.
c
c     10/22/92 - Corrected previous error in the calculation of partials 
c                wrt discontinuity perturbations for anisotropic modes.
c                Note that previous anisotropic kernels calculated from
c                and isotropic starting models are ok.
c 
c     01/24/00 - Updated output file record # 1 to include the record
c                length -- all following programs (partial, calc_dw,
c                idagrn5) must be changed as well to take advantate
c                This avoids the problem that files have no record of
c                the value of knot_t or knot_s that they were written
c                with, which might need to change over time.
c
C23456789112345678921234567893123456789412345678951234567896123456789712
c
c     SEE file FRECHET.DEFS for variable definitions
c
      implicit real*4 (a-h,o-z)
c
      include 'parameter.h'
c                
      integer*4 nqknot
      parameter (nqknot = 1000)
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
c
      data pi/3.14159265350/
      data rn/6371000.0/
      data ric/1221500.0/
      data roc/3480000.0/
      data bigg /6.6732e-11/
      data rhobar/5515.0/
c
c     IF WANT TO FORCE ANISOTROPIC PARTIALS TO BE CALC'ed
cad      data force/.true./
      data force/.false./
c
      do ii = 0, maxn
        do jj = 0, maxl
          old(jj,ii) = 0
        end do
      end do
c
      third = 1.0/3.0
      tthird = 2.0/3.0
      fthird = 4.0/3.0
c
c     scale is the normalizing factor for u, up, v, vp, w, and wp
c     scale2 is the normalizing factor for phi and phip 
c     scale3 is scale2/scale - to avoid numerical instability
c
      scale = 1.0/(rn*sqrt(rn*pi*bigg)*rhobar)
      scale2 = sqrt(pi*bigg/rn)
      scale3 = scale2/scale
c
c     open and read Q file
c
      print*,' enter name of q model'
      read(*,'(a)') q_file
      open(unit=2,file=q_file)
      read(2,*) nq
      if (nq .gt. nqknot) then
        print*,' array limits of q exceeded'
        print*,' reset nqknot'
        stop
      endif
      read(2,*)(rq(i),shearq(i),bulkq(i),i=1,nq)
      close(2)
c
c     convert Q radii to meters
c
      do 10 j=1,nq
        rq(j)=rq(j)*1000.
 10   continue
c
c     get the slopes and intercepts to linearly fit 
c     Q subscript kappa and Q subscript mu
c 
      call interpol(1,nq,rq,bulkq,mbq,bbq)
      call interpol(1,nq,rq,shearq,msq,bsq)
c
c     open branch file
c
      print*,' enter name of branch file'
      read(*,'(256a)') b_file
cad      open(unit=1,file=b_file,form='unformatted',access='sequential',
cad     +    recl=24000)

      open(unit=1,file=b_file,form='unformatted',access='sequential')

      read(1) jcom, nmodes, nnmin, nnmax, llmin, llmax
      read(1) nb, (numb(ii), ii = 1, nb)
      do ii = 1, nb
	read(1) (nnb(jj,ii), llb(jj,ii), 
     &          (stuff(il),il = 1,4),jj = 1,numb(ii))
      end do
      close(1)
c
c     read in output file name -- not opened til later
c
      print*,' enter name of output frechet file'
      read(*,'(256a)') o_file
c
c     open eigenfunction file and read header records 
c     - beginning of loop over eigenfunction files
c
      ifirst = 0
  5   print*,' Enter name of input eigenfunction file'
      read(*,'(a)') m_file 
      if (m_file .eq. ' ') then
        close(2)
        go to 1000
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
      read(2) jcom, wmin, wmax, lmin, lmax, wgrav
      print*, ' mode type = ',jcom,' lmin = ',lmin,' lmax = ', lmax
      print*, ' wmin = ', wmin,' wmax = ', wmax
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
        read(2) ksave,nic,noc,ifanis,tref,(rad(i), i=1,ksave),
     &      (dn(i), i=1,ksave),(alphav(i), i=1,ksave),
     &      (betav(i), i=1,ksave), (muq(i), i=1,ksave),
     &      (kapq(i), i=1,ksave),(alphah(i), i=1,ksave),
     &      (betah(i), i=1,ksave),(eta(i), i=1,ksave)
        print*, ' knots nic, and noc in input= ',ksave,nic,noc
	print*, ' ifanis, force= ',ifanis,force
c
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
	
c ad added this part to check if nrec is large enough
ccc
        write(6,*)'nrec, lmax: ',nrec,lmax
	if(nrec.lt.lmax*4)then	  
	  nrec=lmax*4
	  write(6,"('lmax is large.')")
	  write(6,"('resetting nrec = ',i6)")nrec
	endif

c ad end of addition
ccc

c
c       check to see whether this file exists
c
        inquire(file=o_file,exist=lexist)
c
c       open file output file and write header records
c       pattylin write header record to output file
        if (.not.(lexist)) then
c
          open(unit=3,file=o_file,form='unformatted',access='direct',
     +        recl=nrec,status='new')
          write(3,rec=1)jcom,nmodes,nnmin,nnmax,llmin,llmax,ifanis,nrec
          write(3,rec=2) ksave,nic,noc,ifanis,tref,scale,scale2,ndisc
          write(3,rec=3) (rad(i),i=1,ksave), (kntdsc (i), i = 1, ndisc)
          write(3,rec=4) (dn(i), i=1,ksave)
          write(3,rec=5) (alphav(i), i=1,ksave)
          write(3,rec=6) (betav(i), i=1,ksave)
          if (ifanis.eq.1) then
             write(3,rec=7) (alphah(i), i=1,ksave)
             write(3,rec=8) (betah(i), i=1,ksave)
             write(3,rec=9) (eta(i), i=1,ksave)
             nextrec=10
          else
             nextrec=7
          end if
          print*,'ifanis=',ifanis,' nextrec=',nextrec
          write(3,rec=nextrec) nb, (numb(ii), ii = 1, nb)
          do ii = 1, nb
            irec = nextrec + ii
            write(3,rec=irec) (llb(jj,ii), jj = 1, numb(ii))
          end do
        else 
          if (ifanis.eq.1) then
            nextrec=10
          else
            nextrec=7
          end if
          open(unit=3,file=o_file,form='unformatted',access='direct',
     +        recl=nrec,status='old')
          irec = nextrec + nb
          ic = 0
          do ii = 1, nb
            do im = 1, numb(ii)
              ic = ic + 1
              ni = nnb(im,ii)
              li = llb(im,ii)
              jrec = irec + ic
              read(3,rec=jrec) nn, ll
              if (ni .eq. nn) then
                if (li .eq. ll) then
                  old(ll,nn) = 1
                endif
              endif
            end do
          end do
        endif
      endif
c
c     set up additional constants on a file by file basis
c     expressions using angular order
c
      do j = lmin, lmax
        f = j
        fl1(j) = f * (f + 1.0)
        fl2(j) = sqrt(fl1(j))
        fl3(j) = (f - 1.0) * (f + 2.0)
        fl4(j) = f + 1.0
      end do
c
c     begin reading eigenfunction file
c
c     calculate compressional and shear energy densities and 
c     gravitational potential energy
c
c     NORMALIZATIONS:
c       u, up, v, vp, w and wp should be multiplied by scale;
c       up, vp, and wp are normalized by 1/rn through (rm = (rj/rn)) 
c         in the equations
c       v, vp, w, and wp all carry sqrt(ll(ll+1)) normalization 
c
c       NOTE anisotropic expressions in PREM and A&D do not carry 
c       this normalization thus all the ll(ll+1)) terms in THEIR 
c       expressions are sqrt(ll(ll+1)) terms here.
c
c       phi and phip are normalized by scale2
c       phip is additionally normalized by 1/rn
c
c       unscaling removed 8/1/89 to avoid potential instabilties
c      pylin for spheroidal mode  toridal mode see 847
      if (jcom .eq. 3) then
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
c
c         gravitational energy density from B&G equation 48
c
          rr1 = -(u**2 + v**2)*(w*rj)**2 
          rr2 = 2.0*rj*u*(rm*phip+4.0*pi*dnn(j)*bigg*u*rj-ff*gg(j))
          rr3 = 2.0*rj*fl2(ll)*v*phi
          t = 0.0d0
          do jj = j, knot - 1
            t = t + 0.5d0*(temp(jj) + temp(jj+1))*dr(jj)
          end do
          rr4 = -8.0*pi*bigg*t*rj**2
          rr = rr1 + rr2 + rr3 + rr4
c
c         calculate energy densities in Love's notation--completely
c         general for both the isotropic and anisotropic cases.  
c         Jim's thesis A21 
          xcc = (up*rm)**2.
          xaa = ff**2.
          xff = 2.*rm*up*ff
          xll = (rm*vp + fl2(ll)*u - v)**2.
          xnn = fl3(ll)*v**2. - ff**2.
c
c         kernels used for discontinuity perturbations.  See appendix
c         in Woodhouse and Dziewonski, 1984.
c         Jim's thesis A21 

          xccd = xcc - 2.*(up*rm)**2.
          xlld = xll - 2.*(vp*rm - v + fl2(ll)*u)*vp*rm
          xffd = xff - 2.*ff*up*rm
c
c         now compute the equivilent isotropic energy kernels
c         used even in the anistropic case to compute Q kernel
c         note that mm is the term that is incorrect in PREM
c                                                             
          kk = xaa + xcc + xff
          mm = xll + xnn + tthird*(2.*xaa + 2.*xcc -xff)
c
c         now compute equiv. isotrop. Q kernel
c                                                        
          intg(j,1) = (kappa(j)/kapq(j))*kk + (mu(j)/muq(j))*mm
c
c         if anistropic, %% compute the density and 5 velocity kernels %%
c
          if (ifanis.eq.1) then
c            Jim's thesis A20
c            density kernel
c        
             intg(j,2) = (rr + (xaa + eta(j)*xff)*alphah(j)**2.
     &           + xcc*alphav(j)**2. + (xll-2.*eta(j)*xff)*betav(j)**2.
     &           + xnn*betah(j)**2.)
c
c            SV velocity kernel
c
c             intg(j,3) = 2.*dnn(j)*betav(j)*(xll-2.*eta(j)*xff)

c             intg(j,3) = xll
c
c            PV kernel
c
c             intg(j,4) = 2.*dnn(j)*alphav(j)*xcc
             
c
c            SH kernel
c
c             intg(j,5) = 2.*dnn(j)*betah(j)*xnn
c
c            PH kernel
c
c             intg(j,6) = 2.*dnn(j)*alphah(j)*(xaa + eta(j)*xff)
c
c            eta kernel
c
c             intg(j,7) = dnn(j)*(alphah(j)**2. 
c     &                   - 2.*betav(j)**2.)*xff
c

c            ANISOTROPIC KERNELS (2/16/17 JBR)
c            L kernel (G)
             intg(j,3) = dnn(j)* betav(j)**2 *xll
             
c            C kernel
             intg(j,4) = dnn(j)* alphav(j)**2 *xcc               

c            N kernel
             intg(j,5) = dnn(j)* betah(j)**2 *xnn
             
c            A kernel
             intg(j,6) = dnn(j)* alphah(j)**2 *xaa
             
c            F kernel
             intg(j,7) = dnn(j)* eta(j)*(alphah(j)**2
     &                   - 2.*betav(j)**2) *xff
              

c            discontinuity kernel -- COE-side - surface-side
c
             if (isdisc (j)) then       
               discsv = dnn(j)*(rr + xaa*alphah(j)**2. + 
     &               xccd*alphav(j)**2. + xlld*betav(j)**2. +
     &               xnn*betah(j)**2. + xffd*eta(j)*
     &               (alphah(j)**2. + 2.*betav(j)**2.))
             else if (isdisc (j-1)) then
               dscnum = dscnum + 1                           
               disc (dscnum) = dnn(j)*(rr + xaa*alphah(j)**2. + 
     &               xccd*alphav(j)**2. + xlld*betav(j)**2. +
     &               xnn*betah(j)**2. + xffd*eta(j)*
     &               (alphah(j)**2. + 2.*betav(j)**2.))- discsv
             else
             end if
c
c         otherwise, calculate equivilent isotropic kernels
c
          else
             mvs=mm*betav(j)
c
c            density kernel
c
             intg(j,2) = rr + mvs*betav(j) + kap(j)*kk 
c
c            shear velocity
c
             intg(j,3) = 2.0*dnn(j)*betav(j)*(mm - kk*fthird)
c
c            compressional velocity
c
             intg(j,4) = 2.0*dnn(j)*alphav(j)*kk
c
             if (isdisc (j)) then 
c      
c              on the center_of_the_earth-side of the discontinuity (-)
c
c               kkbar = kk - 2.0 * (up * rm + ff) * up * rm
c               mmbar = mm - tthird * (2.0 * up * rm - ff) * (2 * up * rm)
c     &               - 2.0 * vp * rm * (vp * rm - v + fl2 (ll) * u)
c
               kkbar = xaa + xccd + xffd
               mmbar = xlld + xnn + tthird*(2.*xaa + 2.*xccd -xffd)
               discsv = kappa(j)*kkbar + mu(j)*mmbar + dnn(j)*rr
             else if (isdisc (j-1)) then
c
c              on the surface-side of the discontinuity (+)
c
               kkbar = xaa + xccd + xffd
               mmbar = xlld + xnn + tthird*(2.*xaa + 2.*xccd -xffd)
               dscnum = dscnum + 1                           
               disc (dscnum) = (kappa(j)*kkbar + mu(j)*mmbar + 
     &                       dnn(j)*rr) - discsv
             else
             end if
c
          end if
c        
        end do
c
c       do not overwrite Q if it were calculated by mineos
c        
        if(tref.lt.0.) then
          t = 0.0d0
          do i = 1, knot-1
            t = t + 0.5d0*(intg(i,1)+intg(i+1,1))*dr(i)
          end do
          qq = (1.d0/(t*scale**2))
        else
          qq = q
        end if
c
c       calculate phase velocity
c
        if (ll .ne. 0) then
          cv = w * (rn/1000.) / fl2(ll)
        else
          cv = 99999999.
        endif
c
c       all done for this mode
c       
        if(ifanis.eq.1) then
          nkern=7
          fmt='(f10.1,6f20.15)'
        else 
          nkern=4
          fmt='(f10.1,3f20.15)'
        end if
        write(3,rec=irec) nn, ll, w, qq, gv, cv,
     &                    ((intg(j,ii),j = 1, knot),ii=2,nkern),
     &                    (disc (j), j = 1, ndisc)
         write(*,*) nn, ll, w, qq, gv, cv,knot, nkern, ndisc
c     &                    ((intg(j,ii),j = 1, knot),ii=2,nkern),
c     &                    (disc (j), j = 1, ndisc)    
     
c        if (mod(ll,100).eq.0) then
c          write(6,'(2i4,4f15.5)') nn, ll, w, qq, gv, cv
c          write(6,'(i4,e15.5)') (j,disc(j), j=1,ndisc)
c        end if
c         write(6,fmt)(radius(j),(intg(j,ii)*scale**2,
c     &         ii=2,nkern),j=1,knot)
c
        go to 100
c
c     toroidal modes
c
      elseif (jcom .eq. 2) then
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
c
c         calculate %% energy densities %% in Love's notation--completely
c         general for both the isotropic and anisotropic cases
c
          xll = (rm*wp - wl)**2.
          xnn = fl3(ll) * wl**2.
c
c         kernels used for discontinuity perturbations.  See appendix
c         in Woodhouse and Dziewonski, 1984.
c
          xlld = xll - 2.*wp*rm*(wp*rm-wl)
c
c         now compute the equivilent isotropic energy kernels
c         used even in the anistropic case to compute Q kernel
c                                                             
          mm = xll + xnn
c
c         gravitational energy density
c
          rr = -(wl*w*rj)**2 
c
c         equivilent isotropic Q kernel
c
          intg(j,1) = (mu(j)/muq(j))*mm
c 
c         if model is anisotropic, calculate 2 velocity kernels
c
          if (ifanis.eq.1) then
c
c            density kernel
c        
             intg(j,2) = (rr + xll*betav(noc+j)**2.
     &           + xnn*betah(noc+j)**2.)
c
c            SV velocity kernel
c
c             intg(j,3) = 2.*dnn(j)*betav(noc+j)*xll
c
c            SH kernel
c
c             intg(j,4) = 2.*dnn(j)*betah(noc+j)*xnn
             
c            ANISOTROPIC KERNELS (2/16/17 JBR)
c            L kernel (G)
             intg(j,3) = dnn(j)* betav(noc+j)**2 *xll           

c            N kernel
             intg(j,4) = dnn(j)* betah(noc+j)**2 *xnn
             
                          
c
c            discontinuity kernel -- COE-side - surface-side
c
             if (isdisc (j)) then       
               discsv = dnn(j)*(rr + xlld*betav(noc+j)**2. +
     &               xnn*betah(noc+j)**2.) 
             else if (isdisc (j-1)) then
               dscnum = dscnum + 1                           
               disc (dscnum) = dnn(j)*(rr + xlld*betav(noc+j)**2. +
     &               xnn*betah(j+noc)**2.)- discsv
             else
             end if
c
c         otherwise, calculate equivilent isotropic kernels
c
          else           
             mvs = mm*betav(noc + j)
c
c            density
c
             intg(j,2) = rr + mvs*betav(noc + j)
c
c            shear velocity
c
             intg(j,3) = 2.0*dnn(j)*mvs
c
             if (isdisc (j)) then
c
c              on the center_of_the_earth-side of the discontinuity (-)
c
c               mmbar = mm - 2.0 * wp * rm * (wp * rm - wl)
c
               mmbar = xnn + xlld
               discsv = mu (j) * mmbar + dnn (j) * rr
             else if (isdisc (j-1)) then
c
c              on the surface-side of the discontinuity (+)
c
               mmbar = xnn + xlld
               dscnum = dscnum + 1         
               disc(dscnum) = (mu(j)*mmbar + dnn(j)*rr) - discsv
             else
             end if
c    
          end if 
c
c
        end do
c
c       do not overwrite Q if it were calculated by mineos
c 
        if (tref.lt.0.) then
           t = 0.0d0
           do i = 1, knot-1
             t = t + 0.5d0*(intg(i,1)+intg(i+1,1))*dr(i)
           end do
           qq = (1.d0/(t*scale**2))
        else
           qq=q
        end if
c
c       calculate phase velocity
c
        if (ll .ne. 0) then
          cv = w * (rn/1000.) / fl2(ll)
        else
          cv = 99999999.
        endif
c
c       all done for this mode
c
        if (ifanis.eq.1) then
           nkern=4
           fmt='(f10.1,3f20.15)'
        else 
           nkern=3
           fmt='(f10.1,2f20.15)'
        end if
        write(3,rec=irec) nn, ll, w, qq, gv, cv,
     &                    ((intg(j,ii),j = 1, knot),ii=2,nkern),
     &                    (disc (j), j = 1, ndisc)
c        if (mod(ll,100).eq.0) then
c          write(6,'(2i4,4f15.5)') nn, ll, w, qq, gv, cv
c          write(6,'(f15.5,e15.5)') (j,disc(j), j=1,ndisc)
c        end if
c        write(6,fmt)(radius(j),(intg(j,ii)*scale**2,
c     &         ii=2,nkern),j=1,knot)
c        write(6,'(3f20.15)') ((intg(j,ii)*scale**2,j = 1, knot),ii=2,nkern)
c        write(6,'(3f20.15)') (intg(5,ii)*scale**2,ii=2,nkern)
c
        go to 101
c
c     radial modes--copied down from spheroidal modes section.  Not
c     tested with anisotropic model.  Note that V terms are included,
c     but v and vp are set to zero
c
      elseif (jcom .eq. 1) then
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
c
c         gravitational energy density from B&G equation 48
c
          rr1 = -(u**2 + v**2)*(w*rj)**2 
          rr2 = (2.0*rj*u)*(rm*phip+4.0*pi*dnn(j)*bigg*u*rj-ff*gg(j))
          rr3 = 2.0*rj*fl2(ll)*v*phi
          t = 0.0d0
          do jj = j, knot - 1
            t = t + 0.5d0*(temp(jj) + temp(jj+1))*dr(jj)
          end do
          rr4 = -8.0*pi*bigg*t*rj**2
          rr = rr1 + rr2 + rr3 + rr4
c
c         calculate energy densities in Love's notation--completely
c         general for both the isotropic and anisotropic cases
c         pattylin 201410 note: in Jim's thesis, eqA21  
          xcc = (up*rm)**2.
          xaa = ff**2.
          xff = 2.*rm*up*ff
          xll = (rm*vp + fl2(ll)*u - v)**2.
          xnn = fl3(ll)*v**2. - ff**2.
c
c         kernels used for discontinuity perturbations.  See appendix
c         in Woodhouse and Dziewonski, 1984.
c         pattylin 201410 note: in Jim's thesis, eqA23 
c
          xccd = xcc - 2.*(up*rm)**2.
          xlld = xll - 2.*(vp*rm - v + fl2(ll)*u)*vp*rm
          xffd = xff - 2.*ff*up*rm
c
c         now compute the equivilent isotropic energy kernels
c         used even in the anistropic case to compute Q kernel
c         note that mm is the term that is incorrect in PREM
c                                                             
          kk = xaa + xcc + xff
          mm = xll + xnn + tthird*(2.*xaa + 2.*xcc -xff)
c
c         now compute equiv. isotrop. Q kernel
c                                                        
          intg(j,1) = (kappa(j)/kapq(j))*kk + (mu(j)/muq(j))*mm
c
c         if anistropic, compute the density and 5 velocity kernels
c
          if (ifanis.eq.1) then
c
c            density kernel
c        
             intg(j,2) = (rr + (xaa + eta(j)*xff)*alphah(j)**2.
     &           + xcc*alphav(j)**2. + (xll-2.*eta(j)*xff)*betav(j)**2.
     &           + xnn*betah(j)**2.)
c
c            SV velocity kernel
c
             intg(j,3) = 2.*dnn(j)*betav(j)*(xll-2.*eta(j)*xff)
c
c            PV kernel
c
             intg(j,4) = 2.*dnn(j)*alphav(j)*xcc
c
c            SH kernel
c
             intg(j,5) = 2.*dnn(j)*betah(j)*xnn
c
c            PH kernel
c
             intg(j,6) = 2.*dnn(j)*alphah(j)*(xaa + eta(j)*xff)
c
c            eta kernel
c
             intg(j,7) = dnn(j)*(alphah(j)**2. 
     &                   - 2.*betav(j)**2.)*xff
c
c            discontinuity kernel -- COE-side - surface-side
c
             if (isdisc (j)) then       
               discsv = 2.*dnn(j)*(rr + xaa*alphah(j)**2. + 
     &               xccd*alphav(j)**2. + xlld*betav(j)**2. +
     &               xnn*betah(j)**2. + xffd*eta(j)*
     &               (alphah(j)**2. + 2.*betav(j)**2.))
             else if (isdisc (j-1)) then
               dscnum = dscnum + 1                           
               disc (dscnum) = 2.*dnn(j)*(rr + xaa*alphah(j)**2. + 
     &               xccd*alphav(j)**2. + xlld*betav(j)**2. +
     &               xnn*betah(j)**2. + xffd*eta(j)*
     &               (alphah(j)**2. + 2.*betav(j)**2.))- discsv
             else
             end if
c
c         otherwise, calculate equivilent isotropic kernels
c
          else
             mvs = mm*betav(j)
c
c            density
c
             intg(j,2) = rr + mvs*betav(j) + kap(j)*kk
c
c            shear velocity
c
             intg(j,3) = 2.0*dnn(j)*(mvs - betav(j)*kk*fthird)
c
c            compressional velocity
c
             intg(j,4) = 2.0*dnn(j)*alphav(j)*kk  
c
             if (isdisc (j)) then 
c      
c              on the center_of_the_earth-side of the discontinuity (-)
c
c               kkbar = kk - 2.0 * (up * rm + ff) * up * rm
c               mmbar = mm - tthird * (2.0 * up * rm - ff) * (2 * up * rm)
c     &               - 2.0 * vp * rm * (vp * rm - v + fl2 (ll) * u)
c
               kkbar = xaa + xccd + xffd
               mmbar = xlld + xnn + tthird*(2.*xaa + 2.*xccd -xffd)
               discsv = kappa(j)*kkbar + mu(j)*mmbar + dnn(j)*rr
             else if (isdisc (j-1)) then
c
c              on the surface-side of the discontinuity (+)
c
               kkbar = xaa + xccd + xffd
               mmbar = xlld + xnn + tthird*(2.*xaa + 2.*xccd -xffd)
               dscnum = dscnum + 1                           
               disc (dscnum) = (kappa(j)*kkbar + mu(j)*mmbar + 
     &                       dnn(j)*rr) - discsv
             else
             end if
          end if
c
        end do
c
c       do not overwrite Q if it were calculated by mineos
c       
        if (tref.lt.0.) then
           t = 0.0d0
           do i = 1, knot-1
             t = t + 0.5d0*(intg(i,1)+intg(i+1,1))*dr(i)
           end do
           qq = (1.d0/(t*scale**2))
        else
           qq=q
        end if
c
c       calculate phase velocity
c
        if (ll .ne. 0) then
          cv = w * (rn/1000.) / fl2(ll)
        else
          cv = 99999999.
        endif
c
c       all done for this mode
c
c        write(6,'(2i4,4f15.5)') nn, ll, w, qq, gv, cv
c        write(6,'(3f20.15)') ((intg(j,ii)*scale**2,j = 1, knot),ii=2,4)
        if (ifanis.eq.1) then
           nkern=7
        else
           nkern=4
        end if
		
c		!!!! pattylin OUPUT THINGS is here !!!!
c		!!!! intg(j,ii) jim's A20  length(ii)= 4 if isotropoc, = 7 if anisortropic 
c		!!!  disc is calculated particularly for discontinues. 
        write(3,rec=irec) nn, ll, w, qq, gv, cv,
     &                    ((intg(j,ii),j = 1, knot),ii=2,nkern),
     &                    (disc (j), j = 1, ndisc)
        go to 102
c
      endif
      print*,' i am lost'
c
1000  continue
      ic = 0 
      irec = nextrec + nb
      do ii = 1, nb
        do im = 1, numb(ii)
          ic = ic + 1
          ni = nnb(im,ii)
          li = llb(im,ii)
          jrec = irec + ic
          read(3,rec=jrec) nn, ll
          if ((ni .ne. nn) .or. (li .ne. ll)) then
            print*,' mismatch in branch structure: ', ni, li, nn, ll
          endif
        end do
        print*,' branch completed: ', ii
      end do
c
      close(3)
      stop
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

