      program draw_frechet_gv_perc
      
cad April 8, 2008
cad modified this program to write out gr-velo kernels
cad      
c
c     program to test frechet kernels
c
c     will integrate spheroidal, toroidal, and radial modes
c
c     expects a mineos_frechet outputfile

c       updated 9/91 to handle anisotropic frechet files
c
c       spheroidal and radial modes
c          nn, ll, w, q, gv, cv, kr(1:knot), kb(1:knot), ka(1:knot), kh (1 : ndisc)
c            plus kbh, kah, ke in anisotropic case
c       toroidal modes
c          nn, ll, w, q, gv, cv, kr(1:knot), kb(1:knot), kh (1 : ndisc)
c            plus kbh in anisotropic case 
c
c       where
c          kr is the density kernel
c          kb is the shear velocity kernel (SV in anisotropic case)
c          ka is the compressional velocity kernel (PV   "      " )    
c          kh is the discontinuity kernel
c          kah is the PH velocity kernel
c          kbh is the SH velocity kernel
c          ke is the eta kernel
c
c       changed to read direct access frechet files 9/27/89
c
c       rewritten (again) for the new format frechet files
c                                                                          c
C23456789112345678921234567893123456789412345678951234567896123456789712
C       Modified by
C       Steven S. Shapiro 8 October 1990      
C       to handle perturbations in discontinuity locations.  Also
C       modified to read eigenfunction files instead of mode tables.
C23456789112345678921234567893123456789412345678951234567896123456789712
C       Modified by
C       Josh Russell 8/18/17
c     
C       Writes out kernels for percent perturbation instead of m/s to check
c       against Colleen Dalton's version of kernel codes. Must run this version
c       with frech_cv_perc.
c
c       MULTIPLIES KERNELS BY V
c
c       OUR SCALING :    dc = integral( (c^2/U) * K * dV ) dr
c       THEIR SCALING: dc/c = integral( (c/U) * V * K * dV/V ) dr
c
C       
C23456789112345678921234567893123456789412345678951234567896123456789712
c
      implicit real*4 (a-h,o-z)
c
      include 'parameter.h'
c
      real*4 rad(nknot), buf(nknot6)
c
c     original model parameters
c
      real*4 dn(nknot), alpha(nknot), beta(nknot)
      real*4 alphah(nknot),betah(nknot),eta(nknot)
c
c     perturbed model parameters
c
      real*4 r1(nknot), d1(nknot), a1(nknot), b1(nknot)
      real*4 ah1(nknot), bh1(nknot), et1(nknot)
c
c     slopes and intercepts of interpolated models
c
      real*4 md1(nknot), bd1(nknot), ma1(nknot), ba1(nknot)
      real*4 mah1(nknot),bah1(nknot),mbh1(nknot),bbh1(nknot)
      real*4 met1(nknot), bet1(nknot)
      real*4 mb1(nknot), bb1(nknot)
c
      real*4 interple, scale, scale2, temp(4), kest, dw
c
      real*4                                          
     &      pdr (nknot),
c  difference in radii between successive knots in the perturbed model
     &      dh (maxdisc)
c  differences in discontinuity radii
      real*8 w2
c  eigenfrequency from eigenfunction file
c
      real*8 dbeta(nknot), dalpha(nknot), drho(nknot)
      real*8 dbetah(nknot),dalphah(nknot),deta(nknot)
      real*8 dmu(nknot), dkappa(nknot)
      real*8 dr(nknot), intg(nknot), t, wdb, fthird
c
cad added this 
      real*8 vker(nknot,6)
      parameter (maxlayers=10)
      real*8 tint(maxlayers,7)
      real*8 ave_alphah(maxlayers),ave_alpha(maxlayers)
      real*8 ave_betah(maxlayers),ave_beta(maxlayers)
      real*4 depthtop(maxlayers),depthbottom(maxlayers)
      integer indexradtop(maxlayers),indexradbottom(maxlayers)
      integer nnsave(10000),llsave(10000),irecsave(10000)
      real*4 wsave(10000)
cad end of addition
c
      integer*4 nn, ll, knot, ifirst
      integer*4 nm(nbranch)
      integer*4 nnb(maxl,nbranch), llb(maxl,nbranch)
c
      integer*4                 
     &         maxskp,
c  maximum number of discontinuities allowed to skip
     &         skp,
c  number of discontinuities skipped
     &         k,
c  knot index
     &         ndisc,
c  number of discontinuities in the unperturbed model
     &         npdisc,
c  number of discontinuities in the perturbed model
     &         kntdsc (maxdisc),
c  contains the value of the knot corresponding to the 
c  center_of_the_earth-side of every discontinuity in the unperturbed model
     &         kntpdsc (maxdisc),
c  contains the value of the knot corresponding to the 
c  center_of_the_earth-side of every discontinuity in the perturbed model
     &         n2,
c  n, radial order from eigenfunction file
     &         l2
c  l, angular order from eigenfunction file
c
      logical ldirect
c
      character*1 comp(3), ans, ans2
      character*256 m_file,o_file,b_file,t_file
c
      data comp /'S','T','S'/
      data pi/3.14159265350/
      data rn/6371000./
      data bigg /6.6723e-11/
      data rhobar/5515.0/
c
      fthird = 4.0d0/3.0d0
c
c     open frechet file and read header records 
c     - beginning of loop over frechet files
c
  5   print*,' Enter name of input frechet file'
      read(*,'(a)')m_file 
      if (m_file .eq. ' ') then
        close(2)
        stop
      endif
c
c
      print*,' Enter name of output file'
      read(*,'(a)')b_file 
      open(11,file=b_file) 
c
      index=2
      newvec = index*nknot_t + 6 + maxdisc
      nrec = newvec * 4 
      print*,'first nrec',nrec

8     open(unit=2,file=m_file,form='unformatted',access='direct',
     +      recl=nrec)
      read(2,rec=1) jcom,nmodes,nnmin,nnmax,llmin,llmax,ifanis,nrecw
      print*,'2nd nrec',nrecw
      if (nrec.ne.nrecw) then
         close(2)
         nrec = nrecw
         go to 8
      endif
      read(2,rec=2) ksave,nic,noc,ifanis,tref,scale,scale2,ndisc
      read(2,rec=3) (rad(i),i=1,ksave), (kntdsc (i), i = 1, ndisc)
      read(2,rec=4) (dn(i), i=1,ksave)
      read(2,rec=5) (alpha(i), i=1,ksave)
      read(2,rec=6) (beta(i), i=1,ksave)
      if(ifanis.ne.0) then
        read(2,rec=7) (alphah(i), i=1,ksave)
        read(2,rec=8) (betah(i), i=1,ksave)
        read(2,rec=9) (eta(i), i=1,ksave)
        nextrec = 10
      else
        do i=1,ksave
	  alphah(i)=0.
	  betah(i)=0.
	  eta(i)=0.
	enddo  
	nextrec = 7
      end if 
      read(2,rec=nextrec) nb, (nm(ii), ii = 1, nb)
      print*,' reading branch structure: ', nb
      do ii = 1, nb
        irec = nextrec + ii
        read(2, rec=irec) (llb(jj,ii), jj = 1,nm(ii))
        do jj = 1, nm(ii)
          nnb(jj,ii) = ii - 1
        end do
      end do
      irec = irec+1

c
      if(ifanis.eq.1) then
        if (jcom .eq. 2) then
          index = 3
        else 
          index = 6
        endif
      else
        if (jcom .eq. 2) then
          index = 2
        else 
          index = 3
        endif
      endif
     
      knot = ksave
      nocor = knot - noc
      print*, ' mode type = ',jcom,' lmin = ', llmin,' lmax = ', llmax
      print*, ' scalings: ', scale, scale2
c
c     set up integration parameters and constants
c
      rad(knot+1)=0.0
      rad(1)=1.0
      if (jcom .eq. 2) then
        knot = nocor
        ind = noc
        kind = index*knot + ndisc
      else
        ind = 0
        kind = index*knot + ndisc
      endif
      print*, ' knots = ',knot
      do j = 1, knot
        dr(j) = rad(ind + j + 1) - rad(ind + j)
      end do
      print*, 'Discontinuities in the Frechet file (rad, depth):'
      do i = 1, ndisc
         print*, i, rad (kntdsc (i)) / 1000.0, 
     &     6371.-0.001*rad (kntdsc (i))
      end do
      
      k6 = 6*knot
      k5 = 5*knot
      k4=  4*knot
      k3 = 3*knot
      k2 = 2*knot
      k1 = knot
      k0 = knot - 1
 
ccc find all frequencies
ccc find the mode closest to desired frequency
ccc
c      open(15,file='testfrechet')
c      write(6,*)'got to here'
c      ir=irec
c      do i=1,nb
c        do j=1,nm(i)
c          read(2,rec=ir,err=998) 
c     &              nn,ll,w,qq,gv,cv,(buf(kk), kk = 1, kind)	
c          write(15,*)i,j,nn,ll,w
c	  ir=ir+1
c	enddo
c      enddo  
c      close(15)
c      
     
      
      irecold=irec
      ib=1
      if(nm(ib).gt.10000) stop 'nm(ib) too large!'
      do i=1,nm(ib)
c       print*,'nm(ib)',nm(ib)
c        print*,'irec',irec
        read(2,rec=irec,err=998) 
     &              nn,ll,w,qq,gv,cv,(buf(kk), kk = 1, kind)
        wsave(i)=w
	nnsave(i)=nn
	llsave(i)=ll
	irecsave(i)=irec
	irec=irec+1
	write(6,*)w,nn,ll
      enddo
      irec=irecold
      
      write(6,"('Period of interest (s)?')")
      read(5,*)period
      vomega=2.*pi/period
      imodefnd=0
      diff=99999.
      do i=1,nm(ib)
        dd=abs(vomega-wsave(i))
	if(dd.lt.diff)then
	  diff=dd
	  imodefnd=i
	endif
      enddo
      p=2.*pi/wsave(imodefnd)
      if(imodefnd.eq.0)then
        stop 'did not find a close mode'
      else
        write(6,"('closest mode found: nn,ll,p ',2i6,f10.3)")
     &    nnsave(imodefnd),llsave(imodefnd),p
      endif
 
c     loop over branch -- recall that gv frechet files are only for 1 branch
cad changed end=998 to err=998 in read statement below
cad seems to solve compilation problem

      ib=1
      do i=1,nm(ib)
        read(2,rec=irec,err=998) 
     &              nn,ll,w,qq,gv,cv,(buf(kk), kk = 1, kind)
c
c     integrate the kernels for the perturbation to eigenfrequency
c
        if(i.ne.imodefnd) goto 123
	
	do j = 1, knot
          intg(j) = 0.0d0
        end do
        if (jcom .ne. 2) then
          if(ifanis.eq.0) then
            do j = 1, knot                         
              intg(j) = buf(j)*drho(j) + buf(k1+j)*dbeta(j)
     &             + buf(k2+j)*dalpha(j)
              
	      vker(j,1)=buf(k1+j)
	      vker(j,2)=buf(k2+j)
	      vker(j,3)=buf(j)
	      
              kdscst=k3
            end do
          else
            do j = 1, knot
              intg(j) = buf(j)*drho(j) + buf(k1+j)*dbeta(j)
     &             + buf(k2+j)*dalpha(j) + buf(k3+j)*dbetah(j)
     &             + buf(k4+j)*dalphah(j) + buf(k5+j)*deta(j)

	      vker(j,1)=beta(j) * buf(k1+j)	      	      
	      vker(j,2)=alpha(j) * buf(k2+j)
	      vker(j,3)=betah(j) * buf(k3+j)
	      vker(j,4)=alphah(j) * buf(k4+j)
	      vker(j,5)=eta(j) * buf(k5+j)
	      vker(j,6)=dn(j) * buf(j)	      

              kdscst=k6
            end do
          end if
        else             
          if(ifanis.eq.0) then
            do j = 1, knot              
              intg(j) = buf(j)*drho(ind+j) + buf(k1 + j)*dbeta(ind+j)
              
	      vker(j,1)=buf(k1+j)
	      vker(j,2)=buf(j)
	      
	      kdscst=k2
            end do
          else
            do j = 1, knot
              intg(j) = buf(j)*drho(ind+j) + buf(k1+j)*dbeta(ind+j)
     &             + buf(k2+j)*dbetah(ind+j)
     
              vker(j,1)=beta(ind+j) * buf(k1+j)
	      vker(j,2)=betah(ind+j) * buf(k2+j)
	      vker(j,3)=dn(ind+j) * buf(j)
     
              kdscst=k3
            end do
          endif
        endif
c
ccc write out
ccc spheroidal, no aniso: 1=Vs,2=Vp,3=rho
ccc spheroidal, aniso: 1=Vsv,2=Vpv,3=Vsh,4=Vph,5=eta,6=rho
ccc toroidal, no aniso: 1=Vs,2=rho
ccc toroidal, aniso: 1=Vsv,2=Vsh,3=rho
        do j=1,knot
	  write(11,"(7e15.5)")rad(j+ind),(0.5*scale*scale*vker(j,k),k=1,6)
	enddo
	close(11)
	 
	
c                  
c        dw = 0.5d0*wdb*t*scale*scale
c        w2 = wdb + dw
c        T&D scaling suggests kernels here are factor of 4 smaller than T/D
c		 test that integration works with omega and cv kernels as well
c        gv = w
c        gv = cv
        
cwrong!	dgv = 0.5*w*t*scale*scale
cright!	dgv = 0.5*t*scale*scale
c        gve = gv + dgv
		

123	irec=irec+1
      end do 
998   continue
      close(11)
      close(2)      

111   format(13e16.8)
      
      end
c----------------------------------------------------
c----------------------------------------------------
c----------------------------------------------------
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

