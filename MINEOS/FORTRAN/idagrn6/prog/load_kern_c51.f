      subroutine load_kern_c51(m_file)
c
c     routine to load mode frechet files (from frechet.f)
c
c     returns:		crker(kcrtot,6,maxmodes) -- volumetric kernels
c                       crkers(maxmodes) -- solid surface kernel
c                       crkerm(maxmodes) -- moho kernel
c                       gvi(maxmodes) -- group velocity for each mode
c                       ipnt(maxn,maxl) -- location pointer of each kernel
c                       kcrtot -- # of knots in crust
c                       rmoho -- radius of moho
c                       rsol -- radius of solid surface
c                       rcr,vs,vp,ro -- radii,velocities/densities of
c                            starting model through crust (arrays of 
c                            dimension kcrtot)
c
c     written 9/99 -- JBG
c     01/24/00 -- updated to read in new frechet files, with nrec 
c                 specified in 1st record
c     10/13/04 -- updated to output kernels for two crust layers (seds and solid rock) rather than one
c
c
      implicit real*4 (a-h,o-z)
c
      include 'parameter5.h'
c
      real*4 rad(nknot), buf(nknot6), scale, scale2

      real*4 dn(nknot), alpha(nknot), beta(nknot)
c
      integer*4 nn, ll, knot, kstop,kmoho,ksol
      integer*4 ndisc,kntdsc (50)
c
      character*256 m_file

c     common block for kernel information

      real*4 crker(nlith,6,maxmodes),crkers(maxmodes),
     &       crkerb(maxmodes), crkerm(maxmodes), gvi(maxmodes)
      real*4 vs(nlith),vp(nlith),ro(nlith),rcr(nlith)
      real*4 rsol, rmoho, rbase
      integer*4 kcrtot,ipnt(maxn,maxl),ifanis
      common /kern51/ ifanis,crker,crkers,crkerb,crkerm,gvi,ipnt
      common /omod51/ rsol,rbase,rmoho,vs,vp,ro,rcr,kcrtot
      common/ahdr/ jcom,knots


      print*,' '
      print*,'Loading crust5.1 and/or dbdb5 kernels'
      print*,' '
      do kk = 1,maxmodes
        do jj = 1,6
          do ii = 1,nlith
            crker(ii,jj,kk) = 0.
          enddo
        enddo
        crkers(kk) = 0.
        crkerb(kk) = 0.
        crkerm(kk) = 0.
      enddo
c
C     complex to handle the various record lengths that could be set
c     in frechet.f.  Depends on ifanis and jcom flags, both now in 
c     first record.  To start, assume a short record length (isotropic 
c     toroidal modes) to allow reading of model-type flags.
c
c     no longer true -- now just esimates an nrec, and then reads in
c     the actual one, closing and reopening if necessary
c

      call kblnk(m_file,k)
      m_file = m_file(1:k)
      print*, 'm_file= ', m_file
      index=2
      newvec = index*nknot_t + 6 + maxdisc
      nrec = newvec * 4 
10    open(unit=2,file=m_file,form='unformatted',access='direct',
     +      recl=nrec)
      read(2,rec=1)jcom,nmodes,nnmin,nnmax,llmin,llmax,ifanis,nrecw
      if (nrec.ne.nrecw) then
         close(2)
         nrec = nrecw
         go to 10
      endif
c
      read(2,rec=2) ksave,nic,noc,ifanis,tref,scale,scale2,ndisc
c      print*, ksave,nic,noc,ifanis,tref,scale,scale2,ndisc

      read(2,rec=3) (rad(i),i=1,ksave), (kntdsc (i), i = 1, ndisc)
c      do i=1,ksave
c         print*, i, rad(i)
c      enddo
c      do i=1,ndisc
c         print*, i, kntdsc(i)
c      enddo

      read(2,rec=4) (dn(i), i=1,ksave)
      read(2,rec=5) (alpha(i), i=1,ksave)
      read(2,rec=6) (beta(i), i=1,ksave)
      if(ifanis.ne.0) then
        nextrec = 10
      else
        nextrec = 7
      end if 
      read(2,rec=nextrec) nb
      jrec = nextrec + nb

c     still need value for index

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
     
      rstop = rad(ksave) - 60000.
c
c     Find upper solid boundary, basement, and Moho.  
c       water = vs < 10 m/s
c       sed = 10 < vs < 2000 m/s
c       crust = vs < 4000.0 m/s, vp < 7000.0 km/s

      if (beta(ksave).ge.10.) then
        ksol = ksave
        nsol = -1
        print*,'problem -- no kernel for solid surface'
      else
        do i = ndisc,1,-1
c           print*, 'test: ', ndisc, i, kntdsc(i), beta(kntdsc(i))
         if(beta(kntdsc(i)).ge.10.) then
           ksol = kntdsc(i)
           nsol = i
           rsol = rad(kntdsc(i))
           go to 80
         endif
        enddo
      endif

80    continue

c     arbitrary test to find basement.  vs seems best
  
      do i = nsol,1,-1
c        if(beta(kntdsc(i)).gt.4000..and.alpha(kntdsc(i)).gt.7000.) then

cc    Po Chen
        if(beta(kntdsc(i)).gt.3000.) then
cc    Po Chen
          kbase = kntdsc(i)
          nbase = i 
          rbase = rad(kntdsc(i))
		  print*,'kbase,nbase,rbase,beta',
     &           kbase,nbase,rbase,beta(kntdsc(i))
c                  stop
          go to 85
        endif
      enddo
        
 85   continue
c     arbitrary test to find moho.  vp seems best
  
      do i = nsol,1,-1
c        if(beta(kntdsc(i)).gt.4000..and.alpha(kntdsc(i)).gt.7000.) then
        if(alpha(kntdsc(i)).gt.7000.) then
          kmoho = kntdsc(i)
          nmoho = i 
          rmoho = rad(kntdsc(i))
          go to 90
        endif
      enddo
        
 90   continue

cc     find knot at about 60 km depth - we'll save only up to here
c      NO -- save only up to moho -- only perturb crustal velocities
c
      i = kmoho
      do while (rad(i).gt.rstop)
        i = i-1
      end do
      kstop = i

      kstop = kmoho + 1
      rstop = rad(kstop)

      kcrtot = ksol-kstop 
c      print*, 'kcrtot= ', kcrtot
c     save lithospheric velocities for later use

      j = 1
      do i = kstop,ksol
        vs(j) = beta(i)
        vp(j) = alpha(i)
        ro(j) = dn(i)
        rcr(j) = rad(i)
c        print*, 'in load_kern:', j, rcr(j)
        j=j+1
      enddo
        
c
c     reset start and stop knots for toroidal modes -- buf only holds mantle values
c
      knot = ksave
      nocor = knot - noc
      if (jcom .eq. 2) then
        knot = nocor
        kstop = kstop - noc
      endif
      kind = index*knot + ndisc
      kdscst=index*knot

c
c     begin reading frechet files and load kernels for 0-80 km plus surface and moho
c
      do ii=1,nmodes
        read(2,rec=jrec+ii) nn,ll,w,qq,gv,cv,
     &                      (buf(kk), kk = 1,kind)
        gvi(ii)=gv
c        if(ii.eq.1) print*,nn,ll,w,qq,gv,cv
        ipnt(nn,ll) = ii
        do jj = 1,index
          kk = 0      
          j1 = (jj-1)*knot
          do k = kstop, knot                         
            kk = kk + 1
            crker(kk,jj,ii) = buf(j1+k)
          end do
        end do
        if (nsol.gt.0) then
           crkers(ii) = buf (kdscst + nsol) 
        endif
		crkerb(ii) = buf (kdscst + nbase)
        crkerm(ii) = buf (kdscst + nmoho)

c       DEBUGGING
c
c        if(ll.eq.200 .and. nn.eq.0) then 
c          do kk = 1, knot-kstop + 1
c           k1 = kstop  - 1 + kk
c           k1 = kstop  + noc - 1 + kk
c           print*,rad(k1),beta(k1),(crker(kk,jj,ii), jj=1,index)
c          enddo
c        elseif (ll.eq.54 .and. nn.eq.4) then 
c          do kk = 1, knot-kstop + 1
c           k1 = kstop  - 1 + kk
c           print*,rad(k1),beta(k1),(crker(kk,jj,ii), jj=1,index)
c          enddo
c        endif
c                    
c        if(ll.eq.200 .and. nn.eq.0) then 
c          print*,rad(kntdsc(nmoho)),crkerm(ii)
c          print*,rad(kntdsc(nsol)),crkers(ii)
c         elseif (ll.eq.54 .and. nn.eq.4) then
c           print*,rad(kntdsc(nmoho)),crkerm(ii)
c           print*,rad(kntdsc(nsol)),crkers(ii)
c         endif

      end do

      close(2)
      return
      end

