      program frechet_asc
c
c     dumps ascii version of specified n,l kernels, for plotting nominally 
c
      implicit real*4 (a-h,o-z)
c
      include 'parameter.h'
c
c     horizontal velocities and eta for anisotropic case
c
      real*4 alphah(nknot), betah(nknot), eta(nknot)
c
      real*4 radius(nknot), dn(nknot), alpha(nknot), beta(nknot)
      real*4 abuf(mfrechet), intg(nknot,6), disc(maxdisc)
c
      integer*4 nmt(nbranch)
      integer*4 llt(maxl,nbranch)
      integer*4 kntdsc(maxdisc)
c
      character*256 f_file, p_file
c
c     tricky equivalence statement to read in large blocks of data
c
      common/ablk/n1,l1,w,q1,gv1,cv1,buf(mbuf)
      equivalence(n1,abuf)
c
      data i0 /0/
      data f0 /0.0/
c
c     include some useful numerical constants
c
      include 'numerical.h'
c 
c****************
c     now open direct access frechet file.  Assumes a value
c     for nknot_t  --- new version reads in the actual rec
c     length at end of first rec, and reopens if necessary
c
c
c      Commented out by NJA to fix readonly error
c      print*,' enter name of frechet kernel file'
c      read(*,'(256a)') f_file
c      index=2
c      newvec = index*nknot_t + 6 + maxdisc
c      nrec = newvec * 4 
c10    open(unit=4,file=f_file,form='unformatted',access='direct',
c     +     status='readonly',recl=nrec)
c      read(4,rec=1)jcom2,nmodes,nnmin,nnmax,llmin,llmax,ifanis,nrecw
c      if (nrec.ne.nrecw) then
c         close(4)
c         nrec = nrecw
c         go to 10
c       endif
      print*,' enter name of frechet kernel file'
      read(*,'(256a)') f_file
      index=2
      newvec = index*nknot_t + 6 + maxdisc
      nrec = newvec * 4 
10    open(unit=4,file=f_file,form='unformatted',access='direct',
     +     recl=nrec)
      read(4,rec=1)jcom2,nmodes,nnmin,nnmax,llmin,llmax,ifanis,nrecw
      if (nrec.ne.nrecw) then
         close(4)
         nrec = nrecw
         go to 10
       endif

c
c     if anisotropic, reset record and vector lengths,
c     close, reopen, and reread -- not any more -- still
c     need to set index though
c
c      if (ifanis.eq.1) then
c        if (jcom2 .eq. 2) then
c          index = 3
c          newvec = index*nknot_t + 6 + maxdisc
c        else 
c          index = 6
c          newvec = index*nknot_s + 6 + maxdisc
c        endif
c        nrec = newvec * 4 
c        close(4)
c        open(unit=4,file=f_file,form='unformatted',access='direct',
c     +        status='readonly',recl=nrec)
c        read(4,rec=1) jcom2,nmodes,nnmin,nnmax,llmin,llmax,ifanis 
c        print*,'computing ANISOTROPIC partials'
c
c      else
c        if(jcom2.eq.3) then
c          index=3
c          newvec = index*nknot_s + 6 + maxdisc
cc          nrec = newvec * 4 
c          close(4)
c          open(unit=4,file=f_file,form='unformatted',access='direct',
c     +        status='readonly',recl=nrec)
c          read(4,rec=1) jcom2,nmodes,nnmin,nnmax,llmin,llmax,ifanis
c        end if
c        print*,'computing ISOTROPIC partials'
c      end if
c
c     frechet file should now be open with correct record length
c************
c
c     continue reading header and model records
c
      read(4,rec=2) ksave, nic, noc, ifanis, tref, scale, scale2, ndisc
      read(4,rec=3) (radius(i),i=1,ksave), (kntdsc(i), i=1, ndisc)
      read(4,rec=4) (dn(i), i=1,ksave)
      read(4,rec=5) (alpha(i), i=1,ksave)
      read(4,rec=6) (beta(i), i=1,ksave)
      if (ifanis.eq.1) then 
        read(4,rec=7) (alphah(i), i=1,ksave)
        read(4,rec=8) (betah(i), i=1,ksave)
        read(4,rec=9) (eta(i), i=1,ksave)
        nextrec=10
      else
        nextrec=7
      endif
      read(4,rec=nextrec) nt, (nmt(ii), ii = 1, nt)
      do ii = 1, nt
        irec = nextrec + ii
        read(4, rec=irec) (llt(jj,ii), jj = 1,nmt(ii))
      end do

      if(ifanis.eq.1) then
        if (jcom2 .eq. 2) then
          index = 3
        else 
          index = 6
        endif
      else
        if (jcom2 .eq. 2) then
          index = 2
        else 
          index = 3
        endif
      endif
     
      print*,' jcom = ', jcom2,' # modes = ',nmodes
      print*,' '
      knot = ksave
c
c     now open output file and write a header record
c
      print*,' enter name of output file'
      read(*,'(256a)') p_file
      open(unit=7,file=p_file)
c
      nocor = knot - noc
      if (jcom2 .eq. 2) then
        knot = nocor
      else
        noc = 0
      endif
c
c     nlen is the length of the buf array, nvec is the length of abuf
c
      nlen = index*knot + ndisc
      nvec = 6 + nlen
c
c     ok, now get the mode to read up to
c
2     print*,'Enter n, l value to dump'
      read(*,*)nd,ld
      if(nd .lt.0) go to 999
      write(7,*) nd,ld

      if(jcom2.eq.2) then
        write(7,*) 'ro, sv, sh'
      else
        write(7,*) 'ro, sv, sh, pv, ph, eta'
      endif 
c
      ii=1
4     irec = nextrec + nt + ii
c      read(4,rec=irec,end=999) n,l
      read(4,rec=irec) n,l
      if(l.eq.ld .and. n.eq.nd) then
        read(4,rec=irec) nn, ll, w, qq, gv, cv,
     &                    ((intg(j,ii),j = 1, knot),ii=1,index),
     &                    (disc (j), j = 1, ndisc)
        write(7,*) nn,ll,w*500./pi
        do j = 1,knot
          write(7,1000) radius(noc+j),(intg(j,ii),ii=1,index)
        end do
        write(7,'(a1)') 'discontinuity kernels'
        do j = 1,ndisc
          write(7,'(2e12.4)') radius(kntdsc(j)),disc(j)
        end do
        ii=ii+1
        go to 2
      else
        ii=ii+1
        go to 4
      end if
999   continue
1000  format(7e12.4)
      close(4)
      close(7)
      stop
      end
