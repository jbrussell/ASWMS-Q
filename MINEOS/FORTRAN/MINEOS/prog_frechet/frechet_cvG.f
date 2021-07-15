      program frechet_cvG
c
c     program to calculate frechet kernels for phase velocity for a specified mode branch
c
c     based on mineos_test
c     based on frechet_gv
c     expects a mineos_frechet outputfile
c
c     outputs a sequential access file with a single branch of gv kernels
c     format generally similar to frechet format, accept it is sequential, and only
c     the branch file information for the selected branch is output.  Kernels still
c     need to be scaled by multiplying by scale^2
c     JBG March 2008
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
c
c       Big rewrites by NJA, 2014 to calculate dc/dm
c
C23456789112345678921234567893123456789412345678951234567896123456789712
c
      implicit real*4 (a-h,o-z)
c
      include 'parameter.h'
c
      real*4 rad(nknot), buf(nknot6), buf2(nknot6), buf3(nknot6)
c
c     original model parameters
c
      real*4 dn(nknot), alpha(nknot), beta(nknot)
      real*4 alphah(nknot),betah(nknot),eta(nknot)
c
      real*4 interple, scale, scale2, temp(4), kest, dw
	  real*8 del, diffcv, wdb1, wdb2, delta, derivdbl
c
      integer*4 nn, ll, knot, ifirst
      integer*4 nm(nbranch)
      integer*4 nnb(maxl,nbranch), llb(maxl,nbranch)
c
      integer*4 kntdsc (maxdisc)                 
c
      character*256 m_file,o_file,b_file,t_file
c
c     open frechet file and read header records 
c     - beginning of loop over frechet files
c
      print*,'Calculate phase velcoity kernels for a mode branch'
  5   print*,' Enter name of input frechet file'
      read(*,'(a)')m_file 
      if (m_file .eq. ' ') then
        close(2)
        stop
      endif
c
      print*,' Specify mode branch (0=fundamental)'
      read(*,*) ndob
      idob = ndob+1
c
      print*,' Enter name of output file'
      read(*,'(a)')b_file 
c      open(11,file=b_file,form='unformatted',access='sequential')
c
c    debugging files
c      open(12,file='w_cv_gv.old')
c      open(13,file='w_cv_gv.new')
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
      jrec = irec



c
c     write output header
c
      open(11,file=b_file,form='unformatted',access='direct',
     +      recl=nrec)
      write(11,rec=1) jcom,nmodes,nnmin,nnmax,llmin,llmax,ifanis,nrecw
      write(11,rec=2) ksave,nic,noc,ifanis,tref,scale,scale2,ndisc
      write(11,rec=3) (rad(i),i=1,ksave), (kntdsc (i), i = 1, ndisc)
      write(11,rec=4) (dn(i), i=1,ksave)
      write(11,rec=5) (alpha(i), i=1,ksave)
      write(11,rec=6) (beta(i), i=1,ksave)
      if(ifanis.ne.0) then
        write(11,rec=7) (alphah(i), i=1,ksave)
        write(11,rec=8) (betah(i), i=1,ksave)
        write(11,rec=9) (eta(i), i=1,ksave)
        nextrec=10
      else
        nextrec = 7
      end if 
      numberb = 1
      nmdob = nm(idob)-1
c     output number of branches (always 1), number of modes in that branch (input minus 1),
c     and branch number.  Note that number of modes could be less that nm-1 if -10's are encountered
      write(11,rec=nextrec) numberb, nmdob, ndob
      write(11,rec=nextrec+1) (llb(jj,idob), jj = 1,nmdob)
	  
cc %%%%% 
     

c      write(*,*) jcom,nmodes,nnmin,nnmax,llmin,llmax,ifanis,nrecw
c      write(*,*) ksave,nic,noc,ifanis,tref,scale,scale2,ndisc
c      write(*,*) (rad(i),i=1,ksave), (kntdsc (i), i = 1, ndisc)
c      write(*,*) (dn(i), i=1,ksave)
c      write(*,*) (alpha(i), i=1,ksave)
c      write(*,*) (beta(i), i=1,ksave)
c      if(ifanis.ne.0) then
c        write(101,rec=7) (alphah(i), i=1,ksave)
c        write(101,rec=8) (betah(i), i=1,ksave)
c        write(101,rec=9) (eta(i), i=1,ksave)
c        nextrec=10
c      else
c        nextrec = 7
c      end if 
c      numberb = 1
c      nmdob = nm(idob)-1
cc     output number of branches (always 1), number of modes in that branch (input minus 1),
cc     and branch number.  Note that number of modes could be less that nm-1 if -10's are encountered
c      write(101,rec=nextrec) numberb, nmdob, ndob
c      write(101,rec=nextrec+1) (llb(jj,idob), jj = 1,nmdob)
cc %%%%%%      

      jrec2 = nextrec + 1
	      
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
      print*, ' knots = ',knot
      print*, ' scalings: ', scale, scale2
      print*, 'nocor, ndisc: ',nocor,ndisc
	
      if (jcom .eq. 2) then
        knot = nocor
        ind = noc
        kind = index*knot + ndisc
      else
        ind = 0
        kind = index*knot + ndisc
      endif
c
c     begin reading frechet files -- organized along branches.  Skips over 
c     unwanted records 1 branch at a time
c
      ibskp = 0
      do ii=1,idob-1
        ibskp = ibskp + nm(ii)
      end do

c      open(54,file='gvnew_test')
      print*,ndob,idob,nm(idob)
      do jj=1,nm(idob)-1
        irec = jrec + ibskp + jj
        read(2,rec=irec) nn1,ll1,w1,qq1,gv1,cv1,
     *                   (buf(kk), kk = 1, kind)
c        write(44,*)'1: ',nn1,ll1,w1,irec
	
c	if(ll1.eq.160)then
c	  open(45,file='test_buf_1')
c	  do kk=1,knot
c	    write(45,"(4(e15.5))")rad(kk+noc),buf(kk),buf(kk+knot),
c     *        buf(kk+2*knot)
c	  enddo
c	  close(45)
c	endif
		
		if (nn1 .ne. (idob-1)) then
           	print*,'Confused on branch #:',nn1,idob-1
           	stop
		endif
		
        if (nn1 .ge. 0) then
15        irec = irec + 1	 
          read(2,rec=irec) nn2,ll2,w2,qq2,gv2,cv2,
     *                   (buf2(kk), kk = 1, kind)
c          write(44,*)'2: ',nn2,ll2,w2,irec
	  		if(nn2 .lt. 0) go to 15
          	if(nn2 .ne. nn1) then
             	print*,'nn2 .ne. nn1:',nn2,nn1
             	print*,' only ok if end of branch '			 
             	go to 100
	     endif
		  
c         found two neighboring modes -- interpolate w, cv, gv d(cv)/dm, and then d(gv)/dm
c         based on Rodi et al, 1975, BSSA.  Terms with "3" label are the new interpolated values
          if((ll2-ll1).ne.1) print*,'ll1,ll2,',ll1,ll2
		  
          wdb2 = dble(w2)
          wdb1 = dble(w1)
          del=0.5d0*dlog(wdb2/wdb1)
          w3 = w1*exp(del)
          w3 = (w2+w1)/2.
          delta = (wdb2-wdb1)
          qq3 = (qq2+qq1)/2.
          gv3 = (gv2+gv1)/2.
          cv3 = (cv2+cv1)/2.


	  do kk=1,kind

c              dcv1=((cv1**2.)/(gv1))*buf(kk)
c              dcv2=((cv2**2.)/(gv2))*buf2(kk)
c for G kernels
              dcv1 = buf(kk)
              dcv2 = buf2(kk) 
              dcv3=(dcv1+dcv2)/2.
			  
              buf3(kk)=dcv3  

          end do
          jrec2 = jrec2 + 1

          write(11,rec=jrec2) nn1,ll1,w3,qq3,gv3,cv3,
     &                        (buf3(kk), kk = 1,kind)

c!temporary!!!
c          write(11,rec=jrec2) nn2,ll2,w2,qq2,gv2,cv2,
c     &                        (buf3(kk), kk = 1,kind)

c          write(12,*) nn2,ll2,w2, gv2, cv2
c          write(13,*) nn1,ll1,w3, gv3, cv3
100       continue
        endif
      end do

      close(2)
      close(11)
c      close(12)
c      close(13)
      
c      close(54)
      
      end

