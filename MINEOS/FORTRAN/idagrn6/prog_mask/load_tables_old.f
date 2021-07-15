      subroutine load_source(io,source_file,
     &                       jcom,modes,llmin,llmax)

cc ***************** DOESN'T WORK ****************
cc --- JBR : Fails at command 'read(4) a_desc'

cc
cc loads arrays for the source strains
cc
      include 'parameter5.h'

      integer*4 ip_s(0:maxl,0:maxn),ip_r(0:maxl,0:maxn)
      real*4 wwmin,wwmax,dum1,dum2,dum3,tref
      real*4 f0_temp(maxmodes),f1_temp(maxmodes),f_temp(4,3,maxmodes)
      real*4 x1, x2
      real*4 r(nknot), r1
      character*256 source_file

      common/scl/x1,r1,x2,f
      common/ipoint/ip_s,ip_r
      common/tempor/f_temp,f0_temp,f1_temp

cc...source_strains mode table_hdr
cc   dont read radii,velocities and densities

      print*,' Loading eigenfunctions for source strain'

      call kblnk(source_file,k)

      open(io,form='unformatted',file=source_file(1:k)//'_hdr',
     &     access = 'sequential')
      read(io) jcom, modes, wwmin, wwmax, llmin, llmax
      print*, 'jcom, modes...:', jcom, modes, wwmin, wwmax, 
     &     llmin, llmax
      read(io) dum1, dum2, dum3, ifanis, tref
      read(io) knots, (r(i), i = 1, knots)
      read(io) nrad
      read(io)
      if (modes .gt. maxmodes) then
        print*,' array dimension maxmodes exceeded'
        stop
      endif
      if (llmax .gt. maxl) then
        print*,' array dimension maxl exceeded'
        stop
      endif

cc...initialize pointer array
      do nn = 0, maxn
        do ll = 0, maxl
          ip_s(ll,nn) = 0
        enddo
      enddo

cc...determine modal content of header file
cc   set up pointer
      do ii = 1,modes
          read(io) nn,ll
          ip_s(ll,nn) = ii
      enddo
      close(io)

cc...determine pointers
c
c     determine bracketing depth 
c       if the source lies on a model knot,
c       make it deeper by  .1 km
c
      do 15 i=1,knots
         if (r(i) .eq. r1) then
           r1 = r1 - .0001
           r0 = r0 - .0001d0
           nbrack = i - 1
           go to 16
         elseif (r(i) .gt. r1) then
           nbrack = i - 1
           go to 16
         endif
   15 continue
   16 continue
      x1 = r(nbrack)
      x2 = r(nbrack+1)

c
cc... now open mode table - this assumes one using mineos_format
c     kmode determines the record length: 4 bytes for every mode
c
      kmode = 4*modes
      open(11,form='unformatted',file=source_file,
     &    access='direct',recl=kmode)

      if (jcom .eq. 3) then
        npoint = 4*(nbrack - 1) + 1
      else
        npoint = 2*(nbrack - 1) + 1
      endif
c
c     read all of the modes for the bracketing depths
c
      do 30 j = 1, 3, 2
        read(11, rec=npoint) (f_temp(1,j,k), k = 1, modes)
        npoint = npoint + 1
        read(11, rec=npoint) (f_temp(2,j,k), k = 1, modes)
        if (jcom .eq. 3) then
          npoint = npoint + 1
          read(11, rec=npoint) (f_temp(3,j,k), k = 1, modes)
          npoint = npoint + 1
          read(11, rec=npoint) (f_temp(4,j,k), k = 1, modes)
        endif
        npoint = npoint + 1
   30 continue

      close(11)
      return
      end



      subroutine load_receiver(io,source_file,
     &                         jcom,modes,llmin,llmax)
cc
cc load surface eigenfunction for receiver scalars
cc
      include 'parameter5.h'

      integer*4 ip_s(0:maxl,0:maxn),ip_r(0:maxl,0:maxn)
      real*4 f0_temp(maxmodes),f1_temp(maxmodes),f_temp(4,3,maxmodes)
      real*4 rad(nknot), dn(nknot), pv(nknot), sv(nknot)
      real*4 wwmin,wwmax,dum1,dum2,dum3,tref
      real*4 dt,wmhz
      character*1 lcomp
      character*256 source_file
      logical ocean, ldisp, lida

      common/limits/dt,wmhz,lcomp,ldisp,lida
      common/ipoint/ip_s,ip_r
      common/tempor/f_temp,f0_temp,f1_temp

      print*,' Loading eigenfunctions for receiver scalars'

      call kblnk(source_file,k)

      open(io,form='unformatted',file=source_file(1:k)//'_hdr',
     &     access = 'sequential')
      read(io) jcom, modes, wwmin, wwmax, llmin, llmax
      read(io) dum1, dum2, dum3, ifanis, tref
      read(io) knots
      read(io) nrad, (rad(i), i = 1, nrad), (dn(i), i = 1, nrad),
     &         (pv(i), i=1,nrad), (sv(i), i=1,nrad)
      read(io)
c
c     if this is an oceanic model we need to know the number of 
c     knots under water (obs)
c
        if (sv(nrad) .eq. 0.e0) then
          do i = nrad,1,-1
            if (sv(i) .ne. 0.e0) then
              obs = nrad - i
              go to 5
            endif
          enddo
    5     ocean = .true.
        endif
        if (modes .gt. maxmodes) then
          print*,' array dimension maxmodes exceeded'
          stop
        endif
        if (llmax .gt. maxl) then
          print*,' array dimension maxl exceeded'
          stop
        endif

cc...initialize pointers arrays
      do nn = 0, maxn
        do ll = 0, maxl
          ip_r(ll,nn) = 0
        enddo
      enddo

cc...determine modal content of header file
      do ii = 1,modes
          read(io) nn,ll
          ip_r(ll,nn) = ii
      enddo
      close(io)


cc...set up pointers
c
      if (ocean) then
        print*,' ocean model - seismometer placed on ocean floor'
        if (jcom .eq. 3) then
          npoint = 4*(knots - obs - 1) + 1
        else
          npoint = 2*(knots - obs - 1) + 1
        endif
      else        
        if (jcom .eq. 3) then
          npoint = 4*(knots - 1) + 1
        else
          npoint = 2*(knots - 1) + 1
        endif
      endif

	write(*,*) 'mmmm ',obs,npoint,modes,jcom
c
cc... now open mode table - this assumes one using mineos_format
c     kmode determines the record length: 4 bytes for every mode
c
      kmode = 4*modes
      open(11,form='unformatted',file=source_file,
     &    access='direct',recl=kmode)

      read(11, rec=npoint) (f0_temp(k), k = 1, modes)
c
c     read v at the surface if these are spheroidal modes 
c     for horizontal polariztion
c
      if ((jcom .eq. 3) .and. (lcomp .ne. 'Z')) then
        npoint = npoint + 2
        read(11, rec=npoint) (f1_temp(k), k = 1, modes)
      endif

      close(11)
      return
      end



      subroutine load_bm(modes,jcom2,sort,lpart)
cc
cc reads mask or branch files
cc
      include 'parameter5.h'

      integer*4 modes
      integer*4 nnb(maxl,nbranch),llb(maxl,nbranch),numb(nbranch)
      real*4 sort(0:maxl,0:maxn),stuff(4)
      character*256 m_file, a_desc, b_file
      logical lmask,lbranch,lpart

      common/mask/ m_file, lmask
      common/branch/ b_file, lbranch, nbr


c MASK AND BRANCH FILES
c     using the mask file is the preferred way to limit summation
c       this array now stored as a real array, to allow for amplitude
c
c     it is assumed that the user is using either a mask file
c       or a branch file - not both!
c
      if (lmask) then
        do ii = 0, maxn
          do jj = 0, maxl
            sort(jj,ii) = 0.
          end do
        end do
        open(4,file=m_file,form='unformatted',access='sequential')
        read(4) jcom2, nmodes, nlow, nup, llow, lup        
        read(4) a_desc
        read(4)
        read(4)
        if (nmodes .ne. modes) then
          print*,' mismatch in mode count: ', modes, nmodes
        endif
        print*,' '
        print*, a_desc
        print*,' '
        do ii = nlow, nup
          read(4) kk,(sort(jj,kk), jj = llow, lup)
          if (kk .ne. ii) then
            print*,' potential error in reading mask file'
          endif
        end do
        close(4)
        lpart = .true.
      elseif (lbranch) then
        do ii = 0, maxn
          do jj = 0, maxl
            sort(jj,ii) = 0.
          end do
        end do
        open(4,file=b_file,form='unformatted',access='sequential')
        read(4) jcom2, nmodes, nlow, nup, llow, lup
        read(4) nb, (numb(ii), ii = 1, nb)
        if (nmodes .ne. modes) then
          print*,' mismatch in mode count: ', modes, nmodes
        endif
        print*,' reading information for ', nb,' branches'
        do ii = 1, nb
          read(4) (nnb(jj,ii), llb(jj,ii), 
     &                (stuff(kk),kk = 1,4),jj = 1,numb(ii))
        end do
        close(4)
c
c       nbr picks out a specific branch
c
        if (nbr .ne. 0) then
          print*,' selecting modes for branch: ',nbr
          il = 0
          do im = 1, numb(nbr)
            n1 = nnb(im,nbr)
            l1 = llb(im,nbr)
            sort(l1,n1) = 1.
            il = il + 1
          end do
        else
          print*,' selecting modes for every branch'
          il = 0
          do ib = 1, nb
            do im = 1, numb(ib)
              n1 = nnb(im,ib)
              l1 = llb(im,ib)
              sort(l1,n1) = 1.
              il = il + 1
            end do
          end do  
        endif
        print*,' # of modes selected: ', il
        lpart = .true.
      endif        
      
      return
      end
