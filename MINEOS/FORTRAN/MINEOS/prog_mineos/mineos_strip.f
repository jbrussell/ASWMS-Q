      program mineos_strip
c 
c     mineos_strip strips the upper 800 km of eigenfunctions
c     and dumps the rest, as a disk saving proceedure.  
c     This progam may be run directly after
c     mineos or any time before mineos_format
c
c     01/22/87 radial modes added 
c     01/09/91 check for overlapping modes added
c     04/22/91 bugs cleaned up (concerning anisotropic models, overlapping c      c              eigenfrequency test) pp
c     09/21/91 mineos eigenfunction file input of wd,qd,gvd as real*8 
c              instead of real*4 pp
c
      include 'parameter.h'
c
      real*8 wd,qd,gvd
c
      real*4 abuf,buf,w,gv,q,cv
      real*4 wmin, wmax
      real*4 radius(nknot), r(nknot)
      real*4 dn(nknot),pv(nknot),sv(nknot),ph(nknot),sh(nknot),
     &       qa(nknot),qb(nknot),eta(nknot)
      real*4 wlast
      real*4 old(maxn,maxl)
c
      integer*4 lmin, lmax, nic, noc
      integer*4 nn,ll,i,j,k,l,ifirst
c
      character*256 filen,fileo
c
      logical op
      dimension abuf(maxbyte4+3)
      common/blk/nn,ll,wd,qd,gvd,buf(nknot4)
      equivalence(nn,abuf(1))
c
cc      data old /maxold*0.0/
      data op /.false./
c
c     query the user once about output and other stuff
c
      print*,' Enter the pathname of the stripped mode file.'
      read(*,'(a)')fileo
      open(unit=2,file=fileo,form='unformatted',access='sequential')
c
c     open input mode file - this is a loop to allow many files to be
c     condensed into one
c
      ifirst = 0
 5    print*,' Enter the pathname of the mineos file to be stripped.'
      read(*,'(a)')filen
      if (op) then
        close(1)
      endif
      op = .true.
      if (filen .eq. ' ') then
        go to 100
      else
        ifirst = ifirst + 1
      endif
      open(unit=1,file=filen,form='unformatted',access='sequential')
c
c   read the two header records and define normalization constants
c
      read(1) jcom, wmin, wmax, lmin, lmax, wgrav
c     N.B. tref = [sec/2 pi]
      read(1) knot, nic, noc, ifanis, tref, (radius(i),i=1,knot),
     &        (dn(i), i=1,knot),(pv(i), i=1,knot),(sv(i), i=1,knot),
     &        (qb(i), i=1,knot),(qa(i), i=1,knot),(ph(i), i=1,knot),
     &        (sh(i), i=1,knot),(eta(i), i=1,knot)
      print*, ' mode type = ',jcom,' lmin = ',lmin,' lmax = ', lmax
      print*, ' wmin = ', wmin,' wmax = ', wmax
      print*, ' knots = ',knot, nic, noc
      print*, ' radius = ',radius(knot)
c
c   do the appropriate bookkeeping the first time through
c
      if (ifirst .eq. 1) then
c
c       write out the two record header
c
        write(2) jcom, wmin, wmax, lmin, lmax, wgrav
        write(2) knot, nic, noc, ifanis, tref, (radius(i),i=1,knot),
     &           (dn(i), i=1,knot),(pv(i), i=1,knot),(sv(i), i=1,knot),
     &           (qa(i), i=1,knot),(qb(i), i=1,knot),(ph(i), i=1,knot),
     &           (sh(i), i=1,knot),(eta(i), i=1,knot)
c
        rn = radius(knot)
c
c       strip all the knots deeper than 800 km
c
        kstrip = 0
        do i = 1, knot
          if ((rn - radius(i)) .le. 800000.) then
            kstrip = kstrip + 1
          endif
        end do
c
c       set up bookkeeping devices
c
        if (jcom .eq. 3) then
          nvec = knot*4 + 8
        elseif (jcom .eq. 2) then
          knot = knot - noc
          nvec = knot*2 + 8
        else
          nvec = knot*2 + 8
        endif
        index = knot - kstrip
        print*,' saving ',kstrip,' knots'
      endif
c
c     read modes and rewrite them to unit 2
c       assumes raw output file from mineos
c       for toroidals, read w and wp
c       for radials, read u and up
c       for spheroidals, read u, up, v, vp 
c
c     spheroidal modes
c
ccheck      irec = 0
      if (jcom .eq. 3) then
 20     continue
ccheck        irec = irec + 1
ccheck        print*,'reading record',irec
        read(1,end=5) (abuf(i), i = 1, nvec)
c....   put w,q,gv in real*4 array
        w  = sngl(wd)
        q  = sngl(qd)
        gv = sngl(gvd)
ccheck        print*,'record',irec,' read n=',nn,' l=',ll
        if (nn .lt. 0) then
          print *, ' apparent error in mode calculation'
          print *, nn,' ',ll
          if (nn.lt.-1) then
            print*,'n=-10: skipping next mode as well'
            read(1,end=5)
          end if
          go to 20
        else if (nn .gt. 0) then
          wlast = old(nn,ll+1)
          if (abs(w - wlast) .lt. 1.e-7*w) then
            print*,' skipping overlapping eigenfrequency'
            print *, nn,' ',ll
            go to 20 
          endif
        endif
        if (old(nn+1,ll+1) .ne. w) then
          old(nn+1,ll+1) = w
          write(2) nn,ll,w,q,gv,(buf(index+j),j=1,kstrip),
     &      (buf(knot+index+j),j=1,kstrip),
     &      (buf(2*knot+index+j),j=1,kstrip),
     &      (buf(3*knot+index+j),j=1,kstrip)
        else
          print *, ' mode',nn,' ',ll,' overlaps with previous mode'
        endif
        go to 20
c
      else
c
c     toroidal or radial modes
c
 30     continue
        read(1,end=5) (abuf(i), i = 1, nvec)
c....   put w,q,gv in real*4 array
        w  = sngl(wd)
        q  = sngl(qd)
        gv = sngl(gvd)
        if (nn .lt. 0) then
          print *, ' apparent error in mode calculation'
          print *, nn,' ',ll
          read(1,end=5)
          go to 30
        else if (nn .gt. 0) then
          wlast = old(nn,ll+1)
          if (abs(w - wlast) .lt. 1.e-7*w) then
            print*,' skipping overlapping eigenfrequency'
            print *, nn,' ',ll
            go to 30 
          endif
        endif
        if (old(nn+1,ll+1) .ne. w) then
          old(nn+1,ll+1) = w
          write(2) nn,ll,w,q,gv,(buf(index+j),j=1,kstrip),
     &      (buf(knot+index+j),j=1,kstrip)
        else
          print *, ' mode',nn,' ',ll,' overlaps with previous mode'
        endif
        go to 30
c
      endif
c
c
100   continue
c
c     shut down and go home
c
      close(2)
      stop
      end
  
