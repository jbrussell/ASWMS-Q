      program eig_recover
c 
c     eig_recover allows reading and rewriting of eig files
c     for which mineos was aborted.  Specify last angular
c     order (l) value to recover.
c      jbg 9/30/92
c
      include 'parameter.h'
c
      real*4 wmin, wmax
      real*4 radius(nknot)
      real*4 dn(nknot), pv(nknot), sv(nknot), ph(nknot), sh(nknot)
      real*4 eta(nknot)
      real*4 qalpha(nknot), qbeta(nknot)
      real*4 abuf(maxbyte+3)
c
      real*8 wd,h1d,gvd
c
      character*256 fileo
c
      common/ablk/nn,ll,wd,h1d,gvd,buf(nknot6)
      equivalence(nn,abuf)
c
c     get the name of the eigenfunction file
c
      print*,' Enter the pathname of the eig file to be fixed.'
      read(*,'(a)')fileo
      open(unit=2, file=fileo, form = 'unformatted',
     &   access = 'sequential')
csun   , recl = 28000)
      call kblnk(fileo,k)
      fileo = fileo(1:k)//'_fix'
      open(unit=8, file=fileo, form = 'unformatted',
     &   access = 'sequential')
csun   , recl = 28000)
c
c     read and rewrite top of header
c
      read(2) jcom, wmin, wmax, llmin, llmax, wgrav
      read(2) knot,nic,noc,ifanis,tref,
     &  (radius(i),i=1, knot), (dn(i), i=1,knot),
     &  (pv(i), i=1,knot),(sv(i), i=1,knot),
     &  (qalpha(i),i=1,knot),(qbeta(i),i=1,knot),
     &  (ph(i), i=1,knot),(sh(i), i=1,knot),
     &  (eta(i), i=1,knot)
c
      write(8) jcom, wmin, wmax, llmin, llmax, wgrav
      write(8) knot,nic,noc,ifanis,tref,
     &  (radius(i),i=1, knot), (dn(i), i=1,knot),
     &  (pv(i), i=1,knot),(sv(i), i=1,knot),
     &  (qalpha(i),i=1,knot),(qbeta(i),i=1,knot),
     &  (ph(i), i=1,knot),(sh(i), i=1,knot),
     &  (eta(i), i=1,knot)
c
      if (jcom .eq. 2) then
        knot = knot-noc
      endif
c
      if (jcom .eq. 3) then
        nvec=6*knot+8
      else
        nvec=2*knot+8
      endif
c
c     get n,l values for first mode to fix
c
      print*,'Enter last l value to recover:'
      read(*,*) lf
c
c     read modes and correct as neccessary
c
2     read(2,end=999)(abuf(i),i=1,nvec)
      if(ll.le.lf) then
        write(8)(abuf(i),i=1,nvec)
      else
        print*,'Halted prior to ',nn,ll
        go to 999
      end if
      go to 2
999   continue
      close(2)
      close(8)
      stop
      end
  
