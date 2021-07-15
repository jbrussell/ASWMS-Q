      program mineos_table
c 
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c     mineos_table restructures the formatted output
c     file from the alliant - that is it assumes mineos_q has
c     been run already.  It writes intermediate files to the disk
c     which should be automatically deleted when the program
c     terminates
c
c     this version assumes that mineos_strip has been run - 
c     and that all the eigenfunctions deeper than 800 km have
c     been thrown out
c
c     this version has been adapted to incorporate radial modes as well
c
c     this version keeps the model parameters in the header, 
c     and includes the q model parameters as well
c
c     7/2/90 - file structure of mode header file changed to allow
c              for anisotropic models and Q corrections in mineos
c              record length changed from 24000 to 28000
c     04/22/91 minor bug fixed (concerning overlapping frequencies,
c              output of rq) pp
c     11/17/91 input of Q-model from mineos_q corrected pp
c     11/12/96 skips modes with Q=0 -- usually stonely modes that
c              mineos got stuck on.  jbg
c
      include 'parameter.h'
c
      real*8 pi,bigg,rhobar,rad, fl
c
      real*4 abuf,buf,w,gv,q,cv,dum1
      real*4 wmin, wmax, wwmin, wwmax
      real*4 radius(nknot), r(nknot), sql(maxmodes)
      real*4 t1(maxmodes), t2(maxmodes), t3(maxmodes), t4(maxmodes)
      real*4 dn(nknot) ,pv(nknot), ph(nknot), sv(nknot), sh(nknot)
      real*4 qa(nknot), qb(nknot), eta(nknot)
      real*4 rq(nknot), qalpha(nknot), qbeta(nknot)
      real*4 old(maxn,maxll)
c
      integer*4 lmin, lmax, nic, noc, llmin, llmax
      integer*4 point(nknot), jump, rec_num, rlen
      integer*4 nn,ll,i,j,k,ifirst
c
      character*256 filen,fileo,q_file
c
      logical op
      dimension abuf(maxbyte4)
      common/blk/nn,ll,w,q,dum1,buf(nknot4)
      equivalence(nn,abuf(1))
c
      data old /maxold*0./
      data llmin, llmax /100000,0/
      data wwmin, wwmax /1000.,0./
      data bigg,rhobar/6.6723d-11,5515.0d0/
      data pi /3.14159265358979d0/
      data tol /1.e-5/
      data op /.false./
c
c     query the user once about output and other stuff
c
      print*,' Enter the pathname of the output mode table.'
      read(*,'(a)')fileo
      print*,' Enter the best guess for the max number of modes'
      print*,'  please over-estimate'
      read*, nmode
      if (nmode .gt. maxmodes) then
        print*,' your best guess exceeds the array dimensions'
        print*,' please increase maxmodes'
        stop
      endif
c
c     establish min & max limits for modes in frequency & angular order
c
      print*,' Enter min and max frequencies to be reformatted (mhz)'
      read*, frmin, frmax
      frmin = frmin*pi/500.
      frmax = frmax*pi/500.
      print*,' Enter min and max angular orders to be reformatted '
      read*, lllmin, lllmax
c
c     query user (once only) about the mineos_q file & read the qmodel
c
      print*,'mineos_q output file for this mode set'
      read(*,'(a)') q_file
      open(9, file=q_file)
      read(9,129) nqk
      read(9,130) (rq(i),qbeta(i),qalpha(i),i=1,nqk)
129   format(i4)
130   format(3f11.0)
c
c     open input mineos_strip file
c     loop to allow many files to be condensed into one
c
      modes = 0
      ifirst = 0
 5    print*,'mineos_strip output file to be reformatted.'
      read(*,'(a)') filen
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
c     read the two header records and define normalization constants
c
      read(1) jcom, wmin, wmax, lmin, lmax, wgrav
c     N.B. tref = [sec/2 pi]
      read(1) knot, nic, noc, ifanis, tref, (radius(i), i=1,knot),
     &        (dn(i), i=1,knot),(pv(i), i=1,knot),(sv(i), i=1,knot),
     &        (qa(i), i=1,knot),(qb(i), i=1,knot),(ph(i), i=1,knot),
     &        (sh(i), i=1,knot),(eta(i), i=1,knot)
c.... q calculation performed in mineos <-> tref > 0
      if (tref .gt. 0.) then
        do ii = 1, knot
          rq(ii) = radius(ii)
          qalpha(ii) = qa(ii)
          qbeta(ii) = qb(ii)
        end do
        nqk = knot
      endif
      print*, ' mode type = ',jcom,' lmin = ',lmin,' lmax = ', lmax
      print*, ' wmin = ', wmin,' wmax = ', wmax
      print*, ' knots = ',knot
c
c     do the appropriate bookkeeping the first time through
c
      if (ifirst .eq. 1) then
        rn = radius(knot)
        rad = dble(rn)
        wn = sngl(dsqrt(pi*bigg*rhobar))
        accn = sngl(1.0d+20/(rhobar*rad**4))
c
c       strip all the knots deeper than 800 km
c
        kstrip = 0
        do i = 1, knot
          if ((rn - radius(i)) .le. 800000.) then
            kstrip = kstrip + 1
            r(kstrip) = radius(i)/rn
            point(kstrip) = (kstrip-1)*nmode + 1
          endif
        end do
c
c       open temporary files for storage -- now put in same 
c       directory as output, avoiding filling home space and
c       allowing for simultaneos runs.
c
        call kblnk(fileo,k)
        filen = fileo(1:k)//'_scratch.h'
c        open(3,file=filen,access='sequential',form='unformatted')
        open(3,status='scratch',access='sequential',form='unformatted')
        if (jcom .eq. 3) then
           rlen = 16
        else
           rlen = 8
        endif
        filen = fileo(1:k)//'_scratch.e'
c        open(4,access='direct',form='unformatted',file=filen,recl=rlen)
        open(4,access='direct',form='unformatted',status='scratch',recl=rlen)
c
c       set up bookkeeping devices
c         because the eigenfunctions have been stripped already,
c         the number of knots = kstrip
c
        print*,' stripped knots = ',kstrip
        if (jcom .eq. 3) then
          nvec = kstrip*4 + 5
        else
          nvec = kstrip*2 + 5
        endif
      endif
      k1 = kstrip
      k2 = kstrip*2
      k3 = kstrip*3
c
      wwmin = amin1(wmin, wwmin)
      wwmax = amax1(wmax, wwmax)
c
c     read modes and write them to temporary storage
c       assumes mineos_strip output file
c       for toroidals, read w and wp
c       for radials, read u and up
c       for spheroidals, read u, up, v, vp 
c
 20   continue
      read(1,end=5) (abuf(i), i = 1, nvec)
c
c     look for spurious mode
c
      if (nn .lt. 0) then
        print *, ' apparent error in mode calculation'
        print *, nn,' ',ll
        if (nn.lt.-1 .and. nn.gt.-11) then
          print*,'n=-10: skipping next mode as well'
          read(1,end=5)
        end if
        go to 20
      else if (nn .gt. 0) then
        wlast = old(nn,ll+1)
        if (abs(w-wlast) .lt. 1.e-7*w) then
          print*,' skipping overlapping eigenfrequency'
          print *, nn,' ',ll
          go to 20 
        endif
      endif
      if (ll .ge. lllmin .and. ll .le. lllmax) then
        if (w .ge. frmin .and. w .le. frmax) then
          if (old(nn+1,ll+1) .ne. w) then
            old(nn+1,ll+1) = w
            read(9,*) nq, lq, wq, qq, phiq, cv, gv
            if (nn .ne. nq .and. ll .ne. lq) then
              print*,' error in q file match', nn, ll, nq, lq
              stop
            endif
c
c           skip modes with Q = 0 (i.e. Stonely)
c
            if (qq.lt.tol) then
              print *, ' Mode Q = 0, Skipping ',nn,' ',ll
              go to 20
            end if
c
            llmin = min0(ll, llmin)
            llmax = max0(ll, llmax)
            write(3) nn, ll, w, qq, phiq, cv, gv
            if (jcom .eq. 3) then
              do j = 1, kstrip
                index = j
                jump = point(j) + modes
                write(4, rec = jump) buf(index), buf(index+k1),
     &             buf(index+k2), buf(index+k3)
              end do
            else
              do j = 1, kstrip
                index = j
                jump = point(j) + modes
                write(4, rec = jump) buf(index), buf(index+k1)
              end do
            endif
            modes = modes + 1
            if (modes .ge. nmode) then
              print*,' est. number of modes exceeded'
              print*,'modes > nmode',modes,nmode
              stop
            endif
          else
            print *, ' mode',nn,' ',ll,' overlaps with previous mode'
            read(9,*) nq, lq, wq, qq, phiq, cv, gv
            if (nn .ne. nq .or. ll .ne. lq) then
              backspace(9)
            endif
          endif
        else
          print*,'Frequency out of range of input bounds: ',w*500/pi,
     &           nn,' ',ll
          read(9,*) nq, lq, wq, qq, phiq, cv, gv
          if (nn .ne. nq .or. ll .ne. lq) then
            backspace(9)
          endif
        endif
      endif
      go to 20
c
c
100   continue
c
c     rewind temporary files
c
      close(9)
      rewind(3)
c
c     now that we are finished reading in mineos_strip files,
c       we are going to condense the output into the 
c       mineos_table output files
c
      print*,' now writing mode table for ',modes,' modes'
c
      call kblnk(fileo,k)
      filen = fileo(1:k)//'_hdr'
      rlen = modes*4
      open(unit=2,file=filen,form='unformatted',access='sequential')
      open(unit=7,file=fileo,form='unformatted',access='direct',
     &     recl=rlen)
     
c JBR 9/23/15
ccc JBR      open(unit=8,file='table.log',form='formatted',access='sequential')
      
c
c   write some useful stuff to the header file and set up l(l+1) norm
c
      frminout = frmin*500./pi
      frmaxout = frmax*500./pi
c
      write(2) jcom,modes,frminout,frmaxout,llmin,llmax
      write(2) rn, wn, accn, ifanis, tref
      write(2) kstrip, (r(i), i = 1, kstrip)
      write(2) knot,(radius(i),i=1, knot), (dn(i), i=1,knot),
     &         (pv(i), i=1,knot), (sv(i), i=1,knot),
     &         (ph(i), i=1,knot), (sh(i), i=1,knot),
     &         (eta(i), i=1,knot)
      write(2)nqk,(rq(i),i=1,nqk),(qalpha(i),i=1,nqk),(qbeta(i),i=1,nqk)
      if (jcom .ne. 1) then
        do i = 1, modes
          read(3) nn, ll, w, q, phi, cv, gv
          alpha = 0.5*w/q
          write(2) nn, ll, w, q, alpha, phi, cv, gv
          fl = dble(ll)
          sql(i) = sngl(1.d0/dsqrt(fl*(fl + 1.d0)))
        end do
      else
        do i = 1, modes
          read(3) nn, ll, w, q, phi, cv, gv
          alpha = 0.5*w/q
          write(2) nn, ll, w, q, alpha, phi, cv, gv
        end do
      endif
      close(2)
      close(3,status='delete')
c
c   loop over knots, and within each knot loop over modes
c     t1 - w(toroidal) or u(spheroidal or radial)
c     t2 - wp(toroidal) or up(spheroidal or radial)
c     t3 - v(spheroidal)
c     t4 - vp(spheroidal)
c
c   in the mode table, the v, vp, w, and wp eigenfunctions are 
c     normalized by 1/sqrt(ll(ll+1))
c     this is the only normalization performed by mineos_table
c     the eigenfunctions are further normalized in the summation routine
c
      if (jcom .eq. 3) then
        do k = 1, kstrip
          npoint = 4*(k-1) + 1
          rec_num = point(k)
          do n = 1, modes
            read(4, rec = rec_num) t1(n),t2(n),t3(n),t4(n)
            t3(n) = t3(n)*sql(n)
            t4(n) = t4(n)*sql(n)
            rec_num = rec_num + 1
          end do
          write(7, rec = npoint) (t1(m), m=1, modes)
          npoint = npoint + 1
          write(7, rec = npoint) (t2(m), m=1, modes)
          npoint = npoint + 1
          write(7, rec = npoint) (t3(m), m=1, modes)
          npoint = npoint + 1
          write(7, rec = npoint) (t4(m), m=1, modes)
          print*,' radius ',r(k), ' written to output'
        end do
      elseif (jcom .eq. 2) then 
        do k = 1, kstrip
          npoint = 2*(k-1) + 1
          rec_num = point(k)
          do n = 1, modes
            read(4, rec = rec_num) t1(n),t2(n)
            t1(n) = t1(n)*sql(n)
            t2(n) = t2(n)*sql(n)
            rec_num = rec_num + 1
ccc JBR            write(8, *) t1(n), t2(n)
          end do
          write(7, rec = npoint) (t1(m), m=1, modes)
          npoint = npoint + 1
          write(7, rec = npoint) (t2(m), m=1, modes)
          print*,' radius ',r(k), ' written to output'
        end do
      else
        do k = 1, kstrip
          npoint = 2*(k-1) + 1
          rec_num = point(k)
          do n = 1, modes
            read(4, rec = rec_num) t1(n),t2(n)
            rec_num = rec_num + 1
          end do
          write(7, rec = npoint) (t1(m), m=1, modes)
          npoint = npoint + 1
          write(7, rec = npoint) (t2(m), m=1, modes)
          print*,' radius ',r(k), ' written to output'
        end do
      endif
c
c   all done - close up shop and leave the country
c
      close(4,status='delete')
      close(7)
      close(8)
      stop
      end
  
