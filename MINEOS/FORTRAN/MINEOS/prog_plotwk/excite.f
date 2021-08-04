c
c
c
      subroutine excite(e_file)
c
c     subroutine to open mode excitation file
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop)
      real*4 fill(0:maxmodes), phase(0:maxmodes)
      real*4 group(0:maxmodes), source(0:maxmodes), instr(0:maxmodes)
      real*4 delta, mhz, rad, pi, tpi, xkm, exc, ak, phik, phip
      real*4 dw, dist, t1, t2, w1, w2, tg, lam, angle
      real*4 amp(0:maxmodes), dep
c
      integer*4 ifill(0:maxmodes),cfill(0:maxmodes)
      integer*4 inx(0:maxl,0:maxn),inx2(0:maxmodes,2)
      integer*4 nnb(0:maxl,nbranch),llb(0:maxl,nbranch),nm(nbranch)
      integer*4 e_time(5), nn, ll
      integer*4 jcom, nmodes, nlow, nup, llow, lup
c
      character*4 statn, compn
      character*256 e_file
c
      logical lsummary, lbranch, lonce
c
      common /mode/ inx, modes, inx2
      common /fill/ lsummary, fill, ifill, cfill
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /c_excite/ dist, dep
      common /nparam/ nparam, nextra
      common /branch/ lbranch, nb, nm, nnb, llb
c
      data ra /57.2957795/
      data tpi /6.2831853071796/
      data xkm /111.2/
      data lonce /.true./
c
      mhz = 1000./tpi
      rad = 1.0/mhz  
      pi = tpi/2.
      pi4 = tpi/8.
c
c     open and read header records
c
      open(unit=3,file=e_file,form='unformatted',access='sequential')
      read(3) jcom, num_modes, nup2, nlow2, lup2, llow2
      read(3) (e_time(jj),jj=1,5), statn, compn, delta, dep
      print*,' loading excitation kernels: ',(e_time(ii), ii = 1, 5),
     +       ' ', statn, ' ', compn, ' ', delta, dep
      if (num_modes .ne. nmodes) then
        print*, num_modes, ' modes read in'
      endif
c
      angle = delta/ra
      dist = delta*xkm
c
c     read excitation table
c
      do ii = 1, num_modes
        read(3,end = 25) nn, ll, ak, phik, x1, exc
        index = inx(ll,nn)
        phase(index) = phik
        amp(index) = ak
        fill(index) = exc
        ifill(index) = 0
        modes(index,27) = 0.0
        modes(index,28) = 0.0
c
c       optional check to amplitude
c
c        lam = float(ll) + 0.5
c        ae = ak*cos(lam*angle - pi4 - phik)
c        if (exc .ne. 0.0) then
c          write(6,'(2i4,8f15.8)') nn,ll,ak,phik,phase(index),ae,
c     &          exc,(ae-exc)*100./exc
c        endif
c
c       end optional check
c
      end do
25    continue
      close(3)
c
      if (.not.(lbranch)) then
        print*,' branch information required '
        return
      endif
      do ii = 1, nb
c
c       unwind phase (eliminate 2pi jumps) and fix amplitude errors by interpolation
c
        do im = nm(ii), 2, -1
          n1 = nnb(im,ii)
          l1 = llb(im,ii)
          n2 = nnb(im-1,ii)
          l2 = llb(im-1,ii)
          icx1 = inx(l1,n1)
          icx2 = inx(l2,n2)
          t1 = phase(icx1)
          t2 = phase(icx2)
          nest = nint((t1 - t2)/tpi)
          phase(icx2) = t2 + nest*tpi
        end do
c
c       now calculate phase center derivative and form total group time
c       phip is the center derivative of source phase
c
        do im = 2, nm(ii) - 1
          n0 = nnb(im,ii)
          l0 = llb(im,ii)
          icx0 = inx(l0,n0)
          n1 = nnb(im-1,ii)
          l1 = llb(im-1,ii)
          icx1 = inx(l1,n1)
          n2 = nnb(im+1,ii)
          l2 = llb(im+1,ii)
          icx2 = inx(l2,n2)
          w1 = modes(icx1,1)
          w2 = modes(icx2,1)
          t0 = phase(icx0)
          t1 = phase(icx1)
          t2 = phase(icx2)
          a0 = amp(icx0)
          a1 = amp(icx1)
          a2 = amp(icx2)
c
          if (a2 .eq. -9999.) then
            if (a0 .eq. -9999.) then
              modes(icx0,28) = 0.0
            else
              modes(icx0,28) = modes(icx1,28)
            endif
          elseif (a1 .eq. -9999.) then
            modes(icx0,28) = 0.0
          else  
            dw21 = (w2 - w1)*rad
            phip = (t2 - t1)/dw21
            modes(icx0,28) = phip
          endif
          modes(icx0,27) = dist/modes(icx0,4)
          ifill(icx0) = 1
ccc
ccc       added check to see if the phase suddenly drops to zero
ccc
          if (abs(t1 - t2) .gt. pi/2.) then
            if (abs(t2) .lt. 0.000001) then
              print*,' zeroing phase derivative:', n0, l0, a0, t0, 
     &                 modes(icx0,28)
              phip = 0.0
              modes(icx,28) = 0.0
            endif
          endif
ccc
ccc
ccc
c
c         optional check on phase
c
c         write(6,'(2i4,8f15.8)') n0,l0,a0,t0,modes(icx0,28)
c
c         end optional check
c
        end do
c
c       do the endpoints
c
        n1 = nnb(1,ii)
        l1 = llb(1,ii)
        icx1 = inx(l1,n1)
        n2 = nnb(2,ii)
        l2 = llb(2,ii)
        icx2 = inx(l2,n2)
        modes(icx1,27) = dist/modes(icx1,4)
        modes(icx1,28) = modes(icx2,28)
        ifill(icx1) = 1
c
        n1 = nnb(nm(ii),ii)
        l1 = llb(nm(ii),ii)
        icx1 = inx(l1,n1)
        n2 = nnb(nm(ii)-1,ii)
        l2 = llb(nm(ii)-1,ii)
        icx2 = inx(l2,n2)
        modes(icx1,27) = dist/modes(icx1,4)
        modes(icx1,28) = modes(icx2,28)
        ifill(icx1) = 1
      end do
c
      if (lonce) then
        nextra = nextra + 2
        lonce = .false.
      endif
      call color(nmodes)
c
      return
      end
