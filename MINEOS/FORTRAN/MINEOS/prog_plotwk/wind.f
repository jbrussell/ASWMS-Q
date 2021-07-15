c
c
c
      subroutine wind(gv0,s0,np)
c
c     subroutine to window about group velocity
c
      include 'parameter.h'
c
      character*1 wtype
      real*4 modes(0:maxmodes,nprop)
      real*4 fill(0:maxmodes),  dist, dep
      real*4 gv0, s0, t0, tgv, arg, amp
      integer*4 ifill(0:maxmodes)
      integer*4 cfill(0:maxmodes)
      integer*4 inx(0:maxl,0:maxn)
      integer*4 inx2(0:maxmodes,2)
      logical lsummary
c
      common /mode/ inx, modes, inx2
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /fill/ lsummary, fill, ifill, cfill
      common /c_excite/ dist, dep
c
c     sign convention here is that source (modes(ii,28)) and instrument delay (modes(ii,29))
c       have minus signs
c
      t0 = dist/gv0
      do ii = 1, nmodes
        if (ifill(ii) .ne. 0) then
          if (np .eq. 0) then
            tgv = modes(ii,27)
          elseif (np .eq. 1) then
            tgv = modes(ii,27) - modes(ii,28)
          elseif (np .eq. 2) then
            tgv = modes(ii,27) - modes(ii,28) - modes(ii,29)
          endif
          arg = s0*(tgv - t0)
          amp = exp(-0.5*arg**2)
          fill(ii) = amp*fill(ii)
        endif
      end do
      call color(nmodes)
      return
      end
