      subroutine response(r_file)
c
c     subroutine to add instrument group delay
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop)
      real*4 mhz, rad, pi, tpi, dist
      real*4 omega(5000), gi(5000), mgi(5000), bgi(5000)
      real*4 interple
c
      integer*4 inx(0:maxl,0:maxn),inx2(0:maxmodes,2)
c
      character*256 r_file
c
      common /mode/ inx, modes, inx2
      common /nparam/ nparam, nextra
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /limits2/ wmin, wmax, pmin, pmax, gmin, gmax
c
      data tpi /6.2831853071796/
c
      mhz = 1000./tpi
      rad = 1.0/mhz  
      pi = tpi/2.
c
c     open and read file
c
      open(unit=3,file=r_file,form='unformatted',access='sequential')
c     &     status='readonly')
      read(3) ic
      do ii = 1, ic
        read(3) w, t1, t2, t3, dg
        omega(ii) = w*mhz
        gi(ii) = dg
      end do
      close(3)
c
c     check to find the index which is greater than the maximum frequency
c
      w = 0.0
      if = 0
      do while (w .lt. wmax)
        if = if + 1
        w = omega(if)
      end do
c
c     interpolate these values
c 
      call interpol(1,if,omega,gi,mgi,bgi)
c
c     now loop over the modes and add the instrument group delay to the total group delay
c     
      wl = 0.
      do ii = 1, nmodes
        w = modes(ii,1)
        wf = 0.0
        lf = 0
        do while (wf .lt. w)
          lf = lf + 1
          wf = omega(lf)
        end do
        lf = lf - 1
        if (lf .le. 0) then
          lf = 1
        endif
        dg = interple(lf,if,omega,w,wl,gi,mgi,bgi)
        modes(ii,29) = dg
      end do
      nextra = nextra + 1
c
c     all done
c
      return
      end
