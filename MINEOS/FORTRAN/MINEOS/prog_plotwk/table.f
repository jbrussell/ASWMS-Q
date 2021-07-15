      subroutine table(m_file)
c
c     subroutine to open mode table file
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop)
      real*4 wmin, wmax, pmin, pmax, gmin, gmax
      integer*4 inx(0:maxl,0:maxn), inx2(0:maxmodes,2)
      integer*4 jcom, nmodes, nlow, nup, llow, lup
      character*256 m_file
c
      common /mode/ inx, modes, inx2
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /limits2/ wmin, wmax, pmin, pmax, gmin, gmax
      common /nparam/ nparam, nextra
c
      data tpi / 6.2831853071796/
c
c     open and read header records
c
      open(unit=3,file=m_file,form='unformatted',access='sequential')
      read(3) jcom, nmodes, wmin, wmax, llmin, llmax
      read(3) dum
      read(3) idum
      read(3) idum
      read(3) idum
c
c     intialize limits
c
      wmin = 999999.
      wmax = 0.
      pmin = 999999.
      pmax = 0.
      gmin = 999999.
      gmax = 0.
      nlow = 999999
      nup  = 0
      llow = 999999
      lup  = 0
c
c     now read mode table
c
      do ii = 1, nmodes
        read(3,end = 25) nn, ll, w, qq, alpha, phi, cv, gv
        inx(ll,nn) = ii
        inx2(ii,1) = ll
        inx2(ii,2) = nn
        modes(ii,1) = w*1000./tpi
        modes(ii,2) = qq
        modes(ii,3) = cv
        modes(ii,4) = gv
        do kk = 5, nparam
          modes(ii,kk) = 0.0
        end do
        wmin = amin1(wmin,modes(ii,1))
        wmax = amax1(wmax,modes(ii,1))
        pmin = amin1(pmin,modes(ii,3))
        pmax = amax1(pmax,modes(ii,3))
        gmin = amin1(gmin,modes(ii,4))
        gmax = amax1(gmax,modes(ii,4))
        nlow = min0(nlow,nn)
        nup  = max0(nup,nn)
        llow = min0(llow,ll)
        lup  = max0(lup,ll)
      end do
25    continue
      close(3)
      print*, ' jcom = ', jcom,' # modes = ',nmodes
      print*, ' wmin = ', wmin,' wmax = ', wmax
      print*, ' nmin = ', nlow,' nmax = ', nup
      print*, ' lmin = ', llow,' lmax = ', lup
c
c     check over table for lost modes
c
      print*, ' '      
      do ii = nlow, nup
        do jj = llow, lup
          if (inx(jj,ii) .eq. 0) then
            do kk = jj, lup
              if (inx(kk,ii) .ne. 0) then
                go to 30
              end if
            end do
          go to 35
30        continue
          print*,' mode not found: ', ii, jj
          end if
        end do
35      continue
      end do
c
      call color(nmodes)
c
      return
      end
