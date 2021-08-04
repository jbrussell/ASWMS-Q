      subroutine summary(m_file)
c
c     subroutine to open summary file
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop)
      real*4 wmin, wmax, pmin, pmax, gmin, gmax
      real*4 temp(nprop)
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
      read(3) jcom, nmodes, nnmin, nnmax, llmin, llmax, wgrav
      read(3) 
      read(3) 
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
c     now read mode summary - this version skips the
c       gravitational potential
c
      do ii = 1, nmodes
        read(3,end = 15) nn,ll,(temp(kk), kk = 1, nparam)
        inx(ll,nn) = ii
        inx2(ii,1) = ll
        inx2(ii,2) = nn
        do kk = 1, 10
          modes(ii,kk) = temp(kk)
        end do
        if (nparam .gt. 13) then
c
c         reorganize summary table
c
c         ce - ic
          modes(ii,5)  = temp(5)
c         ce - oc
          modes(ii,6)  = temp(6)
c         ce - lm
          modes(ii,7)  = temp(17)
c         ce - tz
          modes(ii,8)  = temp(18)
c         ce - um
          modes(ii,9)  = temp(19)
c         ce - m
          modes(ii,10) = temp(7)
c         ce - c
          modes(ii,11) = temp(20)
c         se - ic
          modes(ii,12) = temp(8)
c         se - oc
          modes(ii,13) = temp(9)
c         se - lm
          modes(ii,14) = temp(21)
c         se - tz
          modes(ii,15) = temp(22)
c         se - um
          modes(ii,16) = temp(23)
c         se - m
          modes(ii,17) = temp(10)
c         se - c
          modes(ii,18) = temp(24)
c
c         add the ratios of shear/comp 
c
          p1 = temp(14)
          s1 = temp(15)
c         
c         se/ce - ic
c
          if ((modes(ii,12) .eq. 0.0) .and. (modes(ii,5) .eq. 0)) then
            modes(ii,19) = 0
          else
            modes(ii,19) = modes(ii,12)*s1/
     &                   (modes(ii,5)*p1+modes(ii,12)*s1)
          endif
c
c         se/ce - oc
c
          if ((modes(ii,13) .eq. 0.0) .and. (modes(ii,6) .eq. 0)) then
            modes(ii,20) = 0
          else
            modes(ii,20) = modes(ii,13)*s1/
     &                   (modes(ii,6)*p1+modes(ii,13)*s1)
          endif
c
c         se/ce - lm
c
          if ((modes(ii,14) .eq. 0.0) .and. (modes(ii,7) .eq. 0)) then
            modes(ii,21) = 0
          else
            modes(ii,21) = modes(ii,14)*s1/
     &                   (modes(ii,7)*p1+modes(ii,14)*s1)
          endif
c
c         se/ce - tz
c
          if ((modes(ii,15) .eq. 0.0) .and. (modes(ii,8) .eq. 0)) then
            modes(ii,22) = 0
          else
            modes(ii,22) = modes(ii,15)*s1/
     &                   (modes(ii,8)*p1+modes(ii,15)*s1)
          endif
c
c         se/ce - um
c
          if ((modes(ii,16) .eq. 0.0) .and. (modes(ii,9) .eq. 0)) then
            modes(ii,23) = 0
          else
            modes(ii,23) = modes(ii,16)*s1/
     &                   (modes(ii,9)*p1+modes(ii,16)*s1)
          endif
c
c         se/ce - m
c
          if ((modes(ii,17) .eq. 0.0) .and. (modes(ii,10) .eq. 0)) then
            modes(ii,24) = 0
          else
            modes(ii,24) = modes(ii,17)*s1/
     &                   (modes(ii,10)*p1+modes(ii,17)*s1)
          endif
c
c         se/ce - c
c
          if ((modes(ii,18) .eq. 0.0) .and. (modes(ii,11) .eq. 0)) then
            modes(ii,25) = 0
          else
            modes(ii,25) = modes(ii,17)*s1/
     &                   (modes(ii,10)*p1+modes(ii,17)*s1)
          endif
c
c         se/ce - total
c
          if ((s1 .eq. 0.0) .and. (p1 .eq. 0)) then
            modes(ii,26) = 0
          else
            modes(ii,26) = s1/(p1 + s1)
          endif
        endif         
        modes(ii,1) = modes(ii,1)*1000./tpi
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
15    continue
      close(3)
c
      print*, ' jcom = ', jcom,' # modes = ',nmodes
      print*, ' wmin = ', wmin,' wmax = ', wmax
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
