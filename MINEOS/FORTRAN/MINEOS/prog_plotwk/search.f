c
c
c
      subroutine search(nm,lselect)
c
c     subroutine to search over mode table
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop)
      real*4 value(nprop,2)
      real*4 fill(0:maxmodes)
      integer*4 inx(0:maxl,0:maxn)
      integer*4 inx2(0:maxmodes,2)
      integer*4 ifill(0:maxmodes)
      integer*4 cfill(0:maxmodes)
      integer*4 index(nprop), ida
      character*80 cval(nprop)
      character*256 m_file
      logical land, lsummary, lselect
c
      common /mode/ inx, modes, inx2
      common /fill/ lsummary, fill, ifill, cfill
      common /nparam/ nparam, nextra
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /limits2/ wmin, wmax, pmin, pmax, gmin, gmax
      common /value/ cval   
      common /c_search/ icount, index, value
c
c
      data tpi / 6.2831853071796/
c
      CVAL(1) = '  1-W (mHz)                                     '
      CVAL(2) = '  2-Q                                           '
      CVAL(3) = '  3-PHASE VELOCITY (km/s)                       '
      CVAL(4) = '  4-GROUP VELOCITY (km/s)                       '
      CVAL(5) = '  5-COMPRESSIONAL ENERGY DENSITY - INNER CORE   '
      CVAL(6) = '  6-COMPRESSIONAL ENERGY DENSITY - OUTER CORE   '
      CVAL(7) = '  7-COMPRESSIONAL ENERGY DENSITY - LOWER MANTLE '
      CVAL(8) = '  8-COMPRESSIONAL ENERGY DENSITY - TRANSITION   '
      CVAL(9) = '  9-COMPRESSIONAL ENERGY DENSITY - UPPER MANTLE '
      CVAL(10)= ' 10-COMPRESSIONAL ENERGY DENSITY - MANTLE       '
      CVAL(11)= ' 11-COMPRESSIONAL ENERGY DENSITY - CRUST        '
      CVAL(12)= ' 12-SHEAR ENERGY DENSITY - INNER CORE           '
      CVAL(13)= ' 13-SHEAR ENERGY DENSITY - OUTER CORE           '
      CVAL(14)= ' 14-SHEAR ENERGY DENSITY - LOWER MANTLE         '
      CVAL(15)= ' 15-SHEAR ENERGY DENSITY - TRANSITION           '
      CVAL(16)= ' 16-SHEAR ENERGY DENSITY - UPPER MANTLE         '
      CVAL(17)= ' 17-SHEAR ENERGY DENSITY - MANTLE               '
      CVAL(18)= ' 18-SHEAR ENERGY DENSITY - CRUST                '
      CVAL(19)= ' 19-SHEAR/TOTAL - INNER CORE                    '
      CVAL(20)= ' 20-SHEAR/TOTAL - OUTER CORE                    '
      CVAL(21)= ' 21-SHEAR/TOTAL - LOWER MANTLE                  '
      CVAL(22)= ' 22-SHEAR/TOTAL - TRANSITION                    '
      CVAL(23)= ' 23-SHEAR/TOTAL - UPPER MANTLE                  '
      CVAL(24)= ' 24-SHEAR/TOTAL - MANTLE                        '
      CVAL(25)= ' 25-SHEAR/TOTAL - CRUST                         '
      CVAL(26)= ' 26-SHEAR/TOTAL - TOTAL                         '
      CVAL(27)= ' 27-ELASTIC GROUP DELAY (s)                     '
      CVAL(28)= ' 28-SOURCE GROUP DELAY  (s)                     '
      CVAL(29)= ' 29-INSTRUMENT GROUP DELAY (s)                  '

c

      ida = 0
      fmax = -99999999999.
      do kk = 1, nprop
        index(kk) = 0
        value(kk,1) = 0
        value(kk,2) = 0
      end do
c
c     print out existing table
c
      print*,' '
      do ii = 1, 4
        call klen(cval(ii),kk)
        print*, cval(ii)(1:kk)
      end do
      if (lsummary) then
        do ii = 5, nprop
          call klen(cval(ii),kk)
          print*, cval(ii)(1:kk)
        end do
      endif
      print*,' '
      print*,' Enter selection number with min and max values'
      print*,' Complete selection with one of the following and "0 0"'
      print*,' '
      print*,' 99 - exit loop with "and" search'
      print*,'100 - exit loop with "or" search'
      print*,'999 - quit'
      print*,' '
      icount = 0
35    continue
      read(*,*) ind, vmin, vmax
c
c     check for exit command
c
      if (ind .ge. 99) then 
        if (ind .eq. 999) then
          return
        elseif (ind .eq. 99) then
          land = .true.
        else
          land = .false.
        endif
        go to 36
      else
        icount = icount + 1
        index(icount) = ind
        value(ind,1) = vmin
        value(ind,2) = vmax
        go to 35
      endif
36    continue
      print*,' '
      do ii = 1, 4
        call klen(cval(ii),kk)
        print*, cval(ii)(1:kk), value(ii,1), value(ii,2)
      end do
      print*,' '
      if (lsummary) then
        do ii = 5, nprop
          call klen(cval(ii),kk)
          print*, cval(ii)(1:kk), value(ii,1), value(ii,2)
        end do
      endif
c      print*,' '
c      print*,'Amplitudes: '
c      print*,'  Enter 0 for no scale; id# for amplitude'
c      read(*,*) ida
      nm = 0
      if (.not.(lselect)) then
      if (land) then
        do in = 1, nmodes
          jj = inx2(in,1)
          ii = inx2(in,2)
          fill(in) = 0.
          ifill(in) = 0
          do kk = 1, icount
            ind = index(kk)
            val = modes(in,ind)
            if((val.lt.value(ind,1)).or.(val.gt.value(ind,2)))then
              go to 37
            endif
          end do
          nm = nm + 1
          if (ida .ne. 0) then
            fill(in) = modes(in,ida)
          else
            fill(in) = 1.0
          endif
          ifill(in) = 1
37        continue
        end do
      else
c
c     'or' search - must meet any criteria
c
        do in = 1, nmodes
          jj = inx2(in,1)
          ii = inx2(in,2)
          fill(in) = 0.
          ifill(in) = 0
          do kk = 1, icount
            ind = index(kk)
            val = modes(in,ind)
            if((val.ge.value(ind,1)).and.(val.le.value(ind,2)))then
              nm = nm + 1
              if (ida .ne. 0) then
                fill(in) = modes(in,ida)
              else
                fill(in) = 1.0
              endif
              ifill(in) = 1
              go to 38
            endif
38          continue
          end do
        end do
      endif
c
c     only search over the modes which have been selected already
c
      else
      if (land) then
        do in = 1, nmodes
          jj = inx2(in,1)
          ii = inx2(in,2)
          if (ifill(in) .ne. 0) then
            fill(in) = 0.
            ifill(in) = 0
            do kk = 1, icount
              ind = index(kk)
              val = modes(in,ind)
              if((val.lt.value(ind,1)).or.(val.gt.value(ind,2)))then
                go to 39
              endif
            end do
            nm = nm + 1
            if (ida .ne. 0) then
              fill(in) = modes(in,ida)
            else
              fill(in) = 1.0
            endif
            ifill(in) = 1
39          continue
          endif
        end do
      else
c
c     'or' search - must meet any criteria
c
        do in = 1, nmodes
          jj = inx2(in,1)
          ii = inx2(in,2)
          if (ifill(in) .ne. 0) then
            fill(in) = 0.
            ifill(in) = 0
            do kk = 1, icount
              ind = index(kk)
              val = modes(in,ind)
              if((val.ge.value(ind,1)).and.(val.le.value(ind,2)))then
                nm = nm + 1
                if (ida .ne. 0) then
                  fill(in) = modes(in,ida)
                else
                  fill(in) = 1.0
               endif
                ifill(in) = 1
                go to 40
              endif
40            continue
            end do
          endif
        end do
      endif
      endif
c
      call color(nmodes)
c
      lselect = .true.
c
c     all done
c
      return
      end
