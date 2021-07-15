      subroutine sadd(dateo,rdiff)
c
c     subtracts rdiff seconds to dateo
c
      integer*4 dateo(6)
      real*4 rdiff
c
c     print*, dateo, rdiff
c
      iyear = dateo(1)
      iday = dateo(2)
      ihour = dateo(3)
      imin = dateo(4)
      isec = dateo(5) - int(rdiff)
      msec = dateo(6) - (rdiff - int(rdiff))*1000.
c
      ly = lpyr(iyear)
      if (ly .eq. 0) then
        jd = 365
      else
        jd = 366 
      endif
c
c     print*, iyear, iday, ihour, imin, isec, msec
c
      if (abs(msec) .ge. 1000) then
        irmsec = msec/1000
        msec = msec - irmsec * 1000
        isec = isec + irmsec
      endif
      if (msec .lt. 0) then
        isec = isec - 1
        msec = 1000 + msec
      endif
      dateo(6) = msec
c
      if (abs(isec) .ge. 60) then
        irisec = isec/60
        isec = isec - irisec * 60
        imin = imin + irisec
      endif
      if (isec .lt. 0) then
        imin = imin - 1
        isec = 60 + isec
      endif
      dateo(5) = isec
c
      if (abs(imin) .ge. 60) then
        irimin = imin/60
        imin = imin - irimin * 60
        ihour = ihour + irimin
      endif
      if (imin .lt. 0) then
        ihour = ihour - 1
        imin = 60 + imin
      endif
      dateo(4) = imin
c
      if (abs(ihour) .ge. 24) then
        irihour = ihour/24
        ihour = ihour - irihour * 24
        iday = iday + irihour
      endif
      if (ihour .lt. 0) then
        iday = iday - 1
        ihour = 24 + ihour
      endif
      dateo(3) = ihour
c
      if (abs(iday) .ge. jd) then
        iriday = iday/jd
        iday = iday - iriday * jd
        iyear = iyear + iriday
      endif
      if (iday .lt. 0) then
        iyear = iyear - 1
        iday = jd + iday
      endif
      dateo(2) = iday
      dateo(1) = iyear
c
c     print*, dateo
c
      return
      end
