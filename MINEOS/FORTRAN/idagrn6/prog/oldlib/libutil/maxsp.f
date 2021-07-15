c
c
c
      subroutine maxsp(x,low,high,max,maxloc,min1,minl1,min2,minl2,idum)
c
c     maxsp locates the absolute maximum time sample in a time series
c
c     for idum = 1, it searches for the max positive amplitude
c     for idum = -1, it searchs for the max absolute amplitude
c                                  
      integer*4 low, high, maxloc, minl1, minl2, isg, idum
      real*4 x(*), amax, min1, min2, max
c
      amax = -99999999999.0
c
c     seek max positive amplitude
c
      if (idum .gt. 0) then
        do j = low, high
          if (x(j) .gt. amax) then
            max = x(j)
            maxloc = j
            amax = max
          end if
        end do
      else
c
c     seek max absolute amplitude
c
        do j = low, high
          if (abs(x(j)) .gt. amax) then
            max = x(j)
            maxloc = j
            amax = abs(max)
          end if
        end do
      end if
c
      if (max .le. 0.) then
        isg = -1
      else
        isg = 1
      endif
c
c     search forward and back for the associated minimums
c
      if (isg .ne. -1) then
        do j = maxloc - 1, low, -1
          if (x(j-1) .gt. x(j)) then
            min1 = x(j)
            minl1 = j
            go to 15
          endif
        end do
  15    continue
        do j = maxloc + 1, high, 1
          if (x(j+1) .gt. x(j)) then
            min2 = x(j)
            minl2 = j
            go to 20
          endif
        end do
  20    continue
      else
        do j = maxloc - 1, low, -1
          if (x(j-1) .lt. x(j)) then
            min1 = x(j)
            minl1 = j
            go to 25
          endif
        end do
  25    continue
        do j = maxloc + 1, high, 1
          if (x(j+1) .lt. x(j)) then
            min2 = x(j)
            minl2 = j
            go to 30
          endif
        end do
  30    continue
      endif
      return
      end
