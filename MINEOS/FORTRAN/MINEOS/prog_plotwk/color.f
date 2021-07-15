c
c
c
      subroutine color(nmodes)
c
c     subroutine to scale color
c
      include 'parameter.h'
c
      real*4 fill(0:maxmodes)
      real*4 fmax, fmid, fmin, ftotal, fint
      integer*4 ifill(0:maxmodes)
      integer*4 cfill(0:maxmodes)
      logical lsummary, lcolor
c
      common /color2/ fmin, fmid, fmax
      common /color3/ lcolor, clmin, clmax
      common /fill/ lsummary, fill, ifill, cfill
c
c     use hardwired limits or loop over modes and determine limits of 'fill'
c
      if (lcolor) then
        fmin = clmin
        fmax = clmax
      else
        fmax = -999999999999999999.
        fmin =  999999999999999999.
        do ii = 1, nmodes
          if (ifill(ii) .ne. 0) then
            fmax = amax1(fmax,fill(ii))
            fmin = amin1(fmin,fill(ii))
          endif
        end do
      endif
c     print*,' color limits: ', fmin, fmax
      ftotal = fmax - fmin
      fmid = (fmax + fmin)/2.
      fint = abs(ftotal)/4.
      f1 = fmin + fint
      f2 = f1 + fint
      f3 = f2 + fint
      do ii = 1, nmodes
        if (ifill(ii) .ne. 0) then
          amp = fill(ii)
          if ((amp .ge. fmin) .and. (amp .lt. f1)) then 
            cfill(ii) = 4
          elseif ((amp .ge. f1) .and. (amp .lt. f2)) then 
            cfill(ii) = 5
          elseif ((amp .ge. f2) .and. (amp .lt. f3)) then 
            cfill(ii) = 7
          elseif ((amp .ge. f3) .and. (amp .le. fmax)) then 
            cfill(ii) = 8
          else
            cfill(ii) = 2
          endif
        else
          cfill(ii) = 2
        endif
      end do
      return
      end
