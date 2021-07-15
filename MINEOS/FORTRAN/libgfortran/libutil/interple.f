c
c
c
      real function interple(n1, n2, x, dx, xlast, y, m, b)
c
c     given the coefficients for linear interpolation
c     this routine calculates y for an input x
c
c     inputs:
c       n1:      lower bound
c       n2:      upper bound
c       x(n):    array of x-values
c       dx:      point a which the function is to be evaluated
c       y(n):    function to be interpolated
c       m(n-1):  slopes
c       b(n-1):  intercepts
c     returned
c       y:       interpolated value
c
      real x(*), dx, y(*)
      real b(*), m(*), xlast
c
      do i = n1, n2
        if (dx .eq. x(i)) then
          if (dx .eq. x(i+1)) then
            if (xlast .eq. 0.) then
              interple = y(i+1)
              return
            elseif (xlast .lt. x(i)) then
              interple = y(i)
              return
            else
              interple = y(i+1)
              return
            endif
          else
            interple = y(i)
            return
          endif
        elseif ((dx .gt. x(i)) .and. (dx .lt. x(i+1))) then
          if (m(i) .ge. 999.0) then
            if (xlast .lt. dx) then
              interple = y(i)
            else
              interple = y(i+1)
            endif
          else
            interple = m(i)*dx + b(i)
          endif
          return
        endif
      end do
20    continue
c
c     outside array bounds - extrapolate
c
      if (dx .lt. x(n1)) then
        interple = m(n1)*dx + b(n1)
      elseif (dx .gt. x(n2)) then
        interple = m(n2)*dx + b(n2)
      else
        print*,' error in interpolation'
      endif
      return
      end
