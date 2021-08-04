c
c
c
      subroutine interpol(n1, n2, x, y, m, b)
c
c     computes the coefficients for linear interpolation
c     y = mx + b
c
c     inputs:
c       n1:      lower bound for interpolation
c       n2:      upper bound for interpolation
c       x(n):    points at which the function is evaluated
c       y(n):    function to be interpolated
c     outputs:
c       m(n):    slopes of lines
c       b(n):    intercepts
c
      save
      parameter (n=5000)
      real x(n), y(n)
      real b(n), m(n)
c
      if ((n2-n1) .gt. n) then
        print*,'array limits exceeded in interpol',n2-n1,n
        stop
      endif
      do i = n1, n2-1
        dx = x(i+1) - x(i)
        dy = y(i+1) - y(i)
        if (dx .eq. 0.) then
          m(i) = 999.0
        else
          m(i) = dy/dx
        endif
        b(i) = y(i) - m(i)*x(i)
      end do
      return
      end
