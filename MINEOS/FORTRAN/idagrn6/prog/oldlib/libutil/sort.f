c
c
c
      subroutine sort(n, array, nindex)
c
c     subroutine to sort array in increasing order
c     and store the index of order in nindex
c
c     taken from Numerical Recipes, Press et al., p. 233
c
      real*4 array(*)
      integer*4 nindex(*)
c
c     check for a single number
c
      if (n .eq. 1) then
        nindex(1) = 1
        return
      endif
c
      do j = 1, n
        nindex(j) = j
      end do
      l = (n/2) + 1
      ir = n
  10  continue
      if (l .gt. 1) then
        l = l - 1
        ndt = nindex(l)
        q = array(ndt)
      else
        ndt = nindex(ir)
        q = array(ndt)
        nindex(ir) = nindex(1)
        ir = ir - 1
        if (ir .eq. 1) then
          nindex(1) = ndt
          return
        endif
      endif
      i = l
      j = l + l
  20  if (j .le. ir) then
        if (j .lt. ir) then
          if(array(nindex(j)) .lt. array(nindex(j+1))) j = j + 1
        endif
        if (q .lt. array(nindex(j))) then
          nindex(i) = nindex(j)
          i = j
          j = j + j
        else
          j = ir + 1
        endif
        go to 20
      endif
      nindex(i) = ndt
      go to 10
      end
