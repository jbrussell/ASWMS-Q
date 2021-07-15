c
c
c
      integer function ival(ifunc)
c
c     returns 1 ifunc .gt. 0
c             0 otherwise
c
      integer ifunc
c
      if (ifunc .gt. 0) then
        ival = 1
      else
        ival = 0
      endif
c
      return
      end     
