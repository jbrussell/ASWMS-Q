      subroutine seek(a,i1,i2,value,index)
c
      include 'parameter.h'
c
      real*4 a(nknot)
c
      dv = 999999999
c
      do ii = i1, i2
        if (dv .gt. abs(a(ii) - value)) then
          dv = abs(a(ii) - value)
          index = ii
        endif
      end do                  
      return
      end
