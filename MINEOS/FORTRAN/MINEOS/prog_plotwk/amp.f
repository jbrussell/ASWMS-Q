c
c
c
      subroutine amp(ida,lselect,title)
c
c     subroutine to search over mode table
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop)
      real*4 value(nprop,2)
      real*4 fill(0:maxmodes)
      integer*4 inx(0:maxl,0:maxn),ifill(0:maxmodes),cfill(0:maxmodes)
      integer*4 inx2(0:maxmodes,2)
      integer*4 index(nprop), ida
      character*80 cval(nprop)
      character*256 title
      logical lsummary, lselect
c
      common /mode/ inx, modes, inx2
      common /fill/ lsummary, fill, ifill, cfill
      common /nparam/ nparam, nextra
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /value/ cval
c
c     print out parameter
c
      if ((ida .le. 0) .or. (ida .gt. 29)) then
        print*,' error selecting parameter'
        print*, ida
        ida = 0
        return
      else
c        print*, cval(ida)
      endif
      if (lselect) then
        do in = 1, nmodes
          if (ifill(in) .ne. 0) then
            fill(in) = modes(in,ida)
          endif
        end do
      else
        do in = 1, nmodes
          fill(in) = modes(in,ida)
          ifill(in) = 1
       end do
      endif
      title = cval(ida)(5:80)
      call color(nmodes)
      return
      end
