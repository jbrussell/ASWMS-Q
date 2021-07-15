      subroutine get_unit(typbuf,ibp,sub,index,scale)
c
c     return scale factor depending on units
c     index = 1 for time units
c     index = 2 for frequency units
c
      character*(*) typbuf
      character*(*) sub
      integer*4 ibp, index
      real*4 scale
c
      include 'numerical.h'
c
c     time units - convert everything to seconds
c
      if (index .eq. 1) then
        call ciptok(typbuf,ibp,sub)
        if (sub .eq. 'eoi') then 
          scale = spmin
        elseif (sub .eq. 'x') then 
          scale = spmin
        elseif (sub(1:1) .eq. 's') then
          scale = 1.0
        elseif (sub(1:1) .eq. 'm') then
          scale = spmin
        elseif (sub(1:1) .eq. 'h') then
          scale = sphr  
        elseif (sub(1:1) .eq. 'd') then
          scale = spday
        else
          scale = spmin
        endif
c
c     frequency units - convert everything to mhz
c     -- note that this implies scale/value for "s"
c
      elseif (index .eq. 2) then
        call ciptok(typbuf,ibp,sub)
        if (sub .eq. 'eoi') then
          scale = 1.0
        elseif (sub .eq. 'x') then 
          scale = 1.0
        elseif (sub(1:1) .eq. 's') then
          scale = 1000.
        elseif (sub(1:1) .eq. 'h') then
          scale = 1000.
        elseif (sub(1:1) .eq. 'r') then
          scale = rad
        else
          scale = 1.0
        endif
      end if
      return
      end
