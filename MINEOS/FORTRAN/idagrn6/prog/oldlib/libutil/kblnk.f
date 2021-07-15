      subroutine kblnk(string,k)
c
c     returns the position of the first blank in a string
c
      character*(*) string
      character*1 blank
      data blank/' '/
c
      k=0
      l=len(string)
      do i=1,l
        if(string(i:i) .eq. blank) go to 2
        k=i
      end do
    2 return
      end
