	subroutine klen(string,k)
c
c     returns effect length of string - based on longue
c
	character*(*) string
      character*1 blank
      data blank/' '/
c
      k=0
      l=len(string)
      do i=l,1,-1
        if (string(i:i) .ne. blank) go to 2
      end do
	i=0
    2 continue
      k=i
	return
	end
