      subroutine kblnk(string,k)
c
c  returns the number of non-blank characters in string
c
      character*(*) string
      character*1 blank
      data blank/' '/
      k=0
      do 1 i=1,300
      if(string(i:i).eq.blank) go to 2
      k=i
    1 continue
    2 return
      end

