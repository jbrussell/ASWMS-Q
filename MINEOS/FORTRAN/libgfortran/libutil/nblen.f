      integer function nblen (string)
c
c     given a character string, nblen returns the length of the string
c     to the last non-blank character, presuming the string is left-
c     justified, i.e. if string = '   xs  j   ', nblen = 8.
c
c     called non-library routines: none
c     language: standard fortran 77
c                       
      integer ls, i
      character*(*) string, blank*1, null*1
      data blank /' '/
c
      null = char(0)
      nblen = 0
      ls = len(string)
      if (ls .eq. 0) return
      do 1 i = ls, 1, -1
         if (string(i:i) .ne. blank .and. string(i:i) .ne. null) go to 2
    1    continue
      return
    2 nblen = i
      return
      end
