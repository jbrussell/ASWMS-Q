
      function lnblnk(string)

c**** find the # length of a text string, defined as the
c**** position of the last non-blank character. 

c**** this is a UTX/32 intrinsic
c**** may 88

      character*(*) string

c**** get # bytes in string

      nbytes = len(string)

      lnb = 0

      do i = nbytes,1,-1
         if(string(i:i) .ne. ' ') then
            lnb = i
            go to 200
         end if
      end do

200   continue
      lnblnk = lnb

      return
      end
