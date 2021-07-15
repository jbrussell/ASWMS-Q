      subroutine cipnum( token,rnum,iflag )
c
c   purpose:
c      Interpret a token as a number. A real number is returned if
c      numeric interpretation is possible.
c
c   args:
c
c      token       character string to be interpreted
c      rnum        result if numeric - else = 0
c      iflag       =0 if o.k. else = ichar(first non-numeric)
c
c   versions and revisions:
c      for UNIVAC R.Goff Sept. 1981
c      for VAX/UNIX R.Goff Jan. 1981
c
      character*(*) token
      character*15 numstr
c
c   find length of token
c
      lt  =  index( token,' ' ) - 1
      if( lt .le. 0 ) lt = len( token )
c
c   scan token for bum chars
c
      do 100 i = 1,lt
      if( lge(token(i:i),'0') .and. lle(token(i:i),'9') ) go to 100
      if( token(i:i) .eq. 'e' ) go to 100
      if( token(i:i) .eq. '.' ) go to 100
      if( token(i:i) .eq. '+' ) go to 100
      if( token(i:i) .eq. '-' ) go to 100
c
c   return bum char in iflag
c
      rnum  =  0.
      iflag  =  ichar( token(i:i) )
      return
c
100   continue
c
c   move it to numstr right justified
c
      numstr  =  ' '
      numstr(16-lt:15) = token(1:lt)
c
c   use an internal read to decode it
c
      read( numstr,200 ) rnum
200   format(e15.0)
c
c   success
c
      iflag  =  0
      return
      end
