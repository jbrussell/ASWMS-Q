      subroutine cipstr( typbuf,ibp,string )
c
c   purpose:
c      To retrieve that portion of a typed-in line to the right
c      of the buffer pointer ( including imbedded blanks or
c      commas ) for use in titles or headings.
c
c   args:
c      typbuf      typein buffer returned by ciptyp
c      ibp         buffer pointer ( same as in ciptok )
c      string      character string returned
c
c   versions and revisions:
c      for UNIVAC R.Goff Sept. 1981
c      for VAX/UNIX R.Goff jan. 1982
c
      character*(*) typbuf,string
c
c   find length of typein buffer
c
      l=len (typbuf)
c
c   see if buffer pointer is in bounds
c
      if (ibp.lt.0) then
      ibp=0
      string=' '
      return
      end if
      if (ibp.ge.l)  then
      ibp=l
      string='eoi'
      return
      end if
c
c   strip off leading blanks
c
      ist=ibp
100   ist=ist+1
      if (typbuf(ist:ist).ne.' ')  go to 200
      if (ist.lt.l)  go to 100
c
c   ran out of buffer
c
      ibp=l
      string='eoi'
      return
200   string=typbuf(ist:l)
      ibp=l
      return
      end
