      subroutine cipget( prompt,token )
c
c   purpose:
c      This routine performs the routine parsing of all input
c      into tokens. It maintains it's own typein buffer and
c      simply returns a new token each time it is called
c   args:
c      prompt      string printed when new input is needed
c      token       token returned to calling program
c
c   versions and revisions:
c      for UNIX R.Goff Jan. 1982
c
      character*(*) prompt,token
      character*80 typbuf
      data typbuf/' '/
c
c   see if a buffer flush is desired
c
      if( prompt .eq. 'flush' ) then
        typbuf=' '
        ibp = 0
        return
      end if
c
c   see if no re-prompting is desired
c
      if( prompt .eq. 'noprompt' ) then
        call ciptok( typbuf,ibp,token )
        return
      end if
      go to 200
c
c   prompt for new input
c
100   call ciptyp( prompt,typbuf )
      ibp = 0
      if( typbuf .eq. ' ' ) then
        token = ' '
        return
      end if
c
c   get a token from input buffer
c
200   call ciptok( typbuf,ibp,token )
      if( token .eq. 'eoi' ) go to 100
c
c   return token to caller
c
      return
      end
