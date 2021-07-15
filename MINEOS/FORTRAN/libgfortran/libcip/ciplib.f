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
      subroutine cipped( nr,rnam,rval,ni,inam,ival )
c
c   purpose:
c      To facilitate the editing of arrays of numbers. This routine
c      allows the user to display a table of run-time
c      parameters and update them as needed. They are referred to by
c      names supplied by the calling program.
c
c   args:
c      nr          number of parameters in real table
c      rnam        array of character variables containing name of
c                  individual parameters in the same order as rval
c      rval        array of real valued parameters
c      ni          number of integer parameters
c      inam        array of character variables containing names of
c                  integer valued parameters
c      ival        array of integer valued parameters
c
      character*(*) rnam(*),inam(*)
      character*20 token
      character*80 typbuf
      dimension rval (*),ival (*)
c
c   get an editor command
c
100   call cipget('pedit:',typbuf)
      ibp=0
      call ciptok(typbuf,ibp,token)
105   call cipscn( token,'l l! u',2,igo)
      if( igo .gt. 0 ) go to 200
      if( igo .eq. 0 ) then
        print *,token,'   ???'
        call cipget('flush',token)
        go to 100
      end if
      go to ( 110,120,130,140,140,140 ),-igo
c
c   give them some help with the commands
c
110   print *
      print *,'parameter editing commands are:'
      print *
      print *,'l      list by name'
      print *,'l!     list entire table'
      print *,'u      update by name'
      print *
      go to 100
c
c   they asked to quit this nonesense
c
120   return
c
c   no menu for now
c
130   print *,'no menu available'
      go to 100
c
c   don't recognize yes no or blank
c
140   print *,token,'  ???'
      go to 100
c
c   they have specified a legal command
c
200   go to ( 300,400,500 ),igo
c
c   list by name
c
300   call cipget('name:',token)
305   do 310 ir = 1,nr
      if( token .eq. rnam(ir) ) go to 320
310   continue
      do 315 ii = 1,ni
      if( token .eq. inam(ii) ) go to 340
315   continue
      go to 105
c
c   found it
c
320   print 330,rnam(ir),'=',rval(ir)
330   format(a12,a,f10.3)
      go to 390
340   print 350,inam(ii),'=',ival(ii)
350   format(a12,a,i6)
390   call cipget('noprompt',token)
      if( token .ne. 'eoi' ) go to 305
      go to 100
c
c   list all parameters in table
c
400   print 410,(rnam(i),rval(i),i=1,nr)
410   format(4(a12,f10.3,' '))
      print 420,(inam(i),ival(i),i=1,ni)
420   format(4(a12,i6,'     '))
      go to 100
c
c   update by name
c
500   call ciptok(typbuf,ibp,token)
      if (token .eq. 'eoi') then
        call cipget('name:',token)
      endif
505   do 510 ir = 1,nr
      if( token .eq. rnam(ir) ) go to 520
510   continue
      do 515 ii = 1,ni
      if( token .eq. inam(ii) ) go to 540
515   continue
      go to 105
c
c   found it
c
520   call ciptok(typbuf,ibp,token)
      if (token .eq. 'eoi') then
         call cipget('real value:',token)
      endif
      call cipnum(token,rnum,iflag )
      if( iflag .ne. 0 ) go to 600
      rval(ir) = rnum
      go to 590
540   call ciptok(typbuf,ibp,token)
      if (token .eq. 'eoi') then
         call cipget('integer value:',token)
      endif
      call cipnum(token,rnum,iflag)
      if( iflag .ne. 0 ) go to 600
      ival(ii) = nint(rnum)
590   call cipget('noprompt',token)
      if( token .ne. 'eoi' ) go to 505
      go to 100
600   print *,token,' is not a legal numeric'
      call cipget('flush',token)
      go to 100
c
      end
      subroutine cipscn( token,clist,ncm,igo )
c
c   purpose:
c      To scan a list ( clist ) of possible commands for the
c      occurence of token. The list of possible commands is a
c      single character variable with the commands strung together
c      and delimited by blanks or commas just as a typed in line
c      would appear. A default list of commands is also scanned
c      if token is not found in the user supplied list. Shortened
c      typeins of at least ncm characters will cause a match.
c
c   args:
c      token      the typed in command ( see ciptok )
c      clist      the list of possible commands
c      ncm        integer indicating how many characters to match
c      igo        integer returned for use in a computed go to
c                 igo = 1   first command matched etc.
c                 igo = 0   no match made
c                 igo < 0   default command matched
c
c   versions and revisions:
c      for UNIVAC R.Goff Sept. 1981
c      for VAX/UNIX R.Goff Dec. 1981
c
      character*(*) token,clist
      character*20 comand
      character*80 dclist
      data dclist/'help,quit,menu,yes,no,stop,,'/
c
c   use ciptok to parse clist
c
      lt  =  index( token,' ' ) -1
      if( lt .le. 0 ) lt = len( token )
      icp  =  0
      igo  =  0
100   call ciptok( clist,icp,comand )
      if( comand .eq. 'eoi' ) go to 200
      igo  =  igo + 1
      lc  =  index( comand,' ' ) -1
      if( lc .le. 0 ) lc = len( comand )
      lm  =  min( lc,ncm )
      if( lt .lt. lm ) go to 100
      if( token(1:lt) .eq. comand(1:lt) ) return
      go to 100
c
c   search the default list now
c
200   igo  =  0
      icp  =  0
300   call ciptok( dclist,icp,comand )
      if( comand .eq. 'eoi' ) go to 400
      igo  =  igo - 1
      lc  =  index( comand,' ' ) -1
      if( lc .le. 0 ) lc = len( comand )
      lm  =  min( lc,ncm )
      if( lt .lt. lm ) go to 300
      if( token(1:lt) .eq. comand(1:lt) ) return
      go to 300
c
c   take care of 'eoi' and 'eot'
c
400   if( token .eq. 'eoi' ) then
        igo = -99
        return
      end if
      if( token .eq. 'eot' ) then
        igo  =  -999
        return
      end if
c
c   not in there either
c
      igo  =  0
      return
      end
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
      subroutine cipsub( skell,subs,string )
c
c   purpose:
c      Build a string by taking skeleton string and substituting
c      words into it where indicated by substitution metasequences
c      similar to those used in shell scripts. This routine is
c   args:
c      skell     skeleton string
c      subs      words to use in substitutions
c      string    output string which is built
c
c   example:
c      skell  =  '/u1/sss/minster/$1/$2.dat'
c      subs   =  'sssdata event1'
c      string =  'u1/sss/minster/sssdata/event1.dat'
c   author:
c      R. Goff  S-Cubed  Jan. 1982
c
      character*(*) skell,subs,string
      character*20 words(9)
      character*20 sword
      dimension iend(20),iwrd(20)
c
c   first tear the subs string into words
c
      ibp = 0
      do 100 iwd = 1,9
100   call ciptok( subs,ibp,words(iwd) )
c
c   find places where substitution is indicated
c
      nsub = 0
      do 200 i = 1,len(skell)
      if( skell(i:i) .ne. '$' ) go to 200
      iarg = ichar( skell(i+1:i+1) ) - ichar( '0' )
      if( iarg .lt. 0 .or. iarg .gt. 9 ) go to 200
      nsub = nsub + 1
      if( nsub .gt. 20 ) go to 300
      iend(nsub) = i - 1
      iwrd(nsub) = iarg
200   continue
c
c   now build the string
c
300   string = ' '
      ie = 0
      last = 1
      do 400 isub = 1,nsub
      is = ie + 1
      ie = is + iend(isub) - last
      string(is:ie) = skell(last:iend(isub))
      last = iend(isub) + 3
      sword = words(iwrd(isub))
      if( sword.eq.'eoi' ) go to 400
      lw = index(sword,' ') - 1
      if( lw .le. 0 ) lw = len(sword)
      is = ie + 1
      ie = is + lw - 1
      string(is:ie) = sword(1:lw)
400   continue
      is = ie + 1
      string(is:len(string) ) = skell(last:len(skell) )
      return
      end
      subroutine ciptok( typbuf,ibp,token )
c
c   purpose:
c      To take a typed in string of characters and break it
c      up into smaller strings called tokens. Tokens are
c      fields delimited by blanks or commas. Multiple or
c      redundant blanks are treated as one blank however
c      two successive commas will result in a blank token
c      being returned to the calling program whether are not
c      there are blanks in between them.
c
c   args:
c      typbuf      buffer containing the string to be parsed
c      ibp         buffer pointer ( 0 or -> delimiter )
c      token       token found ( or eoi for end of input )
c
c   versions and revisions:
c      for UNIVAC R.Goff Sept. 1981
c      for VAX/UNIX R.Goff Dec. 1981
c
      character*(*) typbuf,token
      character*1 olddel
c
c   see that buffer pointer is in bounds
c
      l  =  len( typbuf )
      if( ibp .lt. 0 ) then
        ibp  =  0
        token  =  ' '
        return
      end if
      if( ibp .ge. l ) then
        ibp  =  l
        token  =  'eoi'
        return
      end if
c
c   pick up old delimiter
c
      olddel  =  ' '
      if( ibp .gt. 0 ) olddel  =  typbuf( ibp:ibp )
c
c   skip white space
c
100   ibp  =  ibp + 1
      if( typbuf(ibp:ibp) .ne. ' ' ) go to 200
      if( ibp .lt. l ) go to 100
c
c   end of buffer - token = 'eoi'
c
150   token  =  'eoi'
      return
c
c   check for two commas
c
200   if( typbuf(ibp:ibp) .ne. ',' ) go to 300
      if( olddel .eq. ',' ) then
        token  =  ' '
        return
      end if
      olddel  =  ','
      if( ibp .ge. l ) go to 150
      go to 100
c
c   copy token to output
c
300   ist  =  ibp
400   ibp  =  ibp + 1
      if( typbuf(ibp:ibp) .eq. ' ' ) go to 500
      if( typbuf(ibp:ibp) .eq. ',' ) go to 500
      if( ibp .lt. l ) go to 400
      token  =  typbuf( ist:l )
      return
500   token  =  typbuf( ist:ibp-1 )
      return
      end
      subroutine ciptyp( insol,typbuf )
c
c   purpose:
c       To solicit and accept typeins from the keyboard. This routine
c       prints a prompt on the screen leaving the cursor just to the
c       left of it waiting for the user to type the answer. Since no
c       facility for suppresion of cr/lf is provided in the ANSI
c       standard for FORTRAN 77 this routine uses some machine dependant
c       extentions to the language.
c
c   args:
c      insol      character string used for the prompt
c      typbuf     character variable which will receive the typein
c
c   versions and revisions:
c      for UNIVAC R.Goff Sept 1981
c      for VAX/UNIX R.Goff Dec. 1981
c
      character*(*) insol,typbuf
c
c   print the prompt
c
      print 100,' ',insol
100   format(a)
c
c   retieve the typein
c
      read (*,200,end=300)typbuf
200   format( a )
      return
300   typbuf  =  'eot'
      end
