      subroutine cipscn(token,clist,ncm,igo )
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
c
c     changing this line to allow for matches when
c       a longer token is typed
c
c     if( token(1:lt) .eq. comand(1:lt) ) return
c
      if( token(1:lm) .eq. comand(1:lm) ) return
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
