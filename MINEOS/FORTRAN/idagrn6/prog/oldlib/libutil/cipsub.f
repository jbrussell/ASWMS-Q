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
      character*80 words(9)
      character*80 sword
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
