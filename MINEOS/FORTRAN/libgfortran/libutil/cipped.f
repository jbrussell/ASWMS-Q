      subroutine cipped(nr,rnam,rval,ni,inam,ival)
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
c
c changed 5/13/92 lsg
c
c     character*80 typbuf
      character*256 typbuf
c
      real*4 rval(*)
      integer*4 ival(*)
c
c   get an editor command
c
100   call cipget('pedit:',typbuf)
      ibp=0
      call ciptok(typbuf,ibp,token)
105   call cipscn( token,'la l u',2,igo)
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
      print *,'la     list all'
      print *,'u      update by name'
      print *,'quit   quit editor'
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
200   go to ( 400,300,500 ),igo
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
400   continue
      print*,' ' 
      print 420,(inam(i),ival(i),i=1,ni)
      print*,' ' 
      print 410,(rnam(i),rval(i),i=1,nr)
      print*,' ' 
410   format(4(a12,f10.3,' '))
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
