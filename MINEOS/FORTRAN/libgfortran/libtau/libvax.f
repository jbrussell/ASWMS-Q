      subroutine abort(msg)
c
c $$$$$ calls VMS library routines $$$$$
c
c   Subroutine abort issues the VMS style fatal error message:
c
c      USRMSG-F-ABORT, <string>
c
c   accompanied by a standard FORTRAN traceback, where <string> is the
c   string in character variable msg.  Subroutine warn issues the VMS
c   style warning:
c
c      USRMSG-W-WARN, <string>
c
c   also accompanied by a FORTRAN traceback.  Calling abort terminates
c   the calling program.  However, warn returns control to the caller
c   after issuing the traceback.  Programmed on 15 June 1983 by
c   R. Buland.
c
      character*(*) msg
      external usrmsg_abort,usrmsg_warn
      call lib$stop(usrmsg_abort,%val(2),%val(len(msg)),%ref(msg))
      return
c   Subroutine warn.
      entry warn(msg)
      call lib$signal(usrmsg_warn,%val(2),%val(len(msg)),%ref(msg))
      return
      end
      subroutine dasign(lu,mode,ia,len)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine dasign opens (connects) logical unit lu to the disk file
c   named by the character string ia with mode mode.  If iabs(mode) = 1,
c   then open the file for reading.  If iabs(mode) = 2, then open the
c   file for writing.  If iabs(mode) = 3, then open a scratch file for
c   writing.  If mode > 0, then the file is formatted.  If mode < 0,
c   then the file is unformatted.  All files opened by dasign are
c   assumed to be direct access.  Programmed on 3 December 1979 by
c   R. Buland.
c
      character*(*) ia
      if(mode.ge.0) nf=1
      if(mode.lt.0) nf=2
      ns=iabs(mode)
      if(ns.le.0.or.ns.gt.4) ns=3
      go to (1,2),nf
 1    go to (11,12,13,14),ns
 11   open(lu,file=ia,status='old',form='formatted',carriagecontrol=
     1 'list',access='direct',organization='relative',recl=len,
     2 readonly,shared)
      return
 12   open(lu,file=ia,status='new',form='formatted',carriagecontrol=
     1 'list',access='direct',organization='relative',recl=len)
      return
 13   open(lu,file=ia,status='scratch',form='formatted',carriagecontrol=
     1 'list',access='direct',organization='relative',recl=len)
      return
 14   open(lu,file=ia,status='old',form='formatted',carriagecontrol=
     1 'list',access='direct',organization='relative',recl=len)
      return
 2    nn=(len+3)/4
      go to (21,22,23,24),ns
 21   open(lu,file=ia,status='old',form='unformatted',access='direct',
     1 organization='relative',recl=nn,readonly,shared)
      return
 22   open(lu,file=ia,status='new',form='unformatted',access='direct',
     1 organization='relative',recl=nn)
      return
 23   open(lu,file=ia,status='scratch',form='unformatted',access=
     1 'direct',organization='relative',recl=nn)
      return
 24   open(lu,file=ia,status='old',form='unformatted',access='direct',
     1 organization='relative',recl=nn)
      return
      end
      subroutine getarg(n,ia)
c
c $$$$$ calls getcl and nxarg $$$$$
c
c   Getarg returns the n th command line argument in character
c   variable ia.  For this to work you must run your job 'foreign'.
c   For example, suppose your executable file is named xxx.exe and
c   is in directory [zzz].  Make a global symbol as follows:
c
c   yyy:==$dba0:[zzz]xxx
c
c   Your program is executed by typing:
c
c   yyy arg1 arg2 ... argn
c
c   Arguments must be seperated from the global symbol and from each
c   other by one or more blanks.  Characters other than alphanumeric
c   ones usually give strange results.  Most special characters will
c   work if preceeded by an alphanumeric.  Alphabetic characters
c   are converted to upper case.  Note that ia will be blank if the
c   n th argument doesn't exist.  The argument will be truncated
c   if longer than ia and ia will be blank padded to the right if
c   longer than the argument.  Programmed on 28 January 1981 by
c   R. Buland.
c
      character*(*) ia
      character*100 ib
      common/gtarc1/nsw,narg,lb
      common/gtarc2/ib
c   Read the command line into ib if it hasn't already been done.
      if(nsw.lt.0) call getcl
c   Check for a valid argument.
      if(n.ge.0.and.n.le.narg) go to 1
c   Return a blank field if not.
      ia=' '
      return
 1    k=0
c   Skip to the n th argument.
      do 2 i=1,n
 2    j=nxarg(k)
      la=len(ia)
c   Copy the argument into ia.
      do 3 i=1,la
      ia(i:i)=ib(k:k)
      k=k+1
      if(k.gt.lb.or.ib(k:k).eq.' ') go to 4
 3    continue
      return
c   Blank pad if necessary.
 4    ia(i+1:la)=' '
      return
      end
      subroutine getcl
c
c $$$$$ calls only VMS library routines $$$$$
c
c   Getcl reads in the command line of a job run 'foreign' and
c   counts the number of arguments (delimited by blanks).
c   Programmed on 28 January 1981 by R. Buland.
c
      character*100 ib
      common/gtarc1/nsw,narg,lb
      common/gtarc2/ib
      data nsw/-1/
c   Flag that the command line has been read.
      nsw=1
c   Read the command line.
      call lib$get_foreign(ib,,lb)
c   Count the arguments.
      narg=0
      k=0
 1    if(nxarg(k).le.0) go to 13
      narg=narg+1
      go to 1
 13   return
      end
      function iargc(n)
c
c $$$$$ calls getcl $$$$$
c
c   Iargc returns the number of command line arguments present.
c   To use command line arguments see the comments in getarg.
c   Argument n of iargc is a dummy and is never used.
c   Programmed on 28 January 1981 by R. Buland.
c
      character*100 ib
      common/gtarc1/nsw,narg,lb
      common/gtarc2/ib
c   Read in the command line into ib if it hasn't already been done.
      if(nsw.lt.0) call getcl
      iargc=narg
      return
      end
      function nxarg(k)
c
c $$$$$ calls no other routine $$$$$
c
c   Nxarg assumes that the k th character of the command line
c   buffer, ib, is the first character of an argument or k=0.
c   It skips the remainder of the argument and the blanks
c   delimiting the next argument.  On normal return k will
c   point to the first character of the next argument and
c   nxarg will be set to 1.  If there is no next argument,
c   nxarg will be set to 0.  Programmed on 28 January 1981
c   by R. Buland.
c
      character*100 ib
      common/gtarc1/nsw,narg,lb
      common/gtarc2/ib
c   Exit if k already points beyond the end of the buffer.
      if(k.ge.lb) go to 13
c   Don't skip over the argument if we haven't found the first
c   one yet.
      if(k.le.0) go to 3
      k=k+1
c   Skip the remainder of the current argument.
 4    if(ib(k:k).eq.' ') go to 3
      k=k+1
      if(k-lb)4,4,13
 3    k=k+1
c   Skip blanks until the next argument or end of buffer.
 5    if(ib(k:k).ne.' ') go to 2
      k=k+1
      if(k-lb)5,5,13
c   Flag normal termination.
 2    nxarg=1
      return
c   Flag no next argument.
 13   nxarg=0
      return
      end
      subroutine tnoua(ia)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine tnoua writes the character string ia to the standard
c   output without the trailing newline (allowing user input on the
c   same line).  Programmed on 17 September 1980 by R. Buland.
c
      save
      character*(*) ia
      write(*,100)ia
 100  format(1x,a,$)
      return
      end
	.facility	usrmsg,1
	.severity	fatal
	abort	<!AD> /fao_count=1
	.severity	warning
	warn	<!AD> /fao_count=1
	.end
      subroutine vexit(ierr)
      if(ierr.eq.0) call exit
      if(ierr.eq.1) call abort('Terminated with errors.')
      call exit(ierr)
      end
