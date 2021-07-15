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
      print 100,insol
100   format( a )
c
c   retieve the typein
c
      read (*,200,end=300)typbuf
200   format( a )
      return
300   typbuf  =  'eot'
      end
