MYBIN = /Users/russell/Lamont/GITHUB/ASWMS-ani/MINEOS/FORTRAN/bin
FC = gfortran
#
PROG= mineos_nohang
SUBS= 
OBJS= mineos_bran_cadmod_jbrmod.o $(SUBS:.f=.o)

.f.o:
	$(FC) $(FFLAGS) -c $*.f

#----------------------------------------------------------------------------------

$(PROG): $(OBJS) 
	$(FC) $(FFLAGS) $(LFLAGS) -o $(MYBIN)/$@ $(OBJS)
	
clean: 
	rm -rf *.o
