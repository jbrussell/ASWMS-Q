MYBIN = /Users/russell/Lamont/GITHUB/ASWMS-ani/MINEOS/FORTRAN/bin
FC = gfortran
FFLAGS=-ffixed-line-length-none
#
PROG= mineos_table
SUBS= kblnk.f
OBJS= $(PROG).o $(SUBS:.f=.o)

.f.o:
	$(FC) $(FFLAGS) -c $*.f

#----------------------------------------------------------------------------------

$(PROG): $(OBJS) 
	$(FC) $(FFLAGS) -o $(MYBIN)/$@ $(OBJS)

# check object files for dependency on .h files
$(OBJS): parameter.h
	$(FC) $(FFLAGS) -c $*.f
