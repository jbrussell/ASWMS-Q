MYBIN = /Users/russell/Lamont/GITHUB/MINEOS_synthetics/FORTRAN/bin
FC = gfortran
FFLAGS=-ffixed-line-length-none
#
PROG= mineos_strip
SUBS=
OBJS= $(PROG).o $(SUBS:.f=.o)

.f.o:
	$(FC) $(FFLAGS) -c $*.f

#----------------------------------------------------------------------------------

$(PROG): $(OBJS) 
	$(FC) $(FFLAGS)  -o $(MYBIN)/$@ $(OBJS)

# check object files for dependency on .h files
#$(OBJS): parameter.h
#	$(FC) $(FFLAGS) -c $*.f
