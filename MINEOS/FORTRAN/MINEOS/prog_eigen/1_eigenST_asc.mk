FBIN = /Users/russell/Lamont/GITHUB/ASWMS-ani/MINEOS/FORTRAN/bin
FC = gfortran
FFLAGS=-ffixed-line-length-none 
#-L/usr/local/include 
#FFLAGS2=-march=x86_64

all:  $(FBIN)/eigenST_asc

.f.o: 
	$(FC) $(FFLAGS) $(FFLAGS2) -c $*.f

#----------------------------------

$(FBIN)/eigenST_asc: eigenST_asc.f
	$(FC) $(FFLAGS) -o $(FBIN)/eigenST_asc eigenST_asc.f

clean: 
	rm -rf *.o $(FBIN)/eigenST_asc
