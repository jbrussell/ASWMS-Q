#FBIN = ~/bin
FBIN = .
FC=gfortran
FFLAGS=-ffixed-line-length-none 
#-L/usr/local/include 
#FFLAGS2=-march=x86_64

all:  $(FBIN)/eigenSTscaled_asc

.f.o: 
	$(FC) $(FFLAGS) $(FFLAGS2) -c $*.f

#----------------------------------

$(FBIN)/frechet: frechet.f
	$(FC) $(FFLAGS) -o $(FBIN)/eigenST_asc eigenST_asc.f

clean: 
	rm -rf *o $(FBIN)/eigenST_asc