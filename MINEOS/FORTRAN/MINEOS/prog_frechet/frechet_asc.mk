#
# current libraries
#
# MYLIB = /quake/data1/gsdf.local/lib
# FFLAGS= -O 
MYLIB = ../lib
MYBIN = /Users/naccardo/Unix/MINEOS/bin
FFLAGS= $(MYFFLAGS)
MYBIN= ../bin

#
frechet_asc:    frechet_asc.o
	$(FC) frechet_asc.o -o $(MYBIN)/frechet_asc
