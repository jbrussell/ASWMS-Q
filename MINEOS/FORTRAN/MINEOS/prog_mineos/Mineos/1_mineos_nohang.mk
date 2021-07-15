MYBIN = /Users/russell/Lamont/PROJ_LoveOvertones/FORTRAN/bin
FC = gfortran
#
PROG= mineos_nohang
SUBS= baylis.f bfs.f dermf.f derms.f detqn_nohang.f drspln.f dsplin.f eifout.f entry.f\
      fprop.f fprpmn.f fpsm.f fsbdry.f fsbm.f gauslv.f grav.f intgds.f match.f\
      model.f modout.f ortho.f remedy_nohang.f rkdot.f rotspl_nohang.f rprop.f rps.f sdepth.f\
      sfbdry.f sfbm.f sprop.f sprpmn.f spsm.f startl.f steps.f svd.f tprop.f\
      tps.f trknt.f whead.f wtable.f zknt.f
OBJS= mineos.o $(SUBS:.f=.o)

.f.o:
	$(FC) $(FFLAGS) -c $*.f

#----------------------------------------------------------------------------------

$(PROG): $(OBJS) 
	$(FC) $(FFLAGS) $(LFLAGS) -o $(MYBIN)/$@ $(OBJS)

# check object files for dependency on .h files
$(OBJS): parameter.h
	$(FC) $(FFLAGS) -c $*.f
	
	
clean: 
	rm -rf *.o

