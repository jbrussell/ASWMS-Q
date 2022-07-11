MYBIN = /Users/russell/Lamont/GITHUB/MINEOS_synthetics/FORTRAN/bin
FC = gfortran
MYLIB = /Users/russell/Lamont/GITHUB/MINEOS_synthetics/FORTRAN/libgfortran
#
PROG= plot_wk
SUBS= amp.f branch_sort.f class.f color.f cvtaper.f excite.f fix_class_c.f \
      fix_class_k.f fix_class_p.f fix_class_r.f fix_class_v.f interple.f \
      interpol.f response.f search.o seek.f summary.f table.f wind.f
OBJS= $(PROG).o $(SUBS:.f=.o)

.f.o:
	$(FC) $(FFLAGS) -c $*.f

#----------------------------------------------------------------------------------


$(PROG): $(OBJS) 
	$(FC) $(FFLAGS) $(LFLAGS) -o $(MYBIN)/$@ $(OBJS) \
        $(MYLIB)/libcip.a \
        $(MYLIB)/libutil.a $(MYLIB)/libtau.a 

# clean up huge .o files
	rm plot_wk.o branch_sort.o
# check object files for dependency on .h files
$(OBJS): parameter.h
	$(FC) $(FFLAGS) -c $*.f
	
clean: 
	rm -rf *.o
