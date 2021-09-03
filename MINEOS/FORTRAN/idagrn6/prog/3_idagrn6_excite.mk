MYBIN = /Users/russell/Lamont/GITHUB/MINEOS_synthetics/FORTRAN/bin
FC = gfortran
MYLIB = /Users/russell/Lamont/GITHUB/MINEOS_synthetics/FORTRAN/libgfortran
#
PROG= idagrn6_excite
SUBS= spheroidal5.f toroidal5_excite.f spheroidal_exc5_jbr.f toroidal_exc5_jbr.f \
      curse4.f source_strains5_jbr.f cubic.f zfcns.f load_tables.f dwdd_h.f \
      readint_h.f load_crust20.f load_kern_c51.f avg_crust20.f avg_dbdb5.f \
      dwdd_c51.f dwdd_db.f

OBJS= $(PROG).o $(SUBS:.f=.o)

.f.o:
	$(FC) $(FFLAGS) -c $*.f

#----------------------------------------------------------------------------------

$(PROG): $(OBJS) 
	$(FC) $(FFLAGS) $(LFLAGS) -o $(MYBIN)/$@ $(OBJS) \
	$(MYLIB)/libcip.a \
	/opt/local/sac/lib/sacio.a \
	$(MYLIB)/libharm.a \
	$(MYLIB)/libutil.a

# check object files for dependency on .h files
$(OBJS): parameter5.h
	$(FC) $(FFLAGS) -c $*.f
