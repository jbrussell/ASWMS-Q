FFLAGS= $(MYFFLAGS)
LIBNAM= $(MYLIB)/libharm.a
#
.f.a:
	$(FC) $(FFLAGS) -c  $<
	$(AR) $@ $*.o
	rm -f $*.o
#
#  List all the target objects
#
$(LIBNAM): \
	$(LIBNAM)(angles.subs.f) \
	$(LIBNAM)(avgspher.subs.f)
		ranlib $(LIBNAM)

#
#   Check for dependency on .h files
#
$(LIBNAM)(avgspher.subs.o): harm_param.h

$(LIBNAM): ; $(RANLIB) $(LIBNAM)
