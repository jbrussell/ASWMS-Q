FFLAGS= $(MYFLAGS)
LIBNAM= $(USRLIB)/libcip.a
.f.a:
	$(FC) $(FFLAGS) -c  $<
	ar rv $@ $*.o
	rm -f $*.o
#
#  List all the target objects
#
$(LIBNAM): \
	$(LIBNAM)(ciplib.o)
		ranlib $(LIBNAM)
