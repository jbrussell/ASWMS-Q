#
#  Compiler options.
#
FFLAGS= $(MYFFLAGS) 
LIBNAM= ../libcip.a
#
#  Library pathname.
#
#
#  Compile, archive and clean.
#
.f.a:
	gfortran $(FFLAGS) -c  $<
	ar rv $@ $*.o
	rm -f $*.o
#
#  List all the target objects.
#
$(LIBNAM): \
	$(LIBNAM)(ciplib.o) 
#
#  Set index.
#
$(LIBNAM): ; ranlib $(LIBNAM)
