#
#  Compiler options.
#
FFLAGS= $(MYFFLAGS) 
LIBNAM= $(MYLIB)/libcip.a
#
#  Library pathname.
#
#
#  Compile, archive and clean.
#
.f.a:
	f77 $(FFLAGS) -c  $<
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
