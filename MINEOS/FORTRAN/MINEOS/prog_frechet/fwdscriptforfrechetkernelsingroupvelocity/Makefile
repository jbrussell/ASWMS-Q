FFLAGS=-ffixed-line-length-none

.f.o: 
	g77 $(FFLAGS) -c $*.f

#----------------------------------

frechet: $(MYBIN)/frechet

$(MYBIN)/frechet: $(MYOBJ)/frechet.o
	g77 $(FFLAGS) -o $(MYBIN)/frechet \
	$(MYOBJ)/frechet.o 

$(MYOBJ)/frechet.o: frechet.f
	g77 $(FFLAGS) -c -o $(MYOBJ)/frechet.o frechet.f

#----------------------------------

frechet_gv: $(MYBIN)/frechet_gv

$(MYBIN)/frechet_gv: $(MYOBJ)/frechet_gv.o
	g77 $(FFLAGS) -o $(MYBIN)/frechet_gv \
	$(MYOBJ)/frechet_gv.o 

$(MYOBJ)/frechet_gv.o: frechet_gv.f
	g77 $(FFLAGS) -c -o $(MYOBJ)/frechet_gv.o frechet_gv.f

#----------------------------------

frechet_cv: $(MYBIN)/frechet_cv

$(MYBIN)/frechet_cv: $(MYOBJ)/frechet_cv.o
	g77 $(FFLAGS) -o $(MYBIN)/frechet_cv \
	$(MYOBJ)/frechet_cv.o 

$(MYOBJ)/frechet_cv.o: frechet_cv.f
	g77 $(FFLAGS) -c -o $(MYOBJ)/frechet_cv.o frechet_cv.f

#----------------------------------
#------------------------------------------------------

draw_frechet_gv: $(MYBIN)/draw_frechet_gv

$(MYBIN)/draw_frechet_gv: $(MYOBJ)/draw_frechet_gv.o
	g77 $(FFLAGS) -o $(MYBIN)/draw_frechet_gv \
	$(MYOBJ)/draw_frechet_gv.o 

$(MYOBJ)/draw_frechet_gv.o: draw_frechet_gv.f
	g77 $(FFLAGS) -c -o $(MYOBJ)/draw_frechet_gv.o draw_frechet_gv.f


#----------------------------------

frechet_writejunk: $(MYBIN)/frechet_writejunk

$(MYBIN)/frechet_writejunk: $(MYOBJ)/frechet_writejunk.o
	g77 $(FFLAGS) -o $(MYBIN)/frechet_writejunk \
	$(MYOBJ)/frechet_writejunk.o 

$(MYOBJ)/frechet_writejunk.o: frechet_writejunk.f
	g77 $(FFLAGS) -c -o $(MYOBJ)/frechet_writejunk.o frechet_writejunk.f
