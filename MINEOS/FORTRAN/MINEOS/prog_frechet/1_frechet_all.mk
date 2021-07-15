FBIN = /Users/russell/Lamont/PROJ_LoveOvertones/FORTRAN/bin
FC = gfortran
FFLAGS=-ffixed-line-length-none 
#-L/usr/local/include 
#FFLAGS2=-march=x86_64

all:  $(FBIN)/frechet  $(FBIN)/frechet_cv $(FBIN)/frechet_cv_ms $(FBIN)/frechet_cv_perc $(FBIN)/frechet_gv  $(FBIN)/draw_frechet_gv $(FBIN)/draw_frechet_gv_perc  $(FBIN)/frechet_psi $(FBIN)/frechet_cvG

.f.o: 
	$(FC) $(FFLAGS) $(FFLAGS2) -c $*.f

#----------------------------------

$(FBIN)/frechet: frechet.f
	$(FC) $(FFLAGS) -o $(FBIN)/frechet frechet.f
$(FBIN)/frechet_psi: frechet_psi.f
	$(FC) $(FFLAGS) -o $(FBIN)/frechet_psi frechet_psi.f
$(FBIN)/frechet_cv: frechet_cv.f
	$(FC) $(FFLAGS) -o $(FBIN)/frechet_cv frechet_cv.f
$(FBIN)/frechet_cv_ms: frechet_cv_ms.f
	$(FC) $(FFLAGS) -o $(FBIN)/frechet_cv_ms frechet_cv_ms.f
$(FBIN)/frechet_cv_perc: frechet_cv_perc.f
	$(FC) $(FFLAGS) -o $(FBIN)/frechet_cv_perc frechet_cv_perc.f

$(FBIN)/frechet_cvG: frechet_cvG.f
	$(FC) $(FFLAGS) -o $(FBIN)/frechet_cvG frechet_cvG.f
$(FBIN)/frechet_gv: frechet_gv.f
	$(FC) $(FFLAGS) -o $(FBIN)/frechet_gv frechet_gv.f

$(FBIN)/draw_frechet_gv: draw_frechet_gv.f
	$(FC) $(FFLAGS)  -o $(FBIN)/draw_frechet_gv draw_frechet_gv.f
$(FBIN)/draw_frechet_gv_perc: draw_frechet_gv_perc.f
	$(FC) $(FFLAGS)  -o $(FBIN)/draw_frechet_gv_perc draw_frechet_gv_perc.f

clean: 
	rm -rf *.o $(FBIN)/frechet $(FBIN)/frechet_cv $(FBIN)/frechet_cv_ms $(FBIN)/frechet_cv_perc $(FBIN)/frechet_gv $(FBIN)/draw_frechet_gv $(FBIN)/draw_frechet_gv_perc $(FBIN)/frechet_psi $(FBIN)/frechet_cvG
