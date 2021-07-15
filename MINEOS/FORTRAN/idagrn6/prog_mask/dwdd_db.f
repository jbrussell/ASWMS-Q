      subroutine dw_dbd(nn,ll,w,ista,dw)
c
c     program to calculate eigenfrequency perturbations for 
c     dbdb5 bathymetry models.  Integrates crustal model 
c     perturbation against mode-specific frechet kernels to get dw for each 
c     mode.  
c
c     written 9/99 -- JBG
c
C23456789112345678921234567893123456789412345678951234567896123456789712
c
      implicit real*4 (a-h,o-z)
c
      include 'parameter5.h'

      real*8 w,dw,scale
c
c common block for kernel information

      real*4 crker(nlith,6,maxmodes),crkers(maxmodes),
     &       crkerb(maxmodes),crkerm(maxmodes), gvi(maxmodes)
      real*4 vs(nlith),vp(nlith),ro(nlith),rcr(nlith)
      real*4 rsol, rbase, rmoho
      integer*4 ifanis,kcrtot,ipnt(maxn,maxl)
      common /kern51/ ifanis,crker,crkers,crkerb,crkerm,gvi,ipnt
      common /omod51/ rsol,rbase,rmoho,vs,vp,ro,rcr,kcrtot

c common block for station by station crustal averages

      real*4 cvp1gc(maxstat),cvs1gc(maxstat),crho1gc(maxstat)
      real*4 cvp2gc(maxstat),cvs2gc(maxstat),crho2gc(maxstat)
      real*4 ctopgc(maxstat),cbasgc(maxstat),cmohgc(maxstat)
      real*4 cvp1la(maxstat),cvs1la(maxstat),crho1la(maxstat)
      real*4 cvp2la(maxstat),cvs2la(maxstat),crho2la(maxstat)
      real*4 ctopla(maxstat),cbasla(maxstat),cmohla(maxstat)
      common /savg51/ cvp1gc,cvs1gc,crho1gc,cvp2gc,cvs2gc,
     &  crho2gc,ctopgc,cbasgc,cmohgc,cvp1la,
     &  cvs1la,crho1la,cvp2la,cvs2la,crho2la,
     &  ctopla,cbasla,cmohla


      data pi/3.14159265350/
      data rn/6371000./
      data bigg /6.6723e-11/
      data rhobar/5515.0/
c
      scale = dble(1.0/(rn*sqrt(rn*pi*bigg)*rhobar))
c
c     set up integration parameters and constants
c


      dsol = (rn + ctopgc(ista)) - rsol 
      dmoh = (rn + cmohgc(ista)) - rmoho 

c      if(nn.eq.0.and.mod(ll,100).eq.0)print*,'top,moh,rs,rm,ds,dm',
c     &         rn,ctopgc(ista),cmohgc(ista),rsol,rmoho,dsol,dmoh
c
      t = 0.0d0
c                    
      t = t - crkers(ipnt(nn,ll))*dsol - crkerm(ipnt(nn,ll))*dmoh
c                  
      dw = 0.5d0*w*dble(t)*scale*scale

      return
      end


      subroutine dd_dbd(nn,ll,w,delt,ista,dd)
c
c     program to calculate distance perturbations for
c     dbdb5 bathymetry models.  Integrates crustal model 
c     perturbation against mode-specific frechet kernels to get dd for each 
c     mode.  
c
c     written 9/99 -- JBG
c
C23456789112345678921234567893123456789412345678951234567896123456789712
c
      implicit real*4 (a-h,o-z)
c
      include 'parameter5.h'

      real*8 w,dw,scale,fact,delt,dd
c
c common block for kernel information

c common block for kernel information

      real*4 crker(nlith,6,maxmodes),crkers(maxmodes),
     &       crkerb(maxmodes),crkerm(maxmodes), gvi(maxmodes)
      real*4 vs(nlith),vp(nlith),ro(nlith),rcr(nlith)
      real*4 rsol, rbase, rmoho
      integer*4 ifanis,kcrtot,ipnt(maxn,maxl)
      common /kern51/ ifanis,crker,crkers,crkerb,crkerm,gvi,ipnt
      common /omod51/ rsol,rbase,rmoho,vs,vp,ro,rcr,kcrtot

c common block for station by station crustal averages

      real*4 cvp1gc(maxstat),cvs1gc(maxstat),crho1gc(maxstat)
      real*4 cvp2gc(maxstat),cvs2gc(maxstat),crho2gc(maxstat)
      real*4 ctopgc(maxstat),cbasgc(maxstat),cmohgc(maxstat)
      real*4 cvp1la(maxstat),cvs1la(maxstat),crho1la(maxstat)
      real*4 cvp2la(maxstat),cvs2la(maxstat),crho2la(maxstat)
      real*4 ctopla(maxstat),cbasla(maxstat),cmohla(maxstat)
      common /savg51/ cvp1gc,cvs1gc,crho1gc,cvp2gc,cvs2gc,
     &  crho2gc,ctopgc,cbasgc,cmohgc,cvp1la,
     &  cvs1la,crho1la,cvp2la,cvs2la,crho2la,
     &  ctopla,cbasla,cmohla


      data pi/3.14159265350/
      data rn/6371000./
      data bigg /6.6723e-11/
      data rhobar/5515.0/
c
      scale = dble(1.0/(rn*sqrt(rn*pi*bigg)*rhobar))
      rnd = dble(rn)
      gv = gvi(ipnt(nn,ll))*1000.

      dsol = (rn + ctopla(ista)) - rsol 
      dmoh = (rn + cmohla(ista)) - rmoho 

c     integrate the kernels for the perturbation to eigenfrequency

      wdb = dble(w)

      t = crkers(ipnt(nn,ll))*dsol - crkerm(ipnt(nn,ll))*dmoh
c                  
      fact = rnd/((dble(ll)+0.5d0)*dble(gv))

c     assume that dd is the same as dw*fact*delt (taken from dwdd_h analogy)

      dw = 0.5d0*w*dble(t)*scale*scale
      dd = dw*fact*delt
cc...ellipticity correction
c      fact2 = -3.d0*fact*dsin(delt)*dcos(2.d0*azim)*dsin(thet)**2
c      dd = dd + wa_t(ip)*fact2

      return
      end
