      subroutine dw_c51(nn,ll,w,ista,dw)
c
c     program to calculate eigenfrequency perturbations for crust5.1
c     and/or dbdb5 crustal/bathymetry models.  Integrates crustal model 
c     perturbation against mode-specific frechet kernels to get dw for each 
c     mode.  
c
c     written 9/99 -- JBG
c
c     altered for 2 layer crust
c
C23456789112345678921234567893123456789412345678951234567896123456789712
c
      implicit real*4 (a-h,o-z)
c
      include 'parameter5.h'

      real*8 w,dw,scale
      real*4 dr(nlith),dalpha(nlith),dbeta(nlith),drho(nlith)
      real*4 intg(nlith)
c
c common block for kernel information

      real*4 crker(nlith,6,maxmodes),crkers(maxmodes),
     &       crkerb(maxmodes),crkerm(maxmodes), gvi(maxmodes)
      real*4 vs(nlith),vp(nlith),ro(nlith),rcr(nlith)
      real*4 rsol, rmoho, rbase
      integer*4 jcom,ifanis,kcrtot,ipnt(maxn,maxl)
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
     & crho2gc,ctopgc,cbasgc,cmohgc,cvp1la,
     & cvs1la,crho1la,cvp2la,cvs2la,crho2la,
     & ctopla,cbasla,cmohla
	 

      common/ahdr/ jcom,knots

      data pi/3.14159265350/
      data rn/6371000./
      data bigg /6.6723e-11/
      data rhobar/5515.0/
c
      scale = dble(1.0/(rn*sqrt(rn*pi*bigg)*rhobar))
c
c     set up integration parameters and constants.  Complex if loops
c     ensure that radii spanned by discontinuity perturbations do
c     not also have volumetric perturbations calculated
c

      do j = 1, kcrtot-1
c         print*, rcr(j), rbase
        dr(j) = rcr(j + 1) - rcr(j)
	if (rcr(j).eq.rbase) then
	  ibasec = j
	endif
      end do
      ibasec = ibasec - 1
c      print*, 'ibasec= ', ibasec

      rmohnew = rn + cmohgc(ista)
      rbasnew = rn + cbasgc(ista)
      rsolnew = rn + ctopgc(ista)
      dsol = rsolnew - rsol 
      dbas = rbasnew - rbase
      dmoh = rmohnew - rmoho 
	  

c     only do volumetric perturbations within ORIGINAL layering, not new perturbed layering.
c     zeros out perturbations if that layer is undefined in crust2.0 in that particular station path
c
      if (crho2gc(ista).lt.-900.) then
c         print*, 'here in a...'
         do i = 1,ibasec
           dalpha(i) = 0.
           dbeta(i) = 0.
           drho(i) = 0.
         end do
      else
c         print*, 'here in b...', ibasec
         do i = 1,ibasec
           dalpha(i) = cvp2gc(ista) - vp(i)
           dbeta(i) = cvs2gc(ista) - vs(i)
           drho(i) = crho2gc(ista) - ro(i)
         end do
      end if
      if (crho1gc(ista).lt.-900.) then
        do i = ibasec+1, kcrtot
           dalpha(i) = 0.
           dbeta(i) = 0.
           drho(i) = 0.
         end do
      else
         do i = ibasec+1, kcrtot
           dalpha(i) = cvp1gc(ista) - vp(i)
           dbeta(i) = cvs1gc(ista) - vs(i)
           drho(i) = crho1gc(ista) - ro(i)
         end do
      end if
c      print*, 'here in c...'
c      if(nn.eq.0.and.ll.eq.100) then
c         print*,'dsol,dbas,dmoh',dsol,dbas,dmoh
c        do i=1,kcrtot
c          write(*,'(4f12.3)')rcr(i),dalpha(i),dbeta(i),drho(i)
c        end do
c      end if

c      if(nn.eq.0.and.ll.eq.100)print*,'jcom,ifanis',jcom,ifanis
c     integrate the kernels for the perturbation to eigenfrequency

      if (jcom .ne. 2) then
        if(ifanis.eq.0) then
          do j = 1, kcrtot                         
            intg(j) = crker(j,1,ipnt(nn,ll))*drho(j) 
     &              + crker(j,2,ipnt(nn,ll))*dbeta(j)
     &              + crker(j,3,ipnt(nn,ll))*dalpha(j)
          end do
        else
          do j = 1, kcrtot
            intg(j) = crker(j,1,ipnt(nn,ll))*drho(j) 
     &              + crker(j,2,ipnt(nn,ll))*dbeta(j)
     &              + crker(j,3,ipnt(nn,ll))*dalpha(j)
     &              + crker(j,4,ipnt(nn,ll))*dbeta(j) 
     &              + crker(j,5,ipnt(nn,ll))*dalpha(j)
          end do
        end if
      else             
        if(ifanis.eq.0) then
          do j = 1, kcrtot              
            intg(j) = crker(j,1,ipnt(nn,ll))*drho(j) 
     &              + crker(j,2,ipnt(nn,ll))*dbeta(j)
          end do
        else
          do j = 1, kcrtot
c            if(nn.eq.0.and.ll.eq.100) then
c                print*,'rho betav betah kernels',
c     &          (crker(j,jjj,ipnt(nn,ll)),jjj=1,3)
c                print*,'dro,dbeta',rcr(j),drho(j),dbeta(j)
c            endif

            intg(j) = crker(j,1,ipnt(nn,ll))*drho(j) 
     &              + crker(j,2,ipnt(nn,ll))*dbeta(j)
     &              + crker(j,3,ipnt(nn,ll))*dbeta(j)
          end do
        endif
      endif
c
      t = 0.0d0
      do j = 1, kcrtot-1
        t = t + 0.5d0*(intg(j) + intg(j+1))*dr(j)
      end do
c                    
c      if(nn.eq.0.and.ll.eq.100) 
c     &  print*,'surf, moho kerns',crkers(ipnt(nn,ll)),
c     &  crkerm(ipnt(nn,ll))
      t = t - crkers(ipnt(nn,ll))*dsol - crkerm(ipnt(nn,ll))*dmoh
      t = t - crkerb(ipnt(nn,ll))*dbas
c                  
      dw = 0.5d0*w*dble(t)*scale*scale

c               if(nn.eq.0 .and. mod(ll,100).eq.0) then
c                 print*,nn,ll,w,dw
c               endif

      return
      end


      subroutine dd_c51(nn,ll,w,delt,ista,dd)
c
c     program to calculate distance perturbations for crust5.1
c     and/or dbdb5 crustal/bathymetry models.  Integrates crustal model 
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
      real*4 dr(nlith),dalpha(nlith),dbeta(nlith),drho(nlith)
      real*4 intg(nlith)
c
c common block for kernel information

      real*4 crker(nlith,6,maxmodes),crkers(maxmodes),
     &       crkerb(maxmodes),crkerm(maxmodes), gvi(maxmodes)
      real*4 vs(nlith),vp(nlith),ro(nlith),rcr(nlith)
      real*4 rsol, rmoho, rbase
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
     & crho2gc,ctopgc,cbasgc,cmohgc,cvp1la,
     & cvs1la,crho1la,cvp2la,cvs2la,crho2la,
     & ctopla,cbasla,cmohla

      data pi/3.14159265350/
      data rn/6371000./
      data bigg /6.6723e-11/
      data rhobar/5515.0/
c
      scale = dble(1.0/(rn*sqrt(rn*pi*bigg)*rhobar))
      rnd = dble(rn)
      gv = gvi(ipnt(nn,ll))*1000.
c
c     set up integration parameters and constants
c

      do j = 1, kcrtot-1
        dr(j) = rcr(j + 1) - rcr(j)
		if (rcr(i).eq.rbase) then
		  ibasec = i
		endif
      end do
      ibasec = ibasec - 1

      rmohnew = rn + cmohgc(ista)
	  rbasnew = rn + cbasgc(ista)
      rsolnew = rn + ctopgc(ista)
      dsol = rsolnew - rsol 
	  dbas = rbasnew - rbase
      dmoh = rmohnew - rmoho 

c
      if (crho2la(ista).lt.-900.) then
         do i = 1,ibasec
           dalpha(i) = 0.
           dbeta(i) = 0.
           drho(i) = 0.
         end do
      else
         do i = 1,ibasec
           dalpha(i) = cvp2la(ista) - vp(i)
           dbeta(i) = cvs2la(ista) - vs(i)
           drho(i) = crho2la(ista) - ro(i)
         end do
      end if
      if (crho1la(ista).lt.-900.) then
        do i = ibasec+1, kcrtot
           dalpha(i) = 0.
           dbeta(i) = 0.
           drho(i) = 0.
         end do
      else
         do i = ibasec+1, kcrtot
           dalpha(i) = cvp1la(ista) - vp(i)
           dbeta(i) = cvs1la(ista) - vs(i)
           drho(i) = crho1la(ista) - ro(i)
         end do
      end if

c     integrate the kernels for the perturbation to eigenfrequency

      wdb = dble(w)

      if (jcom .ne. 2) then
        if(ifanis.eq.0) then
          do j = 1, kcrtot                         
            intg(j) = crker(j,1,ipnt(nn,ll))*drho(j) 
     &              + crker(j,2,ipnt(nn,ll))*dbeta(j)
     &              + crker(j,3,ipnt(nn,ll))*dalpha(j)
          end do
        else
          do j = 1, kcrtot
            intg(j) = crker(j,1,ipnt(nn,ll))*drho(j) 
     &              + crker(j,2,ipnt(nn,ll))*dbeta(j)
     &              + crker(j,3,ipnt(nn,ll))*dalpha(j)
     &              + crker(j,4,ipnt(nn,ll))*dbeta(j) 
     &              + crker(j,5,ipnt(nn,ll))*dalpha(j)
          end do
        end if
      else             
        if(ifanis.eq.0) then
          do j = 1, kcrtot              
            intg(j) = crker(j,1,ipnt(nn,ll))*drho(j) 
     &              + crker(j,2,ipnt(nn,ll))*dbeta(j)
          end do
        else
          do j = 1, kcrtot
            intg(j) = crker(j,1,ipnt(nn,ll))*drho(j) 
     &              + crker(j,2,ipnt(nn,ll))*dbeta(j)
     &              + crker(j,3,ipnt(nn,ll))*dbeta(j)
          end do
        endif
      endif
c
      t = 0.0d0
      do j = 1, kcrtot-1
        t = t + 0.5d0*(intg(j) + intg(j+1))*dr(j)
      end do
c                    
      t = t - crkers(ipnt(nn,ll))*dsol - crkerm(ipnt(nn,ll))*dmoh
	  t = t - crkerb(ipnt(nn,ll))*dbas

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
