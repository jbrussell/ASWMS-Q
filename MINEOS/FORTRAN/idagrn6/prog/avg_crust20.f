      subroutine avg_crust51(elat,elon,delta,azi,ldbdb5)
c bastardized version of program getCN2point (bottom half), from Gabi Laske
c given source and receiver location, calculates an average
c crustal model along the path and along the total great circle
c
c version 2 -- calculates separate averages and depths for sediments and basemen
c              in certain unusual cases, a layer may not exist for a particular path
c              it is HOPED that the boundary perturbation will still be correct.  In
c              this case, returned velocities are -999.
c JBG 10/04

c based on load_crust51, rather than getCN2point directl
c
c retains crust51 as the name for this subroutine (internally),
c the common blocks, and the variable names within the common.
c This allows the routine to be updated without changing the main
c program.  The file name records the CRUST X.X version number.
c
c input:  crust2.0 params passed via crust51 common block
c         elat,elon -- event location
c         dist,azi -- distance and azimuth from event to station
c         ldbdb5 -- option to read and average dbdb5 bathymetry,
c                   which has much higher resolution that crust51  
c
c outputs:  For both short-arc AND total great-circle
c             avg bathym/topo (negative if below sea level)
c             avg moho depth (negative if below sea level)
c             avg P, S velocity, avg density
c
c JBG 6/02
c        
c layer one and two flipped, after the read statement!
c layer 1: water
c layer 2: ice


      parameter(maxdel = 360)
      parameter(nla=90,nlo=180)
      real*4 delt(maxdel*24),gcloc(2,maxdel*24)
      real*4 elat,elon,delta,azi
      logical ldbdb5

c     common block of crustal parameters to input

      real*4 amapvp(8,nlo,nla),amaprho(8,nlo,nla),
     +          amapvs(8,nlo,nla),amapthi(7,nlo,nla),
     +          amapele(nlo,nla)

      common /crust51/ amapvp,amaprho,amapvs,amapthi,
     +                 amapele

c     common block of crustal params to return
c     new params -- vvp1, vvp2 are P velocity for seds and solid crust, etc.  
c     vbas is basement topography
      real*4 vtopsa,vmohsa,vbassa
      real*4 vvp1sa,vvs1sa,vrho1sa,vvp2sa,vvs2sa,vrho2sa
      real*4 vtopgc,vmohgc,vbasgc
      real*4 vvp1gc,vvs1gc,vrho1gc,vvp2gc,vvs2gc,vrho2gc
      common /avg51/ vvp1sa,vvs1sa,vrho1sa,vvp2sa,vvs2sa,vrho2sa,
     &               vtopsa,vbassa,vmohsa,vvp1gc,vvs1gc,vrho1gc,
     &                vvp2gc,vvs2gc,vrho2gc,vtopgc,vbasgc,vmohgc
  
      include 'numerical.h'


c     calculate points along entire great circle.  Assumes 1 degree sampling

      ndelgc = 360
      do i = 1,ndelgc
        delt(i) = float(i)*dkm
      end do
      call gcpath(elat,elon,azi,delt,ndelgc,gcloc)

c     loop over array of points -- internal to loop, once flat and
c     flon are set, it's mostly Laske code
c
c     new counter set to make sure that only layers with non-zero thicknesses enter the average
     
      vtopgc = 0.
      vmohgc = 0.
	  vbasgc = 0.
      vvp1gc = 0.
      vvs1gc = 0.
      vrho1gc = 0.
      vvp2gc = 0.
      vvs2gc = 0.
      vrho2gc = 0.
      ndelsa = int(delta)
      dx=360./nlo

      print*,'Averaging over ',delta,' sa, from ',elat,elon
      icount1 = 0
      icount2 = 0
      do ii = 1,ndelsa
        flat = gcloc(1,ii)
        flon = gcloc(2,ii)
        if(flon.gt.180.)flon=flon-360.
        cola=90.-flat
        ilat=int(cola/dx)+1
        ilon=int((flon+180.)/dx)+1
c
c       sediment layer

        vthi=0.
        vvp=0.
        vvs=0.
        vrho=0.
        do i=3,4
c          print*,ii,flat,flon,i,amapthi(i,ilon,ilat),
c     &           amapvp(i,ilon,ilat),
c     &           amapvs(i,ilon,ilat),amaprho(i,ilon,ilat)
          vthi=vthi+amapthi(i,ilon,ilat)
          vvp=vvp+amapthi(i,ilon,ilat)/amapvp(i,ilon,ilat)
          vvs=vvs+amapthi(i,ilon,ilat)/amapvs(i,ilon,ilat)
          vrho=vrho+amapthi(i,ilon,ilat)*amaprho(i,ilon,ilat)
        end do
c        print*,'elev: ',flat,flon,amapele(ilon,ilat)
        if (vthi .gt.0.) then
          vvp1gc= vvp1gc + vthi/vvp
          vvs1gc= vvs1gc + vthi/vvs
          vrho1gc= vrho1gc + vrho/vthi
          icount1=icount1+1
        endif
        vtopgc = vtopgc + amapele(ilon,ilat)
        vbasgc = vbasgc + (vthi - amapele(ilon,ilat)/1000.)
c
c       bedrock crust

        vthi1 = vthi
        vthi=0.
        vvp=0.
        vvs=0.
        vrho=0.
        do i=5,7
c          print*,ii,flat,flon,i,amapthi(i,ilon,ilat),
c     &           amapvp(i,ilon,ilat),
c     &           amapvs(i,ilon,ilat),amaprho(i,ilon,ilat)
          vthi=vthi+amapthi(i,ilon,ilat)
          vvp=vvp+amapthi(i,ilon,ilat)/amapvp(i,ilon,ilat)
          vvs=vvs+amapthi(i,ilon,ilat)/amapvs(i,ilon,ilat)
          vrho=vrho+amapthi(i,ilon,ilat)*amaprho(i,ilon,ilat)
        end do
c        print*,'elev: ',flat,flon,amapele(ilon,ilat)
        if (vthi.gt.0.) then
          vvp2gc= vvp2gc + vthi/vvp
          vvs2gc= vvs2gc + vthi/vvs
          vrho2gc= vrho2gc + vrho/vthi
          icount2=icount2+1
        endif
        vmohgc = vmohgc + (vthi+vthi1 - amapele(ilon,ilat)/1000.)
      end do

      if (icount1.gt.0) then
        vvp1sa = vvp1gc/icount1*1000.
        vvs1sa = vvs1gc/icount1*1000.
        vrho1sa = vrho1gc/icount1*1000.
      else
        vvp1sa = -999.
        vvs1sa = -999.
        vrho1sa = -999.
      endif
      if (icount2.gt.0) then
        vvp2sa = vvp2gc/icount2*1000.
        vvs2sa = vvs2gc/icount2*1000.
        vrho2sa = vrho2gc/icount2*1000.
      else
        vvp2sa = -999.
        vvs2sa = -999.
        vrho2sa = -999.
      endif
      vtopsa = vtopgc/ndelsa
      vbassa = -1.*vbasgc/ndelsa*1000.
      vmohsa = -1.*vmohgc/ndelsa*1000.

c     continue on around the great circle.  vvpgc, etc retain their values

      do ii = ndelsa+1, ndelgc
        flat = gcloc(1,ii)
        flon = gcloc(2,ii)
        if(flon.gt.180.)flon=flon-360.
        cola=90.-flat
        ilat=int(cola/dx)+1
        ilon=int((flon+180.)/dx)+1
c
c       sediment layer

        vthi=0.
        vvp=0.
        vvs=0.
        vrho=0.
        do i=3,4
          vthi=vthi+amapthi(i,ilon,ilat)
          vvp=vvp+amapthi(i,ilon,ilat)/amapvp(i,ilon,ilat)
          vvs=vvs+amapthi(i,ilon,ilat)/amapvs(i,ilon,ilat)
          vrho=vrho+amapthi(i,ilon,ilat)*amaprho(i,ilon,ilat)
        end do
c        print*,'elev: ',flat,flon,amapele(ilon,ilat)
        if (vthi.gt.0.) then
          vvp1gc= vvp1gc + vthi/vvp
          vvs1gc= vvs1gc + vthi/vvs
          vrho1gc= vrho1gc + vrho/vthi
          icount1=icount1+1
        endif
        vtopgc = vtopgc + amapele(ilon,ilat)
        vbasgc = vbasgc + (vthi - amapele(ilon,ilat)/1000.)
c
c       bedrock crust

        vthi1 = vthi
        vthi=0.
        vvp=0.
        vvs=0.
        vrho=0.
        do i=5,7
          vthi=vthi+amapthi(i,ilon,ilat)
          vvp=vvp+amapthi(i,ilon,ilat)/amapvp(i,ilon,ilat)
          vvs=vvs+amapthi(i,ilon,ilat)/amapvs(i,ilon,ilat)
          vrho=vrho+amapthi(i,ilon,ilat)*amaprho(i,ilon,ilat)
        end do
c        print*,'elev: ',flat,flon,amapele(ilon,ilat)
        if (vthi.gt.0.) then
          vvp2gc= vvp2gc + vthi/vvp
          vvs2gc= vvs2gc + vthi/vvs
          vrho2gc= vrho2gc + vrho/vthi
          icount2=icount2+1
        endif
        vmohgc = vmohgc + (vthi+vthi1 - amapele(ilon,ilat)/1000.)
      end do
	  
      if (icount1.gt.0) then
        vvp1gc = vvp1gc/icount1*1000.
        vvs1gc = vvs1gc/icount1*1000.
        vrho1gc = vrho1gc/icount1*1000.
      else
        vvp1gc = -999.
        vvs1gc = -999.
        vrho1gc = -999.
      endif
      if (icount2.gt.0) then
        vvp2gc = vvp2gc/icount2*1000.
        vvs2gc = vvs2gc/icount2*1000.
        vrho2gc = vrho2gc/icount2*1000.
      else
        vvp2gc = -999.
        vvs2gc = -999.
        vrho2gc = -999.
      endif
      vtopgc = vtopgc/ndelgc
      vbasgc = -1.*vbasgc/ndelgc*1000.
      vmohgc = -1.*vmohgc/ndelgc*1000.

      write(*,'(a)')'short arc: elev, z(base),z(moho):'
      write(*,'(3f12.1)') vtopsa,vbassa,vmohsa                              
      write(*,'(a)')'short arc: vp1,vs1,rho1,vp2,vs2,rho2:'
      write(*,'(6f12.1)')vvp1sa,vvs1sa,vrho1sa,vvp2sa,vvs2sa,vrho2sa
      write(*,'(a)')'great circle: elev, z(base),z(moho):'
      write(*,'(3f12.1)') vtopgc,vbasgc,vmohgc                              
      write(*,'(a)')'great circle: vp1,vs1,rho1,vp2,vs2,rho2:'
      write(*,'(6f12.1)')vvp1gc,vvs1gc,vrho1gc,vvp2gc,vvs2gc,vrho2gc
 794  format(5f12.1)


c     section to replace topography values with better averages from
c     etopo5/dbdb5 itself, rather than the 2x2 degree representation
c     This section is likely to be slow, may want to rewrite if using
c     extensively in the future. -- ONLY DONE FOR SHORT-ARC
c     nsample is the number of samples per degree -- dbdb5 is defined at 12/deg

      if(ldbdb5) then
        nsample = 24
        ndelgcb = ndelgc*nsample
        ndelsab = int(delta*nsample)
        do ii = 1,ndelgcb
          delt(ii) = float(ii)*dkm/float(nsample)
        end do

        call gcpath(elat,elon,azi,delt,ndelgcb,gcloc)

        dbgc = 0

        do ii = 1,ndelsab
         call get_bath(gcloc(1,ii),gcloc(2,ii),db)
         dbgc = dbgc + db
        enddo

        dbsa = dbgc / float(ndelsab)
        print*,'Crust5.1 SA elevation ',vtopsa,'replaced by ',dbsa
        vtopsa = dbsa

c too slow, and 2x2 is quite good for great circle -- skip this
c        do ii = ndelsab+1,ndelgcb
c           call get_bath(gcloc(1,ii),gcloc(2,ii),db)
c         dbgc = dbgc + db
c        enddo
c
c        dbgc = dbgc / float(ndelgcb)
c        print*,'Crust5.1 GC elevation ',vtopgc,'replaced by ',dbgc
c        vtopgc = dbgc

      endif

 99   continue
   
      return

      end
