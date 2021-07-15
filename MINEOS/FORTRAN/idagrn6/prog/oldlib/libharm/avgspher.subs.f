**********************************************************************
C
c                       MAIN SUBROUTINE
C
C
C     AVERAGE OVER MINOR AND MAJOR PATHS, AND OF GREAT CIRCLE OF 
C                  COMPLEX SPHERICAL HARMONICS                         
c
c----------------------------------------------------------------------
c
c      subroutine spher_harm(ismin,ismax,lat,lon,spher)
c      subroutine avg_shm(type,rlat1,rlon1,rlat2,rlon2,ismax,avg_spher)
c      subroutine avg_great_pole(polelat,polelon,ismax,avg_spher)
c
c      subroutine minor_arc(elat,elon,slat,slon,polelat,polelon,
c     &                   npts,alpha,olat,olon,ds)
c      subroutine major_arc(elat,elon,slat,slon,polelat,polelon,
c     &                   npts,alpha,olat,olon,ds)
c
c      subroutine midpt(elat,elon,slat,slon,dlat,dlon)
c
c      subroutine rot_matrices(lat,lon,a,b)
c      subroutine rotate(c,olat,olon,nlat,nlon)
c      subroutine pole(elat,elon,slat,slon,polelat,polelon)
c
c      double precision FUNCTION DFACTORL(N)
c      double precision FUNCTION DGAMMALN(XX)
c      double precision FUNCTION DPLGNDRE(L,M,X)
c
c**********************************************************************

      subroutine spher_harm(ismin,ismax,lat,lon,spher)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     compute complex spherical harmonics
c----------------------------------------------------------------------

      include 'harm_param.h'

      integer*4 ismin,ismax
      real*8 dplgndre,dfactorl
      real*8 lat,lon,xxx,yyy,alon,aa,bb,radcon
      complex*16 coeff(0:nord,-nord:nord),spher(0:nord,-nord:nord)     

      pi = 4.d0*datan(1.d0)
      radcon = pi/180.d0

      do iss = ismin,ismax
        do itt = 0,iss
          coeff(iss,itt) = dsqrt(dble(2*iss+1)*DFACTORL(iss-itt)
     &                  /(4.d0*pi*DFACTORL(iss+itt)))
        enddo
      enddo

      xxx  = (90.d0 - lat)*radcon
      xxx  = dcos(xxx)
      yyy  = lon*radcon
      do iss = ismin,ismax
        do itt = 0,iss
          alon = yyy*dble(itt)
          aa = coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dcos(alon)
          bb = coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dsin(alon)
          spher(iss,itt)  = dcmplx(aa, bb)
          spher(iss,-itt) = dconjg(spher(iss,itt))*(-1.d0)**(itt)
        enddo
      enddo
      return
      end

      subroutine avg_shm(type,rlat1,rlon1,rlat2,rlon2,ismax,avg_spher)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c input :
c     type ---> 'great ','minor ','major '
c     rlat1,rlon1 --> event location
c     rlat2,rlon2 --> station location
c     ismax ---> order of harmonic expansion
c output :
c     avg_spher --> average over spherical harmonics
c
c     tested: 12/5/91 PI
c
c     sampling for numerical averages : ---> 1 degrees
c     precision expected from comparison between analytical
c     and numerical integration for great circle averages
c                    up to degree 12  : ---> 1e-6
c
c     warning: the function plgndr(l,m,x) has a term (-1)**m
c     already included (see Numerical Recipes 1986).This is 
c     different from Edmonds (1960) convention, but eventually 
c     the result is the same
c
c     structure :                                     
c
c     - computes pole to great circle connecting event and station
c       event-station-pole are a right coordinate system
c     - compute analytical average over great circle path
c     - compute minor arc average
c              - compute equidistant points on minor arc in subroutine
c              - do average in main program
c     - compute major arc average
c              - compute equidistant points on major arc in subroutine
c              - do average in main program
c
c     no common blocks, double precision
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      include 'harm_param.h'

      integer*4 ismax
      
      real*8 dplgndre,dfactorl
      real*8 aa,bb,xxx,yyy
      real*8 bigthe,bigphi
      real*8 pi,alon,fact,alpha

      real*8 rlat1,rlon1,rlat2,rlon2
      real*8 elat,elon,slat,slon,polelat,polelon
      real*8 min_lat(nsamp),min_lon(nsamp),min_ds(nsamp)
      real*8 maj_lat(nsamp),maj_lon(nsamp),maj_ds(nsamp)
      real*8 coeff(0:nord,-nord:nord)

      complex*16 avg_spher(0:nord,-nord:nord)
      complex*16 sum(0:nord,-nord:nord)
      complex*16 spher(0:nord,-nord:nord,nsamp)
      
      character*6 type

      common /dummy/coeff,spher

c----------------------------------------------------------------------
      pi = 4.d0*datan(1.d0)
      radcon = pi/180.d0
c----------------------------------------------------------------------

      elat = rlat1 
      elon = rlon1 
      slat = rlat2 
      slon = rlon2

c----------------------------------------------------------------------
c     compute pole to great circle trough event and station location
c----------------------------------------------------------------------
      call pole(elat,elon,slat,slon,polelat,polelon)

      bigthe = (90.d0 - polelat)*pi/180.d0
      bigphi = polelon*pi/180.d0           

c----------------------------------------------------------------------
c     compute coefficient of complex spherical harmonics
c----------------------------------------------------------------------
      do iss = 0,ismax
        do itt = 0,iss
          coeff(iss,itt) = dsqrt(dble(2*iss+1)*DFACTORL(iss-itt)
     &                  /(4.d0*pi*DFACTORL(iss+itt)))
        enddo
      enddo
c----------------------------------------------------------------------
c     GREAT circle , analytical average over the
c     path with pole (bigthe,bigphi) (Backus, 1964, eq.44)
c     tested against numerical integration  12/3/91   PI
c----------------------------------------------------------------------
      if (type.eq.'great ') then                    

        xxx = dcos(bigthe)  
        do iss = 0,ismax
          fact = DPLGNDRE(iss,0,0.d0)/DPLGNDRE(iss,0,1.d0)
          do itt = 0,iss
            alon = bigphi*dble(itt)
            aa = coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dcos(alon)
            bb = coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dsin(alon)
            avg_spher(iss,-itt) = fact*dcmplx(aa,-bb)*(-1.d0)**(itt)
            avg_spher(iss,itt)  = fact*dcmplx(aa, bb)
          enddo
        enddo

c----------------------------------------------------------------------
c     MINOR arc, numerical average
c----------------------------------------------------------------------
      elseif (type.eq.'minor ') then 

        call minor_arc(elat,elon,slat,slon,polelat,polelon,
     &                  numpts,alpha,min_lat,min_lon,min_ds)

        do i = 1,numpts+1
          xxx  = (90.d0 - min_lat(i))*radcon
          xxx  = dcos(xxx)
          yyy  = min_lon(i)*radcon
          do iss = 0,ismax
            do itt = 0,iss
              alon = yyy*dble(itt)
              aa = coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dcos(alon)
              bb = coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dsin(alon)
              spher(iss,-itt,i) = dcmplx(aa,-bb)*(-1.d0)**(itt)
              spher(iss,itt, i) = dcmplx(aa, bb)
            enddo
          enddo
        enddo

        do iss = 0,ismax
          do itt = -iss,iss
            sum(iss,itt) = (0.d0,0.d0)
            do i = 1, numpts
              sum(iss,itt) = sum(iss,itt)+ 0.5d0*(spher(iss,itt,i)
     &                       +spher(iss,itt,i+1))*min_ds(i)*radcon
            enddo
            avg_spher(iss,itt) = sum(iss,itt)/(alpha*radcon)
          enddo
        enddo                                

c----------------------------------------------------------------------
c     MAJOR arc, numerical average
c----------------------------------------------------------------------
      elseif (type.eq.'major ') then 

        call major_arc(elat,elon,slat,slon,polelat,polelon,
     &                  numpts,alpha,maj_lat,maj_lon,maj_ds)

        do i = 1,numpts+1
          xxx  = (90.d0 - maj_lat(i))*radcon
          xxx  = dcos(xxx)
          yyy  = maj_lon(i)*radcon
          do iss = 0,ismax
            do itt = 0,iss
              alon = yyy*dble(itt)
              aa=coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dcos(alon)
              bb=coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dsin(alon)
              spher(iss,-itt,i) = dcmplx(aa,-bb)*(-1.d0)**(itt)
              spher(iss,itt, i) = dcmplx(aa, bb)
            enddo
          enddo
        enddo

        do iss = 0,ismax
          do itt = -iss,iss
            sum(iss,itt) = (0.d0,0.d0)
            do i = 1, numpts
              sum(iss,itt) = sum(iss,itt)+ 0.5d0*(spher(iss,itt,i)+
     &                       spher(iss,itt,i+1))*maj_ds(i)*radcon
            enddo
            avg_spher(iss,itt) = sum(iss,itt)/(alpha*radcon)
          enddo
        enddo

      endif
      return
      end 

c----------------------------------------------------------------------
      subroutine avg_great_pole(polelat,polelon,ismax,avg_spher)
c----------------------------------------------------------------------
c same as before,only great circle average
c----------------------------------------------------------------------
      include 'harm_param.h'

      integer*4 ismax
      
      real*8 dplgndre,dfactorl
      real*8 aa,bb,xxx
      real*8 bigthe,bigphi
      real*8 pi,alon,fact

      real*8 polelat,polelon
      real*8 coeff(0:nord,-nord:nord)

      complex*16 avg_spher(0:nord,-nord:nord)
      
c----------------------------------------------------------------------
      pi = 4.d0*datan(1.d0)
c----------------------------------------------------------------------

      bigthe = (90.d0 - polelat)*pi/180.d0
      bigphi = polelon*pi/180.d0           

c----------------------------------------------------------------------
c     compute coefficient of complex spherical harmonics
c----------------------------------------------------------------------
      do iss = 0,ismax
        do itt = 0,iss
          coeff(iss,itt) = dsqrt(dble(2*iss+1)*DFACTORL(iss-itt)
     &                  /(4.d0*pi*DFACTORL(iss+itt)))
        enddo
      enddo
c----------------------------------------------------------------------
c     GREAT circle , analytical average over the
c     path with pole (bigthe,bigphi) (Backus, 1964, eq.44)
c     tested against numerical integration  12/3/91   PI
c----------------------------------------------------------------------

      xxx = dcos(bigthe)  
      do iss = 0,ismax
        fact = DPLGNDRE(iss,0,0.d0)/DPLGNDRE(iss,0,1.d0)
        do itt = 0,iss
          alon = bigphi*dble(itt)
          aa = coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dcos(alon)
          bb = coeff(iss,itt)*DPLGNDRE(iss,itt,xxx)*dsin(alon)
          avg_spher(iss,-itt) = fact*dcmplx(aa,-bb)*(-1.d0)**(itt)
          avg_spher(iss,itt)  = fact*dcmplx(aa, bb)
        enddo
      enddo

      return
      end 

c**********************************************************************
C
c                RELATED SUBROUTINES AND FUNCTION
c
c**********************************************************************
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine minor_arc(elat,elon,slat,slon,polelat,polelon,
     &                   npts,alpha,olat,olon,ds)
c----------------------------------------------------------------------
c     computes location at regular interval on the minor arc from 
c     event to station
c----------------------------------------------------------------------
      include 'harm_param.h'

      integer*4 npts 
      real*8 samp,alpha,rest
      real*8 nlat(nsamp),nlon(nsamp),olat(nsamp),olon(nsamp),ds(nsamp)
      real*8 elat,elon,enlat,enlon
      real*8 slat,slon,snlat,snlon
      real*8 polelat,polelon
      real*8 a(3,3),b(3,3)
                                         
c...computes the elements of the rotation matrices
c   "b" rotates great circle to equatorial great circle
c   "a" rotates equatorial great circle to great circle with pole
c    polelat,polelon

      call rot_matrices(polelat,polelon,a,b)

c...rotate event location towards equatorial circle
c...rotate station location towards equatorial circle

      call rotate(b,elat,elon,enlat,enlon)
      call rotate(b,slat,slon,snlat,snlon)

c...computes the position of the minor arc on equatorial circle
c...between rotated event and station locations
c...event/station/nord-pole form a right-coordinate system
c...the minor arc is in anticlockwise direction starting
c...from the rotated event location                                                
c...samp degrees sampling
c...find number of points

      samp  = 1.d0
      alpha = dabs(snlon-enlon)                          
      if (alpha.eq.180.d0) then
        print*,' WARNING non-uniqueness : station is anti-podal'
        stop
      endif
      if (alpha.gt.180.d0) alpha = 360.d0-alpha
      rest = dmod(alpha,samp)
      npts = idnint((alpha-rest)/samp)+1

c...computes longitude points and spacing

      do i = 1,npts
        nlat(i) = 0.d0
        nlon(i) = enlon + dble(i-1)*samp
        ds(i)   = samp
      enddo
      nlon(npts+1) = nlon(npts)+rest
      ds(npts) = rest

c...rotate back to great circle defined by event and station location

      do i = 1,npts+1
        call rotate(a,nlat(i),nlon(i),olat(i),olon(i))
      enddo

      return
      end

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine major_arc(elat,elon,slat,slon,polelat,polelon,
     &                   npts,alpha,olat,olon,ds)
c----------------------------------------------------------------------
c     computes location at regular interval on the MAJOR arc from 
c     event to station
c----------------------------------------------------------------------
      include 'harm_param.h'

      integer*4 npts 
      real*8 samp,alpha,rest
      real*8 nlat(nsamp),nlon(nsamp),olat(nsamp),olon(nsamp),ds(nsamp)
      real*8 elat,elon,enlat,enlon
      real*8 slat,slon,snlat,snlon
      real*8 polelat,polelon
      real*8 a(3,3),b(3,3)
                                         
c...computes the elements of the rotation matrices
c   "b" rotates great circle to equatorial great circle
c   "a" rotates equatorial great circle to great circle with pole
c    polelat,polelon

      call rot_matrices(polelat,polelon,a,b)

c...rotate event location towards equatorial circle
c...rotate station location towards equatorial circle

      call rotate(b,elat,elon,enlat,enlon)
      call rotate(b,slat,slon,snlat,snlon)

c...computes the position of the minor arc on equatorial circle
c...between rotated event and station locations
c...event/station/nord-pole form a right-coordinate system
c...the major arc is in anticlockwise direction starting
c...from the rotated station location                                                
c...samp degrees sampling
c...find number of points

      samp  = 1.d0
      alpha = dabs(snlon-enlon)                          
      if (alpha.eq.180.d0) then
        print*,' WARNING non-uniqueness : station is anti-podal'
        stop
      endif
      if (alpha.lt.180.d0) alpha = 360.d0-alpha
      rest = dmod(alpha,samp)
      npts = idnint((alpha-rest)/samp)+1

c...computes longitude points and spacing

      do i = 1,npts
        nlat(i) = 0.d0
        nlon(i) = snlon + dble(i-1)*samp
        ds(i)   = samp
      enddo
      nlon(npts+1) = nlon(npts)+rest
      ds(npts) = rest

c...rotate back to great circle defined by event and station location

      do i = 1,npts+1
        call rotate(a,nlat(i),nlon(i),olat(i),olon(i))
      enddo

      return
      end

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine midpt(elat,elon,slat,slon,dlat,dlon)
c----------------------------------------------------------------------
c     computes the midpoint on small arc between two points
c     event to station
c----------------------------------------------------------------------

      real*8 alpha,alpha2
      real*8 elat,elon,enlat,enlon
      real*8 slat,slon,snlat,snlon
      real*8 dlat,dlon,dnlat,dnlon
      real*8 polelat,polelon
      real*8 a(3,3),b(3,3)

      call pole(elat,elon,slat,slon,polelat,polelon)
                                         
c...computes the elements of the rotation matrices
c   "b" rotates great circle to equatorial great circle
c   "a" rotates equatorial great circle to great circle with pole
c    polelat,polelon

      call rot_matrices(polelat,polelon,a,b)

c...rotate event location towards equatorial circle
c...rotate station location towards equatorial circle

      call rotate(b,elat,elon,enlat,enlon)
      call rotate(b,slat,slon,snlat,snlon)

c...computes the position of the minor arc on equatorial circle
c...between rotated event and station locations
c...event/station/nord-pole form a right-coordinate system
c...the minor arc is in anticlockwise direction starting
c...from the rotated event location                                                
c...samp degrees sampling
c...find number of points

      alpha = dabs(snlon-enlon)                          
      if (alpha.eq.180.d0) then
        print*,' WARNING non-uniqueness : station is anti-podal'
        stop
      endif
      if (alpha.gt.180.d0) alpha = 360.d0-alpha
      alpha2 = 0.5d0*alpha
      dnlat = 0.d0
      dnlon = enlon+alpha2
c...rotate back to great circle with pole(polelat,polelon)
      call rotate(a,dnlat,dnlon,dlat,dlon)
      return
      end

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine rot_matrices(lat,lon,a,b)
c----------------------------------------------------------------------
c    lat lon --> position of rotation pole (in degrees)
c    a,b     --> elements of rotation matrix
c----------------------------------------------------------------------

      real*8 pi,lat,colat,lon,long,a(3,3),b(3,3)
       
      pi    = 4.d0*datan(1.d0)
      colat = (90.d0 - lat)*pi/180.d0
      long  = lon*pi/180.d0   
c elements of rotation matrix
c---> equatorial circle to great circle with pole lat,lon
      a(1,1) =  dcos(long)*dcos(colat)
      a(1,2) = -dsin(long)
      a(1,3) =  dcos(long)*dsin(colat)
      a(2,1) =  dsin(long)*dcos(colat)
      a(2,2) =  dcos(long)
      a(2,3) =  dsin(long)*dsin(colat)
      a(3,1) = -dsin(colat)
      a(3,2) =  0.d0
      a(3,3) =  dcos(colat)
c---> great circle to equatorial circle 
      b(1,1) =  dcos(long)*dcos(colat)
      b(1,2) =  dsin(long)*dcos(colat)
      b(1,3) = -dsin(colat)
      b(2,1) = -dsin(long)
      b(2,2) =  dcos(long)
      b(2,3) =  0.d0
      b(3,1) =  dcos(long)*dsin(colat)
      b(3,2) =  dsin(long)*dsin(colat)
      b(3,3) =  dcos(colat)

      return
      end

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine rotate(c,olat,olon,nlat,nlon)
c----------------------------------------------------------------------
c     c : elements of the 3x3 rotation matrix
c     olat,olon : old lat and lon (in degrees)
c     nlat,nlon : new lat and lon (in degrees)
c----------------------------------------------------------------------

      real*8 olat,olon,nlat,nlon,c(3,3)
      real*8 x,y,z,xn,yn,zn
      real*8 colat,lon,pi
                             
      pi = 4.d0*datan(1.d0)

      colat = (90.d0-olat)*pi/180.d0
      lon   = olon*pi/180.d0

      x = dcos(lon)*dsin(colat)
      y = dsin(lon)*dsin(colat)
      z = dcos(colat)

      xn = x*c(1,1) + y*c(1,2) + z*c(1,3)
      yn = x*c(2,1) + y*c(2,2) + z*c(2,3)
      zn = x*c(3,1) + y*c(3,2) + z*c(3,3)

      nlat = 90.d0 - dacos(zn)*180.d0/pi
      if ((xn.eq.0.d0).and.(yn.eq.0.d0)) then
        nlon = 0.d0
      else 
        nlon = datan2(yn,xn)*180.d0/pi
      endif                   

      return
      end

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine pole(elat,elon,slat,slon,polelat,polelon)
c----------------------------------------------------------------------
c input  : event lat long  --> elat elon (in degrees)
c          station         --> slat slon (    "     )
c output : pole  lat long  --> polelat polelon (in degrees)
c----------------------------------------------------------------------
      real*8 elat,elon,slat,slon,etheta,stheta,ephi,sphi
      real*8 polelat, polelon, exyz(3), rxyz(3), polexyz(3)
      real*8 temp, radcon, pi

      pi     = 4.d0*datan(1.d0)
      radcon = pi / 180.d0
c conversion to radians
      etheta = (90.d0 - elat)* radcon
      ephi   = elon          * radcon
      stheta = (90.d0 - slat)* radcon
      sphi   = slon          * radcon
c cartesian coordinates
      exyz(1) = dsin(etheta)*dcos(ephi)
      exyz(2) = dsin(etheta)*dsin(ephi)
      exyz(3) = dcos(etheta)
      rxyz(1) = dsin(stheta)*dcos(sphi)
      rxyz(2) = dsin(stheta)*dsin(sphi)
      rxyz(3) = dcos(stheta)    
c cross product of exyx and rxyz
      polexyz(1) = exyz(2) * rxyz(3) - exyz(3) * rxyz(2)
      polexyz(2) = exyz(3) * rxyz(1) - exyz(1) * rxyz(3)
      polexyz(3) = exyz(1) * rxyz(2) - exyz(2) * rxyz(1)
c polar coordinates
      temp    = dsqrt(polexyz(1)**2+polexyz(2)**2 +polexyz(3)**2)
      polelat = (pi/2.d0 - dacos(polexyz(3)/temp))/radcon
      polelon = datan2(polexyz(2),polexyz(1))/radcon
      
      return
      end

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     calculate n!
c----------------------------------------------------------------------
      double precision FUNCTION DFACTORL(N)

      REAL*8 A(33),dgammaln

      DATA NTOP,A(1)/0,1.d0/

      IF (N.LT.0) THEN
        PAUSE 'negative factorial'
      ELSE IF (N.LE.NTOP) THEN
        DFACTORL=A(N+1)
      ELSE IF (N.LE.32) THEN
        DO 41 J=NTOP+1,N
          A(J+1)=dble(J)*A(J)
 41      CONTINUE
        NTOP=N
        DFACTORL=A(N+1)
      ELSE
        DFACTORL=DEXP(DGAMMALN(N+1.))
      ENDIF
      RETURN
      END 

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     compute the gamma function
c----------------------------------------------------------------------
      double precision FUNCTION DGAMMALN(XX)

      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER

      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/

      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*DLOG(TMP)-TMP
      SER=ONE
      DO 51 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
 51   CONTINUE
      DGAMMALN=TMP+DLOG(STP*SER)
      RETURN
      END
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     compute legendre polynomials
c----------------------------------------------------------------------
      double precision FUNCTION DPLGNDRE(L,M,X)

      real*8 x,somx2,fact,pmm,pmmp1,pll

      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.)PAUSE 'bad arguments'
      PMM=1.d0
      IF(M.GT.0) THEN
        SOMX2=SQRT((1.d0-X)*(1.d0+X))
        FACT=1.d0
        DO 61 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.d0
 61     CONTINUE
      ENDIF
      IF(L.EQ.M) THEN
        DPLGNDRE=PMM
      ELSE
        PMMP1=X*(dble(2*M)+dble(1))*PMM
        IF(L.EQ.M+1) THEN
          DPLGNDRE=PMMP1
        ELSE
          DO 62 LL=M+2,L
            PLL=(X*dble(2*LL-1)*PMMP1-dble(LL+M-1)*PMM)/dble(LL-M)
            PMM=PMMP1
            PMMP1=PLL
 62       CONTINUE
          DPLGNDRE=PLL
        ENDIF
      ENDIF
      RETURN
      END


