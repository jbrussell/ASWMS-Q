c** study the inversion grid
      character*70 freggrid
      print *,'enter grid file'
      read(*,'(a)') freggrid
      call genreggrid(freggrid,ntot,nxpt,dx,dy)
      print *,ntot,nxpt,dx,dy       
      end

      subroutine genreggrid(freggrid,ntot,nxpt,dx,dy)
!***************************************************************************
!*   Create regularly spaced grid points in an area with center point (slat,slon). *
!*   In order to avoid distortion at high latitude, change coordinates so   *
!*   center of area lies along the equator of coordinate system and there is 
!*   compensation for curvature.                                           *
!*************************************************************************** 
!  WARNING   checkout scheme by running separately first gridreg < brwestreginp.dat
       parameter (nxptmx=100, nyptmx=100)
       parameter (maxnodes = 2000)

       real*4 dazim(nxptmx),ddelt(nyptmx),glat(nxptmx,nyptmx)
       real*4 glon(nxptmx,nyptmx)
      real*4 boxlat(4), boxlon(4)
      real*4 nodelat(maxnodes),nodelon(maxnodes)
       character*70 freggrid
       common /gengrid/ nodelat,nodelon,boxlat,boxlon
       open(15, file = freggrid, status = 'unknown')
       radius = 6371.
       pi = 3.14158
       convdeg = 2.0*pi/360.
!   'Enter the center point (slat, slon):'
       read(15,*) slat,slon
       sdelta=90.
!      sazimz is azimuth from the center point that will tilt grid relative to North
       read(15,*) sazimz
!  delx is increment in degrees at equator - the true increment will vary with latitude from
!  equator of projection
       read(15,*) nxpt, delx,begx
! to follow our traditional convention of listing points from right to left and bottom to top.
!  the degree increment for distance from the pole, dely, should be negative and the
!  beginning point should be farthest in degrees from projection pole (assuming sazimz is
!  northwards), (e.g., begy = 6 (degrees south of center point) is 96 degrees from pole)
!  begx would be negative and increase with positive delx
       read(15,*) nypt,dely, begy
! read in corners for wave intercepts to pass along to main program
       do i = 1, 4
         read(15, *) boxlat(i), boxlon(i)       
         write(57,*) boxlat(i), boxlon(i)
       enddo
       degdist = 2.0*pi*radius/360.
       dx = delx*degdist
       dy = dely*degdist
       call gohead(slat,slon,sdelta,sazimz,plat,plon)
       call disthead(plat,plon,slat,slon,delta,azimz)
       ntot = nxpt*nypt
       do 150 j=1,nypt
         delyinc = (j-1)*dely +begy     
         delt = delta + delyinc 
         do 100 i=1,nxpt
           azim = azimz + (begx + (i-1)*delx)/cos(convdeg*delyinc)
           call gohead(plat,plon,delt,azim,glat(i,j),glon(i,j))
           if (glon(i,j).gt.180.0) then
             glon(i,j) = glon(i,j) - 360.0
           endif
  100    enddo
  150  enddo
!  arrange points in standard order in single dimensioned arrays
       nodecnt = 0
       do i = 1, nxpt
         do j = 1, nypt
           nodecnt = nodecnt + 1
           nodelat(nodecnt) = glat(i,j)
           nodelon(nodecnt) = glon(i,j)
           write(56,200) glat(i,j),glon(i,j)
         enddo
       enddo
       close(15)
  200  format(f7.3,f10.3)
       return
       end


      subroutine disthead(slat,slon,flat,flon,delta,azim)
!  Calculates distance and azimuth on sphere from starting point s 
!  to finishing point f
      dtor= 3.1415928/180.
      slt = slat*dtor
      sln = slon*dtor
      flt = flat*dtor
      fln = flon*dtor
      delta = acos(sin(slt)*sin(flt)+cos(slt)*cos(flt)*cos(fln-sln))
      azim = atan2(sin(fln-sln)*cos(flt),
     1   sin(flt)*cos(slt) - cos(fln-sln)*cos(flt)*sin(slt))
      delta = delta/dtor
      azim = azim/dtor
      return
      end

      subroutine gohead(slat,slon,delta,azim,flat,flon)
! Calculates final latitude and longitude f when starting at point s
! traveling delta degrees on sphere at initial heading azim
      dtor= 3.1415928/180.
      slt = slat*dtor
      dlta = delta*dtor
      azm = azim*dtor
      flat = asin(cos(slt)*sin(dlta)*cos(azm) + sin(slt)*cos(dlta))
      flon = atan2(sin(dlta)*sin(azm),
     1    cos(slt)*cos(dlta) - sin(slt)*sin(dlta)*cos(azm))
      flat = flat/dtor
      flon = slon + flon/dtor
      if (flon.gt.360.) flon = flon - 360.
      if (flon.lt.-360.) flon = flon + 360.
      return
      end 
