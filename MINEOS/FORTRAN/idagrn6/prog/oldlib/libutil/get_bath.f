      subroutine get_bath(elat,elon,fbath)
c
c     given elat and elon of a point, return value of bathymetry
c
      integer*2 ibath1(4320),ibath2(4320)
      character*256 filen
c
      filen='/seismo/data1/datalib/dbdb5/dbdb5.data.sun'
      lu = 173
c
c     open DBDB5 bathymetry file and set up some parameters
c     - 5 min by 5 min
c
      open(lu,file=filen,form='unformatted',recl=8640,access='direct') 
      spd = 1.0/12.0
      xlatm = 90.0
      xlonm = 0.0
      xlat = elat
      if (elon .lt. 0.0) then
        xlon = elon + 360
      else
        xlon = elon
      end if
c
c     determine the latitude records
c
      rlat = (xlatm - xlat)/spd + 1.0
      i1 = int(rlat)
      i2 = i1 + 1
      xlati = xlatm  - xlat - float(i1 - 1)*spd
      read(lu,rec=i1) (ibath1(i),i=1,4320)
      read(lu,rec=i2) (ibath2(i),i=1,4320)
c
c     determine the longitude elements
c
      rlon = (xlon - xlonm)/spd + 1.0
      i1 = int(rlon)
      i2 = i1 + 1
      xloni = xlon - xlonm  - float(i1 - 1)*spd
c
c     we have bounded the value - interpolate
c
      f1 = ibath1(i1)
      f2 = ibath2(i1)
      f3 = ibath1(i2)
      f4 = ibath2(i2)
      f5 = f1 + (xloni/spd)*(f3 - f1)
      f6 = f2 + (xloni/spd)*(f4 - f2)
      fbath = f5 + (xlati/spd)*(f6 - f5)
c     
      close(lu)
      return
      end     
