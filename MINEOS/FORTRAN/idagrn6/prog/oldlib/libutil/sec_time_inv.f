      subroutine sec_time_inv(tsec,ibyr,iyr,idoy,ihr,imin,sec)
c
c     subroutine to convert time in sec relative to a base year
c     to time in y,d,h,m,s.  Inverts the action of sec_tim.f.
c
c     input tsec and ibyr (reference year for t=0)
c     outut time in yr,d,h,m,s
c
c     NOTE use dsec_time and dsec_time_inv for more accuracy
c  
c
      integer*4 iyr, idoy, ihr, imin, ibyr, lpyr
      real*4 sec,tsec
c
      include 'numerical.h'
c
c     first account for year offset

      ierr = lpyr(ibyr)
      if (ierr .eq. 1) then
        iday = 366
      else
        iday = 365
      end if

      spyr = real(iday)*spday
c      nyr = int(tsec/spyr)
c      print*,nyr

c     add 1 to base year, and subtract correct number of seconds,
c     for each whole year of offset.  Requires checking for leap
c     year and calculating spyr with each step.  Loop is not 
c     entered if tsec is less than a year.  

      iyr = ibyr
      do while (tsec.ge.spyr)
c          i = 1,nyr
        tsec = tsec - spyr
        iyr = iyr + 1
        ierr = lpyr(iyr)
        if (ierr .eq. 1) then
          iday = 366
        else
          iday = 365
        end if
        spyr = real(iday)*spday
      end do

c     year should now be correct, and tsec is less than 1 year.  Note
c     that idoy starts at 1 since there is no day 0.

      idoy = 1 + int(tsec/spday)
      tsec = tsec - real(idoy-1)*spday 
      ihr =  int(tsec/sphr)
      tsec = tsec - real(ihr)*sphr
      imin = int(tsec/spmin)
      sec =  real(tsec) - real(imin)*spmin

      if(idoy.gt.iday .or. ihr.gt.24 .or. imin.gt.60 
     &    .or. sec.gt.60. .or. iday.lt.0 .or. ihr.lt.0 .or.
     &     imin.lt.0 .or. sec.lt.0.) then
        print*,'**********'
        print*,'ERROR in sec_time_inv -- failed to reduce'
        print*,'**********'
      end if        

c
      return
      end
