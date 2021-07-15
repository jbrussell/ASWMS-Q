      subroutine sec_time(iyr,idoy,ihr,imin,sec,ibyr,tsec)
c
c     subroutine to convert time to sec relative to a base year
c     input time in gfs format (yr,day,hr,min,sec) and a reference
c     year for t=0
c     returns tsec, time in sec
c
c     NOTE tsec should be double precision in calling program
c     to ensure accuracy.
c  
      integer*4 iyr, idoy, ihr, imin, ibyr
      real*4 sec,tsec
      real*8 ttmp
c
      include 'numerical.h'
c
c     first account for day, hour, minutes, and seconds
c     subtract 1 from idoy since idoy starts at 1, not zero
c
      ttmp = dble(idoy-1)*spday + dble(ihr)*sphr
     .   + dble(imin)*spmin
     .   + sec
c
c     now calculate the year offset
c
      icount = 0
      do ii = 1, iyr-ibyr
        ierr = lpyr(ibyr + ii - 1)
        if (ierr .eq. 1) then
          iday = 366
        else
          iday = 365
        end if
        icount = icount + iday
c         print*,ibyr + ii - 1,iday,icount
      end do
      ttmp = ttmp + dble(icount)*spday
      tsec = sngl(ttmp)
c
      return
      end
