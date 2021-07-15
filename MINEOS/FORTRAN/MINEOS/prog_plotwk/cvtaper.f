c
c
c
cSSS  Steven S. Shapiro
cSSS  4 January 1991                
cSSS  Subroutine "cvtaper" performs a taper in phase velocity space.
cSSS  Function:
cSSS           Tapers phase velocity from cv0 to cv1 using Hanning window
cSSS           Sets phase velocities between cv1 and cv2 inclusive = 1
cSSS           Tapers phase velocity from cv2 to cv3 using Hanning window
cSSS
      subroutine cvtaper (cv0, cv1, cv2, cv3)
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop)
      real*4 fill(0:maxmodes)
      real*4
     &      tcv,
c  the phase velocity
     &      cv0,
     &      cv1,
     &      cv2,
     &      cv3,
     &      pi
      integer*4 ifill(0:maxmodes)
      integer*4 cfill(0:maxmodes)
      integer*4 inx(0:maxl,0:maxn)
      integer*4 inx2(0:maxmodes,2)
      logical lsummary
c
      common /mode/ inx, modes, inx2
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /fill/ lsummary, fill, ifill, cfill
c
      data pi / 3.141592654/                         
      if ((cv0 .gt. cv1) .or. (cv1 .gt. cv2) .or. (cv2 .gt. cv3)) then
         print*, 'Error in "cvtaper": values are not in ascending order.
     &.'
         print*, 'Redo command.'
         return
      else
      end if
      do ii = 1, nmodes
         if (ifill(ii) .ne. 0) then
c  Mode has been selected
            tcv = modes(ii,3)
            if ((tcv .ge. cv1) .and. (tcv .le. cv2)) then
               fill(ii) = 1.0
            else if ((tcv .lt. cv1) .and. (tcv .ge. cv0)) then
               fill(ii) = 0.5 * (1 - cos (pi * (tcv - cv0)
     &                     / (cv1 - cv0)))
            else if ((tcv .gt. cv2) .and. (tcv .le. cv3)) then
               fill(ii) = 0.5 * (1 - cos (pi * (tcv - cv3)
     &                     / (cv2 - cv3)))
            else  
               fill(ii) = 0.0
            end if
         else
         endif
      end do
      call color(nmodes)
      return
      end
