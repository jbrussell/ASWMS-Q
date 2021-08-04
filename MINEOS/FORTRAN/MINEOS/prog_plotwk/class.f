      subroutine class(type,icolor_fill)
c
c     subroutine to classify spheroidal modes 
c     loosely based on the criteria of Okal (1978)
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop)
      real*4 wmin, wmax, pmin, pmax, gmin, gmax
      real*4 fill(0:maxmodes)
      integer*4 ifill(0:maxmodes),cfill(0:maxmodes)
      integer*4 inx(0:maxl,0:maxn), inx2(0:maxmodes,2)
      integer*4 jcom, nmodes, nlow, nup, llow, lup
      integer*4 icolor_fill
      logical lsummary
      character*1 type
c
      common /mode/ inx, modes, inx2
      common /fill/ lsummary, fill, ifill, cfill
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
c
c
      if ((type .eq. 'k') .or. (type .eq. 'K')) then
        do ii = 1, nmodes
          if (cfill(ii) .le. 2) then
            qq = modes(ii,2)
            cv = modes(ii,3)
            gv = modes(ii,4)
            if (cv .ge. 200.0) then
              if (gv .ge. 22.0) then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              endif
            elseif (cv .ge. 50.0) then
              if (gv .ge. 20.0) then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              endif
            elseif ((cv .ge. 16.5) .and. (cv .lt. 50.)) then
              if (gv .ge. 16.0) then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              endif
            elseif ((cv .gt. 11.) .and. (cv .lt. 16.5)) then
              if (gv .ge. 14.)  then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              endif
            elseif (cv .le. 11.) then
              if ((gv .ge. 8.) .and. (gv .lt. 10.)) then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              elseif ((gv .ge. 14.) .and. (gv .lt. 16.)) then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              endif
            endif         
          endif
        end do
        icolor_fill = 4
      elseif ((type .eq. 'p') .or. (type .eq. 'P')) then
        do ii = 1, nmodes
          if (cfill(ii) .le. 2) then
            qq = modes(ii,2)
            cv = modes(ii,3)
            gv = modes(ii,4)
            cs = modes(ii,24)
            ce = modes(ii,6)
            if ((cv .gt. 14.76) .and. (cv .le. 25.08)) then
              if (gv .le. 14.) then
                if (cs .le. 0.80) then
                  if (ce .le. 0.5) then
                    ifill(ii) = 1
                    fill(ii)  = 1.
                    cfill(ii) = 7
                  endif
                endif
              endif
            elseif ((cv .gt. 11.94) .and. (cv .le. 14.76)) then
              if (gv .le. 14.) then
                if (cs .le. 0.80) then
                  if (ce .le. 0.5) then
                    ifill(ii) = 1
                    fill(ii)  = 1.
                    cfill(ii) = 7
                  endif
                endif
              endif
            elseif ((cv .gt. 8.23) .and. (cv .le. 11.94)) then
              if (gv .le. 8.) then
                if (cs .le. 0.86) then
                  if (ce .le. 0.5) then
                    ifill(ii) = 1
                    fill(ii)  = 1.
                    cfill(ii) = 7
                  endif
                endif
              elseif ((gv .ge. 10.) .and. (gv .le. 14.)) then
                if (cs .le. 0.86) then
                  if (ce .le. 0.5) then
                    ifill(ii) = 1
                    fill(ii)  = 1.
                    cfill(ii) = 7
                  endif
                endif
              endif
            endif
          endif
        end do
        icolor_fill = 7
      elseif ((type .eq. 'v') .or. (type .eq. 'V')) then
        do ii = 1, nmodes
          if (cfill(ii) .le. 2) then
            qq = modes(ii,2)
            cv = modes(ii,3)
            gv = modes(ii,4)
            if (cv .ge. 200.0) then
              if ((gv .lt. 22.0) .and. (gv .gt. 8))then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              endif
            elseif (cv .ge. 50.0) then
              if ((gv .lt. 20.0) .and. (gv .gt. 8)) then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              endif
            elseif ((cv .ge. 25.08) .and. (cv .lt. 50.)) then
              if ((gv .lt. 16.0) .and. (gv .gt. 8.0)) then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              endif
            endif         
          endif
        end do
        icolor_fill = 5
      elseif ((type .eq. 'c') .or. (type .eq. 'C')) then
        do ii = 1, nmodes
          if (cfill(ii) .le. 2) then
            qq = modes(ii,2)
            cv = modes(ii,3)
            gv = modes(ii,4)
            cs = modes(ii,21)
            if (cv .gt. 25.07) then
              if (gv .lt. 8.0) then
                ifill(ii) = 1
                fill(ii)  = 1.
                cfill(ii) = 7
              endif
            endif         
          endif
        end do
        icolor_fill = 3
      elseif ((type .eq. 'r') .or. (type .eq. 'R')) then
        do ii = 1, nmodes
          if (cfill(ii) .le. 2) then
            qq = modes(ii,2)
            cv = modes(ii,3)
            gv = modes(ii,4)
            cs = modes(ii,24)
            ce = modes(ii,6)
            if ((cv .gt. 14.76) .and. (cv .lt. 25.08)) then
              if (gv .le. 14.) then
c                if (cs .gt. 0.80) then
c                  if (ce .le. 0.5) then
                    ifill(ii) = 1
                    fill(ii)  = 1.
                    cfill(ii) = 7
c                  endif
c                endif
              endif
            elseif ((cv .gt. 11.94) .and. (cv .le. 14.76)) then
              if (gv .le. 14.) then
c                if (cs .gt. 0.80) then
c                  if (ce .le. 0.5) then
                    ifill(ii) = 1
                    fill(ii)  = 1.
                    cfill(ii) = 7
c                  endif
c                endif
              endif
            elseif (cv .le. 11.94) then
              if (gv .le. 8.) then
c                if (cs .gt. 0.86) then
c                  if (ce .le. 0.5) then
                    ifill(ii) = 1
                    fill(ii)  = 1.
                    cfill(ii) = 7
c                  endif
c                endif
              elseif ((gv .ge. 10.) .and. (gv .le. 14.)) then
c                if (cs .gt. 0.86) then
c                  if (ce .le. 0.5) then
                    ifill(ii) = 1
                    fill(ii)  = 1.
                    cfill(ii) = 7
c                  endif
c                endif
              endif
            endif
          endif
        end do
        icolor_fill = 6
      elseif ((type .eq. 'a') .or. (type .eq. 'A')) then
        print*,' type A not recognized'
      endif
c
      return
      end
