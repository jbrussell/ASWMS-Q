c
c
c
      subroutine fix_class_c (b_file,icolor_fill)
c
c     subroutine to fix classification of spheroidal modes
c     attempt to make clever branch checking
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop)
      real*4 wmin, wmax, pmin, pmax, gmin, gmax
      real*4 fill(0:maxmodes), temp(4)
      real*4 tpi
c
      integer*4 ifill(0:maxmodes),cfill(0:maxmodes)
      integer*4 inx(0:maxl,0:maxn), inx2(0:maxmodes,2)
      integer*4 index(0:maxl,0:maxn), lost(0:maxl,0:maxn)
      integer*4 nnb(0:maxl,nbranch), llb(0:maxl,nbranch), numb(nbranch)
      integer*4 nt(0:maxl,nbranch), lt(0:maxl,nbranch)
      integer*4 jcom, nmodes, nlow, nup, llow, lup
      integer*4 icolor_fill
c
      logical lsummary, lread, lbranch
c
      character*1 type, ans
      character*256 b_file
c
      common /mode/ inx, modes, inx2
      common /fill/ lsummary, fill, ifill, cfill
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /branch/ lbranch, nb, numb, nnb, llb
      common /lost/ lost
c
      data tpi /6.2831853071796/
c
      lbranch = .true.
c
c     first, zero any previous branch info
c 
      do ii = 1, nbranch
        do jj = 0, maxl
          nnb(jj,ii) = 0
          llb(jj,ii) = 0
          nt(jj,ii)  = 0
          lt(jj,ii)  = 0
          index(jj,ii) = 0
          lost(jj,ii) = 0
        end do
      end do
      ib = 0
c
      do ii = 1, nmodes
        if ((cfill(ii) .eq. 7) .or. (cfill(ii) .eq. 8)) then
          index(inx2(ii,1),inx2(ii,2)) = 6
        else
          index(inx2(ii,1),inx2(ii,2)) = 0
        end if
        cfill(ii) = 0
      end do
c
c     sorting - first cut - order the branches
c
      nb = 0
      do il = llow, lup
        ib = 0
        do in = nlow, nup
          if (index(il,in) .eq. 6) then
            ib = ib + 1
          endif
        end do
        nb = max0(ib, nb)
      end do
c
      do il = llow, lup
        ib = 0
        do in = nup, nlow, -1
          if (index(il,in) .eq. 6) then
            ib = ib + 1
            nt(il,nb + 1 - ib) = in
            lt(il,nb + 1 - ib) = il
          endif
        end do
      end do
c
c     everything is fixed for the first branch
c
      do ib = 1, nb
        nnb(1,ib) = nt(llow,ib)
        llb(1,ib) = lt(llow,ib)
      end do 
      do il = llow + 1, lup
        do ib = 1, nb
          do ii = il - 1, 1, -1
            n1 = nnb(ii, ib)
            l1 = llb(ii, ib)
            if (n1 .ne. maxn) then
              wl = modes(inx(l1, n1), 1)
              gl = modes(inx(l1, n1), 4)
              fl = float(l1)
              slope = gl*co*((2.*fl+1.)/sqrt(fl*(fl + 1.)))
              be = wl - slope*fl
              we = slope*float(il) + be
              go to 101
            endif
          end do
101       continue
          dwa = 99999.
          do jj = 1, nb
            n2 = nt(il,jj)
            l2 = lt(il,jj)
            wc = modes(inx(l2,n2),1)
            dw = abs(wc - we)
            if (dw .lt. dwa) then
              dwa = dw
              jsave = jj
            endif
          end do
          if (dwa .lt. 0.25) then
            nnb(il,ib) = nt(il,jsave)
            llb(il,ib) = lt(il,jsave)
            nt(il, jsave) = 0
            lt(il, jsave) = 0
          else
            nnb(il,ib) = maxn        
            llb(il,ib) = il
          endif
        end do
      end do
c
      do ib = 1, nb
        do il = lup, llow, -1
          if (nnb(il,ib) .ne. maxn) then
            go to 102
          else
            nnb(il,ib) = 0
            llb(il,ib) = 0
          end if
        end do
102     continue
        numb(ib) = il
      end do
cccc
cccc  debugging
cccc
c      print*,' '
c      do ii = 1, nb
c        do il = 1, numb(ii)
c          nn = nnb(il,ii)
c          ll = llb(il,ii)
c          if (nn .ne. maxn) then
c            ind = inx(ll,nn)
c            write(6,'(4i5,4f15.5)') ii,il,nn,ll,(modes(ind,jj),jj = 1,4)
c          else
c            write(6,'(4i5)') ii,il,nn,ll
c          endif
c        end do
c        print*,' '
c      end do
cccc
cccc
cccc
c
c     find the lost modes
c
      do ii = 1, nmodes
        nn = inx2(ii,2)
        ll = inx2(ii,1)
        if (index(ll,nn) .eq. 6) then
          do ib = 1, nb
            do im = 1, numb(ib)
              nc = nnb(im, ib)
              lc = llb(im, ib)
              if ((nc .eq. nn) .and. (lc .eq. ll)) then
                 go to 998
              endif
            end do
          end do
          lost(ll,nn) = 1
        endif
998     continue         
      enddo
c
c     try to correlate skipped modes with lost modes
c
      print*,' '
      print*,' checking for lost modes'
      do ib = 1, nb
        im = 0
        do il = 1, lup
          nc = nnb(il, ib)
          lc = llb(il, ib)
c
c         this is a skipped mode
c
          if (nc .eq. maxn) then
            do ii = il - 1, 1, -1
              n1 = nnb(ii, ib)
              l1 = llb(ii, ib)
              if (n1 .ne. maxn) then
                wl = modes(inx(l1, n1), 1)
                gl = modes(inx(l1, n1), 4)
                fl = float(l1)
                slope = gl*co*((2.*fl+1.)/sqrt(fl*(fl + 1.)))
                be = wl - slope*fl
                we = slope*float(lc) + be
                dwl = slope
                go to 1000
              endif
            end do
1000        continue
c               
            do ii = nlow, nup
              if (lost(lc,ii) .gt. 0) then
                wc = modes(inx(lc,ii),1)
                if (abs(we - wc) .lt. 0.25) then
                  nnb(il, ib) = ii
                  print*,' correlate: ', ii, il, wl, wc, wh
                  lost(il, ii) = 0
                  go to 999                   
                else
c                 print*,' not quite: ', we, wc
                endif
              endif
            end do
          endif
999       continue         
          if (lc .ne. 0) then
            im = im + 1
          endif
        enddo
        numb(ib) = im
      enddo
cccc
      print*,' '
      do il = llow, lup
        do ib = nlow, nup
          if (lost(il,ib) .gt. 0) then
            print*,' lost mode: ', ib, il, modes(inx(il,ib),1)
          endif
        end do
      end do
cccc
cccc  debugging
cccc
c      print*,' '
c      do ii = 1, nb
c        do il = 1, numb(ii)
c          nn = nnb(il,ii)
c          ll = llb(il,ii)
c          if (nn .ne. maxn) then
c            ind = inx(ll,nn)
c            write(6,'(4i5,4f15.5)') ii,il,nn,ll,(modes(ind,jj),jj = 1,4)
c          else
c            write(6,'(4i5)') ii,il,nn,ll
c          endif
c        end do
c        print*,' '
c      end do
cccc
cccc
cccc
c
c     loop over branches and assign rotating color values
c
      do ii = 1, nb
        do il = 1, numb(ii)
          nn = nnb(il,ii)
          ll = llb(il,ii)
          if (nn .ne. maxn) then
            ind = inx(ll,nn)
            cfill(ind) = 3 + (ii - (ii/3)*3)
          endif
        end do
      end do
c
c     fill the lost modes in green
c
      do ii = nlow, nup
        do il = llow, lup
          if (lost(il,ii) .gt. 0) then
            ind = inx(il,ii)
            cfill(ind) = 9
          endif
        end do
      end do

c
c     now write output file containing branch information
c       two header records, followed by one record for each branch
c       each branch record contains a 6-tuple of info about each mode - nn, ll, w, q, cv, gv
c
      open(unit=4,file=b_file,form='unformatted',access='sequential')
      write(4) jcom, nmodes, nlow, nup, llow, lup
      write(4) nb, (numb(ii), ii = 1, nb)
      do ii = 1, nb
        write(4) (nnb(jj,ii), llb(jj,ii), 
     &    (modes(inx(llb(jj,ii),nnb(jj,ii)),kk),kk=1,4),jj=1,numb(ii))
      end do
      close(4)
c
      return
      end
