c
c
c
      subroutine fix_class_p (b_file,icolor_fill)
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
      integer*4 icolor_fill, nl(0:maxl)
c
      logical lsummary, lbranch
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
      const = (1000./tpi)*(1.0/(2.*6371.))
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
      lm = 999
      do il = llow, lup
        ib = 0
        do in = nlow, nup
          if (index(il,in) .eq. 6) then
            lm = min0(lm, il)
            ib = ib + 1
          endif
        end do
   5    continue
        nb = max0(ib, nb)
        nl(il) = nlow
      end do
c
      do ib = 1, nb
        im = 0
        do il = lm, lup
          do in = nl(il), nup
            if (index(il,in) .eq. 6) then
              go to 10
            end if
          end do
          if (im .ne. 0) then
            im = im + 1   
            nnb(im,ib) = maxn
            llb(im,ib) = il
            nl(il) = nup
          endif
          go to 50
  10      continue
          if (im .eq. 0) then            
            wl = modes(inx(il,in),1)
            if (ib .gt. 1) then
              wt = modes(inx(llb(1,ib-1),nnb(1,ib-1)),1)
              if (wl .lt. wt) then
                go to 50
              endif
            endif
            gl = modes(inx(il,in),4)
            fl = float(il)
            slope = gl*((2.*fl+1.)/sqrt(fl*(fl + 1.)))*const
            be = wl - slope*fl
            im = im + 1
            nnb(im,ib) = in
            llb(im,ib) = il
            nl(il) = in + 1
          else
            fl = float(il)
            we = slope*fl + be
            wc = modes(inx(il,in),1)
            gc = modes(inx(il,in),4)
            do ii = in + 1, nup
              if (index(il,ii) .eq. 6) then
                go to 20
              end if
            end do
  20        continue
            wf = modes(inx(il,ii),1)
            gf = modes(inx(il,ii),4)
            if (abs(we - wc) .lt. 0.75) then
              dg = abs(gc - gl) 
              if (abs(we - wf) .lt. 0.75) then
                df = abs(gf - gl) 
                if (dg .lt. df) then
                  im = im + 1
                  nnb(im,ib) = in
                  llb(im,ib) = il
                  wl = wc
                  gl = gc
                  fl = float(il)
                  slope = gl*((2.*fl+1.)/sqrt(fl*(fl + 1.)))*const
                  be = wl - slope*fl
                  nl(il) = ii + 1
                else
                  im = im + 1
                  nnb(im,ib) = ii
                  llb(im,ib) = il
                  wl = wf
                  gl = gf
                  fl = float(il)
                  slope = gl*((2.*fl+1.)/sqrt(fl*(fl + 1.)))*const
                  be = wl - slope*fl
                  nl(il) = ii + 1
                endif
              else
                im = im + 1
                nnb(im,ib) = in
                llb(im,ib) = il
                wl = wc
                gl = gc
                fl = float(il)
                slope = gl*((2.*fl+1.)/sqrt(fl*(fl + 1.)))*const
                be = wl - slope*fl
                nl(il) = in + 1
              endif
            elseif (abs(we-wf) .lt. 0.75) then
              wc = wf
              gc = gf
              in = ii
              do ii = in + 1, nup
                if (index(il,ii) .eq. 6) then
                  go to 30
                end if
              end do
  30          continue
              wf = modes(inx(il,ii),1)
              gf = modes(inx(il,ii),4)
              if (abs(we - wc) .lt. 0.75) then
                dg = abs(gc - gl) 
                if (abs(we - wf) .lt. 0.75) then
                  df = abs(gf - gl) 
                  if (dg .lt. df) then
                    im = im + 1
                    nnb(im,ib) = in
                    llb(im,ib) = il
                    wl = wc
                    gl = gc
                    fl = float(il)
                    slope = gl*((2.*fl+1.)/sqrt(fl*(fl + 1.)))*const
                    be = wl - slope*fl
                    nl(il) = in + 1
                  else
                    im = im + 1
                    nnb(im,ib) = ii
                    llb(im,ib) = il
                    wl = wf
                    gl = gf
                    fl = float(il)
                    slope = gl*((2.*fl+1.)/sqrt(fl*(fl + 1.)))*const
                    be = wl - slope*fl
                    nl(il) = ii + 1
                  endif
                else
                  im = im + 1
                  nnb(im,ib) = in
                  llb(im,ib) = il
                  wl = wc
                  gl = gc
                  fl = float(il)
                  slope = gl*((2.*fl+1.)/sqrt(fl*(fl + 1.)))*const
                  be = wl - slope*fl
                  nl(il) = in + 1
                endif
              else
                im = im + 1
                nnb(im,ib) = in
                llb(im,ib) = il
                wl = wc
                gl = gc
                fl = float(il)
                slope = gl*((2.*fl+1.)/sqrt(fl*(fl + 1.)))*const
                be = wl - slope*fl
                nl(il) = in + 1
              endif
            else
              im = im + 1   
              nnb(im,ib) = maxn
              llb(im,ib) = il
              nl(il) = in
            endif
          endif
  50      continue
        end do
      end do
c
      do ib = 1, nb
        if ((nnb(1,ib) .eq. 0) .and. (llb(1,ib) .eq. 0)) then
          go to 103
        endif
        do il = lup, llow, -1
          nn = nnb(il,ib)
          ll = llb(il,ib)
          if (nn .ne. maxn) then
            if ((nn .ne. 0) .and. (ll .ne. 0)) then
              go to 102
            endif
          else
            nnb(il,ib) = 0
            llb(il,ib) = 0
          end if
        end do
102     continue
        numb(ib) = il
      end do
103   continue
      nb = ib - 1
cccc
cccc  debugging
cccc
      print*,' '
      do ii = 1, nb
        do il = 1, numb(ii)
          nn = nnb(il,ii)
          ll = llb(il,ii)
          if (nn .ne. maxn) then
            ind = inx(ll,nn)
            write(6,'(4i5,4f15.5)') ii,il,nn,ll,(modes(ind,jj),jj = 1,4)
          else
            write(6,'(4i5)') ii,il,nn,ll
          endif
        end do
        print*,' '
      end do
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
                if (abs(we - wc) .lt. 0.50) then
                  nnb(il, ib) = ii
                  lost(il, ii) = 0
                  go to 999                   
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
