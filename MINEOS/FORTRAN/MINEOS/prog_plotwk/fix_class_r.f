c
c
c
      subroutine fix_class_r (b_file,icolor_fill)
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
            nt(il,ib) = in
            lt(il,ib) = il
          endif
        end do
        nb = max0(ib, nb)
      end do
c
      do ib = 1, nb
        im = 0
        do il = llow, lup
          nn = nt(il,ib)
          ll = lt(il,ib)
          if (ll .ne. 0) then
            im = im + 1
c
c           believe the first few
c
            if (im .le. 2) then
              nnb(im,ib) = nn
              llb(im,ib) = ll
              wl = modes(inx(ll,nn),1)
              cl = modes(inx(ll,nn),3)
              gl = modes(inx(ll,nn),4)
              if (im .ge. 2) then
                dwl = abs(wl - modes(inx(llb(im-1,ib),nnb(im-1,ib)),1))
              endif
            else
              wc = modes(inx(ll,nn),1)
              cc = modes(inx(ll,nn),3)
              gc = modes(inx(ll,nn),4)
              dw = abs(wc - wl)
              lf = lt(il+1,ib)
              nf = nt(il+1,ib)
              wf = modes(inx(lf,nf),1)
              cf = modes(inx(lf,nf),3)
              gf = modes(inx(lf,nf),4)
c
              fl = float(ll-1)
              slope = gl*(1000./tpi)*((2.*fl+1.)/sqrt(fl*(fl + 1.)))
     &                    *(1.0/(2.*6371.))
              be = wl - slope*fl
              we = slope*(fl+2.0) + be
c              
              if (abs(wf - we) .gt. 0.85) then
                wf = we
              endif
              dw2 = abs(wl - wf)/2.
c
c             have we skipped an l?
c
              if (ll .ne. (llb(im-1,ib) + 1)) then
                do ik = llb(im-1,ib) + 1, ll - 1
                  nnb(im,ib) = maxn
                  llb(im,ib) = ik
                  wl = wl + dwl
                  im = im +  1
                end do
              endif
c
c             sorting - second-order problem - skipped mode
c
              if((abs(dw-dwl).gt.0.30).and.(abs(dw-dw2).gt.0.30))then
                nnb(im,ib) = maxn
                llb(im,ib) = ll
                do ii = nb + 1, ib + 1, -1
                  i1 = ii - 1
                  if (index(lt(il, i1),nt(il, i1)) .eq. 6) then
                    nt(il, ii) = nt(il, i1)
                    lt(il, ii) = lt(il, i1)
                  endif
                end do
                wl = wl + dwl
c
c             ok
c 
              else
                nnb(im,ib) = nn
                llb(im,ib) = ll
                wl = wc
                cl = cc
                gl = gc
                dwl = dw
              endif
            endif
          endif
        end do
        numb(ib) = im
      end do
c
c
c
      do ib = 1, nb
        if ((nnb(1,ib) .eq. 0) .and. (llb(1,ib) .eq. 0)) then
          go to 103
        endif
        do il = lup, llow, -1
          nn = nnb(il,ib)
          ll = llb(il,ib)
          if (nn .ne. maxn) then
            if ((nn .ne. 0) .or. (ll .ne. 0)) then
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
ccc
      do ib = 1, nb
        do il = 1, numb(ib)
          nn = nnb(il,ib)
          ll = llb(il,ib)
          ind = inx(ll,nn)          
          if (cfill(ind) .ne. 0) then
            nnb(il,ib) = maxn
          else
            cfill(ind) = 1
          endif
        end do
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
              if (nc .eq. nn) then
                if (lc .eq. ll) then
                   go to 998
                endif
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
                go to 1000
              endif
            end do
1000        continue
            do ii = il + 1, numb(ib)
              n2 = nnb(ii, ib)
              l2 = llb(ii, ib)
              if (n2 .ne. maxn) then
                wh = modes(inx(l2, n2), 1)
                gh = modes(inx(l2, n2), 4)
                go to 1001
              endif
            end do
            wh = 50.5
            gh = 0.0
1001        continue
c               
            do ii = nlow, nup
              if (lost(lc,ii) .gt. 0) then
                wc = modes(inx(lc,ii),1)
c               if ((wc .gt. wl) .and. (wc .lt. wh)) then
                  fl = float(l1)
                  slope = gl*(1000./tpi)*((2.*fl+1.)/sqrt(fl*(fl + 1.)))
     &                    *(1.0/(2.*6371.))
                  be = wl - slope*fl
                  we = slope*lc + be
                  if (abs(we - wc) .lt. 0.75) then
                    wf = modes(inx(lc,ii+1),1)
                    if (abs(we-wf) .lt. abs(we-wc)) then
                      nnb(il, ib) = ii + 1
c                     print*,' correlate: ', ii+1, il, wl, wc, wh
                      lost(il, ii+1) = 0
                      go to 999                   
                    else 
                      nnb(il, ib) = ii
c                     print*,' correlate: ', ii, il, wl, wc, wh
                      lost(il, ii) = 0
                      go to 999                   
                    endif
                  endif
c                endif
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
