      subroutine fix_class_v (b_file,icolor_fill)
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
      co = (1000./tpi)*(1.0/(2.*6371.))
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
      ib = 0
      do in = nlow, nup
        if (index(1,in) .eq. 6) then
          ib = ib + 1
        end if
      end do   
      nb = ib
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
cc      do ib = 1, nb
cc        nnb(1,ib) = nt(llow,ib)
cc        llb(1,ib) = lt(llow,ib)
cc      end do 
cc      do il = llow + 1, lup
cc        do ib = 1, nb
cc          do ii = il - 1, 1, -1
cc            n1 = nnb(ii, ib)
cc            l1 = llb(ii, ib)
cc            if (n1 .ne. maxn) then
cc              wl = modes(inx(l1, n1), 1)
cc              gl = modes(inx(l1, n1), 4)
cc              fl = float(l1)
cc              slope = gl*co*((2.*fl+1.)/sqrt(fl*(fl + 1.)))
cc              be = wl - slope*fl
cc              we = slope*float(il) + be
cc              go to 101
cc            endif
cc          end do
cc101       continue
cc          dwa = 99999.
cc          do jj = 1, nb
cc            n2 = nt(il,jj)
cc            l2 = lt(il,jj)
cc            wc = modes(inx(l2,n2),1)
cc            dw = abs(wc - we)
cc            if (dw .lt. dwa) then
cc              dwa = dw
cc              jsave = jj
cc            endif
cc          end do
cc          if (dwa .lt. 0.25) then
cc            nnb(il,ib) = nt(il,jsave)
cc            llb(il,ib) = lt(il,jsave)
cc            nt(il, jsave) = 0
cc            lt(il, jsave) = 0
cc          else
cc            nnb(il,ib) = maxn        
cc            llb(il,ib) = il
cc          endif
cc        end do
cc      end do
ccc
cc      do ib = 1, nb
cc        do il = lup, llow, -1
cc          if (nnb(il,ib) .ne. maxn) then
cc            go to 102
cc          else
cc            nnb(il,ib) = 0
cc            llb(il,ib) = 0
cc          end if
cc        end do
cc102     continue
cc        numb(ib) = il
cc      end do
cccc
      do ib = nb, 1, -1
        im = 0
        do il = llow, lup
          nn = nt(il,ib)
          ll = lt(il,ib)
          if (ll .ne. 0) then
            im = im + 1
c
c           believe the first one
c
            if (im .eq. 1) then
              nnb(im,ib) = nn
              llb(im,ib) = ll
              wl = modes(inx(ll,nn),1)
              cl = modes(inx(ll,nn),3)
              fl = float(ll)
              gl = modes(inx(ll,nn),4)
              dwl = gl*co*((2.*fl+1.)/sqrt(fl*(fl + 1.)))
            else
              wc = modes(inx(ll,nn),1)
              cc = modes(inx(ll,nn),3)
              gc = modes(inx(ll,nn),4)
              dw = abs(wc - wl)
              lf = lt(il+1,ib)
              nf = nt(il+1,ib)
              wf = modes(inx(lf,nf),1)
c
              fl = float(ll-1)
              slope = gl*co*((2.*fl+1.)/sqrt(fl*(fl + 1.)))
              be = wl - slope*fl
              we = slope*(fl+2.0) + be
              if (abs(wf - we) .gt. 0.5) then
                wf = we
              endif
              dw2 = abs(wl - wf)/2.
c
c             have we skipped an l?
c
              if (ll .ne. (llb(im-1,ib) + 1)) then
                nnb(im,ib) = maxn
                llb(im,ib) = llb(im-1,ib) + 1
                wl = wl + dwl
                im = im +  1
              endif
c
c             sorting 
c
              if (wc .lt. wl) then
                nnb(im,ib) = maxn
                llb(im,ib) = llb(im-1,ib) + 1
                wl = wl + dwl
                do ii = 1, ib - 1
                  i1 = ii + 1
                  nt(il, ii) = nt(il, i1)
                  lt(il, ii) = lt(il, i1)
                end do
              elseif ((abs(dw-dwl).gt..25).and.(abs(dw-dw2).gt..25))then
                if (ib .eq. 1) then
                  nnb(im,ib) = maxn
                  llb(im,ib) = llb(im-1,ib) + 1
                  wl = wl + dwl
                else 
                  nn = nt(il,ib-1)
                  ll = lt(il,ib-1)
                  wc = modes(inx(ll,nn),1)
                  cc = modes(inx(ll,nn),3)
                  gc = modes(inx(ll,nn),4)
                  dw = abs(wc - wl)
                  if (wc .lt. wl) then
                    nnb(im,ib) = maxn
                    llb(im,ib) = llb(im-1,ib) + 1
                    wl = wl + dwl
                  elseif((abs(dw-dwl).gt..25).and.
     &                   (abs(dw-dw2).gt..25))then
                    nnb(im,ib) = maxn
                    llb(im,ib) = llb(im-1,ib) + 1
                    wl = wl + dwl
                  else
                    nnb(im,ib) = nn
                    llb(im,ib) = ll
                    do ii = ib, 2, -1
                       nt(il, ii) = nt(il, ii - 1)
                       lt(il, ii) = lt(il, ii - 1)
                    end do
                    nt(il, 1) = 0
                    lt(il, 1) = 0
                    wl = wc
                    cl = cc
                    gl = gc
                    dwl = dw
                  endif
                endif 
              else
c
c               check the next branch down for a better fit
c
                if (ib .ne. 1) then
                  n1 = nt(il,ib-1)
                  l1 = lt(il,ib-1)
                  wn = modes(inx(l1,n1),1)
                  dwn = abs(wn - wl)
                  if (wn .lt. wl) then
                    nnb(im,ib) = nn
                    llb(im,ib) = ll
                    wl = wc
                    cl = cc
                    gl = gc
                    dwl = dw
                  elseif((abs(dwn-dwl).gt..25).and.
     &                   (abs(dwn-dw2).gt..25))then
                    nnb(im,ib) = nn
                    llb(im,ib) = ll
                    wl = wc
                    cl = cc
                    gl = gc
                    dwl = dw
                  else
                    if (abs(dwn-dwl) .lt. abs(dwc-dwl)) then
                      nnb(im,ib) = n1
                      llb(im,ib) = l1
                      wl = wn
                      cl = modes(inx(l1,n1),3)
                      gl = modes(inx(l1,n1),4)
                      dwl = dwn
                      do ii = ib, 2, -1
                         nt(il, ii) = nt(il, ii - 1)
                         lt(il, ii) = lt(il, ii - 1)
                      end do
                      nt(il, 1) = 0
                      lt(il, 1) = 0
                    else
                      nnb(im,ib) = nn
                      llb(im,ib) = ll
                      wl = wc
                      cl = cc
                      gl = gc
                      dwl = dw
                    endif
                  endif
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
c
c           check to make sure that we aren't going off into infinity

            fl = llb(im,ib)
            slope = gl*co*((2.*fl+1.)/sqrt(fl*(fl + 1.)))
            be = wl - slope*fl
            we = slope*(fl+1.0) + be
            if (we .ge. 50.0) then
              do ii = il + 1, lup
                do in = 1, ib - 1
                  nt(ii,in) = nt(ii,in+1)
                  lt(ii,in) = lt(ii,in+1)
                end do
              end do
              go to 888
            endif
          endif
        end do
888     continue
        numb(ib) = im
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
      do ib = nb, 1, -1
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
            print*,' checking ', ib, il, wl, wh
            do ii = nup, nlow, -1
              if (lost(lc,ii) .gt. 0) then
                wc = modes(inx(lc,ii),1)
                if ((wc .gt. wl) .and. (wc .lt. wh)) then
                  dw = abs(we - wc)
                  if (dw .lt. slope) then
                    nnb(il, ib) = ii
                    print*,' correlate: ', ii, il, wl, wc, wh, dw, slope
                    lost(il, ii) = 0
                    go to 999                   
                  else
                    print*,' not quite: ', we, wc, dw, slope
                  endif
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
      print*,' '
      do ii = 1, nmodes
        nn = inx2(ii,2)
        ll = inx2(ii,1)
        if (lost(ll,nn) .gt. 0) then
          print*,' lost mode: ', nn, ll, modes(ii,1)
        endif
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
      do ii = 1, nmodes
        nn = inx2(ii,2)
        ll = inx2(ii,1)
        if (lost(ll,nn) .gt. 0) then
          cfill(ii) = 9
        endif
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
