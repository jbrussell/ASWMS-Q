      subroutine branch_sort(b_file)
c
c     subroutine to sort branches by nn
c
      include 'parameter.h'
c
      real*4 modes(0:maxmodes,nprop), wmin, wmax
      real*4 fill(0:maxmodes), temp(4)
      integer*4 jcom, nmodes
      integer*4 ifill(0:maxmodes),cfill(0:maxmodes)
      integer*4 inx(0:maxl,0:maxn),inx2(0:maxmodes,2)
      integer*4 nnb(0:maxl,nbranch), llb(0:maxl,nbranch), numb(nbranch)
      character*256 b_file
      logical lsummary, lbranch
c
      common /mode/ inx, modes, inx2
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /limits2/ wmin, wmax, pmin, pmax, gmin, gmax
      common /fill/ lsummary, fill, ifill, cfill
      common /branch/ lbranch, ib, numb, nnb, llb
c
      data lbranch /.false./
c
c     first, zero any previous branch info
c 
      do ii = 1, nbranch
        do jj = 0, maxl
          nnb(jj,ii) = 0
          llb(jj,ii) = 0
        end do
      end do
      ib = 0
c
      lbranch = .true.
c
c     sort these branches by n
c
      ib = 0
      do ii = nlow, nup
        ib = ib + 1
        im = 0
        do jj = llow, lup
          ic1 = inx(jj,ii)
          if ((modes(ic1,1) .ne. 0.) .and. (cfill(ic1) .gt. 3)) then
            ifill(ic1) = 1
            cfill(ic1) = 3 + (ib - (ib/3)*3)
            im = im + 1
            nnb(im,ib) = ii
            llb(im,ib) = jj            
          endif
        end do
        numb(ib) = im
      end do
c
c     now write output file containing branch information
c       two header records, followed by one record for each branch
c       each branch record contains a 6-tuple of info about each mode - nn, ll, w, q, cv, gv
c
      open(unit=4,file=b_file,form='unformatted',access='sequential')
      write(4) jcom, nmodes, nlow, nup, llow, lup
      write(4) ib, (numb(ii), ii = 1, ib)
      do ii = 1, ib
        write(4) (nnb(jj,ii), llb(jj,ii), 
     &    (modes(inx(llb(jj,ii),nnb(jj,ii)),kk),kk=1,4),jj=1,numb(ii))
      end do
      close(4)
      return
      end
