      subroutine mtimes(depth,delta,plist,nnip,ttip,ppip,dtdhh,
     & dddpp,ppic)
c
c     subroutine for calculating multiple travel times
c       using IASPEI 91 tables
c     very crude attempt to deal with phases like SSS and ScS4
c
c     calls ttimes
c
c     depth - event depth in km
c     delta - epicentral distance in degrees
c     plist - character string describing phases
c     nnip - number of travel times returned
c     ttip - array of nnip travel times
c     ppip - array of nnip ray parameters - d time /d distance values
c     dtdhh - array of nnip d time /d depth values
c     dddpp - array of nnip d distane /d ray parameter values
c     ppic - array of nnip phase names
c
      include 'parameter.f'
c
      real*4 depth, delta
      real*4 ppip(nphase), ttip(nphase), tip(nphase), pip(nphase)
      real*4 pips(nphase), tips(nphase), dddps(nphase)
      real*4 dtdh(nphase), dddp(nphase), slow(2), sls(2), slr(2)
      real*4 dtdhh(nphase), dddpp(nphase)
      logical lmatch
      character*1 ans
      character*8 ppic(nphase), phip(nphase), phips(nphase)
      character*8 plist, slist, phcsd, phcs0
c
      common /t_slow/ slow
      common /s_slow/ sls, slr
c
      data xmatch /0.0000025/
      data maxmult /6/
c
      nnip = 0
      iz = ichar(plist(2:2))
      call klen(plist,kl)
      plist = plist(2:kl)
c
c     downward going phases
c
      if (iz .ge. 65 .and. iz .lt. 97) then
c
c     find the 0th multiple - such as ScS or PcP
c
c       if multiple arrivals, find the right one
c
        call ttimes(depth,delta,plist,nip,tip,pip,dtdh,dddp,phip)
c
c       save the source slowness
c
        do ii = 1, 2
          sls(ii) = slow(ii)
        end do
c
        do is = 1, nip
          if (phip(is) .eq. plist) then
            nnip = nnip  + 1
            ttip(nnip) = tip(is)
            ppic(nnip) = plist
            ppip(nnip) = pip(is)
            dtdhh(nnip) = dtdh(is)
            dddpp(nnip) = dddp(is)
          endif
        end do
c
c       now look for multiples - up to maxmult
c
        xd = 0.0
        zdep = 0.0
        do kk = 2, maxmult
c
c       figure about approximate distances
c
          zdeltad = delta/float(kk)
          zdelta0 = (delta - zdeltad)/float(kk - 1)
c
c         there is a a problem here with multiple distances
c           which needs to be examined in greater detail
c           in order to do S2,3,4, and P2,3,4, successfully.
c
          iloop = 0
4810      continue
          iloop = iloop + 1
c
c         check for a lost cause
c
          if (iloop .gt. 200) then
            print*,' exiting travel time loop ',plist,' ',kk
            go to 4820
          endif
c
c         calculate the 'depth' ray
c
          call ttimes(depth,zdeltad,plist,nip,tip,pip,dtdh,dddp,phip)
          do is = 1, nip
            if (phip(is) .eq. plist) then
              tips(is) = tip(is)
              pips(is) = pip(is)
              phips(is) = phip(is)
              dddps(is) = dddp(is)
            end if
          end do
          nips = nip
c
c         calculate the 'surface' ray
c
          call ttimes(zdep,zdelta0,plist,nip,tip,pip,dtdh,dddp,phip)
c
c         save the receiver slowness
c
          if (iloop .eq. 1) then
            do ii = 1, 2
              slr(ii) = slow(ii)
            end do
          endif
c
c         now begin double loop - look for a match - any match
c
          lmatch = .false.
          xdm = 10000.0
          do is = 1, nips
            tscsd = tips(is)
            pscsd = pips(is)
            phcsd = phips(is)
            dddpsd = dddps(is)
            do id = 1, nip
              tscs0 = tip(id)
              pscs0 = pip(id)
              phcs0 = phip(id)
              dddps0 = dddp(id)
c
c             attempt to avoid serious mistakes
c
              if (phcsd .eq. phcs0) then
c
c               try to match them up by matching ray parameter
c
                xd = abs(pscsd - pscs0)
                if (dddpsd*dddps0 .ge. 0) then
c
c                 check that we are looking at the same type of phase
c
                  if (xd .le. xmatch) then
c
c                   ray parameter match !
c
                    nnip = nnip  + 1
                    ttip(nnip) = tscsd + tscs0*float(kk-1)
                    ppip(nnip) = pip(id)
                    dtdhh(nnip) = dtdh(id)
                    dddpp(nnip) = dddp(id)*(delta/zdeltad)
                    write(ans, '(i1)') kk
                    call kblnk(plist,kl)
                    ppic(nnip) = plist(1:kl)//ans
                    lmatch = .true.
                  else
                    if (xd .lt. xdm) then
                      xdm = xd
                      pscsdm = pscsd
                      pscs0m = pscs0
                      dddpdm = dddpsd
                      dddp0m = dddps0
                    end if
                  end if
                end if
              end if
            end do
          end do
          if (.not.(lmatch)) then
c
c           close, but not quite
c
c           this is correct for a back branch
c
            if (dddpdm .gt. 0.0) then
              if (pscsdm .lt. pscs0m) then
                zdeltad = zdeltad  + 2.0*xdm
              else
                zdeltad = zdeltad  - 2.0*xdm
              endif
            else
c
c           this should be correct for a forward branch
c
              if (pscsdm .lt. pscs0m) then
                zdeltad = zdeltad  - 2.0*xdm
              else
                zdeltad = zdeltad  + 2.0*xdm
              endif
            endif
            zdelta0 = (delta - zdeltad)/float(kk - 1)
            if (zdelta0 .lt. 0) then
              print*,' error in iterating to solution'
              go to 4820
            end if
c
c           write(6,'(i4,7f9.5)')
c    &        kk,zdeltad,zdelta0,xdm,pscsdm,pscs0m,dddpdm,dddp0m
c
            go to 4810
          end if
4820      continue
        end do
      elseif (iz .ge. 97 .and. iz .lt. 123) then
c
c       no sScS or pPcP in the tables - they should be shot!
c       I will need to do something clever here
c
        print*,' upgoing phase search not implemented'
c
c       form a table of upward going phases
c       --this must be done in two steps
c         
c       slist = plist(1:2)
c       zdep = 0.0
c
c       start from the current distance and work back
c
c       do xd = delta, 2.0, -2.0
c
c       calculate the 'sS' or the 'pP' or whatever time
c  
c          call ttimes(depth,xd,slist,nip,tip,pip,dtdh,dddp,phip)
c          do is = 1, nip
c            if (phip(is) .eq. slist) then
c              tscsd = tip(is)
c              pscsd = pip(is)
c              go to 4815
c            endif
c          end do
c4815      continue
c
c         now try to find the 0 depth S or P wave which corresponds 
c
c          zdelta = 
c4820      continue
c
      endif
c
c     I guess we can return now
c
      return
      end
