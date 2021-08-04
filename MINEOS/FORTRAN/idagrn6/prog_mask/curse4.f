c***********************************************************************
c                 CURSE
c
c contains the modified subroutine curse 
c implements a given laterally heterogenous model,i.e
c (source-station sensitive)
c
c      subroutine curse_h2(num_modes,ista,npts,index,dat)
c
c***********************************************************************

      subroutine curse_h2(num_modes,ista,npts,index,dat)
c
c     cosine recursion subroutine - for radially symmetric earth
c
c     now let's perform a little vector multiplication
c     we have a(j,nentry) - a function of modes and stations
c     we will now loop over time points, and within each time
c     point form the cosine recurison for all the modes
c
c     the loop structure here is about as good as it gets
c        the inner loop on modes goes vector
c        the outer loop on stations goes concurrent
c        the time loop has a data dependency - and we can't do anything about it
c
c     this is an expansion of
c     (cos(wt) + (1/2Q)*sin(wt))*exp(-alpha*t)
c     where each time point is expanded by t(n+1) = t(n) + dt
c     the term exp(-alpha*t(n)) is neglected
c
c     stolen directly from jsrc
c
cc     pfi 1/21/92 changed common /mode/ to /mode_h/ to allow for station by
cc           station recursions necessary to deal with laterally 
cc           heterogeneous models
cc           variables dr,di,ct depend now on maxstat, 
cc           loop structure essentially not changed
c
      implicit real*8 (a-h,o-z)
c
      include 'parameter5.h'
c
      integer*4 nptsar(maxstat)
      real*4 digar(maxstat)

      real*8 ct(maxmodes), a
      real*8 dat(maxtime,maxcomp), qinv, ddr, ddi
      real*8 dr(maxmodes), di(maxmodes), lsin, lcos
      real*8 tst,temp
      real*8 tst1,tst2,tst3,tst4,tst5,tst6
c
      integer*4 num_modes, npts, index
c
      common/mode_h/ddr(maxmodes),ddi(maxmodes),
     &                  qinv(maxmodes)
      common/aa/ a(maxmodes,maxcomp)
      common /stasen/digar,nptsar
c
      print*,' calculating recursion for ',num_modes,' modes'
      if (npts.ne.0) then
      print*,'                       for ',npts,' time points'
      else
      print*,'                       for variable time points'
      endif
      print*,'                   and for ',index,' components'

cc...initialize the recursion variables

      do j = 1, num_modes
        di(j) = 0.d0
        dr(j) = 1.d0
        ct(j) = 1.d0
      enddo

cc...test for index here

      if (index .gt. 1) then

cc...do the first time point

        do k = 1, index
          tst = 0.d0
          do j = 1, num_modes
            tst = tst + a(j,k)
          end do
          dat(1,k) = tst
        end do
c
          npts = nptsar(ista)
          do l = 2, npts
            tst1 = 0.d0
            tst2 = 0.d0
            tst3 = 0.d0
            tst4 = 0.d0
            tst5 = 0.d0
            tst6 = 0.d0
            do j = 1, num_modes
              lcos = dr(j)
              lsin = di(j)
              dr(j) = -lsin * ddi(j) + lcos * ddr(j)
              di(j) =  lsin * ddr(j) + lcos * ddi(j)
              temp    =  dr(j) + di(j) * qinv(j)
              tst1 = tst1 + a(j,1)*temp
              tst2 = tst2 + a(j,2)*temp
              tst3 = tst3 + a(j,3)*temp
              tst4 = tst4 + a(j,4)*temp
              tst5 = tst5 + a(j,5)*temp
              tst6 = tst6 + a(j,6)*temp
            enddo
            dat(l, 1) = tst1
            dat(l, 2) = tst2
            dat(l, 3) = tst3
            dat(l, 4) = tst4
            dat(l, 5) = tst5
            dat(l, 6) = tst6
          enddo

      else
c
c       do the first time point
c
        tst = 0.d0
        do j = 1, num_modes
          tst = tst + a(j,1)
        end do
        dat(1, 1) = tst
c
          npts = nptsar(ista)
          do l = 2, npts
            tst = 0.d0
            do j = 1, num_modes
              lcos = dr(j)
              lsin = di(j)
              dr(j) = -lsin * ddi(j) + lcos * ddr(j)
              di(j) =  lsin * ddr(j) + lcos * ddi(j)
              tst = tst + a(j,1)*(dr(j)+di(j)*qinv(j))
            enddo
            dat(l, 1) = tst
          enddo

      endif   
c
c     all done
c
      return
      end
