      subroutine avg_dbdb5(elat,elon,delta,azi,dbsa,dbgc)
c given source and receiver location, calculates an average
c value of bathymetry using dbdb5 and the libutil call get_bath
c
c input:  dbdb5 via get_bath (location of data file is hardwired there)
c         elat,elon -- event location
c         dist,azi -- distance and azimuth from event to station

c
c outputs:  For both short-arc AND total great-circle
c             avg bathym/topo 
c
c JBG 9/99
c        


      parameter(maxdel = 360)
      real*4 delt(maxdel*24),gcloc(2,maxdel*24)
      real*4 elat,elon,delta,azi,dbsa,dbgc



c     common block of crustal params includes bathym values to return

      real*4 vvpsa,vvssa,vrhosa,vtopsa,vmohsa,
     &        vvpgc,vvsgc,vrhogc,vtopgc,vmohgc
      common /avg51/ vvpsa,vvssa,vrhosa,vtopsa,vmohsa,
     &                vvpgc,vvsgc,vrhogc,vtopgc,vmohgc
  
      include 'numerical.h'


      vtopgc = 0.

c     nsample is the number of samples per degree -- dbdb5 is defined at 12/deg
c        print*,elat,elon,delta
        nsample = 24
        ndelgc = 360
        ndelgcb = ndelgc*nsample
        ndelsab = int(delta*nsample)
        do ii = 1,ndelgcb
          delt(ii) = float(ii)*dkm/float(nsample)
        end do

        call gcpath(elat,elon,azi,delt,ndelgcb,gcloc)

        dbgc = 0

        do ii = 1,ndelsab
         call get_bath(gcloc(1,ii),gcloc(2,ii),db)
c          print*,gcloc(1,ii),gcloc(2,ii),delt(ii),db
         dbgc = dbgc + db
        enddo

        dbsa = dbgc / float(ndelsab)

        do ii = ndelsab+1,ndelgcb
           call get_bath(gcloc(1,ii),gcloc(2,ii),db)
         dbgc = dbgc + db
        enddo

        dbgc = dbgc / float(ndelgcb)

c        print*,'In avg_dbdb5, dbsa, dbgc: ',dbsa,dbgc

 99   continue
   
      return

      end
