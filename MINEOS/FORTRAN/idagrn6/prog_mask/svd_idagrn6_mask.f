      program idagrn6_sac

c idagrn6 that outputs sac files
c based on idagrn6.f
c how to build: see Makefile
c usage: edit idagrn.in in current directory then run idagrn6_sac

c Po Chen
c 2006-07-20

      implicit real*8 (a-h,o-z)
c
      include 'parameter5.h'
c     
      real*8 dat(maxtime,maxcomp)
c     
      real*4 data(maxtime), rnum, rdump
      real*4 dt, hours
      real*4 wref
      real*4 del(maxstat), azm(maxstat), dbgc,dbsa
      real*4 utime1(2), utime2(2),etime
c     
      character*1 lcomp
      character*10 sname(maxstat), stname(maxstat), comps
      character*4 compt(6)
      character*12 modeid, cmodel
      character*80 filep, fileq
      character*80 m_file, b_file, adump
c
      integer*4 num(maxstat)
      integer*4 npts, num_modes, idump
      integer*4 evt(5)
      integer*4 index2, nbr
      integer*4 hnum
c************************************************************
c declaration for laterally heterogeneous model
c
      parameter(hnum=5)

      integer*4 ind
      real*8 slat(maxstat),slon(maxstat),sele(maxstat)
      character*4 hmod(hnum)
      character*80 hname,h_file

c declaration for angular parameters and averages over complex
c spherical harmonics

      real*8      dist(maxstat),azmp(maxstat),theta(maxstat)
      complex*16  avg_g(-nord:nord,0:nord),avg_d(-nord:nord,0:nord)
      complex*16  avg_great(-nord:nord,0:nord,maxstat)
      complex*16  avg_diff(-nord:nord,0:nord,maxstat)

      common /ang_h/avg_great,avg_diff,dist,azmp,theta

c declaration for radial integrals and ellipticity correction factor
c spheroidal as well as toroidal (part of it)

      integer*4    jcom_h,ifanis_h,imode
      real*4       tref_h,rn_h
      character*80 model_h

      common /integ_const_h/model_h,tref_h,rn_h,jcom_h,ifanis_h,imode
c
c station sensitive (in source_strain and curse)

      integer*4 nptsar(maxstat)
      real*4    digar(maxstat)
      logical   lstasen,lequal
    
      common /stasen/digar,nptsar

c common block for crust/dbdb options

      logical lcrust51, ldbdb5, ldb1st, lcr1st
      character c51opt*2,cdbopt*2,c_file*80,d_file*80

      common /coption/ lcrust51,ldbdb5,c51opt,cdbopt
       
c common block of crust5.1 model parameters


      real*4 amapvp(8,72,36),amaprho(8,72,36),
     +          amapvs(8,72,36),amapthi(7,72,36),
     +          amapele(72,36)
      common /crust51/ amapvp,amaprho,amapvs,amapthi,
     +                 amapele

c common block of crust5.1 calculated average values

      real*4 vtopsa,vmohsa,vbassa
      real*4 vvp1sa,vvs1sa,vrho1sa,vvp2sa,vvs2sa,vrho2sa
      real*4 vtopgc,vmohgc,vbasgc
      real*4 vvp1gc,vvs1gc,vrho1gc,vvp2gc,vvs2gc,vrho2gc
      common /avg51/ vvp1sa,vvs1sa,vrho1sa,vvp2sa,vvs2sa,vrho2sa,
     &               vtopsa,vbassa,vmohsa,vvp1gc,vvs1gc,vrho1gc,
     &                vvp2gc,vvs2gc,vrho2gc,vtopgc,vbasgc,vmohgc
  

c common block for crust5.1/dbdb5 kernel information

      real*4 crker(nlith,6,maxmodes),crkers(maxmodes),
     &       crkerb(maxmodes),crkerm(maxmodes), gvi(maxmodes)
      real*4 vs(nlith),vp(nlith),ro(nlith),rcr(nlith)
      real*4 rsol, rmoho, rbase
      integer*4 ifanis,kcrtot,ipnt(maxn,maxl)
      common /kern51/ ifanis,crker,crkers,crkerb,crkerm,gvi,ipnt
      common /omod51/ rsol,rbase,rmoho,vs,vp,ro,rcr,kcrtot

c common block for station by station crustal averages cbasegc

      real*4 cvp1gc(maxstat),cvs1gc(maxstat),crho1gc(maxstat)
      real*4 cvp2gc(maxstat),cvs2gc(maxstat),crho2gc(maxstat)
      real*4 ctopgc(maxstat),cbasgc(maxstat),cmohgc(maxstat)
      real*4 cvp1la(maxstat),cvs1la(maxstat),crho1la(maxstat)
      real*4 cvp2la(maxstat),cvs2la(maxstat),crho2la(maxstat)
      real*4 ctopla(maxstat),cbasla(maxstat),cmohla(maxstat)
      common /savg51/ cvp1gc,cvs1gc,crho1gc,cvp2gc,cvs2gc,
     &  crho2gc,ctopgc,cbasgc,cmohgc,cvp1la,
     &  cvs1la,crho1la,cvp2la,cvs2la,crho2la,
     &  ctopla,cbasla,cmohla

c************************************************************

      logical gfs, table, lstat, lmask, lover, lexcite,lread
      logical ldisp,lida,lsource,lbranch,lhetero,lstrain,lreceiv

c
c     sac header commons
c
      integer*4 num_pts, data_time(5), event_time(5)
      integer*4 djulian, ejulian, sum, nerr
c
      real*4 dig, elat, elon, rlat, rlon, relev, dsec, esec
      real*4 m0, moment(6), d0, delta, azim, bazim, dep
      real*4 hdr(200), hdum1, hdum2, hdum3, hdum4
c
      character statn*4, stype*4, comp*4, catch*8, model*12, sourid*24,
     & synid*12, comment*80, sacfnam*80
c
      common /head/ num_pts, statn, stype, comp, catch, data_time,
     &               dsec, djulian, dig, elat, elon, d0, hdum1(2),
     &               rlat, rlon, relev, event_time, esec, ejulian,
     &               hdum2(63), moment, hdum3(10), model, sourid, 
     &               synid, hdum4(60), comment
c
      equivalence (hdr, num_pts)
c
      common/ahdr/ jcom,knots
      common/limits/ dt,wref,lcomp,ldisp,lida
      common/mask/ m_file, lmask
      common/branch/ b_file, lbranch, nbr
      common/exc/ evt, filep, modeid, stname, comps, dep
c
c     command processor declarations
c
      character*256 typbuf, token, sub
      character*80 source_file(3)
      character*256 cmd
c
      data rad/57.29578d0/,icomp/0/
      data rn /6371000.d0/
      data in /8/
      data iout /9/
      data itype /1/
      data hours /1./
      data compt /'zgrr', 'zgtt', 'zgpp', 'zgrt', 'zgrp', 'zgtp'/
      data typbuf, token, sub /3*' '/
      data fileq, source_file /4*' '/
      data lstat, nstat /.false., 0/
      data gfs, table /2*.false./
      data lover, lexcite /2*.false./
      data lsource /.false./
      data lstasen,lequal /.false.,.true./
      data lhetero,lstrain,lreceiv /3*.false./
      data lcr1st,ldb1st /2*.false./
      data lread /.false./
      data hname /' '/
      data hmod /'M84A','M84C','SH8W','SH8U','SH12'/ 
      data cmodel /'prem'/

c     variables in common blocks cannot be initialized with data statements utime1
      wref=4.0
      dt=1.
      lcomp='Z'
      filep=' '
      modeid='mode syn'
c      lmask=.false.
      lida=.false.
      ldisp=.false.
      lcrust51=.false.
      ldbdb5=.false.

c     read the input file for idagrn
      open(11,file='idagrn.in',form='formatted')
      read(11,'(a)') source_file(2) 
      table = .true.
      read(11,*) hours, dt         
c      print*, 'hours, dt', hours, dt
      if (dt.le.0.0) lstasen = .true.
      read(11,*) wref               
      read(11,'(a)') cmodel        
      read(11,'(a)') lcomp          
      read(11,*) ldisp          
c      read(11,'(a)') modeid    

      read(11,*) lcrust51
      if(lcrust51) then
         read(11,'(a)') c51opt 
         read(11,'(a)') c_file 
         lcr1st=.true.
      else
         read(11,'(a)') adump
         read(11,'(a)') adump
      endif

      read(11,'(a)') sourid
      read(11,*) elat, elon, d0 
      read(11,*) m0
      do ii=1,6
         read(11,*) moment(ii)
         moment(ii)=moment(ii)*m0
      enddo

      read(11,*) nstat    
      if (nstat .gt. maxstat) then
        print*,' parameter maxstat exceeded'
        stop
      endif
      do ii=1,nstat
         read(11,*) stname(ii), slat(ii), 
     &        slon(ii), sele(ii)
         print*, stname(ii), slat(ii), slon(ii), sele(ii)
      enddo

      lstat=.true.

c      section below was uncommented to allow masking of branches

      read(11,*) lmask
      if(lmask) then
      	 print*,'Looking for mask file'
	 read(11,'(a)') m_file  
         read(11,*) idump
         lbranch=.false.
      else
         read(11,'(a)') b_file  
	 print*,'Branch File',b_file
         read(11,*) nbr
         lbranch=.true.
      endif

c      read(11,*) lexcite 

c      read(11,*) lida    
c      read(11,*) lsource 

c      read(11,*) lhetero 
c      if(lhetero) then
c         read(11,'(a)') hanme
c         read(11,'(a)') h_file
c         cmodel(1:4)='prem'
c      else
c         read(11,'(a)') adump
c         read(11,'(a)') adump
c      endif
      

c      read(11,*) ldbdb5
c      if(ldbdb5) then
c         read(11,'(a)') cdbopt 
c         read(11,'(a)') d_file
c         ldb1st=.true.
c      else
c         read(11,'(a)') adump
c         read(11,'(a)') adump
c      endif
      
c      read(11,'(a)') source_file(1) 
c      lstrain=.true.
c      read(11,'(a)') source_file(3)
c      lreceiv=.true.
      
      close(11)
c      print*, 'end reading input file idagrn.in'
c     end reading the input file for idagrn

      len1 = lnblnk(filep)
      len2 = lnblnk(fileq)
      len3 = lnblnk(m_file)
      len4 = lnblnk(b_file)
      len5 = lnblnk(source_file(1))
      len6 = lnblnk(source_file(2))
      len7 = lnblnk(source_file(3))

      if (table.and..not.lstrain.and..not.lreceiv) then
         print*,' Mode table s_p_r  :  ', source_file(2)(1:len6)
      elseif (lstrain.and..not.table.and..not.lreceiv) then
         print*,' Mode table s_p_r  :  ', source_file(1)(1:len5)
      elseif (lreceiv.and..not.table.and..not.lstrain) then
         print*,' Mode table s_p_r  :  ', source_file(3)(1:len7)
      elseif (table.and.lstrain.and.lreceiv) then
         print*,' Mode table s_     :  ', source_file(1)(1:len5)
         print*,' Mode table   p_   :  ', source_file(2)(1:len6)
         print*,' Mode table     r  :  ', source_file(3)(1:len7)
      elseif (table.and.lstrain) then
         print*,' Mode table   p_r  :  ', source_file(2)(1:len6)
         print*,' Mode table s_     :  ', source_file(1)(1:len5)
      elseif (table.and.lreceiv) then
         print*,' Mode table s_p    :  ', source_file(2)(1:len6)
         print*,' Mode table     r  :  ', source_file(3)(1:len7)
      elseif (lstrain.and.lreceiv) then
         print*,' Mode table s_p_   :  ', source_file(1)(1:len5)
         print*,' Mode table     r  :  ', source_file(3)(1:len7)
      endif    
      print*,' Radial Earth model:  ', cmodel
      if (lhetero) then
         print*,' Aspherical model  :  ', hname(1:4)
      endif                  
      print*,' Component type    :  ', lcomp,'  Override: ',lover
      if (lstasen) then
         print*,' Length(hours), dt :  ', hours,' ',dt,' <- from header'
      else
         print*,' Length(hours), dt :  ', hours,' ',dt,' <- fixed'
      endif
      print*,' Ref. frequency    :  ', wref
      print*,' Catch             :  ', catch
      if (lmask) then
         print*,' Masking           :  ', m_file(1:len3)
      endif
      if (lbranch) then
         print*,' Branches          :  ', b_file(1:len4)
         if (nbr .ne. 0) then
            print*,' Branch number     :  ',nbr
         endif
      endif
      if (ldisp) then
         print*,' Displacement records  '
      else
         print*,' Acceleration records  '
      endif
      print*,' '
      if (lsource) then
         print*,' Source problem : 6 Greens functions'
      else
         print*,' Structure problem : Single seismogram'
      endif
      print*,' '
      if (lida) then
         print*,'  Z-down  R-away from the source  T-westerly'
      else
         print*,'  Z-up    R-toward the source     T-easterly'
      endif
      print*,' '
      if (lstat) then
         print*,' Synthetics calculated for these stations:'
         do ii = 1, nstat
c            print*,stname(ii)
         end do
      else
         print*,' Synthetics calculated for every station'
      endif
      

cc
cc  load name of source_file correctly
cc  if the strain_mode_table has not been specified
cc
      if (table.and..not.lstrain.and..not.lreceiv) then
         source_file(1) = source_file(2)
         source_file(3) = source_file(2)
      elseif (lstrain.and..not.table.and..not.lreceiv) then
         source_file(2) = source_file(1)
         source_file(3) = source_file(1)
      elseif (lreceiv.and..not.table.and..not.lstrain) then
         source_file(1) = source_file(3)
         source_file(2) = source_file(3)
      elseif (table.and.lstrain.and..not.lreceiv) then
         source_file(3) = source_file(2)
      elseif (table.and.lreceiv.and..not.lstrain) then
         source_file(1) = source_file(2)
      elseif (lstrain.and.lreceiv.and..not.table) then
         source_file(2) = source_file(1)
      endif
      len5 = lnblnk(source_file(1))
      len6 = lnblnk(source_file(2))
      len7 = lnblnk(source_file(3))
      print*,' Mode table s_     :  ', source_file(1)(1:len5)
      print*,' Mode table   p_   :  ', source_file(2)(1:len6)
      print*,' Mode table     r  :  ', source_file(3)(1:len7)
      
cc
cc  load big kernel and database files for crustal and/or dbdb5
cc  corrections.  Only need to do it for first comp/stat for
cc  each mode type.  dbdb5 (if done separately from crust5.1)
cc  requires kernels, but not the datafile.

c      if(lcrust51.and.lcr1st) then
c         call load_kern_c51(c_file)
c         call load_crust51
c      elseif(ldbdb5.and.ldb1st) then
c         call load_kern_c51(c_file)
c      endif
      

cc------------------------------------------------------------------
cc   PREPARE FOR SUMMATION, READ HEADER INFO
cc   
cc------------------------------------------------------------------
      nentry=0
      do is=1,nstat
         nentry = nentry + 1
         num(nentry) = is
         call azimth(elat,elon,real(slat(is)),real(slon(is))
     &        ,delta,azim,bazim)
         
         del(nentry) = delta
         azm(nentry) = azim
c         print*, is, elat, slat(is), elon, slon(is), del(is)
      enddo
c      stop
cc...set index2 for source or structure problem

      if (lsource)      index2 = 6
      if (.not.lsource) index2 = 1

      if (.not.lstasen) then
         npts = int(hours*3600./dt)
         if (npts .gt. maxtime) then
            print*,' array dimension maxtime exceeded'
            stop
         endif
         do i = 1, nentry
           digar(i) = dt
           nptsar(i) = npts
         enddo
      endif

cc...implement aspherical model in source_strains
cc
cc   compute first a number of angular parameter and averages over 
cc   complex spherical harmonics and load  common /ang_h/
cc   and reads integral file loading common /integ_h

      if (lhetero) then
         if (hname(1:4).eq.'M84A') ismax=8
         if (hname(1:4).eq.'M84C') ismax=8
         if (hname(1:4).eq.'SH8W') ismax=8
         if (hname(1:4).eq.'SH8U') ismax=8
         if (hname(1:4).eq.'SH12') ismax=12
         do ii = 1,nentry
            call ang_param(dble(elat),dble(elon),slat(ii),slon(ii),
     &           dist(ii),azmp(ii),theta(ii),avg_g,avg_d,ismax)
            do iss = 0,ismax
               do itt=-iss,iss
                  avg_great(itt,iss,ii)=avg_g(itt,iss)
                  avg_diff(itt,iss,ii)=avg_d(itt,iss)
               enddo
            enddo
         enddo
         call readint_h(hname,h_file)
      endif
      
cc...implement crust5.1 (Mooney et al) crustal corrections
c    This requires the crust5.1 model, as well as partial
c    derivative kernels for this mode table in the output
c    format of frechet.f

c      print*, 'a lcrust51=', lcrust51
      if (lcrust51) then
c         print*, 'b lcrust51=', lcrust51
         call load_kern_c51(c_file)
         call load_crust51
         if (c51opt(1:2).eq.'sa') then
            if (ldbdb5) then 
               modeid = 'Cr5.1/db5 sa'
            else
               modeid = 'Crust5.1 sa'
            endif
            do ii = 1,nentry
               print*,'Calling avg crust for; ',elat,elon,del(ii)
               call avg_crust51(elat,elon,del(ii),azm(ii),ldbdb5)
               cvp1gc(ii)=vvp1sa
               cvs1gc(ii)=vvs1sa
               crho1gc(ii)=vrho1sa
               cvp2gc(ii)=vvp2sa
               cvs2gc(ii)=vvs2sa
               crho2gc(ii)=vrho2sa
               ctopgc(ii)=vtopsa
               cbasgc(ii)=vbassa
               cmohgc(ii)=vmohsa
            end do
         else
            if (ldbdb5) then 
               modeid = 'Cr5.1/db5 gc'
            else
               modeid = 'Crust5.1 gc'
            endif
            do ii = 1,nentry
               call avg_crust51(elat,elon,del(ii),azm(ii),ldbdb5)
               fsa = dist(ii)/360.
               cvp1gc(ii)=vvp1gc
               cvs1gc(ii)=vvs1gc
               crho1gc(ii)=vrho1gc
               cvp2gc(ii)=vvp2gc
               cvs2gc(ii)=vvs2gc
               crho2gc(ii)=vrho2gc
               ctopgc(ii)=vtopgc
               cbasgc(ii)=vbasgc
               cmohgc(ii)=vmohgc
               cvp1la(ii)=(vvp1gc-fsa*vvp1sa)/(1.-fsa)
               cvs1la(ii)=(vvs1gc-fsa*vvs1sa)/(1.-fsa)
               crho1la(ii)=(vrho1gc-fsa*vrho1sa)/(1.-fsa)
               cvp2la(ii)=(vvp2gc-fsa*vvp2sa)/(1.-fsa)
               cvs2la(ii)=(vvs2gc-fsa*vvs2sa)/(1.-fsa)
               crho2la(ii)=(vrho2gc-fsa*vrho2sa)/(1.-fsa)
               ctopla(ii)=(vtopgc-fsa*vtopsa)/(1.-fsa)
               cbasla(ii)=(vbasgc-fsa*vbassa)/(1.-fsa)
               cmohla(ii)=(vmohgc-fsa*vmohsa)/(1.-fsa)
            enddo
         endif
      endif


cc...implement dbdb5 bathymetry corrections
c    This currently uses as slow subroutine call.  Needs partial
c    derivative kernels for this mode table in the output
c    format of frechet.f

      if (ldbdb5 .and. .not.(lcrust51)) then
         call load_kern_c51(d_file)
         if (cdbopt(1:2).eq.'sa') then
            modeid = 'dbdb5 sa'
            do ii = 1,nentry
c     print*,'Calling avg dbdb5 for; ',elat,elon,del(ii)
               call avg_dbdb5(elat,elon,del(ii),azm(ii),dbsa,dbgc)
c     print*,elat,elon,del(ii),azm(ii),dbsa,dbgc
               cvp1gc(ii)=0.
               cvs1gc(ii)=0.
               crho1gc(ii)=0.
               cvp2gc(ii)=0.
               cvs2gc(ii)=0.
               crho2gc(ii)=0.
               ctopgc(ii)=dbsa
               diff = rsol-ctopgc(ii)
               cbasgc(ii)=rbase-diff
               cmohgc(ii)=rmoho-diff
               print*,'Crustal thickness preserved:'
               print*,'Topo: ',dbsa,' Base: ',cbasgc(ii),
     &              ' Moho: ',cmohgc(ii)
            end do
         else
            modeid = 'dbdb5 gc'
            do ii = 1,nentry
               call avg_dbdb5(elat,elon,del(ii),azm(ii),dbsa,dbgc)
               fsa = dist(ii)/360.
               cvp1gc(ii)=0.
               cvs1gc(ii)=0.
               crho1gc(ii)=0.
               cvp2gc(ii)=0.
               cvs2gc(ii)=0.
               crho2gc(ii)=0.
               ctopgc(ii)=dbgc
               diff = rsol-ctopgc(ii)
               cbasgc(ii)=rbase-diff
               cmohgc(ii)=rmoho-diff
               print*,'Crustal thickness preserved:'
               print*,'Topo: ',dbgc,' Moho: ',cmohgc(ii)
               cvp1la(ii)=0.
               cvs1la(ii)=0.
               crho1la(ii)=0.
               cvp2la(ii)=0.
               cvs2la(ii)=0.
               crho2la(ii)=0.
               ctopla(ii)=(ctopgc(ii)-fsa*dbsa)/(1.-fsa)
               diff = rsol - dbsa
               vmohsa = rmoho-diff
               cmohla(ii)=(cmohgc(ii)-fsa*vmohsa)/(1.-fsa)
               vbassa = rbase-diff
               cbasla(ii)=(cbasgc(ii)-fsa*vbassa)/(1.-fsa)
            enddo
         endif
      endif

cc------------------------------------------------------------------
cc   LOOP OVER STATIONS
cc   
cc------------------------------------------------------------------

      do  ista = 1, nentry
c         print*, 'here 1.....'
cc....................................................
cc   SOURCE STRAINS
cc
cc   get mode file for source scalars and calculate the 
cc   appropriate strains also calculate acceleration scalar 
cc   and load w and alpha for cosine recursion
cc   the mode table is read only once (when lread=.false.)
cc
cc....................................................
 
c      total = etime(utime1)    

         call source_strains_h2(source_file,d0,num_modes,ista,
     &        lread,lhetero)


cc... check if mode table and integral table are consistent
         print*, 'here a.....'
         if (lhetero.and..not.lread) then
            if (jcom.ne.jcom_h) then
               print*,' Mode table    for jcom : ',jcom
               print*,' Integral file for jcom : ',jcom_h
               print*,' INCONSISTENT'
               stop
            endif
         endif
         
cc....................................................
cc    EXCITATION ONLY
cc
cc    test for mode type and form excitation kernels
cc    a(num_modes,maxcomp)
cc
cc    for source phase - do not calculate seismograms
cc....................................................

         if (lexcite) then
            do ii = 1, 5
               evt(ii) = event_time(ii)
            end do
            dep = d0
            if (jcom .ne. 2) then
               if (lcomp .eq. 'T') then
                  print*,' no excitation for T comp and 
     &                     spheroidal modes'
               else
                  call spheroidal_exc_h2(lcomp,num_modes,ista,
     &                 del,azm,moment,lhetero)
               endif
            else
               if (lcomp .eq. 'R') then
                  print*,' no excitation for R comp and toroidal modes'
               else
                  call toroidal_exc_h2(lcomp,num_modes,ista,
     &                 del,azm,moment,lhetero)
               endif
            endif
         endif

cc....................................................
cc   SUMMATION
cc
cc   for summation - do calculate seismograms
cc
cc....................................................

         if (.not.lexcite) then
            if (jcom .ne. 2) then
               call spheroidal_h2(lcomp,num_modes,ista,
     &              del,azm,moment,index2,lhetero)
            else
               call toroidal_h2(lcomp,num_modes,ista,
     &              del,azm,moment,index2,lhetero)
            endif
            
             print*,'back from summation'
cc...do the cosine recursion and matrix 
cc   multiplication to form the time series
c            print*, 'here 2.....'
            call curse_h2(num_modes,ista,npts,index2,dat)

            do kk = 1, index2
               if (index2.eq.6) then
                  comp=compt(kk)
                  do ll = 1, npts
                     data(ll) = dat(ll,kk)
                  end do
               else
                  do ll = 1, npts
                     data(ll) = dat(ll,kk)
                  end do
               endif
            enddo
         endif 

         print*,' elapsed time for this station: '


cc....................................................
cc   WRITE SAC_FILE
cc
cc   if summation
cc....................................................
         len1=lnblnk(sourid)
         len2=lnblnk(stname(ista))
c         sacfnam=sourid(1:len1)//"."//stname(ista)(1:len2)
c     &        //"."//lcomp//".sac"
         sacfnam=stname(ista)(1:len2)//"."//lcomp//".sac"
c         open(12,file=stname(ista)(1:len2)//".asc",form='formatted')
c         do ll=1,npts
c            write(12,*) (ll-1)*dt, data(ll)
c         enddo
c         close(12)

         print*, 'stname: ', stname(ista)(1:len2)
         print*, 'stla= ', slat(ista)
         print*, 'stlo= ', slon(ista)
         print*, 'dist= ', del(ista)
         print*, 'azim= ', azm(ista)
         call newhdr
         call setnhv('NPTS',npts,nerr)
         call setfhv('DELTA',real(dt),nerr)
         call setfhv('B', 0.0, nerr)
         call setfhv('EVLA',real(elat),nerr)
         call setfhv('EVLO',real(elon),nerr)
         call setfhv('EVDP',real(d0),nerr)
         call setfhv('STLA',real(slat(ista)),nerr)
         call setfhv('STLO',real(slon(ista)),nerr)
c         call setfhv('DIST',real(del(ista)),nerr)
c         call setfhv('AZ',real(azm(ista)),nerr)
         call wsac0(sacfnam,rdump,data,nerr)


      enddo
      
      end program
