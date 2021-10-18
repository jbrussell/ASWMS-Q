*  modified rdsetupsimulf from DWF.
 
c  rdsetupsimulf.2.JdF.f
c  This program is the part of simannerr that reads in the data and does the
c  Fourier analysis.  It is separated so that simannerr# can run on katmai or
c  other machines that do not have sac.  Writes out the data ready for use in
c  simannerr#  
c  Pipe in data, e.g.,   rdsetupsimul < eqlistper50
c
c  Current version does make station corrections.
      parameter (maxnfreq=20, maxnsta=300, 
     1     maxpts = 60000, maxnodes = 1000, maxevnts = 500)

      real*4 staph(maxevnts,maxnsta,maxnfreq)
      real*4 staamp(maxevnts,maxnsta,maxnfreq)
      real*4 freq(maxnfreq)
      real*4 stadist(maxevnts,maxnsta), staazi(maxevnts,maxnsta)
      real*4 stacor(maxnsta), tempcor(maxnsta),geomsprd(maxnsta)
      real*4 bazi(maxevnts,maxnsta), stadelt(maxevnts,maxnsta)
      real*4 tdata(maxpts),beg(maxevnts,maxnsta),delt(maxnsta)
      real*4 stalat(maxevnts,maxnsta),stalon(maxevnts,maxnsta)
      real*4 attnfac(maxnfreq),fattnfac(10,2)
      integer*4 nsta(maxevnts), dot, idot
      integer*4 nfreq,nstapts(maxnsta)
      integer*4 istanum(maxnsta),istacor(maxnsta), nevents
      character*70 foutput, fsummary, fstalist
c   modify here       
      character*125 fn(maxevnts,maxnsta)
      character*70 finvrsnodes, fftoutput, fvariance,fmaxavamp
      character*70 ftemp,fvelarea,sensfn,fvelout,startvel
      character*70 fresdiag,fendvel,fresmat      
      character*60 staampcor
      character*60 dummy1
      character*3 idnum(maxnsta)
      character*2 nettemp
      character*4 statemp,staname(maxnsta)
      character*8 statemp2
      pi = 3.1415928
      convdeg = 3.1415928/180.
      circ = 6371.*3.1415928/180.
      twopi = 3.1415928*2.
c  **************************************************
c  WARNING:  The following two statements need to be switched depending on which
c  compiler is used
c  ************************************************

200         format(a75, i2)
201         format(a70, a2)
c202         format(a60, a4)
202         format(a85, a3)
c202         format(a39,a3)
c202         format(a41,a3)
c210	    format(a49, i2)

c	cattnfac=0.25e-3
      cattnfac=5e-5

c  read list of files and frequencies to be analyzed and file to output results
c  Usually will pipe in data from some file like fasearrayinp
      read(*,*) nevents
      nobs = 0
      write(*,*) nevents
      do iev = 1, nevents
        read(*,*) nsta(iev),idnum(iev)
c	write(*,*) nsta(iev),idnum(iev)
        nobs = nobs+ 2*nsta(iev)
        read(*,'(a)') (fn(iev,i), i=1,nsta(iev))
c	write(14,'(a)') (fn(iev,i), i=1,nsta(iev))
c        write(*,*) iev, nsta(iev),idnum(iev)
      enddo
c  a lot of the following input is unneeded, just here to parallel simannerr#
c  input, so can use same fasearrayinp# or eqlistper# file for both
      read(*,*) nfreq
      read(*,*) (freq(j), j= 1, nfreq)
      read(*,*) foutput
      read(*,*) fsummary
      read(*,*) finvrsnodes
      read(*,*) fftoutput
      read(*,*) fvariance
      read(*,*) fmaxavamp
      read(*,*) ftemp
      read(*,*) fvelarea
      read(*,*) fstalist
      read(*,*) unifvel
c      read(*,*) iterlimit, dampvel,dampaniso, divfac
      read(*,'(1x)')
      read(*,'(a)') sensfn
      read(*,*) fvelout
      read(*,*) startvel
      read(*,*) fresdiag
      read(*,*) fendvel
      read(*,*) fresmat
      
      read(*,*) staampcor


c      idot=dot(foutput)
c      write(*,*) idot
c      staampcor = 'staampcor'//foutput(idot+1:idot+4)
c      write(14,*) staampcor
c      staampcor = 'staampcorBR3.dat' 	
      open(11, file = fftoutput)
      open(12, file = staampcor, status='old')
      open(13, file = ftemp)
      open(14, file = 'followit13')
      open(18, file = fstalist)


c  first, read in master list of stations    
      do ista2 = 1, maxnsta
	read(18,*) staname(ista2)
	stacor(ista2) = 1.0
	if (staname(ista2).eq.'nope') then
	  mxnsta = maxnsta -1
	  go to 1111
	endif
      enddo
1111  continue

c  fetch station corrections for amplitudes and assign to correction array
      read(12,*) nstacor
      do i = 1,nstacor
        read(12,*) istacor(i),tempcor(i)
      enddo      
      do i = 1, nstacor
        stacor(istacor(i)) = tempcor(i)
      enddo

c  start input loop over events
      do iev = 1, nevents
        write(14,*) idnum(iev), nsta(iev)
        do ista = 1, nsta(iev)
          call rsac1(fn(iev,ista),tdata,nstapts(ista),beg(iev,ista),
     1       delt(ista),maxpts,nerr)
          call getfhv('DIST',stadist(iev,ista),nerr)
          call getfhv('AZ', staazi(iev,ista),nerr)
          call getfhv('BAZ', bazi(iev,ista),nerr)
          call getfhv('GCARC', stadelt(iev,ista), nerr)
          call getfhv('STLA', stalat(iev,ista),nerr)
          call getfhv('STLO', stalon(iev,ista),nerr)
	  call getkhv('KSTNM', statemp2,nerr)
c	  rembazi = bazi(iev,1)
c	  remdist = stadist(iev,1)
c	  if ((abs(rembazi- bazi(iev,ista)).gt.4.0).or.
c     1       (abs(remdist-stadist(iev,ista)).gt.600.)) then    
c	  endif
	  
          geomsprd(ista) = sqrt(sin(stadelt(iev,ista)*convdeg))
c  do following silly steps to extract station number from filename
          write(13,'(a)') fn(iev,ista)
c	  write(*,*) fn(iev,ista)
          rewind (13)
      
c ********************************modify here***********************************
	  read (13,202) dummy1, statemp
	  rewind (13)
	  write(14,*) statemp,statemp2
	  istanum(ista) = 0
	  do ista2 = 1, mxnsta
	     if (statemp.eq.staname(ista2)) then
	       istanum(ista) = ista2
c	       write(*,2070)stalat(iev,ista), stalon(iev,ista),statemp,ista2
c2070   format(f10.4, f10.4,2x, a4, 2x, i3)
	     endif
	  enddo
	  if (istanum(ista).eq.0) then
	    write(*,*) statemp, ista
	    write(*,*) 'WARNING ', statemp,' not in station list'
	  endif
c          write(13,'(a)') fn(iev,ista)
c          rewind (13)
	  
c ********************************modify here***********************************
	  
c
c  get phases and amplitudes at desired frequencies
c
         
          do ifreq = 1, nfreq
            call frt(tdata, freq(ifreq), staamp(iev,ista,ifreq),
     1                      staph(iev,ista,ifreq),
     1                      nstapts(ista),delt(ista))
c  correct amplitudes for geometrical spreading, attenuation and station factor
            attneffect= exp(cattnfac*(stadist(iev,ista)-
     1                                         stadist(iev,1)))
c            attneffect= exp(attnfac(ifreq)*(stadist(iev,ista)-
c     1                                         stadist(iev,1)))
            staamp(iev,ista,ifreq)= staamp(iev,ista,ifreq)*attneffect
     1         *(geomsprd(ista)/stacor(istanum(ista)))
c            staamp(iev,ista,ifreq)= staamp(iev,ista,ifreq)*attneffect
c     1         *(geomsprd(ista))

          enddo  
        enddo
      enddo
c  end of input loop over events
c  output loop over events
      do iev = 1, nevents
        write (11,*) iev, idnum(iev)
        do ista = 1, nsta(iev)
          write(11,*) beg(iev,ista)
       write (11,*) stadist(iev,ista),staazi(iev,ista), bazi(iev,ista),
     1    stadelt(iev,ista), stalat(iev,ista), stalon(iev,ista)
          write(11,*) (staamp(iev,ista,ifreq),staph(iev,ista,ifreq),
     1                               ifreq=1,nfreq)
        enddo
      enddo
      close(unit = 13,status='delete')
      close(unit = 11)
      close(unit = 12)
      close(unit = 14)
      close(unit = 18)
      end

c---positive fourier transform
c 
      SUBROUTINE FRT(UP,FR,ARZ,PRZ,NALL,DELT)
      dimension UP(40000)
      DIMENSION W(40000)
      THETA=6.283185*FR*DELT
      C=COS(THETA)*2.0
      NR1=1
      NR2=NALL
      NDR1=NR2-1
      W(1)=UP(NR2)
      W(2)=C*W(1)+UP(NDR1)
      NDATA=NR2-NR1+1
      NTR1=NDATA-1
      NTR2=NDATA-2
      DO 97 I=3,NDATA
      I1=I-1
      I2=I-2
      NDRI=NR2-I+1
      W(I)=C*W(I1)-W(I2)+UP(NDRI)
97    CONTINUE
      ZRE=(W(NDATA)-W(NTR2)+UP(NR1))*DELT/2.0
      ZIM=W(NTR1)*SIN(THETA)*DELT
      CALL PHASE(ZRE,ZIM,PRZ)
      ARZ=SQRT(ZRE*ZRE+ZIM*ZIM)
      RETURN
      END

      SUBROUTINE PHASE(X,Y,PHI)
      IF(X) 21,20,22
20    IF(Y) 23,24,25
23    PHI=1.5*3.141592
      GO TO 28
24    PHI=0.0
      GO TO 28
25    PHI=0.5*3.141592
      GO TO 28
21    PHI=ATAN(Y/X) +3.141592
      GO TO 28
22    IF(Y) 26,27,27
27    PHI=ATAN(Y/X)
      GO TO 28
26    PHI=ATAN(Y/X)+2.0*3.141592
      GO TO 28
28    CONTINUE
      PHI=PHI/6.283184
      PHI=PHI-AINT(PHI)
      RETURN
      END

      integer function dot(file)
      character file*70
      do 50 i=1,70
      if(file(i:i).ne.'.') goto 50
      dot=i-1
      return
50     continue
      write(1,100) file
100   format(' no dots found in ',a70)
      dot = 0
      return
      end
