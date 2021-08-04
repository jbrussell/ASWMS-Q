      SUBROUTINE WTABLE(IIN,IOUT,IOEIG)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** MAKES UP TABLE OF FREQUENCIES ***
c
c    calculate eigenfrequencies, eigenfunctions
c    write earth model, mode file info to binary eig-file
c    write eigenfrequencies, eigenfunctions to binary eig-file
c
c    version that recalculates wmax to account for an input upper
c    phase velocity cutoff.  Makes use of the parameter WMAXcv, 
c    which is WMAX for each l.
c
c    calls: DETQN, ENTRY, ROTSPL, STEPS, whead
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
c
      include 'parameter.h'
c
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      COMMON/SHANKS/B(46),C(10),DX,STEP(8),STEPF,IN,MAXO
      COMMON/MTAB/WE(nbranch2),DE(nbranch2),KE(nbranch2)
      DIMENSION WT(nbranch2)
      DATA NMX/nbranch/,INSS/5/
      DATA B/46*0.D0/,C/10*0.D0/,STEP/8*0.D0/
      DATA WE/nbranch2*0.D0/,DE/nbranch2*0.D0/,KE/nbranch2*0/,
     &     WT/nbranch2*0.D0/
      data nexit /0/
      data rn /6371000.d0/
c
      STEPF=1.D0
      READ(IIN,*) EPS,EPS1,EPS2,WGRAV
      WRITE(IOUT,100) EPS,EPS1,WGRAV
  100 FORMAT(/,'INTEGRATION PRECISION =',G12.4,'  ROOT PRECISION =',
     +   G12.4,'  GRAVITY CUT OFF =',G12.4,' RAD/S',///,6X,'MODE',
     +   8X,'W(RAD/S)',7X,'W(MHZ)',10X,'T(SECS)',6X,'GRP VEL(KM/S)',
     +   8X,'Q',13X,'RAYLQUO',/)
c
      CALL STEPS(EPS)
c
    5 READ(IIN,*) JCOM
      IF(JCOM.LE.0) RETURN
c
      READ(IIN,*) LMIN,LMAX,WMIN,WMAX,NBRAN,cvmax
      IF(NBRAN.LE.0.OR.NBRAN.GT.NMX) NBRAN=NMX
      IF(LMIN.LE.0) LMIN=1
      NEV=2
c
c     write header records to output file
c
      call whead(ioeig, jcom, wmin, wmax, lmin, lmax, wgrav)
c
      wmin=2.d-3*pi*wmin
      wmax=2.d-3*pi*wmax
      WT(1)=WMIN
      WT(2)=WMAX
      IF(JCOM.NE.1) GOTO 6
      LMIN=0
      LMAX=0
    6 L=LMIN-1
c
c TOP of big loop over l
c
   10 L=L+1
c
c     determine the number of modes (radial orders) for this frequency interval and angular order
c
      IF(L.GT.LMAX.OR.WT(1).GE.WMAX) GOTO 5
      KNSW=1
      MAXO=INSS
      FL=L
      FL1=FL+1.D0
      FL2=FL+FL1
      FL3=FL*FL1
      SFL3=DSQRT(FL3)
cc** Reset wmax, wt(nev) (returned from rotspl) according to cvmax
      if (cvmax.gt.0.) then
         wcv = cvmax*sfl3/rn
         if (wcv.lt.WMAX) then
           WT(NEV) = wcv
           WMAXcv = wcv
         else
           WT(NEV) = WMAX
           WMAXcv = WMAX
         end if
      else
         WMAXcv = WMAX
      end if
C*** DETERMINE MODE COUNT *** 
      WE(1)=WT(1)
c
      CALL DETQN(WE(1),KE(1),DE(1),0)
c
      IMAX=2*NMX
      DO 15 I=2,IMAX
      WE(I)=WMAXcv
   15 KE(I)=-10
      kmod=1
      DO 20 I=2,NEV
      kmod=kmod+1
c
      CALL ENTRY(WT(I),IMAX,KEI)
c
      NMODE=KEI-KE(1)
      IF(NMODE.GE.NBRAN) GOTO 25
   20 CONTINUE
   25 continue
      PRINT 900,WT(1),WT(kmod),NMODE,NBRAN
  900 FORMAT(' COUNT BETWEEN ',G14.6,' AND',G14.6,' =',I8,' : MAX =',I4)
      IF(NMODE.LE.0) then
        if(cvmax.le.0 .or. WMAXcv.ge.WMAX) then
          GOTO 5
        else
          GOTO 10
        endif
      ENDIF
      IF(NMODE.GT.NBRAN) then
        NMODE=NBRAN
        print*, ' max radial order of modes is set to',nbranch
      endif
      IMAX=2*NMODE
C*** FILL UP TABLE USING BISECTION ***
      INDX=2
      ICHK=KE(1)+1
      icount = 0
   35 IF(KE(INDX).NE.ICHK) GOTO 40
      INDX=INDX+2
      IF(INDX.GT.IMAX) GOTO 60
      ICHK=ICHK+1
      GOTO 35
   40 I1=INDX-1
   45 INDX=INDX+2
      IF(INDX.GE.IMAX) GOTO 46
      ICHK=ICHK+1
      IF(KE(INDX).NE.ICHK) GOTO 45
   46 I2=MIN0(INDX,IMAX)
      icount = icount + 1
c
c     debugging problem where it get's stuck here
c     next section disabled
c
c      go to 4444
      if (icount .gt. 10*nmode) then
        print*,' possible infinite loop'
        print*,' reinitializing arrays'
c
c       resetting values to re-initialize loop
c
        l = l - 1
        ke(1) = 0
        de(1) = 0.d0
        nev = 2
        wt(1) = wmin
        wt(2) = wmax
        nexit = nexit + 1
        if (nexit .gt. 50) then
          print*,' i quit'
          stop
        endif
        go to 10
      endif
4444  continue
      WTST=0.5D0*(WE(I2)+WE(I1))
      IF((WE(I2)-WE(I1))/WTST.LT.EPS2) GOTO 50
c
      CALL ENTRY(WTST,IMAX,KTST)
c
      INDX=I1+1
      ICHK=KE(I1)+1
      GOTO 35
   50 PRINT 901,WE(I1),KE(I1),WE(I2),KE(I2)
  901 FORMAT(' PROBLEM IN TABLE : ',2(G16.8,I5))
      J1=I1+1
      J2=I2-1
      DO 55 I=J1,J2
   55 DE(I)=1.D0
      INDX=I2+2
      IF(INDX.GE.IMAX) GOTO 60
      ICHK=KE(I2)+1
      GOTO 35
   60 NEV=IMAX
ccccccccccccccccccccccccccccccccccc
c     Check omega boundaries that each root is within.
c       Second boundary value MUST be larger than first.
c       If this isn't true, reset 2nd b.p. by looking
c       at previous and/or following pairs.
c       JBG 1.91
c
cc      do 70 jg=2,nev,+2
cc       if (we(jg) .lt. we(jg-1)) then
cc       print*,'problems with bounds, resetting (Alg. 1); l,jg= ',l,jg
cc         if (nev.gt.2) then
cc           if (jg.eq.nev) then
cc             we(jg)=(we(jg-2)-we(jg-3)) + we(jg-1)
cc           else
cc             if (jg.gt.2) then
cc               we(jg)=(we(jg+2)-we(jg+1)+we(jg-2)-we(jg-3))/2.d0 
cc     &               + we(jg-1)
cc             else
cc               we(jg)=(we(jg+2)-we(jg+1)) + we(jg-1)
cc             endif
cc           endif
cc         else
cc           we(jg)=wmax
cc         endif
cc       endif
cc70    continue
c
c     The algorithm above fails when the frequency spacing of the
c     preceeding & following pairs differs a lot.  This may lead to
c     an upper bound that encompasses more than one mode which than
c     can lead to eigenfrequencies not being found.  Therefore, we
c     use a different algorithm by resetting the upper bound to the
c     frequency of the following mode's lower bound * (1-eps2)
      do jg=2,nev,+2
        if (we(jg) .lt. we(jg-1)) then
          print*,'problems with bounds, resetting (Alg. 2); l,jg= ',l,jg
          if (jg.eq.nev) then
            we(jg)=wmax
          else
            we(jg)=we(jg+1)-we(jg+1)*eps2
          endif
        endif
      enddo
ccccccccccccccccccccccccccccccccccc	   
      PRINT 902,(I,WE(I),KE(I),DE(I),I=1,NEV)
  902 FORMAT(I5,1PD17.7,I6,1PD17.7)
C*** FIND ROOTS ***
      KNSW=0
      MAXO=8
c
      CALL ROTSPL(NEV,EPS1,WMIN,WMAXcv,WT,IOUT,IOEIG)
c
      GOTO 10
c
      END
