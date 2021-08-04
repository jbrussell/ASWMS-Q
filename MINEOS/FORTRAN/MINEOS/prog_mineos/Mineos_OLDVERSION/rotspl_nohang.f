      SUBROUTINE ROTSPL(NEV,EPS1,WMIN,WMAX,WT,IOUT,IOEIG)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** FIND ROOTS BY SPLINE INTERPOLATION ***
c
c    calls: DETQN, DSPLIN, MODOUT
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
c
      include 'parameter.h'
c
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      COMMON/MTAB/WE(nbranch2),DE(nbranch2),KE(nbranch2)
c      DIMENSION X(30),DET(30),QX(3,30),WRK(90),WT(nbranch2),ICHAR(4)
      DIMENSION X(30),DET(30),QX(3,30),WRK(90),WT(nbranch2)
      character*2 ichar(4)
      DATA X/30*0.D0/,DET/30*0.D0/,QX/90*0.D0/,WRK/90*0.D0/
      DATA TOL/1.D-9/,ITMAX/30/,ICHAR/' S',' T',' S',' C'/
      KNEV=1
      NMO=NEV/2
c
c     1/9/91 attempt to check for overlapping eigenfrequencies
c
      blast = 0.d0
c
      DO 100 I=1,NMO
      iflag1 = 0
      iflag2 = 0
5      K1=2*I
      K=K1-1
      FC=DE(K)
      FB=DE(K1)
      IF(FC*FB.GE.0.D0) GOTO 100
      DET(1)=FC
      DET(2)=FB
c     set intial guess of root.  if already tried and failed, do
c      a restart with initial guess based on last root
c
      if (iflag1.gt.0 .and. iflag2.gt.0) then
	x(1)=blast+dfloat(iflag1)*3.142d-4
        iflag2=0
        print*,'1st guess, ', x(1),' fc= ',fc,' fb= ',fb
      else
        X(1)=WE(K)
      endif
      X(2)=WE(K1)
      C=X(1)
      IF(DABS(FC).LT.DABS(FB)) GOTO 10
      C=X(2)
      FC=DET(2)
   10 J=1
      M=2
      NTRY=2
      B=0.5D0*(X(J)+X(M))
c
   15 T=DABS(B*EPS1)
c
c     1/9/91 test whether this eigenfrequency differs from the last value
c       if not, set the initial eigenfreq a little larger than previous mode,
c       and restart the search.  Note that first flag prevents infinite loops
c       for truly stuck modes, and second flag prevents reset on every
c       iteration
c
      IF (DABS(B-C).LT.T) then
        if (dabs(b-blast) .gt. t) then
          GOTO 65
        else
          print*,' overlapping eigenfrequency! ', b, blast
          iflag1=iflag1+1
          iflag2=1
          if (iflag1.gt.5) then
            write(iout,905) b, blast
            print*,'i quit'
            stop
          else
	    go to 5
	  endif
        endif
      endif
905   format('stuck on eigenfrequency, ',2f8.4)
c
      CALL DETQN(B,KNT,FB,0)
c
      PRINT 900,B,FB
  900 FORMAT(2G20.12)
      icrap=1
      DO 20 M=2,NTRY
      icrap=icrap+1
   20 IF(B.LT.X(M)) GOTO 25
   25 NTRY=NTRY+1
      m=icrap
      J2=NTRY
   30 J1=J2-1
      X(J2)=X(J1)
      DET(J2)=DET(J1)
      IF(J1.EQ.M) GOTO 35
      J2=J2-1
      GOTO 30
   35 X(M)=B
      DET(M)=FB
      DO 40 M=2,NTRY
      J=M-1
      FC=DET(J)
      FB=DET(M)
   40 IF(FC*FB.LT.0.D0) GOTO 45
   45 IND=j+1
      IF(DABS(FC).LT.DABS(FB)) IND=J
      C=X(IND)
      FC=DET(IND)
c
      CALL DSPLIN(NTRY,X,DET,QX,WRK)
c
      DEL=-DET(IND)/QX(1,IND)
   50 DELX=-DET(IND)/(QX(1,IND)+DEL*QX(2,IND))
      IF(DABS(DELX-DEL).LT.TOL) GOTO 55
      IF(DEL*DELX.LT.0.D0) GOTO 60
      DEL=DELX
      GOTO 50
   55 B=C+DELX
      IF(NTRY.GE.ITMAX) GOTO 60
      IF(B.GT.X(J).AND.B.LT.X(M)) GOTO 15
   60 X(1)=X(J)
      X(2)=X(M)
      DET(1)=DET(J)
      DET(2)=DET(M)
      GOTO 10
c
c     found an eigenfrequency
c
C*** WRITE OUT FREQUENCIES ***
c
   65 NORD=KE(K1)
      IF(L.EQ.1) NORD=NORD+1
c
      CALL DETQN(B,KNT,FB,1)
c
      PRINT 900,B,FB
      WDIFF=(B-WRAY*WN)/B
      TCOM=2.D0*PI/B
      WMHZ=1000.D0/TCOM
      GCOM=VN*CG/1000.D0
      QMOD=0.D0
      IF(QINV.GT.0.D0) QMOD=1.D0/QINV
c
c     set n=-1 if eigenfunction was skipped in detqn due to
c     likelihood of hanging -- occasional for "tough" S modes
        
      nhold = NORD
      if(cg.lt.1.d-4.and.qinv.lt.1.d-4.and.wray.lt.1.d-4) then
        NORD = -1
      endif
c
c     1/9/91 save this eigenfrequency for future comparisions
c
      blast = b
      ib = 0
c
      PRINT 200,NORD,ICHAR(JCOM),L,B,WMHZ,TCOM,GCOM,QMOD,WDIFF
      WRITE(IOUT,200) NORD,ICHAR(JCOM),L,B,WMHZ,TCOM,GCOM,QMOD,WDIFF
  200 FORMAT(I5,A2,I7,6G16.7)
c
      CALL MODOUT(B,QMOD,GCOM,IOEIG)

      NORD = nhold
    
c      print*,'B',B
      print*,'CG',CG
      print*,'WN',WN
c JBR - 2/8/16 - avoid NaN values
      IF(CG.NE.CG) THEN
        print*,'ERROR rotspl: CG = NaN'
        GOTO 100
c        return
      ENDIF
c      print*,'WN',WN
c
      BP1=B+2.D0*CG*WN
c
      if ((nord .eq. 0) .or. ((i .eq. 1) .and. (nord .ne. 1))) then
        wmin = 0.95*b
        go to 70
      endif      
      KNEV=KNEV+1
      WT(KNEV)=B
   70 IF(BP1.GE.WMAX) GOTO 100
      KNEV=KNEV+1
      WT(KNEV)=BP1
c  JBR    print*,'WT(KNEV) = ',WT(KNEV)
c  JBR    print*,'KNEV = ',KNEV
c  JBR    print*,'BP1 = ',BP1
  100 CONTINUE
      NEV=KNEV+1
      WT(1)=WMIN
      WT(NEV)=WMAX
c
      RETURN
      END
