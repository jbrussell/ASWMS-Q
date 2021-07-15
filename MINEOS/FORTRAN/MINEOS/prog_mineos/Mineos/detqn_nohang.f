      SUBROUTINE DETQN(WDIM,KNT,DET,IFEIF)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C**** SUPEVISES THE INTEGRATION OF THE EQUATIONS,IT RETURNS THE VALUE
C**** OF THE SECULAR DETERMINANT AS DET AND THE COUNT OF ZERO CROSSINGS.
c
c     when the secular determinant is too large (>0.1) do not calculate
c     eigenfunctions, etc. to avoid that the code crashes.  However,
c     the nature of these undetermined modes should be investigated
c     a posteriori, to decide if another attempt needs to be made to
c     recover these modes or if they can be safely ignored.
c     PP 11/01/96
c
c     This previous attempt discounts too many modes in some cases.
c     Many of these do not crash the code, i.e. "REMEDY" works 
c     properly (apparently).  Now a check is done for problem modes --
c     those in which the "afr" variable inside of "remedy" is zero
c     In this case, a flag is sent back and these modes are ignored.
c     The caveat on checking the importantce of these ignored
c     modes still holds.
c     JBG 01/19/00
c
c     calls: EIFOUT, FPRPMN, FPSM, FSBM, GAUSLV, REMEDY, RPROP, RPS, 
c            SDEPTH, SFBM, SPRPMN, SPSM, STARTL, TPROP, TPS
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
c
      include 'parameter.h'
c
      REAL*8 LCON,NCON,LSPL,NSPL
      COMMON R(nknot),FMU(nknot),FLAM(nknot),QSHEAR(nknot),
     &       QKAPPA(nknot),XA2(nknot),XLAM(nknot),RHO(nknot),
     &       QRO(3,nknot),G(nknot),QG(3,nknot),FCON(nknot),
     &       FSPL(3,nknot),LCON(nknot),LSPL(3,nknot),NCON(nknot),
     &       NSPL(3,nknot),CCON(nknot),CSPL(3,nknot),ACON(nknot),
     &       ASPL(3,nknot)
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      COMMON/EIFX/A(14,nknot),adum(nknot)
      COMMON/RINDX/NIC,NOC,NSL,NICP1,NOCP1,NSLP1,N
      DIMENSION ASS(14),VF(nknot),ZI(4)
      DATA ASS/14*0.D0/,VF/nknot*0.D0/,ZI/4*0.D0/
c
      IBACK=0
      W=WDIM/WN
c JBR     print*,'WDIM = ',WDIM
      WSQ=W*W
      IEXP=0
      KOUNT=0
      KG=0
      FCT=0.D0
      IF(TREF.GT.0.D0) FCT=2.D0*DLOG(TREF*WDIM)/PI
      GOTO (2,3,1,3),JCOM
c
c.... spheroidal modes
    1 IF(WDIM.LE.WGRAV) KG=1
      NVEFM=2+KG*3
      NVESM=5+KG*9
c
      CALL SDEPTH(WDIM,LS)
c
      IF(LS.GT.NOCP1) GOTO 25
      IF(LS.GT.NICP1) GOTO 20
      IF(LS.GT.2) GOTO 15
      R10=4.5D-4*(FL+.5D0)/WDIM
      IF(R10.GE.R(2)) GOTO 15
      R(1)=R10
      G(1)=RHO(1)*R(1)*1.333333333333333D0
      LS=1
c
   15 CALL SPSM(LS,NVESM,ASS)
c
C*** PROPAGATE THROUGH INNER CORE ***
c
      CALL SPRPMN(LS,NIC,ASS,VF,NVESM,IEXP)
c
      R(1)=0.D0
      G(1)=0.D0
c
      CALL SFBM(ASS,KG,IBACK)
c
   20 IS=MAX0(LS,NICP1)
c
      IF(IS.EQ.LS) CALL FPSM(LS,NVEFM,ASS)
c
C*** PROPAGATE THROUGH OUTER CORE ***
c
      CALL FPRPMN(IS,NOC,ASS,VF,NVEFM,IEXP)
c
      CALL FSBM(ASS,KG,IBACK)
c
   25 IS=MAX0(LS,NOCP1)
c
      IF(IS.EQ.LS) CALL SPSM(LS,NVESM,ASS)
c
C*** PROPAGATE THROUGH MANTLE ***
c
      CALL SPRPMN(IS,NSL,ASS,VF,NVESM,IEXP)
c
      IF(NSL.NE.N) GOTO 40
      DNORM=A(1,NSL)*A(1,NSL)
      DO 26 I=2,NVESM
   26 DNORM=DNORM+A(I,NSL)*A(I,NSL)
      DET=A(5,NSL)/DSQRT(DNORM)
      GOTO 45
c
   40 CALL SFBM(ASS,KG,IBACK)
c
C*** PROPAGATE THROUGH OCEAN ***
c
      CALL FPRPMN(NSLP1,N,ASS,VF,NVEFM,IEXP)
c
      IF(KG.EQ.0) DET=A(2,N)/DSQRT(A(1,N)*A(1,N)+A(2,N)*A(2,N))
      IF(KG.NE.0) DET=A(5,N)/DSQRT(A(1,N)**2+A(2,N)**2+A(3,N)**2+
     +   A(4,N)**2+A(5,N)**2)
   45 IF(LS.GT.NOC) DET=-DET
      IF(KNSW.NE.1) GOTO 50
      IF(LS.GT.NOC) KOUNT=KOUNT-2
      IREM=MOD(KOUNT,2)
      IF(IREM.EQ.0.AND.DET.LT.0.D0) KOUNT=KOUNT+1
      IF(IREM.NE.0.AND.DET.GT.0.D0) KOUNT=KOUNT+1
      KNT=KOUNT
   50 IF(IFEIF.EQ.0) RETURN
c
C*** THIS DOES EIGENFUNCTION CALCULATION FOR SPHEROIDAL MODES ***
      IBACK=1
      JEXP=0
      NBAKF=1+KG*3
      NBAKS=4+KG*10
      DO 55 I=1,NBAKS
   55 ASS(I)=0.D0
      IF(N.EQ.NSL) GOTO 65
      IF(KG.NE.0) GOTO 75
      ASS(1)=DSIGN(1.D0,A(1,N))
      GOTO 80
   65 IF(KG.EQ.0) GOTO 75
      ASI1=A(3,N)*A(3,N)+A(12,N)*A(12,N)
      ASI2=A(4,N)*A(4,N)+A(11,N)*A(11,N)
      IF(ASI2.LE.ASI1) ASS(1)=DSIGN(1.D0,A(3,N))
      IF(ASI2.GT.ASI1) ASS(2)=DSIGN(1.D0,A(2,N))
      GOTO 85
   75 ASI1=A(3,N)*A(3,N)
      ASI2=A(4,N)*A(4,N)
      IF(ASI2.LE.ASI1) ASS(1)=DSIGN(1.D0,A(3,N))
      IF(ASI2.GT.ASI1) ASS(2)=DSIGN(1.D0,A(2,N))
      IF(N.EQ.NSL) GOTO 85
c
   80 CALL FPRPMN(N,NSLP1,ASS,VF,NBAKF,JEXP)
c
      CALL FSBM(ASS,KG,IBACK)
c
   85 NTO=MAX0(LS,NOCP1)
c
      CALL SPRPMN(NSL,NTO,ASS,VF,NBAKS,JEXP)
c
      IF(NTO.EQ.LS) GOTO 90
c
      CALL SFBM(ASS,KG,IBACK)
c
      NTO=MAX0(LS,NICP1)
c
      CALL FPRPMN(NOC,NTO,ASS,VF,NBAKF,JEXP)
c
      IF(NTO.EQ.LS) GOTO 90
c
      CALL FSBM(ASS,KG,IBACK)
c
      NTO=MAX0(LS,2)
c
      CALL SPRPMN(NIC,NTO,ASS,VF,NBAKS,JEXP)
c
cjg >>
   90 continue
c      print*,'det = ',DABS(DET)
      IF (DABS(DET).le.1.D-4) then
c
        CALL EIFOUT(LS)
c
      else 
        print*,'DETQN, det = ',DABS(DET)
        lshold = LS
        CALL REMEDY(LS)

        if (LS.eq.-999) then
c....     do not calculate eigenfunctions
        print*,'detqn: secular det. too large -> no efun calc.',det
        print*,' i quit'
        stop
          cg = 0.d0
          wray = 0.d0
          qinv = 0.d0
          do ii = 1,n
            a(1,ii) = 0.d0
            a(2,ii) = 0.d0
            a(3,ii) = 0.d0
            a(4,ii) = 0.d0
            a(5,ii) = 0.d0
            a(6,ii) = 0.d0
          enddo

        else

          LS = lshold
          CALL EIFOUT(LS)
        endif
c
      endif
cjg <<
c
      RETURN
c
C*** RADIAL MODES ***
    2 LS=2
c
      CALL RPS(LS,ASS)
c
      CALL RPROP(LS,N,ASS)
c
      DET=A(2,N)/DSQRT(A(1,N)*A(1,N)+A(2,N)*A(2,N))
      KNT=KOUNT-1
      IF(IFEIF.EQ.0) RETURN
c
      A(1,1)=0.D0
      A(2,1)=0.D0
      DO 205 I=LS,N
      FF=FCON(I)*(1.D0+XLAM(I)*FCT)
      CC=CCON(I)*(1.D0+XA2(I)*FCT)
  205 A(2,I)=(A(2,I)-2.D0*FF*A(1,I)/R(I))/CC
      ZI(1)=0.D0
      ZI(2)=0.D0
      ZI(3)=0.D0
      DO 210 I=LS,N
      IM=I-1
c
  210 IF(R(I).NE.R(IM)) CALL GAUSLV(R(IM),R(I),IM,ZI,3)
c
      RNRM=1.D0/(W*DSQRT(ZI(1)))
      CG=0.D0
      QINV=ZI(2)/(WSQ*ZI(1))
      WRAY=DSQRT(ZI(3)/ZI(1))
      DO 215 J=1,2
      DO 215 I=2,N
  215 A(J,I)=A(J,I)*RNRM
      RETURN
c
C*** TOROIDAL MODES ***
    3 NB=NOCP1
      N2=NSL
      ASS(1)=1.D0
      ASS(2)=0.D0
      IF(JCOM.EQ.2) goto 300
      NB=2
      A(1,1)=0.D0
      A(2,1)=0.D0
      N2=NIC
  300 Q=0.D0
      LS=NB
c
      CALL STARTL(LS,N2,FMU,LSJSR,Q)
c
      LS=LSJSR
c
      IF(LS.NE.NOCP1) CALL TPS(LS,ASS)
c
      CALL TPROP(LS,N2,ASS)
c
      DET=A(2,N2)/DSQRT(A(1,N2)*A(1,N2)+A(2,N2)*A(2,N2))
      IF(IFEIF.EQ.0) GOTO 335
      DO 305 I=LS,N2
  305 A(2,I)=A(1,I)/R(I)+A(2,I)/(LCON(I)*(1.D0+QSHEAR(I)*FCT))
      IF(LS.EQ.NB) GOTO 315
      LS1=LS-1
      DO 310 I=NB,LS1
      A(1,I)=0.D0
  310 A(2,I)=0.D0
  315 DO 320 I=1,4
  320 ZI(I)=0.D0
      DO 325 I=LS,N2
      IM=I-1
c
  325 IF(R(I).NE.R(IM)) CALL GAUSLV(R(IM),R(I),IM,ZI,4)
c
      RNRM=1.D0/(W*DSQRT(ZI(1)))
      CG=(FL+0.5D0)*ZI(2)/(W*ZI(1))
      QINV=ZI(3)/(WSQ*ZI(1))
      WRAY=DSQRT(ZI(4)/ZI(1))
      DO 330 J=1,2
      DO 330 I=LS,N2
  330 A(J,I)=A(J,I)*RNRM
c
      RETURN
c
  335 IF(KNSW.NE.1) RETURN
c
      KNT=KOUNT-1
c
      IF(JCOM.EQ.4.OR.L.EQ.1) RETURN
c
      IREM=MOD(KNT,2)
c
      IF(IREM.EQ.0.AND.DET.LT.0.D0) RETURN
c
      IF(IREM.NE.0.AND.DET.GT.0.D0) RETURN
c
      KNT=KNT+1
c
      RETURN
      END
