      SUBROUTINE SPROP(LI,JF,JL,F,IEXP)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C    SPROP PROPAGATES THE FUNDAMENTAL MATRIX F FROM JF TO JL (A SOLID REGION)
C    IF IORTH=1 THE COLUMNS OF F ARE ORTHOGONALIZED AT EACH LEVEL
C    EXCEPT IN REGIONS OF OSCILLATORY P AND S.
c
c    calls: BAYLIS, ORTHO
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
c
      include 'parameter.h'
c
      REAL*8 LCON,NCON,LSPL,NSPL,NN,LL
      COMMON R(nknot),FMU(nknot),FLAM(nknot),QSHEAR(nknot),
     &       QKAPPA(nknot),XA2(nknot),XLAM(nknot),RHO(nknot),
     &       QRO(3,nknot),G(nknot),QG(3,nknot),FCON(nknot),
     &       FSPL(3,nknot),LCON(nknot),LSPL(3,nknot),NCON(nknot),
     &       NSPL(3,nknot),CCON(nknot),CSPL(3,nknot),ACON(nknot),
     &       ASPL(3,nknot)
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      COMMON/EIFX/AR(14,nknot),INORM(nknot),idum2(nknot)
      COMMON/AREM/A(6,3,nknot)
      COMMON/RINDX/NIC,NOC,NSL,NICP1,NOCP1,NSLP1,N
      COMMON/SHANKS/B(46),C(10),DX,STEP(8),SDUM,IN,IDUM
      DIMENSION F(6,3),S(6,3),H(6,3,10)
      DATA ECONST/1048576.D0/
      DATA S/18*0.D0/,H/180*0.D0/
c
      KK=KG+2
      JJ=2*KK
      JUD=1
      IF(JL.LT.JF) JUD=-1
      Y=R(JF)
      I=JF
      GO TO 80
   10 X=Y
      Y=R(I)
      IF(X.EQ.Y) GOTO 80
      IQ=MIN0(I,I-JUD)
      QFF=1.D0+XLAM(IQ)*FCT
      QLL=1.D0+QSHEAR(IQ)*FCT
      QAA=1.D0+XA2(IQ)*FCT
      ZS=DMIN1(X,Y)
      XI=G(I)/Y
      VPSQ=(FLAM(I)+2.D0*FMU(I))/RHO(I)
      VSSQ=FMU(I)/RHO(I)
      ALFSQ=(WSQ+4.D0*RHO(I)+XI)/VPSQ
      BETASQ=WSQ/VSSQ
      DELSQ=DSQRT((BETASQ-ALFSQ)**2+4.D0*FL3*XI*XI/(VSSQ*VPSQ))
      FKSQ=.5D0*(ALFSQ+BETASQ+DELSQ)
      AL=FL3/(X*X)
      JORTH=1
      AQ=FKSQ-DELSQ-AL
      IF(AQ.GT.0.D0) JORTH=0
      QS=DSQRT(DABS(FKSQ-AL))+1.D0/ZS
      QF=DSQRT(DABS(AQ))+1.D0/ZS
      Q=DMAX1(SFL3/X,QS,QF)
      DEL=JUD*STEP(8)/Q
      DXS=0.D0
   15 Y=X+DEL
      IF(FLOAT(JUD)*(Y-R(I)).GT.0.D0) Y=R(I)
      DX=Y-X
c
      IF(DX.NE.DXS) CALL BAYLIS(Q,7)
c
      DXS=DX
      DO 20 K=1,KK
      DO 20 J=1,JJ
   20 S(J,K)=F(J,K)
      DO 50 NI=1,IN
      Z=X+C(NI)
      T=Z-ZS
      RO=RHO(IQ)+T*(QRO(1,IQ)+T*(QRO(2,IQ)+T*QRO(3,IQ)))
      GR=G(IQ)+T*(QG(1,IQ)+T*(QG(2,IQ)+T*QG(3,IQ)))
      FF=(FCON(IQ)+T*(FSPL(1,IQ)+T*(FSPL(2,IQ)+T*FSPL(3,IQ))))*QFF
      LL=(LCON(IQ)+T*(LSPL(1,IQ)+T*(LSPL(2,IQ)+T*LSPL(3,IQ))))*QLL
      IF(IFANIS.NE.0) GOTO 25
      NN=LL
      CC=FF+LL+LL
      AA=CC
      GOTO 30
   25 NN=(NCON(IQ)+T*(NSPL(1,IQ)+T*(NSPL(2,IQ)+T*NSPL(3,IQ))))*QLL
      CC=(CCON(IQ)+T*(CSPL(1,IQ)+T*(CSPL(2,IQ)+T*CSPL(3,IQ))))*QAA
      AA=(ACON(IQ)+T*(ASPL(1,IQ)+T*(ASPL(2,IQ)+T*ASPL(3,IQ))))*QAA
   30 ZR=1.D0/Z
      SFL3Z=SFL3*ZR
      ROGR=RO*GR
      C11=1.D0/CC
      C22=1.D0/LL
      DMG=AA-NN-FF*FF*C11
      ZDMG=ZR*DMG
      T11=-2.D0*FF*ZR*C11+ZR
      T12=SFL3Z*FF*C11
      T21=-SFL3Z
      T22=ZR+ZR
      S22=-RO*WSQ
      S11=S22+4.D0*ZR*(ZDMG-ROGR)
      S22=S22+ZR*ZR*(FL3*(DMG+NN)-NN-NN)
      S12=SFL3Z*(ROGR-ZDMG-ZDMG)
      IF(KG.EQ.0) S11=S11+4.D0*RO*RO
      IF(KG.EQ.0) GOTO 35
      T31=-4.D0*RO
      T33=-FL*ZR
      S13=-FL1*ZR*RO
      S23=RO*SFL3Z
   35 DO 70 K=1,KK
      IF(KG.EQ.1) GOTO 40
      H(1,K,NI)=T11*F(1,K)+T12*F(2,K)+C11*F(3,K)
      H(2,K,NI)=T21*F(1,K)+T22*F(2,K)+C22*F(4,K)
      H(3,K,NI)=S11*F(1,K)+S12*F(2,K)-T11*F(3,K)-T21*F(4,K)
      H(4,K,NI)=S12*F(1,K)+S22*F(2,K)-T12*F(3,K)-T22*F(4,K)
      GOTO 45
   40 H(1,K,NI)=T11*F(1,K)+T12*F(2,K)+C11*F(4,K)
      H(2,K,NI)=T21*F(1,K)+T22*F(2,K)+C22*F(5,K)
      H(3,K,NI)=T31*F(1,K)+T33*F(3,K)+4.D0*F(6,K)
      H(4,K,NI)=S11*F(1,K)+S12*F(2,K)+S13*F(3,K)-T11*F(4,K)-T21*F(5,K)
     +    -T31*F(6,K)
      H(5,K,NI)=S12*F(1,K)+S22*F(2,K)+S23*F(3,K)-T12*F(4,K)-T22*F(5,K)
      H(6,K,NI)=S13*F(1,K)+S23*F(2,K)-T33*F(6,K)
   45 GO TO (701,702,703,704,705,706,707,708,709,710),NI
  701 DO 7001 J=1,JJ
 7001 F(J,K)=S(J,K)+B(1)*H(J,K,1)
      GO TO 70
  702 DO 7002 J=1,JJ
 7002 F(J,K)=S(J,K)+B(2)*(H(J,K,1)+B(3)*H(J,K,2))
      GO TO 70
  703 DO 7003 J=1,JJ
 7003 F(J,K)=S(J,K)+B(4)*(H(J,K,1)+B(5)*H(J,K,2)+B(6)*H(J,K,3))
      GO TO 70
  704 DO 7004 J=1,JJ
 7004 F(J,K)=S(J,K)+B(7)*(H(J,K,1)+B(8)*H(J,K,2)+B(9)*H(J,K,3)+
     +B(10)*H(J,K,4))
      GO TO 70
  705 DO 7005 J=1,JJ
 7005 F(J,K)=S(J,K)+B(11)*(H(J,K,1)+B(12)*H(J,K,2)+B(13)*H(J,K,3)+
     +B(14)*H(J,K,4)+B(15)*H(J,K,5))
      GO TO 70
  706 DO 7006 J=1,JJ
 7006 F(J,K)=S(J,K)+B(16)*(H(J,K,1)+B(17)*H(J,K,2)+B(18)*H(J,K,3)+
     +B(19)*H(J,K,4)+B(20)*H(J,K,5)+B(21)*H(J,K,6))
      GO TO 70
  707 DO 7007 J=1,JJ
 7007 F(J,K)=S(J,K)+B(22)*(H(J,K,1)+B(23)*H(J,K,3)+B(24)*H(J,K,4)+
     +B(25)*H(J,K,5)+B(26)*H(J,K,6)+B(27)*H(J,K,7))
      GO TO 70
  708 DO 7008 J=1,JJ
 7008 F(J,K)=S(J,K)+B(28)*(H(J,K,1)+B(29)*H(J,K,3)+B(30)*H(J,K,4)+
     +B(31)*H(J,K,5)+B(32)*H(J,K,6)+B(33)*H(J,K,7)+B(34)*H(J,K,8))
      GO TO 70
  709 DO 7009 J=1,JJ
 7009 F(J,K)=S(J,K)+B(35)*(H(J,K,1)+B(36)*H(J,K,3)+B(37)*H(J,K,4)+
     +B(38)*H(J,K,5)+B(39)*H(J,K,6)+B(40)*H(J,K,7)+B(41)*H(J,K,8)+
     +B(42)*H(J,K,9))
      GO TO 70
  710 DO 7010 J=1,JJ
 7010 F(J,K)=S(J,K)+B(43)*(H(J,K,1)+H(J,K,10)+B(45)*H(J,K,5)+
     +B(44)*(H(J,K,4)+H(J,K,6))+B(46)*(H(J,K,7)+H(J,K,9)))
   70 CONTINUE
   50 CONTINUE
      X=Y
      IF(Y.NE.R(I)) GO TO 15
   80 SIZE=0.D0
      DO 81 K=1,KK
      DO 81 J=1,JJ
   81 SIZE=DMAX1(SIZE,DABS(F(J,K)))
   82 IF(SIZE.LT.1024.D0) GOTO 84
      DO 83 K=1,KK
      DO 83 J=1,JJ
   83 F(J,K)=F(J,K)/ECONST
      SIZE=SIZE/ECONST
      IEXP=IEXP+20
      GOTO 82
   84 INORM(I)=IEXP
      DO 85 K=1,KK
      DO 85 J=1,JJ
   85 A(J,K,I)=F(J,K)
c
      IF(JORTH.EQ.1) CALL ORTHO(LI,I,F,KG)
c
      IF(I.EQ.JL) RETURN
c
      I=I+JUD
      GO TO 10
c
      END
