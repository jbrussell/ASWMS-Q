      SUBROUTINE FPROP(JF,JL,F,IEXP)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C    FPROP PROPAGATES THE FUNDAMENTAL MATRIX F FROM JF TO JL (A FLUID REGION)
c
c    calls: BAYLIS
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
      COMMON/EIFX/AR(14,nknot),INORM(nknot),idum2(nknot)
      COMMON/AREM/A(6,3,nknot)
      COMMON/RINDX/NIC,NOC,NSL,NICP1,NOCP1,NSLP1,N
      COMMON/SHANKS/B(46),C(10),DX,STEP(8),SDUM,IN,IDUM
      DIMENSION F(4,2),S(4,2),H(4,2,10)
      DATA ECONST/1048576.D0/
      DATA S/8*0.D0/,H/80*0.D0/
c
      KK=KG+1
      JJ=2*KK
      JUD=1
      IF(JL.LT.JF) JUD=-1
      Y=R(JF)
      I=JF
      GO TO 80
   10 X=Y
      Y=R(I)
      IF(Y.EQ.X) GOTO 80
      IQ=MIN0(I,I-JUD)
      QFF=1.D0+XLAM(IQ)*FCT
      ZS=DMIN1(X,Y)
      XI=G(I)/Y
      ALFSQ=(WSQ+4.D0*RHO(I)+XI-FL3*XI*XI/WSQ)*RHO(I)/FLAM(I)
      Q=DMAX1(SFL3/X,DSQRT(DABS(ALFSQ-FL3/(X*X)))+1.D0/ZS)
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
      D=FL3/WSQ
      DO 40 NI=1,IN
      Z=X+C(NI)
      T=Z-ZS
      ZR=1.D0/Z
      RO=RHO(IQ)+T*(QRO(1,IQ)+T*(QRO(2,IQ)+T*QRO(3,IQ)))
      FLU=(FCON(IQ)+T*(FSPL(1,IQ)+T*(FSPL(2,IQ)+T*FSPL(3,IQ))))*QFF
      GR=(G(IQ)+T*(QG(1,IQ)+T*(QG(2,IQ)+T*QG(3,IQ))))*ZR
      T21=-4.D0*RO
      T12=D*ZR*ZR
      T11=(GR*D-1.D0)*ZR
      S11=-RO*(WSQ+4.D0*GR-GR*GR*D)
      C11=-T12/RO+1.D0/FLU
      IF(KG.EQ.0) S11=S11-T21*RO
      IF(KG.EQ.0) GOTO 25
      T22=-FL*ZR
      S22=RO*T12
      S12=RO*(T11+T22)
   25 DO 70 K=1,KK
      IF(KG.NE.0) GOTO 30
      H(1,K,NI)=T11*F(1,K)+C11*F(2,K)
      H(2,K,NI)=S11*F(1,K)-T11*F(2,K)
      GOTO 35
   30 H(1,K,NI)=T11*F(1,K)+T12*F(2,K)+C11*F(3,K)
      H(2,K,NI)=T21*F(1,K)+T22*F(2,K)+4.D0*F(4,K)
      H(3,K,NI)=S11*F(1,K)+S12*F(2,K)-T11*F(3,K)-T21*F(4,K)
      H(4,K,NI)=S12*F(1,K)+S22*F(2,K)-T12*F(3,K)-T22*F(4,K)
   35 GO TO (701,702,703,704,705,706,707,708,709,710),NI
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
   40 CONTINUE
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
      IF(I.EQ.JL) RETURN
c
      I=I+JUD
      GO TO 10
c
      END
