      SUBROUTINE DERMS(IQ,Z,F,FP,IKNT,QFF,QLL,QAA)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** CALCULATES MINOR VECTOR DERIVATIVE (FP) IN A SOLID ***
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'parameter.h'
c
      REAL*8 NN,LL,LCON,NCON,LSPL,NSPL
      COMMON R(nknot),FMU(nknot),FLAM(nknot),QSHEAR(nknot),
     &       QKAPPA(nknot),XA2(nknot),XLAM(nknot),RHO(nknot),
     &       QRO(3,nknot),G(nknot),QG(3,nknot),FCON(nknot),
     &       FSPL(3,nknot),LCON(nknot),LSPL(3,nknot),NCON(nknot),
     &       NSPL(3,nknot),CCON(nknot),CSPL(3,nknot),ACON(nknot),
     &       ASPL(3,nknot)
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      DIMENSION F(1),FP(1)
c
      IF(IKNT.NE.0) GOTO 19
      T=Z-R(IQ)
      IF(T.NE.0.D0) GOTO 5
      RO=RHO(IQ)
      GR=G(IQ)
      FF=FCON(IQ)*QFF
      LL=LCON(IQ)*QLL
      NN=NCON(IQ)*QLL
      CC=CCON(IQ)*QAA
      AA=ACON(IQ)*QAA
      GOTO 15
    5 RO=RHO(IQ)+T*(QRO(1,IQ)+T*(QRO(2,IQ)+T*QRO(3,IQ)))
      GR=G(IQ)+T*(QG(1,IQ)+T*(QG(2,IQ)+T*QG(3,IQ)))
      FF=(FCON(IQ)+T*(FSPL(1,IQ)+T*(FSPL(2,IQ)+T*FSPL(3,IQ))))*QFF
      LL=(LCON(IQ)+T*(LSPL(1,IQ)+T*(LSPL(2,IQ)+T*LSPL(3,IQ))))*QLL
      IF(IFANIS.NE.0) GOTO 10
      NN=LL
      CC=FF+LL+LL
      AA=CC
      GOTO 15
   10 NN=(NCON(IQ)+T*(NSPL(1,IQ)+T*(NSPL(2,IQ)+T*NSPL(3,IQ))))*QLL
      CC=(CCON(IQ)+T*(CSPL(1,IQ)+T*(CSPL(2,IQ)+T*CSPL(3,IQ))))*QAA
      AA=(ACON(IQ)+T*(ASPL(1,IQ)+T*(ASPL(2,IQ)+T*ASPL(3,IQ))))*QAA
   15 ZR=1.D0/Z
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
      IF(KG.NE.0) GOTO 25
      S11=S11+4.D0*RO*RO
      IF(IBACK.EQ.1) GOTO 20
      B11=T11+T22
      B33=T11-T22
      FP(1)=B11*F(1)+C22*F(3)-C11*F(4)
      FP(2)=S12*F(1)-T21*F(3)+T12*F(4)
   19 IF(KG.NE.0) GOTO 29
      FP(3)=S22*F(1)-2.D0*T12*F(2)+B33*F(3)+C11*F(5)
      FP(4)=-S11*F(1)+2.D0*T21*F(2)-B33*F(4)-C22*F(5)
      FP(5)=-2.D0*S12*F(2)+S11*F(3)-S22*F(4)-B11*F(5)
      RETURN
c
   20 FP(1)=T22*F(1)-T21*F(2)-C22*F(3)
      FP(2)=-T12*F(1)+T11*F(2)-C11*F(4)
      FP(3)=-S22*F(1)+S12*F(2)-T22*F(3)+T12*F(4)
      FP(4)=S12*F(1)-S11*F(2)+T21*F(3)-T11*F(4)
      RETURN
c
   25 T31=-4.D0*RO
      T33=-FL*ZR
      S13=-FL1*ZR*RO
      S23=RO*SFL3Z
      IF(IBACK.EQ.1) GOTO 30
      B11=T11+T22-T33
      B33=T11-T22-T33
      B44=T22-T11-T33
      B55=-T11-T22-T33
      B32=-T12-T12
      B42=T21+T21
      B52=-S12-S12
      B313=-S23-S23
      B414=S13+S13
      B914=T31+T31
      FP(1)=B11*F(1)+C22*F(3)-C11*F(4)
      FP(2)=S12*F(1)-T33*F(2)-T21*F(3)+T12*F(4)-S13*F(13)-S23*F(14)
      FP(6)=4.D0*F(1)-B55*F(6)+C22*F(8)-C11*F(9)
      FP(7)=4.D0*F(2)+S12*F(6)+T33*F(7)-T21*F(8)+T12*F(9)-T31*F(13)
      FP(8)=4.D0*F(3)+S22*F(6)+B32*F(7)-B44*F(8)+C11*F(10)
      FP(9)=4.D0*F(4)-S11*F(6)+B42*F(7)-B33*F(9)-C22*F(10)+B914*F(14)
      FP(10)=4.D0*F(5)+B52*F(7)+S11*F(8)-S22*F(9)-B11*F(10)+B914*F(12)
      FP(11)=-T31*F(2)+S13*F(7)+S23*F(9)-T11*F(11)+T21*F(12)
     +      -S11*F(13)+S12*F(14)
      FP(12)=T31*F(3)+S23*F(7)-S13*F(8)+T12*F(11)-T22*F(12)
     +      +S12*F(13)-S22*F(14)
      FP(13)=S23*F(6)-C11*F(11)+T11*F(13)-T12*F(14)
      FP(14)=-T31*F(1)+S13*F(6)-C22*F(12)-T21*F(13)+T22*F(14)
   29 FP(3)=S22*F(1)+B32*F(2)+B33*F(3)+C11*F(5)+B313*F(13)
      FP(4)=-S11*F(1)+B42*F(2)+B44*F(4)-C22*F(5)+B414*F(14)
      FP(5)=B52*F(2)+S11*F(3)-S22*F(4)+B55*F(5)-B313*F(11)+B414*F(12)
      RETURN
c
   30 B11=T22+T33
      B22=T11+T33
      B33=T11+T22
      B55=T22-T33
      B66=T11-T33
      B99=T11-T22
      T4=F(4)+F(8)
      T5=T4+F(8)
      T4=T4+F(4)
      FP(1)=B11*F(1)-T21*F(2)-T31*F(3)-4.D0*F(5)+C22*F(7)
      FP(2)=-T12*F(1)+B22*F(2)-4.D0*F(6)+C11*F(11)
      FP(3)=B33*F(3)-C22*F(9)+C11*F(12)
      FP(4)=-S23*F(1)+S13*F(2)+T31*F(6)
      FP(5)=S13*F(3)+B55*F(5)-T21*F(6)-C22*F(10)
      FP(6)=S23*F(3)-T12*F(5)+B66*F(6)-C11*F(13)
      FP(7)=S22*F(1)-S12*F(2)-B55*F(7)+T31*F(9)+4.D0*F(10)+T12*F(11)
      FP(8)=S23*F(1)-S12*F(3)-T21*F(9)+T12*F(12)
      FP(9)=S23*F(2)-S22*F(3)-T12*T5+B99*F(9)-C11*F(14)
      FP(10)=S23*(F(4)-F(8))-S22*F(5)+S12*F(6)+S13*F(9)-B11*F(10)
     1      +T12*F(13)
      FP(11)=-S12*F(1)+S11*F(2)-T4*T31+T21*F(7)-B66*F(11)+4.D0*F(13)
      FP(12)=-S13*F(1)+S11*F(3)+T21*T5-T31*F(5)-B99*F(12)+C22*F(14)
      FP(13)=-T4*S13+S12*F(5)-S11*F(6)+T21*F(10)-S23*F(12)-B22*F(13)
      FP(14)=S12*T5-S13*F(7)-S11*F(9)+T31*F(10)-S23*F(11)+S22*F(12)
     1      -B33*F(14)
c
      RETURN
      END
