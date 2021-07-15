      SUBROUTINE DERMF(IQ,Z,F,FP,IKNT,QFF)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** CALCULATES MINOR VECTOR DERIVATIVE (FP) IN A FLUID ***
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
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
      DIMENSION F(1),FP(1)
c
      IF(IKNT.NE.0) GOTO 14
      T=Z-R(IQ)
      IF(T.NE.0.D0) GOTO 5
      RO=RHO(IQ)
      FLU=FCON(IQ)*QFF
      GR=G(IQ)
      GOTO 10
    5 RO=RHO(IQ)+T*(QRO(1,IQ)+T*(QRO(2,IQ)+T*QRO(3,IQ)))
      FLU=(FCON(IQ)+T*(FSPL(1,IQ)+T*(FSPL(2,IQ)+T*FSPL(3,IQ))))*QFF
      GR=G(IQ)+T*(QG(1,IQ)+T*(QG(2,IQ)+T*QG(3,IQ)))
   10 T21=-4.D0*RO
      ZR=1.D0/Z
      T12=FL3*ZR*ZR/WSQ
      T11=GR*T12-ZR
      S11=RO*(GR*GR*T12-WSQ)+T21*GR*ZR
      C11=-T12/RO+1.D0/FLU
   14 IF(KG.NE.0) GOTO 15
      FP(1)=T11*F(1)+C11*F(2)
      FP(2)=(S11-T21*RO)*F(1)-T11*F(2)
      RETURN
c
   15 IF(IKNT.NE.0) GOTO 19
      T22=-FL*ZR
      S22=RO*T12
      B11=T11+T22
      S12=RO*B11
      IF(IBACK.EQ.1) GOTO 20
      B33=T11-T22
      FP(1)=B11*F(1)+4.D0*F(3)-C11*F(4)
      FP(2)=S12*F(1)-T21*F(3)+T12*F(4)
   19 FP(3)=S22*F(1)-(T12+T12)*F(2)+B33*F(3)+C11*F(5)
      FP(4)=-S11*F(1)+(T21+T21)*F(2)-B33*F(4)-4.D0*F(5)
      FP(5)=-(S12+S12)*F(2)+S11*F(3)-S22*F(4)-B11*F(5)
      RETURN
c
   20 FP(1)=T22*F(1)-T21*F(2)-4.D0*F(3)
      FP(2)=-T12*F(1)+T11*F(2)-C11*F(4)
      FP(3)=-S22*F(1)+S12*F(2)-T22*F(3)+T12*F(4)
      FP(4)=S12*F(1)-S11*F(2)+T21*F(3)-T11*F(4)
c
      RETURN
      END
