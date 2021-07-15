      SUBROUTINE INTGDS(RR,IQ,VALS)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** INTERPOLATES INTEGRANDS FOR NORMALISATION,CG,Q ETC..FOR USE WITH GAUSLV.
c
c    calls no other routines
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
      COMMON/EIFX/AR(14,nknot),adum(nknot)
      DIMENSION Q(3),QP(3),VALS(1)
      DATA D1,D2,D3,D4,D5,D6,D7/.111111111111111D0,
     + 0.066666666666667D0,0.666666666666667D0,1.333333333333333D0,
     + 2.666666666666667D0,3.333333333333333D0,5.333333333333333D0/
      DATA Q/3*0.D0/,QP/3*0.D0/
c
      T=RR-R(IQ)
      HN=1.D0/(R(IQ+1)-R(IQ))
      HSQ=HN*HN
      QFF=1.D0+XLAM(IQ)*FCT
      QLL=1.D0+QSHEAR(IQ)*FCT
      IQ1=IQ+1
      IFUN=3
      IF(JCOM.NE.3) IFUN=1
      DO 10 I=1,IFUN
      I2=2*I
      I1=I2-1
      A=((AR(I2,IQ)+AR(I2,IQ1))+2.D0*HN*(AR(I1,IQ)-AR(I1,IQ1)))*HSQ
      B=-(2.D0*AR(I2,IQ)+AR(I2,IQ1))*HN-3.D0*(AR(I1,IQ)-AR(I1,IQ1))*HSQ
      Q(I)=(AR(I1,IQ)+T*(AR(I2,IQ)+T*(B+T*A)))/RR
   10 QP(I)=AR(I2,IQ)+T*(2.D0*B+T*3.D0*A)
      RRO=(RHO(IQ)+T*(QRO(1,IQ)+T*(QRO(2,IQ)+T*QRO(3,IQ))))*RR
      GR=G(IQ)+T*(QG(1,IQ)+T*(QG(2,IQ)+T*QG(3,IQ)))
      FF=(FCON(IQ)+T*(FSPL(1,IQ)+T*(FSPL(2,IQ)+T*FSPL(3,IQ))))*QFF
      LL=(LCON(IQ)+T*(LSPL(1,IQ)+T*(LSPL(2,IQ)+T*LSPL(3,IQ))))*QLL
      IF(IFANIS.NE.0) GOTO 15
      NN=LL
      CC=FF+LL+LL
      AA=CC
      GOTO 20
   15 QAA=1.D0+XA2(IQ)*FCT
      NN=(NCON(IQ)+T*(NSPL(1,IQ)+T*(NSPL(2,IQ)+T*NSPL(3,IQ))))*QLL
      CC=(CCON(IQ)+T*(CSPL(1,IQ)+T*(CSPL(2,IQ)+T*CSPL(3,IQ))))*QAA
      AA=(ACON(IQ)+T*(ASPL(1,IQ)+T*(ASPL(2,IQ)+T*ASPL(3,IQ))))*QAA
   20 QRKA=D1*(4.D0*(AA+FF-NN)+CC)
     1     *(QKAPPA(IQ)+T*HN*(QKAPPA(IQ1)-QKAPPA(IQ)))
      QRMU=D2*(AA+CC-2.D0*FF+5.D0*NN+6.D0*LL)
     1     *(QSHEAR(IQ)+T*HN*(QSHEAR(IQ1)-QSHEAR(IQ)))
      IF(JCOM.NE.3) GOTO 25
      Q1SQ=Q(1)*Q(1)
      Q2SQ=Q(2)*Q(2)
      VALS(1)=RR*RRO*(Q1SQ+Q2SQ)
      FAC=(FL+.5D0)/SFL3
      VALS(2)=(SFL3*(LL*Q1SQ+AA*Q2SQ)+Q(2)*((RRO*GR+2.D0*(NN-AA-LL)+FF)
     +   *Q(1)+RRO*Q(3)-FF*QP(1))+LL*QP(2)*Q(1))*FAC
     +   +.25D0*Q(3)*(QP(3)+FL*Q(3))
      T2=QRKA+D7*QRMU
      T3=QRKA+D4*QRMU
      T4=QRKA+D6*QRMU
      T5=QRKA-D5*QRMU
      T6=QRKA-D3*QRMU
      VALS(3)=.5D0*((FL3*QRMU+T2)*Q1SQ+(2.D0*QRMU+FL3*T3)*Q2SQ)
     1 -Q(1)*SFL3*T4*Q(2)+Q(1)*(T5*QP(1)+SFL3*QRMU*QP(2))+Q(2)*(-2.D0*
     2 QRMU*QP(2)-SFL3*T6*QP(1))+.5D0*(T3*QP(1)*QP(1)+QRMU*QP(2)*QP(2))
      VALS(4)=.5D0*((FL3*LL+4.D0*(RRO*(RRO-GR)+AA-NN-FF)+CC)*Q1SQ+
     +(4.D0*LL-NN-NN+FL3*AA)*Q2SQ +FL*FL*.25D0*Q(3)*Q(3)+CC*QP(1)*QP(1)+
     +LL*QP(2)*QP(2)+.25D0*QP(3)*QP(3))+Q(3)*(RRO*SFL3*Q(2)+FL*.25D0*QP
     +(3))+Q(1)*(SFL3*(RRO*GR+2.D0*(NN-AA-LL)+FF)*Q(2)+RRO*(QP(3)-Q(3))+
     +(FF+FF-CC)*QP(1)+SFL3*LL*QP(2))-Q(2)*(SFL3*FF*QP(1)+(LL+LL)*QP(2))
      RETURN
c
   25 Q(1)=Q(1)*RR
      VALS(1)=RR*RRO*Q(1)*Q(1)
      IF(JCOM.EQ.1) GOTO 30
      VALS(2)=NN*Q(1)*Q(1)
      T1=(RR*QP(1)-Q(1))**2
      T2=(FL3-2.D0)*Q(1)*Q(1)
      VALS(3)=(T1+T2)*QRMU
      VALS(4)=T1*LL+T2*NN
      RETURN
c
   30 T1=(RR*QP(1)+2.D0*Q(1))**2
      T2=D4*(RR*QP(1)-Q(1))**2
      VALS(2)=T1*QRKA+T2*QRMU
      VALS(3)=RR*QP(1)*(CC*RR*QP(1)+4.D0*FF*Q(1))+4.D0*Q(1)*Q(1)
     +    *(AA-NN-RRO*GR)
c
      RETURN
      END
