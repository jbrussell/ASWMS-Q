      SUBROUTINE TPROP(JF,JL,F)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** PROPAGATES F FROM JF TO JL - TOROIDAL MODES ***
c
c    calls: BAYLIS, RKDOT, TRKNT
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
      COMMON/EIFX/A(14,nknot),adum(nknot)
      COMMON/SHANKS/B(46),C(10),DX,STEP(8),STEPF,IN,MAXO
      DIMENSION H(2,10),S(2),F(2)
      DATA H/20*0.D0/,S/2*0.D0/
c
      FL3M2=FL3-2.D0
      MAXO1=MAXO-1
      Y=R(JF)
      VY=FMU(JF)/RHO(JF)
      I=JF
      GO TO 50
   10 IQ=I
      I=I+1
      X=Y
      Y=R(I)
      IF(Y.EQ.X) GOTO 50
      QLL=1.D0+QSHEAR(IQ)*FCT
      VX=VY
      VY=FMU(I)/RHO(I)
      QX=1.D0/X+DSQRT(DABS(WSQ/(VX)-FL3/(X*X)))
      QY=1.D0/Y+DSQRT(DABS(WSQ/(VY)-FL3/(Y*Y)))
      Q=DMAX1(QX,QY)
      DEL=STEP(MAXO)/Q
      DXS=0.D0
   15 Y=X+DEL
      IF(Y.GT.R(I)) Y=R(I)
      DX=Y-X
c
      IF(DX.NE.DXS) CALL BAYLIS(Q,MAXO1)
c
      DXS=DX
      S(1)=F(1)
      S(2)=F(2)
      DO 40 NI=1,IN
      Z=X+C(NI)
      T=Z-R(IQ)
      RO=RHO(IQ)+T*(QRO(1,IQ)+T*(QRO(2,IQ)+T*QRO(3,IQ)))
      LL=(LCON(IQ)+T*(LSPL(1,IQ)+T*(LSPL(2,IQ)+T*LSPL(3,IQ))))*QLL
      NN=LL
      IF(IFANIS.NE.0) NN=(NCON(IQ)+
     +    T*(NSPL(1,IQ)+T*(NSPL(2,IQ)+T*NSPL(3,IQ))))*QLL
      Z=1.D0/Z
      H(1,NI)=Z*F(1)+F(2)/LL
      H(2,NI)=(NN*FL3M2*Z*Z-RO*WSQ)*F(1)-3.D0*Z*F(2)
c
   40 CALL RKDOT(F,S,H,2,NI)
c
      IF(KNSW.NE.1) GOTO 45
      FP=(NN*FL3M2*Z*Z-RO*WSQ)*F(1)-3.D0*Z*F(2)
c
      CALL TRKNT(S(2),H(2,1),F(2),FP,X,Y)
c
   45 X=Y
      IF(Y.NE.R(I)) GOTO 15
   50 A(1,I)=F(1)
      A(2,I)=F(2)
      IF(I.NE.JL) GO TO 10
c
      RETURN
      END
