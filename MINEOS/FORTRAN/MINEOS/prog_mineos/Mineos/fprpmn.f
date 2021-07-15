      SUBROUTINE FPRPMN(JF,JL,F,H,NVEFM,IEXP)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** PROPAGATE THE MINOR VECTOR IN A FLUID REGION FROM LEVEL JF TO JL ***
c
c    calls: BAYLIS, DERMF, RKDOT, ZKNT
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
      COMMON/EIFX/AR(14,nknot),INORM(nknot),idum(nknot)
      COMMON/SHANKS/B(46),C(10),DX,STEP(8),STEPF,IN,MAXO
      DIMENSION F(1),H(NVEFM,1),S(5),FP(5)
      DATA ECONST/1048576.D0/
      DATA S/5*0.D0/,FP/5*0.D0/
c
      IF(NVEFM.EQ.1) GOTO 85
      MAXO1=MAXO-1
      JUD=1
      IF(JL.LT.JF) JUD=-1
      Y=R(JF)
      I=JF
      GO TO 45
   10 X=Y
      Y=R(I)
      IF(Y.EQ.X) GOTO 45
      IQ=MIN0(I,I-JUD)
      QFF=1.D0+XLAM(IQ)*FCT
      ZS=DMIN1(X,Y)
      XI=G(I)/Y
      ALFSQ=(WSQ+4.D0*RHO(I)+XI-FL3*XI*XI/WSQ)*RHO(I)/FLAM(I)
      Q=(DSQRT(DABS(ALFSQ-FL3/(X*X)))+1.D0/ZS+FLOAT(KG)*SFL3/X)/STEPF
      DEL=FLOAT(JUD)*STEP(MAXO)/Q
      DXS=0.D0
   15 Y=X+DEL
      IF(FLOAT(JUD)*(Y-R(I)).GT.0.D0) Y=R(I)
      DX=Y-X
c
      IF(DX.NE.DXS) CALL BAYLIS(Q,MAXO1)
c
      DXS=DX
      DO 30 J=1,NVEFM
   30 S(J)=F(J)
      DO 35 NI=1,IN
      Z=X+C(NI)
c
      CALL DERMF(IQ,Z,F,H(1,NI),0,QFF)
c
   35 CALL RKDOT(F,S,H,NVEFM,NI)
c
      IF(KNSW.NE.1) GOTO 40
c
      CALL DERMF(IQ,Y,F,FP,1,QFF)
c
      CALL ZKNT(S,H,F,FP,X,Y,0)
c
   40 X=Y
      IF(Y.NE.R(I)) GO TO 15
   45 SIZE=DABS(F(1))
      DO 50 J=2,NVEFM
   50 SIZE=DMAX1(SIZE,DABS(F(J)))
   55 IF(SIZE.LT.1024.D0) GOTO 65
      DO 60 J=1,NVEFM
   60 F(J)=F(J)/ECONST
      SIZE=SIZE/ECONST
      IEXP=IEXP+20
      GOTO 55
   65 IF(IBACK.EQ.0) GOTO 70
      INORM(I)=INORM(I)+IEXP
      RNE2   =-AR(1,I)*F(4)+AR(4,I)*F(2)+AR(2,I)*F(1)
      AR(1,I)=-AR(1,I)*F(3)+AR(2,I)*F(2)-AR(3,I)*F(1)
      RNE3   =-AR(2,I)*F(4)+AR(4,I)*F(3)-AR(5,I)*F(1)
      AR(4,I)=-AR(3,I)*F(4)-AR(2,I)*F(3)-AR(5,I)*F(2)
      AR(2,I)=RNE2
      AR(3,I)=RNE3
      GOTO 80
   70 INORM(I)=IEXP
      DO 75 J=1,NVEFM
   75 AR(J,I)=F(J)
   80 IF(I.EQ.JL) RETURN

      I=I+JUD
      GO TO 10
   85 DO 90 I=JL,JF
      INORM(I)=INORM(I)+IEXP
   90 CONTINUE
      DO 91 J=1,2
      DO 91 I=JL,JF
      AR(J,I)=AR(J,I)*F(1)
   91 CONTINUE

      RETURN
      END
