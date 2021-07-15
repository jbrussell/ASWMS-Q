      SUBROUTINE SPRPMN(JF,JL,F,H,NVESM,IEXP)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** PROPAGATE A MINOR VECTOR IN A SOLID REGION FROM LEVEL JF TO JL ***
c
c    calls: BAYLIS, DERMS, RKDOT, ZKNT
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
      DIMENSION F(1),H(NVESM,1),S(14),FP(14),RNE(6)
      DATA ECONST/1048576.D0/
      DATA S/14*0.D0/,FP/14*0.D0/,RNE/6*0.D0/
c
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
      QLL=1.D0+QSHEAR(IQ)*FCT
      QAA=1.D0+XA2(IQ)*FCT
      ZS=DMIN1(X,Y)
      XI=G(I)/Y
      VPSQ=(FLAM(I)+2.D0*FMU(I))/RHO(I)
      VSSQ=FMU(I)/RHO(I)
      ALFSQ=(WSQ+4.D0*RHO(I)+XI)/VPSQ
      BETASQ=WSQ/VSSQ
      DELSQ=DSQRT((BETASQ-ALFSQ)**2+4.D0*FL3*XI*XI/(VSSQ*VPSQ))
      FKSQ=.5D0*(ALFSQ+BETASQ+DELSQ)-FL3/(X*X)
      QT=DSQRT(DABS(FKSQ))+DSQRT(DABS(FKSQ-DELSQ))+2.D0/ZS
      Q=(QT+FLOAT(KG)*SFL3/X)/STEPF
      DEL=FLOAT(JUD)*STEP(MAXO)/Q
      DXS=0.D0
   15 Y=X+DEL
      IF(FLOAT(JUD)*(Y-R(I)).GT.0.D0) Y=R(I)
      DX=Y-X
c
      IF(DX.NE.DXS) CALL BAYLIS(Q,MAXO1)
c
      DXS=DX
      DO 30 J=1,NVESM
   30 S(J)=F(J)
      DO 35 NI=1,IN
      Z=X+C(NI)
c
      CALL DERMS(IQ,Z,F,H(1,NI),0,QFF,QLL,QAA)
c
   35 CALL RKDOT(F,S,H,NVESM,NI)
c
      IF(KNSW.NE.1) GOTO 40
c
      CALL DERMS(IQ,Y,F,FP,1,QFF,QLL,QAA)
c
      CALL ZKNT(S,H,F,FP,X,Y,1)
c
   40 X=Y
      IF(Y.NE.R(I)) GOTO 15
   45 SIZE=DABS(F(1))
      DO 50 J=2,NVESM
   50 SIZE=DMAX1(SIZE,DABS(F(J)))
   55 IF(SIZE.LT.1024.D0) GOTO 65
c JBR - Check that SIZE /= NaN      
      IF(SIZE.NE.SIZE) then
        print*,'SIZE = ',SIZE
        print*,'ERROR sprpmn.f: SIZE = NaN'
        print*,'check spsm.f : QSQ < 0'
c        print*,' i quit'
        stop
      ENDIF
      DO 60 J=1,NVESM
   60 F(J)=F(J)/ECONST
      SIZE=SIZE/ECONST
      IEXP=IEXP+20
      GOTO 55
   65 IF(IBACK.EQ.0) GOTO 85
      INORM(I)=INORM(I)+IEXP
      IF(KG.EQ.0) GOTO 70
      T1=F(4)+F(8)
      T2=T1+F(4)
      T1=T1+F(8)
      T3=F(8)-F(4)
      RNE(1)=AR(6,I)*F(10)-AR(14,I)*F(9)+AR(13,I)*T3
     1      -AR(1,I)*F(7)-AR(7,I)*F(6)+AR(8,I)*F(5)
     2      +AR(12,I)*F(3)-AR(2,I)*F(2)+AR(3,I)*F(1)
      RNE(2)=AR(6,I)*F(13)+AR(14,I)*T2+AR(13,I)*F(12)
     1      -AR(1,I)*F(11)-AR(9,I)*F(6)-AR(7,I)*F(5)
     2      +AR(11,I)*F(3)-AR(4,I)*F(2)-AR(2,I)*F(1)
      RNE(3)=AR(6,I)*F(14)-AR(7,I)*T1-AR(8,I)*F(12)
     1      +AR(13,I)*F(11)-AR(9,I)*F(9)+AR(14,I)*F(7)
     2      +AR(10,I)*F(3)+AR(11,I)*F(2)+AR(12,I)*F(1)
      RNE(4)=AR(14,I)*F(14)+AR(7,I)*F(13)+AR(12,I)*F(12)
     1      -AR(2,I)*F(11)-AR(9,I)*F(10)-AR(11,I)*T3
     2      +AR(4,I)*F(7)+AR(10,I)*F(5)+AR(5,I)*F(1)
      RNE(5)=AR(13,I)*F(14)+AR(8,I)*F(13)-AR(12,I)*T2
     1      -AR(3,I)*F(11)+AR(7,I)*F(10)-AR(11,I)*F(9)
     2      -AR(2,I)*F(7)+AR(10,I)*F(6)+AR(5,I)*F(2)
      RNE(6)=AR(1,I)*F(14)+AR(13,I)*F(13)-AR(2,I)*T1
     1      -AR(3,I)*F(12)+AR(14,I)*F(10)-AR(4,I)*F(9)
     2      -AR(11,I)*F(6)-AR(12,I)*F(5)+AR(5,I)*F(3)
      GOTO 75
   70 RNE(1)=-AR(1,I)*F(3)+AR(2,I)*F(2)-AR(3,I)*F(1)
      RNE(2)=-AR(1,I)*F(4)+AR(4,I)*F(2)+AR(2,I)*F(1)
      RNE(3)=-AR(2,I)*F(4)+AR(4,I)*F(3)-AR(5,I)*F(1)
      RNE(4)=-AR(3,I)*F(4)-AR(2,I)*F(3)-AR(5,I)*F(2)
   75 DO 80 JJ=1,6
   80 AR(JJ,I)=RNE(JJ)
      GOTO 95
   85 INORM(I)=IEXP
      DO 90 J=1,NVESM
   90 AR(J,I)=F(J)
c
   95 IF(I.EQ.JL) RETURN
c
      I=I+JUD
      GO TO 10
c
      END
