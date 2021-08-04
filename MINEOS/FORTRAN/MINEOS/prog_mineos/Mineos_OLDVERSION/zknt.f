      SUBROUTINE ZKNT(S,SP,F,FP,X,Y,IFSOL)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** GIVEN MINOR VECTOR AND DERIVS,CONSTRUCTS MODE COUNT ***
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      DIMENSION S(1),SP(1),F(1),FP(1),XS(4),VAL(4)
      DATA XS/4*0.D0/,VAL/4*0.D0/
c
      IF(IFSOL.EQ.0.AND.KG.EQ.0) GOTO 5
      Y1=S(5)
      Y2=F(5)
      Y1P=SP(5)
      Y2P=FP(5)
      T1=S(3)-S(4)
      T2=F(3)-F(4)
      T1P=SP(3)-SP(4)
      T2P=FP(3)-FP(4)
      GOTO 10
    5 Y1=S(2)
      Y2=F(2)
      Y1P=SP(2)
      Y2P=FP(2)
      T1=S(1)
      T2=F(1)
      T1P=SP(1)
      T2P=FP(1)
   10 H=Y-X
      NS=0
      IF(KOUNT.NE.0) GOTO 15
      A1=Y2-Y1
      A2=0.D0
      A3=0.D0
      A22=0.D0
      A33=0.D0
      GOTO 50
   15 A1=H*Y1P
      A2=-H*(2.D0*Y1P+Y2P)+3.D0*(Y2-Y1)
      A3=H*(Y1P+Y2P)-2.D0*(Y2-Y1)
      A33=3.D0*A3
      A22=2.D0*A2
      IF(A3.NE.0.D0) GOTO 20
      IF(A2.EQ.0.D0) GOTO 50
      XS(2)=-A1/A22
      IF(XS(2).GE.0.D0.AND.XS(2).LE.1.D0) NS=1
      GOTO 50
   20 DISC=A2*A2-A1*A33
      IF(DISC) 50,25,30
   25 XS(2)=-A2/A33
      IF(XS(2).GE.0.D0.AND.XS(2).LE.1.D0) NS=1
      GOTO 50
   30 DISC=DSQRT(DISC)
      TR1=(-A2+DISC)/A33
      TR2=(-A2-DISC)/A33
      IF(DABS(A33).GT.DABS(A1)) GOTO 35
      FAC=A1/A33
      TR1=FAC/TR1
      TR2=FAC/TR2
   35 IF(TR1.LT.0.D0.OR.TR1.GT.1.D0) GOTO 40
      XS(2)=TR1
      NS=1
   40 IF(TR2.LT.0.D0.OR.TR2.GT.1.D0) GOTO 50
      NS=NS+1
      XS(NS+1)=TR2
      IF(NS.LT.2) GOTO 50
      IF(TR2.GE.TR1) GOTO 50
      XS(2)=TR2
      XS(3)=TR1
   50 VAL(1)=Y1
      XS(1)=0.D0
      NS2=NS+2
      VAL(NS2)=Y2
      XS(NS2)=1.D0
      IF(NS.EQ.0) GOTO 60
      NS1=NS+1
      DO 55 J=2,NS1
      T=XS(J)
   55 VAL(J)=Y1+T*(A1+T*(A2+T*A3))
   60 IFT=0
      DO 100 J=2,NS2
      IF(VAL(J-1)*VAL(J).GT.0.D0) GOTO 100
      IF(VAL(J-1).NE.0.D0) GOTO 65
      TES=T1*A1
      GOTO 90
   65 RT1=0.5D0*(XS(J-1)+XS(J))
      RT=RT1
      DO 70 I=1,5
      V=Y1+RT*(A1+RT*(A2+RT*A3))
      VP=A1+RT*(A22+RT*A33)
      ADD=-V/VP
      RT=RT+ADD
      IF(DABS(ADD).LT.1.D-5) GOTO 75
      IF(DABS(RT-RT1).LE..5D0) GOTO 70
      RT=RT1
      GOTO 75
   70 CONTINUE
   75 IF(IFT.NE.0) GOTO 85
      IF(KOUNT.NE.0) GOTO 80
      B1=T2-T1
      B2=0.D0
      B3=0.D0
      GOTO 85
   80 B1=H*T1P
      B2=-H*(2.D0*T1P+T2P)+3.D0*(T2-T1)
      B3=H*(T1P+T2P)-2.D0*(T2-T1)
      IFT=1
   85 TES=T1+RT*(B1+RT*(B2+RT*B3))
      VP=A1+RT*(A22+RT*A33)
      TES=TES*VP
   90 IF(TES.LT.0.D0) KOUNT=1+KOUNT
      IF(TES.GT.0.D0) KOUNT=KOUNT-1
  100 CONTINUE
c
      RETURN
      END
