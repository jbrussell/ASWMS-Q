      SUBROUTINE EIFOUT(LSMIN)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** MASSAGES SPHEROIDAL MODE EIGENFUNCTIONS BEFORE OUTPUT ***
c
c    calls: GAUSLV
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
c
      include 'parameter.h'
c
      REAL*8 LL,LCON,NCON,LSPL,NSPL
      COMMON R(nknot),FMU(nknot),FLAM(nknot),QSHEAR(nknot),
     &       QKAPPA(nknot),XA2(nknot),XLAM(nknot),RHO(nknot),
     &       QRO(3,nknot),G(nknot),QG(3,nknot),FCON(nknot),
     &       FSPL(3,nknot),LCON(nknot),LSPL(3,nknot),NCON(nknot),
     &       NSPL(3,nknot),CCON(nknot),CSPL(3,nknot),ACON(nknot),
     &       ASPL(3,nknot)
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      COMMON/EIFX/A(14,nknot),INORM(nknot),idum(nknot)
      COMMON/RINDX/NIC,NOC,NSL,NICP1,NOCP1,NSLP1,N
      DIMENSION ZI(4)
      DATA ZI/4*0.D0/
c
      I1=MIN0(NIC,MAX0(2,LSMIN))
      I2=NIC
    5 IF(I1.EQ.I2) GOTO 20
      DO 10 IQ=I1,I2
      FF=FCON(IQ)*(1.D0+XLAM(IQ)*FCT)
      LL=LCON(IQ)*(1.D0+QSHEAR(IQ)*FCT)
      ZR=1.D0/R(IQ)
      SFL3Z=SFL3*ZR
      D=1.D0/(CCON(IQ)*(1.D0+XA2(IQ)*FCT))
      V=A(2,IQ)
      IF(KG.NE.0) GOTO 15
      A(2,IQ)=(ZR-2.D0*FF*D*ZR)*A(1,IQ)+SFL3Z*FF*D*V+D*A(3,IQ)
      A(4,IQ)=-SFL3Z*A(1,IQ)+(ZR+ZR)*V+A(4,IQ)/LL
      A(5,IQ)=0.D0
      A(6,IQ)=0.D0
      GOTO 10
   15 A(2,IQ)=(ZR-2.D0*FF*D*ZR)*A(1,IQ)+SFL3Z*FF*D*V+D*A(4,IQ)
      A(4,IQ)=-SFL3Z*A(1,IQ)+(ZR+ZR)*V+A(5,IQ)/LL
      A(5,IQ)=A(3,IQ)
      A(6,IQ)=4.D0*(A(6,IQ)-RHO(IQ)*A(1,IQ))-FL*ZR*A(5,IQ)
   10 A(3,IQ)=V
   20 IF(I2.EQ.NSL) GOTO 25
      I1=MIN0(NSL,MAX0(LSMIN,NOCP1))
      I2=NSL
      GOTO 5
   25 I1=MIN0(NOC,MAX0(LSMIN,NICP1))
      I2=NOC
   30 IF(I1.EQ.I2) GOTO 50
      DO 35 IQ=I1,I2
      ZR=1.D0/R(IQ)
      SFL3Z=SFL3*ZR
      FFI=1.D0/(FLAM(IQ)*(1.D0+XLAM(IQ)*FCT))
      IF(KG.NE.0) GOTO 40
      P=A(2,IQ)
      A(5,IQ)=0.D0
      A(6,IQ)=0.D0
      GOTO 45
   40 P=A(3,IQ)
      A(5,IQ)=A(2,IQ)
      A(6,IQ)=4.D0*(A(4,IQ)-RHO(IQ)*A(1,IQ))-FL*ZR*A(5,IQ)
   45 A(3,IQ)=SFL3Z*(G(IQ)*A(1,IQ)-P/RHO(IQ)+A(5,IQ))/WSQ
      A(2,IQ)=SFL3Z*A(3,IQ)-A(1,IQ)*ZR+P*FFI
   35 A(4,IQ)=SFL3Z*(A(1,IQ)+P*(QRO(1,IQ)/(RHO(IQ)**2)+G(IQ)*FFI)/WSQ)
   50 IF(N.EQ.NSL.OR.I2.EQ.N) GOTO 55
      I1=NSLP1
      I2=N
      GOTO 30
   55 IMAX=0
      DO 60 IQ=LSMIN,N
   60 IMAX=MAX0(INORM(IQ),IMAX)
      DO 65 IQ=LSMIN,N
      IEXP=INORM(IQ)-IMAX
      AL=0.D0
      IF(IEXP.GE.-80) AL=2.D0**IEXP
      DO 65 J=1,6
   65 A(J,IQ)=A(J,IQ)*AL
      LSM1=MAX0(1,LSMIN-1)
      DO 70 I=1,LSM1
      DO 70 J=1,6
   70 A(J,I)=0.D0
      IF(L.GT.1.OR.LSMIN.GT.2) GOTO 75
      A(2,1)=1.5D0*A(1,2)/R(2)-.5D0*A(2,2)
      A(4,1)=1.5D0*A(3,2)/R(2)-.5D0*A(4,2)
   75 DO 80 J=1,4
   80 ZI(J)=0.D0
      I1=MAX0(LSMIN,2)
      DO 85 IQ=I1,N
      IP=IQ-1
c
   85 IF(R(IQ).NE.R(IP)) CALL GAUSLV(R(IP),R(IQ),IP,ZI,4)
c
c     print*, (zi(i),i=1,4)      
      CG=ZI(2)/(W*ZI(1))
      WRAY=DSQRT(2.D0*ZI(4)/ZI(1))
      QINV=2.D0*ZI(3)/(WSQ*ZI(1))
      RNORM=1.D0/(W*DSQRT(ZI(1)))
      DO 90 IQ=I1,N
      ZR=1.D0/R(IQ)
      A(1,IQ)=A(1,IQ)*ZR
      A(2,IQ)=(A(2,IQ)-A(1,IQ))*ZR
      A(3,IQ)=A(3,IQ)*ZR
      A(4,IQ)=(A(4,IQ)-A(3,IQ))*ZR
      A(5,IQ)=A(5,IQ)*ZR
      A(6,IQ)=(A(6,IQ)-A(5,IQ))*ZR
      A(1,IQ)=A(1,IQ)*RNORM
      A(2,IQ)=A(2,IQ)*RNORM
      A(3,IQ)=A(3,IQ)*RNORM
      A(4,IQ)=A(4,IQ)*RNORM
      A(5,IQ)=A(5,IQ)*RNORM
   90 A(6,IQ)=A(6,IQ)*RNORM
      IF(LSMIN.GT.2.OR.L.GT.2) RETURN
c
      IF(L.EQ.2) GOTO 95
      A(1,1)=A(1,2)-.5D0*A(2,2)*R(2)
      A(2,1)=0.D0
      A(3,1)=A(3,2)-.5D0*A(4,2)*R(2)
      A(4,1)=0.D0
      A(6,1)=1.5D0*A(5,2)/R(2)-.5D0*A(6,2)
      RETURN
c
   95 A(2,1)=1.5D0*A(1,2)/R(2)-.5D0*A(2,2)
      A(4,1)=1.5D0*A(3,2)/R(2)-.5D0*A(4,2)
c
      RETURN
      END
