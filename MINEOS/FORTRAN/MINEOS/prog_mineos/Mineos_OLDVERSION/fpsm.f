      SUBROUTINE FPSM(LS,NVEFM,ASS)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** COMPUTES SPHEROIDAL MODE START SOLUTION IN A FLUID REGION USING POWER
C*** SERIES OR SPH. BESSEL FNS. IF THE ARGUMENT IS TOO LARGE.
c
c    calls: BFS
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
      DIMENSION ASS(1)
c
      X=R(LS)
      FLA=FLAM(LS)*(1.D0+XLAM(LS)*FCT)
      VPSQ=FLA/RHO(LS)
      ZETA=4.D0*RHO(LS)
      XI=G(LS)/X
      QSQ=(WSQ+FLOAT(KG)*ZETA+XI-FL3*XI*XI/WSQ)/VPSQ
      XSQ=X*X
      ZSQ=QSQ*XSQ
      IF(ZSQ/(4.D0*FL+6.D0).GT.0.1D0) GO TO 10
      F=1.D0
      DT=1.D0
      U=0.D0
      S2=0.D0
      C=0.D0
      D1=FL+FL3*XI/WSQ
      IF(KG.EQ.0) U=-D1/QSQ
      D3=FL2-D1
    5 C=C+2.D0
      C2=1.D0/(C*(FL2+C))
      D1=D1+2.D0
      S=DT*XSQ*C2
      U=U+D1*S
      S2=S2+S
      IF(DABS(DT/F).LT.EPS) GOTO 6
      DT=-DT*ZSQ*C2
      F=F+DT
      GOTO 5
    6 P=ZETA*S2-VPSQ
      S=ZETA*S2*D3-VPSQ*FL2
      GOTO 15
   10 Z=DSQRT(ZSQ)
c      IF(ZSQ.LT.0) then
c        print*,'******** ZSQ NEGATIVE... SKIPPING ********'
c        return
c      ENDIF
c
      CALL BFS(L,Z,EPS,F,FP,FPP)
c
      P=-F*VPSQ
      S=FL2*P
      U=(F*FL-Z*FP)/QSQ
      IF(KG.EQ.0) U=-(FL3*XI*F/WSQ+Z*FP)/QSQ
   15 IF(KG.EQ.0) GOTO 20
      C1=FL*G(LS)-WSQ*X
      C2=FL2*C1*0.25D0/X-RHO(LS)*FL
      ASS(1)=X*FL*P-C1*U
      ASS(2)=-X*FL*F*FLA
      ASS(3)=FL*S*0.25D0-U*C2
      ASS(4)=X*F*FLA*C1
      ASS(5)=-X*F*FLA*C2
      GOTO 25
   20 ASS(1)=U
      ASS(2)=X*F*FLA
   25 SUM=ASS(1)*ASS(1)
      DO 30 I=2,NVEFM
   30 SUM=SUM+ASS(I)*ASS(I)
      SUM=1.D0/DSQRT(SUM)
      IF(ASS(NVEFM).LT.0.D0) SUM=-SUM
      DO 35 I=1,NVEFM
   35 ASS(I)=ASS(I)*SUM
c
      RETURN
      END
