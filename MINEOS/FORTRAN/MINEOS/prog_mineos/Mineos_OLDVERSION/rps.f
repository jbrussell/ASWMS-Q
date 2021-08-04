      SUBROUTINE RPS(I,AA)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** COMPUTES RADIAL MODE START SOLN USING POWER SERIES OR SPH BESSEL
C*** FNS IF ARGUMENT IS TOO LARGE.
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
      DIMENSION AA(2)
c
      X=R(I)
      FLA=FLAM(I)*(1.D0+XLAM(I)*FCT)
      FU=FMU(I)*(1.D0+QSHEAR(I)*FCT)
      AA(1)=X
      Z=RHO(I)*X*X*(WSQ+4.D0*RHO(I)+G(I)/X)/(FLA+2.D0*FU)
      IF(Z.GT.0.1D0) GO TO 2
      D=1.D0
      U=0.D0
      UP=0.D0
      DC=1.D0
      A=0.D0
      B=1.D0
    1 C=B
      A=A+2.D0
      B=B+2.D0
      UT=DC/B
      U=U+UT
      UP=UP+C*UT
      C=1.D0/(A*B)
      DC=-DC*C*Z
      D=D+DC
      IF(DABS(DC/D).GE.EPS) GO TO 1
      AA(2)=(FLA*D+2.D0*FU*UP)/U
      RETURN
c
    2 C=DSQRT(Z)
c
      CALL BFS(L,C,EPS,D,A,B)
c
      U=-A/C
      AA(2)=(FLA*D-2.D0*FU*B)/U
c
      RETURN
      END
