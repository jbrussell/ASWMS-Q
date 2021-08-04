      SUBROUTINE TPS(I,A)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** COMPUTES TOROIDAL MODE START SOLN USING POWER SERIES OR SPH BESSEL
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
      DIMENSION A(2)
c
      A(1)=R(I)
      FU=FMU(I)*(1.D0+QSHEAR(I)*FCT)
      X2=WSQ*R(I)*R(I)*RHO(I)/FU
      TE=X2/(4.D0*FL+6.D0)
      IF(TE.GT.0.1D0) GO TO 10
      C=0.D0
      D=1.D0
      E=FL-1.D0
      F=D
      H=E
    5 C=C+2.D0
      D=-D*X2/(C*(FL2+C))
      E=E+2.D0
      F=F+D
      H=H+D*E
      IF(DABS(D/F).GT.EPS) GO TO 5
      A(2)=FU*H/F
      RETURN
c
   10 X=DSQRT(X2)
c
      CALL BFS(L,X,EPS,F,H,D)
c
      A(2)=FU*(X*H-F)/F
c
      RETURN
      END
