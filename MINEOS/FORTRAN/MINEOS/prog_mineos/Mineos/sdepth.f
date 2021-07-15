      SUBROUTINE SDEPTH(WDIM,LS)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** FINDS STARTING LEVEL,LS, FOR A GIVEN L AND W ***
c
c    calls: STARTL
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
      COMMON/RINDX/NIC,NOC,NSL,NICP1,NOCP1,NSLP1,N
      DATA AW,BW,DW/-2.D-3,2.25D-3,1.28D-3/
c
      Q=0.D0
      W=WDIM/WN
      WSOC=AW+DW*FL
      IF(WDIM.GT.WSOC) GOTO 10
c
      CALL STARTL(NOCP1,NSL,FMU,LS,Q)
c
      IF(LS.EQ.NSL) LS=LS-1
      IF(LS.GT.NOCP1) RETURN
c
   10 WSIC=AW+BW*FL
      IF(WDIM.GT.WSIC) GOTO 20
c
      CALL STARTL(NICP1,NOC,FLAM,LS,Q)
c
      IF(LS.EQ.NOC) LS=LS-1
      IF(LS.GT.NICP1) RETURN
c
   20 CALL STARTL(2,NIC,FMU,LS,Q)
c
      IF(LS.EQ.NIC) LS=LS-1
c
      RETURN
      END
