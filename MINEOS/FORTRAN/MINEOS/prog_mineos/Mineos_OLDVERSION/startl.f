      SUBROUTINE STARTL(JF,JL,V,LS,Q)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** FINDS START LEVEL BETWEEN JF AND JL USING VELOCITYV AND ANG. ORD. L.
C*** UPON ENTRY Q IS THE VALUE OF THE EXPONENT AT R(JF) OR AT THE TURNING
C*** POINT(Q=0) DEPENDING ON PREVIOUS CALLS TO STARTL. UPON EXIT Q IS THE
C*** VALUE OF THE EXPONENT AT THE STARTING LEVEL LS.
c
c    calls no other routines
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
      DIMENSION RRLOG(nknot),P(nknot),V(1)
      DATA IFIRST/1/
      DATA RRLOG/nknot*0.D0/,P/nknot*0.D0/
c
      IF(IFIRST.NE.1) GOTO 5
      IFIRST=0
      VERTNO=-DLOG(EPS)
      DO 1 I=3,N
    1 RRLOG(I)=.5D0*DLOG(R(I)/R(I-1))
    5 DO 10 J=JF,JL
      PP=FL3-WSQ*R(J)*R(J)*RHO(J)/V(J)
      IF(PP.LE.0.D0) GOTO 15
   10 P(J)=DSQRT(PP)
   15 P(J)=0.D0
   20 K=J
      J=J-1
      IF(J.LE.JF) GO TO 25
      Q=Q+RRLOG(K)*(P(J)+P(K))
      IF(Q.LT.VERTNO) GO TO 20
      LS=J
      RETURN
c
   25 LS=JF
c
      RETURN
      END
