      SUBROUTINE MODOUT(WCOM,QMOD,GCOM,IOEIG)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c  writes output eigenfrequencies, eigenfunctions to mode file
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
c
      include 'parameter.h'
c
      REAL*4 ABUF,BUF
c            ,WW,QQ,GC
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      COMMON/RINDX/NIC,NOC,NSL,NICP1,NOCP1,NSLP1,N
      COMMON/EIFX/A(14,nknot),DUM(nknot)
      COMMON/BUFCOM/NN,LL,WW,QQ,GC,BUF(nknot6)
      DIMENSION ABUF(maxbyte+3)
      EQUIVALENCE (NN,ABUF)
      DATA A/nknot14*0.D0/,DUM/nknot*0.D0/,BUF/nknot6*0.D0/
c
      NN=NORD
      LL=L
      WW=WCOM
      QQ=QMOD
      GC=GCOM
      if(JCOM.ne.2)goto 5
      nocor=N-NOC
      do 20, i=1,nocor
      buf(i)=A(1,NOC+i)
      j=i+nocor
      buf(j)=A(2,NOC+i)
   20 continue
      nvec=2*nocor+8
      goto 15
    5 NVEC=2*N+8
      IF (JCOM .EQ. 3) NVEC=6*N+8
      DO 10 I=1,N
      BUF(I)=A(1,I)
      J=I+N
      BUF(J)=A(2,I)
      IF(JCOM.NE.3) GOTO 10
      J=J+N
      BUF(J)=A(3,I)
      J=J+N
      BUF(J)=A(4,I)
      J=J+N
      BUF(J)=A(5,I)
      J=J+N
      BUF(J)=A(6,I)
   10 CONTINUE
   15 write(ioeig) (abuf(i),i=1,nvec)
      print *,'!!!!!!!!!!!!!!! nvec = ',nvec
c
      RETURN
      END
