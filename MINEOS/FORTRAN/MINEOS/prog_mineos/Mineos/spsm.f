      SUBROUTINE SPSM(LS,NVESM,ASS)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** COMPUTES SPHEROIDAL MODE START SOLNS IN A SOLID REGION USING POWER
C*** SERIES OR SPH. BESSEL FNS. IF THE ARGUMENT IS TOO LARGE.
c
c    calls: BFS
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
      DIMENSION A(6,2),E(15),ASS(1)
      DATA A/12*0.D0/,E/15*0.D0/
c
      X=R(LS)
      XSQ=X*X
      RO=RHO(LS)
      FU=FMU(LS)*(1.D0+QSHEAR(LS)*FCT)
      FLU=FLAM(LS)*(1.D0+XLAM(LS)*FCT)+2.D0*FU
      VSSQ=FU/RO
      VPSQ=FLU/RO
      ZETA=4.D0*RO
      XI=G(LS)/X
      ALFSQ=(WSQ+FLOAT(KG)*ZETA+XI)/VPSQ
      BETASQ=WSQ/VSSQ
      IF(XSQ*BETASQ/(4.D0*FL+6.D0).LT.0.1D0) GOTO 5
      DELSQ=DSQRT((BETASQ-ALFSQ)**2+4.D0*FL3*XI*XI/(VPSQ*VSSQ))
      FKSQ=.5D0*(ALFSQ+BETASQ+DELSQ)
      QSQ=FKSQ-DELSQ
c      print*,'QSQ = ',QSQ
c      print*,'WSQ = ',WSQ
      IF(QSQ.LT.0) then
        print*,'******** QSQ NEGATIVE... SKIPPING ********'
c        QSQ = QSQ*(-1)
        return
      ENDIF
c      print*,'FKSQ ',FKSQ
c      print*,'DELSQ ',DELSQ
      ARG=DSQRT(QSQ)*X
c      print*,'DSQRT(QSQ) ',DSQRT(QSQ)
c      print*,'QSQ ',QSQ
      B=XI/(VSSQ*(BETASQ-QSQ))
      K=1
c
    4 CALL BFS(L,ARG,EPS,F,FP,FPP)
c
      C=ARG*FP/F
      A(1,K)=FL3*B+C
      A(2,K)=1.D0+B+B*C
      A(3,K)=-ARG*ARG
      A(4,K)=B*A(3,K)
      A(5,K)=1.D0
      A(6,K)=FL1-FL3*B
c      print*,'B ',B
c      print*,'C ',C
c      print*,'FL1 ',FL1
c      print*,'FL3 ',FL3
c      print*,'ARG ',ARG
c      print*,'A(1,K) ',A(1,K)
c      print*,'A(2,K) ',A(2,K)
c      print*,'A(3,K) ',A(3,K)
c      print*,'A(4,K) ',A(4,K)
c      print*,'A(5,K) ',A(5,K)
c      print*,'A(6,K) ',A(6,K)
      IF(K.EQ.2) GOTO 25
      ARG=DSQRT(FKSQ)*X
      B=-FLU/(FU*FL3*B)
      K=2
      GOTO 4
    5 D0=1.D0
      H0=-VPSQ/(FL1*VSSQ)
      DO 10 K=1,2
      C=2.D0
      C2=1.D0/(C*(FL2+C))
      D1=C2*(FL3*XI*H0/VPSQ-ALFSQ*D0)*XSQ
      H1=C2*(XI*D0/VSSQ-BETASQ*H0)*XSQ
      U=C2*(FL3*H0+(FL+C)*D0)
      V=C2*(D0+(FL1+C)*H0)
      IF(K.EQ.1.OR.KG.NE.0) GOTO 15
      TE=FL1*VSSQ/(XSQ*(WSQ-XI*FL))
      U=U+FL*TE
      V=V+TE
   15 P=C2*D0
      S=(FL2+C)*P-U
      H=H0+H1
      D=D0+D1
      TE=1.D0
   20 C=C+2.D0
      C2=1.D0/(C*(FL2+C))
      UN=C2*(FL3*H1+(FL+C)*D1)
      U=U+UN
      V=V+C2*(D1+(FL1+C)*H1)
      PN=C2*D1
      P=P+PN
      S=S+(FL2+C)*PN-UN
      IF(TE.LT.EPS) GOTO 21
      D2=C2*(FL3*XI*H1/VPSQ-ALFSQ*D1)*XSQ
      H1=C2*(XI*D1/VSSQ-BETASQ*H1)*XSQ
      D1=D2
      D=D+D1
      H=H+H1
      TE=DMAX1(DABS(D1/D),DABS(H1/H))
      GOTO 20
   21 A(1,K)=U
      A(2,K)=V
      A(3,K)=D
      A(4,K)=H
      A(5,K)=P
      A(6,K)=S
      D0=0.D0
   10 H0=-1.D0
      TE=FL1*VSSQ/(ZETA*XSQ)
      A(5,2)=A(5,2)+TE
      A(6,2)=A(6,2)+TE*FL2
   25 JJ=3+2*KG
      KK=JJ+1
      LL=0
      DO 26 I=1,JJ
      I1=I+1
      DO 26 J=I1,KK
      LL=LL+1
c      print*,'A(I,1) ',A(I,1)
c      print*,'A(I,2) ',A(I,2)
c      print*,'A(J,2) ',A(J,2)
c      print*,'A(J,1) ',A(J,1)
   26 E(LL)=A(I,1)*A(J,2)-A(J,1)*A(I,2)
      IF(KG.NE.0) GOTO 27
      ASS(1)=XSQ*E(1)
      ASS(2)=FU*X*SFL3*(2.D0*E(1)-E(5))
      ASS(3)=FU*X*(E(3)-2.D0*E(1))
      ASS(4)=X*(FLU*E(4)+4.D0*FU*E(1))
      ASS(5)=FU*(FLU*(E(6)+2.D0*E(4))+4.D0*FU*(FL3*(E(5)-E(1))
     +     -E(3)+2.D0*E(1)))
      GOTO 28
   27 C0=WSQ-XI*FL
      C1=RO*FL+0.25D0*FL2*C0
      C2=2.D0*FU/X
      C3=C2*(FL-1.D0)
      ASS(6)=XSQ*(C0*E(1)-ZETA*(FL*E(8)-E(4)))
      ASS(14)=FLU*(FL*E(6)-E(2))
      ASS(13)=FU*SFL3*(FL*E(7)-E(3))
      ASS(1)=X*(C1*E(1)-RO*(FL*E(9)-E(5)))
c      print*,'X ',X
c      print*,'E(1) ',E(1)
c      print*,'C1 ',C1
c      print*,'RO ',RO
c      print*,'FL ',FL
c      print*,'E(9) ',E(9)
c      print*,'E(5) ',E(5)
      ASS(7)=X*FLU*(C0*E(2)-ZETA*FL*E(11))/SFL3+C2*SFL3*ASS(6)
      ASS(8)=X*FU*(C0*E(3)-ZETA*FL*E(13))-C2*ASS(6)
      ASS(12)=(FLU*FL*E(10)+2.D0*(ASS(14)+SFL3*ASS(13)))*FU/X
      ASS(2)=FLU*(C1*E(2)-RO*FL*E(12))/SFL3+C2*SFL3*ASS(1)
      ASS(3)=FU*(C1*E(3)-RO*FL*E(14))-C2*ASS(1)
      ASS(9)=(X*C0*ASS(14)+SFL3*ASS(7)-C3*FL*ASS(6))/FL
      ASS(11)=(SFL3*ASS(12)+C3*(SFL3*ASS(14)-FL*ASS(13)))/FL
      ASS(4)=(C1*ASS(14)+SFL3*ASS(2)-C3*FL*ASS(1))/FL
      ASS(10)=(X*C0*ASS(11)-C3*(SFL3*ASS(9)+FL*ASS(7)))/SFL3
      ASS(5)=(C1*ASS(11)-C3*(SFL3*ASS(4)+FL*ASS(2)))/SFL3
   28 SUM=ASS(1)*ASS(1)
      DO 29 I=2,NVESM
   29 SUM=SUM+ASS(I)*ASS(I)
      SUM=1.D0/DSQRT(SUM)
      IF(ASS(5).LT.0.D0) SUM=-SUM
      DO 30 I=1,NVESM
   30 ASS(I)=ASS(I)*SUM
c
      RETURN
      END
