      SUBROUTINE DSPLIN(N,X,Y,Q,F)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
      DIMENSION X(1),Y(1),Q(3,1),F(3,1),YY(3)
      EQUIVALENCE (YY(1),Y0)
      DATA YY/3*0.D0/
c
      A0=0.D0
      J2=N-2
      H=X(2)-X(1)
      H2=X(3)-X(1)
      Y0=H*H2*(H2-H)
      H=H*H
      H2=H2*H2
      B0=(Y(1)*(H-H2)+Y(2)*H2-Y(3)*H)/Y0
      B1=B0
      DO 5 I=1,J2
      H=X(I+1)-X(I)
      Y0=Y(I+1)-Y(I)
      H2=H*H
      HA=H-A0
      H2A=H-2.D0*A0
      H3A=2.D0*H-3.*A0
      H2B=H2*B0
      Q(1,I)=H2/HA
      Q(2,I)=-HA/(H2A*H2)
      Q(3,I)=-H*H2A/H3A
      F(1,I)=(Y0-H*B0)/(H*HA)
      F(2,I)=(H2B-Y0*(2.D0*H-A0))/(H*H2*H2A)
      F(3,I)=-(H2B-3.D0*Y0*HA)/(H*H3A)
      A0=Q(3,I)
    5 B0=F(3,I)
      I=J2+1
      H=X(I+1)-X(I)
      Y0=Y(I+1)-Y(I)
      H2=H*H
      HA=H-A0
      H2A=H*HA
      H2B=H2*B0-Y0*(2.D0*H-A0)
      Q(1,I)=H2/HA
      F(1,I)=(Y0-H*B0)/H2A
      HA=X(J2)-X(I+1)
      Y0=-H*HA*(HA+H)
      HA=HA*HA
      Y0=(Y(I+1)*(H2-HA)+Y(I)*HA-Y(J2)*H2)/Y0
      Q(3,I)=(Y0*H2A+H2B)/(H*H2*(H-2.D0*A0))
      Q(2,I)=F(1,I)-Q(1,I)*Q(3,I)
      DO 10 J=1,J2
      K=I-1
      Q(1,I)=F(3,K)-Q(3,K)*Q(2,I)
      Q(3,K)=F(2,K)-Q(2,K)*Q(1,I)
      Q(2,K)=F(1,K)-Q(1,K)*Q(3,K)
   10 I=K
      Q(1,I)=B1
      DO 15 J=1,3
   15 Q(J,N)=YY(J)
c
      RETURN
      END
