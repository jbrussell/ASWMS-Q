      SUBROUTINE SVD(A,MROW,NCOL)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C    SECTION I CHAPTER 10 WILKENSON AND REINSCH (1971 ,SPRINGER).
C    THE MATRIX A IS OVERWRITTEN WITH V(NCOL,NCOL), THE RIGHT SIDE ORTHOGONAL
C    MATRIX IN THE SVD DECOMPOSITION. FOR USE ONLY IN EOS SUBS AS ,TO REDUCE
C    BRANCHING POINTS, I HAVE USED THE FACT THAT NCOL IS LT MROW.
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
      DIMENSION A(6,1),E(3),Q(3)
      DATA E/3*0.D0/,Q/3*0.D0/
c
      EPS=1.5D-14
      TOL=1.D-293
      G=0.D0
      X=0.D0
      DO 60 I=1,NCOL
      L=I+1      
      E(I)=G
      S=0.D0
      DO 10 J=I,MROW
   10 S=S+A(J,I)*A(J,I)
      IF(S.GT.TOL) GO TO 15
      Q(I)=0.D0
      IF(L.GT.NCOL) GOTO 60
      GO TO 30
   15 Q(I)=DSIGN(DSQRT(S),-A(I,I))
      H=A(I,I)*Q(I)-S
      A(I,I)=A(I,I)-Q(I)
      IF(L.GT.NCOL) GO TO 60
      DO 25 J=L,NCOL
      S=0.D0
      DO 20 K=I,MROW
   20 S=S+A(K,I)*A(K,J)
      F=S/H
      DO 25 K=I,MROW
   25 A(K,J)=A(K,J)+F*A(K,I)
   30 S=0.D0
      DO 35 J=L,NCOL
   35 S=S+A(I,J)*A(I,J)
      IF(S.GE.TOL)GO TO 40
      G=0.D0
      GO TO 60
   40 G=DSIGN(DSQRT(S),-A(I,L))
      H=A(I,L)*G-S
      A(I,L)=A(I,L)-G
      DO 45 J=L,NCOL
   45 E(J)=A(I,J)/H
      DO 55 J=L,MROW
      S=0.D0
      DO 50 K=L,NCOL
   50 S=S+A(J,K)*A(I,K)
      DO 55 K=L,NCOL
   55 A(J,K)=A(J,K)+S*E(K)
   60 X=DMAX1(DABS(Q(I))+DABS(E(I)),X)
      GOTO 100
   75 IF(G.EQ.0.D0)GO TO 91
      H=A(I,L)*G
      DO 80 J=L,NCOL
   80 A(J,I)=A(I,J)/H
      DO 90 J=L,NCOL
      S=0.D0
      DO 85 K=L,NCOL
   85 S=S+A(I,K)*A(K,J)
      DO 90 K=L,NCOL
   90 A(K,J)=A(K,J)+S*A(K,I)
   91 DO 95 J=L,NCOL
      A(I,J)=0.D0
   95 CONTINUE
      DO 96 J=L,NCOL
      A(J,I)=0.D0
   96 CONTINUE
c
c correction to prevent accessing arrays out of bounds
c
  100 if (i .eq. (ncol+1)) then
        i = ncol
      endif
c
      A(I,I)=1.D0
      G=E(I)
      L=I
      I=I-1
      IF(I.GE.1)GO TO 75
      EP=EPS*X
      K=NCOL
  105 L=K
  110 IF(DABS(E(L)).LE.EP)GO TO 125
      IF(DABS(Q(L-1)).LE.EP) GO TO 115
      L=L-1
      IF(L.GE.1)GO TO 110
  115 C=0.D0
      S=1.D0
      DO 120 I=L,K
      F=S*E(I)
      E(I)=C*E(I)
      IF(DABS(F).LE.EP)GO TO 125
      G=Q(I)
      H=DSQRT(F*F+G*G)
      C=G/H
      S=-F/H
  120 Q(I)=H
  125 Z=Q(K)
      IF(L.EQ.K)GO TO 145
      X=Q(L)
      Y=Q(K-1)
      G=E(K-1)
      H=E(K)
      F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.D0*H*Y)
      G=DSQRT(F*F+1.D0)
      F=((X-Z)*(X+Z)+H*(Y/(F+DSIGN(G,F))-H))/X
      C=1.D0
      S=1.D0
      LP1=L+1
      DO 140 I=LP1,K
      G=E(I)
      Y=Q(I)
      H=S*G
      G=C*G
      Z=DSQRT(F*F+H*H)
      IM1=I-1
      E(IM1)=Z
      C=F/Z
      S=H/Z
      F=S*G+C*X
      G=C*G-S*X
      H=S*Y
      Y=C*Y
      DO 130 J=1,NCOL
      X=A(J,IM1)
      Z=A(J,I)
      A(J,IM1)=C*X+S*Z
  130 A(J,I)=C*Z-S*X
      Z=DSQRT(F*F+H*H)
      Q(IM1)=Z
      C=F/Z
      S=H/Z
      F=S*Y+C*G
  140 X=C*Y-S*G
      E(L)=0.D0
      E(K)=F
      Q(K)=X
      GO TO 105
  145 IF(Z.GE.0.D0)GO TO 155
      Q(K)=-Z
      DO 150 J=1,NCOL
  150 A(J,K)=-A(J,K)
  155 K=K-1
      IF(K.GE.1)GO TO 105
c
      RETURN
      END
