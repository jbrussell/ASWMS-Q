      SUBROUTINE SFBDRY(JF,JL,AS,AF,KG)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** THE TANGENTIAL TRACTION SCALAR IS FORCED TO VANISH AT THE SOLID
C*** SIDE OF A S/F BOUNDARY(LEVEL JL).A(J,3,I) IS ELLIMINATED FOR
C*** I=JF...JL AND AF IS LOADED FROM A AT LEVEL JL.
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'parameter.h'
c
      COMMON/AREM/A(6,3,nknot)
      DIMENSION AS(6,1),AF(4,1)
c
      N1=MIN0(JF,JL)
      N2=MAX0(JF,JL)
      IF(KG.NE.0) GOTO 25
      I1=1
      I2=2
      IF(DABS(AS(4,2)).GT.DABS(AS(4,1))) GOTO 10
      I1=2
      I2=1
   10 RAT=-AS(4,I1)/AS(4,I2)
      DO 15 J=1,4
      DO 15 I=N1,N2
   15 A(J,1,I)=A(J,I1,I)+RAT*A(J,I2,I)
      AF(1,1)=A(1,1,JL)
      AF(2,1)=A(3,1,JL)
      RETURN
c
   25 AB53=DABS(AS(5,3))
      DO 30 K=1,2
      I1=K
      I2=3
      IF(AB53.GT.DABS(AS(5,K))) GOTO 35
      I1=3
      I2=K
   35 RAT=-AS(5,I1)/AS(5,I2)
      DO 40 I=N1,N2
      DO 40 J=1,6
   40 A(J,K,I)=A(J,I1,I)+RAT*A(J,I2,I)
      AF(1,K)=A(1,K,JL)
      AF(2,K)=A(3,K,JL)
      AF(3,K)=A(4,K,JL)
   30 AF(4,K)=A(6,K,JL)
c
      RETURN
      END
