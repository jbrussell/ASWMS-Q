      SUBROUTINE MATCH(N,J,KG,AF,AFR)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'parameter.h'
c
      COMMON/EIFX/AR(14,nknot),INORM(nknot),idum2(nknot)
      COMMON/AREM/A(6,3,nknot)
      DIMENSION AF(4,1),AFR(1)
c
      PRINT 999,'match 0',AF(1,1),AF(2,1),AFR(1),AFR(2)

      K=J+2
      IF(KG.EQ.1) GO TO 20
      C=(AF(1,1)*AFR(1)+AF(2,1)*AFR(2))/(AF(1,1)**2+AF(2,1)**2)
      AF1=AF(1,1)*C
      AF2=AF(2,1)*C
      PRINT 999,'match 1',AF1,AF2,AFR(1),AFR(2)
      IDIFF=INORM(J)-INORM(J+1)
      INORM(J+1)=INORM(J)
  999 FORMAT(a8,4G20.10)
      DO 10 I=K,N
      INORM(I)=INORM(I)+IDIFF
   10 CONTINUE
      DO 11 J=1,4
      DO 11 I=K,N
      AR(J,I)=C*A(J,1,I)
   11 CONTINUE
      RETURN
c
   20 A2=(AF(3,1)*AFR(1)-AF(1,1)*AFR(3))/(AF(1,2)*AF(3,1)-AF(1,1)
     +   *AF(3,2))
      A1=(AF(3,2)*AFR(1)-AF(1,2)*AFR(3))/(AF(1,1)*AF(3,2)-AF(3,1)
     +   *AF(1,2))
      AF1=A1*AF(1,1)+A2*AF(1,2)
      AF2=A1*AF(2,1)+A2*AF(2,2)
      AF3=A1*AF(3,1)+A2*AF(3,2)
      AF4=A1*AF(4,1)+A2*AF(4,2)
      PRINT 999,'match 2',AF1,AF2,AF3,AF4
      PRINT 999,'match 3',(AFR(I),I=1,4)
      IDIFF=INORM(J)-INORM(J+1)
      INORM(J+1)=INORM(J)
      DO 25 I=K,N
      INORM(I)=INORM(I)+IDIFF
      DO 25 J=1,6
   25 AR(J,I)=A1*A(J,1,I)+A2*A(J,2,I)
c
      RETURN
      END
