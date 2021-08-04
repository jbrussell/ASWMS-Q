      SUBROUTINE FSBM(ASS,KG,IBACK)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** CONVERT MINOR VECTOR AT A FLUID/SOLID BOUNDARY ***
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      DIMENSION ASS(14),AS(14)
c
      DO 10 J=1,14
      AS(J)=ASS(J)
   10 ASS(J)=0.D0
      IF(IBACK.EQ.1) GOTO 30
      IF(KG.NE.0) GOTO 20
      ASS(1)=AS(1)
      ASS(4)=-AS(2)
      RETURN
c
   20 ASS(6)=AS(1)
      ASS(14)=AS(2)
      ASS(1)=AS(3)
      ASS(9)=AS(4)
      ASS(4)=-AS(5)
      RETURN
c
   30 IF(KG.NE.0) GOTO 40
      ASS(1)=-AS(1)
      RETURN
c
   40 ASS(1)=-AS(1)
      ASS(3)=-AS(2)
      ASS(5)=-AS(3)
      ASS(12)=AS(4)
c
      RETURN
      END
