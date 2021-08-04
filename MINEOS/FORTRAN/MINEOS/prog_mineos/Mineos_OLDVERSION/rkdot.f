       SUBROUTINE RKDOT(F,S,H,NVEC,NI)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C*** PERFORMS DOT PRODUCT WITH RKS COEFFICIENTS ***
c
c    calls no other routines
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      COMMON/SHANKS/B(46),C(10),DX,STEP(8),STEPF,IN,MAXO
      DIMENSION S(1),F(1),H(NVEC,1)
c
      GOTO (1,2,3,4,5,6,7,8,9,10),NI
    1 DO 21 J=1,NVEC
   21 F(J)=S(J)+B(1)*H(J,1)
      RETURN
c
    2 DO 22 J=1,NVEC
   22 F(J)=S(J)+B(2)*(H(J,1)+B(3)*H(J,2))
      RETURN
c
    3 DO 23 J=1,NVEC
   23 F(J)=S(J)+B(4)*(H(J,1)+B(5)*H(J,2)+B(6)*H(J,3))
      RETURN
c
    4 DO 24 J=1,NVEC
   24 F(J)=S(J)+B(7)*(H(J,1)+B(8)*H(J,2)+B(9)*H(J,3)+B(10)*H(J,4))
      RETURN
c
    5 DO 25 J=1,NVEC
   25 F(J)=S(J)+B(11)*(H(J,1)+B(12)*H(J,2)+B(13)*H(J,3)+B(14)*H(J,4)+
     +B(15)*H(J,5))
      RETURN
c
    6 DO 26 J=1,NVEC
   26 F(J)=S(J)+B(16)*(H(J,1)+B(17)*H(J,2)+B(18)*H(J,3)+B(19)*H(J,4)+
     +B(20)*H(J,5)+B(21)*H(J,6))
      RETURN
c
    7 DO 27 J=1,NVEC
   27 F(J)=S(J)+B(22)*(H(J,1)+B(23)*H(J,3)+B(24)*H(J,4)+B(25)*H(J,5)+
     +B(26)*H(J,6)+B(27)*H(J,7))
      RETURN
c
    8 DO 28 J=1,NVEC
   28 F(J)=S(J)+B(28)*(H(J,1)+B(29)*H(J,3)+B(30)*H(J,4)+B(31)*H(J,5)+
     +B(32)*H(J,6)+B(33)*H(J,7)+B(34)*H(J,8))
      RETURN
c
    9 DO 29 J=1,NVEC
   29 F(J)=S(J)+B(35)*(H(J,1)+B(36)*H(J,3)+B(37)*H(J,4)+B(38)*H(J,5)+
     +B(39)*H(J,6)+B(40)*H(J,7)+B(41)*H(J,8)+B(42)*H(J,9))
      RETURN
c
   10 DO 30 J=1,NVEC
   30 F(J)=S(J)+B(43)*(H(J,1)+H(J,10)+B(44)*(H(J,4)+H(J,6))+
     +B(45)*H(J,5)+B(46)*(H(J,7)+H(J,9)))
c
      RETURN
      END
