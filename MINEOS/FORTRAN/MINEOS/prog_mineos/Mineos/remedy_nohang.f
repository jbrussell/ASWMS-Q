      SUBROUTINE REMEDY(LS)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C    OBTAINS THE EIGENFUNCTION OF AN AWKWARD SPHEROIDAL MODE BY
C    INTEGRATING TO THE ICB OR THE MCB.
c
c    version to avoid hangups on tough modes -- see bottom of
c    routine -- jbg 1/00
c
c    calls: FPROP, FSBDRY, MATCH, ORTHO, SFBDRY, SPROP
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      SAVE
c
      include 'parameter.h'
c
      COMMON/BITS/PI,RN,VN,WN,W,WSQ,WRAY,QINV,CG,WGRAV,TREF,FCT,EPS,FL,
     +  FL1,FL2,FL3,SFL3,JCOM,NORD,L,KG,KOUNT,KNSW,IFANIS,IBACK
      COMMON/EIFX/AR(14,nknot),INORM(nknot),idum2(nknot)
      COMMON/AREM/A(6,3,nknot)
      COMMON/RINDX/NIC,NOC,NSL,NICP1,NOCP1,NSLP1,N
      DIMENSION AF(4,2),AS(6,3),AFR(4)
      DATA AF/8*0.D0/,AS/18*0.D0/,AFR/4*0.D0/,A/nknot18*0.D0/
c
      IEXP=0
      DO 10 K=1,2
      DO 10 J=1,4
   10 AF(J,K)=0.D0
      AF(1,1)=1.D0
      IF(KG.EQ.1) AF(2,2)=1.D0
      IF(NSL.EQ.N) GOTO 5
      DO 6 I=NSLP1,N
      DO 6 K=1,3
      DO 6 J=1,6
    6 A(J,K,I)=0.D0
c
      CALL FPROP(N,NSLP1,AF,IEXP)
c
    5 CALL FSBDRY(AF,AS,KG)
c
      DO 7 K=1,3
      DO 7 J=1,6
    7 A(J,K,NSL)=AS(J,K)
c
      IF(N.NE.NSL) CALL ORTHO(N,NSL,AS,KG)
c
      CALL SPROP(N,NSL,NOCP1,AS,IEXP)
c
      CALL SFBDRY(N,NOCP1,AS,AF,KG)
c
      IMTCH=NOC
      DO 11 I=1,4
   11 AFR(I)=AR(I,NOC)
      IF(LS.GT.NIC) GOTO 15
c
      CALL FPROP(NOC,NICP1,AF,IEXP)
c
      IMTCH=NIC
      DO 12 I=1,4
   12 AFR(I)=AR(I,NICP1)
c
c
c     check to make sure first two elements of afr are
c     not zero -- if they are, then match hangs.  Changes
c     LS to an error flag that is used back in detqn -- jbg 1/00 

   15 if(DABS(AFR(1)).lt.1.D-8 .and. DABS(AFR(2)).lt.1.D-8) then
         LS = -999
      else
         CALL MATCH(N,IMTCH,KG,AF,AFR)
      endif
c
      RETURN
      END
