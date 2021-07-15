      INTEGER FUNCTION LPYR(YEAR)
C
C Function lpyr determines if year
C is a leap year.
C
C This function uses the intrinsic
C function mod. If your machine
C does not supply this function,
C make one -
C mod(i,j) = iabs(i - (i/j)*j)
C
c returns 1 for leap year, 0 otherwise -- jbg 6/97
C
C Calls:
C   mod - intrinsic funtion
C
C      Programmed by Madeleine Zirbes
C         September 15,1980
C
C YEAR - INPUT
      INTEGER YEAR
      IF (.NOT.(MOD(YEAR, 400) .EQ. 0)) GOTO 6100
        LPYR = (1)
        RETURN
6100  CONTINUE
      IF (.NOT.(MOD(YEAR, 4) .NE. 0)) GOTO 6120
        LPYR = (0)
        RETURN
6120  CONTINUE
      IF (.NOT.(MOD(YEAR, 100) .EQ. 0)) GOTO 6140
        LPYR = (0)
        RETURN
6140  CONTINUE
      LPYR = (1)
      RETURN
      END
