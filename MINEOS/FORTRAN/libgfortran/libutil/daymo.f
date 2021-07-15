      INTEGER FUNCTION DAYMO(DOFY, MONTH, DAY, YEAR)
C
C Function daymo determines the month and day
C of the month,given the year and day of year.
C It returns 1 if it was successful,0 otherwise.
C If dofy is not within legal limits,month and
C day will be returned as zero.
C
C
C Calls:
C   lpyr
C
C      Programmed by Madeleine Zirbes
C         September 15,1980
C
C DAY OF YEAR - INPUT
      INTEGER DOFY
C MONTH - OUTPUT
      INTEGER MONTH
C DAY OF MONTH - OUTPUT
      INTEGER DAY
C YEAR - INPUT
      INTEGER YEAR
C
C DAY OF YEAR
      INTEGER IDAY
C FUNCTION
      INTEGER LPYR
C NUMBER OF DAYS IN MONTH
      INTEGER MDAYS(12)
      DATA MDAYS/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
C
      IDAY = DOFY
      IF (.NOT.(IDAY .LT. 1)) GOTO 2060
        MONTH = 0
        DAY = 0
        DAYMO = (0)
        RETURN
C
2060  CONTINUE
      IF (.NOT.(LPYR(YEAR) .EQ. 1)) GOTO 2080
        MDAYS(2) = 29
        GOTO 2090
2080  CONTINUE
        MDAYS(2) = 28
C
2090  CONTINUE
      DO 2100 MONTH = 1, 12
        DAY = IDAY
        IDAY = IDAY - MDAYS(MONTH)
        IF (.NOT.(IDAY .LE. 0)) GOTO 2120
          DAYMO = (1)
          RETURN
2120    CONTINUE
C
2100  CONTINUE
      MONTH = 0
      DAY = 0
      DAYMO = (0)
      RETURN
C
      END
