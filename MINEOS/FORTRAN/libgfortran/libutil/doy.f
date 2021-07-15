      INTEGER FUNCTION DOY(MONTH, DAY, YEAR)
C
C Function doy determines the day of the
C year,given the month,day and year.
C If month or day are illegal,the return
C value of the function is zero.
C
C
C Calls:
C   lpyr
C
C      Programmed by Madeleine Zirbes
C         September 15,1980
C
C MONTH - INPUT
      INTEGER MONTH
C DAY OF MONTH - INPUT
      INTEGER DAY
C YEAR - INPUT
      INTEGER YEAR
C FUNCTION
      INTEGER LPYR
      INTEGER INC
      INTEGER NDAYS(12)
      DATA NDAYS/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
C
      IF (.NOT.(MONTH .LT. 1 .OR. MONTH .GT. 12)) GOTO 2140
        DOY = (0)
        RETURN
2140  CONTINUE
      IF (.NOT.(DAY .LT. 1 .OR. DAY .GT. 31)) GOTO 2160
        DOY = (0)
        RETURN
2160  CONTINUE
      IF (.NOT.(LPYR(YEAR) .EQ. 1 .AND. MONTH .GT. 2)) GOTO 2180
        INC = 1
        GOTO 2190
2180  CONTINUE
        INC = 0
2190  CONTINUE
      DOY = (NDAYS(MONTH) + DAY + INC)
      RETURN
      END     
