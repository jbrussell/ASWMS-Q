C**** PROGRAM MINEOS ****
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
C  THIS PROGRAM USES TWO INPUT AND TWO OUTPUT FILES
C  THE INPUT FILES ARE 1) A MODEL FILE AND 2) A CONTROL FILE FOR SPECIFYING
C  WHICH MODES ARE TO BE CALCULATED AND TO WHAT ACCURACY ETC.
C  THE OUTPUT FILES ARE 1) A MODEL LISTING + A SUMMARY OF MODE PROPERTIES
C  AND 2) A FILE FOR THE EIGENFUNCTIONS (DISK OR TAPE)
C
C                     INPUT
C                    *******
C
C 1) MODEL FILE
C  CARD 1  :   TITLE (80 CHARS)                (20A4 FORMAT)
C  CARD 2  :   IFANIS,TREF,IFDECK              (UNFORMATTED)
C           IFANIS=1 FOR ANISOTROPIC MODEL, 0 FOR ISOTROPIC
C           TREF=REF PERIOD(SECS) OF MODEL FOR DISPERSION CORRECTION.
C           IF TREF IS .LE. 0. NO CORRECTION IS MADE.
C           IFDECK=1 FOR CARD DECK MODEL, 0 FOR POLYNOMIAL MODEL.
C             *** CARD DECK MODEL ***
C  CARD 3  :   N,NIC,NOC                       (UNFORMATTED)
C           N=NO OF LEVELS,NIC=INDEX OF SOLID SIDE OF ICB,NOC=INDEX OF
C           FLUID SIDE OF MCB.NOTE THAT N MUST BE .LE. nknot.
C  CARD 4-N:  R,RHO,VPV,VSV,QKAPPA,QSHEAR,VPH,VSH,ETA (F8.0,3F9.2,2F9.1,2F9.2,
C                                                           F9.5)
C           IF ISOTROPIC VP=VPV,VS=VSV AND VPH,VSH,ETA ARE UNSPECIFIED.
C           IF THE Q MODEL IS NOT SPECIFIED THEN NO DISP. CORRECTION IS MADE.
C           S.I. UNITS ARE USED,E.G. RADIUS IN METRES AND VELOCITIES IN M/S.
C             *** POLYNOMIAL MODEL ***
C  CARD 3  :   NREG,NIC,NOC,RX                 (UNFORMATTED)
C           NREG IS THE NUMBER OF REGIONS IN THE MODEL,NIC AND NOC AS BEFORE
C           AND RX IS THE NORMALISING RADIUS FOR THE POLYNOMIALS.
C           RX IS GIVEN IN KMS AND IS USUALLY 6371.
C  CARD 4  :   NLAY,R1,R2                      (UNFORMATTED)
C           NLAY IS THE NUMBER OF LEVELS TO BE USED IN THE REGION EXTENDING
C           FROM RADIUS R1 TO R2 (IN KMS).
C  CARD 5-9(ISO) OR CARD 5-12(ANI) : COEFS     (5F9.5)
C           5 SETS OF COEFFICIENTS ARE REQUIRED FOR EACH REGION OF AN ISOTROPIC
C           MODEL AND 8 SETS FOR AN ANISOTROPIC MODEL.EACH POLYNOMIAL CAN  BE
C           UP TO A QUARTIC IN R/RX ( 5 COEFFICIENTS) AND ARE ORDERED THUSLY:
C           RHO,VPV,VSV,QKAPPA,QSHEAR,VPH,VSH,ETA. THE COEFFS ARE GIVEN IN THE
C           USUAL MIXED SEISMOLOGICAL UNITS (RHO IN G/CC,VEL IN KM/S ETC.)
C           CONVERSION TO S.I. IS DONE BY THE PROGRAM.
C           CARDS 4-9(ISO) OR 4-12(ANI) ARE REPEATED FOR EACH REGION OF THE
C           MODEL
C
C 2) MODE FILE
C   CARD 1  :  EPS,EPS1,EPS2,WGRAV             (UNFORMATTED)
C           EPS CONTROLS THE ACCURACY OF THE RUNGE-KUTTA INTEGRATION SCHEME.
C           THE RELATIVE ACCURACY OF AN EIGEN FREQUENCY WILL BE 2-3 X EPS.
C           EPS1 CONTROLS THE PRECISION WITH WHICH A ROOT IS FOUND.
C           EPS2 IS THE MINIMUM RELATIVE SEPARATION OF TWO ROOTS WITH THE
C           SAME ANGULAR ORDER.IT IS SAFE TO SET EPS=EPS1=EPS2=1.D-7.
C           WGRAV IS THE FREQUENCY IN RAD/S ABOVE WHICH GRAVITATIONAL TERMS
C           ARE NEGLECTED-THIS GIVES ABOUT A FACTOR OF 3 INCREASE IN SPEED.
C   CARD 2  :  JCOM                            (UNFORMATTED)
C           JCOM=1 RADIAL MODES, 2 TOROIDAL MODES, 3 SPHEROIDAL MODES,
C           4 INNER CORE TOROIDAL MODES. IF JCOM IS .LE. 0 THE PROGRAM QUITS.
C   CARD 3  :  LMIN,LMAX,WMIN,WMAX             (UNFORMATTED)
C           LMIN - LMAX DEFINES THE RANGE OF ANGULAR ORDERS TO BE COMPUTED.
C           IF JCOM=1 THIS IS IGNORED. WMIN - WMAX DEFINES THE FREQUENCY RANGE
C           TO BE COMPUTED (IN mHz)
C           CARDS 2-3 ARE REPEATED UNTIL JCOM .LE. 0 IS ENCOUNTERED.
C
C                     OUTPUT
C                    ********
C     
c      05/22/91 - pp
c      common block /SHANKS/ reorganized to avoid performance degradation
c      common block /EIFX/   has the same length everywhere now (15*nknot r*8)
c
c      09/21/91 - pp
c      mode file output of ww,qq,gc (in modout) as real*8 instead of real*4
c
c      03/27/96 - pp
c      subroutines separated into individual files
c      dimensions changed: 1000 -> nknot; 215 -> nbranch; 430 -> nbranch2
c      include file parameter.h added
c
c      calls: MODEL, WTABLE
c     
      character name*300
      integer*4 iret
      SAVE
      IIN=2
c
      print *,' input model file : '
      read(*,101) name
101   format(a300)
      open(iin,file=name,status='old',
     + form='formatted',iostat=iret)
      if(iret.ne.0) then
         print *,' error opening model file',name
         stop
      end if
      iout=16
csun      iout=6
      print *,' output ascii file : '
      read(*,101) name
      open(iout,file=name,form='formatted',
     + iostat=iret)
      if(iret.ne.0) then
         print *,' error opening ascii output file',name
         stop
      end if
      IOEIG=3
c
      CALL MODEL(IIN,IOUT)
c
      close(iin)
      print *,' eigenfunction file name : '
      read(*,101) name
      open(ioeig,file=name,form='unformatted',iostat=iret)
      if(iret.ne.0) then
         print *,' error opening eigenfunction file',name
         stop
      end if
      print *,' mode file name : '
      read(*,101) name
      open(iin,file=name,status='old',
     + form='formatted',iostat=iret)
      if(iret.ne.0) then
         print *,' error opening mode file : ',name
         stop
      end if
c
      CALL WTABLE(IIN,IOUT,IOEIG)
c
      close(iin)
      close(iout)
      close(ioeig)
c
      print *,'All Finished'
      stop
      END
