      Block Data model_init
      REAL*8 LCON(nknot),NCON(nknot),LSPL(3,nknot),NSPL(3,nknot)
      real*4 tdum(nknot,9)
      real*4 R(nknot),FMU(nknot),FLAM(nknot),QSHEAR(nknot),
     &       QKAPPA(nknot),XA2(nknot),XLAM(nknot),RHO(nknot),
     &       QRO(3,nknot),G(nknot),QG(3,nknot),FCON(nknot),
     &       FSPL(3,nknot),CCON(nknot),CSPL(3,nknot),ACON(nknot),
     &       ASPL(3,nknot)
      COMMON R,FMU,FLAM,QSHEAR,
     &       QKAPPA,XA2,XLAM,RHO,
     &       QRO,G,QG,FCON,
     &       FSPL,LCON,LSPL,NCON,
     &       NSPL,CCON,CSPL,ACON,
     &       ASPL
      common /head/ tdum
      DATA R/nknot*0.D0/,FMU/nknot*0.D0/,FLAM/nknot*0.D0/,
     &     QSHEAR/nknot*0.D0/,QKAPPA/nknot*0.D0/,XA2/nknot*0.D0/,
     &     XLAM/nknot*0.D0/,RHO/nknot*0.D0/,QRO/nknot3*0.D0/,
     &     G/nknot*0.D0/,QG/nknot3*0.D0/,FCON/nknot*0.D0/,
     &     FSPL/nknot3*0.D0/,LCON/nknot*0.D0/,LSPL/nknot3*0.D0/,
     &     NCON/nknot*0.D0/,NSPL/nknot3*0.D0/,CCON/nknot*0.D0/,
     &     CSPL/nknot3*0.D0/,ACON/nknot*0.D0/,ASPL/nknot3*0.D0/,
     &     tdum/nknot9*0.0/
      END Block Data model_init

