C****************************************************************
C****************************************************************
C**REACTION PATH (AB INITIO POTENTIAL....NON-ECKART)
C****************************************************************
C****************************************************************
      SUBROUTINE REACT(NATOM,NMODE,XM,X0,OMEGA,XL,XLP,XX,XXP,RR,XR,XP,
     1TEMP,BB,BBS,B,BS,BSS,IREACT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ECKART
C****************************************************************
C****************************************************************
C**ALLOW ONE POINT FOR EQUILIBRIUM (1)
C**ALLOWS FOR ONE POINT EVERY HALF-DEGREE (721)
C**ALLOW EXTRA POINT FOR EQUILIBRIUM 2ND TIME (722)
C**
C**WORK IN DEGREES
C**
C****************************************************************
C****************************************************************
      DIMENSION XM(NATOM),X0(NATOM,3),XL(NATOM,NMODE,3,722)
      DIMENSION XLP(NATOM,NMODE,3,722)
      DIMENSION XR(NATOM,3),OMEGA(NMODE,722),RR(NATOM,NATOM)
      DIMENSION XX(NATOM,3,722),XXP(NATOM,3,722)
      DIMENSION XP(3*NATOM,3*NATOM,722),TEMP(3*NATOM,3*NATOM)
      DIMENSION BB(NMODE,NMODE,3,722),BBS(NMODE,NMODE,722)
      DIMENSION B(NMODE,3,3,722),BS(NMODE,3,722),BSS(NMODE,722)
      DIMENSION XNORM(722)

      DIMENSION SUP1(3,3),SUP2(3,3),SUP3(10),SUP4(4,4),SUP5(4,4)
      DIMENSION SUP6(20)

      COMMON/SCOORD/DSTAU(722),SPLUS(722),SMINUS(722)
      COMMON/ECKCON/ECKART
      COMMON/MOMI/XK(3,3),XMU(3,3)
      COMMON/REACTN/IDUMMY,MMTAU,INIT,INCTAU

C**SPECIFIC TO METHANOL
      DIMENSION VVV(362)
      COMMON/DAVE/VE,V(722),WE(12),DOMEGA(11,722),DERIV(11,722)
C**SPECIFIC TO METHANOL

C**SPECIFIC TO ABINITIO
      DIMENSION DERIVS(NATOM,3,722),PAROT(3,3,722)
C**SPECIFIC TO ABINITIO

      COMMON/RETURN/IRET
      COMMON/CHECK/MCHECK
      COMMON/PLMNMX/MXCD(3)
      COMMON/PATH/ISCFCI
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
220   FORMAT(/,1X,'MODE',I4,'  OMEGA = ',F20.12,/)
245   FORMAT(//,1X,'DISPLACEMENTS OF NORMAL MODES, L(alpha i,k)')
250   FORMAT(/,1X,'MODE(k) = ',I4)
255   FORMAT(1X,'ATOM(i) ',I2,':  x =',F12.6,'  y =',F12.6,'  z =',
     1F12.6)
260   FORMAT(//,1X,'ROTATED GEOMETRY',/)
265   FORMAT(1X,'ATOM ',I2,':  X0 =',F12.7,'  Y0 =',F12.7,'  Z0 =',
     1F12.7)
C*******************************************************************
C*******************************************************************
      ECKART=.FALSE.
      DO I=1,3
        MXCD(I)=0
      END DO

C**OPEN ALL DISC FILES BEFORE START OF LOOP
C**
C**EQUILIBRIUM ENERGY
      open(3,file='ve.dat')
      read(3,*) VE
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'EQUILIBRIUM POTENTIAL = ',VE
      END IF
C**EQUILIBRIUM OMEGAs
      open(4,file='we.dat')
      read(4,*) (WE(I), I=1,12)
      DO I=1,12
        WE(I)=WE(I)/WAVENM
      END DO
C**REACTION PATH ENERGIES
      open(5,file='v.dat')
      read(5,*) VVV
      DO ITAU=1,362
C**USE SYMMETRY AT 'TRANS' (300 DEGREES)
        IF(ITAU.GT.242)GO TO 5
C**USE SYMMETRY AT 'CIS' (240 DEGREES)
        IF(ITAU.GT.122)GO TO 6
        V(ITAU)=VVV(ITAU)
        GO TO 4
6       CONTINUE
        JNCR=ITAU-122
        IPLUS=122+JNCR
        IMINUS=122-JNCR
        V(IPLUS)=V(IMINUS)
        GO TO 4
5       CONTINUE
        JNCR=ITAU-242
        IPLUS=242+JNCR
        IMINUS=242-JNCR
        V(IPLUS)=V(IMINUS)
4       CONTINUE
      END DO
C**USE SYMMETRY AT 'CIS' (360 DEGREES)
      DO ITAU=1,360
        V(362+ITAU)=V(362-ITAU)
      END DO
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)'REACTION PATH POTENTIAL'
        DO ITAU=2,722
          WRITE(IOUT,*)(V(ITAU)-VE)*WAVENM
        END DO
      END IF
C**GEOMETRIES
      open(16,file='geom.dat')
C**FORCE CONSTANTS (SECOND DERIVATIVES)
      open(17,file='fcm.dat')
C**FIRST DERIVATIVES
      open(18,file='deriv.dat')

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'****************EQUILIBRIUM COORDINATE PARAMETERS'
        WRITE(IOUT,*)
C**ORIGINAL EQUILIBRIUM FROM DATA (FORT.1) FILE
        DO I=1,NATOM
          WRITE(IOUT,*)(X0(I,K),K=1,3)
        END DO
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**READ SINGLE 'EQU' GEOMETRY INTO XR(NATOM,3) FROM CADPAC
      read(16,*) ((XR(I,J), J=1,3),I=1,NATOM)

C**RETURN MARKER FOR CHECKM
      IRET=0
C**MOVE TO C. of M. AND ROTATE TO PRINCIPAL AXES
      CALL CHECKM(XM,XR,XX,RR,NATOM)
C**  XR CONTAINS NEW (TRANSLATED, ROTATED) GEOMETRY
      DO I=1,NATOM
        DO K=1,3
          XX(I,K,1)=XR(I,K)
        END DO
      END DO
C**RETURN MARKER FOR CHECKM
      IRET=-1
C*******************************************************************
C*******************************************************************
C**IN TERMS OF CURVELINEAR COORDINATE (TAU), EQUILIBRIUM IS TAU=0
      TAU=0
C**GET ITAU=1 GEOMETRY
      TAUE=0
C**DEFINE FIRST POINT 'TRANS'
      INCR=180
C*********************
C**GET INCREMENT TO NEXT POINT
      TAUINC=0.0D0
C**INIT WILL POINT AT TRANS (180 DEGREES)
      INIT=0
C**JNIT WILL POINT AT CIS (240 DEGREES)
      JNIT=0
C**KNIT WILL POINT AT TRANS (300 DEGREES)
      KNIT=0
C**LNIT WILL POINT AT CIS (360 DEGREES)
      LNIT=0

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'INCR,TAUINC',INCR,TAUINC
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************START REACTION COORDINATE ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      DO ITAU=1,361
C**ONLY GO AS FAR AS 'CIS' (360 DEGREES)
        IF(LNIT.NE.0)GO TO 9999

        IF(ITAU.EQ.1)THEN
          TAU=TAU+TAUINC/RAD
        ELSE
          TAU=TAU+0.5D0/RAD
        END IF

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)
          WRITE(IOUT,*)'****************NEXT POINT********************'
          WRITE(IOUT,*)
          WRITE(IOUT,*)'ITAU = ',ITAU+1,'  TAU = ',(TAU+TAUE)*RAD
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**USE SYMMETRY AT 'TRANS' (300 DEGREES)
        IF(ITAU+1.GT.KNIT.AND.KNIT.NE.0)GO TO 5555
C**USE SYMMETRY AT 'CIS' (240 DEGREES)
        IF(ITAU+1.GT.JNIT.AND.JNIT.NE.0)GO TO 6666

C
C

C**DATA FOR EACH POINT GOES INTO POSITION ITAU+1.....STARTS AT 2
C**READ SINGLE GEOMETRY INTO XX(NATOM,3,ITAU+1) FROM CADPAC
        read(16,*) ((XR(I,J), J=1,3),I=1,NATOM)

C
C

C**SET TEMPORARY EQUILIBRIUM FROM PREVIOUS POINT
        DO I=1,NATOM
          DO K=1,3
            X0(I,K)=XX(I,K,ITAU)
          END DO
        END DO
C**POINT TO 180 DEGREES
        IF(INIT.EQ.0)THEN
          IF(INCR+(ITAU-1)/2.EQ.180)INIT=ITAU+1
          INITCD=1
          JX=1
          DO J=1,NATOM
            DO K=1,3
              IF(DABS(XX(J,K,INIT-1)).LT.DABS(XX(JX,INITCD,INIT-1)))
     1        THEN
                INITCD=K
                JX=J
              END IF
            END DO
          END DO
          IF(INIT.NE.0)THEN
            IF(JPRINT.GT.0)WRITE(IOUT,*)'180 DEGREES IS AT POINT ',INIT
            IF(ITAU.NE.1)INITCD=MXCD(1)
            IF(JPRINT.GT.0)WRITE(IOUT,*)'ZERO COORDINATE IS ',INITCD
          END IF
        END IF
C**POINT TO 240 DEGREES
        IF(JNIT.EQ.0)THEN
          IF(INCR+(ITAU-1)/2.EQ.240)JNIT=ITAU+1
          JNITCD=1
          JX=1
          DO J=1,NATOM
            DO K=1,3
              IF(DABS(XX(J,K,JNIT-1)).LT.DABS(XX(JX,JNITCD,JNIT-1)))
     1        THEN
                JNITCD=K
                JX=J
              END IF
            END DO
          END DO
          IF(JNIT.NE.0)THEN
            IF(JPRINT.GT.0)WRITE(IOUT,*)'240 DEGREES IS AT POINT ',JNIT
            IF(ITAU.NE.1)JNITCD=MXCD(1)
            IF(JPRINT.GT.0)WRITE(IOUT,*)'ZERO COORDINATE IS ',JNITCD
          END IF
        END IF
C**NOW MOVE TO CENTRE OF MASS AND ROTATE TO PRINCIPAL AXES
        CALL GETSAY(NATOM,XX(1,1,ITAU),XM,XX(1,1,ITAU+1),XR,RR)
C**SAVE ROTATION MATRIX
        DO I=1,3
          DO J=1,3
            PAROT(I,J,ITAU+1)=XK(I,J)
            IF(ITAU.EQ.1)PAROT(I,J,ITAU)=XK(I,J)
          END DO
        END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)
          WRITE(IOUT,*)'OLD ROTATION MATRIX'
          DO I=1,3
            WRITE(IOUT,*)(PAROT(I,J,ITAU+1),J=1,3)
          END DO
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**READ FIRST DERIVATIVE MATRIX HERE INTO DERIVS
C**Atom 1: x,y,z ... Atom 2: x,y,z , etc.
        READ(18,*)((DERIVS(I,J,ITAU+1),J=1,3),I=1,NATOM)
C**GET PHASE RIGHT
        ICH=0
        IX=1
        DO I=2,NATOM
          IF(DABS(XX(I,1,ITAU+1)).GT.DABS(XX(IX,1,ITAU+1)))IX=I
        END DO
        IF(XX(IX,1,ITAU+1).GE.0.AND.XX(IX,1,ITAU).LT.0)ICH=1
        IF(XX(IX,1,ITAU+1).LE.0.AND.XX(IX,1,ITAU).GT.0)ICH=1
        IF(ICH.NE.0)THEN
          DO I=1,NATOM
            XX(I,1,ITAU+1)=-XX(I,1,ITAU+1)
          END DO
          DO J=1,3
            PAROT(J,1,ITAU+1)=-PAROT(J,1,ITAU+1)
          END DO
        END IF
        ICH=0
        IY=1
        DO I=2,NATOM
          IF(DABS(XX(I,2,ITAU+1)).GT.DABS(XX(IY,2,ITAU+1)))IY=I
        END DO
        IF(XX(IY,2,ITAU+1).GE.0.AND.XX(IY,2,ITAU).LT.0)ICH=1
        IF(XX(IY,2,ITAU+1).LE.0.AND.XX(IY,2,ITAU).GT.0)ICH=1
        IF(ICH.NE.0)THEN
          DO I=1,NATOM
            XX(I,2,ITAU+1)=-XX(I,2,ITAU+1)
          END DO
          DO J=1,3
            PAROT(J,2,ITAU+1)=-PAROT(J,2,ITAU+1)
          END DO
        END IF
        ICH=0
        IZ=1
        DO I=2,NATOM
          IF(DABS(XX(I,3,ITAU+1)).GT.DABS(XX(IZ,3,ITAU+1)))IZ=I
        END DO
        IF(XX(IZ,3,ITAU+1).GE.0.AND.XX(IZ,3,ITAU).LT.0)ICH=1
        IF(XX(IZ,3,ITAU+1).LE.0.AND.XX(IZ,3,ITAU).GT.0)ICH=1
        IF(ICH.NE.0)THEN
          DO I=1,NATOM
            XX(I,3,ITAU+1)=-XX(I,3,ITAU+1)
          END DO
          DO J=1,3
            PAROT(J,3,ITAU+1)=-PAROT(J,3,ITAU+1)
          END DO
        END IF
C**USE ROTATION OF ORIGINAL GEOMETRY TO GET CORRECT 1ST DERIVATIVES
        DO I=1,NATOM
          DO J=1,3
            RR(I,J)=0
            DO K=1,3
              RR(I,J)=RR(I,J)+PAROT(K,J,ITAU+1)*DERIVS(I,K,ITAU+1)
            END DO
          END DO
        END DO
        DO I=1,NATOM
          DO K=1,3
            DERIVS(I,K,ITAU+1)=RR(I,K)
          END DO
        END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)
          WRITE(IOUT,*)'NEW ROTATION MATRIX'
          DO I=1,3
            WRITE(IOUT,*)(PAROT(I,J,ITAU+1),J=1,3)
          END DO
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        Z=0
        DO I=1,NATOM
          DO K=1,3
            Z=Z+XX(I,K,ITAU)*XX(I,K,ITAU+1)
          END DO
        END DO
        Z=SQRT(Z)
        IF(Z.LT.1.D-03)THEN
          IF(JPRINT.GT.0)WRITE(IOUT,*)'FAILURE!!!  Z = ',Z
        ELSE
          IF(JPRINT.GT.0)WRITE(IOUT,*)'SUCCESS!!!  Z = ',Z
        END IF
        GO TO 4444
6666    CONTINUE
C**POINT TO 300 DEGREES
        IF(KNIT.EQ.0)THEN
          IF(INCR+(ITAU-1)/2.EQ.300)KNIT=ITAU+1
          KNITCD=1
          JX=1
          DO J=1,NATOM
            DO K=1,3
              IF(DABS(XX(J,K,KNIT-1)).LT.DABS(XX(JX,KNITCD,KNIT-1)))
     1        THEN
                KNITCD=K
                JX=J
              END IF
            END DO
          END DO
          IF(KNIT.NE.0)THEN
            IF(JPRINT.GT.0)WRITE(IOUT,*)'300 DEGREES IS AT POINT ',KNIT
            IF(ITAU.NE.1)KNITCD=MXCD(1)
            IF(JPRINT.GT.0)WRITE(IOUT,*)'ZERO COORDINATE IS ',KNITCD
          END IF
        END IF
C**************************************
C**************************************
C**SET GEOMETRIES CORRECTLY ABOUT '240'
C**************************************
C**************************************
        JNCR=ITAU-121
        IPLUS=122+JNCR
        IMINUS=122-JNCR
        DO J=1,NATOM
          DO K=1,3
            XX(J,K,IPLUS)=XX(J,K,IMINUS)
            DERIVS(J,K,IPLUS)=DERIVS(J,K,IMINUS)
            IF(K.EQ.JNITCD)THEN
              XX(J,K,IPLUS)=-XX(J,K,IMINUS)
              DERIVS(J,K,IPLUS)=-DERIVS(J,K,IMINUS)
            END IF
          END DO
        END DO
        DO K=1,3
          XX(3,K,IPLUS)=XX(4,K,IMINUS)
          XX(4,K,IPLUS)=XX(3,K,IMINUS)
          DERIVS(3,K,IPLUS)=DERIVS(4,K,IMINUS)
          DERIVS(4,K,IPLUS)=DERIVS(3,K,IMINUS)
          IF(K.EQ.JNITCD)THEN
            XX(3,K,IPLUS)=-XX(4,K,IMINUS)
            XX(4,K,IPLUS)=-XX(3,K,IMINUS)
            DERIVS(3,K,IPLUS)=-DERIVS(4,K,IMINUS)
            DERIVS(4,K,IPLUS)=-DERIVS(3,K,IMINUS)
          END IF
        END DO
        GO TO 4444
5555    CONTINUE
C**POINT TO 360 DEGREES
        IF(LNIT.EQ.0)THEN
          IF(INCR+(ITAU-1)/2.EQ.360)LNIT=ITAU+1
          LNITCD=1
          JX=1
          DO J=1,NATOM
            DO K=1,3
              IF(DABS(XX(J,K,LNIT-1)).LT.DABS(XX(JX,LNITCD,LNIT-1)))
     1        THEN
                LNITCD=K
                JX=J
              END IF
            END DO
          END DO
          IF(LNIT.NE.0)THEN
            IF(JPRINT.GT.0)WRITE(IOUT,*)'360 DEGREES IS AT POINT ',LNIT
            IF(ITAU.NE.1)LNITCD=MXCD(1)
            IF(JPRINT.GT.0)WRITE(IOUT,*)'ZERO COORDINATE IS ',LNITCD
          END IF
        END IF
C**************************************
C**************************************
C**SET GEOMETRIES CORRECTLY ABOUT '300'
C**************************************
C**************************************
        JNCR=ITAU-241
        IPLUS=242+JNCR
        IMINUS=242-JNCR
        DO J=1,NATOM
          DO K=1,3
            XX(J,K,IPLUS)=XX(J,K,IMINUS)
            DERIVS(J,K,IPLUS)=DERIVS(J,K,IMINUS)
            IF(K.EQ.KNITCD)THEN
              XX(J,K,IPLUS)=-XX(J,K,IMINUS)
              DERIVS(J,K,IPLUS)=-DERIVS(J,K,IMINUS)
            END IF
          END DO
        END DO
        DO K=1,3
          XX(3,K,IPLUS)=XX(5,K,IMINUS)
          XX(5,K,IPLUS)=XX(3,K,IMINUS)
          DERIVS(3,K,IPLUS)=DERIVS(5,K,IMINUS)
          DERIVS(5,K,IPLUS)=DERIVS(3,K,IMINUS)
          IF(K.EQ.KNITCD)THEN
            XX(3,K,IPLUS)=-XX(5,K,IMINUS)
            XX(5,K,IPLUS)=-XX(3,K,IMINUS)
            DERIVS(3,K,IPLUS)=-DERIVS(5,K,IMINUS)
            DERIVS(5,K,IPLUS)=-DERIVS(3,K,IMINUS)
          END IF
        END DO
4444    CONTINUE

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C**INTERNAL DISPLACEMENT COORDINATES
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,260)
          DO I=1,NATOM
            WRITE(IOUT,265)I,(XX(I,J,ITAU+1),J=1,3)
          END DO
          WRITE(IOUT,*)
          WRITE(IOUT,*)
          WRITE(IOUT,*)'*********************END ROTATION OF GEOMETRY'
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

9999    CONTINUE
      END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************END REACTION COORDINATE ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**************************************
C**************************************
C**SET GEOMETRIES CORRECTLY ABOUT '360'
C**************************************
C**************************************
      IPLUS=LNIT
      IMINUS=LNIT
      DO I=1,360
        IPLUS=IPLUS+1
        IF(IPLUS.GT.722)IPLUS=3
        IMINUS=IMINUS-1
        IF(IMINUS.LT.2)IMINUS=721
        DO J=1,NATOM
          DO K=1,3
            XX(J,K,IPLUS)=XX(J,K,IMINUS)
            DERIVS(J,K,IPLUS)=DERIVS(J,K,IMINUS)
            IF(K.EQ.LNITCD)THEN
              XX(J,K,IPLUS)=-XX(J,K,IMINUS)
              DERIVS(J,K,IPLUS)=-DERIVS(J,K,IMINUS)
            END IF
          END DO
        END DO
        DO K=1,3
          XX(4,K,IPLUS)=XX(5,K,IMINUS)
          XX(5,K,IPLUS)=XX(4,K,IMINUS)
          DERIVS(4,K,IPLUS)=DERIVS(5,K,IMINUS)
          DERIVS(5,K,IPLUS)=DERIVS(4,K,IMINUS)
          IF(K.EQ.LNITCD)THEN
            XX(4,K,IPLUS)=-XX(5,K,IMINUS)
            XX(5,K,IPLUS)=-XX(4,K,IMINUS)
            DERIVS(4,K,IPLUS)=-DERIVS(5,K,IMINUS)
            DERIVS(5,K,IPLUS)=-DERIVS(4,K,IMINUS)
          END IF
        END DO
      END DO
C*******************************************
C**SET TOTAL MASS
      XMTOT=0
      DO I=1,NATOM
C**SUM MASSES
        XMTOT=XMTOT+XM(I)
        DO K=1,3
C**RESET TRUE EQUILIBRIUM
          X0(I,K)=XX(I,K,1)
        END DO
      END DO
C**WHEN HAPPY FORM MASS-WEIGHTED COORDINATES
      DO ITAU=2,722
        DO K=1,3
          DO I=1,NATOM
            XX(I,K,ITAU)=XX(I,K,ITAU)*SQRT(XM(I))
          END DO
        END DO
      END DO
C*******************************************************************
C*******************************************************************

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'***********START COORDINATE DERIVATIVE ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      DO ITAU=2,722
C**CALCULATE GRADIENTS OF COORDINATES AT TAU = ITAU
C**************************************************
C**ITAU-1 AND ITAU+1 DIFFER BY ONE DEGREE
        DO I=1,NATOM
          DO K=1,3
            XXP(I,K,ITAU)=(XX(I,K,ITAU+1)-XX(I,K,ITAU-1))
            IF(ITAU.EQ.2)
     1      XXP(I,K,ITAU)=(XX(I,K,ITAU+1)-XX(I,K,ITAU))*2
            IF(ITAU.EQ.722)
     1      XXP(I,K,ITAU)=(XX(I,K,ITAU)-XX(I,K,ITAU-1))*2
          END DO
        END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'****************NEXT POINT********************'
          WRITE(IOUT,*)
          WRITE(IOUT,*)'ITAU= ',ITAU
          WRITE(IOUT,*)
          WRITE(IOUT,*)'GEOMETRY AT POINT ITAU = ',ITAU
          WRITE(IOUT,*)
          DO I=1,NATOM
            WRITE(IOUT,265)I,(XX(I,J,ITAU)/SQRT(XM(I)),J=1,3)
          END DO
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C********************************************
C**FORM REACTION COORDINATE (S) AT ITAU+1, ITAU-1
        X=0
        DO I=1,NATOM
          IF(ITAU.NE.722.AND.ITAU.NE.INIT-1.AND.ITAU.NE.LNIT-1)THEN
            DO K=1,3
              X=X+(XX(I,K,ITAU+1)-XX(I,K,ITAU))**2
            END DO
          ELSE
            DO K=1,3
              X=X+(XX(I,K,ITAU)-XX(I,K,ITAU-1))**2
            END DO
          END IF
        END DO
        SPLUS(ITAU)=SQRT(X)
        X=0
        DO I=1,NATOM
          IF(ITAU.NE.2.AND.ITAU.NE.INIT+1.AND.ITAU.NE.LNIT+1)THEN
            DO K=1,3
              X=X+(XX(I,K,ITAU)-XX(I,K,ITAU-1))**2
            END DO
          ELSE
            DO K=1,3
              X=X+(XX(I,K,ITAU+1)-XX(I,K,ITAU))**2
            END DO
          END IF
        END DO
        SMINUS(ITAU)=SQRT(X)
C********************************************
C**GET ARC-LENGTH/RADIAN JACOBIAN JUST FOR KICKS!
        DSTAU(ITAU)=(SPLUS(ITAU)+SMINUS(ITAU))*RAD

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)WRITE(IOUT,*)'REACTION PATH = ',DSTAU(ITAU)
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**NORMALISE
        Z=0
        DO I=1,NATOM
          DO K=1,3
            Z=Z+XXP(I,K,ITAU)**2
          END DO
        END DO
        Z=SQRT(Z)
        IF(Z.LT.1.D-03)THEN
          IF(JPRINT.GT.0)WRITE(IOUT,*)'BEWARE!!!  Z = ',Z
          DO I=1,NATOM
            DO K=1,3
              XXP(I,K,ITAU)=XXP(I,K,ITAU-1)
            END DO
          END DO
        END IF
        IF(Z.GE.1.D-03)THEN
          DO I=1,NATOM
            DO K=1,3
              XXP(I,K,ITAU)=XXP(I,K,ITAU)/Z
            END DO
          END DO
        END IF
C**SAVE NORMALISATION CONSTANT
        XNORM(ITAU)=Z

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'XNORM = ',Z
          WRITE(IOUT,*)
          WRITE(IOUT,*)'DERIVATIVES AT POINT ITAU = ',ITAU
          WRITE(IOUT,*)
          DO I=1,NATOM
            WRITE(IOUT,265)I,(XXP(I,K,ITAU),K=1,3)
          END DO
        END IF
      END DO
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************END COORDINATE DERIVATIVE ALGORITHM'
        WRITE(IOUT,*)

C**PRINT EQUILIBRIUM VALUES NORMAL MODES FROM FORT.1 FILE
        WRITE(IOUT,*)
        WRITE(IOUT,*)
        WRITE(IOUT,*)'EQUILIBRIUM 3N-7 NORMAL MODES'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      J0=0
      DO J=1,NMODE-1
        J0=J0+1
C**MODE FOR REACTION PATH WILL BE MISSING, SO SKIP
        IF(J0.EQ.IREACT)J0=J0+1

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)'OMEGA ',J,': ',OMEGA(J0,1)*WAVENM
          DO I=1,NATOM
            WRITE(IOUT,*)(XL(I,J0,K,1),K=1,3)
          END DO
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      END DO
C*******************************************************************
C*******************************************************************

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************START PROJECTION ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      DO ITAU=1,362
C**SET UP MOMENT OF INERTIA MATRIX
C*********************************

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'****************NEXT POINT********************'
          WRITE(IOUT,*)
          WRITE(IOUT,*)'ITAU= ',ITAU
          WRITE(IOUT,*)
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        IF(ITAU.EQ.1)GO TO 7777

C**USE SYMMETRY AT 'TRANS' (300 DEGREES)
        IF(ITAU.GT.KNIT)GO TO 555
C**USE SYMMETRY AT 'CIS' (240 DEGREES)
        IF(ITAU.GT.JNIT)GO TO 666

        DO IX=1,3
          DO IY=1,IX
            XK(IY,IX)=0
          END DO
        END DO
        DO IX=1,3
          DO I=1,NATOM
            X=0
            DO J=1,3
              IF(J.NE.IX)X=X+XX(I,J,ITAU)*XX(I,J,ITAU)
            END DO
            XK(IX,IX)=XK(IX,IX)+X
          END DO
        END DO
        DO IX=1,3
          DO IY=1,IX-1
            DO I=1,NATOM
              XK(IY,IX)=XK(IY,IX)-XX(I,IY,ITAU)*XX(I,IX,ITAU)
            END DO
          END DO
        END DO
        DO IX=1,3
          DO IY=1,IX
            XK(IX,IY)=XK(IY,IX)
          END DO
        END DO
C**CALCULATE INVERSE MOMENT OF INERTIA
C*************************************
        DO I=1,3
          DO J=1,3
            SUP1(I,J)=XK(I,J)
          END DO
        END DO
        IFAIL=1
        CALL MATINV(SUP1,3,3,SUP2,SUP3,IFAIL)
        IF(IFAIL.NE.0)THEN
          IF(JPRINT.GT.0)WRITE(IOUT,*)'IFAIL = ',IFAIL,' IN XMU MATINV'
          STOP 'ERROR IN MATINV'
        END IF

C**TEMPORARY ??
        DO I=1,3
          DO J=1,3
            IF(SUP2(I,J).LT.0)SUP2(I,J)=0
          END DO
        END DO
C**TEMPORARY ??

        DO I=1,3
          DO J=1,3
            XMU(I,J)=SQRT(SUP2(I,J))
          END DO
        END DO
C**FORM PROJECTION MATRIX
C************************
        DO I=1,NATOM
          DO K=1,3
            XL(I,4,K,ITAU)=XXP(I,K,ITAU)
            XL(I,K,1,ITAU)=(XMU(K,2)*XX(I,3,ITAU)-
     1                      XMU(K,3)*XX(I,2,ITAU))
            XL(I,K,2,ITAU)=(XMU(K,3)*XX(I,1,ITAU)-
     1                      XMU(K,1)*XX(I,3,ITAU))
            XL(I,K,3,ITAU)=(XMU(K,1)*XX(I,2,ITAU)-
     1                      XMU(K,2)*XX(I,1,ITAU))
          END DO
        END DO
C**CALCULATE OVERLAP MATRIX
        DO KP=1,4
          DO KQ=1,4
            SUP5(KP,KQ)=0
            DO I=1,NATOM
              DO K=1,3
                SUP5(KP,KQ)=SUP5(KP,KQ)+XL(I,KP,K,ITAU)*
     1                                  XL(I,KQ,K,ITAU)
              END DO
            END DO
          END DO
        END DO
C**CALCULATE INVERSE OF OVERLAP
C*************************************
        IFAIL=1
        CALL MATINV(SUP5,4,4,SUP4,SUP6,IFAIL)
        IF(IFAIL.NE.0)THEN
          IF(JPRINT.GT.0)WRITE(IOUT,*)'IFAIL = ',IFAIL,' IN S MATINV'
          STOP 'ERROR IN MATINV'
        END IF
        ILHS=0
        DO I=1,NATOM
C**KI LABELS GAMMA (LHS)
          DO KI=1,3
            ILHS=ILHS+1
            IRHS=0
            DO J=1,NATOM
C**KJ LABELS GAMMA (RHS)
              DO KJ=1,3
                IRHS=IRHS+1
                XP(ILHS,IRHS,ITAU)=0
                DO KP=1,4
                  DO KQ=1,4
                    XP(ILHS,IRHS,ITAU)=XP(ILHS,IRHS,ITAU)+
     1              XL(I,KP,KI,ITAU)*XL(J,KQ,KJ,ITAU)*SUP4(KP,KQ)
                  END DO
                END DO
                IF(KI.EQ.KJ)XP(ILHS,IRHS,ITAU)=XP(ILHS,IRHS,ITAU)+
     1          SQRT(XM(I)*XM(J))/XMTOT
              END DO
            END DO
          END DO
        END DO
7777    CONTINUE
C**FORM THE PROJECTED FORCE CONSTANTS AND 3N-7 NORMAL COORDINATES
C****************************************************************
        CALL PROJEC(NATOM,NMODE,XM,XX(1,1,ITAU),OMEGA(1,ITAU),
     1  XL(1,1,1,ITAU),XR,RR,XP(1,1,ITAU),TEMP,PAROT(1,1,ITAU),ITAU)

C**FIRST TIME NORMALS RUN...NO PROJECTION
        IF(ITAU.EQ.1)THEN

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)THEN
            WRITE(IOUT,*)
            WRITE(IOUT,*)
            WRITE(IOUT,*)'EQUILIBRIUM 3N-7 NORMAL MODES'
            WRITE(IOUT,*)
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

          J0=0
          DO J=1,NMODE-1
            J0=J0+1
C**MODE FOR REACTION PATH WILL BE MISSING, SO SKIP
            IF(J0.EQ.IREACT)J0=J0+1

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)THEN
              WRITE(IOUT,*)'OMEGA ',J,': ',OMEGA(J0,1)*WAVENM
              DO I=1,NATOM
                WRITE(IOUT,*)(XL(I,J0,K,1),K=1,3)
              END DO
            END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

          END DO
          GO TO 8888
        END IF
        J0=0
        DO J=1,NMODE-1
          J0=J0+1
C**MODE FOR REACTION PATH WILL BE MISSING, SO SKIP FOR ITAU=1
          IF(J0.EQ.IREACT)J0=J0+1
          IF(ITAU.GT.2)J0=J
C**STANDARDISE THE VECTORS W.R.T. PREVIOUS VECTORS
          IX=1
          IK=1
          DO I=1,NATOM
            DO K=1,3
              IF(DABS(XL(I,J,K,ITAU)).GT.DABS(XL(IX,J,IK,ITAU)))THEN
                IX=I
                IK=K
              END IF
            END DO
          END DO
          IF((XL(IX,J,IK,ITAU).GE.0.AND.XL(IX,J0,IK,ITAU-1).LT.0).OR.
     1      (XL(IX,J,IK,ITAU).LT.0.AND.XL(IX,J0,IK,ITAU-1).GE.0))THEN
            DO I=1,NATOM
              DO K=1,3
                XL(I,J,K,ITAU)=-XL(I,J,K,ITAU)
              END DO
            END DO
          END IF
C**TEST TO SEE IF MODES SWITCHED
          SUM=0
          DO I=1,NATOM
            DO K=1,3
              SUM=SUM+XL(I,J,K,ITAU)*XL(I,J0,K,ITAU-1)
            END DO
          END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)THEN
            WRITE(IOUT,*)'SUM SHOULD BE 1',SUM
            IF(DABS(SUM).LT.0.99D0)WRITE(IOUT,*)
     1      'MODES ',J,J+1,'  SWITCHED ? - SUM',SUM
            IF(DABS(SUM).GT.1.01D0)WRITE(IOUT,*)
     1      'MODES ',J,J+1,'  SWITCHED ? - SUM',SUM
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**ARBITRARY GUESS OF 1.D-03 DEFINES ORTHOGONAL VECTORS
CC        IF(DABS(SUM).LT.1.D-03)THEN
C**TEMPORARY - DIFFICULT VECTOR
          IF(DABS(SUM).LT.0.5D0)THEN
C**VECTORS HAVE SWITCHED

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)WRITE(IOUT,*)
     1      'SWITCHING MODES ',J,' AND ',J+1
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

            SAVE=OMEGA(J,ITAU)
            OMEGA(J,ITAU)=OMEGA(J+1,ITAU)
            OMEGA(J+1,ITAU)=SAVE
            DO I=1,NATOM
              DO K=1,3
                SAVE=XL(I,J,K,ITAU)
                XL(I,J,K,ITAU)=XL(I,J+1,K,ITAU)
                XL(I,J+1,K,ITAU)=SAVE
              END DO
            END DO
C**MUST NOW STANDARDISE NEW (SWITCHED) VECTOR
            IX=1
            IK=1
            DO I=1,NATOM
              DO K=1,3
                IF(DABS(XL(I,J,K,ITAU)).GT.DABS(XL(IX,J,IK,ITAU)))THEN
                  IX=I
                  IK=K
                END IF
              END DO
            END DO
            IF((XL(IX,J,IK,ITAU).GE.0.AND.XL(IX,J0,IK,ITAU-1).LT.0).OR.
     1        (XL(IX,J,IK,ITAU).LT.0.AND.XL(IX,J0,IK,ITAU-1).GE.0))THEN
              DO I=1,NATOM
                DO K=1,3
                  XL(I,J,K,ITAU)=-XL(I,J,K,ITAU)
                END DO
              END DO
            END IF
C**TEST TO SEE IF MODES SWITCHED
            SUM=0
            DO I=1,NATOM
              DO K=1,3
                SUM=SUM+XL(I,J,K,ITAU)*XL(I,J0,K,ITAU-1)
              END DO
            END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)THEN
              IF(DABS(SUM).LT.0.99D0)WRITE(IOUT,*)'SUM SHOULD BE 1',SUM
            END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

          END IF
C**CHECK ALL SIGNS RELATIVE TO PREVIOUS VECTOR
          ICLEAR=0
          DO I=1,NATOM
            DO K=1,3
              IF(XL(I,J,K,ITAU).LT.0.AND.XL(I,J0,K,ITAU-1).GE.0)ICLEAR=1
              IF(XL(I,J,K,ITAU).GE.0.AND.XL(I,J0,K,ITAU-1).LT.0)ICLEAR=1
            END DO
          END DO
          IF(ICLEAR.NE.0)THEN

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)THEN
              WRITE(IOUT,*)'*******************************'
              WRITE(IOUT,*)'FOLLOWING VECTOR LOOKS WRONG FOR MODE ',J
              DO I=1,NATOM
                WRITE(IOUT,*)(XL(I,J,K,ITAU),K=1,3)
              END DO
              WRITE(IOUT,*)
              DO I=1,NATOM
                WRITE(IOUT,*)(XL(I,J0,K,ITAU-1),K=1,3)
              END DO
            END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

          END IF
        END DO
        GO TO 444
666     CONTINUE
C**************************************
C**************************************
C**SET VECTORS CORRECTLY ABOUT '240'
C**************************************
C**************************************
        JNCR=ITAU-122
        IPLUS=122+JNCR
        IMINUS=122-JNCR
        V(IPLUS)=V(IMINUS)
        DO L=1,NMODE-1
          OMEGA(L,IPLUS)=OMEGA(L,IMINUS)
          DOMEGA(L,IPLUS)=DOMEGA(L,IMINUS)
C**FOR THIS MODE, DETERMINE IF IN-PLANE AT '180'
          JX=1
          MINVEC=1
          DO J=1,NATOM
            DO K=1,3
              IF(DABS(XL(J,L,K,INIT)).LT.DABS(XL(JX,L,MINVEC,INIT)))
     1        THEN
                MINVEC=K
                JX=J
              END IF
            END DO
          END DO
          DO J=1,NATOM
            DO K=1,3
              XL(J,L,K,IPLUS)=XL(J,L,K,IMINUS)
              IF(MINVEC.EQ.INITCD)THEN
C**IN-PLANE VECTOR
                IF(K.EQ.INITCD)XL(J,L,K,IPLUS)=-XL(J,L,K,IMINUS)
              ELSE
C**OUT-OF-PLANE VECTOR
                IF(K.NE.INITCD)XL(J,L,K,IPLUS)=-XL(J,L,K,IMINUS)
              END IF
            END DO
          END DO
          DO K=1,3
            XL(3,L,K,IPLUS)=XL(4,L,K,IMINUS)
            XL(4,L,K,IPLUS)=XL(3,L,K,IMINUS)
            IF(MINVEC.EQ.INITCD)THEN
C**IN-PLANE VECTOR
              IF(K.EQ.INITCD)THEN
                XL(3,L,K,IPLUS)=-XL(4,L,K,IMINUS)
                XL(4,L,K,IPLUS)=-XL(3,L,K,IMINUS)
              END IF
            ELSE
C**OUT-OF-PLANE VECTOR
              IF(K.NE.INITCD)THEN
                XL(3,L,K,IPLUS)=-XL(4,L,K,IMINUS)
                XL(4,L,K,IPLUS)=-XL(3,L,K,IMINUS)
              END IF
            END IF
          END DO
        END DO
        GO TO 444
555     CONTINUE
C**************************************
C**************************************
C**SET VECTORS CORRECTLY ABOUT '300'
C**************************************
C**************************************
        JNCR=ITAU-242
        IPLUS=242+JNCR
        IMINUS=242-JNCR
        V(IPLUS)=V(IMINUS)
        DO L=1,NMODE-1
          OMEGA(L,IPLUS)=OMEGA(L,IMINUS)
          DOMEGA(L,IPLUS)=DOMEGA(L,IMINUS)
C**FOR THIS MODE, DETERMINE IF IN-PLANE AT '180'
          JX=1
          MINVEC=1
          DO J=1,NATOM
            DO K=1,3
              IF(DABS(XL(J,L,K,INIT)).LT.DABS(XL(JX,L,MINVEC,INIT)))
     1        THEN
                MINVEC=K
                JX=J
              END IF
            END DO
          END DO
          DO J=1,NATOM
            DO K=1,3
              XL(J,L,K,IPLUS)=XL(J,L,K,IMINUS)
              IF(MINVEC.EQ.INITCD)THEN
C**IN-PLANE VECTOR
                IF(K.EQ.INITCD)XL(J,L,K,IPLUS)=-XL(J,L,K,IMINUS)
              ELSE
C**OUT-OF-PLANE VECTOR
                IF(K.NE.INITCD)XL(J,L,K,IPLUS)=-XL(J,L,K,IMINUS)
              END IF
            END DO
          END DO
          DO K=1,3
            XL(3,L,K,IPLUS)=XL(5,L,K,IMINUS)
            XL(5,L,K,IPLUS)=XL(3,L,K,IMINUS)
            IF(MINVEC.EQ.INITCD)THEN
C**IN-PLANE VECTOR
              IF(K.EQ.INITCD)THEN
                XL(3,L,K,IPLUS)=-XL(5,L,K,IMINUS)
                XL(5,L,K,IPLUS)=-XL(3,L,K,IMINUS)
              END IF
            ELSE
C**OUT-OF-PLANE VECTOR
              IF(K.NE.INITCD)THEN
                XL(3,L,K,IPLUS)=-XL(5,L,K,IMINUS)
                XL(5,L,K,IPLUS)=-XL(3,L,K,IMINUS)
              END IF
            END IF
          END DO
        END DO
444     CONTINUE

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
C**FINAL VECTORS
          WRITE(IOUT,*)
          WRITE(IOUT,*)'FINAL VECTORS'
          DO J=1,NMODE-1
            WRITE(IOUT,*)'OMEGA ',J,': ',OMEGA(J,ITAU)*WAVENM
            DO I=1,NATOM
              WRITE(IOUT,*)(XL(I,J,K,ITAU),K=1,3)
            END DO
          END DO
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        DO J=1,NMODE-1
          DOMEGA(J,ITAU)=OMEGA(J,ITAU)
        END DO
8888    CONTINUE
      END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************END PROJECTION ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**************************************
C**************************************
C**SET VECTORS CORRECTLY ABOUT '360'
C**************************************
C**************************************
      IPLUS=LNIT
      IMINUS=LNIT
      DO I=1,360
        IPLUS=IPLUS+1
        IF(IPLUS.GT.722)IPLUS=3
        IMINUS=IMINUS-1
        IF(IMINUS.LT.2)IMINUS=721
        V(IPLUS)=V(IMINUS)
        DO L=1,NMODE-1
          OMEGA(L,IPLUS)=OMEGA(L,IMINUS)
          DOMEGA(L,IPLUS)=DOMEGA(L,IMINUS)
C**FOR THIS MODE, DETERMINE IF IN-PLANE AT '180'
          JX=1
          MINVEC=1
          DO J=1,NATOM
            DO K=1,3
              IF(DABS(XL(J,L,K,INIT)).LT.DABS(XL(JX,L,MINVEC,INIT)))
     1        THEN
                MINVEC=K
                JX=J
              END IF
            END DO
          END DO
          DO J=1,NATOM
            DO K=1,3
              XL(J,L,K,IPLUS)=XL(J,L,K,IMINUS)
              IF(MINVEC.EQ.INITCD)THEN
C**IN-PLANE VECTOR
                IF(K.EQ.INITCD)XL(J,L,K,IPLUS)=-XL(J,L,K,IMINUS)
              ELSE
C**OUT-OF-PLANE VECTOR
                IF(K.NE.INITCD)XL(J,L,K,IPLUS)=-XL(J,L,K,IMINUS)
              END IF
            END DO
          END DO
          DO K=1,3
            XL(4,L,K,IPLUS)=XL(5,L,K,IMINUS)
            XL(5,L,K,IPLUS)=XL(4,L,K,IMINUS)
            IF(MINVEC.EQ.INITCD)THEN
C**IN-PLANE VECTOR
              IF(K.EQ.INITCD)THEN
                XL(4,L,K,IPLUS)=-XL(5,L,K,IMINUS)
                XL(5,L,K,IPLUS)=-XL(4,L,K,IMINUS)
              END IF
            ELSE
C**OUT-OF-PLANE VECTOR
              IF(K.NE.INITCD)THEN
                XL(4,L,K,IPLUS)=-XL(5,L,K,IMINUS)
                XL(5,L,K,IPLUS)=-XL(4,L,K,IMINUS)
              END IF
            END IF
          END DO
        END DO

C**USE VECTORS OF PROJECTED MODES TO TRANSFORM DERIVS WRT X,Y,Z 
C**TO DERIVS WRT Q IN DERIV(NMODE-1,1)
        DO L=1,NMODE-1
          DERIV(L,IPLUS)=0
          DERIV(L,IMINUS)=0
          DO J=1,NATOM
            DO K=1,3
              DERIV(L,IPLUS)=DERIV(L,IPLUS)+
     1        DERIVS(J,K,IPLUS)*XL(J,L,K,IPLUS)/SQRT(XM(J))
              DERIV(L,IMINUS)=DERIV(L,IMINUS)+
     1        DERIVS(J,K,IMINUS)*XL(J,L,K,IMINUS)/SQRT(XM(J))
            END DO
          END DO
        END DO
      END DO
C*******************************************************************
C*******************************************************************

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************START VECTOR DERIVATIVE ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      DO ITAU=2,722
C********************************************
C**NON-ECKART PATH IN DEGREES (1/2 DEGREE)
        SPLUS(ITAU)=0.5D0
        SMINUS(ITAU)=0.5D0
C**THEREFORE JACOBIAN = RAD
        DSTAU(ITAU)=(SPLUS(ITAU)+SMINUS(ITAU))*RAD
C********************************************

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'****************NEXT POINT********************'
          WRITE(IOUT,*)
          WRITE(IOUT,*)'ITAU= ',ITAU
          WRITE(IOUT,*)
          WRITE(IOUT,*)'GEOMETRY AT POINT ITAU = ',ITAU
          DO I=1,NATOM
            WRITE(IOUT,*)(XX(I,K,ITAU)/SQRT(XM(I)),K=1,3)
          END DO
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**UN-NORMALISE XXP
        DO I=1,NATOM
          DO K=1,3
            XXP(I,K,ITAU)=XXP(I,K,ITAU)*XNORM(ITAU)
          END DO
        END DO
C**CALCULATE GRADIENT OF VECTOR AT TAU = ITAU
C********************************************
        DO J=1,NMODE-1
          DO I=1,NATOM
            DO K=1,3
              XLP(I,J,K,ITAU)=(XL(I,J,K,ITAU+1)-XL(I,J,K,ITAU-1))

              IF(ITAU.EQ.2)
     1        XLP(I,J,K,ITAU)=(XL(I,J,K,ITAU+1)-XL(I,J,K,ITAU))*2
              IF(ITAU.EQ.722)
     1        XLP(I,J,K,ITAU)=(XL(I,J,K,ITAU)-XL(I,J,K,ITAU-1))*2

              XLP(I,J,K,ITAU)=XLP(I,J,K,ITAU)/
     1        (SPLUS(ITAU)+SMINUS(ITAU))
            END DO
          END DO
C**NORMALISE ???
          Z=0
          DO I=1,NATOM
            DO K=1,3
              Z=Z+XLP(I,J,K,ITAU)**2
            END DO
          END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)THEN
            WRITE(IOUT,*)
            Z=SQRT(Z)
            IF(Z.LT.1.D-03)THEN
              WRITE(IOUT,*)'BEWARE!!!  NORM = ',Z
            ELSE
              WRITE(IOUT,*)'POSSIBLE!  NORM = ',Z
            END IF
            WRITE(IOUT,*)'OMEGA',OMEGA(J,ITAU)*WAVENM
            WRITE(IOUT,*)'VECTORS'
            DO I=1,NATOM
              WRITE(IOUT,*)(XL(I,J,K,ITAU),K=1,3)
            END DO
            WRITE(IOUT,*)'VECTOR GRADS'
            DO I=1,NATOM
              WRITE(IOUT,*)(XLP(I,J,K,ITAU),K=1,3)
            END DO
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C
          SUM1=0
          DO I=1,NATOM
            SUM1=SUM1+XX(I,2,ITAU)*XL(I,J,3,ITAU)-XX(I,3,ITAU)*
     1      XL(I,J,2,ITAU)
          END DO
          SUM2=0
          DO I=1,NATOM
            SUM2=SUM2+XX(I,3,ITAU)*XL(I,J,1,ITAU)-XX(I,1,ITAU)*
     1      XL(I,J,3,ITAU)
          END DO
          SUM3=0
          DO I=1,NATOM
            SUM3=SUM3+XX(I,1,ITAU)*XL(I,J,2,ITAU)-XX(I,2,ITAU)*
     1      XL(I,J,1,ITAU)
          END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)
     1    WRITE(IOUT,*)'VECTOR PRODUCT XX*XL ',SUM1,SUM2,SUM3
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C
          SUM1=0
          DO I=1,NATOM
            SUM1=SUM1+XXP(I,2,ITAU)*XL(I,J,3,ITAU)-XXP(I,3,ITAU)*
     1      XL(I,J,2,ITAU)
          END DO
          SUM2=0
          DO I=1,NATOM
            SUM2=SUM2+XXP(I,3,ITAU)*XL(I,J,1,ITAU)-XXP(I,1,ITAU)*
     1      XL(I,J,3,ITAU)
          END DO
          SUM3=0
          DO I=1,NATOM
            SUM3=SUM3+XXP(I,1,ITAU)*XL(I,J,2,ITAU)-XXP(I,2,ITAU)*
     1      XL(I,J,1,ITAU)
          END DO
          SUM1=SUM1
          SUM2=SUM2
          SUM3=SUM3

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)
     1    WRITE(IOUT,*)'VECTOR PRODUCT XXP*XL ',SUM1,SUM2,SUM3
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C
          SUM1P=0
          DO I=1,NATOM
            SUM1P=SUM1P+XX(I,2,ITAU)*XLP(I,J,3,ITAU)-XX(I,3,ITAU)*
     1      XLP(I,J,2,ITAU)
          END DO
          SUM2P=0
          DO I=1,NATOM
            SUM2P=SUM2P+XX(I,3,ITAU)*XLP(I,J,1,ITAU)-XX(I,1,ITAU)*
     1      XLP(I,J,3,ITAU)
          END DO
          SUM3P=0
          DO I=1,NATOM
            SUM3P=SUM3P+XX(I,1,ITAU)*XLP(I,J,2,ITAU)-XX(I,2,ITAU)*
     1      XLP(I,J,1,ITAU)
          END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)
     1    WRITE(IOUT,*)'VECTOR PRODUCT XX*XLP ',SUM1P,SUM2P,SUM3P
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

          IF(JPRINT.GT.0)THEN
            IF(DABS(SUM1)-DABS(SUM1P).GT.1.D-3)WRITE(IOUT,*)'SUSPECT 1'
            IF(DABS(SUM2)-DABS(SUM2P).GT.1.D-3)WRITE(IOUT,*)'SUSPECT 2'
            IF(DABS(SUM3)-DABS(SUM3P).GT.1.D-3)WRITE(IOUT,*)'SUSPECT 3'
          END IF
        END DO
      END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************END VECTOR DERIVATIVE ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**************************************
C**************************************
C**FORM WORKING QUANTITIES
C**************************************
C**************************************
      DO ITAU=2,722
        DO J=1,NMODE-1
          BSS(J,ITAU)=0
          DO K=1,3
            BS(J,K,ITAU)=0
            DO L=1,3
              B(J,K,L,ITAU)=0
            END DO
            DO L=1,NMODE-1
              BB(L,J,K,ITAU)=0
            END DO
          END DO
          DO L=1,NMODE-1
            BBS(L,J,ITAU)=0
          END DO
        END DO
        DO K=1,NMODE-1
          DO L=1,NMODE-1
            DO I=1,NATOM
              BB(K,L,1,ITAU)=BB(K,L,1,ITAU)+XL(I,K,2,ITAU)*
     1        XL(I,L,3,ITAU)-XL(I,K,3,ITAU)*XL(I,L,2,ITAU)
              BB(K,L,2,ITAU)=BB(K,L,2,ITAU)+XL(I,K,3,ITAU)*
     1        XL(I,L,1,ITAU)-XL(I,K,1,ITAU)*XL(I,L,3,ITAU)
              BB(K,L,3,ITAU)=BB(K,L,3,ITAU)+XL(I,K,1,ITAU)*
     1        XL(I,L,2,ITAU)-XL(I,K,2,ITAU)*XL(I,L,1,ITAU)
              DO J=1,3
                BBS(K,L,ITAU)=BBS(K,L,ITAU)+XLP(I,K,J,ITAU)*
     1          XL(I,L,J,ITAU)
              END DO
            END DO
          END DO
          DO J=1,3
            DO L=1,3
              DO I=1,NATOM
                B(K,J,J,ITAU)=B(K,J,J,ITAU)+XL(I,K,L,ITAU)*XX(I,L,ITAU)
                B(K,J,L,ITAU)=B(K,J,L,ITAU)-XL(I,K,J,ITAU)*XX(I,L,ITAU)
              END DO
            END DO
          END DO
          DO I=1,NATOM
            BS(K,1,ITAU)=BS(K,1,ITAU)+XXP(I,3,ITAU)*XL(I,K,2,ITAU)-
     1      XXP(I,2,ITAU)*XL(I,K,3,ITAU)
            BS(K,2,ITAU)=BS(K,2,ITAU)+XXP(I,1,ITAU)*XL(I,K,3,ITAU)-
     1      XXP(I,3,ITAU)*XL(I,K,1,ITAU)
            BS(K,3,ITAU)=BS(K,3,ITAU)+XXP(I,2,ITAU)*XL(I,K,1,ITAU)-
     1      XXP(I,1,ITAU)*XL(I,K,2,ITAU)
            DO L=1,3
              BSS(K,ITAU)=BSS(K,ITAU)+XLP(I,K,L,ITAU)*XXP(I,L,ITAU)
            END DO
          END DO
        END DO

C**NCH TEST
        IF(ITAU.GT.3.AND.ITAU.LT.721)THEN
          DO I=1,NMODE-1
            DO J=1,I
              IF(I.NE.J)THEN
                IF(DABS(BBS(I,J,ITAU)+BBS(J,I,ITAU)).GT.0.01D0)THEN
                  IF(JPRINT.GT.0)THEN
                    WRITE(IOUT,*)'THE FOLLOWING SHOULD ADD TO ZERO'
                    WRITE(IOUT,*)ITAU,I,J,BBS(I,J,ITAU),BBS(J,I,ITAU)
                  END IF
                END IF
              END IF
            END DO
          END DO
        END IF
C**NCH TEST

      END DO
C**************************************
C**************************************
      IF(ISCFCI.LT.0)STOP 'CURVELINEAR'
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETSAY(NATOM,X0,XM,XX,XR,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XM(NATOM),XX(NATOM,3),X0(NATOM,3)
      DIMENSION XR(NATOM,3),RR(NATOM,NATOM)
      COMMON/RETURN/IRET
      COMMON/PLMNMX/MXCD(3)
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
C*****************************************************
200   FORMAT(//,1X,'DISPLACED INTERNAL COORDINATES',/)
210   FORMAT(//,1X,'ECKART INTERNAL COORDINATES',/)
220   FORMAT(//,1X,'DISTORTED INTERNAL COORDINATES',/)
240   FORMAT(1X,'ATOM ',I2,':  X0 =',F12.7,'  Y0 =',F12.7,'  Z0 =',
     1F12.7)
245   FORMAT(//,60(1H*),/,1X,'DISPLACED GEOMETRY',2X,I4/)
250   FORMAT(//,1X,'ECKART GEOMETRY',/)
255   FORMAT(//,1X,'ALPHA,BETA,GAMMA',/)
260   FORMAT(//,1X,'PREVIOUS GEOMETRY',/)
C*****************************************************

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'*********************START ROTATION OF GEOMETRY'
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C***********************************
C**XR CONTAINS DISTORTED GEOMETRY
      CALL CHECKM(XM,XR,XX,RR,NATOM)
C**XX CONTAINS DISTORTED GEOMETRY AT CENTRE OF MASS
C**XR CONTAINS NEW (TRANSLATED, ROTATED) GEOMETRY
C***********************************
C**XR CONTAINS DISTORTED GEOMETRY
      DO I=1,NATOM
        DO K=1,3
          IF(IRET.NE.0)XX(I,K)=XR(I,K)
        END DO
      END DO
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,260)
        DO I=1,NATOM
          WRITE(IOUT,240)I,(X0(I,J),J=1,3)
        END DO
      END IF
C**LOOK FOR LARGEST COORDINATE (3)
      XMAXC=0
      DO J=1,3
        DO I=1,NATOM
          IF(DABS(XX(I,J)).GT.XMAXC)THEN
            XMAXC=DABS(XX(I,J))
            MXCD(3)=J
          END IF
        END DO
      END DO
C**LOOK FOR SMALLEST COORDINATE (1)
      XMINC=1000000
      DO J=1,3
        DO I=1,NATOM
          IF(DABS(XX(I,J)).LT.XMINC)THEN
            XMINC=DABS(XX(I,J))
            MXCD(1)=J
          END IF
        END DO
      END DO
C**SET REMAINING COORDINATE
      IF(MXCD(1).EQ.1)THEN
        IF(MXCD(3).EQ.2)THEN
          MXCD(2)=3
        ELSE
          MXCD(2)=2
        END IF
      END IF
      IF(MXCD(1).EQ.2)THEN
        IF(MXCD(3).EQ.3)THEN
          MXCD(2)=1
        ELSE
          MXCD(2)=3
        END IF
      END IF
      IF(MXCD(1).EQ.3)THEN
        IF(MXCD(3).EQ.1)THEN
          MXCD(2)=2
        ELSE
          MXCD(2)=1
        END IF
      END IF
      IF(JPRINT.GT.0)WRITE(IOUT,*)'MXCD ',(MXCD(I),I=1,3)
      RETURN
      END
