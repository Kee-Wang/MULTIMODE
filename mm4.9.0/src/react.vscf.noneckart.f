C****************************************************************
C****************************************************************
C**REACTION PATH (ANALYTICAL POTENTIAL....NON-ECKART)
C****************************************************************
C****************************************************************
      SUBROUTINE REACT(NATOM,NMODE,XM,X0,OMEGA,XL,XLP,XX,XXP,RR,XR,XP,
     1TEMP,BB,BBS,B,BS,BSS,IREACT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ECKART
C***************************************************************
C***************************************************************
C**REACTION PATH (NON-ECKART VARIATION) TAYLORED TO HYDROGENPEROXIDE
C**ALLOW ONE POINT FOR EQUILIBRIUM (1)
C**ALLOWS FOR ONE POINT EVERY HALF DEGREE (721)
C**ALLOW EXTRA POINT FOR EQUILIBRIUM 2ND TIME (722)
C**
C**WORK IN RADIANS
C**
C***************************************************************
C***************************************************************
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
C**
      DIMENSION IGOT(30,3),IPLAT(30)
C**
      COMMON/RETURN/IRET
      COMMON/CHECK/MCHECK
      COMMON/PLMNMX/MXCD(3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
220   FORMAT(/,1X,'MODE',I4,'  OMEGA = ',F20.12,/)
245   FORMAT(//,1X,'DISPLACEMENTS OF NORMAL MODES, L(alpha i,k)')
255   FORMAT(1X,'ATOM(i) ',I2,':  x =',F12.6,'  y =',F12.6,'  z =',
     1F12.6)
260   FORMAT(/,1X,'PRINCIPAL AXIS GEOMETRY',/)
265   FORMAT(1X,'ATOM ',I2,':  X0 =',F12.7,'  Y0 =',F12.7,'  Z0 =',
     1F12.7)
270   FORMAT(/,1X,'PREVIOUS GEOMETRY',/)
C*******************************************************************
      ECKART=.FALSE.

      DO I=1,3
        MXCD(I)=0
      END DO
      IRET=-1
C*******************************************************************
C*******************************************************************
C**MAKE SURE TRUE EQUILIBRIUM ORIENTATED CORRECTLY IF REQUIRED
C**X0 CONTAINS EQUILIBRIUM GEOMETRY WHICH SATISFIES EVERYTHING
      IF(MCHECK.NE.0)CALL CHECKM(XM,X0,XX,RR,NATOM)

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'****************EQUILIBRIUM COORDINATE PARAMETERS'
        WRITE(IOUT,*)
        CALL GETPOT(V,NATOM,X0,RR)
        DO I=1,NATOM
          WRITE(IOUT,265)I,(X0(I,K),K=1,3)
        END DO
        WRITE(IOUT,*)
        WRITE(IOUT,*)'INITIAL V =',V*WAVENM
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**IN TERMS OF CURVELINEAR COORDINATE (TAU), EQUILIBRIUM IS AT TAU=0
      TAU=0
C**GET EQUILIBRIUM TAUE
      CALL MINPOT(TAU,NATOM,XR,RR,TAUE)
C*******************************************************************
C*******************************************************************
C**GET INCREMENT TO TRANS
      INCR=180
      TAUINC=180-TAUE*RAD
C**JNIT WILL POINT AT 'TRANS'
      JNIT=0
C**INIT WILL POINT AT 'CIS'
      INIT=0

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'TAUINC',TAUINC
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************START REACTION COORDINATE ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**FIRST LOOP FROM 0 TO 360 DEGREES, INCREMENTS OF 1/2 DEGREE
      VMIN=1.D+20
      DO ITAU=0,721
C**ONLY GO AS FAR AS 'CIS'
        IF(INIT.NE.0)GO TO 9999

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'**********************************'
          WRITE(IOUT,*)'GEOMETRIES AT POINT ITAU = ',ITAU+1
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        IF(ITAU.EQ.0)THEN
          TAU=TAU+TAUINC/RAD
        ELSE
          IF(ITAU.GT.1)TAU=TAU+0.5D0/RAD
        END IF
C**GET (UNROTATED) CARTESIAN COORDINATES FOR NEXT VALUE OF TAU IN XR
        CALL MINPOT(TAU,NATOM,XR,RR,TAUX)
        CALL GETPOT(V,NATOM,XR,RR)
C**FIND ITAU CORRESPONDING TO MINIMUM ENERGY
        IF(V.LT.VMIN)THEN
          VMIN=V
          ITAUMN=ITAU+1
        END IF
C**NOW MOVE TO CENTRE OF MASS AND ROTATE TO PRINCIPAL AXES
        CALL GETSAY(NATOM,XM,XX(1,1,ITAU+1),XR,RR)

C**NORMALISE
        IF(ITAU.EQ.0)THEN
          Z=0
          DO I=1,NATOM
            DO K=1,3
              Z=Z+XX(I,K,1)*XX(I,K,1)
            END DO
          END DO
          XNORM(1)=SQRT(Z)
          GO TO 9898
        END IF
C**NORMALISE

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
        END IF

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C**PREVIOUS INTERNAL DISPLACEMENT COORDINATES
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,270)
          DO I=1,NATOM
            WRITE(IOUT,265)I,(XX(I,K,ITAU),K=1,3)
          END DO
C**NEW INTERNAL DISPLACEMENT COORDINATES ALONG PRINCIPAL AXES
          WRITE(IOUT,260)
          DO I=1,NATOM
            WRITE(IOUT,265)I,(XX(I,K,ITAU+1),K=1,3)
          END DO
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**NORMALISED TO PREVIOUS?
        Z=0
        DO I=1,NATOM
          DO K=1,3
            Z=Z+XX(I,K,ITAU)*XX(I,K,ITAU+1)
          END DO
        END DO
        Z=SQRT(Z)/XNORM(ITAU)
        IF(DABS(Z-1.D0).GT.1.D-02)THEN
          IF(JPRINT.GT.0)WRITE(IOUT,*)'FAILURE!!!  Z = ',Z
C**TRY TO SORT IT OUT
C********************
C**FIRST CHECK 'X','Y','Z' COORDS...IND's ARE ZERO IF OK
          CALL NORXYZ(XX,NATOM,ITAU,INDX,1)
          CALL NORXYZ(XX,NATOM,ITAU,INDY,2)
          CALL NORXYZ(XX,NATOM,ITAU,INDZ,3)
          IF(JPRINT.GT.0)WRITE(IOUT,*)'INDX,INDY,INDZ = ',INDX,INDY,INDZ
          IF(INDX.EQ.0.AND.INDY.EQ.0.AND.INDZ.EQ.0)THEN
            IF(JPRINT.GT.0)WRITE(IOUT,*)'DISASTER - ALL SEEMS FINE!'
          END IF
          IF(INDX.NE.0.AND.INDY.NE.0.AND.INDZ.NE.0)THEN
            IF(JPRINT.GT.0)WRITE(IOUT,*)'DISASTER - ALL SEEMS WRONG!'
          END IF
C**INTERCHANGE REQUIRED COORDINATES
C**CURRENT 'X' IS CORRECT
          IF(INDX.EQ.0)THEN

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)WRITE(IOUT,*)'INTERCHANGE Y AND Z'
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

            DO I=1,NATOM
              TTTT=XX(I,2,ITAU+1)
              XX(I,2,ITAU+1)=XX(I,3,ITAU+1)
              XX(I,3,ITAU+1)=TTTT
            END DO
          END IF
C**CURRENT 'Y' IS CORRECT
          IF(INDY.EQ.0)THEN

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)WRITE(IOUT,*)'INTERCHANGE X AND Z'
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

            DO I=1,NATOM
              TTTT=XX(I,1,ITAU+1)
              XX(I,1,ITAU+1)=XX(I,3,ITAU+1)
              XX(I,3,ITAU+1)=TTTT
            END DO
          END IF
C**CURRENT 'Z' IS CORRECT
          IF(INDZ.EQ.0)THEN

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)WRITE(IOUT,*)'INTERCHANGE X AND Y'
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

            DO I=1,NATOM
              TTTT=XX(I,1,ITAU+1)
              XX(I,1,ITAU+1)=XX(I,2,ITAU+1)
              XX(I,2,ITAU+1)=TTTT
            END DO
          END IF
C**CHECK PHASES AGAIN
C**LOOK FOR BIGGEST X-COORDINATE
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
          END IF
C**LOOK FOR BIGGEST Y-COORDINATE
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
          END IF
C**LOOK FOR BIGGEST Z-COORDINATE
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
          END IF
C**NEW INTERNAL DISPLACEMENT COORDINATES ALONG PRINCIPAL AXES

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)THEN
            WRITE(IOUT,260)
            DO I=1,NATOM
              WRITE(IOUT,265)I,(XX(I,K,ITAU+1),K=1,3)
            END DO
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**NORMALISED TO PREVIOUS?
          Z=0
          DO I=1,NATOM
            DO K=1,3
              Z=Z+XX(I,K,ITAU)*XX(I,K,ITAU+1)
            END DO
          END DO
          Z=SQRT(Z)/XNORM(ITAU)
          IF(DABS(Z-1.D0).GT.1.D-02)THEN
            IF(JPRINT.GT.0)
     1      WRITE(IOUT,*)'DISASTER - STILL SEEMS WRONG!  Z = ',Z
          ELSE
            IF(JPRINT.GT.0)WRITE(IOUT,*)'SUCCESS!!!  Z = ',Z
            GO TO 9998
          END IF
C********************
C**SECOND, CHECK ATOMS...IPLAT's ARE ZERO IF OK
          DO IABC=1,NATOM
            CALL NOATOM(XX,NATOM,ITAU,IPLAT,IABC)
          END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)WRITE(IOUT,*)'IPLAT = ',(IPLAT(I),I=1,NATOM)
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

          IND=0
          DO I=1,NATOM
            IF(IPLAT(I).NE.0)IND=IND+1
          END DO
C**ALL ATOMS APPEAR NORMALISED
          IF(IND.EQ.0)THEN
            IF(JPRINT.GT.0)
     1      WRITE(IOUT,*)'DOUBLE DISASTER - ALL SEEMS FINE!'
          END IF
C**NO ATOM APPEARS NORMALISED
          IF(IND.EQ.NATOM)THEN
            IF(JPRINT.GT.0)
     1      WRITE(IOUT,*)'DOUBLE DISASTER - ALL SEEMS WRONG!'
          END IF
C**FIRST MAKE SURE PHASES OF CORRECT ATOMS ARE THE SAME
C**LOOK FOR BIGGEST X-COORDINATE
          ICH=0
          IXX=0
          DO I=1,NATOM
            IF(IPLAT(I).EQ.0.AND.IXX.EQ.0)IXX=I
          END DO
          IX=IXX
          DO I=IXX+1,NATOM
            IF(IPLAT(I).EQ.0.AND.(DABS(XX(I,1,ITAU+1)).GT.
     1      DABS(XX(IX,1,ITAU+1))))IX=I
          END DO
          IF(XX(IX,1,ITAU+1).GE.0.AND.XX(IX,1,ITAU).LT.0)ICH=1
          IF(XX(IX,1,ITAU+1).LE.0.AND.XX(IX,1,ITAU).GT.0)ICH=1
          IF(ICH.NE.0)THEN
            DO I=1,NATOM
              XX(I,1,ITAU+1)=-XX(I,1,ITAU+1)
            END DO
          END IF
C**LOOK FOR BIGGEST Y-COORDINATE
          ICH=0
          IYY=0
          DO I=1,NATOM
            IF(IPLAT(I).EQ.0.AND.IYY.EQ.0)IYY=I
          END DO
          IY=IYY
          DO I=IYY+1,NATOM
            IF(IPLAT(I).EQ.0.AND.(DABS(XX(I,2,ITAU+1)).GT.
     1      DABS(XX(IY,2,ITAU+1))))IY=I
          END DO
          IF(XX(IY,2,ITAU+1).GE.0.AND.XX(IY,2,ITAU).LT.0)ICH=1
          IF(XX(IY,2,ITAU+1).LE.0.AND.XX(IY,2,ITAU).GT.0)ICH=1
          IF(ICH.NE.0)THEN
            DO I=1,NATOM
              XX(I,2,ITAU+1)=-XX(I,2,ITAU+1)
            END DO
          END IF
C**LOOK FOR BIGGEST Z-COORDINATE
          ICH=0
          IZZ=0
          DO I=1,NATOM
            IF(IPLAT(I).EQ.0.AND.IZZ.EQ.0)IZZ=I
          END DO
          IZ=IZZ
          DO I=IZZ+1,NATOM
            IF(IPLAT(I).EQ.0.AND.(DABS(XX(I,3,ITAU+1)).GT.
     1      DABS(XX(IZ,3,ITAU+1))))IZ=I
          END DO
          IF(XX(IZ,3,ITAU+1).GE.0.AND.XX(IZ,3,ITAU).LT.0)ICH=1
          IF(XX(IZ,3,ITAU+1).LE.0.AND.XX(IZ,3,ITAU).GT.0)ICH=1
          IF(ICH.NE.0)THEN
            DO I=1,NATOM
              XX(I,3,ITAU+1)=-XX(I,3,ITAU+1)
            END DO
          END IF
C**CHECK FOR NON-PLANAR COORDINATES
        DO I=1,NATOM
          DO K=1,3
            IGOT(I,K)=0
          END DO
        END DO
        DO I=1,NATOM
          IF(IPLAT(I).EQ.0)GO TO 9991
          DO J=1,NATOM
            IF((I.EQ.J).OR.(IPLAT(J).EQ.0))GO TO 9990

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)WRITE(IOUT,*)'CHECKING ATOMS ',I,' AND ',J
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

            DO K=1,3
C**CAN'T USE ATOM J
              IF(IGOT(J,K).NE.0)GO TO 9990
C**CAN'T USE ATOM I
              IF(IGOT(I,K).NE.0)GO TO 9991
C**ITAU=2 DEFINES FIRST SYMMETRIC STRUCTURE (ALWAYS!!!)
              IF(DABS(XX(I,K,2)-XX(J,K,2)).LT.1.D-6)THEN
                IGOT(I,K)=1
                IGOT(J,K)=1
              END IF
              IF(DABS(XX(I,K,2)+XX(J,K,2)).LT.1.D-6)THEN
                IGOT(I,K)=-1
                IGOT(J,K)=-1
              END IF
9989          CONTINUE
            END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)THEN
              WRITE(IOUT,*)'ATOM ',I,' IGOT = ',(IGOT(I,K),K=1,3)
              WRITE(IOUT,*)'ATOM ',J,' IGOT = ',(IGOT(J,K),K=1,3)
            END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**EVEN A SINGLE NON-QUALIFYING COORD. CONSTITUTES TEST-FAILURE
            DO K=1,3
              IF(IGOT(J,K).EQ.0)GO TO 9990
            END DO
C**MUST BE C2-EQUIVALENT (TWO COORDINATES SAME, ONE OPPOSITE)
            ISUM=0
            DO K=1,3
              ISUM=ISUM+IGOT(J,K)
            END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)
     1      WRITE(IOUT,*)'ISUM = ',ISUM,' FOR ATOMS ',I,' AND ',J
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**NO GOOD IF FAILS TEST - RESET IGOT TO ZERO
            IF(ISUM.LE.0)THEN
              DO K=1,3
                IGOT(I,K)=0
                IGOT(J,K)=0
              END DO
              GO TO 9990
            END IF
C**INTERCHANGE REQUIRED ATOMS

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            IF(JPRINT.GT.0)WRITE(IOUT,*)'INTERCHANGE ATOMS ',I,' AND ',J
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

            DO K=1,3
              TTTT=XX(I,K,ITAU+1)
              XX(I,K,ITAU+1)=XX(J,K,ITAU+1)
              XX(J,K,ITAU+1)=TTTT
            END DO
9990        CONTINUE
          END DO
9991      CONTINUE
        END DO
C**CHECK PHASES AGAIN
C**LOOK FOR BIGGEST X-COORDINATE
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
          END IF
C**LOOK FOR BIGGEST Y-COORDINATE
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
          END IF
C**LOOK FOR BIGGEST Z-COORDINATE
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
          END IF
C**NEW INTERNAL DISPLACEMENT COORDINATES ALONG PRINCIPAL AXES

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)THEN
            WRITE(IOUT,260)
            DO I=1,NATOM
              WRITE(IOUT,265)I,(XX(I,K,ITAU+1),K=1,3)
            END DO
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**NORMALISED TO PREVIOUS?
          Z=0
          DO I=1,NATOM
            DO K=1,3
              Z=Z+XX(I,K,ITAU)*XX(I,K,ITAU+1)
            END DO
          END DO
          Z=SQRT(Z)/XNORM(ITAU)
          IF(DABS(Z-1.D0).GT.1.D-02)THEN
            IF(JPRINT.GT.0)
     1      WRITE(IOUT,*)'DOUBLE DISASTER - STILL SEEMS WRONG!  Z = ',Z
          ELSE
            IF(JPRINT.GT.0)WRITE(IOUT,*)'SUCCESS!!!  Z = ',Z
          END IF
        ELSE
          IF(JPRINT.GT.0)WRITE(IOUT,*)'SUCCESS!!!  Z = ',Z
        END IF
9998    CONTINUE
C**POINT TO 180 DEGREES 
        IF(JNIT.EQ.0.AND.ITAU.GT.0)THEN
          IF(INCR+(ITAU-1)/2.EQ.180)JNIT=ITAU+1
          IF(JNIT.NE.0)THEN
            IF(JPRINT.GT.0)WRITE(IOUT,*)'180 DEGREES IS AT POINT ',JNIT
            JNITCD=MXCD(1)
            IF(JPRINT.GT.0)WRITE(IOUT,*)'ZERO COORDINATE IS ',JNITCD
          END IF
        END IF
C**POINT TO ZERO DEGREES (360 DEGREES)
        IF(INIT.EQ.0)THEN
          IF(INCR+(ITAU-1)/2.EQ.360)INIT=ITAU+1
          IF(INIT.NE.0)THEN
            IF(JPRINT.GT.0)WRITE(IOUT,*)'360 DEGREES IS AT POINT ',INIT
            INITCD=MXCD(1)
            IF(JPRINT.GT.0)WRITE(IOUT,*)'ZERO COORDINATE IS ',INITCD
          END IF
        END IF
C**NOW FORM NORMALISATION CONSTANT
        Z=0
        DO I=1,NATOM
          DO K=1,3
            Z=Z+XX(I,K,ITAU+1)*XX(I,K,ITAU+1)
          END DO
        END DO
        XNORM(ITAU+1)=SQRT(Z)

9898    CONTINUE
C*********************************************************
C*********************************************************
C**LOOK FOR LARGEST COORDINATE (3)
        XMAXC=0
        DO J=1,3
          DO I=1,NATOM
            IF(DABS(XX(I,J,ITAU+1)).GT.XMAXC)THEN
              XMAXC=DABS(XX(I,J,ITAU+1))
              MXCD(3)=J
            END IF
          END DO
        END DO
C**LOOK FOR SMALLEST COORDINATE (1)
        XMINC=1000000
        DO J=1,3
          DO I=1,NATOM
            IF(DABS(XX(I,J,ITAU+1)).LT.XMINC)THEN
              XMINC=DABS(XX(I,J,ITAU+1))
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

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'MXCD ',(MXCD(I),I=1,3)
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
9999    CONTINUE
      END DO
C**AT END OF LOOP, CURRENT POINT SHOULD BE EQUILIBRIUM
C**XX(2) CORRESPONDS TO XX(722), BUT BEWARE ROTATION OF PRINCIPAL AXES

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
C**SET GEOMETRIES CORRECTLY ABOUT 'CIS'
C**************************************
C**************************************
      IPLUS=INIT
      IMINUS=INIT
      DO ITAU=1,360
        IPLUS=IPLUS+1
        IF(IPLUS.GT.722)IPLUS=3
        DO I=1,NATOM
          DO K=1,3
            XX(I,K,IPLUS)=0
          END DO
        END DO
        IMINUS=IMINUS-1
        IF(IMINUS.LT.2)IMINUS=721
        IPLAN=0
        DO I=1,30
          IPLAT(I)=0
          DO K=1,3
            IGOT(I,K)=0
          END DO
        END DO
C**FIRST CHECK PLANAR COORDINATES
        DO I=1,NATOM
          IF(DABS(XX(I,1,INIT)).LT.1.D-6)THEN
            IF(DABS(XX(I,1,INIT-1)).GT.1.D-6)THEN
C**PLANAR COORDINATE IS 'X'
              IPLAN=1
              IPLAT(I)=I
            END IF
            XX(I,1,IPLUS)=-XX(I,1,IMINUS)
            XX(I,2,IPLUS)=XX(I,2,IMINUS)
            XX(I,3,IPLUS)=XX(I,3,IMINUS)
            DO K=1,3
              IGOT(I,K)=K
            END DO
          END IF
          IF(DABS(XX(I,2,INIT)).LT.1.D-6)THEN
            IF(DABS(XX(I,2,INIT-1)).GT.1.D-6)THEN
C**PLANAR COORDINATE IS 'Y'
              IPLAN=2
              IPLAT(I)=I
            END IF
            XX(I,2,IPLUS)=-XX(I,2,IMINUS)
            XX(I,1,IPLUS)=XX(I,1,IMINUS)
            XX(I,3,IPLUS)=XX(I,3,IMINUS)
            DO K=1,3
              IGOT(I,K)=K
            END DO
          END IF
          IF(DABS(XX(I,3,INIT)).LT.1.D-6)THEN
            IF(DABS(XX(I,3,INIT-1)).GT.1.D-6)THEN
C**PLANAR COORDINATE IS 'Z'
              IPLAN=3
              IPLAT(I)=I
            END IF
            XX(I,3,IPLUS)=-XX(I,3,IMINUS)
            XX(I,1,IPLUS)=XX(I,1,IMINUS)
            XX(I,2,IPLUS)=XX(I,2,IMINUS)
            DO K=1,3
              IGOT(I,K)=K
            END DO
          END IF
        END DO
C**COULD BE 'CIS'
        IF(IPLAN.NE.0)THEN
          DO I=1,NATOM
            IF(DABS(XX(I,IPLAN,INIT)).LT.1.D-6.AND.
     1      DABS(XX(I,IPLAN,INIT-1)).LT.1.D-6)IPLAT(I)=I
          END DO
        END IF
C**NOW CHECK NON-PLANAR COORDINATES
        DO I=1,NATOM
          DO J=1,NATOM
            IF(I.EQ.J)GO TO 9994
            DO K=1,3
              IF(IGOT(I,K).NE.0.OR.IGOT(J,K).NE.0)GO TO 9992
              IF((DABS(XX(I,K,INIT)-XX(J,K,INIT)).LT.1.D-6).AND.
     1        (XX(I,K,IMINUS).NE.XX(J,K,IMINUS)))THEN
                XX(I,K,IPLUS)=XX(J,K,IMINUS)
                XX(J,K,IPLUS)=XX(I,K,IMINUS)
                IGOT(I,K)=K
                IGOT(J,K)=K
              END IF
9992          CONTINUE
            END DO
            DO K=1,3
              IF(IGOT(I,K).NE.0.OR.IGOT(J,K).NE.0)GO TO 9993
              IF((DABS(XX(I,K,INIT)+XX(J,K,INIT)).LT.1.D-6).AND.
     1        (XX(I,K,IMINUS).NE.-XX(J,K,IMINUS)))THEN
                XX(I,K,IPLUS)=-XX(J,K,IMINUS)
                XX(J,K,IPLUS)=-XX(I,K,IMINUS)
                IGOT(I,K)=K
                IGOT(J,K)=K
              END IF
9993          CONTINUE
            END DO
9994        CONTINUE
          END DO
        END DO
      END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)'PLANAR COORDINATE AT CIS: ',IPLAN
        WRITE(IOUT,*)'PLANAR ATOMS'
        DO I=1,NATOM
          WRITE(IOUT,*)IPLAT(I)
        END DO
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C*******************************************************************
C*******************************************************************
      XMTOT=0
      DO I=1,NATOM
C**SUM MASSES
        XMTOT=XMTOT+XM(I)
        DO K=1,3
C**RESET TRUE EQUILIBRIUM
          X0(I,K)=XX(I,K,ITAUMN)
        END DO
      END DO
C*******************************************
C**WHEN HAPPY FORM MASS-WEIGHTED COORDINATES
C*******************************************
      DO ITAU=1,722
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
C**RELY ON REFLECTED GEOMETRIES TO GET DERIVATIVES RIGHT

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'**********************************'
          WRITE(IOUT,*)'DERIVATIVES AT POINT ITAU = ',ITAU
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**CALCULATE GRADIENTS OF COORDINATES AT TAU = ITAU
C**************************************************
        DO I=1,NATOM
          DO K=1,3
            XXP(I,K,ITAU)=(XX(I,K,ITAU+1)-XX(I,K,ITAU-1))*RAD
            IF(ITAU.EQ.2)
     1      XXP(I,K,ITAU)=(XX(I,K,ITAU+1)-XX(I,K,ITAU))*2*RAD
            IF(ITAU.EQ.722)
     1      XXP(I,K,ITAU)=(XX(I,K,ITAU)-XX(I,K,ITAU-1))*2*RAD
          END DO
        END DO
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
          WRITE(IOUT,*)'XNORM = ',Z
          WRITE(IOUT,*)
          WRITE(IOUT,*)'GEOMETRIES AT POINT ITAU = ',ITAU
          DO I=1,NATOM
            WRITE(IOUT,*)(XX(I,K,ITAU)/SQRT(XM(I)),K=1,3)
          END DO
          WRITE(IOUT,*)
          WRITE(IOUT,*)'DERIVATIVES AT POINT ITAU = ',ITAU
          DO I=1,NATOM
            WRITE(IOUT,*)(XXP(I,K,ITAU),K=1,3)
          END DO
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

8888  CONTINUE
      END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************END COORDINATE DERIVATIVE ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C*******************************************************************
C**PRINT EQUILIBRIUM VALUES NORMAL MODES FROM FORT.1 FILE

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)
        WRITE(IOUT,*)'EQUILIBRIUM 3N-6 NORMAL MODES IN DATA FILE'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      J0=0
      DO J=1,NMODE
        J0=J0+1
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)'OMEGA ',J,': ',OMEGA(J0,1)*WAVENM
          DO I=1,NATOM
            WRITE(IOUT,*)(XL(I,J0,K,1),K=1,3)
          END DO
        END IF
      END DO
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

      DO ITAU=1,722
C**ONLY GO AS FAR AS 'CIS'
        IF(ITAU.GT.INIT)GO TO 7777
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'**********************************'
          WRITE(IOUT,*)'NORMAL COORDS.  AT POINT ITAU = ',ITAU
        END IF
        IF(ITAU.EQ.1)GO TO 7776
C**SET UP MOMENT OF INERTIA MATRIX
C*********************************
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
          IF(JPRINT.GT.0)WRITE(IOUT,*)'ERROR IN MATINV'
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

7776    CONTINUE
C**FORM THE PROJECTED FORCE CONSTANTS AND 3N-7 NORMAL COORDINATES
C****************************************************************
        CALL PROJEC(NATOM,NMODE,XM,XX(1,1,ITAU),OMEGA(1,ITAU),
     1  XL(1,1,1,ITAU),XR,RR,XP(1,1,ITAU),TEMP,DUMMY,ITAU)

C**FIRST TIME NORMALS RUN...NO PROJECTION
        IF(ITAU.EQ.1)THEN

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)THEN
            WRITE(IOUT,*)
            WRITE(IOUT,*)
            WRITE(IOUT,*)'REFERENCE 3N-6 NORMAL MODES'
            WRITE(IOUT,*)
            J0=0
            DO J=1,NMODE
              J0=J0+1
              WRITE(IOUT,*)'OMEGA ',J,': ',OMEGA(J0,1)*WAVENM
              DO I=1,NATOM
                WRITE(IOUT,*)(XL(I,J0,K,1),K=1,3)
              END DO
            END DO
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

          GO TO 7777
        END IF

C**ON RETURN (NMODE-1) VALUES OF OMEGA AND XL FOR ITAU=2,722
        J0=0
        DO J=1,NMODE-1
          J0=J0+1
C**MODE FOR REACTION PATH WILL BE MISSING, SO SKIP
          IF(J0.EQ.IREACT)J0=J0+1
          IF(ITAU.GT.2)J0=J

C****************************************************************
C**SET PLANAR VECTOR TO ZERO
          IF(ITAU.EQ.JNIT)THEN
            DO I=1,NATOM
              XL(I,J,JNITCD,ITAU)=0
            END DO
C**IN CASE IT IS SWITCHED
            IF(J.NE.NMODE-1)THEN
              DO I=1,NATOM
                XL(I,J+1,JNITCD,ITAU)=0
              END DO
            END IF
          END IF
          IF(ITAU.EQ.INIT)THEN
            DO I=1,NATOM
              XL(I,J,INITCD,ITAU)=0
            END DO
C**IN CASE IT IS SWITCHED
            IF(J.NE.NMODE-1)THEN
              DO I=1,NATOM
                XL(I,J+1,INITCD,ITAU)=0
              END DO
            END IF
          END IF
C****************************************************************

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
          IF(ITAU.EQ.2)GO TO 5555
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
            IF(DABS(SUM).LT.0.99D0)WRITE(IOUT,*)
     1      'MODES ',J,J+1,'  SWITCHED ? - SUM',SUM
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**ARBITRARY GUESS OF 1.D-03 DEFINES ORTHOGONAL VECTORS
          IF(DABS(SUM).LT.1.D-03)THEN
C**VECTORS HAVE SWITCHED
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
C**STANDARDISE WHOLE VECTOR
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

C**CHECK ALL SIGNS RELATIVE TO PREVIOUS VECTOR
          ICLEAR=0
          DO I=1,NATOM
            DO K=1,3
              IF(XL(I,J,K,ITAU).LT.0.AND.XL(I,J0,K,ITAU-1).GE.0)ICLEAR=1
              IF(XL(I,J,K,ITAU).GE.0.AND.XL(I,J0,K,ITAU-1).LT.0)ICLEAR=1
            END DO
          END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)THEN
            IF(ICLEAR.NE.0)THEN
              WRITE(IOUT,*)
              WRITE(IOUT,*)'FOR ITAU = ',ITAU
              WRITE(IOUT,*)'FOLLOWING VECTOR LOOKS WRONG FOR MODE ',J
              WRITE(IOUT,*)'MODE ',J,' OMEGA = ',OMEGA(J,ITAU)*WAVENM
              DO I=1,NATOM
                WRITE(IOUT,*)(XL(I,J,K,ITAU),K=1,3)
              END DO
              WRITE(IOUT,*)'PREVIOUS VECTOR'
              DO I=1,NATOM
                WRITE(IOUT,*)(XL(I,J0,K,ITAU-1),K=1,3)
              END DO
            END IF
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

5555      CONTINUE
        END DO

C**SET MINIMUM OMEGAS AND VECTORS FOR SELECTION OF BASIS
        IF(ITAU.EQ.ITAUMN)THEN
          J0=0
          DO J=1,NMODE-1
            J0=J0+1
            IF(J0.EQ.IREACT)J0=J0+1
C**SET MINIMUM OMEGAS FOR SELECTION OF BASIS
C**USE MINIMUM-ENERGY VECTORS FOR BASIS
            OMEGA(J0,1)=DABS(OMEGA(J,ITAU))
            DO I=1,NATOM
              DO K=1,3
                XL(I,J0,K,1)=XL(I,J,K,ITAU)
              END DO
            END DO
          END DO
        END IF

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'FINAL NORMAL COORDINATES'
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        DO J=1,NMODE-1

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)THEN
            WRITE(IOUT,*)'MODE = ',J,' OMEGA = ',OMEGA(J,ITAU)*WAVENM
            DO I=1,NATOM
              WRITE(IOUT,*)(XL(I,J,K,ITAU),K=1,3)
            END DO
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        END DO

7777    CONTINUE
      END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************END PROJECTION ALGORITHM'
        WRITE(IOUT,*)
        WRITE(IOUT,*)
        WRITE(IOUT,*)'GEOMETRY TO DEFINE PRIMITIVE BASIS'
        WRITE(IOUT,*)
        DO I=1,NATOM
          WRITE(IOUT,265)I,(X0(I,K),K=1,3)
        END DO
        WRITE(IOUT,*)
        WRITE(IOUT,*)'OMEGAs TO DEFINE PRIMITIVE BASIS'
        WRITE(IOUT,*)
        J0=0
        DO J=1,NMODE-1
          J0=J0+1
          IF(J0.EQ.IREACT)J0=J0+1
          WRITE(IOUT,*)'FOR MODE ',J0,' OMEGA = ',OMEGA(J0,1)*
     1    WAVENM
        END DO
        WRITE(IOUT,*)
        WRITE(IOUT,*)'VECTORS TO DEFINE PRIMITIVE BASIS'
        WRITE(IOUT,*)
        J0=0
        DO J=1,NMODE-1
          J0=J0+1
          IF(J0.EQ.IREACT)J0=J0+1
          WRITE(IOUT,*)'FOR MODE ',J0
          DO I=1,NATOM
            WRITE(IOUT,*)(XL(I,J0,K,1),K=1,3)
          END DO
        END DO
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**************************************
C**************************************
C**SET VECTORS CORRECTLY ABOUT 'CIS'
C**************************************
C**************************************
      IPLUS=INIT
      IMINUS=INIT
      DO I=1,360
        IPLUS=IPLUS+1
        IF(IPLUS.GT.722)IPLUS=3
        IMINUS=IMINUS-1
        IF(IMINUS.LT.2)IMINUS=721
        DO L=1,NMODE-1
C**FOR THIS MODE, DETERMINE IF IN-PLANE AT 'CIS'
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
C**NO OUT-OF-PLANE (3N-7) VECTORS FOR HOOH
CTMP          IF(MINVEC.EQ.INITCD)THEN
C**IN-PLANE VECTOR
                IF(K.EQ.INITCD)XL(J,L,K,IPLUS)=-XL(J,L,K,IMINUS)
CTMP          ELSE
C**OUT-OF-PLANE VECTOR
CTMP            IF(K.NE.INITCD)XL(J,L,K,IPLUS)=-XL(J,L,K,IMINUS)
CTMP          END IF
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
        WRITE(IOUT,*)'**************START VECTOR DERIVATIVES ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      DO ITAU=2,722

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'**********************************'
          WRITE(IOUT,*)'VECTOR DERIVATIVES  AT POINT ITAU = ',ITAU
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**CALCULATE GRADIENT OF VECTOR AT TAU = ITAU
C********************************************
C**FORM REACTION COORDINATE (S) AT ITAU+1, ITAU-1
        X=0
        DO I=1,NATOM
          IF(ITAU.NE.722)THEN
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
          IF(ITAU.NE.2)THEN
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
        DSTAU(ITAU)=(SPLUS(ITAU)+SMINUS(ITAU))*RAD

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C**WRITE ARC-LENGTH/RADIAN JACOBIAN JUST FOR FANCY!
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'DSTAU = ',DSTAU(ITAU)
          WRITE(IOUT,*)
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C********************************************
C********************************************
C**NON-ECKART PATH IN RADIANS (1/2 DEGREES)
        SPLUS(ITAU)=0.5D0/RAD
        SMINUS(ITAU)=0.5D0/RAD
C**THEREFORE JACOBIAN = 1
        DSTAU(ITAU)=1.0D0
C********************************************
C********************************************
C**UN-NORMALISE XXP
        DO I=1,NATOM
          DO K=1,3
            XXP(I,K,ITAU)=XXP(I,K,ITAU)*XNORM(ITAU)
          END DO
        END DO
C********************************************
C********************************************
        DO J=1,NMODE-1

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
          IF(JPRINT.GT.0)THEN
            WRITE(IOUT,*)
            WRITE(IOUT,*)'MODE ',J
            WRITE(IOUT,*)
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

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
          SUM1=SUM1*XNORM(ITAU)
          SUM2=SUM2*XNORM(ITAU)
          SUM3=SUM3*XNORM(ITAU)

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
          IF(JPRINT.GT.0)THEN
            WRITE(IOUT,*)'VECTOR PRODUCT XX*XLP ',SUM1P,SUM2P,SUM3P
            IF(DABS(SUM1)-DABS(SUM1P).GT.1.D-3)WRITE(IOUT,*)'SUSPECT'
            IF(DABS(SUM2)-DABS(SUM2P).GT.1.D-3)WRITE(IOUT,*)'SUSPECT'
            IF(DABS(SUM3)-DABS(SUM3P).GT.1.D-3)WRITE(IOUT,*)'SUSPECT'
          END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        END DO
6666    CONTINUE
      END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'**************END VECTOR DERIVATIVES ALGORITHM'
        WRITE(IOUT,*)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C**FORM WORKING QUANTITIES
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

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                  IF(JPRINT.GT.0)THEN
                    WRITE(IOUT,*)'THE FOLLOWING SHOULD ADD TO ZERO'
                    WRITE(IOUT,*)ITAU,I,J,BBS(I,J,ITAU),BBS(J,I,ITAU)
                  END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

                END IF
              END IF
            END DO
          END DO
        END IF
C**NCH TEST

      END DO
C*******************************************************************
C*******************************************************************
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETSAY(NATOM,XM,XX,XR,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XM(NATOM),XX(NATOM,3)
      DIMENSION XR(NATOM,3),RR(NATOM,NATOM)
      COMMON/PLMNMX/MXCD(3)
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
200   FORMAT(//,1X,'DISPLACED INTERNAL COORDINATES',/)
220   FORMAT(//,1X,'DISTORTED INTERNAL COORDINATES',/)
240   FORMAT(1X,'ATOM ',I2,':  X0 =',F12.7,'  Y0 =',F12.7,'  Z0 =',
     1F12.7)
245   FORMAT(//,60(1H*),/,1X,'DISPLACED GEOMETRY',2X,I4/)
255   FORMAT(//,1X,'ALPHA,BETA,GAMMA',/)
260   FORMAT(//,1X,'PREVIOUS GEOMETRY',/)
C*****************************************************
C**XR CONTAINS DISTORTED GEOMETRY
      CALL CHECKM(XM,XR,XX,RR,NATOM)
C**XR CONTAINS DISTORTED GEOMETRY ALONG P.A.
C*****************************************************

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)WRITE(IOUT,200)
      DO I=1,NATOM
        DO K=1,3
          XX(I,K)=XR(I,K)
        END DO
        IF(JPRINT.GT.0)WRITE(IOUT,240)I,(XX(I,J),J=1,3)
      END DO
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE NORXYZ(XX,NATOM,ITAU,IND,IABC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3,722)
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT
C**NORMALISE PREVIOUS 'COORD'
      IND=0
      Z=0
      DO I=1,NATOM
        Z=Z+XX(I,IABC,ITAU)*XX(I,IABC,ITAU)
      END DO
      Z=SQRT(Z)
C**IS CURRENT 'COORD' NORMALISED TO PREVIOUS?
      ZZ=0
      DO I=1,NATOM
        ZZ=ZZ+XX(I,IABC,ITAU)*XX(I,IABC,ITAU+1)
      END DO
      ZZ=SQRT(DABS(ZZ))/Z

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)WRITE(IOUT,*)'X-NORM = ',ZZ
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IF(DABS(ZZ-1.D0).GT.1.D-02)IND=1
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE NOATOM(XX,NATOM,ITAU,IPLAT,IABC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3,722),IPLAT(NATOM)
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT
C**NORMALISE PREVIOUS 'ATOM'
      IPLAT(IABC)=0
      Z=0
      DO I=1,3
        Z=Z+XX(IABC,I,ITAU)*XX(IABC,I,ITAU)
      END DO
      Z=SQRT(Z)
C**IS CURRENT 'ATOM' NORMALISED TO PREVIOUS?
      ZZ=0
      DO I=1,3
        ZZ=ZZ+XX(IABC,I,ITAU)*XX(IABC,I,ITAU+1)
      END DO
      ZZ=SQRT(DABS(ZZ))/Z

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)WRITE(IOUT,*)'X-NORM = ',ZZ
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IF(DABS(ZZ-1.D0).GT.1.D-02)IPLAT(IABC)=1
      RETURN
      END
