C****************************************************************
C****************************************************************
C**REACTION PATH (ANALYTICAL POTENTIAL....ECKART)
C****************************************************************
C****************************************************************
      SUBROUTINE REACT(NATOM,NMODE,XM,X0,OMEGA,XL,XLP,XX,XXP,RR,XR,XP,
     1TEMP,BB,BBS,B,BS,BSS,IREACT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PLANAR(3),PLNMOL,ECKART
C***************************************************************
C***************************************************************
C**REACTION PATH (ECKART VARIATION) TAYLORED TO HYDROGENPEROXIDE
C**ALLOW ONE POINT FOR EQUILIBRIUM (1)
C**ALLOWS FOR ONE POINT EVERY HALF DEGREE (721)
C**ALLOW EXTRA POINT FOR EQUILIBRIUM 2ND TIME (722)
C**
C**WORK IN ARC-LENGTH
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
      DIMENSION SUP1(3,3),SUP2(3,3),SUP3(10)
      DIMENSION XNORM(722)
      COMMON/SCOORD/DSTAU(722),SPLUS(722),SMINUS(722)
      COMMON/MOMI/XK(3,3),XMU(3,3)
      COMMON/CHECK/MCHECK
      COMMON/ROTIND/ANGSAV(3),INDROT
      COMMON/PLANE/PLNMOL,PLANAR
      COMMON/PLMNMX/MXCD(3)
      COMMON/ECKCON/ECKART
      COMMON/ECKIND/INDECK
      COMMON/REACTN/IDUMMY,MMTAU,INIT,INCTAU
      COMMON/RETURN/IRET
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
220   FORMAT(/,1X,'MODE',I4,'  OMEGA = ',F20.12,/)
245   FORMAT(//,1X,'DISPLACEMENTS OF NORMAL MODES, L(alpha i,k)')
250   FORMAT(/,1X,'MODE(k) = ',I4)
255   FORMAT(1X,'ATOM(i) ',I2,':  x =',F12.6,'  y =',F12.6,'  z =',
     1F12.6)
260   FORMAT(/,1X,'ECKART GEOMETRY',/)
265   FORMAT(1X,'ATOM ',I2,':  X0 =',F12.7,'  Y0 =',F12.7,'  Z0 =',
     1F12.7)
270   FORMAT(/,1X,'PREVIOUS GEOMETRY',/)
C*******************************************************************
      ECKART=.TRUE.

      DO I=1,3
        MXCD(I)=0
      END DO
      IRET=-1
C*******************************************************************
C*******************************************************************
C**MAKE SURE TRUE EQUILIBRIUM ORIENTATED CORRECTLY IF REQUIRED
C**X0 CONTAINS EQUILIBRIUM GEOMETRY WHICH SATISFIES EVERYTHING
      IF(MCHECK.NE.0)CALL CHECKM(XM,X0,XX,RR,NATOM)
C**XX(1) CONTAINS EQUILIBRIUM
      DO I=1,NATOM
        DO K=1,3
          XX(I,K,1)=X0(I,K)
        END DO
      END DO

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
C**GET INCREMENT TO NEXT INTEGRAL TAU
      INCR=TAUE*RAD+1
      TAUINC=INCR-TAUE*RAD
C**JNIT WILL POINT AT 'TRANS'
      JNIT=0
C**INIT WILL POINT AT 'CIS'
      INIT=0
C**INDECK SET WHEN PLANAR FOR NEXT POINT
      INDECK=0
C**INDROT SET AFTER FIRST CALL TO GETSAY
      INDROT=0

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

C**FIRST LOOP FROM 0 TO 360 DEGREES, INCREMENTS OF 1/2 DEGREE
      VMIN=1.D+20
      DO ITAU=1,721
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

C**SET TEMPORARY EQUILIBRIUM
        DO I=1,NATOM
          DO K=1,3
            X0(I,K)=XX(I,K,ITAU)
          END DO
        END DO
        IF(ITAU.EQ.1)THEN
          TAU=TAU+TAUINC/RAD
        ELSE
          TAU=TAU+0.5D0/RAD
        END IF
C**POINT TO 180 DEGREES 
        IF(JNIT.EQ.0)THEN
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
C**GET (UNROTATED) CARTESIAN COORDINATES FOR NEXT VALUE OF TAU IN XR
        CALL MINPOT(TAU,NATOM,XR,RR,TAUX)
        CALL GETPOT(V,NATOM,XR,RR)
C**FIND ITAU CORRESPONDING TO MINIMUM ENERGY
        IF(V.LT.VMIN)THEN
          VMIN=V
          ITAUMN=ITAU+1
        END IF
C**NOW MOVE TO CENTRE OF MASS AND ROTATE
        CALL GETSAY(NATOM,XX(1,1,ITAU-1),XM,XX(1,1,ITAU+1),XR,RR)
C**GET PHASE RIGHT
        IF(PLNMOL.OR.INDECK.EQ.0)THEN
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
        END IF

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C**PREVIOUS INTERNAL DISPLACEMENT COORDINATES
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,270)
          DO I=1,NATOM
            WRITE(IOUT,265)I,(XX(I,K,ITAU),K=1,3)
          END DO
C**ECKART INTERNAL DISPLACEMENT COORDINATES
          WRITE(IOUT,260)
          DO I=1,NATOM
            WRITE(IOUT,265)I,(XX(I,K,ITAU+1),K=1,3)
          END DO
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        SUM1=0
        DO I=1,NATOM
          SUM1=SUM1+(XX(I,2,ITAU)*XX(I,3,ITAU+1)-XX(I,3,ITAU)*
     1    XX(I,2,ITAU+1))*XM(I)
        END DO
        SUM2=0
        DO I=1,NATOM
          SUM2=SUM2+(XX(I,3,ITAU)*XX(I,1,ITAU+1)-XX(I,1,ITAU)*
     1    XX(I,3,ITAU+1))*XM(I)
        END DO
        SUM3=0
        DO I=1,NATOM
          SUM3=SUM3+(XX(I,1,ITAU)*XX(I,2,ITAU+1)-XX(I,2,ITAU)*
     1    XX(I,1,ITAU+1))*XM(I)
        END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)WRITE(IOUT,*)'ECKART PRODUCT XX*XX ',
     1  SUM1,SUM2,SUM3
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        IF(INDECK.NE.0.AND..NOT.PLNMOL)INDECK=0
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
      DO I=1,359
        IPLUS=IPLUS+1
        IF(IPLUS.GT.722)IPLUS=3
        IMINUS=IMINUS-1
        IF(IMINUS.LT.2)IMINUS=721
        DO J=1,NATOM
          DO K=1,3
            XX(J,K,IPLUS)=XX(J,K,IMINUS)
            IF(K.EQ.INITCD)XX(J,K,IPLUS)=-XX(J,K,IMINUS)
          END DO
        END DO
      END DO
C**GET PHASE ITAU=2 CORRECT (THIS HAS BEEN MOVED FROM
C**TAU TO 4*PI - TAU
      ICH=0
      IX=1
      DO I=2,NATOM
        IF(DABS(XX(I,1,3)).GT.DABS(XX(IX,1,3)))IX=I
      END DO
      IF(XX(IX,1,3).GE.0.AND.XX(IX,1,2).LT.0)ICH=1
      IF(XX(IX,1,3).LE.0.AND.XX(IX,1,2).GT.0)ICH=1
      IF(ICH.NE.0)THEN
        DO I=1,NATOM
          XX(I,1,2)=-XX(I,1,2)
        END DO
      END IF
      ICH=0
      IY=1
      DO I=2,NATOM
        IF(DABS(XX(I,2,3)).GT.DABS(XX(IY,2,3)))IY=I
      END DO
      IF(XX(IY,2,3).GE.0.AND.XX(IY,2,2).LT.0)ICH=1
      IF(XX(IY,2,3).LE.0.AND.XX(IY,2,2).GT.0)ICH=1
      IF(ICH.NE.0)THEN
        DO I=1,NATOM
          XX(I,2,2)=-XX(I,2,2)
        END DO
      END IF
      ICH=0
      IZ=1
      DO I=2,NATOM
        IF(DABS(XX(I,3,3)).GT.DABS(XX(IZ,3,3)))IZ=I
      END DO
      IF(XX(IZ,3,3).GE.0.AND.XX(IZ,3,2).LT.0)ICH=1
      IF(XX(IZ,3,3).LE.0.AND.XX(IZ,3,2).GT.0)ICH=1
      IF(ICH.NE.0)THEN
        DO I=1,NATOM
          XX(I,3,2)=-XX(I,3,2)
        END DO
      END IF
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
C**ONLY GO AS FAR AS 'CIS'
        IF(ITAU.GT.INIT)GO TO 8888

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
            XXP(I,K,ITAU)=(XX(I,K,ITAU+1)-XX(I,K,ITAU-1))
C**PATH CONTINUOUS THROUGH ITAU = ..., 721, 722(2), 3, ...
            IF(ITAU.EQ.2)
     1      XXP(I,K,ITAU)=(XX(I,K,ITAU+1)-XX(I,K,721))
            IF(ITAU.EQ.722)
     1      XXP(I,K,ITAU)=(XX(I,K,3)-XX(I,K,ITAU-1))
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

        SUM1=0
        DO I=1,NATOM
          SUM1=SUM1+XX(I,2,ITAU)*XXP(I,3,ITAU)-XX(I,3,ITAU)*
     1    XXP(I,2,ITAU)
        END DO
        SUM2=0
        DO I=1,NATOM
          SUM2=SUM2+XX(I,3,ITAU)*XXP(I,1,ITAU)-XX(I,1,ITAU)*
     1    XXP(I,3,ITAU)
        END DO
        SUM3=0
        DO I=1,NATOM
          SUM3=SUM3+XX(I,1,ITAU)*XXP(I,2,ITAU)-XX(I,2,ITAU)*
     1    XXP(I,1,ITAU)
        END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)WRITE(IOUT,*)'ECKART PRODUCT XX*XXP ',
     1  SUM1,SUM2,SUM3
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

C**************************************
C**************************************
C**SET DERIVATIVES CORRECTLY ABOUT 'CIS'
C**************************************
C**************************************
      IPLUS=INIT
      IMINUS=INIT
      DO I=1,359
        IPLUS=IPLUS+1
        IF(IPLUS.GT.722)IPLUS=3
        IMINUS=IMINUS-1
        IF(IMINUS.LT.2)IMINUS=721
        DO J=1,NATOM
          DO K=1,3
            XXP(J,K,IPLUS)=XXP(J,K,IMINUS)
            IF(K.NE.INITCD)XXP(J,K,IPLUS)=-XXP(J,K,IMINUS)
          END DO
        END DO
      END DO

C**RESET PLANAR GEOMETRIES AND DERIVATIVES (SPECIAL ACTION)
      DO J=1,NATOM
        DO K=1,3
C**TEMPORARY STORAGE UNTIL LATER FOR SAFETY IN PROJECTION SCHEME
          BB(1,J,K,JNIT)=XX(J,K,JNIT)
          BB(2,J,K,JNIT)=XXP(J,K,JNIT)
          XX(J,K,JNIT)=XX(J,K,JNIT+1)
          XXP(J,K,JNIT)=XXP(J,K,JNIT+1)
          B(J,1,K,INIT)=XX(J,K,INIT)
          B(J,2,K,INIT)=XXP(J,K,INIT)
          XX(J,K,INIT)=XX(J,K,INIT-1)
          XXP(J,K,INIT)=XXP(J,K,INIT-1)
        END DO
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

      write(*,*) omega(1:6,1)
      DO ITAU=2,722
C**ONLY GO AS FAR AS 'CIS'
        IF(ITAU.GT.INIT)GO TO 7777

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IF(JPRINT.GT.0)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)'**********************************'
          WRITE(IOUT,*)'NORMAL COORDS.  AT POINT ITAU = ',ITAU
        END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

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
        DO I=1,3
          DO J=1,3
            XMU(I,J)=SUP2(I,J)
          END DO
        END DO
C**FORM PROJECTION MATRIX
C************************
        ILHS=0
        DO I=1,NATOM
C**KI LABELS GAMMA (LHS)
          DO KI=1,3
            ILHS=ILHS+1
            IRHS=0
            DO J=1,NATOM
              IF(KI.EQ.1)THEN
                SUP3(1)=XX(I,3,ITAU)*
     2          (XMU(2,2)*XX(J,3,ITAU)-XMU(2,3)*XX(J,2,ITAU))
     3                 -XX(I,2,ITAU)*
     4          (XMU(3,2)*XX(J,3,ITAU)-XMU(3,3)*XX(J,2,ITAU))
                SUP3(2)=XX(I,3,ITAU)*
     2          (XMU(2,3)*XX(J,1,ITAU)-XMU(2,1)*XX(J,3,ITAU))
     3                 -XX(I,2,ITAU)*
     4          (XMU(3,3)*XX(J,1,ITAU)-XMU(3,1)*XX(J,3,ITAU))
                SUP3(3)=XX(I,3,ITAU)*
     2          (XMU(2,1)*XX(J,2,ITAU)-XMU(2,2)*XX(J,1,ITAU))
     3                 -XX(I,2,ITAU)*
     4          (XMU(3,1)*XX(J,2,ITAU)-XMU(3,2)*XX(J,1,ITAU))
              END IF
              IF(KI.EQ.2)THEN
                SUP3(1)=XX(I,1,ITAU)*
     2          (XMU(3,2)*XX(J,3,ITAU)-XMU(3,3)*XX(J,2,ITAU))
     3                 -XX(I,3,ITAU)*
     4          (XMU(1,2)*XX(J,3,ITAU)-XMU(1,3)*XX(J,2,ITAU))
                SUP3(2)=XX(I,1,ITAU)*
     2          (XMU(3,3)*XX(J,1,ITAU)-XMU(3,1)*XX(J,3,ITAU))
     3                 -XX(I,3,ITAU)*
     4          (XMU(1,3)*XX(J,1,ITAU)-XMU(1,1)*XX(J,3,ITAU))
                SUP3(3)=XX(I,1,ITAU)*
     2          (XMU(3,1)*XX(J,2,ITAU)-XMU(3,2)*XX(J,1,ITAU))
     3                 -XX(I,3,ITAU)*
     4          (XMU(1,1)*XX(J,2,ITAU)-XMU(1,2)*XX(J,1,ITAU))
              END IF
              IF(KI.EQ.3)THEN
                SUP3(1)=XX(I,2,ITAU)*
     2          (XMU(1,2)*XX(J,3,ITAU)-XMU(1,3)*XX(J,2,ITAU))
     3                 -XX(I,1,ITAU)*
     4          (XMU(2,2)*XX(J,3,ITAU)-XMU(2,3)*XX(J,2,ITAU))
                SUP3(2)=XX(I,2,ITAU)*
     2          (XMU(1,3)*XX(J,1,ITAU)-XMU(1,1)*XX(J,3,ITAU))
     3                 -XX(I,1,ITAU)*
     4          (XMU(2,3)*XX(J,1,ITAU)-XMU(2,1)*XX(J,3,ITAU))
                SUP3(3)=XX(I,2,ITAU)*
     2          (XMU(1,1)*XX(J,2,ITAU)-XMU(1,2)*XX(J,1,ITAU))
     3                 -XX(I,1,ITAU)*
     4          (XMU(2,1)*XX(J,2,ITAU)-XMU(2,2)*XX(J,1,ITAU))
              END IF
C**KJ LABELS GAMMA (RHS)
              DO KJ=1,3
                IRHS=IRHS+1
                XP(ILHS,IRHS,ITAU)=XXP(I,KI,ITAU)*XXP(J,KJ,ITAU)
                IF(KI.EQ.KJ)XP(ILHS,IRHS,ITAU)=XP(ILHS,IRHS,ITAU)+
     1          SQRT(XM(I)*XM(J))/XMTOT
                XP(ILHS,IRHS,ITAU)=XP(ILHS,IRHS,ITAU)+SUP3(KJ)
              END DO
            END DO
          END DO
        END DO
C**FORM THE PROJECTED FORCE CONSTANTS AND 3N-7 NORMAL COORDINATES
C****************************************************************
        CALL PROJEC(NATOM,NMODE,XM,XX(1,1,ITAU),OMEGA(1,ITAU),
     1  XL(1,1,1,ITAU),XR,RR,XP(1,1,ITAU),TEMP,DUMMY,ITAU)
C**ON RETURN (NMODE-1) VALUES OF OMEGA AND XL FOR ITAU=2,722
        J0=0
        DO J=1,NMODE-1
          J0=J0+1
C**MODE FOR REACTION PATH WILL BE MISSING, SO SKIP
          IF(J0.EQ.IREACT)J0=J0+1
          IF(ITAU.GT.2)J0=J

C****************************************************************
C**SET PLANAR VECTOR TO ZERO (SPECIAL ACTION)
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
                WRITE(IOUT,*)(XL(I,J,K,ITAU-1),K=1,3)
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
      write(*,*) omega(1:6,1)

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
C**RESET PLANAR GEOMETRIES AND DERIVATIVES
      DO J=1,NATOM
        DO K=1,3
C**TEMPORARY STORAGE UNTIL LATER
          XX(J,K,JNIT)=BB(1,J,K,JNIT)
          XXP(J,K,JNIT)=BB(2,J,K,JNIT)
          XX(J,K,INIT)=B(J,1,K,INIT)
          XXP(J,K,INIT)=B(J,2,K,INIT)
        END DO
      END DO
      IPLUS=INIT
      IMINUS=INIT
      DO I=1,359
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
C**ONLY GO AS FAR AS 'CIS'
        IF(ITAU.GT.INIT)GO TO 6666

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
          IF(ITAU.NE.722.AND.ITAU.NE.INIT-1.AND.ITAU.NE.JNIT-1)THEN
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
          IF(ITAU.NE.2.AND.ITAU.NE.INIT+1.AND.ITAU.NE.JNIT+1)THEN
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
C**FORM JACOBIAN BETWEEN ARC-LENGTH AND RADIANS
        DSTAU(ITAU)=(SPLUS(ITAU)+SMINUS(ITAU))*RAD
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
     1        XLP(I,J,K,ITAU)=(XL(I,J,K,ITAU+1)-XL(I,J,K,721))
              IF(ITAU.EQ.722)
     1        XLP(I,J,K,ITAU)=(XL(I,J,K,3)-XL(I,J,K,ITAU-1))
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
          IF(JPRINT.GT.0)WRITE(IOUT,*)'VECTOR PRODUCT XX*XL ',
     1    SUM1,SUM2,SUM3
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
          IF(JPRINT.GT.0)WRITE(IOUT,*)'VECTOR PRODUCT XXP*XL ',
     1    SUM1,SUM2,SUM3
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

C**************************************
C**************************************
C**SET VECTOR DERIVATIVES CORRECTLY ABOUT 'CIS'
C**************************************
C**************************************
      IPLUS=INIT
      IMINUS=INIT
      DO I=1,359
        IPLUS=IPLUS+1
        IF(IPLUS.GT.722)IPLUS=3
        IMINUS=IMINUS-1
        IF(IMINUS.LT.2)IMINUS=721
        DSTAU(IPLUS)=DSTAU(IMINUS)
        SPLUS(IPLUS)=SPLUS(IMINUS)
        SMINUS(IPLUS)=SMINUS(IMINUS)
        DO L=1,NMODE-1
C**FOR THIS MODE, DETERMINE IF IN-PLANE AT 'CIS'
          JX=1
          MINVEC=1
          DO J=1,NATOM
            DO K=1,3
              IF(DABS(XL(J,L,K,INIT)).LT.DABS(XL(J,L,MINVEC,INIT)))
     1        THEN
                MINVEC=K
                JX=J
              END IF
            END DO
          END DO
          DO J=1,NATOM
            DO K=1,3
              XLP(J,L,K,IPLUS)=XLP(J,L,K,IMINUS)
C**NO OUT-OF-PLANE (3N-7) VECTORS FOR HOOH
CTMP          IF(MINVEC.EQ.INITCD)THEN
C**IN-PLANE VECTOR
                IF(K.NE.INITCD)XLP(J,L,K,IPLUS)=-XLP(J,L,K,IMINUS)
CTMP          ELSE
C**OUT-OF-PLANE VECTOR
CTMP            IF(K.EQ.INITCD)XLP(J,L,K,IPLUS)=-XLP(J,L,K,IMINUS)
CTMP          END IF
            END DO
          END DO
        END DO
      END DO

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
      SUBROUTINE GETSAY(NATOM,X0,XM,XX,XR,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PLANAR(3),PLNMOL
      DIMENSION XM(NATOM),XX(NATOM,3),X0(NATOM,3)
      DIMENSION XR(NATOM,3),RR(NATOM,NATOM)
      DIMENSION R(3,3)
      DIMENSION WRK(100),SOL(3),F(3)
      COMMON/ROTIND/ANGSAV(3),INDROT
      COMMON/XACC/XACC
      COMMON/PLMNMX/MXCD(3)
      COMMON/FILASS/IOUT,INP
      COMMON/ECKIND/INDECK
      COMMON/PLANE/PLNMOL,PLANAR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/PRINT/IPRINT,JPRINT
C*****************************************************
      EXTERNAL GUNCT,MONIT
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
C**XR CONTAINS DISTORTED GEOMETRY
      CALL CHECKM(XM,XR,XX,RR,NATOM)
C**XX CONTAINS DISTORTED GEOMETRY AT CENTRE OF MASS
C*****************************************************
C**SEE IF PLANAR (X-COORDINATE ZERO)
      PLANAR(1)=.TRUE.
      DO I=1,NATOM
        PLANAR(1)=(PLANAR(1).AND.DABS(XX(I,1)).LT.1.D-6)
      END DO
C**SEE IF PLANAR (Y-COORDINATE ZERO)
      PLANAR(2)=.TRUE.
      DO I=1,NATOM
        PLANAR(2)=(PLANAR(2).AND.DABS(XX(I,2)).LT.1.D-6)
      END DO
C**SEE IF PLANAR (Z-COORDINATE ZERO)
      PLANAR(3)=.TRUE.
      DO I=1,NATOM
        PLANAR(3)=(PLANAR(3).AND.DABS(XX(I,3)).LT.1.D-6)
      END DO
      IF(PLANAR(1).AND.PLANAR(2))STOP 'PLANAR IN X AND Y'
      IF(PLANAR(1).AND.PLANAR(3))STOP 'PLANAR IN X AND Z'
      IF(PLANAR(2).AND.PLANAR(3))STOP 'PLANAR IN Y AND Z'
      PLNMOL=(PLANAR(1).OR.PLANAR(2).OR.PLANAR(3))

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        IF(PLANAR(1))WRITE(IOUT,*)'PLANAR IN X'
        IF(PLANAR(2))WRITE(IOUT,*)'PLANAR IN Y'
        IF(PLANAR(3))WRITE(IOUT,*)'PLANAR IN Z'
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C*********************************************************
C*********************************************************
      IF(INDECK.NE.0)GO TO 9999
C**XR CONTAINS DISTORTED GEOMETRY
      DO I=1,NATOM
        DO K=1,3
          XR(I,K)=XX(I,K)
        END DO
      END DO

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,200)
        DO I=1,NATOM
          WRITE(IOUT,240)I,(XX(I,J),J=1,3)
        END DO
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IFAIL=1
      NFIT=3
      STEP=.1D0
      XTOL=1.D-9
      FTOL=1.D-15
      MFIT=3
      STEP=1.D0
      XACC=FTOL
      MAXCAL=1000
      KPRINT=1000
      IWXY=2*MFIT*(NFIT+MFIT)+2*NFIT+5*MFIT
      IF(IWXY.GT.100)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'WRK TOO SMALL'
        STOP 'WRK TOO SMALL IN GETSAY'
      END IF
      IF(INDROT.EQ.0)THEN
C**ANSWER NEAR ZERO, SO INITIAL GUESS SET TO ZERO
        SOL(1)=0.D0
        SOL(2)=0.D0
        SOL(3)=0.D0
        SOL(1)=90.D0/RAD
        SOL(2)=90.D0/RAD
        SOL(3)=-90.D0/RAD
      ELSE
C**SET INITIAL GUESS TO PREVIOUS RESULT
        SOL(1)=ANGSAV(1)
        SOL(2)=ANGSAV(2)
        SOL(3)=ANGSAV(3)
      END IF
      CALL E04FBF(NFIT,MFIT,SOL,F,SUMSQ,FTOL,XTOL,STEP,WRK,100,
     1GUNCT,MONIT,KPRINT,MAXCAL,IFAIL)
      IF(IFAIL.NE.0.AND.IFAIL.NE.3)THEN
        STOP 'ERROR IN E04FBF'
      END IF
      AL=SOL(1)
      BE=SOL(2)
      GA=SOL(3)

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,255)
        WRITE(IOUT,*)AL*RAD,BE*RAD,GA*RAD
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IF(.NOT.PLNMOL)THEN
        ANGSAV(1)=SOL(1)
        ANGSAV(2)=SOL(2)
        ANGSAV(3)=SOL(3)
      END IF
      SA=SIN(AL)
      CA=COS(AL)
      SB=SIN(BE)
      CB=COS(BE)
      SG=SIN(GA)
      CG=COS(GA)
      R(1,1)=CB*CA*CG-SA*SG
      R(1,2)=CB*SA*CG+CA*SG
      R(1,3)=-SB*CG
      R(2,1)=-CB*CA*SG-SA*CG
      R(2,2)=-CB*SA*SG+CA*CG
      R(2,3)=SB*SG
      R(3,1)=SB*CA
      R(3,2)=SB*SA
      R(3,3)=CB
      DO I=1,NATOM
        DO K=1,3
          XX(I,K)=0
          DO J=1,3
            XX(I,K)=XX(I,K)+R(K,J)*XR(I,J)
          END DO
        END DO
      END DO
      XMAXC=0
      DO J=1,3
        DO I=1,NATOM
          IF(DABS(XX(I,J)).GT.XMAXC)THEN
            XMAXC=DABS(XX(I,J))
            MAXCOD=J
          END IF
        END DO
      END DO
C*********************************************************
C*********************************************************
9999  CONTINUE
      IF(PLNMOL.OR.INDECK.NE.0)THEN
        IF(PLNMOL)THEN
C**SET CORRECT COORD TO ZERO FOR PLANAR GEOMETRY
          K=MXCD(1)
          DO I=1,NATOM
            XX(I,K)=0
          END DO
        END IF
        IF(INDECK.NE.0)THEN
C**SUBSTITUTE PREVIOUS GEOMETRY MINUS 1
          DO I=1,NATOM
            DO K=1,3
              XX(I,K)=X0(I,K)
            END DO
          END DO
C**INTERCHANGE SIGN (PLANAR) COORDINATE
          K=MXCD(1)
          DO I=1,NATOM
            XX(I,K)=-XX(I,K)
          END DO
        ELSE
C**SET INDECK FOR NEXT TIME IF PLANAR THIS TIME
          IF(MXCD(1).NE.0)INDECK=1
        END IF
      END IF
C*********************************************************
C*********************************************************
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

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'MXCD ',(MXCD(I),I=1,3)
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      IF(PLNMOL)THEN
        DO K=1,3
          IF(MXCD(1).NE.0)PLANAR(K)=(MXCD(1).EQ.K)
        END DO
      END IF

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IF(JPRINT.GT.0)THEN
        IF(PLANAR(1))WRITE(IOUT,*)'PLANAR IN X'
        IF(PLANAR(2))WRITE(IOUT,*)'PLANAR IN Y'
        IF(PLANAR(3))WRITE(IOUT,*)'PLANAR IN Z'
      END IF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      INDROT=INDROT+1
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GUNCT(M,N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PLANAR(3),PLNMOL
      DIMENSION X(N),F(M),R(3,3)
      COMMON/MODES/NMODE,NATOM
      COMMON/PLANE/PLNMOL,PLANAR
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/CWORK/W(1)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF
C**************************************************ASSIGN TOTAL STORAGE
      AL=X(1)
      BE=X(2)
      GA=X(3)
      SA=SIN(AL)
      CA=COS(AL)
      SB=SIN(BE)
      CB=COS(BE)
      SG=SIN(GA)
      CG=COS(GA)
      R(1,1)=CB*CA*CG-SA*SG
      R(1,2)=CB*SA*CG+CA*SG
      R(1,3)=-SB*CG
      R(2,1)=-CB*CA*SG-SA*CG
      R(2,2)=-CB*SA*SG+CA*CG
      R(2,3)=SB*SG
      R(3,1)=SB*CA
      R(3,2)=SB*SA
      R(3,3)=CB
      CALL GETZM(NATOM,R,W(LYZ),W(LYZ+3*NATOM),W(LXM),W(LX0),F(1))
      RETURN
      END
