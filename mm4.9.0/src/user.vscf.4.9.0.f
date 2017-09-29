C***********************************************************
C***********************************************************
C**USER-SUPPLIED
C***********************************************************
C***********************************************************
      SUBROUTINE ROTATE(NATOM,X0,OMEGA,NMODE,XL,WAVENM,NXMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X0(NATOM,3),OMEGA(NMODE),XL(NATOM,NMODE,3)
      COMMON/FILASS/IOUT
C******************************************************************
C**ROUTINE TO UNSCRAMBLE RPH INTERSECTING NORMAL COORDINATE VECTORS
C**SPECIFIC TO H5O2+
C******************************************************************
        DO J=1,NXMODE
          IF(J.LT.NXMODE)THEN
            IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.5.0D0)THEN
C**ROTATION ANGLE ZERO
              IF(DABS(XL(5,J,1)-XL(6,J,1)).LT.1.D-5)GO TO 100
C**ROTATION ANGLE INFINITE
              IF(DABS(XL(5,J+1,1)-XL(6,J+1,1)).LT.1.D-5)GO TO 100
              TH=(XL(5,J,1)-XL(6,J,1))/
     1           (-XL(5,J+1,1)+XL(6,J+1,1))
              TH=DATAN(TH)
              CTH=DCOS(TH)
              STH=DSIN(TH)
              DO I=1,NATOM
                DO K=1,3
                  TEMP1=CTH*XL(I,J,K)+STH*XL(I,J+1,K)
                  TEMP2=-STH*XL(I,J,K)+CTH*XL(I,J+1,K)
                  XL(I,J,K)=TEMP1
                  XL(I,J+1,K)=TEMP2
                END DO
              END DO
100           CONTINUE
            END IF
          END IF
        END DO
C******************************************************************
C**ROUTINE TO UNSCRAMBLE E AND T SYMMETRY NORMAL COORDINATE VECTORS
C**SPECIFIC TO CH4
C******************************************************************
CCCC  IF(DABS(X0(1,1))-DABS(X0(1,2)).LT.1.D-6)THEN
C**FIND E SYMMETRIES
CC      DO J=1,NXMODE
CC        IF(J.LT.NXMODE)THEN
CC          IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.1.D-3)THEN
C**ROTATION ANGLE ZERO
CC            IF(DABS(XL(2,J,1)+XL(3,J,1)).LT.1.D-5)GO TO 101
C**ROTATION ANGLE INFINITE
CC            IF(DABS(XL(2,J+1,2)+XL(3,J+1,2)).LT.1.D-5)GO TO 101
CC            TH=(XL(2,J,1)+XL(3,J,1))/
CC   1           (-(XL(2,J+1,2)+XL(3,J+1,2)))
CC            TH=DATAN(TH)
CC            CTH=DCOS(TH)
CC            STH=DSIN(TH)
CC            DO I=1,NATOM
CC              DO K=1,3
CC                TEMP1=CTH*XL(I,J,K)+STH*XL(I,J+1,K)
CC                TEMP2=-STH*XL(I,J,K)+CTH*XL(I,J+1,K)
CC                XL(I,J,K)=TEMP1
CC                XL(I,J+1,K)=TEMP2
CC              END DO
CC            END DO
101           CONTINUE
CC          END IF
CC        END IF
CC      END DO
C**********************************************
C**FIND T SYMMETRIES
CCC     DO J=1,NXMODE
CCC       IF(J.LT.NXMODE-1)THEN
CCC         IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.1.D-3)THEN
CCC         IF(DABS(OMEGA(J)-OMEGA(J+2))*WAVENM.LT.1.D-3)THEN
C**ROTATION ANGLE ZERO
CCC           IF(DABS(XL(2,J,3)-XL(3,J,3)).LT.1.D-5)GO TO 200
C**ROTATION ANGLE INFINITE
CCC           IF(DABS(XL(2,J+1,3)-XL(3,J+1,3)).LT.1.D-5)GO TO 200
CCC           TH=(XL(2,J,3)-XL(3,J,3))/
CCC  1           (-XL(2,J+1,3)+XL(3,J+1,3))
CCC           TH=DATAN(TH)
CCC           CTH=DCOS(TH)
CCC           STH=DSIN(TH)
CCC           DO I=1,NATOM
CCC             DO K=1,3
CCC               TEMP1=CTH*XL(I,J,K)+STH*XL(I,J+1,K)
CCC               TEMP2=-STH*XL(I,J,K)+CTH*XL(I,J+1,K)
CCC               XL(I,J,K)=TEMP1
CCC               XL(I,J+1,K)=TEMP2
CCC             END DO
CCC           END DO
200           CONTINUE
CCC         END IF
CCC         END IF
CCC       END IF
CCC     END DO
CCCC  END IF
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE RTGEOM(NATOM,X0,RR,E,WAVENM,INDRTG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X0(NATOM,3),E(3),RR(NATOM,NATOM)
      COMMON/MOMI/XK(3,3),XMU(3,3)
C*************************************************
C**ROUTINE TO ROTATE VECTORS OF PRINCIPAL AXIS SYSTEM
C**SPECIFIC TO CH4
C*************************************************
C**TEMPORARY
      RETURN
C**TEMPORARY
C**FIND E SYMMETRIES
CCC   IF(DABS(X0(4,1))-DABS(X0(4,2)).LT.1.D-6)THEN
        DO J=1,3
          IF(J.LT.3)THEN
            IF(DABS(1/E(J)-1/E(J+1))*WAVENM.LT.1.D-1)THEN
      CALL BONDS(NATOM,RR,X0)
      R1=RR(1,3)
      R2=RR(1,4)
      R3=RR(1,2)
      R4=RR(2,3)
      R5=RR(2,4)
      R6=RR(3,4)
      CTH13=(R1*R1+R3*R3-R4*R4)/(2*R1*R3)
      IF(CTH13.GE.1)CTH13=1
      IF(CTH13.LE.-1)CTH13=-1
      CTH23=(R2*R2+R3*R3-R5*R5)/(2*R2*R3)
      IF(CTH23.GE.1)CTH23=1
      IF(CTH23.LE.-1)CTH23=-1
      CTH12=(R1*R1+R2*R2-R6*R6)/(2*R1*R2)
      IF(CTH12.GE.1)CTH12=1
      IF(CTH12.LE.-1)CTH12=-1
      TH13=DACOS(CTH13)
      TH23=DACOS(CTH23)
      TH12=DACOS(CTH12)
      STH132=SIN(TH13/2)
      STH232=SIN(TH23/2)
      STH122=SIN(TH12/2)
      CTH12=(STH132*STH132+STH122*STH122-STH232*STH232)/
     1(2*STH122*STH132)
      IF(CTH12.GE.1)CTH12=1
      IF(CTH12.LE.-1)CTH12=-1
      TH12=DACOS(CTH12)
      TH1=2*TH12
      CTH1=COS(TH1)
      STH1=SIN(TH1)
      CTH32=(STH232*STH232+STH132*STH132-STH122*STH122)/
     1(2*STH132*STH232)
      IF(CTH32.GE.1)CTH32=1
      IF(CTH32.LE.-1)CTH32=-1
      TH32=DACOS(CTH32)
      TH3=2*TH32
      CTH3=COS(TH3)
      STH3=SIN(TH3)
      CBESQ=(SIN(TH32)*SIN(TH32)-STH122*STH122)/
     1(SIN(TH32)*SIN(TH32))
      IF(CBESQ.LT.0)CBESQ=0.D0
      COSBET=SQRT(CBESQ)
      IF(COSBET.GE.1)COSBET=1
      IF(COSBET.LE.-1)COSBET=-1
      BET=DACOS(COSBET)
      SINBET=SIN(BET)
      RX1=R1*SINBET
      RX2=R2*SINBET
      RX3=R3*SINBET
      Z1=R1*COSBET
      Z2=R2*COSBET
      Z3=R3*COSBET
      DO I=1,3
        X0(1,I)=0
      END DO
      X0(2,1)=RX3*CTH1
      X0(2,2)=RX3*STH1
      X0(2,3)=Z3
      X0(3,1)=RX1*CTH3
      X0(3,2)=-RX1*STH3
      X0(3,3)=Z1
      X0(4,1)=RX2
      X0(4,2)=0
      X0(4,3)=Z2
            INDRTG=1
            END IF
          END IF
        END DO
CCC   END IF
      RETURN
      END
