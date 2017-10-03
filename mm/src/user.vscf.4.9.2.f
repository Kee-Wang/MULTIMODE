C***********************************************************
C***********************************************************
C**USER-SUPPLIED
C***********************************************************
C***********************************************************
      SUBROUTINE ROTATE(NATOM,X0,OMEGA,NMODE,XL,WAVENM,NXMODE,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     DIMENSION X0(NATOM,3),OMEGA(NMODE),XL(NATOM,NMODE,3)
      DIMENSION X0(NATOM,3),OMEGA(NMODE),XL(NATOM,NXMODE,3)
      COMMON/FILASS/IOUT
C**TEMPORARY
      RETURN
C**TEMPORARY
      PI=DACOS(-1.D0)
      PI2=PI/2
      PI4=PI/4
      GO TO(1,2,3),IND
1     CONTINUE
C******************************************************************
C**ROUTINE TO UNSCRAMBLE E AND T SYMMETRY NORMAL COORDINATE VECTORS
C**SPECIFIC TO CH4
C******************************************************************
CCCC  IF(DABS(X0(1,1))-DABS(X0(1,2)).LT.1.D-6)THEN
C**FIND E SYMMETRIES
        DO J=1,NXMODE
          IF(J.LT.NXMODE)THEN
            IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.1.D-1)THEN
              IF(DABS(OMEGA(J)-OMEGA(J-1))*WAVENM.LT.1.D-1)GO TO 100
              IF(J.LT.NXMODE-1)THEN
                IF(DABS(OMEGA(J)-OMEGA(J+2))*WAVENM.LT.1.D-1)GO TO 100
              END IF
C             CALL ROTATV(NATOM,NXMODE,XL,J,1)
C**ROTATION ANGLE ZERO
              IF(DABS(XL(2,J,1)+XL(3,J,3)).LT.1.D-5)GO TO 100
C**ROTATION ANGLE INFINITE
              IF(DABS(XL(2,J+1,1)+XL(3,J+1,3)).LT.1.D-5)GO TO 100
              TH=(XL(2,J,1)+XL(3,J,3))/
     1           (-(XL(2,J+1,1)+XL(3,J+1,3)))
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
        RETURN
2     CONTINUE
C**********************************************
C**FIND T SYMMETRIES
        DO J=1,NXMODE
          IF(J.LT.NXMODE-1)THEN
            IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.1.D-1)THEN
              IF(DABS(OMEGA(J)-OMEGA(J+2))*WAVENM.LT.1.D-1)THEN
                CALL ROTATV(NATOM,NXMODE,XL,J,3)
              END IF
            END IF
          END IF
        END DO
CCCC  END IF
      RETURN
3     CONTINUE
C******************************************************************
C**ROUTINE TO UNSCRAMBLE RPH INTERSECTING NORMAL COORDINATE VECTORS
C**SPECIFIC TO H5O2+ OR H3O2- OR H5+
C******************************************************************
      DO J=1,NXMODE
        IF(J.LT.NXMODE)THEN
          IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.5.0D0)THEN
C**ROTATION ANGLE ZERO
C**H5O2+
CC          IF(DABS(XL(5,J,1)-XL(6,J,1)).LT.1.D-5)GO TO 300
C**H3O2-
CC          IF(DABS(XL(3,J,1)-XL(5,J,1)).LT.1.D-5)GO TO 300
C**H5+
            IF((DABS(XL(3,J,2)-XL(5,J,2)).LT.1.D-5).AND.
     1      (DABS(XL(3,J,2)-XL(4,J,2)).LT.1.D-5))GO TO 300
C**ROTATION ANGLE INFINITE
C**H5O2+
CC          IF(DABS(XL(5,J+1,1)-XL(6,J+1,1)).LT.1.D-5)GO TO 300
CC          TH=(XL(5,J,1)-XL(6,J,1))/
CC   1         (-XL(5,J+1,1)+XL(6,J+1,1))
C**H3O2-
CC          IF(DABS(XL(3,J+1,1)-XL(5,J+1,1)).LT.1.D-5)GO TO 300
CC          TH=(XL(3,J,1)-XL(5,J,1))/
CC   1         (-XL(3,J+1,1)+XL(5,J+1,1))
C**H5+
            IF(DABS(XL(3,J+1,2)-XL(5,J+1,2)).LT.1.D-5)GO TO 301
            TH=(XL(3,J,2)-XL(5,J,2))/
     1         (-XL(3,J+1,2)+XL(5,J+1,2))
            GO TO 302
301         CONTINUE
            IF(DABS(XL(3,J+1,2)-XL(4,J+1,2)).LT.1.D-5)GO TO 300
            TH=(XL(3,J,2)-XL(4,J,2))/
     1         (-XL(3,J+1,2)+XL(4,J+1,2))
302         CONTINUE
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
300         CONTINUE
          END IF
        END IF
      END DO
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE ROTATV(NATOM,NMODE,XL,J,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XL(NATOM,NMODE,3)
      DIMENSION F(15),SOL(3),WRK(1000)
      COMMON/FILASS/IOUT
      COMMON/VECTR/XLCURR(5,9,3),XLWORK(5,9,3),
     1ITHIS,ITHAT,INEXT
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/PRINT/IPRINT,JPRINT
C*****************************************************
      EXTERNAL VDERV1,VDERV2,MONIT
C*****************************************************
      ITHIS=J
      ITHAT=J+1
      INEXT=J+2
C**TEMPORARY
C       IF(IND.EQ.1)
C    1  WRITE(IOUT,*)'ITHIS = ',ITHIS
        IF(IND.EQ.1)
     1  WRITE(IOUT,*)'ITHIS,ITHAT = ',ITHIS,ITHAT
        IF(IND.EQ.3)
     1  WRITE(IOUT,*)'ITHIS,ITHAT,INEXT = ',ITHIS,ITHAT,INEXT
        CALL FLUSH(IOUT)
C**TEMPORARY
      DO K=1,3
        DO I=1,NATOM
          XLCURR(I,ITHIS,K)=XL(I,ITHIS,K)
          XLCURR(I,ITHAT,K)=XL(I,ITHAT,K)
          IF(IND.GT.2)THEN
            XLCURR(I,INEXT,K)=XL(I,INEXT,K)
          END IF
        END DO
      END DO
C**
      DO MODE=1,INEXT
        WRITE(IOUT,*)'MODE = ',MODE
        DO I=1,NATOM
          WRITE(IOUT,100)(XL(I,MODE,K),K=1,3)
        END DO
      END DO
      CALL FLUSH(IOUT)
100   FORMAT(1X,3(F12.8,1X))
C**
      IFAIL=1
      MFIT=IND*3
      IF(IND.EQ.1)MFIT=2*MFIT
      NFIT=IND
      STEP=1/RAD
      XTOL=1.D-8
      FTOL=1.D-15
      MAXCAL=2000
      KPRINT=2000
      IWXY=2*NFIT*(NFIT+MFIT)+2*MFIT+5*NFIT
      IF(IWXY.GT.1000)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'WRK TOO SMALL'
        STOP 'WRK TOO SMALL IN ROTATV'
      END IF
C**INITIAL ANGLES
      SOL(1)=0
      SOL(2)=0
      SOL(3)=0
      CALL E04FBF(MFIT,NFIT,SOL,F,SUMSQ,FTOL,XTOL,STEP,WRK,1000,
     1VDERV1,MONIT,KPRINT,MAXCAL,IFAIL)
C     IF(JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'ROTATV IFAIL = ',IFAIL
        WRITE(IOUT,*)
C     END IF
      IF(IFAIL.NE.0.AND.IFAIL.NE.3)THEN
        STOP 'ERROR IN E04FBF'
      END IF
C**RELOAD ROTATED VECTORS
      DO K=1,3
        DO I=1,NATOM
          XL(I,ITHIS,K)=XLWORK(I,ITHIS,K)
          XL(I,ITHAT,K)=XLWORK(I,ITHAT,K)
          IF(IND.GT.2)THEN
            XL(I,INEXT,K)=XLWORK(I,INEXT,K)
          END IF
        END DO
      END DO
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE VDERV1(M,N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),F(M)
      DIMENSION XLPREV(5,9,3)
      COMMON/VECTR/XLCURR(5,9,3),XL(5,9,3),ITHIS,ITHAT,INEXT
      COMMON/FILASS/IOUT
      IF(N.EQ.2)GO TO 1002
      IF(N.EQ.3)GO TO 1003
      TH=X(1)
C**FORM TWO NEW VECTORS FROM "ITHIS", "ITHAT"
      CTH=DCOS(TH)
      STH=DSIN(TH)
      DO I=1,5
        DO K=1,3
          XL(I,ITHIS,K)=CTH*XLCURR(I,ITHIS,K)+STH*XLCURR(I,ITHAT,K)
          XL(I,ITHAT,K)=-STH*XLCURR(I,ITHIS,K)+CTH*XLCURR(I,ITHAT,K)
        END DO
      END DO
C**EACH NEW VECTOR SHOULD HAVE C2(X),C2(Y),C2(Z)
C**   C2(X)
CC    SUM1=0
      DO K=1,2
CC      SUM1=SUM1+DABS(XL(1,ITHIS+K-1,2))-DABS(XL(2,ITHIS+K-1,2))
        F(K)=DABS(XL(1,ITHIS+K-1,2))-DABS(XL(2,ITHIS+K-1,2))
      END DO
C**   C2(Y)
CC    SUM2=0
      DO K=1,2
CC      SUM2=SUM2+DABS(XL(1,ITHIS+K-1,3))-DABS(XL(3,ITHIS+K-1,3))
        F(K+2)=DABS(XL(1,ITHIS+K-1,3))-DABS(XL(3,ITHIS+K-1,3))
      END DO
C**   C2(Z)
CC    SUM3=0
      DO K=1,2
CC      SUM3=SUM3+DABS(XL(1,ITHIS+K-1,1))-DABS(XL(4,ITHIS+K-1,1))
        F(K+4)=DABS(XL(1,ITHIS+K-1,1))-DABS(XL(4,ITHIS+K-1,1))
      END DO
      WRITE(IOUT,*)'F6: ',F(1),F(2),F(3),F(4),F(5),F(6)
      RETURN

1002  CONTINUE
      TH=X(1)
      PH=X(2)
C**FORM THREE NEW VECTORS FROM "ITHIS", "ITHAT", "INEXT"
      CTH=DCOS(TH)
      STH=DSIN(TH)
      CPH=DCOS(PH)
      SPH=DSIN(PH)
      DO I=1,9
        DO K=1,3
          XL(I,INEXT,K)=CTH*XLCURR(I,INEXT,K)+
     1    STH*SPH*XLCURR(I,ITHIS,K)-
     2    STH*CPH*XLCURR(I,ITHAT,K)
          XL(I,ITHIS,K)=
     1    CPH*XLCURR(I,ITHIS,K)+
     2    SPH*XLCURR(I,ITHAT,K)
          XL(I,ITHAT,K)=STH*XLCURR(I,INEXT,K)-
     1    CTH*SPH*XLCURR(I,ITHIS,K)+
     2    CTH*CPH*XLCURR(I,ITHAT,K)
        END DO
      END DO
C**DOT PRODUCT BETWEEN "ITHIS" VECTORS
      SUM1=0
      DO I=1,9
        DO K=1,3
          SUM1=SUM1+XL(I,ITHIS,K)*XLPREV(I,ITHIS,K)
        END DO
      END DO
C**DOT PRODUCT BETWEEN "ITHAT" VECTORS
      SUM2=0
      DO I=1,9
        DO K=1,3
          SUM2=SUM2+XL(I,ITHAT,K)*XLPREV(I,ITHAT,K)
        END DO
      END DO
C**SHOULD ALL BE 1
      F(1)=1-DABS(SUM1)
      F(2)=1-DABS(SUM2)
      WRITE(IOUT,*)'F2: ',F(1),F(2)
      RETURN

1003  CONTINUE
      TH=X(1)
      PH=X(2)
      XH=X(3)
C**FORM THREE NEW VECTORS FROM "ITHIS", "ITHAT", "INEXT"
      CTH=DCOS(TH)
      STH=DSIN(TH)
      CPH=DCOS(PH)
      SPH=DSIN(PH)
      CXH=DCOS(XH)
      SXH=DSIN(XH)
      DO I=1,5
        DO K=1,3
          XL(I,INEXT,K)=(CTH*CXH-STH*CPH*SXH)*XLCURR(I,INEXT,K)+
     1    SPH*SXH*XLCURR(I,ITHIS,K)-
     2    (STH*CXH+CTH*CPH*SXH)*XLCURR(I,ITHAT,K)
          XL(I,ITHIS,K)=STH*SPH*XLCURR(I,INEXT,K)+
     1    CPH*XLCURR(I,ITHIS,K)+
     2    CTH*SPH*XLCURR(I,ITHAT,K)
          XL(I,ITHAT,K)=(CTH*SXH+STH*CPH*CXH)*XLCURR(I,INEXT,K)-
     1    SPH*CXH*XLCURR(I,ITHIS,K)-
     2    (STH*SXH-CTH*CPH*CXH)*XLCURR(I,ITHAT,K)
        END DO
      END DO
C**EACH NEW VECTOR SHOULD HAVE C2(X),C2(Y),C2(Z)
C**   C2(X)
CC    SUM1=0
      DO K=1,3
CC      SUM1=SUM1+DABS(XL(1,ITHIS+K-1,2))-DABS(XL(2,ITHIS+K-1,2))
        F(K)=DABS(XL(1,ITHIS+K-1,2))-DABS(XL(2,ITHIS+K-1,2))
      END DO
C**   C2(Y)
CC    SUM2=0
      DO K=1,3
CC      SUM2=SUM2+DABS(XL(1,ITHIS+K-1,3))-DABS(XL(3,ITHIS+K-1,3))
        F(K+3)=DABS(XL(1,ITHIS+K-1,3))-DABS(XL(3,ITHIS+K-1,3))
      END DO
C**   C2(Z)
CC    SUM3=0
      DO K=1,3
CC      SUM3=SUM3+DABS(XL(1,ITHIS+K-1,1))-DABS(XL(4,ITHIS+K-1,1))
        F(K+6)=DABS(XL(1,ITHIS+K-1,1))-DABS(XL(4,ITHIS+K-1,1))
      END DO
C**RELOAD ROTATED VECTORS
C     DO K=1,3
C       DO I=1,NATOM
C         XLCURR(I,ITHIS,K)=XL(I,ITHIS,K)
C         XLCURR(I,ITHAT,K)=XL(I,ITHAT,K)
C         IF(IND.GT.1)THEN
C           XLCURR(I,INEXT,K)=XL(I,INEXT,K)
C         END IF
C       END DO
C     END DO
CC    F(1)=SUM1
CC    F(2)=SUM2
CC    F(3)=SUM3
CC    WRITE(IOUT,*)'F3: ',F(1),F(2),F(3)
      WRITE(IOUT,*)'F9: ',F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8),F(9)
      CALL FLUSH(IOUT)
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE VDERV2(M,N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),F(M)
      DIMENSION XLPREV(5,9,3)
      COMMON/VECTR/XLCURR(5,9,3),XL(5,9,3),ITHIS,ITHAT,INEXT
      COMMON/FILASS/IOUT
      IF(N.EQ.2)GO TO 1002
      IF(N.EQ.3)GO TO 1003
      TH=X(1)
C**FORM TWO NEW VECTORS FROM "ITHIS", "ITHAT"
      CTH=DCOS(TH)
      STH=DSIN(TH)
      DO I=1,5
        DO K=1,3
          XL(I,ITHIS,K)=CTH*XLCURR(I,ITHIS,K)+STH*XLCURR(I,ITHAT,K)
          XL(I,ITHAT,K)=-STH*XLCURR(I,ITHIS,K)+CTH*XLCURR(I,ITHAT,K)
        END DO
      END DO
C**DOT PRODUCT BETWEEN "ITHIS" VECTORS
      SUM1=0
      DO I=1,5
        DO K=1,3
          SUM1=SUM1+XL(I,ITHIS,K)*XLPREV(I,ITHIS,K)
        END DO
      END DO
C**DOT PRODUCT BETWEEN "ITHAT" VECTORS
      SUM2=0
      DO I=1,5
        DO K=1,3
          SUM2=SUM2+XL(I,ITHAT,K)*XLPREV(I,ITHAT,K)
        END DO
      END DO
C**TOTAL
      SUM=DABS(SUM1)+DABS(SUM2)
C**SHOULD BE 2
C     F(1)=2-SUM
      F(1)=1-SUM1
CCCC  WRITE(IOUT,*)'F1: ',F(1)
      RETURN

1002  CONTINUE
      TH=X(1)
      PH=X(2)
C**FORM THREE NEW VECTORS FROM "ITHIS", "ITHAT", "INEXT"
      CTH=DCOS(TH)
      STH=DSIN(TH)
      CPH=DCOS(PH)
      SPH=DSIN(PH)
      DO I=1,5
        DO K=1,3
          XL(I,INEXT,K)=CTH*XLCURR(I,INEXT,K)+
     1    STH*SPH*XLCURR(I,ITHIS,K)-
     2    STH*CPH*XLCURR(I,ITHAT,K)
          XL(I,ITHIS,K)=
     1    CPH*XLCURR(I,ITHIS,K)+
     2    SPH*XLCURR(I,ITHAT,K)
          XL(I,ITHAT,K)=STH*XLCURR(I,INEXT,K)-
     1    CTH*SPH*XLCURR(I,ITHIS,K)+
     2    CTH*CPH*XLCURR(I,ITHAT,K)
        END DO
      END DO
C**DOT PRODUCT BETWEEN "ITHIS" VECTORS
      SUM1=0
      DO I=1,5
        DO K=1,3
          SUM1=SUM1+XL(I,ITHIS,K)*XLPREV(I,ITHIS,K)
        END DO
      END DO
C**DOT PRODUCT BETWEEN "ITHAT" VECTORS
      SUM2=0
      DO I=1,5
        DO K=1,3
          SUM2=SUM2+XL(I,ITHAT,K)*XLPREV(I,ITHAT,K)
        END DO
      END DO
C**SHOULD ALL BE 1
      F(1)=1-DABS(SUM1)
      F(2)=1-DABS(SUM2)
      WRITE(IOUT,*)'F2: ',F(1),F(2)
      RETURN

1003  CONTINUE
      TH=X(1)
      PH=X(2)
      XH=X(3)
C**FORM THREE NEW VECTORS FROM "ITHIS", "ITHAT", "INEXT"
      CTH=DCOS(TH)
      STH=DSIN(TH)
      CPH=DCOS(PH)
      SPH=DSIN(PH)
      CXH=DCOS(XH)
      SXH=DSIN(XH)
      DO I=1,5
        DO K=1,3
          XL(I,INEXT,K)=(CTH*CXH-STH*CPH*SXH)*XLCURR(I,INEXT,K)+
     1    SPH*SXH*XLCURR(I,ITHIS,K)-
     2    (STH*CXH+CTH*CPH*SXH)*XLCURR(I,ITHAT,K)
          XL(I,ITHIS,K)=STH*SPH*XLCURR(I,INEXT,K)+
     1    CPH*XLCURR(I,ITHIS,K)+
     2    CTH*SPH*XLCURR(I,ITHAT,K)
          XL(I,ITHAT,K)=(CTH*SXH+STH*CPH*CXH)*XLCURR(I,INEXT,K)-
     1    SPH*CXH*XLCURR(I,ITHIS,K)-
     2    (STH*SXH-CTH*CPH*CXH)*XLCURR(I,ITHAT,K)
        END DO
      END DO
C**EACH NEW VECTOR SHOULD HAVE C2(X),C2(Y),C2(Z)
C**   C2(X)
CC    SUM1=0
      DO K=1,3
CC      SUM1=SUM1+DABS(XL(1,ITHIS+K-1,2))-DABS(XL(2,ITHIS+K-1,2))
        F(K)=DABS(XL(1,ITHIS+K-1,2))-DABS(XL(2,ITHIS+K-1,2))
      END DO
C**   C2(Y)
CC    SUM2=0
      DO K=1,3
CC      SUM2=SUM2+DABS(XL(1,ITHIS+K-1,3))-DABS(XL(3,ITHIS+K-1,3))
        F(K+3)=DABS(XL(1,ITHIS+K-1,3))-DABS(XL(3,ITHIS+K-1,3))
      END DO
C**   C2(Z)
CC    SUM3=0
      DO K=1,3
CC      SUM3=SUM3+DABS(XL(1,ITHIS+K-1,1))-DABS(XL(4,ITHIS+K-1,1))
        F(K+6)=DABS(XL(1,ITHIS+K-1,1))-DABS(XL(4,ITHIS+K-1,1))
      END DO
      F(10)=DABS(XL(1,ITHIS,1))-DABS(XL(1,ITHAT,1))
      F(11)=DABS(XL(1,ITHIS,1))-DABS(XL(1,INEXT,1))
      F(12)=DABS(XL(1,ITHIS,2))-DABS(XL(1,ITHAT,2))
      F(13)=DABS(XL(1,ITHIS,2))-DABS(XL(1,INEXT,2))
      F(14)=DABS(XL(1,ITHIS,3))-DABS(XL(1,ITHAT,3))
      F(15)=DABS(XL(1,ITHIS,3))-DABS(XL(1,INEXT,3))
C**RELOAD ROTATED VECTORS
C     DO K=1,3
C       DO I=1,NATOM
C         XLCURR(I,ITHIS,K)=XL(I,ITHIS,K)
C         XLCURR(I,ITHAT,K)=XL(I,ITHAT,K)
C         IF(IND.GT.1)THEN
C           XLCURR(I,INEXT,K)=XL(I,INEXT,K)
C         END IF
C       END DO
C     END DO
CC    F(1)=SUM1
CC    F(2)=SUM2
CC    F(3)=SUM3
CC    WRITE(IOUT,*)'F3: ',F(1),F(2),F(3)
      WRITE(IOUT,*)'F15: ',F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8),
     1F(9),F(10),F(11),F(12),F(13),F(14),F(15)
      CALL FLUSH(IOUT)
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE RTGEOM(NATOM,X0,RR,E,WAVENM,INDRTG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X0(NATOM,3),E(3),RR(NATOM,NATOM)
      COMMON/MOMI/XK(3,3),XMU(3,3)
C*************************************************
C**ROUTINE TO ROTATE PRINCIPAL AXIS SYSTEM IF E-MODES
C**SPECIFIC TO NH3
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
     1        (2*STH122*STH132)
              IF(CTH12.GE.1)CTH12=1
              IF(CTH12.LE.-1)CTH12=-1
              TH12=DACOS(CTH12)
              TH1=2*TH12
              CTH1=COS(TH1)
              STH1=SIN(TH1)
              CTH32=(STH232*STH232+STH132*STH132-STH122*STH122)/
     1        (2*STH132*STH232)
              IF(CTH32.GE.1)CTH32=1
              IF(CTH32.LE.-1)CTH32=-1
              TH32=DACOS(CTH32)
              TH3=2*TH32
              CTH3=COS(TH3)
              STH3=SIN(TH3)
      CTH3=(CTH1+CTH3)/2
      CTH1=CTH3
      STH3=SQRT(1.D0-CTH3*CTH3)
      STH1=STH3
              CBESQ=(SIN(TH32)*SIN(TH32)-STH122*STH122)/
     1        (SIN(TH32)*SIN(TH32))
              IF(CBESQ.LT.0)CBESQ=0.D0
              COSBET=SQRT(CBESQ)
      COSBET=0.D0
              IF(COSBET.GE.1)COSBET=1
              IF(COSBET.LE.-1)COSBET=-1
              BET=DACOS(COSBET)
              SINBET=SIN(BET)
      SINBET=1.D0
              RX1=R1*SINBET
      RX1=(R1+R3)/2
              RX2=R2*SINBET
              RX3=R3*SINBET
      RX3=RX1
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
