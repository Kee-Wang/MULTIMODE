C**************************************************************
C**************************************************************
      SUBROUTINE USERIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FILASS/IOUT,INP
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/DUMP/JDUMP,IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/ROTPAR/RX1,RX2,RX3,THX1,THX2,TAUX
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(200),NTERM,
     1KTERM(200,6),ITERM(200,6),JTERM(200,6),
     2NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3),CDIP(100,3)
200   FORMAT(//,1X,'EQUILIBRIUM GEOMETRY',/)
205   FORMAT(1X,'R1 = ',F10.6,'  R2 = ',F10.6,'  R3 = ',F10.6,/)
210   FORMAT(1X,'TH1 = ',F10.6,'  TH2 = ',F10.6,'  PHI = ',F10.6,/)
215   FORMAT(/,1X,'POTENTIAL TERMS: ITERM, JTERM, CTERM',/)
220   FORMAT(1X,6I2,1X,6I2,1X,D20.12)
225   FORMAT(/,1X,'A1-COMPONENT',/)
230   FORMAT(/,1X,'B2-COMPONENT',/)
235   FORMAT(/,1X,'B1-COMPONENT',/)
240   FORMAT(/,1X,'DIPOLE TERMS: ITERM, JTERM, CTERM',/)
      WRITE(IOUT,200)
C**REFERENCE GEOMETRY IN ANGSTROMS + DEGREES
      READ(INP,*)
      READ(INP,*)RE1,RE2,RE3
      WRITE(IOUT,205)RE1,RE2,RE3
C**CONVERT ANGSTROMS TO BOHR
      RE1=RE1/BOHR
      RE2=RE2/BOHR
      RE3=RE3/BOHR
      RX1=RE1
      RX2=RE2
      RX3=RE3
      READ(INP,*)
      READ(INP,*)THE1,THE2,PHIE
      WRITE(IOUT,210)THE1,THE2,PHIE
C**CONVERT DEGREES TO RADIANS
      THE1=THE1/RAD
      THE2=THE2/RAD
      PHIE=PHIE/RAD
      THX1=THE1
      THX2=THE2
      WRITE(IOUT,200)
      WRITE(IOUT,205)RE1,RE2,RE3
      WRITE(IOUT,210)THE1,THE2,PHIE
      READ(INP,*)
      READ(INP,*)NTERM
      DO I=1,NTERM
        READ(INP,*)(KTERM(I,J),J=1,6),CTERM(I)
      END DO
      DO J=1,6
        DO I=1,NTERM
          ITERM(I,J)=0
          JTERM(I,J)=0
        END DO
      END DO
      DO I=1,NTERM
        K=1
        DO J=1,6
          IF(KTERM(I,J).NE.0)THEN
            ITERM(I,K)=J
            JTERM(I,K)=KTERM(I,J)
            K=K+1
          END IF
        END DO
      END DO
      WRITE(IOUT,215)
      DO I=1,NTERM
        WRITE(IOUT,220)(ITERM(I,J),J=1,6),(JTERM(I,J),J=1,6),CTERM(I)
      END DO
      IF(LDUMP.EQ.0)RETURN
C*************************
      DO L=1,3
        IF(L.EQ.1)WRITE(IOUT,225)
        IF(L.EQ.2)WRITE(IOUT,230)
        IF(L.EQ.3)WRITE(IOUT,235)
        READ(INP,*)
        READ(INP,*)NDIP(L)
        DO I=1,NDIP(L)
          READ(INP,*)(KDIP(I,J,L),J=1,6),CDIP(I,L)
        END DO
        DO J=1,6
          DO I=1,NDIP(L)
            IDIP(I,J,L)=0
            JDIP(I,J,L)=0
          END DO
        END DO
        DO I=1,NDIP(L)
          K=1
          DO J=1,6
            IF(KDIP(I,J,L).NE.0)THEN
              IDIP(I,K,L)=J
              JDIP(I,K,L)=KDIP(I,J,L)
              K=K+1
            END IF
          END DO
        END DO
        WRITE(IOUT,240)
        DO I=1,NDIP(L)
          WRITE(IOUT,220)(IDIP(I,J,L),J=1,6),(JDIP(I,J,L),J=1,6),
     1    CDIP(I,L)
        END DO
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETPOT(V,NATOM,XX,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM)
      DIMENSION QQ(6)
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(200),NTERM,
     1KTERM(200,6),ITERM(200,6),JTERM(200,6),
     2NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3),CDIP(100,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/ECKCNT/ICNT,INTC
      COMMON/FILASS/IOUT,INP
200   FORMAT(1X,6F12.6)
      CALL BONDS(NATOM,RR,XX)
      R1=RR(1,3)
      R2=RR(2,4)
      R3=RR(3,4)
      R4=RR(1,4)
      R5=RR(2,3)
      R6=RR(1,2)
      TH1=DACOS((R1*R1+R3*R3-R4*R4)/(2*R1*R3))
      TH2=DACOS((R2*R2+R3*R3-R5*R5)/(2*R2*R3))
      X1=R1*DSIN(TH1)
      X2=R2*DSIN(TH2)
      X=R3-R1*DCOS(TH1)-R2*DCOS(TH2)
      DISC=R6*R6-X*X
      IF(DISC.LE.0.D0)THEN
        R=0.D0
      ELSE
        R=DSQRT(R6*R6-X*X)
      END IF
      COSPHI=(X1*X1+X2*X2-R*R)/(2*X1*X2)
      IF(COSPHI.GE.1.D0)COSPHI=1.D0
      IF(COSPHI.LE.-1.D0)COSPHI=-1.D0
      PHI=DACOS(COSPHI)
      IF(INTC.NE.0)THEN
        IF(INTC.GT.0)THEN
          RE1=R1
          RE2=R2
          RE3=R3
          THE1=TH1
          THE2=TH2
          PHIE=PHI
        ELSE
          R1=R1-RE1
          R2=R2-RE2
          R3=R3-RE3
          TH1=TH1-THE1
          TH2=TH2-THE2
          PHI=PHI-PHIE
        END IF
        WRITE(IOUT,200)R1,R2,R3,TH1*RAD,TH2*RAD,PHI*RAD
        RETURN
      END IF
C**POTENTIAL USES BOHR + RADIANS
      QQ(1)=(R1-RE1)/R1
      QQ(2)=(R2-RE2)/R2
      QQ(3)=(R3-RE3)/R3
      QQ(4)=TH1-THE1
      QQ(5)=TH2-THE2
      QQ(6)=PHI
      V=0.D0
      DO I=1,NTERM
        TERM=1
C**MAXIMUM OF 6 MODES COUPLED
        DO J=1,6
          IF(ITERM(I,J).NE.0)THEN
            K=ITERM(I,J)
            L=JTERM(I,J)
            FACT=QQ(K)**L
            IF(K.EQ.6)FACT=DCOS(L*PHI)
            TERM=TERM*FACT
          END IF
        END DO
        V=V+CTERM(I)*TERM
      END DO
      V=V+876.831891D0/WAVENM
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETQPT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETAPT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDIP
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE GETQDT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE MINPOT(TAU,NATOM,XR,RR,TAUE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XR(NATOM,3),RR(NATOM,NATOM)
      DIMENSION WRK(200),SOL(5),F(5)
      COMMON/ROTPAR/RX1,RX2,RX3,THX1,THX2,TAUX
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(200),NTERM,
     1KTERM(200,6),ITERM(200,6),JTERM(200,6),
     2NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3),CDIP(100,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
C**RETURNS WITH CARTESIAN COORDINATES IN XR CORRESPONDING TO
C**MINIMUM ENERGY AT FIXED TAU
C*****************************************************
      EXTERNAL DERIV,MONIT
C*****************************************************
      IFAIL=1
      MFIT=5
C**TEMPORARY
C     STEP=10.D0
      STEP=5.D0
      XTOL=1.D-8
      FTOL=1.D-15
      MAXCAL=1000
      KPRINT=1000
C**TEMPORARY
C     WRITE(IOUT,*)'*******************************'
      IWXY=2*MFIT*(MFIT+MFIT)+2*MFIT+5*MFIT
      IF(IWXY.GT.200)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'WRK TOO SMALL'
        STOP 'WRK TOO SMALL IN MINPOT'
      END IF
C**SET CURRENT TAU
      TAUX=TAU+PHIE
      TAUE=TAUX
C**ANSWER NEAR LAST STRUCTURE, SO INITIAL GUESS SET TO IT
      SOL(1)=RX1
      SOL(2)=RX2
      SOL(3)=RX3
      SOL(4)=THX1
      SOL(5)=THX2
      CALL E04FBF(MFIT,MFIT,SOL,F,SUMSQ,FTOL,XTOL,STEP,WRK,200,
     1DERIV,MONIT,KPRINT,MAXCAL,IFAIL)
C     WRITE(IOUT,*)
C     WRITE(IOUT,*)'MINPOT IFAIL = ',IFAIL
      IF(IFAIL.NE.0.AND.IFAIL.NE.3)THEN
        STOP 'ERROR IN E04FBF'
      END IF
C**SAVE THIS STRUCTURE FOR NEXT TIME
      RX1=SOL(1)
      RX2=SOL(2)
      RX3=SOL(3)
      THX1=SOL(4)
      THX2=SOL(5)
C**TEMPORARY
C     WRITE(IOUT,*)'TAU: ',TAUX*RAD
C     WRITE(IOUT,*)'GEOM: ',RX1*BOHR,RX2*BOHR,RX3*BOHR,
C    1THX1*RAD,THX2*RAD
C**TEMPORARY
C**SET CARTESIAN COORDINATES FOR THIS STRUCTURE
      CT=DCOS(TAUX/2)
      IF(DABS(CT).LT.1.D-10)CT=0.D0
      ST=DSIN(TAUX/2)
      IF(DABS(ST).LT.1.D-10)ST=0.D0
      XR(1,1)=RX1*DSIN(THX1)*CT
      XR(1,2)=RX1*DSIN(THX1)*ST
      XR(1,3)=-RX3/2+RX1*DCOS(THX1)
      XR(2,1)=RX2*DSIN(THX2)*CT
      XR(2,2)=-RX2*DSIN(THX2)*ST
      XR(2,3)=RX3/2-RX2*DCOS(THX2)
      XR(3,1)=0.D0
      XR(3,2)=0.D0
      XR(3,3)=-RX3/2
      XR(4,1)=0.D0
      XR(4,2)=0.D0
      XR(4,3)=RX3/2
      CALL GETPOT(V,NATOM,XR,RR)
C     WRITE(IOUT,*)'V = ',V*WAVENM
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DERIV(M,N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),F(M)
      DIMENSION QQ(6),TERM(6)
      COMMON/ROTPAR/RX1,RX2,RX3,THX1,THX2,TAUX
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(200),NTERM,
     1KTERM(200,6),ITERM(200,6),JTERM(200,6),
     2NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3),CDIP(100,3)
      COMMON/FILASS/IOUT,INP
C**RESET INTERNAL COORDINATES
      R1=X(1)
      R2=X(2)
      R3=X(3)
      TH1=X(4)
      TH2=X(5)
      QQ(1)=(R1-RE1)/R1
      QQ(2)=(R2-RE2)/R2
      QQ(3)=(R3-RE3)/R3
      QQ(4)=TH1-THE1
      QQ(5)=TH2-THE2
      QQ(6)=TAUX
C**FIRST DERIVATIVES
      DO I=1,5
        F(I)=0.D0
      END DO
      DO I=1,NTERM
        DO ID=1,5
C**GET DERIVATIVE WRT COORDINATE ID
          FACT=0.D0
          DO K=1,6
            TERM(K)=1.D0
          END DO
C**MAXIMUM OF 6 MODES COUPLED
          DO J=1,6
            IF(ITERM(I,J).NE.0)THEN
              K=ITERM(I,J)
              L=JTERM(I,J)
              IF(K.EQ.ID)THEN
                TERM(K)=CTERM(I)
                FACT=L*QQ(K)**(L-1)
                IF(K.EQ.1)FACT=FACT*RE1/(R1*R1)
                IF(K.EQ.2)FACT=FACT*RE2/(R2*R2)
                IF(K.EQ.3)FACT=FACT*RE3/(R3*R3)
              ELSE
                TERM(K)=QQ(K)**L
                IF(K.EQ.6)TERM(K)=DCOS(L*TAUX)
              END IF
            END IF
          END DO
          Z=1.D0
          DO K=1,6
            Z=Z*TERM(K)
          END DO
          F(ID)=F(ID)+Z*FACT
        END DO
      END DO
CCCC  WRITE(IOUT,*)(F(I),I=1,5)
      RETURN
      END

