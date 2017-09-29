C**************************************************************
C**************************************************************
      SUBROUTINE USERIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FILASS/IOUT,INP
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/DUMP/JDUMP,IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(148),NTERM,
     1KTERM(148,6),ITERM(148,6),JTERM(148,6),NDUMMY,
     2CDIP(100,3),NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3)
      COMMON/ALPHA/ALPHA1,ALPHA2,ALPHA3
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
      READ(INP,*)
      READ(INP,*)THE1,THE2,PHIE
      WRITE(IOUT,210)THE1,THE2,PHIE
C**CONVERT DEGREES TO RADIANS
      THE1=THE1/RAD
      THE2=THE2/RAD
      PHIE=PHIE/RAD
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
      READ(INP,*)
      READ(INP,*)ALPHA1,ALPHA2,ALPHA3
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
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(148),NTERM,
     1KTERM(148,6),ITERM(148,6),JTERM(148,6),NDUMMY,
     2CDIP(100,3),NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/ALPHA/ALPHA1,ALPHA2,ALPHA3
      COMMON/ECKCNT/ICNT,INTC
      COMMON/FILASS/IOUT,INP
200   FORMAT(1X,6F12.6)
      CALL BONDS(NATOM,RR,XX)
      R1=RR(1,3)
      R2=RR(1,4)
      R3=RR(1,2)
      R4=RR(2,3)
      R5=RR(2,4)
      R6=RR(3,4)
      TH1=DACOS((R1*R1+R3*R3-R4*R4)/(2*R1*R3))
      TH2=DACOS((R2*R2+R3*R3-R5*R5)/(2*R2*R3))
      ALPHA=DACOS((R1*R1+R2*R2-R6*R6)/(2*R1*R2))
      COSPHI=((COS(ALPHA)-COS(TH1)*COS(TH2))/(SIN(TH1)*SIN(TH2)))
      IF(COSPHI.GE.1)COSPHI=1
      IF(COSPHI.LE.-1)COSPHI=-1
      PHI=DACOS(COSPHI)
C     IF(INTC.NE.0)THEN
C       IF(INTC.GT.0)THEN
C         RE1=R1
C         RE2=R2
C         RE3=R3
C         THE1=TH1
C         THE2=TH2
C         PHIE=PHI
C       ELSE
C         R1=R1-RE1
C         R2=R2-RE2
C         R3=R3-RE3
C         TH1=TH1-THE1
C         TH2=TH2-THE2
C         PHI=PHI-PHIE
C       END IF
C       WRITE(IOUT,200)R1,R2,R3,TH1*RAD,TH2*RAD,PHI*RAD
C       RETURN
C     END IF
C**POTENTIAL USES BOHR + RADIANS
      QQ(1)=(1-EXP(-ALPHA1*(R1-RE1)))
      QQ(2)=(1-EXP(-ALPHA2*(R2-RE2)))
      QQ(3)=(1-EXP(-ALPHA3*(R3-RE3)))
      QQ(4)=TH1-THE1
      QQ(5)=TH2-THE2
      QQ(6)=PHI-PHIE
      V=0
      DO I=1,NTERM
        TERM=1
C**MAXIMUM OF 6 MODES COUPLED
        DO J=1,6
          IF(ITERM(I,J).NE.0)THEN
            K=ITERM(I,J)
            L=JTERM(I,J)
            TERM=TERM*QQ(K)**L
          END IF
        END DO
        V=V+CTERM(I)*TERM
      END DO
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE GETQPT
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE GETAPT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDIP(V,NATOM,XX,RR,LDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM)
      DIMENSION QQ(6)
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(148),NTERM,
     1KTERM(148,6),ITERM(148,6),JTERM(148,6),NDUMMY,
     2CDIP(100,3),NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/ALPHA/ALPHA1,ALPHA2,ALPHA3
      COMMON/FILASS/IOUT,INP
200   FORMAT(1X,6F12.6)
      IF(XX(1,1).GE.0)IPOS=1
      IF(XX(1,1).LT.0)IPOS=-1
      CALL BONDS(NATOM,RR,XX)
      R1=RR(1,3)
      R2=RR(1,4)
      R3=RR(1,2)
      R4=RR(2,3)
      R5=RR(2,4)
      R6=RR(3,4)
      TH1=DACOS((R1*R1+R3*R3-R4*R4)/(2*R1*R3))
      TH2=DACOS((R2*R2+R3*R3-R5*R5)/(2*R2*R3))
      ALPHA=DACOS((R1*R1+R2*R2-R6*R6)/(2*R1*R2))
      COSPHI=((COS(ALPHA)-COS(TH1)*COS(TH2))/(SIN(TH1)*SIN(TH2)))
      IF(COSPHI.GE.1)COSPHI=1
      IF(COSPHI.LE.-1)COSPHI=-1
      PHI=DACOS(COSPHI)
C**DIPOLES USE ANGSTROMS + RADIANS
      QQ(1)=(R1-RE1)*BOHR
      QQ(2)=(R2-RE2)*BOHR
      QQ(3)=(R3-RE3)*BOHR
      QQ(4)=TH1-THE1
      QQ(5)=TH2-THE2
      QQ(6)=(PHI-PHIE)*IPOS
      V=0
      DO I=1,NDIP(LDIP)
        TERM=1
C**MAXIMUM OF 6 MODES COUPLED
        DO J=1,6
          IF(IDIP(I,J,LDIP).NE.0)THEN
            K=IDIP(I,J,LDIP)
            L=JDIP(I,J,LDIP)
            TERM=TERM*QQ(K)**L
          END IF
        END DO
        V=V+CDIP(I,LDIP)*TERM
      END DO
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE GETQDT
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE MINPOT
      RETURN
      END
