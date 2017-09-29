C***********************************************************************
C subroutine for furan B97-1/TZ2P; Taylor coordinates; with dipoles
C***********************************************************************
      SUBROUTINE USERIN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/STUART/PHI2,PHI3,PHI4
      PARAMETER (NM=21)
C** NM=# of modes, NLx=# of lines in FCMQx_au
      DIMENSION PHI2(NM),PHI3(NM,NM,NM),PHI4(NM,NM,NM)

      COMMON/RUDOLF/DPI0,DPI1,DPI2,M2,DPI3,M3
      PARAMETER (NMD=21, ND2=231, ND3=441, NAXIS=3)
C** NMD=# of modes, NDx=# of terms in DIPQx_A1
C** NAXIS=# of non-zero cartesian components of dipole moment 
      DIMENSION DPI0(NAXIS),DPI1(NMD,NAXIS),DPI2(ND2,NAXIS)
      DIMENSION DPI3(ND3,NAXIS)
      DIMENSION M2(ND2,NMD,NAXIS),M3(ND3,NMD,NAXIS)

      COMMON/DUMP/JDUMP,IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/FILASS/IOUT,INP
C     COMMON/STUART/ALPHA(3),BETA(3),PHI2(3),PHI3(6),N3(6,3),
C    1PHI4(9),N4(9,3)
C     COMMON/RUDOLF/DPI0(3),DPI1(3,3),DPI2(6,3),M2(6,3,3),DPI3(9,3),
C    1M3(9,3,3)

C
C** initialise the arrays
      DO 10, K=1, NM
         PHI2(K)=0.0D0
         DO 10, J=1,NM
            DO 10, I=1,NM
               PHI3(I,J,K)=0.0D0
               PHI4(I,J,K)=0.0D0
 10   CONTINUE

C** read the PHI2(I)=d^2 V/dQ_i^2 from fcmq2_au.
C
      READ(INP,*)N2
      DO 20, L=1, N2
         READ(INP,*) TEMP, I
         PHI2(I)=TEMP
 20   CONTINUE
C

C** read the PHI3(I,J,K)=d^3 V/dQ_i dQ_j dQ_k   from fcmq3_au.
C
      READ(INP,*)N3
      DO 30, L=1, N3
         READ(INP,*) TEMP, I, J, K
         PHI3(I,J,K)=TEMP
 30   CONTINUE
C

C** read the PHI4(I,J,K)=d^4 V/dQ_i dQ_j dQ_k dQ_k   from fcmq4_au.
C
      READ(INP,*)N4
      DO 40, L=1, N4
         READ(INP,*) TEMP, I, J, K
         PHI4(I,J,K)=TEMP
 40   CONTINUE
C
      IF(LDUMP.EQ.0)RETURN

C*****************************************
C**READ DIPOLE MOMENTS (A1, B2, B1)
C*****************************************

      DO IDIP=1,NAXIS

C**READ A1,B2,B1 DIPOLES

C** read the DPI1(I)=d V/dQ_i from DIPQ1
C
        READ(INP,*)
        READ(INP,*)DPI0(IDIP)
C
        READ(INP,*)N1
        DO 50, L=1, N1
           READ(INP,*) DPI1(L,IDIP)
 50     CONTINUE
C

C** read the PHI2(I)=d^2 V/dQ_i dQ_j from DIPQ2
        DO I=1,ND2
          DO J=1,NMD
            M2(I,J,IDIP)=0
          END DO
        END DO
C
        READ(INP,*)N2
        DO 60, L=1, N2
           READ(INP,*) TEMP, I, J
           DPI2(L,IDIP)=TEMP
           M2(L,I,IDIP)=M2(L,I,IDIP)+1
           M2(L,J,IDIP)=M2(L,J,IDIP)+1
 60     CONTINUE
C

C** read the PHI3(I,J,K)=d^3 V/dQ_i dQ_j dQ_k from DIPQ3
        DO I=1,ND3
          DO J=1,NMD
            M3(I,J,IDIP)=0
          END DO
        END DO
C
        READ(INP,*)N3
        DO 70, L=1, N3
CC         READ(13,*) TEMP, I, J, K
           READ(INP,*) TEMP, I, J, K
           DPI3(L,IDIP)=TEMP
           M3(L,I,IDIP)=M3(L,I,IDIP)+1
           M3(L,J,IDIP)=M3(L,J,IDIP)+1
           M3(L,K,IDIP)=M3(L,K,IDIP)+1
 70     CONTINUE
C

      END DO

      RETURN
      END
CC****************************************************************
CC****************************************************************
      SUBROUTINE GETPOT(V,NATOM,XX,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM)
      DIMENSION QQ(30)
      COMMON/POTPAR/RE1,RE2,PHIE,CTERM(6500),NTERM,KTERM(6500,4),
     1ITERM(148,6),JTERM(148,6),NDIP(3),KDIP(100,6,3),
     2IDIP(100,6,3),JDIP(100,6,3)
      COMMON/POTPAR1/CDIP(100,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/ECKCNT/ICNT,INTC
      COMMON/FILASS/IOUT,INP
      COMMON/CART/XE(12),YE(12),ZE(12)
      COMMON/TENL/X(36,12),Y(36,12),Z(36,12) 
200   FORMAT(1X,6F12.6)
C
C**POTENTIAL IN CM-1
      DO I=1,30
      QQ(I)=0.0D0
      DO J=1,12
      QQ(I)=QQ(I)+(XX(J,1)-XE(J))*X(I,J)
      QQ(I)=QQ(I)+(XX(J,2)-YE(J))*Y(I,J)
      QQ(I)=QQ(I)+(XX(J,3)-ZE(J))*Z(I,J)
      END DO
      QQ(I)=QQ(I)*DSQRT(CTERM(I))
      END DO
      V=0
      DO I=1,NTERM
        TERM=1
      IF(KTERM(I,3).EQ.0.AND.KTERM(I,4).EQ.0)FAC=0.5D0
      IF(KTERM(I,3).NE.0.AND.KTERM(I,4).EQ.0)FAC=1.0D0/6.0D0
      IF(KTERM(I,4).NE.0)FAC=1.0D0/24.0D0
        DO J=1,4
          IF(KTERM(I,J).NE.0)THEN
            K=KTERM(I,J)
            TERM=TERM*QQ(K)
          END IF
        END DO
        V=V+FAC*CTERM(I)*TERM
      END DO
      V=V/WAVENM
      RETURN
      END
C**************************************************
C**************************************************
      SUBROUTINE GETQPT(V,NMODE,QQ,XTANPM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),XTANPM(NMODE)

      PARAMETER (HALF = 0.5D0, SIXTH = 1.0D0/6.0D0,
     &A24TH = 1.0D0/24.0D0, QUART = 0.25, THIRD = 1.0D0/3.0D0)

      COMMON/STUART/PHI2,PHI3,PHI4
      PARAMETER (NM=21)
C** NM=# of modes
      DIMENSION PHI2(NM),PHI3(NM,NM,NM),PHI4(NM,NM,NM)


C     COMMON/STUART/ALPHA(3),BETA(3),PHI2(3),PHI3(6),N3(6,3),
C    1PHI4(9),N4(9,3)

      V=0.0D0

C** coupled potential with the correct coefficients
C**   (due to the permutations of the indices).
      DO 30, I=1, NMODE
C         get the diagonal terms
         TEMP=QUART*PHI4(I,I,I)*QQ(I)+PHI3(I,I,I)
         TEMP=THIRD*TEMP*QQ(I)+PHI2(I)
         V=V+HALF*TEMP*QQ(I)*QQ(I)
C         now get the phi3(aab) and phi4(aaab) terms.
         DO 40, J=1, NMODE
            IF(J .EQ. I) GO TO 40
            TEMP=THIRD*PHI4(J,I,I)*QQ(I)+PHI3(J,I,I)
            V=V+HALF*TEMP*QQ(I)*QQ(I)*QQ(J)
 40      CONTINUE
C         now get the phi4(aabb) terms
         DO 50, J=I+1, NMODE
            TEMP=QQ(I)*QQ(J)
            V=V+QUART*PHI4(J,J,I)*TEMP*TEMP
 50      CONTINUE
C         now get phi4(aabc) terms
         DO 60, J=1, NMODE
            IF(J .EQ. I) GO TO 60
            DO 70, K=J+1, NMODE
               IF(K .EQ. I) GO TO 70
               V=V+HALF*PHI4(K,J,I)*QQ(I)*QQ(I)*QQ(J)*QQ(K)
 70         CONTINUE
 60      CONTINUE
C         and finally the phi3(abc) terms
         DO 80, J=I+1, NMODE
            DO 90, K=J+1, NMODE
               V=V+PHI3(K,J,I)*QQ(I)*QQ(J)*QQ(K)
 90         CONTINUE
 80      CONTINUE
 30   CONTINUE

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
C**************************************************
C**************************************************
      SUBROUTINE GETQDT(V,NMODE,QQ,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE)

      COMMON/RUDOLF/DPI0,DPI1,DPI2,M2,DPI3,M3
      PARAMETER (NMD=21, ND2=231, ND3=441, NAXIS=3)
C** NMD=# of modes, NDx=# of lines in DIPQx_A1
C** NAXIS=# of non-zero cartesian components of dipole moment 
      DIMENSION DPI0(NAXIS),DPI1(NMD,NAXIS),DPI2(ND2,NAXIS)
      DIMENSION DPI3(ND3,NAXIS)
      DIMENSION M2(ND2,NMD,NAXIS),M3(ND3,NMD,NAXIS)

C     COMMON/RUDOLF/DPI0(3),DPI1(3,3),DPI2(6,3),M2(6,3,3),DPI3(9,3),
C    1M3(9,3,3)

      DIMENSION Y(NMD)

      DO I=1,NMD
        Y(I)=QQ(I)
      END DO

C**CONSTANT
      V=DPI0(IDIP)

C**LINEAR
      DO I=1,NMD
        V=V+DPI1(I,IDIP)*Y(I)
      END DO

C**QUADRATICS
      DO I=1,ND2
        TERM=1
        DO J=1,NMD
          IF(M2(I,J,IDIP).NE.0)TERM=TERM*Y(J)**M2(I,J,IDIP)
        END DO
        V=V+DPI2(I,IDIP)*TERM
      END DO

C**CUBICS
      DO I=1,ND3
        TERM=1
        DO J=1,NMD
          IF(M3(I,J,IDIP).NE.0)TERM=TERM*Y(J)**M3(I,J,IDIP)
        END DO
        V=V+DPI3(I,IDIP)*TERM
      END DO

      RETURN
      END
C**************************************************
C**************************************************
      SUBROUTINE MINPOT
      RETURN
      END
