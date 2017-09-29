C**************************************************************
C**************************************************************
C**ECKART CONDITIONS
C**************************************************************
C**************************************************************
      SUBROUTINE ECKART(NATOM,NMODE,X0,XL,XM,XX,XR,RR,MBF,XQ,QQ,NNMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XM(NATOM),XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3)
      DIMENSION XR(NATOM,3),RR(NATOM,NATOM),MBF(NMODE),XQ(1),QQ(NMODE)
      COMMON/COUPLE/ICOUPL
      COMMON/FILASS/IOUT,INP
      COMMON/ECKCNT/ICNT,INTC
      COMMON/MOLPRO/MOLPRO,MOUT
200   FORMAT(//,1X,'EQUILIBRIUM BOND DISTANCES',/)
      REWIND MOUT
C**EQUILIBRIUM INTERNAL COORDINATES
      WRITE(IOUT,200)
      CALL BONDS(NATOM,RR,X0)
      DO I=1,NATOM
        WRITE(IOUT,*)(RR(I,J),J=I+1,NATOM)
      END DO
C**DISPLACED CARTESIAN COORDINATES
      K3=0
      DO K=1,NNMODE
        IK2=MBF(K)
        CALL DUMEC1(XQ(1+K3),IK2,NMODE,NATOM,QQ,XM,XR,X0,XL,K,XX,RR)
        IF(ICOUPL.LE.1)GO TO 5001
        L3=0
        DO L=1,K-1
          IL2=MBF(L)
          CALL DUMEC2(XQ(1+K3),XQ(1+L3),IK2,IL2,NMODE,NATOM,
     1    QQ,XM,XR,X0,XL,K,L,XX,RR)
          IF(ICOUPL.EQ.2)GO TO 5002
          N3=0
          DO N=1,L-1
            IN2=MBF(N)
            CALL DUMEC3(XQ(1+K3),XQ(1+L3),XQ(1+N3),IK2,IL2,IN2,
     1      NMODE,NATOM,QQ,XM,XR,X0,XL,K,L,N,XX,RR)
            IF(ICOUPL.EQ.3)GO TO 5003
            M3=0
            DO M=1,N-1
              IM2=MBF(M)
              CALL DUMEC4(XQ(1+K3),XQ(1+L3),XQ(1+N3),XQ(1+M3),IK2,
     1        IL2,IN2,IM2,NMODE,NATOM,QQ,XM,XR,X0,XL,K,L,N,M,XX,RR)
              IF(ICOUPL.EQ.4)GO TO 5004
C*****************************************
5004          CONTINUE
              M3=M3+IM2
            END DO
5003        CONTINUE
            N3=N3+IN2
          END DO
5002      CONTINUE
          L3=L3+IL2
        END DO
5001    CONTINUE
        K3=K3+IK2
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMEC1(XQ,MM,NMODE,NATOM,QQ,XM,XR,X0,XL,MODE,XX,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),SPACE,CHSYM(4)
      DIMENSION XQ(MM),QQ(NMODE),X0(NATOM,3),XM(NATOM)
      DIMENSION XR(NATOM,3),XL(NATOM,NMODE,3)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM)
      COMMON/DUMP/JDUMP,IDUMP,KDUMP,MDUMP
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
204   FORMAT(I4)
205   FORMAT(1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2)
206   FORMAT(A2,1X,3F20.10)
208   FORMAT(2F20.10)
C**ONLY DUMP IF UNIQUE SYMMETRY
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE)MSYM=I
        END DO
      END DO
      ISKIP=0
      IF(MSYM.NE.1)ISKIP=1

      DO K=1,NMODE
        QQ(K)=0
      END DO
      MDUM=MDUMP
      IF(ISKIP.NE.0)MDUM=MDUMP/2
      DO M=1,MDUM
C       N=MM/2+1+((-1)**M)*(MM/2+1-(M+2)/2)
        N=M
        QQ(MODE)=XQ(N)
        DO I=1,NATOM
          DO J=1,3
            XR(I,J)=X0(I,J)+XL(I,MODE,J)*QQ(MODE)/SQRT(XM(I))
          END DO
        END DO
        CALL GETECK(NATOM,X0,XM,XX,XR,RR)
        WRITE(MOUT,204)NATOM
        WRITE(MOUT,205)MODE,CHSYM(MSYM),M
        WRITE(MOUT,*)'**********************************'
        DO I=1,NATOM
          WRITE(MOUT,206)SYMBOL(I),(XX(I,K)*BOHR,K=1,3)
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMEC2(XQ1,XQ2,MM1,MM2,NMODE,NATOM,QQ,XM,XR,X0,XL,
     1MODE1,MODE2,XX,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),SPACE,CHSYM(4)
      DIMENSION XQ1(MM1),XQ2(MM2),QQ(NMODE),X0(NATOM,3),XM(NATOM)
      DIMENSION XR(NATOM,3),XL(NATOM,NMODE,3)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM)
      COMMON/DUMP/JDUMP,IDUMP,KDUMP,MDUMP
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
204   FORMAT(I4)
205   FORMAT(1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2,
     1       1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2)
206   FORMAT(A2,1X,3F20.10)
208   FORMAT(3F20.10)
C**ONLY DUMP IF UNIQUE SYMMETRY
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE1)MSYM1=I
        END DO
      END DO
      ISKIP1=0
      IF(MSYM1.NE.1)ISKIP1=1
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE2)MSYM2=I
        END DO
      END DO
      ISKIP2=0
      IF(MSYM2.NE.1)ISKIP2=1
      N12=ISYMP(MSYM1,MSYM2)
C**SPECIAL CASE B2 WITH B2, B1 WITH B1, A2 WITH A2
      IF(MSYM1.NE.1.AND.(N12.EQ.1))ISKIP1=0

      DO K=1,NMODE
        QQ(K)=0
      END DO
      MDUM1=MDUMP
      IF(ISKIP1.NE.0)MDUM1=MDUMP/2
      DO M1=1,MDUM1
C       N1=MM1/2+1+((-1)**M1)*(MM1/2+1-(M1+2)/2)
        N1=M1
        QQ(MODE1)=XQ1(N1)
        MDUM2=MDUMP
        IF(ISKIP2.NE.0)MDUM2=MDUMP/2
        DO M2=1,MDUM2
C         N2=MM2/2+1+((-1)**M2)*(MM2/2+1-(M2+2)/2)
          N2=M2
          QQ(MODE2)=XQ2(N2)
          DO I=1,NATOM
            DO J=1,3
              XR(I,J)=X0(I,J)+XL(I,MODE1,J)*QQ(MODE1)/SQRT(XM(I))
              XR(I,J)=XR(I,J)+XL(I,MODE2,J)*QQ(MODE2)/SQRT(XM(I))
            END DO
          END DO
          CALL GETECK(NATOM,X0,XM,XX,XR,RR)
          WRITE(MOUT,204)NATOM
          WRITE(MOUT,205)MODE1,CHSYM(MSYM1),M1,MODE2,CHSYM(MSYM2),M2
          WRITE(MOUT,*)'**********************************'
          DO I=1,NATOM
            WRITE(MOUT,206)SYMBOL(I),(XX(I,K)*BOHR,K=1,3)
          END DO
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMEC3(XQ1,XQ2,XQ3,MM1,MM2,MM3,NMODE,NATOM,QQ,XM,
     1XR,X0,XL,MODE1,MODE2,MODE3,XX,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),SPACE,CHSYM(4)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION QQ(NMODE),X0(NATOM,3),XM(NATOM)
      DIMENSION XR(NATOM,3),XL(NATOM,NMODE,3)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM)
      COMMON/DUMP/JDUMP,IDUMP,KDUMP,MDUMP
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
204   FORMAT(I4)
205   FORMAT(1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2,
     1       1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2,
     2       1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2)
206   FORMAT(A2,1X,3F20.10)
208   FORMAT(4F20.10)
C**ONLY DUMP IF UNIQUE SYMMETRY
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE1)MSYM1=I
        END DO
      END DO
      ISKIP1=0
      IF(MSYM1.NE.1)ISKIP1=1
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE2)MSYM2=I
        END DO
      END DO
      ISKIP2=0
      IF(MSYM2.NE.1)ISKIP2=1
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE3)MSYM3=I
        END DO
      END DO
      ISKIP3=0
      IF(MSYM3.NE.1)ISKIP3=1
C**SPECIAL CASE B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2, ETC.
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM2))ISKIP1=0
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM3))ISKIP1=0
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM3))ISKIP2=0

C**SPECIAL CASE S2*S3 = S1
      N23=ISYMP(MSYM2,MSYM3)
      IF(N23.EQ.MSYM1)ISKIP1=0

      DO K=1,NMODE
        QQ(K)=0
      END DO
      MDUM1=MDUMP
      IF(ISKIP1.NE.0)MDUM1=MDUMP/2
      DO M1=1,MDUM1
C       N1=MM1/2+1+((-1)**M1)*(MM1/2+1-(M1+2)/2)
        N1=M1
        QQ(MODE1)=XQ1(N1)
        MDUM2=MDUMP
        IF(ISKIP2.NE.0)MDUM2=MDUMP/2
        DO M2=1,MDUM2
C         N2=MM2/2+1+((-1)**M2)*(MM2/2+1-(M2+2)/2)
          N2=M2
          QQ(MODE2)=XQ2(N2)
          MDUM3=MDUMP
          IF(ISKIP3.NE.0)MDUM3=MDUMP/2
          DO M3=1,MDUM3
C           N3=MM3/2+1+((-1)**M3)*(MM3/2+1-(M3+2)/2)
            N3=M3
            QQ(MODE3)=XQ3(N3)
            DO I=1,NATOM
              DO J=1,3
                XR(I,J)=X0(I,J)+XL(I,MODE1,J)*QQ(MODE1)/SQRT(XM(I))
                XR(I,J)=XR(I,J)+XL(I,MODE2,J)*QQ(MODE2)/SQRT(XM(I))
                XR(I,J)=XR(I,J)+XL(I,MODE3,J)*QQ(MODE3)/SQRT(XM(I))
              END DO
            END DO
            CALL GETECK(NATOM,X0,XM,XX,XR,RR)
            WRITE(MOUT,204)NATOM
            WRITE(MOUT,205)MODE1,CHSYM(MSYM1),M1,MODE2,CHSYM(MSYM2),M2,
     1      MODE3,CHSYM(MSYM2),M3
            WRITE(MOUT,*)'**********************************'
            DO I=1,NATOM
              WRITE(MOUT,206)SYMBOL(I),(XX(I,K)*BOHR,K=1,3)
            END DO
          END DO
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMEC4(XQ1,XQ2,XQ3,XQ4,MM1,MM2,MM3,MM4,NMODE,NATOM,
     1QQ,XM,XR,X0,XL,MODE1,MODE2,MODE3,MODE4,XX,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),SPACE,CHSYM(4)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION QQ(NMODE),X0(NATOM,3),XM(NATOM)
      DIMENSION XR(NATOM,3),XL(NATOM,NMODE,3)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM)
      COMMON/DUMP/JDUMP,IDUMP,KDUMP,MDUMP
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
204   FORMAT(I4)
205   FORMAT(1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2,
     1       1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2,
     2       1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2,
     3       1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2)
206   FORMAT(A2,1X,3F20.10)
208   FORMAT(5F20.10)
C**ONLY DUMP POTENTIAL DATA IF UNIQUE SYMMETRY CONFIGURATION
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE1)MSYM1=I
        END DO
      END DO
      ISKIP1=0
      IF(MSYM1.NE.1)ISKIP1=1
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE2)MSYM2=I
        END DO
      END DO
      ISKIP2=0
      IF(MSYM2.NE.1)ISKIP2=1
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE3)MSYM3=I
        END DO
      END DO
      ISKIP3=0
      IF(MSYM3.NE.1)ISKIP3=1
      DO I=1,NWSYM
        DO K=1,NSYM(I)
         IF(ISYM(I,K).EQ.MODE4)MSYM4=I
        END DO
      END DO
      ISKIP4=0
      IF(MSYM4.NE.1)ISKIP4=1
C**SPECIAL CASE B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2, ETC.
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM2))ISKIP1=0
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM3))ISKIP1=0
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM4))ISKIP1=0
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM3))ISKIP2=0
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM4))ISKIP2=0
      IF(MSYM3.NE.1.AND.(MSYM3.EQ.MSYM4))ISKIP3=0

C**SPECIAL CASE S2*S3 = S1
      N23=ISYMP(MSYM2,MSYM3)
      IF(N23.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S2*S4 = S1
      N24=ISYMP(MSYM2,MSYM4)
      IF(N24.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S4 = S1
      N34=ISYMP(MSYM3,MSYM4)
      IF(N34.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S4 = S2
      IF(N34.EQ.MSYM2)ISKIP2=0

C**SPECIAL CASE S2*S3*S4 = S1
      N234=ISYMP(MSYM2,N34)
      IF(N234.EQ.MSYM1)ISKIP1=0

      DO K=1,NMODE
        QQ(K)=0
      END DO
      MDUM1=MDUMP
      IF(ISKIP1.NE.0)MDUM1=MDUMP/2
      DO M1=1,MDUM1
C       N1=MM1/2+1+((-1)**M1)*(MM1/2+1-(M1+2)/2)
        N1=M1
        QQ(MODE1)=XQ1(N1)
        MDUM2=MDUMP
        IF(ISKIP2.NE.0)MDUM2=MDUMP/2
        DO M2=1,MDUM2
C         N2=MM2/2+1+((-1)**M2)*(MM2/2+1-(M2+2)/2)
          N2=M2
          QQ(MODE2)=XQ2(N2)
          MDUM3=MDUMP
          IF(ISKIP3.NE.0)MDUM3=MDUMP/2
          DO M3=1,MDUM3
C           N3=MM3/2+1+((-1)**M3)*(MM3/2+1-(M3+2)/2)
            N3=M3
            QQ(MODE3)=XQ3(N3)
            MDUM4=MDUMP
            IF(ISKIP4.NE.0)MDUM4=MDUMP/2
            DO M4=1,MDUM4
C             N4=MM4/2+1+((-1)**M4)*(MM4/2+1-(M4+2)/2)
              N4=M4
              QQ(MODE4)=XQ4(N4)
              DO I=1,NATOM
                DO J=1,3
                  XR(I,J)=X0(I,J)+XL(I,MODE1,J)*QQ(MODE1)/SQRT(XM(I))
                  XR(I,J)=XR(I,J)+XL(I,MODE2,J)*QQ(MODE2)/SQRT(XM(I))
                  XR(I,J)=XR(I,J)+XL(I,MODE3,J)*QQ(MODE3)/SQRT(XM(I))
                  XR(I,J)=XR(I,J)+XL(I,MODE4,J)*QQ(MODE4)/SQRT(XM(I))
                END DO
              END DO
              CALL GETECK(NATOM,X0,XM,XX,XR,RR)
              WRITE(MOUT,204)NATOM
              WRITE(MOUT,205)MODE1,CHSYM(MSYM1),M1,MODE2,CHSYM(MSYM2),
     1        M2,MODE3,CHSYM(MSYM3),M3,MODE4,CHSYM(MSYM4),M4
              WRITE(MOUT,*)'**********************************'
              DO I=1,NATOM
                WRITE(MOUT,206)SYMBOL(I),(XX(I,K)*BOHR,K=1,3)
              END DO
            END DO
          END DO
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETECK(NATOM,X0,XM,XX,XR,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PLANAR(3),PLNMOL
      DIMENSION XM(NATOM),XX(NATOM,3),X0(NATOM,3)
      DIMENSION XR(NATOM,3),RR(NATOM,NATOM)
      DIMENSION R(3,3)
      DIMENSION WRK(100),SOL(3),F(3)
      COMMON/RETURN/IRET
      COMMON/ROTIND/ANGSAV(3),INDROT
      COMMON/XACC/XACC
      COMMON/PLMNMX/MXCD(3)
      COMMON/FILASS/IOUT,INP
      COMMON/ECKCNT/ICNT,INTC
      COMMON/ECKIND/INDECK
      COMMON/PLANE/PLNMOL,PLANAR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
C*****************************************************
      EXTERNAL FUNCT,MONIT
C*****************************************************
200   FORMAT(//,1X,'DISPLACED INTERNAL COORDINATES',/)
210   FORMAT(//,1X,'ECKART BOND DISTANCES',/)
220   FORMAT(//,1X,'DISTORTED INTERNAL COORDINATES',/)
240   FORMAT(1X,'ATOM ',I2,':  X0 =',F12.7,'  Y0 =',F12.7,'  Z0 =',
     1F12.7)
245   FORMAT(//,60(1H*),/,1X,'DISPLACED GEOMETRY',2X,I4/)
250   FORMAT(//,1X,'ECKART GEOMETRY',/)
255   FORMAT(/,1X,'ALPHA,BETA,GAMMA',/)
260   FORMAT(//,1X,'PREVIOUS GEOMETRY',/)
265   FORMAT(//,1X,'REVISED ECKART GEOMETRY',/)
C*****************************************************
      IF(IRET.NE.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'*********************START ECKART'
C***********************************
C**XR CONTAINS DISTORTED GEOMETRY
        CALL CHECKM(XM,XR,XX,RR,NATOM)
C**XX CONTAINS DISTORTED GEOMETRY AT CENTRE OF MASS
C***********************************
        ICNT=ICNT+1
        WRITE(IOUT,245)ICNT
        DO I=1,NATOM
          WRITE(IOUT,240)I,(XX(I,J),J=1,3)
        END DO
      END IF
C***********************************
C********************************************************
C********************************************************
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
      IF(PLANAR(1))WRITE(IOUT,*)'PLANAR IN X'
      IF(PLANAR(2))WRITE(IOUT,*)'PLANAR IN Y'
      IF(PLANAR(3))WRITE(IOUT,*)'PLANAR IN Z'
C*********************************************************
C*********************************************************
C**PUT PLANAR MOLECULE INTO CORRECT PLANE
      IF(PLANAR(1))THEN
C**************
C**X COORD ZERO
C**************
        IF(MXCD(1).EQ.2)THEN
C**PREVIOUS Y COORD APPROACHING ZERO
          DO I=1,NATOM
            SAVE=XX(I,1)
            XX(I,1)=XX(I,2)
            XX(I,2)=SAVE
          END DO
        END IF
        IF(MXCD(1).EQ.3)THEN
C**PREVIOUS Z COORD APPROACHING ZERO
          DO I=1,NATOM
            SAVE=XX(I,1)
            XX(I,1)=XX(I,3)
            XX(I,3)=SAVE
          END DO
        END IF
      END IF
      IF(PLANAR(2))THEN
C**************
C**Y COORD ZERO
C**************
        IF(MXCD(1).EQ.1)THEN
C**PREVIOUS X COORD APPROACHING ZERO
          DO I=1,NATOM
            SAVE=XX(I,2)
            XX(I,2)=XX(I,1)
            XX(I,1)=SAVE
          END DO
        END IF
        IF(MXCD(1).EQ.3)THEN
C**PREVIOUS Z COORD APPROACHING ZERO
          DO I=1,NATOM
            SAVE=XX(I,2)
            XX(I,2)=XX(I,3)
            XX(I,3)=SAVE
          END DO
        END IF
      END IF
      IF(PLANAR(3))THEN
C**************
C**Z COORD ZERO
C**************
        IF(MXCD(1).EQ.1)THEN
C**PREVIOUS X COORD APPROACHING ZERO
          DO I=1,NATOM
            SAVE=XX(I,3)
            XX(I,3)=XX(I,1)
            XX(I,1)=SAVE
          END DO
        END IF
        IF(MXCD(1).EQ.2)THEN
C**PREVIOUS Y COORD APPROACHING ZERO
          DO I=1,NATOM
            SAVE=XX(I,3)
            XX(I,3)=XX(I,2)
            XX(I,2)=SAVE
          END DO
        END IF
      END IF
      IF(PLNMOL)THEN
        DO K=1,3
          PLANAR(K)=(MXCD(1).EQ.K)
        END DO
      END IF
C*********************************************************
C*********************************************************
      IF(PLANAR(1))WRITE(IOUT,*)'PREVIOUSLY PLANAR IN X'
      IF(PLANAR(2))WRITE(IOUT,*)'PREVIOUSLY PLANAR IN Y'
      IF(PLANAR(3))WRITE(IOUT,*)'PREVIOUSLY PLANAR IN Z'
C********************************************************
C********************************************************
      IF(IRET.NE.0)THEN
C**XR CONTAINS DISTORTED GEOMETRY
        DO I=1,NATOM
          DO K=1,3
            XR(I,K)=XX(I,K)
          END DO
        END DO
        WRITE(IOUT,200)
        DO I=1,NATOM
          WRITE(IOUT,240)I,(XX(I,J),J=1,3)
        END DO
      ELSE
C**RR CONTAINS DISTORTED GEOMETRY
        DO I=1,NATOM
          DO K=1,3
            RR(I,K)=XX(I,K)
          END DO
        END DO
      END IF
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
        STOP 'WRK TOO SMALL IN GETECK'
      END IF
C**ANSWER NEAR ZERO, SO INITIAL GUESS SET TO ZERO
      SOL(1)=0
      SOL(2)=0
      SOL(3)=0
      CALL E04FBF(NFIT,MFIT,SOL,F,SUMSQ,FTOL,XTOL,STEP,WRK,100,
     1FUNCT,MONIT,KPRINT,MAXCAL,IFAIL)
      WRITE(IOUT,*)
      WRITE(IOUT,*)'GETECK IFAIL = ',IFAIL
      WRITE(IOUT,*)
      WRITE(IOUT,*)'ZEROs ',F(1),F(2),F(3)
      IF(IFAIL.NE.0.AND.IFAIL.NE.3)THEN
        STOP 'ERROR IN E04FBF'
      END IF
      AL=SOL(1)
      BE=SOL(2)
      GA=SOL(3)
      WRITE(IOUT,255)
      WRITE(IOUT,*)AL,BE,GA
      WRITE(IOUT,*)
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
            IF(IRET.NE.0)ZZZZ=XR(I,J)
            IF(IRET.EQ.0)ZZZZ=RR(I,J)
            XX(I,K)=XX(I,K)+R(K,J)*ZZZZ
          END DO
        END DO
      END DO
C**ECKART INTERNAL DISPLACEMENT COORDINATES
      WRITE(IOUT,250)
      DO I=1,NATOM
        WRITE(IOUT,240)I,(XX(I,J),J=1,3)
      END DO
C**LOOK FOR LARGEST COORDINATE
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
C**PUT PLANAR MOLECULE INTO CORRECT PLANE
      IF(PLANAR(1))THEN
C**************
C**X COORD ZERO
C**************
        IF(MXCD(3).EQ.2.AND.MAXCOD.EQ.3)THEN
C**PREVIOUS Y COORD WAS MAXIMUM
          DO I=1,NATOM
            SAVE=XX(I,2)
            XX(I,2)=XX(I,3)
            XX(I,3)=SAVE
          END DO
        END IF
        IF(MXCD(3).EQ.3.AND.MAXCOD.EQ.2)THEN
C**PREVIOUS Z COORD WAS MAXIMUM
          DO I=1,NATOM
            SAVE=XX(I,2)
            XX(I,2)=XX(I,3)
            XX(I,3)=SAVE
          END DO
        END IF
      END IF
      IF(PLANAR(2))THEN
C**************
C**Y COORD ZERO
C**************
        IF(MXCD(3).EQ.1.AND.MAXCOD.EQ.3)THEN
C**PREVIOUS X COORD WAS MAXIMUM
          DO I=1,NATOM
            SAVE=XX(I,1)
            XX(I,1)=XX(I,3)
            XX(I,3)=SAVE
          END DO
        END IF
        IF(MXCD(3).EQ.3.AND.MAXCOD.EQ.1)THEN
C**PREVIOUS Z COORD WAS MAXIMUM
          DO I=1,NATOM
            SAVE=XX(I,1)
            XX(I,1)=XX(I,3)
            XX(I,3)=SAVE
          END DO
        END IF
      END IF
      IF(PLANAR(3))THEN
C**************
C**Z COORD ZERO
C**************
        IF(MXCD(3).EQ.1.AND.MAXCOD.EQ.2)THEN
C**PREVIOUS X COORD WAS MAXIMUM
          DO I=1,NATOM
            SAVE=XX(I,1)
            XX(I,1)=XX(I,2)
            XX(I,2)=SAVE
          END DO
        END IF
        IF(MXCD(3).EQ.2.AND.MAXCOD.EQ.1)THEN
C**PREVIOUS Y COORD WAS MAXIMUM
          DO I=1,NATOM
            SAVE=XX(I,1)
            XX(I,1)=XX(I,2)
            XX(I,2)=SAVE
          END DO
        END IF
      END IF
C*********************************************************
C*********************************************************
9999  CONTINUE
      IF(PLNMOL)THEN
        IF(PLANAR(1))WRITE(IOUT,*)'PLANAR IN X'
        IF(PLANAR(2))WRITE(IOUT,*)'PLANAR IN Y'
        IF(PLANAR(3))WRITE(IOUT,*)'PLANAR IN Z'
        WRITE(IOUT,265)
        DO I=1,NATOM
          WRITE(IOUT,240)I,(XX(I,J),J=1,3)
        END DO
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
      WRITE(IOUT,*)'MXCD ',(MXCD(I),I=1,3)
      INTC=-1
      WRITE(IOUT,210)
      CALL BONDS(NATOM,RR,XX)
      DO I=1,NATOM
        WRITE(IOUT,*)(RR(I,J),J=I+1,NATOM)
      END DO
      INTC=0
      IF(IRET.NE.0)THEN
        WRITE(IOUT,*)'*********************END ECKART'
        WRITE(IOUT,*)
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE FUNCT(M,N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PLANAR(3),PLNMOL
      DIMENSION X(N),F(M),R(3,3)
      COMMON/RETURN/IRET
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
      IF(IRET.NE.0)
     1CALL GETZM(NATOM,R,W(LYZ),W(LYZ+3*NATOM),W(LXM),W(LX0),F(1))
      IF(IRET.EQ.0)
     1CALL GETMZ(NATOM,R,W(LYZ),W(LYZ+3*NATOM),W(LXM),W(LX0),F(1),
     2W(LRR))
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETMZ(NATOM,R,XR,XRP,XM,X0,F,XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PLANAR(3),PLNMOL
      DIMENSION R(3,3),XR(NATOM,3),XRP(NATOM,3),X0(NATOM,3),F(3)
      DIMENSION XM(NATOM),XX(NATOM,3)
      COMMON/PLANE/PLNMOL,PLANAR
      DO I=1,NATOM
        DO K=1,3
          XRP(I,K)=0
          DO J=1,3
            XRP(I,K)=XRP(I,K)+R(K,J)*XX(I,J)
          END DO
        END DO
      END DO
      F(1)=0
      F(2)=0
      F(3)=0
      IF(PLNMOL)THEN
        DO I=1,NATOM
          IF(PLANAR(1))
     1    F(1)=F(1)-XM(I)*(X0(I,2)*XRP(I,3)-X0(I,3)*XRP(I,2))
          IF(PLANAR(2))
     1    F(1)=F(1)-XM(I)*(X0(I,1)*XRP(I,3)-X0(I,3)*XRP(I,1))
          IF(PLANAR(3))
     1    F(1)=F(1)-XM(I)*(X0(I,1)*XRP(I,2)-X0(I,2)*XRP(I,1))
        END DO
        DO I=1,NATOM
          IF(PLANAR(1))
     1    F(1)=F(1)+XM(I)*(XR(I,2)*XRP(I,3)-XR(I,3)*XRP(I,2))
          IF(PLANAR(2))
     1    F(1)=F(1)+XM(I)*(XR(I,1)*XRP(I,3)-XR(I,3)*XRP(I,1))
          IF(PLANAR(3))
     1    F(1)=F(1)+XM(I)*(XR(I,1)*XRP(I,2)-XR(I,2)*XRP(I,1))
        END DO
      ELSE
        DO I=1,NATOM
          F(1)=F(1)-XM(I)*(X0(I,2)*XRP(I,3)-X0(I,3)*XRP(I,2))
          F(2)=F(2)-XM(I)*(X0(I,1)*XRP(I,3)-X0(I,3)*XRP(I,1))
          F(3)=F(3)-XM(I)*(X0(I,1)*XRP(I,2)-X0(I,2)*XRP(I,1))
        END DO
        DO I=1,NATOM
          F(1)=F(1)+XM(I)*(XR(I,2)*XRP(I,3)-XR(I,3)*XRP(I,2))
          F(2)=F(2)+XM(I)*(XR(I,1)*XRP(I,3)-XR(I,3)*XRP(I,1))
          F(3)=F(3)+XM(I)*(XR(I,1)*XRP(I,2)-XR(I,2)*XRP(I,1))
        END DO
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETZM(NATOM,R,XR,XRP,XM,X0,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(3,3),XR(NATOM,3),XRP(NATOM,3),X0(NATOM,3),F(3)
      DIMENSION XM(NATOM)
      COMMON/FILASS/IOUT,INP
      DO I=1,NATOM
        DO K=1,3
          XRP(I,K)=0
          DO J=1,3
            XRP(I,K)=XRP(I,K)+R(K,J)*XR(I,J)
          END DO
        END DO
      END DO
      F(1)=0
      F(2)=0
      F(3)=0
      DO I=1,NATOM
        F(1)=F(1)+XM(I)*(X0(I,2)*XRP(I,3)-X0(I,3)*XRP(I,2))
        F(2)=F(2)+XM(I)*(X0(I,1)*XRP(I,3)-X0(I,3)*XRP(I,1))
        F(3)=F(3)+XM(I)*(X0(I,1)*XRP(I,2)-X0(I,2)*XRP(I,1))
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE MONIT(M,N,X,F,S,NCALL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**MONITOR PROGRESS OF E04FBF IF REQUIRED
      DIMENSION X(N),F(M)
      RETURN
      END
