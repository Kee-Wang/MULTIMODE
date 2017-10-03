C****************************************************************
C****************************************************************
C**ROT
C****************************************************************
C****************************************************************
      SUBROUTINE ROTC(XX,XM,NATOM,A,B,C,BBAR,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),XM(NATOM),WR(3),E(3)
      COMMON/MOMI/XK(3,3),XMU(3,3)

      IF(IND.NE.0)THEN
C**SET UP MOMENT OF INERTIA MATRIX
        DO IX=1,3
          DO IY=1,3  ! chen
            XK(IY,IX)=0
          END DO
        END DO

        DO IX=1,3
          DO I=1,NATOM
            X=0
            DO J=1,3
              X=X+XX(I,J)*XX(I,J)  ! chen
            END DO
            XK(IX,IX)=XK(IX,IX)+XM(I)*X
          END DO
        END DO

        DO IX=1,3
          DO IY=1,IX  ! chen
            DO I=1,NATOM
              XK(IY,IX)=XK(IY,IX)-XM(I)*XX(I,IY)*XX(I,IX)
            END DO
          END DO
        END DO

        DO IX=1,3
          DO IY=1,IX-1  ! chen
            XK(IX,IY)=XK(IY,IX)
          END DO
        END DO
      END IF
C**DIAGONALISE MOMENT OF INERTIA MATRIX
      CALL DIAG(XK,XK,3,3,-1,WR,E,WR,3,3,XK,3,3,XK,XK,XK,IDUM,IDUM,
     1IDUM)
C**E ARE IN ASCENDING ORDER
      IF(E(1).GT.1.D-35)THEN
        A=1/(2*E(1))
      ELSE
        A=1.D+35
      END IF
      IF(E(2).GT.1.D-35)THEN
        B=1/(2*E(2))
      ELSE
        B=1.D+35
      END IF
      IF(E(3).GT.1.D-35)THEN
        C=1/(2*E(3))
      ELSE
        C=1.D+35
      END IF
      BBAR=(B+C)/2
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI0(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,V,
     1VR,S,J21,V0,W0,E0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21)
      REAL*4 VR(J21)
      DIMENSION QQ(NMODE)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION XZ(NMODE,NMODE,3),AB(NMODE,3),B(NMODE,NMODE)
      DIMENSION AA(NMODE,3,3),BB(NMODE)
      DIMENSION S(J21,9,J21)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
100   FORMAT(//,1X,'ROTATIONAL CONSTANTS AT EQUILIBRIUM ',3F10.5,///)
      CALL ROTC(X0,XM,NATOM,AROT,BROT,CROT,BBAR,1)
      WRITE(IOUT,100)AROT*WAVENM,BROT*WAVENM,CROT*WAVENM
      MX=1
      MY=2
      MZ=3
      DO K=1,NMODE
        QQ(K)=0
        MULT(K)=0
      END DO
      CALL CORIOL(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,ZZ)
      IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J)=(S(J,1,J)*XMU(MX,MX)+S(J,2,J)*XMU(MY,MY)+
     1      S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        ELSE
          DO J=1,J21
            VR(J)=(S(J,1,J)*XMU(MX,MX)+S(J,2,J)*XMU(MY,MY)+
     1      S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        END IF
      ELSE
C**ADIABATIC ROTATION
        CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,0)
        DO I=1,J21
          DO J=1,J21
            V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
          END DO
        END DO
        CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,V0,V0,
     1  V0,IDUM,IDUM,IDUM)
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J)=E0(J)
          END DO
        ELSE
          DO J=1,J21
            VR(J)=E0(J)
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMV0(NMODE,NATOM,QQ,BB,BS,BT,BBS,B,BBT,
     1BSS,XX,XXP,X0,XL,XM,V,VR,S,J21,V0,W0,E0,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21)
      REAL*4 VR(J21)
      DIMENSION MODINT(NMODE),XM(NATOM)
CCCC  DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,362)
CCCC  DIMENSION BBS(NMODE,NMODE,362),BSS(NMODE,362),BBT(NMODE)
CCCC  DIMENSION BS(NMODE,3,362),BT(NMODE,NMODE),B(NMODE,3,3,362)
CCCC  DIMENSION XX(NATOM,3,362),X0(NATOM,3),XL(NATOM,NMODE,3,362)
      DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,722)
      DIMENSION BBS(NMODE,NMODE,722),BSS(NMODE,722),BBT(NMODE)
      DIMENSION BS(NMODE,3,722),BT(NMODE,NMODE),B(NMODE,3,3,722)
      DIMENSION XX(NATOM,3,722),X0(NATOM,3),XL(NATOM,NMODE,3,722)
      DIMENSION XXP(NATOM,3,722)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION S(J21,9,J21)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
C**ASSUME STARTING VALUE TAU: INIT
C**ASSUME NEXT VALUE TAU: INIT+INCTAU
C**ITAU=2 EQUIVALENT TO ITAU=722
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/VIBMOD/NSMODE,NVMODE
      ITAU=INIT-INCTAU
      MX=1
      MY=2
      MZ=3
      MDT=MODINT(NSMODE)
C**LOOP ROUND TAU
      DO MTAU=1,MMTAU/MDT
        ITAU=ITAU+INCTAU
CCCC    IF(ITAU.GT.362)ITAU=ITAU-360
        IF(ITAU.GT.722)ITAU=ITAU-720
        DO K=1,NMODE
          QQ(K)=0
          MULT(K)=0
        END DO
        CALL MILLER(NMODE,NATOM,QQ,BB,BS,BBS,B,BSS,XX,XL,XXP,XM,
     1  BBT,BT,ZZ,ITAU,1)
        IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
          IF(JCOUPC.GE.0)THEN
            DO J=1,J21
              V(J)=(S(J,1,J)*XMU0(MX,MX)+S(J,2,J)*XMU0(MY,MY)+
     1        S(J,3,J)*XMU0(MZ,MZ))/2
            END DO
          ELSE
            DO J=1,J21
              VR(J)=(S(J,1,J)*XMU0(MX,MX)+S(J,2,J)*XMU0(MY,MY)+
     1        S(J,3,J)*XMU0(MZ,MZ))/2
            END DO
          END IF
        ELSE
C**ADIABATIC ROTATION - BEWARE XMU
          CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,1)
          DO I=1,J21
            DO J=1,J21
              V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
            END DO
          END DO
          CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,V0,V0,
     1    V0,IDUM,IDUM,IDUM)
          IF(JCOUPC.GE.0)THEN
            DO J=1,J21
              V(J)=E0(J)
            END DO
          ELSE
            DO J=1,J21
              VR(J)=E0(J)
            END DO
          END IF
        END IF
        IF(JCOUPC.GE.0)THEN
          WRITE(61)V
        ELSE
          WRITE(61)VR
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI1(XQ,MM,NMODE,NATOM,QQ,XZ,AB,B,AA,BB,
     1XX,X0,XL,XM,MODE,V,VR,S,J21,V0,W0,E0,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM)
      REAL*4 VR(J21,MM)
      DIMENSION MODINT(NMODE)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION XZ(NMODE,NMODE,3),AB(NMODE,3),B(NMODE,NMODE)
      DIMENSION AA(NMODE,3,3),BB(NMODE)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ(MM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/TORS/QTAU
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
      MX=1
      MY=2
      MZ=3
      DO K=1,NMODE
        QQ(K)=0
        MULT(K)=0
      END DO
      MULT(MODE)=1
      MD=MODINT(MODE)
      DO M=1,MM/MD
        QQ(MODE)=XQ(M)
        CALL CORIOL(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,ZZ)
        IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
          CALL SET61(V,VR,MM/MD,E0,J21,M,S,MX,MY,MZ,0)
C         IF(JCOUPC.GT.0)THEN
C           DO J=1,J21
C             V(J,M)=(S(J,1,J)*XMU(MX,MX)+S(J,2,J)*XMU(MY,MY)+
C    1        S(J,3,J)*XMU(MZ,MZ))/2
C           END DO
C         ELSE
C           DO J=1,J21
C             VR(J,M)=(S(J,1,J)*XMU(MX,MX)+S(J,2,J)*XMU(MY,MY)+
C    1        S(J,3,J)*XMU(MZ,MZ))/2
C           END DO
C         END IF
        ELSE
C**ADIABATIC ROTATION
          CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,0)
          DO I=1,J21
            DO J=1,J21
              V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
            END DO
          END DO
          CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,V0,V0,
     1    V0,IDUM,IDUM,IDUM)
          CALL SET61(V,VR,MM/MD,E0,J21,M,S,MX,MY,MZ,1)
C         IF(JCOUPC.GT.0)THEN
C           DO J=1,J21
C             V(J,M)=E0(J)
C           END DO
C         ELSE
C           DO J=1,J21
C             VR(J,M)=E0(J)
C           END DO
C         END IF
        END IF
      END DO
      CALL OUT61(V,VR,MM/MD,J21)
C     IF(JCOUPC.GT.0)THEN
C       WRITE(I61)V
C     ELSE
C       WRITE(I61)VR
C     END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMV1(XQ,MM,NMODE,NATOM,QQ,BB,BS,BT,BBS,B,BBT,
     1BSS,XX,XXP,X0,XL,XM,MODE,V,VR,S,J21,V0,W0,E0,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM)
      REAL*4 VR(J21,MM)
      DIMENSION MODINT(NMODE),XM(NATOM)
CCCC  DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,362)
CCCC  DIMENSION BBS(NMODE,NMODE,362),BSS(NMODE,362),BBT(NMODE)
CCCC  DIMENSION BS(NMODE,3,362),BT(NMODE,NMODE),B(NMODE,3,3,362)
CCCC  DIMENSION XX(NATOM,3,362),X0(NATOM,3),XL(NATOM,NMODE,3,362)
      DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,722)
      DIMENSION BBS(NMODE,NMODE,722),BSS(NMODE,722),BBT(NMODE)
      DIMENSION BS(NMODE,3,362),BT(NMODE,NMODE),B(NMODE,3,3,722)
      DIMENSION XX(NATOM,3,722),X0(NATOM,3),XL(NATOM,NMODE,3,722)
      DIMENSION XXP(NATOM,3,722)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ(MM)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
C**ASSUME STARTING VALUE TAU: INIT
C**ASSUME NEXT VALUE TAU: INIT+INCTAU
C**ITAU=2 EQUIVALENT TO ITAU=722
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/VIBMOD/NSMODE,NVMODE
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      ITAU=INIT-INCTAU
      MX=1
      MY=2
      MZ=3
      MDT=MODINT(NSMODE)
      MD1=MODINT(MODE)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MODE.EQ.ISYM(I,J))N1=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
C**LOOP ROUND TAU
      DO MTAU=1,MMTAU/MDT
        ITAU=ITAU+INCTAU
CCCC    IF(ITAU.GT.362)ITAU=ITAU-360
        IF(ITAU.GT.722)ITAU=ITAU-720
        DO K=1,NMODE
          QQ(K)=0
          MULT(K)=0
        END DO
        MULT(MODE)=1
        DO M=1,MM/MD1
          QQ(MODE)=XQ(M)
          CALL MILLER(NMODE,NATOM,QQ,BB,BS,BBS,B,BSS,XX,XL,XXP,XM,
     1    BBT,BT,ZZ,ITAU,1)
          IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
            CALL VSET61(V,VR,MM/MD1,E0,J21,M,S,MX,MY,MZ,0)
C           IF(JCOUPC.GT.0)THEN
C             DO J=1,J21
C               V(J,M)=(S(J,1,J)*XMU0(MX,MX)+S(J,2,J)*XMU0(MY,MY)+
C    1          S(J,3,J)*XMU0(MZ,MZ))/2
C             END DO
C           ELSE
C             DO J=1,J21
C               VR(J,M)=(S(J,1,J)*XMU0(MX,MX)+S(J,2,J)*XMU0(MY,MY)+
C    1          S(J,3,J)*XMU0(MZ,MZ))/2
C             END DO
C           END IF
          ELSE
C**ADIABATIC ROTATION - BEWARE XMU
            CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,1)
            DO I=1,J21
              DO J=1,J21
                V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
              END DO
            END DO
            CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,V0,
     1      V0,V0,IDUM,IDUM,IDUM)
            CALL VSET61(V,VR,MM/MD1,E0,J21,M,S,MX,MY,MZ,1)
C           IF(JCOUPC.GT.0)THEN
C             DO J=1,J21
C               V(J,M)=E0(J)
C             END DO
C           ELSE
C             DO J=1,J21
C               VR(J,M)=E0(J)
C             END DO
C           END IF
          END IF
        END DO
        CALL OUT61(V,VR,MM/MD1,J21)
C       IF(JCOUPC.GT.0)THEN
C         WRITE(I61)V
C       ELSE
C         WRITE(I61)VR
C       END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI2(XQ1,XQ2,MM1,MM2,NMODE,NATOM,QQ,XZ,AB,B,AA,BB,
     1XX,X0,XL,XM,MODE1,MODE2,V,VR,S,J21,V0,W0,E0,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM2,MM1)
      REAL*4 VR(J21,MM2,MM1)
      DIMENSION MODINT(NMODE)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION XZ(NMODE,NMODE,3),AB(NMODE,3),B(NMODE,NMODE)
      DIMENSION AA(NMODE,3,3),BB(NMODE)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ1(MM1),XQ2(MM2),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/TORS/QTAU
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      MX=1
      MY=2
      MZ=3
      DO K=1,NMODE
        QQ(K)=0
        MULT(K)=0
      END DO
      MULT(MODE1)=1
      MULT(MODE2)=1
      MD1=MODINT(MODE1)
      MD2=MODINT(MODE2)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MODE1.EQ.ISYM(I,J))N1=I
          IF(MODE2.EQ.ISYM(I,J))N2=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      DO M1=1,MM1/MD1
        QQ(MODE1)=XQ1(M1)
        DO M2=1,MM2/MD2
          QQ(MODE2)=XQ2(M2)
          CALL CORIOL(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,ZZ)
          IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
            CALL SET62(V,VR,MM1/MD1,MM2/MD2,E0,J21,M1,M2,S,MX,MY,
     1      MZ,0)
C           IF(JCOUPC.GT.0)THEN
C             DO J=1,J21
C               V(J,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
C    1          S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
C             END DO
C           ELSE
C             DO J=1,J21
C               VR(J,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
C    1          S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
C             END DO
C           END IF
          ELSE
C**ADIABATIC ROTATION
            CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,0)
            DO I=1,J21
              DO J=1,J21
                V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
              END DO
            END DO
            CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,V0,
     1      V0,V0,IDUM,IDUM,IDUM)
            CALL SET62(V,VR,MM1/MD1,MM2/MD2,E0,J21,M1,M2,S,MX,MY,
     1      MZ,1)
C           IF(JCOUPC.GT.0)THEN
C             DO J=1,J21
C               V(J,M2,M1)=E0(J)
C             END DO
C           ELSE
C             DO J=1,J21
C               VR(J,M2,M1)=E0(J)
C             END DO
C           END IF
          END IF
        END DO
      END DO
      CALL OUT62(V,VR,MM1/MD1,MM2/MD2,J21)
C     IF(JCOUPC.GT.0)THEN
C       WRITE(I62)V
C     ELSE
C       WRITE(I62)VR
C     END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMV2(XQ1,XQ2,MM1,MM2,NMODE,NATOM,QQ,BB,BS,BT,BBS,B,
     1BBT,BSS,XX,XXP,X0,XL,XM,MODE1,MODE2,V,VR,S,J21,V0,W0,E0,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM2,MM1)
      REAL*4 VR(J21,MM2,MM1)
      DIMENSION MODINT(NMODE),XM(NATOM)
CCCC  DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,362)
CCCC  DIMENSION BBS(NMODE,NMODE,362),BSS(NMODE,362),BBT(NMODE)
CCCC  DIMENSION BS(NMODE,3,362),BT(NMODE,NMODE),B(NMODE,3,3,362)
CCCC  DIMENSION XX(NATOM,3,362),X0(NATOM,3),XL(NATOM,NMODE,3,362)
      DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,722)
      DIMENSION BBS(NMODE,NMODE,722),BSS(NMODE,722),BBT(NMODE)
      DIMENSION BS(NMODE,3,722),BT(NMODE,NMODE),B(NMODE,3,3,722)
      DIMENSION XX(NATOM,3,722),X0(NATOM,3),XL(NATOM,NMODE,3,722)
      DIMENSION XXP(NATOM,3,722)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ1(MM1),XQ2(MM2)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
C**ASSUME STARTING VALUE TAU: INIT
C**ASSUME NEXT VALUE TAU: INIT+INCTAU
C**ITAU=2 EQUIVALENT TO ITAU=722
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/VIBMOD/NSMODE,NVMODE
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      ITAU=INIT-INCTAU
      MX=1
      MY=2
      MZ=3
      MDT=MODINT(NSMODE)
      MD1=MODINT(MODE1)
      MD2=MODINT(MODE2)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MODE1.EQ.ISYM(I,J))N1=I
          IF(MODE2.EQ.ISYM(I,J))N2=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      IF(N2.EQ.NT.AND.MDT.GT.1)MD2=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.GT.1)MD2=1
C**LOOP ROUND TAU
      DO MTAU=1,MMTAU/MDT
        ITAU=ITAU+INCTAU
CCCC    IF(ITAU.GT.362)ITAU=ITAU-360
        IF(ITAU.GT.722)ITAU=ITAU-720
        DO K=1,NMODE
          QQ(K)=0
          MULT(K)=0
        END DO
        MULT(MODE1)=1
        MULT(MODE2)=1
        DO M1=1,MM1/MD1
          QQ(MODE1)=XQ1(M1)
          DO M2=1,MM2/MD2
            QQ(MODE2)=XQ2(M2)
            CALL MILLER(NMODE,NATOM,QQ,BB,BS,BBS,B,BSS,XX,XL,XXP,
     1      XM,BBT,BT,ZZ,ITAU,1)
            IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
              CALL VSET62(V,VR,MM1/MD1,MM2/MD2,E0,J21,M1,M2,S,MX,
     1        MY,MZ,0)
C             IF(JCOUPC.GT.0)THEN
C               DO J=1,J21
C                 V(J,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
C    1            S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
C               END DO
C             ELSE
C               DO J=1,J21
C                 VR(J,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
C    1            S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
C               END DO
C             END IF
            ELSE
C**ADIABATIC ROTATION - BEWARE XMU
              CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,1)
              DO I=1,J21
                DO J=1,J21
                  V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
                END DO
              END DO
              CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,V0,
     1        V0,V0,IDUM,IDUM,IDUM)
              CALL VSET62(V,VR,MM1/MD1,MM2/MD2,E0,J21,M1,M2,S,MX,
     1        MY,MZ,1)
C             IF(JCOUPC.GT.0)THEN
C               DO J=1,J21
C                 V(J,M2,M1)=E0(J)
C               END DO
C             ELSE
C               DO J=1,J21
C                 VR(J,M2,M1)=E0(J)
C               END DO
C             END IF
            END IF
          END DO
        END DO
        CALL OUT62(V,VR,MM1/MD1,MM2/MD2,J21)
C       IF(JCOUPC.GT.0)THEN
C         WRITE(I62)V
C       ELSE
C         WRITE(I62)VR
C       END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI3(XQ1,XQ2,XQ3,MM1,MM2,MM3,NMODE,NATOM,QQ,XZ,AB,B,
     1AA,BB,XX,X0,XL,XM,MODE1,MODE2,MODE3,V,VR,S,J21,V0,W0,E0,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM3,MM2,MM1)
      REAL*4 VR(J21,MM3,MM2,MM1)
      DIMENSION MODINT(NMODE)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION XZ(NMODE,NMODE,3),AB(NMODE,3),B(NMODE,NMODE)
      DIMENSION AA(NMODE,3,3),BB(NMODE)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/TORS/QTAU
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      MX=1
      MY=2
      MZ=3
      DO K=1,NMODE
        QQ(K)=0
        MULT(K)=0
      END DO
      MULT(MODE1)=1
      MULT(MODE2)=1
      MULT(MODE3)=1
      MD1=MODINT(MODE1)
      MD2=MODINT(MODE2)
      MD3=MODINT(MODE3)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MODE1.EQ.ISYM(I,J))N1=I
          IF(MODE2.EQ.ISYM(I,J))N2=I
          IF(MODE3.EQ.ISYM(I,J))N3=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N2.EQ.N3)MD3=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      DO M1=1,MM1/MD1
        QQ(MODE1)=XQ1(M1)
        DO M2=1,MM2/MD2
          QQ(MODE2)=XQ2(M2)
          DO M3=1,MM3/MD3
            QQ(MODE3)=XQ3(M3)
            CALL CORIOL(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,ZZ)
            IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
              CALL SET63(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,E0,J21,M1,M2,M3,
     1        S,MX,MY,MZ,0)
C             IF(JCOUPC.GT.0)THEN
C               DO J=1,J21
C                 V(J,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
C    1            S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
C               END DO
C             ELSE
C               DO J=1,J21
C                 VR(J,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
C    1            S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
C               END DO
C             END IF
            ELSE
C**ADIABATIC ROTATION
              CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,0)
              DO I=1,J21
                DO J=1,J21
                  V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
                END DO
              END DO
              CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,
     1        V0,V0,V0,IDUM,IDUM,IDUM)
              CALL SET63(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,E0,J21,M1,M2,M3,
     1        S,MX,MY,MZ,1)
C             IF(JCOUPC.GT.0)THEN
C               DO J=1,J21
C                 V(J,M3,M2,M1)=E0(J)
C               END DO
C             ELSE
C               DO J=1,J21
C                 VR(J,M3,M2,M1)=E0(J)
C               END DO
C             END IF
            END IF
          END DO
        END DO
      END DO
      CALL OUT63(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,J21)
C     IF(JCOUPC.GT.0)THEN
C       WRITE(I63)V
C     ELSE
C       WRITE(I63)VR
C     END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMV3(XQ1,XQ2,XQ3,MM1,MM2,MM3,NMODE,NATOM,QQ,BB,BS,
     1BT,BBS,B,BBT,BSS,XX,XXP,X0,XL,XM,MODE1,MODE2,MODE3,V,VR,S,J21,
     2V0,W0,E0,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM3,MM2,MM1)
      REAL*4 VR(J21,MM3,MM2,MM1)
      DIMENSION MODINT(NMODE),XM(NATOM)
CCCC  DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,362)
CCCC  DIMENSION BBS(NMODE,NMODE,362),BSS(NMODE,362),BBT(NMODE)
CCCC  DIMENSION BS(NMODE,3,362),BT(NMODE,NMODE),B(NMODE,3,3,362)
CCCC  DIMENSION XX(NATOM,3,362),X0(NATOM,3),XL(NATOM,NMODE,3,362)
      DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,722)
      DIMENSION BBS(NMODE,NMODE,722),BSS(NMODE,722),BBT(NMODE)
      DIMENSION BS(NMODE,3,722),BT(NMODE,NMODE),B(NMODE,3,3,722)
      DIMENSION XX(NATOM,3,722),X0(NATOM,3),XL(NATOM,NMODE,3,722)
      DIMENSION XXP(NATOM,3,722)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
C**ASSUME STARTING VALUE TAU: INIT
C**ASSUME NEXT VALUE TAU: INIT+INCTAU
C**ITAU=2 EQUIVALENT TO ITAU=722
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/VIBMOD/NSMODE,NVMODE
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      ITAU=INIT-INCTAU
      MX=1
      MY=2
      MZ=3
      MDT=MODINT(NSMODE)
      MD1=MODINT(MODE1)
      MD2=MODINT(MODE2)
      MD3=MODINT(MODE3)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MODE1.EQ.ISYM(I,J))N1=I
          IF(MODE2.EQ.ISYM(I,J))N2=I
          IF(MODE3.EQ.ISYM(I,J))N3=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N2.EQ.N3)MD3=1
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      IF(N2.EQ.NT.AND.MDT.GT.1)MD2=1
      IF(N3.EQ.NT.AND.MDT.GT.1)MD3=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.GT.1)MD2=1
      IF(N1T.EQ.N3.AND.MDT.GT.1)MD3=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.GT.1)MD3=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.GT.1)MD3=1
C**LOOP ROUND TAU
      DO MTAU=1,MMTAU/MDT
        ITAU=ITAU+INCTAU
        IF(ITAU.GT.722)ITAU=ITAU-720
        DO K=1,NMODE
          QQ(K)=0
          MULT(K)=0
        END DO
        MULT(MODE1)=1
        MULT(MODE2)=1
        MULT(MODE3)=1
        DO M1=1,MM1/MD1
          QQ(MODE1)=XQ1(M1)
          DO M2=1,MM2/MD2
            QQ(MODE2)=XQ2(M2)
            DO M3=1,MM3/MD3
              QQ(MODE3)=XQ3(M3)
              CALL MILLER(NMODE,NATOM,QQ,BB,BS,BBS,B,BSS,XX,XL,XXP,
     1        XM,BBT,BT,ZZ,ITAU,1)
              IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
                CALL VSET63(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,E0,J21,M1,M2,
     1          M3,S,MX,MY,MZ,0)
C               IF(JCOUPC.GT.0)THEN
C                 DO J=1,J21
C                   V(J,M3,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
C    1              S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
C                 END DO
C               ELSE
C                 DO J=1,J21
C                   VR(J,M3,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
C    1              S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
C                 END DO
C               END IF
              ELSE
C**ADIABATIC ROTATION - BEWARE XMU
                CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,1)
                DO I=1,J21
                  DO J=1,J21
                    V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
                  END DO
                END DO
                CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,
     1          V0,V0,V0,IDUM,IDUM,IDUM)
                CALL VSET63(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,E0,J21,M1,M2,
     1          M3,S,MX,MY,MZ,1)
C               IF(JCOUPC.GT.0)THEN
C                 DO J=1,J21
C                   V(J,M3,M2,M1)=E0(J)
C                 END DO
C               ELSE
C                 DO J=1,J21
C                   VR(J,M3,M2,M1)=E0(J)
C                 END DO
C               END IF
              END IF
            END DO
          END DO
        END DO
        CALL OUT63(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,J21)
C       IF(JCOUPC.GT.0)THEN
C         WRITE(I63)V
C       ELSE
C         WRITE(I63)VR
C       END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI4(XQ1,XQ2,XQ3,XQ4,MM1,MM2,MM3,MM4,NMODE,NATOM,QQ,
     1XZ,AB,B,AA,BB,XX,X0,XL,XM,MODE1,MODE2,MODE3,MODE4,V,VR,S,J21,V0,
     2W0,E0,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM4,MM3,MM2,MM1)
      REAL*4 VR(J21,MM4,MM3,MM2,MM1)
      DIMENSION MODINT(NMODE)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION XZ(NMODE,NMODE,3),AB(NMODE,3),B(NMODE,NMODE)
      DIMENSION AA(NMODE,3,3),BB(NMODE)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/TORS/QTAU
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      MX=1
      MY=2
      MZ=3
      DO K=1,NMODE
        QQ(K)=0
        MULT(K)=0
      END DO
      MULT(MODE1)=1
      MULT(MODE2)=1
      MULT(MODE3)=1
      MULT(MODE4)=1
      MD1=MODINT(MODE1)
      MD2=MODINT(MODE2)
      MD3=MODINT(MODE3)
      MD4=MODINT(MODE4)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MODE1.EQ.ISYM(I,J))N1=I
          IF(MODE2.EQ.ISYM(I,J))N2=I
          IF(MODE3.EQ.ISYM(I,J))N3=I
          IF(MODE4.EQ.ISYM(I,J))N4=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N3.EQ.N4)MD4=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      DO M1=1,MM1/MD1
        QQ(MODE1)=XQ1(M1)
        DO M2=1,MM2/MD2
          QQ(MODE2)=XQ2(M2)
          DO M3=1,MM3/MD3
            QQ(MODE3)=XQ3(M3)
            DO M4=1,MM4/MD4
              QQ(MODE4)=XQ4(M4)
              CALL CORIOL(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,ZZ)
              IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
                CALL SET64(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,E0,J21,
     1          M1,M2,M3,M4,S,MX,MY,MZ,0)
C               IF(JCOUPC.GT.0)THEN
C                 DO J=1,J21
C                   V(J,M4,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
C    1              S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
C                 END DO
C               ELSE
C                 DO J=1,J21
C                   VR(J,M4,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
C    1              S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
C                 END DO
C               END IF
              ELSE
C**ADIABATIC ROTATION
                CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,0)
                DO I=1,J21
                  DO J=1,J21
                    V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
                  END DO
                END DO
                CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,
     1          V0,V0,V0,IDUM,IDUM,IDUM)
                CALL SET64(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,E0,J21,
     1          M1,M2,M3,M4,S,MX,MY,MZ,1)
C               IF(JCOUPC.GT.0)THEN
C                 DO J=1,J21
C                   V(J,M4,M3,M2,M1)=E0(J)
C                 END DO
C               ELSE
C                 DO J=1,J21
C                   VR(J,M4,M3,M2,M1)=E0(J)
C                 END DO
C               END IF
              END IF
            END DO
          END DO
        END DO
      END DO
      CALL OUT64(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,J21)
C     IF(JCOUPC.GT.0)THEN
C       WRITE(I64)V
C     ELSE
C       WRITE(I64)VR
C     END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMV4(XQ1,XQ2,XQ3,XQ4,MM1,MM2,MM3,MM4,NMODE,NATOM,QQ,
     1BB,BS,BT,BBS,B,BBT,BSS,XX,XXP,X0,XL,XM,MODE1,MODE2,MODE3,MODE4,V,
     2VR,S,J21,V0,W0,E0,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM4,MM3,MM2,MM1)
      REAL*4 VR(J21,MM4,MM3,MM2,MM1)
      DIMENSION MODINT(NMODE),XM(NATOM)
CCCC  DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,362)
CCCC  DIMENSION BBS(NMODE,NMODE,362),BSS(NMODE,362),BBT(NMODE)
CCCC  DIMENSION BS(NMODE,3,362),BT(NMODE,NMODE),B(NMODE,3,3,362)
CCCC  DIMENSION XX(NATOM,3,362),X0(NATOM,3),XL(NATOM,NMODE,3,362)
      DIMENSION QQ(NMODE),BB(NMODE,NMODE,3,722)
      DIMENSION BBS(NMODE,NMODE,722),BSS(NMODE,722),BBT(NMODE)
      DIMENSION BS(NMODE,3,722),BT(NMODE,NMODE),B(NMODE,3,3,722)
      DIMENSION XX(NATOM,3,722),X0(NATOM,3),XL(NATOM,NMODE,3,722)
      DIMENSION XXP(NATOM,3,722)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
C**ASSUME STARTING VALUE TAU: INIT
C**ASSUME NEXT VALUE TAU: INIT+INCTAU
C**ITAU=2 EQUIVALENT TO ITAU=722
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/VIBMOD/NSMODE,NVMODE
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      ITAU=INIT-INCTAU
      MX=1
      MY=2
      MZ=3
      MDT=MODINT(NSMODE)
      MD1=MODINT(MODE1)
      MD2=MODINT(MODE2)
      MD3=MODINT(MODE3)
      MD4=MODINT(MODE4)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MODE1.EQ.ISYM(I,J))N1=I
          IF(MODE2.EQ.ISYM(I,J))N2=I
          IF(MODE3.EQ.ISYM(I,J))N3=I
          IF(MODE4.EQ.ISYM(I,J))N4=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N3.EQ.N4)MD4=1
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      IF(N2.EQ.NT.AND.MDT.GT.1)MD2=1
      IF(N3.EQ.NT.AND.MDT.GT.1)MD3=1
      IF(N4.EQ.NT.AND.MDT.GT.1)MD4=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.GT.1)MD2=1
      IF(N1T.EQ.N3.AND.MDT.GT.1)MD3=1
      IF(N1T.EQ.N4.AND.MDT.GT.1)MD4=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.GT.1)MD3=1
      IF(N2T.EQ.N4.AND.MDT.GT.1)MD4=1
      N3T=ISYMP(N3,NT)
      IF(N3T.EQ.N4.AND.MDT.GT.1)MD4=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.GT.1)MD3=1
      IF(N12T.EQ.N4.AND.MDT.GT.1)MD4=1
      N13T=ISYMP(N13,NT)
      IF(N13T.EQ.N4.AND.MDT.GT.1)MD4=1
      N23T=ISYMP(N23,NT)
      IF(N23T.EQ.N4.AND.MDT.GT.1)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      N123T=ISYMP(N123,NT)
      IF(N123T.EQ.N4.AND.MDT.GT.1)MD4=1
C**LOOP ROUND TAU
      DO MTAU=1,MMTAU/MDT
        ITAU=ITAU+INCTAU
        IF(ITAU.GT.722)ITAU=ITAU-720
        DO K=1,NMODE
          QQ(K)=0
          MULT(K)=0
        END DO
        MULT(MODE1)=1
        MULT(MODE2)=1
        MULT(MODE3)=1
        MULT(MODE4)=1
        DO M1=1,MM1/MD1
          QQ(MODE1)=XQ1(M1)
          DO M2=1,MM2/MD2
            QQ(MODE2)=XQ2(M2)
            DO M3=1,MM3/MD3
              QQ(MODE3)=XQ3(M3)
              DO M4=1,MM4/MD4
                QQ(MODE4)=XQ4(M4)
                CALL MILLER(NMODE,NATOM,QQ,BB,BS,BBS,B,BSS,XX,XL,XXP,
     1          XM,BBT,BT,ZZ,ITAU,1)
                IF(KMAX.GT.0.AND.ICI.LT.0)THEN
C**FULL ROTATION
                  CALL VSET64(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,E0,
     1            J21,M1,M2,M3,M4,S,MX,MY,MZ,0)
C                 IF(JCOUPC.GT.0)THEN
C                   DO J=1,J21
C                     V(J,M4,M3,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
C    1                S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
C                   END DO
C                 ELSE
C                   DO J=1,J21
C                     VR(J,M4,M3,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
C    1                S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
C                   END DO
C                 END IF
                ELSE
C**ADIABATIC ROTATION - BEWARE XMU
                  CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,1)
                  DO I=1,J21
                    DO J=1,J21
                      V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
                    END DO
                  END DO
                  CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,
     1            J21,V0,V0,V0,IDUM,IDUM,IDUM)
                  CALL VSET64(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,E0,
     1            J21,M1,M2,M3,M4,S,MX,MY,MZ,1)
C                 IF(JCOUPC.GT.0)THEN
C                   DO J=1,J21
C                     V(J,M4,M3,M2,M1)=E0(J)
C                   END DO
C                 ELSE
C                   DO J=1,J21
C                     VR(J,M4,M3,M2,M1)=E0(J)
C                   END DO
C                 END IF
                END IF
              END DO
            END DO
          END DO
        END DO
        CALL OUT64(V,VR,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,J21)
C       IF(JCOUPC.GT.0)THEN
C         WRITE(I64)V
C       ELSE
C         WRITE(I64)VR
C       END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE SET61(V,VR,MM,E0,J21,M,S,MX,MY,MZ,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM)
      REAL*4 VR(J21,MM)
      DIMENSION S(J21,9,J21)
      DIMENSION E0(J21)
      COMMON/MOMI/XK(3,3),XMU(3,3)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      IF(ICOUPC.LT.1)RETURN
      IF(IND.EQ.0)THEN
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M)=(S(J,1,J)*XMU(MX,MX)+S(J,2,J)*XMU(MY,MY)+
     1      S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        ELSE
          DO J=1,J21
            VR(J,M)=(S(J,1,J)*XMU(MX,MX)+S(J,2,J)*XMU(MY,MY)+
     1      S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        END IF
      ELSE
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M)=E0(J)
          END DO
        ELSE
          DO J=1,J21
            VR(J,M)=E0(J)
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VSET61(V,VR,MM,E0,J21,M,S,MX,MY,MZ,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM)
      REAL*4 VR(J21,MM)
      DIMENSION S(J21,9,J21)
      DIMENSION E0(J21)
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      IF(ICOUPC.LT.1)RETURN
      IF(IND.EQ.0)THEN
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M)=(S(J,1,J)*XMU0(MX,MX)+S(J,2,J)*XMU0(MY,MY)+
     1      S(J,3,J)*XMU0(MZ,MZ))/2
          END DO
        ELSE
          DO J=1,J21
            VR(J,M)=(S(J,1,J)*XMU0(MX,MX)+S(J,2,J)*XMU0(MY,MY)+
     1      S(J,3,J)*XMU0(MZ,MZ))/2
          END DO
        END IF
      ELSE
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M)=E0(J)
          END DO
        ELSE
          DO J=1,J21
            VR(J,M)=E0(J)
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE OUT61(V,VR,MM,J21)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM)
      REAL*4 VR(J21,MM)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.1)WRITE(I61)V
      ELSE
        IF(ICOUPC.GE.1)WRITE(I61)VR
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE SET62(V,VR,MM1,MM2,E0,J21,M1,M2,S,MX,MY,MZ,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM2,MM1)
      REAL*4 VR(J21,MM2,MM1)
      DIMENSION S(J21,9,J21)
      DIMENSION E0(J21)
      COMMON/MOMI/XK(3,3),XMU(3,3)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      IF(ICOUPC.LT.2)RETURN
      IF(IND.EQ.0)THEN
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1      S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        ELSE
          DO J=1,J21
            VR(J,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1      S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        END IF
      ELSE
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M2,M1)=E0(J)
          END DO
        ELSE
          DO J=1,J21
            VR(J,M2,M1)=E0(J)
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VSET62(V,VR,MM1,MM2,E0,J21,M1,M2,S,MX,MY,MZ,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM2,MM1)
      REAL*4 VR(J21,MM2,MM1)
      DIMENSION S(J21,9,J21)
      DIMENSION E0(J21)
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      IF(ICOUPC.LT.2)RETURN
      IF(IND.EQ.0)THEN
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
     1      S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
          END DO
        ELSE
          DO J=1,J21
            VR(J,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
     1      S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
          END DO
        END IF
      ELSE
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M2,M1)=E0(J)
          END DO
        ELSE
          DO J=1,J21
            VR(J,M2,M1)=E0(J)
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE OUT62(V,VR,MM1,MM2,J21)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM2,MM1)
      REAL*4 VR(J21,MM2,MM1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.2)WRITE(I62)V
      ELSE
        IF(ICOUPC.GE.2)WRITE(I62)VR
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE SET63(V,VR,MM1,MM2,MM3,E0,J21,M1,M2,M3,S,MX,MY,MZ,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM3,MM2,MM1)
      REAL*4 VR(J21,MM3,MM2,MM1)
      DIMENSION S(J21,9,J21)
      DIMENSION E0(J21)
      COMMON/MOMI/XK(3,3),XMU(3,3)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FILASS/IOUT
      IF(ICOUPC.LT.3)RETURN
      IF(IND.EQ.0)THEN
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1      S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        ELSE
          DO J=1,J21
            VR(J,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1      S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        END IF
      ELSE
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M3,M2,M1)=E0(J)
          END DO
        ELSE
          DO J=1,J21
            VR(J,M3,M2,M1)=E0(J)
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VSET63(V,VR,MM1,MM2,MM3,E0,J21,M1,M2,M3,S,MX,MY,MZ,
     1IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM3,MM2,MM1)
      REAL*4 VR(J21,MM3,MM2,MM1)
      DIMENSION S(J21,9,J21)
      DIMENSION E0(J21)
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      IF(ICOUPC.LT.3)RETURN
      IF(IND.EQ.0)THEN
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M3,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
     1      S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
          END DO
        ELSE
          DO J=1,J21
            VR(J,M3,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
     1      S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
          END DO
        END IF
      ELSE
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M3,M2,M1)=E0(J)
          END DO
        ELSE
          DO J=1,J21
            VR(J,M3,M2,M1)=E0(J)
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE OUT63(V,VR,MM1,MM2,MM3,J21)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM3,MM2,MM1)
      REAL*4 VR(J21,MM3,MM2,MM1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.3)WRITE(I63)V
      ELSE
        IF(ICOUPC.GE.3)WRITE(I63)VR
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE SET64(V,VR,MM1,MM2,MM3,MM4,E0,J21,M1,M2,M3,M4,S,MX,MY,
     1MZ,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM4,MM3,MM2,MM1)
      REAL*4 VR(J21,MM4,MM3,MM2,MM1)
      DIMENSION S(J21,9,J21)
      DIMENSION E0(J21)
      COMMON/MOMI/XK(3,3),XMU(3,3)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      IF(ICOUPC.LT.4)RETURN
      IF(IND.EQ.0)THEN
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M4,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1      S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        ELSE
          DO J=1,J21
            VR(J,M4,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1      S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
          END DO
        END IF
      ELSE
        IF(JCOUPC.GE.0)THEN
          DO J=1,J21
            V(J,M4,M3,M2,M1)=E0(J)
          END DO
        ELSE
          DO J=1,J21
            VR(J,M4,M3,M2,M1)=E0(J)
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VSET64(V,VR,MM1,MM2,MM3,MM4,E0,J21,M1,M2,M3,M4,S,MX,
     1MY,MZ,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM4,MM3,MM2,MM1)
      REAL*4 VR(J21,MM4,MM3,MM2,MM1)
      DIMENSION S(J21,9,J21)
      DIMENSION E0(J21)
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      IF(ICOUPC.LT.4)RETURN
      IF(IND.EQ.0)THEN
        IF(JCOUPL.GT.0)THEN
          DO J=1,J21
            V(J,M4,M3,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
     1      S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
          END DO
        ELSE
          DO J=1,J21
            VR(J,M4,M3,M2,M1)=(S(J,1,J)*XMU0(MX,MX)+
     1      S(J,2,J)*XMU0(MY,MY)+S(J,3,J)*XMU0(MZ,MZ))/2
          END DO
        END IF
      ELSE
        IF(JCOUPL.GT.0)THEN
          DO J=1,J21
            V(J,M4,M3,M2,M1)=E0(J)
          END DO
        ELSE
          DO J=1,J21
            VR(J,M4,M3,M2,M1)=E0(J)
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE OUT64(V,VR,MM1,MM2,MM3,MM4,J21)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM4,MM3,MM2,MM1)
      REAL*4 VR(J21,MM4,MM3,MM2,MM1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.4)WRITE(I64)V
      ELSE
        IF(ICOUPC.GE.4)WRITE(I64)VR
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MV0(NMODE,MODE,HTAU,XQTAU,XA,NSIZE,NNTAU,MMTAU,
     1IP,ISIZE,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(12)
      REAL*4 VMR(12)
      DIMENSION MODINT(NMODE)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU),XA(1)
      DIMENSION IP(ISIZE,NMODE)
      DIMENSION TEMP1(2),TEMP(2)
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP

      IF(ITIM.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MV0'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACTC=JNTFAC(NMODE,ICOUPC,0)

      KAL=KROTL/2
      INCTL=MOD(IFLAUD,2)*MOD(KAL,2)
      LMAXL=IFLAUD-(IFLAUD-1)*MOD(KAL+1,2)
      KAR=KROTR/2
      INCTR=MOD(IFLAUD,2)*MOD(KAR,2)
      LMAXR=IFLAUD-(IFLAUD-1)*MOD(KAR+1,2)
      FACTRC=0
      IF(J21.GT.1)FACTRC=IFACTC
      MDT=MODINT(NSMODE)

C***********************************************************

      IOFF=(IABC-7)*2
C**LOOP ROUND TAU
      ITAU=INIT-INCTAU
      DO MTAU=1,MMTAU/MDT
        ITAU=ITAU+INCTAU
CCCC    IF(ITAU.GT.362)ITAU=ITAU-360
        IF(ITAU.GT.722)ITAU=ITAU-720

C***********************************************************

        IF(JCOUPC.GE.0)THEN
          IF(ICOUPC.GE.0)READ(I91)VM
        ELSE
          IF(ICOUPC.GE.0)READ(I91)VMR
        END IF

C***********************************************************

        IF(IABC.LT.7)THEN
          IF(JCOUPC.GE.0)THEN
            TEMP(1)=VM(IABC)*MDT*IFACTC
          ELSE
            TEMP(1)=VMR(IABC)*MDT*IFACTC
          END IF
          TEMP(2)=0
        ELSE
          IF(JCOUPC.GE.0)THEN
            DO K=1,2
              TEMP1(K)=VM(6+K+IOFF)*MDT*IFACTC
            END DO
          ELSE
            DO K=1,2
              TEMP1(K)=VMR(6+K+IOFF)*MDT*IFACTC
            END DO
          END IF
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          IRTAU=IP(IRHS,MODE)+1-MOD(KAR,2)+INCTR
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            ILTAU=IP(ILHS,MODE)+1-MOD(KAL,2)+INCTL
C**FORM REQUIRED COMBINATIONS
            TEMP(1)=TEMP1(1)*HTAU(ILTAU,MTAU,1,LMAXL)+
     1              TEMP1(2)*HTAU(ILTAU,MTAU,2,LMAXL)
            TEMP(2)=-TEMP1(2)*HTAU(ILTAU,MTAU,1,LMAXL)
            DO K=1,2
              XA(ILHS+J0)=XA(ILHS+J0)+TEMP(K)*HTAU(IRTAU,MTAU,K,LMAXR)*
     1        DSTAU(ITAU)
            END DO
          END DO
        END DO
      END DO
      IF(ITIM.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMI0(XA,ISIZEL,ISIZER,IABC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V,VM(6)
      REAL*4 VR,VMR(6)
      DIMENSION XA(ISIZEL,ISIZER)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MUREF/V,VM,VR,VMR
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6),FACTS(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP

      FACT=1
      DO I=1,ICOUPC
        FACT=FACT-FACTOR(I)*FACTS(I)
      END DO

      DO IRHS=1,ISIZER
        ILHS=IRHS
        IF(JCOUPC.GE.0)
     1  XA(ILHS,IRHS)=XA(ILHS,IRHS)+VM(IABC)*FACT*IHERM
        IF(JCOUPC.LT.0)
     1  XA(ILHS,IRHS)=XA(ILHS,IRHS)+VMR(IABC)*FACT*IHERM
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMV0(NMODE,XA,XK,NSIZE,IPL,IPR,ISIZMX,ISIZEL,ISIZER,
     1IP1,ISIZE1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE)
      DIMENSION XA(ISIZEL,ISIZER),XK(1),IP1(ISIZE1,1)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMV0'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NRTAU=IPR(IRHS,NMODE)
C**FIND RHS INDEX (TRIVIAL)
        DO IR=1,NSIZE
          IF(NRTAU.EQ.IP1(IR,1))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.NMODE.AND.(IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NLTAU=IPL(ILHS,NMODE)
C**FIND LHS INDEX (TRIVIAL)
          DO IL=1,NSIZE
            IF(NLTAU.EQ.IP1(IL,1))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MI1(NMODE,MODE,MOD1,HL,HR,XQ,XA,NSIZE,NN,MM,IP,
     1ISIZE,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM,7)
      REAL*4 VMR(MM,7)
      DIMENSION MODINT(NMODE)
      DIMENSION HL(NN,MM,3),HR(NN,MM,3),XQ(MM),XA(1)
      DIMENSION IP(ISIZE,NMODE)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP
      MD=MODINT(MOD1)
      CALL VDMI1(NMODE,MODE,MOD1,HL,HR,XQ,XA,NSIZE,NN,MM,MM/MD,IP,
     1ISIZE,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDMI1(NMODE,MODE,MOD1,HL,HR,XQ,XA,NSIZE,NN,MH,MM,IP,
     1ISIZE,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM,7)
      REAL*4 VMR(MM,7)
      DIMENSION MODINT(NMODE)
      DIMENSION HL(NN,MH,3),HR(NN,MH,3),XQ(MM),XA(1)
      DIMENSION IP(ISIZE,NMODE)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP

      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MI1'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACTC=INTFAC(NMODE,ICOUPC,1)

C**TEMPORARY-LIN
C     LMAXL=1
C     LMAXR=1
C     IF(MOD1.GT.NONLIN)THEN
C       KAL=KROTL/2
C       LMAXL=1+MOD(KAL,2)
C       KAR=KROTR/2
C       LMAXR=1+MOD(KAR,2)
C     END IF
C**TEMPORARY-LIN

C***********************************************************

      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GT.0)READ(I91)VM
      ELSE
        IF(ICOUPC.GT.0)READ(I91)VMR
      END IF

C***********************************************************

      MD=MODINT(MOD1)

C***********************************************************
CCCC  DO M=1,MM/MD
      DO M=1,MM
        IF(JCOUPC.GE.0)THEN
          TERM=VM(M,IABC)*IFACTC
        ELSE
          TERM=VMR(M,IABC)*IFACTC
        END IF
C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          NR=IP(IRHS,MODE)
          X=TERM*HR(NR,M,1)*MD
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP(ILHS,MODE)
            Y=HL(NL,M,1)
            XA(ILHS+J0)=XA(ILHS+J0)+Y*X
          END DO
        END DO
      END DO
      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM1A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MV1(NMODE,MODE1,MODE2,MOD1,H1,XQ1,HTAU,XQTAU,NN1,
     1MM1,NNTAU,MMTAU,XA,NSIZE,IP,ISIZE,TEMP,JCI1,JCIM,X0,T0,VM,VMR,
     2J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      REAL*8 VM(MM1,15)
      REAL*4 VMR(MM1,15)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),HTAU(NNTAU,MMTAU,3,1)
      DIMENSION IP(ISIZE,NMODE)
      DIMENSION TEMP(5,JCI1,JCI1)
      DIMENSION X0(5,JCI1,JCI1,MM1),T0(5,JCIM,JCIM,MMTAU)
      DIMENSION X(5),Y(5),C(5)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQTAU(MMTAU)
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      CALL VDMV1(NMODE,MODE1,MODE2,MOD1,H1,XQ1,HTAU,XQTAU,NN1,MM1,
     1MM1/MD1,NNTAU,MMTAU,XA,NSIZE,IP,ISIZE,TEMP,JCI1,JCIM,X0,T0,VM,
     2VMR,J21,KROTL,KROTR,IABC,MODINT)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDMV1(NMODE,MODE1,MODE2,MOD1,H1,XQ1,HTAU,XQTAU,NN1,
     1MH1,MM1,NNTAU,MMTAU,XA,NSIZE,IP,ISIZE,TEMP,JCI1,JCIM,X0,T0,VM,
     2VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      REAL*8 VM(MM1,15)
      REAL*4 VMR(MM1,15)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MH1,3),HTAU(NNTAU,MMTAU,3,1)
      DIMENSION IP(ISIZE,NMODE)
      DIMENSION TEMP(5,JCI1,JCI1)
      DIMENSION X0(5,JCI1,JCI1,MM1),T0(5,JCIM,JCIM,MMTAU)
      DIMENSION X(5),Y(5),C(5)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQTAU(MMTAU)
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP

      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MV1'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACTC=JNTFAC(NMODE,ICOUPC,1)

      KAL=KROTL/2
      INCTL=MOD(IFLAUD,2)*MOD(KAL,2)
      LMAXL=IFLAUD-(IFLAUD-1)*MOD(KAL+1,2)
      KAR=KROTR/2
      INCTR=MOD(IFLAUD,2)*MOD(KAR,2)
      LMAXR=IFLAUD-(IFLAUD-1)*MOD(KAR+1,2)
      FACTRC=0
      IF(J21.GT.1)FACTRC=IFACTC
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      MD=MD1*MDT

C***********************************************************

      IOFF=(IABC-7)*3
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,1)*MD*IFACTC
          X(2)=X(1)
          X(3)=H1(NR1,M1,2)*MD*IFACTC
          X(4)=X(1)
          X(5)=X(1)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,1)
            Y(2)=H1(NL1,M1,2)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            DO K=1,5
              X0(K,NL1,NR1,M1)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
      DO MTAU=1,MMTAU/MDT
        DO NRR=1,JCIM
          NR=NRR+1-MOD(KAR,2)+INCTR
          X(1)=HTAU(NR,MTAU,1,LMAXR)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=HTAU(NR,MTAU,2,LMAXR)
          DO NLL=1,JCIM
            NL=NLL+1-MOD(KAL,2)+INCTL
            Y(1)=HTAU(NL,MTAU,1,LMAXL)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=HTAU(NL,MTAU,2,LMAXL)
            Y(5)=Y(1)
            DO K=1,5
              T0(K,NLL,NRR,MTAU)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**LOOP ROUND TAU (START 2-MODE INTEGRATION)
      ITAU=INIT-INCTAU
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          ITAU=ITAU+INCTAU
CCCC      IF(ITAU.GT.362)ITAU=ITAU-360
          IF(ITAU.GT.722)ITAU=ITAU-720
        ELSE
CCCC      QTAU=XQTAU(MTAU)
C**DELTA4 AND DELTA5 FROM TORSION (ARBITRARY SET EULER 'GAMMA' TO ZERO)
CCCC      DELTA(4)=+QTAU
CCCC      DELTA(5)=-QTAU
        END IF

C***********************************************************

        IF(JCOUPC.GE.0)THEN
          IF(ICOUPC.GE.1)READ(I91)VM
        ELSE
          IF(ICOUPC.GE.1)READ(I91)VMR
        END IF

C***********************************************************

        DO IRHS1=1,JCI1
          DO ILHS1=1,JCI1
            DO K=1,5
              TEMP(K,ILHS1,IRHS1)=0
            END DO
          END DO
        END DO

        IF(IABC.LT.7)THEN
C**START 1-MODE INTEGRATION
          IF(JCOUPL.GT.0)THEN
CCCC        DO M1=1,MM1/MD1
            DO M1=1,MM1
              DO IRHS1=1,JCI1
                DO ILHS1=1,JCI1
                  TEMP(1,ILHS1,IRHS1)=TEMP(1,ILHS1,IRHS1)+
     1            X0(1,ILHS1,IRHS1,M1)*VM(M1,IABC)
                END DO
              END DO
            END DO
          ELSE
CCCC        DO M1=1,MM1/MD1
            DO M1=1,MM1
              DO IRHS1=1,JCI1
                DO ILHS1=1,JCI1
                  TEMP(1,ILHS1,IRHS1)=TEMP(1,ILHS1,IRHS1)+
     1            X0(1,ILHS1,IRHS1,M1)*VMR(M1,IABC)
                END DO
              END DO
            END DO
          END IF
C**END 1-MODE INTEGRATION
        ELSE
C**START 1-MODE INTEGRATION
CCCC      DO M1=1,MM1/MD1
          DO M1=1,MM1
            IF(JCOUPL.GT.0)THEN
              C(1)=VM(M1,7+IOFF)
              C(2)=VM(M1,8+IOFF)
              C(3)=-VM(M1,8+IOFF)
              C(4)=VM(M1,9+IOFF)
              C(5)=-VM(M1,9+IOFF)
            ELSE
              C(1)=VMR(M1,7+IOFF)
              C(2)=VMR(M1,8+IOFF)
              C(3)=-VMR(M1,8+IOFF)
              C(4)=VMR(M1,9+IOFF)
              C(5)=-VMR(M1,9+IOFF)
            END IF
            DO IRHS1=1,JCI1
              DO ILHS1=1,JCI1
                DO K=1,5
                  TEMP(K,ILHS1,IRHS1)=TEMP(K,ILHS1,IRHS1)+
     1            X0(K,ILHS1,IRHS1,M1)*C(K)
                END DO
              END DO
            END DO
          END DO
        END IF
C**END 1-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
        DO IRHS=1,NSIZE
          NR=IP(IRHS,MODE1)
          IRTAU=IP(IRHS,MODE2)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP(ILHS,MODE1)
            ILTAU=IP(ILHS,MODE2)
            DO K=1,5
              XA(ILHS+J0)=XA(ILHS+J0)+TEMP(K,NL,NR)*
     1        T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
            END DO
          END DO
        END DO
      END DO
C**END 2-MODE INTEGRATION
      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM1A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMI1(NMODE,MODE,XA,XK,NSIZE,IPL,IPR,ISIZMX,ISIZEL,
     1ISIZER,IP1,ISIZE1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(ISIZEL,ISIZER),XK(1)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE),IP1(ISIZE1,1)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMI1'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR=IPR(IRHS,MODE)
C**FIND RHS INDEX (TRIVIAL CASE)
        DO IR=1,NSIZE
          IF(NR.EQ.IP1(IR,1))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE.AND.(IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL=IPL(ILHS,MODE)
C**FIND LHS INDEX (TRIVIAL CASE)
          DO IL=1,NSIZE
            IF(NL.EQ.IP1(IL,1))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM1A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMV1(NMODE,MODE,XA,XK,NSIZE,IPL,IPR,ISIZMX,ISIZEL,
     1ISIZER,IP2,ISIZE2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE)
      DIMENSION XA(ISIZEL,ISIZER),XK(1),IP2(ISIZE2,2)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMV1'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR=IPR(IRHS,MODE)
        NRTAU=IPR(IRHS,NMODE)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR.EQ.IP2(IR,1).AND.NRTAU.EQ.IP2(IR,2))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE.AND.K.NE.NMODE.AND.(IPR(IRHS,K).NE.
     1      IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL=IPL(ILHS,MODE)
          NLTAU=IPL(ILHS,NMODE)
C**FIND LHS INDEX
          DO IL=1,NSIZE
            IF(NL.EQ.IP2(IL,1).AND.NLTAU.EQ.IP2(IL,2))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM1A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MI2(NMODE,MODE1,MODE2,MOD1,MOD2,H1L,H1R,XQ1,H2L,H2R,
     1XQ2,NN1,MM1,NN2,MM2,XA,NSIZE,IP,ISIZE,TEMP,TEMP1,JCI1,JCI2,X0,Y0,
     2VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM2,MM1,13)
      REAL*4 VMR(MM2,MM1,13)
      DIMENSION MODINT(NMODE)
      DIMENSION H1L(NN1,MM1,3),H2L(NN2,MM2,3)
      DIMENSION H1R(NN1,MM1,3),H2R(NN2,MM2,3)
      DIMENSION IP(ISIZE,NMODE)
C     DIMENSION TEMP1(4,JCI2,JCI2)
      DIMENSION TEMP(5,JCI2,JCI2)
      DIMENSION X0(5,JCI1,JCI1,MM1),Y0(5,JCI2,JCI2,MM2)
      DIMENSION X(2),Y(2)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      CALL VDMI2(NMODE,MODE1,MODE2,MOD1,MOD2,H1L,H1R,XQ1,H2L,H2R,XQ2,
     1NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,XA,NSIZE,IP,ISIZE,TEMP,TEMP1,
     2JCI1,JCI2,X0,Y0,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDMI2(NMODE,MODE1,MODE2,MOD1,MOD2,H1L,H1R,XQ1,H2L,H2R,
     1XQ2,NN1,MH1,MM1,NN2,MH2,MM2,XA,NSIZE,IP,ISIZE,TEMP,TEMP1,JCI1,
     2JCI2,X0,Y0,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM2,MM1,13)
      REAL*4 VMR(MM2,MM1,13)
      DIMENSION MODINT(NMODE)
      DIMENSION H1L(NN1,MH1,3),H2L(NN2,MH2,3)
      DIMENSION H1R(NN1,MH1,3),H2R(NN2,MH2,3)
      DIMENSION IP(ISIZE,NMODE)
C     DIMENSION TEMP1(4,JCI2,JCI2)
      DIMENSION TEMP(5,JCI2,JCI2)
      DIMENSION X0(5,JCI1,JCI1,MM1),Y0(5,JCI2,JCI2,MM2)
      DIMENSION X(2),Y(2)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP

      IF(ITIM2A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MI2'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACTC=INTFAC(NMODE,ICOUPC,2)

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(ICOUPC.GE.2)READ(I92)VM
      ELSE
        IF(ICOUPC.GE.2)READ(I92)VMR
      END IF

C**TEMPORARY-LIN
      KAL=KROTL/2
      KAR=KROTR/2
      KKAR=MOD(KAR,2)
      KKAL=MOD(KAL,2)
      LMAXL1=1
      LMAXR1=1
      IF(MOD1.GT.NONLIN)THEN
        LMAXL1=1+MOD(KAL,2)
        LMAXR1=1+MOD(KAR,2)
      END IF
      IFACT1=IABS(LMAXL1-LMAXR1)
      LMAXL2=1
      LMAXR2=1
      IF(MOD2.GT.NONLIN)THEN
        LMAXL2=1+MOD(KAL,2)
        LMAXR2=1+MOD(KAR,2)
      END IF
      IFACT2=IABS(LMAXL2-LMAXR2)
C**TEMPORARY-LIN

C***********************************************************

      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      MD=MD1*MD2

C***********************************************************

      IOFF=(IABC-7)*2
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          F=H1R(NR1,M1,1)*MD
          DO NL1=1,JCI1
            X0(1,NL1,NR1,M1)=F*H1L(NL1,M1,1)
          END DO
C**DIFFERENTIATE TO THE RIGHT
          DO K=1,2
            X(K)=H1R(NR1,M1,3-K)*MD
          END DO
          DO NL1=1,JCI1
            F=H1L(NL1,M1,1)
            DO K=1,2
              X0(K+1,NL1,NR1,M1)=F*X(K)
            END DO
          END DO
C**DIFFERENTIATE TO THE LEFT
          F=H1R(NR1,M1,1)*MD
          DO NL1=1,JCI1
            DO K=1,2
              X0(K+3,NL1,NR1,M1)=F*H1L(NL1,M1,3-K)
            END DO
          END DO
        END DO
      END DO

CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          F=H2R(NR2,M2,1)*IFACTC
          DO NL2=1,JCI2
            Y0(1,NL2,NR2,M2)=F*H2L(NL2,M2,1)
          END DO
C**DIFFERENTIATE TO THE RIGHT
          DO K=1,2
            Y(K)=H2R(NR2,M2,K)*IFACTC
          END DO
          DO NL2=1,JCI2
            F=H2L(NL2,M2,1)
            DO K=1,2
              Y0(K+1,NL2,NR2,M2)=F*Y(K)
            END DO
          END DO
C**DIFFERENTIATE TO THE LEFT
          F=H2R(NR2,M2,1)*IFACTC
          DO NL2=1,JCI2
            DO K=1,2
              Y0(K+3,NL2,NR2,M2)=F*H2L(NL2,M2,K)
            END DO
          END DO
        END DO
      END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**START 2-MODE INTEGRATION
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO IRHS2=1,JCI2
          DO ILHS2=1,JCI2
            DO K=1,5
              TEMP(K,ILHS2,IRHS2)=0
            END DO
          END DO
        END DO
        IF(IABC.LT.7)THEN
C**START 1-MODE INTEGRATION
          IF(JCOUPL.GT.0)THEN
CCCC        DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS2=1,JCI2
                DO ILHS2=1,JCI2
                  TEMP(1,ILHS2,IRHS2)=TEMP(1,ILHS2,IRHS2)+
     1            Y0(1,ILHS2,IRHS2,M2)*VM(M2,M1,IABC)
                END DO
              END DO
            END DO
          ELSE
CCCC        DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS2=1,JCI2
                DO ILHS2=1,JCI2
                  TEMP(1,ILHS2,IRHS2)=TEMP(1,ILHS2,IRHS2)+
     1            Y0(1,ILHS2,IRHS2,M2)*VMR(M2,M1,IABC)
                END DO
              END DO
            END DO
          END IF
C**END 1-MODE INTEGRATION
        ELSE
C**START 1-MODE INTEGRATION
          IF(JCOUPL.GT.0)THEN
CCCC        DO M2=1,MM2/MD2
            DO M2=1,MM2
              ZZ=VM(M2,M1,13)
              DO IRHS2=1,JCI2
                DO ILHS2=1,JCI2
C**DIFFERENTIATE TO RIGHT
                  DO K=1,2
                    TEMP(K+1,ILHS2,IRHS2)=TEMP(K+1,ILHS2,IRHS2)+
     1              Y0(K+1,ILHS2,IRHS2,M2)*VM(M2,M1,6+K+IOFF)
                  END DO
C**DIFFERENTIATE TO LEFT
                  DO K=1,2
                    TEMP(K+3,ILHS2,IRHS2)=TEMP(K+3,ILHS2,IRHS2)-
     1              Y0(K+3,ILHS2,IRHS2,M2)*VM(M2,M1,6+K+IOFF)
                  END DO
C**EXTRA TERM FOR LINEAR
                  TEMP(1,ILHS2,IRHS2)=TEMP(1,ILHS2,IRHS2)+
     1            ZZ*Y0(1,ILHS2,IRHS2,M2)*KKAR*(KKAR-KKAL)*
     2            (IFACT1+IFACT2)
                END DO
              END DO
            END DO
          ELSE
CCCC        DO M2=1,MM2/MD2
            DO M2=1,MM2
              ZZ=VMR(M2,M1,13)
              DO IRHS2=1,JCI2
                DO ILHS2=1,JCI2
C**DIFFERENTIATE TO RIGHT
                  DO K=1,2
                    TEMP(K+1,ILHS2,IRHS2)=TEMP(K+1,ILHS2,IRHS2)+
     1              Y0(K+1,ILHS2,IRHS2,M2)*VMR(M2,M1,6+K+IOFF)
                  END DO
C**DIFFERENTIATE TO LEFT
                  DO K=1,2
                    TEMP(K+3,ILHS2,IRHS2)=TEMP(K+3,ILHS2,IRHS2)-
     1              Y0(K+3,ILHS2,IRHS2,M2)*VMR(M2,M1,6+K+IOFF)
                  END DO
C**EXTRA TERM FOR LINEAR
                  TEMP(1,ILHS2,IRHS2)=TEMP(1,ILHS2,IRHS2)+
     1            ZZ*Y0(1,ILHS2,IRHS2,M2)*KKAR*(KKAR-KKAL)*
     2            (IFACT1+IFACT2)
                END DO
              END DO
            END DO
          END IF
C**END 1-MODE INTEGRATION
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
        DO IRHS=1,NSIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            DO K=1,5
              XA(ILHS+J0)=XA(ILHS+J0)+
     1        TEMP(K,NL2,NR2)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO

C**END 2-MODE INTEGRATION
      IF(ITIM2A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM2A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MV2(NMODE,MODE1,MODE2,MODE3,MOD1,MOD2,H1,XQ1,H2,XQ2,
     1HTAU,XQTAU,NN1,MM1,NN2,MM2,NNTAU,MMTAU,XA,NSIZE,IP,ISIZE,TEMP,
     2TEMP1,JCI1,JCI2,JCIM,X0,Y0,T0,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      REAL*8 VM(MM2,MM1,18)
      REAL*4 VMR(MM2,MM1,18)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),HTAU(NNTAU,MMTAU,3,1)
      DIMENSION IP(ISIZE,NMODE)
      DIMENSION TEMP1(7,JCI2,JCI2)
      DIMENSION TEMP(7,JCI1,JCI2,JCI1,JCI2)
      DIMENSION T0(7,JCIM,JCIM,MMTAU),X0(7,JCI1,JCI1,MM1)
      DIMENSION Y0(7,JCI2,JCI2,MM2)
      DIMENSION X(7),Y(7),C(7)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQTAU(MMTAU)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N2.EQ.NT.AND.MDT.GT.1)MD2=1
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.GT.1)MD2=1
      CALL VDMV2(NMODE,MODE1,MODE2,MODE3,MOD1,MOD2,H1,XQ1,H2,XQ2,
     1HTAU,XQTAU,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,NNTAU,MMTAU,XA,NSIZE,
     2IP,ISIZE,TEMP,TEMP1,JCI1,JCI2,JCIM,X0,Y0,T0,VM,VMR,J21,KROTL,
     3KROTR,IABC,MODINT)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDMV2(NMODE,MODE1,MODE2,MODE3,MOD1,MOD2,H1,XQ1,H2,XQ2,
     1HTAU,XQTAU,NN1,MH1,MM1,NN2,MH2,MM2,NNTAU,MMTAU,XA,NSIZE,IP,ISIZE,
     2TEMP,TEMP1,JCI1,JCI2,JCIM,X0,Y0,T0,VM,VMR,J21,KROTL,KROTR,IABC,
     3MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      REAL*8 VM(MM2,MM1,18)
      REAL*4 VMR(MM2,MM1,18)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),HTAU(NNTAU,MMTAU,3,1)
      DIMENSION IP(ISIZE,NMODE)
      DIMENSION TEMP1(7,JCI2,JCI2)
      DIMENSION TEMP(7,JCI1,JCI2,JCI1,JCI2)
      DIMENSION T0(7,JCIM,JCIM,MMTAU),X0(7,JCI1,JCI1,MM1)
      DIMENSION Y0(7,JCI2,JCI2,MM2)
      DIMENSION X(7),Y(7),C(7)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQTAU(MMTAU)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP

      IF(ITIM2A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MV2'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACTC=JNTFAC(NMODE,ICOUPC,2)

      KAL=KROTL/2
      INCTL=MOD(IFLAUD,2)*MOD(KAL,2)
      LMAXL=IFLAUD-(IFLAUD-1)*MOD(KAL+1,2)
      KAR=KROTR/2
      INCTR=MOD(IFLAUD,2)*MOD(KAR,2)
      LMAXR=IFLAUD-(IFLAUD-1)*MOD(KAR+1,2)
      FACTRC=0
      IF(J21.GT.1)FACTRC=IFACTC
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N2.EQ.NT.AND.MDT.GT.1)MD2=1
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.GT.1)MD2=1
      MD=MD1*MD2*MDT

C***********************************************************

      IOFF=(IABC-7)*4
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,1)*MD*IFACTC
          X(2)=X(1)
          X(3)=H1(NR1,M1,2)*MD*IFACTC
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,1)
            Y(2)=H1(NL1,M1,2)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            DO K=1,7
              X0(K,NL1,NR1,M1)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=H2(NR2,M2,2)
          X(6)=X(1)
          X(7)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=H2(NL2,M2,2)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            DO K=1,7
              Y0(K,NL2,NR2,M2)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
      DO MTAU=1,MMTAU/MDT
        DO NRR=1,JCIM
          NR=NRR+1-MOD(KAR,2)+INCTR
          X(1)=HTAU(NR,MTAU,1,LMAXR)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=HTAU(NR,MTAU,2,LMAXR)
          DO NLL=1,JCIM
            NL=NLL+1-MOD(KAL,2)+INCTL
            Y(1)=HTAU(NL,MTAU,1,LMAXL)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=HTAU(NL,MTAU,2,LMAXL)
            Y(7)=Y(1)
            DO K=1,7
              T0(K,NLL,NRR,MTAU)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**LOOP ROUND TAU (START 3-MODE INTEGRATION)
      ITAU=INIT-INCTAU
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          ITAU=ITAU+INCTAU
CCCC      IF(ITAU.GT.362)ITAU=ITAU-360
          IF(ITAU.GT.722)ITAU=ITAU-720
        ELSE
CCCC      QTAU=XQTAU(MTAU)
C**DELTA4 AND DELTA5 FROM TORSION (ARBITRARY SET EULER 'GAMMA' TO ZERO)
CCCC      DELTA(4)=+QTAU
CCCC      DELTA(5)=-QTAU
        END IF

C***********************************************************

        IF(JCOUPC.GE.0)THEN
          IF(ICOUPC.GE.2)READ(I92)VM
        ELSE
          IF(ICOUPC.GE.2)READ(I92)VMR
        END IF

C***********************************************************

        DO IRHS2=1,JCI2
          DO IRHS1=1,JCI1
            DO ILHS2=1,JCI2
              DO ILHS1=1,JCI1
                DO K=1,7
                  TEMP(K,ILHS1,ILHS2,IRHS1,IRHS2)=0
                END DO
              END DO
            END DO
          END DO
        END DO

        IF(IABC.LT.7)THEN
C**START 2-MODE INTEGRATION
CCCC      DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS2=1,JCI2
              DO ILHS2=1,JCI2
                TEMP1(1,ILHS2,IRHS2)=0
              END DO
            END DO
C**START 1-MODE INTEGRATION
            IF(JCOUPL.GT.0)THEN
CCCC          DO M2=1,MM2/MD2
              DO M2=1,MM2
                DO IRHS2=1,JCI2
                  DO ILHS2=1,JCI2
                    TEMP1(1,ILHS2,IRHS2)=TEMP1(1,ILHS2,IRHS2)+
     1              Y0(1,ILHS2,IRHS2,M2)*VM(M2,M1,IABC)
                  END DO
                END DO
              END DO
            ELSE
CCCC          DO M2=1,MM2/MD2
              DO M2=1,MM2
                DO IRHS2=1,JCI2
                  DO ILHS2=1,JCI2
                    TEMP1(1,ILHS2,IRHS2)=TEMP1(1,ILHS2,IRHS2)+
     1              Y0(1,ILHS2,IRHS2,M2)*VMR(M2,M1,IABC)
                  END DO
                END DO
              END DO
            END IF
C**END 1-MODE INTEGRATION
              DO IRHS2=1,JCI2
                DO IRHS1=1,JCI1
                  DO ILHS2=1,JCI2
                    DO ILHS1=1,JCI1
                    TEMP(1,ILHS1,ILHS2,IRHS1,IRHS2)=
     1              TEMP(1,ILHS1,ILHS2,IRHS1,IRHS2)+
     2              X0(1,ILHS1,IRHS1,M1)*TEMP1(1,ILHS2,IRHS2)
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 2-MODE INTEGRATION
        ELSE
C**START 2-MODE INTEGRATION
CCCC      DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS2=1,JCI2
              DO ILHS2=1,JCI2
                DO K=1,7
                  TEMP1(K,ILHS2,IRHS2)=0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
CCCC        DO M2=1,MM2/MD2
            DO M2=1,MM2
              IF(JCOUPL.GT.0)THEN
                C(1)=VM(M2,M1,7+IOFF)
                C(2)=VM(M2,M1,8+IOFF)
                C(3)=-VM(M2,M1,8+IOFF)
                C(4)=VM(M2,M1,9+IOFF)
                C(5)=-VM(M2,M1,9+IOFF)
                C(6)=VM(M2,M1,10+IOFF)
                C(7)=-VM(M2,M1,10+IOFF)
              ELSE
                C(1)=VMR(M2,M1,7+IOFF)
                C(2)=VMR(M2,M1,8+IOFF)
                C(3)=-VMR(M2,M1,8+IOFF)
                C(4)=VMR(M2,M1,9+IOFF)
                C(5)=-VMR(M2,M1,9+IOFF)
                C(6)=VMR(M2,M1,10+IOFF)
                C(7)=-VMR(M2,M1,10+IOFF)
              END IF
              DO IRHS2=1,JCI2
                DO ILHS2=1,JCI2
                  DO K=1,7
                    TEMP1(K,ILHS2,IRHS2)=TEMP1(K,ILHS2,IRHS2)+
     1              Y0(K,ILHS2,IRHS2,M2)*C(K)
                  END DO
                END DO
              END DO
            END DO
C**END 1-MODE INTEGRATION
            DO IRHS2=1,JCI2
              DO IRHS1=1,JCI1
                DO ILHS2=1,JCI2
                  DO ILHS1=1,JCI1
                    DO K=1,7
                      TEMP(K,ILHS1,ILHS2,IRHS1,IRHS2)=
     1                TEMP(K,ILHS1,ILHS2,IRHS1,IRHS2)+
     2                X0(K,ILHS1,IRHS1,M1)*TEMP1(K,ILHS2,IRHS2)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 2-MODE INTEGRATION
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
        DO IRHS=1,NSIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          IRTAU=IP(IRHS,MODE3)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            ILTAU=IP(ILHS,MODE3)
            DO K=1,7
              XA(ILHS+J0)=XA(ILHS+J0)+TEMP(K,NL1,NL2,NR1,NR2)*
     1        T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
            END DO
          END DO
        END DO
      END DO
C**END 3-MODE INTEGRATION
      IF(ITIM2A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM2A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMI2(NMODE,MODE1,MODE2,XA,XK,NSIZE,IPL,IPR,ISIZMX,
     1ISIZEL,ISIZER,IP2,ISIZE2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(ISIZEL,ISIZER),XK(1),
     1IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE),IP2(ISIZE2,2)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM2A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMI2'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IPR(IRHS,MODE1)
        NR2=IPR(IRHS,MODE2)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.(IPR(IRHS,K).NE.
     1      IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPL(ILHS,MODE1)
          NL2=IPL(ILHS,MODE2)
C**FIND LHS INDEX
          DO IL=1,NSIZE
            IF(NL1.EQ.IP2(IL,1).AND.NL2.EQ.IP2(IL,2))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM2A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM2A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMV2(NMODE,MODE1,MODE2,XA,XK,NSIZE,IPL,IPR,ISIZMX,
     1ISIZEL,ISIZER,IP3,ISIZE3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE)
      DIMENSION XA(ISIZEL,ISIZER),XK(1),IP3(ISIZE3,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM2A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMV2'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IPR(IRHS,MODE1)
        NR2=IPR(IRHS,MODE2)
        NRTAU=IPR(IRHS,NMODE)
C**FIND RHS INDEX
C*************************************************************
C**
C**                     NOTE*NOTE*NOTE
C**
C**FOR STRETCH-BEND K WILL DENOTE BEND, L WILL DENOTE STRETCH
C**            K IS INDEX 1        L IS INDEX 2
C*************************************************************
        DO IR=1,NSIZE
          IF(NR1.EQ.IP3(IR,1).AND.NR2.EQ.IP3(IR,2).AND.
     1       NRTAU.EQ.IP3(IR,3))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.NMODE.AND.(
     1      IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPL(ILHS,MODE1)
          NL2=IPL(ILHS,MODE2)
          NLTAU=IPL(ILHS,NMODE)
C**FIND LHS INDEX
          DO IL=1,NSIZE
            IF(NL1.EQ.IP3(IL,1).AND.NL2.EQ.IP3(IL,2).AND.
     1         NLTAU.EQ.IP3(IL,3))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM2A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM2A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MI3(NMODE,MODE1,MODE2,MODE3,MOD1,MOD2,MOD3,H1L,H1R,
     1XQ1,H2L,H2R,XQ2,H3L,H3R,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,XA,NSIZE,IP,
     2ISIZE,TEMP,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,Z0,VM,VMR,J21,KROTL,
     3KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM3,MM2,MM1,16)
      REAL*4 VMR(MM3,MM2,MM1,16)
      DIMENSION MODINT(NMODE)
      DIMENSION H1L(NN1,MM1,3),H2L(NN2,MM2,3),H3L(NN3,MM3,3)
      DIMENSION H1R(NN1,MM1,3),H2R(NN2,MM2,3),H3R(NN3,MM3,3)
      DIMENSION IP(ISIZE,NMODE)
      DIMENSION TEMP1(7,JCI3,JCI3)
C     DIMENSION TEMP2(3,JCI2,JCI3,JCI2,JCI3)
      DIMENSION TEMP(7,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(7,JCI1,JCI1,MM1),Y0(7,JCI2,JCI2,MM2)
      DIMENSION Z0(7,JCI3,JCI3,MM3)
      DIMENSION X(2),Y(2),Z(2)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
        END DO
      END DO
      IF(N2.EQ.N3)MD3=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N2)MD2=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      CALL VDMI3(NMODE,MODE1,MODE2,MODE3,MOD1,MOD2,MOD3,H1L,H1R,XQ1,
     1H2L,H2R,XQ2,H3L,H3R,XQ3,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,
     2MM3/MD3,XA,NSIZE,IP,ISIZE,TEMP,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,
     3Z0,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDMI3(NMODE,MODE1,MODE2,MODE3,MOD1,MOD2,MOD3,H1L,H1R,
     1XQ1,H2L,H2R,XQ2,H3L,H3R,XQ3,NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,
     2XA,NSIZE,IP,ISIZE,TEMP,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,Z0,VM,
     3VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM3,MM2,MM1,16)
      REAL*4 VMR(MM3,MM2,MM1,16)
      DIMENSION MODINT(NMODE)
      DIMENSION H1L(NN1,MH1,3),H2L(NN2,MH2,3),H3L(NN3,MH3,3)
      DIMENSION H1R(NN1,MH1,3),H2R(NN2,MH2,3),H3R(NN3,MH3,3)
      DIMENSION IP(ISIZE,NMODE)
      DIMENSION TEMP1(7,JCI3,JCI3)
C     DIMENSION TEMP2(3,JCI2,JCI3,JCI2,JCI3)
      DIMENSION TEMP(7,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(7,JCI1,JCI1,MM1),Y0(7,JCI2,JCI2,MM2)
      DIMENSION Z0(7,JCI3,JCI3,MM3)
      DIMENSION X(2),Y(2),Z(2)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP

      IF(ITIM3A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MI3'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACTC=INTFAC(NMODE,ICOUPC,3)

C**TEMPORARY-LIN
      KAL=KROTL/2
      KAR=KROTR/2
      KKAR=MOD(KAR,2)
      KKAL=MOD(KAL,2)
      LMAXL1=1
      LMAXR1=1
      IF(MOD1.GT.NONLIN)THEN
        LMAXL1=1+MOD(KAL,2)
        LMAXR1=1+MOD(KAR,2)
      END IF
      IFACT1=IABS(LMAXL1-LMAXR1)
      LMAXL2=1
      LMAXR2=1
      IF(MOD2.GT.NONLIN)THEN
        LMAXL2=1+MOD(KAL,2)
        LMAXR2=1+MOD(KAR,2)
      END IF
      IFACT2=IABS(LMAXL2-LMAXR2)
      LMAXL3=1
      LMAXR3=1
      IF(MOD3.GT.NONLIN)THEN
        LMAXL3=1+MOD(KAL,2)
        LMAXR3=1+MOD(KAR,2)
      END IF
      IFACT3=IABS(LMAXL3-LMAXR3)
C**TEMPORARY-LIN

C***********************************************************

      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.3)READ(I93)VM
      ELSE
        IF(ICOUPC.GE.3)READ(I93)VMR
      END IF

C***********************************************************

      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
        END DO
      END DO
      IF(N2.EQ.N3)MD3=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N2)MD2=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      MD=MD1*MD2*MD3

C***********************************************************

      IOFF=(IABC-7)*3
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          F=H1R(NR1,M1,1)*MD
          DO NL1=1,JCI1
            X0(1,NL1,NR1,M1)=F*H1L(NL1,M1,1)
          END DO
C**DIFFERENTIATE TO THE RIGHT
          DO K=1,2
            X(K)=H1R(NR1,M1,3-K)*MD
          END DO
          DO NL1=1,JCI1
            F=H1L(NL1,M1,1)
            DO K=1,2
              X0(K+1,NL1,NR1,M1)=F*X(K)
            END DO
            X0(4,NL1,NR1,M1)=X0(1,NL1,NR1,M1)
          END DO
C**DIFFERENTIATE TO THE LEFT
          F=H1R(NR1,M1,1)*MD
          DO NL1=1,JCI1
            DO K=1,2
              X0(K+4,NL1,NR1,M1)=F*H1L(NL1,M1,3-K)
            END DO
            X0(7,NL1,NR1,M1)=X0(1,NL1,NR1,M1)
          END DO
        END DO
      END DO

CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          F=H2R(NR2,M2,1)*IFACTC
          DO NL2=1,JCI2
            Y0(1,NL2,NR2,M2)=F*H2L(NL2,M2,1)
          END DO
C**DIFFERENTIATE TO THE RIGHT
          DO K=1,2
            Y(K)=H2R(NR2,M2,K)*IFACTC
          END DO
          DO NL2=1,JCI2
            F=H2L(NL2,M2,1)
            DO K=1,2
              Y0(K+1,NL2,NR2,M2)=F*Y(K)
            END DO
            Y0(4,NL2,NR2,M2)=Y0(1,NL2,NR2,M2)
          END DO
C**DIFFERENTIATE TO THE LEFT
          F=H2R(NR2,M2,1)*IFACTC
          DO NL2=1,JCI2
            DO K=1,2
              Y0(K+4,NL2,NR2,M2)=F*H2L(NL2,M2,K)
            END DO
            Y0(7,NL2,NR2,M2)=Y0(1,NL2,NR2,M2)
          END DO
        END DO
      END DO

CCCC  DO M3=1,MM3/MD3
      DO M3=1,MM3
        DO NR3=1,JCI3
          F=H3R(NR3,M3,1)
          DO NL3=1,JCI3
            Z0(1,NL3,NR3,M3)=F*H3L(NL3,M3,1)
          END DO
C**DIFFERENTIATE TO THE RIGHT
          DO K=1,2
            Z(K)=H3R(NR3,M3,K)
          END DO
          DO NL3=1,JCI3
            Z0(2,NL3,NR3,M3)=Z0(1,NL3,NR3,M3)
            F=H3L(NL3,M3,1)
            DO K=1,2
              Z0(K+2,NL3,NR3,M3)=F*Z(K)
            END DO
          END DO
C**DIFFERENTIATE TO THE LEFT
          F=H3R(NR3,M3,1)
          DO NL3=1,JCI3
            Z0(5,NL3,NR3,M3)=Z0(1,NL3,NR3,M3)
            DO K=1,2
              Z0(K+5,NL3,NR3,M3)=F*H3L(NL3,M3,K)
            END DO
          END DO
        END DO
      END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**START 3-MODE INTEGRATION
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO IRHS3=1,JCI3
          DO IRHS2=1,JCI2
            DO ILHS3=1,JCI3
              DO ILHS2=1,JCI2
                DO K=1,7
                  TEMP(K,ILHS2,ILHS3,IRHS2,IRHS3)=0
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(IABC.LT.7)THEN
C**START 2-MODE INTEGRATION
CCCC      DO M2=1,MM2/MD2
          DO M2=1,MM2
            DO IRHS3=1,JCI3
              DO ILHS3=1,JCI3
                TEMP1(1,ILHS3,IRHS3)=0
              END DO
            END DO
C**START 1-MODE INTEGRATION
            IF(JCOUPL.GT.0)THEN
CCCC          DO M3=1,MM3/MD3
              DO M3=1,MM3
                DO IRHS3=1,JCI3
                  DO ILHS3=1,JCI3
                    TEMP1(1,ILHS3,IRHS3)=TEMP1(1,ILHS3,IRHS3)+
     1              Z0(1,ILHS3,IRHS3,M3)*VM(M3,M2,M1,IABC)
                  END DO
                END DO
              END DO
            ELSE
CCCC          DO M3=1,MM3/MD3
              DO M3=1,MM3
                DO IRHS3=1,JCI3
                  DO ILHS3=1,JCI3
                    TEMP1(1,ILHS3,IRHS3)=TEMP1(1,ILHS3,IRHS3)+
     1              Z0(1,ILHS3,IRHS3,M3)*VMR(M3,M2,M1,IABC)
                  END DO
                END DO
              END DO
            END IF
C**END 1-MODE INTEGRATION
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    TEMP(1,ILHS2,ILHS3,IRHS2,IRHS3)=
     1              TEMP(1,ILHS2,ILHS3,IRHS2,IRHS3)+
     2              Y0(1,ILHS2,IRHS2,M2)*TEMP1(1,ILHS3,IRHS3)
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 2-MODE INTEGRATION
        ELSE
C**START 2-MODE INTEGRATION
CCCC      DO M2=1,MM2/MD2
          DO M2=1,MM2
            DO IRHS3=1,JCI3
              DO ILHS3=1,JCI3
                DO K=1,7
                  TEMP1(K,ILHS3,IRHS3)=0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
            IF(JCOUPL.GT.0)THEN
CCCC          DO M3=1,MM3/MD3
              DO M3=1,MM3
                ZZ=VM(M3,M2,M1,16)
                DO IRHS3=1,JCI3
                  DO ILHS3=1,JCI3
C**DIFFERENTIATE TO RIGHT
                    DO K=1,3
                      TEMP1(K+1,ILHS3,IRHS3)=TEMP1(K+1,ILHS3,IRHS3)+
     1                Z0(K+1,ILHS3,IRHS3,M3)*VM(M3,M2,M1,6+K+IOFF)
                    END DO
C**DIFFERENTIATE TO LEFT
                    DO K=1,3
                      TEMP1(K+4,ILHS3,IRHS3)=TEMP1(K+4,ILHS3,IRHS3)-
     1                Z0(K+4,ILHS3,IRHS3,M3)*VM(M3,M2,M1,6+K+IOFF)
                    END DO
C**EXTRA TERM FOR LINEAR
                    TEMP1(1,ILHS3,IRHS3)=TEMP1(1,ILHS3,IRHS3)+
     1              ZZ*Z0(1,ILHS3,IRHS3,M3)*KKAR*(KKAR-KKAL)*(IFACT1+
     2              IFACT2+IFACT3)
                  END DO
                END DO
              END DO
            ELSE
CCCC          DO M3=1,MM3/MD3
              DO M3=1,MM3
                ZZ=VMR(M3,M2,M1,16)
                DO IRHS3=1,JCI3
                  DO ILHS3=1,JCI3
C**DIFFERENTIATE TO RIGHT
                    DO K=1,3
                      TEMP1(K+1,ILHS3,IRHS3)=TEMP1(K+1,ILHS3,IRHS3)+
     1                Z0(K+1,ILHS3,IRHS3,M3)*VMR(M3,M2,M1,6+K+IOFF)
                    END DO
C**DIFFERENTIATE TO LEFT
                    DO K=1,3
                      TEMP1(K+4,ILHS3,IRHS3)=TEMP1(K+4,ILHS3,IRHS3)-
     1                Z0(K+4,ILHS3,IRHS3,M3)*VMR(M3,M2,M1,6+K+IOFF)
                    END DO
C**EXTRA TERM FOR LINEAR
                    TEMP1(1,ILHS3,IRHS3)=TEMP1(1,ILHS3,IRHS3)+
     1              ZZ*Z0(1,ILHS3,IRHS3,M3)*KKAR*(KKAR-KKAL)*(IFACT1+
     2              IFACT2+IFACT3)
                  END DO
                END DO
              END DO
            END IF
C**END 1-MODE INTEGRATION
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO K=1,7
                      TEMP(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1                TEMP(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2                Y0(K,ILHS2,IRHS2,M2)*TEMP1(K,ILHS3,IRHS3)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 2-MODE INTEGRATION
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
        DO IRHS=1,NSIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          NR3=IP(IRHS,MODE3)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            NL3=IP(ILHS,MODE3)
            DO K=1,7
              XA(ILHS+J0)=XA(ILHS+J0)+
     1        TEMP(K,NL2,NL3,NR2,NR3)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 3-MODE INTEGRATION
      IF(ITIM3A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM3A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MV3(NMODE,MODE1,MODE2,MODE3,MODE4,MOD1,MOD2,MOD3,
     1H1,XQ1,H2,XQ2,H3,XQ3,HTAU,XQTAU,NN1,MM1,NN2,MM2,NN3,MM3,NNTAU,
     2MMTAU,XA,NSIZE,IP,ISIZE,TEMP,TEMP1,TEMP2,JCI1,JCI2,JCI3,JCIM,X0,
     3Y0,Z0,T0,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      REAL*8 VM(MM3,MM2,MM1,21)
      REAL*4 VMR(MM3,MM2,MM1,21)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),IP(ISIZE,NMODE)
      DIMENSION TEMP1(9,JCI3,JCI3),TEMP2(9,JCI2,JCI3,JCI2,JCI3)
      DIMENSION TEMP(9,JCI1,JCI2,JCI3,JCI1,JCI2,JCI3)
      DIMENSION T0(9,JCIM,JCIM,MMTAU),X0(9,JCI1,JCI1,MM1)
      DIMENSION Y0(9,JCI2,JCI2,MM2),Z0(9,JCI3,JCI3,MM3)
      DIMENSION X(9),Y(9),C(9)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQTAU(MMTAU)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N2.EQ.N3)MD3=1
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      IF(N2.EQ.NT.AND.MDT.GT.1)MD2=1
      IF(N3.EQ.NT.AND.MDT.GT.1)MD3=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.GT.1)MD2=1
      IF(N1T.EQ.N3.AND.MDT.GT.1)MD3=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.GT.1)MD3=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.GT.1)MD3=1
      CALL VDMV3(NMODE,MODE1,MODE2,MODE3,MODE4,MOD1,MOD2,MOD3,H1,XQ1,
     1H2,XQ2,H3,XQ3,HTAU,XQTAU,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,
     2MM3/MD3,NNTAU,MMTAU,XA,NSIZE,IP,ISIZE,TEMP,TEMP1,TEMP2,JCI1,JCI2,
     3JCI3,JCIM,X0,Y0,Z0,T0,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDMV3(NMODE,MODE1,MODE2,MODE3,MODE4,MOD1,MOD2,MOD3,
     1H1,XQ1,H2,XQ2,H3,XQ3,HTAU,XQTAU,NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,
     2MM3,NNTAU,MMTAU,XA,NSIZE,IP,ISIZE,TEMP,TEMP1,TEMP2,JCI1,JCI2,
     3JCI3,JCIM,X0,Y0,Z0,T0,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      REAL*8 VM(MM3,MM2,MM1,21)
      REAL*4 VMR(MM3,MM2,MM1,21)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION IP(ISIZE,NMODE),HTAU(NNTAU,MMTAU,3,1)
      DIMENSION TEMP1(9,JCI3,JCI3),TEMP2(9,JCI2,JCI3,JCI2,JCI3)
      DIMENSION TEMP(9,JCI1,JCI2,JCI3,JCI1,JCI2,JCI3)
      DIMENSION T0(9,JCIM,JCIM,MMTAU),X0(9,JCI1,JCI1,MM1)
      DIMENSION Y0(9,JCI2,JCI2,MM2),Z0(9,JCI3,JCI3,MM3)
      DIMENSION X(9),Y(9),C(9)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQTAU(MMTAU)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP

      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MV3'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACTC=JNTFAC(NMODE,ICOUPC,3)

      KAL=KROTL/2
      INCTL=MOD(IFLAUD,2)*MOD(KAL,2)
      LMAXL=IFLAUD-(IFLAUD-1)*MOD(KAL+1,2)
      KAR=KROTR/2
      INCTR=MOD(IFLAUD,2)*MOD(KAR,2)
      LMAXR=IFLAUD-(IFLAUD-1)*MOD(KAR+1,2)
      FACTRC=0
      IF(J21.GT.1)FACTRC=IFACTC
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N2.EQ.N3)MD3=1
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      IF(N2.EQ.NT.AND.MDT.GT.1)MD2=1
      IF(N3.EQ.NT.AND.MDT.GT.1)MD3=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.GT.1)MD2=1
      IF(N1T.EQ.N3.AND.MDT.GT.1)MD3=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.GT.1)MD3=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.GT.1)MD3=1
      MD=MD1*MD2*MD3*MDT

C***********************************************************

      IOFF=(IABC-7)*5
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,1)*MD*IFACTC
          X(2)=X(1)
          X(3)=H1(NR1,M1,2)*MD*IFACTC
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(1)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,1)
            Y(2)=H1(NL1,M1,2)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            DO K=1,9
              X0(K,NL1,NR1,M1)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=H2(NR2,M2,2)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=H2(NL2,M2,2)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            DO K=1,9
              Y0(K,NL2,NR2,M2)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M3=1,MM3/MD3
      DO M3=1,MM3
        DO NR3=1,JCI3
          X(1)=H3(NR3,M3,1)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=H3(NR3,M3,2)
          X(8)=X(1)
          X(9)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=H3(NL3,M3,2)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            DO K=1,9
              Z0(K,NL3,NR3,M3)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
      DO MTAU=1,MMTAU/MDT
        DO NRR=1,JCIM
          NR=NRR+1-MOD(KAR,2)+INCTR
          X(1)=HTAU(NR,MTAU,1,LMAXR)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=HTAU(NR,MTAU,2,LMAXR)
          DO NLL=1,JCIM
            NL=NLL+1-MOD(KAL,2)+INCTL
            Y(1)=HTAU(NL,MTAU,1,LMAXL)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=HTAU(NL,MTAU,2,LMAXL)
            Y(9)=Y(1)
            DO K=1,9
              T0(K,NLL,NRR,MTAU)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**LOOP ROUND TAU (START 4-MODE INTEGRATION)
      ITAU=INIT-INCTAU
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          ITAU=ITAU+INCTAU
CCCC      IF(ITAU.GT.362)ITAU=ITAU-360
          IF(ITAU.GT.722)ITAU=ITAU-720
        ELSE
CCCC      QTAU=XQTAU(MTAU)
C**DELTA4 AND DELTA5 FROM TORSION (ARBITRARY SET EULER 'GAMMA' TO ZERO)
CCCC      DELTA(4)=+QTAU
CCCC      DELTA(5)=-QTAU
        END IF

C***********************************************************

        IF(JCOUPC.GE.0)THEN
          IF(ICOUPC.GE.3)READ(I93)VM
        ELSE
          IF(ICOUPC.GE.3)READ(I93)VMR
        END IF

C***********************************************************

        DO IRHS3=1,JCI3
          DO IRHS2=1,JCI2
            DO IRHS1=1,JCI1
              DO ILHS3=1,JCI3
                DO ILHS2=1,JCI2
                  DO ILHS1=1,JCI1
                    DO K=1,9
                      TEMP(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=0
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

        IF(IABC.LT.7)THEN
C**START 3-MODE INTEGRATION
CCCC      DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    TEMP2(1,ILHS2,ILHS3,IRHS2,IRHS3)=0
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
CCCC        DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS3=1,JCI3
                DO ILHS3=1,JCI3
                  TEMP1(1,ILHS3,IRHS3)=0
                END DO
              END DO
C**START 1-MODE INTEGRATION
              IF(JCOUPL.GT.0)THEN
CCCC            DO M3=1,MM3/MD3
                DO M3=1,MM3
                  DO IRHS3=1,JCI3
                    DO ILHS3=1,JCI3
                      TEMP1(1,ILHS3,IRHS3)=TEMP1(1,ILHS3,IRHS3)+
     1                Z0(1,ILHS3,IRHS3,M3)*VM(M3,M2,M1,IABC)
                    END DO
                  END DO
                END DO
              ELSE
CCCC            DO M3=1,MM3/MD3
                DO M3=1,MM3
                  DO IRHS3=1,JCI3
                    DO ILHS3=1,JCI3
                      TEMP1(1,ILHS3,IRHS3)=TEMP1(1,ILHS3,IRHS3)+
     1                Z0(1,ILHS3,IRHS3,M3)*VMR(M3,M2,M1,IABC)
                    END DO
                  END DO
                END DO
              END IF
C**END 1-MODE INTEGRATION
              DO IRHS3=1,JCI3
                DO IRHS2=1,JCI2
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      TEMP2(1,ILHS2,ILHS3,IRHS2,IRHS3)=
     1                TEMP2(1,ILHS2,ILHS3,IRHS2,IRHS3)+
     2                Y0(1,ILHS2,IRHS2,M2)*TEMP1(1,ILHS3,IRHS3)
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 2-MODE INTEGRATION
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO IRHS1=1,JCI1
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO ILHS1=1,JCI1
                        TEMP(1,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=
     1                  TEMP(1,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)+
     2                  X0(1,ILHS1,IRHS1,M1)*
     3                  TEMP2(1,ILHS2,ILHS3,IRHS2,IRHS3)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 3-MODE INTEGRATION
        ELSE
C**START 3-MODE INTEGRATION
CCCC      DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO K=1,9
                      TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=0
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
CCCC        DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS3=1,JCI3
                DO ILHS3=1,JCI3
                  DO K=1,9
                    TEMP1(K,ILHS3,IRHS3)=0
                  END DO
                END DO
              END DO
C**START 1-MODE INTEGRATION
CCCC          DO M3=1,MM3/MD3
              DO M3=1,MM3
                IF(JCOUPL.GT.0)THEN
                  C(1)=VM(M3,M2,M1,7+IOFF)
                  C(2)=VM(M3,M2,M1,8+IOFF)
                  C(3)=-VM(M3,M2,M1,8+IOFF)
                  C(4)=VM(M3,M2,M1,9+IOFF)
                  C(5)=-VM(M3,M2,M1,9+IOFF)
                  C(6)=VM(M3,M2,M1,10+IOFF)
                  C(7)=-VM(M3,M2,M1,10+IOFF)
                  C(8)=VM(M3,M2,M1,11+IOFF)
                  C(9)=-VM(M3,M2,M1,11+IOFF)
                ELSE
                  C(1)=VMR(M3,M2,M1,7+IOFF)
                  C(2)=VMR(M3,M2,M1,8+IOFF)
                  C(3)=-VMR(M3,M2,M1,8+IOFF)
                  C(4)=VMR(M3,M2,M1,9+IOFF)
                  C(5)=-VMR(M3,M2,M1,9+IOFF)
                  C(6)=VMR(M3,M2,M1,10+IOFF)
                  C(7)=-VMR(M3,M2,M1,10+IOFF)
                  C(8)=VMR(M3,M2,M1,11+IOFF)
                  C(9)=-VMR(M3,M2,M1,11+IOFF)
                END IF
                DO IRHS3=1,JCI3
                  DO ILHS3=1,JCI3
                    DO K=1,9
                      TEMP1(K,ILHS3,IRHS3)=TEMP1(K,ILHS3,IRHS3)+
     1                Z0(K,ILHS3,IRHS3,M3)*C(K)
                    END DO
                  END DO
                END DO
              END DO
C**END 1-MODE INTEGRATION
              DO IRHS3=1,JCI3
                DO IRHS2=1,JCI2
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO K=1,9
                        TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1                  TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2                  Y0(K,ILHS2,IRHS2,M2)*TEMP1(K,ILHS3,IRHS3)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 2-MODE INTEGRATION
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO IRHS1=1,JCI1
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO ILHS1=1,JCI1
                        DO K=1,9
                          TEMP(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=
     1                    TEMP(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)+
     2                    X0(K,ILHS1,IRHS1,M1)*
     3                    TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 3-MODE INTEGRATION
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (4-DIM)
        DO IRHS=1,NSIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          NR3=IP(IRHS,MODE3)
          IRTAU=IP(IRHS,MODE4)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            NL3=IP(ILHS,MODE3)
            ILTAU=IP(ILHS,MODE4)
            DO K=1,9
              XA(ILHS+J0)=XA(ILHS+J0)+DSTAU(ITAU)*
     1        TEMP(K,NL1,NL2,NL3,NR1,NR2,NR3)*T0(K,ILTAU,IRTAU,MTAU)
            END DO
          END DO
        END DO
      END DO
C**END 4-MODE INTEGRATION
      IF(ITIM4A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM4A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMI3(NMODE,MODE1,MODE2,MODE3,XA,XK,NSIZE,IPL,IPR,
     1ISIZMX,ISIZEL,ISIZER,IP3,ISIZE3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE),IP3(ISIZE3,3)
      DIMENSION XA(ISIZEL,ISIZER),XK(1)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM3A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMI3'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IPR(IRHS,MODE1)
        NR2=IPR(IRHS,MODE2)
        NR3=IPR(IRHS,MODE3)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP3(IR,1).AND.NR2.EQ.IP3(IR,2).AND.
     1       NR3.EQ.IP3(IR,3))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.(
     1      IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPL(ILHS,MODE1)
          NL2=IPL(ILHS,MODE2)
          NL3=IPL(ILHS,MODE3)
C**FIND LHS INDEX
          DO IL=1,NSIZE
            IF(NL1.EQ.IP3(IL,1).AND.NL2.EQ.IP3(IL,2).AND.
     1         NL3.EQ.IP3(IL,3))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM3A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM3A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMV3(NMODE,MODE1,MODE2,MODE3,XA,XK,NSIZE,IPL,IPR,
     1ISIZMX,ISIZEL,ISIZER,IP4,ISIZE4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE)
      DIMENSION XA(ISIZEL,ISIZER),XK(1),IP4(ISIZE4,4)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM3A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMV3'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IPR(IRHS,MODE1)
        NR2=IPR(IRHS,MODE2)
        NR3=IPR(IRHS,MODE3)
        NRTAU=IPR(IRHS,NMODE)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP4(IR,1).AND.NR2.EQ.IP4(IR,2).AND.
     1       NR3.EQ.IP4(IR,3).AND.NRTAU.EQ.IP4(IR,4))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1      NMODE.AND.(IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPL(ILHS,MODE1)
          NL2=IPL(ILHS,MODE2)
          NL3=IPL(ILHS,MODE3)
          NLTAU=IPL(ILHS,NMODE)
C**FIND LHS INDEX
          DO IL=1,NSIZE
            IF(NL1.EQ.IP4(IL,1).AND.NL2.EQ.IP4(IL,2).AND.
     1         NL3.EQ.IP4(IL,3).AND.NLTAU.EQ.IP4(IL,4))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM3A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM3A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MI4(NMODE,MODE1,MODE2,MODE3,MODE4,MOD1,MOD2,MOD3,
     1MOD4,H1L,H1R,XQ1,H2L,H2R,XQ2,H3L,H3R,XQ3,H4L,H4R,XQ4,NN1,MM1,NN2,
     2MM2,NN3,MM3,NN4,MM4,XA,NSIZE,IP,ISIZE,TEMP,TEMP1,TEMP2,TEMP3,
     3JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VM,VMR,J21,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM4,MM3,MM2,MM1,19)
      REAL*4 VMR(MM4,MM3,MM2,MM1,19)
      DIMENSION MODINT(NMODE)
      DIMENSION H1L(NN1,MM1,3),H2L(NN2,MM2,3),H3L(NN3,MM3,3),
     1H4L(NN4,MM4,3)
      DIMENSION H1R(NN1,MM1,3),H2R(NN2,MM2,3),H3R(NN3,MM3,3),
     1H4R(NN4,MM4,3)
      DIMENSION IP(ISIZE,NMODE)
      DIMENSION TEMP1(9,JCI4,JCI4),TEMP2(9,JCI3,JCI4,JCI3,JCI4)
C     DIMENSION TEMP3(4,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION TEMP(9,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(9,JCI1,JCI1,MM1),Y0(9,JCI2,JCI2,MM2)
      DIMENSION Z0(9,JCI3,JCI3,MM3),W0(9,JCI4,JCI4,MM4)
      DIMENSION X(2),Y(2),Z(2),W(2)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
        END DO
      END DO
      IF(N3.EQ.N4)MD4=1
      IF(N2.EQ.N4)MD4=1
      IF(N1.EQ.N4)MD4=1
      IF(N2.EQ.N3)MD3=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N2)MD2=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      CALL VDMI4(NMODE,MODE1,MODE2,MODE3,MODE4,MOD1,MOD2,MOD3,MOD4,H1L,
     1H1R,XQ1,H2L,H2R,XQ2,H3L,H3R,XQ3,H4L,H4R,XQ4,NN1,MM1,MM1/MD1,NN2,
     2MM2,MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,XA,NSIZE,IP,ISIZE,
     3TEMP,TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VM,VMR,
     4J21,IABC,MODINT)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDMI4(NMODE,MODE1,MODE2,MODE3,MODE4,MOD1,MOD2,MOD3,
     1MOD4,H1L,H1R,XQ1,H2L,H2R,XQ2,H3L,H3R,XQ3,H4L,H4R,XQ4,NN1,MH1,MM1,
     2NN2,MH2,MM2,NN3,MH3,MM3,NN4,MH4,MM4,XA,NSIZE,IP,ISIZE,TEMP,TEMP1,
     3TEMP2,TEMP3,JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VM,VMR,J21,IABC,
     4MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM4,MM3,MM2,MM1,19)
      REAL*4 VMR(MM4,MM3,MM2,MM1,19)
      DIMENSION MODINT(NMODE)
      DIMENSION H1L(NN1,MH1,3),H2L(NN2,MH2,3),H3L(NN3,MH3,3),
     1H4L(NN4,MH4,3)
      DIMENSION H1R(NN1,MH1,3),H2R(NN2,MH2,3),H3R(NN3,MH3,3),
     1H4R(NN4,MH4,3)
      DIMENSION IP(ISIZE,NMODE)
      DIMENSION TEMP1(9,JCI4,JCI4),TEMP2(9,JCI3,JCI4,JCI3,JCI4)
C     DIMENSION TEMP3(4,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION TEMP(9,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(9,JCI1,JCI1,MM1),Y0(9,JCI2,JCI2,MM2)
      DIMENSION Z0(9,JCI3,JCI3,MM3),W0(9,JCI4,JCI4,MM4)
      DIMENSION X(2),Y(2),Z(2),W(2)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP

      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MI4'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACTC=INTFAC(NMODE,ICOUPC,4)

C***********************************************************

      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.4)READ(I94)VM
      ELSE
        IF(ICOUPC.GE.4)READ(I94)VMR
      END IF

C***********************************************************

      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
        END DO
      END DO
      IF(N3.EQ.N4)MD4=1
      IF(N2.EQ.N4)MD4=1
      IF(N1.EQ.N4)MD4=1
      IF(N2.EQ.N3)MD3=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N2)MD2=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      MD=MD1*MD2*MD3*MD4

C***********************************************************

      IOFF=(IABC-7)*4
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          F=H1R(NR1,M1,1)*MD
          DO NL1=1,JCI1
            X0(1,NL1,NR1,M1)=F*H1L(NL1,M1,1)
          END DO
C**DIFFERENTIATE TO THE RIGHT
          DO K=1,2
            X(K)=H1R(NR1,M1,3-K)*MD
          END DO
          DO NL1=1,JCI1
            F=H1L(NL1,M1,1)
            DO K=1,2
              X0(K+1,NL1,NR1,M1)=F*X(K)
            END DO
            X0(4,NL1,NR1,M1)=X0(1,NL1,NR1,M1)
            X0(5,NL1,NR1,M1)=X0(1,NL1,NR1,M1)
          END DO
C**DIFFERENTIATE TO THE LEFT
          F=H1R(NR1,M1,1)*MD
          DO NL1=1,JCI1
            DO K=1,2
              X0(K+5,NL1,NR1,M1)=F*H1L(NL1,M1,3-K)
            END DO
            X0(8,NL1,NR1,M1)=X0(1,NL1,NR1,M1)
            X0(9,NL1,NR1,M1)=X0(1,NL1,NR1,M1)
          END DO
        END DO
      END DO

CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          F=H2R(NR2,M2,1)*IFACTC
          DO NL2=1,JCI2
            Y0(1,NL2,NR2,M2)=F*H2L(NL2,M2,1)
          END DO
C**DIFFERENTIATE TO THE RIGHT
          DO K=1,2
            Y(K)=H2R(NR2,M2,K)*IFACTC
          END DO
          DO NL2=1,JCI2
            F=H2L(NL2,M2,1)
            DO K=1,2
              Y0(K+1,NL2,NR2,M2)=F*Y(K)
            END DO
            Y0(4,NL2,NR2,M2)=Y0(1,NL2,NR2,M2)
            Y0(5,NL2,NR2,M2)=Y0(1,NL2,NR2,M2)
          END DO
C**DIFFERENTIATE TO THE LEFT
          F=H2R(NR2,M2,1)*IFACTC
          DO NL2=1,JCI2
            DO K=1,2
              Y0(K+5,NL2,NR2,M2)=F*H2L(NL2,M2,K)
            END DO
            Y0(8,NL2,NR2,M2)=Y0(1,NL2,NR2,M2)
            Y0(9,NL2,NR2,M2)=Y0(1,NL2,NR2,M2)
          END DO
        END DO
      END DO

CCCC  DO M3=1,MM3/MD3
      DO M3=1,MM3
        DO NR3=1,JCI3
          F=H3R(NR3,M3,1)
          DO NL3=1,JCI3
            Z0(1,NL3,NR3,M3)=F*H3L(NL3,M3,1)
          END DO
C**DIFFERENTIATE TO THE RIGHT
          DO K=1,2
            Z(K)=H3R(NR3,M3,K)
          END DO
          DO NL3=1,JCI3
            Z0(2,NL3,NR3,M3)=Z0(1,NL3,NR3,M3)
            F=H3L(NL3,M3,1)
            DO K=1,2
              Z0(K+2,NL3,NR3,M3)=F*Z(K)
            END DO
            Z0(5,NL3,NR3,M3)=Z0(1,NL3,NR3,M3)
          END DO
C**DIFFERENTIATE TO THE LEFT
          F=H3R(NR3,M3,1)
          DO NL3=1,JCI3
            Z0(6,NL3,NR3,M3)=Z0(1,NL3,NR3,M3)
            DO K=1,2
              Z0(K+6,NL3,NR3,M3)=F*H3L(NL3,M3,K)
            END DO
            Z0(9,NL3,NR3,M3)=Z0(1,NL3,NR3,M3)
          END DO
        END DO
      END DO

CCCC  DO M4=1,MM4/MD4
      DO M4=1,MM4
        DO NR4=1,JCI4
          F=H4R(NR4,M4,1)
          DO NL4=1,JCI4
            W0(1,NL4,NR4,M4)=F*H4L(NL4,M4,1)
          END DO
C**DIFFERENTIATE TO THE RIGHT
          DO K=1,2
            W(K)=H4R(NR4,M4,K)
          END DO
          DO NL4=1,JCI4
            W0(2,NL4,NR4,M4)=W0(1,NL4,NR4,M4)
            W0(3,NL4,NR4,M4)=W0(1,NL4,NR4,M4)
            F=H4L(NL4,M4,1)
            DO K=1,2
              W0(K+3,NL4,NR4,M4)=F*W(K)
            END DO
          END DO
C**DIFFERENTIATE TO THE LEFT
          F=H4R(NR4,M4,1)
          DO NL4=1,JCI4
            W0(6,NL4,NR4,M4)=W0(1,NL4,NR4,M4)
            W0(7,NL4,NR4,M4)=W0(1,NL4,NR4,M4)
            DO K=1,2
              W0(K+7,NL4,NR4,M4)=F*H4L(NL4,M4,K)
            END DO
          END DO

        END DO
      END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**START 4-MODE INTEGRATION
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO IRHS4=1,JCI4
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO ILHS4=1,JCI4
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO K=1,9
                      TEMP(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=0
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(IABC.LT.7)THEN
C**START 3-MODE INTEGRATION
CCCC      DO M2=1,MM2/MD2
          DO M2=1,MM2
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    TEMP2(1,ILHS3,ILHS4,IRHS3,IRHS4)=0
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
CCCC        DO M3=1,MM3/MD3
            DO M3=1,MM3
              DO IRHS4=1,JCI4
                DO ILHS4=1,JCI4
                  TEMP1(1,ILHS4,IRHS4)=0
                END DO
              END DO
C**START 1-MODE INTEGRATION
              IF(JCOUPL.GT.0)THEN
CCCC            DO M4=1,MM4/MD4
                DO M4=1,MM4
                  DO IRHS4=1,JCI4
                    DO ILHS4=1,JCI4
                      TEMP1(1,ILHS4,IRHS4)=TEMP1(1,ILHS4,IRHS4)+
     1                W0(1,ILHS4,IRHS4,M4)*VM(M4,M3,M2,M1,IABC)
                    END DO
                  END DO
                END DO
              ELSE
CCCC            DO M4=1,MM4/MD4
                DO M4=1,MM4
                  DO IRHS4=1,JCI4
                    DO ILHS4=1,JCI4
                      TEMP1(1,ILHS4,IRHS4)=TEMP1(1,ILHS4,IRHS4)+
     1                W0(1,ILHS4,IRHS4,M4)*VMR(M4,M3,M2,M1,IABC)
                    END DO
                  END DO
                END DO
              END IF
C**END 1-MODE INTEGRATION
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS4=1,JCI4
                    DO ILHS3=1,JCI3
                      TEMP2(1,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                TEMP2(1,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                Z0(1,ILHS3,IRHS3,M3)*TEMP1(1,ILHS4,IRHS4)
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 2-MODE INTEGRATION
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO IRHS2=1,JCI2
                  DO ILHS4=1,JCI4
                    DO ILHS3=1,JCI3
                      DO ILHS2=1,JCI2
                        TEMP(1,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=
     1                  TEMP(1,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)+
     2                  Y0(1,ILHS2,IRHS2,M2)*
     3                  TEMP2(1,ILHS3,ILHS4,IRHS3,IRHS4)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 3-MODE INTEGRATION
        ELSE
C**START 3-MODE INTEGRATION
CCCC      DO M2=1,MM2/MD2
          DO M2=1,MM2
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO K=1,9
                      TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=0
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
CCCC        DO M3=1,MM3/MD3
            DO M3=1,MM3
              DO IRHS4=1,JCI4
                DO ILHS4=1,JCI4
                  DO K=1,9
                    TEMP1(K,ILHS4,IRHS4)=0
                  END DO
                END DO
              END DO
C**START 1-MODE INTEGRATION
              IF(JCOUPL.GT.0)THEN
CCCC            DO M4=1,MM4/MD4
                DO M4=1,MM4
                  DO IRHS4=1,JCI4
                    DO ILHS4=1,JCI4
C**DIFFERENTIATE TO RIGHT
                      DO K=1,4
                        TEMP1(K+1,ILHS4,IRHS4)=TEMP1(K+1,ILHS4,IRHS4)+
     1                  W0(K+1,ILHS4,IRHS4,M4)*VM(M4,M3,M2,M1,6+K+IOFF)
                      END DO
C**DIFFERENTIATE TO LEFT
                      DO K=1,4
                        TEMP1(K+5,ILHS4,IRHS4)=TEMP1(K+5,ILHS4,IRHS4)-
     1                  W0(K+5,ILHS4,IRHS4,M4)*VM(M4,M3,M2,M1,6+K+IOFF)
                      END DO
                    END DO
                  END DO
                END DO
              ELSE
CCCC            DO M4=1,MM4/MD4
                DO M4=1,MM4
                  DO IRHS4=1,JCI4
                    DO ILHS4=1,JCI4
C**DIFFERENTIATE TO RIGHT
                      DO K=1,4
                        TEMP1(K+1,ILHS4,IRHS4)=TEMP1(K+1,ILHS4,IRHS4)+
     1                  W0(K+1,ILHS4,IRHS4,M4)*VMR(M4,M3,M2,M1,6+K+IOFF)
                      END DO
C**DIFFERENTIATE TO LEFT
                      DO K=1,4
                        TEMP1(K+5,ILHS4,IRHS4)=TEMP1(K+5,ILHS4,IRHS4)-
     1                  W0(K+5,ILHS4,IRHS4,M4)*VMR(M4,M3,M2,M1,6+K+IOFF)
                      END DO
                    END DO
                  END DO
                END DO
              END IF
C**END 1-MODE INTEGRATION
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS4=1,JCI4
                    DO ILHS3=1,JCI3
                      DO K=1,9
                        TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                  TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                  Z0(K,ILHS3,IRHS3,M3)*TEMP1(K,ILHS4,IRHS4)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 2-MODE INTEGRATION
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO IRHS2=1,JCI2
                  DO ILHS4=1,JCI4
                    DO ILHS3=1,JCI3
                      DO ILHS2=1,JCI2
                        DO K=1,9
                          TEMP(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=
     1                    TEMP(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)+
     2                    Y0(K,ILHS2,IRHS2,M2)*
     3                    TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 3-MODE INTEGRATION
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (4-DIM)
        DO IRHS=1,NSIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          NR3=IP(IRHS,MODE3)
          NR4=IP(IRHS,MODE4)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            NL3=IP(ILHS,MODE3)
            NL4=IP(ILHS,MODE4)
            DO K=1,9
              XA(ILHS+J0)=XA(ILHS+J0)+
     1        TEMP(K,NL2,NL3,NL4,NR2,NR3,NR4)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 4-MODE INTEGRATION
      IF(ITIM4A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM4A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MV4(NMODE,MODE1,MODE2,MODE3,MODE4,MODE5,MOD1,MOD2,
     1MOD3,MOD4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,HTAU,XQTAU,NN1,MM1,NN2,MM2,
     2NN3,MM3,NN4,MM4,NNTAU,MMTAU,XA,NSIZE,IP,ISIZE,TEMP1,TEMP2,TEMP3,
     3JCI1,JCI2,JCI3,JCI4,JCIM,X0,Y0,Z0,W0,T0,VM,VMR,J21,KROTL,KROTR,
     4IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      REAL*8 VM(MM4,MM3,MM2,MM1,24)
      REAL*4 VMR(MM4,MM3,MM2,MM1,24)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),IP(ISIZE,NMODE)
      DIMENSION TEMP1(11,JCI4,JCI4),TEMP2(11,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(11,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION T0(11,JCIM,JCIM,MMTAU)
      DIMENSION X0(11,JCI1,JCI1,MM1),Y0(11,JCI2,JCI2,MM2)
      DIMENSION Z0(11,JCI3,JCI3,MM3),W0(11,JCI4,JCI4,MM4)
      DIMENSION X(11),Y(11),C(11)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQTAU(MMTAU)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N3.EQ.N4)MD4=1
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      IF(N2.EQ.NT.AND.MDT.GT.1)MD2=1
      IF(N3.EQ.NT.AND.MDT.GT.1)MD3=1
      IF(N4.EQ.NT.AND.MDT.GT.1)MD4=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.GT.1)MD2=1
      IF(N1T.EQ.N3.AND.MDT.GT.1)MD3=1
      IF(N1T.EQ.N4.AND.MDT.GT.1)MD4=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.GT.1)MD3=1
      IF(N2T.EQ.N4.AND.MDT.GT.1)MD4=1
      N3T=ISYMP(N3,NT)
      IF(N3T.EQ.N4.AND.MDT.GT.1)MD4=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.GT.1)MD3=1
      IF(N12T.EQ.N4.AND.MDT.GT.1)MD4=1
      N13T=ISYMP(N13,NT)
      IF(N13T.EQ.N4.AND.MDT.GT.1)MD4=1
      N23T=ISYMP(N23,NT)
      IF(N23T.EQ.N4.AND.MDT.GT.1)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      N123T=ISYMP(N123,NT)
      IF(N123T.EQ.N4.AND.MDT.GT.1)MD4=1
      CALL VDMV4(NMODE,MODE1,MODE2,MODE3,MODE4,MODE5,MOD1,MOD2,MOD3,
     1MOD4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,HTAU,XQTAU,NN1,MM1,MM1/MD1,NN2,
     2MM2,MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,NNTAU,MMTAU,XA,NSIZE,
     3IP,ISIZE,TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,JCI4,JCIM,X0,Y0,Z0,W0,
     4T0,VM,VMR,J21,KROTL,KROTR,IABC,MODINT)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDMV4(NMODE,MODE1,MODE2,MODE3,MODE4,MODE5,MOD1,MOD2,
     1MOD3,MOD4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,HTAU,XQTAU,NN1,MH1,MM1,NN2,
     2MH2,MM2,NN3,MH3,MM3,NN4,MH4,MM4,NNTAU,MMTAU,XA,NSIZE,IP,ISIZE,
     3TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,JCI4,JCIM,X0,Y0,Z0,W0,T0,VM,VMR,
     4J21,KROTL,KROTR,IABC,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      REAL*8 VM(MM4,MM3,MM2,MM1,24)
      REAL*4 VMR(MM4,MM3,MM2,MM1,24)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),IP(ISIZE,NMODE)
      DIMENSION TEMP1(11,JCI4,JCI4),TEMP2(11,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(11,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION T0(11,JCIM,JCIM,MMTAU)
      DIMENSION X0(11,JCI1,JCI1,MM1),Y0(11,JCI2,JCI2,MM2)
      DIMENSION Z0(11,JCI3,JCI3,MM3),W0(11,JCI4,JCI4,MM4)
      DIMENSION X(11),Y(11),C(11)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQTAU(MMTAU)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/FILASS/IOUT,INP

      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MV4'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACTC=JNTFAC(NMODE,ICOUPC,4)

      KAL=KROTL/2
      INCTL=MOD(IFLAUD,2)*MOD(KAL,2)
      LMAXL=IFLAUD-(IFLAUD-1)*MOD(KAL+1,2)
      KAR=KROTR/2
      INCTR=MOD(IFLAUD,2)*MOD(KAR,2)
      LMAXR=IFLAUD-(IFLAUD-1)*MOD(KAR+1,2)
      FACTRC=0
      IF(J21.GT.1)FACTRC=IFACTC
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N3.EQ.N4)MD4=1
      IF(N1.EQ.NT.AND.MDT.GT.1)MD1=1
      IF(N2.EQ.NT.AND.MDT.GT.1)MD2=1
      IF(N3.EQ.NT.AND.MDT.GT.1)MD3=1
      IF(N4.EQ.NT.AND.MDT.GT.1)MD4=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.GT.1)MD2=1
      IF(N1T.EQ.N3.AND.MDT.GT.1)MD3=1
      IF(N1T.EQ.N4.AND.MDT.GT.1)MD4=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.GT.1)MD3=1
      IF(N2T.EQ.N4.AND.MDT.GT.1)MD4=1
      N3T=ISYMP(N3,NT)
      IF(N3T.EQ.N4.AND.MDT.GT.1)MD4=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.GT.1)MD3=1
      IF(N12T.EQ.N4.AND.MDT.GT.1)MD4=1
      N13T=ISYMP(N13,NT)
      IF(N13T.EQ.N4.AND.MDT.GT.1)MD4=1
      N23T=ISYMP(N23,NT)
      IF(N23T.EQ.N4.AND.MDT.GT.1)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      N123T=ISYMP(N123,NT)
      IF(N123T.EQ.N4.AND.MDT.GT.1)MD4=1
      MD=MD1*MD2*MD3*MD4*MDT

C***********************************************************

      IOFF=(IABC-7)*6
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,1)*MD*IFACTC
          X(2)=X(1)
          X(3)=H1(NR1,M1,2)*MD*IFACTC
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(1)
          X(10)=X(1)
          X(11)=X(1)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,1)
            Y(2)=H1(NL1,M1,2)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            Y(11)=Y(1)
            DO K=1,11
              X0(K,NL1,NR1,M1)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=H2(NR2,M2,2)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(1)
          X(10)=X(1)
          X(11)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=H2(NL2,M2,2)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            Y(11)=Y(1)
            DO K=1,11
              Y0(K,NL2,NR2,M2)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M3=1,MM3/MD3
      DO M3=1,MM3
        DO NR3=1,JCI3
          X(1)=H3(NR3,M3,1)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=H3(NR3,M3,2)
          X(8)=X(1)
          X(9)=X(1)
          X(10)=X(1)
          X(11)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=H3(NL3,M3,2)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            Y(11)=Y(1)
            DO K=1,11
              Z0(K,NL3,NR3,M3)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M4=1,MM4/MD4
      DO M4=1,MM4
        DO NR4=1,JCI4
          X(1)=H4(NR4,M4,1)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=H4(NR4,M4,2)
          X(10)=X(1)
          X(11)=X(1)
          DO NL4=1,JCI4
            Y(1)=H4(NL4,M4,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=H4(NL4,M4,2)
            Y(9)=Y(1)
            Y(10)=Y(1)
            Y(11)=Y(1)
            DO K=1,11
              W0(K,NL4,NR4,M4)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
      DO MTAU=1,MMTAU/MDT
        DO NRR=1,JCIM
          NR=NRR+1-MOD(KAR,2)+INCTR
          X(1)=HTAU(NR,MTAU,1,LMAXR)
          X(2)=X(1)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(1)
          X(10)=X(1)
          X(11)=HTAU(NR,MTAU,2,LMAXR)
          DO NLL=1,JCIM
            NL=NLL+1-MOD(KAL,2)+INCTL
            Y(1)=HTAU(NL,MTAU,1,LMAXL)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=HTAU(NL,MTAU,2,LMAXL)
            Y(11)=Y(1)
            DO K=1,11
              T0(K,NLL,NRR,MTAU)=Y(K)*X(K)
            END DO
          END DO
        END DO
      END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**LOOP ROUND TAU (START 5-MODE INTEGRATION)
      ITAU=INIT-INCTAU
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          ITAU=ITAU+INCTAU
CCCC      IF(ITAU.GT.362)ITAU=ITAU-360
          IF(ITAU.GT.722)ITAU=ITAU-720
        ELSE
CCCC      QTAU=XQTAU(MTAU)
C**DELTA4 AND DELTA5 FROM TORSION (ARBITRARY SET EULER 'GAMMA' TO ZERO)
CCCC      DELTA(4)=+QTAU
CCCC      DELTA(5)=-QTAU
        END IF

C***********************************************************

        IF(JCOUPC.GE.0)THEN
          IF(ICOUPC.GE.4)READ(I94)VM
        ELSE
          IF(ICOUPC.GE.4)READ(I94)VMR
        END IF

C***********************************************************

C**START 4-MODE INTEGRATION
CCCC    DO M1=1,MM1/MD1
        DO M1=1,MM1
          DO IRHS4=1,JCI4
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO K=1,11
                        TEMP3(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
          IF(IABC.LT.7)THEN
C**START 3-MODE INTEGRATION
CCCC        DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS4=1,JCI4
                    DO ILHS3=1,JCI3
                      TEMP2(1,ILHS3,ILHS4,IRHS3,IRHS4)=0
                    END DO
                  END DO
                END DO
              END DO
C**START 2-MODE INTEGRATION
CCCC          DO M3=1,MM3/MD3
              DO M3=1,MM3
                DO IRHS4=1,JCI4
                  DO ILHS4=1,JCI4
                    TEMP1(1,ILHS4,IRHS4)=0
                  END DO
                END DO
C**START 1-MODE INTEGRATION
                IF(JCOUPL.GT.0)THEN
CCCC              DO M4=1,MM4/MD4
                  DO M4=1,MM4
                    DO IRHS4=1,JCI4
                      DO ILHS4=1,JCI4
                        TEMP1(1,ILHS4,IRHS4)=TEMP1(1,ILHS4,IRHS4)+
     1                  W0(1,ILHS4,IRHS4,M4)*VM(M4,M3,M2,M1,IABC)
                      END DO
                    END DO
                  END DO
                ELSE
CCCC              DO M4=1,MM4/MD4
                  DO M4=1,MM4
                    DO IRHS4=1,JCI4
                      DO ILHS4=1,JCI4
                        TEMP1(1,ILHS4,IRHS4)=TEMP1(1,ILHS4,IRHS4)+
     1                  W0(1,ILHS4,IRHS4,M4)*VMR(M4,M3,M2,M1,IABC)
                      END DO
                    END DO
                  END DO
                END IF
C**END 1-MODE INTEGRATION
                DO IRHS4=1,JCI4
                  DO IRHS3=1,JCI3
                    DO ILHS4=1,JCI4
                      DO ILHS3=1,JCI3
                        TEMP2(1,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                  TEMP2(1,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                  Z0(1,ILHS3,IRHS3,M3)*TEMP1(1,ILHS4,IRHS4)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
C**END 2-MODE INTEGRATION
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO IRHS2=1,JCI2
                    DO ILHS4=1,JCI4
                      DO ILHS3=1,JCI3
                        DO ILHS2=1,JCI2
                          TEMP3(1,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=
     1                    TEMP3(1,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)+
     2                    Y0(1,ILHS2,IRHS2,M2)*
     3                    TEMP2(1,ILHS3,ILHS4,IRHS3,IRHS4)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 3-MODE INTEGRATION
          ELSE
C**START 3-MODE INTEGRATION
CCCC        DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS4=1,JCI4
                    DO ILHS3=1,JCI3
                      DO K=1,11
                        TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
C**START 2-MODE INTEGRATION
CCCC          DO M3=1,MM3/MD3
              DO M3=1,MM3
                DO IRHS4=1,JCI4
                  DO ILHS4=1,JCI4
                    DO K=1,11
                      TEMP1(K,ILHS4,IRHS4)=0
                    END DO
                  END DO
                END DO
C**START 1-MODE INTEGRATION
CCCC            DO M4=1,MM4/MD4
                DO M4=1,MM4
                  IF(JCOUPL.GT.0)THEN
                    C(1)=VM(M4,M3,M2,M1,7+IOFF)
                    C(2)=VM(M4,M3,M2,M1,8+IOFF)
                    C(3)=-VM(M4,M3,M2,M1,8+IOFF)
                    C(4)=VM(M4,M3,M2,M1,9+IOFF)
                    C(5)=-VM(M4,M3,M2,M1,9+IOFF)
                    C(6)=VM(M4,M3,M2,M1,10+IOFF)
                    C(7)=-VM(M4,M3,M2,M1,10+IOFF)
                    C(8)=VM(M4,M3,M2,M1,11+IOFF)
                    C(9)=-VM(M4,M3,M2,M1,11+IOFF)
                    C(10)=VM(M4,M3,M2,M1,12+IOFF)
                    C(11)=-VM(M4,M3,M2,M1,12+IOFF)
                  ELSE
                    C(1)=VMR(M4,M3,M2,M1,7+IOFF)
                    C(2)=VMR(M4,M3,M2,M1,8+IOFF)
                    C(3)=-VMR(M4,M3,M2,M1,8+IOFF)
                    C(4)=VMR(M4,M3,M2,M1,9+IOFF)
                    C(5)=-VMR(M4,M3,M2,M1,9+IOFF)
                    C(6)=VMR(M4,M3,M2,M1,10+IOFF)
                    C(7)=-VMR(M4,M3,M2,M1,10+IOFF)
                    C(8)=VMR(M4,M3,M2,M1,11+IOFF)
                    C(9)=-VMR(M4,M3,M2,M1,11+IOFF)
                    C(10)=VMR(M4,M3,M2,M1,12+IOFF)
                    C(11)=-VMR(M4,M3,M2,M1,12+IOFF)
                  END IF
                  DO IRHS4=1,JCI4
                    DO ILHS4=1,JCI4
                      DO K=1,11
                        TEMP1(K,ILHS4,IRHS4)=TEMP1(K,ILHS4,IRHS4)+
     1                  W0(K,ILHS4,IRHS4,M4)*C(K)
                      END DO
                    END DO
                  END DO
                END DO
C**END 1-MODE INTEGRATION
                DO IRHS4=1,JCI4
                  DO IRHS3=1,JCI3
                    DO ILHS4=1,JCI4
                      DO ILHS3=1,JCI3
                        DO K=1,11
                          TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                    TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                    Z0(K,ILHS3,IRHS3,M3)*TEMP1(K,ILHS4,IRHS4)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
C**END 2-MODE INTEGRATION
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO IRHS2=1,JCI2
                    DO ILHS4=1,JCI4
                      DO ILHS3=1,JCI3
                        DO ILHS2=1,JCI2
                          DO K=1,11
                        TEMP3(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=
     1                  TEMP3(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)+
     2                  Y0(K,ILHS2,IRHS2,M2)*
     3                  TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)
                          END DO
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 3-MODE INTEGRATION
          END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (5-DIM)
          DO IRHS=1,NSIZE
            NR1=IP(IRHS,MODE1)
            NR2=IP(IRHS,MODE2)
            NR3=IP(IRHS,MODE3)
            NR4=IP(IRHS,MODE4)
            IRTAU=IP(IRHS,MODE5)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP(ILHS,MODE1)
              NL2=IP(ILHS,MODE2)
              NL3=IP(ILHS,MODE3)
              NL4=IP(ILHS,MODE4)
              ILTAU=IP(ILHS,MODE5)
              DO K=1,11
                XA(ILHS+J0)=XA(ILHS+J0)+
     1          TEMP3(K,NL2,NL3,NL4,NR2,NR3,NR4)*X0(K,NL1,NR1,M1)*
     2          T0(K,ILTAU,IRTAU,MTAU)
              END DO
            END DO
          END DO
        END DO
C**END 4-MODE INTEGRATION
      END DO
C**END 5-MODE INTEGRATION
      IF(ITIM4A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM4A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMI4(NMODE,MODE1,MODE2,MODE3,MODE4,XA,XK,NSIZE,IPL,
     1IPR,ISIZMX,ISIZEL,ISIZER,IP4,ISIZE4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE),IP4(ISIZE4,4)
      DIMENSION XA(ISIZEL,ISIZER),XK(1)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMI4'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IPR(IRHS,MODE1)
        NR2=IPR(IRHS,MODE2)
        NR3=IPR(IRHS,MODE3)
        NR4=IPR(IRHS,MODE4)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP4(IR,1).AND.NR2.EQ.IP4(IR,2).AND.
     1       NR3.EQ.IP4(IR,3).AND.NR4.EQ.IP4(IR,4))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1      MODE4.AND.(IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPL(ILHS,MODE1)
          NL2=IPL(ILHS,MODE2)
          NL3=IPL(ILHS,MODE3)
          NL4=IPL(ILHS,MODE4)
C**FIND LHS INDEX
          DO IL=1,NSIZE
            IF(NL1.EQ.IP4(IL,1).AND.NL2.EQ.IP4(IL,2).AND.
     1         NL3.EQ.IP4(IL,3).AND.NL4.EQ.IP4(IL,4))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM4A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM4A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMV4(NMODE,MODE1,MODE2,MODE3,MODE4,XA,XK,NSIZE,IPL,
     1IPR,ISIZMX,ISIZEL,ISIZER,IP5,ISIZE5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE)
      DIMENSION XA(ISIZEL,ISIZER),XK(1),IP5(ISIZE5,5)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/HERM/IHERM
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMV4'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IPR(IRHS,MODE1)
        NR2=IPR(IRHS,MODE2)
        NR3=IPR(IRHS,MODE3)
        NR4=IPR(IRHS,MODE4)
        NRTAU=IPR(IRHS,NMODE)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP5(IR,1).AND.NR2.EQ.IP5(IR,2).AND.
     1       NR3.EQ.IP5(IR,3).AND.NR4.EQ.IP5(IR,4).AND.
     2       NRTAU.EQ.IP5(IR,5))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1      MODE4.AND.K.NE.NMODE.AND.(IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPL(ILHS,MODE1)
          NL2=IPL(ILHS,MODE2)
          NL3=IPL(ILHS,MODE3)
          NL4=IPL(ILHS,MODE4)
          NLTAU=IPL(ILHS,NMODE)
C**FIND LHS INDEX
          DO IL=1,NSIZE
            IF(NL1.EQ.IP5(IL,1).AND.NL2.EQ.IP5(IL,2).AND.
     1         NL3.EQ.IP5(IL,3).AND.NL4.EQ.IP5(IL,4).AND.
     2         NLTAU.EQ.IP5(IL,5))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM4A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM4A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE ROTMAT(XA,KXA,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(KXA)
      IF(IND.EQ.0)THEN
        READ(59)XA
      ELSE
        IF(IND.EQ.1)THEN
          WRITE(58)XA
        ELSE
          READ(58)XA
          WRITE(59)XA
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DUMVM(XK,ISIZEL,ISIZER,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XK(ISIZEL,ISIZER)
      COMMON/FILASS/IOUT,INP
      DO IX=1,ISIZER
        WRITE(IND)(XK(IY,IX),IY=1,ISIZEL)
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DIAGEL(ISIZE,KROT,XA,EVCI,NVALCF,J21,JROT,NS,ISTART,
     1KSTART,IOFF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ
      DIMENSION XA(1),EVCI(NVALCF,J21,1)
      COMMON/LANCZO/LANCZ
      COMMON/MATRIX/NVAL,NVALR,KSTEP,KSIGN
      COMMON/FILASS/IOUT,INP
      IEL=(JROT-1)/KSTEP+1
      IF(LANCZ)THEN
        KSIZE=(KROT-KSTART)*NVAL
        JSIZE=ISIZE-KSIZE-(KSTART-1)*NVAL
        L0=KSIZE*JSIZE+KSIZE*(KSIZE+1)/2
        DO ILHS=1,NVAL
          XA(ILHS+L0)=EVCI(ILHS,IEL,NS)
          L0=L0+JSIZE-ILHS
        END DO
      ELSE
C       IOFF=(KROT-1)*NVAL
        DO I=1,NVAL
          JOFF=IOFF+I-1
          J0=JOFF*(JOFF+1)/2+IOFF
          XA(J0+I)=EVCI(I,IEL,NS)
        END DO
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE OFFDEL(ISIZE,XA,CFSL,CFSR,ISIZMX,TEMP,ISIZEL,
     1ISIZER,NVALCF,YK,S,J21,KROTL,KROTR,JROTL,JROTR,IELX,IERX,
     2NSL,NSR,ISTART,KSTART,IABC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LINEAR
      DIMENSION XA(1)
      DIMENSION CFSL(ISIZMX,NVALCF,1),CFSR(ISIZMX,NVALCF,1)
      DIMENSION TEMP(ISIZMX,1)
      DIMENSION YK(ISIZMX,ISIZMX),S(J21,9,J21)
      COMMON/LANCZO/LANCZ
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10),NSYMSL(10),MSYMSL(10,100),
     1NSYMSR(10,100),MSYMSR(10,100,10)
      COMMON/FILASS/IOUT,INP
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU,IFLAUD
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/REACTL/JREACT
      COMMON/TYPE/LINEAR
      COMMON/MATRIX/NVAL,NVALR,KSTEP,KSIGN
      IF(ITIM.EQ.0)THEN
        WRITE(IOUT,*)'Calculating OFFDEL'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      KTHISR=1
      KTHISL=1
C**TEMPORARY-LIN
CC    IF(JREACT.GT.0)THEN
C     IF(JREACT.GT.0.OR.LINEAR)THEN
      IF((JREACT.GT.0.AND.IFLAUD.EQ.2).OR.LINEAR)THEN
C**TEMPORARY-LIN
        KTHISR=1+MOD(JROTR/2,2)
        KTHISL=1+MOD(JROTL/2,2)
      END IF
      IELR=(IERX-1)/KSTEP+1
      IELL=(IELX-1)/KSTEP+1
      ILOFF=NVAL*(KROTL-1)
      IROFF=NVAL*(KROTR-1)
      IF(LANCZ)THEN
        KL0=KROTL
        KSIZE=(KROTL-KSTART)*NVAL
        JSIZE=ISIZE-KSIZE-(KSTART-1)*NVAL
        L0=KSIZE*JSIZE+KSIZE*(KSIZE+1)/2
      END IF
C     IND=33
      IROT=IABC
C**READ BASIC INTEGRALS
C     REWIND IND
C**K=0/K=0
      REWIND 33
C**K=1/K=0
      REWIND 34
C**K=0/K=1
      REWIND 35
C**K=1/K=1
      REWIND 36
      LMAXR=1
      LMAXL=1
C**TEMPORARY-LIN
C     IF(JREACT.GT.0)THEN
      IF(JREACT.GT.0.OR.LINEAR)THEN
C**TEMPORARY-LIN
C       LMAXR=2
C       LMAXL=2
        LMAXR=IFLAUD
        LMAXL=IFLAUD
      END IF
      MSL=NSL
      IF(NSR.LT.NSL)MSL=NSR
      MSR=NSR
      IF(NSR.LT.NSL)MSR=NSL
      DO LTHISR=1,LMAXR
        DO LTHISL=1,LMAXL
CC        MSLL=NVSYM
CC        IF(LTHISL.EQ.KTHISL.AND.LTHISR.EQ.KTHISR)MSLL=MSL
CC        DO KSL=1,MSLL
          IF(KTHISR.EQ.LTHISR.AND.KTHISL.EQ.LTHISL)THEN
            IND=MOD(LTHISR-1,2)*LMAXR+LTHISL+32
            DO KSL=1,NVSYM
              KSIZEL=NTOT(KSL)
              IF(KSIZEL.EQ.0)GO TO 999
              DO KSR=KSL,NVSYM
                KSIZER=NTOT(KSR)
                IF(KSIZER.EQ.0)GO TO 888
C**
C**SEARCH FOR MATCH FOR SYMMETRIES NSL,NSR IN RO-VIB BLOCKS
      IGOT=0
      DO NROT=1,NRSYM
        NL=NSYMSL(NROT)
C**NO LHS SYMMETRIES
        IF(NL.EQ.0)GO TO 888
        DO NLS=1,NL
          NSSL=MSYMSL(NROT,NLS)
          NR=NSYMSR(NROT,NLS)
C**LHS SYMMETRY WRONG OR NO CORRESPONDING RHS SYMMETRIES
          IF(NSSL.NE.KSL.OR.NR.EQ.0)GO TO 5504
          DO NRS=1,NR
            NSSR=MSYMSR(NROT,NLS,NR)
C**FOUND A MATCH FOR THIS NSL,NSR
            IF(KSR.EQ.NSSR)IGOT=1
          END DO
C**FOUND MATCH - PROCESS IT
          IF(IGOT.NE.0)GO TO 5505
5504      CONTINUE
        END DO
      END DO
C**IF NO MATCH, DON'T NEED TO DO THIS NSL,NSR
      IF(IGOT.EQ.0)GO TO 888
5505  CONTINUE
C**
                IF(NSL.LE.NSR)THEN
                  DO IX=1,KSIZER
                    READ(IND)(YK(IY,IX),IY=1,KSIZEL)
                  END DO
                ELSE
                  DO IX=1,KSIZER
                    READ(IND)(YK(IX,IY),IY=1,KSIZEL)
                  END DO
                END IF
CC              IF(LTHISL.EQ.KTHISL.AND.LTHISR.EQ.KTHISR.AND.KSL.EQ.MSL.
                IF(KSL.EQ.MSL.AND.KSR.EQ.MSR)GO TO 1000
888             CONTINUE
              END DO
999           CONTINUE
            END DO
          END IF
        END DO
      END DO
1000  CONTINUE
C**RHS
C**CALL MATRIX MULT. ROUTINE MXMA TO SET UP TEMP(IY,I2)
C     CALL MXMA(YK(1,1),1,ISIZMX,CFSR(1,1,IELR),1,ISIZMX,
C    &      TEMP(1,1),1,ISIZMX,ISIZEL,ISIZER,NVAL)
      CALL DGEMM('N','N',ISIZEL,NVAL,ISIZER,1.0D0,YK(1,1),ISIZMX,
     &       CFSR(1,1,IELR),ISIZMX,0.0D0,TEMP,ISIZMX)
      DO IX=1,NVAL
        DO IY=1,NVAL
          YK(IY,IX)=0
        END DO
      END DO
C**LHS
C**CALL MXMA TO MULT. TEMP() BY LHS CFS
C     CALL MXMB(CFSL(1,1,IELL),ISIZMX,1,TEMP(1,1),1,ISIZMX,
C    &      YK(1,1),1,ISIZMX,NVAL,ISIZEL,NVAL)
      CALL DGEMM('T','N',NVAL,NVAL,ISIZEL,1.0D0,CFSL(1,1,IELL),
     &  ISIZMX,TEMP,ISIZMX,1.0D0,YK(1,1),ISIZMX)
C**MULTIPLY OFF-DIAGONAL BLOCK BY RELEVANT S-MATRIX ELEMENT
      IF(KROTL.EQ.KROTR)THEN
        IF(LANCZ)THEN
          IY=0
          JY=0
          IX0=0
          DO ILHS=1,NVAL
            IY=IY+1
            KR0=KL0
            IX=IX0
            DO IRHS=ILHS,JSIZE
              JRHS=IRHS-ILHS+1+JY
              IX=IX+1
              IF(KL0.EQ.KROTL.AND.KR0.EQ.KROTR)THEN
                XA(IRHS+L0)=XA(IRHS+L0)+YK(IY,IX)*S(JROTL,IROT,JROTR)*
     1          KSIGN
              END IF
              IF(MOD(JRHS,NVAL).EQ.0)THEN
                KR0=KR0+1
                IX=IX0
              END IF
            END DO
            L0=L0+JSIZE-ILHS
            JY=JY+1
            IX0=IX0+1
          END DO
        ELSE
          DO IX=1,NVAL
            JOFF=IROFF+IX-1
            J0=JOFF*(JOFF+1)/2+ILOFF
            DO IY=1,IX
              XA(J0+IY)=XA(J0+IY)+YK(IY,IX)*S(JROTL,IROT,JROTR)*KSIGN
            END DO
          END DO
        END IF
      ELSE
        IF(LANCZ)THEN
          IY=0
          JY=0
          DO ILHS=1,NVAL
            IY=IY+1
            KR0=KL0
            IX=0
            DO IRHS=ILHS,JSIZE
              JRHS=IRHS-ILHS+1+JY
              IX=IX+1
              IF(KL0.EQ.KROTL.AND.KR0.EQ.KROTR)THEN
                XA(IRHS+L0)=XA(IRHS+L0)+YK(IY,IX)*S(JROTL,IROT,JROTR)*
     1          KSIGN
              END IF
              IF(MOD(JRHS,NVAL).EQ.0)THEN
                KR0=KR0+1
                IX=0
              END IF
            END DO
            L0=L0+JSIZE-ILHS
            JY=JY+1
          END DO
        ELSE
          DO IX=1,NVAL
            JOFF=IROFF+IX-1
            J0=JOFF*(JOFF+1)/2+ILOFF
            DO IY=1,NVAL
              XA(J0+IY)=XA(J0+IY)+YK(IY,IX)*S(JROTL,IROT,JROTR)*KSIGN
            END DO
          END DO
        END IF
      END IF
      IF(ITIM.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM=ITIM+1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE ROTEL(MAXJ2,S,SX,JM,N,M1,M2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SX(JM,9,JM),S(JM,9,JM)
      COMMON/REACTL/JREACT,IFITRP
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      DO K=1,JM
        DO M=1,9
          DO I=1,JM
            S(I,M,K)=0
          END DO
        END DO
      END DO
C**ON ENTRY MAXJ2 = 2*J
C**         JM = 2*J + 1
      JK=MAXJ2+1
      J=MAXJ2/2
      DO I2=1,JM
        DO K=1,9
          DO I1=1,JM
            SX(I1,K,I2)=0.D0
          END DO
        END DO
      END DO
C*********************************************
C     L1=L2
      L=-J-1
      DO I=1,JK
        L=L+1
C**NB....FUNCTION FOR MZ IS EXP(+I.K.GAMMA),NB.....NOT BRINK & SATCHLER
C**Jx**2
        SX(I,1,I)=0.5D0*(J*(J+1)-L*L)
C**Jy**2
        SX(I,2,I)=SX(I,1,I)
C**Jz**2
        SX(I,3,I)=L*L
C**Jz TAKE ACCOUNT OF 'i'
        SX(I,9,I)=L
      END DO
C*********************************************
C     L1=L2+1,  L1=L2-1
      L1=-J-1
      JK1=JK-1
      DO I1=1,JK1
        I2=I1+1
        L1=L1+1
        L2=L1+1
C**Jx ALWAYS NEGATIVE
        SX(I1,7,I2)=-0.5D0*SQRT((J+L2)*(J-L2+1.D0))
        SX(I2,7,I1)=+0.5D0*SQRT((J-L1)*(J+L1+1.D0))
C**Jy INITIALLY POSITIVE
        SX(I1,8,I2)=-SX(I1,7,I2)
        SX(I2,8,I1)=SX(I2,7,I1)
C**JxJz + JzJx INITIALLY LIKE POSITIVE Jx
        SX(I1,4,I2)=-(2*L2-1)*SX(I1,7,I2)
        SX(I2,4,I1)=-(2*L1+1)*SX(I2,7,I1)
C**JyJz + JzJy ALWAYS LIKE POSITIVE Jy
        SX(I1,5,I2)=(2*L2-1)*SX(I1,8,I2)
        SX(I2,5,I1)=(2*L1+1)*SX(I2,8,I1)
      END DO
C*********************************************
C     L1=L2+2,   L1=L2-2
      L1=-J-1
      JK2=JK-2
      IF(JK2.GT.0)THEN
        DO I1=1,JK2
          I2=I1+2
          L1=L1+1
          L2=L1+2
C**Jx**2
          SX(I1,1,I2)=-0.25D0*SQRT((J-L2+1.D0)*(J-L2+2.D0)*(J+L2)*
     1    (J+L2-1.D0))
          SX(I2,1,I1)=-0.25D0*SQRT((J+L1+1.D0)*(J+L1+2.D0)*(J-L1)*
     1    (J-L1-1.D0))
C**Jy**2
          SX(I1,2,I2)=-SX(I1,1,I2)
          SX(I2,2,I1)=-SX(I2,1,I1)
C**JxJy + JyJx INITIALLY POSITIVE
          SX(I1,6,I2)=0.5D0*SQRT((J+L2-1.D0)*(J-L2+2.D0)*(J+L2)*
     1    (J-L2+1.D0))
          SX(I2,6,I1)=-0.5D0*SQRT((J-L1-1.D0)*(J+L1+2.D0)*(J-L1)*
     1    (J+L1+1.D0))
        END DO
      END IF
C***********************************************************
C***********************************************************
      DO I=1,9
C**SET ALTERNATING SIGNS FOR K=0
C**Jy AND [JxJz]+ TAKE ACCOUNT OF 'i'
        IS48=1
        IF(I.EQ.4.OR.I.EQ.8)IS48=-1
C**[JxJy]+ TAKE ACCOUNT OF 'i'
        IS6=1
        IF(I.EQ.6)IS6=-1
C**TAKE ACCOUNT OF ANTI-COMMUTATION
C       IF(JREACT.EQ.0)THEN
C       IS789=1
C       IF(I.EQ.7.OR.I.EQ.8.OR.I.EQ.9)IS789=-1
C       END IF
C*******************************
C**INDEX FOR D(0,0)J
        J0=J+1
C*******************
        I1=1
        SQ2=SQRT(2.D0)
        I2=I1
C**<D(0,0)J/O/D(0,0)J>
        S(I1,I,I2)=SX(J0,I,J0)*M1
        IF(J.EQ.0)GO TO 5000
C**<D(0,0)J/O/D(0,K)J>
        DO L=1,J
          IF(L.GT.2)GO TO 2000
          I2=I2+1
          IS=(-1)**(N+L)
          J1=J0+L
          J2=J0-L
C**Jy AND [JxJz]+ TAKE ACCOUNT OF 'i'
C**[JxJy]+ TAKE ACCOUNT OF 'i'
          S(I1,I,I2)=IS48*IS6*(SX(J0,I,J1)+SX(J0,I,J2)*IS)/SQ2
          I2=I2+1
          S(I1,I,I2)=(SX(J0,I,J1)-SX(J0,I,J2)*IS)/SQ2
2000      CONTINUE
        END DO
        I1=2
        ISTART=I1
C**<D(0,K+)J/O/D(0,K+)J>
C**<D(0,K+)J/O/D(0,K-)J>
C**<D(0,K-)J/O/D(0,K+)J>
C**<D(0,K-)J/O/D(0,K-)J>
        DO L=1,J
          I2=I1
          IS=(-1)**(N+L)
          J1=J0+L
          J2=J0-L
          S(I1,I,I2)=(M1*(SX(J1,I,J1)+SX(J2,I,J2))+IS*M2*(SX(J1,I,J2)
     1    +SX(J2,I,J1)))/2.D0
          I2=I2+1
          S(I1,I,I2)=(SX(J1,I,J1)-SX(J2,I,J2)-IS*(SX(J1,I,J2)-
     1    SX(J2,I,J1)))/2.D0
          S(I2,I,I1)=(SX(J1,I,J1)-SX(J2,I,J2)+IS*(SX(J1,I,J2)-
     1    SX(J2,I,J1)))/2.D0
          I1=I1+1
          I2=I1
          S(I1,I,I2)=(M1*(SX(J1,I,J1)+SX(J2,I,J2))-IS*M2*(SX(J1,I,J2)
     1    +SX(J2,I,J1)))/2.D0
          I1=I1+1
        END DO
        I1=ISTART
        L1=1
3000    I2=I1+2
C**<D(0,K)J/O/D(0,K+1)J>
        L2=L1+1
        IF(L2.GT.J)GO TO 5000
        IS1=(-1)**(N+L1)
        IS2=(-1)**(N+L2)
        J1=J0+L1
        J2=J0-L1
        K1=J0+L2
        K2=J0-L2
        S(I1,I,I2)=(SX(J1,I,K1)+IS1*IS2*SX(J2,I,K2)+IS2*SX(J1,I,K2)+
     1  IS1*SX(J2,I,K1))/2.D0
        I2=I2+1
        S(I1,I,I2)=(SX(J1,I,K1)-IS1*IS2*SX(J2,I,K2)-IS2*SX(J1,I,K2)+
     1  IS1*SX(J2,I,K1))/2.D0
        I1=I1+1
        I2=I2-1
C**Jy AND [JxJz]+ TAKE ACCOUNT OF 'i'
        S(I1,I,I2)=(SX(J1,I,K1)-IS1*IS2*SX(J2,I,K2)+IS2*SX(J1,I,K2)-
     1  IS1*SX(J2,I,K1))*IS48/2.D0
        I2=I2+1
        S(I1,I,I2)=(SX(J1,I,K1)+IS1*IS2*SX(J2,I,K2)-IS2*SX(J1,I,K2)-
     1  IS1*SX(J2,I,K1))/2.D0
C**<D(0,K)J/O/D(0,K+2)J>
        I1=I1-1
        I2=I2+1
        L2=L2+1
        IF(L2.GT.J)GO TO 4000
        IS2=(-1)**(N+L2)
        K1=J0+L2
        K2=J0-L2
        S(I1,I,I2)=(SX(J1,I,K1)+IS1*IS2*SX(J2,I,K2)+IS2*SX(J1,I,K2)+
     1  IS1*SX(J2,I,K1))/2.D0
        I2=I2+1
        S(I1,I,I2)=(SX(J1,I,K1)-IS1*IS2*SX(J2,I,K2)-IS2*SX(J1,I,K2)+
     1  IS1*SX(J2,I,K1))/2.D0
        I1=I1+1
        I2=I2-1
C**[JxJy]+ TAKE ACCOUNT OF 'i'
        S(I1,I,I2)=(SX(J1,I,K1)-IS1*IS2*SX(J2,I,K2)+IS2*SX(J1,I,K2)-
     1  IS1*SX(J2,I,K1))*IS6/2.D0
        I2=I2+1
        S(I1,I,I2)=(SX(J1,I,K1)+IS1*IS2*SX(J2,I,K2)-IS2*SX(J1,I,K2)-
     1  IS1*SX(J2,I,K1))/2.D0
4000    CONTINUE
        I1=I1+1
        L1=L1+1
        GO TO 3000
5000    CONTINUE
C**TAKE ACCOUNT OF ANTI-COMMUTATION
C       IF(JREACT.EQ.0)THEN
C       DO IX=1,JM
C       DO JX=1,JM
C         S(JX,I,IX)=S(JX,I,IX)*IS789
C       END DO
C       END DO
C       END IF
      END DO
C**TEMPORARY
C     IF(JREACT.NE.0)THEN
      DO IX=1,JM
      DO JX=1,JM
C**INTERCHANGE X^2 AND Y^2
        X=S(JX,1,IX)
        S(JX,1,IX)=S(JX,2,IX)
        S(JX,2,IX)=X
C**INTERCHANGE XZ AND YZ
        X=S(JX,4,IX)
        S(JX,4,IX)=S(JX,5,IX)
        S(JX,5,IX)=X
C**INTERCHANGE X AND Y
        X=S(JX,7,IX)
        S(JX,7,IX)=S(JX,8,IX)
        S(JX,8,IX)=X
      END DO
      END DO
C     END IF
C**TEMPORARY
      IF(IPRINT.GT.2)WRITE(IOUT,*)'ROTATION MATRIX'
      DO I=1,JM
      DO J=1,JM
      IF(IPRINT.GT.2)WRITE(IOUT,500)I,J,(S(I,K,J),K=1,9)
500   FORMAT(/,1X,2I3,6(5X,F12.8))
      END DO
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE OUTCIR(ISIZE,XK,X,NDUMP,MSYM,LDUMP,IDUMP,WRK,NMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**NEED TO WRITE ROT FUNCTION IN K-GROUPS (K=0,...,JTHIS) OF NVAL EACH
C**NB FOR K>0, EACH K HAS +/- COMPONENTS, UNLESS TRIATOMIC
      LOGICAL LANCZ,LANZA,LANZB,TRIAT
      DIMENSION NDUMP(IDUMP,1),XK(ISIZE,1),X(ISIZE),WRK(1)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      DIMENSION IASS(60)
      COMMON/MATRIX/NVALV,NVALR,KSTEP,KSIGN,NVALCF
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/TRIATO/TRIAT
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/EVL/EVL
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7
      REWIND 60
C     IF(TRIAT)THEN
C       IF(MSYM.LE.NVSYM)KSIZE=JTHIS+1
C       IF(MSYM.GT.NVSYM)KSIZE=JTHIS
C     ELSE
C       KSIZE=2*JTHIS+1
C     END IF
      KSIZE=ISIZE/NVALV
      WRITE(INP5)NVALV,KSIZE*NVALV,IDUMP
      DO J=1,IDUMP
        N=J
        IF(LDUMP.LT.0)THEN
          N=NDUMP(J,MSYM)
        END IF
        WRITE(INP5)N
        IF(N.GE.0)THEN
          WRITE(INP5)WRK(N)-EVL/WAVENM
1         READ(60)JX,KA,KC,(IASS(MODE),MODE=1,NMODE)
          IF(N.NE.JX)GO TO 1
          WRITE(INP5)(IASS(MODE),MODE=1,NMODE),KA,KC
          ISTART=0
          DO K=1,KSIZE
C           IF(LANZA)THEN
            IF(LANCZ)THEN
              REWIND 54
              DO M=1,N
                READ(54)X
              END DO
              WRITE(INP5)(X(I+ISTART),I=1,NVALV)
            ELSE
              WRITE(INP5)(XK(I+ISTART,N),I=1,NVALV)
            END IF
            ISTART=ISTART+NVALV
          END DO
        END IF
      END DO
      RETURN
      END

