C****************************************************************
C****************************************************************
C**CURVELINEAR
C****************************************************************
C****************************************************************
      SUBROUTINE FLIP(NATOM,NMODE,NBF,MBF,NVF,OMEGA,XL,MODINT,XLS,
     1MAXBAS,MEFF,MAXJ,ICI,IREACT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NBF(NMODE),MBF(NMODE),NVF(NMODE),MODINT(NMODE)
      DIMENSION MEFF(NMODE)
      DIMENSION XLS(NATOM,3),MAXBAS(NMODE,MAXJ),MAXS(4)
C**EQUILIBRIUM VALUES ONLY
      DIMENSION OMEGA(NMODE),XL(NATOM,NMODE,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
C**TEMPORARY (DIMENSIONS)
      COMMON/DAVE/VE,V(722),WE(12),DOMEGA(11,722)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100)
C************************
      WRITE(IOUT,*)'FLIP IREACT = ',IREACT
      N0=0
      DO N=1,NMODE-1
        N0=N0+1
        IF(N.EQ.IREACT)THEN
          N0=N0+1
C**N0 IS CURRENT LOCATION - N IS NEW LOCATION
          NBFS=NBF(IREACT)
          MBFS=MBF(IREACT)
          NVFS=NVF(IREACT)
          IF(MEFF(IREACT).GT.IREACT)THEN
            MEFFS=MEFF(IREACT)-1
          ELSE
            MEFFS=MEFF(IREACT)
          END IF
          MODNTS=MODINT(IREACT)
          OMEGAS=OMEGA(IREACT)
          WES=WE(IREACT)
          IF(ICI.LT.0)THEN
            DO I=1,MAXJ
              MAXS(I)=MAXBAS(IREACT,I)
            END DO
          END IF
          DO I=1,NATOM
            DO K=1,3
              XLS(I,K)=XL(I,IREACT,K)
            END DO
          END DO
        END IF
        NBF(N)=NBF(N0)
        MBF(N)=MBF(N0)
        NVF(N)=NVF(N0)
        IF(MEFF(N0).EQ.IREACT)THEN
          MEFF(N)=NMODE
        ELSE
          IF(MEFF(N0).GT.IREACT)THEN
            MEFF(N)=MEFF(N0)-1
          ELSE
            MEFF(N)=MEFF(N0)
          END IF
        END IF
        MODINT(N)=MODINT(N0)
        OMEGA(N)=OMEGA(N0)
        write(*,*) N, OMEGA(N)
        WE(N)=WE(N0)
        IF(ICI.LT.0)THEN
          DO I=1,MAXJ
            MAXBAS(N,I)=MAXBAS(N0,I)
          END DO
        END IF
        DO I=1,NATOM
          DO K=1,3
            XL(I,N,K)=XL(I,N0,K)
          END DO
        END DO
      END DO
      IF(IREACT.LT.NMODE)THEN
        NBF(NMODE)=NBFS
        MBF(NMODE)=MBFS
        NVF(NMODE)=NVFS
        MEFF(NMODE)=MEFFS
        MODINT(NMODE)=MODNTS
        OMEGA(NMODE)=OMEGAS
        WE(NMODE)=WES
        IF(ICI.LT.0)THEN
          DO I=1,MAXJ
            MAXBAS(NMODE,I)=MAXS(I)
          END DO
        END IF
        DO I=1,NATOM
          DO K=1,3
            XL(I,NMODE,K)=XLS(I,K)
          END DO
        END DO
C**SYMMETRY DEFINED BY EQUILIBRIUM VECTORS
        DO I=1,NWSYM
          DO J=1,NSYM(I)
            IF(IREACT.EQ.ISYM(I,J))THEN
              IR=I
              JR=J
            END IF
          END DO
        END DO
        DO N=IREACT+1,NMODE
          DO I=1,NWSYM
            DO J=1,NSYM(I)
              IF(N.EQ.ISYM(I,J))THEN
                ISYM(I,J)=N-1
                IF(I.EQ.IR.AND.J.GT.JR)ISYM(I,J-1)=ISYM(I,J)
              END IF
            END DO
          END DO
        END DO
        ISYM(IR,NSYM(IR))=NMODE
C**CONTRACTION SCHEMES DEFINED BY EQUILIBRIUM VECTORS
        DO I=1,NCONT
          DO J=1,ICONT(I)
            IF(IREACT.EQ.JCONT(I,J))THEN
              IR=I
              JR=J
            END IF
          END DO
        END DO
        DO N=IREACT+1,NMODE
          DO I=1,NCONT
            DO J=1,ICONT(I)
              IF(N.EQ.JCONT(I,J))THEN
                JCONT(I,J)=N-1
                IF(I.EQ.IR.AND.J.GT.JR)JCONT(I,J-1)=JCONT(I,J)
              END IF
            END DO
          END DO
        END DO
        JCONT(IR,ICONT(IR))=NMODE
      END IF
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,*)'RE-ORDERING OF NORMAL COORDINATES'
      WRITE(IOUT,*)
      DO I=1,NMODE
        WRITE(IOUT,*)'MODE ',I
        WRITE(IOUT,*)'NBF ',NBF(I)
        WRITE(IOUT,*)'MBF ',MBF(I)
        WRITE(IOUT,*)'NVF ',NVF(I)
        WRITE(IOUT,*)'MEFF ',MEFF(I)
        WRITE(IOUT,*)'MODINT ',MODINT(I)
        WRITE(IOUT,*)'OMEGA ',OMEGA(I)*WAVENM
        WRITE(IOUT,*)'WE ',WE(I)*WAVENM
        IF(ICI.LT.0)THEN
          WRITE(IOUT,*)'MAXBAS'
          DO K=1,MAXJ
            WRITE(IOUT,*)(MAXBAS(J,K),J=1,NMODE)
          END DO
        END IF
        WRITE(IOUT,*)'XL'
        DO J=1,NATOM
          WRITE(IOUT,*)(XL(J,I,K),K=1,3)
        END DO
      END DO
      WRITE(IOUT,*)'SYM'
      DO I=1,NWSYM
        WRITE(IOUT,*)(ISYM(I,J),J=1,NSYM(I))
      END DO
      WRITE(IOUT,*)'CONT'
      DO I=1,NCONT
        WRITE(IOUT,*)(JCONT(I,J),J=1,ICONT(I))
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE TORSPT(M,H,XQ,XJ,LV,INDT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**   GET INTEGRAL COS(L.PHI) AND SIN(L.PHI)
      LOGICAL LINEAR
      DIMENSION H(2*LV,M,3,1),XQ(M),XJ(M)
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
CCCC  COMMON/SCOORD/DSTAU(362),SPLUS(362),SMINUS(362)
      COMMON/SCOORD/DSTAU(722),SPLUS(722),SMINUS(722)
      COMMON/ROTS/JMAX
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/TYPE/LINEAR
      COMMON/VIBMOD/NSMODE,NVMODE
500   FORMAT(//,1X,50(1H*),/,
     11X,'TORSION INTEGRATION POINTS AND WEIGHTS FOR MODE ',I2,/)
501   FORMAT(8X,'TAU(RAD)',10X,'GAUSS WEIGHT (DTAU = 2*PI/M)',3X,
     1'POTENTIAL(CM-1)',/)
502   FORMAT(2X,F20.12,2X,D25.12)
      IF(INDT.EQ.0)THEN
        IF(MOD(M,2).NE.0)STOP 'TAU POINTS MUST BE EVEN'
CCCC  IF(MOD(180,M).NE.0)STOP '180/M MUST BE INTEGRAL'
        IF(MOD(360,M).NE.0)STOP '360/M MUST BE INTEGRAL'
C**ON ENTRY, INIT IS ITAU INDEX FOR ZERO DEGREES
C**FIRST POINT IS THIS PLUS 180/M
CCCC    INIT=INIT+180/M
        INIT=INIT+360/M
C**INCREMENT WILL BE TWICE THIS
CCCC    INCTAU=360/M
        INCTAU=720/M
        RETURN
      END IF
      IF(IPRINT.GT.0)THEN
        WRITE(IOUT,500)NSMODE
        WRITE(IOUT,501)
        CALL FLUSH(IOUT)
      END IF
      PI2=2*DACOS(-1.0D0)
      MM=2*M
      SP=PI2/MM
      XQ(1)=SP
      WRITE(IOUT,*)'NUMBER OF POINTS: ',M
      WRITE(IOUT,*)'FIRST POINT: ',INIT,XQ(1)*RAD,' DEGREES'
      DO I=2,M
        XQ(I)=XQ(I-1)+2*SP
      END DO
      ITAU=INIT-INCTAU
      DO J=1,M
        ITAU=ITAU+INCTAU
        IF(ITAU.GT.722)ITAU=ITAU-720
      END DO
      WRITE(IOUT,*)'LAST POINT: ',ITAU,XQ(M)*RAD,' DEGREES'
      DO I=1,M
        XJ(I)=PI2/M
      END DO
      DO J=1,M
        IF(IPRINT.GT.0)WRITE(IOUT,502)XQ(J),XJ(J)
      END DO
      PI=PI2/2
      PIQ=1.D0/SQRT(PI)
      PI2Q=1.D0/SQRT(PI2)
      L=2*LV
C**LOOP ROUND TAU
      ITAU=INIT-INCTAU
      DO 55 J=1,M
      ITAU=ITAU+INCTAU
CCCC  IF(ITAU.GT.362)ITAU=ITAU-360
      IF(ITAU.GT.722)ITAU=ITAU-720
      XJS=SQRT(XJ(J))

      H(1,J,1,1)=0
      H(1,J,2,1)=0
      H(1,J,3,1)=0
      H(2,J,1,1)=PI2Q
      H(2,J,2,1)=0.D0
      H(2,J,3,1)=0.D0
      IX=3
      LMAX=1
      IF(JMAX.GT.0)LMAX=2
      IND=0
      DO LL=LV-2+1,LV-2+LMAX
        IND=IND+1
        I=IX
        IF(LL.EQ.0)GO TO 3
        DO II=1,LL
          IF(IND.EQ.1)T=II
          IF(IND.EQ.2)T=0.5D0*(2*II-1)
          H(I,J,1,IND)=SIN(T*XQ(J))
          IF(.NOT.LINEAR)THEN
            H(I,J,2,IND)=T*COS(T*XQ(J))/DSTAU(ITAU)
          ELSE
            H(I,J,2,IND)=T*COS(T*XQ(J))
          END IF
          H(I,J,3,IND)=-T*T*SIN(T*XQ(J))
          I=I+1
          H(I,J,1,IND)=COS(T*XQ(J))
          IF(.NOT.LINEAR)THEN
            H(I,J,2,IND)=-T*SIN(T*XQ(J))/DSTAU(ITAU)
          ELSE
            H(I,J,2,IND)=-T*SIN(T*XQ(J))
          END IF
          H(I,J,3,IND)=-T*T*COS(T*XQ(J))
          I=I+1
        END DO
3       CONTINUE
        IX=1
      END DO

      IX=3
      DO IND=1,LMAX
        IT=IX
        DO I=1,L
          DO K=1,3
            H(I,J,K,IND)=H(I,J,K,IND)*XJS
          END DO
        END DO
        DO I=IT,L
          DO K=1,3
            H(I,J,K,IND)=H(I,J,K,IND)*PIQ
          END DO
        END DO
        IX=1
      END DO
55    CONTINUE
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE MILLER(NMODE,NATOM,QQ,BB,BS,BBS,B,BSS,XX,XL,XXP,
     1XM,BBT,BT,ZZ,ITAU,JMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ECKART
      DIMENSION QQ(NMODE),XM(NATOM),BBT(NMODE),BT(NMODE,NMODE)
CCCC  DIMENSION BB(NMODE,NMODE,3,362),BBS(NMODE,NMODE,362)
CCCC  DIMENSION B(NMODE,3,3,362),BS(NMODE,3,362),BSS(NMODE,362)
CCCC  DIMENSION XX(NATOM,3,362),XL(NATOM,NMODE,3,362)
      DIMENSION BB(NMODE,NMODE,3,722),BBS(NMODE,NMODE,722)
      DIMENSION B(NMODE,3,3,722),BS(NMODE,3,722),BSS(NMODE,722)
      DIMENSION XX(NATOM,3,722),XL(NATOM,NMODE,3,722)
      DIMENSION XXP(NATOM,3,722)
      DIMENSION SUP1(4,4),SUP2(4,4),SUP3(20),SUP4(4),SUP5(4)
      COMMON/FILASS/IOUT,INP
      COMMON/ECKCON/ECKART
      COMMON/MOMI/XIN(3,3),XMU(3,3)
      COMMON/MOMI0/XMU0(4,4),XIN0(4,4),XM0(4,4),SST,DETMU
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
CCCC  COMMON/SCOORD/DSTAU(362),SPLUS(362),SMINUS(362)
      COMMON/SCOORD/DSTAU(722),SPLUS(722),SMINUS(722)
      COMMON/DETS/DETM,DET(2)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
C****************************************************************
C****************************************************************
C**NEED FIRST PART 3 TIMES:
C**     AT ITAU+1
C**     AT ITAU-1
C**     AT ITAU
      SST=0
      DO IS=-1,3,2
        INCR=MOD(IS,3)+ITAU
CCCC    IF(IS.EQ.-1.AND.ITAU.EQ.2)INCR=361
CCCC    IF(IS.EQ.1.AND.ITAU.EQ.362)INCR=3
        IF(IS.EQ.-1.AND.ITAU.EQ.2)INCR=721
        IF(IS.EQ.1.AND.ITAU.EQ.722)INCR=3
C*********************************
C**SET UP MOMENT OF INERTIA MATRIX
C*********************************
        DO IX=1,4
          DO IY=1,IX
            XMU0(IY,IX)=0
            XIN0(IY,IX)=0
          END DO
        END DO
C**FIRST THE MATRIX B (XMU0)
        DO K=1,NMODE-1
          DO J=1,3
            DO I=1,J
              XMU0(I,J)=XMU0(I,J)+QQ(K)*B(K,I,J,INCR)
            END DO
            XMU0(J,4)=XMU0(J,4)+QQ(K)*BS(K,J,INCR)
          END DO
          XMU0(4,4)=XMU0(4,4)+QQ(K)*BSS(K,INCR)
        END DO
C**NOW THE MATRIX I0 (XIN0)
        DO IX=1,3
          DO I=1,NATOM
            X=0
            DO J=1,3
              IF(J.NE.IX)X=X+XX(I,J,INCR)*XX(I,J,INCR)
            END DO
            XIN0(IX,IX)=XIN0(IX,IX)+X
          END DO
        END DO
        DO IX=1,3
          DO IY=1,IX-1
            DO I=1,NATOM
              XIN0(IY,IX)=XIN0(IY,IX)-XX(I,IY,INCR)*XX(I,IX,INCR)
            END DO
          END DO
        END DO

        IF(ECKART)THEN
          XIN0(4,4)=1
        ELSE
          DO I=1,NATOM
            XIN0(1,4)=XIN0(1,4)+
     &                 XX(I,2,INCR)*XXP(I,3,INCR)-
     &                 XX(I,3,INCR)*XXP(I,2,INCR)
            XIN0(2,4)=XIN0(2,4)+
     &                 XX(I,3,INCR)*XXP(I,1,INCR)-
     &                 XX(I,1,INCR)*XXP(I,3,INCR)
            XIN0(3,4)=XIN0(3,4)+
     &                 XX(I,1,INCR)*XXP(I,2,INCR)-
     &                 XX(I,2,INCR)*XXP(I,1,INCR)
          ENDDO
          DO I=1,NATOM
            DO J=1,3
              XIN0(4,4)=XIN0(4,4)+XXP(I,J,INCR)**2
            END DO
          END DO
        END IF

        DO IX=1,4
          DO IY=1,IX
            XMU0(IX,IY)=XMU0(IY,IX)
            XIN0(IX,IY)=XIN0(IY,IX)
          END DO
        END DO
C**FORM I0 + B = M (XM0)
        DO I=1,4
          DO J=1,4
            XM0(I,J)=XIN0(I,J)+XMU0(I,J)
C**KEEP M
            IF(I.LE.3.AND.J.LE.3)XIN(I,J)=XM0(I,J)
          END DO
        END DO
C**SAVE M
        DO I=1,4
          DO J=1,4
            SUP1(I,J)=XM0(I,J)
          END DO
        END DO

        IF(.NOT.ECKART)THEN
C         IF(SUP1(4,4).LT.0.5d0*XIN0(4,4)) SUP1(4,4) = 0.5d0*XIN0(4,4)
C         IF(SUP1(4,4).GT.1.5d0*XIN0(4,4)) SUP1(4,4) = 1.5d0*XIN0(4,4)
        END IF

C**INVERT M
        IFAIL=1
        CALL MATINV(SUP1,4,4,SUP2,SUP3,IFAIL)
        IF(IFAIL.NE.0)THEN
          WRITE(IOUT,*)'ERROR IN MATINV'
          STOP 'ERROR IN MATINV'
        END IF
C**FORM MU (XMU0)
        CALL DGEMM('N','N',4,4,4,1.0D0,XIN0,
     &    4,SUP2,4,0.0D0,SUP1,4)
        CALL DGEMM('N','N',4,4,4,1.0D0,SUP2,
     &    4,SUP1,4,0.0D0,XMU0,4)
C**KEEP MU
        DO I=1,3
          DO J=1,3
            XMU(I,J)=XMU0(I,J)
          END DO
        END DO
C**FORM DETERMINANT OF MU (DETMU)
        DETMU=DET4(XMU0,4)
C**FORM d/dTAU(DETMU**1/4)  (SST)
C*********************************************************************
C** N.B.  NEGLECT TERM IN MU**(-3/4) BECAUSE IT WILL BE MADE TO CANCEL
C*********************************************************************
        IF(IS.NE.3)THEN
          IF(IS.EQ.-1.AND.ITAU.NE.2)SST=SST+IS*DETMU/4
          IF(IS.EQ.1.AND.ITAU.NE.722)SST=SST+IS*DETMU/4
        END IF
C**SET CURRENT GEOMETRY
        DO I=1,NATOM
          DO K=1,3
            XX(I,K,1)=XX(I,K,INCR)/SQRT(XM(I))
          END DO
        END DO
      END DO
      STOTAL=SPLUS(ITAU)+SMINUS(ITAU)
      IF(ITAU.EQ.2)STOTAL=SPLUS(ITAU)
      IF(ITAU.EQ.722)STOTAL=SMINUS(ITAU)
      SST=SST/STOTAL
C**AFTER LOOP, CURRENT TAU IS AT ITAU
C****************************************************************
C****************************************************************
C**FORM DETERMINANT OF M (DETM)
      DETM=DET4(XM0,4)
C**FORM COFACTOR OF M (XIN0)
      DO I=1,4
        DO J=1,4
          K=0
          DO KK=1,4
            IF(KK.NE.I)THEN
              K=K+1
              L=0
              DO LL=1,4
                IF(LL.NE.J)THEN
                  L=L+1
                  SUP2(K,L)=XM0(KK,LL)
                END IF
              END DO
            END IF
          END DO
C**FORM DETERMINANT OF M(-ROW I,-COLUMN J)
          XIN0(I,J)=((-1)**(I+J))*DET3(SUP2,4)
        END DO
      END DO
C**FORM d/dQ(DETMU**1/4)  (BBT)
      DO K=1,NMODE-1
C**FORM Bk(s)  (XM0)
        DO I=1,3
          DO J=1,3
            XM0(I,J)=B(K,I,J,ITAU)
          END DO
          XM0(I,4)=BS(K,I,ITAU)
          XM0(4,I)=XM0(I,4)
        END DO
        XM0(4,4)=BSS(K,ITAU)
        X=0
        DO I=1,4
          DO J=1,4
            X=X+XM0(I,J)*XIN0(I,J)
          END DO
        END DO
C******************************************************************
C** N.B.  NEGLECT TERM IN MU**1/4 BECAUSE IT WILL BE MADE TO CANCEL
C******************************************************************
        BBT(K)=-MULT(K)*X/(2*DETM)
      END DO
      IF(JMAX.NE.0)RETURN
C*************************
C**NOW GET WHAT'S REQUIRED
C*************************
      DO I=1,3
        SUP3(I)=0
        DO K=1,NMODE-1
          DO L=1,NMODE-1
            SUP3(I)=SUP3(I)+QQ(K)*MULT(L)*BB(K,L,I,ITAU)*BBT(L)
          END DO
        END DO
      END DO
C**CONSTANT
C**********
      ZZ=0
      DO I=1,3
        DO J=1,3
          ZZ=ZZ+SUP3(I)*XMU0(I,J)*SUP3(J)
        END DO
      END DO
      XYZ=0
      DO K=1,NMODE-1
        XYZ=XYZ+MULT(K)*BBT(K)*BBT(K)
      END DO
      ZZ=ZZ+XYZ
      Z=0
      DO K=1,NMODE-1
        DO L=1,NMODE-1
          Z=Z+QQ(K)*MULT(L)*BBS(K,L,ITAU)*BBT(L)
        END DO
      END DO
      ZZ=ZZ+Z*XMU0(4,4)*Z
      DO I=1,3
        ZZ=ZZ+SUP3(I)*XMU0(I,4)*Z+Z*XMU0(4,I)*SUP3(I)
        ZZ=ZZ-SUP3(I)*XMU0(I,4)*SST/DETMU-SST*XMU0(4,I)*SUP3(I)/DETMU
      END DO
      ZZ=ZZ-SST*XMU0(4,4)*Z/DETMU-Z*XMU0(4,4)*SST/DETMU
C*************************
      DO L=1,NMODE-1
        DO I=1,3
          SUP4(I)=0
          DO K=1,NMODE-1
            SUP4(I)=SUP4(I)+QQ(K)*MULT(L)*BB(K,L,I,ITAU)
          END DO
        END DO
        X=0
        DO K=1,NMODE-1
          X=X+QQ(K)*MULT(L)*BBS(K,L,ITAU)
        END DO
C**FIRST DERIVS
C**************
        BBT(L)=BBT(L)+X*XMU0(4,4)*Z
        DO I=1,3
          DO J=1,3
            BBT(L)=BBT(L)+SUP3(I)*XMU0(I,J)*SUP4(J)
          END DO
          BBT(L)=BBT(L)+SUP4(I)*XMU0(I,4)*Z+SUP3(I)*XMU0(I,4)*X
          BBT(L)=BBT(L)-SUP4(I)*XMU0(I,4)*SST/DETMU
        END DO
        BBT(L)=BBT(L)-SST*XMU0(4,4)*X/DETMU
C**SECOND DERIVS
C***************
        DO N=1,L
          BT(N,L)=0
        END DO
        BT(L,L)=1
C*************************
        DO N=1,L
          DO J=1,3
            SUP5(J)=0
            DO M=1,NMODE-1
              SUP5(J)=SUP5(J)+QQ(M)*MULT(N)*BB(M,N,J,ITAU)
            END DO
          END DO
          Y=0
          DO K=1,NMODE-1
            Y=Y+QQ(K)*MULT(N)*BBS(K,N,ITAU)
          END DO
          BT(N,L)=BT(N,L)+Y*XMU0(4,4)*X
          DO I=1,3
            DO J=1,3
              BT(N,L)=BT(N,L)+SUP4(I)*XMU0(I,J)*SUP5(J)
            END DO
            BT(N,L)=BT(N,L)+X*XMU0(4,I)*SUP5(I)
            BT(N,L)=BT(N,L)+Y*XMU0(I,4)*SUP4(I)
          END DO
        END DO
        BT(L,NMODE)=-XMU0(4,4)*X
        DO I=1,3
          BT(L,NMODE)=BT(L,NMODE)-SUP4(I)*XMU0(I,4)
        END DO
      END DO
C*************************
C*************************
C**CONSTANT
C**********
      ZZ=ZZ+SST*XMU0(4,4)*SST/(DETMU*DETMU)
C**FIRST DERIVS
C**************
      BBT(NMODE)=SST*XMU0(4,4)/DETMU
      BBT(NMODE)=BBT(NMODE)-XMU0(4,4)*Z
      DO I=1,3
        BBT(NMODE)=BBT(NMODE)-SUP3(I)*XMU0(I,4)
      END DO
C**SECOND DERIVS
C***************
      BT(NMODE,NMODE)=XMU0(4,4)
C**NOW GET CORRECT SIGN
      ZZ=+ZZ/2
      DO I=1,NMODE
        BBT(I)=+BBT(I)/2
        DO J=1,I
          BT(J,I)=+BT(J,I)/2
        END DO
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      DOUBLE PRECISION FUNCTION DET3(X,IX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(IX,3)
      DET3=X(1,1)*(X(2,2)*X(3,3)-X(3,2)*X(2,3))
     1    -X(1,2)*(X(2,1)*X(3,3)-X(3,1)*X(2,3))
     2    +X(1,3)*(X(2,1)*X(3,2)-X(3,1)*X(2,2))
      RETURN
      END
C****************************************************************
C****************************************************************
      DOUBLE PRECISION FUNCTION DET4(X,IX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(IX,4),Z(3,3)
      DET4=0
      I=1
      DO J=1,4
        K=0
        DO KK=1,4
          IF(KK.NE.I)THEN
            K=K+1
            L=0
            DO LL=1,4
              IF(LL.NE.J)THEN
                L=L+1
                Z(K,L)=X(KK,LL)
              END IF
            END DO
          END IF
        END DO
        DET4=DET4+((-1)**(I+J))*DET3(Z,3)*X(I,J)
      END DO
      RETURN
      END
