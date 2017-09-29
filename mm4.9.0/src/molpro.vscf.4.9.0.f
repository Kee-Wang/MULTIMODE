C****************************************************************
C****************************************************************
C**LINK TO MOLPRO
C****************************************************************
C****************************************************************
      SUBROUTINE GETQ1(XK,EVAL,SOL,N,XQ,XV,M,MODE,INDEX,XTANPM,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XK(N,N),EVAL(N,3),SOL(N,1),XQ(M),XV(M)
      DIMENSION XTANPM(1)

      DIMENSION D(N,N)
      COMMON/VMIN/VMIN,VMAX

      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
100   FORMAT(/,1X,'DGESV IFAIL = ',I3)
101   FORMAT(D20.10,5X,I3)
102   FORMAT(5X,'SHOULD NOT OCCUR')
103   FORMAT(1X,'SUMSQ = ',D20.10,/)
105   FORMAT(/,7X,'COEFFICIENT',7X,'POWER',/)
106   FORMAT(/,1X,'MODE: ',I3,/)
107   FORMAT(/,1X,'VALUE OF Q AT MINIMUM: ',D20.10,/)
108   FORMAT(/,1X,'FITTING OF POTENTIAL FOR MODE: ',I3)
109   FORMAT(//,4X,'OBSERVED',7X,'CALCULATED',/)
110   FORMAT(1X,F12.4,5X,F12.4)
111   FORMAT(1X,'SYMMETRY ',A2)
112   FORMAT(/,1X,'RE-FIT POTENTIAL FOR MODE: ',I3)
113   FORMAT(//,4X,'LIMIT(-)',7X,'LIMIT(+)',7X,'XTANH PARAM = ',F12.8,
     1/)
114   FORMAT(//,1X,'ENTERING LEAST SQUARES FIT ROUTINE',/,
     11X,'MFIT = ',I4,' NFIT = ',I4,/)
115   FORMAT(//,1X,'ENTERING EXACT FIT ROUTINE',/)
116   FORMAT(D20.10)
C*********************************************************
      IF(M.LT.4)STOP 'TOO FEW POINTS'
      IF(M.LT.N)RETURN
      WRITE(IOUT,108)MODE
      IF(MOLPRO.EQ.2)WRITE(IOUT,111)CHSYM(INDEX)
      IF(MOLPRO.EQ.0)INDEX=0
      ISTART=(M-N)/2
      NM=N
      MN=1
      IF(MODE.GT.NONLIN)MN=2
      IF(MOLPRO.EQ.2.AND.INDEX.NE.1)THEN
        NM=N/2
        MN=2
      END IF
1000  CONTINUE
      DO I=1,NM
        EVAL(I,1)=XV(ISTART+I)
        DO J=1,NM
          XK(I,J)=XQ(ISTART+I)**(MN*(J-1))
          IF(MOLPRO.EQ.2)THEN
            XK(I,J)=XTANH(XTANPM(MODE)*XQ(ISTART+I))**(MN*(J-1))
          END IF
        END DO
      END DO
C***********************************
      WRITE(IOUT,115)
      IFAIL=1
      CALL DGESV(NM,1,XK,N,SOL,EVAL,N,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C***********************************
      DO I=1,NM
        SOL(I,1)=EVAL(I,1)
      END DO
      WRITE(IOUT,105)
      WRITE(MOUT,116)XTANPM(MODE)
      WRITE(MOUT,*)NM
      DO J=1,NM
        WRITE(IOUT,101)SOL(J,1),MN*(J-1)
        WRITE(MOUT,101)SOL(J,1),MN*(J-1)
      END DO
      CALL FLUSH(IOUT)

      WRITE(IOUT,109)
      DO I=1,NM
        EVAL(I,1)=0
        DO J=1,NM
          EVAL(I,1)=EVAL(I,1)+SOL(J,1)*
     1    XTANH(XTANPM(MODE)*XQ(ISTART+I))**(MN*(J-1))
        END DO
        WRITE(IOUT,110)XV(ISTART+I)*WAVENM,EVAL(I,1)*WAVENM
      END DO
      WRITE(IOUT,113)XTANPM(MODE)
      VNEG=0
      QNEG=-1
      DO J=1,NM
C**TEMPORARY
        VNEG=VNEG+SOL(J,1)*
     1  QNEG**(MN*(J-1))
C**TEMPORARY
      END DO
      VPOS=0
      QPOS=+1
      DO J=1,NM
C**TEMPORARY
        VPOS=VPOS+SOL(J,1)*
     1  QPOS**(MN*(J-1))
C**TEMPORARY
      END DO
      WRITE(IOUT,110)WAVENM*VNEG,WAVENM*VPOS
      RETURN
C******************************
C**LEAST SQUARES
C******************************
      NFIT=NM/2
      MFIT=NM
      WRITE(IOUT,114)MFIT,NFIT
C**TEMPORARY
      WRITE(IOUT,*)'VMIN = ',VMIN
C**TEMPORARY


C**A-MATRIX
      DO I=1,NM
        EVAL(I,1)=XV(ISTART+I)
        EVAL(I,2)=EVAL(I,1)
C**WEIGHT (MINIMUM SCALED TO UNITY)
        VSCAL=WAVENM*(EVAL(I,1)-VMIN)+1
        EVAL(I,3)=1/VSCAL
        DO J=1,NFIT
          XK(I,J)=XQ(ISTART+I)**(MN*(J-1))
          IF(MOLPRO.EQ.2)THEN
            XK(I,J)=XTANH(XTANPM(MODE)*XQ(ISTART+I))**(MN*(J-1))
          END IF
        END DO
      END DO
C**D-MATRIX
      DO J=1,NFIT
        DO K=1,NFIT
          X=0
          DO I=1,MFIT
            X=X+EVAL(I,3)*XK(I,K)*XK(I,J)
          END DO
          D(K,J)=X
        END DO
      END DO
C**E-MATRIX
      DO J=1,NFIT
        X=0
        DO I=1,MFIT
          X=X+EVAL(I,3)*XK(I,J)*EVAL(I,2)
        END DO
        EVAL(J,1)=X
      END DO

C******************************
      IFAIL=1
      CALL DGESV(NFIT,1,D,N,SOL,EVAL,N,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C******************************
      DO I=1,NM
        SOL(I,1)=EVAL(I,1)
      END DO

      WRITE(IOUT,105)
      WRITE(MOUT,*)NFIT
      DO J=1,NFIT
        WRITE(IOUT,101)SOL(J,1),MN*(J-1)
        WRITE(MOUT,101)SOL(J,1),MN*(J-1)
      END DO
      CALL FLUSH(IOUT)

      WRITE(IOUT,109)

      RMS=0

      DO I=1,MFIT
        EVAL(I,1)=0
        DO J=1,NFIT
          EVAL(I,1)=EVAL(I,1)+SOL(J,1)*
     1    XTANH(XTANPM(MODE)*XQ(ISTART+I))**(MN*(J-1))
        END DO
        WRITE(IOUT,110)XV(ISTART+I)*WAVENM,EVAL(I,1)*WAVENM

C       RMS=RMS+((XV(ISTART+I)-EVAL(I,1))*EVAL(I,3))**2
        RMS=RMS+(XV(ISTART+I)-EVAL(I,1))**2

      END DO

      RMS=WAVENM*DSQRT(RMS/MFIT)
      WRITE(IOUT,*)'ROOT MEAN SQUARE (CM-1) = ',RMS

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETH1(XQ,M,XV,MODE,MSYM,XTANPM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XQ(M),XV(M,2)
      DIMENSION XTANPM(1)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FILASS/IOUT,INP
103   FORMAT(5F15.8)
106   FORMAT(I5)
C**************************************************************
      MID=M/2
      IF(MSYM.NE.1)THEN
        DO I=1,MID
          K=MID+I
          L=MID+1-I
          XV(K,1)=XV(L,1)
          XV(K,2)=-XV(L,2)
        END DO
      END IF
      WRITE(MOUT,106)M
      WRITE(MOUT,103)(XQ(I),I=1,M)
      WRITE(MOUT,103)(XV(I,1),I=1,M)
      WRITE(MOUT,103)(XV(I,2),I=1,M)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE MOLDP2(XQ1,XQ2,MM1,MM2,NMODE,NATOM,QQ,XX,X0,XL,XM,
     1MODE1,MODE2,ENERGY,MSYM1,MSYM2,XTANPM,WK,XJ1,XJ2,NP1,CP1,IP1,
     2VP1,DP1,NTOT1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4),SYMBAD
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION VP1(MMX1,NTOT1),DP1(MMX1,NTOT1)
      DIMENSION XQ1(MM1),XQ2(MM2),QQ(NMODE),ENERGY(2+MM2,2+MM1,3)
      DIMENSION XTANPM(NMODE),WK(2,NMODE),XJ1(MM1),XJ2(MM2)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION DQ(NMODE)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/INCREM/INCR
      COMMON/WHICH/IWHICH

      COMMON/VMIN/VMIN,VMAX

      COMMON/FITTER/MFIT
      COMMON/FILASS/IOUT,INP
204   FORMAT(I4)
205   FORMAT(1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2,
     1       1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2)
206   FORMAT(A2,1X,3F20.10)
208   FORMAT(3F20.10)

      VMAX=0

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
      MMM1=MM1
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
      IF(MOLPRO.LT.0)MMM1=MM1+2
C**EVEN NUMBER OF POINTS
      MID1=MMM1/2
      MMM2=MM2
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
      IF(MOLPRO.LT.0)MMM2=MM2+2
      MID2=MMM2/2

      MH1=0
      MJ1=0
C************************************************INCREMENT
      INCR=MOLINC
      INCR1=0
C************************************************INCREMENT
      INC1=1
      IF(MOLPRO.GT.0.AND.MOD(MID1,2).EQ.0)INC1=MOLINC
      DO M1=1,MMM1
        INCR1=INCR1+1
        IF(INCR1.EQ.INC1)THEN
C**NEED THIS POINT - RESET COUNT
          INCR1=0
        ELSE
C**SKIP THIS POINT
          GO TO 1000
        END IF
C**SET DEFAULT INCREMENT
        INC1=INCR
        IF(MOLPRO.LT.0)THEN
C**ONE AFTER FIRST GAUSS NEXT
          IF(M1.EQ.1)INC1=1
C**SECOND GAUSS NEXT
          IF(M1+1.EQ.MMM1)INC1=1
C**ONLY IF MID1 ODD
          IF(MOD(MID1,2).NE.0)THEN
            IF(M1.EQ.2)INC1=1
            IF(M1+2.EQ.MMM1)INC1=1
          END IF
        END IF
C**ONE AFTER MIDDLE ONE NEXT
        IF(M1.EQ.MID1)INC1=1
        IF(MOLPRO.LT.0)THEN
          IF(M1.GT.1.AND.M1.LT.MMM1)THEN
            MJ1=MJ1+1
            XJ1(MJ1)=XQ1(M1-1)
          END IF
          MX1=M1-1
        ELSE
          MX1=M1
        END IF
        IF(ISKIP1.NE.0.AND.M1.GT.MID1)GO TO 1000
        MH1=MH1+1
        QQ(MODE1)=XQ1(MX1)
        IF(MOLPRO.LT.0)THEN
          IF(M1.EQ.1)QQ(MODE1)=WK(1,MODE1)
          IF(M1.EQ.MMM1)QQ(MODE1)=WK(2,MODE1)
        END IF

        MH2=0
        MJ2=0
        INCR2=0
C************************************************INCREMENT
        INC2=1
        IF(MOLPRO.GT.0.AND.MOD(MID2,2).EQ.0)INC2=MOLINC
        DO M2=1,MMM2
          INCR2=INCR2+1
          IF(INCR2.EQ.INC2)THEN
            INCR2=0
          ELSE
            GO TO 2000
          END IF
          INC2=INCR
          IF(MOLPRO.LT.0)THEN
            IF(M2.EQ.1)INC2=1
            IF(M2+1.EQ.MMM2)INC2=1
            IF(MOD(MID2,2).NE.0)THEN
              IF(M2.EQ.2)INC2=1
              IF(M2+2.EQ.MMM2)INC2=1
            END IF
          END IF
          IF(M2.EQ.MID2)INC2=1
          IF(MOLPRO.LT.0)THEN
            IF(M2.GT.1.AND.M2.LT.MMM2)THEN
              MJ2=MJ2+1
              XJ2(MJ2)=XQ2(M2-1)
            END IF
            MX2=M2-1
          ELSE
            MX2=M2
          END IF
          IF(ISKIP2.NE.0.AND.M2.GT.MID2)GO TO 2000
          MH2=MH2+1
          QQ(MODE2)=XQ2(MX2)
          IF(MOLPRO.LT.0)THEN
            IF(M2.EQ.1)QQ(MODE2)=WK(1,MODE2)
            IF(M2.EQ.MMM2)QQ(MODE2)=WK(2,MODE2)
          END IF

          IF(IABS(MOLPRO).EQ.3)THEN
            IF(IWHICH.GT.0)THEN
C**FITTING GLOBAL POTENTIAL
C
C**FIRST GET POTENTIAL AT M1,M2
              DO I=1,NATOM
                DO K=1,3
                  XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1            SQRT(XM(I))
                  XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1            SQRT(XM(I))
                END DO
              END DO
              CALL GETPOT(VVV,NATOM,XX,ENERGY)
              IF(MOLPRO.LT.0)THEN
C**SECOND GET DERIVATIVE WRT MODE1
                IF(M1.EQ.1.OR.M1.LT.MMM1)THEN
                  QQQ=QQ(MODE1)+1.D-3
                  DO I=1,NATOM
                    DO K=1,3
                      XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQQ/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                SQRT(XM(I))
                    END DO
                  END DO
                  CALL GETPOT(VVP1,NATOM,XX,ENERGY)
                END IF
                IF(M1.EQ.1)DD1=(VVP1-VVV)*1.D3
                IF(M1.EQ.MMM1.OR.M1.GT.1)THEN
                  QQQ=QQ(MODE1)-1.D-3
                  DO I=1,NATOM
                    DO K=1,3
                      XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQQ/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                SQRT(XM(I))
                    END DO
                  END DO
                  CALL GETPOT(VVM1,NATOM,XX,ENERGY)
                END IF
                IF(M1.EQ.MMM1)DD1=(VVV-VVM1)*1.D3
                IF(M1.GT.1.AND.M1.LT.MMM1)DD1=(VVP1-VVM1)*1.D3/2
C**THIRD GET DERIVATIVE WRT MODE2
                IF(M2.EQ.1.OR.M2.LT.MMM2)THEN
                  QQQ=QQ(MODE2)+1.D-3
                  DO I=1,NATOM
                    DO K=1,3
                      XX(I,K)=X0(I,K)+XL(I,MODE2,K)*QQQ/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                SQRT(XM(I))
                    END DO
                  END DO
                  CALL GETPOT(VVP2,NATOM,XX,ENERGY)
                END IF
                IF(M2.EQ.1)DD2=(VVP2-VVV)*1.D3
                IF(M2.EQ.MMM2.OR.M2.GT.1)THEN
                  QQQ=QQ(MODE2)-1.D-3
                  DO I=1,NATOM
                    DO K=1,3
                      XX(I,K)=X0(I,K)+XL(I,MODE2,K)*QQQ/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                SQRT(XM(I))
                    END DO
                  END DO
                  CALL GETPOT(VVM2,NATOM,XX,ENERGY)
                END IF
                IF(M2.EQ.MMM2)DD2=(VVV-VVM2)*1.D3
                IF(M2.GT.1.AND.M2.LT.MMM2)DD2=(VVP2-VVM2)*1.D3/2
              END IF
            END IF

            WRITE(MOUT,204)NATOM
            WRITE(MOUT,205)MODE1,CHSYM(MSYM1),M1,MODE2,CHSYM(MSYM2),M2

            IF(IWHICH.GT.0)THEN
C**FITTING GLOBAL POTENTIAL
              IF(MOLPRO.GT.0)WRITE(MOUT,208)VVV
              IF(MOLPRO.LT.0)WRITE(MOUT,208)VVV,DD1,DD2
            ELSE
C**FITTING ABINITIO POTENTIAL
              WRITE(MOUT,*)'**********************************'
            END IF
          END IF 

          IF(MOLPRO.EQ.4)THEN
C**FIT ABINITIO POTENTIAL
            READ(MINP,204)NATOM
            READ(MINP,205)MODE1,CHSYM(MSYM1),MDUM1,MODE2,CHSYM(MSYM2),
     1      MDUM2
            READ(MINP,*)ENERGY(MH2,MH1,1)

            IF(ENERGY(MH2,MH1,1).GT.VMAX)VMAX=ENERGY(MH2,MH1,1)

          END IF
          IF(MOLPRO.EQ.-4)THEN
C**INTERPOLATE ABINITIO POTENTIAL
            READ(MINP,204)NATOM
            READ(MINP,205)MODE1,CHSYM(MSYM1),MDUM1,MODE2,CHSYM(MSYM2),
     1      MDUM2
            READ(MINP,*)ENERGY(MH2,MH1,1),ENERGY(MH2,MH1,2),
     1      ENERGY(MH2,MH1,3)
          END IF
          DO I=1,NATOM
            IF(IABS(MOLPRO).EQ.3)THEN
              DO K=1,3
                XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1          SQRT(XM(I))
                XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1          SQRT(XM(I))
              END DO
              WRITE(MOUT,206)SYMBOL(I),(XX(I,K)*BOHR,K=1,3)
            END IF
            IF(IABS(MOLPRO).EQ.4)THEN
              READ(MINP,206)SYMBAD,(XX(I,K),K=1,3)
              DO K=1,3
                XX(I,K)=XX(I,K)/BOHR
              END DO
            END IF
          END DO
C**TAKE OUT 1-DIM POTENTIALS + CONSTANT
          IF(MOLPRO.EQ.4)THEN
            CALL GETPQ1(VQ1,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1)
            ENERGY(MH2,MH1,2)=ENERGY(MH2,MH1,1)-VQ1
          END IF
          IF(MOLPRO.EQ.-4)THEN
          END IF
2000      CONTINUE
        END DO
1000    CONTINUE
      END DO
      MFIT=MH1*MH2
      IF(MOLPRO.EQ.-4)CALL GETH2(XJ1,XJ2,MM1,MM2,ENERGY,WK,MODE1,
     1MODE2,MSYM1,MSYM2)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETQ2(XK,EVAL,SOL,MM1,XQ1,MM2,XQ2,XV,MODE1,MODE2,
     1MSYM1,MSYM2,NMODE,QQ,XTANPM,WRK,MDIM,VINF,D,NP1,CP1,IP1,NTOT1,
     2MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION XK(MDIM,MDIM),EVAL(MDIM,4),SOL(MDIM,1)
      DIMENSION XQ1(MM1),XQ2(MM2),XV(2+MM2,2+MM1,2),QQ(NMODE)
      DIMENSION XTANPM(NMODE),WRK(1),VINF(NMODE)


      DIMENSION D(MDIM,MDIM)
      COMMON/VMIN/VMIN,VMAX


      COMMON/SYMMP/ISYMP(10,10)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FITTER/MFIT,MN1,MMM1,ISYM1,ISKIP1,XTAN1,MN2,MMM2,ISYM2,
     1ISKIP2,XTAN2
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/RMS/RMS
      COMMON/FILASS/IOUT,INP
100   FORMAT(/,1X,'DGESV IFAIL = ',I3)
101   FORMAT(D20.10,5X,I3,2X,I3)
102   FORMAT(//,1X,'ENTERING EXACT FIT ROUTINE',/)
103   FORMAT(//,1X,'ENTERING LEAST SQUARES FIT ROUTINE',/,
     11X,'MFIT = ',I4,' NFIT = ',I4,/)
105   FORMAT(/,7X,'COEFFICIENT',9X,'POWERS',/)
108   FORMAT(/,1X,'FITTING OF POTENTIAL FOR MODES: ',I3,2X,I3)
109   FORMAT(//,7X,'OBSERVED',7X,'CALCULATED',/)
110   FORMAT(1X,F12.4,5X,F12.4)
111   FORMAT(1X,'SYMMETRIES ',A2,2X,A2)
112   FORMAT(//,7X,'OBSERVED',7X,'CALCULATED',9X,'TOTAL',/)
113   FORMAT(1X,F12.4,5X,F12.4,5X,F12.4)
C***********************************************
      DO K=1,NMODE
        QQ(K)=0
      END DO
      WRITE(IOUT,108)MODE1,MODE2
      WRITE(IOUT,111)CHSYM(MSYM1),CHSYM(MSYM2)
      MMM1=MM1
      MMM2=MM2
C************************************************INCREMENT
      MX1=MM1/MOLINC
      MX2=MM2/MOLINC
      MMM1=MX1+MOD(MX1,MOLINC)
      MMM2=MX2+MOD(MX2,MOLINC)
      ISYM1=MSYM1
      ISYM2=MSYM2
      XTAN1=XTANPM(MODE1)
      XTAN2=XTANPM(MODE2)

      NDIM=MM1*MM2
      N=MMM1*MMM2
      N12=ISYMP(MSYM1,MSYM2)
C**SIZE OF POLYNOMIAL DEPENDS ON SYMMETRY
      NM=N
      IF(MSYM1.NE.1)NM=N/2
      IF(MSYM2.NE.1)NM=N/2
      IF(MSYM1.NE.1.AND.MSYM2.NE.1.AND.N12.NE.1)NM=N/4

C**POWERS OF POLYNOMIAL DEPENDS ON SYMMETRY
      MN1=1
      MN2=1
      IF(MSYM1.NE.1)MN1=2
      IF(MSYM2.NE.1)MN2=2
C**SPECIAL CASE B2 WITH B2, B1 WITH B1, A2 WITH A2
      IF(MSYM1.NE.1.AND.(N12.EQ.1))MN1=1

C**ONLY NEED POTENTIAL DATA IF UNIQUE SYMMETRY CONFIGURATION
      ISKIP1=0
      IF(MSYM1.NE.1)ISKIP1=1
      ISKIP2=0
      IF(MSYM2.NE.1)ISKIP2=1
C**SPECIAL CASE B2 WITH B2, B1 WITH B1, A2 WITH A2
      IF(MSYM1.NE.1.AND.(N12.EQ.1))ISKIP1=0

      IOFF2=0
      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 1000
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 2000
          I=I+1
          EVAL(I,1)=XV(MH2,MH1,1)
          EVAL(I,2)=XV(MH2,MH1,2)
C**WEIGHT (MINIMUM SCALED TO UNITY)
          VSCAL=WAVENM*(EVAL(I,1)-VMIN)+1
          EVAL(I,3)=1/VSCAL

          J=0
          DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
            IF(MSYM1.NE.1.AND.(N12.EQ.1))
     1      IOFF2=MOD(J1,2)+1-MN2
            DO J2=1,MMM2/MN2
              J=J+1
              XK(I,J)=(XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     1        (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))
            END DO
          END DO
2000      CONTINUE
        END DO
1000    CONTINUE
      END DO
C******************************
      WRITE(IOUT,102)
      IFAIL=1
      CALL DGESV(NM,1,XK,MDIM,SOL,EVAL,NM,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C******************************
      DO I=1,NM
        SOL(I,1)=EVAL(I,1)
      END DO
      WRITE(IOUT,105)
      J=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.NE.1.AND.(N12.EQ.1))
     1  IOFF2=MOD(J1,2)+1-MN2
        DO J2=1,MMM2/MN2
          J=J+1
          WRITE(IOUT,101)SOL(J,1),MN1*(J1-1),MN2*(J2-1)-IOFF2
        END DO
      END DO
      CALL FLUSH(IOUT)

      WRITE(IOUT,109)
      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 3000
        QQ(MODE1)=XQ1(MX1)
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 4000
          QQ(MODE2)=XQ2(MX2)
          I=I+1
          EVAL(I,1)=0
          J=0
          DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
            IF(MSYM1.NE.1.AND.(N12.EQ.1))
     1      IOFF2=MOD(J1,2)+1-MN2
            DO J2=1,MMM2/MN2
              J=J+1
              EVAL(I,1)=EVAL(I,1)+SOL(J,1)*
     1        (XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     2        (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))
            END DO
          END DO
          WRITE(IOUT,110)WAVENM*XV(MH2,MH1,1),WAVENM*EVAL(I,1)
4000      CONTINUE
        END DO
3000    CONTINUE
      END DO

C******************************
C**LEAST SQUARES
C******************************
      NFIT=0
      J=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.NE.1.AND.(N12.EQ.1))
     1  IOFF2=MOD(J1,2)+1-MN2
        DO J2=1,MMM2/MN2
          J=J+1
          IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0)THEN
            NFIT=NFIT+1
C**COPY ONLY RELEVANT COEFFICIENTS FOR Vij(2)
            SOL(NFIT,1)=SOL(J,1)
          END IF
        END DO
      END DO
      WRITE(IOUT,103)MFIT,NFIT

C**A-MATRIX
      I=0
      MX1=MM1/2
      IS1=-1
      INC1=0
      DO M1=1,MMM1
        MX1=MX1+INC1*IS1
        IS1=-IS1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 1111
        MX2=MM2/2
        IS2=-1
        INC2=0
        DO M2=1,MMM2
          MX2=MX2+INC2*IS2
          IS2=-IS2
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 2222
C**LEFT-HAND SIDE
          I=I+1
          J=0
          DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
            IF(MSYM1.NE.1.AND.(N12.EQ.1))
     1      IOFF2=MOD(J1,2)+1-MN2
            DO J2=1,MMM2/MN2
              IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0)THEN
C**RIGHT-HAND SIDE
                J=J+1
                XK(I,J)=(XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     1          (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))
              END IF
            END DO
          END DO
2222      CONTINUE
        END DO
1111    CONTINUE
      END DO
C**D-MATRIX
      DO J=1,NFIT
        DO K=1,NFIT
          X=0
          DO I=1,MFIT
            X=X+EVAL(I,3)*XK(I,K)*XK(I,J)
          END DO
          D(K,J)=X
        END DO
      END DO
C**E-MATRIX
      DO J=1,NFIT
        X=0
        DO I=1,MFIT
          X=X+EVAL(I,3)*XK(I,J)*EVAL(I,2)
        END DO
        EVAL(J,1)=X
      END DO

C******************************
      IFAIL=1
      CALL DGESV(NFIT,1,D,MDIM,SOL,EVAL,NFIT,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C******************************
      DO I=1,NFIT
        SOL(I,1)=EVAL(I,1)
      END DO

      WRITE(IOUT,105)
      WRITE(MOUT,*)NFIT
      NFIT=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.NE.1.AND.(N12.EQ.1))
     1  IOFF2=MOD(J1,2)+1-MN2
        DO J2=1,MMM2/MN2
          IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0)THEN
            NFIT=NFIT+1
            WRITE(IOUT,101)SOL(NFIT,1),MN1*(J1-1),MN2*(J2-1)-IOFF2
            WRITE(MOUT,101)SOL(NFIT,1),MN1*(J1-1),MN2*(J2-1)-IOFF2
          END IF
        END DO
      END DO
      WRITE(IOUT,112)
      CALL FLUSH(IOUT)

      RMS=0

      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 5000
        QQ(MODE1)=XQ1(MX1)
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 6000
          QQ(MODE2)=XQ2(MX2)
          I=I+1
          EVAL(I,1)=0
          NFIT=0
          DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
            IF(MSYM1.NE.1.AND.(N12.EQ.1))
     1      IOFF2=MOD(J1,2)+1-MN2
            DO J2=1,MMM2/MN2
              IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0)THEN
                NFIT=NFIT+1
                EVAL(I,1)=EVAL(I,1)+SOL(NFIT,1)*
     1          (XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     2          (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))
              END IF
            END DO
          END DO
C**ADD IN 1-DIM POTENTIALS + CONSTANT
          CALL GETPQ1(VQ1,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1)
          WRITE(IOUT,113)WAVENM*XV(MH2,MH1,2),WAVENM*EVAL(I,1),
     1    WAVENM*XV(MH2,MH1,1)

          RMS=RMS+(XV(MH2,MH1,2)-EVAL(I,1))**2

6000      CONTINUE
        END DO
5000    CONTINUE
      END DO

      RMS=WAVENM*DSQRT(RMS/MFIT)
      WRITE(IOUT,*)'ROOT MEAN SQUARE (CM-1) = ',RMS


      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETH2(XQ1,XQ2,MM1,MM2,ENERGY,WK,MODE1,MODE2,MSYM1,
     1MSYM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XQ1(MM1),XQ2(MM2),ENERGY(2+MM2,2+MM1,3),WK(2,1)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/INCREM/INCR
      COMMON/FILASS/IOUT,INP
103   FORMAT(5F15.8)
106   FORMAT(2I5)
C******************************
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS
      MMM1=MM1+2
C**EVEN NUMBER OF POINTS
      MID1=MMM1/2
      MMM2=MM2+2
      MID2=MMM2/2
C************************************************INCREMENT
C**NUMBER HEG EACH SIDE
      M1HALF=MM1/2
C**RESET MID1
      MID1=M1HALF/INCR + INCR
C**RESET MMM1 
      MMM1=MID1*2
      M2HALF=MM2/2
      MID2=M2HALF/INCR + INCR
      MMM2=MID2*2
C************************************************INCREMENT
      MH1=MMM1-2
      MX1=MMM1
      IF(MSYM1.NE.1)THEN
        MX1=MID1
        IF(MSYM1.EQ.MSYM2)MX1=MMM1
      END IF
      MH2=MMM2-2
      MX2=MMM2
      IF(MSYM2.NE.1)MX2=MID2
C**
      IF(MX1.EQ.MID1)THEN
        DO M1=1,MID1
          K1=MID1+M1
          L1=MID1+1-M1
          IF(MX2.EQ.MID2)THEN
C**B2 WITH B1, B2 WITH A2, B1 WITH A2
            DO M2=1,MID2
              ENERGY(M2,K1,1)=ENERGY(M2,L1,1)
              ENERGY(M2,K1,2)=-ENERGY(M2,L1,2)
              ENERGY(M2,K1,3)=ENERGY(M2,L1,3)
              K2=MID2+M2
              L2=MID2+1-M2
              ENERGY(K2,M1,1)=ENERGY(L2,M1,1)
              ENERGY(K2,M1,2)=ENERGY(L2,M1,2)
              ENERGY(K2,M1,3)=-ENERGY(L2,M1,3)
              ENERGY(K2,K1,1)=ENERGY(L2,L1,1)
              ENERGY(K2,K1,2)=-ENERGY(L2,L1,2)
              ENERGY(K2,K1,3)=-ENERGY(L2,L1,3)
            END DO
          ELSE
C**B2 WITH A1, B1 WITH A1, A2 WITH A1
            DO M2=1,MMM2
              ENERGY(M2,K1,1)=ENERGY(M2,L1,1)
              ENERGY(M2,K1,2)=-ENERGY(M2,L1,2)
              ENERGY(M2,K1,3)=ENERGY(M2,L1,3)
            END DO
          END IF
        END DO
      ELSE
        DO M1=1,MMM1
          IF(MX2.EQ.MID2)THEN
            IF(MSYM1.EQ.MSYM2)THEN
C**B2 WITH B2, B1 WITH B1, A2 WITH A2
              K1=MMM1+1-M1
              DO M2=1,MID2
                K2=MMM2+1-M2
                ENERGY(K2,K1,1)=ENERGY(M2,M1,1)
                ENERGY(K2,K1,2)=-ENERGY(M2,M1,2)
                ENERGY(K2,K1,3)=-ENERGY(M2,M1,3)
              END DO
            ELSE
C**A1 WITH B2, A1 WITH B1, A1 WITH A2
              DO M2=1,MID2
                K2=MID2+M2
                L2=MID2+1-M2
                ENERGY(K2,M1,1)=ENERGY(L2,M1,1)
                ENERGY(K2,M1,2)=ENERGY(L2,M1,2)
                ENERGY(K2,M1,3)=-ENERGY(L2,M1,3)
              END DO
            END IF
          END IF
        END DO
      END IF
      WRITE(MOUT,106)MMM1,MMM2
      WRITE(MOUT,103)WK(1,MODE1),(XQ1(M),M=1,MH1),WK(2,MODE1)
      WRITE(MOUT,103)WK(1,MODE2),(XQ2(M),M=1,MH2),WK(2,MODE2)
C**ENERGIES AT M2 FOR EACH M1
      DO M1=1,MMM1
        WRITE(MOUT,103)(ENERGY(M2,M1,1),M2=1,MMM2)
      END DO
C**DERIVATIVES OF M1 FOR EACH M2 (MODE1 > MODE2)
      DO M2=1,MMM2
        WRITE(MOUT,103)(ENERGY(M2,M1,2),M1=1,MMM1)
      END DO
C**DERIVATIVES OF M2 FOR EACH M1 (MODE2 < MODE1)
      DO M1=1,MMM1
        WRITE(MOUT,103)(ENERGY(M2,M1,3),M2=1,MMM2)
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE MOLDP3(XQ1,XQ2,XQ3,MM1,MM2,MM3,NMODE,NATOM,QQ,XX,X0,
     1XL,XM,MODE1,MODE2,MODE3,ENERGY,MSYM1,MSYM2,MSYM3,XTANPM,WK,XJ1,
     2XJ2,XJ3,NP1,CP1,IP1,VP1,DP1,NTOT1,MMX1,NP2,CP2,IP2,VP2,DP2A,DP2B,
     3NTOT2,
     3MAX2,INDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4),SYMBAD
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),QQ(NMODE)
      DIMENSION XJ1(MM1),XJ2(MM2),XJ3(MM3)
      DIMENSION ENERGY(2+MM3,2+MM2,2+MM1,4)
      DIMENSION XTANPM(NMODE),WK(2,NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION VP1(MMX1,NTOT1),DP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION VP2(MAX2,MAX2,NTOT2),DP2A(MAX2,MAX2,NTOT2)
      DIMENSION DP2B(MAX2,MAX2,NTOT2),INDK(1),DQ(NMODE)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/INCREM/INCR
      COMMON/WHICH/IWHICH
      COMMON/FITTER/MFIT
      COMMON/FILASS/IOUT,INP
204   FORMAT(I4)
205   FORMAT(1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2,
     1       1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2,
     2       1X,'MODE:',I2,' SYMMETRY:',A2,' POINT:',I2)
206   FORMAT(A2,1X,3F20.10)
208   FORMAT(4F20.10)

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
      MMM1=MM1
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
      IF(MOLPRO.LT.0)MMM1=MM1+2
C**EVEN NUMBER OF POINTS
      MID1=MMM1/2
      MMM2=MM2
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
      IF(MOLPRO.LT.0)MMM2=MM2+2
      MID2=MMM2/2
      MMM3=MM3
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
      IF(MOLPRO.LT.0)MMM3=MM3+2
      MID3=MMM3/2

      MH1=0
      MJ1=0
C************************************************INCREMENT
      INCR=MOLINC
      INCR1=0
C************************************************INCREMENT
      INC1=1
      IF(MOLPRO.GT.0.AND.MOD(MID1,2).EQ.0)INC1=MOLINC
      DO M1=1,MMM1
        INCR1=INCR1+1
        IF(INCR1.EQ.INC1)THEN
C**NEED THIS POINT - RESET COUNT
          INCR1=0
        ELSE
C**SKIP THIS POINT
          GO TO 1000
        END IF
C**SET DEFAULT INCREMENT
        INC1=INCR
        IF(MOLPRO.LT.0)THEN
C**ONE AFTER FIRST GAUSS NEXT
          IF(M1.EQ.1)INC1=1
C**SECOND GAUSS NEXT
          IF(M1+1.EQ.MMM1)INC1=1
C**ONLY IF MID1 ODD
          IF(MOD(MID1,2).NE.0)THEN
            IF(M1.EQ.2)INC1=1
            IF(M1+2.EQ.MMM1)INC1=1
          END IF
        END IF
C**ONE AFTER MIDDLE ONE NEXT
        IF(M1.EQ.MID1)INC1=1
        IF(MOLPRO.LT.0)THEN
          IF(M1.GT.1.AND.M1.LT.MMM1)THEN
            MJ1=MJ1+1
            XJ1(MJ1)=XQ1(M1-1)
          END IF
          MX1=M1-1
        ELSE
          MX1=M1
        END IF
        IF(ISKIP1.NE.0.AND.M1.GT.MID1)GO TO 1000
        MH1=MH1+1
        QQ(MODE1)=XQ1(MX1)
        IF(MOLPRO.LT.0)THEN
          IF(M1.EQ.1)QQ(MODE1)=WK(1,MODE1)
          IF(M1.EQ.MMM1)QQ(MODE1)=WK(2,MODE1)
        END IF

        MH2=0
        MJ2=0
        INCR2=0
C************************************************INCREMENT
        INC2=1
        IF(MOLPRO.GT.0.AND.MOD(MID2,2).EQ.0)INC2=MOLINC
        DO M2=1,MMM2
          INCR2=INCR2+1
          IF(INCR2.EQ.INC2)THEN
            INCR2=0
          ELSE
            GO TO 2000
          END IF
          INC2=INCR
          IF(MOLPRO.LT.0)THEN
            IF(M2.EQ.1)INC2=1
            IF(M2+1.EQ.MMM2)INC2=1
            IF(MOD(MID2,2).NE.0)THEN
              IF(M2.EQ.2)INC2=1
              IF(M2+2.EQ.MMM2)INC2=1
            END IF
          END IF
          IF(M2.EQ.MID2)INC2=1
          IF(MOLPRO.LT.0)THEN
            IF(M2.GT.1.AND.M2.LT.MMM2)THEN
              MJ2=MJ2+1
              XJ2(MJ2)=XQ2(M2-1)
            END IF
            MX2=M2-1
          ELSE
            MX2=M2
          END IF
          IF(ISKIP2.NE.0.AND.M2.GT.MID2)GO TO 2000
          MH2=MH2+1
          QQ(MODE2)=XQ2(MX2)
          IF(MOLPRO.LT.0)THEN
            IF(M2.EQ.1)QQ(MODE2)=WK(1,MODE2)
            IF(M2.EQ.MMM2)QQ(MODE2)=WK(2,MODE2)
          END IF

          MH3=0
          MJ3=0
          INCR3=0
C************************************************INCREMENT
          INC3=1
          IF(MOLPRO.GT.0.AND.MOD(MID3,2).EQ.0)INC3=MOLINC
          DO M3=1,MMM3
            INCR3=INCR3+1
            IF(INCR3.EQ.INC3)THEN
              INCR3=0
            ELSE
              GO TO 3000
            END IF
            INC3=INCR
            IF(MOLPRO.LT.0)THEN
              IF(M3.EQ.1)INC3=1
              IF(M3+1.EQ.MMM3)INC3=1
              IF(MOD(MID3,2).NE.0)THEN
                IF(M3.EQ.2)INC3=1
                IF(M3+2.EQ.MMM3)INC3=1
              END IF
            END IF
            IF(M3.EQ.MID3)INC3=1
            IF(MOLPRO.LT.0)THEN
              IF(M3.GT.1.AND.M3.LT.MMM3)THEN
                MJ3=MJ3+1
                XJ3(MJ3)=XQ3(M3-1)
              END IF
              MX3=M3-1
            ELSE
              MX3=M3
            END IF
            IF(ISKIP3.NE.0.AND.M3.GT.MID3)GO TO 3000
            MH3=MH3+1
            QQ(MODE3)=XQ3(MX3)
            IF(MOLPRO.LT.0)THEN
              IF(M3.EQ.1)QQ(MODE3)=WK(1,MODE3)
              IF(M3.EQ.MMM3)QQ(MODE3)=WK(2,MODE3)
            END IF

            IF(IABS(MOLPRO).EQ.5)THEN

              IF(IWHICH.GT.0)THEN
C**FITTING GLOBAL POTENTIAL
C
C**FIRST GET POTENTIAL AT M1,M2,M3
                DO I=1,NATOM
                  DO K=1,3
                    XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1              SQRT(XM(I))
                    XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1              SQRT(XM(I))
                    XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1              SQRT(XM(I))
                  END DO
                END DO
                CALL GETPOT(VVV,NATOM,XX,ENERGY)
                IF(MOLPRO.LT.0)THEN
C**SECOND GET DERIVATIVE WRT MODE1
                  IF(M1.EQ.1.OR.M1.LT.MMM1)THEN
                    QQQ=QQ(MODE1)+1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVP1,NATOM,XX,ENERGY)
                  END IF
                  IF(M1.EQ.1)DD1=(VVP1-VVV)*1.D3
                  IF(M1.EQ.MMM1.OR.M1.GT.1)THEN
                    QQQ=QQ(MODE1)-1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVM1,NATOM,XX,ENERGY)
                  END IF
                  IF(M1.EQ.MMM1)DD1=(VVV-VVM1)*1.D3
                  IF(M1.GT.1.AND.M1.LT.MMM1)DD1=(VVP1-VVM1)*1.D3/2
C**THIRD GET DERIVATIVE WRT MODE2
                  IF(M2.EQ.1.OR.M2.LT.MMM2)THEN
                    QQQ=QQ(MODE2)+1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE2,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVP2,NATOM,XX,ENERGY)
                  END IF
                  IF(M2.EQ.1)DD2=(VVP2-VVV)*1.D3
                  IF(M2.EQ.MMM2.OR.M2.GT.1)THEN
                    QQQ=QQ(MODE2)-1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE2,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVM2,NATOM,XX,ENERGY)
                  END IF
                  IF(M2.EQ.MMM2)DD2=(VVV-VVM2)*1.D3
                  IF(M2.GT.1.AND.M2.LT.MMM2)DD2=(VVP2-VVM2)*1.D3/2
C**FOURTH GET DERIVATIVE WRT MODE3
                  IF(M3.EQ.1.OR.M3.LT.MMM3)THEN
                    QQQ=QQ(MODE3)+1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE3,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVP3,NATOM,XX,ENERGY)
                  END IF
                  IF(M3.EQ.1)DD3=(VVP3-VVV)*1.D3
                  IF(M3.EQ.MMM3.OR.M3.GT.1)THEN
                    QQQ=QQ(MODE3)-1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE3,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVM3,NATOM,XX,ENERGY)
                  END IF
                  IF(M3.EQ.MMM3)DD3=(VVV-VVM3)*1.D3
                  IF(M3.GT.1.AND.M3.LT.MMM3)DD3=(VVP3-VVM3)*1.D3/2
                END IF
C**TO HERE IS FUDGE FOR 3-DIM ABINITIO PACKAGE
              END IF

              WRITE(MOUT,204)NATOM
              WRITE(MOUT,205)MODE1,CHSYM(MSYM1),M1,MODE2,CHSYM(MSYM2),
     1        M2,MODE3,CHSYM(MSYM3),M3

              IF(IWHICH.GT.0)THEN
C**FITTING GLOBAL POTENTIAL
                IF(MOLPRO.GT.0)WRITE(MOUT,208)VVV
                IF(MOLPRO.LT.0)WRITE(MOUT,208)VVV,DD1,DD2,DD3
              ELSE
C**FITTING ABINITIO POTENTIAL
                WRITE(MOUT,*)'**********************************'
              END IF

            END IF 

            IF(MOLPRO.EQ.6)THEN
C**FIT ABINITIO POTENTIAL
              READ(MINP,204)NATOM
              READ(MINP,205)MODE1,CHSYM(MSYM1),MD1,MODE2,CHSYM(MSYM2),
     1        MD2,MODE3,CHSYM(MSYM3),MD3
              READ(MINP,*)ENERGY(MH3,MH2,MH1,1)
            END IF
            IF(MOLPRO.EQ.-6)THEN
C**INTERPOLATE ABINITIO POTENTIAL
              READ(MINP,204)NATOM
              READ(MINP,205)MODE1,CHSYM(MSYM1),MD1,MODE2,CHSYM(MSYM2),
     1        MD2,MODE3,CHSYM(MSYM3),MD3
              READ(MINP,*)ENERGY(MH3,MH2,MH1,1),
     1        ENERGY(MH3,MH2,MH1,2),ENERGY(MH3,MH2,MH1,3),
     2        ENERGY(MH3,MH2,MH1,4)
            END IF
            DO I=1,NATOM
              IF(IABS(MOLPRO).EQ.5)THEN
                DO K=1,3
                  XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1            SQRT(XM(I))
                  XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1            SQRT(XM(I))
                  XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1            SQRT(XM(I))
                END DO
                WRITE(MOUT,206)SYMBOL(I),(XX(I,K)*BOHR,K=1,3)
              END IF
              IF(IABS(MOLPRO).EQ.6)THEN
                READ(MINP,206)SYMBAD,(XX(I,K),K=1,3)
                DO K=1,3
                  XX(I,K)=XX(I,K)/BOHR
                END DO
              END IF
            END DO
C**TAKE OUT 1-DIM, 2-DIM POTENTIALS + CONSTANT
            IF(MOLPRO.EQ.6)THEN
              CALL GETPQ2(VQ2,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1        NP2,CP2,IP2,NTOT2,MAX2,INDK)
              ENERGY(MH3,MH2,MH1,2)=ENERGY(MH3,MH2,MH1,1)-VQ2
            END IF
            IF(MOLPRO.EQ.-6)THEN
            END IF
3000        CONTINUE
          END DO
2000      CONTINUE
        END DO
1000    CONTINUE
      END DO
      MFIT=MH1*MH2*MH3
      IF(MOLPRO.EQ.-6)CALL GETH3(XJ1,XJ2,XJ3,MM1,MM2,MM3,ENERGY,
     1WK,MODE1,MODE2,MODE3,MSYM1,MSYM2,MSYM3)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETQ3(XK,EVAL,SOL,MM1,XQ1,MM2,XQ2,MM3,XQ3,XV,MODE1,
     1MODE2,MODE3,MSYM1,MSYM2,MSYM3,NMODE,QQ,XTANPM,WRK,MDIM,VINF,D,
     2NP1,CP1,IP1,NTOT1,MMX1,NP2,CP2,IP2,NTOT2,MAX2,INDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XK(MDIM,MDIM),EVAL(MDIM,5)
      DIMENSION SOL(MDIM,1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XV(2+MM3,2+MM2,2+MM1,2)
      DIMENSION QQ(NMODE),VINF(NMODE)
      DIMENSION XTANPM(NMODE),WRK(1)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)


      DIMENSION D(MDIM,MDIM)
      COMMON/VMIN/VMIN,VMAX


      COMMON/SYMMP/ISYMP(10,10)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FITTER/MFIT,MN1,MMM1,ISYM1,ISKIP1,XTAN1,MN2,MMM2,ISYM2,
     1ISKIP2,XTAN2,MN3,MMM3,ISYM3,ISKIP3,XTAN3
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/RMS/RMS
      COMMON/FILASS/IOUT,INP
      EXTERNAL FPRINT,FUNC3
100   FORMAT(/,1X,'DGESV IFAIL = ',I3)
101   FORMAT(D20.10,5X,I3,2X,I3,2X,I3)
102   FORMAT(//,1X,'ENTERING EXACT FIT ROUTINE',/)
103   FORMAT(//,1X,'ENTERING LEAST SQUARES FIT ROUTINE',/,
     11X,'MFIT = ',I4,' NFIT = ',I4,/)
105   FORMAT(/,7X,'COEFFICIENT',9X,'POWERS',/)
108   FORMAT(/,1X,'FITTING OF POTENTIAL FOR MODES: ',I3,2X,I3,2X,I3)
109   FORMAT(//,7X,'OBSERVED',7X,'CALCULATED',/)
110   FORMAT(1X,F12.4,5X,F12.4)
111   FORMAT(1X,'SYMMETRIES ',A2,2X,A2,2X,A2)
112   FORMAT(//,7X,'OBSERVED',7X,'CALCULATED',9X,'TOTAL',/)
113   FORMAT(1X,F12.4,5X,F12.4,5X,F12.4)
C***********************************************

CCCC  FACT=WAVENM*(VMAX-VMIN)/10

      DO K=1,NMODE
        QQ(K)=0
      END DO
      WRITE(IOUT,108)MODE1,MODE2,MODE3
      WRITE(IOUT,111)CHSYM(MSYM1),CHSYM(MSYM2),CHSYM(MSYM3)
      MMM1=MM1
      MMM2=MM2
      MMM3=MM3
C************************************************INCREMENT
      MX1=MM1/MOLINC
      MX2=MM2/MOLINC
      MX3=MM3/MOLINC
      MMM1=MX1+MOD(MX1,MOLINC)
      MMM2=MX2+MOD(MX2,MOLINC)
      MMM3=MX3+MOD(MX3,MOLINC)
      ISYM1=MSYM1
      ISYM2=MSYM2
      ISYM3=MSYM3
      N12=ISYMP(MSYM1,MSYM2)
      N13=ISYMP(MSYM1,MSYM3)
      N23=ISYMP(MSYM2,MSYM3)
      XTAN1=XTANPM(MODE1)
      XTAN2=XTANPM(MODE2)
      XTAN3=XTANPM(MODE3)

      NDIM=MM1*MM2*MM3
      N=MMM1*MMM2*MMM3
C**SIZE OF POLYNOMIAL DEPENDS ON SYMMETRY
      NM=N
      IF(MSYM1.NE.1)NM=N/2
      IF(MSYM2.NE.1)NM=N/2
      IF(MSYM3.NE.1)NM=N/2
      IF(MSYM1.NE.1.AND.MSYM2.NE.1.AND.N12.NE.1)NM=N/4
      IF(MSYM1.NE.1.AND.MSYM3.NE.1.AND.N13.NE.1)NM=N/4
      IF(MSYM2.NE.1.AND.MSYM3.NE.1.AND.N23.NE.1)NM=N/4

C**POWERS OF POLYNOMIAL DEPENDS ON SYMMETRY
      MN1=1
      MN2=1
      MN3=1
      IF(MSYM1.NE.1)MN1=2
      IF(MSYM2.NE.1)MN2=2
      IF(MSYM3.NE.1)MN3=2
C**SPECIAL CASE B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2, ETC.
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM2))MN1=1
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM3))MN1=1
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM3))MN2=1

C**SPECIAL CASE S2*S3 = S1
      N23=ISYMP(MSYM2,MSYM3)
      IF(N23.EQ.MSYM1)MN1=1

C**ONLY NEED POTENTIAL DATA IF UNIQUE SYMMETRY CONFIGURATION
      ISKIP1=0
      IF(MSYM1.NE.1)ISKIP1=1
      ISKIP2=0
      IF(MSYM2.NE.1)ISKIP2=1
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

      IOFF3=0
      IOFF2=0
      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 1000
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 2000
          MH3=MMM3/2
          MX3=MM3/2
          IS3=-1
          INH3=0
          INC3=0
          DO M3=1,MMM3
            MH3=MH3+INH3*IS3
            MX3=MX3+INC3*IS3
            IS3=-IS3
            INH3=INH3+1
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 3000
            I=I+1
            EVAL(I,1)=XV(MH3,MH2,MH1,1)
            EVAL(I,2)=XV(MH3,MH2,MH1,2)
C**WEIGHT (MINIMUM SCALED TO UNITY)
            VSCAL=WAVENM*(EVAL(I,1)-VMIN)+1
            EVAL(I,3)=1/VSCAL

CCCC        TERM=WAVENM*(EVAL(I,1)-VMIN)/FACT+1
CCCC        EVAL(I,3)=EVAL(I,3)*TERM

            J=0
            DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
              IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1        IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
              IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1        IOFF3=MOD(J1,2)+1-MN3
              DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     2          IOFF3=MOD(J1+J2,2)-MN3
                DO J3=1,MMM3/MN3
                 J=J+1
                 XK(I,J)=(XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     1           (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     2           (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))
                END DO
              END DO
            END DO
3000        CONTINUE
          END DO
2000      CONTINUE
        END DO
1000    CONTINUE
      END DO
C******************************
      WRITE(IOUT,102)
      IFAIL=1
      CALL DGESV(NM,1,XK,MDIM,SOL,EVAL,NM,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C******************************
      DO I=1,NM
        SOL(I,1)=EVAL(I,1)
      END DO
      WRITE(IOUT,105)
      J=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1  IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1  IOFF3=MOD(J1,2)+1-MN3
        DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1    IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     2    IOFF3=MOD(J1+J2,2)-MN3
          DO J3=1,MMM3/MN3
            J=J+1
            WRITE(IOUT,101)SOL(J,1),MN1*(J1-1),MN2*(J2-1)-IOFF2,
     1      MN3*(J3-1)-IOFF3
          END DO
        END DO
      END DO
      CALL FLUSH(IOUT)

      WRITE(IOUT,109)
      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 4000
        QQ(MODE1)=XQ1(MX1)
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 5000
          QQ(MODE2)=XQ2(MX2)
          MH3=MMM3/2
          MX3=MM3/2
          IS3=-1
          INH3=0
          INC3=0
          DO M3=1,MMM3
            MH3=MH3+INH3*IS3
            MX3=MX3+INC3*IS3
            IS3=-IS3
            INH3=INH3+1
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 6000
            QQ(MODE3)=XQ3(MX3)
            I=I+1
            EVAL(I,1)=0
            J=0
            DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
              IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1        IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
              IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1        IOFF3=MOD(J1,2)+1-MN3
              DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     2          IOFF3=MOD(J1+J2,2)-MN3
                DO J3=1,MMM3/MN3
                  J=J+1
                  EVAL(I,1)=EVAL(I,1)+SOL(J,1)*
     1            (XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     2            (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     3            (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))
                END DO
              END DO
            END DO
            WRITE(IOUT,110)WAVENM*XV(MH3,MH2,MH1,1),WAVENM*EVAL(I,1)
6000        CONTINUE
          END DO
5000      CONTINUE
        END DO
4000    CONTINUE
      END DO

C******************************
C**LEAST SQUARES
C******************************
      NFIT=0
      J=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1  IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1  IOFF3=MOD(J1,2)+1-MN3
        DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1    IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     2    IOFF3=MOD(J1+J2,2)-MN3
          DO J3=1,MMM3/MN3
            J=J+1
            IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1      MN3*(J3-1)-IOFF3.NE.0)THEN
              NFIT=NFIT+1
C**COPY ONLY RELEVANT COEFFICIENTS FOR Vijk(3)
              SOL(NFIT,1)=SOL(J,1)
            END IF
          END DO
        END DO
      END DO
      WRITE(IOUT,103)MFIT,NFIT


C**A-MATRIX
      I=0
      MX1=MM1/2
      IS1=-1
      INC1=0
      DO M1=1,MMM1
        MX1=MX1+INC1*IS1
        IS1=-IS1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 1111
        MX2=MM2/2
        IS2=-1
        INC2=0
        DO M2=1,MMM2
          MX2=MX2+INC2*IS2
          IS2=-IS2
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 2222
          MX3=MM3/2
          IS3=-1
          INC3=0
          DO M3=1,MMM3
            MX3=MX3+INC3*IS3
            IS3=-IS3
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 3333
C**LEFT-HAND SIDE
            I=I+1
            J=0
            DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
              IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1        IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
              IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1        IOFF3=MOD(J1,2)+1-MN3
              DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     2          IOFF3=MOD(J1+J2,2)-MN3
                DO J3=1,MMM3/MN3
                  IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1            MN3*(J3-1)-IOFF3.NE.0)THEN
                 J=J+1
             XK(I,J)=(XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     1       (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     2       (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))
                  END IF
                END DO
              END DO
            END DO
3333        CONTINUE
          END DO
2222      CONTINUE
        END DO
1111    CONTINUE
      END DO
C**D-MATRIX
      DO J=1,NFIT
        DO K=1,NFIT
          X=0
          DO I=1,MFIT
            X=X+EVAL(I,3)*XK(I,K)*XK(I,J)
          END DO
          D(K,J)=X
        END DO
      END DO
C**E-MATRIX
      DO J=1,NFIT
        X=0
        DO I=1,MFIT
          X=X+EVAL(I,3)*XK(I,J)*EVAL(I,2)
        END DO
        EVAL(J,1)=X
      END DO

C******************************
      IFAIL=1
      CALL DGESV(NFIT,1,D,MDIM,SOL,EVAL,NFIT,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C******************************
      DO I=1,NFIT
        SOL(I,1)=EVAL(I,1)
      END DO

      WRITE(IOUT,105)
      WRITE(MOUT,*)NFIT
      NFIT=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1  IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1  IOFF3=MOD(J1,2)+1-MN3
        DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1    IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     2    IOFF3=MOD(J1+J2,2)-MN3
          DO J3=1,MMM3/MN3
            IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1      MN3*(J3-1)-IOFF3.NE.0)THEN
              NFIT=NFIT+1
              WRITE(IOUT,101)SOL(NFIT,1),MN1*(J1-1),MN2*(J2-1)-IOFF2,
     1        MN3*(J3-1)-IOFF3
              WRITE(MOUT,101)SOL(NFIT,1),MN1*(J1-1),MN2*(J2-1)-IOFF2,
     1        MN3*(J3-1)-IOFF3
            END IF
          END DO
        END DO
      END DO
      CALL FLUSH(IOUT)

      WRITE(IOUT,112)


      RMS=0


      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 7000
        QQ(MODE1)=XQ1(MX1)
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 8000
          QQ(MODE2)=XQ2(MX2)
          MH3=MMM3/2
          MX3=MM3/2
          IS3=-1
          INH3=0
          INC3=0
          DO M3=1,MMM3
            MH3=MH3+INH3*IS3
            MX3=MX3+INC3*IS3
            IS3=-IS3
            INH3=INH3+1
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 9000
            QQ(MODE3)=XQ3(MX3)
            I=I+1
            EVAL(I,1)=0
            NFIT=0
            DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
              IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1        IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
              IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1        IOFF3=MOD(J1,2)+1-MN3
              DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     2          IOFF3=MOD(J1+J2,2)-MN3
                DO J3=1,MMM3/MN3
                  IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1            MN3*(J3-1)-IOFF3.NE.0)THEN
                    NFIT=NFIT+1
                    EVAL(I,1)=EVAL(I,1)+SOL(NFIT,1)*
     1           (XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     2           (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     3           (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))
                  END IF
                END DO
              END DO
            END DO
C**ADD IN 1-DIM, 2-DIM POTENTIALS + CONSTANT
            CALL GETPQ2(VQ2,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1      NP2,CP2,IP2,NTOT2,MAX2,INDK)
            WRITE(IOUT,113)WAVENM*XV(MH3,MH2,MH1,2),WAVENM*EVAL(I,1),
     1      WAVENM*XV(MH3,MH2,MH1,1)


          RMS=RMS+(XV(MH3,MH2,MH1,2)-EVAL(I,1))**2


9000        CONTINUE
          END DO
8000      CONTINUE
        END DO
7000    CONTINUE
      END DO


      RMS=WAVENM*DSQRT(RMS/MFIT)
      WRITE(IOUT,*)'ROOT MEAN SQUARE (CM-1) = ',RMS


      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETH3(XQ1,XQ2,XQ3,MM1,MM2,MM3,ENERGY,WK,
     1MODE1,MODE2,MODE3,MSYM1,MSYM2,MSYM3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),ENERGY(2+MM3,2+MM2,2+MM1,4)
      DIMENSION WK(2,1)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/INCREM/INCR
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/FILASS/IOUT,INP
103   FORMAT(5F15.8)
106   FORMAT(3I5)
C******************************
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS
      MMM1=MM1+2
C**EVEN NUMBER OF POINTS
      MID1=MMM1/2
      MMM2=MM2+2
      MID2=MMM2/2
      MMM3=MM3+2
      MID3=MMM3/2
C************************************************INCREMENT
C**NUMBER HEG EACH SIDE
      M1HALF=MM1/2
C**RESET MID1
      MID1=M1HALF/INCR + INCR
C**RESET MMM1
      MMM1=MID1*2
      M2HALF=MM2/2
      MID2=M2HALF/INCR + INCR
      MMM2=MID2*2
      M3HALF=MM3/2
      MID3=M3HALF/INCR + INCR
      MMM3=MID3*2
C************************************************INCREMENT
      MH1=MMM1-2
      MH2=MMM2-2
      MH3=MMM3-2
C**DO NOTHING IF ALL THREE MODES HAVE MSYM = 1
      IF(MSYM1.EQ.1.AND.MSYM2.EQ.1.AND.MSYM3.EQ.1)GO TO 9999

C**ANY TWO MODES HAVE MSYM = 1
      IF(MSYM1.EQ.1.AND.MSYM2.EQ.1)THEN
        DO M1=1,MMM1
          DO M2=1,MMM2
            DO M3=1,MID3
              K3=MID3+M3
              L3=MID3+1-M3
              ENERGY(K3,M2,M1,1)=ENERGY(L3,M2,M1,1)
              ENERGY(K3,M2,M1,2)=ENERGY(L3,M2,M1,2)
              ENERGY(K3,M2,M1,3)=ENERGY(L3,M2,M1,3)
              ENERGY(K3,M2,M1,4)=-ENERGY(L3,M2,M1,4)
            END DO
          END DO
        END DO
        GO TO 9999
      END IF
      IF(MSYM1.EQ.1.AND.MSYM3.EQ.1)THEN
        DO M1=1,MMM1
          DO M2=1,MID2
            K2=MID2+M2
            L2=MID2+1-M2
            DO M3=1,MMM3
              ENERGY(M3,K2,M1,1)=ENERGY(M3,L2,M1,1)
              ENERGY(M3,K2,M1,2)=ENERGY(M3,L2,M1,2)
              ENERGY(M3,K2,M1,3)=-ENERGY(M3,L2,M1,3)
              ENERGY(M3,K2,M1,4)=ENERGY(M3,L2,M1,4)
            END DO
          END DO
        END DO
        GO TO 9999
      END IF
      IF(MSYM2.EQ.1.AND.MSYM3.EQ.1)THEN
        DO M1=1,MID1
          K1=MID1+M1
          L1=MID1+1-M1
          DO M2=1,MMM2
            DO M3=1,MMM3
              ENERGY(M3,M2,K1,1)=ENERGY(M3,M2,L1,1)
              ENERGY(M3,M2,K1,2)=-ENERGY(M3,M2,L1,2)
              ENERGY(M3,M2,K1,3)=ENERGY(M3,M2,L1,3)
              ENERGY(M3,M2,K1,4)=ENERGY(M3,M2,L1,4)
            END DO
          END DO
        END DO
        GO TO 9999
      END IF

C**JUST ONE MODE HAS MSYM = 1
      IF(MSYM1.EQ.1)THEN
        DO M1=1,MMM1
          IF(MSYM2.EQ.MSYM3)THEN
C**REMAINING TWO MODES HAVE SAME SYMMETRY
            DO M2=1,MMM2
              K2=MMM2+1-M2
              DO M3=1,MID3
                K3=MMM3+1-M3
                ENERGY(K3,K2,M1,1)=ENERGY(M3,M2,M1,1)
                ENERGY(K3,K2,M1,2)=ENERGY(M3,M2,M1,2)
                ENERGY(K3,K2,M1,3)=-ENERGY(M3,M2,M1,3)
                ENERGY(K3,K2,M1,4)=-ENERGY(M3,M2,M1,4)
              END DO
            END DO
          ELSE
C**REMAINING TWO MODES HAVE DIFFERENT SYMMETRIES
            DO M2=1,MID2
              K2=MID2+M2
              L2=MID2+1-M2
              DO M3=1,MID3
                ENERGY(M3,K2,M1,1)=ENERGY(M3,L2,M1,1)
                ENERGY(M3,K2,M1,2)=ENERGY(M3,L2,M1,2)
                ENERGY(M3,K2,M1,3)=-ENERGY(M3,L2,M1,3)
                ENERGY(M3,K2,M1,4)=ENERGY(M3,L2,M1,4)
                K3=MID3+M3
                L3=MID3+1-M3
                ENERGY(K3,M2,M1,1)=ENERGY(L3,M2,M1,1)
                ENERGY(K3,M2,M1,2)=ENERGY(L3,M2,M1,2)
                ENERGY(K3,M2,M1,3)=ENERGY(L3,M2,M1,3)
                ENERGY(K3,M2,M1,4)=-ENERGY(L3,M2,M1,4)
                ENERGY(K3,K2,M1,1)=ENERGY(L3,L2,M1,1)
                ENERGY(K3,K2,M1,2)=ENERGY(L3,L2,M1,2)
                ENERGY(K3,K2,M1,3)=-ENERGY(L3,L2,M1,3)
                ENERGY(K3,K2,M1,4)=-ENERGY(L3,L2,M1,4)
              END DO
            END DO
          END IF
        END DO
        GO TO 9999
      END IF
      IF(MSYM2.EQ.1)THEN
        DO M2=1,MMM2
          IF(MSYM1.EQ.MSYM3)THEN
C**REMAINING TWO MODES HAVE SAME SYMMETRY
            DO M1=1,MMM1
              K1=MMM1+1-M1
              DO M3=1,MID3
                K3=MMM3+1-M3
                ENERGY(K3,M2,K1,1)=ENERGY(M3,M2,M1,1)
                ENERGY(K3,M2,K1,2)=-ENERGY(M3,M2,M1,2)
                ENERGY(K3,M2,K1,3)=ENERGY(M3,M2,M1,3)
                ENERGY(K3,M2,K1,4)=-ENERGY(M3,M2,M1,4)
              END DO
            END DO
          ELSE
C**REMAINING TWO MODES HAVE DIFFERENT SYMMETRIES
            DO M1=1,MID1
              K1=MID1+M1
              L1=MID1+1-M1
              DO M3=1,MID3
                ENERGY(M3,M2,K1,1)=ENERGY(M3,M2,L1,1)
                ENERGY(M3,M2,K1,2)=-ENERGY(M3,M2,L1,2)
                ENERGY(M3,M2,K1,3)=ENERGY(M3,M2,L1,3)
                ENERGY(M3,M2,K1,4)=ENERGY(M3,M2,L1,4)
                K3=MID3+M3
                L3=MID3+1-M3
                ENERGY(K3,M2,M1,1)=ENERGY(L3,M2,M1,1)
                ENERGY(K3,M2,M1,2)=ENERGY(L3,M2,M1,2)
                ENERGY(K3,M2,M1,3)=ENERGY(L3,M2,M1,3)
                ENERGY(K3,M2,M1,4)=-ENERGY(L3,M2,M1,4)
                ENERGY(K3,M2,K1,1)=ENERGY(L3,M2,L1,1)
                ENERGY(K3,M2,K1,2)=-ENERGY(L3,M2,L1,2)
                ENERGY(K3,M2,K1,3)=ENERGY(L3,M2,L1,3)
                ENERGY(K3,M2,K1,4)=-ENERGY(L3,M2,L1,4)
              END DO
            END DO
          END IF
        END DO
        GO TO 9999
      END IF
      IF(MSYM3.EQ.1)THEN
        DO M3=1,MMM3
          IF(MSYM1.EQ.MSYM2)THEN
C**REMAINING TWO MODES HAVE SAME SYMMETRY
            DO M1=1,MMM1
              K1=MMM1+1-M1
              DO M2=1,MID2
                K2=MMM2+1-M2
                ENERGY(M3,K2,K1,1)=ENERGY(M3,M2,M1,1)
                ENERGY(M3,K2,K1,2)=-ENERGY(M3,M2,M1,2)
                ENERGY(M3,K2,K1,3)=-ENERGY(M3,M2,M1,3)
                ENERGY(M3,K2,K1,4)=ENERGY(M3,M2,M1,4)
              END DO
            END DO
          ELSE
C**REMAINING TWO MODES HAVE DIFFERENT SYMMETRIES
            DO M1=1,MID1
              K1=MID1+M1
              L1=MID1+1-M1
              DO M2=1,MID2
                ENERGY(M3,M2,K1,1)=ENERGY(M3,M2,L1,1)
                ENERGY(M3,M2,K1,2)=-ENERGY(M3,M2,L1,2)
                ENERGY(M3,M2,K1,3)=ENERGY(M3,M2,L1,3)
                ENERGY(M3,M2,K1,4)=ENERGY(M3,M2,L1,4)
                K2=MID2+M2
                L2=MID2+1-M2
                ENERGY(M3,K2,M1,1)=ENERGY(M3,L2,M1,1)
                ENERGY(M3,K2,M1,2)=ENERGY(M3,L2,M1,2)
                ENERGY(M3,K2,M1,3)=-ENERGY(M3,L2,M1,3)
                ENERGY(M3,K2,M1,4)=ENERGY(M3,L2,M1,4)
                ENERGY(M3,K2,K1,1)=ENERGY(M3,L2,L1,1)
                ENERGY(M3,K2,K1,2)=-ENERGY(M3,L2,L1,2)
                ENERGY(M3,K2,K1,3)=-ENERGY(M3,L2,L1,3)
                ENERGY(M3,K2,K1,4)=ENERGY(M3,L2,L1,4)
              END DO
            END DO
          END IF
        END DO
        GO TO 9999
      END IF

C**********************************************MISSING
C**SPECIAL CASE MSYM1=MSYM2=MSYM3
      IF(MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)THEN
        GO TO 9999
      END IF
C**********************************************MISSING

C**TWO MODES EQUAL
      IF(MSYM1.EQ.MSYM2)THEN
        DO M1=1,MMM1
          DO M2=1,MID2
            K2=MID2+M2
            L2=MID2+1-M2
            DO M3=1,MID3
              ENERGY(M3,K2,M1,1)=ENERGY(M3,L2,M1,1)
              ENERGY(M3,K2,M1,2)=ENERGY(M3,L2,M1,2)
              ENERGY(M3,K2,M1,3)=-ENERGY(M3,L2,M1,3)
              ENERGY(M3,K2,M1,4)=ENERGY(M3,L2,M1,4)
              K3=MID3+M3
              L3=MID3+1-M3
              ENERGY(K3,M2,M1,1)=ENERGY(L3,M2,M1,1)
              ENERGY(K3,M2,M1,2)=ENERGY(L3,M2,M1,2)
              ENERGY(K3,M2,M1,3)=ENERGY(L3,M2,M1,3)
              ENERGY(K3,M2,M1,4)=-ENERGY(L3,M2,M1,4)
              ENERGY(K3,K2,M1,1)=ENERGY(L3,L2,M1,1)
              ENERGY(K3,K2,M1,2)=ENERGY(L3,L2,M1,2)
              ENERGY(K3,K2,M1,3)=-ENERGY(L3,L2,M1,3)
              ENERGY(K3,K2,M1,4)=-ENERGY(L3,L2,M1,4)
            END DO
          END DO
        END DO
        DO M3=1,MMM3
          DO M1=1,MMM1
            K1=MMM1+1-M1
            DO M2=1,MID2
              K2=MMM2+1-M2
              ENERGY(M3,K2,K1,1)=ENERGY(M3,M2,M1,1)
              ENERGY(M3,K2,K1,2)=-ENERGY(M3,M2,M1,2)
              ENERGY(M3,K2,K1,3)=-ENERGY(M3,M2,M1,3)
              ENERGY(M3,K2,K1,4)=ENERGY(M3,M2,M1,4)
            END DO
          END DO
        END DO
        GO TO 9999
      END IF
      IF(MSYM1.EQ.MSYM3)THEN
        DO M1=1,MMM1
          DO M2=1,MID2
            K2=MID2+M2
            L2=MID2+1-M2
            DO M3=1,MID3
              ENERGY(M3,K2,M1,1)=ENERGY(M3,L2,M1,1)
              ENERGY(M3,K2,M1,2)=ENERGY(M3,L2,M1,2)
              ENERGY(M3,K2,M1,3)=-ENERGY(M3,L2,M1,3)
              ENERGY(M3,K2,M1,4)=ENERGY(M3,L2,M1,4)
              K3=MID3+M3
              L3=MID3+1-M3
              ENERGY(K3,M2,M1,1)=ENERGY(L3,M2,M1,1)
              ENERGY(K3,M2,M1,2)=ENERGY(L3,M2,M1,2)
              ENERGY(K3,M2,M1,3)=ENERGY(L3,M2,M1,3)
              ENERGY(K3,M2,M1,4)=-ENERGY(L3,M2,M1,4)
              ENERGY(K3,K2,M1,1)=ENERGY(L3,L2,M1,1)
              ENERGY(K3,K2,M1,2)=ENERGY(L3,L2,M1,2)
              ENERGY(K3,K2,M1,3)=-ENERGY(L3,L2,M1,3)
              ENERGY(K3,K2,M1,4)=-ENERGY(L3,L2,M1,4)
            END DO
          END DO
        END DO
        DO M2=1,MMM2
          DO M1=1,MMM1
            K1=MMM1+1-M1
            DO M3=1,MID3
              K3=MMM3+1-M3
              ENERGY(K3,M2,K1,1)=ENERGY(M3,M2,M1,1)
              ENERGY(K3,M2,K1,2)=-ENERGY(M3,M2,M1,2)
              ENERGY(K3,M2,K1,3)=ENERGY(M3,M2,M1,3)
              ENERGY(K3,M2,K1,4)=-ENERGY(M3,M2,M1,4)
            END DO
          END DO
        END DO
        GO TO 9999
      END IF
      IF(MSYM2.EQ.MSYM3)THEN
        DO M2=1,MMM2
          DO M1=1,MID1
            K1=MID1+M1
            L1=MID1+1-M1
            DO M3=1,MID3
              ENERGY(M3,M2,K1,1)=ENERGY(M3,M2,L1,1)
              ENERGY(M3,M2,K1,2)=-ENERGY(M3,M2,L1,2)
              ENERGY(M3,M2,K1,3)=ENERGY(M3,M2,L1,3)
              ENERGY(M3,M2,K1,4)=ENERGY(M3,M2,L1,4)
              K3=MID3+M3
              L3=MID3+1-M3
              ENERGY(K3,M2,M1,1)=ENERGY(L3,M2,M1,1)
              ENERGY(K3,M2,M1,2)=ENERGY(L3,M2,M1,2)
              ENERGY(K3,M2,M1,3)=ENERGY(L3,M2,M1,3)
              ENERGY(K3,M2,M1,4)=-ENERGY(L3,M2,M1,4)
              ENERGY(K3,M2,K1,1)=ENERGY(L3,M2,L1,1)
              ENERGY(K3,M2,K1,2)=-ENERGY(L3,M2,L1,2)
              ENERGY(K3,M2,K1,3)=ENERGY(L3,M2,L1,3)
              ENERGY(K3,M2,K1,4)=-ENERGY(L3,M2,L1,4)
            END DO
          END DO
        END DO
        DO M1=1,MMM1
          DO M2=1,MMM2
            K2=MMM2+1-M2
            DO M3=1,MID3
              K3=MMM3+1-M3
              ENERGY(K3,K2,M1,1)=ENERGY(M3,M2,M1,1)
              ENERGY(K3,K2,M1,2)=ENERGY(M3,M2,M1,2)
              ENERGY(K3,K2,M1,3)=-ENERGY(M3,M2,M1,3)
              ENERGY(K3,K2,M1,4)=-ENERGY(M3,M2,M1,4)
            END DO
          END DO
        END DO
        GO TO 9999
      END IF
C**REMAINING MUST ALL BE DIFFERENT

C**********************************************MISSING
C**SPECIAL CASE MSYM1*MSYM2 = MSYM3
      N12=ISYMP(MSYM1,MSYM2)
      IF(MSYM3.NE.1.AND.(N12.EQ.MSYM3))THEN
      END IF
C**********************************************MISSING

9999  CONTINUE
      WRITE(MOUT,106)MMM1,MMM2,MMM3
      WRITE(MOUT,103)WK(1,MODE1),(XQ1(M),M=1,MH1),WK(2,MODE1)
      WRITE(MOUT,103)WK(1,MODE2),(XQ2(M),M=1,MH2),WK(2,MODE2)
      WRITE(MOUT,103)WK(1,MODE3),(XQ3(M),M=1,MH3),WK(2,MODE3)
C**ENERGIES AT M3 FOR EACH M2 AND M1
      DO M1=1,MMM1
        DO M2=1,MMM2
          WRITE(MOUT,103)(ENERGY(M3,M2,M1,1),M3=1,MMM3)
        END DO
      END DO
C**DERIVATIVES OF M1 FOR EACH M3 AND M2 (MODE1 > MODE2 > MODE3)
      DO M2=1,MMM2
        DO M3=1,MMM3
          WRITE(MOUT,103)(ENERGY(M3,M2,M1,2),M1=1,MMM1)
        END DO
      END DO
C**DERIVATIVES OF M2 FOR EACH M3 AND M1 (MODE1 > MODE2 > MODE3)
      DO M1=1,MMM1
        DO M3=1,MMM3
          WRITE(MOUT,103)(ENERGY(M3,M2,M1,3),M2=1,MMM2)
        END DO
      END DO
C**DERIVATIVES OF M3 FOR EACH M2 AND M1 (MODE1 > MODE2 > MODE3)
      DO M1=1,MMM1
        DO M2=1,MMM2
          WRITE(MOUT,103)(ENERGY(M3,M2,M1,4),M3=1,MMM3)
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE MOLDP4(XQ1,XQ2,XQ3,XQ4,MM1,MM2,MM3,MM4,NMODE,NATOM,QQ,
     1XX,X0,XL,XM,MODE1,MODE2,MODE3,MODE4,ENERGY,MSYM1,MSYM2,MSYM3,
     2MSYM4,XTANPM,WK,XJ1,XJ2,XJ3,XJ4,NP1,CP1,IP1,VP1,DP1,NTOT1,MMX1,
     3NP2,CP2,IP2,VP2,DP2A,DP2B,NTOT2,MAX2,NP3,CP3,IP3,VP3,DP3A,DP3B,
     4DP3C,NTOT3,MAX3,INDK,INDL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4),SYMBAD
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),QQ(NMODE)
      DIMENSION XJ1(MM1),XJ2(MM2),XJ3(MM3),XJ4(MM4)
      DIMENSION ENERGY(2+MM4,2+MM3,2+MM2,2+MM1,5)
      DIMENSION XTANPM(NMODE),WK(2,NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION VP1(MMX1,NTOT1),DP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION VP2(MAX2,MAX2,NTOT2),DP2A(MAX2,MAX2,NTOT2),
     1DP2B(MAX2,MAX2,NTOT2)
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),IP3(MAX3,NTOT3,3)
      DIMENSION VP3(MAX3*MAX3*MAX3,NTOT3),DP3A(MAX3*MAX3*MAX3,NTOT3),
     1DP3B(MAX3*MAX3*MAX3,NTOT3),DP3C(MAX3*MAX3*MAX3,NTOT3)
      DIMENSION INDK(1),INDL(1),DQ(NMODE)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/INCREM/INCR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FITTER/MFIT
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/FILASS/IOUT,INP
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
      MMM1=MM1
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
      IF(MOLPRO.LT.0)MMM1=MM1+2
C**EVEN NUMBER OF POINTS
      MID1=MMM1/2
      MMM2=MM2
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
      IF(MOLPRO.LT.0)MMM2=MM2+2
      MID2=MMM2/2
      MMM3=MM3
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
      IF(MOLPRO.LT.0)MMM3=MM3+2
      MID3=MMM3/2
      MMM4=MM4
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
      IF(MOLPRO.LT.0)MMM4=MM4+2
      MID4=MMM4/2

      MH1=0
      MJ1=0
C************************************************INCREMENT
      INCR=MOLINC
      INCR1=0
C************************************************INCREMENT
      INC1=1
      IF(MOLPRO.GT.0.AND.MOD(MID1,2).EQ.0)INC1=MOLINC
      DO M1=1,MMM1
        INCR1=INCR1+1
        IF(INCR1.EQ.INC1)THEN
C**NEED THIS POINT - RESET COUNT
          INCR1=0
        ELSE
C**SKIP THIS POINT
          GO TO 1000
        END IF
C**SET DEFAULT INCREMENT
        INC1=INCR
        IF(MOLPRO.LT.0)THEN
C**ONE AFTER FIRST GAUSS NEXT
          IF(M1.EQ.1)INC1=1
C**SECOND GAUSS NEXT
          IF(M1+1.EQ.MMM1)INC1=1
C**ONLY IF MID1 ODD
          IF(MOD(MID1,2).NE.0)THEN
            IF(M1.EQ.2)INC1=1
            IF(M1+2.EQ.MMM1)INC1=1
          END IF
        END IF
C**ONE AFTER MIDDLE ONE NEXT
        IF(M1.EQ.MID1)INC1=1
        IF(MOLPRO.LT.0)THEN
          IF(M1.GT.1.AND.M1.LT.MMM1)THEN
            MJ1=MJ1+1
            XJ1(MJ1)=XQ1(M1-1)
          END IF
          MX1=M1-1
        ELSE
          MX1=M1
        END IF
        IF(ISKIP1.NE.0.AND.M1.GT.MID1)GO TO 1000
        MH1=MH1+1
        QQ(MODE1)=XQ1(MX1)
        IF(MOLPRO.LT.0)THEN
          IF(M1.EQ.1)QQ(MODE1)=WK(1,MODE1)
          IF(M1.EQ.MMM1)QQ(MODE1)=WK(2,MODE1)
        END IF

        MH2=0
        MJ2=0
        INCR2=0
C************************************************INCREMENT
        INC2=1
        IF(MOLPRO.GT.0.AND.MOD(MID2,2).EQ.0)INC2=MOLINC
        DO M2=1,MMM2
          INCR2=INCR2+1
          IF(INCR2.EQ.INC2)THEN
            INCR2=0
          ELSE
            GO TO 2000
          END IF
          INC2=INCR
          IF(MOLPRO.LT.0)THEN
            IF(M2.EQ.1)INC2=1
            IF(M2+1.EQ.MMM2)INC2=1
            IF(MOD(MID2,2).NE.0)THEN
              IF(M2.EQ.2)INC2=1
              IF(M2+2.EQ.MMM2)INC2=1
            END IF
          END IF
          IF(M2.EQ.MID2)INC2=1
          IF(MOLPRO.LT.0)THEN
            IF(M2.GT.1.AND.M2.LT.MMM2)THEN
              MJ2=MJ2+1
              XJ2(MJ2)=XQ2(M2-1)
            END IF
            MX2=M2-1
          ELSE
            MX2=M2
          END IF
          IF(ISKIP2.NE.0.AND.M2.GT.MID2)GO TO 2000
          MH2=MH2+1
          QQ(MODE2)=XQ2(MX2)
          IF(MOLPRO.LT.0)THEN
            IF(M2.EQ.1)QQ(MODE2)=WK(1,MODE2)
            IF(M2.EQ.MMM2)QQ(MODE2)=WK(2,MODE2)
          END IF

          MH3=0
          MJ3=0
          INCR3=0
C************************************************INCREMENT
          INC3=1
          IF(MOLPRO.GT.0.AND.MOD(MID3,2).EQ.0)INC3=MOLINC
          DO M3=1,MMM3
            INCR3=INCR3+1
            IF(INCR3.EQ.INC3)THEN
              INCR3=0
            ELSE
              GO TO 3000
            END IF
            INC3=INCR
            IF(MOLPRO.LT.0)THEN
              IF(M3.EQ.1)INC3=1
              IF(M3+1.EQ.MMM3)INC3=1
              IF(MOD(MID3,2).NE.0)THEN
                IF(M3.EQ.2)INC3=1
                IF(M3+2.EQ.MMM3)INC3=1
              END IF
            END IF
            IF(M3.EQ.MID3)INC3=1
            IF(MOLPRO.LT.0)THEN
              IF(M3.GT.1.AND.M3.LT.MMM3)THEN
                MJ3=MJ3+1
                XJ3(MJ3)=XQ3(M3-1)
              END IF
              MX3=M3-1
            ELSE
              MX3=M3
            END IF
            IF(ISKIP3.NE.0.AND.M3.GT.MID3)GO TO 3000
            MH3=MH3+1
            QQ(MODE3)=XQ3(MX3)
            IF(MOLPRO.LT.0)THEN
              IF(M3.EQ.1)QQ(MODE3)=WK(1,MODE3)
              IF(M3.EQ.MMM3)QQ(MODE3)=WK(2,MODE3)
            END IF

            MH4=0
            MJ4=0
            INCR4=0
C************************************************INCREMENT
            INC4=1
            IF(MOLPRO.GT.0.AND.MOD(MID4,2).EQ.0)INC4=MOLINC
            DO M4=1,MMM4
              INCR4=INCR4+1
              IF(INCR4.EQ.INC4)THEN
                INCR4=0
              ELSE
                GO TO 4000
              END IF
              INC4=INCR
              IF(MOLPRO.LT.0)THEN
                IF(M4.EQ.1)INC4=1
                IF(M4+1.EQ.MMM4)INC4=1
                IF(MOD(MID4,2).NE.0)THEN
                  IF(M4.EQ.2)INC4=1
                  IF(M4+2.EQ.MMM4)INC4=1
                END IF
              END IF
              IF(M4.EQ.MID4)INC4=1
              IF(MOLPRO.LT.0)THEN
                IF(M4.GT.1.AND.M4.LT.MMM4)THEN
                  MJ4=MJ4+1
                  XJ4(MJ4)=XQ4(M4-1)
                END IF
                MX4=M4-1
              ELSE
                MX4=M4
              END IF
              IF(ISKIP4.NE.0.AND.M4.GT.MID4)GO TO 4000
              MH4=MH4+1
              QQ(MODE4)=XQ4(MX4)
              IF(MOLPRO.LT.0)THEN
                IF(M4.EQ.1)QQ(MODE4)=WK(1,MODE4)
                IF(M4.EQ.MMM4)QQ(MODE4)=WK(2,MODE4)
              END IF

              IF(IABS(MOLPRO).EQ.7)THEN
                IF(IWHICH.GT.0)THEN
C**FITTING GLOBAL POTENTIAL
C
C**FIRST GET POTENTIAL AT M1,M2,M3,M4
                  DO I=1,NATOM
                    DO K=1,3
                      XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1                SQRT(XM(I))
                    END DO
                  END DO
                  CALL GETPOT(VVV,NATOM,XX,ENERGY)
C**SECOND GET DERIVATIVE WRT MODE1
                  IF(M1.EQ.1.OR.M1.LT.MMM1)THEN
                    QQQ=QQ(MODE1)+1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVP1,NATOM,XX,ENERGY)
                  END IF
                  IF(M1.EQ.1)DD1=(VVP1-VVV)*1.D3
                  IF(M1.EQ.MMM1.OR.M1.GT.1)THEN
                    QQQ=QQ(MODE1)-1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVM1,NATOM,XX,ENERGY)
                  END IF
                  IF(M1.EQ.MMM1)DD1=(VVV-VVM1)*1.D3
                  IF(M1.GT.1.AND.M1.LT.MMM1)DD1=(VVP1-VVM1)*1.D3/2
C**THIRD GET DERIVATIVE WRT MODE2
                  IF(M2.EQ.1.OR.M2.LT.MMM2)THEN
                    QQQ=QQ(MODE2)+1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE2,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVP2,NATOM,XX,ENERGY)
                  END IF
                  IF(M2.EQ.1)DD2=(VVP2-VVV)*1.D3
                  IF(M2.EQ.MMM2.OR.M2.GT.1)THEN
                    QQQ=QQ(MODE2)-1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE2,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVM2,NATOM,XX,ENERGY)
                  END IF
                  IF(M2.EQ.MMM2)DD2=(VVV-VVM2)*1.D3
                  IF(M2.GT.1.AND.M2.LT.MMM2)DD2=(VVP2-VVM2)*1.D3/2
C**FOURTH GET DERIVATIVE WRT MODE3
                  IF(M3.EQ.1.OR.M3.LT.MMM3)THEN
                    QQQ=QQ(MODE3)+1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE3,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVP3,NATOM,XX,ENERGY)
                  END IF
                  IF(M3.EQ.1)DD3=(VVP3-VVV)*1.D3
                  IF(M3.EQ.MMM3.OR.M3.GT.1)THEN
                    QQQ=QQ(MODE3)-1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE3,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVM3,NATOM,XX,ENERGY)
                  END IF
                  IF(M3.EQ.MMM3)DD3=(VVV-VVM3)*1.D3
                  IF(M3.GT.1.AND.M3.LT.MMM3)DD3=(VVP3-VVM3)*1.D3/2
C**FIFTH GET DERIVATIVE WRT MODE4
                  IF(M4.EQ.1.OR.M4.LT.MMM4)THEN
                    QQQ=QQ(MODE4)+1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE4,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVP4,NATOM,XX,ENERGY)
                  END IF
                  IF(M4.EQ.1)DD4=(VVP4-VVV)*1.D3
                  IF(M4.EQ.MMM4.OR.M4.GT.1)THEN
                    QQQ=QQ(MODE4)-1.D-3
                    DO I=1,NATOM
                      DO K=1,3
                        XX(I,K)=X0(I,K)+XL(I,MODE4,K)*QQQ/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                  SQRT(XM(I))
                        XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                  SQRT(XM(I))
                      END DO
                    END DO
                    CALL GETPOT(VVM4,NATOM,XX,ENERGY)
                  END IF
                  IF(M4.EQ.MMM4)DD4=(VVV-VVM4)*1.D3
                  IF(M4.GT.1.AND.M4.LT.MMM4)DD4=(VVP4-VVM4)*1.D3/2
                END IF

                WRITE(MOUT,204)NATOM
                WRITE(MOUT,205)MODE1,CHSYM(MSYM1),M1,MODE2,
     1          CHSYM(MSYM2),M2,MODE3,CHSYM(MSYM3),M3,MODE4,
     2          CHSYM(MSYM4),M4

                IF(IWHICH.GT.0)THEN
C**FITTING GLOBAL POTENTIAL
                  IF(MOLPRO.GT.0)WRITE(MOUT,208)VVV
                  IF(MOLPRO.LT.0)WRITE(MOUT,208)VVV,DD1,DD2,DD3,DD4
                ELSE
C**FITTING ABINITIO POTENTIAL
                  WRITE(MOUT,*)'**********************************'
                END IF

              END IF 

              IF(MOLPRO.EQ.8)THEN
C**FIT ABINITIO POTENTIAL
                READ(MINP,204)NATOM
                READ(MINP,205)MODE1,CHSYM(MSYM1),MD1,MODE2,
     1          CHSYM(MSYM2),MD2,MODE3,CHSYM(MSYM3),MD3,MODE4,
     2          CHSYM(MSYM4),MD4
                READ(MINP,*)ENERGY(MH4,MH3,MH2,MH1,1)
              END IF
              IF(MOLPRO.EQ.-8)THEN
C**INTERPOLATE ABINITIO POTENTIAL
                READ(MINP,204)NATOM
                READ(MINP,205)MODE1,CHSYM(MSYM1),MD1,MODE2,
     1          CHSYM(MSYM2),MD2,MODE3,CHSYM(MSYM3),MD3,MODE4,
     2          CHSYM(MSYM4),MD4
                READ(MINP,*)ENERGY(MH4,MH3,MH2,MH1,1),
     1          ENERGY(MH4,MH3,MH2,MH1,2),ENERGY(MH4,MH3,MH2,MH1,3),
     2          ENERGY(MH4,MH3,MH2,MH1,4),ENERGY(MH4,MH3,MH2,MH1,5)
              END IF
              DO I=1,NATOM
                IF(IABS(MOLPRO).EQ.7)THEN
                  DO K=1,3
                    XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1              SQRT(XM(I))
                    XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1              SQRT(XM(I))
                    XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1              SQRT(XM(I))
                    XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1              SQRT(XM(I))
                  END DO
                  WRITE(MOUT,206)SYMBOL(I),(XX(I,K)*BOHR,K=1,3)
                END IF
                IF(IABS(MOLPRO).EQ.8)THEN
                  READ(MINP,206)SYMBAD,(XX(I,K),K=1,3)
                  DO K=1,3
                    XX(I,K)=XX(I,K)/BOHR
                  END DO
                END IF
              END DO
C**TAKE OUT 1-DIM, 2-DIM, 3-DIM POTENTIALS + CONSTANT
              IF(MOLPRO.EQ.8)THEN
                CALL GETPQ3(VQ3,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1          NP2,CP2,IP2,NTOT2,MAX2,NP3,CP3,IP3,NTOT3,MAX3,INDK,
     2          INDL)
                ENERGY(MH4,MH3,MH2,MH1,2)=ENERGY(MH4,MH3,MH2,MH1,1)-
     1          VQ3
              END IF
              IF(MOLPRO.EQ.-8)THEN
              END IF
4000          CONTINUE
            END DO
3000        CONTINUE
          END DO
2000      CONTINUE
        END DO
1000    CONTINUE
      END DO
      MFIT=MH1*MH2*MH3*MH4
      IF(MOLPRO.EQ.-8)CALL GETH4(XJ1,XJ2,XJ3,XJ4,MM1,MM2,MM3,MM4,
     1ENERGY,WK,MODE1,MODE2,MODE3,MODE4,MSYM1,MSYM2,MSYM3,MSYM4)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETQ4(XK,EVAL,SOL,MM1,XQ1,MM2,XQ2,MM3,XQ3,MM4,XQ4,XV,
     1MODE1,MODE2,MODE3,MODE4,MSYM1,MSYM2,MSYM3,MSYM4,NMODE,QQ,XTANPM,
     2WRK,MDIM,VINF,D,NP1,CP1,IP1,NTOT1,MMX1,NP2,CP2,IP2,NTOT2,MAX2,
     3NP3,CP3,IP3,NTOT3,MAX3,INDK,INDL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XK(MDIM,MDIM),EVAL(MDIM,6)
      DIMENSION SOL(MDIM,1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION XV(2+MM4,2+MM3,2+MM2,2+MM1,2)
      DIMENSION QQ(NMODE),VINF(NMODE)
      DIMENSION XTANPM(NMODE),WRK(1)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),IP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)


      DIMENSION D(MDIM,MDIM)
      COMMON/VMIN/VMIN,VMAX


      COMMON/SYMMP/ISYMP(10,10)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FITTER/MFIT,MN1,MMM1,ISYM1,ISKIP1,XTAN1,MN2,MMM2,ISYM2,
     1ISKIP2,XTAN2,MN3,MMM3,ISYM3,ISKIP3,XTAN3,MN4,MMM4,ISYM4,ISKIP4,
     2XTAN4
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/RMS/RMS
      COMMON/FILASS/IOUT,INP
100   FORMAT(/,1X,'DGESV IFAIL = ',I3)
101   FORMAT(D20.10,5X,I3,2X,I3,2X,I3,2X,I3)
102   FORMAT(//,1X,'ENTERING EXACT FIT ROUTINE',/)
103   FORMAT(//,1X,'ENTERING LEAST SQUARES FIT ROUTINE',/,
     11X,'MFIT = ',I4,' NFIT = ',I4,/)
105   FORMAT(/,7X,'COEFFICIENT',9X,'POWERS',/)
108   FORMAT(/,1X,'FITTING OF POTENTIAL FOR MODES: ',I3,2X,I3,2X,I3,2X,
     1I3)
109   FORMAT(//,7X,'OBSERVED',7X,'CALCULATED',/)
110   FORMAT(1X,F12.4,5X,F12.4)
111   FORMAT(1X,'SYMMETRIES ',A2,2X,A2,2X,A2,2X,A2)
112   FORMAT(//,7X,'OBSERVED',7X,'CALCULATED',9X,'TOTAL',/)
113   FORMAT(1X,F12.4,5X,F12.4,5X,F12.4)
C***********************************************

      DO K=1,NMODE
        QQ(K)=0
      END DO
      WRITE(IOUT,108)MODE1,MODE2,MODE3,MODE4
      WRITE(IOUT,111)CHSYM(MSYM1),CHSYM(MSYM2),CHSYM(MSYM3),
     1CHSYM(MSYM4)
      MMM1=MM1
      MMM2=MM2
      MMM3=MM3
      MMM4=MM4
C************************************************INCREMENT
      MX1=MM1/MOLINC
      MX2=MM2/MOLINC
      MX3=MM3/MOLINC
      MX4=MM4/MOLINC
      MMM1=MX1+MOD(MX1,MOLINC)
      MMM2=MX2+MOD(MX2,MOLINC)
      MMM3=MX3+MOD(MX3,MOLINC)
      MMM4=MX4+MOD(MX4,MOLINC)
      ISYM1=MSYM1
      ISYM2=MSYM2
      ISYM3=MSYM3
      ISYM4=MSYM4
      N12=ISYMP(MSYM1,MSYM2)
      N13=ISYMP(MSYM1,MSYM3)
      N14=ISYMP(MSYM1,MSYM4)
      N23=ISYMP(MSYM2,MSYM3)
      N24=ISYMP(MSYM2,MSYM4)
      N34=ISYMP(MSYM3,MSYM4)
      XTAN1=XTANPM(MODE1)
      XTAN2=XTANPM(MODE2)
      XTAN3=XTANPM(MODE3)
      XTAN4=XTANPM(MODE4)

      NDIM=MM1*MM2*MM3*MM4
      N=MMM1*MMM2*MMM3*MMM4
C**SIZE OF POLYNOMIAL DEPENDS ON SYMMETRY
      NM=N
      IF(MSYM1.NE.1)NM=N/2
      IF(MSYM2.NE.1)NM=N/2
      IF(MSYM3.NE.1)NM=N/2
      IF(MSYM4.NE.1)NM=N/2
      IF(MSYM1.NE.1.AND.MSYM2.NE.1.AND.N12.NE.1)NM=N/4
      IF(MSYM1.NE.1.AND.MSYM3.NE.1.AND.N13.NE.1)NM=N/4
      IF(MSYM1.NE.1.AND.MSYM4.NE.1.AND.N14.NE.1)NM=N/4
      IF(MSYM2.NE.1.AND.MSYM3.NE.1.AND.N23.NE.1)NM=N/4
      IF(MSYM2.NE.1.AND.MSYM4.NE.1.AND.N24.NE.1)NM=N/4
      IF(MSYM3.NE.1.AND.MSYM4.NE.1.AND.N34.NE.1)NM=N/4

C**POWERS OF POLYNOMIAL DEPENDS ON SYMMETRY
      MN1=1
      MN2=1
      MN3=1
      MN4=1
      IF(MSYM1.NE.1)MN1=2
      IF(MSYM2.NE.1)MN2=2
      IF(MSYM3.NE.1)MN3=2
      IF(MSYM4.NE.1)MN4=2
C**SPECIAL CASE B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2, ETC.
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM2))MN1=1
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM3))MN1=1
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM4))MN1=1
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM3))MN2=1
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM4))MN2=1
      IF(MSYM3.NE.1.AND.(MSYM3.EQ.MSYM4))MN3=1

C**SPECIAL CASE S2*S3 = S1
      N23=ISYMP(MSYM2,MSYM3)
      IF(N23.EQ.MSYM1)MN1=1
C**SPECIAL CASE S2*S4 = S1
      N24=ISYMP(MSYM2,MSYM4)
      IF(N24.EQ.MSYM1)MN1=1
C**SPECIAL CASE S3*S4 = S1
      N34=ISYMP(MSYM3,MSYM4)
      IF(N34.EQ.MSYM1)MN1=1
C**SPECIAL CASE S3*S4 = S2
      IF(N34.EQ.MSYM2)MN2=1

C**SPECIAL CASE S2*S3*S4 = S1
      N234=ISYMP(MSYM2,N34)
      IF(N234.EQ.MSYM1)MN1=1

C**ONLY NEED POTENTIAL DATA IF UNIQUE SYMMETRY CONFIGURATION
      ISKIP1=0
      IF(MSYM1.NE.1)ISKIP1=1
      ISKIP2=0
      IF(MSYM2.NE.1)ISKIP2=1
      ISKIP3=0
      IF(MSYM3.NE.1)ISKIP3=1
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
 
      IOFF4=0
      IOFF3=0
      IOFF2=0
      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 1000
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 2000
          MH3=MMM3/2
          MX3=MM3/2
          IS3=-1
          INH3=0
          INC3=0
          DO M3=1,MMM3
            MH3=MH3+INH3*IS3
            MX3=MX3+INC3*IS3
            IS3=-IS3
            INH3=INH3+1
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 3000
            MH4=MMM4/2
            MX4=MM4/2
            IS4=-1
            INH4=0
            INC4=0
            DO M4=1,MMM4
              MH4=MH4+INH4*IS4
              MX4=MX4+INC4*IS4
              IS4=-IS4
              INH4=INH4+1
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC4=INC4+1
              ELSE
                INC4=2*M4-1
              END IF
              IF(ISKIP4.NE.0.AND.MOD(M4,2).EQ.0)GO TO 4000
              I=I+1
              EVAL(I,1)=XV(MH4,MH3,MH2,MH1,1)
              EVAL(I,2)=XV(MH4,MH3,MH2,MH1,2)
C**WEIGHT (MINIMUM SCALED TO UNITY)
              VSCAL=WAVENM*(EVAL(I,1)-VMIN)+1
              EVAL(I,3)=1/VSCAL

              J=0
              DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1          IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1          IOFF4=MOD(J1,2)+1-MN4
                DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1            IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1            IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1            IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1            IOFF4=MOD(J1+J2,2)-MN4
                  DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1              IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM4)IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     2              IOFF4=MOD(J1+J2+J3,2)-MN4
                    DO J4=1,MMM4/MN4
                 J=J+1
                 XK(I,J)=(XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     1           (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     2           (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))*
     3           (XTANH(XTANPM(MODE4)*XQ4(MX4))**(MN4*(J4-1)-IOFF4))
                    END DO
                  END DO
                END DO
              END DO
4000          CONTINUE
            END DO
3000        CONTINUE
          END DO
2000      CONTINUE
        END DO
1000    CONTINUE
      END DO
C******************************
      WRITE(IOUT,102)
      IFAIL=1
      CALL DGESV(NM,1,XK,MDIM,SOL,EVAL,NM,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C******************************
      DO I=1,NM
        SOL(I,1)=EVAL(I,1)
      END DO
      WRITE(IOUT,105)
      J=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1  IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1  IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1  IOFF4=MOD(J1,2)+1-MN4
        DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1    IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1    IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1    IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1    IOFF4=MOD(J1+J2,2)-MN4
          DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
            IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1      IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     1      IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM4)
     1      IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1      MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J2+J3,2)-MN4
            DO J4=1,MMM4/MN4
              J=J+1
              WRITE(IOUT,101)SOL(J,1),MN1*(J1-1),MN2*(J2-1)-IOFF2,
     1        MN3*(J3-1)-IOFF3,MN4*(J4-1)-IOFF4
            END DO
          END DO
        END DO
      END DO
      CALL FLUSH(IOUT)

      WRITE(IOUT,109)
      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 5000
        QQ(MODE1)=XQ1(MX1)
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 6000
          QQ(MODE2)=XQ2(MX2)
          MH3=MMM3/2
          MX3=MM3/2
          IS3=-1
          INH3=0
          INC3=0
          DO M3=1,MMM3
            MH3=MH3+INH3*IS3
            MX3=MX3+INC3*IS3
            IS3=-IS3
            INH3=INH3+1
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 7000
            QQ(MODE3)=XQ3(MX3)
            MH4=MMM4/2
            MX4=MM4/2
            IS4=-1
            INH4=0
            INC4=0
            DO M4=1,MMM4
              MH4=MH4+INH4*IS4
              MX4=MX4+INC4*IS4
              IS4=-IS4
              INH4=INH4+1
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC4=INC4+1
              ELSE
                INC4=2*M4-1
              END IF
              IF(ISKIP4.NE.0.AND.MOD(M4,2).EQ.0)GO TO 8000
              QQ(MODE4)=XQ4(MX4)
              I=I+1
              EVAL(I,1)=0
              J=0
              DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1          IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1          IOFF4=MOD(J1,2)+1-MN4
                DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1            IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1            IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1            IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1            IOFF4=MOD(J1+J2,2)-MN4
                  DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1              IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM4)IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     2              IOFF4=MOD(J1+J2+J3,2)-MN4
                    DO J4=1,MMM4/MN4
                      J=J+1
                      EVAL(I,1)=EVAL(I,1)+SOL(J,1)*
     1            (XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     2            (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     3            (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))*
     4            (XTANH(XTANPM(MODE4)*XQ4(MX4))**(MN4*(J4-1)-IOFF4))
                    END DO
                  END DO
                END DO
              END DO
C**ADD IN 1-DIM, 2-DIM, 3-DIM POTENTIALS + CONSTANT
              WRITE(IOUT,110)WAVENM*XV(MH4,MH3,MH2,MH1,1),
     1        WAVENM*EVAL(I,1)
8000          CONTINUE
            END DO
7000        CONTINUE
          END DO
6000      CONTINUE
        END DO
5000    CONTINUE
      END DO
C******************************
C**LEAST SQUARES
C******************************
      NFIT=0
      J=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1  IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1  IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1  IOFF4=MOD(J1,2)+1-MN4
        DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1    IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1    IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1    IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1    IOFF4=MOD(J1+J2,2)-MN4
          DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
            IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1      IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     1      IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM4)
     1      IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1      MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J2+J3,2)-MN4
            DO J4=1,MMM4/MN4
              J=J+1
              IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1        MN3*(J3-1)-IOFF3.NE.0.AND.MN4*(J4-1)-IOFF4.NE.0)THEN
                NFIT=NFIT+1
C**COPY ONLY RELEVANT COEFFICIENTS FOR Vijkl(4)
                SOL(NFIT,1)=SOL(J,1)
              END IF
            END DO
          END DO
        END DO
      END DO
      WRITE(IOUT,103)MFIT,NFIT


C**A-MATRIX
      I=0
      MX1=MM1/2
      IS1=-1
      INC1=0
      DO M1=1,MMM1
        MX1=MX1+INC1*IS1
        IS1=-IS1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 1111
        MX2=MM2/2
        IS2=-1
        INC2=0
        DO M2=1,MMM2
          MX2=MX2+INC2*IS2
          IS2=-IS2
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 2222
          MX3=MM3/2
          IS3=-1
          INC3=0
          DO M3=1,MMM3
            MX3=MX3+INC3*IS3
            IS3=-IS3
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 3333
            MX4=MM4/2
            IS4=-1
            INC4=0
            DO M4=1,MMM4
              MX4=MX4+INC4*IS4
              IS4=-IS4
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC4=INC4+1
              ELSE
                INC4=2*M4-1
              END IF
              IF(ISKIP4.NE.0.AND.MOD(M4,2).EQ.0)GO TO 4444
              I=I+1
              J=0
              DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1          IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1          IOFF4=MOD(J1,2)+1-MN4
                DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1            IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1            IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1            IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1            IOFF4=MOD(J1+J2,2)-MN4
                  DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1              IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM4)IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     2              IOFF4=MOD(J1+J2+J3,2)-MN4
                    DO J4=1,MMM4/MN4
                      IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1                MN3*(J3-1)-IOFF3.NE.0.AND.MN4*(J4-1)-IOFF4.NE.0)
     2                THEN
                 J=J+1
                 XK(I,J)=(XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     1           (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     2           (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))*
     3           (XTANH(XTANPM(MODE4)*XQ4(MX4))**(MN4*(J4-1)-IOFF4))
                      END IF
                    END DO
                  END DO
                END DO
              END DO
4444          CONTINUE
            END DO
3333        CONTINUE
          END DO
2222      CONTINUE
        END DO
1111    CONTINUE
      END DO
C**D-MATRIX
      DO J=1,NFIT
        DO K=1,NFIT
          X=0
          DO I=1,MFIT
            X=X+EVAL(I,3)*XK(I,K)*XK(I,J)
          END DO
          D(K,J)=X
        END DO
      END DO
C**E-MATRIX
      DO J=1,NFIT
        X=0
        DO I=1,MFIT
          X=X+EVAL(I,3)*XK(I,J)*EVAL(I,2)
        END DO
        EVAL(J,1)=X
      END DO

C******************************
      IFAIL=1
      CALL DGESV(NFIT,1,D,MDIM,SOL,EVAL,NFIT,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C******************************
      DO I=1,NFIT
        SOL(I,1)=EVAL(I,1)
      END DO


      WRITE(IOUT,105)
      WRITE(MOUT,*)NFIT
      NFIT=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1  IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1  IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1  IOFF4=MOD(J1,2)+1-MN4
        DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1    IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1    IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1    IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1    IOFF4=MOD(J1+J2,2)-MN4
          DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
            IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1      IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     1      IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM4)
     1      IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1      MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J2+J3,2)-MN4
            DO J4=1,MMM4/MN4
              IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1        MN3*(J3-1)-IOFF3.NE.0.AND.MN4*(J4-1)-IOFF4.NE.0)THEN
                NFIT=NFIT+1
                WRITE(IOUT,101)SOL(NFIT,1),MN1*(J1-1),MN2*(J2-1)-IOFF2,
     1          MN3*(J3-1)-IOFF3,MN4*(J4-1)-IOFF4
                WRITE(MOUT,101)SOL(NFIT,1),MN1*(J1-1),MN2*(J2-1)-IOFF2,
     1          MN3*(J3-1)-IOFF3,MN4*(J4-1)-IOFF4
              END IF
            END DO
          END DO
        END DO
      END DO
      CALL FLUSH(IOUT)

      WRITE(IOUT,112)


      RMS=0


      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 9000
        QQ(MODE1)=XQ1(MX1)
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 10000
          QQ(MODE2)=XQ2(MX2)
          MH3=MMM3/2
          MX3=MM3/2
          IS3=-1
          INH3=0
          INC3=0
          DO M3=1,MMM3
            MH3=MH3+INH3*IS3
            MX3=MX3+INC3*IS3
            IS3=-IS3
            INH3=INH3+1
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 11000
            QQ(MODE3)=XQ3(MX3)
            MH4=MMM4/2
            MX4=MM4/2
            IS4=-1
            INH4=0
            INC4=0
            DO M4=1,MMM4
              MH4=MH4+INH4*IS4
              MX4=MX4+INC4*IS4
              IS4=-IS4
              INH4=INH4+1
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC4=INC4+1
              ELSE
                INC4=2*M4-1
              END IF
              IF(ISKIP4.NE.0.AND.MOD(M4,2).EQ.0)GO TO 12000
              QQ(MODE4)=XQ4(MX4)
              I=I+1
              EVAL(I,1)=0
              NFIT=0
              DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1          IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1          IOFF4=MOD(J1,2)+1-MN4
                DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1            IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1            IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1            IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1            IOFF4=MOD(J1+J2,2)-MN4
                  DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1              IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM4)IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     2              IOFF4=MOD(J1+J2+J3,2)-MN4
                    DO J4=1,MMM4/MN4
                      IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1                MN3*(J3-1)-IOFF3.NE.0.AND.MN4*(J4-1)-IOFF4.NE.0)
     2                THEN
                        NFIT=NFIT+1
                        EVAL(I,1)=EVAL(I,1)+SOL(NFIT,1)*
     1             (XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     2             (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     3             (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))*
     4             (XTANH(XTANPM(MODE4)*XQ4(MX4))**(MN4*(J4-1)-IOFF4))
                      END IF
                    END DO
                  END DO
                END DO
              END DO
C**ADD IN 1-DIM, 2-DIM, 3-DIM POTENTIALS + CONSTANT
              CALL GETPQ3(VQ3,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1        NP2,CP2,IP2,NTOT2,MAX2,NP3,CP3,IP3,NTOT3,MAX3,INDK,INDL)
              WRITE(IOUT,113)WAVENM*XV(MH4,MH3,MH2,MH1,2),
     1        WAVENM*EVAL(I,1),WAVENM*XV(MH4,MH3,MH2,MH1,1)


          RMS=RMS+(XV(MH4,MH3,MH2,MH1,2)-EVAL(I,1))**2


12000         CONTINUE
            END DO
11000       CONTINUE
          END DO
10000     CONTINUE
        END DO
9000    CONTINUE
      END DO


      RMS=WAVENM*DSQRT(RMS/MFIT)
      WRITE(IOUT,*)'ROOT MEAN SQUARE (CM-1) = ',RMS


      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETH4(XQ1,XQ2,XQ3,XQ4,MM1,MM2,MM3,MM4,ENERGY,WK,
     1MODE1,MODE2,MODE3,MODE4,MSYM1,MSYM2,MSYM3,MSYM4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION ENERGY(2+MM4,2+MM3,2+MM2,2+MM1,5)
      DIMENSION WK(2,1)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/INCREM/INCR
      COMMON/FILASS/IOUT,INP
103   FORMAT(5F15.8)
106   FORMAT(4I5)
C******************************
C******************************
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS
      MMM1=MM1+2
C**EVEN NUMBER OF POINTS
      MID1=MMM1/2
      MMM2=MM2+2
      MID2=MMM2/2
      MMM3=MM3+2
      MID3=MMM3/2
      MMM4=MM4+2
      MID4=MMM4/2
C************************************************INCREMENT
C**NUMBER HEG EACH SIDE
      M1HALF=MM1/2
C**RESET MID1
      MID1=M1HALF/INCR + INCR
C**RESET MMM1
      MMM1=MID1*2
      M2HALF=MM2/2
      MID2=M2HALF/INCR + INCR
      MMM2=MID2*2
      M3HALF=MM3/2
      MID3=M3HALF/INCR + INCR
      MMM3=MID3*2
      M4HALF=MM4/2
      MID4=M4HALF/INCR + INCR
      MMM4=MID4*2
C************************************************INCREMENT
      MH1=MMM1-2
      MH2=MMM2-2
      MH3=MMM3-2
      MH4=MMM4-2
      WRITE(MOUT,106)MMM1,MMM2,MMM3,MMM4
      WRITE(MOUT,103)WK(1,MODE1),(XQ1(M),M=1,MH1),WK(2,MODE1)
      WRITE(MOUT,103)WK(1,MODE2),(XQ2(M),M=1,MH2),WK(2,MODE2)
      WRITE(MOUT,103)WK(1,MODE3),(XQ3(M),M=1,MH3),WK(2,MODE3)
      WRITE(MOUT,103)WK(1,MODE4),(XQ4(M),M=1,MH4),WK(2,MODE4)
C**ENERGIES AT M4 FOR EACH M3, M2 AND M1
      DO M1=1,MMM1
        DO M2=1,MMM2
          DO M3=1,MMM3
            WRITE(MOUT,103)(ENERGY(M4,M3,M2,M1,1),M4=1,MMM4)
          END DO
        END DO
      END DO
C**DERIVATIVES OF M1 FOR EACH M4, M3 AND M2 
      DO M2=1,MMM2
        DO M3=1,MMM3
          DO M4=1,MMM4
            WRITE(MOUT,103)(ENERGY(M4,M3,M2,M1,2),M1=1,MMM1)
          END DO
        END DO
      END DO
C**DERIVATIVES OF M2 FOR EACH M4, M3 AND M1
      DO M1=1,MMM1
        DO M3=1,MMM3
          DO M4=1,MMM4
            WRITE(MOUT,103)(ENERGY(M4,M3,M2,M1,3),M2=1,MMM2)
          END DO
        END DO
      END DO
C**DERIVATIVES OF M3 FOR EACH M4, M2 AND M1
      DO M1=1,MMM1
        DO M2=1,MMM2
          DO M4=1,MMM4
            WRITE(MOUT,103)(ENERGY(M4,M3,M2,M1,4),M3=1,MMM3)
          END DO
        END DO
      END DO
C**DERIVATIVES OF M4 FOR EACH M3, M2 AND M1
      DO M1=1,MMM1
        DO M2=1,MMM2
          DO M3=1,MMM3
            WRITE(MOUT,103)(ENERGY(M4,M3,M2,M1,5),M4=1,MMM4)
          END DO
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE MOLDP5(XQ1,XQ2,XQ3,XQ4,XQ5,MM1,MM2,MM3,MM4,MM5,NMODE,
     1NATOM,QQ,XX,X0,XL,XM,MODE1,MODE2,MODE3,MODE4,MODE5,ENERGY,MSYM1,
     2MSYM2,MSYM3,MSYM4,MSYM5,XTANPM,WK,XJ1,XJ2,XJ3,XJ4,XJ5,NP1,CP1,
     3IP1,NTOT1,MMX1,NP2,CP2,IP2,NTOT2,MAX2,NP3,CP3,IP3,NTOT3,MAX3,
     4NP4,CP4,IP4,NTOT4,MAX4,INDK,INDL,INDN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4),SYMBAD
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQ5(MM5),QQ(NMODE)
      DIMENSION XJ1(MM1),XJ2(MM2),XJ3(MM3),XJ4(MM4),XJ5(MM5)
      DIMENSION ENERGY(MM5,MM4,MM3,MM2,MM1,2)
      DIMENSION XTANPM(NMODE),WK(2,NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),IP3(MAX3,NTOT3,3)
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),IP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/INCREM/INCR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FITTER/MFIT
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/FILASS/IOUT,INP
204   FORMAT(I4)
205   FORMAT(1X,'MODE:',I2,' SYM:',A2,' POINT:',I2,
     1       1X,'MODE:',I2,' SYM:',A2,' POINT:',I2,
     2       1X,'MODE:',I2,' SYM:',A2,' POINT:',I2,
     3       1X,'MODE:',I2,' SYM:',A2,' POINT:',I2,
     4       1X,'MODE:',I2,' SYM:',A2,' POINT:',I2 )
206   FORMAT(A2,1X,3F20.10)
208   FORMAT(6F15.8)

C**ONLY DUMP POTENTIAL DATA IF UNIQUE SYMMETRY CONFIGURATION
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE1)MSYM1=I
          IF(ISYM(I,K).EQ.MODE2)MSYM2=I
          IF(ISYM(I,K).EQ.MODE3)MSYM3=I
          IF(ISYM(I,K).EQ.MODE4)MSYM4=I
          IF(ISYM(I,K).EQ.MODE5)MSYM5=I
        END DO
      END DO
      ISKIP1=0
      IF(MSYM1.NE.1)ISKIP1=1
      ISKIP2=0
      IF(MSYM2.NE.1)ISKIP2=1
      ISKIP3=0
      IF(MSYM3.NE.1)ISKIP3=1
      ISKIP4=0
      IF(MSYM4.NE.1)ISKIP4=1
      ISKIP5=0
      IF(MSYM5.NE.1)ISKIP5=1
C**SPECIAL CASE B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2, ETC.
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM2))ISKIP1=0
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM3))ISKIP1=0
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM4))ISKIP1=0
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM5))ISKIP1=0
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM3))ISKIP2=0
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM4))ISKIP2=0
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM5))ISKIP2=0
      IF(MSYM3.NE.1.AND.(MSYM3.EQ.MSYM4))ISKIP3=0
      IF(MSYM3.NE.1.AND.(MSYM3.EQ.MSYM5))ISKIP3=0
      IF(MSYM4.NE.1.AND.(MSYM4.EQ.MSYM5))ISKIP4=0
      
      N12=ISYMP(MSYM1,MSYM2)
      N13=ISYMP(MSYM1,MSYM3)
      N14=ISYMP(MSYM1,MSYM4)
      N15=ISYMP(MSYM1,MSYM5)
      N23=ISYMP(MSYM2,MSYM3)
      N24=ISYMP(MSYM2,MSYM4)
      N25=ISYMP(MSYM2,MSYM5)
      N34=ISYMP(MSYM3,MSYM4)
      N35=ISYMP(MSYM3,MSYM5)
      N45=ISYMP(MSYM4,MSYM5)
      
C**SPECIAL CASE S2*S3 = S1
      IF(N23.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S2*S4 = S1
      IF(N24.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S2*S5 = S1
      IF(N25.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S4 = S1
      IF(N34.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S5 = S1
      IF(N35.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S4*S5 = S1
      IF(N45.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S4 = S2
      IF(N34.EQ.MSYM2)ISKIP2=0
C**SPECIAL CASE S3*S5 = S2
      IF(N35.EQ.MSYM2)ISKIP2=0
C**SPECIAL CASE S4*S5 = S2
      IF(N45.EQ.MSYM2)ISKIP2=0
C**SPECIAL CASE S4*S5 = S3
      IF(N45.EQ.MSYM3)ISKIP3=0

      N234=ISYMP(MSYM2,N34)
      N235=ISYMP(MSYM2,N35)
      N245=ISYMP(MSYM2,N45)
      N345=ISYMP(MSYM3,N45)
      
C**SPECIAL CASE S2*S3*S4 = S1
      IF(N234.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S2*S3*S5 = S1
      IF(N235.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S2*S4*S5 = S1
      IF(N245.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S4*S5 = S1
      IF(N345.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S4*S5 = S2
      IF(N345.EQ.MSYM2)ISKIP2=0

      N2345=ISYMP(N23,N45)
      
C**SPECIAL CASE S2*S3*S4*S5 = S1
      IF(N2345.EQ.MSYM1)ISKIP1=0      
     

      DO K=1,NMODE
        QQ(K)=0
      END DO
      MMM1=MM1
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
C     IF(MOLPRO.LT.0)MMM1=MM1+2
C**EVEN NUMBER OF POINTS
      MID1=MMM1/2
      MMM2=MM2
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
C     IF(MOLPRO.LT.0)MMM2=MM2+2
      MID2=MMM2/2
      MMM3=MM3
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
C     IF(MOLPRO.LT.0)MMM3=MM3+2
      MID3=MMM3/2
      MMM4=MM4
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
C     IF(MOLPRO.LT.0)MMM4=MM4+2
      MID4=MMM4/2
      MMM5=MM5
C**MAKE ALLOWANCE FOR EXTREME GAUSS POINTS IF INTERPOLATION
C     IF(MOLPRO.LT.0)MMM5=MM5+2
      MID5=MMM5/2

      MH1=0
      MJ1=0
C************************************************INCREMENT
      INCR=MOLINC
      INCR1=0
C************************************************INCREMENT
      INC1=1
      IF(MOLPRO.GT.0.AND.MOD(MID1,2).EQ.0)INC1=MOLINC
      DO M1=1,MMM1
        INCR1=INCR1+1
        IF(INCR1.EQ.INC1)THEN
C**NEED THIS POINT - RESET COUNT
          INCR1=0
        ELSE
C**SKIP THIS POINT
          GO TO 1000
        END IF
C**SET DEFAULT INCREMENT
        INC1=INCR
C**ONE AFTER MIDDLE ONE NEXT
        IF(M1.EQ.MID1)INC1=1
        MX1=M1
        IF(ISKIP1.NE.0.AND.M1.GT.MID1)GO TO 1000
        MH1=MH1+1
        QQ(MODE1)=XQ1(MX1)

        MH2=0
        MJ2=0
        INCR2=0
C************************************************INCREMENT
        INC2=1
        IF(MOLPRO.GT.0.AND.MOD(MID2,2).EQ.0)INC2=MOLINC
        DO M2=1,MMM2
          INCR2=INCR2+1
          IF(INCR2.EQ.INC2)THEN
            INCR2=0
          ELSE
            GO TO 2000
          END IF
          INC2=INCR
          IF(M2.EQ.MID2)INC2=1
          MX2=M2
          IF(ISKIP2.NE.0.AND.M2.GT.MID2)GO TO 2000
          MH2=MH2+1
          QQ(MODE2)=XQ2(MX2)

          MH3=0
          MJ3=0
          INCR3=0
C************************************************INCREMENT
          INC3=1
          IF(MOLPRO.GT.0.AND.MOD(MID3,2).EQ.0)INC3=MOLINC
          DO M3=1,MMM3
            INCR3=INCR3+1
            IF(INCR3.EQ.INC3)THEN
              INCR3=0
            ELSE
              GO TO 3000
            END IF
            INC3=INCR
            IF(M3.EQ.MID3)INC3=1
            MX3=M3
            IF(ISKIP3.NE.0.AND.M3.GT.MID3)GO TO 3000
            MH3=MH3+1
            QQ(MODE3)=XQ3(MX3)

            MH4=0
            MJ4=0
            INCR4=0
C************************************************INCREMENT
            INC4=1
            IF(MOLPRO.GT.0.AND.MOD(MID4,2).EQ.0)INC4=MOLINC
            DO M4=1,MMM4
              INCR4=INCR4+1
              IF(INCR4.EQ.INC4)THEN
                INCR4=0
              ELSE
                GO TO 4000
              END IF
              INC4=INCR
              IF(M4.EQ.MID4)INC4=1
              MX4=M4
              IF(ISKIP4.NE.0.AND.M4.GT.MID4)GO TO 4000
              MH4=MH4+1
              QQ(MODE4)=XQ4(MX4)

              MH5=0
              MJ5=0 
              INCR5=0
C************************************************INCREMENT
              INC5=1
              IF(MOLPRO.GT.0.AND.MOD(MID5,2).EQ.0)INC5=MOLINC
              DO M5=1,MMM5
                INCR5=INCR5+1
                IF(INCR5.EQ.INC5)THEN
                INCR5=0
              ELSE
                GO TO 5000
              END IF
              INC5=INCR
              IF(M5.EQ.MID5)INC5=1
              MX5=M5
              IF(ISKIP5.NE.0.AND.M5.GT.MID5)GO TO 5000
              MH5=MH5+1
              QQ(MODE5)=XQ5(MX5)
              
              IF(IABS(MOLPRO).EQ.9)THEN

                IF(IWHICH.GT.0)THEN
C**FITTING GLOBAL POTENTIAL
C
C**FIRST GET POTENTIAL AT M1,M2,M3,M4,M5
                  DO I=1,NATOM
                    DO K=1,3
                      XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1                SQRT(XM(I))
                      XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
     1                SQRT(XM(I)) 
                    END DO
                  END DO
                  CALL GETPOT(VVV,NATOM,XX,ENERGY)
C**SECOND GET DERIVATIVE WRT MODE1
C                 IF(M1.EQ.1.OR.M1.LT.MMM1)THEN
C                   QQQ=QQ(MODE1)+1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
C    1                  SQRT(XM(I)) 
C                     END DO
C                   END DO
C                   CALL GETPOT(VVP1,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M1.EQ.1)DD1=(VVP1-VVV)*1.D3
C                 IF(M1.EQ.MMM1.OR.M1.GT.1)THEN
C                   QQQ=QQ(MODE1)-1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
C    1                  SQRT(XM(I)) 
C                     END DO
C                   END DO
C                   CALL GETPOT(VVM1,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M1.EQ.MMM1)DD1=(VVV-VVM1)*1.D3
C                 IF(M1.GT.1.AND.M1.LT.MMM1)DD1=(VVP1-VVM1)*1.D3/2
C**THIRD GET DERIVATIVE WRT MODE2
C                 IF(M2.EQ.1.OR.M2.LT.MMM2)THEN
C                   QQQ=QQ(MODE2)+1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE2,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
C    1                  SQRT(XM(I))
C                     END DO
C                   END DO
C                   CALL GETPOT(VVP2,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M2.EQ.1)DD2=(VVP2-VVV)*1.D3
C                 IF(M2.EQ.MMM2.OR.M2.GT.1)THEN
C                   QQQ=QQ(MODE2)-1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE2,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
C    1                  SQRT(XM(I))
C                     END DO
C                   END DO
C                   CALL GETPOT(VVM2,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M2.EQ.MMM2)DD2=(VVV-VVM2)*1.D3
C                 IF(M2.GT.1.AND.M2.LT.MMM2)DD2=(VVP2-VVM2)*1.D3/2
C**FOURTH GET DERIVATIVE WRT MODE3
C                 IF(M3.EQ.1.OR.M3.LT.MMM3)THEN
C                   QQQ=QQ(MODE3)+1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE3,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
C    1                  SQRT(XM(I))
C                     END DO
C                   END DO
C                   CALL GETPOT(VVP3,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M3.EQ.1)DD3=(VVP3-VVV)*1.D3
C                 IF(M3.EQ.MMM3.OR.M3.GT.1)THEN
C                   QQQ=QQ(MODE3)-1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE3,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
C    1                  SQRT(XM(I))
C                     END DO
C                   END DO
C                   CALL GETPOT(VVM3,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M3.EQ.MMM3)DD3=(VVV-VVM3)*1.D3
C                 IF(M3.GT.1.AND.M3.LT.MMM3)DD3=(VVP3-VVM3)*1.D3/2
C**FIFTH GET DERIVATIVE WRT MODE4
C                 IF(M4.EQ.1.OR.M4.LT.MMM4)THEN
C                   QQQ=QQ(MODE4)+1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE4,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
C    1                  SQRT(XM(I))     
C                     END DO
C                   END DO
C                   CALL GETPOT(VVP4,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M4.EQ.1)DD4=(VVP4-VVV)*1.D3
C                 IF(M4.EQ.MMM4.OR.M4.GT.1)THEN
C                   QQQ=QQ(MODE4)-1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE4,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
C    1                  SQRT(XM(I))     
C                     END DO
C                   END DO
C                   CALL GETPOT(VVM4,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M4.EQ.MMM4)DD4=(VVV-VVM4)*1.D3
C                 IF(M4.GT.1.AND.M4.LT.MMM4)DD4=(VVP4-VVM4)*1.D3/2
C**SIXTH GET DERIVATIVE WRT MODE5
C                 IF(M5.EQ.1.OR.M5.LT.MMM5)THEN
C                   QQQ=QQ(MODE5)+1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE5,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
C    1                  SQRT(XM(I))     
C                     END DO
C                   END DO
C                   CALL GETPOT(VVP5,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M5.EQ.1)DD5=(VVP5-VVV)*1.D3
C                 IF(M5.EQ.MMM5.OR.M5.GT.1)THEN
C                   QQQ=QQ(MODE5)-1.D-3
C                   DO I=1,NATOM
C                     DO K=1,3
C                       XX(I,K)=X0(I,K)+XL(I,MODE5,K)*QQQ/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
C    1                  SQRT(XM(I))
C                       XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
C    1                  SQRT(XM(I))     
C                     END DO
C                   END DO
C                   CALL GETPOT(VVM5,NATOM,XX,ENERGY)
C                 END IF
C                 IF(M5.EQ.MMM5)DD5=(VVV-VVM5)*1.D3
C                 IF(M5.GT.1.AND.M5.LT.MMM5)DD5=(VVP5-VVM5)*1.D3/2
                END IF

                WRITE(MOUT,204)NATOM
                WRITE(MOUT,205)MODE1,CHSYM(MSYM1),M1,MODE2,
     1          CHSYM(MSYM2),M2,MODE3,CHSYM(MSYM3),M3,MODE4,
     2          CHSYM(MSYM4),M4,MODE5,CHSYM(MSYM5),M5

                IF(IWHICH.GT.0)THEN
C**FITTING GLOBAL POTENTIAL
C                 WRITE(MOUT,208)VVV,DD1,DD2,DD3,DD4,DD5
                  WRITE(MOUT,208)VVV
                ELSE
C**FITTING ABINITIO POTENTIAL
                  WRITE(MOUT,*)'**********************************'
                END IF

              END IF 

              IF(MOLPRO.EQ.10)THEN
C**FIT ABINITIO POTENTIAL
                READ(MINP,204)NATOM
                READ(MINP,205)MODE1,CHSYM(MSYM1),MD1,MODE2,
     1          CHSYM(MSYM2),MD2,MODE3,CHSYM(MSYM3),MD3,MODE4,
     2          CHSYM(MSYM4),MD4,MODE5,CHSYM(MSYM5),MD5
                READ(MINP,*)ENERGY(MH5,MH4,MH3,MH2,MH1,1)
              END IF
              DO I=1,NATOM
                DO K=1,3
                  XX(I,K)=X0(I,K)+XL(I,MODE1,K)*QQ(MODE1)/
     1            SQRT(XM(I))
                  XX(I,K)=XX(I,K)+XL(I,MODE2,K)*QQ(MODE2)/
     1            SQRT(XM(I))
                  XX(I,K)=XX(I,K)+XL(I,MODE3,K)*QQ(MODE3)/
     1            SQRT(XM(I))
                  XX(I,K)=XX(I,K)+XL(I,MODE4,K)*QQ(MODE4)/
     1            SQRT(XM(I))
                  XX(I,K)=XX(I,K)+XL(I,MODE5,K)*QQ(MODE5)/
     1            SQRT(XM(I))
                END DO
                IF(IABS(MOLPRO).EQ.9)THEN
                  WRITE(MOUT,206)SYMBOL(I),(XX(I,K)*BOHR,K=1,3)
                END IF
                IF(IABS(MOLPRO).EQ.10)THEN
                  READ(MINP,206)SYMBAD,(XX(I,K),K=1,3)
                  DO K=1,3
                    XX(I,K)=XX(I,K)/BOHR
                  END DO
                END IF
              END DO
C**TAKE OUT 1-DIM, 2-DIM, 3-DIM, 4-DIM POTENTIALS + CONSTANT
              IF(MOLPRO.EQ.10)THEN
                CALL GETPQ4(VQ4,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1          NP2,CP2,IP2,NTOT2,MAX2,NP3,CP3,IP3,NTOT3,MAX3,NP4,CP4,
     2          IP4,NTOT4,MAX4,INDK,INDL,INDN)
                ENERGY(MH5,MH4,MH3,MH2,MH1,2)=ENERGY(MH5,MH4,MH3,MH2,
     1          MH1,1)-VQ4
              END IF
5000            CONTINUE
              END DO
4000          CONTINUE
            END DO
3000        CONTINUE
          END DO
2000      CONTINUE
        END DO
1000    CONTINUE
      END DO
      MFIT=MH1*MH2*MH3*MH4*MH5
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETQ5(XK,EVAL,SOL,MM1,XQ1,MM2,XQ2,MM3,XQ3,MM4,XQ4,
     1MM5,XQ5,XV,MODE1,MODE2,MODE3,MODE4,MODE5,MSYM1,MSYM2,MSYM3,MSYM4,
     2MSYM5,NMODE,QQ,XTANPM,WRK,MDIM,VINF,D,NP1,CP1,IP1,NTOT1,MMX1,
     3NP2,CP2,IP2,NTOT2,MAX2,NP3,CP3,IP3,NTOT3,MAX3,NP4,CP4,IP4,NTOT4,
     4MAX4,INDK,INDL,INDN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XK(MDIM,MDIM),EVAL(MDIM,7)
      DIMENSION SOL(MDIM,1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQ5(MM5)
      DIMENSION XV(MM5,MM4,MM3,MM2,MM1,2)
      DIMENSION QQ(NMODE),VINF(NMODE)
      DIMENSION XTANPM(NMODE),WRK(1)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),IP3(MAX3,NTOT3,3)
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),IP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)


      DIMENSION D(MDIM,MDIM)
      COMMON/VMIN/VMIN,VMAX


      COMMON/SYMMP/ISYMP(10,10)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FITTER/MFIT,MN1,MMM1,ISYM1,ISKIP1,XTAN1,MN2,MMM2,ISYM2,
     1ISKIP2,XTAN2,MN3,MMM3,ISYM3,ISKIP3,XTAN3,MN4,MMM4,ISYM4,ISKIP4,
     2XTAN4,MN5,MMM5,ISYM5,ISKIP5,XTAN5
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/RMS/RMS
      COMMON/FILASS/IOUT,INP
100   FORMAT(/,1X,'DGESV IFAIL = ',I3)
101   FORMAT(D20.10,5X,I3,2X,I3,2X,I3,2X,I3,2X,I3)
102   FORMAT(//,1X,'ENTERING EXACT FIT ROUTINE',/)
103   FORMAT(//,1X,'ENTERING LEAST SQUARES FIT ROUTINE',/,
     11X,'MFIT = ',I4,' NFIT = ',I4,/)
105   FORMAT(/,7X,'COEFFICIENT',9X,'POWERS',/)
108   FORMAT(/,1X,'FITTING OF POTENTIAL FOR MODES: ',I3,2X,I3,2X,I3,2X,
     1I3,2X,I3)
109   FORMAT(//,7X,'OBSERVED',7X,'CALCULATED',/)
110   FORMAT(1X,F12.4,5X,F12.4)
111   FORMAT(1X,'SYMMETRIES ',A2,2X,A2,2X,A2,2X,A2,2X,A2)
112   FORMAT(//,7X,'OBSERVED',7X,'CALCULATED',9X,'TOTAL',/)
113   FORMAT(1X,F12.4,5X,F12.4,5X,F12.4)
C***********************************************
      DO K=1,NMODE
        QQ(K)=0
      END DO
      WRITE(IOUT,108)MODE1,MODE2,MODE3,MODE4,MODE5
      WRITE(IOUT,111)CHSYM(MSYM1),CHSYM(MSYM2),CHSYM(MSYM3),
     1CHSYM(MSYM4),CHSYM(MSYM5)
      MMM1=MM1
      MMM2=MM2
      MMM3=MM3
      MMM4=MM4
      MMM5=MM5
C************************************************INCREMENT
      MX1=MM1/MOLINC
      MX2=MM2/MOLINC
      MX3=MM3/MOLINC
      MX4=MM4/MOLINC
      MX5=MM5/MOLINC
      MMM1=MX1+MOD(MX1,MOLINC)
      MMM2=MX2+MOD(MX2,MOLINC)
      MMM3=MX3+MOD(MX3,MOLINC)
      MMM4=MX4+MOD(MX4,MOLINC)
      MMM5=MX5+MOD(MX5,MOLINC)
      ISYM1=MSYM1
      ISYM2=MSYM2
      ISYM3=MSYM3
      ISYM4=MSYM4
      ISYM5=MSYM5
      XTAN1=XTANPM(MODE1)
      XTAN2=XTANPM(MODE2)
      XTAN3=XTANPM(MODE3)
      XTAN4=XTANPM(MODE4)
      XTAN5=XTANPM(MODE5)

C**POLYNOMIAL DEPENDS ON SYMMETRY 
      NDIM=MM1*MM2*MM3*MM4*MM5
      N=MMM1*MMM2*MMM3*MMM4*MMM5
      NM=N
      IF(MSYM1.NE.1)NM=N/2
      IF(MSYM2.NE.1)NM=N/2
      IF(MSYM3.NE.1)NM=N/2
      IF(MSYM4.NE.1)NM=N/2
      IF(MSYM5.NE.1)NM=N/2
      IF(MSYM1.NE.1.AND.MSYM2.NE.1.AND.MSYM1.NE.MSYM2)NM=N/4
      IF(MSYM1.NE.1.AND.MSYM3.NE.1.AND.MSYM1.NE.MSYM3)NM=N/4
      IF(MSYM1.NE.1.AND.MSYM4.NE.1.AND.MSYM1.NE.MSYM4)NM=N/4
      IF(MSYM1.NE.1.AND.MSYM5.NE.1.AND.MSYM1.NE.MSYM5)NM=N/4
      IF(MSYM2.NE.1.AND.MSYM3.NE.1.AND.MSYM2.NE.MSYM3)NM=N/4
      IF(MSYM2.NE.1.AND.MSYM4.NE.1.AND.MSYM2.NE.MSYM4)NM=N/4
      IF(MSYM2.NE.1.AND.MSYM5.NE.1.AND.MSYM2.NE.MSYM5)NM=N/4
      IF(MSYM3.NE.1.AND.MSYM4.NE.1.AND.MSYM3.NE.MSYM4)NM=N/4
      IF(MSYM3.NE.1.AND.MSYM5.NE.1.AND.MSYM3.NE.MSYM5)NM=N/4
      IF(MSYM4.NE.1.AND.MSYM5.NE.1.AND.MSYM4.NE.MSYM5)NM=N/4

      MN1=1
      MN2=1
      MN3=1
      MN4=1
      MN5=1
      IF(MSYM1.NE.1)MN1=2
      IF(MSYM2.NE.1)MN2=2
      IF(MSYM3.NE.1)MN3=2
      IF(MSYM4.NE.1)MN4=2
      IF(MSYM5.NE.1)MN5=2
C**SPECIAL CASE B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2, ETC.
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM2))MN1=1
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM3))MN1=1
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM4))MN1=1
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM5))MN1=1
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM3))MN2=1
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM4))MN2=1
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM5))MN2=1
      IF(MSYM3.NE.1.AND.(MSYM3.EQ.MSYM4))MN3=1
      IF(MSYM3.NE.1.AND.(MSYM3.EQ.MSYM5))MN3=1
      IF(MSYM4.NE.1.AND.(MSYM4.EQ.MSYM5))MN4=1
      
      N12=ISYMP(MSYM1,MSYM2)
      N13=ISYMP(MSYM1,MSYM3)
      N14=ISYMP(MSYM1,MSYM4)
      N15=ISYMP(MSYM1,MSYM5)
      N23=ISYMP(MSYM2,MSYM3)
      N24=ISYMP(MSYM2,MSYM4)
      N25=ISYMP(MSYM2,MSYM5)
      N34=ISYMP(MSYM3,MSYM4)
      N35=ISYMP(MSYM3,MSYM5)
      N45=ISYMP(MSYM4,MSYM5)
      
C**SPECIAL CASE S2*S3 = S1
      IF(N23.EQ.MSYM1)MN1=1
C**SPECIAL CASE S2*S4 = S1
      IF(N24.EQ.MSYM1)MN1=1
C**SPECIAL CASE S2*S5 = S1
      IF(N25.EQ.MSYM1)MN1=1
C**SPECIAL CASE S3*S4 = S1
      IF(N34.EQ.MSYM1)MN1=1
C**SPECIAL CASE S3*S5 = S1
      IF(N35.EQ.MSYM1)MN1=1
C**SPECIAL CASE S4*S5 = S1
      IF(N45.EQ.MSYM1)MN1=1
C**SPECIAL CASE S3*S4 = S2
      IF(N34.EQ.MSYM2)MN2=1
C**SPECIAL CASE S3*S5 = S2
      IF(N35.EQ.MSYM2)MN2=1
C**SPECIAL CASE S4*S5 = S2
      IF(N45.EQ.MSYM2)MN2=1            
C**SPECIAL CASE S4*S5 = S3
      IF(N45.EQ.MSYM3)MN3=1     
      
      N234=ISYMP(MSYM2,N34)
      N235=ISYMP(MSYM2,N35)
      N245=ISYMP(MSYM2,N45)
      N345=ISYMP(MSYM3,N45)
      
C**SPECIAL CASE S2*S3*S4 = S1
      IF(N234.EQ.MSYM1)MN1=1
C**SPECIAL CASE S2*S3*S5 = S1
      IF(N235.EQ.MSYM1)MN1=1      
C**SPECIAL CASE S2*S4*S5 = S1
      IF(N245.EQ.MSYM1)MN1=1
C**SPECIAL CASE S3*S4*S5 = S1
      IF(N345.EQ.MSYM1)MN1=1      
C**SPECIAL CASE S3*S4*S5 = S2      
      IF(N345.EQ.MSYM2)MN2=1      
      
      N2345=ISYMP(N23,N45)      
      
c**SPECIAL CASE S2*S3*S4*S5 = S1
      IF(N2345.EQ.MSYM1)MN1=1      
      
      ISKIP1=0
      IF(MSYM1.NE.1)ISKIP1=1
      ISKIP2=0
      IF(MSYM2.NE.1)ISKIP2=1
      ISKIP3=0
      IF(MSYM3.NE.1)ISKIP3=1
      ISKIP4=0
      IF(MSYM4.NE.1)ISKIP4=1
      ISKIP5=0
      IF(MSYM5.NE.1)ISKIP5=1
C**SPECIAL CASE B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2, ETC.
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2, ETC.
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM2))ISKIP1=0
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM3))ISKIP1=0
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM4))ISKIP1=0
      IF(MSYM1.NE.1.AND.(MSYM1.EQ.MSYM5))ISKIP1=0
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM3))ISKIP2=0
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM4))ISKIP2=0
      IF(MSYM2.NE.1.AND.(MSYM2.EQ.MSYM5))ISKIP2=0
      IF(MSYM3.NE.1.AND.(MSYM3.EQ.MSYM4))ISKIP3=0
      IF(MSYM3.NE.1.AND.(MSYM3.EQ.MSYM5))ISKIP3=0
      IF(MSYM4.NE.1.AND.(MSYM4.EQ.MSYM5))ISKIP4=0

C**SPECIAL CASE S2*S3 = S1
      IF(N23.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S2*S4 = S1
      IF(N24.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S2*S5 = S1
      IF(N25.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S4 = S1
      IF(N34.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S5 = S1
      IF(N35.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S4*S5 = S1
      IF(N45.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S4 = S2
      IF(N34.EQ.MSYM2)ISKIP2=0
C**SPECIAL CASE S3*S5 = S2
      IF(N35.EQ.MSYM2)ISKIP2=0
C**SPECIAL CASE S4*S5 = S2
      IF(N45.EQ.MSYM2)ISKIP2=0
C**SPECIAL CASE S4*S5 = S3
      IF(N45.EQ.MSYM3)ISKIP3=0      
      
C**SPECIAL CASE S2*S3*S4 = S1
      IF(N234.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S2*S3*S5 = S1
      IF(N235.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S2*S4*S5 = S1
      IF(N245.EQ.MSYM1)ISKIP1=0
C**SPECIAL CASE S3*S4*S5 = S1
      IF(N345.EQ.MSYM1)ISKIP1=0            
C**SPECIAL CASE S3*S4*S5 = S2
      IF(N345.EQ.MSYM2)ISKIP2=0      

C**SPECIAL CASE S2*S3*S4*S5 = S1
      IF(N2345.EQ.MSYM1)ISKIP1=0
 
      IOFF5=0
      IOFF4=0
      IOFF3=0
      IOFF2=0
      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 5001
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 5002
          MH3=MMM3/2
          MX3=MM3/2
          IS3=-1
          INH3=0
          INC3=0
          DO M3=1,MMM3
            MH3=MH3+INH3*IS3
            MX3=MX3+INC3*IS3
            IS3=-IS3
            INH3=INH3+1
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 5003
            MH4=MMM4/2
            MX4=MM4/2
            IS4=-1
            INH4=0
            INC4=0
            DO M4=1,MMM4
              MH4=MH4+INH4*IS4
              MX4=MX4+INC4*IS4
              IS4=-IS4
              INH4=INH4+1
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC4=INC4+1
              ELSE
                INC4=2*M4-1
              END IF
              IF(ISKIP4.NE.0.AND.MOD(M4,2).EQ.0)GO TO 5004

            MH5=MMM5/2
            MX5=MM5/2
            IS5=-1
            INH5=0
            INC5=0
            DO M5=1,MMM5
              MH5=MH5+INH5*IS5
              MX5=MX5+INC5*IS5
              IS5=-IS5
              INH5=INH5+1
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC5=INC5+1
              ELSE
                INC5=2*M5-1
              END IF
              IF(ISKIP5.NE.0.AND.MOD(M5,2).EQ.0)GO TO 5005
              I=I+1
              EVAL(I,1)=XV(MH5,MH4,MH3,MH2,MH1,1)
              EVAL(I,2)=XV(MH5,MH4,MH3,MH2,MH1,2)
C**WEIGHT (MINIMUM SCALED TO UNITY)
              VSCAL=WAVENM*(EVAL(I,1)-VMIN)+1
              EVAL(I,3)=1/VSCAL

              J=0
              DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1          IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1          IOFF4=MOD(J1,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM5.AND.MN5.EQ.2)
     1          IOFF5=MOD(J1,2)+1-MN5     
                DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1            IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1            IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM5.AND.MN5.EQ.2)
     1            IOFF5=MOD(J2,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1            IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1            IOFF4=MOD(J1+J2,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM5)
     1            IOFF5=MOD(J1+J2,2)-MN5
                  DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1              IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM5.AND.MN5.EQ.2)
     1              IOFF5=MOD(J3,2)+1-MN5     
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM4)IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     2              IOFF4=MOD(J1+J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM5)
     2              IOFF5=MOD(J1+J2+J3,2)-MN5
                    DO J4=1,MMM4/MN4
C**SPECIAL CASE B2 WITH B2
                      IF(MSYM4.EQ.MSYM5.AND.MN5.EQ.2)
     1                IOFF5=MOD(J4,2)+1-MN5                     
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM4.AND.
     1                MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM4.AND.
     1                MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM3.NE.1.AND.MSYM3.EQ.MSYM4.AND.
     1                MSYM3.EQ.MSYM5)IOFF5=MOD(J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1                MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2                IOFF5=MOD(J1+J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1                MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2                IOFF5=MOD(J1+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1                MSYM2.EQ.MSYM4.AND.MSYM2.EQ.MSYM5)
     2                IOFF5=MOD(J2+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1                MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4.AND.
     2                MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J3+J4,2)-MN5
                      DO J5=1,MMM5/MN5
                 J=J+1
                 XK(I,J)=(XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     1           (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     2           (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))*
     3           (XTANH(XTANPM(MODE4)*XQ4(MX4))**(MN4*(J4-1)-IOFF4))*
     4           (XTANH(XTANPM(MODE5)*XQ5(MX5))**(MN5*(J5-1)-IOFF5))
                      END DO
                    END DO
                  END DO
                END DO
              END DO
5005          CONTINUE
            ENDDO
5004          CONTINUE
            END DO
5003        CONTINUE
          END DO
5002      CONTINUE
        END DO
5001    CONTINUE
      END DO
C******************************
      WRITE(IOUT,102)
      IFAIL=1
C     CALL DGESV(NM,1,XK,NDIM,SOL,EVAL,N,IFAIL)
      CALL DGESV(NM,1,XK,MDIM,SOL,EVAL,NM,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C******************************
      DO I=1,NM
        SOL(I,1)=EVAL(I,1)
      END DO
      WRITE(IOUT,105)
      J=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1  IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1  IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1  IOFF4=MOD(J1,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM5.AND.MN5.EQ.2)
     1  IOFF5=MOD(J1,2)+1-MN5     
        DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1    IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1    IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM5.AND.MN5.EQ.2)
     1    IOFF5=MOD(J2,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1    IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1    IOFF4=MOD(J1+J2,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM5)
     1    IOFF5=MOD(J1+J2,2)-MN5
          DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
            IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1      IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
            IF(MSYM3.EQ.MSYM5.AND.MN5.EQ.2)
     1      IOFF5=MOD(J3,2)+1-MN5     
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     1      IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM5)
     1      IOFF5=MOD(J1+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM4)
     1      IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM5)
     1      IOFF5=MOD(J2+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1      MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1      MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J3,2)-MN5
            DO J4=1,MMM4/MN4
C**SPECIAL CASE B2 WITH B2
              IF(MSYM4.EQ.MSYM5.AND.MN5.EQ.2)
     1        IOFF5=MOD(J4,2)+1-MN5                     
C**SPECIAL CASE B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     1        IOFF5=MOD(J1+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
              IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM4.AND.MSYM2.EQ.MSYM5)
     1        IOFF5=MOD(J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
              IF(MSYM3.NE.1.AND.MSYM3.EQ.MSYM4.AND.MSYM3.EQ.MSYM5)
     1        IOFF5=MOD(J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4.AND.
     1        MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4.AND.
     1        MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
              IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM4.AND.
     1        MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1        MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2        IOFF5=MOD(J1+J2+J3+J4,2)-MN5
              DO J5=1,MMM5/MN5
                J=J+1
                WRITE(IOUT,101)SOL(J,1),MN1*(J1-1),MN2*(J2-1)-IOFF2,
     1          MN3*(J3-1)-IOFF3,MN4*(J4-1)-IOFF4,MN5*(J5-1)-IOFF5
              END DO
            END DO
          END DO
        END DO
      END DO
      CALL FLUSH(IOUT)

      WRITE(IOUT,109)
      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 6001
        QQ(MODE1)=XQ1(MX1)
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 6002
          QQ(MODE2)=XQ2(MX2)
          MH3=MMM3/2
          MX3=MM3/2
          IS3=-1
          INH3=0
          INC3=0
          DO M3=1,MMM3
            MH3=MH3+INH3*IS3
            MX3=MX3+INC3*IS3
            IS3=-IS3
            INH3=INH3+1
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 6003
            QQ(MODE3)=XQ3(MX3)
            MH4=MMM4/2
            MX4=MM4/2
            IS4=-1
            INH4=0
            INC4=0
            DO M4=1,MMM4
              MH4=MH4+INH4*IS4
              MX4=MX4+INC4*IS4
              IS4=-IS4
              INH4=INH4+1
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC4=INC4+1
              ELSE
                INC4=2*M4-1
              END IF
              IF(ISKIP4.NE.0.AND.MOD(M4,2).EQ.0)GO TO 6004
              QQ(MODE4)=XQ4(MX4)
              MH5=MMM5/2
              MX5=MM5/2
              IS5=-1
              INH5=0
              INC5=0
              DO M5=1,MMM5
                MH5=MH5+INH5*IS5
                MX5=MX5+INC5*IS5
                IS5=-IS5
                INH5=INH5+1
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC5=INC5+1
              ELSE
                INC5=2*M5-1
              END IF
              IF(ISKIP5.NE.0.AND.MOD(M5,2).EQ.0)GO TO 6005
              QQ(MODE5)=XQ5(MX5)
              I=I+1
              EVAL(I,1)=0
              J=0
              DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1          IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1          IOFF4=MOD(J1,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM5.AND.MN5.EQ.2)
     1          IOFF5=MOD(J1,2)+1-MN5
                DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1            IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1            IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM5.AND.MN5.EQ.2)
     1            IOFF5=MOD(J2,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1            IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1            IOFF4=MOD(J1+J2,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM5)
     1            IOFF5=MOD(J1+J2,2)-MN5
                  DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1              IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM5.AND.MN5.EQ.2)
     1              IOFF5=MOD(J3,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM4)IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     2              IOFF4=MOD(J1+J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM5)
     2              IOFF5=MOD(J1+J2+J3,2)-MN5
                    DO J4=1,MMM4/MN4
C**SPECIAL CASE B2 WITH B2
                      IF(MSYM4.EQ.MSYM5.AND.MN5.EQ.2)
     1                IOFF5=MOD(J4,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM4.AND.
     1                MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM4.AND.
     1                MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM3.NE.1.AND.MSYM3.EQ.MSYM4.AND.
     1                MSYM3.EQ.MSYM5)IOFF5=MOD(J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1                MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2                IOFF5=MOD(J1+J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1                MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2                IOFF5=MOD(J1+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1                MSYM2.EQ.MSYM4.AND.MSYM2.EQ.MSYM5)
     2                IOFF5=MOD(J2+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1                MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4.AND.
     2                MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J3+J4,2)-MN5
                      DO J5=1,MMM5/MN5
                        J=J+1
                        EVAL(I,1)=EVAL(I,1)+SOL(J,1)*
     1            (XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     2            (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     3            (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))*
     4            (XTANH(XTANPM(MODE4)*XQ4(MX4))**(MN4*(J4-1)-IOFF4))*
     5            (XTANH(XTANPM(MODE5)*XQ5(MX5))**(MN5*(J5-1)-IOFF5))
     
                      END DO
                    END DO
                  END DO
                END DO
              END DO
C**ADD IN 1-DIM, 2-DIM, 3-DIM, 4-DIM POTENTIALS + CONSTANT
              WRITE(IOUT,110)WAVENM*XV(MH5,MH4,MH3,MH2,MH1,1),
     1        WAVENM*EVAL(I,1)
     
6005            CONTINUE
              END DO
6004          CONTINUE
            END DO
6003        CONTINUE
          END DO
6002      CONTINUE
        END DO
6001    CONTINUE
      END DO

C******************************
C**LEAST SQUARES
C******************************
      NFIT=0
      J=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1  IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1  IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1  IOFF4=MOD(J1,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM5.AND.MN5.EQ.2)
     1  IOFF5=MOD(J1,2)+1-MN5
        DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1    IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1    IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM5.AND.MN5.EQ.2)
     1    IOFF5=MOD(J2,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1    IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1    IOFF4=MOD(J1+J2,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM5)
     1    IOFF5=MOD(J1+J2,2)-MN5
          DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
            IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1      IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
            IF(MSYM3.EQ.MSYM5.AND.MN5.EQ.2)
     1      IOFF5=MOD(J3,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     1      IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM5)
     1      IOFF5=MOD(J1+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM4)
     1      IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM5)
     1      IOFF5=MOD(J2+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1      MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1      MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J3,2)-MN5
            DO J4=1,MMM4/MN4
C**SPECIAL CASE B2 WITH B2
              IF(MSYM4.EQ.MSYM5.AND.MN5.EQ.2)
     1        IOFF5=MOD(J4,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     1        IOFF5=MOD(J1+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
              IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM4.AND.MSYM2.EQ.MSYM5)
     1        IOFF5=MOD(J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
              IF(MSYM3.NE.1.AND.MSYM3.EQ.MSYM4.AND.MSYM3.EQ.MSYM5)
     1        IOFF5=MOD(J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4.AND.
     1        MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4.AND.
     1        MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
              IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM4.AND.
     1        MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1        MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2        IOFF5=MOD(J1+J2+J3+J4,2)-MN5
              DO J5=1,MMM5/MN5
                J=J+1
                IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1          MN3*(J3-1)-IOFF3.NE.0.AND.MN4*(J4-1)-IOFF4.NE.0.AND.
     2          MN5*(J5-1)-IOFF5.NE.0)THEN
                  NFIT=NFIT+1
C**COPY ONLY RELEVANT COEFFICIENTS FOR Vijkl(5)
                  SOL(NFIT,1)=SOL(J,1)
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO
      WRITE(IOUT,103)MFIT,NFIT


C**A-MATRIX
      I=0
C******************************
      MX1=MM1/2
      IS1=-1
      INC1=0
      DO M1=1,MMM1
        MX1=MX1+INC1*IS1
        IS1=-IS1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 1111
        MX2=MM2/2
        IS2=-1
        INC2=0
        DO M2=1,MMM2
          MX2=MX2+INC2*IS2
          IS2=-IS2
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 2222
          MX3=MM3/2
          IS3=-1
          INC3=0
          DO M3=1,MMM3
            MX3=MX3+INC3*IS3
            IS3=-IS3
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 3333
            MX4=MM4/2
            IS4=-1
            INC4=0
            DO M4=1,MMM4
              MX4=MX4+INC4*IS4
              IS4=-IS4
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC4=INC4+1
              ELSE
                INC4=2*M4-1
              END IF
              IF(ISKIP4.NE.0.AND.MOD(M4,2).EQ.0)GO TO 4444
              MX5=MM5/2
              IS5=-1
              INC5=0
              DO M5=1,MMM5
                MX5=MX5+INC5*IS5
                IS5=-IS5
C************************************************INCREMENT
                IF(MOLINC.EQ.1)THEN
                  INC5=INC5+1
                ELSE
                  INC5=2*M5-1
                END IF
                IF(ISKIP5.NE.0.AND.MOD(M5,2).EQ.0)GO TO 5555

              I=I+1
              J=0
              DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2 
                IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1          IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2 
                IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2 
                IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1          IOFF4=MOD(J1,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM5.AND.MN5.EQ.2)
     1          IOFF5=MOD(J1,2)+1-MN5     
                DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1            IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1            IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM5.AND.MN5.EQ.2)
     1            IOFF5=MOD(J2,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1            IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1            IOFF4=MOD(J1+J2,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM5)
     1            IOFF5=MOD(J1+J2,2)-MN5
                  DO J3=1,MMM3/MN3 
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1              IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 
                    IF(MSYM3.EQ.MSYM5.AND.MN5.EQ.2)
     1              IOFF5=MOD(J3,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM4)IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     2              IOFF4=MOD(J1+J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM5)
     2              IOFF5=MOD(J1+J2+J3,2)-MN5
                    DO J4=1,MMM4/MN4
C**SPECIAL CASE B2 WITH B2
                      IF(MSYM4.EQ.MSYM5.AND.MN5.EQ.2)
     1                IOFF5=MOD(J4,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM4.AND.
     1                MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM4.AND.
     1                MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM3.NE.1.AND.MSYM3.EQ.MSYM4.AND.
     1                MSYM3.EQ.MSYM5)IOFF5=MOD(J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1                MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2                IOFF5=MOD(J1+J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1                MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2                IOFF5=MOD(J1+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1                MSYM2.EQ.MSYM4.AND.MSYM2.EQ.MSYM5)
     2                IOFF5=MOD(J2+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1                MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4.AND.
     2                MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J3+J4,2)-MN5
                      DO J5=1,MMM5/MN5
                      IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1                MN3*(J3-1)-IOFF3.NE.0.AND.MN4*(J4-1)-IOFF4.NE.0.
     2                AND.MN5*(J5-1)-IOFF5.NE.0)THEN
                       J=J+1
                 XK(I,J)=(XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     1           (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     2           (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))*
     3           (XTANH(XTANPM(MODE4)*XQ4(MX4))**(MN4*(J4-1)-IOFF4))*
     4           (XTANH(XTANPM(MODE5)*XQ5(MX5))**(MN5*(J5-1)-IOFF5))
                      END IF
                      END DO
                    END DO
                  END DO
                END DO
              END DO
5555          CONTINUE              
            END DO
4444          CONTINUE
            END DO
3333        CONTINUE
          END DO
2222      CONTINUE
        END DO
1111    CONTINUE
      END DO
C**D-MATRIX
      DO J=1,NFIT
        DO K=1,NFIT
          X=0
          DO I=1,MFIT
            X=X+EVAL(I,3)*XK(I,K)*XK(I,J)
          END DO
          D(K,J)=X
        END DO
      END DO
C**E-MATRIX
      DO J=1,NFIT
        X=0
        DO I=1,MFIT
          X=X+EVAL(I,3)*XK(I,J)*EVAL(I,2)
        END DO
        EVAL(J,1)=X
      END DO

C******************************
      IFAIL=1
      CALL DGESV(NFIT,1,D,MDIM,SOL,EVAL,NFIT,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
C******************************
      DO I=1,NFIT
        SOL(I,1)=EVAL(I,1)
      END DO

      WRITE(IOUT,105)
      WRITE(MOUT,*)NFIT
      NFIT=0
      DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1  IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1  IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1  IOFF4=MOD(J1,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
        IF(MSYM1.EQ.MSYM5.AND.MN5.EQ.2)
     1  IOFF5=MOD(J1,2)+1-MN5
        DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1    IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1    IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
          IF(MSYM2.EQ.MSYM5.AND.MN5.EQ.2)
     1    IOFF5=MOD(J2,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1    IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1    IOFF4=MOD(J1+J2,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
          IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM5)
     1    IOFF5=MOD(J1+J2,2)-MN5
          DO J3=1,MMM3/MN3
C**SPECIAL CASE B2 WITH B2
            IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1      IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
            IF(MSYM3.EQ.MSYM5.AND.MN5.EQ.2)
     1      IOFF5=MOD(J3,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     1      IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM5)
     1      IOFF5=MOD(J1+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM4)
     1      IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
            IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM5)
     1      IOFF5=MOD(J2+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1      MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
            IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1      MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J3,2)-MN5
            DO J4=1,MMM4/MN4
C**SPECIAL CASE B2 WITH B2
              IF(MSYM4.EQ.MSYM5.AND.MN5.EQ.2)
     1        IOFF5=MOD(J4,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     1        IOFF5=MOD(J1+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
              IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM4.AND.MSYM2.EQ.MSYM5)
     1        IOFF5=MOD(J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
              IF(MSYM3.NE.1.AND.MSYM3.EQ.MSYM4.AND.MSYM3.EQ.MSYM5)
     1        IOFF5=MOD(J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4.AND.
     1        MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4.AND.
     1        MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
              IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.MSYM2.EQ.MSYM4.AND.
     1        MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2
              IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3.AND.
     1        MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2        IOFF5=MOD(J1+J2+J3+J4,2)-MN5
              DO J5=1,MMM5/MN5
                IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1          MN3*(J3-1)-IOFF3.NE.0.AND.MN4*(J4-1)-IOFF4.NE.0.AND.
     2          MN5*(J5-1)-IOFF5.NE.0)THEN
                NFIT=NFIT+1
                WRITE(IOUT,101)SOL(NFIT,1),MN1*(J1-1),MN2*(J2-1)-IOFF2,
     1          MN3*(J3-1)-IOFF3,MN4*(J4-1)-IOFF4,MN5*(J5-1)-IOFF5
                WRITE(MOUT,101)SOL(NFIT,1),MN1*(J1-1),MN2*(J2-1)-IOFF2,
     1          MN3*(J3-1)-IOFF3,MN4*(J4-1)-IOFF4,MN5*(J5-1)-IOFF5
                END IF
              
              END DO
            END DO
          END DO
        END DO
      END DO
      CALL FLUSH(IOUT)

      WRITE(IOUT,112)


      RMS=0


      I=0
      MH1=MMM1/2
      MX1=MM1/2
      IS1=-1
      INH1=0
      INC1=0
      DO M1=1,MMM1
        MH1=MH1+INH1*IS1
        MX1=MX1+INC1*IS1
        IS1=-IS1
        INH1=INH1+1
C************************************************INCREMENT
        IF(MOLINC.EQ.1)THEN
          INC1=INC1+1
        ELSE
          INC1=2*M1-1
        END IF
        IF(ISKIP1.NE.0.AND.MOD(M1,2).EQ.0)GO TO 7001
        QQ(MODE1)=XQ1(MX1)
        MH2=MMM2/2
        MX2=MM2/2
        IS2=-1
        INH2=0
        INC2=0
        DO M2=1,MMM2
          MH2=MH2+INH2*IS2
          MX2=MX2+INC2*IS2
          IS2=-IS2
          INH2=INH2+1
C************************************************INCREMENT
          IF(MOLINC.EQ.1)THEN
            INC2=INC2+1
          ELSE
            INC2=2*M2-1
          END IF
          IF(ISKIP2.NE.0.AND.MOD(M2,2).EQ.0)GO TO 7002
          QQ(MODE2)=XQ2(MX2)
          MH3=MMM3/2
          MX3=MM3/2
          IS3=-1
          INH3=0
          INC3=0
          DO M3=1,MMM3
            MH3=MH3+INH3*IS3
            MX3=MX3+INC3*IS3
            IS3=-IS3
            INH3=INH3+1
C************************************************INCREMENT
            IF(MOLINC.EQ.1)THEN
              INC3=INC3+1
            ELSE
              INC3=2*M3-1
            END IF
            IF(ISKIP3.NE.0.AND.MOD(M3,2).EQ.0)GO TO 7003
            QQ(MODE3)=XQ3(MX3)
            MH4=MMM4/2
            MX4=MM4/2
            IS4=-1
            INH4=0
            INC4=0
            DO M4=1,MMM4
              MH4=MH4+INH4*IS4
              MX4=MX4+INC4*IS4
              IS4=-IS4
              INH4=INH4+1
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC4=INC4+1
              ELSE
                INC4=2*M4-1
              END IF
              IF(ISKIP4.NE.0.AND.MOD(M4,2).EQ.0)GO TO 7004
              QQ(MODE4)=XQ4(MX4)

            MH5=MMM5/2
            MX5=MM5/2
            IS5=-1
            INH5=0
            INC5=0
            DO M5=1,MMM5
              MH5=MH5+INH5*IS5
              MX5=MX5+INC5*IS5
              IS5=-IS5
              INH5=INH5+1
C************************************************INCREMENT
              IF(MOLINC.EQ.1)THEN
                INC5=INC5+1
              ELSE
                INC5=2*M5-1
              END IF
              IF(ISKIP5.NE.0.AND.MOD(M5,2).EQ.0)GO TO 7005
              QQ(MODE5)=XQ5(MX5)


              I=I+1
              EVAL(I,1)=0
              NFIT=0
              DO J1=1,MMM1/MN1
C**SPECIAL CASE B2 WITH B2 
                IF(MSYM1.EQ.MSYM2.AND.MN2.EQ.2)
     1          IOFF2=MOD(J1,2)+1-MN2
C**SPECIAL CASE B2 WITH B2 
                IF(MSYM1.EQ.MSYM3.AND.MN3.EQ.2)
     1          IOFF3=MOD(J1,2)+1-MN3
C**SPECIAL CASE B2 WITH B2 
                IF(MSYM1.EQ.MSYM4.AND.MN4.EQ.2)
     1          IOFF4=MOD(J1,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                IF(MSYM1.EQ.MSYM5.AND.MN5.EQ.2)
     1          IOFF5=MOD(J1,2)+1-MN5     
                DO J2=1,MMM2/MN2
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM3.AND.MN3.EQ.2)
     1            IOFF3=MOD(J2,2)+1-MN3
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM4.AND.MN4.EQ.2)
     1            IOFF4=MOD(J2,2)+1-MN4
C**SPECIAL CASE B2 WITH B2
                  IF(MSYM2.EQ.MSYM5.AND.MN5.EQ.2)
     1            IOFF5=MOD(J2,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM3)
     1            IOFF3=MOD(J1+J2,2)-MN3
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM4)
     1            IOFF4=MOD(J1+J2,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                  IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.MSYM1.EQ.MSYM5)
     1            IOFF5=MOD(J1+J2,2)-MN5
                  DO J3=1,MMM3/MN3 
C**SPECIAL CASE B2 WITH B2
                    IF(MSYM3.EQ.MSYM4.AND.MN4.EQ.2)
     1              IOFF4=MOD(J3,2)+1-MN4
C**SPECIAL CASE B2 WITH B2 
                    IF(MSYM3.EQ.MSYM5.AND.MN5.EQ.2)
     1              IOFF5=MOD(J3,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM4)IOFF4=MOD(J1+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1              MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM4)IOFF4=MOD(J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2
                    IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1              MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J3,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4)
     2              IOFF4=MOD(J1+J2+J3,2)-MN4
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                    IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1              MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM5)
     2              IOFF5=MOD(J1+J2+J3,2)-MN5
                    DO J4=1,MMM4/MN4
C**SPECIAL CASE B2 WITH B2
                      IF(MSYM4.EQ.MSYM5.AND.MN5.EQ.2)
     1                IOFF5=MOD(J4,2)+1-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM4.AND.
     1                MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM4.AND.
     1                MSYM2.EQ.MSYM5)IOFF5=MOD(J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2
                      IF(MSYM3.NE.1.AND.MSYM3.EQ.MSYM4.AND.
     1                MSYM3.EQ.MSYM5)IOFF5=MOD(J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1                MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2                IOFF5=MOD(J1+J2+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM3.AND.
     1                MSYM1.EQ.MSYM4.AND.MSYM1.EQ.MSYM5)
     2                IOFF5=MOD(J1+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM2.NE.1.AND.MSYM2.EQ.MSYM3.AND.
     1                MSYM2.EQ.MSYM4.AND.MSYM2.EQ.MSYM5)
     2                IOFF5=MOD(J2+J3+J4,2)-MN5
C**SPECIAL CASE B2 WITH B2 WITH B2 WITH B2 WITH B2
                      IF(MSYM1.NE.1.AND.MSYM1.EQ.MSYM2.AND.
     1                MSYM1.EQ.MSYM3.AND.MSYM1.EQ.MSYM4.AND.
     2                MSYM1.EQ.MSYM5)IOFF5=MOD(J1+J2+J3+J4,2)-MN5
                      DO J5=1,MMM5/MN5
                      IF(MN1*(J1-1).NE.0.AND.MN2*(J2-1)-IOFF2.NE.0.AND.
     1                MN3*(J3-1)-IOFF3.NE.0.AND.MN4*(J4-1)-IOFF4.NE.0.
     2                AND.MN5*(J5-1)-IOFF5.NE.0)THEN
                        NFIT=NFIT+1
                        EVAL(I,1)=EVAL(I,1)+SOL(NFIT,1)*
     1             (XTANH(XTANPM(MODE1)*XQ1(MX1))**(MN1*(J1-1)))*
     2             (XTANH(XTANPM(MODE2)*XQ2(MX2))**(MN2*(J2-1)-IOFF2))*
     3             (XTANH(XTANPM(MODE3)*XQ3(MX3))**(MN3*(J3-1)-IOFF3))*
     4             (XTANH(XTANPM(MODE4)*XQ4(MX4))**(MN4*(J4-1)-IOFF4))*
     5             (XTANH(XTANPM(MODE5)*XQ5(MX5))**(MN5*(J5-1)-IOFF5))
                      END IF
                      END DO
                    END DO
                  END DO
                END DO
              END DO
C**ADD IN 1-DIM, 2-DIM, 3-DIM,4-DIM POTENTIALS + CONSTANT
              CALL GETPQ4(VQ4,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1        NP2,CP2,IP2,NTOT2,MAX2,NP3,CP3,IP3,NTOT3,MAX3,NP4,CP4,
     2        IP4,NTOT4,MAX4,INDK,INDL,INDN)
C             WRITE(IOUT,110)WAVENM*XV(MH5,MH4,MH3,MH2,MH1,1),
C    1        WAVENM*EVAL(I,1)
C    2        +WAVENM*VQ4
              WRITE(IOUT,113)WAVENM*XV(MH5,MH4,MH3,MH2,MH1,2),
     1        WAVENM*EVAL(I,1),WAVENM*XV(MH5,MH4,MH3,MH2,MH1,1)


C         RMS=RMS+((XV(MH5,MH4,MH3,MH2,MH1,2)-EVAL(I,1))*EVAL(I,3))**2
          RMS=RMS+(XV(MH5,MH4,MH3,MH2,MH1,2)-EVAL(I,1))**2
          
7005           CONTINUE
              END DO
7004         CONTINUE
            END DO
7003       CONTINUE
          END DO
7002     CONTINUE
        END DO
7001    CONTINUE
      END DO


      RMS=WAVENM*DSQRT(RMS/MFIT)
      WRITE(IOUT,*)'ROOT MEAN SQUARE (CM-1) = ',RMS


      RETURN
      END
C**************************************************************
C**************************************************************
      DOUBLE PRECISION FUNCTION XTANH(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      XPLUS=DEXP(X)
      XMIN=DEXP(-X)
      XTANH=(XPLUS-XMIN)/(XPLUS+XMIN)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE FITSIN(NMODE,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1NP2,CP2,IP2,NTOT2,MAX2,NP3,CP3,IP3,NTOT3,MAX3,NP4,CP4,IP4,NTOT4,
     2MAX4,NP5,CP5,IP5,NTOT5,MAX5,INDK,INDL,INDN,INDM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XTANPM(NMODE),INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),IP3(MAX3,NTOT3,3)
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),IP4(MAX4,NTOT4,4)
      DIMENSION NP5(NTOT5),CP5(MAX5,NTOT5),IP5(MAX5,NTOT5,5)
      COMMON/FILASS/IOUT,INP,MOUTIN
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/ENTER/IENTER,IENTMX(5)
      COMMON/AVCON/AVCON
      COMMON/ANALYT/MAXQU,MAXPOW
C*********************************************************************
      IF(MOLPRO.GT.2.OR.(MOLINC.GT.0.AND.IWHICH.LT.0.AND.ICOUPL.GT.0))
     1THEN
        REWIND MOUTIN
        IENTER=IENTER+1
        GO TO(1,2,12,12,12,12,12,12,12,12)IENTER
C**ONE-DIMENSIONAL FITS
1       CONTINUE
        AV=0
        READ(INP,*)
        IENTMX(1)=0
        DO I=1,NMODE
          READ(INP,*)XTANPM(I)
          READ(INP,*)NP1(I)
          IF(NP1(I).GT.IENTMX(1))IENTMX(1)=NP1(I)
          WRITE(MOUTIN,*)NP1(I)
          DO J=1,NP1(I)
            READ(INP,*)DUMMY1,IDUMMY
            IF(IDUMMY.EQ.0)AV=AV+DUMMY1
            WRITE(MOUTIN,*)DUMMY1,IDUMMY
          END DO
        END DO
        AVCON=AV/NMODE
        GO TO 12
2       DO I=1,NMODE
          READ(MOUTIN,*)NP1(I)
          DO J=1,NP1(I)
            READ(MOUTIN,*)CP1(J,I),IP1(J,I)
            IF(IP1(J,I).GT.MAXPOW)MAXPOW=IP1(J,I)
            IF(IP1(J,I).EQ.0)CP1(J,I)=AVCON/NMODE
          END DO
        END DO
        IF(MOLPRO.EQ.3.OR.MOLPRO.EQ.4.OR.(IWHICH.LT.0.AND.
     1  MOLINC.GT.0.AND.ICOUPL.EQ.1))IENTER=0
12      CONTINUE
      END IF
      IF(IENTER.LT.3)RETURN

      IF(MOLPRO.GT.4.OR.(MOLINC.GT.0.AND.IWHICH.LT.0.AND.ICOUPL.GT.1))
     1THEN
        GO TO(99,99,3,4,34,34,34,34,34,34)IENTER
C**TWO-DIMENSIONAL FITS
3       READ(INP,*)
        IENTMX(2)=0
        IND=0
        DO K=2,NMODE
          DO L=1,K-1
            IND=IND+1
            READ(INP,*)NP2(IND)
            IF(NP2(IND).GT.IENTMX(2))IENTMX(2)=NP2(IND)
            WRITE(MOUTIN,*)NP2(IND)
            DO J=1,NP2(IND)
              READ(INP,*)DUMMY1,DUMMY2,IDUMMY
              WRITE(MOUTIN,*)DUMMY1,DUMMY2,IDUMMY
            END DO
          END DO
        END DO
        GO TO 34
99      WRITE(IOUT,*)'IENTER = ',IENTER
        STOP 'BUG'
4       IND=0
        DO K=2,NMODE
          DO L=1,K-1
            IND=IND+1
            READ(MOUTIN,*)NP2(IND)
            DO J=1,NP2(IND)
              READ(MOUTIN,*)CP2(J,IND),IP2(J,IND,1),
     1        IP2(J,IND,2)
              IF(IP2(J,IND,1).GT.MAXPOW)MAXPOW=IP2(J,IND,1)
              IF(IP2(J,IND,2).GT.MAXPOW)MAXPOW=IP2(J,IND,2)
            END DO
          END DO
        END DO
        IF(MOLPRO.EQ.5.OR.MOLPRO.EQ.6.OR.(IWHICH.LT.0.AND.
     1  MOLINC.GT.0.AND.ICOUPL.EQ.2))IENTER=0
34      CONTINUE
      END IF
      IF(IENTER.LT.5)RETURN

      IF(MOLPRO.GT.6.OR.(MOLINC.GT.0.AND.IWHICH.LT.0.AND.ICOUPL.GT.2))
     1THEN
        GO TO(999,999,999,999,5,6,56,56,56,56)IENTER
C**THREE-DIMENSIONAL FITS
5       READ(INP,*)
        IENTMX(3)=0
        IND=0
        DO K=3,NMODE
          DO L=2,K-1
            DO N=1,L-1
              IND=IND+1
              READ(INP,*)NP3(IND)
              IF(NP3(IND).GT.IENTMX(3))IENTMX(3)=NP3(IND)
              WRITE(MOUTIN,*)NP3(IND)
              DO J=1,NP3(IND)
                READ(INP,*)DUMMY1,DUMMY2,DUMMY3,IDUMMY
                WRITE(MOUTIN,*)DUMMY1,DUMMY2,DUMMY3,IDUMMY
              END DO
            END DO
          END DO
        END DO
        GO TO 56
999     WRITE(IOUT,*)'IENTER = ',IENTER
        STOP 'BUG'
6       IND=0
        DO K=3,NMODE
          DO L=2,K-1
            DO N=1,L-1
              IND=IND+1
              READ(MOUTIN,*)NP3(IND)
              DO J=1,NP3(IND)
                READ(MOUTIN,*)CP3(J,IND),IP3(J,IND,1),IP3(J,IND,2),
     1          IP3(J,IND,3)
                IF(IP3(J,IND,1).GT.MAXPOW)MAXPOW=IP3(J,IND,1)
                IF(IP3(J,IND,2).GT.MAXPOW)MAXPOW=IP3(J,IND,2)
                IF(IP3(J,IND,3).GT.MAXPOW)MAXPOW=IP3(J,IND,3)
              END DO
            END DO
          END DO
        END DO
        IF(MOLPRO.EQ.7.OR.MOLPRO.EQ.8.OR.(IWHICH.LT.0.AND.
     1  MOLINC.GT.0.AND.ICOUPL.EQ.3))IENTER=0
56      CONTINUE
      END IF
      IF(IENTER.LT.7)RETURN

      IF(MOLPRO.GT.8.OR.(MOLINC.GT.0.AND.IWHICH.LT.0.AND.ICOUPL.GT.3))
     1THEN
        GO TO(9999,9999,9999,9999,9999,9999,7,8,78,78)IENTER
C**FOUR-DIMENSIONAL FITS
7       READ(INP,*)
        IENTMX(4)=0
        IND=0
        DO K=4,NMODE
          DO L=3,K-1
            DO N=2,L-1
              DO M=1,N-1
                IND=IND+1
                READ(INP,*)NP4(IND)
                IF(NP4(IND).GT.IENTMX(4))IENTMX(4)=NP4(IND)
                WRITE(MOUTIN,*)NP4(IND)
                DO J=1,NP4(IND)
                  READ(INP,*)DUMMY1,DUMMY2,DUMMY3,DUMMY4,IDUMMY
                  WRITE(MOUTIN,*)DUMMY1,DUMMY2,DUMMY3,DUMMY4,IDUMMY
                END DO
              END DO
            END DO
          END DO
        END DO
        GO TO 78
9999    WRITE(IOUT,*)'IENTER = ',IENTER
        STOP 'BUG'
8       IND=0
        DO K=4,NMODE
          DO L=3,K-1
            DO N=2,L-1
              DO M=1,N-1
                IND=IND+1
                READ(MOUTIN,*)NP4(IND)
                DO J=1,NP4(IND)
                  READ(MOUTIN,*)CP4(J,IND),IP4(J,IND,1),
     1            IP4(J,IND,2),IP4(J,IND,3),IP4(J,IND,4)
                  IF(IP4(J,IND,1).GT.MAXPOW)MAXPOW=IP4(J,IND,1)
                  IF(IP4(J,IND,2).GT.MAXPOW)MAXPOW=IP4(J,IND,2)
                  IF(IP4(J,IND,3).GT.MAXPOW)MAXPOW=IP4(J,IND,3)
                  IF(IP4(J,IND,4).GT.MAXPOW)MAXPOW=IP4(J,IND,4)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(MOLPRO.EQ.9.OR.MOLPRO.EQ.10.OR.(IWHICH.LT.0.AND.
     1  MOLINC.GT.0.AND.ICOUPL.EQ.4))IENTER=0
78      CONTINUE
      END IF
      IF(IENTER.LT.9)RETURN

      IF(MOLPRO.GT.10.OR.(MOLINC.GT.0.AND.IWHICH.LT.0.AND.ICOUPL.GT.4))
     1THEN
        GO TO(99999,99999,99999,99999,99999,99999,99999,99999,9,10)
     1  IENTER
C**FIVE-DIMENSIONAL FITS
9       READ(INP,*)
        IENTMX(5)=0
        IND=0
        DO K=5,NMODE
          DO L=4,K-1
            DO N=3,L-1
              DO M=2,N-1
                DO J=1,M-1
                  IND=IND+1
                  READ(INP,*)NP5(IND)
                  IF(NP5(IND).GT.IENTMX(5))IENTMX(5)=NP5(IND)
                  WRITE(MOUTIN,*)NP5(IND)
                  DO I=1,NP5(IND)
                    READ(INP,*)DUMMY1,DUMMY2,DUMMY3,DUMMY4,DUMMY5,
     1              IDUMMY
                    WRITE(MOUTIN,*)DUMMY1,DUMMY2,DUMMY3,DUMMY4,DUMMY5,
     1              IDUMMY
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        GO TO 90
99999   WRITE(IOUT,*)'IENTER = ',IENTER
        STOP 'BUG'
10      IND=0
        DO K=5,NMODE
          DO L=4,K-1
            DO N=3,L-1
              DO M=2,N-1
                DO J=1,M-1
                  IND=IND+1
                  READ(MOUTIN,*)NP5(IND)
                  DO I=1,NP5(IND)
                    READ(MOUTIN,*)CP5(I,IND),IP5(I,IND,1),
     1              IP5(I,IND,2),IP5(I,IND,3),
     2              IP5(I,IND,4),IP5(I,IND,5)
                    IF(IP5(I,IND,1).GT.MAXPOW)MAXPOW=IP5(I,IND,1)
                    IF(IP5(I,IND,2).GT.MAXPOW)MAXPOW=IP5(I,IND,2)
                    IF(IP5(I,IND,3).GT.MAXPOW)MAXPOW=IP5(I,IND,3)
                    IF(IP5(I,IND,4).GT.MAXPOW)MAXPOW=IP5(I,IND,4)
                    IF(IP5(I,IND,5).GT.MAXPOW)MAXPOW=IP5(I,IND,5)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(MOLPRO.EQ.11.OR.MOLPRO.EQ.12.OR.(IWHICH.LT.0.AND.
     1  MOLINC.GT.0.AND.ICOUPL.EQ.5))IENTER=0
90      CONTINUE
      END IF
      IF(IENTER.LT.11)RETURN
      IF(IENTER.NE.0)THEN
        WRITE(IOUT,*)'IENTER = ',IENTER
        STOP 'BUG'
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETPQ1(V,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),XTANPM(NMODE)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/AVCON/AVCON
      COMMON/FILASS/IOUT,INP
C**CORRECTION FOR CONSTANT + 1-DIM POTENTIALS
      V=AVCON
      DO K=1,NMODE
        IF(QQ(K).NE.0)THEN
          DO I=2,NP1(K)
            V=V+CP1(I,K)*XTANH(XTANPM(K)*QQ(K))**IP1(I,K)
          END DO
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETPQ2(V,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1NP2,CP2,IP2,NTOT2,MAX2,INDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),XTANPM(NMODE)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/AVCON/AVCON
      COMMON/FILASS/IOUT,INP
C**CORRECTION FOR CONSTANT + 1-DIM + 2-DIM POTENTIALS
      V=AVCON
      DO K=1,NMODE
        IF(QQ(K).NE.0)THEN
          DO I=2,NP1(K)
            V=V+CP1(I,K)*XTANH(XTANPM(K)*QQ(K))**IP1(I,K)
          END DO
          DO L=1,K-1
            IF(QQ(L).NE.0)THEN
              IND=INDK(K)+L
              DO I=1,NP2(IND)
                V=V+CP2(I,IND)*XTANH(XTANPM(K)*QQ(K))**IP2(I,IND,1)*
     1                         XTANH(XTANPM(L)*QQ(L))**IP2(I,IND,2)
              END DO
            END IF
          END DO
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETPQ3(V,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1NP2,CP2,IP2,NTOT2,MAX2,NP3,CP3,IP3,NTOT3,MAX3,INDK,INDL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),XTANPM(NMODE)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),IP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/AVCON/AVCON
      COMMON/FILASS/IOUT,INP
C**CORRECTION FOR CONSTANT + 1-DIM + 2-DIM + 3-DIM POTENTIALS
      V=AVCON
      DO K=1,NMODE
        IF(QQ(K).NE.0)THEN
          DO I=2,NP1(K)
            V=V+CP1(I,K)*XTANH(XTANPM(K)*QQ(K))**IP1(I,K)
          END DO
          DO L=1,K-1
            IF(QQ(L).NE.0)THEN
              IND=INDK(K)+L
              DO I=1,NP2(IND)
                V=V+CP2(I,IND)*XTANH(XTANPM(K)*QQ(K))**IP2(I,IND,1)*
     1                         XTANH(XTANPM(L)*QQ(L))**IP2(I,IND,2)
              END DO
              DO N=1,L-1
                IF(QQ(N).NE.0)THEN
                  IND=INDL(K)+INDK(L)+N
                  DO I=1,NP3(IND)
                    V=V+CP3(I,IND)*
     1              XTANH(XTANPM(K)*QQ(K))**IP3(I,IND,1)*
     2              XTANH(XTANPM(L)*QQ(L))**IP3(I,IND,2)*
     3              XTANH(XTANPM(N)*QQ(N))**IP3(I,IND,3)
                  END DO
                END IF
              END DO
            END IF
          END DO
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETPQ4(V,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1,
     1NP2,CP2,IP2,NTOT2,MAX2,NP3,CP3,IP3,NTOT3,MAX3,NP4,CP4,IP4,
     2NTOT4,MAX4,INDK,INDL,INDN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),XTANPM(NMODE)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),IP2(MAX2,NTOT2,2)
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),IP3(MAX3,NTOT3,3)
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),IP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/AVCON/AVCON
      COMMON/FILASS/IOUT,INP
C**CORRECTION FOR CONSTANT + 1-DIM + 2-DIM + 3-DIM + 4-DIM POTENTIALS
      V=AVCON
      DO K=1,NMODE
        IF(QQ(K).NE.0)THEN
          DO I=2,NP1(K)
            V=V+CP1(I,K)*XTANH(XTANPM(K)*QQ(K))**IP1(I,K)
          END DO
          DO L=1,K-1
            IF(QQ(L).NE.0)THEN
              IND=INDK(K)+L
              DO I=1,NP2(IND)
                V=V+CP2(I,IND)*XTANH(XTANPM(K)*QQ(K))**IP2(I,IND,1)*
     1                         XTANH(XTANPM(L)*QQ(L))**IP2(I,IND,2)
              END DO
              DO N=1,L-1
                IF(QQ(N).NE.0)THEN
                  IND=INDL(K)+INDK(L)+N
                  DO I=1,NP3(IND)
                    V=V+CP3(I,IND)*
     1              XTANH(XTANPM(K)*QQ(K))**IP3(I,IND,1)*
     2              XTANH(XTANPM(L)*QQ(L))**IP3(I,IND,2)*
     3              XTANH(XTANPM(N)*QQ(N))**IP3(I,IND,3)
                  END DO
                  DO M=1,N-1
                    IF(QQ(M).NE.0)THEN
                      IND=INDN(K)+INDL(L)+INDK(N)+M
                      DO I=1,NP4(IND)
                        V=V+CP4(I,IND)*
     1                  XTANH(XTANPM(K)*QQ(K))**IP4(I,IND,1)*
     2                  XTANH(XTANPM(L)*QQ(L))**IP4(I,IND,2)*
     3                  XTANH(XTANPM(N)*QQ(N))**IP4(I,IND,3)*
     4                  XTANH(XTANPM(M)*QQ(M))**IP4(I,IND,4)
                      END DO
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END DO
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETPQT(V,NMODE,QQ,XTANPM,NP1,CP1,IP1,NTOT1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),XTANPM(NMODE)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),IP1(MMX1,NTOT1)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/AVCON/AVCON
      COMMON/FILASS/IOUT,INP
C**NORMAL COORDINATE POTENTIAL
      V=0
      DO K=1,NMODE
        IF(QQ(K).NE.0)THEN
          DO I=1,NP1(K)
            V=V+CP1(I,K)*XTANH(XTANPM(K)*QQ(K))**IP1(I,K)
          END DO
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDQT(V,NMODE,QQ,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE)
C**NOT YET IMPLEMENTED
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE HERMIN(NMODE,XTANPM,NP1,GP1,V1,D1,NTOT1,MMX1,NP2,GP2,
     1V2,D2A,D2B,NTOT2,MAX2,NP3,GP3,V3,D3A,D3B,D3C,NTOT3,MAX3,NP4,GP4,
     2V4,D4A,D4B,D4C,D4D,NTOT4,MAX4,INDK,INDL,INDN,INDM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 SYMBOL(100),CHSYM(4)
      DIMENSION XTANPM(NMODE),INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION NP1(NTOT1),GP1(MMX1,NTOT1),V1(MMX1,NTOT1),
     1D1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),GP2(MAX2,NMODE),
     1V2(MAX2,MAX2,NTOT2),D2A(MAX2,MAX2,NTOT2),D2B(MAX2,MAX2,NTOT2)
      DIMENSION NP3(NTOT3),GP3(MAX3,NMODE),
     1V3(MAX3,MAX3,MAX3,NTOT3),D3A(MAX3,MAX3,MAX3,NTOT3),
     2D3B(MAX3,MAX3,MAX3,NTOT3),D3C(MAX3,MAX3,MAX3,NTOT3)
      DIMENSION NP4(NTOT4),GP4(MAX4,NMODE),
     1V4(MAX4,MAX4,MAX4,MAX4,NTOT4),D4A(MAX4,MAX4,MAX4,MAX4,NTOT4),
     2D4B(MAX4,MAX4,MAX4,MAX4,NTOT4),D4C(MAX4,MAX4,MAX4,MAX4,NTOT4),
     3D4D(MAX4,MAX4,MAX4,MAX4,NTOT4)
      DIMENSION DK(100),QQ(100)
      COMMON/FILASS/IOUT,INP,MOUTIN
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/ENTER/IENTER,IENTMX(5)
106   FORMAT(2I5)
107   FORMAT(5F15.8)

      IF(MOLPRO.LT.-2.OR.(MOLINC.LT.0.AND.IWHICH.LT.0.AND.ICOUPL.GT.0))
     1THEN
        REWIND MOUTIN
        IF(IENTER.GE.2)THEN
          IENTER=IENTER+2
        ELSE
          IENTER=IENTER+1
        END IF
        GO TO(1,2,12,12,12,12,12,12,12,12)IENTER
C**ONE-DIMENSIONAL INTERPOLATIONS
1       CONTINUE
        READ(INP,*)
        IENTMX(1)=0
        DO K=1,NMODE
          READ(INP,*)NP1(K)
          NPTSK=NP1(K)
          IF(NPTSK.GT.IENTMX(1))IENTMX(1)=NPTSK
          WRITE(MOUTIN,*)NPTSK
          DO IDUM=1,3
            READ(INP,*)(DK(I),I=1,NPTSK)
            WRITE(MOUTIN,*)(DK(I),I=1,NPTSK)
          END DO
        END DO
        IENTMX(2)=IENTMX(1)
        IENTMX(3)=IENTMX(1)
        IENTMX(4)=IENTMX(1)
        IENTMX(5)=IENTMX(1)
        GO TO 12
2       DO K=1,NMODE
          READ(MOUTIN,*)NP1(K)
          NPTSK=NP1(K)
          READ(MOUTIN,*)(GP1(M,K),M=1,NPTSK)
          GP1(1,K)=GP1(1,K)-1.D-8
          GP1(NPTSK,K)=GP1(NPTSK,K)+1.D-8
C**V1(M,K) ARE ENERGIES AT 'K' POINTS
          READ(MOUTIN,*)(V1(M,K),M=1,NPTSK)
C**D1(M,K) ARE DERIVATIVES WRT K AT 'K' POINTS
          READ(MOUTIN,*)(D1(M,K),M=1,NPTSK)
        END DO
        IF(MOLPRO.EQ.-3.OR.MOLPRO.EQ.-4.OR.(IWHICH.LT.0.AND.
     1  MOLINC.LT.0.AND.ICOUPL.EQ.1))IENTER=0
12      CONTINUE
      END IF
      IF(MOLINC.LT.0.AND.IENTER.EQ.0)GO TO 100
      IF(IENTER.LT.3)RETURN

      IF(MOLPRO.LT.-4.OR.(MOLINC.LT.0.AND.IWHICH.LT.0.AND.ICOUPL.GT.1))
     1THEN
        GO TO(99,99,99,4,34,34,34,34,34,34)IENTER
C**TWO-DIMENSIONAL INTERPOLATIONS
99      WRITE(IOUT,*)'IENTER = ',IENTER
        STOP 'BUG'
4       IND=0
        READ(INP,*)
        DO K=2,NMODE
          DO L=1,K-1
            IND=IND+1
            READ(INP,*)NP2(K),NP2(L)
            NPTSK=NP2(K)
            READ(INP,*)(GP2(M,K),M=1,NPTSK)
            NPTSL=NP2(L)
            READ(INP,*)(GP2(M,L),M=1,NPTSL)
C**V2(ML,MK,IND) ARE ENERGIES AT 'L' POINTS FOR EACH 'K' POINT
            DO MK=1,NPTSK
              READ(INP,*)(V2(ML,MK,IND),ML=1,NPTSL)
            END DO
C**D2(MK,ML,IND) ARE DERIVATIVES WRT K AT 'K' POINTS FOR EACH
C**'L' POINT (K > L)
            DO ML=1,NPTSL
              READ(INP,*)(D2A(MK,ML,IND),MK=1,NPTSK)
            END DO
C**D2(ML,MK,IND) ARE DERIVATIVES WRT L AT 'L' POINTS FOR EACH
C**'K' POINT (L < K)
            DO MK=1,NPTSK
              READ(INP,*)(D2B(ML,MK,IND),ML=1,NPTSL)
            END DO
          END DO
        END DO
        IF(MOLPRO.EQ.-5.OR.MOLPRO.EQ.-6.OR.(IWHICH.LT.0.AND.
     1  MOLINC.LT.0.AND.ICOUPL.EQ.2))IENTER=0
34      CONTINUE
      END IF
      IF(MOLINC.LT.0.AND.IENTER.EQ.0)GO TO 100
      IF(IENTER.LT.5)RETURN

      IF(MOLPRO.LT.-6.OR.(MOLINC.LT.0.AND.IWHICH.LT.0.AND.ICOUPL.GT.2))
     1THEN
        GO TO(999,999,999,999,999,6,56,56,56,56)IENTER
C**THREE-DIMENSIONAL INTERPOLATIONS
999     WRITE(IOUT,*)'IENTER = ',IENTER
        STOP 'BUG'
6       IND=0
        READ(INP,*)
        DO K=3,NMODE
          DO L=2,K-1
            DO N=1,L-1
              IND=IND+1
              READ(INP,*)NP3(K),NP3(L),NP3(N)
              NPTSK=NP3(K)
              READ(INP,*)(GP3(M,K),M=1,NPTSK)
              NPTSL=NP3(L)
              READ(INP,*)(GP3(M,L),M=1,NPTSL)
              NPTSN=NP3(N)
              READ(INP,*)(GP3(M,N),M=1,NPTSN)
C**V3(MN,ML,MK,IND) ARE ENERGIES AT 'N' POINTS FOR EACH 'L','K' POINT
              DO MK=1,NPTSK
                DO ML=1,NPTSL
                  READ(INP,*)(V3(MN,ML,MK,IND),MN=1,NPTSN)
                END DO
              END DO
C**D3(MK,MN,ML,IND) ARE DERIVS. WRT K AT 'K' POINTS FOR EACH
C**'L','N' POINT (K > L > N)
              DO ML=1,NPTSL
                DO MN=1,NPTSN
                  READ(INP,*)(D3A(MK,MN,ML,IND),
     1            MK=1,NPTSK)
                END DO
              END DO
C**D3(ML,MN,MK,IND) ARE DERIVS. WRT L AT 'L' POINTS FOR EACH
C**'N','K' POINT (K > L > N)
              DO MK=1,NPTSK
                DO MN=1,NPTSN
                  READ(INP,*)(D3B(ML,MN,MK,IND),
     1            ML=1,NPTSL)
                END DO
              END DO
C**D3(MN,ML,MK,IND) ARE DERIVS. WRT N AT 'N' POINTS FOR EACH
C**'L','K' POINT (K > L > N)
              DO MK=1,NPTSK
                DO ML=1,NPTSL
                  READ(INP,*)(D3C(MN,ML,MK,IND),
     1            MN=1,NPTSN)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(MOLPRO.EQ.-7.OR.MOLPRO.EQ.-8.OR.(IWHICH.LT.0.AND.
     1  MOLINC.LT.0.AND.ICOUPL.EQ.3))IENTER=0
56      CONTINUE
      END IF
      IF(MOLINC.LT.0.AND.IENTER.EQ.0)GO TO 100
      IF(IENTER.LT.7)RETURN

      IF(MOLPRO.LT.-8.OR.(MOLINC.LT.0.AND.IWHICH.LT.0.AND.ICOUPL.GT.3))
     1THEN
        GO TO(9999,9999,9999,9999,9999,9999,9999,8,78,78)IENTER
C**FOUR-DIMENSIONAL INTERPOLATIONS
9999    WRITE(IOUT,*)'IENTER = ',IENTER
        STOP 'BUG'
8       IND=0
        READ(INP,*)
        DO K=4,NMODE
          DO L=3,K-1
            DO N=2,L-1
              DO M=1,N-1
                IND=IND+1
                READ(INP,*)NP4(K),NP4(L),NP4(N),NP4(M)
                NPTSK=NP4(K)
                READ(INP,*)(GP4(I,K),I=1,NPTSK)
                NPTSL=NP4(L)
                READ(INP,*)(GP4(I,L),I=1,NPTSL)
                NPTSN=NP4(N)
                READ(INP,*)(GP4(I,N),I=1,NPTSN)
                NPTSM=NP4(M)
                READ(INP,*)(GP4(I,M),I=1,NPTSM)
C**V4(MM,MN,ML,MK,IND) ARE ENERGIES AT 'M' POINTS FOR EACH
C**'N','L','K' POINT
                DO MK=1,NPTSK
                  DO ML=1,NPTSL
                    DO MN=1,NPTSN
                      READ(INP,*)(V4(MM,MN,ML,MK,IND),MM=1,NPTSM)
                    END DO
                  END DO
                END DO
C**D4(MK,MM,MN,ML,IND) ARE DERIVS. WRT K AT 'K' POINTS FOR EACH
C**'M','N','L' POINT
                DO ML=1,NPTSL
                  DO MN=1,NPTSN
                    DO MM=1,NPTSM
                      READ(INP,*)(D4A(MM,MN,ML,MK,IND),MK=1,NPTSK)
                    END DO
                  END DO
                END DO
C**D4(ML,MM,MN,MK,IND) ARE DERIVS. WRT L AT 'L' POINTS FOR EACH
C**'M','N','K' POINT
                DO MK=1,NPTSK
                  DO MN=1,NPTSN
                    DO MM=1,NPTSM
                      READ(INP,*)(D4B(MM,MN,ML,MK,IND),ML=1,NPTSL)
                    END DO
                  END DO
                END DO
C**D4(MN,MM,ML,MK,N,L,K) ARE DERIVS. WRT N AT 'N' POINTS FOR EACH
C**'M','L','K' POINT
                DO MK=1,NPTSK
                  DO ML=1,NPTSL
                    DO MM=1,NPTSM
                      READ(INP,*)(D4C(MM,MN,ML,MK,IND),MN=1,NPTSN)
                    END DO
                  END DO
                END DO
C**D4(MM,MN,ML,MK,N,L,K) ARE DERIVS. WRT M AT 'M' POINTS FOR EACH
C**'N','L','K' POINT
                DO MK=1,NPTSK
                  DO ML=1,NPTSL
                    DO MN=1,NPTSN
                      READ(INP,*)(D4D(MM,MN,ML,MK,IND),MM=1,NPTSM)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(MOLPRO.EQ.-9.OR.MOLPRO.EQ.-10.OR.(IWHICH.LT.0.AND.
     1  MOLINC.LT.0.AND.ICOUPL.EQ.4))IENTER=0
78      CONTINUE
      END IF
      IF(MOLINC.LT.0.AND.IENTER.EQ.0)GO TO 100
      IF(IENTER.LT.9)RETURN
C
100   CONTINUE

      IF(ICOUPL.EQ.1)GO TO 101
      IF(ICOUPL.EQ.2)GO TO 102
      IF(ICOUPL.EQ.3)GO TO 103
      IF(ICOUPL.EQ.4)GO TO 104
      STOP 'BUG IN HERMIN'

C**MODIFY 4-D DATA
104   CONTINUE
      WRITE(IOUT,*)'MODIFY 4-D'
C**CORRECTION FOR 1-, 2- AND 3-DIM POTENTIALS
      DO K=4,NMODE
        NPTSK=NP4(K)
        DO MK=1,NPTSK
          QK=GP4(MK,K)
          DO L=3,K-1
            NPTSL=NP4(L)
            DO ML=1,NPTSL
              QL=GP4(ML,L)
              DO N=2,L-1
                NPTSN=NP4(N)
                DO MN=1,NPTSN
                  QN=GP4(MN,N)
                  DO M=1,N-1
                    NPTSM=NP4(M)
                    DO MM=1,NPTSM
                      QM=GP4(MM,M)
                      DO I=1,NMODE
                        QQ(I)=0
                      END DO
                      QQ(K)=QK
                      QQ(L)=QL
                      QQ(N)=QN
                      QQ(M)=QM
                      IND=INDN(K)+INDL(L)+INDK(N)+M
                      CALL GETQP1(V,DK,NMODE,QQ,NP1,GP1,V1,D1,NTOT1,
     1                MMX1)
                      V4(MM,MN,ML,MK,IND)=V4(MM,MN,ML,MK,IND)-V
                      D4A(MM,MN,ML,MK,IND)=D4A(MM,MN,ML,MK,IND)-DK(K)
                      D4B(MM,MN,ML,MK,IND)=D4B(MM,MN,ML,MK,IND)-DK(L)
                      D4C(MM,MN,ML,MK,IND)=D4C(MM,MN,ML,MK,IND)-DK(N)
                      D4D(MM,MN,ML,MK,IND)=D4D(MM,MN,ML,MK,IND)-DK(M)
                      CALL GETQP2(V,DK,NMODE,QQ,NP1,GP1,V1,D1,NTOT1,
     1                MMX1,NP2,GP2,V2,D2A,D2B,NTOT2,MAX2,INDK)
                      V4(MM,MN,ML,MK,IND)=V4(MM,MN,ML,MK,IND)+V
                      D4A(MM,MN,ML,MK,IND)=D4A(MM,MN,ML,MK,IND)+DK(K)
                      D4B(MM,MN,ML,MK,IND)=D4B(MM,MN,ML,MK,IND)+DK(L)
                      D4C(MM,MN,ML,MK,IND)=D4C(MM,MN,ML,MK,IND)+DK(N)
                      D4D(MM,MN,ML,MK,IND)=D4D(MM,MN,ML,MK,IND)+DK(M)
                      CALL GETQP3(V,DK,NMODE,QQ,NP1,GP1,V1,D1,NTOT1,
     1                MMX1,NP2,GP2,V2,D2A,D2B,NTOT2,MAX2,NP3,GP3,V3,
     2                D3A,D3B,D3C,NTOT3,MAX3,INDK,INDL)
                      V4(MM,MN,ML,MK,IND)=V4(MM,MN,ML,MK,IND)-V
                      D4A(MM,MN,ML,MK,IND)=D4A(MM,MN,ML,MK,IND)-DK(K)
                      D4B(MM,MN,ML,MK,IND)=D4B(MM,MN,ML,MK,IND)-DK(L)
                      D4C(MM,MN,ML,MK,IND)=D4C(MM,MN,ML,MK,IND)-DK(N)
                      D4D(MM,MN,ML,MK,IND)=D4D(MM,MN,ML,MK,IND)-DK(M)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

C**MODIFY 3-D DATA
103   CONTINUE
      WRITE(IOUT,*)'MODIFY 3-D'
C**CORRECTION FOR 1- AND 2-DIM POTENTIALS
      DO K=3,NMODE
        NPTSK=NP3(K)
        DO MK=1,NPTSK
          QK=GP3(MK,K)
          DO L=2,K-1
            NPTSL=NP3(L)
            DO ML=1,NPTSL
              QL=GP3(ML,L)
              DO N=1,L-1
                NPTSN=NP3(N)
                DO MN=1,NPTSN
                  QN=GP3(MN,N)
                  DO I=1,NMODE
                    QQ(I)=0
                  END DO
                  QQ(K)=QK
                  QQ(L)=QL
                  QQ(N)=QN
                  IND=INDL(K)+INDK(L)+N
                  CALL GETQP1(V,DK,NMODE,QQ,NP1,GP1,V1,D1,NTOT1,MMX1)
                  V3(MN,ML,MK,IND)=V3(MN,ML,MK,IND)+V
                  D3A(MK,MN,ML,IND)=D3A(MK,MN,ML,IND)+DK(K)
                  D3B(ML,MN,MK,IND)=D3B(ML,MN,MK,IND)+DK(L)
                  D3C(MN,ML,MK,IND)=D3C(MN,ML,MK,IND)+DK(N)
                  CALL GETQP2(V,DK,NMODE,QQ,NP1,GP1,V1,D1,NTOT1,MMX1,
     1            NP2,GP2,V2,D2A,D2B,NTOT2,MAX2,INDK)
                  V3(MN,ML,MK,IND)=V3(MN,ML,MK,IND)-V
                  D3A(MK,MN,ML,IND)=D3A(MK,MN,ML,IND)-DK(K)
                  D3B(ML,MN,MK,IND)=D3B(ML,MN,MK,IND)-DK(L)
                  D3C(MN,ML,MK,IND)=D3C(MN,ML,MK,IND)-DK(N)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

C**MODIFY 2-D DATA
102   CONTINUE
      WRITE(IOUT,*)'MODIFY 2-D'
C**CORRECTION FOR 1-DIM POTENTIALS
      DO K=2,NMODE
        NPTSK=NP2(K)
        DO MK=1,NPTSK
          QK=GP2(MK,K)
          DO L=1,K-1
            NPTSL=NP2(L)
            DO ML=1,NPTSL
              QL=GP2(ML,L)
              DO I=1,NMODE
                QQ(I)=0
              END DO
              QQ(K)=QK
              QQ(L)=QL
              IND=INDK(K)+L
              CALL GETQP1(V,DK,NMODE,QQ,NP1,GP1,V1,D1,NTOT1,MMX1)
              V2(ML,MK,IND)=V2(ML,MK,IND)-V
              D2A(MK,ML,IND)=D2A(MK,ML,IND)-DK(K)
              D2B(ML,MK,IND)=D2B(ML,MK,IND)-DK(L)
            END DO
          END DO
        END DO
      END DO

101   CONTINUE
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETQP1(V,DK,NMODE,QQ,NP1,GP1,VP1,DP1,NTOT1,
     1MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),DK(NMODE)
      DIMENSION NP1(NTOT1),GP1(MMX1,NTOT1),VP1(MMX1,NTOT1),
     1DP1(MMX1,NTOT1)
C**CORRECTION FOR 1-DIM POTENTIALS
      V=0
C**   V=AVCON
      DO K=1,NMODE
        IF(QQ(K).NE.0)THEN
          QK=QQ(K)
          CALL GETHP1(VK,K,QK,NP1,GP1,VP1,DP1,NTOT1,MMX1)
          V=V+VK
C***************
CC    GO TO 1000
C***************
          QKP=QK+1.D-6
          CALL GETHP1(VKP,K,QKP,NP1,GP1,VP1,DP1,NTOT1,MMX1)
          QKM=QK-1.D-6
          CALL GETHP1(VKM,K,QKM,NP1,GP1,VP1,DP1,NTOT1,MMX1)
          DK(K)=(VKP-VKM)*1.D6/2
1000  CONTINUE
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETHP1(V,MODE,QQ,NP1,GP1,V1,D1,NTOT1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NP1(NTOT1),GP1(MMX1,NTOT1),V1(MMX1,NTOT1),D1(MMX1,NTOT1)
C*************************************************
      dimension a(0:1),ax(0:1)
C*************************************************
      COMMON/FILASS/IOUT,INP
      V=0
C**MODE K WILL BE X
      K=MODE
      X=QQ
C**FIND X ORIGIN (K)
      DO M=1,NP1(K)
        IF(X.LE.GP1(M,K))THEN
          MK=M
          GO TO 11
        END IF
      END DO
11    CONTINUE
c  get x relative to x(k-1)
      XIR=X-GP1(MK-1,K)
C**GET INTERVAL IN X
      DX=GP1(MK,K)-GP1(MK-1,K)
      A(0)=V1(MK-1,K)
      A(1)=V1(MK,K)
      AX(0)=D1(MK-1,K)
      AX(1)=D1(MK,K)
      V=V+hermi1(dx,a,ax,xir)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETQP2(V,DK,NMODE,QQ,NP1,GP1,VP1,DP1,NTOT1,MMX1,
     1NP2,GP2,VP2,DP2A,DP2B,NTOT2,MAX2,INDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),DK(NMODE),INDK(1)
      DIMENSION NP1(NTOT1),GP1(MMX1,NTOT1),VP1(MMX1,NTOT1),
     1DP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),GP2(MAX2,NMODE),
     1VP2(MAX2,MAX2,NTOT2),DP2A(MAX2,MAX2,NTOT2),DP2B(MAX2,MAX2,NTOT2)
      V=0
      DO I=1,NMODE
        DK(I)=0
      END DO
C**CORRECTION FOR 1-DIM POTENTIALS
CC    CALL GETQP1(V,DK,NMODE,QQ,NP1,GP1,VP1,DP1,NTOT1,MMX1)
C**CORRECTION FOR 2-DIM POTENTIALS
      DO K=1,NMODE
        IF(QQ(K).NE.0)THEN
          QK=QQ(K)
          QKP=QK+1.D-6
          QKM=QK-1.D-6
          DO L=1,K-1
            IF(QQ(L).NE.0)THEN
              QL=QQ(L)
              QLP=QL+1.D-6
              QLM=QL-1.D-6
              CALL GETHP2(VKL,K,L,QK,QL,NP2,GP2,VP2,DP2A,DP2B,NTOT2,
     1        MAX2,INDK)
              V=V+VKL
C***************
CC    GO TO 1000
C***************
              CALL GETHP2(VKP,K,L,QKP,QL,NP2,GP2,VP2,DP2A,DP2B,NTOT2,
     1        MAX2,INDK)
              CALL GETHP2(VKM,K,L,QKM,QL,NP2,GP2,VP2,DP2A,DP2B,NTOT2,
     1        MAX2,INDK)
              DK(K)=DK(K)+(VKP-VKM)*1.D6/2
              CALL GETHP2(VLP,K,L,QK,QLP,NP2,GP2,VP2,DP2A,DP2B,NTOT2,
     1        MAX2,INDK)
              CALL GETHP2(VLM,K,L,QK,QLM,NP2,GP2,VP2,DP2A,DP2B,NTOT2,
     1        MAX2,INDK)
              DK(L)=DK(L)+(VLP-VLM)*1.D6/2
1000  CONTINUE
            END IF
          END DO
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETHP2(V,MODEK,MODEL,QK,QL,NP2,HP2,V2,D2A,D2B,NTOT2,
     1MAX2,INDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION INDK(1)
      DIMENSION NP2(NTOT2),HP2(MAX2,1),
     1V2(MAX2,MAX2,NTOT2),D2A(MAX2,MAX2,NTOT2),D2B(MAX2,MAX2,NTOT2)
C*************************************************
      dimension f(0:1,0:1),fx(0:1,0:1),fy(0:1,0:1)
C*************************************************
      COMMON/FILASS/IOUT,INP
      V=0
C**MODE K WILL BE X
      K=MODEK
      X=QK
C**MODE L WILL BE Y
      L=MODEL
      Y=QL
C**FIND X ORIGIN (K)
      DO M=1,NP2(K)
        IF(X.LE.HP2(M,K))THEN
          MK=M
          GO TO 21
        END IF
      END DO
21    CONTINUE
C**FIND Y ORIGIN (L)
      DO M=1,NP2(L)
        IF(Y.LE.HP2(M,L))THEN
          ML=M
          GO TO 22
        END IF
      END DO
22    CONTINUE
c  get x relative to x(k-1)
      XIR=X-HP2(MK-1,K)
c  get y relative to y(l-1)
      YIR=Y-HP2(ML-1,L)
C**GET INTERVAL IN X
      DX=HP2(MK,K)-HP2(MK-1,K)
C**GET INTERVAL IN Y
      DY=HP2(ML,L)-HP2(ML-1,L)
C
      IND=INDK(K)+L
      F(0,0)=V2(ML-1,MK-1,IND)
      F(1,0)=V2(ML-1,MK,IND)
      F(0,1)=V2(ML,MK-1,IND)
      F(1,1)=V2(ML,MK,IND)
C
      FX(0,0)=D2A(MK-1,ML-1,IND)
      FX(1,0)=D2A(MK,ML-1,IND)
      FX(0,1)=D2A(MK-1,ML,IND)
      FX(1,1)=D2A(MK,ML,IND)
C
      FY(0,0)=D2B(ML-1,MK-1,IND)
      FY(1,0)=D2B(ML-1,MK,IND)
      FY(0,1)=D2B(ML,MK-1,IND)
      FY(1,1)=D2B(ML,MK,IND)
      V=V+hermi2(dx,dy,f,fx,fy,xir,yir)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETQP3(V,DK,NMODE,QQ,NP1,GP1,VP1,DP1,NTOT1,
     1MMX1,NP2,GP2,VP2,DP2A,DP2B,NTOT2,MAX2,NP3,GP3,VP3,DP3A,DP3B,
     2DP3C,NTOT3,MAX3,INDK,INDL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),INDK(1),INDL(1),DK(NMODE)
      DIMENSION NP1(NTOT1),GP1(MMX1,NTOT1),VP1(MMX1,NTOT1),
     1DP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),GP2(MAX2,NMODE),
     1VP2(MAX2,MAX2,NTOT2),DP2A(MAX2,MAX2,NTOT2),DP2B(MAX2,MAX2,NTOT2)
      DIMENSION NP3(NTOT3),GP3(MAX3,NMODE),
     1VP3(MAX3,MAX3,MAX3,NTOT3),DP3A(MAX3,MAX3,MAX3,NTOT3),
     2DP3B(MAX3,MAX3,MAX3,NTOT3),DP3C(MAX3,MAX3,MAX3,NTOT3)
      V=0
      DO I=1,NMODE
        DK(I)=0
      END DO
C**CORRECTION FOR 1-DIM AND 2-DIM POTENTIALS
CC    CALL GETQP2(V,DK,NMODE,QQ,NP1,GP1,VP1,DP1,NTOT1,MMX1,
CC   1NP2,GP2,VP2,DP2A,DP2B,NTOT2,MAX2,INDK)
C**CORRECTION FOR 3-DIM POTENTIALS
      DO K=1,NMODE
        IF(QQ(K).NE.0)THEN
          QK=QQ(K)
          QKP=QK+1.D-6
          QKM=QK-1.D-6
          DO L=1,K-1
            IF(QQ(L).NE.0)THEN
              QL=QQ(L)
              QLP=QL+1.D-6
              QLM=QL-1.D-6
              DO N=1,L-1
                IF(QQ(N).NE.0)THEN
                  QN=QQ(N)
                  QNP=QN+1.D-6
                  QNM=QN-1.D-6
                  CALL GETHP3(VKLN,K,L,N,QK,QL,QN,NP3,GP3,VP3,DP3A,
     1            DP3B,DP3C,NTOT3,MAX3,INDK,INDL)
                  V=V+VKLN
C***************
CC    GO TO 1000
C***************
                  CALL GETHP3(VKP,K,L,N,QKP,QL,QN,NP3,GP3,VP3,DP3A,
     1            DP3B,DP3C,NTOT3,MAX3,INDK,INDL)
                  CALL GETHP3(VKM,K,L,N,QKM,QL,QN,NP3,GP3,VP3,DP3A,
     1            DP3B,DP3C,NTOT3,MAX3,INDK,INDL)
                  DK(K)=DK(K)+(VKP-VKM)*1.D6/2
                  CALL GETHP3(VLP,K,L,N,QK,QLP,QN,NP3,GP3,VP3,DP3A,
     1            DP3B,DP3C,NTOT3,MAX3,INDK,INDL)
                  CALL GETHP3(VLM,K,L,N,QK,QLM,QN,NP3,GP3,VP3,DP3A,
     1            DP3B,DP3C,NTOT3,MAX3,INDK,INDL)
                  DK(L)=DK(L)+(VLP-VLM)*1.D6/2
                  CALL GETHP3(VNP,K,L,N,QK,QL,QNP,NP3,GP3,VP3,DP3A,
     1            DP3B,DP3C,NTOT3,MAX3,INDK,INDL)
                  CALL GETHP3(VNM,K,L,N,QK,QL,QNM,NP3,GP3,VP3,DP3A,
     1            DP3B,DP3C,NTOT3,MAX3,INDK,INDL)
                  DK(N)=DK(N)+(VNP-VNM)*1.D6/2
1000  CONTINUE
                END IF
              END DO
            END IF
          END DO
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETHP3(V,MODEK,MODEL,MODEN,QK,QL,QN,NP3,HP3,V3,D3A,
     1D3B,D3C,NTOT3,MAX3,INDK,INDL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION INDK(1),INDL(1)
      DIMENSION NP3(NTOT3),HP3(MAX3,1),
     1V3(MAX3,MAX3,MAX3,NTOT3),D3A(MAX3,MAX3,MAX3,NTOT3),
     2D3B(MAX3,MAX3,MAX3,NTOT3),D3C(MAX3,MAX3,MAX3,NTOT3)
C*************************************************
      dimension g(0:1,0:1,0:1),gx(0:1,0:1,0:1),gy(0:1,0:1,0:1)
      dimension gz(0:1,0:1,0:1)
C*************************************************
      COMMON/FILASS/IOUT,INP
      V=0
C**MODE K WILL BE X
      K=MODEK
      X=QK
C**MODE L WILL BE Y
      L=MODEL
      Y=QL
C**MODE N WILL BE Z
      N=MODEN
      Z=QN
C**FIND X ORIGIN (K)
      DO M=1,NP3(K)
        IF(X.LE.HP3(M,K))THEN
          MK=M
          GO TO 31
        END IF
      END DO
31    CONTINUE
C**FIND Y ORIGIN (L)
      DO M=1,NP3(L)
        IF(Y.LE.HP3(M,L))THEN
          ML=M
          GO TO 32
        END IF
      END DO
32    CONTINUE
C**FIND Z ORIGIN (N)
      DO M=1,NP3(N)
        IF(Z.LE.HP3(M,N))THEN
          MN=M
          GO TO 33
        END IF
      END DO
33    CONTINUE
c  get x relative to x(k-1)
      XIR=X-HP3(MK-1,K)
c  get y relative to y(l-1)
      YIR=Y-HP3(ML-1,L)
c  get z relative to z(n-1)
      ZIR=Z-HP3(MN-1,N)
C**GET INTERVAL IN X
      DX=HP3(MK,K)-HP3(MK-1,K)
C**GET INTERVAL IN Y
      DY=HP3(ML,L)-HP3(ML-1,L)
C**GET INTERVAL IN Z
      DZ=HP3(MN,N)-HP3(MN-1,N)
C
      IND=INDL(K)+INDK(L)+N
      G(0,0,0)=V3(MN-1,ML-1,MK-1,IND)
      G(1,0,0)=V3(MN-1,ML-1,MK,IND)
      G(0,1,0)=V3(MN-1,ML,MK-1,IND)
      G(0,0,1)=V3(MN,ML-1,MK-1,IND)
      G(1,1,0)=V3(MN-1,ML,MK,IND)
      G(1,0,1)=V3(MN,ML-1,MK,IND)
      G(0,1,1)=V3(MN,ML,MK-1,IND)
      G(1,1,1)=V3(MN,ML,MK,IND)
C
      GX(0,0,0)=D3A(MK-1,MN-1,ML-1,IND)
      GX(1,0,0)=D3A(MK,MN-1,ML-1,IND)
      GX(0,1,0)=D3A(MK-1,MN-1,ML,IND)
      GX(0,0,1)=D3A(MK-1,MN,ML-1,IND)
      GX(1,1,0)=D3A(MK,MN-1,ML,IND)
      GX(1,0,1)=D3A(MK,MN,ML-1,IND)
      GX(0,1,1)=D3A(MK-1,MN,ML,IND)
      GX(1,1,1)=D3A(MK,MN,ML,IND)
C
      GY(0,0,0)=D3B(ML-1,MN-1,MK-1,IND)
      GY(1,0,0)=D3B(ML-1,MN-1,MK,IND)
      GY(0,1,0)=D3B(ML,MN-1,MK-1,IND)
      GY(0,0,1)=D3B(ML-1,MN,MK-1,IND)
      GY(1,1,0)=D3B(ML,MN-1,MK,IND)
      GY(1,0,1)=D3B(ML-1,MN,MK,IND)
      GY(0,1,1)=D3B(ML,MN,MK-1,IND)
      GY(1,1,1)=D3B(ML,MN,MK,IND)
C
      GZ(0,0,0)=D3C(MN-1,ML-1,MK-1,IND)
      GZ(1,0,0)=D3C(MN-1,ML-1,MK,IND)
      GZ(0,1,0)=D3C(MN-1,ML,MK-1,IND)
      GZ(0,0,1)=D3C(MN,ML-1,MK-1,IND)
      GZ(1,1,0)=D3C(MN-1,ML,MK,IND)
      GZ(1,0,1)=D3C(MN,ML-1,MK,IND)
      GZ(0,1,1)=D3C(MN,ML,MK-1,IND)
      GZ(1,1,1)=D3C(MN,ML,MK,IND)
      V=V+hermi3(dx,dy,dz,g,gx,gy,gz,xir,yir,zir)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETQP4(V,DK,NMODE,QQ,NP1,GP1,VP1,DP1,NTOT1,
     1MMX1,NP2,GP2,VP2,DP2A,DP2B,NTOT2,MAX2,NP3,GP3,VP3,DP3A,DP3B,
     2DP3C,NTOT3,MAX3,NP4,GP4,VP4,DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,
     3INDL,INDN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(NMODE),INDK(1),INDL(1),INDN(1),DK(NMODE)
      DIMENSION NP1(NTOT1),GP1(MMX1,NTOT1),VP1(MMX1,NTOT1),
     1DP1(MMX1,NTOT1)
      DIMENSION NP2(NTOT2),GP2(MAX2,NMODE),
     1VP2(MAX2,MAX2,NTOT2),DP2A(MAX2,MAX2,NTOT2),DP2B(MAX2,MAX2,NTOT2)
      DIMENSION NP3(NTOT3),GP3(MAX3,NMODE),
     1VP3(MAX3,MAX3,MAX3,NTOT3),DP3A(MAX3,MAX3,MAX3,NTOT3),
     2DP3B(MAX3,MAX3,MAX3,NTOT3),DP3C(MAX3,MAX3,MAX3,NTOT3)
      DIMENSION NP4(NTOT4),GP4(MAX4,NMODE),
     1VP4(MAX4,MAX4,MAX4,MAX4,NTOT4),DP4A(MAX4,MAX4,MAX4,MAX4,NTOT4),
     2DP4B(MAX4,MAX4,MAX4,MAX4,NTOT4),DP4C(MAX4,MAX4,MAX4,MAX4,NTOT4),
     3DP4D(MAX4,MAX4,MAX4,MAX4,NTOT4)
      V=0
      DO I=1,NMODE
        DK(I)=0
      END DO
C**CORRECTION FOR 1-DIM, 2-DIM AND 3-DIM POTENTIALS
CC    CALL GETQP3(V,DK,NMODE,QQ,NP1,GP1,VP1,DP1,NTOT1,MMX1,
CC   1NP2,GP2,VP2,DP2A,DP2B,NTOT2,MAX2,NP3,GP3,VP3,DP3A,DP3B,DP3C,
CC   2NTOT3,MAX3,INDK,INDL)
C**CORRECTION FOR 4-DIM POTENTIALS
      DO K=1,NMODE
        IF(QQ(K).NE.0)THEN
          QK=QQ(K)
          QKP=QK+1.D-6
          QKM=QK-1.D-6
          DO L=1,K-1
            IF(QQ(L).NE.0)THEN
              QL=QQ(L)
              QLP=QL+1.D-6
              QLM=QL-1.D-6
              DO N=1,L-1
                IF(QQ(N).NE.0)THEN
                  QN=QQ(N)
                  QNP=QN+1.D-6
                  QNM=QN-1.D-6
                  DO M=1,N-1
                    IF(QQ(M).NE.0)THEN
                      QM=QQ(M)
                      QMP=QM+1.D-6
                      QMM=QM-1.D-6
                      CALL GETHP4(VKLNM,K,L,N,M,QK,QL,QN,QM,NP4,GP4,VP4,
     1                DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,INDL,INDN)
                      V=V+VKLNM
C***************
CC    GO TO 1000
C***************
                      CALL GETHP4(VKP,K,L,N,M,QKP,QL,QN,QM,NP4,GP4,VP4,
     1                DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,INDL,INDN)
                      CALL GETHP4(VKM,K,L,N,M,QKM,QL,QN,QM,NP4,GP4,VP4,
     1                DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,INDL,INDN)
                      DK(K)=DK(K)+(VKP-VKM)*1.D6/2
                      CALL GETHP4(VLP,K,L,N,M,QK,QLP,QN,QM,NP4,GP4,VP4,
     1                DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,INDL,INDN)
                      CALL GETHP4(VLM,K,L,N,M,QK,QLM,QN,QM,NP4,GP4,VP4,
     1                DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,INDL,INDN)
                      DK(L)=DK(L)+(VLP-VLM)*1.D6/2
                      CALL GETHP4(VNP,K,L,N,M,QK,QL,QNP,QM,NP4,GP4,VP4,
     1                DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,INDL,INDN)
                      CALL GETHP4(VNM,K,L,N,M,QK,QL,QNM,QM,NP4,GP4,VP4,
     1                DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,INDL,INDN)
                      DK(N)=DK(N)+(VNP-VNM)*1.D6/2
                      CALL GETHP4(VMP,K,L,N,M,QK,QL,QN,QMP,NP4,GP4,VP4,
     1                DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,INDL,INDN)
                      CALL GETHP4(VMM,K,L,N,M,QK,QL,QN,QMM,NP4,GP4,VP4,
     1                DP4A,DP4B,DP4C,DP4D,NTOT4,MAX4,INDK,INDL,INDN)
                      DK(M)=DK(M)+(VMP-VMM)*1.D6/2
1000  CONTINUE
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END DO
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETHP4(V,MODEK,MODEL,MODEN,MODEM,QK,QL,QN,QM,NP4,HP4,
     1V4,D4A,D4B,D4C,D4D,NTOT4,MAX4,INDK,INDL,INDN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION INDK(1),INDL(1),INDN(1)
      DIMENSION NP4(NTOT4),HP4(MAX4,1),
     1V4(MAX4,MAX4,MAX4,MAX4,NTOT4),D4A(MAX4,MAX4,MAX4,MAX4,NTOT4),
     2D4B(MAX4,MAX4,MAX4,MAX4,NTOT4),D4C(MAX4,MAX4,MAX4,MAX4,NTOT4),
     3D4D(MAX4,MAX4,MAX4,MAX4,NTOT4)
C*************************************************
      dimension h(0:1,0:1,0:1,0:1),hx(0:1,0:1,0:1,0:1)
      dimension hy(0:1,0:1,0:1,0:1),hz(0:1,0:1,0:1,0:1)
      dimension hw(0:1,0:1,0:1,0:1)
C*************************************************
      COMMON/FILASS/IOUT,INP
      V=0
C**MODE K WILL BE X
      K=MODEK
      X=QK
C**MODE L WILL BE Y
      L=MODEL
      Y=QL
C**MODE N WILL BE Z
      N=MODEN
      Z=QN
C**MODE M WILL BE W
      M=MODEM
      W=QM
C**FIND X ORIGIN (K)
      DO I=1,NP4(K)
        IF(X.LE.HP4(I,K))THEN
          MK=I
          GO TO 41
        END IF
      END DO
41    CONTINUE
C**FIND Y ORIGIN (L)
      DO I=1,NP4(L)
        IF(Y.LE.HP4(I,L))THEN
          ML=I
          GO TO 42
        END IF
      END DO
42    CONTINUE
C**FIND Z ORIGIN (N)
      DO I=1,NP4(N)
        IF(Z.LE.HP4(I,N))THEN
          MN=I
          GO TO 43
        END IF
      END DO
43    CONTINUE
C**FIND W ORIGIN (M)
      DO I=1,NP4(M)
        IF(W.LE.HP4(I,M))THEN
          MM=I
          GO TO 44
        END IF
      END DO
44    CONTINUE
c  get x relative to x(k-1)
      XIR=X-HP4(MK-1,K)
c  get y relative to y(l-1)
      YIR=Y-HP4(ML-1,L)
c  get z relative to z(n-1)
      ZIR=Z-HP4(MN-1,N)
c  get w relative to w(m-1)
      WIR=W-HP4(MM-1,M)
C**GET INTERVAL IN X
      DX=HP4(MK,K)-HP4(MK-1,K)
C**GET INTERVAL IN Y
      DY=HP4(ML,L)-HP4(ML-1,L)
C**GET INTERVAL IN Z
      DZ=HP4(MN,N)-HP4(MN-1,N)
C**GET INTERVAL IN W
      DW=HP4(MM,M)-HP4(MM-1,M)
C
      IND=INDN(K)+INDL(L)+INDK(N)+M
      H(0,0,0,0)=V4(MM-1,MN-1,ML-1,MK-1,IND)
      H(1,0,0,0)=V4(MM-1,MN-1,ML-1,MK,IND)
      H(0,1,0,0)=V4(MM-1,MN-1,ML,MK-1,IND)
      H(0,0,1,0)=V4(MM-1,MN,ML-1,MK-1,IND)
      H(0,0,0,1)=V4(MM,MN-1,ML-1,MK-1,IND)
      H(1,1,0,0)=V4(MM-1,MN-1,ML,MK,IND)
      H(1,0,1,0)=V4(MM-1,MN,ML-1,MK,IND)
      H(1,0,0,1)=V4(MM,MN-1,ML-1,MK,IND)
      H(0,1,1,0)=V4(MM-1,MN,ML,MK-1,IND)
      H(0,1,0,1)=V4(MM,MN-1,ML,MK-1,IND)
      H(0,0,1,1)=V4(MM,MN,ML-1,MK-1,IND)
      H(1,1,1,0)=V4(MM-1,MN,ML,MK,IND)
      H(1,1,0,1)=V4(MM,MN-1,ML,MK,IND)
      H(1,0,1,1)=V4(MM,MN,ML-1,MK,IND)
      H(0,1,1,1)=V4(MM,MN,ML,MK-1,IND)
      H(1,1,1,1)=V4(MM,MN,ML,MK,IND)
C
      HX(0,0,0,0)=D4A(MK-1,MM-1,MN-1,ML-1,IND)
      HX(1,0,0,0)=D4A(MK,MM-1,MN-1,ML-1,IND)
      HX(0,1,0,0)=D4A(MK-1,MM-1,MN-1,ML,IND)
      HX(0,0,1,0)=D4A(MK-1,MM-1,MN,ML-1,IND)
      HX(0,0,0,1)=D4A(MK-1,MM,MN-1,ML-1,IND)
      HX(1,1,0,0)=D4A(MK,MM-1,MN-1,ML,IND)
      HX(1,0,1,0)=D4A(MK,MM-1,MN,ML-1,IND)
      HX(1,0,0,1)=D4A(MK,MM,MN-1,ML-1,IND)
      HX(0,1,1,0)=D4A(MK-1,MM-1,MN,ML,IND)
      HX(0,1,0,1)=D4A(MK-1,MM,MN-1,ML,IND)
      HX(0,0,1,1)=D4A(MK-1,MM,MN,ML-1,IND)
      HX(1,1,1,0)=D4A(MK,MM-1,MN,ML,IND)
      HX(1,1,0,1)=D4A(MK,MM,MN-1,ML,IND)
      HX(1,0,1,1)=D4A(MK,MM,MN,ML-1,IND)
      HX(0,1,1,1)=D4A(MK-1,MM,MN,ML,IND)
      HX(1,1,1,1)=D4A(MK,MM,MN,ML,IND)
C
      HY(0,0,0,0)=D4B(ML-1,MK-1,MM-1,MN-1,IND)
      HY(1,0,0,0)=D4B(ML-1,MK,MM-1,MN-1,IND)
      HY(0,1,0,0)=D4B(ML,MK-1,MM-1,MN-1,IND)
      HY(0,0,1,0)=D4B(ML-1,MK-1,MM-1,MN,IND)
      HY(0,0,0,1)=D4B(ML-1,MK-1,MM,MN-1,IND)
      HY(1,1,0,0)=D4B(ML,MK,MM-1,MN-1,IND)
      HY(1,0,1,0)=D4B(ML-1,MK,MM-1,MN,IND)
      HY(1,0,0,1)=D4B(ML-1,MK,MM,MN-1,IND)
      HY(0,1,1,0)=D4B(ML,MK-1,MM-1,MN,IND)
      HY(0,1,0,1)=D4B(ML,MK-1,MM,MN-1,IND)
      HY(0,0,1,1)=D4B(ML-1,MK-1,MM,MN,IND)
      HY(1,1,1,0)=D4B(ML,MK,MM-1,MN,IND)
      HY(1,1,0,1)=D4B(ML,MK,MM,MN-1,IND)
      HY(1,0,1,1)=D4B(ML-1,MK,MM,MN,IND)
      HY(0,1,1,1)=D4B(ML,MK-1,MM,MN,IND)
      HY(1,1,1,1)=D4B(ML,MK,MM,MN,IND)
C
      HZ(0,0,0,0)=D4C(MN-1,ML-1,MK-1,MM-1,IND)
      HZ(1,0,0,0)=D4C(MN-1,ML-1,MK,MM-1,IND)
      HZ(0,1,0,0)=D4C(MN-1,ML,MK-1,MM-1,IND)
      HZ(0,0,1,0)=D4C(MN,ML-1,MK-1,MM-1,IND)
      HZ(0,0,0,1)=D4C(MN-1,ML-1,MK-1,MM,IND)
      HZ(1,1,0,0)=D4C(MN-1,ML,MK,MM-1,IND)
      HZ(1,0,1,0)=D4C(MN,ML-1,MK,MM-1,IND)
      HZ(1,0,0,1)=D4C(MN-1,ML-1,MK,MM,IND)
      HZ(0,1,1,0)=D4C(MN,ML,MK-1,MM-1,IND)
      HZ(0,1,0,1)=D4C(MN-1,ML,MK-1,MM,IND)
      HZ(0,0,1,1)=D4C(MN,ML-1,MK-1,MM,IND)
      HZ(1,1,1,0)=D4C(MN,ML,MK,MM-1,IND)
      HZ(1,1,0,1)=D4C(MN-1,ML,MK,MM,IND)
      HZ(1,0,1,1)=D4C(MN,ML-1,MK,MM,IND)
      HZ(0,1,1,1)=D4C(MN,ML,MK-1,MM,IND)
      HZ(1,1,1,1)=D4C(MN,ML,MK,MM,IND)
C
      HW(0,0,0,0)=D4D(MM-1,MN-1,ML-1,MK-1,IND)
      HW(1,0,0,0)=D4D(MM-1,MN-1,ML-1,MK,IND)
      HW(0,1,0,0)=D4D(MM-1,MN-1,ML,MK-1,IND)
      HW(0,0,1,0)=D4D(MM-1,MN,ML-1,MK-1,IND)
      HW(0,0,0,1)=D4D(MM,MN-1,ML-1,MK-1,IND)
      HW(1,1,0,0)=D4D(MM-1,MN-1,ML,MK,IND)
      HW(1,0,1,0)=D4D(MM-1,MN,ML-1,MK,IND)
      HW(1,0,0,1)=D4D(MM,MN-1,ML-1,MK,IND)
      HW(0,1,1,0)=D4D(MM-1,MN,ML,MK-1,IND)
      HW(0,1,0,1)=D4D(MM,MN-1,ML,MK-1,IND)
      HW(0,0,1,1)=D4D(MM,MN,ML-1,MK-1,IND)
      HW(1,1,1,0)=D4D(MM-1,MN,ML,MK,IND)
      HW(1,1,0,1)=D4D(MM,MN-1,ML,MK,IND)
      HW(1,0,1,1)=D4D(MM,MN,ML-1,MK,IND)
      HW(0,1,1,1)=D4D(MM,MN,ML,MK-1,IND)
      HW(1,1,1,1)=D4D(MM,MN,ML,MK,IND)
C***************************************
      V=V+hermi4(dx,dy,dz,dw,h,hx,hy,hz,hw,xir,yir,zir,wir)
C************************************************
      RETURN
      END
C****************************************************************
C****************************************************************
      double precision function hermi1 (dx, f, fx, x)
      implicit real*8 (a-h,o-z)
      dimension  f(0:1), fx(0:1)
*     ------------------------------------------------------------------
*     Compute a 1-D Hermite interpolant.
*     ------------------------------------------------------------------
*   ..procedures
      gx(i) = dx*fx(i)-(f(1)-f(0))
*     ------------------------------------------------------------------
*   ..multilinear term
      s = ((dx-x)*f(0)+x*f(1))/dx
*   ..multicubic term
      s = s + (dx-x)*x*((dx-x)*gx(0)-x*gx(1))/dx**3
*   ..result
      hermi1 = s
      return
*     ------------------------------------------------------------------
*     B. J. Braams, Courant Institute, NYU.
      end
C****************************************************************
C****************************************************************
      double precision function hermi2 (dx, dy, f, fx, fy, x, y)
      implicit real*8 (a-h,o-z)
      dimension f(0:1,0:1), fx(0:1,0:1), fy(0:1,0:1)
*     ------------------------------------------------------------------
*     Compute a 2-D Hermite interpolant on a rectangle.
*     ------------------------------------------------------------------
*   ..procedures
      gx(i,j) = dx*fx(i,j)-(f(1,j)-f(0,j))
      gy(i,j) = dy*fy(i,j)-(f(i,1)-f(i,0))
*     ------------------------------------------------------------------
*   ..multilinear term
      s = ((dy-y)*((dx-x)*f(0,0)+x*f(1,0))+
     .          y*((dx-x)*f(0,1)+x*f(1,1)))/(dx*dy)
*   ..multicubic (x) term
      s = s + (dx-x)*x*
     .  ((dy-y)*((dx-x)*gx(0,0)-x*gx(1,0))+
     .        y*((dx-x)*gx(0,1)-x*gx(1,1)))/(dx**3*dy)
*   ..multicubic (y) term
      s = s + (dy-y)*y*
     .  ((dy-y)*((dx-x)*gy(0,0)+x*gy(1,0))-
     .        y*((dx-x)*gy(0,1)+x*gy(1,1)))/(dx*dy**3)
*   ..result
      hermi2 = s
      return
*     ------------------------------------------------------------------
*     B. J. Braams, Courant Institute, NYU.
      end       
C****************************************************************
C****************************************************************
      double precision function hermi3 (dx,dy,dz,f,fx,fy,fz,x,y,z)
      implicit real*8 (a-h,o-z)
      dimension f(0:1,0:1,0:1), fx(0:1,0:1,0:1),
     .  fy(0:1,0:1,0:1), fz(0:1,0:1,0:1)
*     ------------------------------------------------------------------
*     Compute a 3-D Hermite interpolant on a box.
*     ------------------------------------------------------------------
*   ..procedures
      gx(i,j,k) = dx*fx(i,j,k)-(f(1,j,k)-f(0,j,k))
      gy(i,j,k) = dy*fy(i,j,k)-(f(i,1,k)-f(i,0,k))
      gz(i,j,k) = dz*fz(i,j,k)-(f(i,j,1)-f(i,j,0))
*     ------------------------------------------------------------------
*   ..multilinear term
      s = ((dz-z)*((dy-y)*((dx-x)*f(0,0,0)+x*f(1,0,0))+
     .                  y*((dx-x)*f(0,1,0)+x*f(1,1,0)))+
     .          z*((dy-y)*((dx-x)*f(0,0,1)+x*f(1,0,1))+
     .                  y*((dx-x)*f(0,1,1)+x*f(1,1,1))))/(dx*dy*dz)
*   ..multicubic (x) term
      s = s + (dx-x)*x*
     .  ((dz-z)*((dy-y)*((dx-x)*gx(0,0,0)-x*gx(1,0,0))+
     .                y*((dx-x)*gx(0,1,0)-x*gx(1,1,0)))+
     .        z*((dy-y)*((dx-x)*gx(0,0,1)-x*gx(1,0,1))+
     .                y*((dx-x)*gx(0,1,1)-x*gx(1,1,1))))/(dx**3*dy*dz)
*   ..multicubic (y) term
      s = s + (dy-y)*y*
     .  ((dz-z)*((dy-y)*((dx-x)*gy(0,0,0)+x*gy(1,0,0))-
     .                y*((dx-x)*gy(0,1,0)+x*gy(1,1,0)))+
     .        z*((dy-y)*((dx-x)*gy(0,0,1)+x*gy(1,0,1))-
     .                y*((dx-x)*gy(0,1,1)+x*gy(1,1,1))))/(dx*dy**3*dz)
*   ..multicubic (z) term
      s = s + (dz-z)*z*
     .  ((dz-z)*((dy-y)*((dx-x)*gz(0,0,0)+x*gz(1,0,0))+
     .                y*((dx-x)*gz(0,1,0)+x*gz(1,1,0)))-
     .        z*((dy-y)*((dx-x)*gz(0,0,1)+x*gz(1,0,1))+
     .                y*((dx-x)*gz(0,1,1)+x*gz(1,1,1))))/(dx*dy*dz**3)
*   ..result
      hermi3 = s
      return
*     ------------------------------------------------------------------
*     B. J. Braams, Courant Institute, NYU.
      end
C****************************************************************
C****************************************************************
      double precision function hermi4 (dx, dy, dz, dw, f, fx, fy, fz, 
     .fw, x, y, z, w)
      implicit double precision (a-h,o-z)
      dimension f(0:1,0:1,0:1,0:1), fx(0:1,0:1,0:1,0:1),
     .  fy(0:1,0:1,0:1,0:1), fz(0:1,0:1,0:1,0:1), fw(0:1,0:1,0:1,0:1)
*     ------------------------------------------------------------------
*     Compute a 4-D Hermite interpolant on a box.
*     ------------------------------------------------------------------
*   ..local variables
      dimension s(0:1), sw(0:1)
*   ..procedures
      gw(i,j,k,l) = dw*fw(i,j,k,l)-(f(i,j,k,1)-f(i,j,k,0))
*     ------------------------------------------------------------------
*   ..interpolants at w=0 and w=dw
      do l = 0, 1
       s(l) = hermi3 (dx, dy, dz, f(0,0,0,l), fx(0,0,0,l), fy(0,0,0,l),
     .   fz(0,0,0,l), x, y, z)
      enddo
*   ..multicubic (w) correction term
      do l = 0, 1
       sw(l) =
     .   ((dz-z)*((dy-y)*((dx-x)*gw(0,0,0,l)+x*gw(1,0,0,l))+
     .                 y*((dx-x)*gw(0,1,0,l)+x*gw(1,1,0,l)))+
     .         z*((dy-y)*((dx-x)*gw(0,0,1,l)+x*gw(1,0,1,l))+
     .                 y*((dx-x)*gw(0,1,1,l)+x*gw(1,1,1,l))))/(dx*dy*dz)
      enddo
*   ..result
      hermi4 = ((dw-w)*s(0)+w*s(1))/dw +
     .  (dw-w)*w*((dw-w)*sw(0)-w*sw(1))/dw**3
      return
*     ------------------------------------------------------------------
*     B. J. Braams, Courant Institute, NYU.
      end
C***********************************************************************
C**STILL TO BE CONVERTED TO DOUBLE PRECISION - AND INCLUDED INTO CODE
C***********************************************************************
      subroutine hgrad1 (nx, lx, dx, f, fx, x, r, rx)
      implicit none
      integer nx, lx
      real dx, f(0:nx), fx(0:nx), x, r, rx
*     ------------------------------------------------------------------
*     Compute a 1-D Hermite interpolant and its gradient.
*     ------------------------------------------------------------------
*     ..local variables
      integer i
      real s, sx, t0
*     ..procedures
      real gx
      gx(i) = dx*fx(i)-(f(lx+1)-f(lx))
*     ------------------------------------------------------------------
*     ..test inputs
      if (.not.(0.le.lx.and.lx.lt.nx)) then
       stop 'hgrad1--bad lx'
      endif
*     ..linear term
      s = ((dx-x)*f(lx)+x*f(lx+1))/dx
      sx = (f(lx+1)-f(lx))/dx
*     ..cubic term
      t0 = ((dx-x)*gx(lx)-x*gx(lx+1))/dx**3
      s = s + (dx-x)*x*t0
      sx = sx + (dx-2*x)*t0-(dx-x)*x*(gx(lx)+gx(lx+1))/dx**3
*     ..result
      r = s
      rx = sx
      return
*     ------------------------------------------------------------------
*     B. J. Braams, Courant Institute, NYU.
      end
      subroutine hgrad2 (nx, ny, lx, ly, dx, dy, f, fx, fy, x, y,
     .  r, rx, ry)
      implicit none
      integer nx, ny, lx, ly
      real dx, dy, f(0:nx,0:ny), fx(0:nx,0:ny), fy(0:nx,0:ny), x, y,
     .  r, rx, ry
*     ------------------------------------------------------------------
*     Compute a 2-D Hermite interpolant and its gradient.
*     ------------------------------------------------------------------
*     ..local variables
      integer i, j
      real s, sx, sy, t0
*     ..procedures
      real gx, gy
      gx(i,j) = dx*fx(i,j)-(f(lx+1,j)-f(lx,j))
      gy(i,j) = dy*fy(i,j)-(f(i,ly+1)-f(i,ly))
*     ------------------------------------------------------------------
*     ..test inputs
      if (.not.(0.le.lx.and.lx.lt.nx.and.0.le.ly.and.ly.lt.ny)) then
       stop 'hgrad2--bad lx, ly'
      endif
*     ..multilinear term
      s = (1-y/dy)*((1-x/dx)*f(lx,ly)+(x/dx)*f(lx+1,ly))+
     .      (y/dy)*((1-x/dx)*f(lx,ly+1)+(x/dx)*f(lx+1,ly+1))
      sx = (1-y/dy)*(f(lx+1,ly)-f(lx,ly))+
     .       (y/dy)*(f(lx+1,ly+1)-f(lx,ly+1))
      sy = (1-x/dx)*(f(lx,ly+1)-f(lx,ly))+
     .       (x/dx)*(f(lx+1,ly+1)-f(lx+1,ly))
*     ..multicubic (x) term
      t0 = (1-y/dy)*((1-x/dx)*gx(lx,ly)-(x/dx)*gx(lx+1,ly))+
     .       (y/dy)*((1-x/dx)*gx(lx,ly+1)-(x/dx)*gx(lx+1,ly+1))
      s = s + (1-x/dx)*(x/dx)*t0
      sx = sx + (1-2*(x/dx))*t0 -
     .  (1-x/dx)*(x/dx)*
     .  ((1-y/dy)*(gx(lx,ly)+gx(lx+1,ly))+
     .     (y/dy)*(gx(lx,ly+1)+gx(lx+1,ly+1)))
      sy = sy + (1-x/dx)*(x/dx)*
     .  ((1-x/dx)*(gx(lx,ly+1)-gx(lx,ly))-
     .     (x/dx)*(gx(lx+1,ly+1)-gx(lx+1,ly)))
*     ..multicubic (y) term
      t0 = ((1-y/dy)*((1-x/dx)*gy(lx,ly)+(x/dx)*gy(lx+1,ly))-
     .        (y/dy)*((1-x/dx)*gy(lx,ly+1)+(x/dx)*gy(lx+1,ly+1)))
      s = s + (1-y/dy)*(y/dy)*t0
      sx = sx + (1-y/dy)*(y/dy)*
     .  ((1-y/dy)*(gy(lx+1,ly)-gy(lx,ly))-
     .     (y/dy)*(gy(lx+1,ly+1)-gy(lx,ly+1)))
      sy = sy + (1-2*(y/dy))*t0 -
     .  (1-y/dy)*(y/dy)*
     .  ((1-x/dx)*(gy(lx,ly)+gy(lx,ly+1))+
     .     (x/dx)*(gy(lx+1,ly)+gy(lx+1,ly+1)))
*     ..result
      r = s
      rx = sx/dx
      ry = sy/dy
      return
*     ------------------------------------------------------------------
*     B. J. Braams, Courant Institute, NYU.
      end
      subroutine hgrad3 (nx, ny, nz, lx, ly, lz, dx, dy, dz,
     .  f, fx, fy, fz, x, y, z, r, rx, ry, rz)
      implicit none
      integer nx, ny, nz, lx, ly, lz
      real dx, dy, dz, f(0:nx,0:ny,0:nz), fx(0:nx,0:ny,0:nz),
     .  fy(0:nx,0:ny,0:nz), fz(0:nx,0:ny,0:nz), x, y, z, r, rx, ry, rz
*     ------------------------------------------------------------------
*     Compute a 3-D Hermite interpolant and its gradient.
*     ------------------------------------------------------------------
*     ..local variables
      integer i, j, k
      real s, sx, sy, sz, f1(0:1,0:1), f1x(0:1,0:1), f1y(0:1,0:1),
     .  f1z(0:1,0:1), f1xz(0:1,0:1), f1yz(0:1,0:1), f2(0:1,0:1),
     .  f2z(0:1,0:1)
*     ..procedures
      real hermi2, gz
      external hgrad2, hermi2
      gz(i,j,k) = dz*fz(i,j,k)-(f(i,j,lz+1)-f(i,j,lz))
*     ------------------------------------------------------------------
*     ..test inputs
      if (.not.(0.le.lx.and.lx.lt.nx.and.0.le.ly.and.ly.lt.ny.and.
     .  0.le.lz.and.lz.lt.nz)) then
       stop 'hgrad3--bad lx, ly, lz'
      endif
*     ..compute auxiliary quantities
*     f1 for terms linear in z
*     f2 for multicubic (z) correction term
      do j = 0, 1
       do i = 0, 1
        f1(i,j) =
     .    (1-z/dz)*f(lx+i,ly+j,lz)+
     .      (z/dz)*f(lx+i,ly+j,lz+1)
        f1x(i,j) =
     .    (1-z/dz)*fx(lx+i,ly+j,lz)+
     .      (z/dz)*fx(lx+i,ly+j,lz+1)
        f1y(i,j) =
     .    (1-z/dz)*fy(lx+i,ly+j,lz)+
     .      (z/dz)*fy(lx+i,ly+j,lz+1)
        f1z(i,j) = (f(lx+i,ly+j,lz+1)-f(lx+i,ly+j,lz))/dz
        f1xz(i,j) = (fx(lx+i,ly+j,lz+1)-fx(lx+i,ly+j,lz))/dz
        f1yz(i,j) = (fy(lx+i,ly+j,lz+1)-fy(lx+i,ly+j,lz))/dz
        f2(i,j) = (1-z/dz)*(z/dz)*
     .    ((1-z/dz)*gz(lx+i,ly+j,lz)-(z/dz)*gz(lx+i,ly+j,lz+1))
        f2z(i,j) =
     .    ((1-z/dz)*(1-3*(z/dz))*gz(lx+i,ly+j,lz)-
     .      (z/dz)*(2-3*(z/dz))*gz(lx+i,ly+j,lz+1))/dz
       enddo
      enddo
*     ..terms linear in z
      call hgrad2 (1, 1, 0, 0, dx, dy, f1, f1x, f1y, x, y, s, sx, sy)
      sz = hermi2 (dx, dy, f1z, f1xz, f1yz, x, y)
*     ..result
      r = s +
     .  (1-y/dy)*((1-x/dx)*f2(0,0)+(x/dx)*f2(1,0))+
     .    (y/dy)*((1-x/dx)*f2(0,1)+(x/dx)*f2(1,1))
      rx = sx +
     .  ((1-y/dy)*(f2(1,0)-f2(0,0))+
     .     (y/dy)*(f2(1,1)-f2(0,1)))/dx
      ry = sy +
     .  ((1-x/dx)*(f2(0,1)-f2(0,0))+
     .     (x/dx)*(f2(1,1)-f2(1,0)))/dy
      rz = sz +
     .  (1-y/dy)*((1-x/dx)*f2z(0,0)+(x/dx)*f2z(1,0))+
     .    (y/dy)*((1-x/dx)*f2z(0,1)+(x/dx)*f2z(1,1))
      return
*     ------------------------------------------------------------------
*     B. J. Braams, Courant Institute, NYU.
      end
      subroutine hgrad4 (nx, ny, nz, nw, lx, ly, lz, lw, dx, dy, dz,
     .  dw, f, fx, fy, fz, fw, x, y, z, w, r, rx, ry, rz, rw)
      implicit none
      integer nx, ny, nz, nw, lx, ly, lz, lw
      real dx, dy, dz, dw, f(0:nx,0:ny,0:nz,0:nw),
     .  fx(0:nx,0:ny,0:nz,0:nw), fy(0:nx,0:ny,0:nz,0:nw),
     .  fz(0:nx,0:ny,0:nz,0:nw), fw(0:nx,0:ny,0:nz,0:nw),
     .  x, y, z, w, r, rx, ry, rz, rw
*     ------------------------------------------------------------------
*     Compute a 4-D Hermite interpolant and its gradient.
*     ------------------------------------------------------------------
*     ..local variables
      integer i, j, k, l
      real s, sx, sy, sz, sw, f1(0:1,0:1,0:1), f1x(0:1,0:1,0:1),
     .  f1y(0:1,0:1,0:1), f1z(0:1,0:1,0:1), f1w(0:1,0:1,0:1),
     .  f1xw(0:1,0:1,0:1), f1yw(0:1,0:1,0:1), f1zw(0:1,0:1,0:1),
     .  f2(0:1,0:1,0:1), f2w(0:1,0:1,0:1)
*     ..procedures
      real hermi3, gw
      external hgrad3, hermi3
      gw(i,j,k,l) = dw*fw(i,j,k,l)-(f(i,j,k,lw+1)-f(i,j,k,lw))
*     ------------------------------------------------------------------
*     ..test inputs
      if (.not.(0.le.lx.and.lx.lt.nx.and.0.le.ly.and.ly.lt.ny.and.
     .  0.le.lz.and.lz.lt.nz.and.0.le.lw.and.lw.lt.nw)) then
       stop 'hgrad4--bad lx, ly, lz, lw'
      endif
*     ..compute auxiliary quantities
*     f1 for terms linear in w
*     f2 for multicubic (w) correction term
      do k = 0, 1
       do j = 0, 1
        do i = 0, 1
         f1(i,j,k) =
     .     (1-w/dw)*f(lx+i,ly+j,lz+k,lw)+
     .       (w/dw)*f(lx+i,ly+j,lz+k,lw+1)
         f1x(i,j,k) =
     .     (1-w/dw)*fx(lx+i,ly+j,lz+k,lw)+
     .       (w/dw)*fx(lx+i,ly+j,lz+k,lw+1)
         f1y(i,j,k) =
     .     (1-w/dw)*fy(lx+i,ly+j,lz+k,lw)+
     .       (w/dw)*fy(lx+i,ly+j,lz+k,lw+1)
         f1z(i,j,k) =
     .     (1-w/dw)*fz(lx+i,ly+j,lz+k,lw)+
     .       (w/dw)*fz(lx+i,ly+j,lz+k,lw+1)
         f1w(i,j,k) =
     .     (f(lx+i,ly+j,lz+k,lw+1)-f(lx+i,ly+j,lz+k,lw))/dw
         f1xw(i,j,k) =
     .     (fx(lx+i,ly+j,lz+k,lw+1)-fx(lx+i,ly+j,lz+k,lw))/dw
         f1yw(i,j,k) =
     .     (fy(lx+i,ly+j,lz+k,lw+1)-fy(lx+i,ly+j,lz+k,lw))/dw
         f1zw(i,j,k) =
     .     (fz(lx+i,ly+j,lz+k,lw+1)-fz(lx+i,ly+j,lz+k,lw))/dw
         f2(i,j,k) = (1-w/dw)*(w/dw)*
     .     ((1-w/dw)*gw(lx+i,ly+j,lz+k,lw)-
     .        (w/dw)*gw(lx+i,ly+j,lz+k,lw+1))
         f2w(i,j,k) =
     .     ((1-w/dw)*(1-3*(w/dw))*gw(lx+i,ly+j,lz+k,lw)-
     .        (w/dw)*(2-3*(w/dw))*gw(lx+i,ly+j,lz+k,lw+1))/dw
        enddo
       enddo
      enddo
*     ..terms linear in w
      call hgrad3 (1, 1, 1, 0, 0, 0, dx, dy, dz, f1, f1x, f1y, f1z,
     .  x, y, z, s, sx, sy, sz)
      sw = hermi3 (dx, dy, dz, f1w, f1xw, f1yw, f1zw, x, y, z)
*     ..result
      r = s +
     .  (1-z/dz)*((1-y/dy)*((1-x/dx)*f2(0,0,0)+(x/dx)*f2(1,0,0))+
     .              (y/dy)*((1-x/dx)*f2(0,1,0)+(x/dx)*f2(1,1,0)))+
     .    (z/dz)*((1-y/dy)*((1-x/dx)*f2(0,0,1)+(x/dx)*f2(1,0,1))+
     .              (y/dy)*((1-x/dx)*f2(0,1,1)+(x/dx)*f2(1,1,1)))
      rx = sx +
     .  ((1-z/dz)*((1-y/dy)*(f2(1,0,0)-f2(0,0,0))+
     .               (y/dy)*(f2(1,1,0)-f2(0,1,0)))+
     .     (z/dz)*((1-y/dy)*(f2(1,0,1)-f2(0,0,1))+
     .               (y/dy)*(f2(1,1,1)-f2(0,1,1))))/dx
      ry = sy +
     .  ((1-z/dz)*((1-x/dx)*(f2(0,1,0)-f2(0,0,0))+
     .               (x/dx)*(f2(1,1,0)-f2(1,0,0)))+
     .     (z/dz)*((1-x/dx)*(f2(0,1,1)-f2(0,0,1))+
     .               (x/dx)*(f2(1,1,1)-f2(1,0,1))))/dy
      rz = sz +
     .  ((1-y/dy)*((1-x/dx)*(f2(0,0,1)-f2(0,0,0))+
     .               (x/dx)*(f2(1,0,1)-f2(1,0,0)))+
     .     (y/dy)*((1-x/dx)*(f2(0,1,1)-f2(0,1,0))+
     .               (x/dx)*(f2(1,1,1)-f2(1,1,0))))/dz
      rw = sw +
     .  (1-z/dz)*((1-y/dy)*((1-x/dx)*f2w(0,0,0)+(x/dx)*f2w(1,0,0))+
     .              (y/dy)*((1-x/dx)*f2w(0,1,0)+(x/dx)*f2w(1,1,0)))+
     .    (z/dz)*((1-y/dy)*((1-x/dx)*f2w(0,0,1)+(x/dx)*f2w(1,0,1))+
     .              (y/dy)*((1-x/dx)*f2w(0,1,1)+(x/dx)*f2w(1,1,1)))
      return
*     ------------------------------------------------------------------
*     B. J. Braams, Courant Institute, NYU.
      end

