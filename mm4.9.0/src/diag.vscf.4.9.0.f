C****************************************************************
C****************************************************************
C**DIAG
C****************************************************************
C****************************************************************
      SUBROUTINE DIAG(XA,XK,NXN,NN,NV,SUP4,EVAL,WRK,NVAL,NVEC,
     1ISTAT,NSTAT,NMODE,OV,VEC,IASSIG,ISIZMX,J21,NS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PRT,LGIV,TRIAT
      DIMENSION SUP4(NN),EVAL(NN),WRK(NN)
      DIMENSION XK(NXN,1),XA(1),OV(NXN,1),VEC(NXN,1)
      DIMENSION ISTAT(NSTAT,NMODE),IASSIG(ISIZMX,J21,3,1)
      DIMENSION JJMAX(3),XMAX(4)
      DIMENSION NS1(3),NS2(3),IOFF1(3),IOFF2(3)
      COMMON/FILASS/IOUT
      COMMON/TRIATO/TRIAT
      COMMON/GIVEN/LGIV
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/EVL/EVL,CUT
      COMMON/RETURN/IRET
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/ROTS/JMAX,KMAX,IDUM,KEL21,KEL
      COMMON/REACTL/JREACT
      COMMON/MATRIX/NVALV,NVALR,KSTEP
      COMMON/JKAKC/JTHIS,KA,KC
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100)
      COMMON/RPHROT/IROTV
C**********************************************************
200   FORMAT(/,1X,'VIBRATIONAL ENERGY LEVELS (RECIPROCAL CMS)',/)
205   FORMAT(5(I4,': ',F10.2))
C**********************************************************
C     TZERO=TIMER(0D0)
C     WRITE(IOUT,*)'CALCULATING DIAG...'
C     CALL TIMIT(1)
      PRT=(IPRINT.GT.0)
C**********************************************************
      IF(LGIV)THEN
        IVECT=1000000
        CALL GIVENS(NN,NVAL,NXN,XA,SUP4,EVAL,XK,WRK,IVECT)
      ELSE
        DO 8 I=1,NN
        DO 8 J=1,I
8       XK(I,J)=XK(J,I)
        IF(NV.LE.0)THEN
          CALL TRID(NN,XK,NXN,SUP4,EVAL)
        ELSE
          DO 11 I=1,NN
          DO 11 J=1,I
11        OV(I,J)=OV(J,I)
          IFAIL=1
          CALL DSYGV(1,'V','U',NN,XK,NXN,OV,NXN,EVAL,VEC,3*NN-1,IFAIL)
          IF(IFAIL.NE.0)THEN
            WRITE(IOUT,*)'IFAIL = ',IFAIL
            WRITE(IOUT,*)'ERROR IN DSYGV'
            STOP 'ERROR IN DSYGV'
          END IF
          IF(IROTV.NE.0)THEN
C**SEARCH FOR DEGENERACIES
            DO J=1,NN
              IF(J.LT.NN)THEN
                IF(DABS(EVAL(J)-EVAL(J+1))*WAVENM.LT.1.D-3)THEN
C**SEE WHICH IS LARGER OF SIN(1) and COS(1) IN SECOND VECTOR
                  I2=2
                  I3=3
                  IF(DABS(XK(3,J+1)).GT.DABS(XK(2,J+1)))THEN
                    I2=3
                    I3=2
                  END IF
C**I2 IS POSITION OF 'ZERO' IN FIRST VECTOR
                  TH=XK(I2,J)/XK(I2,J+1)
                  TH=DATAN(TH)
                  CTH=DCOS(TH)
                  STH=DSIN(TH)
                  DO I=1,NN
                    TEMP1=CTH*XK(I,J)+STH*XK(I,J+1)
                    TEMP2=-STH*XK(I,J)+CTH*XK(I,J+1)
                    XK(I,J)=TEMP1
                    XK(I,J+1)=TEMP2
                  END DO
                END IF
              END IF
            END DO
            DO I=1,NN
              WRITE(IOUT,*)'ENERGY ',EVAL(I)*WAVENM
              WRITE(IOUT,*)'VECTOR ',I
              WRITE(IOUT,*)(XK(J,I),J=1,NN)
              CALL FLUSH(IOUT)
            END DO
          END IF
        END IF
      END IF
C     IF(IRET.NE.0.AND.NV.LE.0)RETURN
      ICUTX=NVAL
      IF(NVAL.EQ.0.OR..NOT.LGIV)ICUTX=NN
      NN1=ICUTX-1
      DO 3 I=1,NN1
      L=I+1
      DO 3 J=L,ICUTX
      IF(EVAL(I)-EVAL(J))3,3,1
1     A1=EVAL(I)
      EVAL(I)=EVAL(J)
      EVAL(J)=A1
      DO 2 K=1,NN
      A1=XK(K,I)
      XK(K,I)=XK(K,J)
2     XK(K,J)=A1
3     CONTINUE
C**********************************************************
C**RETURN FOR HERMPT AND VSCF (NOT FINAL MODE)
C**********************************************************
      IF(NV.LT.0.OR.IRET.NE.0)RETURN
      MT=1
      IF(EVL.NE.0.0D0)GO TO 4
      MT=2
4     CONTINUE
      DO 5 I=1,ICUTX
      WRK(I)=EVAL(I)
5     EVAL(I)=EVAL(I)*WAVENM
      IF(EVL.NE.0.0D0)GO TO 6
      EVL=EVAL(1)
6     CONTINUE
C**********************************************************
      IF(IPRINT.LT.0)THEN
      WRITE(IOUT,200)
      WRITE(IOUT,205)(I-MT+1,EVAL(I),I=1,ICUTX)
      END IF
C**********************************************************
      DO 7 I=1,NN
      EVAL(I)=EVAL(I)-EVL
7     CONTINUE
C**********************************************************
C**RETURN FOR VSCF (FINAL MODE)
C**********************************************************
      IF(NV.EQ.0.AND.ICID.EQ.0)RETURN
      XMAX(1)=2
      NS1(1)=NS
      NS2(1)=NS
      IF(LCOUNT.LT.0)THEN
C**GET SIZES OF TOTAL CONTRACTION SCHEME BASES
        ICSIZ1=0
        DO I=1,NVSYM
          ICSIZ1=ICSIZ1+ISIZC(1,I)
        END DO
        ICSIZ2=0
        DO I=1,NVSYM
          ICSIZ2=ICSIZ2+ISIZC(2,I)
        END DO
      END IF
      IF(ICID.EQ.0)THEN
        NSIZE=NN
C**SIZE OF SINGLE UNIT
CCCC    IF(LCOUNT.LT.0)NSIZE=NSTAT
        KKC=KC
        KKA=KA
      ELSE
C       J21=2*JTHIS+1
        NSIZE=NVALV
      END IF
      ICYCL=1
      IF(JPRINT.LT.-1.AND.ICID.EQ.0)ICYCL=2
      IF(JPRINT.LT.-2.AND.ICID.EQ.0)ICYCL=3
      DO 10 I=1,NN
CCCC  IF(I.GT.NVAL.AND.LCOUNT.LT.0)GO TO 1000
      IF(I.GT.NVAL)GO TO 1000
      IF(EVAL(I).GT.CUT.AND.LCOUNT.LT.0)GO TO 1000
      DO 100 ICYC=1,ICYCL
      NSX=NS
      IOFF=NSX-NVSYM
      IF(IOFF.GT.0)NSX=IOFF
      IF(ICID.NE.0)THEN
        KKC=JTHIS
        KKA=0
        IF(TRIAT.AND.NS.GT.NVSYM)THEN
          KKC=JTHIS-1
          KKA=1
        END IF
      END IF
      JX=I+1-MT
      if (icyc.lt.3) then
c  find two max ci coeffs
        XMAX(ICYC+1)=0.D0
        DO 9 J=1,NN
        IF(XK(J,I).EQ.XMAX(ICYC))GO TO 9
        IF(ABS(XK(J,I)).LE.ABS(XMAX(ICYC+1)))GO TO 9
        XMAX(ICYC+1)=XK(J,I)
        IMAX=J
9       CONTINUE
      else
c  get third biggest coeff
        xmax(4)=0.0
        do 90 j=1,NN
        if((xk(j,i).eq.xmax(2)).or.(xk(j,i).eq.xmax(3))) go to 90
        if (abs(xk(j,i)).lt.abs(xmax(4))) go to 90
        xmax(4)=xk(j,i)
        imax=j
90      continue
      endif
      ICMAX=IMAX
      IIA=0
      IOFF1(ICYC)=0
      IOFF2(ICYC)=0
C*****************************************
      IF(LCOUNT.LT.0)THEN
C**GET STARTING POINT
        IGOT=0
C**A1G
        IF(NS.EQ.1)THEN
          DO J=1,NVSYM
            IF(IGOT.EQ.0)THEN
              IF(ISIZC(1,J).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(J)
                IGOT=1
                IIA=J+1
                NS1(ICYC)=J
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**B2G
        IF(NS.EQ.2)THEN
          IS=1
          DO JJ=1,NVSYM
            J=JJ+IS
            IF(IGOT.EQ.0)THEN
              IS=-IS
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                IIS=IS
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**B1G
        IF(NS.EQ.3)THEN
          IS=1
          DO JJ=1,NVSYM
            J=JJ+2*IS
            IF(IGOT.EQ.0)THEN
              JJJ=MOD(JJ,2)+1
              IS=IS*(-1)**JJJ
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                IIS=IS
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**A2G
        IF(NS.EQ.4)THEN
          DO JJ=1,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
C           J=NVSYM+1-JJ
            J=5+NVSYM*KSR-JJ
            IF(IGOT.EQ.0)THEN
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**A1U
        IF(NS.EQ.5)THEN
          DO JJ=1,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
            J=JJ+4*LSR
            IF(IGOT.EQ.0)THEN
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**B2U
        IF(NS.EQ.6)THEN
          IS=1
          DO JJ=1,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
            J=JJ+IS+4*LSR
            IF(IGOT.EQ.0)THEN
              IS=-IS
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                IIS=IS
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**B1U
        IF(NS.EQ.7)THEN
          IS=1
          DO JJ=1,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
            J=JJ+2*IS+4*LSR
            IF(IGOT.EQ.0)THEN
              JJJ=MOD(JJ,2)+1
              IS=IS*(-1)**JJJ
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                IIS=IS
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**A2U
        IF(NS.EQ.8)THEN
          DO JJ=1,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
            J=5+NVSYM*KSR-JJ+4*LSR
            IF(IGOT.EQ.0)THEN
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
        DO JJJ=1,NS1(ICYC)-1
          IOFF1(ICYC)=IOFF1(ICYC)+ISIZC(1,JJJ)
        END DO
        DO JJJ=1,NS2(ICYC)-1
          IOFF2(ICYC)=IOFF2(ICYC)+ISIZC(2,JJJ)
        END DO
      END IF
C*****************************************

12    CONTINUE
      IF(NSIZE.GE.IMAX)GO TO 13
      IMAX=IMAX-NSIZE
      IF(LCOUNT.LT.0)THEN
        IIB=IIA
        IGOT=0
C**A1G
        IF(NS.EQ.1)THEN
          DO J=IIB,NVSYM
            IF(IGOT.EQ.0)THEN
              IF(ISIZC(1,J).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(J)
                IGOT=1
                IIA=J+1
                NS1(ICYC)=J
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**B2G
        IF(NS.EQ.2)THEN
          IS=IIS
          DO JJ=IIB,NVSYM
            J=JJ+IS
            IF(IGOT.EQ.0)THEN
              IS=-IS
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                IIS=IS
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**B1G
        IF(NS.EQ.3)THEN
          IS=IIS
          DO JJ=IIB,NVSYM
            J=JJ+2*IS
            IF(IGOT.EQ.0)THEN
              JJJ=MOD(JJ,2)+1
              IS=IS*(-1)**JJJ
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                IIS=IS
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**A2G
        IF(NS.EQ.4)THEN
          DO JJ=IIB,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
C           J=NVSYM+1-JJ
            J=5+NVSYM*KSR-JJ
            IF(IGOT.EQ.0)THEN
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**A1U
        IF(NS.EQ.5)THEN
          DO JJ=IIB,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
            J=JJ+4*LSR
            IF(IGOT.EQ.0)THEN
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**B2U
        IF(NS.EQ.6)THEN
          IS=IIS
          DO JJ=IIB,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
            J=JJ+IS+4*LSR
            IF(IGOT.EQ.0)THEN
              IS=-IS
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                IIS=IS
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**B1U
        IF(NS.EQ.7)THEN
          IS=IIS
          DO JJ=IIB,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
            J=JJ+2*IS+4*LSR
            IF(IGOT.EQ.0)THEN
              JJJ=MOD(JJ,2)+1
              IS=IS*(-1)**JJJ
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                IIS=IS
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
C**A2U
        IF(NS.EQ.8)THEN
          DO JJ=IIB,NVSYM
            KSR=JJ/5
            LSR=(-1)**(KSR+2)
            J=5+NVSYM*KSR-JJ+4*LSR
            IF(IGOT.EQ.0)THEN
              IF(ISIZC(1,JJ).NE.0.AND.ISIZC(2,J).NE.0)THEN
                NSIZE=NCSIZE(JJ)
                IGOT=1
                IIA=JJ+1
                NS1(ICYC)=JJ
                NS2(ICYC)=J
              END IF
            END IF
          END DO
        END IF
        IOFF1(ICYC)=0
        DO JJJ=1,NS1(ICYC)-1
          IOFF1(ICYC)=IOFF1(ICYC)+ISIZC(1,JJJ)
        END DO
        IOFF2(ICYC)=0
        DO JJJ=1,NS2(ICYC)-1
          IOFF2(ICYC)=IOFF2(ICYC)+ISIZC(2,JJJ)
        END DO
        GO TO 12
      END IF
      IIA=IIA+1
      IF(TRIAT)THEN
        KKA=KKA+1
        IF(MOD(IIA,2).EQ.0)KKC=KKC-2
      ELSE
        IA=MOD(IIA,2)
        KKA=KKA+IA
        KKC=JTHIS-KKA
        KKC=KKC+MOD(IIA+1,2)
      END IF
C**UPDATE VIBRATIONAL SYMMETRY FOR NEXT TIME IF REQUIRED
      IF((.NOT.TRIAT).AND.NVSYM.EQ.4)THEN
        IF(MOD(IIA,2).NE.0)THEN
          IF(JREACT.LE.0)THEN
            NSX=NSX+2
          ELSE
          END IF
        ELSE
          IF(JREACT.LE.0)THEN
            IF(MOD(IIA,4).NE.0)THEN
              NSX=NSX+(-1)**NS
            ELSE
              NSX=NSX-(-1)**NS
            END IF
          ELSE
            NSX=NSX+2
          END IF
        END IF
        IF(NSX.GT.NVSYM)THEN
          IF(JREACT.LE.0)THEN
            INCR=NSX-NVSYM
            NSX=INCR
          ELSE
            NSX=NSX-NVSYM
          END IF
        END IF
        IF(NSX.EQ.0)NSX=NVSYM
      ELSE
        NSX=NSX+1
        IF(NSX.GT.NVSYM)NSX=1
      END IF
      GO TO 12
13    CONTINUE
      IF(ICID.EQ.0)THEN
C**STORE CI ASSIGNMENTS
        KROT=KEL
        IF(LCOUNT.GE.0)THEN
          IASSIG(I,KROT,ICYC,NSX)=IMAX
          JJMAX(ICYC)=IMAX
        ELSE
          JJMAX(ICYC)=ICMAX
        END IF
      ELSE
C**RECALL CI ASSIGNMENTS
        IOFF=0
        IF(NSX.GT.1)THEN
          DO K=1,NSX-1
            IOFF=IOFF+NTOT(K)
          END DO
        END IF
        IF(TRIAT)THEN
          KROT=2*KKA
          KOFF=1
          IF(NS.GT.NVSYM)KOFF=0
          JROT=KROT+KOFF
        ELSE
          KROT=2*KKA
          IF(KKA.EQ.0.OR.KKA+KKC.NE.JTHIS)KROT=KROT+1
          JROT=KROT
        END IF
        KKEL=(JROT-1)/KSTEP+1
        JJMAX(1)=IOFF+IASSIG(IMAX,KKEL,1,NSX)
        IF(JPRINT.LT.-1)JJMAX(2)=IOFF+IASSIG(IMAX,KKEL,2,NSX)
        IF(JPRINT.LT.-2)JJMAX(3)=IOFF+IASSIG(IMAX,KKEL,3,NSX)
      END IF
100   CONTINUE
      IF(I.GT.NVAL)GO TO 1000
      IF(EVAL(I).GT.CUT)GO TO 1000
      CALL PRCI(ISTAT,NSTAT,NMODE,JJMAX(1),JJMAX(2),JJMAX(3),
     1XMAX(2),XMAX(3),XMAX(4),JX,EVAL(I),EVL,JTHIS,KKA,KKC,45,
     2OV,VEC,ICSIZ1,ICSIZ2,IOFF1,IOFF2,IASSIG,ISIZMX,J21,NV,NS1,NS2)
1000  CONTINUE
10    CONTINUE
C     TIME=TIMER(TZERO)
C     WRITE(IOUT,'(1X,''TIME TAKEN(DIAG)= '',F14.6,'' SECONDS'',/)
C    1')TIME
C     CALL TIMIT(3)
      RETURN
      END
C*******************************************************************
C*******************************************************************
      SUBROUTINE SVASS(IASSIG,ISIZMX,ICASS,ISCFCI,J21,NVSYM,LCONT,NS,
     1NVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IASSIG(ISIZMX,J21,3,1),ICASS(ISCFCI,J21,NVSYM,1)
      COMMON/ROTS/JMAX,KMAX,IDUM,KEL21,KROT
      COMMON/FILASS/IOUT
      DO I=1,NVAL
        ICASS(I,KROT,NS,LCONT)=IASSIG(I,KROT,1,NS)
      END DO
      RETURN
      END
C*******************************************************************
C*******************************************************************
       SUBROUTINE GIVENS (NX,NROOTX,NJX,A,B,ROOT,VECT,IAD,IVECT)
       IMPLICIT REAL*8 (A-H,O-Z)
C      GIVENS  -EIGENVALUES AND EIGENVECTORS BY THE GIVENS METHOD.
C      BY FRANKLIN PROSSER, INDIANA UNIVERSITY.
C      SEPTEMBER, 1967
C
C      THANKS ARE DUE TO F. E. HARRIS (STANFORD UNIVERSITY) AND H. H.
C      MICHELS (UNITED AIRCRAFT RESEARCH LABORATORIES) FOR EXCELLENT
C      WORK ON NUMERICAL DIFFICULTIES WITH EARLIER VERSIONS OF THIS
C      PROGRAM.
C
C      THE PARAMETERS FOR THE ROUTINE ARE...
C          NX     ORDER OF MATRIX
C          NROOTX NUMBER OF ROOTS WANTED.  THE NROOTX SMALLEST (MOST
C                 NEGATIVE) ROOTS WILL BE CALCULATED.  IF NO VECTORS
C                 ARE WANTED, MAKE THIS NUMBER NEGATIVE.
C          NJX    ROW DIMENSION OF VECT ARRAY.  SEE 'VECT' BELOW.
C                 NJX MUST NOT BE LESS THAN NX.
C          A      MATRIX STORED BY COLUMNS IN PACKED UPPER TRIANGULAR
C                FORM, I.E. OCCUPYING NX*(NX+1)/2 CONSECUTIVE  LOCATION
C          B      SCRATCH ARRAY USED BY GIVENS.  MUST BE AT LEAST NX*5
C          ROOT   ARRAY TO HOLD THE EIGENVALUES.  MUST BE AT LEAST
C                 NROOTX CELLS LONG.  THE NROOTX SMALLEST ROOTS ARE
C                 ORDERED LARGEST FIRST IN THIS ARRAY.
C          VECT   EIGENVECTOR ARRAY.  EACH COLUMN WILL HOLD AN
C                 EIGENVECTOR FOR THE CORRESPONDING ROOT.  MUST BE
C                 DIMENSIONED WITH 'NJX' ROWS AND AT LEAST 'NROOTX'
C                 COLUMNS, UNLESS NO VECTORS
C                 ARE REQUESTED (NEGATIVE NROOTX).  IN THIS LATTER
C                 CASE, THE ARGUMENT VECT IS JUST A DUMMY, AND THE
C                 STORAGE IS NOT USED.
C                 THE EIGENVECTORS ARE NORMALIZED TO UNITY.
       DIMENSION B(NX,5),A(1),ROOT(NROOTX),VECT(NJX,NROOTX)
       DIMENSION IAD(1)
       COMMON/FILASS/IOUT
C
C     THE ARRAYS A AND B ARE DESTROYED BY THE COMPUTATION.  THE RESULTS
C      APPEAR IN ROOT AND VECT.
C      FOR PROPER FUNCTIONING OF THIS ROUTINE, THE RESULT OF A FLOATING
C      POINT UNDERFLOW SHOULD BE A ZERO.
C      TO CONVERT THIS ROUTINE TO DOUBLE PRECISION (E.G. ON IBM 360
C      MACHINES), BE SURE THAT ALL REAL VARIABLES AND FUNCTION
C      REFERENCES ARE PROPERLY MADE DOUBLE PRECISION.
C     THE VALUE OF 'ETA' (SEE BELOW) SHOULD ALSO BE CHANGED, TO REFLECT
C      THE INCREASED PRECISION.
C
C      THE ORIGINAL REFERENCE TO THE GIVENS TECHNIQUE IS IN OAK RIDGE
C      REPORT NUMBER ORNL 1574 (PHYSICS), BY WALLACE GIVENS.
C      THE METHOD AS PRESENTED IN THIS PROGRAM CONSISTS OF FOUR STEPS,
C      ALL MODIFICATIONS OF THE ORIGINAL METHOD...
C
C      FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE
C      HOUSEHOLDER TECHNIQUE (J. H. WILKINSON, COMP. J. 3, 23 (1960)).
C      THE ROOTS ARE THEN LOCATED BY THE STURM SEQUENCE METHOD (J. M.
C      ORTEGA (SEE REFERENCE BELOW).  THE VECTORS OF THE TRIDIAGONAL
C     FORM ARE THEN EVALUATED (J. H. WILKINSON, COMP. J. 1, 90 (1958)),
C      AND LAST THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE
C      ORIGINAL ARRAY (FIRST REFERENCE).
C      VECTORS FOR DEGENERATE (OR NEAR-DEGENERATE) ROOTS ARE FORCED
C      TO BE ORTHOGONAL, USING A METHOD SUGGESTED BY B. GARBOW, ARGONNE
C      NATIONAL LABS (PRIVATE COMMUNICATION, 1964).  THE GRAM-SCHMIDT
C      PROCESS IS USED FOR THE ORTHOGONALIZATION.
C
C      AN EXCELLENT PRESENTATION OF THE GIVENS TECHNIQUE IS FOUND IN
C      J. M. ORTEGA'S ARTICLE IN 'MATHEMATICS FOR DIGITAL COMPUTERS,'
C      VOLUME 2, ED. BY RALSTON AND WILF, WILEY (1967), PAGE 94.
C
      CONTINUE
C
C ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C **   USERS PLEASE NOTE...
C **   THE FOLLOWING TWO PARAMETERS, ETA AND THETA, SHOULD BE ADJUSTED
C **   BY THE USER FOR HIS PARTICULAR MACHINE.
C **   ETA IS AN INDICATION OF THE PRECISION OF THE FLOATING POINT
C **   REPRESENTATION ON THE COMPUTER BEING USED (ROUGHLY 10**(-M),
C **   WHERE M IS THE NUMBER OF DECIMALS OF PRECISION ).
C **   THETA IS AN INDICATION OF THE RANGE OF NUMBERS THAT CAN BE
C **   EXPRESSED IN THE FLOATING POINT REPRESENTATION (ROUGHLY THE
C **   LARGEST NUMBER).
C **   SOME RECOMMENDED VALUES FOLLOW.
C **   FOR IBM 7094, UNIVAC 1108, ETC. (27-BIT BINARY FRACTION, 8-BIT
C **   BINARY EXPONENT), ETA=1.E-8, THETA=1.E37.
C **   FOR CONTROL DATA 3600 (36-BIT BINARY FRACTION, 11-BIT BINARY
C **   EXPONENT), ETA=1.E-11, THETA=1.E307.
C **   FOR CONTROL DATA 6600 (48-BIT BINARY FRACTION, 11-BIT BINARY
C **   EXPONENT), ETA=1.E-14, THETA=1.E307.
C **   FOR IBM 360/50 AND 360/65 DOUBLE PRECISION (56-BIT HEXADECIMAL
C **   FRACTION, 7-BIT HEXADECIMAL EXPONENT), ETA=1.E-16, THETA=1.E75.
C **
      ETA = 1.0D-16
      THETA = 1.D+64
C** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C  NCK MAY BE PASSED THROUGH THE CALLING SEQUENCE BUT TO BE COMPATABLE
C  WITH PREVIOUS VERSIONS I MAKE A CHOICE HERE.
      NCK=-1
C  NCK.GE.0 GIVES THE EIGENVALUES IN DECREASING ALGEBRAIC ORDER.
C  NCK.LT.0 GIVES THE EIGENVALUES IN INCREASING ALGEBRAIC ORDER.
C  IN EITHER CASE, THE EIGENVECTORS ARE ORDERED SIMILARLY.
       DEL1 = ETA/100.
       DELTA = ETA**2*100.
       SMALL = ETA**2/100.
       DELBIG = THETA*DELTA/1000.
       THETA1 = 1000./THETA
C      TOLER  IS A FACTOR USED TO DETERMINE IF TWO ROOTS ARE CLOSE
C      ENOUGH TO BE CONSIDERED DEGENERATE FOR PURPOSES OF ORTHOGONALI-
C      ZING THEIR VECTORS.  FOR THE MATRIX NORMED TO UNITY, IF THE
C      DIFFERENCE BETWEEN TWO ROOTS IS LESS THAN TOLER, THEN
C      ORTHOGONALIZATION WILL OCCUR.
       TOLER = ETA*100.
C
C      INITIAL VALUE FOR PSEUDORANDOM NUMBER GENERATOR... (2**23)-3
       RPOWER = 8388608.
       RPOW1 = RPOWER/2.
       RAND1 = RPOWER - 3.
C
      N = NX
      NROOT = IABS(NROOTX)
      IF (NROOT.EQ.0) GO TO 1001
      IF (N-1) 1001,1003,105
1003  ROOT(1) = A(1)
      IF (NROOTX.GT.0) VECT(1,1) = 1.0
      GO TO 1001
105   CONTINUE
C     NSIZE    NUMBER OF ELEMENTS IN THE PACKED ARRAY
      NSIZE = (N*(N+1))/2
      NM1 = N-1
      NM2 = N-2
C
C     SCALE MATRIX TO EUCLIDEAN NORM OF 1.  SCALE FACTOR IS ANORM.
      FACTOR = 0.
      DO 70 I=1,NSIZE
70    FACTOR = DMAX1(FACTOR,DABS(A(I)))
      IF (FACTOR.NE.0.) GO TO 72
C     NULL MATRIX.  FIX UP ROOTS AND VECTORS, THEN EXIT.
      DO 78 I=1,NROOT
      IF (NROOTX.LT.0) GO TO 78
      DO 77 J=1,N
77    VECT(J,I) = 0.
      VECT(I,I) = 1.0
78    ROOT(I) = 0.
      GO TO 1001
C
72    ANORM = 0.
      J = 1
      K = 1
86    DO 80 I=1,NSIZE
      IF (I.NE.J) GO TO 81
      ANORM = ANORM + (A(I)/FACTOR)**2/2.
      K = K+1
      J = J+K
      GO TO 80
81    ANORM = ANORM + (A(I)/FACTOR)**2
80    CONTINUE
83    ANORM = DSQRT(ANORM*2.)*FACTOR
      DO 91 I=1,NSIZE
91    A(I) = A(I)/ANORM
      ALIMIT = 1.0
C     WRITE(IOUT,*)'END SETUP'
C     CALL FLUSH(IOUT)
C     CALL CPUTIM('SETUP   ',3)
C
C      TRIDIA SECTION.
C      TRIDIAGONALIZATION OF SYMMETRIC MATRIX
      DO 92 I=1,N
 92   IAD(I)=I*(I+1)/2
       ID = 0
       IA = 1
       IF (NM2.EQ.0) GO TO 201
       DO 200  J=1,NM2
C      J       COUNTS ROW  OF A-MATRIX TO BE DIAGONALIZED
C      IA      START OF NON-CODIAGONAL ELEMENTS IN THE ROW
C      ID      INDEX OF CODIAGONAL ELEMENT ON ROW BEING CODIAGONALIZED.
       IA = IA+J+2
       ID = ID + J + 1
       JP2 = J+2
C      SUM SQUARES OF NON-CODIAGONAL ELEMENTS IN ROW J
       II = IA
       SUM = 0.0
       DO 100 I=JP2,N
C      SUM = SUM + A(II)**2
C100   II = II + I
 100   SUM=SUM+A(IAD(I)-I+J)**2
       TEMP = A(ID)
       IF (SUM.GT.SMALL) GO TO 110
C      NO TRANSFORMATION NECESSARY IF ALL THE NON-CODIAGONAL
C      ELEMENTS ARE TINY.
120    B(J,1) = TEMP
       A(ID) = 0.
       GO TO 200
C      NOW COMPLETE THE SUM OF OFF-DIAGONAL SQUARES
110    SUM = DSQRT(SUM + TEMP**2)
C      NEW CODIAGONAL ELEMENT
       B(J,1) = -SIGN(SUM,TEMP)
C      FIRST NON-ZERO ELEMENT OF THIS W-VECTOR
       B(J+1,2) = DSQRT((1.0 + DABS(TEMP)/SUM)/2.0)
C      FORM REST OF THE W-VECTOR ELEMENTS
       TEMP = SIGN(0.5/(B(J+1,2)*SUM),TEMP)
       II = IA
       DO 130 I=JP2,N
C      B(I,2) = A(II)*TEMP
C130   II = II + I
 130   B(I,2)=A(IAD(I)-I+J)*TEMP
C      FORM P-VECTOR AND SCALAR.  P-VECTOR = A-MATRIX*W-VECTOR.
C      SCALAR = W-VECTOR*P-VECTOR.
       AK = 0.0
C      IC      LOCATION OF NEXT DIAGONAL ELEMENT
       IC = ID + 1
       J1 = J + 1
       DO 190  I=J1,N
       JJ = IC
       TEMP = 0.
C      DO 180  II=J1,N
C      I       RUNS OVER THE NON-ZERO P-ELEMENTS
C      II      RUNS OVER ELEMENTS OF W-VECTOR
C      TEMP = TEMP + B(II,2)*A(JJ)
C      CHANGE INCREMENTING MODE AT THE DIAGONAL ELEMENTS.
C      IF (II.LT.I) GO TO 210
C140   JJ = JJ + II
C      GO TO 180
C210   JJ = JJ + 1
C180   CONTINUE
       JJ=IC-J1
       DO 180 II=J1,I-1
 180   TEMP=TEMP+B(II,2)*A(JJ+II)
       DO 182 II=I,N
 182   TEMP=TEMP+B(II,2)*A(IAD(II)-II+I)
C      BUILD UP THE K-SCALAR (AK)
       AK = AK + TEMP*B(I,2)
       B(I,1) = TEMP
C      MOVE IC TO TOP OF NEXT A-MATRIX 'ROW'
 190   IC = IC + I
C      FORM THE Q-VECTOR
       DO 150  I=J1,N
150    B(I,1) = B(I,1) - AK*B(I,2)
C      TRANSFORM THE REST OF THE A-MATRIX
C      JJ      START-1 OF THE REST OF THE A-MATRIX
       JJ = ID
C      MOVE W-VECTOR INTO THE OLD A-MATRIX LOCATIONS TO SAVE SPACE
C      I       RUNS OVER THE SIGNIFICANT ELEMENTS OF THE W-VECTOR
       DO 160  I=J1,N
       A(JJ) = B(I,2)
       B1=-2*B(I,1)
       B2=-2*B(I,2)
       JO=JJ-J
       IF(I-J1.GT.IVECT) THEN
       CALL VACCUM(I-J,A(JO+J1),B(J1,2),B1)
       CALL VACCUM(I-J,A(JO+J1),B(J1,1),B2)
       ELSE
       DO 170  II=J1,I
C      JJ = JJ + 1
C170   A(JJ) = A(JJ) - 2.0*(B(I,1)*B(II,2) + B(I,2)*B(II,1))
 170   A(JO+II) = A(JO+II) + B1*B(II,2) + B2*B(II,1)
       ENDIF
160    JJ = JJ + I
200    CONTINUE
C      MOVE LAST CODIAGONAL ELEMENT OUT INTO ITS PROPER PLACE
201    CONTINUE
       B(NM1,1) = A(NSIZE-1)
       A(NSIZE-1) = 0.
C      WRITE(IOUT,*)'END TRIDIAG'
C      CALL FLUSH(IOUT)
C      CALL CPUTIM('TRIDIAG ',3)
C
C     STURM SECTION.
C     STURM SEQUENCE ITERATION TO OBTAIN ROOTS OF TRIDIAGONAL FORM.
C     MOVE DIAGONAL ELEMENTS INTO SECOND N ELEMENTS OF B-VECTOR.
C     THIS IS A MORE CONVENIENT INDEXING POSITION.
C     ALSO, PUT SQUARE OF CODIAGONAL ELEMENTS IN THIRD N ELEMENTS.
      JUMP=1
      DO 320 J=1,N
      B(J,2)=A(JUMP)
      B(J,3) = B(J,1)**2
320   JUMP = JUMP+J+1
      DO 310 I=1,NROOT
310   ROOT(I) = +ALIMIT
      ROOTL = -ALIMIT
C    ISOLATE THE ROOTS.  THE NROOT LOWEST ROOTS ARE FOUND, LOWEST FIRST
      DO 330 I=1,NROOT
C     FIND CURRENT 'BEST' UPPER BOUND
      ROOTX = +ALIMIT
      DO 335 J=I,NROOT
335   ROOTX = DMIN1(ROOTX,ROOT(J))
      ROOT(I) = ROOTX
C     GET IMPROVED TRIAL ROOT
500   TRIAL = (ROOTL + ROOT(I))*0.5
C     IF (TRIAL.EQ.ROOTL.OR.TRIAL.EQ.ROOT(I)) GO TO 330
      IF(DABS(TRIAL-ROOTL).LE.ETA)GO TO 330
      IF(DABS(TRIAL-ROOT(I)).LE.ETA)GO TO 330
C     FORM STURM SEQUENCE RATIOS, USING ORTEGA'S ALGORITHM (MODIFIED).
C     NOMTCH IS THE NUMBER OF ROOTS LESS THAN THE TRIAL VALUE.
350   CONTINUE
      NOMTCH = N
      J = 1
360   F0 = B(J,2) - TRIAL
370   CONTINUE
      IF (DABS(F0).LT.THETA1) GO TO 380
      IF (F0.GE.0.) NOMTCH = NOMTCH - 1
      J = J + 1
      IF (J.GT.N) GO TO 390
C     SINCE MATRIX IS NORMED TO UNITY, MAGNITUDE OF B(J,3) IS LESS THAN
C    ONE, AND SINCE F0 IS GREATER THAN THETA1, OVERFLOW IS NOT POSSIBLE
C     AT THE DIVISION STEP
      F0 = B(J,2) - TRIAL - B(J-1,3)/F0
      GO TO 370
380   J = J + 2
      NOMTCH = NOMTCH - 1
      IF (J.LE.N) GO TO 360
390   CONTINUE
C     FIX NEW BOUNDS ON ROOTS
      IF (NOMTCH.GE.I) GO TO 540
      ROOTL = TRIAL
      GO TO 500
540   ROOT(I) = TRIAL
      NOM = MIN0(NROOT,NOMTCH)
      ROOT(NOM) = TRIAL
      GO TO 500
330   CONTINUE
C          SKIP THIS SECTION IF NCK IS LT 0
      IF(NCK.LT.0) GO TO 807
C     REVERSE THE ORDER OF THE EIGENVALUES, SINCE CUSTOM DICTATES
C     'LARGEST FIRST'.  THIS SECTION MAY BE REMOVED IF DESIRED WITHOUT
C     AFFECTING THE REMAINDER OF THE ROUTINE.
      NRT = NROOT/2
      DO 10 I=1,NRT
      SAVE = ROOT(I)
      NMIP1 = NROOT - I + 1
      ROOT(I) = ROOT(NMIP1)
10    ROOT(NMIP1) = SAVE
C
C     TRIVEC SECTION.
C     EIGENVECTORS OF CODIAGONAL FORM
807   CONTINUE
C     WRITE(IOUT,*)'END STURM'
C     CALL FLUSH(IOUT)
C     CALL CPUTIM('STURM   ',3)
C     QUIT NOW IF NO VECTORS WERE REQUESTED.
      IF (NROOTX.LT.0) GO TO 1002
C     INITIALIZE VECTOR ARRAY.
      DO 15 I=1,N
      DO 15 J=1,NROOT
15    VECT(I,J) = 1.0
      DO 700 I=1,NROOT
      AROOT = ROOT(I)
C     ORTHOGONALIZE IF ROOTS ARE CLOSE.
      IF (I.EQ.1) GO TO 710
C     THE ABSOLUTE VALUE IN THE NEXT TEST IS TO ASSURE THAT THE TRIVEC
C     SECTION IS INDEPENDENT OF THE ORDER OF THE EIGENVALUES.
715   IF (DABS(ROOT(I-1)-AROOT).LT.TOLER) GO TO 720
710   IA = -1
720   IA = IA + 1
      ELIM1 = A(1) - AROOT
      ELIM2 = B(1,1)
      JUMP = 1
      DO 750  J=1,NM1
      JUMP = JUMP+J+1
C     GET THE CORRECT PIVOT EQUATION FOR THIS STEP.
      IF (DABS(ELIM1).LE.DABS(B(J,1))) GO TO 760
C     FIRST (ELIM1) EQUATION IS THE PIVOT THIS TIME.  CASE 1.
      B(J,2) = ELIM1
      B(J,3) = ELIM2
      B(J,4) = 0.
      TEMP = B(J,1)/ELIM1
      ELIM1 = A(JUMP) - AROOT - TEMP*ELIM2
      ELIM2 = B(J+1,1)
      GO TO 755
C     SECOND EQUATION IS THE PIVOT THIS TIME.  CASE 2.
760   B(J,2) = B(J,1)
      B(J,3) = A(JUMP) - AROOT
      B(J,4) = B(J+1,1)
      TEMP = 1.0
      IF (DABS(B(J,1)).GT.THETA1) TEMP = ELIM1/B(J,1)
      ELIM1 = ELIM2 - TEMP*B(J,3)
      ELIM2 = -TEMP*B(J+1,1)
C     SAVE FACTOR FOR THE SECOND ITERATION.
755   B(J,5) = TEMP
750   CONTINUE
      B(N,2) = ELIM1
      B(N,3) = 0.
      B(N,4) = 0.
      B(NM1,4) = 0.
      ITER = 1
      IF (IA.NE.0) GO TO 801
C     BACK SUBSTITUTE TO GET THIS VECTOR.
790   L = N + 1
      DO 780 J=1,N
      L = L - 1
786   CONTINUE
      ELIM1=VECT(L,I)
      IF(L.LT.N)ELIM1=ELIM1-VECT(L+1,I)*B(L,3)
      IF(L.LT.N-1)ELIM1=ELIM1-VECT(L+2,I)*B(L,4)
C     IF OVERFLOW IS CONCEIVABLE, SCALE THE VECTOR DOWN.
C     THIS APPROACH IS USED TO AVOID MACHINE-DEPENDENT AND SYSTEM-
C     DEPENDENT CALLS TO OVERFLOW ROUTINES.
      IF (DABS(ELIM1).GT.DELBIG) GO TO 782
      TEMP = B(L,2)
      IF (DABS(B(L,2)).LT.DELTA) TEMP = DELTA
      VECT(L,I) = ELIM1/TEMP
      GO TO 780
C     VECTOR IS TOO BIG.  SCALE IT DOWN.
782   DO 784 K=1,N
784   VECT(K,I) = VECT(K,I)/DELBIG
      GO TO 786
780   CONTINUE
      GO TO (820,800), ITER
C     SECOND ITERATION.  (BOTH ITERATIONS FOR REPEATED-ROOT VECTORS).
820   ITER = ITER + 1
890   ELIM1 = VECT(1,I)
      DO 830 J=1,NM1
      IF (B(J,2).EQ.B(J,1)) GO TO 840
C     IF(DABS(B(J,2)-B(J,1)).LE.ETA)GO TO 840
C     CASE ONE.
850   VECT(J,I) = ELIM1
      ELIM1 = VECT(J+1,I) - ELIM1*B(J,5)
      GO TO 830
C     CASE TWO.
840   VECT(J,I) = VECT(J+1,I)
      ELIM1 = ELIM1 - VECT(J+1,I)*TEMP
830   CONTINUE
      VECT(N,I) = ELIM1
      GO TO 790
C     PRODUCE A RANDOM VECTOR
801   CONTINUE
      DO 802 J=1,N
C    GENERATE PSEUDORANDOM NUMBERS WITH UNIFORM DISTRIBUTION IN (-1,1).
C     THIS RANDOM NUMBER SCHEME IS OF THE FORM...
C     RAND1 = DMOD((2**12+3)*RAND1,2**23)
C     IT HAS A PERIOD OF 2**21 NUMBERS.
      RAND1 = DMOD(4099.*RAND1,RPOWER)
802   VECT(J,I) = RAND1/RPOW1 - 1.
      GO TO 790
C
C     ORTHOGONALIZE THIS REPEATED-ROOT VECTOR TO OTHERS WITH THIS ROOT.
800   IF (IA.EQ.0) GO TO 885
      DO 860 J1=1,IA
      K = I - J1
      TEMP = 0.
      DO 870 J=1,N
870   TEMP = TEMP + VECT(J,I)*VECT(J,K)
      DO 880 J=1,N
880   VECT(J,I) = VECT(J,I) - TEMP*VECT(J,K)
860   CONTINUE
885   GO TO (890,900), ITER
C     NORMALIZE THE VECTOR
900   ELIM1 = 0.
      DO 904 J=1,N
904   ELIM1 = DMAX1(DABS(VECT(J,I)),ELIM1)
      TEMP=0.
      DO 910 J=1,N
      ELIM2=VECT(J,I)/ELIM1
910   TEMP = TEMP + ELIM2**2
      TEMP=1.0/(DSQRT(TEMP)*ELIM1)
      DO 920 J=1,N
      VECT(J,I) = VECT(J,I)*TEMP
      IF (DABS(VECT(J,I)).LT.DEL1) VECT(J,I) = 0.
920   CONTINUE
700   CONTINUE
C     WRITE(IOUT,*)'END TRIVEC'
C     CALL FLUSH(IOUT)
C     CALL CPUTIM('TRIVEC  ',3)
C
C      SIMVEC SECTION.
C      ROTATE CODIAGONAL VECTORS INTO VECTORS OF ORIGINAL ARRAY
C      LOOP OVER ALL THE TRANSFORMATION VECTORS
       IF (NM2.EQ.0) GO TO 1002
       JUMP = NSIZE - (N+1)
       IM = NM1
       DO 950  I=1,NM2
       J1 = JUMP
C      MOVE A TRANSFORMATION VECTOR OUT INTO BETTER INDEXING POSITION.
       DO 955  J=IM,N
       B(J,2) = A(J1)
955    J1 = J1 + J
C      MODIFY ALL REQUESTED VECTORS.
       DO 960  K=1,NROOT
       TEMP = 0.
C      FORM SCALAR PRODUCT OF TRANSFORMATION VECTOR WITH EIGENVECTOR
       DO 970  J=IM,N
970    TEMP = TEMP + B(J,2)*VECT(J,K)
       TEMP = TEMP + TEMP
       DO 980  J=IM,N
980    VECT(J,K) = VECT(J,K) - TEMP*B(J,2)
960    CONTINUE
       JUMP = JUMP - IM
950    IM = IM - 1
1002   CONTINUE
C      RESTORE ROOTS TO THEIR PROPER SIZE.
       DO 95 I=1,NROOT
95     ROOT(I) = ROOT(I)*ANORM
1001   CONTINUE
C     WRITE(IOUT,*)'END SIMVEC'
C     CALL FLUSH(IOUT)
C      CALL CPUTIM('SIMVEC  ',3)
       RETURN
       END
C*****************************************************************
C*****************************************************************
      SUBROUTINE VACCUM(L,A,B,S)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(L),B(L)
      DO 1 I=1,L
 1    A(I)=A(I)+B(I)*S
      RETURN
      END
C*****************************************************************
C*****************************************************************
      SUBROUTINE TRID(N,Z,NN,E,D)
      IMPLICIT NONE
      INTEGER J,J1,K,IOUT,IPRINT,JPRINT,N,NN,L,I
      DOUBLE PRECISION H,E,D,HH,Z,G,F,TOL,EPS
      DIMENSION E(N),D(N)
      DIMENSION Z(NN,NN)
      COMMON/FILASS/IOUT
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/TOLS/TOL,EPS
      IF((N-1).GT.0) THEN
        DO 111 I=N,2,-1
          L=I-2
          F=Z(I,I-1)
          G=0.D0
          IF((L-1).GE.0) THEN
            DO 3 K=1,L
3             G=G+Z(I,K)*Z(I,K)
          ENDIF
          H=G+F*F
          IF(G.LE.TOL) THEN
            E(I)=F
            H=0.D0
            D(I)=H
          ELSE
            L=L+1
            G=DSQRT(H)
            IF(F.GE.0.D0) G=-G
            E(I)=G
            H=H-F*G
            Z(I,I-1)=F-G
            F=0.D0
            IF((L-1).GE.0) THEN
              DO 10 J=1,L
                Z(J,I)=Z(I,J)/H
                G=0.D0
                DO 8 K=1,J
8                 G=G+Z(J,K)*Z(I,K)
                J1=J+1
                IF((J1-L).LE.0) THEN
                  DO 9 K=J1,L
9                   G=G+Z(K,J)*Z(I,K)
                ENDIF
                E(J)=G/H
                F=F+G*Z(J,I)
10            CONTINUE
            ENDIF
            HH=F/(H+H)
            IF((L-1).GE.0) THEN
              DO 12 J=1,L
                F=Z(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 11 K=1,J
11                Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
12            CONTINUE
              D(I)=H
            ENDIF
          ENDIF
111     CONTINUE
      ENDIF
      D(1)=0.D0
      E(1)=0.D0
      DO 22 I=1,N
        L=I-1
        IF(D(I).NE.0.D0) THEN
          DO 18 J=1,L
            G=0.D0
            DO 16 K=1,L
16            G=G+Z(I,K)*Z(K,J)
            DO 17 K=1,L
17            Z(K,J)=Z(K,J)-G*Z(K,I)
18        CONTINUE
        ENDIF
        D(I)=Z(I,I)
        Z(I,I)=1.D0
        IF((L-1).GE.0) THEN
          DO 21 J=1,L
            Z(I,J)=0.D0
21          Z(J,I)=0.D0
        ENDIF
22    CONTINUE
      IF((N-2).GE.0) THEN
        DO 102 I=2,N
102       E(I-1)=E(I)
      ENDIF
      E(N)=0.D0
      CALL TDIAG(N,Z,NN,E,D)
      RETURN
      END
C*******************************************************************
C*****************************************************************
      SUBROUTINE TDIAG(N,Z,NN,E,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(N),D(N)
      DIMENSION Z(NN,NN)
      COMMON/FILASS/IOUT
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/TOLS/TOL,EPS
      B=0.D0
      F=0.D0
      DO 22 L=1,N
      J=0
      H=EPS*(ABS(D(L))+ABS(E(L)))
      IF(B-H)4,5,5
4     B=H
5     DO 6 M=L,N
      IF(ABS(E(M))-B)7,7,6
6     CONTINUE
7     IF(M-L)8,21,8
8     IF(J-30)9,23,9
9     J=J+1
      P=(D(L+1)-D(L))/(2.D0*E(L))
      R=SQRT(P*P+1)
      IF(P-0.D0)10,11,11
10    V=P-R
      GO TO 12
11    V=P+R
12    H=D(L)-E(L)/V
      DO 13 I=L,N
13    D(I)=D(I)-H
      F=F+H
      P=D(M)
      C=1.D0
      S=0.D0
      IF(L-M)14,20,20
14    I=M-1
25    G=C*E(I)
      H=C*P
      IF(ABS(P)-ABS(E(I)))16,15,15
15    C=E(I)/P
      R=SQRT(C*C+1.D0)
      E(I+1)=S*P*R
      S=C/R
      C=1.D0/R
      GO TO 17
16    C=P/E(I)
      R=SQRT(C*C+1)
      E(I+1)=S*R*E(I)
      S=1.D0/R
      C=C/R
17    P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
      DO 18 K=1,N
      H=Z(K,I+1)
      Z(K,I+1)=S*Z(K,I)+C*H
18    Z(K,I)=C*Z(K,I)-S*H
      IF(I-L)20,20,19
19    I=I-1
      GO TO 25
20    E(L)=S*P
      D(L)=C*P
      IF(ABS(E(L))-B)21,21,8
21    D(L)=D(L)+F
22    CONTINUE
      GO TO 24
23    WRITE(IOUT,1000)
1000  FORMAT(53H0DIAGONALIZATION HAS FAILED - MORE THAN 30 ITERATIONS)
      STOP 'TDIAG'
24    RETURN
      END
