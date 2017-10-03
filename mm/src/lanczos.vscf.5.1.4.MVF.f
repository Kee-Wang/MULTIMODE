C****************************************************************
C****************************************************************
C**LANCZOS
C****************************************************************
C****************************************************************
      SUBROUTINE LANCZB(X,V,Z,W,Y,S,ISIZE,NVAL,ISTAT,NSTAT,NMODE,OV,
     1VEC,IASSIG,ISIZMX,J21,NS,ENERGY,ELAST,NVEC,XA,XACOPY,XK,WK,SUP4,
     2EVAL,NCYCLE,NVS,LANCNT,LANK,LANCYC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRIAT,LGIV
C**STORE X,V,Z ON DISC
      DIMENSION LANCNT(1),LANK(1)
      DIMENSION X(ISIZE,NVAL),V(ISIZE,LANCYC),Z(ISIZE,LANCYC),
     2W(ISIZE),Y(ISIZE),S(ISIZE,LANCYC),OV(ISIZE,1),VEC(ISIZE,1)
      DIMENSION ENERGY(NVAL),ELAST(NVAL),NVEC(NVAL)
      DIMENSION XK(LANCYC*NCYCLE,LANCYC),XA(1),XACOPY(1)
      DIMENSION SUP4(5*LANCYC*NCYCLE),EVAL(LANCYC*NCYCLE)
      DIMENSION WK(LANCYC*NCYCLE)
      DIMENSION ISTAT(NSTAT,NMODE),IASSIG(ISIZMX,J21,3,1)
      DIMENSION JJMAX(3),XMAX(4)
      DIMENSION NS1(3),NS2(3),IOFF1(3),IOFF2(3)
      COMMON/FILASS/IOUT,INP
      COMMON/TRIATO/TRIAT
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/EVL/EVL,CUT
      COMMON/GIVEN/LGIV
      COMMON/LANTOL/TOLLAN
      COMMON/CIDIAG/ICID
      COMMON/ROTS/JMAX,KMAX,IDUM,KEL21,KEL
      COMMON/REACTL/JREACT
      COMMON/MATRIX/NVALV,NVALR,KSTEP
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LAN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TOLCUT,EVLCUT,TOLMAT,LANCCC
      COMMON/JKAKC/JTHIS,KA,KC
C**TEMPORARY (DIMENSIONS)
      COMMON/KAKCAS/IASKA(101,10),IASKC(101,10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
C****************************************************************
100   FORMAT(/,1X,'INITIAL ENERGY = ',D20.12)
105   FORMAT(/,1X,'LANCZOS NOT CONVERGED FOR FUNCTIONS:')
110   FORMAT(/,1X,'LANCZOS FULLY CONVERGED FOR FUNCTIONS ',I4,' -',I4,
     1' IN ',I3,' CYCLES')
115   FORMAT(20I4)
C****************************************************************
      IF(LLAN20.LT.0)THEN
        REWIND 51
        READ(51)Y
      END IF
C**STORE V ON 52
      REWIND 52
      NCONV=0
      LANCON=0
C**STORE X ON 54
      REWIND 54


9999  CONTINUE
C**STORE Z ON 53
      REWIND 53
C**************************************************************
C**************************************************************
C**INITIAL GUESS !!!
C**************************************************************
C**************************************************************
      ILAN=0
      DO K=1,LANCYC+LANCON
        DO I=1,ISIZE
          X(I,K)=0
        END DO
        XLAST=1.D+10
        DO I=1,ISIZE
          IF(Y(I).LT.XLAST)THEN
            IF(K.EQ.1)THEN
              XLAST=Y(I)
              NVEC(K)=I
            ELSE
              IGOT=0
              DO J=1,K-1
                IF(I.EQ.NVEC(J))IGOT=1
              END DO
              IF(IGOT.EQ.0)THEN
                XLAST=Y(I)
                NVEC(K)=I
              END IF
            END IF
          END IF
        END DO
C**INITIAL GUESS !!!
        NV=NVEC(K)
        X(NV,K)=1
      END DO
      DO K=1,LANCYC
        K1=LANCON+K
        NVEC(K)=NVEC(K1)
        DO I=1,ISIZE
          X(I,K)=X(I,K1)
        END DO
      END DO
      
       CALL TIMIT(3)
       CALL FLUSH(IOUT)
       WRITE(IOUT,*)
       WRITE(IOUT,*)'DIAGONAL ELEMENTS'
       DO K=1,LANCYC
         WRITE(IOUT,*)NVEC(K),Y(NVEC(K))*WAVENM
       END DO
       CALL TIMIT(1)
       CALL FLUSH(IOUT)

C**************************************************************
C**************************************************************
C**********INITIAL SET-UP.....M(rs) = <psi(r)/H/psi(s)>
C**ASSUME INITIAL FUNCTIONS X0 ORTHONORMAL (VIRTUAL SCF STATES)
C**SET INITIAL FUNCTION V0 = X0 (NORMALISED)
C**************************************************************
C**************************************************************
      DO K=1,LANCYC
        DO I=1,ISIZE
          V(I,K)=X(I,K)
        END DO
        CALL SCHMDB(V,K,Z,LANCYC,ISIZE,NCONV)
      END DO
C*******************************DISC (52)
      DO KV=1,LANCYC
        CALL LOUT(V(1,KV),ISIZE,52)
      END DO
C**FORM Z0 = M.V0 (UNNORMALISED)
C**ILAN=0 FIRST TIME....NVEC TEST NOT APPLICABLE
      IF(JCOUPL.GE.0)THEN
        CALL MXBP(V,Z,W,W,ISIZE,NVEC,LANCYC,ILAN,LANCNT,LANK)
      ELSE
        CALL MXBM(V,Z,W,W,ISIZE,NVEC,LANCYC,ILAN,LANCNT,LANK)
      END IF
      DO K=1,LANCYC
C**FORM S0 = Z0
        DO I=1,ISIZE
          S(I,K)=Z(I,K)
        END DO
      END DO
C*******************************DISC (53)
      DO KV=1,LANCYC
        CALL LOUT(Z(1,KV),ISIZE,53)
      END DO
C**FORM V0.(M.V0) = V0.Z0
      DO KV=1,LANCYC
        J0=KV*(KV-1)/2
        DO KZ=1,KV
          CALL DOT(V(1,KV),Z(1,KZ),E,ISIZE)
          XA(J0+KZ)=E
          XACOPY(J0+KZ)=E
        END DO
C**V0.(M.V0) = V0.Z0 = <V0/H/V0> = E-VALUE OF V0
        ENERGY(KV+LANCON)=E*WAVENM
        ELAST(KV)=ENERGY(KV+LANCON)
      END DO
C**************************************************************
C**************************************************************
C**FORM NEXT BASIS FUNCTION V1 (UNNORMALISED)
C**************************************************************
C**************************************************************
      DO KV=1,LANCYC
        J0=KV*(KV-1)/2
C**NVEC POINTS TO POSITION OF UNIT VECTOR FIRST TIME
        NV=NVEC(KV)
        V(NV,KV)=0
        DO I=1,ISIZE
          IF(I.NE.NV)THEN
            V(I,KV)=-(S(I,KV)-XA(J0+KV)*X(I,KV))/(Y(I)-XA(J0+KV))
          END IF
        END DO
C**ORTHONORMALISE V1 TO V0
        CALL SCHMDB(V,KV,Z,LANCYC,ISIZE,1+NCONV)
      END DO
C**SCHMIDT POSITIONS DISC 52 FOR NEW V(K)
C*******************************DISC (52)
      DO KV=1,LANCYC
        CALL LOUT(V(1,KV),ISIZE,52)
      END DO
      KGOT=0
      DO K=1,LANCYC
C**NVEC USED AS DUMMY INDICATOR FIRST TIME
        NVEC(K)=100000
        IGOT=0
        DO I=1,ISIZE
          IF(V(I,K).NE.0)IGOT=1
        END DO
        IF(IGOT.EQ.0)THEN
          ELAST(K)=0.D0
          NVEC(K)=2
        END IF
        IF(ELAST(K).NE.0.D0)KGOT=1
      END DO
      IF(KGOT.EQ.0)GO TO 2000
C**TO BE SUCCESSFUL, MUST BE AT LEAST ONE V(I,K) NON-ZERO
      KSIZE=LANCYC
C**************************************************************
C**************************************************************
C**START LANCZOS ITERATIONS
C**************************************************************
C**************************************************************
      CALL TIMIT(3)
      CALL FLUSH(IOUT)
      DO ILAN=2,NCYCLE
        WRITE(IOUT,*)
        WRITE(IOUT,*)'Cycle ',ILAN
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
        ILOFF=0
C*******************************DISC (53)
        REWIND 53
        DO I=1,ILAN-1
          DO KV=1,LANCYC
            CALL LIN(Z(1,KV),ISIZE,53)
          END DO
          KVV=0
          DO KV=1,LANCYC
            IF(NVEC(KV).GT.ILAN)THEN
              KVV=KVV+1
              IROFF=KSIZE+KVV
              J0=ILOFF+IROFF*(IROFF-1)/2
              KZZ=0
              DO KZ=1,LANCYC
                IF(NVEC(KZ).GT.I)THEN
                  KZZ=KZZ+1
                  CALL DOT(V(1,KV),Z(1,KZ),P,ISIZE)
                  XA(J0+KZZ)=P
                  XACOPY(J0+KZZ)=P
                END IF
              END DO
            END IF
          END DO
          ILOFF=ILOFF+LANCYC
          DO K=1,LANCYC
            IF(NVEC(K).LE.I)ILOFF=ILOFF-1
          END DO
        END DO
C**FORM Z = M.V
        IF(JCOUPL.GE.0)THEN
          CALL MXBP(V,Z,W,W,ISIZE,NVEC,LANCYC,ILAN,LANCNT,LANK)
        ELSE
          CALL MXBM(V,Z,W,W,ISIZE,NVEC,LANCYC,ILAN,LANCNT,LANK)
        END IF
C*******************************DISC (53)
        DO KV=1,LANCYC
          CALL LOUT(Z(1,KV),ISIZE,53)
          CALL FLUSH(53)
        END DO
        KVV=0
        DO KV=1,LANCYC
          IF(NVEC(KV).GT.ILAN)THEN
            KVV=KVV+1
            IROFF=KSIZE+KVV
            J0=KSIZE+IROFF*(IROFF-1)/2
            KZZ=0
            DO KZ=1,KV
              IF(NVEC(KZ).GT.ILAN)THEN
                KZZ=KZZ+1
                CALL DOT(V(1,KV),Z(1,KZ),P,ISIZE)
                XA(J0+KZZ)=P
                XACOPY(J0+KZZ)=P
              END IF
            END DO
          END IF
        END DO
C**GET CURRENT MATRIX SIZE FROM PREVIOUS SIZE
        KSIZE=KSIZE+LANCYC
        DO K=1,LANCYC
          IF(NVEC(K).LE.ILAN)KSIZE=KSIZE-1
        END DO
C**SET UP SQUARE MATRIX FOR QL
        IF(.NOT.LGIV)THEN
          CALL QLCOPY(XACOPY,KSIZE,XK,LANCYC*NCYCLE)
        END IF
C**DIAGONALISE MATRIX ORDER KSIZE (NO OVERLAP)
        CALL DIAG(XA,XK,LANCYC*NCYCLE,KSIZE,-1,SUP4,EVAL,WK,
     1  LANCYC,LANCYC,IP,ISIZMX,NMODE,XK,XK,IASSIG,ISIZMX,J21,NS)
C**************************************************************
C**************************************************************
C**FORM NEW FUNCTION X.....LOWEST-ENERGY EIGENFUNCTIONS K (NORMALISED)
C**************************************************************
C**************************************************************
        DO K=1,LANCYC
          DO I=1,ISIZE
            X(I,K)=0
          END DO
        END DO
        JOFF=0
C*******************************DISC (52)
        REWIND 52
C**SKIP CONVERGED FUNCTIONS
        DO J=1,NCONV
          DO KV=1,LANCYC
            CALL LIN(V(1,KV),ISIZE,52)
          END DO
        END DO
        DO J=1,ILAN
C*******************************DISC (52)
          DO KV=1,LANCYC
            CALL LIN(V(1,KV),ISIZE,52)
          END DO
          DO K=1,LANCYC
              KOFF=JOFF
              DO N=1,LANCYC
                IF(NVEC(N).GT.J)THEN
                  KOFF=KOFF+1
                  DO I=1,ISIZE
                    X(I,K)=X(I,K)+XK(KOFF,K)*V(I,N)
                  END DO
                END IF
              END DO
          END DO
          JOFF=KOFF
        END DO
        DO K=1,LANCYC
          IF(ELAST(K).EQ.0.D0)THEN
            DO I=1,ISIZE
              Z(I,K)=X(I,K)
              X(I,K)=Z(I,K)
            END DO
          END IF
        END DO
C**FORM NEW FUNCTION S
        DO K=1,LANCYC
          ENERGY(K+LANCON)=EVAL(K)*WAVENM
          DO I=1,ISIZE
            S(I,K)=0
          END DO
        END DO
        JOFF=0
C*******************************DISC (53)
        REWIND 53
        DO J=1,ILAN
          DO KV=1,LANCYC
            CALL LIN(Z(1,KV),ISIZE,53)
          END DO
          DO K=1,LANCYC
            KOFF=JOFF
            DO N=1,LANCYC
              IF(NVEC(N).GT.J)THEN
                KOFF=KOFF+1
                DO I=1,ISIZE
                  S(I,K)=S(I,K)+XK(KOFF,K)*Z(I,N)
                END DO
              END IF
            END DO
          END DO
          JOFF=KOFF
        END DO
        IF(ILAN.EQ.NCYCLE)GO TO 1000
        KGOT=0
        DO K=1,LANCYC
          IF(ELAST(K).NE.0.D0)THEN
            IF(DABS(ENERGY(K+LANCON)-ELAST(K)).LT.TOLLAN)THEN
              ELAST(K)=0.D0
C**FUNCTION K CONVERGED EARLY
C**NVEC(K) IS FIRST CYCLE IN WHICH FUNCTION 'K' IS NOT INCLUDED
              NVEC(K)=ILAN+1
            ELSE
              ELAST(K)=ENERGY(K+LANCON)
            END IF
          END IF
          IF(ELAST(K).NE.0.D0)KGOT=1
        END DO
        IF(KGOT.EQ.0)GO TO 2000
C**************************************************************
C**************************************************************
C**FORM NEXT BASIS FUNCTION V(ILAN+1) (UNNORMALISED)
C**************************************************************
C**************************************************************
        KGOT=0
        DO K=1,LANCYC
          IF(ELAST(K).NE.0.D0)THEN
            DO I=1,ISIZE
              V(I,K)=-(S(I,K)-EVAL(K)*X(I,K))/
     1        (Y(I)-EVAL(K))
            END DO
C**ORTHONORMALISE V(ILAN+1) TO V(ILAN)
            CALL SCHMDB(V,K,Z,LANCYC,ISIZE,ILAN+NCONV)
            IGOT=0
            DO I=1,ISIZE
              IF(V(I,K).NE.0)IGOT=1
            END DO
            IF(IGOT.EQ.0)THEN
              ELAST(K)=0.D0
C**ZERO FUNCTION K 
C**NVEC(K) IS FIRST CYCLE IN WHICH FUNCTION 'K' IS NOT INCLUDED
              NVEC(K)=ILAN+1
            END IF
          ELSE
            DO I=1,ISIZE
              V(I,K)=0
            END DO
          END IF
          IF(ELAST(K).NE.0.D0)KGOT=1
        END DO
        IF(KGOT.EQ.0)GO TO 2000
C**SCHMIDT POSITIONS DISC 52 FOR NEW V(K)
C**OTHERWISE DISC 52 ALREADY POSITIONED FORMING NEW X
C*******************************DISC (52)
        DO KV=1,LANCYC
          CALL LOUT(V(1,KV),ISIZE,52)
        END DO
C**RE-LOAD MATRIX ORDER ILAN IN PREPARATION FOR MATRIX ORDER ILAN+1
        DO I=1,KSIZE*(KSIZE+1)/2
          XA(I)=XACOPY(I)
        END DO
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
      END DO
      GO TO 1002
C**************************************************************
C**************************************************************
C**END LANCZOS ITERATIONS
C**************************************************************
C**************************************************************
1000  CONTINUE
      CALL TIMIT(3)
      CALL FLUSH(IOUT)
1002  CONTINUE
      WRITE(IOUT,105)
      CALL FLUSH(IOUT)
      N=0
      DO K=1,LANCYC
        IF(ELAST(K).NE.0)THEN
          N=N+1
          NVEC(N)=K
        END IF
      END DO
      WRITE(IOUT,115)(NVEC(K),K=1,N)
      CALL FLUSH(IOUT)
      GO TO 3000
2000  CONTINUE
C**CONVERGED
      WRITE(IOUT,110)LANCON+1,LANCON+LANCYC,ILAN
      CALL FLUSH(IOUT)
3000  CONTINUE
C
C**WRITE THIS BLOCK TO (54)
      DO KV=1,LANCYC
        CALL LOUT(X(1,KV),ISIZE,54)
      END DO
      NCONV=NCONV+1
      REWIND 54
C**READ ALL BLOCKS FROM (54)
      DO KV=1,NCONV*LANCYC
        CALL LIN(X(1,KV),ISIZE,54)
      END DO
      REWIND 52
C**TRANSFER TO (52)
      DO KV=1,NCONV*LANCYC
        CALL LOUT(X(1,KV),ISIZE,52)
      END DO
      LANCON=LANCON+LANCYC
      IF(LANCON.LT.NVAL)GO TO 9999
C**************************************************************
C**************************************************************
C**ASSIGN LEVELS
C**************************************************************
C**************************************************************
      WRITE(IOUT,*)'Start Assignments'
      CALL TIMIT(1)
      CALL FLUSH(IOUT)
      MT=1
      IF(EVL.NE.0.0D0)GO TO 4
      MT=2
4     CONTINUE
      IF(EVL.EQ.0)EVL=ENERGY(1)
      DO 7 I=1,NVAL
      ELAST(I)=ENERGY(I)/WAVENM
      ENERGY(I)=ENERGY(I)-EVL
7     CONTINUE
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
        NSIZE=ISIZE
        KKC=KC
        KKA=KA
      ELSE
        NSIZE=NVALV
      END IF
      ICYCL=1
      IF(JPRINT.LT.-1.AND.ICID.EQ.0)ICYCL=2
      IF(JPRINT.LT.-2.AND.ICID.EQ.0)ICYCL=3
      DO 10 I=1,NVAL
      IF(ENERGY(I).GT.CUT)GO TO 1001
      DO 101 ICYC=1,ICYCL
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
        DO 9 J=1,ISIZE
        IF(X(J,I).EQ.XMAX(ICYC))GO TO 9
        IF(ABS(X(J,I)).LE.ABS(XMAX(ICYC+1)))GO TO 9
        XMAX(ICYC+1)=X(J,I)
        IMAX=J
9       CONTINUE
      else
c  get third biggest coeff
        xmax(4)=0.0
        do 90 j=1,isize
        if((x(j,i).eq.xmax(2)).or.(x(j,i).eq.xmax(3))) go to 90
        if (abs(x(j,i)).lt.abs(xmax(4))) go to 90
        xmax(4)=x(j,i)
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
C       KKC=KKC+MOD(IIA+1,2)
        IF(MOD(KKA,2).EQ.MOD(IIA+1,2))KKC=KKC+1
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
        KKA=IASKA(KKEL,NSX)
        KKC=IASKC(KKEL,NSX)
        JJMAX(1)=IOFF+IASSIG(IMAX,KKEL,1,NSX)
        IF(JPRINT.LT.-1)JJMAX(2)=IOFF+IASSIG(IMAX,KKEL,2,NSX)
        IF(JPRINT.LT.-2)JJMAX(3)=IOFF+IASSIG(IMAX,KKEL,3,NSX)
      END IF
101   CONTINUE
      CALL PRCI(MT,ISTAT,NSTAT,NMODE,JJMAX(1),JJMAX(2),JJMAX(3),
     1XMAX(2),XMAX(3),XMAX(4),JX,ENERGY(I),EVL,JTHIS,KKA,KKC,45,
     2OV,VEC,ICSIZ1,ICSIZ2,IOFF1,IOFF2,IASSIG,ISIZMX,J21,NVS,NS1,NS2)
      CALL FLUSH(IOUT)
1001  CONTINUE
10    CONTINUE
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE SCHMDB(V,KV,WRK,NVAL,ISIZE,NV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(ISIZE,NVAL),WRK(ISIZE,NVAL)
      COMMON/FILASS/IOUT,INP
C**SCHMIDT ORTHOGONALISE
      N=NV+1
C**NEW FUNCTION IS V(ISIZE,KV,N)
C**NORMALISE NEW FUNCTION
      CALL DOT(V(1,KV),V(1,KV),S,ISIZE)
      SQ=SQRT(S)
C**BEWARE!!! IS 10**-6 CORRECT?
C**WE WILL MAKE THE PROGRAM PRINT...
C**EACH TIME THIS IS OBEYED
      IF(SQ.GT.1.D-6)THEN
        S=1/SQ
      ELSE
        WRITE(IOUT,*)'ZERO VECTOR NO. ',KV
        CALL FLUSH(IOUT)
        S=0
C**POSITION DISC 52 FOR V
        REWIND 52
        DO I=1,N
C************************DISC (52)
          IF(I.NE.N)THEN
            DO KW=1,NVAL
              CALL LIN(WRK(1,KW),ISIZE,52)
            END DO
          END IF
        END DO
      END IF
      DO I=1,ISIZE
        V(I,KV)=V(I,KV)*S
      END DO
      IF(S.EQ.0)RETURN
C**IF I<N, OLD FUNCTIONS ARE V(ISIZE,NVAL,NV)
C**IF I=N, OLD FUNCTIONS ALSO INCLUDE V(ISIZE,KV-1,N)
C************************DISC (52)
      REWIND 52
      DO I=1,N
        IF(I.NE.N)THEN
          DO KW=1,NVAL
            CALL LIN(WRK(1,KW),ISIZE,52)
          END DO
        END IF
        NXVAL=NVAL
        IF(I.EQ.N)NXVAL=KV
        DO K=1,NXVAL
          IF(I.EQ.N)THEN
            DO J=1,ISIZE
              WRK(J,K)=V(J,K)
            END DO
          END IF
          CALL DOT(V(1,KV),WRK(1,K),S,ISIZE)
          IF(I.EQ.N.AND.K.EQ.KV)THEN
            SQ=SQRT(S)
C**BEWARE!!! IS 10**-6 CORRECT?
C**WE WILL MAKE THE PROGRAM PRINT...
C**EACH TIME THIS IS OBEYED
            IF(SQ.GT.1.D-6)THEN
              S=1.D0-1.D0/SQ
            ELSE
              WRITE(IOUT,*)'REJECTED VECTOR NO. ',K
              CALL FLUSH(IOUT)
              DO J=1,ISIZE
                V(J,KV)=0
                WRK(J,KV)=0
              END DO
            END IF
          END IF
          DO J=1,ISIZE
            V(J,KV)=V(J,KV)-WRK(J,K)*S
          END DO
        END DO
      END DO
C**THE ABOVE IS SCHMIDT OF W.M.
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE MXBP(X,V,WK,WKR,N,NVEC,NVAL,ILAN,LANCNT,LANK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 WK(N)
      REAL*4 WKR(N)
      DIMENSION LANCNT(1),LANK(1)
      DIMENSION X(N,NVAL),V(N,NVAL),NVEC(NVAL)
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LAN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TOLCUT,EVLCUT,TOLMAT,LANCYC
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/FILASS/IOUT
C**LOOP ROUND COLUMNS OF 'M'
      LAN20=IABS(LLAN20)
      IF(N.GT.LAN20)THEN
        REWIND 46
        REWIND 47
        REWIND 48
        REWIND 49
        REWIND 50
        INP20=46
      ELSE
        REWIND 20
        INP20=20
      END IF
      DO K=1,NVAL
        DO I=1,N
          V(I,K)=0
        END DO
      END DO
      DO I=1,N
        MCOUNT=LANCNT(I)
C**READ Ith COLUMN FROM DISK INTO WK
        READ(INP20)(LANK(JJ),WK(JJ),JJ=1,MCOUNT)
        IF(N.GT.LAN20)THEN
          INP20=INP20+1
          IF(INP20.GT.50)INP20=46
        END IF
        Z=0
C**DIAGONAL ELEMENT
        IF(LANK(MCOUNT).EQ.I)Z=WK(MCOUNT)
        DO K=1,NVAL
          IF(NVEC(K).GT.ILAN)THEN
C**FORM M.X
            DO JJ=1,MCOUNT
              J=LANK(JJ)
              V(J,K)=V(J,K)+WK(JJ)*X(I,K)
              V(I,K)=V(I,K)+WK(JJ)*X(J,K)
            END DO
            V(I,K)=V(I,K)-Z*X(I,K)
          END IF
        END DO
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE MXBM(X,V,WK,WKR,N,NVEC,NVAL,ILAN,LANCNT,LANK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 WK(N)
      REAL*4 WKR(N)
      DIMENSION LANCNT(1),LANK(1)
      DIMENSION X(N,NVAL),V(N,NVAL),NVEC(NVAL)
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LAN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TOLCUT,EVLCUT,TOLMAT,LANCYC
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/FILASS/IOUT
C**LOOP ROUND COLUMNS OF 'M'
      LAN20=IABS(LLAN20)
      IF(N.GT.LAN20)THEN
        REWIND 46
        REWIND 47
        REWIND 48
        REWIND 49
        REWIND 50
        INP20=46
      ELSE
        REWIND 20
        INP20=20
      END IF
      DO K=1,NVAL
        DO I=1,N
          V(I,K)=0
        END DO
      END DO
      DO I=1,N
        MCOUNT=LANCNT(I)
C**READ Ith COLUMN FROM DISK INTO WK
        READ(INP20)(LANK(JJ),WKR(JJ),JJ=1,MCOUNT)
        IF(N.GT.LAN20)THEN
          INP20=INP20+1
          IF(INP20.GT.50)INP20=46
        END IF
        Z=0
        IF(LANK(MCOUNT).EQ.I)Z=WKR(MCOUNT)
        DO K=1,NVAL
          IF(NVEC(K).GT.ILAN)THEN
C**FORM M.X
            DO JJ=1,MCOUNT
              J=LANK(JJ)
              V(J,K)=V(J,K)+WKR(JJ)*X(I,K)
              V(I,K)=V(I,K)+WKR(JJ)*X(J,K)
            END DO
            V(I,K)=V(I,K)-Z*X(I,K)
          END IF
        END DO
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DOT(X,V,W,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),V(N)
      W=0
      DO J=1,N
        W=W+X(J)*V(J)
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE LANOUT(XK,WK,WKR,Y,ZK,ZKR,N,ISTART,IEND,LANCNT,LANK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 WK(N),ZK(N)
      REAL*4 WKR(N),ZKR(N)
      DIMENSION XK(1),Y(N),LANCNT(1),LANK(1)
      LOGICAL LANCZ,LANZA,LANZB
      COMMON/LANCZO/LANCZ,LANZA,LANZB
C**TEMPORARY (DIMENSIONS)
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LAN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TOLCUT,EVLCUT,TOLMAT,LANCYC
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT

      LANCUT=IABS(LLNCUT)
      LAN20=IABS(LLAN20)
      IF(ISTART.EQ.1)THEN
        REWIND 51
        IF(N.GT.LAN20)THEN
          REWIND 46
          REWIND 47
          REWIND 48
          REWIND 49
          REWIND 50
          INP20=46
        ELSE
          REWIND 20
          INP20=20
        END IF
      END IF
      J0=0
      DO I=ISTART,IEND
        LANCNT(I)=0
        Y(I)=XK(J0+I)
        IF(LANCUT.NE.0.AND.LANZB)THEN
          IF(JCOUPL.GT.0)THEN
            K=0
            DO J=1,I
              K=K+1
              WK(K)=XK(J0+J)
            END DO
            WRITE(INP20)(WK(J),J=1,K)
          ELSE
            K=0
            DO J=1,I
              K=K+1
              WKR(K)=XK(J0+J)
            END DO
            WRITE(INP20)(WKR(J),J=1,K)
          END IF
        ELSE
          IF(JCOUPL.GE.0)THEN
            K=0
            DO J=1,I
              K=K+1
              WK(K)=XK(J0+J)
              IF(DABS(WK(K)).GE.TOLMAT)THEN
                LANCNT(I)=LANCNT(I)+1
                LANK(LANCNT(I))=K
              END IF
            END DO
            WRITE(INP20)(LANK(J),WK(LANK(J)),J=1,LANCNT(I))
          ELSE
            K=0
            DO J=1,I
              K=K+1
              WKR(K)=XK(J0+J)
              IF(ABS(WKR(K)).GE.TOLMAT)THEN
                LANCNT(I)=LANCNT(I)+1
                LANK(LANCNT(I))=K
              END IF
            END DO
            WRITE(INP20)(LANK(J),WKR(LANK(J)),J=1,LANCNT(I))
          END IF
        END IF
        IF(N.GT.LAN20)THEN
          INP20=INP20+1
          IF(INP20.GT.50)INP20=46
        END IF
        J0=J0+I
      END DO

      IF(IEND.EQ.N)THEN
        WRITE(51)Y
        CALL FLUSH(51)
        IF(N.GT.LAN20)THEN
          CALL FLUSH(46)
          CALL FLUSH(47)
          CALL FLUSH(48)
          CALL FLUSH(49)
          CALL FLUSH(50)
        ELSE
          CALL FLUSH(20)
        END IF
        IF(LLNCUT.LT.0)STOP 'END LANCZOS SET-UP'
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE QLCOPY(XACOPY,KSIZE,XK,NXK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XACOPY(1),XK(NXK,1)
      L=0
      DO I=1,KSIZE
        DO J=1,I
          L=L+1
          XK(J,I)=XACOPY(L)
          XK(I,J)=XK(J,I)
        END DO
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE LIN(WRK,ISIZE,INP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WRK(ISIZE)
      READ(INP)WRK
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE LOUT(WRK,ISIZE,INP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WRK(ISIZE)
      WRITE(INP)WRK
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE LANPT(XK,EVAL,LAN12D,KFLAG,XL,XD,ISIZE,Z,ZL,ICOUNT,
     1ICOUNTX,JP,JSIZE,IP,ISIZMX,NS,NMODE,LANCUT,ISTART,IEND,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION KFLAG(ISIZE),XD(ISIZE),
     1XL(LAN12D,1),ZL(LAN12D,ICOUNT-ICOUNTX)
      DIMENSION XK(LAN12D,LAN12D),EVAL(LAN12D),Z(LANCUT)
      DIMENSION IP(ISIZMX,NMODE),JP(JSIZE,NMODE)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LLN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TLLCUT,EVLCUT,TOLMAT,LANCYC,LANINC
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/FILASS/IOUT

      WRITE(IOUT,*)'Calculating LANPT'
      CALL FLUSH(IOUT)
      CALL TIMIT(1)

      TOLCUT=DABS(TLLCUT)
      LANCUT=IABS(LLNCUT)

      IF(IND.EQ.1)GO TO 200
      IF(IND.EQ.2)GO TO 300

C**LOOP ROUND REMAINING BASIS
      DO JCOL=ISTART,IEND
        J=JCOL-ISTART+1
C**ZEROISE PERTURBATION TERMS BETWEEN CONTRACTED BASIS
        DO L=1,LANCUT
C**CONVERT BAND ORIGIN IN CM-1 TO ENERGY IN HARTREE
          E=(EVAL(L)+EVLCUT)/WAVENM
          Z(L)=0
          DO I=1,LAN12D
C**LOOP ROUND CONTRACTED BASIS - GET MATRIX ELEMENTS FOR CONTRACTED FNS.
            Z(L)=Z(L)+XK(I,L)*XL(I,J)
          END DO
C**GET PERTURBATION TERM FOR CURRENT COLUMN JCOL
          X=Z(L)*Z(L)/DABS(XD(JCOL)-E)
          IF(X*WAVENM.GT.TOLCUT)THEN
            ICOUNT=ICOUNT+1
            KFLAG(ICOUNT)=JCOL
            GO TO 100
          END IF
        END DO  
100     CONTINUE
      END DO
C**ICOUNT IS TOTAL SIZE FINAL 3-DIM ETC. MATRIX
  
      CALL TIMIT(3)

      RETURN
200   CONTINUE
C**REDUCE SIZE OF XL
      DO I=ICOUNTX+1,ICOUNT
        K=KFLAG(I)-ISTART+1
        DO J=1,LAN12D
          ZL(J,I-ICOUNTX)=XL(J,K)
        END DO
      END DO
      IF(LANINC.GT.1)THEN
        DO I=1,ICOUNT-ICOUNTX
          WRITE(52)(ZL(J,I),J=1,LAN12D)
        END DO
      END IF

      CALL TIMIT(3)

      RETURN
300   CONTINUE
      IF(LANINC.GT.1)THEN
        REWIND 52
        DO I=1,ICOUNT
          READ(52)(ZL(J,I),J=1,LAN12D)
        END DO
        REWIND 52
      END IF
C**SKIP UNWANTED BASES
      JTOT=0
      DO K=1,NS-1
        JTOT=JTOT+NTOT(K)
      END DO
      DO I=1,ICOUNT
        K=KFLAG(I)
        DO M=1,NMODE
          IP(I+LAN12D,M)=JP(JTOT+K,M)
        END DO
      END DO

      CALL TIMIT(3)

      RETURN
      END
