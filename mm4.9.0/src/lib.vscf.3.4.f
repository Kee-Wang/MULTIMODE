C****************************************************************
C****************************************************************
C**DELETED NAG ROUTINES
C****************************************************************
C****************************************************************
      SUBROUTINE E04FBF(M, N, X, F, FMIN, ACC, DSTEP, DMAX, W, IW,
     * FUNCT, MONIT, IPRINT, MAXFUN, IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(M),X(N),W(IW)
      COMMON/XACC/XACC
      COMMON/FILASS/IOUT
      IF (M.GE.1 .AND. N.GE.1 .AND. MAXFUN.GE.1 .AND. ACC.GE.XACC
     * .AND. DSTEP.GE.XACC .AND. DMAX.GE.XACC .AND. (IFAIL.EQ.1
     * .OR. IFAIL.EQ.0) .AND. IW.GE.2*N*(N+M)+2*M+5*N) GO TO 20
      WRITE(IOUT,*)'MFIT,NFIT,MAXCAL,FTOL,XTOL,STEP,IW,XACC',
     1M,N,MAXFUN,ACC,DSTEP,DMAX,IW,XACC
      IFAIL=1
      RETURN
   20 IV = IFAIL
      IFAIL = 0
C     SET VARIOUS PARAMETERS
      MAXC = 0
C     ;MAXC; COUNTS THE NUMBER OF CALLS OF CALFUN
      MPN = M + N
      NT = N + 2
      NTEST = 0
C     ;NT; AND ;NTEST; CAUSE AN ERROR RETURN IF F(X) DOES NOT
C     DECREASE
      DTEST = FLOAT(N+N) - 0.5
C     ;DTEST; IS USED IN A TEST TO MAINTAIN LINEAR INDEPENDENCE
C     PARTITION THE WORKING SPACE ARRAY W
C     THE FIRST PARTITION HOLDS THE JACOBIAN APPROXIMATION
      NWI = M*N
C     THE NEXT PARTITION HOLDS THE GENERALIZED INVERSE
      NWX = NWI + MPN*N
C     THE NEXT PARTITION HOLDS THE BEST VECTOR X
      NWF = NWX + N
C     THE NEXT PARTITION HOLDS THE BEST VECTOR F
      NWC = NWF + M
C     THE NEXT PARTITION HOLDS THE COUNTS OF THE INDEPENDENT
C     DIRECTIONS
      NWD = NWC + N
C     THE NEXT PARTITION HOLDS THE INDEPENDENT DIRECTIONS
      NWW = NWD + N*N
C     THE REMAINDER OF W IS USED FOR SCRATCH VECTORS
      NWV = NWW + N
      NWT = NWV + M
      NWU = NWT + N
      FMIN = -1.
C     USUALLY ;FMIN; IS THE LEAST CALCULATED VALUE OF F(X)
C     JUST AFTER LABEL 132 OF VA05A,IT IS CONVENIENT TO MAKE A
C     DUMMY
C     ASSIGNMENT OF X(1) TO DW, ALTHOUGH DW MAY NOT HAVE BEEN
C     ASSIGNED
C     A VALUE. THIS DOES NOT AFFECT THE RESULTS, BUT MAY CAUSE A
C     FAILURE WITH SOME COMPILERS. SO WE ARBIRTRARILY SET DW TO THE
C     INITIAL VALUE OF 0.0
      DW = 0.0
      DD = 0.
C     USUALLY ;DD; IS THE SQUARE OF THE CURRENT STEP LENGTH
      DSS = DSTEP*DSTEP
      DM = DMAX*DMAX
      PARM = DSQRT(ACC)/DMAX
C     ;PARM; IS THE LEAST VALUE OF THE MARQUARDT PARAMETER
      DPAR = 10.*DM
C     ;DPAR; AND ;NTPAR; ARE USED TO REGULATE THE MARQUARDT
C     PARAMETER
      IS = 4
C     ;IS; CONTROLS A GO TO STATEMENT FOLLOWING A CALL OF CALFUN
      IC = 0
      TINC = 1.
C     ;TINC; IS USED IN THE CRITERION TO INCREASE THE STEP LENGTH
      IPC = 0
      GO TO 80
C     TEST WHETHER THERE HAVE BEEN MAXFUN CALLS OF CALFUN
   40 IF (MAXFUN-MAXC) 60, 60, 80
   60 IFAIL = 2
      GO TO 200
C     CALL THE SUBROUTINE CALFUN
   80 MAXC = MAXC + 1
      CALL FUNCT(M, N, X, F)
C     CALCULATE THE SUM OF SQUARES
      FSQ = 0.
      DO 100 I=1,M
         FSQ = FSQ + F(I)*F(I)
  100 CONTINUE
C     TEST FOR ERROR RETURN BECAUSE F(X) DOES NOT DECREASE
      GO TO (120, 280, 120, 280), IS
  120 IF (FSQ-FMIN) 260, 140, 140
  140 IF (DD-DSS) 160, 160, 280
  160 NTEST = NTEST - 1
      IF (NTEST) 180, 180, 280
  180 IFAIL = 3
  200 DO 220 I=1,N
         NLH = NWX + I
         X(I) = W(NLH)
  220 CONTINUE
      DO 240 I=1,M
         NLH = NWF + I
         F(I) = W(NLH)
  240 CONTINUE
      RETURN
  260 NTEST = NT
  280 IF (IPRINT) 380, 380, 300
  300 IF (IPRINT-1) 360, 360, 320
  320 IPC = IPC - 1
      IF (IPC) 340, 340, 380
  340 IPC = IPRINT
  360 CALL MONIT(M, N, X, F, FSQ, MAXC)
  380 GO TO (1860, 1940, 1940, 400), IS
C     STORE THE INITIAL VECTORS X AND F
  400 IF (IC) 420, 420, 460
  420 DO 440 I=1,N
         NLH = NWX + I
         W(NLH) = X(I)
  440 CONTINUE
      GO TO 540
C     CALCULATE THE INITIAL JACOBIAN APPROXIMATION
  460 K = IC
      DO 480 I=1,M
         NLH = NWF + I
         W(K) = (F(I)-W(NLH))/DSTEP
         K = K + N
  480 CONTINUE
C     TEST WHETHER THE MOST RECENT X IS BEST
      IF (FMIN-FSQ) 500, 500, 520
  500 NLH = NWX + IC
      X(IC) = W(NLH)
      GO TO 580
  520 NLH = NWX + IC
      W(NLH) = X(IC)
  540 DO 560 I=1,M
         NLH = NWF + I
         W(NLH) = F(I)
  560 CONTINUE
      FMIN = FSQ
C     SET X FOR THE NEXT CALL OF CALFUN
  580 IC = IC + 1
      IF (IC-N) 600, 600, 620
  600 NLH = NWX + IC
      X(IC) = W(NLH) + DSTEP
      GO TO 80
C     SET THE DIRECTION MATRIX
  620 K = NWD
      DO 660 I=1,N
         DO 640 J=1,N
            K = K + 1
            W(K) = 0.
  640    CONTINUE
         NLH = K + I - N
         W(NLH) = 1.
         NLH = NWC + I
         MLH = N - I
         W(NLH) = 1. + FLOAT(MLH)
  660 CONTINUE
C     SET THE MARQUARDT PARAMETER TO ITS LEAST VALUE
  680 PAR = PARM
C     COPY THE JACOBIAN AND APPEND THE MARQUARDT MATRIX
  700 PPAR = PAR*PAR
      NTPAR = 0
  720 KK = 0
      K = NWI + NWI
      DO 780 I=1,N
         DO 740 J=1,M
            KK = KK + 1
            NLH = KK + NWI
            W(NLH) = W(KK)
  740    CONTINUE
         DO 760 J=1,N
            K = K + 1
            W(K) = 0.
  760    CONTINUE
         NLH = K + I - N
         W(NLH) = PAR
  780 CONTINUE
C     CALCULATE THE GENERALIZED INVERSE OF J
      IIW = 3*N + M
      CALL E04FBZ(N, MPN, W(NWI+1), N, W(NWW+1), IIW)
C     NOTE THAT THE THIRD AND FIFTH ENTRIES OF THIS ARGUMENT LIST
C     STAND FOR ONE-DIMENSIONAL ARRAYS.
C     START THE ITERATION BY TESTING FMIN
  800 IF (FMIN-ACC) 200, 200, 820
C     NEXT PREDICT THE DESCENT AND MARQUARDT MINIMA
  820 DS = 0.
      DN = 0.
      SP = 0.
      DO 860 I=1,N
         X(I) = 0.
         F(I) = 0.
         K = I
         DO 840 J=1,M
            NLH = NWF + J
            X(I) = X(I) - W(K)*W(NLH)
            MLH = NWI + K
            F(I) = F(I) - W(MLH)*W(NLH)
            K = K + N
  840    CONTINUE
         DS = DS + X(I)*X(I)
         DN = DN + F(I)*F(I)
         SP = SP + X(I)*F(I)
  860 CONTINUE
C     PREDICT THE REDUCTION IN F(X) DUE TO THE MARQUARDT STEP
C     AND ALSO PREDICT THE LENGTH OF THE STEEPEST DESCENT STEP
      PRED = SP + SP
      DMULT = 0.
      K = 0
      DO 900 I=1,M
         AP = 0.
         AD = 0.
         DO 880 J=1,N
            K = K + 1
            AP = AP + W(K)*F(J)
            AD = AD + W(K)*X(J)
  880    CONTINUE
         PRED = PRED - AP*AP
         DMULT = DMULT + AD*AD
  900 CONTINUE
C     TEST FOR CONVERGENCE
      IF (DN-DM) 920, 920, 940
  920 AP = DSQRT(DN)
      IF (PRED+2.*PPAR*AP*(DMAX-AP)-ACC) 200, 200, 960
  940 IF (PRED+PPAR*(DM-DN)-ACC) 200, 200, 960
C     TEST WHETHER TO APPLY THE FULL MARQUARDT CORRECTION
  960 DMULT = DS/DMULT
      DS = DS*DMULT*DMULT
  980 IS = 2
      IF (DN-DD) 1000, 1000, 1060
C     TEST THAT THE MARQUARDT PARAMETER HAS ITS LEAST VALUE
 1000 IF (PAR-PARM) 1020, 1020, 680
1020  DD=DMAX1(DN,DSS)
      DS = 0.25*DN
      TINC = 1.
      IF (DN-DSS) 1040, 1480, 1480
 1040 IS = 3
      GO TO 1820
C     TEST WHETHER TO INCREASE THE MARQUARDT PARAMETER
 1060 IF (DN-DPAR) 1080, 1080, 1100
 1080 NTPAR = 0
      GO TO 1220
 1100 IF (NTPAR) 1120, 1120, 1140
 1120 NTPAR = 1
      PTM = DN
      GO TO 1220
 1140 NTPAR = NTPAR + 1
      PTM = DMIN1(PTM,DN)
      IF (NTPAR-NT) 1220, 1160, 1160
C     SET THE LARGER VALUE OF THE MARQUARDT PARAMETER
 1160 PAR = PAR*(PTM/DM)**0.25
      IF (6.*DD-DM) 1180, 700, 700
 1180 AP = DSQRT(PRED/DN)
      IF (AP-PAR) 700, 700, 1200
 1200 PAR = DMIN1(AP,PAR*(DM/(6.*DD))**0.25)
      GO TO 700
C     TEST WHETHER TO USE THE STEEPEST DESCENT DIRECTION
 1220 IF (DS-DD) 1320, 1240, 1240
C     TEST WHETHER THE INITIAL VALUE OF DD HAS BEEN SET
 1240 IF (DD) 1260, 1260, 1300
 1260 DD = DMIN1(DM,DS)
      IF (DD-DSS) 1280, 1300, 1300
 1280 DD = DSS
      GO TO 980
C     SET THE MULTIPLIER OF THE STEEPEST DESCENT DIRECTION
 1300 ANMULT = 0.
      DMULT = DMULT*DSQRT(DD/DS)
      GO TO 1340
C     INTERPOLATE BETWEEN THE STEEPEST DESCENT AND MARQUARDT
C     DIRECTIONS
 1320 SP = SP*DMULT
      ANMULT = (DD-DS)/((SP-DS)+DSQRT((SP-DD)**2+(DN-DD)*(DD-DS)))
      DMULT = DMULT*(1.-ANMULT)
C     CALCULATE THE CORRECTION TO X, AND ITS ANGLE WITH THE FIRST
C     DIRECTION
 1340 DN = 0.
      SP = 0.
      DO 1360 I=1,N
         F(I) = DMULT*X(I) + ANMULT*F(I)
         DN = DN + F(I)*F(I)
         NLH = NWD + I
         SP = SP + F(I)*W(NLH)
 1360 CONTINUE
      DS = 0.25*DN
C     TEST WHETHER AN EXTRA STEP IS NEEDED FOR INDEPENDENCE
      IF (W(NWC+1)-DTEST) 1480, 1480, 1380
 1380 IF (SP*SP-DS) 1400, 1480, 1480
C     TAKE THE EXTRA STEP AND UPDATE THE DIRECTION MATRIX
 1400 DO 1420 I=1,N
         NLH = NWX + I
         MLH = NWD + I
         X(I) = W(NLH) + DSTEP*W(MLH)
         NLH = NWC + I
         W(NLH) = W(NLH+1) + 1.
 1420 CONTINUE
      W(NWD) = 1.
      DO 1460 I=1,N
         K = NWD + I
         SP = W(K)
         DO 1440 J=2,N
            NLH = K + N
            W(K) = W(NLH)
            K = K + N
 1440    CONTINUE
         W(K) = SP
 1460 CONTINUE
      GO TO 40
C     EXPRESS THE NEW DIRECTION IN TERMS OF THOSE OF THE DIRECTION
C     MATRIX, AND UPDATE THE COUNTS IN W(NWC+1) ETC.
 1480 SP = 0.
      K = NWD
      DO 1600 I=1,N
         X(I) = DW
         DW = 0.
         DO 1500 J=1,N
            K = K + 1
            DW = DW + F(J)*W(K)
 1500    CONTINUE
         GO TO (1560, 1520), IS
 1520    NLH = NWC + I
         W(NLH) = W(NLH) + 1.
         SP = SP + DW*DW
         IF (SP-DS) 1600, 1600, 1540
 1540    IS = 1
         KK = I
         X(1) = DW
         GO TO 1580
 1560    X(I) = DW
 1580    NLH = NWC + I
         W(NLH) = W(NLH+1) + 1.
 1600 CONTINUE
      W(NWD) = 1.
C     REORDER THE DIRECTIONS SO THAT KK IS FIRST
      IF (KK-1) 1680, 1680, 1620
 1620 KS = NWC + KK*N
      DO 1660 I=1,N
         K = KS + I
         SP = W(K)
         DO 1640 J=2,KK
            NLH = K - N
            W(K) = W(NLH)
            K = K - N
 1640    CONTINUE
         W(K) = SP
 1660 CONTINUE
C     GENERATE THE NEW ORTHOGONAL DIRECTION MATRIX
 1680 DO 1700 I=1,N
         NLH = NWW + I
         W(NLH) = 0.
 1700 CONTINUE
      SP = X(1)*X(1)
      K = NWD
      DO 1740 I=2,N
         DS = DSQRT(SP*(SP+X(I)*X(I)))
         DW = SP/DS
         DS = X(I)/DS
         SP = SP + X(I)*X(I)
         DO 1720 J=1,N
            K = K + 1
            NLH = NWW + J
            W(NLH) = W(NLH) + X(I-1)*W(K)
            MLH = K + N
            W(K) = DW*W(MLH) - DS*W(NLH)
 1720    CONTINUE
 1740 CONTINUE
      SP = 1./DSQRT(DN)
      DO 1760 I=1,N
         K = K + 1
         W(K) = SP*F(I)
 1760 CONTINUE
C     PREDICT THE NEW RIGHT HAND SIDES
      FNP = 0.
      K = 0
      DO 1800 I=1,M
         NLH = NWW + I
         MLH = NWF + I
         W(NLH) = W(MLH)
         DO 1780 J=1,N
            K = K + 1
            W(NLH) = W(NLH) + W(K)*F(J)
 1780    CONTINUE
         FNP = FNP + W(NLH)**2
 1800 CONTINUE
C     CALCULATE THE NEXT VECTOR X, AND THEN CALL CALFUN
 1820 DO 1840 I=1,N
         NLH = NWX + I
         X(I) = W(NLH) + F(I)
 1840 CONTINUE
      GO TO 40
C     UPDATE THE STEP SIZE
 1860 DMULT = 0.9*FMIN + 0.1*FNP - FSQ
      IF (DMULT) 1880, 1900, 1900
 1880 DD = DMAX1(DSS,0.25*DD)
      TINC = 1.
      IF (FSQ-FMIN) 1960, 2060, 2060
C     TRY THE TEST TO DECIDE WHETHER TO INCREASE THE STEP LENGTH
 1900 SP = 0.
      SS = 0.
      DO 1920 I=1,M
         NLH = NWW + I
         SP = SP + DABS(F(I)*(F(I)-W(NLH)))
         SS = SS + (F(I)-W(NLH))**2
 1920 CONTINUE
      PJ = 1. + DMULT/(SP+DSQRT(SP*SP+DMULT*SS))
      SP = DMIN1(4.D0,TINC,PJ)
      TINC = PJ/SP
      DD = DMIN1(DM,SP*DD)
      GO TO 1960
C     IF F(X) IMPROVES STORE THE NEW VALUE OF X
 1940 IF (FSQ-FMIN) 1960, 2020, 2020
 1960 FMIN = FSQ
      DO 1980 I=1,N
         SP = X(I)
         NLH = NWX + I
         X(I) = W(NLH)
         W(NLH) = SP
 1980 CONTINUE
      DO 2000 I=1,M
         SP = F(I)
         NLH = NWF + I
         F(I) = W(NLH)
         W(NLH) = SP
 2000 CONTINUE
 2020 GO TO (2060, 2060, 2040), IS
 2040 IS = 2
      IF (FMIN-ACC) 200, 200, 1400
C     CALCULATE THE CHANGES IN X AND IN F
 2060 DS = 0.
      DO 2080 I=1,N
         NLH = NWX + I
         X(I) = X(I) - W(NLH)
         DS = DS + X(I)*X(I)
 2080 CONTINUE
      DO 2100 I=1,M
         NLH = NWF + I
         F(I) = F(I) - W(NLH)
 2100 CONTINUE
C     CALCULATE THE GENERALIZED INVERSE TIMES THE CHANGE IN X
      K = NWI
      SS = 0.
      DO 2140 I=1,MPN
         SP = 0.
         DO 2120 J=1,N
            K = K + 1
            SP = SP + W(K)*X(J)
 2120    CONTINUE
         NLH = NWV + I
         W(NLH) = SP
         SS = SS + SP*SP
 2140 CONTINUE
C     CALCULATE J TIMES THE CHANGE IN F
C     ALSO APPLY PROJECTION TO THE GENERALIZED INVERSE
      DO 2220 I=1,N
         ST = 0.
         K = NWI + I
         DO 2160 J=1,MPN
            NLH = J + NWV
            ST = ST + W(K)*W(NLH)
            K = K + N
 2160    CONTINUE
         ST = ST/SS
         K = NWI + I
         DO 2180 J=1,MPN
            NLH = J + NWV
            W(K) = W(K) - ST*W(NLH)
            K = K + N
 2180    CONTINUE
         ST = PPAR*X(I)
         K = I
         DO 2200 J=1,M
            ST = ST + W(K)*F(J)
            K = K + N
 2200    CONTINUE
         NLH = NWW + I
         W(NLH) = ST
 2220 CONTINUE
C     REVISE J AND CALCULATE ROW VECTOR FOR CORRECTION TO INVERSE
      IC = 0
      K = 0
      KK = NWI
      SP = 0.
      SPP = 0.
      DO 2280 I=1,M
         SS = F(I)
         ST = F(I)
         DO 2240 J=1,N
            IC = IC + 1
            KK = KK + 1
            SS = SS - W(IC)*X(J)
            NLH = NWW + J
            ST = ST - W(KK)*W(NLH)
 2240    CONTINUE
         SS = SS/DS
         NLH = NWV + I
         W(NLH) = ST
         SP = SP + F(I)*ST
         SPP = SPP + ST*ST
         DO 2260 J=1,N
            K = K + 1
            W(K) = W(K) + SS*X(J)
 2260    CONTINUE
 2280 CONTINUE
      DO 2320 I=1,N
         ST = PAR*X(I)
         DO 2300 J=1,N
            KK = KK + 1
            NLH = NWW + J
            ST = ST - W(KK)*W(NLH)
 2300    CONTINUE
         NLH = NWT + I
         W(NLH) = ST
         SP = SP + PAR*X(I)*ST
         SPP = SPP + ST*ST
 2320 CONTINUE
C     TEST THAT THE SCALAR PRODUCT IS SUFFICIENTLY ACCURATE
      IF (0.01*SPP-DABS(SP-SPP)) 720, 720, 2340
C     CALCULATE THE NEW GENERALIZED INVERSE
 2340 DO 2420 I=1,N
         K = NWI + I
         ST = X(I)
         DO 2360 J=1,M
            ST = ST - W(K)*F(J)
            K = K + N
 2360    CONTINUE
         SS = 0.
         DO 2380 J=1,N
            SS = SS + W(K)*X(J)
            K = K + N
 2380    CONTINUE
         ST = (ST-PAR*SS)/SP
         K = NWI + I
         DO 2400 J=1,MPN
            NLH = NWV + J
            W(K) = W(K) + ST*W(NLH)
            K = K + N
 2400    CONTINUE
 2420 CONTINUE
      GO TO 800
      END
      SUBROUTINE E04FBZ(M, N, A, IA, W, IW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,N),W(IW)
      NRW = M
C     THE SECOND PARTITION RECORDS ROW INTERCHANGES
      NCW = M + M
C     THE THIRD PARTITION RECORDS COLUMN INTERCHANGES
C     SET THE INITIAL RECORDS OF ROW AND COLUMN INTERCHANGES
      DO 20 I=1,M
         NLH = NRW + I
         W(NLH) = 0.5 + FLOAT(I)
   20 CONTINUE
      DO 40 I=1,N
         NLH = NCW + I
         W(NLH) = 0.5 + FLOAT(I)
   40 CONTINUE
C     ;KK; COUNTS THE SEPARATE ELEMENTARY TRANSFORMATIONS
      KK = 1
C     FIND LARGEST ROW AND MAKE ROW INTERCHANGES
   60 RMAX = 0.
      DO 120 I=KK,M
         SUM = 0.
         DO 80 J=KK,N
            SUM = SUM + A(I,J)**2
   80    CONTINUE
         IF (RMAX-SUM) 100, 120, 120
  100    RMAX = SUM
         IR = I
  120 CONTINUE
      IF (IR-KK) 180, 180, 140
  140 NLH = NRW + KK
      SUM = W(NLH)
      MLH = NRW + IR
      W(NLH) = W(MLH)
      W(MLH) = SUM
      DO 160 J=1,N
         SUM = A(KK,J)
         A(KK,J) = A(IR,J)
         A(IR,J) = SUM
  160 CONTINUE
C     FIND LARGEST ELEMENT OF PIVOTAL ROW, AND MAKE COLUMN
C     INTERCHANGES
  180 RMAX = 0.
      SUM = 0.
      DO 220 J=KK,N
         SUM = SUM + A(KK,J)**2
         IF (RMAX-DABS(A(KK,J))) 200, 220, 220
  200    RMAX = DABS(A(KK,J))
         IR = J
  220 CONTINUE
      IF (IR-KK) 280, 280, 240
  240 NLH = NCW + KK
      RMAX = W(NLH)
      MLH = NCW + IR
      W(NLH) = W(MLH)
      W(MLH) = RMAX
      DO 260 I=1,M
         RMAX = A(I,KK)
         A(I,KK) = A(I,IR)
         A(I,IR) = RMAX
  260 CONTINUE
C     REPLACE THE PIVOTAL ROW BY THE VECTOR OF THE TRANSFORMATION
  280 SIGMA = DSQRT(SUM)
      BSQ = DSQRT(SUM+SIGMA*DABS(A(KK,KK)))
      W(KK) = DSIGN(SIGMA+DABS(A(KK,KK)),A(KK,KK))/BSQ
      A(KK,KK) = -DSIGN(SIGMA,A(KK,KK))
      KP = KK + 1
      IF (KP-N) 300, 300, 420
  300 DO 320 J=KP,N
         A(KK,J) = A(KK,J)/BSQ
  320 CONTINUE
C     APPLY THE TRANSFORMATION TO THE REMAINING ROWS OF A
      IF (KP-M) 340, 340, 420
  340 DO 400 I=KP,M
         SUM = W(KK)*A(I,KK)
         DO 360 J=KP,N
            SUM = SUM + A(KK,J)*A(I,J)
  360    CONTINUE
         A(I,KK) = A(I,KK) - SUM*W(KK)
         DO 380 J=KP,N
            A(I,J) = A(I,J) - SUM*A(KK,J)
  380    CONTINUE
  400 CONTINUE
      KK = KP
      GO TO 60
C     AT THIS STAGE THE REDUCTION OF A IS COMPLETE
C     NOW WE BUILD UP THE GENERALIZED INVERSE
C     APPLY THE FIRST ELEMENTARY TRANSFORMATION
  420 KK = M
      KP = M + 1
      SUM = W(M)/A(M,M)
      IF (N-M) 480, 480, 440
  440 DO 460 J=KP,N
         A(M,J) = -SUM*A(M,J)
  460 CONTINUE
  480 A(M,M) = 1./A(M,M) - SUM*W(M)
C     NOW APPLY THE OTHER (M-1) TRANSFORMATIONS
  500 KP = KK
      KK = KP - 1
      IF (KK) 660, 660, 520
C     FIRST TRANSFORM THE LAST (M-KK) ROWS
  520 DO 580 I=KP,M
         SUM = 0.
         DO 540 J=KP,N
            SUM = SUM + A(KK,J)*A(I,J)
  540    CONTINUE
         DO 560 J=KP,N
            A(I,J) = A(I,J) - SUM*A(KK,J)
  560    CONTINUE
         W(I) = -SUM*W(KK)
  580 CONTINUE
C     THEN CALCULATE THE NEW ROW IN POSITION KK
      DO 620 J=KP,N
         SUM = -W(KK)*A(KK,J)
         DO 600 I=KP,M
            SUM = SUM - A(I,KK)*A(I,J)
  600    CONTINUE
         A(KK,J) = SUM/A(KK,KK)
  620 CONTINUE
C     AND REVISE THE COLUMN IN POSITION KK
      SUM = 1. - W(KK)**2
      DO 640 I=KP,M
         SUM = SUM - A(I,KK)*W(I)
         A(I,KK) = W(I)
  640 CONTINUE
      A(KK,KK) = SUM/A(KK,KK)
      GO TO 500
C     RESTORE THE ROW INTERCHANGES
  660 DO 740 I=1,M
  680    NLH = NRW + I
         IR = IFIX(SNGL(W(NLH)))
         IF (I-IR) 700, 740, 740
  700    SUM = W(NLH)
         MLH = NRW + IR
         W(NLH) = W(MLH)
         W(MLH) = SUM
         DO 720 J=1,N
            SUM = A(I,J)
            A(I,J) = A(IR,J)
            A(IR,J) = SUM
  720    CONTINUE
         GO TO 680
  740 CONTINUE
C     RESTORE THE COLUMN INTERCHANGES
      DO 820 J=1,N
  760    NLH = NCW + J
         IR = IFIX(SNGL(W(NLH)))
         IF (J-IR) 780, 820, 820
  780    SUM = W(NLH)
         MLH = NCW + IR
         W(NLH) = W(MLH)
         W(MLH) = SUM
         DO 800 I=1,M
            SUM = A(I,J)
            A(I,J) = A(I,IR)
            A(I,IR) = SUM
  800    CONTINUE
         GO TO 760
  820 CONTINUE
      RETURN
      END
