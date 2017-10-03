C**************************************************************
C**************************************************************
C**CONTRACTION SCHEMES
C**************************************************************
C**************************************************************
      SUBROUTINE GETCP0(JP2,JSIZE2,JCI1,JCI2,ITOT2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JP2(JSIZE2,2)
      DO I1=1,JCI1
        DO I2=1,JCI2
          ITOT2=ITOT2+1
          JP2(ITOT2,1)=I1
          JP2(ITOT2,2)=I2
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETCPS(JP,JSIZE,IP,ISIZE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JP(JSIZE,2),IP(ISIZE,2)
      DO K=1,2
        DO I=1,ISIZE
          IP(I,K)=JP(I,K)
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETCP2(IP2,ISIZE2,NMODE,KMODE,LMODE,MAXBAS,NSIZE2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NMODE,2),IP2(ISIZE2,NMODE)
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/FILASS/IOUT,INP
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      DO I=1,ICONT(1)
        IF(KMODE.EQ.JCONT(1,I))THEN
          IS1=1
          JS1=KMODE
          KS1=MAXBFN(MAXBAS,NNMODE,KMODE,1)
        END IF
        IF(LMODE.EQ.JCONT(1,I))THEN
          IS1=2
          JS1=LMODE
          KS1=MAXBFN(MAXBAS,NNMODE,LMODE,1)
        END IF
      END DO
      DO I=1,ICONT(2)
        IF(KMODE.EQ.JCONT(2,I))THEN
          IS2=1
          JS2=KMODE
          KS2=MAXBFN(MAXBAS,NNMODE,KMODE,1)
        END IF
        IF(LMODE.EQ.JCONT(2,I))THEN
          IS2=2
          JS2=LMODE
          KS2=MAXBFN(MAXBAS,NNMODE,LMODE,1)
        END IF
      END DO
C**NOW PUT IN SCHEME 1 BASIS
      NSZ1=1
      NMAX=MAXSUM(1,1)+1
      CALL GETJC1(IP2,ISIZE2,NMODE,IS1,KS1,NMAX,NSZ1,
     1MAXBAS,2,JS1,2)
C**NOW GET SIZE OF SCHEME 2 BASIS
      NSZ2=1
      NMAX=MAXSUM(1,2)+1
      CALL GETJC1(IP2,ISIZE2,NMODE,IS2,KS2,NMAX,NSZ2,
     1MAXBAS,-2,JS2,2)
C**TOTAL SIZE IS DIRECT PRODUCT
      NSIZE2=NSZ1*NSZ2
C**DUPLICATE SCHEME 1 BASIS FOR ALL SCHEME 2
      K=NSIZE2+1
      DO I=1,NSZ1
        J=NSZ1+1-I
        DO M=1,NSZ2
          K=K-1
          IP2(K,IS1)=IP2(J,IS1)
        END DO
      END DO
C**NOW FILL IN SCHEME 2 BASIS
      NSZ2=0
      DO I=1,NSZ1
        NMAX=MAXSUM(1,2)+1
        CALL GETJC1(IP2,ISIZE2,NMODE,IS2,KS2,NMAX,NSZ2,
     1  MAXBAS,0,JS2,1)
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETCP3(IP3,ISIZE3,NNMODE,KMODE,LMODE,NMODE,MAXBAS,
     1NSIZE3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NNMODE,3),IP3(ISIZE3,NNMODE)
      DIMENSION IS1(2),IS2(2),JS1(2),JS2(2),KS1(2),KS2(2)
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/FILASS/IOUT,INP
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
C**SEPARATE MODES INTO SCHEMES
      NS1=0
      DO I=1,ICONT(1)
        IF(KMODE.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=1
          JS1(NS1)=KMODE
          KS1(NS1)=MAXBFN(MAXBAS,NNMODE,KMODE,1)
        END IF
        IF(LMODE.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=2
          JS1(NS1)=LMODE
          KS1(NS1)=MAXBFN(MAXBAS,NNMODE,LMODE,1)
        END IF
        IF(NMODE.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=3
          JS1(NS1)=NMODE
          KS1(NS1)=MAXBFN(MAXBAS,NNMODE,NMODE,1)
        END IF
      END DO
      NS2=0
      DO I=1,ICONT(2)
        IF(KMODE.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=1
          JS2(NS2)=KMODE
          KS2(NS2)=MAXBFN(MAXBAS,NNMODE,KMODE,1)
        END IF
        IF(LMODE.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=2
          JS2(NS2)=LMODE
          KS2(NS2)=MAXBFN(MAXBAS,NNMODE,LMODE,1)
        END IF
        IF(NMODE.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=3
          JS2(NS2)=NMODE
          KS2(NS2)=MAXBFN(MAXBAS,NNMODE,NMODE,1)
        END IF
      END DO
C**NOW PUT IN SCHEME 1 BASIS
      NSIZE1=1
      NMAX=MAXSUM(1,1)+1
      DO I=1,NS1
        IND=3
        IF(I.GT.1)IND=0
        CALL GETJC1(IP3,ISIZE3,NNMODE,IS1(I),KS1(I),NMAX,NSIZE1,
     1  MAXBAS,IND,JS1(I),2)
      END DO
      IF(NS1.EQ.1)GO TO 7777
      NMAX=MAXSUM(2,1)+2
      DO I=1,NS1
        DO J=1,I-1
          CALL GETJC2(IP3,ISIZE3,NNMODE,IS1(I),IS1(J),KS1(I),KS1(J),
     1    NMAX,NSIZE1,MAXBAS,0,JS1(I),JS1(J))
        END DO
      END DO
7777  CONTINUE
C**NOW GET SIZE OF SCHEME 2 BASIS
      NSIZE2=1
      NMAX=MAXSUM(1,2)+1
      DO I=1,NS2
        CALL GETJC1(IP3,ISIZE3,NNMODE,IS2(I),KS2(I),NMAX,NSIZE2,
     1  MAXBAS,-3,JS2(I),2)
      END DO
      IF(NS2.EQ.1)GO TO 8888
      NMAX=MAXSUM(2,2)+2
      DO I=1,NS2
        DO J=1,I-1
          CALL GETJC2(IP3,ISIZE3,NNMODE,IS2(I),IS2(J),KS2(I),KS2(J),
     1    NMAX,NSIZE2,MAXBAS,-3,JS2(I),JS2(J))
        END DO
      END DO
8888  CONTINUE
C**TOTAL SIZE IS DIRECT PRODUCT
      NSIZE3=NSIZE1*NSIZE2
C**DUPLICATE SCHEME 1 BASIS FOR ALL SCHEME 2
      K=NSIZE3+1
      DO I=1,NSIZE1
        J=NSIZE1+1-I
        DO M=1,NSIZE2
          K=K-1
          DO N=1,NS1
          IP3(K,IS1(N))=IP3(J,IS1(N))
          END DO
        END DO
      END DO
C**NOW FILL IN SCHEME 2 BASIS
      NSIZE2=0
      DO N=1,NSIZE1
        NMAX=MAXSUM(1,2)+1
        DO I=1,NS2
          IND=1
          IF(I.GT.1)IND=2
          CALL GETJC1(IP3,ISIZE3,NNMODE,IS2(I),KS2(I),NMAX,NSIZE2,
     1    MAXBAS,0,JS2(I),IND)
        END DO
        IF(NS2.EQ.1)GO TO 9999
        NMAX=MAXSUM(2,2)+2
        DO I=1,NS2
          DO J=1,I-1
            CALL GETJC2(IP3,ISIZE3,NNMODE,IS2(I),IS2(J),KS2(I),KS2(J),
     1      NMAX,NSIZE2,MAXBAS,0,JS2(I),JS2(J))
          END DO
        END DO
9999    CONTINUE
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETCP4(IP4,ISIZE4,NMODE,K,L,N,M,MAXBAS,NSIZE4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NMODE,4),IP4(ISIZE4,NMODE)
      DIMENSION IS1(3),IS2(3),JS1(3),JS2(3),KS1(3),KS2(3)
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/FILASS/IOUT,INP
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
C**SEPARATE MODES INTO SCHEMES
      NS1=0
      DO I=1,ICONT(1)
        IF(K.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=1
          JS1(NS1)=K
          KS1(NS1)=MAXBFN(MAXBAS,NMODE,K,1)
        END IF
        IF(L.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=2
          JS1(NS1)=L
          KS1(NS1)=MAXBFN(MAXBAS,NMODE,L,1)
        END IF
        IF(N.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=3
          JS1(NS1)=N
          KS1(NS1)=MAXBFN(MAXBAS,NMODE,N,1)
        END IF
        IF(M.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=4
          JS1(NS1)=M
          KS1(NS1)=MAXBFN(MAXBAS,NMODE,M,1)
        END IF
      END DO
      NS2=0
      DO I=1,ICONT(2)
        IF(K.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=1
          JS2(NS2)=K
          KS2(NS2)=MAXBFN(MAXBAS,NMODE,K,1)
        END IF
        IF(L.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=2
          JS2(NS2)=L
          KS2(NS2)=MAXBFN(MAXBAS,NMODE,L,1)
        END IF
        IF(N.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=3
          JS2(NS2)=N
          KS2(NS2)=MAXBFN(MAXBAS,NMODE,N,1)
        END IF
        IF(M.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=4
          JS2(NS2)=M
          KS2(NS2)=MAXBFN(MAXBAS,NMODE,M,1)
        END IF
      END DO
C**NOW PUT IN SCHEME 1 BASIS
      NSIZE1=1
      NMAX=MAXSUM(1,1)+1
      DO I=1,NS1
        IND=4
        IF(I.GT.1)IND=0
        CALL GETJC1(IP4,ISIZE4,NMODE,IS1(I),KS1(I),NMAX,NSIZE1,
     1  MAXBAS,IND,JS1(I),2)
      END DO
      IF(NS1.EQ.1)GO TO 7777
      NMAX=MAXSUM(2,1)+2
      DO I=1,NS1
        DO J=1,I-1
          CALL GETJC2(IP4,ISIZE4,NMODE,IS1(I),IS1(J),KS1(I),KS1(J),
     1    NMAX,NSIZE1,MAXBAS,0,JS1(I),JS1(J))
        END DO
      END DO
      IF(NS1.EQ.2)GO TO 7777
      NMAX=MAXSUM(3,1)+3
      DO I=1,NS1
        DO J=1,I-1
          DO KK=1,J-1
            CALL GETJC3(IP4,ISIZE4,NMODE,IS1(I),IS1(J),IS1(KK),
     1      KS1(I),KS1(J),KS1(KK),NMAX,NSIZE1,MAXBAS,0,JS1(I),JS1(J),
     2      JS1(KK))
          END DO
        END DO
      END DO
7777  CONTINUE
C**NOW GET SIZE OF SCHEME 2 BASIS
      NSIZE2=1
      NMAX=MAXSUM(1,2)+1
      DO I=1,NS2
        CALL GETJC1(IP4,ISIZE4,NMODE,IS2(I),KS2(I),NMAX,NSIZE2,
     1  MAXBAS,-4,JS2(I),2)
      END DO
      IF(NS2.EQ.1)GO TO 8888
      NMAX=MAXSUM(2,2)+2
      DO I=1,NS2
        DO J=1,I-1
          CALL GETJC2(IP4,ISIZE4,NMODE,IS2(I),IS2(J),KS2(I),KS2(J),
     1    NMAX,NSIZE2,MAXBAS,-4,JS2(I),JS2(J))
        END DO
      END DO
      IF(NS2.EQ.2)GO TO 8888
      NMAX=MAXSUM(3,2)+3
      DO I=1,NS2
        DO J=1,I-1
          DO KK=1,J-1
            CALL GETJC3(IP4,ISIZE4,NMODE,IS2(I),IS2(J),IS2(KK),KS2(I),
     1      KS2(J),KS2(KK),NMAX,NSIZE2,MAXBAS,-4,JS2(I),JS2(J),JS2(KK))
          END DO
        END DO
      END DO
8888  CONTINUE
C**TOTAL SIZE IS DIRECT PRODUCT
      NSIZE4=NSIZE1*NSIZE2
C**DUPLICATE SCHEME 1 BASIS FOR ALL SCHEME 2
      KK=NSIZE4+1
      DO I=1,NSIZE1
        J=NSIZE1+1-I
        DO MM=1,NSIZE2
          KK=KK-1
          DO NN=1,NS1
          IP4(KK,IS1(NN))=IP4(J,IS1(NN))
          END DO
        END DO
      END DO
C**NOW FILL IN SCHEME 2 BASIS
      NSIZE2=0
      DO NN=1,NSIZE1
        NMAX=MAXSUM(1,2)+1
        DO I=1,NS2
          IND=1
          IF(I.GT.1)IND=2
          CALL GETJC1(IP4,ISIZE4,NMODE,IS2(I),KS2(I),NMAX,NSIZE2,
     1    MAXBAS,0,JS2(I),IND)
        END DO
        IF(NS2.EQ.1)GO TO 9999
        NMAX=MAXSUM(2,2)+2
        DO I=1,NS2
          DO J=1,I-1
            CALL GETJC2(IP4,ISIZE4,NMODE,IS2(I),IS2(J),KS2(I),KS2(J),
     1      NMAX,NSIZE2,MAXBAS,0,JS2(I),JS2(J))
          END DO
        END DO
        IF(NS2.EQ.2)GO TO 9999
        NMAX=MAXSUM(3,2)+3
        DO I=1,NS2
          DO J=1,I-1
            DO KK=1,J-1
              CALL GETJC3(IP4,ISIZE4,NMODE,IS2(I),IS2(J),IS2(KK),
     1        KS2(I),KS2(J),KS2(KK),NMAX,NSIZE2,MAXBAS,0,JS2(I),JS2(J),
     2        JS2(KK))
            END DO
          END DO
        END DO
9999    CONTINUE
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETCP5(IP5,ISIZE5,NMODE,K,L,N,M,IT,MAXBAS,NSIZE5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NMODE,5),IP5(ISIZE5,NMODE)
      DIMENSION IS1(4),IS2(4),JS1(4),JS2(4),KS1(4),KS2(4)
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/FILASS/IOUT,INP
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
C**SEPARATE MODES INTO SCHEMES
      NS1=0
      DO I=1,ICONT(1)
        IF(K.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=1
          JS1(NS1)=K
          KS1(NS1)=MAXBFN(MAXBAS,NMODE,K,1)
        END IF
        IF(L.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=2
          JS1(NS1)=L
          KS1(NS1)=MAXBFN(MAXBAS,NMODE,L,1)
        END IF
        IF(N.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=3
          JS1(NS1)=N
          KS1(NS1)=MAXBFN(MAXBAS,NMODE,N,1)
        END IF
        IF(M.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=4
          JS1(NS1)=M
          KS1(NS1)=MAXBFN(MAXBAS,NMODE,M,1)
        END IF
        IF(IT.EQ.JCONT(1,I))THEN
          NS1=NS1+1
          IS1(NS1)=5
          JS1(NS1)=IT
          KS1(NS1)=MAXBFN(MAXBAS,NMODE,IT,1)
        END IF
      END DO
      NS2=0
      DO I=1,ICONT(2)
        IF(K.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=1
          JS2(NS2)=K
          KS2(NS2)=MAXBFN(MAXBAS,NMODE,K,1)
        END IF
        IF(L.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=2
          JS2(NS2)=L
          KS2(NS2)=MAXBFN(MAXBAS,NMODE,L,1)
        END IF
        IF(N.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=3
          JS2(NS2)=N
          KS2(NS2)=MAXBFN(MAXBAS,NMODE,N,1)
        END IF
        IF(M.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=4
          JS2(NS2)=M
          KS2(NS2)=MAXBFN(MAXBAS,NMODE,M,1)
        END IF
        IF(IT.EQ.JCONT(2,I))THEN
          NS2=NS2+1
          IS2(NS2)=5
          JS2(NS2)=IT
          KS2(NS2)=MAXBFN(MAXBAS,NMODE,IT,1)
        END IF
      END DO
C**NOW PUT IN SCHEME 1 BASIS
      NSIZE1=1
      NMAX=MAXSUM(1,1)+1
      DO I=1,NS1
        IND=5
        IF(I.GT.1)IND=0
        CALL GETJC1(IP5,ISIZE5,NMODE,IS1(I),KS1(I),NMAX,NSIZE1,
     1  MAXBAS,IND,JS1(I),2)
      END DO
      IF(NS1.EQ.1)GO TO 7777
      NMAX=MAXSUM(2,1)+2
      DO I=1,NS1
        DO J=1,I-1
          CALL GETJC2(IP5,ISIZE5,NMODE,IS1(I),IS1(J),KS1(I),KS1(J),
     1    NMAX,NSIZE1,MAXBAS,0,JS1(I),JS1(J))
        END DO
      END DO
      IF(NS1.EQ.2)GO TO 7777
      NMAX=MAXSUM(3,1)+3
      DO I=1,NS1
        DO J=1,I-1
          DO KK=1,J-1
            CALL GETJC3(IP5,ISIZE5,NMODE,IS1(I),IS1(J),IS1(KK),
     1      KS1(I),KS1(J),KS1(KK),NMAX,NSIZE1,MAXBAS,0,JS1(I),JS1(J),
     2      JS1(KK))
          END DO
        END DO
      END DO
      IF(NS1.EQ.3)GO TO 7777
      NMAX=MAXSUM(4,1)+4
      DO I=1,NS1
        DO J=1,I-1
          DO KK=1,J-1
            DO LL=1,KK-1
              CALL GETJC4(IP5,ISIZE5,NMODE,IS1(I),IS1(J),IS1(KK),
     1        IS1(LL),KS1(I),KS1(J),KS1(KK),KS1(LL),NMAX,NSIZE1,MAXBAS,
     2        0,JS1(I),JS1(J),JS1(KK),JS1(LL))
            END DO
          END DO
        END DO
      END DO
7777  CONTINUE
C**NOW GET SIZE OF SCHEME 2 BASIS
      NSIZE2=1
      NMAX=MAXSUM(1,2)+1
      DO I=1,NS2
        CALL GETJC1(IP5,ISIZE5,NMODE,IS2(I),KS2(I),NMAX,NSIZE2,
     1  MAXBAS,-5,JS2(I),2)
      END DO
      IF(NS2.EQ.1)GO TO 8888
      NMAX=MAXSUM(2,2)+2
      DO I=1,NS2
        DO J=1,I-1
          CALL GETJC2(IP5,ISIZE5,NMODE,IS2(I),IS2(J),KS2(I),KS2(J),
     1    NMAX,NSIZE2,MAXBAS,-5,JS2(I),JS2(J))
        END DO
      END DO
      IF(NS2.EQ.2)GO TO 8888
      NMAX=MAXSUM(3,2)+3
      DO I=1,NS2
        DO J=1,I-1
          DO KK=1,J-1
            CALL GETJC3(IP5,ISIZE5,NMODE,IS2(I),IS2(J),IS2(KK),KS2(I),
     1      KS2(J),KS2(KK),NMAX,NSIZE2,MAXBAS,-5,JS2(I),JS2(J),JS2(KK))
          END DO
        END DO
      END DO
      IF(NS2.EQ.3)GO TO 8888
      NMAX=MAXSUM(4,2)+4
      DO I=1,NS2
        DO J=1,I-1
          DO KK=1,J-1
            DO LL=1,KK-1
              CALL GETJC4(IP5,ISIZE5,NMODE,IS2(I),IS2(J),IS2(KK),
     1        IS2(LL),KS2(I),KS2(J),KS2(KK),KS2(LL),NMAX,NSIZE2,MAXBAS,
     2        -5,JS2(I),JS2(J),JS2(KK),JS2(LL))
            END DO
          END DO
        END DO
      END DO
8888  CONTINUE
C**TOTAL SIZE IS DIRECT PRODUCT
      NSIZE5=NSIZE1*NSIZE2
C**DUPLICATE SCHEME 1 BASIS FOR ALL SCHEME 2
      KK=NSIZE5+1
      DO I=1,NSIZE1
        J=NSIZE1+1-I
        DO MM=1,NSIZE2
          KK=KK-1
          DO NN=1,NS1
          IP5(KK,IS1(NN))=IP5(J,IS1(NN))
          END DO
        END DO
      END DO
C**NOW FILL IN SCHEME 2 BASIS
      NSIZE2=0
      DO NN=1,NSIZE1
        NMAX=MAXSUM(1,2)+1
        DO I=1,NS2
          IND=1
          IF(I.GT.1)IND=2
          CALL GETJC1(IP5,ISIZE5,NMODE,IS2(I),KS2(I),NMAX,NSIZE2,
     1    MAXBAS,0,JS2(I),IND)
        END DO
        IF(NS2.EQ.1)GO TO 9999
        NMAX=MAXSUM(2,2)+2
        DO I=1,NS2
          DO J=1,I-1
            CALL GETJC2(IP5,ISIZE5,NMODE,IS2(I),IS2(J),KS2(I),KS2(J),
     1      NMAX,NSIZE2,MAXBAS,0,JS2(I),JS2(J))
          END DO
        END DO
        IF(NS2.EQ.2)GO TO 9999
        NMAX=MAXSUM(3,2)+3
        DO I=1,NS2
          DO J=1,I-1
            DO KK=1,J-1
              CALL GETJC3(IP5,ISIZE5,NMODE,IS2(I),IS2(J),IS2(KK),
     1        KS2(I),KS2(J),KS2(KK),NMAX,NSIZE2,MAXBAS,0,JS2(I),JS2(J),
     2        JS2(KK))
            END DO
          END DO
        END DO
        IF(NS2.EQ.3)GO TO 9999
        NMAX=MAXSUM(4,2)+4
        DO I=1,NS2
          DO J=1,I-1
            DO KK=1,J-1
              DO LL=1,KK-1
                CALL GETJC4(IP5,ISIZE5,NMODE,IS2(I),IS2(J),IS2(KK),
     1          IS2(LL),KS2(I),KS2(J),KS2(KK),KS2(LL),NMAX,NSIZE2,
     2          MAXBAS,0,JS2(I),JS2(J),JS2(KK),JS2(LL))
              END DO
            END DO
          END DO
        END DO
9999    CONTINUE
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETJC1(JP1,JSIZE1,NMODE,NVMOD1,JCI,NMAX,ITOT1,
     1MAXBAS,NVMODE,KK,ISTART)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JP1(JSIZE1,NVMODE),MAXBAS(NMODE,1)
      COMMON/FILASS/IOUT,INP
C     ISTART=2
      IF(NVMODE.GT.0)THEN
C**ZERO POINT BASIS
        DO J=1,NVMODE
          DO I=1,JSIZE1
            JP1(I,J)=1
          END DO
        END DO
      ELSE
C       IF(NVMODE.EQ.0)ISTART=1
      END IF
      DO I1=ISTART,JCI
        IF(I1.GT.MAXBAS(KK,1))GO TO 100
        IF(I1.GT.NMAX)GO TO 100
        ITOT1=ITOT1+1
        IF(NVMODE.GE.0)JP1(ITOT1,NVMOD1)=I1
100     CONTINUE
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETJC2(JP2,JSIZE2,NMODE,NVMOD1,NVMOD2,JCI1,JCI2,
     1NMAX,ITOT2,MAXBAS,NVMODE,KK,LL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JP2(JSIZE2,NVMODE),MAXBAS(NMODE,2)
      DO I1=2,JCI1
        IF(I1.GT.MAXBAS(KK,2))GO TO 101
        DO I2=2,JCI2
          IF(I2.GT.MAXBAS(LL,2))GO TO 100
          IF(I1+I2.GT.NMAX)GO TO 100
          ITOT2=ITOT2+1
          IF(NVMODE.GE.0)THEN
            JP2(ITOT2,NVMOD1)=I1
            JP2(ITOT2,NVMOD2)=I2
          END IF
100       CONTINUE
        END DO
101     CONTINUE
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETJC3(JP3,JSIZE3,NMODE,NVMOD1,NVMOD2,NVMOD3,
     1JCI1,JCI2,JCI3,NMAX,ITOT3,MAXBAS,NVMODE,KK,LL,NN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JP3(JSIZE3,NVMODE),MAXBAS(NMODE,3)
      DO I1=2,JCI1
        IF(I1.GT.MAXBAS(KK,3))GO TO 102
        DO I2=2,JCI2
          IF(I2.GT.MAXBAS(LL,3))GO TO 101
          DO I3=2,JCI3
            IF(I3.GT.MAXBAS(NN,3))GO TO 100
            IF(I1+I2+I3.GT.NMAX)GO TO 100
            ITOT3=ITOT3+1
            IF(NVMODE.GE.0)THEN
              JP3(ITOT3,NVMOD1)=I1
              JP3(ITOT3,NVMOD2)=I2
              JP3(ITOT3,NVMOD3)=I3
            END IF
100         CONTINUE
          END DO
101       CONTINUE
        END DO
102     CONTINUE
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETJC4(JP4,JSIZE4,NMODE,NVMOD1,NVMOD2,NVMOD3,
     1NVMOD4,JCI1,JCI2,JCI3,JCI4,NMAX,ITOT4,MAXBAS,NVMODE,KK,LL,NN,MM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JP4(JSIZE4,NVMODE),MAXBAS(NMODE,4)
      DO I1=2,JCI1
        IF(I1.GT.MAXBAS(KK,4))GO TO 103
        DO I2=2,JCI2
          IF(I2.GT.MAXBAS(LL,4))GO TO 102
          DO I3=2,JCI3
            IF(I3.GT.MAXBAS(NN,4))GO TO 101
            DO I4=2,JCI4
              IF(I4.GT.MAXBAS(MM,4))GO TO 100
              IF(I1+I2+I3+I4.GT.NMAX)GO TO 100
              ITOT4=ITOT4+1
              IF(NVMODE.GE.0)THEN
                JP4(ITOT4,NVMOD1)=I1
                JP4(ITOT4,NVMOD2)=I2
                JP4(ITOT4,NVMOD3)=I3
                JP4(ITOT4,NVMOD4)=I4
              END IF
100           CONTINUE
            END DO
101         CONTINUE
          END DO
102       CONTINUE
        END DO
103     CONTINUE
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETCP(IP,ISIZE,JP1,JP2,JP3,JP4,JP5,JP6,JSIZE1,JSIZE2,
     1JSIZE3,JSIZE4,JSIZE5,JSIZE6,JP,JSIZE,ITOT,JTOT1,JTOT2,JTOT3,
     2JTOT4,JTOT5,JTOT6,JTOT,ICI,JSTART,INC,KBAS,MAXBAS,NMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IP(ISIZE,1),JP1(JSIZE1,1),JP2(JSIZE2,1),JP3(JSIZE3,1)
      DIMENSION JP4(JSIZE4,1),JP5(JSIZE5,1),JP6(JSIZE6,1),JP(JSIZE,1)
      DIMENSION MAXBAS(NMODE,1)
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
C**TEMPORARY (DIMENSIONS)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
C************************
      COMMON/FILASS/IOUT,INP
C**WE PUT IN (NMOD)-MODE FOR SCHEME (2)
      NMOD=INC+1
      IF(ICI.LT.-INC)THEN
C**ADD IN 0-MODE FOR SCHEME (1) WITH ALL (NMOD) OF SCHEME (2)
C**TOTAL OF 0+NMOD MODES THIS PART (ALWAYS FITS)
        NMOD0=NMOD
        NMAX1=MAXSUM(NMOD0,1)+NVMODE
        NMAX2=MAXSUM(NMOD0,2)+NVMODE
        NMAX=MAX0(NMAX1,NMAX2)
        DO J=1,JTOT
          NSUM=0
          DO I=1,ICONT(2)
            IF(JP(JSTART+J,I)-1.GT.NBAS(NMOD0,2))GO TO 10
            IF(JP(JSTART+J,I).GT.MAXBAS(JCONT(2,I),NMOD0))GO TO 10
            NSUM=NSUM+JP(JSTART+J,I)
          END DO
          IF(NSUM.LE.NMAX)THEN
            ITOT=ITOT+1
            IF(KBAS.GT.0)THEN
              DO I=1,ICONT(2)
                K=JCONT(2,I)
                IP(ITOT,K)=JP(JSTART+J,I)
              END DO
            END IF
          END IF
10        CONTINUE
        END DO
        IF(ICI.LT.-1-INC)THEN
C**ADD IN 1-MODE FOR SCHEME (1) WITH ALL (NMOD) OF SCHEME (2)
C**TOTAL OF 1+NMOD MODES THIS PART
          NMOD1=NMOD+1
          NMAX1=MAXSUM(NMOD1,1)+NVMODE
          NMAX2=MAXSUM(NMOD1,2)+NVMODE
          NMAX=MAX0(NMAX1,NMAX2)
          DO M=2,JTOT1
            DO J=1,JTOT
              NSUM=0
              DO I=1,ICONT(1)
                IF(JP1(M,I)-1.GT.NBAS(NMOD1,1))GO TO 20
                IF(JP1(M,I).GT.MAXBAS(JCONT(1,I),NMOD1))GO TO 20
                NSUM=NSUM+JP1(M,I)
              END DO
              DO I=1,ICONT(2)
                IF(JP(JSTART+J,I)-1.GT.NBAS(NMOD1,2))GO TO 20
                IF(JP(JSTART+J,I).GT.MAXBAS(JCONT(2,I),NMOD1))GO TO 20
                NSUM=NSUM+JP(JSTART+J,I)
              END DO
              IF(NSUM.LE.NMAX)THEN
                ITOT=ITOT+1
                IF(KBAS.GT.0)THEN
                  DO I=1,ICONT(1)
                    K=JCONT(1,I)
                    IP(ITOT,K)=JP1(M,I)
                  END DO
                  DO I=1,ICONT(2)
                    K=JCONT(2,I)
                    IP(ITOT,K)=JP(JSTART+J,I)
                  END DO
                END IF
              END IF
20            CONTINUE
            END DO
          END DO
          IF(ICI.LT.-2-INC)THEN
C**ADD IN 2-MODE FOR SCHEME (1) WITH ALL (NMOD) OF SCHEME (2)
C**TOTAL OF 2+NMOD MODES THIS PART
            NMOD2=NMOD+2
            NMAX1=MAXSUM(NMOD2,1)+NVMODE
            NMAX2=MAXSUM(NMOD2,2)+NVMODE
            NMAX=MAX0(NMAX1,NMAX2)
            DO M=1,JTOT2
              DO J=1,JTOT
                NSUM=0
                DO I=1,ICONT(1)
                  IF(JP2(M,I)-1.GT.NBAS(NMOD2,1))GO TO 30
                  IF(JP2(M,I).GT.MAXBAS(JCONT(1,I),NMOD2))GO TO 30
                  NSUM=NSUM+JP2(M,I)
                END DO
                DO I=1,ICONT(2)
                  IF(JP(JSTART+J,I)-1.GT.NBAS(NMOD2,2))GO TO 30
                  IF(JP(JSTART+J,I).GT.MAXBAS(JCONT(2,I),NMOD2))
     1            GO TO 30
                  NSUM=NSUM+JP(JSTART+J,I)
                END DO
                IF(NSUM.LE.NMAX)THEN
                  ITOT=ITOT+1
                  IF(KBAS.GT.0)THEN
                    DO I=1,ICONT(1)
                      K=JCONT(1,I)
                      IP(ITOT,K)=JP2(M,I)
                    END DO
                    DO I=1,ICONT(2)
                      K=JCONT(2,I)
                      IP(ITOT,K)=JP(JSTART+J,I)
                    END DO
                  END IF
                END IF
30              CONTINUE
              END DO
            END DO
            IF(ICI.LT.-3-INC)THEN
C**ADD IN 3-MODE FOR SCHEME (1) WITH ALL (NMOD) OF SCHEME (2)
C**TOTAL OF 3+NMOD MODES THIS PART
              NMOD3=NMOD+3
              NMAX1=MAXSUM(NMOD3,1)+NVMODE
              NMAX2=MAXSUM(NMOD3,2)+NVMODE
              NMAX=MAX0(NMAX1,NMAX2)
              DO M=1,JTOT3
                DO J=1,JTOT
                  NSUM=0
                  DO I=1,ICONT(1)
                    IF(JP3(M,I)-1.GT.NBAS(NMOD3,1))GO TO 40
                    IF(JP3(M,I).GT.MAXBAS(JCONT(1,I),NMOD3))GO TO 40
                    NSUM=NSUM+JP3(M,I)
                  END DO
                  DO I=1,ICONT(2)
                    IF(JP(JSTART+J,I)-1.GT.NBAS(NMOD3,2))GO TO 40
                    IF(JP(JSTART+J,I).GT.MAXBAS(JCONT(2,I),NMOD3))
     1              GO TO 40
                    NSUM=NSUM+JP(JSTART+J,I)
                  END DO
                  IF(NSUM.LE.NMAX)THEN
                    ITOT=ITOT+1
                    IF(KBAS.GT.0)THEN
                      DO I=1,ICONT(1)
                        K=JCONT(1,I)
                        IP(ITOT,K)=JP3(M,I)
                      END DO
                      DO I=1,ICONT(2)
                        K=JCONT(2,I)
                        IP(ITOT,K)=JP(JSTART+J,I)
                      END DO
                    END IF
                  END IF
40                CONTINUE
                END DO
              END DO
              IF(ICI.LT.-4-INC)THEN
C**ADD IN 4-MODE FOR SCHEME (1) WITH ALL (NMOD) OF SCHEME (2)
C**TOTAL OF 4+NMOD MODES THIS PART
                NMOD4=NMOD+4
                NMAX1=MAXSUM(NMOD4,1)+NVMODE
                NMAX2=MAXSUM(NMOD4,2)+NVMODE
                NMAX=MAX0(NMAX1,NMAX2)
                DO M=1,JTOT4
                  DO J=1,JTOT
                    NSUM=0
                    DO I=1,ICONT(1)
                      IF(JP4(M,I)-1.GT.NBAS(NMOD4,1))GO TO 50
                      IF(JP4(M,I).GT.MAXBAS(JCONT(1,I),NMOD4))GO TO 50
                      NSUM=NSUM+JP4(M,I)
                    END DO
                    DO I=1,ICONT(2)
                      IF(JP(JSTART+J,I)-1.GT.NBAS(NMOD4,2))GO TO 50
                      IF(JP(JSTART+J,I).GT.MAXBAS(JCONT(2,I),NMOD4))
     1                GO TO 50
                      NSUM=NSUM+JP(JSTART+J,I)
                    END DO
                    IF(NSUM.LE.NMAX)THEN
                      ITOT=ITOT+1
                      IF(KBAS.GT.0)THEN
                        DO I=1,ICONT(1)
                          K=JCONT(1,I)
                          IP(ITOT,K)=JP4(M,I)
                        END DO
                        DO I=1,ICONT(2)
                          K=JCONT(2,I)
                          IP(ITOT,K)=JP(JSTART+J,I)
                        END DO
                      END IF
                    END IF
50                  CONTINUE
                  END DO
                END DO
                IF(ICI.LT.-5-INC)THEN
C**ADD IN 5-MODE FOR SCHEME (1) WITH ALL (NMOD) OF SCHEME (2)
C**TOTAL OF 5+NMOD MODES THIS PART
                  NMOD5=NMOD+5
                  NMAX1=MAXSUM(NMOD5,1)+NVMODE
                  NMAX2=MAXSUM(NMOD5,2)+NVMODE
                  NMAX=MAX0(NMAX1,NMAX2)
                  DO M=1,JTOT5
                    DO J=1,JTOT
                      NSUM=0
                      DO I=1,ICONT(1)
                        IF(JP5(M,I)-1.GT.NBAS(NMOD5,1))GO TO 60
                        IF(JP5(M,I).GT.MAXBAS(JCONT(1,I),NMOD5))
     1                  GO TO 60
                        NSUM=NSUM+JP5(M,I)
                      END DO
                      DO I=1,ICONT(2)
                        IF(JP(JSTART+J,I)-1.GT.NBAS(NMOD5,2))GO TO 60
                        IF(JP(JSTART+J,I).GT.MAXBAS(JCONT(2,I),NMOD5))
     1                  GO TO 60
                        NSUM=NSUM+JP(JSTART+J,I)
                      END DO
                      IF(NSUM.LE.NMAX)THEN
                        ITOT=ITOT+1
                        IF(KBAS.GT.0)THEN
                          DO I=1,ICONT(1)
                            K=JCONT(1,I)
                            IP(ITOT,K)=JP5(M,I)
                          END DO
                          DO I=1,ICONT(2)
                            K=JCONT(2,I)
                            IP(ITOT,K)=JP(JSTART+J,I)
                          END DO
                        END IF
                      END IF
60                    CONTINUE
                    END DO
                  END DO
                  IF(ICI.LT.-6-INC)THEN
C**ADD IN 6-MODE FOR SCHEME (1) WITH ALL (NMOD) OF SCHEME (2)
C**TOTAL OF 6+NMOD MODES THIS PART
                    NMOD6=NMOD+6
                    NMAX1=MAXSUM(NMOD6,1)+NVMODE
                    NMAX2=MAXSUM(NMOD6,2)+NVMODE
                    NMAX=MAX0(NMAX1,NMAX2)
                    DO M=1,JTOT6
                      DO J=1,JTOT
                        NSUM=0
                        DO I=1,ICONT(1)
                          IF(JP6(M,I)-1.GT.NBAS(NMOD6,1))GO TO 70
                          IF(JP6(M,I).GT.MAXBAS(JCONT(1,I),NMOD6))
     1                    GO TO 70
                          NSUM=NSUM+JP6(M,I)
                        END DO
                        DO I=1,ICONT(2)
                          IF(JP(JSTART+J,I)-1.GT.NBAS(NMOD6,2))GO TO 70
                          IF(JP(JSTART+J,I).GT.
     1                    MAXBAS(JCONT(2,I),NMOD6))GO TO 70
                          NSUM=NSUM+JP(JSTART+J,I)
                        END DO
                        IF(NSUM.LE.NMAX)THEN
                          ITOT=ITOT+1
                          IF(KBAS.GT.0)THEN
                            DO I=1,ICONT(1)
                              K=JCONT(1,I)
                              IP(ITOT,K)=JP6(M,I)
                            END DO
                            DO I=1,ICONT(2)
                              K=JCONT(2,I)
                              IP(ITOT,K)=JP(JSTART+J,I)
                            END DO
                          END IF
                        END IF
70                      CONTINUE
                      END DO
                    END DO
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE VCCI0(IP,ISIZMX,XA,KEL,NS,EVCI,NVALX,KEL21,NVSYMX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      DIMENSION IP(ISIZMX,2)
      DIMENSION XA(1)
      DIMENSION EVCI(NVALX,KEL21,NVSYMX,2)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFF=IOFF+NCSIZE(I)
      END DO
      J0R=IOFF*(IOFF+1)/2
C**TEMPORARY-LAN
      JS=1
C**TEMPORARY-LAN
C     DO 9999 ISM1=1,NVSYM
      DO 9999 ISM1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
C**TEMPORARY-LAN
        J0=(IROFF)*(IROFF-1)/2+IOFF-J0R
C**TEMPORARY-LAN
C**RHS INDEX FOR SCHEME 1
        NR=IP(IROFF,1)
C**RHS INDEX FOR SCHEME 2
        MR=IP(IROFF,2)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**LHS INDEX FOR SCHEME 1
          NL=IP(ILOFF,1)
C**LHS INDEX FOR SCHEME 2
          ML=IP(ILOFF,2)
C**OVERLAP OF CONTRACTION FUNCTIONS IN SCHEME 1
          IF(NR.NE.NL)GO TO 4000
C**OVERLAP OF CONTRACTION FUNCTIONS IN SCHEME 2
          IF(MR.NE.ML)GO TO 4000
C**GET MATRIX ELEMENT
          XYZ=EVCI(NL,KEL,ISM1,1)+EVCI(ML,KEL,ISM2,2)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE VCCIKE(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MM,IP,ISIZMX,
     1IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,NSIZE,XK,
     2TEMP,XCON,NVAL,ISTART,IEND,MODINT,OMEGA,KEL,NS,CFS,ISIZXX,NVALX,
     3KEL21,NVSYMX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*4 XRA1(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H(NN,MM,3),XQ(MM),IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP1(ISIZE1,1),TEMP(ICSIZE,NVAL)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA1(1)
      DIMENSION XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
C**ZERO ORDER KINETIC ENERGY TERM
      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VCCIKE'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1'
C**IF IT IS IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....IT IS NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0CIKE
      IF(ICONDP.NE.0)GO TO 6666
C**ZERO ORDER KINETIC ENERGY INTEGRAL FOR 'ACTIVE' SCHEME
      MD=MODINT(MOD1)
      DO M=1,MM/MD
        Q=XQ(M)
        IF(IWHICH.EQ.0)THEN
          VHARM=OMEGA*OMEGA*Q*Q/2
        ELSE
          VHARM=0.D0
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          NR=IP1(IRHS,1)
          X=(-H(NR,M,3)/2+VHARM*H(NR,M,1))*MD
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP1(ILHS,1)
C           MULT=1-MOD(IABS(NR-NL),MD)
            Y=H(NL,M,1)
C           XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X*MULT
            XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X
          END DO
        END DO
      END DO
      CALL MATOUT(XA1,XRA1,NSIZE,21)
      GO TO 7777
6666  CALL MATIN(XA1,XRA1,NSIZE,21)
7777  CONTINUE

C***********************************ALGORITHM FROM VCIKE

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR=IPC(IROFF,MODE1)
C**FIND RHS INDEX (TRIVIAL CASE)
        DO IR=1,NSIZE
          IF(NR.EQ.IP1(IR,1))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.(IPC(IROFF,K).NE.IPC(ILOFF,K)))
     1      IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL=IPC(ILOFF,MODE1)
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
          XYZ=XA1(I)
3000      CONTINUE
          XK(ILHS,IRHS)=XYZ*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
        J0=(IROFF)*(IROFF-1)/2+IOFF
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IP(IROFF,KCONT)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IP(IROFF,MCONT).NE.IP(ILOFF,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IP(ILOFF,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE VCCIKC(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MM,IP,ISIZMX,
     1IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,NSIZE,XK,
     2TEMP,IRL,XCON,NVAL,ISTART,IEND,MODINT,OMEGA,KEL,NS,CFS,ISIZXX,
     3NVALX,KEL21,NVSYMX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*4 XRA1(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H(NN,MM,3),XQ(MM),IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP1(ISIZE1,1),TEMP(ICSIZE,NVAL)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),IRL(ICSIZE),XA1(1)
      DIMENSION XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
C**ZERO ORDER KINETIC ENERGY TERM
      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'Calculating VCCIKC'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1'
C**IF IT IS IN SCHEME 1.....
      NKMODE=ICONT(1)
      NLMODE=ICONT(2)
      MCONT=2
C**....IT IS NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        NLMODE=ICONT(1)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0CIKE
      IF(ICONDP.NE.0)GO TO 6666
C**ZERO ORDER KINETIC ENERGY INTEGRAL FOR 'ACTIVE' SCHEME
      MD=MODINT(MOD1)
      DO M=1,MM/MD
        Q=XQ(M)
        IF(IWHICH.EQ.0)THEN
          VHARM=OMEGA*OMEGA*Q*Q/2
        ELSE
          VHARM=0.D0
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          NR=IP1(IRHS,1)
          X=(-H(NR,M,3)/2+VHARM*H(NR,M,1))*MD
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP1(ILHS,1)
C           MULT=1-MOD(IABS(NR-NL),MD)
            Y=H(NL,M,1)
C           XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X*MULT
            XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X
          END DO
        END DO
      END DO
      CALL MATOUT(XA1,XRA1,NSIZE,21)
      GO TO 7777
6666  CALL MATIN(XA1,XRA1,NSIZE,21)
7777  CONTINUE

C***********************************ALGORITHM FROM VCIKE

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        KCOFF=ICOFF1
        LCOFF=ICOFF2
        KCSZ=ICSZ1
        LCSZ=ICSZ2
        ICNTZ=ICCONT(1)
        NKVAL=NCVAL(1,ISM1)
        NLVAL=NCVAL(2,ISM2)
      ELSE
        ISM=ISM2
        KCOFF=ICOFF2
        LCOFF=ICOFF1
        KCSZ=ICSZ2
        LCSZ=ICSZ1
        ICNTZ=ICCONT(2)
        NKVAL=NCVAL(2,ISM2)
        NLVAL=NCVAL(1,ISM1)
      END IF

C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

C
      DO IKRL=1,KCSZ
        IOFFRL=IKRL+KCOFF
        NRL=IPC(IOFFRL,MODE1)
        DO IR=1,NSIZE
          IF(NRL.EQ.IP1(IR,1))GO TO 1000
        END DO
1000    IRL(IKRL)=IR
      END DO
C
C**CASE 1..INTEGRATION MODES IN CONTRACTED SCHEME.
C**      ..OVERLAP OF NON-CONTRACTED BASIS.
      IF(ICNTZ.GT.0)THEN
        DO IRHS=1,KCSZ
          IOFFR=IRHS+KCOFF
          IR=IRL(IRHS)
          DO ILHS=1,IRHS
            IOFFL=ILHS+KCOFF
C**OVERLAP OF REMAINING STATES
            IS=1
            DO K=1,NKMODE
              IF(IS.EQ.0)GO TO 2000
              IF(K.NE.MODE1.AND.(IPC(IOFFR,K).NE.IPC(IOFFL,K)))
     1        IS=0
            END DO
C**OVERLAP OF REMAINING STATES
            IL=IRL(ILHS)
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA1(I)
2000        CONTINUE
            XK(ILHS,IRHS)=XYZ*IS
            XK(IRHS,ILHS)=XK(ILHS,IRHS)
          END DO
        END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',KCSZ,NKVAL,KCSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
        CALL DGEMM('T','N',NKVAL,NKVAL,KCSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

        IRHS=0
        DO IKR=1,NKVAL
          DO ILR=1,LCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            DO IKL=1,IKR
              ILHS=(IKL-1)*LCSZ+ILR
              XYZ=XCON(IKL,IKR)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ
            END DO
          END DO
        END DO

C**CASE 2..INTEGRATION MODES IN NON-CONTRACTED SCHEME.
C**      ..OVERLAP OF CONTRACTED BASIS.
      ELSE
        IRHS=0
        DO ILR=1,NLVAL
          DO IKR=1,KCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            IOFFR=IKR+KCOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            IR=IRL(IKR)
            DO IKL=1,IKR
              ILHS=(ILR-1)*KCSZ+IKL
              IOFFL=IKL+KCOFF
C**OVERLAP OF REMAINING STATES
              IS=1
              DO K=1,NKMODE
                IF(IS.EQ.0)GO TO 3000
                IF(K.NE.MODE1.AND.(IPC(IOFFR,K).NE.IPC(IOFFL,K)))
     1          IS=0
              END DO
C**OVERLAP OF REMAINING STATES
              IL=IRL(IKL)
C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
              XYZ=XA1(I)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ*IS
3000          CONTINUE
            END DO
          END DO
        END DO
      END IF

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCV0(NMODE,NNMODE,ITMODE,HTAU,XQTAU,NNTAU,MMTAU,IP,
     1ISIZMX,IPC,ICSIZE,IPSIZE,ITCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,
     2NSIZE,XK,
     2TEMP,XCON,NVAL,ISTART,IEND,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,
     3CFS,ISIZXX,NVALX,KEL21,NVSYMX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP,VC(3),VR(J21)
      REAL*4 VPR,VCR(3),VRR(J21)
      REAL*4 XRA1(1)
      DIMENSION X(3)
      DIMENSION IPC(IPSIZE,1),IP1(ISIZE1,1),TEMP(ICSIZE,NVAL)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA1(1),XCON(NVAL,NVAL)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU),MODINT(NMODE)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      CALL VDCV0(NMODE,NNMODE,ITMODE,HTAU,XQTAU,NNTAU,MMTAU,
     1IP,ISIZMX,IPC,ICSIZE,IPSIZE,ITCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,
     2NSIZE,
     2XK,TEMP,XCON,NVAL,ISTART,IEND,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,
     3NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCV0(NMODE,NNMODE,ITMODE,HTAU,XQTAU,NNTAU,
     1MMTAU,IP,ISIZMX,IPC,ICSIZE,IPSIZE,ITCONT,XA,ISIZE,IP1,ISIZE1,XA1,
     2XRA1,
     2NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,VC,VCR,VR,VRR,J21,KROT,
     3MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP,VC(3),VR(J21)
      REAL*4 VPR,VCR(3),VRR(J21)
      REAL*4 XRA1(1)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      DIMENSION X(3)
      DIMENSION IPC(IPSIZE,1),IP1(ISIZE1,1),TEMP(ICSIZE,NVAL)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA1(1),XCON(NVAL,NVAL)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU),MODINT(NMODE)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VCCV0'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'TAU'
C**IF IT IS IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....IT IS NOT IN SCHEME 2
      IF(ITCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0CV0
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=JNTFAC(NMODE,ICOUPL,0)
      IFACTC=JNTFAC(NMODE,ICOUPC,0)
C**INTRINSIC
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1

      KA=KROT/2
      LMAX=1+MOD(KA,2)
      FACTRC=0.D0
      IF(J21.GT.1)FACTRC=IFACTC
      MDT=MODINT(NSMODE)
C**LOOP ROUND TAU
      ITAU=INIT-INCTAU
      DO MTAU=1,MMTAU/MDT
        ITAU=ITAU+INCTAU
CCCC    IF(ITAU.GT.362)ITAU=ITAU-360
        IF(ITAU.GT.722)ITAU=ITAU-720

C***********************************************************

        IF(JCOUPC.GE.0)THEN
          IF(J21.GT.1.AND.ICOUPC.GE.0)READ(61)VR
          IF(ICOUPC.GE.0)READ(81)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPC.GE.0)READ(61)VRR
          IF(ICOUPC.GE.0)READ(81)VCR
        END IF
        IF(JCOUPL.GE.0)THEN
          READ(71)VP
        ELSE
          READ(71)VPR
        END IF

C***********************************************************

        IF(JCOUPL.GE.0)THEN
          VV=VP*IFACTL*DSTAU(ITAU)
          IF(ICOUPC.GE.0)THEN
            VV=VV+VR(KROT)*FACTRC*DSTAU(ITAU)
            DO I=1,3
              X(I)=VC(I)*IFACTC*DSTAU(ITAU)
            END DO
          END IF
          X(3)=X(3)+VV
        ELSE
          VV=VPR*IFACTL*DSTAU(ITAU)
          IF(ICOUPC.GE.0)THEN
            V=VV+VRR(KROT)*FACTRC*DSTAU(ITAU)
            DO I=1,3
              X(I)=VCR(I)*IFACTC*DSTAU(ITAU)
            END DO
          END IF
          X(3)=X(3)+VV
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
C**EVEN TERMS
        DO IRHS=1,NSIZE
          IR=IP1(IRHS,1)+1-MOD(KA,2)
          X1E=HTAU(IR,MTAU,1,LMAX)*X(3)*MDT
          X2E=HTAU(IR,MTAU,2,LMAX)*X(2)*MDT
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            IL=IP1(ILHS,1)+1-MOD(KA,2)
            MULTM=1-MOD(IABS(IR-IL),MDT)
CCCC        MULTM=1-MOD(IABS(IR-IL),2)
            Y1=HTAU(IL,MTAU,1,LMAX)
            Y2=HTAU(IL,MTAU,2,LMAX)
            XA1(ILHS+J0)=XA1(ILHS+J0)+(Y1*X1E+Y2*X2E)*MULTM
          END DO
        END DO
C**ODD TERMS
        DO IRHS=1,NSIZE
          IR=IP1(IRHS,1)+1-MOD(KA,2)
          X1O=HTAU(IR,MTAU,2,LMAX)*X(1)*MDT
          X2O=HTAU(IR,MTAU,1,LMAX)*X(1)*MDT
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            IL=IP1(ILHS,1)-MOD(KA,2)
            MULTM=1-MOD(IABS(IR-IL),MDT)
CCCC        MULTM=1-MOD(IABS(IR-IL),2)
            Y1=HTAU(IL,MTAU,1,LMAX)
            Y2=HTAU(IL,MTAU,2,LMAX)
            XA1(ILHS+J0)=XA1(ILHS+J0)+(Y1*X1O+Y2*X2O)*MULTM
          END DO
        END DO
      END DO
      CALL MATOUT(XA1,XRA1,NSIZE,21)
      GO TO 7777
6666  CALL MATIN(XA1,XRA1,NSIZE,21)
7777  CONTINUE

C***********************************ALGORITHM FROM VCV0

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(ITCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
        NLVAL=NCVAL(2,ISM2)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
        NLVAL=NCVAL(1,ISM1)
      END IF

C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NRTAU=IPC(IROFF,ITMODE)
C**FIND RHS INDEX (TRIVIAL)
        DO IR=1,NSIZE
          IF(NRTAU.EQ.IP1(IR,1))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.ITMODE.AND.(IPC(IROFF,K).NE.IPC(ILOFF,K)))
     1      IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NLTAU=IPC(ILOFF,ITMODE)
C**FIND LHS INDEX
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
          XYZ=XA1(I)
3000      CONTINUE
          XK(ILHS,IRHS)=XYZ*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      IF(ITCONT.EQ.1)THEN
C**SCHEME 1 ACTIVE (K POINTS TO 1; L POINTS TO 2)
        KOFF=0
        DO IR=1,NLVAL
C**DIAGONAL BLOCKS IR ALL HAVE NON-ZERO OVERLAP IN L
          LROFF=0
          DO IRHS=1,NKVAL
            IROFF=1+LROFF+KOFF+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            NR=IP(IROFF,ITCONT)
            LLOFF=0
            DO ILHS=1,IRHS
              ILOFF=1+LLOFF+KOFF
              NL=IP(ILOFF+IOFF,ITCONT)
              XYZ=XCON(NL,NR)
              XA(ILOFF+J0)=XA(ILOFF+J0)+XYZ
              LLOFF=LLOFF+NLVAL
            END DO
            LROFF=LROFF+NLVAL
          END DO
          KOFF=KOFF+1
        END DO
      ELSE
C**SCHEME 2 ACTIVE (L POINTS TO 1; K POINTS TO 2)
        LOFF=0
        DO IR=1,NLVAL
C**DIAGONAL BLOCKS IR ALL HAVE NON-ZERO OVERLAP IN L
          DO IRHS=1,NKVAL
            IROFF=IRHS+LOFF+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            NR=IP(IROFF,ITCONT)
            DO ILHS=1,IRHS
              ILOFF=ILHS+LOFF
              NL=IP(ILOFF+IOFF,ITCONT)
              XYZ=XCON(NL,NR)
              XA(ILOFF+J0)=XA(ILOFF+J0)+XYZ
            END DO
          END DO
          LOFF=LOFF+NKVAL
        END DO
      END IF

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI1(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MM,IP,ISIZMX,
     1IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,NSIZE,XK,
     2TEMP,XCON,NVAL,ISTART,IEND,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     3OMEGA,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,XKAN,MAXQU,MAXPOW,
     4NP1,CP1,JP1,NTOT1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM),VC(MM),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM),VRR(J21,MM)
      REAL*4 XRA1(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H(NN,MM,3),XQ(MM),IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP1(ISIZE1,1),TEMP(ICSIZE,NVAL)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA1(1)
      DIMENSION XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),JP1(MMX1,NTOT1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      MD=MODINT(MOD1)
      CALL VDCCI1(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MM,MM/MD,IP,ISIZMX,
     1IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,NSIZE,XK,
     2TEMP,XCON,NVAL,ISTART,IEND,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     3OMEGA,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,XKAN,MAXQU,MAXPOW,
     4NP1,CP1,JP1,NTOT1,MMX1)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCCI1(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MH,MM,IP,
     1ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,
     2NSIZE,XK,
     2TEMP,XCON,NVAL,ISTART,IEND,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     3OMEGA,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,XKAN,MAXQU,MAXPOW,
     4NP1,CP1,JP1,NTOT1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM),VC(MM),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM),VRR(J21,MM)
      REAL*4 XRA1(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H(NN,MM,3),XQ(MM),IP(ISIZMX,NNMODE)
      DIMENSION H(NN,MH,3),XQ(MM),IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP1(ISIZE1,1),TEMP(ICSIZE,NVAL)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA1(1)
      DIMENSION XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),JP1(MMX1,NTOT1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM1A.EQ.1)THEN
        WRITE(IOUT,*)'Calculating VCCI1'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1'
C**IF IT IS IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....IT IS NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0CI1
C**TEMPORARY
      IF(ICONDP.NE.0)THEN
        IF(MOLINC.GT.0)GO TO 99999
        GO TO 6666
      END IF
C**TEMPORARY

      IFACTC=INTFAC(NMODE,ICOUPC,1)
      IFACTL=INTFAC(NMODE,ICOUPL,1)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFACTL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(71)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(71)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.0)READ(61)VR
        IF(ICOUPC.GT.0)READ(81)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.0)READ(61)VRR
        IF(ICOUPC.GT.0)READ(81)VCR
      END IF

      IF(MOLINC.GT.0)GO TO 99999

C***********************************************************

C**ONE-MODE COUPLING INTEGRAL FOR 'ACTIVE' SCHEME
      MD=MODINT(MOD1)
CCCC  DO M=1,MM/MD
      DO M=1,MM
        VV=0
        IF(JCOUPL.GT.0)THEN
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)VV=VP(M)*IFACTL
C**WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
          IF(ICOUPC.GT.0)THEN
            VV=VV+VC(M)*IFACTC
            IF(J21.GT.1)TERM=VR(KROT,M)*IFACTC
          END IF
        ELSE
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)VV=VPR(M)*IFACTL
C**WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
          IF(ICOUPC.GT.0)THEN
            VV=VV+VCR(M)*IFACTC
            IF(J21.GT.1)TERM=VRR(KROT,M)*IFACTC
          END IF
        END IF
C**********************************************************
C**KINETIC ENERGY
        Q=XQ(M)
        IF(IWHICH.EQ.0)THEN
          VHARM=OMEGA*OMEGA*Q*Q/2
        ELSE
          VHARM=0.D0
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          NR=IP1(IRHS,1)
          X=(-H(NR,M,3)/2+VHARM*H(NR,M,1))*MD
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP1(ILHS,1)
C           MULT=1-MOD(IABS(NR-NL),MD)
            Y=H(NL,M,1)
C           XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X*MULT
            XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X
          END DO
        END DO
C**********************************************************

C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          NR=IP1(IRHS,1)
          X=(VV+TERM)*H(NR,M,1)*MD
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP1(ILHS,1)
C           MULT=1-MOD(IABS(NR-NL),MD)
            Y=H(NL,M,1)
C           XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X*MULT
            XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X
          END DO
        END DO
      END DO
      CALL MATOUT(XA1,XRA1,NSIZE,21)
      GO TO 7777
6666  CALL MATIN(XA1,XRA1,NSIZE,21)
7777  CONTINUE

C     IF(MOLINC.GT.0)RETURN

C***********************************ALGORITHM FROM VCI1

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFF=IOFF+NCSIZE(I)
      END DO
      J0R=IOFF*(IOFF+1)/2
C**TEMPORARY-LAN
      JS=1
C**TEMPORARY-LAN
C     DO 9999 ISM1=1,NVSYM
      DO 9999 ISM1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR=IPC(IROFF,MODE1)
C**FIND RHS INDEX (TRIVIAL CASE)
        DO IR=1,NSIZE
          IF(NR.EQ.IP1(IR,1))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.(IPC(IROFF,K).NE.IPC(ILOFF,K)))
     1      IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL=IPC(ILOFF,MODE1)
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
          XYZ=XA1(I)
          ZYX=0
          IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
            DO I=1,NP1(MOD1)
              K=JP1(I,MOD1)+1
              ZYX=ZYX+CP1(I,MOD1)*XKAN(NL,NR,K,MOD1)
            END DO
C**ANALYTIC
          END IF
3000      CONTINUE
          XK(ILHS,IRHS)=(XYZ+ZYX)*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
C**TEMPORARY-LAN
        J0=(IROFF)*(IROFF-1)/2+IOFF-J0R
C**TEMPORARY-LAN
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IP(IROFF,KCONT)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IP(IROFF,MCONT).NE.IP(ILOFF,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IP(ILOFF,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
99999 CONTINUE

      IF(ITIM1A.EQ.1)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM1A=2
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI1C(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MM,IP,ISIZMX,
     1IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,NSIZE,XK,
     2TEMP,IRL,XCON,NVAL,ISTART,IEND,VP,VPR,VC,VCR,VR,VRR,J21,KROT,
     3MODINT,OMEGA,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,XKAN,MAXQU,
     4MAXPOW,NP1,CP1,JP1,NTOT1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM),VC(MM),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM),VRR(J21,MM)
      REAL*4 XRA1(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H(NN,MM,3),XQ(MM),IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP1(ISIZE1,1),TEMP(ICSIZE,NVAL)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),IRL(ICSIZE),XA1(1)
      DIMENSION XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),JP1(MMX1,NTOT1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      MD=MODINT(MOD1)
      CALL VDCI1C(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MM,MM/MD,IP,ISIZMX,
     1IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,NSIZE,XK,
     2TEMP,IRL,XCON,NVAL,ISTART,IEND,VP,VPR,VC,VCR,VR,VRR,J21,KROT,
     3MODINT,OMEGA,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,XKAN,MAXQU,
     4MAXPOW,NP1,CP1,JP1,NTOT1,MMX1)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI1C(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MH,MM,IP,
     1ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP1,ISIZE1,XA1,XRA1,
     2NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,IEND,VP,VPR,VC,VCR,VR,VRR,J21,
     3KROT,MODINT,OMEGA,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,XKAN,
     4MAXQU,MAXPOW,NP1,CP1,JP1,NTOT1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM),VC(MM),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM),VRR(J21,MM)
      REAL*4 XRA1(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H(NN,MM,3),XQ(MM),IP(ISIZMX,NNMODE)
      DIMENSION H(NN,MH,3),XQ(MM),IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP1(ISIZE1,1),TEMP(ICSIZE,NVAL)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),IRL(ICSIZE),XA1(1)
      DIMENSION XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP1(NTOT1),CP1(MMX1,NTOT1),JP1(MMX1,NTOT1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VCCI1C'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1'
C**IF IT IS IN SCHEME 1.....
      NKMODE=ICONT(1)
      NLMODE=ICONT(2)
      MCONT=2
C**....IT IS NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        NLMODE=ICONT(1)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0CI1
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,1)
      IFACTL=INTFAC(NMODE,ICOUPL,1)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(71)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(71)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.0)READ(61)VR
        IF(ICOUPC.GT.0)READ(81)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.0)READ(61)VRR
        IF(ICOUPC.GT.0)READ(81)VCR
      END IF

C***********************************************************

C**ONE-MODE COUPLING INTEGRAL FOR 'ACTIVE' SCHEME
      MD=MODINT(MOD1)
CCCC  DO M=1,MM/MD
      DO M=1,MM
        VV=0
        IF(JCOUPL.GT.0)THEN
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)VV=VP(M)*IFACTL
C**WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
          IF(ICOUPC.GT.0)THEN
            VV=VV+VC(M)*IFACTC
            IF(J21.GT.1)TERM=VR(KROT,M)*IFACTC
          END IF
        ELSE
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)VV=VPR(M)*IFACTL
C**WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
          IF(ICOUPC.GT.0)THEN
            VV=VV+VCR(M)*IFACTC
            IF(J21.GT.1)TERM=VRR(KROT,M)*IFACTC
          END IF
        END IF
C**********************************************************
C**KINETIC ENERGY
        Q=XQ(M)
        IF(IWHICH.EQ.0)THEN
          VHARM=OMEGA*OMEGA*Q*Q/2
        ELSE
          VHARM=0.D0
        END IF

C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          NR=IP1(IRHS,1)
          X=(-H(NR,M,3)/2+VHARM*H(NR,M,1))*MD
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP1(ILHS,1)
C           MULT=1-MOD(IABS(NR-NL),MD)
            Y=H(NL,M,1)
C           XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X*MULT
            XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X
          END DO
        END DO
C**********************************************************

C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          NR=IP1(IRHS,1)
          X=(VV+TERM)*H(NR,M,1)*MD
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP1(ILHS,1)
C           MULT=1-MOD(IABS(NR-NL),MD)
            Y=H(NL,M,1)
C           XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X*MULT
            XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X
          END DO
        END DO
      END DO
      CALL MATOUT(XA1,XRA1,NSIZE,21)
      GO TO 7777
6666  CALL MATIN(XA1,XRA1,NSIZE,21)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI1

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        KCOFF=ICOFF1
        LCOFF=ICOFF2
        KCSZ=ICSZ1
        LCSZ=ICSZ2
        ICNTZ=ICCONT(1)
        NKVAL=NCVAL(1,ISM1)
        NLVAL=NCVAL(2,ISM2)
      ELSE
        ISM=ISM2
        KCOFF=ICOFF2
        LCOFF=ICOFF1
        KCSZ=ICSZ2
        LCSZ=ICSZ1
        ICNTZ=ICCONT(2)
        NKVAL=NCVAL(2,ISM2)
        NLVAL=NCVAL(1,ISM1)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

C
      DO IKRL=1,KCSZ
        IOFFRL=IKRL+KCOFF
        NRL=IPC(IOFFRL,MODE1)
        DO IR=1,NSIZE
          IF(NRL.EQ.IP1(IR,1))GO TO 1000
        END DO
1000    IRL(IKRL)=IR
      END DO
C
C**CASE 1..INTEGRATION MODES IN CONTRACTED SCHEME.
C**      ..OVERLAP OF NON-CONTRACTED BASIS.
      IF(ICNTZ.GT.0)THEN
        DO IRHS=1,KCSZ
          IOFFR=IRHS+KCOFF
          IR=IRL(IRHS)
          DO ILHS=1,IRHS
            IOFFL=ILHS+KCOFF
C**OVERLAP OF REMAINING STATES
            IS=1
            DO K=1,NKMODE
              IF(IS.EQ.0)GO TO 2000
              IF(K.NE.MODE1.AND.(IPC(IOFFR,K).NE.IPC(IOFFL,K)))
     1        IS=0
            END DO
C**OVERLAP OF REMAINING STATES
            IL=IRL(ILHS)
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA1(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              DO I=1,NP1(MOD1)
                K=JP1(I,MOD1)+1
                ZYX=ZYX+CP1(I,MOD1)*XKAN(NL,NR,K,MOD1)
              END DO
C**ANALYTIC
            END IF
2000        CONTINUE
            XK(ILHS,IRHS)=(XYZ+ZYX)*IS
            XK(IRHS,ILHS)=XK(ILHS,IRHS)
          END DO
        END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',KCSZ,NKVAL,KCSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,KCSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

        IRHS=0
        DO IKR=1,NKVAL
          DO ILR=1,LCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            DO IKL=1,IKR
              ILHS=(IKL-1)*LCSZ+ILR
              XYZ=XCON(IKL,IKR)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ
            END DO
          END DO
        END DO

C**CASE 2..INTEGRATION MODES IN NON-CONTRACTED SCHEME.
C**      ..OVERLAP OF CONTRACTED BASIS.
      ELSE
        IRHS=0
        DO ILR=1,NLVAL
          DO IKR=1,KCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            IOFFR=IKR+KCOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            IR=IRL(IKR)
            DO IKL=1,IKR
              ILHS=(ILR-1)*KCSZ+IKL
              IOFFL=IKL+KCOFF
C**OVERLAP OF REMAINING STATES
              IS=1
              DO K=1,NKMODE
                IF(IS.EQ.0)GO TO 3000
                IF(K.NE.MODE1.AND.(IPC(IOFFR,K).NE.IPC(IOFFL,K)))
     1          IS=0
              END DO
C**OVERLAP OF REMAINING STATES
              IL=IRL(IKL)
C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
              XYZ=XA1(I)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ*IS
3000          CONTINUE
            END DO
          END DO
        END DO
      END IF

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM1A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCV1A(NMODE,NNMODE,MOD1,MODE1,ITMODE,H,
     1XQ,HTAU,XQTAU,NN,MM,NNTAU,MMTAU,IP,ISIZMX,IPC,ICSIZE,IPSIZE,
     2KCONT,XA,
     2ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,
     3TEMP1,JCI,JCIM,X0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,OMEGA,
     4KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I9)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM),VC(MM,6),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM,6),VRR(J21,MM)
      REAL*4 XRA2(1)
      DIMENSION X(9),Y(9),C(9)
      DIMENSION MODINT(NMODE)
      DIMENSION H(NN,MM,3,1),XQ(MM)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC(IPSIZE,1),IP2(ISIZE2,2)
      DIMENSION XA(1),XA2(1),XK(ICSIZE,ICSIZE)
      DIMENSION TEMP1(I9,JCI,JCI),TEMP(ICSIZE,NVAL),XCON(NVAL,NVAL)
      DIMENSION X0(I9,JCI,JCI,MM),T0(I9,JCIM,JCIM,MMTAU)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      CALL VDCV1A(NMODE,NNMODE,MOD1,MODE1,ITMODE,H,
     1XQ,HTAU,XQTAU,NN,MM,MM/MD1,NNTAU,MMTAU,IP,ISIZMX,IPC,ICSIZE,
     2IPSIZE,KCONT,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,XK,TEMP,XCON,
     3NVAL,ISTART,IEND,TEMP1,JCI,JCIM,X0,T0,VP,VPR,VC,VCR,VR,VRR,J21,
     4KROT,MODINT,OMEGA,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I9)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCV1A(NMODE,NNMODE,MOD1,MODE1,ITMODE,H,XQ,HTAU,XQTAU,
     1NN,MH,MM,NNTAU,MMTAU,IP,ISIZMX,IPC,ICSIZE,IPSIZE,
     2KCONT,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,XK,TEMP,XCON,NVAL,
     3ISTART,IEND,TEMP1,JCI,JCIM,X0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,
     4MODINT,OMEGA,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I9)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM),VC(MM,6),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM,6),VRR(J21,MM)
      REAL*4 XRA2(1)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      DIMENSION X(9),Y(9),C(9)
      DIMENSION MODINT(NMODE)
C     DIMENSION H(NN,MM,3,1),XQ(MM)
      DIMENSION H(NN,MH,3,1),XQ(MM)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC(IPSIZE,1),IP2(ISIZE2,2)
      DIMENSION XA(1),XA2(1),XK(ICSIZE,ICSIZE)
      DIMENSION TEMP1(I9,JCI,JCI),TEMP(ICSIZE,NVAL),XCON(NVAL,NVAL)
      DIMENSION X0(I9,JCI,JCI,MM),T0(I9,JCIM,JCIM,MMTAU)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VCCV1A'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' AND 'TAU'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0CV1
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=JNTFAC(NMODE,ICOUPL,1)
      IFACTC=JNTFAC(NMODE,ICOUPC,1)
C**INTRINSIC
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1

      KA=KROT/2
      LMAX=1+MOD(KA,2)
      FACTRC=0.D0
      IF(J21.GT.1)FACTRC=IFACTC
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      MD=MD1*MDT

C**FORM INDIVIDUAL INTEGRATION TERMS (START)
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          DO NRR=1,JCIM
            NR=NRR+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            X(2)=X(1)
            X(3)=HTAU(NR,MTAU,2,LMAX)*MD
            X(4)=X(1)
            X(5)=X(1)
            X(6)=X(3)
            X(7)=X(4)
            X(8)=HTAU(NR,MTAU,2,LMAX)*MD
            X(9)=X(1)
            DO NLL=1,JCIM
              NL=NLL+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              Y(2)=Y(1)
              Y(3)=Y(1)
              Y(4)=HTAU(NL,MTAU,2,LMAX)
              Y(5)=Y(1)
              Y(6)=Y(1)
              Y(7)=Y(4)
              Y(8)=Y(4)
              Y(9)=Y(1)
              DO K=1,I9
                T0(K,NLL,NRR,MTAU)=Y(10-K)*X(10-K)
              END DO
            END DO
          END DO
        ELSE
          DO NRR=1,JCIM
            NR=2*NRR-MOD(NRR,2)+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            X(2)=HTAU(NR,MTAU,2,LMAX)*MD
            X(3)=X(1)
            X(4)=X(2)
            X(5)=HTAU(NR,MTAU,3,LMAX)*MD
            X(6)=X(1)
            X(7)=HTAU(NR,MTAU,1,LMAX)*FACTRC*MD
            DO NLL=1,JCIM
              NL=2*NLL-MOD(NLL,2)+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO K=1,7
                T0(K,NLL,NRR,MTAU)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M=1,MM/MD1
      DO M=1,MM
        IF(.NOT.LINEAR)THEN
          DO NR1=1,JCI
            X(1)=H(NR1,M,2,1)
            X(2)=H(NR1,M,1,1)
            X(3)=X(2)
            X(4)=X(2)
            X(5)=X(1)
            X(6)=X(2)
            X(7)=X(1)
            X(8)=X(2)
            X(9)=X(2)
            DO NL1=1,JCI
              Y(1)=H(NL1,M,1,1)
              Y(2)=H(NL1,M,2,1)
              Y(3)=Y(1)
              Y(4)=Y(1)
              Y(5)=Y(2)
              Y(6)=Y(2)
              Y(7)=Y(1)
              Y(8)=Y(1)
              Y(9)=Y(1)
              DO K=1,I9
                X0(K,NL1,NR1,M)=Y(10-K)*X(10-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR1=1,JCI
            X(1)=H(NR1,M,2,1)
            X(2)=H(NR1,M,1,1)
            X(3)=H(NR1,M,3,1)
            X(4)=X(1)
            X(5)=X(2)
            X(6)=X(2)
            X(7)=X(2)
            DO NL1=1,JCI
              Y(1)=H(NL1,M,1,1)
              DO K=1,7
                X0(K,NL1,NR1,M)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
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
          IF(J21.GT.1.AND.ICOUPC.GE.1)READ(61)VR
          IF(ICOUPC.GE.1)READ(81)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPC.GE.1)READ(61)VRR
          IF(ICOUPC.GE.1)READ(81)VCR
        END IF
        IF(JCOUPL.GT.0)THEN
          READ(71)VP
        ELSE
          READ(71)VPR
        END IF

C***********************************************************

        IF(.NOT.LINEAR)THEN
          DO IRHS1=1,JCI
            DO ILHS1=1,JCI
              DO K=1,I9
                TEMP1(K,ILHS1,IRHS1)=0.D0
              END DO
            END DO
          END DO
C**START 1-MODE INTEGRATION
C         DO M=1,MM/MD1
          DO M=1,MM
            DO I=1,9
              C(I)=0
            END DO
            IF(JCOUPL.GT.0)THEN
              IF(ICOUPC.GE.1)THEN
C**NO WATSON TERM IF RPH
                C(1)=VC(M,1)*IFACTC
                C(2)=VC(M,1)*IFACTC
                C(3)=VC(M,2)*IFACTC
                C(4)=VC(M,2)*IFACTC
                C(5)=VC(M,3)*IFACTC
                C(6)=VC(M,4)*IFACTC
                C(7)=VC(M,4)*IFACTC
                C(8)=VC(M,5)*IFACTC
                C(9)=VC(M,6)*IFACTC
                IF(J21.GT.1)C(9)=C(9)+VR(KROT,M)*IFACTC
              END IF
              C(9)=C(9)+VP(M)*IFACTL
            ELSE
              IF(ICOUPC.GE.1)THEN
C**NO WATSON TERM IF RPH
                C(1)=VCR(M,1)*IFACTC
                C(2)=VCR(M,1)*IFACTC
                C(3)=VCR(M,2)*IFACTC
                C(4)=VCR(M,2)*IFACTC
                C(5)=VCR(M,3)*IFACTC
                C(6)=VCR(M,4)*IFACTC
                C(7)=VCR(M,4)*IFACTC
                C(8)=VCR(M,5)*IFACTC
                C(9)=VCR(M,6)*IFACTC
                IF(J21.GT.1)C(9)=C(9)+VRR(KROT,M)*IFACTC
              END IF
              C(9)=C(9)+VPR(M)*IFACTL
            END IF
            DO IRHS1=1,JCI
              DO ILHS1=1,JCI
                DO K=1,I9
                  TEMP1(K,ILHS1,IRHS1)=TEMP1(K,ILHS1,IRHS1)+
     1            X0(K,ILHS1,IRHS1,M)*C(10-K)
                END DO
              END DO
            END DO
          END DO
C**END 1-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
          DO IRHS=1,NSIZE
            NR=IP2(IRHS,1)
            IRTAU=IP2(IRHS,2)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL=IP2(ILHS,1)
              ILTAU=IP2(ILHS,2)
              DO K=1,I9
                XA2(ILHS+J0)=XA2(ILHS+J0)+TEMP1(K,NL,NR)*
     1          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
              END DO
            END DO
          END DO
        ELSE
          DO IRHS1=1,JCI
            DO ILHS1=1,JCI
              DO K=1,7
                TEMP1(K,ILHS1,IRHS1)=0.D0
              END DO
            END DO
          END DO
C**START 1-MODE INTEGRATION
C         DO M=1,MM/MD1
          DO M=1,MM
            Q=XQ(M)
            IF(JCOUPL.GT.0)THEN
              DO IRHS1=1,JCI
                DO ILHS1=1,JCI
                  IF(ICOUPC.GE.1)THEN
C**(6) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                    DO K=1,5
                      TEMP1(K,ILHS1,IRHS1)=TEMP1(K,ILHS1,IRHS1)-
     1                X0(K,ILHS1,IRHS1,M)*VC(M,K)*IFACTC
                    END DO
                    TEMP1(4,ILHS1,IRHS1)=TEMP1(4,ILHS1,IRHS1)-
     1              X0(4,ILHS1,IRHS1,M)*VC(M,4)*IFACTC
                    TEMP1(7,ILHS1,IRHS1)=TEMP1(7,ILHS1,IRHS1)+
     1              X0(7,ILHS1,IRHS1,M)*VR(KROT,M)
                  END IF
C**INCLUDE K.E. ONCE FOR TAU
                  TEMP1(5,ILHS1,IRHS1)=TEMP1(5,ILHS1,IRHS1)-
     1            X0(5,ILHS1,IRHS1,M)/(8*Q*Q)
C***************************
                  TEMP1(6,ILHS1,IRHS1)=TEMP1(6,ILHS1,IRHS1)+
     1            X0(6,ILHS1,IRHS1,M)*VP(M)*IFACTL
                END DO
              END DO
            ELSE
              DO IRHS1=1,JCI
                DO ILHS1=1,JCI
                  IF(ICOUPC.GE.1)THEN
C**(6) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                    DO K=1,5
                      TEMP1(K,ILHS1,IRHS1)=TEMP1(K,ILHS1,IRHS1)-
     1                X0(K,ILHS1,IRHS1,M)*VCR(M,K)*IFACTC
                    END DO
                    TEMP1(4,ILHS1,IRHS1)=TEMP1(4,ILHS1,IRHS1)-
     1              X0(4,ILHS1,IRHS1,M)*VCR(M,4)*IFACTC
                    TEMP1(7,ILHS1,IRHS1)=TEMP1(7,ILHS1,IRHS1)+
     2              X0(7,ILHS1,IRHS1,M)*VRR(KROT,M)
                  END IF
C**INCLUDE K.E. ONCE FOR TAU
                  TEMP1(5,ILHS1,IRHS1)=TEMP1(5,ILHS1,IRHS1)-
     1            X0(5,ILHS1,IRHS1,M)/(8*Q*Q)
C***************************
                  TEMP1(6,ILHS1,IRHS1)=TEMP1(6,ILHS1,IRHS1)+
     1            X0(6,ILHS1,IRHS1,M)*VPR(M)*IFACTL
                END DO
              END DO
            END IF

C**********************************************************
C**INCLUDE K.E. ONCE FOR Q
            IF(IWHICH.EQ.0)THEN
              VHARM=OMEGA*OMEGA*Q*Q/2
            ELSE
              VHARM=0.D0
            END IF
            DO IRHS1=1,JCI
              DO ILHS1=1,JCI
                TEMP1(3,ILHS1,IRHS1)=TEMP1(3,ILHS1,IRHS1)-
     1          X0(3,ILHS1,IRHS1,M)/2
                TEMP1(6,ILHS1,IRHS1)=TEMP1(6,ILHS1,IRHS1)+
     1          X0(6,ILHS1,IRHS1,M)*VHARM
              END DO
            END DO
C**********************************************************

          END DO
C**END 1-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
          DO IRHS=1,NSIZE
            NR=IP2(IRHS,1)
            IRRTAU=IP2(IRHS,2)
            IRTAU=IRRTAU/2+MOD(IRRTAU,2)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL=IP2(ILHS,1)
              ILLTAU=IP2(ILHS,2)
              ILTAU=ILLTAU/2+MOD(ILLTAU,2)
              DO K=1,7
                XA2(ILHS+J0)=XA2(ILHS+J0)+TEMP1(K,NL,NR)*
     1          T0(K,ILTAU,IRTAU,MTAU)
              END DO
            END DO
          END DO
        END IF
      END DO
C**END TAU LOOP (2-MODE INTEGRATION)
      CALL MATOUT(XA2,XRA2,NSIZE,21)
      GO TO 7777
6666  CALL MATIN(XA2,XRA2,NSIZE,21)
7777  CONTINUE

C***********************************ALGORITHM FROM VCV1

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
        NLVAL=NCVAL(2,ISM2)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
        NLVAL=NCVAL(1,ISM1)
      END IF

C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR=IPC(IROFF,MODE1)
        NRTAU=IPC(IROFF,ITMODE)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR.EQ.IP2(IR,1).AND.NRTAU.EQ.IP2(IR,2))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.ITMODE.AND.(IPC(IROFF,K).NE.
     1      IPC(ILOFF,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL=IPC(ILOFF,MODE1)
          NLTAU=IPC(ILOFF,ITMODE)
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
          XYZ=XA2(I)
3000      CONTINUE
          XK(ILHS,IRHS)=XYZ*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      IF(KCONT.EQ.1)THEN
C**SCHEME 1 ACTIVE (K POINTS TO 1; L POINTS TO 2)
        KOFF=0
        DO IR=1,NLVAL
C**DIAGONAL BLOCKS IR ALL HAVE NON-ZERO OVERLAP IN L
          LROFF=0
          DO IRHS=1,NKVAL
            IROFF=1+LROFF+KOFF+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            NR=IP(IROFF,KCONT)
            LLOFF=0
            DO ILHS=1,IRHS
              ILOFF=1+LLOFF+KOFF
              NL=IP(ILOFF+IOFF,KCONT)
              XYZ=XCON(NL,NR)
              XA(ILOFF+J0)=XA(ILOFF+J0)+XYZ
              LLOFF=LLOFF+NLVAL
            END DO
            LROFF=LROFF+NLVAL
          END DO
          KOFF=KOFF+1
        END DO
      ELSE
C**SCHEME 2 ACTIVE (L POINTS TO 1; K POINTS TO 2)
        LOFF=0
        DO IR=1,NLVAL
C**DIAGONAL BLOCKS IR ALL HAVE NON-ZERO OVERLAP IN L
          DO IRHS=1,NKVAL
            IROFF=IRHS+LOFF+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            NR=IP(IROFF,KCONT)
            DO ILHS=1,IRHS
              ILOFF=ILHS+LOFF
              NL=IP(ILOFF+IOFF,KCONT)
              XYZ=XCON(NL,NR)
              XA(ILOFF+J0)=XA(ILOFF+J0)+XYZ
            END DO
          END DO
          LOFF=LOFF+NKVAL
        END DO
      END IF

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM1A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCV1B(NMODE,NNMODE,MOD1,MODE1,MODE2,H,XQ,HTAU,XQTAU,
     1NN,MM,NNTAU,MMTAU,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,
     2IPSIZ2,KCONT1,
     2KCONT2,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,IRL,XK2,TEMP,KTEMP,
     3XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,JCI,JCIM,X0,T0,VP,VPR,VC,VCR,
     4VR,VRR,J21,KROT,MODINT,OMEGA,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,
     5NVSYMX,IND,NONZC1,NONZC2,I9,MAXVAL,MAXQ2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM),VC(MM,6),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM,6),VRR(J21,MM)
      REAL*4 XRA2(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H(NN,MM,3,1),XQ(MM)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP2(ISIZE2,2),TEMP(KTEMP,1)
      DIMENSION TEMP1(I9,JCI,JCI)
      DIMENSION X0(I9,JCI,JCI,MM),T0(I9,JCIM,JCIM,MMTAU)
      DIMENSION X(9),Y(9),C(9)
C******************************************************************
      DIMENSION XA(1),XA2(1)
      DIMENSION IRL(MAXQ2,MAXQ2),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C******************************************************************
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      CALL VDCV1B(NMODE,NNMODE,MOD1,MODE1,MODE2,H,XQ,HTAU,XQTAU,
     1NN,MM,MM/MD1,NNTAU,MMTAU,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,
     2ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,
     3IRL,XK2,
     3TEMP,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,JCI,JCIM,X0,T0,VP,
     4VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,OMEGA,KEL,NS,CFS1,CFS2,ISIZXX,
     5NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I9,MAXVAL,MAXQ2)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCV1B(NMODE,NNMODE,MOD1,MODE1,MODE2,H,XQ,HTAU,XQTAU,
     1NN,MH,MM,NNTAU,MMTAU,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,
     2IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,IRL,XK2,
     3TEMP,
     3KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,JCI,JCIM,X0,T0,VP,VPR,
     4VC,VCR,VR,VRR,J21,KROT,MODINT,OMEGA,KEL,NS,CFS1,CFS2,ISIZXX,
     5NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I9,MAXVAL,MAXQ2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM),VC(MM,6),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM,6),VRR(J21,MM)
      REAL*4 XRA2(1)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD

C**TEMPORARY
      DIMENSION IRL(MAXQ2,MAXQ2)
C**TEMPORARY

      DIMENSION MODINT(NMODE)
C     DIMENSION H(NN,MM,3,1),XQ(MM)
      DIMENSION H(NN,MH,3,1),XQ(MM)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP2(ISIZE2,2),TEMP(KTEMP,1)
      DIMENSION TEMP1(I9,JCI,JCI)
      DIMENSION X0(I9,JCI,JCI,MM),T0(I9,JCIM,JCIM,MMTAU)
      DIMENSION X(9),Y(9),C(9)
C*****************************************************************
      DIMENSION XA(1),XA2(1)
      DIMENSION XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C*****************************************************************
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM1B.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCV1B'
          CALL TIMIT(1)
          CALL FLUSH(IOUT)
        END IF
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'TAU'
C**IF IT IS IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF

      IF(IND.NE.0)GO TO 7777

C**FIND MAX QUANTUM
      MAXQ2=0
      DO IR=1,NSIZE
        NR1=IP2(IR,1)
        NR2=IP2(IR,2)
        IF(NR1.GT.MAXQ2)MAXQ2=NR1
        IF(NR2.GT.MAXQ2)MAXQ2=NR2
      END DO

C***********************************ALGORITHM FROM V0CV1
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=JNTFAC(NMODE,ICOUPL,1)
      IFACTC=JNTFAC(NMODE,ICOUPC,1)
C**INTRINSIC
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1

      KA=KROT/2
      LMAX=1+MOD(KA,2)
      FACTRC=0.D0
      IF(J21.GT.1)FACTRC=IFACTC
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      MD=MD1*MDT

C**FORM INDIVIDUAL INTEGRATION TERMS (START)
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          DO NRR=1,JCIM
            NR=NRR+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            X(2)=X(1)
            X(3)=HTAU(NR,MTAU,2,LMAX)*MD
            X(4)=X(1)
            X(5)=X(1)
            X(6)=X(3)
            X(7)=X(4)
            X(8)=HTAU(NR,MTAU,2,LMAX)*MD
            X(9)=X(1)
            DO NLL=1,JCIM
              NL=NLL+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              Y(2)=Y(1)
              Y(3)=Y(1)
              Y(4)=HTAU(NL,MTAU,2,LMAX)
              Y(5)=Y(1)
              Y(6)=Y(1)
              Y(7)=Y(4)
              Y(8)=Y(4)
              Y(9)=Y(1)
              DO K=1,I9
                T0(K,NLL,NRR,MTAU)=Y(10-K)*X(10-K)
              END DO
            END DO
          END DO
        ELSE
          DO NRR=1,JCIM
            NR=2*NRR-MOD(NRR,2)+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            X(2)=HTAU(NR,MTAU,2,LMAX)*MD
            X(3)=X(1)
            X(4)=X(2)
            X(5)=HTAU(NR,MTAU,3,LMAX)*MD
            X(6)=X(1)
            X(7)=HTAU(NR,MTAU,1,LMAX)*FACTRC*MD
            DO NLL=1,JCIM
              NL=2*NLL-MOD(NLL,2)+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO K=1,7
                T0(K,NLL,NRR,MTAU)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M=1,MM/MD1
      DO M=1,MM
        IF(.NOT.LINEAR)THEN
          DO NR1=1,JCI
            X(1)=H(NR1,M,2,1)
            X(2)=H(NR1,M,1,1)
            X(3)=X(2)
            X(4)=X(2)
            X(5)=X(1)
            X(6)=X(2)
            X(7)=X(1)
            X(8)=X(2)
            X(9)=X(2)
            DO NL1=1,JCI
              Y(1)=H(NL1,M,1,1)
              Y(2)=H(NL1,M,2,1)
              Y(3)=Y(1)
              Y(4)=Y(1)
              Y(5)=Y(2)
              Y(6)=Y(2)
              Y(7)=Y(1)
              Y(8)=Y(1)
              Y(9)=Y(1)
              DO K=1,I9
                X0(K,NL1,NR1,M)=Y(10-K)*X(10-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR1=1,JCI
            X(1)=H(NR1,M,2,1)
            X(2)=H(NR1,M,1,1)
            X(3)=H(NR1,M,3,1)
            X(4)=X(1)
            X(5)=X(2)
            X(6)=X(2)
            X(7)=X(2)
            DO NL1=1,JCI
              Y(1)=H(NL1,M,1,1)
              DO K=1,7
                X0(K,NL1,NR1,M)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
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
          IF(J21.GT.1.AND.ICOUPC.GE.1)READ(61)VR
          IF(ICOUPC.GE.1)READ(81)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPC.GE.1)READ(61)VRR
          IF(ICOUPC.GE.1)READ(81)VCR
        END IF
        IF(JCOUPL.GT.0)THEN
          READ(71)VP
        ELSE
          READ(71)VPR
        END IF

C***********************************************************

C       IF(.NOT.LINEAR)THEN
          DO IRHS1=1,JCI
            DO ILHS1=1,JCI
              DO K=1,I9
                TEMP1(K,ILHS1,IRHS1)=0.D0
              END DO
            END DO
          END DO
C**START 1-MODE INTEGRATION
C         DO M=1,MM/MD1
          DO M=1,MM
            DO I=1,9
              C(I)=0
            END DO
            IF(JCOUPL.GT.0)THEN
              IF(ICOUPC.GE.1)THEN
C**NO WATSON TERM IF RPH
                C(1)=VC(M,1)*IFACTC
                C(2)=VC(M,1)*IFACTC
                C(3)=VC(M,2)*IFACTC
                C(4)=VC(M,2)*IFACTC
                C(5)=VC(M,3)*IFACTC
                C(6)=VC(M,4)*IFACTC
                C(7)=VC(M,4)*IFACTC
                C(8)=VC(M,5)*IFACTC
                C(9)=VC(M,6)*IFACTC
                IF(J21.GT.1)C(9)=C(9)+VR(KROT,M)*IFACTC
              END IF
              C(9)=C(9)+VP(M)*IFACTL
            ELSE
              IF(ICOUPC.GE.1)THEN
C**NO WATSON TERM IF RPH
                C(1)=VCR(M,1)*IFACTC
                C(2)=VCR(M,1)*IFACTC
                C(3)=VCR(M,2)*IFACTC
                C(4)=VCR(M,2)*IFACTC
                C(5)=VCR(M,3)*IFACTC
                C(6)=VCR(M,4)*IFACTC
                C(7)=VCR(M,4)*IFACTC
                C(8)=VCR(M,5)*IFACTC
                C(9)=VCR(M,6)*IFACTC
                IF(J21.GT.1)C(9)=C(9)+VRR(KROT,M)*IFACTC
              END IF
              C(9)=C(9)+VPR(M)*IFACTL
            END IF
            DO IRHS1=1,JCI
              DO ILHS1=1,JCI
                DO K=1,I9
                  TEMP1(K,ILHS1,IRHS1)=TEMP1(K,ILHS1,IRHS1)+
     1            X0(K,ILHS1,IRHS1,M)*C(10-K)
                END DO
              END DO
            END DO
          END DO
C**END 1-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
          DO IRHS=1,NSIZE
            NR=IP2(IRHS,1)
            IRTAU=IP2(IRHS,2)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL=IP2(ILHS,1)
              ILTAU=IP2(ILHS,2)
              DO K=1,I9
                XA2(ILHS+J0)=XA2(ILHS+J0)+TEMP1(K,NL,NR)*
     1          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
              END DO
            END DO
          END DO
C       ELSE
C         DO IRHS1=1,JCI
C           DO ILHS1=1,JCI
C             DO K=1,7
C               TEMP1(K,ILHS1,IRHS1)=0.D0
C             END DO
C           END DO
C         END DO
C**START 1-MODE INTEGRATION
C         DO M=1,MM/MD1
C         DO M=1,MM
C           Q=XQ(M)
C           IF(JCOUPL.GT.0)THEN
C             DO IRHS1=1,JCI
C               DO ILHS1=1,JCI
C                 IF(ICOUPC.GE.1)THEN
C**(6) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
C                   DO K=1,5
C                     TEMP1(K,ILHS1,IRHS1)=TEMP1(K,ILHS1,IRHS1)-
C    1                X0(K,ILHS1,IRHS1,M)*VC(M,K)*IFACTC
C                   END DO
C                   TEMP1(4,ILHS1,IRHS1)=TEMP1(4,ILHS1,IRHS1)-
C    1              X0(4,ILHS1,IRHS1,M)*VC(M,4)*IFACTC
C                   TEMP1(7,ILHS1,IRHS1)=TEMP1(7,ILHS1,IRHS1)+
C    1              X0(7,ILHS1,IRHS1,M)*VR(KROT,M)
C                 END IF
C**INCLUDE K.E. ONCE FOR TAU
C                 TEMP1(5,ILHS1,IRHS1)=TEMP1(5,ILHS1,IRHS1)-
C    1            X0(5,ILHS1,IRHS1,M)/(8*Q*Q)
C***************************
C                 TEMP1(6,ILHS1,IRHS1)=TEMP1(6,ILHS1,IRHS1)+
C    1            X0(6,ILHS1,IRHS1,M)*VP(M)*IFACTL
C               END DO
C             END DO
C           ELSE
C             DO IRHS1=1,JCI
C               DO ILHS1=1,JCI
C                 IF(ICOUPC.GE.1)THEN
C**(6) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
C                   DO K=1,5
C                     TEMP1(K,ILHS1,IRHS1)=TEMP1(K,ILHS1,IRHS1)-
C    1                X0(K,ILHS1,IRHS1,M)*VCR(M,K)*IFACTC
C                   END DO
C                   TEMP1(4,ILHS1,IRHS1)=TEMP1(4,ILHS1,IRHS1)-
C    1              X0(4,ILHS1,IRHS1,M)*VCR(M,4)*IFACTC
C                   TEMP1(7,ILHS1,IRHS1)=TEMP1(7,ILHS1,IRHS1)+
C    2              X0(7,ILHS1,IRHS1,M)*VRR(KROT,M)
C                 END IF
C**INCLUDE K.E. ONCE FOR TAU
C                 TEMP1(5,ILHS1,IRHS1)=TEMP1(5,ILHS1,IRHS1)-
C    1            X0(5,ILHS1,IRHS1,M)/(8*Q*Q)
C***************************
C                 TEMP1(6,ILHS1,IRHS1)=TEMP1(6,ILHS1,IRHS1)+
C    1            X0(6,ILHS1,IRHS1,M)*VPR(M)*IFACTL
C               END DO
C             END DO
C           END IF

C**********************************************************
C**INCLUDE K.E. ONCE FOR Q
C           IF(IWHICH.EQ.0)THEN
C             VHARM=OMEGA*OMEGA*Q*Q/2
C           ELSE
C             VHARM=0.D0
C           END IF
C           DO IRHS1=1,JCI
C             DO ILHS1=1,JCI
C               TEMP1(3,ILHS1,IRHS1)=TEMP1(3,ILHS1,IRHS1)-
C    1          X0(3,ILHS1,IRHS1,M)/2
C               TEMP1(6,ILHS1,IRHS1)=TEMP1(6,ILHS1,IRHS1)+
C    1          X0(6,ILHS1,IRHS1,M)*VHARM
C             END DO
C           END DO
C**********************************************************

C         END DO
C**END 1-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
C         DO IRHS=1,NSIZE
C           NR=IP2(IRHS,1)
C           IRRTAU=IP2(IRHS,2)
C           IRTAU=IRRTAU/2+MOD(IRRTAU,2)
C           J0=IRHS*(IRHS-1)/2
C           DO ILHS=1,IRHS
C             NL=IP2(ILHS,1)
C             ILLTAU=IP2(ILHS,2)
C             ILTAU=ILLTAU/2+MOD(ILLTAU,2)
C             DO K=1,7
C               XA2(ILHS+J0)=XA2(ILHS+J0)+TEMP1(K,NL,NR)*
C    1          T0(K,ILTAU,IRTAU,MTAU)
C             END DO
C           END DO
C         END DO
C       END IF
      END DO
C**END TAU LOOP (2-MODE INTEGRATION)
      CALL MATOUT(XA2,XRA2,NSIZE,21)
      GO TO 7777
6666  CALL MATIN(XA2,XRA2,NSIZE,21)
7777  CONTINUE

      IF(IND.NE.0)THEN
        DO IR=1,NSIZE
          NR1=IP2(IR,1)
          NR2=IP2(IR,2)
          IRL(NR1,NR2)=IR
        END DO
      END IF

C***********************************ALGORITHM FROM VCV1

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFFR=IOFFR+NCSIZE(I)
      END DO
      J0R=IOFFR*(IOFFR+1)/2
C**TEMPORARY-LAN
      JSR=1
C**TEMPORARY-LAN
C     DO 9999 ISMR1=1,NVSYM
      DO 9999 ISMR1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IF(K.NE.MODE1.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE

        IF(NONC1.GT.NONC2)THEN
        NONZ1=0
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          NR1=IPC1(IROFF,MODE1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          KCCCL=KCSZL
CC        IF(NONC1.GT.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             KCCCL=KCSZL
C           ELSE
C             KCCCL=IRH1
C           END IF
CC        END IF
          DO ILH1=1,KCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IF(K.NE.MODE1.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            NL1=IPC1(ILOFF,MODE1)
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=NR1
            NONZC1(2,NONZ1)=NL1
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
        IF(NONZ1.EQ.0)GO TO 5556
        NONZC1(3,NONZ1+1)=0
        END IF

      END IF

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN TAU
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME INVOLVING 'MOD2' (TAU)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IF(K.NE.MODE2.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE

        IF(NONC1.LE.NONC2)THEN
        NONZ2=0
C**CONTRACTION SCHEME INVOLVING 'MOD2' (TAU)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          NR2=IPC2(JROFF,MODE2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          LCCCL=LCSZL
CC        IF(NONC1.LE.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             LCCCL=LCSZL
C           ELSE
C             LCCCL=IRH2
C           END IF
CC        END IF
          DO ILH2=1,LCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IF(K.NE.MODE2.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            NL2=IPC2(JLOFF,MODE2)
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=NR2
            NONZC2(2,NONZ2)=NL2
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
        IF(NONZ2.EQ.0)GO TO 5555
        NONZC2(3,NONZ2+1)=0
        END IF

      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    IF(ITIM1B.EQ.0)THEN
C       WRITE(IOUT,*)'VCCV1B'
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2,
C    1  'NONC1,NONC2 = ',NONC1,NONC2
C       CALL FLUSH(IOUT)
CC    END IF
      IF(NONC1.LE.NONC2)GO TO 3003
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************
C
C**CONTRACTION SCHEME INVOLVING 'MOD2' (TAU)
C
      DO IRH2=1,LCSZR
        JROFF=IRH2+LCOFFR
        NR2=IPC2(JROFF,MODE2)
        DO I=1,NKVALR
          DO J=1,NKVALL
            DO K=1,NLVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          LCCCL=LCSZL
        ELSE
          LCCCL=IRH2
        END IF
        DO ILH2=1,LCCCL
          JLOFF=ILH2+LCOFFL
          NL2=IPC2(JLOFF,MODE2)
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
          DO K=1,NKMOD2
            IF(K.NE.MODE2.AND.
     1      (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 6001
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH2.EQ.IRH2)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NKVALR
            DO J=1,NKVALL
              XCON2(J,I)=0
            END DO
          END DO
C
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
C
          IPREV=0
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            INEXT=NONZC1(3,INH1+1)
            ILH1=NONZC1(4,INH1)
            NR1=NONZC1(1,INH1)
            NL1=NONZC1(2,INH1)
            IF(IPREV.NE.IRH1)THEN
C**NEW RHS
              DO J=1,NKVALL
                TEMP(J,1)=0
              END DO
            END IF

C**TEMPORARY
      IR=IRL(NR1,NR2)
      IL=IRL(NL1,NL2)
C**TEMPORARY

C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XXX=XA2(I)

            DO J=1,NKVALL
              TEMP(J,1)=TEMP(J,1)+XXX*CFS1(ILH1,J,KEL,KSML)
            END DO

            IF(IRH1.NE.INEXT)THEN

            DO I=1,NKVALR
              DO J=1,NKVALL
                XCON2(J,I)=XCON2(J,I)+TEMP(J,1)*CFS1(IRH1,I,KEL,KSMR)
              END DO
            END DO

            END IF
            IPREV=IRH1
          END DO
CC
          DO I=1,NKVALR
            DO J=1,NKVALL
              DO K=1,NLVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6001      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILL,IKR,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKR,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILR,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILR,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO
C*************************************

      GO TO 5555

3003  CONTINUE

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************
C
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
C
C**TEMPORARY
C     WRITE(IOUT,*)'POINT 1'
C     WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C     WRITE(IOUT,*)'KCOFFR,KCOFFL = ',KCOFFR,KCOFFL
C     WRITE(IOUT,*)'MODE1,MODE2 = ',MODE1,MODE2
C     WRITE(IOUT,*)'NKVALR,NKVALL = ',NKVALR,NKVALL
C     WRITE(IOUT,*)'NLVALR,NLVALL = ',NLVALR,NLVALL
C     WRITE(IOUT,*)'ISMR1,ISML1 = ',ISMR1,ISML1
C     CALL FLUSH(IOUT)
C**TEMPORARY

      DO IRH1=1,KCSZR
        IROFF=IRH1+KCOFFR
        NR1=IPC1(IROFF,MODE1)
        DO I=1,NLVALR
          DO J=1,NLVALL
            DO K=1,NKVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          KCCCL=KCSZL
        ELSE
          KCCCL=IRH1
        END IF
C**TEMPORARY
C     WRITE(IOUT,*)'POINT 2, IRH1 = ',IRH1
C     CALL FLUSH(IOUT)
C**TEMPORARY
        DO ILH1=1,KCCCL
          ILOFF=ILH1+KCOFFL
          NL1=IPC1(ILOFF,MODE1)
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
          DO K=1,NKMOD1
            IF(K.NE.MODE1.AND.
     1      (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 6002
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH1.EQ.IRH1)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NLVALR
            DO J=1,NLVALL
              XCON2(J,I)=0
            END DO
          END DO
C
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
C
C**TEMPORARY
C     WRITE(IOUT,*)'POINT 2A, ILH1 = ',ILH1
C     CALL FLUSH(IOUT)
C**TEMPORARY
          JPREV=0
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            JNEXT=NONZC2(3,INH2+1)
            ILH2=NONZC2(4,INH2)
            NR2=NONZC2(1,INH2)
            NL2=NONZC2(2,INH2)
            IF(JPREV.NE.IRH2)THEN
C**NEW RHS
              DO J=1,NLVALL
                TEMP(J,1)=0
              END DO
            END IF

C**TEMPORARY
      IR=IRL(NR1,NR2)
      IL=IRL(NL1,NL2)
C**TEMPORARY

C**TEMPORARY
C     WRITE(IOUT,*)'POINT 2B, INH2,NR1,NR2,NL1,NL2,IL,IR = ',
C    1INH2,NR1,NR2,NL1,NL2,IL,IR
C     CALL FLUSH(IOUT)
C**TEMPORARY

C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XXX=XA2(I)

            DO J=1,NLVALL
              TEMP(J,1)=TEMP(J,1)+XXX*CFS2(ILH2,J,KEL,LSML)
            END DO
5002        CONTINUE

            IF(IRH2.NE.JNEXT)THEN

            DO I=1,NLVALR
              DO J=1,NLVALL
                XCON2(J,I)=XCON2(J,I)+TEMP(J,1)*CFS2(IRH2,I,KEL,LSMR)
              END DO
            END DO

            END IF
            JPREV=IRH2
          END DO
CC
C**TEMPORARY
C     WRITE(IOUT,*)'POINT 2C'
C     CALL FLUSH(IOUT)
C**TEMPORARY
          DO I=1,NLVALR
            DO J=1,NLVALL
              DO K=1,NKVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6002      CONTINUE
C**TEMPORARY
C     WRITE(IOUT,*)'POINT 2D'
C     CALL FLUSH(IOUT)
C**TEMPORARY
        END DO

C**TEMPORARY
C     WRITE(IOUT,*)'POINT 3'
C     CALL FLUSH(IOUT)
C**TEMPORARY

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**TEMPORARY
C     WRITE(IOUT,*)'POINT 4'
C     CALL FLUSH(IOUT)
C**TEMPORARY
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKR,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKL,ILR,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILR,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO
C*************************************

5555  CONTINUE
5556  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM1B.EQ.0)THEN
        IF(IND.NE.0)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM1B=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI2A(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,
     1H1,XQ1,H2,XQ2,NN1,MM1,NN2,MM2,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,
     2XA,
     2ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,
     3TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     4KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I5,XKAN,MAXQU,MAXPOW,NP2,CP2,
     5JP2,NTOT2,MAX2,INDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,6),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,6),VRR(J21,MM2,MM1)
      REAL*4 XRA2(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP2(ISIZE2,2),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I5,JCI2,JCI2)
      DIMENSION X0(I5,JCI1,JCI1,MM1),Y0(I5,JCI2,JCI2,MM2)
      DIMENSION X(5),Y(5),C(5)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA2(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),JP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI2A(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,
     1H1,XQ1,H2,XQ2,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,IP,ISIZMX,IPC,
     2ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,XK,TEMP,
     3XCON,
     3NVAL,ISTART,IEND,TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,
     4VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I5,
     5XKAN,MAXQU,MAXPOW,NP2,CP2,JP2,NTOT2,MAX2,INDK)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI2A(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,
     1H1,XQ1,H2,XQ2,NN1,MH1,MM1,NN2,MH2,MM2,IP,ISIZMX,IPC,ICSIZE,
     2IPSIZE,KCONT,
     2XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,
     3TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     4KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I5,XKAN,MAXQU,MAXPOW,NP2,CP2,
     5JP2,NTOT2,MAX2,INDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,6),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,6),VRR(J21,MM2,MM1)
      REAL*4 XRA2(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP2(ISIZE2,2),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I5,JCI2,JCI2)
      DIMENSION X0(I5,JCI1,JCI1,MM1),Y0(I5,JCI2,JCI2,MM2)
      DIMENSION X(5),Y(5),C(5)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA2(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),JP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM2A.EQ.1)THEN
        WRITE(IOUT,*)'Calculating VCCI2A'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' AND 'MOD2'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C**ANALYTIC
      IND=INDK(MOD1)+MOD2
C**ANALYTIC

C***********************************ALGORITHM FROM V0CI2
C**TEMPORARY
      IF(ICONDP.NE.0)THEN
        IF(MOLINC.GT.0)GO TO 99999
        GO TO 6666
      END IF
C**TEMPORARY

      IFACTC=INTFAC(NMODE,ICOUPC,2)
      IFACTL=INTFAC(NMODE,ICOUPL,2)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(72)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(72)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.1)READ(62)VR
        IF(ICOUPC.GT.1)READ(82)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.1)READ(62)VRR
        IF(ICOUPC.GT.1)READ(82)VCR
      END IF

      IF(MOLINC.GT.0)GO TO 99999

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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(3)
            Y(5)=Y(3)
            DO K=1,I5
              X0(K,NL1,NR1,M1)=Y(6-K)*X(6-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(3)
            Y(5)=Y(1)
            DO K=1,I5
              Y0(K,NL2,NR2,M2)=Y(6-K)*X(6-K)
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
            DO K=1,I5
              TEMP1(K,ILHS2,IRHS2)=0.D0
            END DO
          END DO
        END DO
C**START 1-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO I=1,5
            C(I)=0
          END DO
          IF(JCOUPL.GT.0)THEN
            IF(ICOUPC.GT.1)THEN
              C(1)=VC(M2,M1,3)*IFACTC
              C(2)=VC(M2,M1,4)*IFACTC
              C(3)=VC(M2,M1,4)*IFACTC
              C(4)=VC(M2,M1,5)*IFACTC
              C(5)=VC(M2,M1,6)*IFACTC
              IF(J21.GT.1)C(5)=C(5)+VR(KROT,M2,M1)*IFACTC
            END IF
            IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(5)=C(5)+VP(M2,M1)*IFACTL
          ELSE
            IF(ICOUPC.GT.1)THEN
              C(1)=VCR(M2,M1,3)*IFACTC
              C(2)=VCR(M2,M1,4)*IFACTC
              C(3)=VCR(M2,M1,4)*IFACTC
              C(4)=VCR(M2,M1,5)*IFACTC
              C(5)=VCR(M2,M1,6)*IFACTC
              IF(J21.GT.1)C(5)=C(5)+VRR(KROT,M2,M1)*IFACTC
            END IF
            IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(5)=C(5)+VPR(M2,M1)*IFACTL
          END IF
          DO IRHS2=1,JCI2
            DO ILHS2=1,JCI2
              DO K=1,I5
                TEMP1(K,ILHS2,IRHS2)=TEMP1(K,ILHS2,IRHS2)+
     1          Y0(K,ILHS2,IRHS2,M2)*C(6-K)
              END DO
            END DO
          END DO
        END DO
C**END 1-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
        DO IRHS=1,NSIZE
          NR1=IP2(IRHS,1)
          NR2=IP2(IRHS,2)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP2(ILHS,1)
            NL2=IP2(ILHS,2)
            DO K=1,I5
              XA2(ILHS+J0)=XA2(ILHS+J0)+
     1        TEMP1(K,NL2,NR2)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 2-MODE INTEGRATION
      CALL MATOUT(XA2,XRA2,NSIZE,22)
      GO TO 7777
6666  CALL MATIN(XA2,XRA2,NSIZE,22)
7777  CONTINUE

C     IF(MOLINC.GT.0)RETURN

C***********************************ALGORITHM FROM VCI2

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFF=IOFF+NCSIZE(I)
      END DO
      J0R=IOFF*(IOFF+1)/2
C**TEMPORARY-LAN
      JS=1
C**TEMPORARY-LAN
C     DO 9999 ISM1=1,NVSYM
      DO 9999 ISM1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
        NLVAL=NCVAL(2,ISM2)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
        NLVAL=NCVAL(1,ISM1)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR1=IPC(IROFF,MODE1)
        NR2=IPC(IROFF,MODE2)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.(IPC(IROFF,K).NE.
     1      IPC(ILOFF,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPC(ILOFF,MODE1)
          NL2=IPC(ILOFF,MODE2)
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
          XYZ=XA2(I)
          ZYX=0
          IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
            DO I=1,NP2(IND)
              K=JP2(I,IND,1)+1
              L=JP2(I,IND,2)+1
              ZYZ=ZYX+CP2(I,IND)*XKAN(NL1,NR1,K,MOD1)*
     1        XKAN(NL2,NR2,L,MOD2)
            END DO
C**ANALYTIC
          END IF
3000      CONTINUE
          XK(ILHS,IRHS)=(XYZ+ZYX)*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

C     IF(KCONT.EQ.1)THEN
C**SCHEME 1 ACTIVE (K POINTS TO 1; L POINTS TO 2)
C       KOFF=0
C       DO IR=1,NLVAL
C**DIAGONAL BLOCKS IR ALL HAVE NON-ZERO OVERLAP IN L
C         LROFF=0
C         DO IRHS=1,NKVAL
C           IROFF=1+LROFF+KOFF+IOFF
C**TEMPORARY-LAN
C           J0=IROFF*(IROFF-1)/2+IOFF-J0R
C**TEMPORARY-LAN
C           NR=IP(IROFF,KCONT)
C           LLOFF=0
C           DO ILHS=1,IRHS
C             ILOFF=1+LLOFF+KOFF
C             NL=IP(ILOFF+IOFF,KCONT)
C             XYZ=XCON(NL,NR)
C             XA(ILOFF+J0)=XA(ILOFF+J0)+XYZ
C             LLOFF=LLOFF+NLVAL
C           END DO
C           LROFF=LROFF+NLVAL
C         END DO
C         KOFF=KOFF+1
C       END DO
C     ELSE
C**SCHEME 2 ACTIVE (L POINTS TO 1; K POINTS TO 2)
C       LOFF=0
C       DO IR=1,NLVAL
C**DIAGONAL BLOCKS IR ALL HAVE NON-ZERO OVERLAP IN L
C         DO IRHS=1,NKVAL
C           IROFF=IRHS+LOFF+IOFF
C**TEMPORARY-LAN
C           J0=IROFF*(IROFF-1)/2+IOFF-J0R
C**TEMPORARY-LAN
C           NR=IP(IROFF,KCONT)
C           DO ILHS=1,IRHS
C             ILOFF=ILHS+LOFF
C             NL=IP(ILOFF+IOFF,KCONT)
C             XYZ=XCON(NL,NR)
C             XA(ILOFF+J0)=XA(ILOFF+J0)+XYZ
C           END DO
C         END DO
C         LOFF=LOFF+NKVAL
C       END DO
C     END IF

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
C**TEMPORARY-LAN
        J0=(IROFF)*(IROFF-1)/2+IOFF-J0R
C**TEMPORARY-LAN
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IP(IROFF,KCONT)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IP(IROFF,MCONT).NE.IP(ILOFF,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IP(ILOFF,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO


C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
99999 CONTINUE

      IF(ITIM2A.EQ.1)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM2A=2
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI2B(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,H1,XQ1,H2,
     1XQ2,NN1,MM1,NN2,MM2,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,
     2IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,IRL,
     3XK3,TEMP,KTEMP,XCON1,NVAL1,NVAL2,ISTART,IEND,TEMP1,
     4JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,
     5CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I5,XKAN,MAXQU,
     6MAXPOW,NP2,CP2,JP2,NTOT2,MAX2,INDK,MAXVAL,MAXQ2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,6),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,6),VRR(J21,MM2,MM1)
      REAL*4 XRA2(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP2(ISIZE2,2),TEMP(KTEMP,1)
      DIMENSION TEMP1(I5,JCI2,JCI2)
      DIMENSION X0(I5,JCI1,JCI1,MM1),Y0(I5,JCI2,JCI2,MM2)
      DIMENSION X(5),Y(5),C(5)
C***************************************************************
      DIMENSION XA(1),XA2(1)
      DIMENSION IRL(MAXQ2,MAXQ2),XCON1(MAXVAL,MAXVAL)
      DIMENSION XK3(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C***************************************************************
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),JP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI2B(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,H1,XQ1,H2,
     1XQ2,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,
     2IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,
     3NSIZE,IRL,XK3,TEMP,KTEMP,XCON1,NVAL1,NVAL2,ISTART,
     4IEND,TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     5KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I5,
     6XKAN,MAXQU,MAXPOW,NP2,CP2,JP2,NTOT2,MAX2,INDK,MAXVAL,MAXQ2)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI2B(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,H1,XQ1,H2,
     1XQ2,NN1,MH1,MM1,NN2,MH2,MM2,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,
     2ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,
     3IRL,XK3,TEMP,KTEMP,XCON1,NVAL1,NVAL2,ISTART,IEND,
     4TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,
     5NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I5,XKAN,
     6MAXQU,MAXPOW,NP2,CP2,JP2,NTOT2,MAX2,INDK,MAXVAL,MAXQ2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,6),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,6),VRR(J21,MM2,MM1)
      REAL*4 XRA2(1)

C**TEMPORARY
      DIMENSION IRL(MAXQ2,MAXQ2)
C**TEMPORARY

      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP2(ISIZE2,2),TEMP(KTEMP,1)
      DIMENSION TEMP1(I5,JCI2,JCI2)
      DIMENSION X0(I5,JCI1,JCI1,MM1),Y0(I5,JCI2,JCI2,MM2)
      DIMENSION X(5),Y(5),C(5)
C***************************************************************
      DIMENSION XA(1),XA2(1)
      DIMENSION XCON1(MAXVAL,MAXVAL)
      DIMENSION XK3(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C***************************************************************
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),JP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)
C**ANALYTIC
      COMMON/NONZ/NONZ1,NONZ2
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
C***************************************************************
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C***************************************************************
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM2A.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCI2B'
          CALL TIMIT(1)
          CALL FLUSH(IOUT)
        END IF
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L)
C**IF IT IS IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
      IF(IND.NE.0)GO TO 7777

C**FIND MAX QUANTUM
      MAXQ2=0
      DO IR=1,NSIZE
        NR1=IP2(IR,1)
        NR2=IP2(IR,2)
        IF(NR1.GT.MAXQ2)MAXQ2=NR1
        IF(NR2.GT.MAXQ2)MAXQ2=NR2
      END DO

C**ANALYTIC
      IND2=INDK(MOD1)+MOD2
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI2
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,2)
      IFACTL=INTFAC(NMODE,ICOUPL,2)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GE.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(72)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(72)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.1)READ(62)VR
        IF(ICOUPC.GT.1)READ(82)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.1)READ(62)VRR
        IF(ICOUPC.GT.1)READ(82)VCR
      END IF

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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(3)
            Y(5)=Y(3)
            DO K=1,I5
              X0(K,NL1,NR1,M1)=Y(6-K)*X(6-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(3)
            Y(5)=Y(1)
            DO K=1,I5
              Y0(K,NL2,NR2,M2)=Y(6-K)*X(6-K)
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
            DO K=1,I5
              TEMP1(K,ILHS2,IRHS2)=0
            END DO
          END DO
        END DO
C**START 1-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO I=1,5
            C(I)=0
          END DO
          IF(JCOUPL.GE.0)THEN
            IF(ICOUPC.GT.1)THEN
              C(1)=VC(M2,M1,3)*IFACTC
              C(2)=VC(M2,M1,4)*IFACTC
              C(3)=VC(M2,M1,4)*IFACTC
              C(4)=VC(M2,M1,5)*IFACTC
              C(5)=VC(M2,M1,6)*IFACTC
              IF(J21.GT.1)C(5)=C(5)+VR(KROT,M2,M1)*IFACTC
            END IF
            IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(5)=C(5)+VP(M2,M1)*IFACTL
          ELSE
            IF(ICOUPC.GT.1)THEN
              C(1)=VCR(M2,M1,3)*IFACTC
              C(2)=VCR(M2,M1,4)*IFACTC
              C(3)=VCR(M2,M1,4)*IFACTC
              C(4)=VCR(M2,M1,5)*IFACTC
              C(5)=VCR(M2,M1,6)*IFACTC
              IF(J21.GT.1)C(5)=C(5)+VRR(KROT,M2,M1)*IFACTC
            END IF
            IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(5)=C(5)+VPR(M2,M1)*IFACTL
          END IF
          DO IRHS2=1,JCI2
            DO ILHS2=1,JCI2
              DO K=1,I5
                TEMP1(K,ILHS2,IRHS2)=TEMP1(K,ILHS2,IRHS2)+
     1          Y0(K,ILHS2,IRHS2,M2)*C(6-K)
              END DO
            END DO
          END DO
        END DO
C**END 1-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
        DO IRHS=1,NSIZE
          NR1=IP2(IRHS,1)
          NR2=IP2(IRHS,2)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP2(ILHS,1)
            NL2=IP2(ILHS,2)
            DO K=1,I5
              XA2(ILHS+J0)=XA2(ILHS+J0)+
     1        TEMP1(K,NL2,NR2)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 2-MODE INTEGRATION
      CALL MATOUT(XA2,XRA2,NSIZE,22)
      GO TO 7777
6666  CALL MATIN(XA2,XRA2,NSIZE,22)
7777  CONTINUE

      IF(IND.NE.0)THEN
        DO IR=1,NSIZE
          NR1=IP2(IR,1)
          NR2=IP2(IR,2)
          IRL(NR1,NR2)=IR
        END DO
      END IF

C***********************************ALGORITHM FROM VCI2

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFFR=IOFFR+NCSIZE(I)
      END DO
      J0R=IOFFR*(IOFFR+1)/2
C**TEMPORARY-LAN
      JSR=1
C**TEMPORARY-LAN
C     DO 9999 ISMR1=1,NVSYM
      DO 9999 ISMR1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IF(K.NE.MODE1.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE

        IF(NONC1.GT.NONC2)THEN
        NONZ1=0
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          NR1=IPC1(IROFF,MODE1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          KCCCL=KCSZL
CC        IF(NONC1.GT.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             KCCCL=KCSZL
C           ELSE
C             KCCCL=IRH1
C           END IF
CC        END IF
          DO ILH1=1,KCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IF(K.NE.MODE1.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            NL1=IPC1(ILOFF,MODE1)
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=NR1
            NONZC1(2,NONZ1)=NL1
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
        IF(NONZ1.EQ.0)GO TO 5556
        NONZC1(3,NONZ1+1)=0
        END IF

      END IF

CCCC  IF(NONZ1.EQ.0)GO TO 5556
CCCC  NONZC1(3,NONZ1+1)=0

4444  CONTINUE

C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IF(K.NE.MODE2.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE

        IF(NONC1.LE.NONC2)THEN
        NONZ2=0
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          NR2=IPC2(JROFF,MODE2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          LCCCL=LCSZL
CC        IF(NONC1.LE.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             LCCCL=LCSZL
C           ELSE
C             LCCCL=IRH2
C           END IF
CC        END IF
          DO ILH2=1,LCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IF(K.NE.MODE2.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            NL2=IPC2(JLOFF,MODE2)
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=NR2
            NONZC2(2,NONZ2)=NL2
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
        IF(NONZ2.EQ.0)GO TO 5555
        NONZC2(3,NONZ2+1)=0
        END IF

      END IF

CCCC  IF(NONZ2.EQ.0)GO TO 5555
CCCC  NONZC2(3,NONZ2+1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM2A.EQ.0)THEN
C       WRITE(IOUT,*)'VCCI2B - ISMR1,ISML1 = ',ISMR1,ISML1
C       WRITE(IOUT,*)'ICONDP = ',ICONDP
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C       WRITE(IOUT,*)'NONC1,NONC2 = ',NONC1,NONC2
C       WRITE(IOUT,*)
C     END IF

      IF(NONC1.LE.NONC2)GO TO 3003
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************
C
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
C
      DO IRH2=1,LCSZR
        JROFF=IRH2+LCOFFR
        NR2=IPC2(JROFF,MODE2)
        DO I=1,NKVALR
          DO J=1,NKVALL
            DO K=1,NLVALL
              XK3(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK3(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          LCCCL=LCSZL
        ELSE
          LCCCL=IRH2
        END IF
        DO ILH2=1,LCCCL
          JLOFF=ILH2+LCOFFL
          NL2=IPC2(JLOFF,MODE2)
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
          DO K=1,NKMOD2
            IF(K.NE.MODE2.AND.
     1      (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 6001
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH2.EQ.IRH2)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NKVALR
            DO J=1,NKVALL
              XCON1(J,I)=0
            END DO
          END DO
C
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
C
C         DO IRH1=1,KCSZR
C           IROFF=IRH1+KCOFFR
C           NR1=IPC1(IROFF,MODE1)
C           DO J=1,NKVALL
C             TEMP(J,1)=0
C           END DO
C           DO ILH1=1,KCSZL
C             ILOFF=ILH1+KCOFFL
C             NL1=IPC1(ILOFF,MODE1)
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
CC            DO K=1,NKMOD1
CC              IF(K.NE.MODE1.AND.
CC   1          (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 5001
CC            END DO
C             DO K=1,MODE1-1
C               IF(IPC1(IROFF,K).NE.IPC1(ILOFF,K))GO TO 5001
C             END DO
C             DO K=MODE1+1,NKMOD1
C               IF(IPC1(IROFF,K).NE.IPC1(ILOFF,K))GO TO 5001
C             END DO

          IPREV=0
          DO INH1=1,NONZ1
              IRH1=NONZC1(3,INH1)
              INEXT=NONZC1(3,INH1+1)
              ILH1=NONZC1(4,INH1)
              NR1=NONZC1(1,INH1)
              NL1=NONZC1(2,INH1)
              IF(IPREV.NE.IRH1)THEN
C**NEW RHS
                DO J=1,NKVALL
                  TEMP(J,1)=0
                END DO
              END IF

C**FIND RHS INDEX
C             DO IR=1,NSIZE
C               IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1000
C             END DO
1000          CONTINUE
C**FIND LHS INDEX
C             DO IL=1,NSIZE
C               IF(NL1.EQ.IP2(IL,1).AND.NL2.EQ.IP2(IL,2))GO TO 2000
C             END DO
2000          CONTINUE

C**TEMPORARY
      IR=IRL(NR1,NR2)
      IL=IRL(NL1,NL2)
C**TEMPORARY

C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
C             XYZ=XA2(I)
C             ZYX=0
C             IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
C               DO I=1,NP2(IND2)
C                 K=JP2(I,IND2,1)+1
C                 L=JP2(I,IND2,2)+1
C                 ZYX=ZYX+CP2(I,IND2)*XKAN(NL1,NR1,K,MOD1)*
C    1            XKAN(NL2,NR2,L,MOD2)
C               END DO
C**ANALYTIC
C             END IF
C             XXX=XYZ+ZYX
              XXX=XA2(I)

              DO J=1,NKVALL
                TEMP(J,1)=TEMP(J,1)+XXX*CFS1(ILH1,J,KEL,KSML)
              END DO
5001          CONTINUE

C           END DO
            IF(IRH1.NE.INEXT)THEN

            DO I=1,NKVALR
              DO J=1,NKVALL
                XCON1(J,I)=XCON1(J,I)+TEMP(J,1)*CFS1(IRH1,I,KEL,KSMR)
              END DO
            END DO

            END IF
            IPREV=IRH1
          END DO
CC
          DO I=1,NKVALR
            DO J=1,NKVALL
              DO K=1,NLVALL
                XK3(K,J,I,1)=XK3(K,J,I,1)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON1(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK3(K,J,I,2)=XK3(K,J,I,2)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON1(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6001      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK3(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK3(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK3(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK3(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK3(ILL,IKR,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILL,KEL,LSMR)*XK3(ILR,IKR,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK3(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK3(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK3(ILR,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILR,KEL,LSMR)*XK3(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END
 
        END IF
      END DO
C*************************************

      GO TO 5555 

3003  CONTINUE

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************
C
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
C
      DO IRH1=1,KCSZR
        IROFF=IRH1+KCOFFR
        NR1=IPC1(IROFF,MODE1)
        DO I=1,NLVALR
          DO J=1,NLVALL
            DO K=1,NKVALL
              XK3(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK3(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          KCCCL=KCSZL
        ELSE
          KCCCL=IRH1
        END IF
        DO ILH1=1,KCCCL
          ILOFF=ILH1+KCOFFL
          NL1=IPC1(ILOFF,MODE1)
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
          DO K=1,NKMOD1
            IF(K.NE.MODE1.AND.
     1      (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 6002
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH1.EQ.IRH1)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NLVALR
            DO J=1,NLVALL
              XCON1(J,I)=0
            END DO
          END DO
C
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
C
C         DO IRH2=1,LCSZR
C           JROFF=IRH2+LCOFFR
C           NR2=IPC2(JROFF,MODE2)
C           DO J=1,NLVALL
C             TEMP(J,1)=0
C           END DO
C           DO ILH2=1,LCSZL
C             JLOFF=ILH2+LCOFFL
C             NL2=IPC2(JLOFF,MODE2)
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
CC            DO K=1,NKMOD2
CC              IF(K.NE.MODE2.AND.
CC   1          (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 5002
CC            END DO
C             DO K=1,MODE2-1
C               IF(IPC2(JROFF,K).NE.IPC2(JLOFF,K))GO TO 5002
C             END DO
C             DO K=MODE2+1,NKMOD2
C               IF(IPC2(JROFF,K).NE.IPC2(JLOFF,K))GO TO 5002
C             END DO

          JPREV=0
          DO INH2=1,NONZ2
              IRH2=NONZC2(3,INH2)
              JNEXT=NONZC2(3,INH2+1)
              ILH2=NONZC2(4,INH2)
              NR2=NONZC2(1,INH2)
              NL2=NONZC2(2,INH2)
              IF(JPREV.NE.IRH2)THEN
C**NEW RHS
                DO J=1,NLVALL
                  TEMP(J,1)=0
                END DO
              END IF

C**FIND RHS INDEX
C             DO IR=1,NSIZE
C               IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1001
C             END DO
1001          CONTINUE
C**FIND LHS INDEX
C             DO IL=1,NSIZE
C               IF(NL1.EQ.IP2(IL,1).AND.NL2.EQ.IP2(IL,2))GO TO 2001
C             END DO
2001          CONTINUE

C**TEMPORARY
      IR=IRL(NR1,NR2)
      IL=IRL(NL1,NL2)
C**TEMPORARY

C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
C             XYZ=XA2(I)
C             ZYX=0
C             IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
C               DO I=1,NP2(IND2)
C                 K=JP2(I,IND2,1)+1
C                 L=JP2(I,IND2,2)+1
C                 ZYX=ZYX+CP2(I,IND2)*XKAN(NL1,NR1,K,MOD1)*
C    1            XKAN(NL2,NR2,L,MOD2)
C               END DO
C**ANALYTIC
C             END IF
C             XXX=XYZ+ZYX
              XXX=XA2(I)

              DO J=1,NLVALL
                TEMP(J,1)=TEMP(J,1)+XXX*CFS2(ILH2,J,KEL,LSML)
              END DO
5002          CONTINUE

C           END DO
            IF(IRH2.NE.JNEXT)THEN

            DO I=1,NLVALR
              DO J=1,NLVALL
                XCON1(J,I)=XCON1(J,I)+TEMP(J,1)*CFS2(IRH2,I,KEL,LSMR)
              END DO
            END DO

            END IF
            JPREV=IRH2
          END DO
CC
          DO I=1,NLVALR
            DO J=1,NLVALL
              DO K=1,NKVALL
                XK3(K,J,I,1)=XK3(K,J,I,1)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON1(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK3(K,J,I,2)=XK3(K,J,I,2)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON1(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6002      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CFS1(IRH1,IKR,KEL,KSMR)*XK3(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CFS1(IRH1,IKR,KEL,KSMR)*XK3(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CFS1(IRH1,IKR,KEL,KSMR)*XK3(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK3(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CFS1(IRH1,IKR,KEL,KSMR)*XK3(IKR,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKR,KEL,KSMR)*XK3(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CFS1(IRH1,IKR,KEL,KSMR)*XK3(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK3(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CFS1(IRH1,IKR,KEL,KSMR)*XK3(IKL,ILR,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKL,KEL,KSMR)*XK3(IKR,ILR,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END
 
        END IF
      END DO
C*************************************

5555  CONTINUE
5556  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM2A.EQ.0)THEN
        IF(IND.NE.0)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM2A=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI2C(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,H1,XQ1,H2,
     1XQ2,NN1,MM1,NN2,MM2,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,
     2IP2,ISIZE2,XA2,XRA2,NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,IEND,
     3TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,
     4NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I5,XKAN,MAXQU,MAXPOW,NP2,CP2,
     5JP2,NTOT2,MAX2,INDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,6),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,6),VRR(J21,MM2,MM1)
      REAL*4 XRA2(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP2(ISIZE2,2),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I5,JCI2,JCI2),IRL(ICSIZE)
      DIMENSION X0(I5,JCI1,JCI1,MM1),Y0(I5,JCI2,JCI2,MM2)
      DIMENSION X(5),Y(5),C(5)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA2(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),JP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI2C(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,H1,XQ1,H2,XQ2,NN1,
     1MM1,MM1/MD1,NN2,MM2,MM2/MD2,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,
     2ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,
     3IEND,TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     4KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I5,XKAN,MAXQU,MAXPOW,NP2,
     5CP2,JP2,NTOT2,MAX2,INDK)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI2C(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,H1,XQ1,H2,
     1XQ2,NN1,MH1,MM1,NN2,MH2,MM2,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,
     2ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,
     3IEND,TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     4KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I5,XKAN,MAXQU,MAXPOW,NP2,CP2,
     5JP2,NTOT2,MAX2,INDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,6),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,6),VRR(J21,MM2,MM1)
      REAL*4 XRA2(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP2(ISIZE2,2),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I5,JCI2,JCI2),IRL(ICSIZE)
      DIMENSION X0(I5,JCI1,JCI1,MM1),Y0(I5,JCI2,JCI2,MM2)
      DIMENSION X(5),Y(5),C(5)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA2(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),JP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM2A.EQ.1)THEN
        WRITE(IOUT,*)'Calculating VCCI2C'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' AND 'MOD2'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      NLMODE=ICONT(2)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        NLMODE=ICONT(1)
        MCONT=1
      END IF

C**ANALYTIC
      IND=INDK(MOD1)+MOD2
C**ANALYTIC

C***********************************ALGORITHM FROM V0CI2
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,2)
      IFACTL=INTFAC(NMODE,ICOUPL,2)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(72)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(72)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.1)READ(62)VR
        IF(ICOUPC.GT.1)READ(82)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.1)READ(62)VRR
        IF(ICOUPC.GT.1)READ(82)VCR
      END IF

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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(3)
            Y(5)=Y(3)
            DO K=1,I5
              X0(K,NL1,NR1,M1)=Y(6-K)*X(6-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(3)
            Y(5)=Y(1)
            DO K=1,I5
              Y0(K,NL2,NR2,M2)=Y(6-K)*X(6-K)
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
            DO K=1,I5
              TEMP1(K,ILHS2,IRHS2)=0.D0
            END DO
          END DO
        END DO
C**START 1-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO I=1,5
            C(I)=0
          END DO
          IF(JCOUPL.GT.0)THEN
            IF(ICOUPC.GT.1)THEN
              C(1)=VC(M2,M1,3)*IFACTC
              C(2)=VC(M2,M1,4)*IFACTC
              C(3)=VC(M2,M1,4)*IFACTC
              C(4)=VC(M2,M1,5)*IFACTC
              C(5)=VC(M2,M1,6)*IFACTC
              IF(J21.GT.1)C(5)=C(5)+VR(KROT,M2,M1)*IFACTC
            END IF
            IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(5)=C(5)+VP(M2,M1)*IFACTL
          ELSE
            IF(ICOUPC.GT.1)THEN
              C(1)=VCR(M2,M1,3)*IFACTC
              C(2)=VCR(M2,M1,4)*IFACTC
              C(3)=VCR(M2,M1,4)*IFACTC
              C(4)=VCR(M2,M1,5)*IFACTC
              C(5)=VCR(M2,M1,6)*IFACTC
              IF(J21.GT.1)C(5)=C(5)+VRR(KROT,M2,M1)*IFACTC
            END IF
            IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(5)=C(5)+VPR(M2,M1)*IFACTL
          END IF
          DO IRHS2=1,JCI2
            DO ILHS2=1,JCI2
              DO K=1,I5
                TEMP1(K,ILHS2,IRHS2)=TEMP1(K,ILHS2,IRHS2)+
     1          Y0(K,ILHS2,IRHS2,M2)*C(6-K)
              END DO
            END DO
          END DO
        END DO
C**END 1-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
        DO IRHS=1,NSIZE
          NR1=IP2(IRHS,1)
          NR2=IP2(IRHS,2)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP2(ILHS,1)
            NL2=IP2(ILHS,2)
            DO K=1,I5
              XA2(ILHS+J0)=XA2(ILHS+J0)+
     1        TEMP1(K,NL2,NR2)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 2-MODE INTEGRATION
      CALL MATOUT(XA2,XRA2,NSIZE,22)
      GO TO 7777
6666  CALL MATIN(XA2,XRA2,NSIZE,22)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI2

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        KCOFF=ICOFF1
        LCOFF=ICOFF2
        KCSZ=ICSZ1
        LCSZ=ICSZ2
        ICNTZ=ICCONT(1)
        NKVAL=NCVAL(1,ISM1)
        NLVAL=NCVAL(2,ISM2)
      ELSE
        ISM=ISM2
        KCOFF=ICOFF2
        LCOFF=ICOFF1
        KCSZ=ICSZ2
        LCSZ=ICSZ1
        ICNTZ=ICCONT(2)
        NKVAL=NCVAL(2,ISM2)
        NLVAL=NCVAL(1,ISM1)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

C
      DO IKRL=1,KCSZ
        IOFFRL=IKRL+KCOFF
        NR1=IPC(IOFFRL,MODE1)
        NR2=IPC(IOFFRL,MODE2)
        DO IR=1,NSIZE
          IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1000
        END DO
1000    IRL(IKRL)=IR
      END DO
C
C**CASE 1..INTEGRATION MODES IN CONTRACTED SCHEME.
C**      ..OVERLAP OF NON-CONTRACTED BASIS.
      IF(ICNTZ.GT.0)THEN
        DO IRHS=1,KCSZ
          IOFFR=IRHS+KCOFF
          IR=IRL(IRHS)
          DO ILHS=1,IRHS
            IOFFL=ILHS+KCOFF
C**OVERLAP OF REMAINING STATES
            IS=1
            DO K=1,NKMODE
              IF(IS.EQ.0)GO TO 2000
              IF(K.NE.MODE1.AND.K.NE.MODE2.AND.(IPC(IOFFR,K).NE.
     1        IPC(IOFFL,K)))IS=0
            END DO
C**OVERLAP OF REMAINING STATES
            IL=IRL(ILHS)
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA2(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              DO I=1,NP2(IND)
                K=JP2(I,IND,1)+1
                L=JP2(I,IND,2)+1
                ZYZ=ZYX+CP2(I,IND)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)
              END DO
C**ANALYTIC
            END IF
2000        CONTINUE
            XK(ILHS,IRHS)=(XYZ+ZYX)*IS
            XK(IRHS,ILHS)=XK(ILHS,IRHS)
          END DO
        END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',KCSZ,NKVAL,KCSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
        CALL DGEMM('T','N',NKVAL,NKVAL,KCSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

        IRHS=0
        DO IKR=1,NKVAL
          DO ILR=1,LCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            DO IKL=1,IKR
              ILHS=(IKL-1)*LCSZ+ILR
              XYZ=XCON(IKL,IKR)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ
            END DO
          END DO
        END DO

C**CASE 2..INTEGRATION MODES IN NON-CONTRACTED SCHEME.
C**      ..OVERLAP OF CONTRACTED BASIS.
      ELSE
        IRHS=0
        DO ILR=1,NLVAL
          DO IKR=1,KCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            IOFFR=IKR+KCOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            IR=IRL(IKR)
            DO IKL=1,IKR
              ILHS=(ILR-1)*KCSZ+IKL
              IOFFL=IKL+KCOFF
C**OVERLAP OF REMAINING STATES
              IS=1
              DO K=1,NKMODE
                IF(IS.EQ.0)GO TO 3000
                IF(K.NE.MODE1.AND.K.NE.MODE2.AND.(IPC(IOFFR,K).NE.
     1          IPC(IOFFL,K)))IS=0
              END DO
C**OVERLAP OF REMAINING STATES
              IL=IRL(IKL)
C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
              XYZ=XA2(I)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ*IS
3000          CONTINUE
            END DO
          END DO
        END DO
      END IF

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM2A.EQ.1)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM2A=2
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI2D(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,H1,XQ1,H2,
     1XQ2,NN1,MM1,NN2,MM2,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,
     2IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,XK1,
     3XK3,TEMP,KTEMP,XCON1,NVAL1,NVAL2,ISTART,IEND,TEMP1,
     4JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,
     5CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I5,XKAN,MAXQU,
     6MAXPOW,NP2,CP2,JP2,NTOT2,MAX2,INDK,MAXVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,6),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,6),VRR(J21,MM2,MM1)
      REAL*4 XRA2(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP2(ISIZE2,2),TEMP(KTEMP,1)
      DIMENSION TEMP1(I5,JCI2,JCI2)
      DIMENSION X0(I5,JCI1,JCI1,MM1),Y0(I5,JCI2,JCI2,MM2)
      DIMENSION X(5),Y(5),C(5)
C***************************************************************
      DIMENSION XA(1),XA2(1)
      DIMENSION XK1(KTEMP,KTEMP),XCON1(MAXVAL,MAXVAL)
      DIMENSION XK3(MAXVAL,1)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C***************************************************************
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),JP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI2D(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,H1,XQ1,H2,
     1XQ2,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,
     2IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,
     3NSIZE,XK1,XK3,TEMP,KTEMP,XCON1,NVAL1,NVAL2,ISTART,
     4IEND,TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     5KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I5,
     6XKAN,MAXQU,MAXPOW,NP2,CP2,JP2,NTOT2,MAX2,INDK,MAXVAL)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI2D(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,H1,XQ1,H2,
     1XQ2,NN1,MH1,MM1,NN2,MH2,MM2,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,
     2ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP2,ISIZE2,XA2,XRA2,NSIZE,
     3XK1,XK3,TEMP,KTEMP,XCON1,NVAL1,NVAL2,ISTART,IEND,
     4TEMP1,JCI1,JCI2,X0,Y0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,
     5NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I5,XKAN,
     6MAXQU,MAXPOW,NP2,CP2,JP2,NTOT2,MAX2,INDK,MAXVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,6),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,6),VRR(J21,MM2,MM1)
      REAL*4 XRA2(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP2(ISIZE2,2),TEMP(KTEMP,1)
      DIMENSION TEMP1(I5,JCI2,JCI2)
      DIMENSION X0(I5,JCI1,JCI1,MM1),Y0(I5,JCI2,JCI2,MM2)
      DIMENSION X(5),Y(5),C(5)
C***************************************************************
      DIMENSION XA(1),XA2(1)
      DIMENSION XK1(KTEMP,KTEMP),XCON1(MAXVAL,MAXVAL)
      DIMENSION XK3(MAXVAL,1)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C***************************************************************
C**ANALYTIC
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
      DIMENSION NP2(NTOT2),CP2(MAX2,NTOT2),JP2(MAX2,NTOT2,2)
      DIMENSION INDK(1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
C***************************************************************
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C***************************************************************
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM2A.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCI2D'
          CALL TIMIT(1)
          CALL FLUSH(IOUT)
        END IF
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L)
C**IF IT IS IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
      IF(IND.NE.0)GO TO 7777

C**ANALYTIC
      IND2=INDK(MOD1)+MOD2
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI2
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,2)
      IFACTL=INTFAC(NMODE,ICOUPL,2)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(72)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(72)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.1)READ(62)VR
        IF(ICOUPC.GT.1)READ(82)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.1)READ(62)VRR
        IF(ICOUPC.GT.1)READ(82)VCR
      END IF

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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(3)
            Y(5)=Y(3)
            DO K=1,I5
              X0(K,NL1,NR1,M1)=Y(6-K)*X(6-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(3)
            Y(5)=Y(1)
            DO K=1,I5
              Y0(K,NL2,NR2,M2)=Y(6-K)*X(6-K)
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
            DO K=1,I5
              TEMP1(K,ILHS2,IRHS2)=0
            END DO
          END DO
        END DO
C**START 1-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO I=1,5
            C(I)=0
          END DO
          IF(JCOUPL.GT.0)THEN
            IF(ICOUPC.GT.1)THEN
              C(1)=VC(M2,M1,3)*IFACTC
              C(2)=VC(M2,M1,4)*IFACTC
              C(3)=VC(M2,M1,4)*IFACTC
              C(4)=VC(M2,M1,5)*IFACTC
              C(5)=VC(M2,M1,6)*IFACTC
              IF(J21.GT.1)C(5)=C(5)+VR(KROT,M2,M1)*IFACTC
            END IF
            IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(5)=C(5)+VP(M2,M1)*IFACTL
          ELSE
            IF(ICOUPC.GT.1)THEN
              C(1)=VCR(M2,M1,3)*IFACTC
              C(2)=VCR(M2,M1,4)*IFACTC
              C(3)=VCR(M2,M1,4)*IFACTC
              C(4)=VCR(M2,M1,5)*IFACTC
              C(5)=VCR(M2,M1,6)*IFACTC
              IF(J21.GT.1)C(5)=C(5)+VRR(KROT,M2,M1)*IFACTC
            END IF
            IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(5)=C(5)+VPR(M2,M1)*IFACTL
          END IF
          DO IRHS2=1,JCI2
            DO ILHS2=1,JCI2
              DO K=1,I5
                TEMP1(K,ILHS2,IRHS2)=TEMP1(K,ILHS2,IRHS2)+
     1          Y0(K,ILHS2,IRHS2,M2)*C(6-K)
              END DO
            END DO
          END DO
        END DO
C**END 1-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
        DO IRHS=1,NSIZE
          NR1=IP2(IRHS,1)
          NR2=IP2(IRHS,2)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP2(ILHS,1)
            NL2=IP2(ILHS,2)
            DO K=1,I5
              XA2(ILHS+J0)=XA2(ILHS+J0)+
     1        TEMP1(K,NL2,NR2)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 2-MODE INTEGRATION
      CALL MATOUT(XA2,XRA2,NSIZE,22)
      GO TO 7777
6666  CALL MATIN(XA2,XRA2,NSIZE,22)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI2

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
      JSR=1
      DO 9999 ISMR1=1,NVSYM
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IF(K.NE.MODE1.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE
        NONZ1=0
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          NR1=IPC1(IROFF,MODE1)
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IF(K.NE.MODE1.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR 'MOD1'
            NL1=IPC1(ILOFF,MODE1)
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=NR1
            NONZC1(2,NONZ1)=NL1
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
      END IF
      IF(NONZ1.EQ.0)GO TO 5556
      NONZC1(3,NONZ1+1)=0

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IF(K.NE.MODE2.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE
        NONZ2=0
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          NR2=IPC2(JROFF,MODE2)
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IF(K.NE.MODE2.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR 'MOD2'
            NL2=IPC2(JLOFF,MODE2)
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=NR2
            NONZC2(2,NONZ2)=NL2
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
      END IF
      IF(NONZ2.EQ.0)GO TO 5555
      NONZC2(3,NONZ2+1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM2A.EQ.0)THEN
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
        WRITE(IOUT,*)'NONC1,NONC2 = ',NONC1,NONC2
C     END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CASE 1..K IS IN CONTRACTED SCHEME, L IS IN UNCONTRACTED SCHEME
      IF(ICCONT(1).GT.0)THEN

        IF(NONC1.LE.NONC2)GO TO 3003
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        IPREV=0
        DO INH1=1,NONZ1
          IRH1=NONZC1(3,INH1)
C         IF(INH1.NE.NONZ1)THEN
            INEXT=NONZC1(3,INH1+1)
C         ELSE
C           INEXT=0
C         END IF
          ILH1=NONZC1(4,INH1)
          IF(IPREV.NE.IRH1)THEN
            DO INH2=1,NONZ2
              DO IKL=1,NKVALL
                XK3(IKL,INH2)=0
              END DO
            END DO
          END IF
          NR1=NONZC1(1,INH1)
          NL1=NONZC1(2,INH1)

C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            ILH2=NONZC2(4,INH2)
            NR2=NONZC2(1,INH2)
            NL2=NONZC2(2,INH2)
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1003
            END DO
1003        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IF(NL1.EQ.IP2(IL,1).AND.NL2.EQ.IP2(IL,2))GO TO 2003
            END DO
2003        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA2(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              DO I=1,NP2(IND2)
                K=JP2(I,IND2,1)+1
                L=JP2(I,IND2,2)+1
                ZYX=ZYX+CP2(I,IND2)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)
              END DO
C**ANALYTIC
            END IF
            DO IKL=1,NKVALL
              XK3(IKL,INH2)=XK3(IKL,INH2)+(XYZ+ZYX)*
     1        CFS1(ILH1,IKL,KEL,KSML)
            END DO
          END DO
          IF(IRH1.NE.INEXT)THEN
            DO INH2=1,NONZ2
              IRH2=NONZC2(3,INH2)
              ILH2=NONZC2(4,INH2)
              IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
                DO IKR=1,NKVALR
                  CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                  IRHS=(IKR-1)*LCSZR+IRH2
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO IKL=1,NKVALL
                    ILHS=(IKL-1)*LCSZL+ILH2
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK3(IKL,INH2)
                  END DO
                END DO
C**OFF-DIAGONAL BLOCK - END
              ELSE
C**DIAGONAL BLOCK - START
                DO IKR=1,NKVALR
                  CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                  IRHS=(IKR-1)*LCSZR+IRH2
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO IKL=1,IKR
                    IF(IKL.EQ.IKR.AND.ILH2.GT.IRH2)GO TO 1004
                    ILHS=(IKL-1)*LCSZL+ILH2
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK3(IKL,INH2)
1004                CONTINUE
                  END DO
                END DO
C**DIAGONAL BLOCK - END
              END IF
            END DO
          END IF
          IPREV=IRH1
        END DO
        GO TO 5555

3003    CONTINUE
        DO IRH1=1,KCSZR
          DO ILH1=1,KCSZL
            XK1(ILH1,IRH1)=0
          END DO
        END DO

C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
        DO INH2=1,NONZ2
          IRH2=NONZC2(3,INH2)
          ILH2=NONZC2(4,INH2)
          NR2=NONZC2(1,INH2)
          NL2=NONZC2(2,INH2)

C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            ILH1=NONZC1(4,INH1)
            NR1=NONZC1(1,INH1)
            NL1=NONZC1(2,INH1)
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1000
            END DO
1000        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IF(NL1.EQ.IP2(IL,1).AND.NL2.EQ.IP2(IL,2))GO TO 2000
            END DO
2000        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA2(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              DO I=1,NP2(IND2)
                K=JP2(I,IND2,1)+1
                L=JP2(I,IND2,2)+1
                ZYX=ZYX+CP2(I,IND2)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)
              END DO
C**ANALYTIC
            END IF
            XK1(ILH1,IRH1)=XYZ+ZYX
          END DO

C************************************DGEMM (RHS)
          CALL DGEMM('N','N',KCSZL,NKVALR,KCSZR,1.0D0,XK1(1,1),
     &    KTEMP,CFS1(1,1,KEL,KSMR),ISIZXX,0.0D0,TEMP(1,1),KTEMP)

C************************************DGEMM (LHS)
          CALL DGEMM('T','N',NKVALL,NKVALR,KCSZL,1.D0,
     &    CFS1(1,1,KEL,KSML),ISIZXX,TEMP,KTEMP,0.0D0,XCON1(1,1),MAXVAL)

          IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
            DO IKR=1,NKVALR
              IRHS=(IKR-1)*LCSZR+IRH2
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO IKL=1,NKVALL
                ILHS=(IKL-1)*LCSZL+ILH2
                XA(ILHS+J0)=XA(ILHS+J0)+XCON1(IKL,IKR)
              END DO
            END DO
C**OFF-DIAGONAL BLOCK - END
          ELSE
C**DIAGONAL BLOCK - START
            DO IKR=1,NKVALR
              IRHS=(IKR-1)*LCSZR+IRH2
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO IKL=1,IKR
                IF(IKL.EQ.IKR.AND.ILH2.GT.IRH2)GO TO 1002
                ILHS=(IKL-1)*LCSZL+ILH2
                XA(ILHS+J0)=XA(ILHS+J0)+XCON1(IKL,IKR)
1002            CONTINUE
              END DO
            END DO
C**DIAGONAL BLOCK - END
          END IF
        END DO

C**CASE 2..L IS IN CONTRACTED SCHEME, K IS IN UNCONTRACTED SCHEME
      ELSE

        IF(NONC1.GT.NONC2)GO TO 4004
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
        IPREV=0
        DO INH2=1,NONZ2
          IRH2=NONZC2(3,INH2)
C         IF(INH2.NE.NONZ2)THEN
            INEXT=NONZC2(3,INH2+1)
C         ELSE
C           INEXT=0
C         END IF
          ILH2=NONZC2(4,INH2)
          IF(IPREV.NE.IRH2)THEN
            DO INH1=1,NONZ1
              DO ILL=1,NLVALL
                XK3(ILL,INH1)=0
              END DO
            END DO
          END IF
          NR2=NONZC2(1,INH2)
          NL2=NONZC2(2,INH2)

C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            ILH1=NONZC1(4,INH1)
            NR1=NONZC1(1,INH1)
            NL1=NONZC1(2,INH1)
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1005
            END DO
1005        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IF(NL1.EQ.IP2(IL,1).AND.NL2.EQ.IP2(IL,2))GO TO 2005
            END DO
2005        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA2(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              DO I=1,NP2(IND2)
                K=JP2(I,IND2,1)+1
                L=JP2(I,IND2,2)+1
                ZYX=ZYX+CP2(I,IND2)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)
              END DO
C**ANALYTIC
            END IF
            DO ILL=1,NLVALL
              XK3(ILL,INH1)=XK3(ILL,INH1)+(XYZ+ZYX)*
     1        CFS2(ILH2,ILL,KEL,LSML)
            END DO
          END DO
          IF(IRH2.NE.INEXT)THEN
            DO INH1=1,NONZ1
              IRH1=NONZC1(3,INH1)
              ILH1=NONZC1(4,INH1)
              IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
                DO ILR=1,NLVALR
                  CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                  IRHS=(ILR-1)*KCSZR+IRH1
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO ILL=1,NLVALL
                    ILHS=(ILL-1)*KCSZL+ILH1
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK3(ILL,INH1)
                  END DO
                END DO
C**OFF-DIAGONAL BLOCK - END
              ELSE
C**DIAGONAL BLOCK - START
                DO ILR=1,NLVALR
                  CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                  IRHS=(ILR-1)*KCSZR+IRH1
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO ILL=1,ILR
                    IF(ILL.EQ.ILR.AND.ILH1.GT.IRH1)GO TO 2004
                    ILHS=(ILL-1)*KCSZL+ILH1
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK3(ILL,INH1)
2004                CONTINUE
                  END DO
                END DO
C**DIAGONAL BLOCK - END
              END IF
            END DO
          END IF
          IPREV=IRH2
        END DO
        GO TO 5555

4004    CONTINUE
        DO IRH2=1,LCSZR
          DO ILH2=1,LCSZL
            XK1(ILH2,IRH2)=0
          END DO
        END DO

C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        DO INH1=1,NONZ1
          IRH1=NONZC1(3,INH1)
          ILH1=NONZC1(4,INH1)
          NR1=NONZC1(1,INH1)
          NL1=NONZC1(2,INH1)

C**CONTRACTION SCHEME INVOLVING 'MOD2' (L)
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            ILH2=NONZC2(4,INH2)
            NR2=NONZC2(1,INH2)
            NL2=NONZC2(2,INH2)
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1001
            END DO
1001        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IF(NL1.EQ.IP2(IL,1).AND.NL2.EQ.IP2(IL,2))GO TO 2001
            END DO
2001        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA2(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              DO I=1,NP2(IND2)
                K=JP2(I,IND2,1)+1
                L=JP2(I,IND2,2)+1
                ZYX=ZYX+CP2(I,IND2)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)
              END DO
C**ANALYTIC
            END IF
            XK1(ILH2,IRH2)=XYZ+ZYX
          END DO

C************************************DGEMM (RHS)
          CALL DGEMM('N','N',LCSZL,NLVALR,LCSZR,1.0D0,XK1(1,1),
     &    KTEMP,CFS2(1,1,KEL,LSMR),ISIZXX,0.0D0,TEMP(1,1),KTEMP)

C************************************DGEMM (LHS)
          CALL DGEMM('T','N',NLVALL,NLVALR,LCSZL,1.D0,
     &    CFS2(1,1,KEL,LSML),ISIZXX,TEMP,KTEMP,0.0D0,XCON1(1,1),MAXVAL)

          IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
            DO ILR=1,NLVALR
              IRHS=(ILR-1)*KCSZR+IRH1
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO ILL=1,NLVALL
                ILHS=(ILL-1)*KCSZL+ILH1
                XA(ILHS+J0)=XA(ILHS+J0)+XCON1(ILL,ILR)
              END DO
            END DO
C**OFF-DIAGONAL BLOCK - END
          ELSE
C**DIAGONAL BLOCK - START
            DO ILR=1,NLVALR
              IRHS=(ILR-1)*KCSZR+IRH1
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO ILL=1,ILR
                IF(ILL.EQ.ILR.AND.ILH1.GT.IRH1)GO TO 2002
                ILHS=(ILL-1)*KCSZL+ILH1
                XA(ILHS+J0)=XA(ILHS+J0)+XCON1(ILL,ILR)
2002            CONTINUE
              END DO
            END DO
C**DIAGONAL BLOCK - END
          END IF
        END DO
      END IF

5555  CONTINUE
5556  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM2A.EQ.0)THEN
        IF(IND.NE.0)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM2A=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCV2A(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,ITMODE,H1,
     1XQ1,H2,XQ2,HTAU,XQTAU,NN1,MM1,NN2,MM2,NNTAU,MMTAU,IP,ISIZMX,IPC,
     2ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP3,ISIZE3,XA3,XRA3,NSIZE,XK,TEMP,
     3XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCIM,X0,Y0,T0,VP,VPR,
     4VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,
     5NVSYMX,I16)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,10),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,10),VRR(J21,MM2,MM1)
      REAL*4 XRA3(1)
      DIMENSION X(16),Y(16),C(16)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC(IPSIZE,1),IP3(ISIZE3,3)
      DIMENSION XA(1),XA3(1),XK(ICSIZE,ICSIZE)
      DIMENSION TEMP(ICSIZE,NVAL),XCON(NVAL,NVAL)
      DIMENSION TEMP1(I16,JCI1,JCI2,JCI1,JCI2),TEMP2(I16,JCI2,JCI2)
      DIMENSION X0(I16,JCI1,JCI1,MM1),Y0(I16,JCI2,JCI2,MM2)
      DIMENSION T0(I16,JCIM,JCIM,MMTAU)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      CALL VDCV2A(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,ITMODE,H1,XQ1,H2,
     2XQ2,HTAU,XQTAU,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,NNTAU,MMTAU,IP,
     3ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP3,ISIZE3,XA3,XRA3,
     4NSIZE,XK,
     4TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCIM,X0,Y0,T0,
     5VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,
     6KEL21,NVSYMX,I16)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCV2A(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,ITMODE,H1,
     1XQ1,H2,XQ2,HTAU,XQTAU,NN1,MH1,MM1,NN2,MH2,MM2,NNTAU,MMTAU,IP,
     2ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP3,ISIZE3,XA3,XRA3,
     3NSIZE,XK,
     3TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCIM,X0,Y0,T0,
     4VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,
     5KEL21,NVSYMX,I16)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,10),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,10),VRR(J21,MM2,MM1)
      REAL*4 XRA3(1)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      DIMENSION X(16),Y(16),C(16)
      DIMENSION MODINT(NMODE)
C     DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC(IPSIZE,1),IP3(ISIZE3,3)
      DIMENSION XA(1),XA3(1),XK(ICSIZE,ICSIZE)
      DIMENSION TEMP(ICSIZE,NVAL),XCON(NVAL,NVAL)
      DIMENSION TEMP1(I16,JCI1,JCI2,JCI1,JCI2),TEMP2(I16,JCI2,JCI2)
      DIMENSION X0(I16,JCI1,JCI1,MM1),Y0(I16,JCI2,JCI2,MM2)
      DIMENSION T0(I16,JCIM,JCIM,MMTAU)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM2A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VCCV2A'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1', 'MOD2'
C** AND 'TAU'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0CV2
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=JNTFAC(NMODE,ICOUPL,2)
      IFACTC=JNTFAC(NMODE,ICOUPC,2)
C**INTRINSIC
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1

      KA=KROT/2
      LMAX=1+MOD(KA,2)
      FACTRC=0.D0
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
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      MD=MD1*MD2*MDT

C**FORM INDIVIDUAL INTEGRATION TERMS (START)
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          DO NRR=1,JCIM
            NR=NRR+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,16
              X(I)=X(1)
            END DO
            X(5)=HTAU(NR,MTAU,2,LMAX)*MD
            X(10)=X(5)
            X(13)=X(5)
            X(15)=X(5)
            DO NLL=1,JCIM
              NL=NLL+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO I=2,16
                Y(I)=Y(1)
              END DO
              Y(6)=HTAU(NL,MTAU,2,LMAX)
              Y(11)=Y(6)
              Y(14)=Y(6)
              Y(15)=Y(6)
              DO K=1,I16
                T0(K,NLL,NRR,MTAU)=Y(17-K)*X(17-K)
              END DO
            END DO
          END DO
        ELSE
          DO NRR=1,JCIM
            NR=2*NRR-MOD(NRR,2)+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,10
              X(I)=X(1)
            END DO
            X(3)=HTAU(NR,MTAU,2,LMAX)*MD
            X(6)=X(3)
            X(8)=X(3)
            X(9)=HTAU(NR,MTAU,3,LMAX)*MD
            X(11)=HTAU(NR,MTAU,1,LMAX)*FACTRC*MD
            DO NLL=1,JCIM
              NL=2*NLL-MOD(NLL,2)+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO K=1,11
                T0(K,NLL,NRR,MTAU)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M1=1,MM1/MD1
      DO M1=1,MM1
        IF(.NOT.LINEAR)THEN
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,16
              X(I)=X(2)
            END DO
            X(7)=X(1)
            X(9)=X(1)
            X(11)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              Y(2)=H1(NL1,M1,2)
              DO I=3,16
                Y(I)=Y(1)
              END DO
              Y(7)=Y(2)
              Y(8)=Y(2)
              Y(10)=Y(2)
              DO K=1,I16
                X0(K,NL1,NR1,M1)=Y(17-K)*X(17-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,11
              X(I)=X(2)
            END DO
            X(4)=H1(NR1,M1,3)
            X(5)=X(1)
            X(6)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              DO K=1,11
                X0(K,NL1,NR1,M1)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M2=1,MM2/MD2
      DO M2=1,MM2
        IF(.NOT.LINEAR)THEN
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            DO I=2,16
              X(I)=X(1)
            END DO
            X(3)=H2(NR2,M2,2)
            X(8)=X(3)
            X(12)=X(3)
            X(14)=X(3)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO I=2,16
                Y(I)=Y(1)
              END DO
              Y(4)=H2(NL2,M2,2)
              Y(9)=Y(4)
              Y(12)=Y(4)
              Y(13)=Y(4)
              DO K=1,I16
                Y0(K,NL2,NR2,M2)=Y(17-K)*X(17-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            X(2)=H2(NR2,M2,2)
            DO I=3,11
              X(I)=X(1)
            END DO
            X(5)=X(2)
            X(7)=H2(NR2,M2,3)
            X(8)=X(2)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO K=1,11
                Y0(K,NL2,NR2,M2)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
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
          IF(J21.GT.1.AND.ICOUPC.GE.2)READ(62)VR
          IF(ICOUPC.GE.2)READ(82)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPC.GE.2)READ(62)VRR
          IF(ICOUPC.GE.2)READ(82)VCR
        END IF
        IF(JCOUPL.GT.0)THEN
          READ(72)VP
        ELSE
          READ(72)VPR
        END IF

C***********************************************************

        IF(.NOT.LINEAR)THEN
          DO IRHS2=1,JCI2
            DO IRHS1=1,JCI1
              DO ILHS2=1,JCI2
                DO ILHS1=1,JCI1
                  DO K=1,I16
                    TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)=0.D0
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 2-MODE INTEGRATION
C         DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS2=1,JCI2
              DO ILHS2=1,JCI2
                DO K=1,I16
                  TEMP2(K,ILHS2,IRHS2)=0.D0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
C           DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO I=1,16
                C(I)=0
              END DO
              IF(JCOUPL.GT.0)THEN
                IF(ICOUPC.GE.2)THEN
C**NO WATSON TERM IF RPH
                  C(1)=VC(M2,M1,1)*IFACTC
                  C(2)=VC(M2,M1,1)*IFACTC
                  C(3)=VC(M2,M1,2)*IFACTC
                  C(4)=VC(M2,M1,2)*IFACTC
                  C(5)=VC(M2,M1,3)*IFACTC
                  C(6)=VC(M2,M1,3)*IFACTC
                  C(7)=VC(M2,M1,4)*IFACTC
                  C(8)=VC(M2,M1,5)*IFACTC
                  C(9)=VC(M2,M1,5)*IFACTC
                  C(10)=VC(M2,M1,6)*IFACTC
                  C(11)=VC(M2,M1,6)*IFACTC
                  C(12)=VC(M2,M1,7)*IFACTC
                  C(13)=VC(M2,M1,8)*IFACTC
                  C(14)=VC(M2,M1,8)*IFACTC
                  C(15)=VC(M2,M1,9)*IFACTC
                  C(16)=VC(M2,M1,10)*IFACTC
                  IF(J21.GT.1)C(16)=C(16)+VR(KROT,M2,M1)*IFACTC
                END IF
                C(16)=C(16)+VP(M2,M1)*IFACTL
              ELSE
                IF(ICOUPC.GE.2)THEN
C**NO WATSON TERM IF RPH
                  C(1)=VCR(M2,M1,1)*IFACTC
                  C(2)=VCR(M2,M1,1)*IFACTC
                  C(3)=VCR(M2,M1,2)*IFACTC
                  C(4)=VCR(M2,M1,2)*IFACTC
                  C(5)=VCR(M2,M1,3)*IFACTC
                  C(6)=VCR(M2,M1,3)*IFACTC
                  C(7)=VCR(M2,M1,4)*IFACTC
                  C(8)=VCR(M2,M1,5)*IFACTC
                  C(9)=VCR(M2,M1,5)*IFACTC
                  C(10)=VCR(M2,M1,6)*IFACTC
                  C(11)=VCR(M2,M1,6)*IFACTC
                  C(12)=VCR(M2,M1,7)*IFACTC
                  C(13)=VCR(M2,M1,8)*IFACTC
                  C(14)=VCR(M2,M1,8)*IFACTC
                  C(15)=VCR(M2,M1,9)*IFACTC
                  C(16)=VCR(M2,M1,10)*IFACTC
                  IF(J21.GT.1)C(16)=C(16)+VRR(KROT,M2,M1)*IFACTC
                END IF
                C(16)=C(16)+VPR(M2,M1)*IFACTL
              END IF
              DO IRHS2=1,JCI2
                DO ILHS2=1,JCI2
                  DO K=1,I16
                    TEMP2(K,ILHS2,IRHS2)=TEMP2(K,ILHS2,IRHS2)+
     1              Y0(K,ILHS2,IRHS2,M2)*C(17-K)
                  END DO
                END DO
              END DO
            END DO
C**END 1-MODE INTEGRATION
            DO IRHS2=1,JCI2
              DO IRHS1=1,JCI1
                DO ILHS2=1,JCI2
                  DO ILHS1=1,JCI1
                    DO K=1,I16
                      TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)=
     1                TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)+
     2                X0(K,ILHS1,IRHS1,M1)*TEMP2(K,ILHS2,IRHS2)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 2-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
          DO IRHS=1,NSIZE
            NR1=IP3(IRHS,1)
            NR2=IP3(IRHS,2)
            IRTAU=IP3(IRHS,3)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP3(ILHS,1)
              NL2=IP3(ILHS,2)
              ILTAU=IP3(ILHS,3)
              DO K=1,I16
                XA3(ILHS+J0)=XA3(ILHS+J0)+TEMP1(K,NL1,NL2,NR1,NR2)*
     1          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
              END DO
            END DO
          END DO
        ELSE
          DO IRHS2=1,JCI2
            DO IRHS1=1,JCI1
              DO ILHS2=1,JCI2
                DO ILHS1=1,JCI1
                  DO K=1,11
                    TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)=0.D0
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 2-MODE INTEGRATION
C         DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS2=1,JCI2
              DO ILHS2=1,JCI2
                DO K=1,11
                  TEMP2(K,ILHS2,IRHS2)=0.D0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
C           DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO I=1,9
                C(I)=0
              END DO
              IF(JCOUPL.GT.0)THEN
                IF(ICOUPC.GE.2)THEN
C**(10) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                  C(1)=VC(M2,M1,1)
                  C(2)=VC(M2,M1,2)
                  C(3)=VC(M2,M1,3)
                  C(4)=VC(M2,M1,4)
                  C(5)=2*VC(M2,M1,5)
                  C(6)=2*VC(M2,M1,6)
                  C(7)=VC(M2,M1,7)
                  C(8)=2*VC(M2,M1,8)
                  C(9)=VC(M2,M1,9)
                END IF
                DO IRHS2=1,JCI2
                  DO ILHS2=1,JCI2
                    DO K=1,9
                      TEMP2(K,ILHS2,IRHS2)=TEMP2(K,ILHS2,IRHS2)-
     1                Y0(K,ILHS2,IRHS2,M2)*C(K)*IFACTC
                    END DO
                    TEMP2(10,ILHS2,IRHS2)=TEMP2(10,ILHS2,IRHS2)+
     1              Y0(10,ILHS2,IRHS2,M2)*VP(M2,M1)*IFACTL
                    TEMP2(11,ILHS2,IRHS2)=TEMP2(11,ILHS2,IRHS2)+
     1              Y0(11,ILHS2,IRHS2,M2)*VR(KROT,M2,M1)
                  END DO
                END DO
              ELSE
                IF(ICOUPC.GE.2)THEN
C**(10) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                  C(1)=VCR(M2,M1,1)
                  C(2)=VCR(M2,M1,2)
                  C(3)=VCR(M2,M1,3)
                  C(4)=VCR(M2,M1,4)
                  C(5)=2*VCR(M2,M1,5)
                  C(6)=2*VCR(M2,M1,6)
                  C(7)=VCR(M2,M1,7)
                  C(8)=2*VCR(M2,M1,8)
                  C(9)=VCR(M2,M1,9)
                END IF
                DO IRHS2=1,JCI2
                  DO ILHS2=1,JCI2
                    DO K=1,9
                      TEMP2(K,ILHS2,IRHS2)=TEMP2(K,ILHS2,IRHS2)-
     1                Y0(K,ILHS2,IRHS2,M2)*C(K)*IFACTC
                    END DO
                    TEMP2(10,ILHS2,IRHS2)=TEMP2(10,ILHS2,IRHS2)+
     1              Y0(10,ILHS2,IRHS2,M2)*VPR(M2,M1)*IFACTL
                    TEMP2(11,ILHS2,IRHS2)=TEMP2(11,ILHS2,IRHS2)+
     1              Y0(11,ILHS2,IRHS2,M2)*VRR(KROT,M2,M1)
                  END DO
                END DO
              END IF
            END DO
C**END 1-MODE INTEGRATION
            DO IRHS2=1,JCI2
              DO IRHS1=1,JCI1
                DO ILHS2=1,JCI2
                  DO ILHS1=1,JCI1
                    DO K=1,11
                      TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)=
     1                TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)+
     2                X0(K,ILHS1,IRHS1,M1)*TEMP2(K,ILHS2,IRHS2)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 2-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
C*************************************************************
C**
C**                     NOTE*NOTE*NOTE
C**
C**FOR STRETCH-BEND K WILL DENOTE BEND, L WILL DENOTE STRETCH
C**            K IS INDEX 1        L IS INDEX 2
C*************************************************************
          DO IRHS=1,NSIZE
            NR1=IP3(IRHS,1)
            NR2=IP3(IRHS,2)
            IRRTAU=IP3(IRHS,3)
            IRTAU=IRRTAU/2+MOD(IRRTAU,2)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP3(ILHS,1)
              NL2=IP3(ILHS,2)
              ILLTAU=IP3(ILHS,3)
              ILTAU=ILLTAU/2+MOD(ILLTAU,2)
              DO K=1,11
                XA3(ILHS+J0)=XA3(ILHS+J0)+TEMP1(K,NL1,NL2,NR1,NR2)*
     1          T0(K,ILTAU,IRTAU,MTAU)
              END DO
            END DO
          END DO
        END IF
      END DO
C**END TAU LOOP (3-MODE INTEGRATION)

      CALL MATOUT(XA3,XRA3,NSIZE,22)
      GO TO 7777
6666  CALL MATIN(XA3,XRA3,NSIZE,22)
7777  CONTINUE

C*************************************************************
C**
C**                     NOTE*NOTE*NOTE
C**
C**FOR STRETCH-BEND K WILL DENOTE BEND, L WILL DENOTE STRETCH
C**            K IS INDEX 1        L IS INDEX 2
C*************************************************************
C***********************************ALGORITHM FROM VCV2

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR1=IPC(IROFF,MODE1)
        NR2=IPC(IROFF,MODE2)
        NRTAU=IPC(IROFF,ITMODE)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP3(IR,1).AND.NR2.EQ.IP3(IR,2).AND.
     1       NRTAU.EQ.IP3(IR,3))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.ITMODE.AND.(
     1      IPC(IROFF,K).NE.IPC(ILOFF,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPC(ILOFF,MODE1)
          NL2=IPC(ILOFF,MODE2)
          NLTAU=IPC(ILOFF,ITMODE)
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
          XYZ=XA3(I)
3000      CONTINUE
          XK(ILHS,IRHS)=XYZ*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
        J0=(IROFF)*(IROFF-1)/2+IOFF
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IP(IROFF,KCONT)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IP(IROFF,MCONT).NE.IP(ILOFF,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IP(ILOFF,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE

C**TEMPORARY
C     WRITE(IOUT,*)'VCCV2A MATRIX'
C     WRITE(IOUT,*)(XA(I),I=1,20)
C     WRITE(IOUT,*)
C     WRITE(IOUT,*)(XA(ISIZE-21+I),I=1,20)
C**TEMPORARY

      IF(ITIM2A.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM2A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCV2B(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,MODE3,H1,
     1XQ1,H2,XQ2,HTAU,XQTAU,NN1,MM1,NN2,MM2,NNTAU,MMTAU,IP,ISIZMX,IPC1,
     2ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP3,
     3ISIZE3,XA3,XRA3,
     3NSIZE,IRL,XK2,TEMP,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,
     4TEMP2,JCI1,JCI2,JCIM,X0,Y0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,
     5MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,
     6NONZC1,NONZC2,I16,MAXVAL,MAXQ3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,10),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,10),VRR(J21,MM2,MM1)
      REAL*4 XRA3(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP3(ISIZE3,3),TEMP(KTEMP,1)
      DIMENSION TEMP1(I16,JCI1,JCI2,JCI1,JCI2),TEMP2(I16,JCI2,JCI2)
      DIMENSION X0(I16,JCI1,JCI1,MM1),Y0(I16,JCI2,JCI2,MM2)
      DIMENSION T0(I16,JCIM,JCIM,MMTAU)
      DIMENSION X(16),Y(16),C(16)
C**************************************************************
      DIMENSION XA(1),XA3(1)
      DIMENSION IRL(MAXQ3,MAXQ3,MAXQ3),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C**************************************************************
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      CALL VDCV2B(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,MODE3,H1,XQ1,H2,
     1XQ2,HTAU,XQTAU,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,NNTAU,MMTAU,IP,
     2ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,
     3ISIZE,IP3,ISIZE3,
     3XA3,XRA3,NSIZE,IRL,XK2,TEMP,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,
     4TEMP1,TEMP2,JCI1,JCI2,JCIM,X0,Y0,T0,VP,VPR,VC,VCR,VR,VRR,J21,
     5KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,
     6NONZC1,NONZC2,I16,MAXVAL,MAXQ3)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCV2B(NMODE,NNMODE,MOD1,MOD2,MODE1,MODE2,MODE3,H1,
     1XQ1,H2,XQ2,HTAU,XQTAU,NN1,MH1,MM1,NN2,MH2,MM2,NNTAU,MMTAU,IP,
     2ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,
     3ISIZE,IP3,ISIZE3,
     3XA3,XRA3,NSIZE,IRL,XK2,TEMP,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,
     4TEMP1,TEMP2,JCI1,JCI2,JCIM,X0,Y0,T0,VP,VPR,VC,VCR,VR,VRR,J21,
     5KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,
     6NONZC1,NONZC2,I16,MAXVAL,MAXQ3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,10),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,10),VRR(J21,MM2,MM1)
      REAL*4 XRA3(1)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD

C**TEMPORARY
      DIMENSION IRL(MAXQ3,MAXQ3,MAXQ3),NX3(3),MX3(3)
C**TEMPORARY

      DIMENSION MODINT(NMODE)
C     DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP3(ISIZE3,3),TEMP(KTEMP,1)
      DIMENSION TEMP1(I16,JCI1,JCI2,JCI1,JCI2),TEMP2(I16,JCI2,JCI2)
      DIMENSION X0(I16,JCI1,JCI1,MM1),Y0(I16,JCI2,JCI2,MM2)
      DIMENSION T0(I16,JCIM,JCIM,MMTAU)
      DIMENSION X(16),Y(16),C(16)
C*****************************************************************
      DIMENSION XA(1),XA3(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C*****************************************************************
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM2B.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCV2B'
          CALL TIMIT(1)
          CALL FLUSH(IOUT)
        END IF
      END IF

C**************************************************************
C**************************************************************
C**ONE CONTRACTION SCHEME HAS A SINGLE MODE
C**ONE CONTRACTION SCHEME HAS TWO MODES
C**CHECK SCHEME '1'
      I1=KCONT1
      NC1=0
      DO NN=1,ICONT(I1)
        IF(MOD1.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE1
          JC1(NC1)=1
        END IF
        IF(MOD2.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE2
          JC1(NC1)=2
        END IF
        IF(NNMODE.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE3
          JC1(NC1)=3
        END IF
      END DO
C**CHECK SCHEME '2'
      I2=KCONT2
      NC2=0
      DO NN=1,ICONT(I2)
        IF(MOD1.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE1
          JC2(NC2)=1
        END IF
        IF(MOD2.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE2
          JC2(NC2)=2
        END IF
        IF(NNMODE.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE3
          JC2(NC2)=3
        END IF
      END DO
C**NKMOD1 IS TOTAL NUMBER OF MODES FOR SCHEME INVOLVING 'MOD1' (K)
C**NKMOD2 IS TOTAL NUMBER OF MODES FOR THE OTHER SCHEME (ONE OR BOTH
C**OF 'MOD2' (L) AND 'MOD3' (N))
C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L) OR
C**'MOD3' (N)
C**IF THEY ARE IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
C**************************************************************
      IF(IND.NE.0)GO TO 7777
C**************************************************************

C**FIND MAX QUANTUM
      MAXQ3=0
      DO IR=1,NSIZE
        NR1=IP3(IR,1)
        NR2=IP3(IR,2)
        NR3=IP3(IR,3)
        IF(NR1.GT.MAXQ3)MAXQ3=NR1
        IF(NR2.GT.MAXQ3)MAXQ3=NR2
        IF(NR3.GT.MAXQ3)MAXQ3=NR3
      END DO

C***********************************ALGORITHM FROM V0CV2
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=JNTFAC(NMODE,ICOUPL,2)
      IFACTC=JNTFAC(NMODE,ICOUPC,2)
C**INTRINSIC
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1

      KA=KROT/2
      LMAX=1+MOD(KA,2)
      FACTRC=0.D0
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      MD=MD1*MD2*MDT

C**FORM INDIVIDUAL INTEGRATION TERMS (START)
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          DO NRR=1,JCIM
            NR=NRR+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,16
              X(I)=X(1)
            END DO
            X(5)=HTAU(NR,MTAU,2,LMAX)*MD
            X(10)=X(5)
            X(13)=X(5)
            X(15)=X(5)
            DO NLL=1,JCIM
              NL=NLL+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO I=2,16
                Y(I)=Y(1)
              END DO
              Y(6)=HTAU(NL,MTAU,2,LMAX)
              Y(11)=Y(6)
              Y(14)=Y(6)
              Y(15)=Y(6)
              DO K=1,I16
                T0(K,NLL,NRR,MTAU)=Y(17-K)*X(17-K)
              END DO
            END DO
          END DO
        ELSE
          DO NRR=1,JCIM
            NR=2*NRR-MOD(NRR,2)+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,10
              X(I)=X(1)
            END DO
            X(3)=HTAU(NR,MTAU,2,LMAX)*MD
            X(6)=X(3)
            X(8)=X(3)
            X(9)=HTAU(NR,MTAU,3,LMAX)*MD
            X(11)=HTAU(NR,MTAU,1,LMAX)*FACTRC*MD
            DO NLL=1,JCIM
              NL=2*NLL-MOD(NLL,2)+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO K=1,11
                T0(K,NLL,NRR,MTAU)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M1=1,MM1/MD1
      DO M1=1,MM1
        IF(.NOT.LINEAR)THEN
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,16
              X(I)=X(2)
            END DO
            X(7)=X(1)
            X(9)=X(1)
            X(11)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              Y(2)=H1(NL1,M1,2)
              DO I=3,16
                Y(I)=Y(1)
              END DO
              Y(7)=Y(2)
              Y(8)=Y(2)
              Y(10)=Y(2)
              DO K=1,I16
                X0(K,NL1,NR1,M1)=Y(17-K)*X(17-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,11
              X(I)=X(2)
            END DO
            X(4)=H1(NR1,M1,3)
            X(5)=X(1)
            X(6)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              DO K=1,11
                X0(K,NL1,NR1,M1)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M2=1,MM2/MD2
      DO M2=1,MM2
        IF(.NOT.LINEAR)THEN
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            DO I=2,16
              X(I)=X(1)
            END DO
            X(3)=H2(NR2,M2,2)
            X(8)=X(3)
            X(12)=X(3)
            X(14)=X(3)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO I=2,16
                Y(I)=Y(1)
              END DO
              Y(4)=H2(NL2,M2,2)
              Y(9)=Y(4)
              Y(12)=Y(4)
              Y(13)=Y(4)
              DO K=1,I16
                Y0(K,NL2,NR2,M2)=Y(17-K)*X(17-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            X(2)=H2(NR2,M2,2)
            DO I=3,11
              X(I)=X(1)
            END DO
            X(5)=X(2)
            X(7)=H2(NR2,M2,3)
            X(8)=X(2)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO K=1,11
                Y0(K,NL2,NR2,M2)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
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
          IF(J21.GT.1.AND.ICOUPL.GE.2)READ(62)VR
          IF(ICOUPC.GE.2)READ(82)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPL.GE.2)READ(62)VRR
          IF(ICOUPC.GE.2)READ(82)VCR
        END IF
        IF(JCOUPL.GT.0)THEN
          READ(72)VP
        ELSE
          READ(72)VPR
        END IF

C***********************************************************

        IF(.NOT.LINEAR)THEN
          DO IRHS2=1,JCI2
            DO IRHS1=1,JCI1
              DO ILHS2=1,JCI2
                DO ILHS1=1,JCI1
                  DO K=1,I16
                    TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)=0.D0
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 2-MODE INTEGRATION
C         DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS2=1,JCI2
              DO ILHS2=1,JCI2
                DO K=1,I16
                  TEMP2(K,ILHS2,IRHS2)=0.D0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
C           DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO I=1,16
                C(I)=0
              END DO
              IF(JCOUPL.GT.0)THEN
                IF(ICOUPC.GE.2)THEN
C**NO WATSON TERM IF RPH
                  C(1)=VC(M2,M1,1)*IFACTC
                  C(2)=VC(M2,M1,1)*IFACTC
                  C(3)=VC(M2,M1,2)*IFACTC
                  C(4)=VC(M2,M1,2)*IFACTC
                  C(5)=VC(M2,M1,3)*IFACTC
                  C(6)=VC(M2,M1,3)*IFACTC
                  C(7)=VC(M2,M1,4)*IFACTC
                  C(8)=VC(M2,M1,5)*IFACTC
                  C(9)=VC(M2,M1,5)*IFACTC
                  C(10)=VC(M2,M1,6)*IFACTC
                  C(11)=VC(M2,M1,6)*IFACTC
                  C(12)=VC(M2,M1,7)*IFACTC
                  C(13)=VC(M2,M1,8)*IFACTC
                  C(14)=VC(M2,M1,8)*IFACTC
                  C(15)=VC(M2,M1,9)*IFACTC
                  C(16)=VC(M2,M1,10)*IFACTC
                  IF(J21.GT.1)C(16)=C(16)+VR(KROT,M2,M1)*IFACTC
                END IF
                C(16)=C(16)+VP(M2,M1)*IFACTL
              ELSE
                IF(ICOUPC.GE.2)THEN
C**NO WATSON TERM IF RPH
                  C(1)=VCR(M2,M1,1)*IFACTC
                  C(2)=VCR(M2,M1,1)*IFACTC
                  C(3)=VCR(M2,M1,2)*IFACTC
                  C(4)=VCR(M2,M1,2)*IFACTC
                  C(5)=VCR(M2,M1,3)*IFACTC
                  C(6)=VCR(M2,M1,3)*IFACTC
                  C(7)=VCR(M2,M1,4)*IFACTC
                  C(8)=VCR(M2,M1,5)*IFACTC
                  C(9)=VCR(M2,M1,5)*IFACTC
                  C(10)=VCR(M2,M1,6)*IFACTC
                  C(11)=VCR(M2,M1,6)*IFACTC
                  C(12)=VCR(M2,M1,7)*IFACTC
                  C(13)=VCR(M2,M1,8)*IFACTC
                  C(14)=VCR(M2,M1,8)*IFACTC
                  C(15)=VCR(M2,M1,9)*IFACTC
                  C(16)=VCR(M2,M1,10)*IFACTC
                  IF(J21.GT.1)C(16)=C(16)+VRR(KROT,M2,M1)*IFACTC
                END IF
                C(16)=C(16)+VPR(M2,M1)*IFACTL
              END IF
              DO IRHS2=1,JCI2
                DO ILHS2=1,JCI2
                  DO K=1,I16
                    TEMP2(K,ILHS2,IRHS2)=TEMP2(K,ILHS2,IRHS2)+
     1              Y0(K,ILHS2,IRHS2,M2)*C(17-K)
                  END DO
                END DO
              END DO
            END DO
C**END 1-MODE INTEGRATION
            DO IRHS2=1,JCI2
              DO IRHS1=1,JCI1
                DO ILHS2=1,JCI2
                  DO ILHS1=1,JCI1
                    DO K=1,I16
                      TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)=
     1                TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)+
     2                X0(K,ILHS1,IRHS1,M1)*TEMP2(K,ILHS2,IRHS2)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 2-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
          DO IRHS=1,NSIZE
            NR1=IP3(IRHS,1)
            NR2=IP3(IRHS,2)
            IRTAU=IP3(IRHS,3)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP3(ILHS,1)
              NL2=IP3(ILHS,2)
              ILTAU=IP3(ILHS,3)
              DO K=1,I16
                XA3(ILHS+J0)=XA3(ILHS+J0)+TEMP1(K,NL1,NL2,NR1,NR2)*
     1          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
              END DO
            END DO
          END DO
        ELSE
          DO IRHS2=1,JCI2
            DO IRHS1=1,JCI1
              DO ILHS2=1,JCI2
                DO ILHS1=1,JCI1
                  DO K=1,11
                    TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)=0.D0
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 2-MODE INTEGRATION
C         DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS2=1,JCI2
              DO ILHS2=1,JCI2
                DO K=1,11
                  TEMP2(K,ILHS2,IRHS2)=0.D0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
C           DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO I=1,9
                C(I)=0
              END DO
              IF(JCOUPL.GT.0)THEN
                IF(ICOUPC.GE.2)THEN
C**(10) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                  C(1)=VC(M2,M1,1)
                  C(2)=VC(M2,M1,2)
                  C(3)=VC(M2,M1,3)
                  C(4)=VC(M2,M1,4)
                  C(5)=2*VC(M2,M1,5)
                  C(6)=2*VC(M2,M1,6)
                  C(7)=VC(M2,M1,7)
                  C(8)=2*VC(M2,M1,8)
                  C(9)=VC(M2,M1,9)
                END IF
                DO IRHS2=1,JCI2
                  DO ILHS2=1,JCI2
                    DO K=1,9
                      TEMP2(K,ILHS2,IRHS2)=TEMP2(K,ILHS2,IRHS2)-
     1                Y0(K,ILHS2,IRHS2,M2)*C(K)*IFACTC
                    END DO
                    TEMP2(10,ILHS2,IRHS2)=TEMP2(10,ILHS2,IRHS2)+
     1              Y0(10,ILHS2,IRHS2,M2)*VP(M2,M1)*IFACTL
                    TEMP2(11,ILHS2,IRHS2)=TEMP2(11,ILHS2,IRHS2)+
     1              Y0(11,ILHS2,IRHS2,M2)*VR(KROT,M2,M1)
                  END DO
                END DO
              ELSE
                IF(ICOUPC.GE.2)THEN
C**(10) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                  C(1)=VCR(M2,M1,1)
                  C(2)=VCR(M2,M1,2)
                  C(3)=VCR(M2,M1,3)
                  C(4)=VCR(M2,M1,4)
                  C(5)=2*VCR(M2,M1,5)
                  C(6)=2*VCR(M2,M1,6)
                  C(7)=VCR(M2,M1,7)
                  C(8)=2*VCR(M2,M1,8)
                  C(9)=VCR(M2,M1,9)
                END IF
                DO IRHS2=1,JCI2
                  DO ILHS2=1,JCI2
                    DO K=1,9
                      TEMP2(K,ILHS2,IRHS2)=TEMP2(K,ILHS2,IRHS2)-
     1                Y0(K,ILHS2,IRHS2,M2)*C(K)*IFACTC
                    END DO
                    TEMP2(10,ILHS2,IRHS2)=TEMP2(10,ILHS2,IRHS2)+
     1              Y0(10,ILHS2,IRHS2,M2)*VPR(M2,M1)*IFACTL
                    TEMP2(11,ILHS2,IRHS2)=TEMP2(11,ILHS2,IRHS2)+
     1              Y0(11,ILHS2,IRHS2,M2)*VRR(KROT,M2,M1)
                  END DO
                END DO
              END IF
            END DO
C**END 1-MODE INTEGRATION
            DO IRHS2=1,JCI2
              DO IRHS1=1,JCI1
                DO ILHS2=1,JCI2
                  DO ILHS1=1,JCI1
                    DO K=1,11
                      TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)=
     1                TEMP1(K,ILHS1,ILHS2,IRHS1,IRHS2)+
     2                X0(K,ILHS1,IRHS1,M1)*TEMP2(K,ILHS2,IRHS2)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 2-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
C*************************************************************
C**
C**                     NOTE*NOTE*NOTE
C**
C**FOR STRETCH-BEND K WILL DENOTE BEND, L WILL DENOTE STRETCH
C**            K IS INDEX 1        L IS INDEX 2
C*************************************************************
          DO IRHS=1,NSIZE
            NR1=IP3(IRHS,1)
            NR2=IP3(IRHS,2)
            IRRTAU=IP3(IRHS,3)
            IRTAU=IRRTAU/2+MOD(IRRTAU,2)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP3(ILHS,1)
              NL2=IP3(ILHS,2)
              ILLTAU=IP3(ILHS,3)
              ILTAU=ILLTAU/2+MOD(ILLTAU,2)
              DO K=1,11
                XA3(ILHS+J0)=XA3(ILHS+J0)+TEMP1(K,NL1,NL2,NR1,NR2)*
     1          T0(K,ILTAU,IRTAU,MTAU)
              END DO
            END DO
          END DO
        END IF
      END DO
C**END TAU LOOP (3-MODE INTEGRATION)

      CALL MATOUT(XA3,XRA3,NSIZE,22)
      GO TO 7777
6666  CALL MATIN(XA3,XRA3,NSIZE,22)
7777  CONTINUE

      IF(IND.NE.0)THEN
        DO IR=1,NSIZE
          NR=0
          DO J=1,NC1
            J1=JC1(J)
            NR1=IP3(IR,J1)
            NR=NR+1
            NX3(NR)=NR1
          END DO
          DO J=1,NC2
            J2=JC2(J)
            NR2=IP3(IR,J2)
            NR=NR+1
            NX3(NR)=NR2
          END DO
          IRL(NX3(1),NX3(2),NX3(3))=IR
        END DO
      END IF

C*************************************************************
C**
C**                     NOTE*NOTE*NOTE
C**
C**FOR STRETCH-BEND K WILL DENOTE BEND, L WILL DENOTE STRETCH
C**            K IS INDEX 1        L IS INDEX 2
C*************************************************************
C***********************************ALGORITHM FROM VCV2

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFFR=IOFFR+NCSIZE(I)
      END DO
      J0R=IOFFR*(IOFFR+1)/2
      CALL FLUSH(IOUT)
C**TEMPORARY-LAN
      JSR=1
C**TEMPORARY-LAN
C     DO 9999 ISMR1=1,NVSYM
      DO 9999 ISMR1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE

        IF(NONC1.GT.NONC2)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          KCCCL=KCSZL
CC        IF(NONC1.GT.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             KCCCL=KCSZL
C           ELSE
C             KCCCL=IRH1
C           END IF
CC        END IF
          DO ILH1=1,KCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=IROFF
            NONZC1(2,NONZ1)=ILOFF
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
        IF(NONZ1.EQ.0)GO TO 5556
        NONZC1(3,NONZ1+1)=0
        END IF

      END IF

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR TAU)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE

        IF(NONC1.LE.NONC2)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR TAU)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          LCCCL=LCSZL
CC        IF(NONC1.LE.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             LCCCL=LCSZL
C           ELSE
C             LCCCL=IRH2
C           END IF
CC        END IF
          DO ILH2=1,LCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=JROFF
            NONZC2(2,NONZ2)=JLOFF
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
        IF(NONZ2.EQ.0)GO TO 5555
        NONZC2(3,NONZ2+1)=0
        END IF

      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM2B.EQ.0)THEN
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C     END IF

      IF(NONC1.LE.NONC2)GO TO 3003
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '2' (L AND/OR TAU)

      DO IRH2=1,LCSZR
        JROFF=IRH2+LCOFFR
        DO I=1,NKVALR
          DO J=1,NKVALL
            DO K=1,NLVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          LCCCL=LCSZL
        ELSE
          LCCCL=IRH2
        END IF
        DO ILH2=1,LCCCL
          JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
          DO K=1,NKMOD2
            IGOT=0
            DO I=1,NC2
              IF(K.EQ.MC2(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 6001
          END DO
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH2.EQ.IRH2)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NKVALR
            DO J=1,NKVALL
              XCON2(J,I)=0
            END DO
          END DO

          IPREV=0
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            INEXT=NONZC1(3,INH1+1)
            ILH1=NONZC1(4,INH1)
            IROFF=NONZC1(1,INH1)
            ILOFF=NONZC1(2,INH1)
            IF(IPREV.NE.IRH1)THEN
C**NEW RHS
              DO J=1,NKVALL
                TEMP(J,1)=0
              END DO
            END IF

            DO I=1,NC1
              NCR1(I)=IPC1(IROFF,MC1(I))
              NCL1(I)=IPC1(ILOFF,MC1(I))
            END DO

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX3(NR)=NR1
        MX3(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX3(NR)=NR2
        MX3(NR)=NL2
      END DO
      IR=IRL(NX3(1),NX3(2),NX3(3))
      IL=IRL(MX3(1),MX3(2),MX3(3))
C**TEMPORARY

C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XXX=XA3(I)

            DO K=1,NKVALL
              TEMP(K,1)=TEMP(K,1)+XXX*CFS1(ILH1,K,KEL,KSML)
            END DO
5001        CONTINUE

            IF(IRH1.NE.INEXT)THEN

            DO K=1,NKVALR
              DO L=1,NKVALL
                XCON2(L,K)=XCON2(L,K)+TEMP(L,1)*CFS1(IRH1,K,KEL,KSMR)
              END DO
            END DO

            END IF
            IPREV=IRH1
          END DO
CC
          DO I=1,NKVALR
            DO J=1,NKVALL
              DO K=1,NLVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6001      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILL,IKR,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKR,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILR,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILR,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO
C*************************************

      GO TO 5555

3003  CONTINUE

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '1' (K)

      DO IRH1=1,KCSZR
        IROFF=IRH1+KCOFFR
        DO I=1,NLVALR
          DO J=1,NLVALL
            DO K=1,NKVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          KCCCL=KCSZL
        ELSE
          KCCCL=IRH1
        END IF
        DO ILH1=1,KCCCL
          ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
          DO K=1,NKMOD1
            IGOT=0
            DO I=1,NC1
              IF(K.EQ.MC1(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 6002
          END DO
          DO I=1,NC1
            NCR1(I)=IPC1(IROFF,MC1(I))
            NCL1(I)=IPC1(ILOFF,MC1(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH1.EQ.IRH1)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NLVALR
            DO J=1,NLVALL
              XCON2(J,I)=0
            END DO
          END DO

C**CONTRACTION SCHEME  '2' (L AND/OR N)

          JPREV=0
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            JNEXT=NONZC2(3,INH2+1)
            ILH2=NONZC2(4,INH2)
            JROFF=NONZC2(1,INH2)
            JLOFF=NONZC2(2,INH2)
            IF(JPREV.NE.IRH2)THEN
C**NEW RHS
              DO J=1,NLVALL
                TEMP(J,1)=0
              END DO
            END IF

            DO I=1,NC2
              NCR2(I)=IPC2(JROFF,MC2(I))
              NCL2(I)=IPC2(JLOFF,MC2(I))
            END DO

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX3(NR)=NR1
        MX3(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX3(NR)=NR2
        MX3(NR)=NL2
      END DO
      IR=IRL(NX3(1),NX3(2),NX3(3))
      IL=IRL(MX3(1),MX3(2),MX3(3))
C**TEMPORARY

C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XXX=XA3(I)

            DO J=1,NLVALL
              TEMP(J,1)=TEMP(J,1)+XXX*CFS2(ILH2,J,KEL,LSML)
            END DO
5002        CONTINUE

            IF(IRH2.NE.JNEXT)THEN

            DO I=1,NLVALR
              DO J=1,NLVALL
                XCON2(J,I)=XCON2(J,I)+TEMP(J,1)*CFS2(IRH2,I,KEL,LSMR)
              END DO
            END DO

            END IF
            JPREV=IRH2
          END DO
CC
          DO I=1,NLVALR
            DO J=1,NLVALL
              DO K=1,NKVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6002      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKR,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKL,ILR,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILR,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO

5555  CONTINUE
5556  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE
C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM2B.EQ.0)THEN
        IF(IND.NE.0)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM2B=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI3A(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,IP,ISIZMX,IPC,
     2ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP3,ISIZE3,XA3,XRA3,NSIZE,XK,TEMP,
     3XCON,
     3NVAL,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,Z0,VP,
     4VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,
     5NVSYMX,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,JP3,NTOT3,MAX3,INDK,INDL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,10),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,10),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA3(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP3(ISIZE3,3),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I10,JCI3,JCI3),TEMP2(I10,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(I10,JCI1,JCI1,MM1),Y0(I10,JCI2,JCI2,MM2)
      DIMENSION Z0(I10,JCI3,JCI3,MM3)
      DIMENSION X(10),Y(10),C(10)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA3(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),JP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
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
      CALL VDCI3A(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,
     2MM3/MD3,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP3,ISIZE3,
     3XA3,XRA3,
     3NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,
     4JCI3,X0,Y0,Z0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,
     5ISIZXX,NVALX,KEL21,NVSYMX,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,JP3,NTOT3,
     6MAX3,INDK,INDL)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI3A(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,IP,
     2ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP3,ISIZE3,XA3,XRA3,
     3NSIZE,XK,
     3TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,
     4Y0,Z0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,
     5NVALX,KEL21,NVSYMX,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,JP3,NTOT3,MAX3,
     6INDK,INDL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,10),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,10),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA3(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP3(ISIZE3,3),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I10,JCI3,JCI3),TEMP2(I10,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(I10,JCI1,JCI1,MM1),Y0(I10,JCI2,JCI2,MM2)
      DIMENSION Z0(I10,JCI3,JCI3,MM3)
      DIMENSION X(10),Y(10),C(10)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA3(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),JP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP

      IF(ITIM3A.EQ.1)THEN
        WRITE(IOUT,*)'Calculating VCCI3A'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1', 'MOD2'
C**AND 'MOD3'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C**ANALYTIC
      IND=INDL(MOD1)+INDK(MOD2)+MOD3
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI3
C**TEMPORARY
      IF(ICONDP.NE.0)THEN
        IF(MOLINC.GT.0)GO TO 99999
        GO TO 6666
      END IF
C**TEMPORARY

      IFACTC=INTFAC(NMODE,ICOUPC,3)
      IFACTL=INTFAC(NMODE,ICOUPL,3)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(73)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(73)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.2)READ(63)VR
        IF(ICOUPC.GT.2)READ(83)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.2)READ(63)VRR
        IF(ICOUPC.GT.2)READ(83)VCR
      END IF

      IF(MOLINC.GT.0)GO TO 99999

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

C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(2)
          X(8)=X(2)
          X(9)=X(2)
          X(10)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(1)
            Y(5)=Y(3)
            Y(6)=Y(3)
            Y(7)=Y(3)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(3)
            DO K=1,I10
              X0(K,NL1,NR1,M1)=Y(11-K)*X(11-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(1)
          X(10)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(3)
            Y(7)=Y(3)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            DO K=1,I10
              Y0(K,NL2,NR2,M2)=Y(11-K)*X(11-K)
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
          X(4)=H3(NR3,M3,2)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(4)
          X(8)=X(1)
          X(9)=X(4)
          X(10)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=H3(NL3,M3,2)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(5)
            Y(9)=Y(5)
            Y(10)=Y(1)
            DO K=1,I10
              Z0(K,NL3,NR3,M3)=Y(11-K)*X(11-K)
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
                DO K=1,I10
                  TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=0.D0
                END DO
              END DO
            END DO
          END DO
        END DO
C**START 2-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO IRHS3=1,JCI3
            DO ILHS3=1,JCI3
              DO K=1,I10
                TEMP1(K,ILHS3,IRHS3)=0.D0
              END DO
            END DO
          END DO
C**START 1-MODE INTEGRATION
CCCC      DO M3=1,MM3/MD3
          DO M3=1,MM3
            DO I=1,10
              C(I)=0
            END DO
            IF(JCOUPL.GT.0)THEN
              IF(ICOUPC.GT.2)THEN
                C(1)=VC(M3,M2,M1,4)*IFACTC
                C(2)=VC(M3,M2,M1,5)*IFACTC
                C(3)=VC(M3,M2,M1,5)*IFACTC
                C(4)=VC(M3,M2,M1,6)*IFACTC
                C(5)=VC(M3,M2,M1,6)*IFACTC
                C(6)=VC(M3,M2,M1,7)*IFACTC
                C(7)=VC(M3,M2,M1,8)*IFACTC
                C(8)=VC(M3,M2,M1,8)*IFACTC
                C(9)=VC(M3,M2,M1,9)*IFACTC
                C(10)=VC(M3,M2,M1,10)*IFACTC
                IF(J21.GT.1)C(10)=C(10)+VR(KROT,M3,M2,M1)*IFACTC
              END IF
              IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(10)=C(10)+
     1        VP(M3,M2,M1)*IFACTL
            ELSE
              IF(ICOUPC.GT.2)THEN
                C(1)=VCR(M3,M2,M1,4)*IFACTC
                C(2)=VCR(M3,M2,M1,5)*IFACTC
                C(3)=VCR(M3,M2,M1,5)*IFACTC
                C(4)=VCR(M3,M2,M1,6)*IFACTC
                C(5)=VCR(M3,M2,M1,6)*IFACTC
                C(6)=VCR(M3,M2,M1,7)*IFACTC
                C(7)=VCR(M3,M2,M1,8)*IFACTC
                C(8)=VCR(M3,M2,M1,8)*IFACTC
                C(9)=VCR(M3,M2,M1,9)*IFACTC
                C(10)=VCR(M3,M2,M1,10)*IFACTC
                IF(J21.GT.1)C(10)=C(10)+VRR(KROT,M3,M2,M1)*IFACTC
              END IF
              IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(10)=C(10)+
     1        VPR(M3,M2,M1)*IFACTL
            END IF
            DO IRHS3=1,JCI3
              DO ILHS3=1,JCI3
                DO K=1,I10
                  TEMP1(K,ILHS3,IRHS3)=TEMP1(K,ILHS3,IRHS3)+
     1            Z0(K,ILHS3,IRHS3,M3)*C(11-K)
                END DO
              END DO
            END DO
          END DO
C**END 1-MODE INTEGRATION
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO ILHS3=1,JCI3
                DO ILHS2=1,JCI2
                  DO K=1,I10
                    TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1              TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2              Y0(K,ILHS2,IRHS2,M2)*TEMP1(K,ILHS3,IRHS3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
C**END 2-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
        DO IRHS=1,NSIZE
          NR1=IP3(IRHS,1)
          NR2=IP3(IRHS,2)
          NR3=IP3(IRHS,3)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP3(ILHS,1)
            NL2=IP3(ILHS,2)
            NL3=IP3(ILHS,3)
            DO K=1,I10
              XA3(ILHS+J0)=XA3(ILHS+J0)+
     1        TEMP2(K,NL2,NL3,NR2,NR3)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 3-MODE INTEGRATION
      CALL MATOUT(XA3,XRA3,NSIZE,23)
      GO TO 7777
6666  CALL MATIN(XA3,XRA3,NSIZE,23)
7777  CONTINUE

C     IF(MOLINC.GT.0)RETURN

C***********************************ALGORITHM FROM VCI3

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA)
      IOFF=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFF=IOFF+NCSIZE(I)
      END DO
      J0R=IOFF*(IOFF+1)/2
C**TEMPORARY-LAN
      JS=1
C**TEMPORARY-LAN
C     DO 9999 ISM1=1,NVSYM
      DO 9999 ISM1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR1=IPC(IROFF,MODE1)
        NR2=IPC(IROFF,MODE2)
        NR3=IPC(IROFF,MODE3)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP3(IR,1).AND.NR2.EQ.IP3(IR,2).AND.
     1       NR3.EQ.IP3(IR,3))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.(
     1      IPC(IROFF,K).NE.IPC(ILOFF,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPC(ILOFF,MODE1)
          NL2=IPC(ILOFF,MODE2)
          NL3=IPC(ILOFF,MODE3)
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
          XYZ=XA3(I)
          ZYX=0
          IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
            DO I=1,NP3(IND)
              K=JP3(I,IND,1)+1
              L=JP3(I,IND,2)+1
              N=JP3(I,IND,3)+1
              ZYX=ZYX+CP3(I,IND)*XKAN(NL1,NR1,K,MOD1)*
     1        XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)
            END DO
C**ANALYTIC
          END IF
3000      CONTINUE
          XK(ILHS,IRHS)=(XYZ+ZYX)*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
C**TEMPORARY-LAN
        J0=IROFF*(IROFF-1)/2+IOFF-J0R
C**TEMPORARY-LAN
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IP(IROFF,KCONT)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IP(IROFF,MCONT).NE.IP(ILOFF,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IP(ILOFF,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
99999 CONTINUE

      IF(ITIM3A.EQ.1)THEN
        CALL TIMIT(3)
        ITIM3A=2
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI3B(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,IP,ISIZMX,IPC1,
     2ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP3,
     3ISIZE3,XA3,XRA3,NSIZE,IRL,XK2,XK1,KTEMP,XCON2,NVAL1,NVAL2,
     4ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,Z0,VP,VPR,VC,VCR,VR,
     5VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,
     6IND,NONZC1,NONZC2,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,JP3,NTOT3,MAX3,
     7INDK,INDL,MAXVAL,MAXQ3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,10),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,10),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA3(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP3(ISIZE3,3),IRL(MAXQ3,MAXQ3,MAXQ3)
      DIMENSION TEMP1(I10,JCI3,JCI3),TEMP2(I10,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(I10,JCI1,JCI1,MM1),Y0(I10,JCI2,JCI2,MM2)
      DIMENSION Z0(I10,JCI3,JCI3,MM3)
      DIMENSION X(10),Y(10),C(10)
C*****************************************************************
      DIMENSION XA(1),XK1(KTEMP,1),XA3(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C*****************************************************************
C**ANALYTIC
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),JP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI3B(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,
     2MM3/MD3,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,
     3KCONT2,XA,ISIZE,IP3,ISIZE3,XA3,XRA3,NSIZE,IRL,XK2,XK1,KTEMP,
     4XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,
     5Z0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,
     6NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I10,XKAN,MAXQU,MAXPOW,NP3,
     7CP3,JP3,NTOT3,MAX3,INDK,INDL,MAXVAL,MAXQ3)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI3B(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,IP,
     2ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,
     3ISIZE,IP3,ISIZE3,XA3,XRA3,NSIZE,IRL,XK2,XK1,KTEMP,XCON2,NVAL1,
     4NVAL2,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,Z0,VP,VPR,VC,
     5VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,
     6NVSYMX,IND,NONZC1,NONZC2,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,JP3,NTOT3,
     7MAX3,INDK,INDL,MAXVAL,MAXQ3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,10),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,10),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA3(1)

C**TEMPORARY
      DIMENSION IRL(MAXQ3,MAXQ3,MAXQ3),NX3(3),MX3(3)
C**TEMPORARY

      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP3(ISIZE3,3)
      DIMENSION TEMP1(I10,JCI3,JCI3),TEMP2(I10,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(I10,JCI1,JCI1,MM1),Y0(I10,JCI2,JCI2,MM2)
      DIMENSION Z0(I10,JCI3,JCI3,MM3)
      DIMENSION X(10),Y(10),C(10)
C*****************************************************************
      DIMENSION XA(1),XK1(KTEMP,1),XA3(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C*****************************************************************
C**ANALYTIC
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),JP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM3A.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCI3B'
          CALL TIMIT(1)
          CALL FLUSH(IOUT)
        END IF
      END IF

C**************************************************************
C**************************************************************
C**ONE CONTRACTION SCHEME HAS A SINGLE MODE
C**ONE CONTRACTION SCHEME HAS TWO MODES
C**CHECK SCHEME '1'
      I1=KCONT1
      NC1=0
      DO NN=1,ICONT(I1)
        IF(MOD1.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE1
          JC1(NC1)=1
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD2.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE2
          JC1(NC1)=2
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD3.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE3
          JC1(NC1)=3
        END IF
      END DO
C**CHECK SCHEME '2'
      I2=KCONT2
      NC2=0
      DO NN=1,ICONT(I2)
        IF(MOD1.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE1
          JC2(NC2)=1
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD2.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE2
          JC2(NC2)=2
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD3.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE3
          JC2(NC2)=3
        END IF
      END DO
C**NKMOD1 IS TOTAL NUMBER OF MODES FOR SCHEME INVOLVING 'MOD1' (K)
C**NKMOD2 IS TOTAL NUMBER OF MODES FOR THE OTHER SCHEME (ONE OR BOTH 
C**OF 'MOD2' (L) AND 'MOD3' (N))
C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L) OR 
C**'MOD3' (N)
C**IF THEY ARE IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
C**************************************************************
      IF(IND.NE.0)GO TO 7777
C**************************************************************

C**FIND MAX QUANTUM
      MAXQ3=0
      DO IR=1,NSIZE
        NR1=IP3(IR,1)
        NR2=IP3(IR,2)
        NR3=IP3(IR,3)
        IF(NR1.GT.MAXQ3)MAXQ3=NR1
        IF(NR2.GT.MAXQ3)MAXQ3=NR2
        IF(NR3.GT.MAXQ3)MAXQ3=NR3
      END DO

C**ANALYTIC
      IND3=INDL(MOD1)+INDK(MOD2)+MOD3
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI3
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,3)
      IFACTL=INTFAC(NMODE,ICOUPL,3)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(73)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(73)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.2)READ(63)VR
        IF(ICOUPC.GT.2)READ(83)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.2)READ(63)VRR
        IF(ICOUPC.GT.2)READ(83)VCR
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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(2)
          X(8)=X(2)
          X(9)=X(2)
          X(10)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(1)
            Y(5)=Y(3)
            Y(6)=Y(3)
            Y(7)=Y(3)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(3)
            DO K=1,I10
              X0(K,NL1,NR1,M1)=Y(11-K)*X(11-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(1)
          X(10)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(3)
            Y(7)=Y(3)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            DO K=1,I10
              Y0(K,NL2,NR2,M2)=Y(11-K)*X(11-K)
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
          X(4)=H3(NR3,M3,2)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(4)
          X(8)=X(1)
          X(9)=X(4)
          X(10)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=H3(NL3,M3,2)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(5)
            Y(9)=Y(5)
            Y(10)=Y(1)
            DO K=1,I10
              Z0(K,NL3,NR3,M3)=Y(11-K)*X(11-K)
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
                DO K=1,I10
                  TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=0
                END DO
              END DO
            END DO
          END DO
        END DO
C**START 2-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO IRHS3=1,JCI3
            DO ILHS3=1,JCI3
              DO K=1,I10
                TEMP1(K,ILHS3,IRHS3)=0
              END DO
            END DO
          END DO
C**START 1-MODE INTEGRATION
CCCC      DO M3=1,MM3/MD3
          DO M3=1,MM3
            DO I=1,10
              C(I)=0
            END DO
            IF(JCOUPL.GT.0)THEN
              IF(ICOUPC.GT.2)THEN
                C(1)=VC(M3,M2,M1,4)*IFACTC
                C(2)=VC(M3,M2,M1,5)*IFACTC
                C(3)=VC(M3,M2,M1,5)*IFACTC
                C(4)=VC(M3,M2,M1,6)*IFACTC
                C(5)=VC(M3,M2,M1,6)*IFACTC
                C(6)=VC(M3,M2,M1,7)*IFACTC
                C(7)=VC(M3,M2,M1,8)*IFACTC
                C(8)=VC(M3,M2,M1,8)*IFACTC
                C(9)=VC(M3,M2,M1,9)*IFACTC
                C(10)=VC(M3,M2,M1,10)*IFACTC
                IF(J21.GT.1)C(10)=C(10)+VR(KROT,M3,M2,M1)*IFACTC
              END IF
              IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(10)=C(10)+
     1        VP(M3,M2,M1)*IFACTL
            ELSE
              IF(ICOUPC.GT.2)THEN
                C(1)=VCR(M3,M2,M1,4)*IFACTC
                C(2)=VCR(M3,M2,M1,5)*IFACTC
                C(3)=VCR(M3,M2,M1,5)*IFACTC
                C(4)=VCR(M3,M2,M1,6)*IFACTC
                C(5)=VCR(M3,M2,M1,6)*IFACTC
                C(6)=VCR(M3,M2,M1,7)*IFACTC
                C(7)=VCR(M3,M2,M1,8)*IFACTC
                C(8)=VCR(M3,M2,M1,8)*IFACTC
                C(9)=VCR(M3,M2,M1,9)*IFACTC
                C(10)=VCR(M3,M2,M1,10)*IFACTC
                IF(J21.GT.1)C(10)=C(10)+VRR(KROT,M3,M2,M1)*IFACTC
              END IF
              IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(10)=C(10)+
     1        VPR(M3,M2,M1)*IFACTL
            END IF
            DO IRHS3=1,JCI3
              DO ILHS3=1,JCI3
                DO K=1,I10
                  TEMP1(K,ILHS3,IRHS3)=TEMP1(K,ILHS3,IRHS3)+
     1            Z0(K,ILHS3,IRHS3,M3)*C(11-K)
                END DO
              END DO
            END DO
          END DO
C**END 1-MODE INTEGRATION
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO ILHS3=1,JCI3
                DO ILHS2=1,JCI2
                  DO K=1,I10
                    TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1              TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2              Y0(K,ILHS2,IRHS2,M2)*TEMP1(K,ILHS3,IRHS3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
C**END 2-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
        DO IRHS=1,NSIZE
          NR1=IP3(IRHS,1)
          NR2=IP3(IRHS,2)
          NR3=IP3(IRHS,3)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP3(ILHS,1)
            NL2=IP3(ILHS,2)
            NL3=IP3(ILHS,3)
            DO K=1,I10
              XA3(ILHS+J0)=XA3(ILHS+J0)+
     1        TEMP2(K,NL2,NL3,NR2,NR3)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 3-MODE INTEGRATION
      CALL MATOUT(XA3,XRA3,NSIZE,23)
      GO TO 7777
6666  CALL MATIN(XA3,XRA3,NSIZE,23)
7777  CONTINUE

      IF(IND.NE.0)THEN
        DO IR=1,NSIZE
          NR=0
          DO J=1,NC1
            J1=JC1(J)
            NR1=IP3(IR,J1)
            NR=NR+1
            NX3(NR)=NR1
          END DO
          DO J=1,NC2
            J2=JC2(J)
            NR2=IP3(IR,J2)
            NR=NR+1
            NX3(NR)=NR2
          END DO
          IRL(NX3(1),NX3(2),NX3(3))=IR
        END DO
      END IF

C***********************************ALGORITHM FROM VCI3

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFFR=IOFFR+NCSIZE(I)
      END DO
      J0R=IOFFR*(IOFFR+1)/2
      CALL FLUSH(IOUT)
C**TEMPORARY-LAN
      JSR=1
C**TEMPORARY-LAN
C     DO 9999 ISMR1=1,NVSYM
      DO 9999 ISMR1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE

        IF(NONC1.GT.NONC2)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          KCCCL=KCSZL
CC        IF(NONC1.GT.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             KCCCL=KCSZL
C           ELSE
C             KCCCL=IRH1
C           END IF
CC        END IF
          DO ILH1=1,KCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=IROFF
            NONZC1(2,NONZ1)=ILOFF
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
        IF(NONZ1.EQ.0)GO TO 5556
        NONZC1(3,NONZ1+1)=0
        END IF

      END IF
 
CCCC  IF(NONZ1.EQ.0)GO TO 5556
CCCC  NONZC1(3,NONZ1+1)=0

4444  CONTINUE

C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE
   
        IF(NONC1.LE.NONC2)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          LCCCL=LCSZL
CC        IF(NONC1.LE.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             LCCCL=LCSZL
C           ELSE
C             LCCCL=IRH2
C           END IF
CC        END IF
          DO ILH2=1,LCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=JROFF
            NONZC2(2,NONZ2)=JLOFF
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
        IF(NONZ2.EQ.0)GO TO 5555
        NONZC2(3,NONZ2+1)=0
        END IF

      END IF

CCCC  IF(NONZ2.EQ.0)GO TO 5555
CCCC  NONZC2(3,NONZ2+1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM3A.EQ.0)THEN
C       WRITE(IOUT,*)'VCCI3B - ISMR1,ISML1 = ',ISMR1,ISML1
C       WRITE(IOUT,*)'ICONDP = ',ICONDP
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C       WRITE(IOUT,*)'NONC1,NONC2 = ',NONC1,NONC2
C     END IF

      IF(NONC1.LE.NONC2)GO TO 3003
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '2' (L AND/OR N)

      DO IRH2=1,LCSZR
        JROFF=IRH2+LCOFFR
        DO I=1,NKVALR
          DO J=1,NKVALL
            DO K=1,NLVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          LCCCL=LCSZL
        ELSE
          LCCCL=IRH2
        END IF
        DO ILH2=1,LCCCL
          JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
          DO K=1,NKMOD2
            IGOT=0
            DO I=1,NC2
              IF(K.EQ.MC2(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 6001
          END DO
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH2.EQ.IRH2)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NKVALR
            DO J=1,NKVALL
              XCON2(J,I)=0
            END DO
          END DO
C
C**CONTRACTION SCHEME  '1' (K)
C
C         DO IRH1=1,KCSZR
C           IROFF=IRH1+KCOFFR
C           DO J=1,NKVALL
C             XK1(J,1)=0
C           END DO
C           DO ILH1=1,KCSZL
C             ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
C             DO K=1,NKMOD1
C               IGOT=0
C               DO I=1,NC1
C                 IF(K.EQ.MC1(I))IGOT=1
C               END DO
C               IF(IGOT.EQ.0.AND.
C    1          (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 5001
C             END DO

          IPREV=0
          DO INH1=1,NONZ1
              IRH1=NONZC1(3,INH1)
              INEXT=NONZC1(3,INH1+1)
              ILH1=NONZC1(4,INH1)
              IROFF=NONZC1(1,INH1)
              ILOFF=NONZC1(2,INH1)
              IF(IPREV.NE.IRH1)THEN
C**NEW RHS
                DO J=1,NKVALL
                  XK1(J,1)=0
                END DO
              END IF

              DO I=1,NC1
                NCR1(I)=IPC1(IROFF,MC1(I))
                NCL1(I)=IPC1(ILOFF,MC1(I))
              END DO

C**FIND RHS INDEX
C             DO IR=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCR1(I).NE.IP3(IR,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCR2(I).NE.IP3(IR,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 1000
C             END DO
1000          CONTINUE
C**FIND LHS INDEX
C             DO IL=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCL1(I).NE.IP3(IL,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCL2(I).NE.IP3(IL,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 2000
C             END DO
2000          CONTINUE

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX3(NR)=NR1
        MX3(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX3(NR)=NR2
        MX3(NR)=NL2
      END DO
      IR=IRL(NX3(1),NX3(2),NX3(3))
      IL=IRL(MX3(1),MX3(2),MX3(3))
C**TEMPORARY

C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
C             XYZ=XA3(I)
C             ZYX=0
C             IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
C               NR1=IP3(MR,1)
C               NL1=IP3(ML,1)
C               NR2=IP3(MR,2)
C               NL2=IP3(ML,2)
C               NR3=IP3(MR,3)
C               NL3=IP3(ML,3)
C               DO I=1,NP3(IND3)
C                 K=JP3(I,IND3,1)+1
C                 L=JP3(I,IND3,2)+1
C                 N=JP3(I,IND3,3)+1
C                 ZYX=ZYX+CP3(I,IND3)*XKAN(NL1,NR1,K,MOD1)*
C    1            XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)
C               END DO
C**ANALYTIC
C             END IF
C             XXX=XYZ+ZYX
              XXX=XA3(I)

              DO K=1,NKVALL
                XK1(K,1)=XK1(K,1)+XXX*CFS1(ILH1,K,KEL,KSML)
              END DO
5001          CONTINUE

C           END DO
            IF(IRH1.NE.INEXT)THEN

            DO K=1,NKVALR
              DO L=1,NKVALL
                XCON2(L,K)=XCON2(L,K)+XK1(L,1)*CFS1(IRH1,K,KEL,KSMR)
              END DO
            END DO

            END IF
            IPREV=IRH1
          END DO
CC
          DO I=1,NKVALR
            DO J=1,NKVALL
              DO K=1,NLVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6001      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILL,IKR,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKR,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILR,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILR,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO
C*************************************

      GO TO 5555

3003  CONTINUE

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '1' (K)

      DO IRH1=1,KCSZR
        IROFF=IRH1+KCOFFR
        DO I=1,NLVALR
          DO J=1,NLVALL
            DO K=1,NKVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          KCCCL=KCSZL
        ELSE
          KCCCL=IRH1
        END IF
        DO ILH1=1,KCCCL
          ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
          DO K=1,NKMOD1
            IGOT=0
            DO I=1,NC1
              IF(K.EQ.MC1(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 6002
          END DO
          DO I=1,NC1
            NCR1(I)=IPC1(IROFF,MC1(I))
            NCL1(I)=IPC1(ILOFF,MC1(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH1.EQ.IRH1)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NLVALR
            DO J=1,NLVALL
              XCON2(J,I)=0
            END DO
          END DO

C**CONTRACTION SCHEME  '2' (L AND/OR N)

C         DO IRH2=1,LCSZR
C           JROFF=IRH2+LCOFFR
C           DO J=1,NLVALL
C             XK1(J,1)=0
C           END DO
C           DO ILH2=1,LCSZL
C             JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
C             DO K=1,NKMOD2
C               IGOT=0
C               DO I=1,NC2
C                 IF(K.EQ.MC2(I))IGOT=1
C               END DO
C               IF(IGOT.EQ.0.AND.
C    1          (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 5002
C             END DO

          JPREV=0
          DO INH2=1,NONZ2
              IRH2=NONZC2(3,INH2)
              JNEXT=NONZC2(3,INH2+1)
              ILH2=NONZC2(4,INH2)
              JROFF=NONZC2(1,INH2)
              JLOFF=NONZC2(2,INH2)
              IF(JPREV.NE.IRH2)THEN
C**NEW RHS
                DO J=1,NLVALL
                  XK1(J,1)=0
                END DO
              END IF

              DO I=1,NC2
                NCR2(I)=IPC2(JROFF,MC2(I))
                NCL2(I)=IPC2(JLOFF,MC2(I))
              END DO

C**FIND RHS INDEX
C             DO IR=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCR1(I).NE.IP3(IR,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCR2(I).NE.IP3(IR,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 1001
C             END DO
1001          CONTINUE
C**FIND LHS INDEX
C             DO IL=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCL1(I).NE.IP3(IL,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCL2(I).NE.IP3(IL,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 2001
C             END DO
2001          CONTINUE

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX3(NR)=NR1
        MX3(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX3(NR)=NR2
        MX3(NR)=NL2
      END DO
      IR=IRL(NX3(1),NX3(2),NX3(3))
      IL=IRL(MX3(1),MX3(2),MX3(3))
C**TEMPORARY

C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
C             XYZ=XA3(I)
C             ZYX=0
C             IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
C               NR1=IP3(MR,1)
C               NL1=IP3(ML,1)
C               NR2=IP3(MR,2)
C               NL2=IP3(ML,2)
C               NR3=IP3(MR,3)
C               NL3=IP3(ML,3)
C               DO I=1,NP3(IND3)
C                 K=JP3(I,IND3,1)+1
C                 L=JP3(I,IND3,2)+1
C                 N=JP3(I,IND3,3)+1
C                 ZYX=ZYX+CP3(I,IND3)*XKAN(NL1,NR1,K,MOD1)*
C    1            XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)
C               END DO
C**ANALYTIC
C             END IF
C             XK1(ILH2,IRH2)=XYZ+ZYX
C             XXX=XYZ+ZYX
              XXX=XA3(I)
 
              DO J=1,NLVALL
                XK1(J,1)=XK1(J,1)+XXX*CFS2(ILH2,J,KEL,LSML)
              END DO
5002          CONTINUE

C           END DO
            IF(IRH2.NE.JNEXT)THEN

            DO I=1,NLVALR
              DO J=1,NLVALL
                XCON2(J,I)=XCON2(J,I)+XK1(J,1)*CFS2(IRH2,I,KEL,LSMR)
              END DO
            END DO

            END IF
            JPREV=IRH2
          END DO
CC
          DO I=1,NLVALR
            DO J=1,NLVALL
              DO K=1,NKVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6002      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKR,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKL,ILR,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILR,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO

5555  CONTINUE
5556  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM3A.EQ.0)THEN
        IF(IND.NE.0)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM3A=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI3C(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,IP,ISIZMX,IPC,
     2ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP3,ISIZE3,XA3,XRA3,NSIZE,XK,TEMP,
     3IRL,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,Z0,VP,
     4VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,
     5NVSYMX,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,JP3,NTOT3,MAX3,INDK,INDL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,10),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,10),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA3(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION IP(ISIZMX,NNMODE),IRL(ICSIZE)
      DIMENSION IPC(IPSIZE,1),IP3(ISIZE3,3),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I10,JCI3,JCI3),TEMP2(I10,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(I10,JCI1,JCI1,MM1),Y0(I10,JCI2,JCI2,MM2)
      DIMENSION Z0(I10,JCI3,JCI3,MM3)
      DIMENSION X(10),Y(10),C(10)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA3(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),JP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
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
      CALL VDCI3C(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,
     2MM3/MD3,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP3,ISIZE3,
     3XA3,XRA3,NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,
     4JCI1,JCI2,JCI3,X0,Y0,Z0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,
     5NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,
     6JP3,NTOT3,MAX3,INDK,INDL)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI3C(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,IP,
     2ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP3,ISIZE3,XA3,XRA3,
     3NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,
     4JCI3,X0,Y0,Z0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,
     5ISIZXX,NVALX,KEL21,NVSYMX,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,JP3,
     6NTOT3,MAX3,INDK,INDL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,10),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,10),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA3(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION IP(ISIZMX,NNMODE),IRL(ICSIZE)
      DIMENSION IPC(IPSIZE,1),IP3(ISIZE3,3),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I10,JCI3,JCI3),TEMP2(I10,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(I10,JCI1,JCI1,MM1),Y0(I10,JCI2,JCI2,MM2)
      DIMENSION Z0(I10,JCI3,JCI3,MM3)
      DIMENSION X(10),Y(10),C(10)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA3(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),JP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP

      IF(ITIM3A.EQ.1)THEN
        WRITE(IOUT,*)'Calculating VCCI3C'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1', 'MOD2'
C**AND 'MOD3'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      NLMODE=ICONT(2)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        NLMODE=ICONT(1)
        MCONT=1
      END IF

C**ANALYTIC
      IND=INDL(MOD1)+INDK(MOD2)+MOD3
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI3
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,3)
      IFACTL=INTFAC(NMODE,ICOUPL,3)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(73)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(73)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.2)READ(63)VR
        IF(ICOUPC.GT.2)READ(83)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.2)READ(63)VRR
        IF(ICOUPC.GT.2)READ(83)VCR
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

C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(2)
          X(8)=X(2)
          X(9)=X(2)
          X(10)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(1)
            Y(5)=Y(3)
            Y(6)=Y(3)
            Y(7)=Y(3)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(3)
            DO K=1,I10
              X0(K,NL1,NR1,M1)=Y(11-K)*X(11-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(1)
          X(10)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(3)
            Y(7)=Y(3)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            DO K=1,I10
              Y0(K,NL2,NR2,M2)=Y(11-K)*X(11-K)
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
          X(4)=H3(NR3,M3,2)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(4)
          X(8)=X(1)
          X(9)=X(4)
          X(10)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=H3(NL3,M3,2)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(5)
            Y(9)=Y(5)
            Y(10)=Y(1)
            DO K=1,I10
              Z0(K,NL3,NR3,M3)=Y(11-K)*X(11-K)
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
                DO K=1,I10
                  TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=0.D0
                END DO
              END DO
            END DO
          END DO
        END DO
C**START 2-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO IRHS3=1,JCI3
            DO ILHS3=1,JCI3
              DO K=1,I10
                TEMP1(K,ILHS3,IRHS3)=0.D0
              END DO
            END DO
          END DO
C**START 1-MODE INTEGRATION
CCCC      DO M3=1,MM3/MD3
          DO M3=1,MM3
            DO I=1,10
              C(I)=0
            END DO
            IF(JCOUPL.GT.0)THEN
              IF(ICOUPC.GT.2)THEN
                C(1)=VC(M3,M2,M1,4)*IFACTC
                C(2)=VC(M3,M2,M1,5)*IFACTC
                C(3)=VC(M3,M2,M1,5)*IFACTC
                C(4)=VC(M3,M2,M1,6)*IFACTC
                C(5)=VC(M3,M2,M1,6)*IFACTC
                C(6)=VC(M3,M2,M1,7)*IFACTC
                C(7)=VC(M3,M2,M1,8)*IFACTC
                C(8)=VC(M3,M2,M1,8)*IFACTC
                C(9)=VC(M3,M2,M1,9)*IFACTC
                C(10)=VC(M3,M2,M1,10)*IFACTC
                IF(J21.GT.1)C(10)=C(10)+VR(KROT,M3,M2,M1)*IFACTC
              END IF
              IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(10)=C(10)+
     1        VP(M3,M2,M1)*IFACTL
            ELSE
              IF(ICOUPC.GT.2)THEN
                C(1)=VCR(M3,M2,M1,4)*IFACTC
                C(2)=VCR(M3,M2,M1,5)*IFACTC
                C(3)=VCR(M3,M2,M1,5)*IFACTC
                C(4)=VCR(M3,M2,M1,6)*IFACTC
                C(5)=VCR(M3,M2,M1,6)*IFACTC
                C(6)=VCR(M3,M2,M1,7)*IFACTC
                C(7)=VCR(M3,M2,M1,8)*IFACTC
                C(8)=VCR(M3,M2,M1,8)*IFACTC
                C(9)=VCR(M3,M2,M1,9)*IFACTC
                C(10)=VCR(M3,M2,M1,10)*IFACTC
                IF(J21.GT.1)C(10)=C(10)+VRR(KROT,M3,M2,M1)*IFACTC
              END IF
              IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(10)=C(10)+
     1        VPR(M3,M2,M1)*IFACTL
            END IF
            DO IRHS3=1,JCI3
              DO ILHS3=1,JCI3
                DO K=1,I10
                  TEMP1(K,ILHS3,IRHS3)=TEMP1(K,ILHS3,IRHS3)+
     1            Z0(K,ILHS3,IRHS3,M3)*C(11-K)
                END DO
              END DO
            END DO
          END DO
C**END 1-MODE INTEGRATION
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO ILHS3=1,JCI3
                DO ILHS2=1,JCI2
                  DO K=1,I10
                    TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1              TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2              Y0(K,ILHS2,IRHS2,M2)*TEMP1(K,ILHS3,IRHS3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
C**END 2-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
        DO IRHS=1,NSIZE
          NR1=IP3(IRHS,1)
          NR2=IP3(IRHS,2)
          NR3=IP3(IRHS,3)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP3(ILHS,1)
            NL2=IP3(ILHS,2)
            NL3=IP3(ILHS,3)
            DO K=1,I10
              XA3(ILHS+J0)=XA3(ILHS+J0)+
     1        TEMP2(K,NL2,NL3,NR2,NR3)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 3-MODE INTEGRATION
      CALL MATOUT(XA3,XRA3,NSIZE,23)
      GO TO 7777
6666  CALL MATIN(XA3,XRA3,NSIZE,23)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI3

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        KCOFF=ICOFF1
        LCOFF=ICOFF2
        KCSZ=ICSZ1
        LCSZ=ICSZ2
        ICNTZ=ICCONT(1)
        NKVAL=NCVAL(1,ISM1)
        NLVAL=NCVAL(2,ISM2)
      ELSE
        ISM=ISM2
        KCOFF=ICOFF2
        LCOFF=ICOFF1
        KCSZ=ICSZ2
        LCSZ=ICSZ1
        ICNTZ=ICCONT(2)
        NKVAL=NCVAL(2,ISM2)
        NLVAL=NCVAL(1,ISM1)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

C
      DO IKRL=1,KCSZ
        IOFFRL=IKRL+KCOFF
        NR1=IPC(IOFFRL,MODE1)
        NR2=IPC(IOFFRL,MODE2)
        NR3=IPC(IOFFRL,MODE3)
        DO IR=1,NSIZE
          IF(NR1.EQ.IP3(IR,1).AND.NR2.EQ.IP3(IR,2).AND.
     1    NR3.EQ.IP3(IR,3))GO TO 1000
        END DO
1000    IRL(IKRL)=IR
      END DO
C
C**CASE 1..INTEGRATION MODES IN CONTRACTED SCHEME.
C**      ..OVERLAP OF NON-CONTRACTED BASIS.
      IF(ICNTZ.GT.0)THEN
        DO IRHS=1,KCSZ
          IOFFR=IRHS+KCOFF
          IR=IRL(IRHS)
          DO ILHS=1,IRHS
            IOFFL=ILHS+KCOFF
C**OVERLAP OF REMAINING STATES
            IS=1
            DO K=1,NKMODE
              IF(IS.EQ.0)GO TO 2000
              IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.(
     1        IPC(IOFFR,K).NE.IPC(IOFFL,K)))IS=0
            END DO
C**OVERLAP OF REMAINING STATES
            IL=IRL(ILHS)
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA3(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              DO I=1,NP3(IND)
                K=JP3(I,IND,1)+1
                L=JP3(I,IND,2)+1
                N=JP3(I,IND,3)+1
                ZYX=ZYX+CP3(I,IND)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)
              END DO
C**ANALYTIC
            END IF
2000        CONTINUE
            XK(ILHS,IRHS)=(XYZ+ZYX)*IS
            XK(IRHS,ILHS)=XK(ILHS,IRHS)
          END DO
        END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',KCSZ,NKVAL,KCSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
        CALL DGEMM('T','N',NKVAL,NKVAL,KCSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

        IRHS=0
        DO IKR=1,NKVAL
          DO ILR=1,LCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            DO IKL=1,IKR
              ILHS=(IKL-1)*LCSZ+ILR
              XYZ=XCON(IKL,IKR)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ
            END DO
          END DO
        END DO

C**CASE 2..INTEGRATION MODES IN NON-CONTRACTED SCHEME.
C**      ..OVERLAP OF CONTRACTED BASIS.
      ELSE
        IRHS=0
        DO ILR=1,NLVAL
          DO IKR=1,KCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            IOFFR=IKR+KCOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            IR=IRL(IKR)
            DO IKL=1,IKR
              ILHS=(ILR-1)*KCSZ+IKL
              IOFFL=IKL+KCOFF
C**OVERLAP OF REMAINING STATES
              IS=1
              DO K=1,NKMODE
                IF(IS.EQ.0)GO TO 3000
                IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.(
     1          IPC(IOFFR,K).NE.IPC(IOFFL,K)))IS=0
              END DO
C**OVERLAP OF REMAINING STATES
              IL=IRL(IKL)
C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
              XYZ=XA3(I)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ*IS
3000          CONTINUE
            END DO
          END DO
        END DO
      END IF

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
      IF(ITIM3A.EQ.1)THEN
        CALL TIMIT(3)
        ITIM3A=2
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI3D(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,IP,ISIZMX,IPC1,
     2ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP3,
     3ISIZE3,XA3,XRA3,NSIZE,XK1,XK2,TEMP,KTEMP,XCON2,NVAL1,NVAL2,
     4ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,Z0,VP,VPR,VC,VCR,VR,
     5VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,
     6IND,NONZC1,NONZC2,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,JP3,NTOT3,MAX3,
     7INDK,INDL,MAXVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,10),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,10),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA3(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP3(ISIZE3,3),TEMP(KTEMP,1)
      DIMENSION TEMP1(I10,JCI3,JCI3),TEMP2(I10,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(I10,JCI1,JCI1,MM1),Y0(I10,JCI2,JCI2,MM2)
      DIMENSION Z0(I10,JCI3,JCI3,MM3)
      DIMENSION X(10),Y(10),C(10)
C**************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA3(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,1)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C**************************************************************
C**ANALYTIC
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),JP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI3D(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,
     2MM3/MD3,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,
     3KCONT2,XA,ISIZE,IP3,ISIZE3,XA3,XRA3,NSIZE,XK1,XK2,TEMP,KTEMP,
     4XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,
     5Z0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,
     6NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I10,XKAN,MAXQU,MAXPOW,NP3,
     7CP3,JP3,NTOT3,MAX3,INDK,INDL,MAXVAL)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI3D(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1H1,XQ1,H2,XQ2,H3,XQ3,NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,IP,
     2ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,
     3ISIZE,IP3,ISIZE3,XA3,XRA3,NSIZE,XK1,XK2,TEMP,KTEMP,XCON2,
     4NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,JCI1,JCI2,JCI3,X0,Y0,Z0,VP,
     5VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,
     6KEL21,NVSYMX,IND,NONZC1,NONZC2,I10,XKAN,MAXQU,MAXPOW,NP3,CP3,JP3,
     7NTOT3,MAX3,INDK,INDL,MAXVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,10),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,10),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA3(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP3(ISIZE3,3),TEMP(KTEMP,1)
      DIMENSION TEMP1(I10,JCI3,JCI3),TEMP2(I10,JCI2,JCI3,JCI2,JCI3)
      DIMENSION X0(I10,JCI1,JCI1,MM1),Y0(I10,JCI2,JCI2,MM2)
      DIMENSION Z0(I10,JCI3,JCI3,MM3)
      DIMENSION X(10),Y(10),C(10)
C**************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA3(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,1)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C**************************************************************
C**ANALYTIC
      DIMENSION NP3(NTOT3),CP3(MAX3,NTOT3),JP3(MAX3,NTOT3,3)
      DIMENSION INDK(1),INDL(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM3A.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCI3D'
          CALL TIMIT(1)
          CALL FLUSH(IOUT)
        END IF
      END IF

C**************************************************************
C**************************************************************
C**ONE CONTRACTION SCHEME HAS A SINGLE MODE
C**ONE CONTRACTION SCHEME HAS TWO MODES
C**CHECK SCHEME '1'
      I1=KCONT1
      NC1=0
      DO NN=1,ICONT(I1)
        IF(MOD1.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE1
          JC1(NC1)=1
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD2.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE2
          JC1(NC1)=2
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD3.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE3
          JC1(NC1)=3
        END IF
      END DO
C**CHECK SCHEME '2'
      I2=KCONT2
      NC2=0
      DO NN=1,ICONT(I2)
        IF(MOD1.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE1
          JC2(NC2)=1
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD2.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE2
          JC2(NC2)=2
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD3.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE3
          JC2(NC2)=3
        END IF
      END DO
C**NKMOD1 IS TOTAL NUMBER OF MODES FOR SCHEME INVOLVING 'MOD1' (K)
C**NKMOD2 IS TOTAL NUMBER OF MODES FOR THE OTHER SCHEME (ONE OR BOTH 
C**OF 'MOD2' (L) AND 'MOD3' (N))
C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L) OR 
C**'MOD3' (N)
C**IF THEY ARE IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
C**************************************************************
      IF(IND.NE.0)GO TO 7777
C**************************************************************

C**ANALYTIC
      IND3=INDL(MOD1)+INDK(MOD2)+MOD3
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI3
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,3)
      IFACTL=INTFAC(NMODE,ICOUPL,3)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(73)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(73)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.2)READ(63)VR
        IF(ICOUPC.GT.2)READ(83)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.2)READ(63)VRR
        IF(ICOUPC.GT.2)READ(83)VCR
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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(2)
          X(8)=X(2)
          X(9)=X(2)
          X(10)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(1)
            Y(5)=Y(3)
            Y(6)=Y(3)
            Y(7)=Y(3)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(3)
            DO K=1,I10
              X0(K,NL1,NR1,M1)=Y(11-K)*X(11-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(1)
          X(10)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(3)
            Y(7)=Y(3)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            DO K=1,I10
              Y0(K,NL2,NR2,M2)=Y(11-K)*X(11-K)
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
          X(4)=H3(NR3,M3,2)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(4)
          X(8)=X(1)
          X(9)=X(4)
          X(10)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=H3(NL3,M3,2)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(5)
            Y(9)=Y(5)
            Y(10)=Y(1)
            DO K=1,I10
              Z0(K,NL3,NR3,M3)=Y(11-K)*X(11-K)
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
                DO K=1,I10
                  TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=0
                END DO
              END DO
            END DO
          END DO
        END DO
C**START 2-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO IRHS3=1,JCI3
            DO ILHS3=1,JCI3
              DO K=1,I10
                TEMP1(K,ILHS3,IRHS3)=0
              END DO
            END DO
          END DO
C**START 1-MODE INTEGRATION
CCCC      DO M3=1,MM3/MD3
          DO M3=1,MM3
            IF(JCOUPL.GT.0)THEN
              IF(ICOUPC.GT.2)THEN
                C(1)=VC(M3,M2,M1,4)*IFACTC
                C(2)=VC(M3,M2,M1,5)*IFACTC
                C(3)=VC(M3,M2,M1,5)*IFACTC
                C(4)=VC(M3,M2,M1,6)*IFACTC
                C(5)=VC(M3,M2,M1,6)*IFACTC
                C(6)=VC(M3,M2,M1,7)*IFACTC
                C(7)=VC(M3,M2,M1,8)*IFACTC
                C(8)=VC(M3,M2,M1,8)*IFACTC
                C(9)=VC(M3,M2,M1,9)*IFACTC
                C(10)=VC(M3,M2,M1,10)*IFACTC
                IF(J21.GT.1)C(10)=C(10)+VR(KROT,M3,M2,M1)*IFACTC
              END IF
              IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(10)=C(10)+
     1        VP(M3,M2,M1)*IFACTL
            ELSE
              IF(ICOUPC.GT.2)THEN
                C(1)=VCR(M3,M2,M1,4)*IFACTC
                C(2)=VCR(M3,M2,M1,5)*IFACTC
                C(3)=VCR(M3,M2,M1,5)*IFACTC
                C(4)=VCR(M3,M2,M1,6)*IFACTC
                C(5)=VCR(M3,M2,M1,6)*IFACTC
                C(6)=VCR(M3,M2,M1,7)*IFACTC
                C(7)=VCR(M3,M2,M1,8)*IFACTC
                C(8)=VCR(M3,M2,M1,8)*IFACTC
                C(9)=VCR(M3,M2,M1,9)*IFACTC
                C(10)=VCR(M3,M2,M1,10)*IFACTC
                IF(J21.GT.1)C(10)=C(10)+VRR(KROT,M3,M2,M1)*IFACTC
              END IF
              IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(10)=C(10)+
     1        VPR(M3,M2,M1)*IFACTL
            END IF
            DO IRHS3=1,JCI3
              DO ILHS3=1,JCI3
                DO K=1,I10
                  TEMP1(K,ILHS3,IRHS3)=TEMP1(K,ILHS3,IRHS3)+
     1            Z0(K,ILHS3,IRHS3,M3)*C(11-K)
                END DO
              END DO
            END DO
          END DO
C**END 1-MODE INTEGRATION
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO ILHS3=1,JCI3
                DO ILHS2=1,JCI2
                  DO K=1,I10
                    TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1              TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2              Y0(K,ILHS2,IRHS2,M2)*TEMP1(K,ILHS3,IRHS3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
C**END 2-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
        DO IRHS=1,NSIZE
          NR1=IP3(IRHS,1)
          NR2=IP3(IRHS,2)
          NR3=IP3(IRHS,3)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP3(ILHS,1)
            NL2=IP3(ILHS,2)
            NL3=IP3(ILHS,3)
            DO K=1,I10
              XA3(ILHS+J0)=XA3(ILHS+J0)+
     1        TEMP2(K,NL2,NL3,NR2,NR3)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 3-MODE INTEGRATION
      CALL MATOUT(XA3,XRA3,NSIZE,23)
      GO TO 7777
6666  CALL MATIN(XA3,XRA3,NSIZE,23)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI3

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
      JSR=1
      DO 9999 ISMR1=1,NVSYM
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=IROFF
            NONZC1(2,NONZ1)=ILOFF
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
      END IF
      IF(NONZ1.EQ.0)GO TO 5555
      NONZC1(3,NONZ1+1)=0

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=JROFF
            NONZC2(2,NONZ2)=JLOFF
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
      END IF
      IF(NONZ2.EQ.0)GO TO 5555
      NONZC2(3,NONZ2+1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM3A.EQ.0)THEN
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C     END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CASE 1..K IS IN CONTRACTED SCHEME, L OR N IS IN UNCONTRACTED SCHEME
      IF(ICCONT(1).GT.0)THEN

        IF(NONC1.LE.NONC2)GO TO 3003
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        IPREV=0
        DO INH1=1,NONZ1
          IRH1=NONZC1(3,INH1)
C         IF(INH1.NE.NONZ1)THEN
            INEXT=NONZC1(3,INH1+1)
C         ELSE
C           INEXT=0
C         END IF
          ILH1=NONZC1(4,INH1)
          IF(IPREV.NE.IRH1)THEN
            DO INH2=1,NONZ2
              DO IKL=1,NKVALL
                XK2(IKL,INH2)=0
              END DO
            END DO
          END IF
          IROFF=NONZC1(1,INH1)
          ILOFF=NONZC1(2,INH1)
          DO I=1,NC1
            NCR1(I)=IPC1(IROFF,MC1(I))
            NCL1(I)=IPC1(ILOFF,MC1(I))
          END DO

C**CONTRACTION SCHEME INVOLVING 'MOD2' (L AND/OR N)
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            ILH2=NONZC2(4,INH2)
            JROFF=NONZC2(1,INH2)
            JLOFF=NONZC2(2,INH2)
            DO I=1,NC2
              NCR2(I)=IPC2(JROFF,MC2(I))
              NCL2(I)=IPC2(JLOFF,MC2(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP3(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP3(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1003
            END DO
1003        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP3(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP3(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2003
            END DO
2003        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA3(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP3(MR,1)
              NL1=IP3(ML,1)
              NR2=IP3(MR,2)
              NL2=IP3(ML,2)
              NR3=IP3(MR,3)
              NL3=IP3(ML,3)
              DO I=1,NP3(IND3)
                K=JP3(I,IND3,1)+1
                L=JP3(I,IND3,2)+1
                N=JP3(I,IND3,3)+1
                ZYX=ZYX+CP3(I,IND3)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)
              END DO
C**ANALYTIC
            END IF
            DO IKL=1,NKVALL
              XK2(IKL,INH2)=XK2(IKL,INH2)+(XYZ+ZYX)*
     1        CFS1(ILH1,IKL,KEL,KSML)
            END DO
          END DO
          IF(IRH1.NE.INEXT)THEN
            DO INH2=1,NONZ2
              IRH2=NONZC2(3,INH2)
              ILH2=NONZC2(4,INH2)
              IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
                DO IKR=1,NKVALR
                  CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                  IRHS=(IKR-1)*LCSZR+IRH2
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO IKL=1,NKVALL
                    ILHS=(IKL-1)*LCSZL+ILH2
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(IKL,INH2)
                  END DO
                END DO
C**OFF-DIAGONAL BLOCK - END
              ELSE
C**DIAGONAL BLOCK - START
                DO IKR=1,NKVALR
                  CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                  IRHS=(IKR-1)*LCSZR+IRH2
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO IKL=1,IKR
                    IF(IKL.EQ.IKR.AND.ILH2.GT.IRH2)GO TO 1004
                    ILHS=(IKL-1)*LCSZL+ILH2
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(IKL,INH2)
1004                CONTINUE
                  END DO
                END DO
C**DIAGONAL BLOCK - END
              END IF
            END DO
          END IF
          IPREV=IRH1
        END DO
        GO TO 5555

3003    CONTINUE
        DO IRH1=1,KCSZR
          DO ILH1=1,KCSZL
            XK1(ILH1,IRH1)=0
          END DO
        END DO
 
C**CONTRACTION SCHEME '2' (L AND/OR N)
        DO INH2=1,NONZ2
          IRH2=NONZC2(3,INH2)
          ILH2=NONZC2(4,INH2)
          JROFF=NONZC2(1,INH2)
          JLOFF=NONZC2(2,INH2)
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

C**CONTRACTION SCHEME  '1' (K)
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            ILH1=NONZC1(4,INH1)
            IROFF=NONZC1(1,INH1)
            ILOFF=NONZC1(2,INH1)
            DO I=1,NC1
              NCR1(I)=IPC1(IROFF,MC1(I))
              NCL1(I)=IPC1(ILOFF,MC1(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP3(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP3(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1000
            END DO
1000        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP3(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP3(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2000
            END DO
2000        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA3(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP3(MR,1)
              NL1=IP3(ML,1)
              NR2=IP3(MR,2)
              NL2=IP3(ML,2)
              NR3=IP3(MR,3)
              NL3=IP3(ML,3)
              DO I=1,NP3(IND3)
                K=JP3(I,IND3,1)+1
                L=JP3(I,IND3,2)+1
                N=JP3(I,IND3,3)+1
                ZYX=ZYX+CP3(I,IND3)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)
              END DO
C**ANALYTIC
            END IF
            XK1(ILH1,IRH1)=XYZ+ZYX
          END DO

C************************************DGEMM (RHS)
          CALL DGEMM('N','N',KCSZL,NKVALR,KCSZR,1.0D0,XK1(1,1),
     &    KTEMP,CFS1(1,1,KEL,KSMR),ISIZXX,0.0D0,TEMP(1,1),KTEMP)

C************************************DGEMM (LHS)
          CALL DGEMM('T','N',NKVALL,NKVALR,KCSZL,1.0D0,
     &    CFS1(1,1,KEL,KSML),ISIZXX,TEMP,KTEMP,0.0D0,XCON2(1,1),MAXVAL)

          IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
            DO IKR=1,NKVALR
              IRHS=(IKR-1)*LCSZR+IRH2
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO IKL=1,NKVALL
                ILHS=(IKL-1)*LCSZL+ILH2
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(IKL,IKR)
              END DO
            END DO
C**OFF-DIAGONAL BLOCK - END
          ELSE
C**DIAGONAL BLOCK - START
            DO IKR=1,NKVALR
              IRHS=(IKR-1)*LCSZR+IRH2
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO IKL=1,IKR
                IF(IKL.EQ.IKR.AND.ILH2.GT.IRH2)GO TO 1002
                ILHS=(IKL-1)*LCSZL+ILH2
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(IKL,IKR)
1002            CONTINUE
              END DO
            END DO
C**DIAGONAL BLOCK - END
          END IF
        END DO

C**CASE 2..L OR N IS IN CONTRACTED SCHEME, K IS IN UNCONTRACTED SCHEME
      ELSE

        IF(NONC1.GT.NONC2)GO TO 4004
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L AND/OR N)
        IPREV=0
        DO INH2=1,NONZ2
          IRH2=NONZC2(3,INH2)
C         IF(INH2.NE.NONZ2)THEN
            INEXT=NONZC2(3,INH2+1)
C         ELSE
C           INEXT=0
C         END IF
          ILH2=NONZC2(4,INH2)
          IF(IPREV.NE.IRH2)THEN
            DO INH1=1,NONZ1
              DO ILL=1,NLVALL
                XK2(ILL,INH1)=0
              END DO
            END DO
          END IF
          JROFF=NONZC2(1,INH2)
          JLOFF=NONZC2(2,INH2)
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            ILH1=NONZC1(4,INH1)
            IROFF=NONZC1(1,INH1)
            ILOFF=NONZC1(2,INH1)
            DO I=1,NC1
              NCR1(I)=IPC1(IROFF,MC1(I))
              NCL1(I)=IPC1(ILOFF,MC1(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP3(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP3(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1005
            END DO
1005        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP3(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP3(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2005
            END DO
2005        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA3(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP3(MR,1)
              NL1=IP3(ML,1)
              NR2=IP3(MR,2)
              NL2=IP3(ML,2)
              NR3=IP3(MR,3)
              NL3=IP3(ML,3)
              DO I=1,NP3(IND3)
                K=JP3(I,IND3,1)+1
                L=JP3(I,IND3,2)+1
                N=JP3(I,IND3,3)+1
                ZYX=ZYX+CP3(I,IND3)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)
              END DO
C**ANALYTIC
            END IF
            DO ILL=1,NLVALL
              XK2(ILL,INH1)=XK2(ILL,INH1)+(XYZ+ZYX)*
     1        CFS2(ILH2,ILL,KEL,LSML)
            END DO
          END DO
          IF(IRH2.NE.INEXT)THEN
            DO INH1=1,NONZ1
              IRH1=NONZC1(3,INH1)
              ILH1=NONZC1(4,INH1)
              IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
                DO ILR=1,NLVALR
                  CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                  IRHS=(ILR-1)*KCSZR+IRH1
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO ILL=1,NLVALL
                    ILHS=(ILL-1)*KCSZL+ILH1
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(ILL,INH1)
                  END DO
                END DO
C**OFF-DIAGONAL BLOCK - END
              ELSE
C**DIAGONAL BLOCK - START
                DO ILR=1,NLVALR
                  CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                  IRHS=(ILR-1)*KCSZR+IRH1
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO ILL=1,ILR
                    IF(ILL.EQ.ILR.AND.ILH1.GT.IRH1)GO TO 2004
                    ILHS=(ILL-1)*KCSZL+ILH1
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(ILL,INH1)
2004                CONTINUE
                  END DO
                END DO
C**DIAGONAL BLOCK - END
              END IF
            END DO
          END IF
          IPREV=IRH2
        END DO
        GO TO 5555

4004    CONTINUE
        DO IRH2=1,LCSZR
          DO ILH2=1,LCSZL
            XK1(ILH2,IRH2)=0
          END DO
        END DO
 
C**CONTRACTION SCHEME '1' (K)
        DO INH1=1,NONZ1
          IRH1=NONZC1(3,INH1)
          ILH1=NONZC1(4,INH1)
          JROFF=NONZC1(1,INH1)
          JLOFF=NONZC1(2,INH1)
          DO I=1,NC1
            NCR1(I)=IPC1(JROFF,MC1(I))
            NCL1(I)=IPC1(JLOFF,MC1(I))
          END DO

C**CONTRACTION SCHEME  '2' (L AND/OR N)
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            ILH2=NONZC2(4,INH2)
            IROFF=NONZC2(1,INH2)
            ILOFF=NONZC2(2,INH2)
            DO I=1,NC2
              NCR2(I)=IPC2(IROFF,MC2(I))
              NCL2(I)=IPC2(ILOFF,MC2(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP3(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP3(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1001
            END DO
1001        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP3(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP3(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2001
            END DO
2001        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA3(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP3(MR,1)
              NL1=IP3(ML,1)
              NR2=IP3(MR,2)
              NL2=IP3(ML,2)
              NR3=IP3(MR,3)
              NL3=IP3(ML,3)
              DO I=1,NP3(IND3)
                K=JP3(I,IND3,1)+1
                L=JP3(I,IND3,2)+1
                N=JP3(I,IND3,3)+1
                ZYX=ZYX+CP3(I,IND3)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)
              END DO
C**ANALYTIC
            END IF
            XK1(ILH2,IRH2)=XYZ+ZYX
          END DO

C************************************DGEMM (RHS)
          CALL DGEMM('N','N',LCSZL,NLVALR,LCSZR,1.0D0,XK1(1,1),
     &    KTEMP,CFS2(1,1,KEL,LSMR),ISIZXX,0.0D0,TEMP(1,1),KTEMP)

C************************************DGEMM (LHS)
          CALL DGEMM('T','N',NLVALL,NLVALR,LCSZL,1.0D0,
     &    CFS2(1,1,KEL,LSML),ISIZXX,TEMP,KTEMP,0.0D0,XCON2(1,1),MAXVAL)

          IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
            DO ILR=1,NLVALR
              IRHS=(ILR-1)*KCSZR+IRH1
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO ILL=1,NLVALL
                ILHS=(ILL-1)*KCSZL+ILH1
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(ILL,ILR)
              END DO
            END DO
C**OFF-DIAGONAL BLOCK - END
          ELSE
C**DIAGONAL BLOCK - START
            DO ILR=1,NLVALR
              IRHS=(ILR-1)*KCSZR+IRH1
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO ILL=1,ILR
                IF(ILL.EQ.ILR.AND.ILH1.GT.IRH1)GO TO 2002
                ILHS=(ILL-1)*KCSZL+ILH1
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(ILL,ILR)
2002            CONTINUE
              END DO
            END DO
C**DIAGONAL BLOCK - END
          END IF
        END DO
      END IF
C*************************************

5555  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE
      IF(ITIM3A.EQ.0)THEN
        IF(IND.NE.0)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM3A=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCV3A(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1ITMODE,H1,XQ1,H2,XQ2,H3,XQ3,HTAU,XQTAU,NN1,MM1,NN2,MM2,NN3,MM3,
     2NNTAU,MMTAU,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP4,
     3ISIZE4,XA4,
     3XRA4,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,
     4JCI2,JCI3,JCIM,X0,Y0,Z0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     5KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I25)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,15),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,15),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      DIMENSION X(25),Y(25),C(25)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC(IPSIZE,1),IP4(ISIZE4,4)
      DIMENSION XA(1),XA4(1),XK(ICSIZE,ICSIZE)
      DIMENSION TEMP(ICSIZE,NVAL),XCON(NVAL,NVAL)
      DIMENSION TEMP1(I25,JCI1,JCI2,JCI3,JCI1,JCI2,JCI3)
      DIMENSION TEMP2(I25,JCI2,JCI3,JCI2,JCI3),TEMP3(I25,JCI3,JCI3)
      DIMENSION X0(I25,JCI1,JCI1,MM1),Y0(I25,JCI2,JCI2,MM2)
      DIMENSION Z0(I25,JCI3,JCI3,MM3),T0(I25,JCIM,JCIM,MMTAU)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      IF(N3.EQ.NT.AND.MDT.EQ.2)MD3=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      IF(N1T.EQ.N3.AND.MDT.EQ.2)MD3=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.EQ.2)MD3=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.EQ.2)MD3=1
      CALL VDCV3A(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,ITMODE,
     1H1,XQ1,H2,XQ2,H3,XQ3,HTAU,XQTAU,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,
     2NN3,MM3,MM3/MD3,NNTAU,MMTAU,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,
     3ISIZE,
     3IP4,ISIZE4,XA4,XRA4,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,TEMP1,
     4TEMP2,TEMP3,JCI1,JCI2,JCI3,JCIM,X0,Y0,Z0,T0,VP,VPR,VC,VCR,VR,VRR,
     5J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I25)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCV3A(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1ITMODE,H1,XQ1,H2,XQ2,H3,XQ3,HTAU,XQTAU,NN1,MH1,MM1,NN2,MH2,MM2,
     2NN3,MH3,MM3,NNTAU,MMTAU,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,
     3ISIZE,IP4,
     3ISIZE4,XA4,XRA4,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,
     4TEMP3,JCI1,JCI2,JCI3,JCIM,X0,Y0,Z0,T0,VP,VPR,VC,VCR,VR,VRR,J21,
     5KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I25)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,15),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,15),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      DIMENSION X(25),Y(25),C(25)
      DIMENSION MODINT(NMODE)
C     DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC(IPSIZE,1),IP4(ISIZE4,4)
      DIMENSION XA(1),XA4(1),XK(ICSIZE,ICSIZE)
      DIMENSION TEMP(ICSIZE,NVAL),XCON(NVAL,NVAL)
      DIMENSION TEMP1(I25,JCI1,JCI2,JCI3,JCI1,JCI2,JCI3)
      DIMENSION TEMP2(I25,JCI2,JCI3,JCI2,JCI3),TEMP3(I25,JCI3,JCI3)
      DIMENSION X0(I25,JCI1,JCI1,MM1),Y0(I25,JCI2,JCI2,MM2)
      DIMENSION Z0(I25,JCI3,JCI3,MM3),T0(I25,JCIM,JCIM,MMTAU)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM3A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VCCV3A'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1', 'MOD2',
C**'MOD3' AND 'TAU'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0CV3
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=JNTFAC(NMODE,ICOUPL,3)
      IFACTC=JNTFAC(NMODE,ICOUPC,3)
C**INTRINSIC
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1

      KA=KROT/2
      LMAX=1+MOD(KA,2)
      FACTRC=0.D0
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      IF(N3.EQ.NT.AND.MDT.EQ.2)MD3=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      IF(N1T.EQ.N3.AND.MDT.EQ.2)MD3=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.EQ.2)MD3=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.EQ.2)MD3=1
      MD=MD1*MD2*MD3*MDT
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          DO NRR=1,JCIM
            NR=NRR+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,25
              X(I)=X(1)
            END DO
            X(7)=HTAU(NR,MTAU,2,LMAX)*MD
            X(14)=X(7)
            X(19)=X(7)
            X(22)=X(7)
            X(24)=X(7)
            DO NLL=1,JCIM
              NL=NLL+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO I=2,25
                Y(I)=Y(1)
              END DO
              Y(8)=HTAU(NL,MTAU,2,LMAX)
              Y(15)=Y(8)
              Y(20)=Y(8)
              Y(23)=Y(8)
              Y(24)=Y(8)
              DO K=1,I25
                T0(K,NLL,NRR,MTAU)=Y(26-K)*X(26-K)
              END DO
            END DO
          END DO
        ELSE
          DO NRR=1,JCIM
            NR=2*NRR-MOD(NRR,2)+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,15
              X(I)=X(1)
            END DO
            X(4)=HTAU(NR,MTAU,2,LMAX)*MD
            X(8)=X(4)
            X(11)=X(4)
            X(13)=X(4)
            X(14)=HTAU(NR,MTAU,3,LMAX)*MD
            X(16)=HTAU(NR,MTAU,1,LMAX)*FACTRC*MD
            DO NLL=1,JCIM
              NL=2*NLL-MOD(NLL,2)+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO K=1,16
                T0(K,NLL,NRR,MTAU)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M1=1,MM1/MD1
      DO M1=1,MM1
        IF(.NOT.LINEAR)THEN
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,25
              X(I)=X(2)
            END DO
            X(9)=X(1)
            X(11)=X(1)
            X(13)=X(1)
            X(15)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              Y(2)=H1(NL1,M1,2)
              DO I=3,25
                Y(I)=Y(1)
              END DO
              Y(9)=Y(2)
              Y(10)=Y(2)
              Y(12)=Y(2)
              Y(14)=Y(2)
              DO K=1,I25
                X0(K,NL1,NR1,M1)=Y(26-K)*X(26-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,16
              X(I)=X(2)
            END DO
            X(5)=H1(NR1,M1,3)
            X(6)=X(1)
            X(7)=X(1)
            X(8)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              DO K=1,16
                X0(K,NL1,NR1,M1)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M2=1,MM2/MD2
      DO M2=1,MM2
        IF(.NOT.LINEAR)THEN
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            DO I=2,25
              X(I)=X(1)
            END DO
            X(3)=H2(NR2,M2,2)
            X(10)=X(3)
            X(16)=X(3)
            X(18)=X(3)
            X(20)=X(3)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO I=2,25
                Y(I)=Y(1)
              END DO
              Y(4)=H2(NL2,M2,2)
              Y(11)=Y(4)
              Y(16)=Y(4)
              Y(17)=Y(4)
              Y(19)=Y(4)
              DO K=1,I25
                Y0(K,NL2,NR2,M2)=Y(26-K)*X(26-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            X(2)=H2(NR2,M2,2)
            DO I=3,16
              X(I)=X(1)
            END DO
            X(6)=X(2)
            X(9)=H2(NR2,M2,3)
            X(10)=X(2)
            X(11)=X(2)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO K=1,16
                Y0(K,NL2,NR2,M2)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M3=1,MM3/MD3
      DO M3=1,MM3
        IF(.NOT.LINEAR)THEN
          DO NR3=1,JCI3
            X(1)=H3(NR3,M3,1)
            DO I=2,25
              X(I)=X(1)
            END DO
            X(5)=H3(NR3,M3,2)
            X(12)=X(5)
            X(17)=X(5)
            X(21)=X(5)
            X(23)=X(5)
            DO NL3=1,JCI3
              Y(1)=H3(NL3,M3,1)
              DO I=2,25
                Y(I)=Y(1)
              END DO
              Y(6)=H3(NL3,M3,2)
              Y(13)=Y(6)
              Y(18)=Y(6)
              Y(21)=Y(6)
              Y(22)=Y(6)
              DO K=1,I25
                Z0(K,NL3,NR3,M3)=Y(26-K)*X(26-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR3=1,JCI3
            X(1)=H3(NR3,M3,1)
            DO I=2,16
              X(I)=X(1)
            END DO
            X(3)=H3(NR3,M3,2)
            X(7)=X(3)
            X(10)=X(3)
            X(12)=H3(NR3,M3,3)
            X(13)=X(3)
            DO NL3=1,JCI3
              Y(1)=H3(NL3,M3,1)
              DO K=1,16
                Z0(K,NL3,NR3,M3)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
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
          IF(J21.GT.1.AND.ICOUPC.GE.3)READ(63)VR
          IF(ICOUPC.GE.3)READ(83)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPC.GE.3)READ(63)VRR
          IF(ICOUPC.GE.3)READ(83)VCR
        END IF
        IF(JCOUPL.GT.0)THEN
          READ(73)VP
        ELSE
          READ(73)VPR
        END IF

C***********************************************************

        IF(.NOT.LINEAR)THEN
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO IRHS1=1,JCI1
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO ILHS1=1,JCI1
                      DO K=1,I25
                      TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=0.D0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 3-MODE INTEGRATION
C         DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO K=1,I25
                      TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=0.D0
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
C           DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS3=1,JCI3
                DO ILHS3=1,JCI3
                  DO K=1,I25
                    TEMP3(K,ILHS3,IRHS3)=0.D0
                  END DO
                END DO
              END DO
C**START 1-MODE INTEGRATION
C             DO M3=1,MM3/MD3
              DO M3=1,MM3
                DO I=1,25
                  C(I)=0
                END DO
                IF(JCOUPL.GT.0)THEN
                  IF(ICOUPC.GE.3)THEN
C**NO WATSON TERM IF RPH
                    C(1)=VC(M3,M2,M1,1)*IFACTC
                    C(2)=VC(M3,M2,M1,1)*IFACTC
                    C(3)=VC(M3,M2,M1,2)*IFACTC
                    C(4)=VC(M3,M2,M1,2)*IFACTC
                    C(5)=VC(M3,M2,M1,3)*IFACTC
                    C(6)=VC(M3,M2,M1,3)*IFACTC
                    C(7)=VC(M3,M2,M1,4)*IFACTC
                    C(8)=VC(M3,M2,M1,4)*IFACTC
                    C(9)=VC(M3,M2,M1,5)*IFACTC
                    C(10)=VC(M3,M2,M1,6)*IFACTC
                    C(11)=VC(M3,M2,M1,6)*IFACTC
                    C(12)=VC(M3,M2,M1,7)*IFACTC
                    C(13)=VC(M3,M2,M1,7)*IFACTC
                    C(14)=VC(M3,M2,M1,8)*IFACTC
                    C(15)=VC(M3,M2,M1,8)*IFACTC
                    C(16)=VC(M3,M2,M1,9)*IFACTC
                    C(17)=VC(M3,M2,M1,10)*IFACTC
                    C(18)=VC(M3,M2,M1,10)*IFACTC
                    C(19)=VC(M3,M2,M1,11)*IFACTC
                    C(20)=VC(M3,M2,M1,11)*IFACTC
                    C(21)=VC(M3,M2,M1,12)*IFACTC
                    C(22)=VC(M3,M2,M1,13)*IFACTC
                    C(23)=VC(M3,M2,M1,13)*IFACTC
                    C(24)=VC(M3,M2,M1,14)*IFACTC
                    C(25)=VC(M3,M2,M1,15)*IFACTC
                    IF(J21.GT.1)C(25)=C(25)+VR(KROT,M3,M2,M1)*IFACTC
                  END IF
                  C(25)=C(25)+VP(M3,M2,M1)*IFACTL
                ELSE
                  IF(ICOUPC.GE.3)THEN
C**NO WATSON TERM IF RPH
                    C(1)=VCR(M3,M2,M1,1)*IFACTC
                    C(2)=VCR(M3,M2,M1,1)*IFACTC
                    C(3)=VCR(M3,M2,M1,2)*IFACTC
                    C(4)=VCR(M3,M2,M1,2)*IFACTC
                    C(5)=VCR(M3,M2,M1,3)*IFACTC
                    C(6)=VCR(M3,M2,M1,3)*IFACTC
                    C(7)=VCR(M3,M2,M1,4)*IFACTC
                    C(8)=VCR(M3,M2,M1,4)*IFACTC
                    C(9)=VCR(M3,M2,M1,5)*IFACTC
                    C(10)=VCR(M3,M2,M1,6)*IFACTC
                    C(11)=VCR(M3,M2,M1,6)*IFACTC
                    C(12)=VCR(M3,M2,M1,7)*IFACTC
                    C(13)=VCR(M3,M2,M1,7)*IFACTC
                    C(14)=VCR(M3,M2,M1,8)*IFACTC
                    C(15)=VCR(M3,M2,M1,8)*IFACTC
                    C(16)=VCR(M3,M2,M1,9)*IFACTC
                    C(17)=VCR(M3,M2,M1,10)*IFACTC
                    C(18)=VCR(M3,M2,M1,10)*IFACTC
                    C(19)=VCR(M3,M2,M1,11)*IFACTC
                    C(20)=VCR(M3,M2,M1,11)*IFACTC
                    C(21)=VCR(M3,M2,M1,12)*IFACTC
                    C(22)=VCR(M3,M2,M1,13)*IFACTC
                    C(23)=VCR(M3,M2,M1,13)*IFACTC
                    C(24)=VCR(M3,M2,M1,14)*IFACTC
                    C(25)=VCR(M3,M2,M1,15)*IFACTC
                    IF(J21.GT.1)C(25)=C(25)+VRR(KROT,M3,M2,M1)*IFACTC
                  END IF
                  C(25)=C(25)+VPR(M3,M2,M1)*IFACTL
                END IF
                DO IRHS3=1,JCI3
                  DO ILHS3=1,JCI3
                    DO K=1,I25
                      TEMP3(K,ILHS3,IRHS3)=TEMP3(K,ILHS3,IRHS3)+
     1                Z0(K,ILHS3,IRHS3,M3)*C(26-K)
                    END DO
                  END DO
                END DO
              END DO
C**END 1-MODE INTEGRATION
              DO IRHS3=1,JCI3
                DO IRHS2=1,JCI2
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO K=1,I25
                        TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1                  TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2                  Y0(K,ILHS2,IRHS2,M2)*TEMP3(K,ILHS3,IRHS3)
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
                        DO K=1,I25
                          TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=
     1                    TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)+
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

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
          DO IRHS=1,NSIZE
            NR1=IP4(IRHS,1)
            NR2=IP4(IRHS,2)
            NR3=IP4(IRHS,3)
            IRTAU=IP4(IRHS,4)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP4(ILHS,1)
              NL2=IP4(ILHS,2)
              NL3=IP4(ILHS,3)
              ILTAU=IP4(ILHS,4)
              DO K=1,I25
                XA4(ILHS+J0)=XA4(ILHS+J0)+
     1          TEMP1(K,NL1,NL2,NL3,NR1,NR2,NR3)*
     2          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
              END DO
            END DO
          END DO
        ELSE
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO IRHS1=1,JCI1
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO ILHS1=1,JCI1
                      DO K=1,16
                      TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=0.D0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 3-MODE INTEGRATION
C         DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO K=1,16
                      TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=0.D0
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
C           DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS3=1,JCI3
                DO ILHS3=1,JCI3
                  DO K=1,16
                    TEMP3(K,ILHS3,IRHS3)=0.D0
                  END DO
                END DO
              END DO
C**START 1-MODE INTEGRATION
C             DO M3=1,MM3/MD3
              DO M3=1,MM3
                DO I=1,14
                  C(I)=0
                END DO
                IF(JCOUPL.GT.0)THEN
                  IF(ICOUPC.GE.3)THEN
C**(15) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                    DO K=1,14
                      C(K)=VC(M3,M2,M1,K)
                    END DO
                    C(6)=2*C(6)
                    C(7)=2*C(7)
                    C(8)=2*C(8)
                    C(10)=2*C(10)
                    C(11)=2*C(11)
                    C(13)=2*C(13)
                  END IF
                  DO IRHS3=1,JCI3
                    DO ILHS3=1,JCI3
                      DO K=1,14
                        TEMP3(K,ILHS3,IRHS3)=TEMP3(K,ILHS3,IRHS3)-
     1                  Z0(K,ILHS3,IRHS3,M3)*C(K)*IFACTC
                      END DO
                      TEMP3(15,ILHS3,IRHS3)=TEMP3(15,ILHS3,IRHS3)+
     1                Z0(15,ILHS3,IRHS3,M3)*VP(M3,M2,M1)*IFACTL
                      TEMP3(16,ILHS3,IRHS3)=TEMP3(16,ILHS3,IRHS3)+
     1                Z0(16,ILHS3,IRHS3,M3)*VR(KROT,M3,M2,M1)
                    END DO
                  END DO
                ELSE
                  IF(ICOUPC.GE.3)THEN
C**(15) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                    DO K=1,14
                      C(K)=VCR(M3,M2,M1,K)
                    END DO
                    C(6)=2*C(6)
                    C(7)=2*C(7)
                    C(8)=2*C(8)
                    C(10)=2*C(10)
                    C(11)=2*C(11)
                    C(13)=2*C(13)
                  END IF
                  DO IRHS3=1,JCI3
                    DO ILHS3=1,JCI3
                      DO K=1,14
                        TEMP3(K,ILHS3,IRHS3)=TEMP3(K,ILHS3,IRHS3)-
     1                  Z0(K,ILHS3,IRHS3,M3)*C(K)*IFACTC
                      END DO
                      TEMP3(15,ILHS3,IRHS3)=TEMP3(15,ILHS3,IRHS3)+
     1                Z0(15,ILHS3,IRHS3,M3)*VPR(M3,M2,M1)*IFACTL
                      TEMP3(16,ILHS3,IRHS3)=TEMP3(16,ILHS3,IRHS3)+
     1                Z0(16,ILHS3,IRHS3,M3)*VRR(KROT,M3,M2,M1)
                    END DO
                  END DO
                END IF
              END DO
C**END 1-MODE INTEGRATION
              DO IRHS3=1,JCI3
                DO IRHS2=1,JCI2
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO K=1,16
                        TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1                  TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2                  Y0(K,ILHS2,IRHS2,M2)*TEMP3(K,ILHS3,IRHS3)
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
                        DO K=1,16
                          TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=
     1                    TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)+
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

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
          DO IRHS=1,NSIZE
            NR1=IP4(IRHS,1)
            NR2=IP4(IRHS,2)
            NR3=IP4(IRHS,3)
            IRRTAU=IP4(IRHS,4)
            IRTAU=IRRTAU/2+MOD(IRRTAU,2)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP4(ILHS,1)
              NL2=IP4(ILHS,2)
              NL3=IP4(ILHS,3)
              ILLTAU=IP4(ILHS,4)
              ILTAU=ILLTAU/2+MOD(ILLTAU,2)
              DO K=1,16
                XA4(ILHS+J0)=XA4(ILHS+J0)+
     1          TEMP1(K,NL1,NL2,NL3,NR1,NR2,NR3)*T0(K,ILTAU,IRTAU,MTAU)
              END DO
            END DO
          END DO
        END IF
      END DO
C**END TAU LOOP (4-MODE INTEGRATION)
      CALL MATOUT(XA4,XRA4,NSIZE,23)
      GO TO 7777
6666  CALL MATIN(XA4,XRA4,NSIZE,23)
7777  CONTINUE

C***********************************ALGORITHM FROM VCV3

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR1=IPC(IROFF,MODE1)
        NR2=IPC(IROFF,MODE2)
        NR3=IPC(IROFF,MODE3)
        NRTAU=IPC(IROFF,ITMODE)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP4(IR,1).AND.NR2.EQ.IP4(IR,2).AND.
     1       NR3.EQ.IP4(IR,3).AND.NRTAU.EQ.IP4(IR,4))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1      ITMODE.AND.(IPC(IROFF,K).NE.IPC(ILOFF,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPC(ILOFF,MODE1)
          NL2=IPC(ILOFF,MODE2)
          NL3=IPC(ILOFF,MODE3)
          NLTAU=IPC(ILOFF,ITMODE)
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
          XYZ=XA4(I)
3000      CONTINUE
          XK(ILHS,IRHS)=XYZ*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
        J0=(IROFF)*(IROFF-1)/2+IOFF
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IP(IROFF,KCONT)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IP(IROFF,MCONT).NE.IP(ILOFF,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IP(ILOFF,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM3A.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM3A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCV3B(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1MODE4,H1,XQ1,H2,XQ2,H3,XQ3,HTAU,XQTAU,NN1,MM1,NN2,MM2,NN3,MM3,
     2NNTAU,MMTAU,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,
     3KCONT1,KCONT2,XA,
     3ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK1,XK2,IRL,KTEMP,XCON2,NVAL1,
     4NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,JCIM,X0,Y0,Z0,
     5T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,
     6NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I25,MAXVAL,MAXQ4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,15),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,15),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP4(ISIZE4,4),IRL(MAXQ4,MAXQ4,MAXQ4,MAXQ4)
      DIMENSION TEMP1(I25,JCI1,JCI2,JCI3,JCI1,JCI2,JCI3)
      DIMENSION TEMP2(I25,JCI2,JCI3,JCI2,JCI3),TEMP3(I25,JCI3,JCI3)
      DIMENSION X0(I25,JCI1,JCI1,MM1),Y0(I25,JCI2,JCI2,MM2)
      DIMENSION Z0(I25,JCI3,JCI3,MM3),T0(I25,JCIM,JCIM,MMTAU)
      DIMENSION X(25),Y(25),C(25)
C****************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA4(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C****************************************************************
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      IF(N3.EQ.NT.AND.MDT.EQ.2)MD3=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      IF(N1T.EQ.N3.AND.MDT.EQ.2)MD3=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.EQ.2)MD3=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.EQ.2)MD3=1
      CALL VDCV3B(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,MODE4,
     1H1,XQ1,H2,XQ2,H3,XQ3,HTAU,XQTAU,NN1,MM1,MM1/MD1,NN2,MM2,MM2/MD2,
     2NN3,MM3,MM3/MD3,NNTAU,MMTAU,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,
     3ICSIZ2,IPSIZ2,
     3KCONT1,KCONT2,XA,ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK1,XK2,IRL,
     4KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,JCI2,
     5JCI3,JCIM,X0,Y0,Z0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,
     6NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I25,
     7MAXVAL,MAXQ4)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCV3B(NMODE,NNMODE,MOD1,MOD2,MOD3,MODE1,MODE2,MODE3,
     1MODE4,H1,XQ1,H2,XQ2,H3,XQ3,HTAU,XQTAU,NN1,MH1,MM1,NN2,MH2,MM2,
     2NN3,MH3,MM3,NNTAU,MMTAU,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,
     3IPSIZ2,KCONT1,
     3KCONT2,XA,ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK1,XK2,IRL,KTEMP,
     4XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,
     5JCIM,X0,Y0,Z0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,
     6CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I25,MAXVAL,
     7MAXQ4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,15),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,15),VRR(J21,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD

C**TEMPORARY
      DIMENSION IRL(MAXQ4,MAXQ4,MAXQ4,MAXQ4),NX4(4),MX4(4)
C**TEMPORARY

      DIMENSION MODINT(NMODE)
C     DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP4(ISIZE4,4),TEMP(KTEMP,1)
      DIMENSION TEMP1(I25,JCI1,JCI2,JCI3,JCI1,JCI2,JCI3)
      DIMENSION TEMP2(I25,JCI2,JCI3,JCI2,JCI3),TEMP3(I25,JCI3,JCI3)
      DIMENSION X0(I25,JCI1,JCI1,MM1),Y0(I25,JCI2,JCI2,MM2)
      DIMENSION Z0(I25,JCI3,JCI3,MM3),T0(I25,JCIM,JCIM,MMTAU)
      DIMENSION X(25),Y(25),C(25)
C*****************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA4(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C*****************************************************************
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM3B.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCV3B'
          CALL TIMIT(1)
          CALL FLUSH(IOUT)
        END IF
      END IF

C**************************************************************
C**************************************************************
C**BOTH CONTRACTION SCHEMES HAVE AT LEAST ONE MODE
C**CHECK SCHEME '1'
      I1=KCONT1
      NC1=0
      DO NN=1,ICONT(I1)
        IF(MOD1.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE1
          JC1(NC1)=1
        END IF
        IF(MOD2.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE2
          JC1(NC1)=2
        END IF
        IF(MOD3.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE3
          JC1(NC1)=3
        END IF
        IF(NNMODE.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE4
          JC1(NC1)=4
        END IF
      END DO
C**CHECK SCHEME '2'
      I2=KCONT2
      NC2=0
      DO NN=1,ICONT(I2)
        IF(MOD1.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE1
          JC2(NC2)=1
        END IF
        IF(MOD2.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE2
          JC2(NC2)=2
        END IF
        IF(MOD3.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE3
          JC2(NC2)=3
        END IF
        IF(NNMODE.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE4
          JC2(NC2)=4
        END IF
      END DO
C**NKMOD1 IS TOTAL NUMBER OF MODES FOR SCHEME INVOLVING 'MOD1' (K)
C**NKMOD2 IS TOTAL NUMBER OF MODES FOR THE OTHER SCHEME (ONE OR ALL
C**OF 'MOD2' (L), 'MOD3' (N) AND 'MOD4' (M))
C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L),
C**'MOD3' (N) OR 'MOD4' (M)
C**IF THEY ARE IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
C**************************************************************
      IF(IND.NE.0)GO TO 7777
C**************************************************************

C**FIND MAX QUANTUM
      MAXQ4=0
      DO IR=1,NSIZE
        NR1=IP4(IR,1)
        NR2=IP4(IR,2)
        NR3=IP4(IR,3)
        NR4=IP4(IR,4)
        IF(NR1.GT.MAXQ4)MAXQ4=NR1
        IF(NR2.GT.MAXQ4)MAXQ4=NR2
        IF(NR3.GT.MAXQ4)MAXQ4=NR3
        IF(NR4.GT.MAXQ4)MAXQ4=NR4
      END DO

C***********************************ALGORITHM FROM V0CV3
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=JNTFAC(NMODE,ICOUPL,3)
      IFACTC=JNTFAC(NMODE,ICOUPC,3)
C**INTRINSIC
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1

      KA=KROT/2
      LMAX=1+MOD(KA,2)
      FACTRC=0.D0
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      IF(N3.EQ.NT.AND.MDT.EQ.2)MD3=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      IF(N1T.EQ.N3.AND.MDT.EQ.2)MD3=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.EQ.2)MD3=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.EQ.2)MD3=1
      MD=MD1*MD2*MD3*MDT
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          DO NRR=1,JCIM
            NR=NRR+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,25
              X(I)=X(1)
            END DO
            X(7)=HTAU(NR,MTAU,2,LMAX)*MD
            X(14)=X(7)
            X(19)=X(7)
            X(22)=X(7)
            X(24)=X(7)
            DO NLL=1,JCIM
              NL=NLL+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO I=2,25
                Y(I)=Y(1)
              END DO
              Y(8)=HTAU(NL,MTAU,2,LMAX)
              Y(15)=Y(8)
              Y(20)=Y(8)
              Y(23)=Y(8)
              Y(24)=Y(8)
              DO K=1,I25
                T0(K,NLL,NRR,MTAU)=Y(26-K)*X(26-K)
              END DO
            END DO
          END DO
        ELSE
          DO NRR=1,JCIM
            NR=2*NRR-MOD(NRR,2)+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,15
              X(I)=X(1)
            END DO
            X(4)=HTAU(NR,MTAU,2,LMAX)*MD
            X(8)=X(4)
            X(11)=X(4)
            X(13)=X(4)
            X(14)=HTAU(NR,MTAU,3,LMAX)*MD
            X(16)=HTAU(NR,MTAU,1,LMAX)*FACTRC*MD
            DO NLL=1,JCIM
              NL=2*NLL-MOD(NLL,2)+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO K=1,16
                T0(K,NLL,NRR,MTAU)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M1=1,MM1/MD1
      DO M1=1,MM1
        IF(.NOT.LINEAR)THEN
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,25
              X(I)=X(2)
            END DO
            X(9)=X(1)
            X(11)=X(1)
            X(13)=X(1)
            X(15)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              Y(2)=H1(NL1,M1,2)
              DO I=3,25
                Y(I)=Y(1)
              END DO
              Y(9)=Y(2)
              Y(10)=Y(2)
              Y(12)=Y(2)
              Y(14)=Y(2)
              DO K=1,I25
                X0(K,NL1,NR1,M1)=Y(26-K)*X(26-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,16
              X(I)=X(2)
            END DO
            X(5)=H1(NR1,M1,3)
            X(6)=X(1)
            X(7)=X(1)
            X(8)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              DO K=1,16
                X0(K,NL1,NR1,M1)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M2=1,MM2/MD2
      DO M2=1,MM2
        IF(.NOT.LINEAR)THEN
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            DO I=2,25
              X(I)=X(1)
            END DO
            X(3)=H2(NR2,M2,2)
            X(10)=X(3)
            X(16)=X(3)
            X(18)=X(3)
            X(20)=X(3)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO I=2,25
                Y(I)=Y(1)
              END DO
              Y(4)=H2(NL2,M2,2)
              Y(11)=Y(4)
              Y(16)=Y(4)
              Y(17)=Y(4)
              Y(19)=Y(4)
              DO K=1,I25
                Y0(K,NL2,NR2,M2)=Y(26-K)*X(26-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            X(2)=H2(NR2,M2,2)
            DO I=3,16
              X(I)=X(1)
            END DO
            X(6)=X(2)
            X(9)=H2(NR2,M2,3)
            X(10)=X(2)
            X(11)=X(2)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO K=1,16
                Y0(K,NL2,NR2,M2)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M3=1,MM3/MD3
      DO M3=1,MM3
        IF(.NOT.LINEAR)THEN
          DO NR3=1,JCI3
            X(1)=H3(NR3,M3,1)
            DO I=2,25
              X(I)=X(1)
            END DO
            X(5)=H3(NR3,M3,2)
            X(12)=X(5)
            X(17)=X(5)
            X(21)=X(5)
            X(23)=X(5)
            DO NL3=1,JCI3
              Y(1)=H3(NL3,M3,1)
              DO I=2,25
                Y(I)=Y(1)
              END DO
              Y(6)=H3(NL3,M3,2)
              Y(13)=Y(6)
              Y(18)=Y(6)
              Y(21)=Y(6)
              Y(22)=Y(6)
              DO K=1,I25
                Z0(K,NL3,NR3,M3)=Y(26-K)*X(26-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR3=1,JCI3
            X(1)=H3(NR3,M3,1)
            DO I=2,16
              X(I)=X(1)
            END DO
            X(3)=H3(NR3,M3,2)
            X(7)=X(3)
            X(10)=X(3)
            X(12)=H3(NR3,M3,3)
            X(13)=X(3)
            DO NL3=1,JCI3
              Y(1)=H3(NL3,M3,1)
              DO K=1,16
                Z0(K,NL3,NR3,M3)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
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
          IF(J21.GT.1.AND.ICOUPC.GE.3)READ(63)VR
          IF(ICOUPC.GE.3)READ(83)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPC.GE.3)READ(63)VRR
          IF(ICOUPC.GE.3)READ(83)VCR
        END IF
        IF(JCOUPL.GT.0)THEN
          READ(73)VP
        ELSE
          READ(73)VPR
        END IF

C***********************************************************

        IF(.NOT.LINEAR)THEN
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO IRHS1=1,JCI1
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO ILHS1=1,JCI1
                      DO K=1,I25
                      TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=0.D0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 3-MODE INTEGRATION
C         DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO K=1,I25
                      TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=0.D0
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
C           DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS3=1,JCI3
                DO ILHS3=1,JCI3
                  DO K=1,I25
                    TEMP3(K,ILHS3,IRHS3)=0.D0
                  END DO
                END DO
              END DO
C**START 1-MODE INTEGRATION
C             DO M3=1,MM3/MD3
              DO M3=1,MM3
                DO I=1,25
                  C(I)=0
                END DO
                IF(JCOUPL.GT.0)THEN
                  IF(ICOUPC.GE.3)THEN
C**NO WATSON TERM IF RPH
                    C(1)=VC(M3,M2,M1,1)*IFACTC
                    C(2)=VC(M3,M2,M1,1)*IFACTC
                    C(3)=VC(M3,M2,M1,2)*IFACTC
                    C(4)=VC(M3,M2,M1,2)*IFACTC
                    C(5)=VC(M3,M2,M1,3)*IFACTC
                    C(6)=VC(M3,M2,M1,3)*IFACTC
                    C(7)=VC(M3,M2,M1,4)*IFACTC
                    C(8)=VC(M3,M2,M1,4)*IFACTC
                    C(9)=VC(M3,M2,M1,5)*IFACTC
                    C(10)=VC(M3,M2,M1,6)*IFACTC
                    C(11)=VC(M3,M2,M1,6)*IFACTC
                    C(12)=VC(M3,M2,M1,7)*IFACTC
                    C(13)=VC(M3,M2,M1,7)*IFACTC
                    C(14)=VC(M3,M2,M1,8)*IFACTC
                    C(15)=VC(M3,M2,M1,8)*IFACTC
                    C(16)=VC(M3,M2,M1,9)*IFACTC
                    C(17)=VC(M3,M2,M1,10)*IFACTC
                    C(18)=VC(M3,M2,M1,10)*IFACTC
                    C(19)=VC(M3,M2,M1,11)*IFACTC
                    C(20)=VC(M3,M2,M1,11)*IFACTC
                    C(21)=VC(M3,M2,M1,12)*IFACTC
                    C(22)=VC(M3,M2,M1,13)*IFACTC
                    C(23)=VC(M3,M2,M1,13)*IFACTC
                    C(24)=VC(M3,M2,M1,14)*IFACTC
                    C(25)=VC(M3,M2,M1,15)*IFACTC
                    IF(J21.GT.1)C(25)=C(25)+VR(KROT,M3,M2,M1)*IFACTC
                  END IF
                  C(25)=C(25)+VP(M3,M2,M1)*IFACTL
                ELSE
                  IF(ICOUPC.GE.3)THEN
C**NO WATSON TERM IF RPH
                    C(1)=VCR(M3,M2,M1,1)*IFACTC
                    C(2)=VCR(M3,M2,M1,1)*IFACTC
                    C(3)=VCR(M3,M2,M1,2)*IFACTC
                    C(4)=VCR(M3,M2,M1,2)*IFACTC
                    C(5)=VCR(M3,M2,M1,3)*IFACTC
                    C(6)=VCR(M3,M2,M1,3)*IFACTC
                    C(7)=VCR(M3,M2,M1,4)*IFACTC
                    C(8)=VCR(M3,M2,M1,4)*IFACTC
                    C(9)=VCR(M3,M2,M1,5)*IFACTC
                    C(10)=VCR(M3,M2,M1,6)*IFACTC
                    C(11)=VCR(M3,M2,M1,6)*IFACTC
                    C(12)=VCR(M3,M2,M1,7)*IFACTC
                    C(13)=VCR(M3,M2,M1,7)*IFACTC
                    C(14)=VCR(M3,M2,M1,8)*IFACTC
                    C(15)=VCR(M3,M2,M1,8)*IFACTC
                    C(16)=VCR(M3,M2,M1,9)*IFACTC
                    C(17)=VCR(M3,M2,M1,10)*IFACTC
                    C(18)=VCR(M3,M2,M1,10)*IFACTC
                    C(19)=VCR(M3,M2,M1,11)*IFACTC
                    C(20)=VCR(M3,M2,M1,11)*IFACTC
                    C(21)=VCR(M3,M2,M1,12)*IFACTC
                    C(22)=VCR(M3,M2,M1,13)*IFACTC
                    C(23)=VCR(M3,M2,M1,13)*IFACTC
                    C(24)=VCR(M3,M2,M1,14)*IFACTC
                    C(25)=VCR(M3,M2,M1,15)*IFACTC
                    IF(J21.GT.1)C(25)=C(25)+VRR(KROT,M3,M2,M1)*IFACTC
                  END IF
                  C(25)=C(25)+VPR(M3,M2,M1)*IFACTL
                END IF
                DO IRHS3=1,JCI3
                  DO ILHS3=1,JCI3
                    DO K=1,I25
                      TEMP3(K,ILHS3,IRHS3)=TEMP3(K,ILHS3,IRHS3)+
     1                Z0(K,ILHS3,IRHS3,M3)*C(26-K)
                    END DO
                  END DO
                END DO
              END DO
C**END 1-MODE INTEGRATION
              DO IRHS3=1,JCI3
                DO IRHS2=1,JCI2
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO K=1,I25
                        TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1                  TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2                  Y0(K,ILHS2,IRHS2,M2)*TEMP3(K,ILHS3,IRHS3)
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
                        DO K=1,I25
                          TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=
     1                    TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)+
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

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
          DO IRHS=1,NSIZE
            NR1=IP4(IRHS,1)
            NR2=IP4(IRHS,2)
            NR3=IP4(IRHS,3)
            IRTAU=IP4(IRHS,4)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP4(ILHS,1)
              NL2=IP4(ILHS,2)
              NL3=IP4(ILHS,3)
              ILTAU=IP4(ILHS,4)
              DO K=1,I25
                XA4(ILHS+J0)=XA4(ILHS+J0)+
     1          TEMP1(K,NL1,NL2,NL3,NR1,NR2,NR3)*
     2          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
              END DO
            END DO
          END DO
        ELSE
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO IRHS1=1,JCI1
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO ILHS1=1,JCI1
                      DO K=1,16
                      TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=0.D0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 3-MODE INTEGRATION
C         DO M1=1,MM1/MD1
          DO M1=1,MM1
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    DO K=1,16
                      TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=0.D0
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
C           DO M2=1,MM2/MD2
            DO M2=1,MM2
              DO IRHS3=1,JCI3
                DO ILHS3=1,JCI3
                  DO K=1,16
                    TEMP3(K,ILHS3,IRHS3)=0.D0
                  END DO
                END DO
              END DO
C**START 1-MODE INTEGRATION
C             DO M3=1,MM3/MD3
              DO M3=1,MM3
                DO I=1,14
                  C(I)=0
                END DO
                IF(JCOUPL.GT.0)THEN
                  IF(ICOUPC.GE.3)THEN
C**(15) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                    DO K=1,14
                      C(K)=VC(M3,M2,M1,K)
                    END DO
                    C(6)=2*C(6)
                    C(7)=2*C(7)
                    C(8)=2*C(8)
                    C(10)=2*C(10)
                    C(11)=2*C(11)
                    C(13)=2*C(13)
                  END IF
                  DO IRHS3=1,JCI3
                    DO ILHS3=1,JCI3
                      DO K=1,14
                        TEMP3(K,ILHS3,IRHS3)=TEMP3(K,ILHS3,IRHS3)-
     1                  Z0(K,ILHS3,IRHS3,M3)*C(K)*ICOUPC
                      END DO
                      TEMP3(15,ILHS3,IRHS3)=TEMP3(15,ILHS3,IRHS3)+
     1                Z0(15,ILHS3,IRHS3,M3)*VP(M3,M2,M1)*ICOUPL
                      TEMP3(16,ILHS3,IRHS3)=TEMP3(16,ILHS3,IRHS3)+
     1                Z0(16,ILHS3,IRHS3,M3)*VR(KROT,M3,M2,M1)
                    END DO
                  END DO
                ELSE
                  IF(ICOUPC.GE.3)THEN
C**(15) IS WATSON CORRECTION TERM ONLY - ZERO IF LINEAR (SEE CORIOL)
                    DO K=1,14
                      C(K)=VCR(M3,M2,M1,K)
                    END DO
                    C(6)=2*C(6)
                    C(7)=2*C(7)
                    C(8)=2*C(8)
                    C(10)=2*C(10)
                    C(11)=2*C(11)
                    C(13)=2*C(13)
                  END IF
                  DO IRHS3=1,JCI3
                    DO ILHS3=1,JCI3
                      DO K=1,14
                        TEMP3(K,ILHS3,IRHS3)=TEMP3(K,ILHS3,IRHS3)-
     1                  Z0(K,ILHS3,IRHS3,M3)*C(K)*ICOUPC
                      END DO
                      TEMP3(15,ILHS3,IRHS3)=TEMP3(15,ILHS3,IRHS3)+
     1                Z0(15,ILHS3,IRHS3,M3)*VPR(M3,M2,M1)*ICOUPL
                      TEMP3(16,ILHS3,IRHS3)=TEMP3(16,ILHS3,IRHS3)+
     1                Z0(16,ILHS3,IRHS3,M3)*VRR(KROT,M3,M2,M1)
                    END DO
                  END DO
                END IF
              END DO
C**END 1-MODE INTEGRATION
              DO IRHS3=1,JCI3
                DO IRHS2=1,JCI2
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO K=1,16
                        TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)=
     1                  TEMP2(K,ILHS2,ILHS3,IRHS2,IRHS3)+
     2                  Y0(K,ILHS2,IRHS2,M2)*TEMP3(K,ILHS3,IRHS3)
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
                        DO K=1,16
                          TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)=
     1                    TEMP1(K,ILHS1,ILHS2,ILHS3,IRHS1,IRHS2,IRHS3)+
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

C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
          DO IRHS=1,NSIZE
            NR1=IP4(IRHS,1)
            NR2=IP4(IRHS,2)
            NR3=IP4(IRHS,3)
            IRRTAU=IP4(IRHS,4)
            IRTAU=IRRTAU/2+MOD(IRRTAU,2)
            J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP4(ILHS,1)
              NL2=IP4(ILHS,2)
              NL3=IP4(ILHS,3)
              ILLTAU=IP4(ILHS,4)
              ILTAU=ILLTAU/2+MOD(ILLTAU,2)
              DO K=1,16
                XA4(ILHS+J0)=XA4(ILHS+J0)+
     1          TEMP1(K,NL1,NL2,NL3,NR1,NR2,NR3)*
     2          T0(K,ILTAU,IRTAU,MTAU)
              END DO
            END DO
          END DO
        END IF
      END DO
C**END TAU LOOP (4-MODE INTEGRATION)

      CALL MATOUT(XA4,XRA4,NSIZE,23)
      GO TO 7777
6666  CALL MATIN(XA4,XRA4,NSIZE,23)
7777  CONTINUE

      IF(IND.NE.0)THEN
        DO IR=1,NSIZE
          NR=0
          DO J=1,NC1
            J1=JC1(J)
            NR1=IP4(IR,J1)
            NR=NR+1
            NX4(NR)=NR1
          END DO
          DO J=1,NC2
            J2=JC2(J)
            NR2=IP4(IR,J2)
            NR=NR+1
            NX4(NR)=NR2
          END DO
          IRL(NX4(1),NX4(2),NX4(3),NX4(4))=IR
        END DO
      END IF

C***********************************ALGORITHM FROM VCV3

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFFR=IOFFR+NCSIZE(I)
      END DO
      J0R=IOFFR*(IOFFR+1)/2
C**TEMPORARY-LAN
      JSR=1
C**TEMPORARY-LAN
C     DO 9999 ISMR1=1,NVSYM
      DO 9999 ISMR1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE

        IF(NONC1.GT.NONC2)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          KCCCL=KCSZL
CC        IF(NONC1.GT.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             KCCCL=KCSZL
C           ELSE
C             KCCCL=IRH1
C           END IF
CC        END IF
          DO ILH1=1,KCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C         DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=IROFF
            NONZC1(2,NONZ1)=ILOFF
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
        IF(NONZ1.EQ.0)GO TO 5556
        NONZC1(3,NONZ1+1)=0
        END IF

      END IF

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR TAU)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE

        IF(NONC1.LE.NONC2)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR TAU)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          LCCCL=LCSZL
CC        IF(NONC1.LE.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             LCCCL=LCSZL
C           ELSE
C             LCCCL=IRH2
C           END IF
CC        END IF
          DO ILH2=1,LCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=JROFF
            NONZC2(2,NONZ2)=JLOFF
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
        IF(NONZ2.EQ.0)GO TO 5555
        NONZC2(3,NONZ2+1)=0
        END IF

      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM3B.EQ.0)THEN
C       WRITE(IOUT,*)'NKVALL,NKVALR = ',NKVALL,NKVALR
C       WRITE(IOUT,*)'NLVALL,NLVALR = ',NLVALL,NLVALR
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C       WRITE(IOUT,*)'NONC1,NONC2 = ',NONC1,NONC2
C     END IF

      IF(NONC1.LE.NONC2)GO TO 3003
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR/ TAU)

      DO IRH2=1,LCSZR
        JROFF=IRH2+LCOFFR
        DO I=1,NKVALR
          DO J=1,NKVALL
            DO K=1,NLVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          LCCCL=LCSZL
        ELSE
          LCCCL=IRH2
        END IF
        DO ILH2=1,LCCCL
          JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
          DO K=1,NKMOD2
            IGOT=0
            DO I=1,NC2
              IF(K.EQ.MC2(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 6001
          END DO
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH2.EQ.IRH2)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NKVALR
            DO J=1,NKVALL
              XCON2(J,I)=0
            END DO
          END DO

          IPREV=0
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            INEXT=NONZC1(3,INH1+1)
            ILH1=NONZC1(4,INH1)
            IROFF=NONZC1(1,INH1)
            ILOFF=NONZC1(2,INH1)
            IF(IPREV.NE.IRH1)THEN
C**NEW RHS
              DO J=1,NKVALL
                XK1(J,1)=0
              END DO
            END IF

            DO I=1,NC1
              NCR1(I)=IPC1(IROFF,MC1(I))
              NCL1(I)=IPC1(ILOFF,MC1(I))
            END DO

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX4(NR)=NR1
        MX4(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX4(NR)=NR2
        MX4(NR)=NL2
      END DO
      IR=IRL(NX4(1),NX4(2),NX4(3),NX4(4))
      IL=IRL(MX4(1),MX4(2),MX4(3),MX4(4))
C**TEMPORARY

C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XXX=XA4(I)

            DO K=1,NKVALL
              XK1(K,1)=XK1(K,1)+XXX*CFS1(ILH1,K,KEL,KSML)
            END DO
5001        CONTINUE

            IF(IRH1.NE.INEXT)THEN

            DO K=1,NKVALR
              DO L=1,NKVALL
                XCON2(L,K)=XCON2(L,K)+XK1(L,1)*CFS1(IRH1,K,KEL,KSMR)
              END DO
            END DO

            END IF
            IPREV=IRH1
          END DO
CC
          DO I=1,NKVALR
            DO J=1,NKVALL
              DO K=1,NLVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6001      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILL,IKR,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKR,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILR,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILR,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO
C*************************************

      GO TO 5555

3003  CONTINUE

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '1' (K)

      DO IRH1=1,KCSZR
        IROFF=IRH1+KCOFFR
        DO I=1,NLVALR
          DO J=1,NLVALL
            DO K=1,NKVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          KCCCL=KCSZL
        ELSE
          KCCCL=IRH1
        END IF
        DO ILH1=1,KCCCL
          ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
          DO K=1,NKMOD1
            IGOT=0
            DO I=1,NC1
              IF(K.EQ.MC1(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 6002
          END DO
          DO I=1,NC1
            NCR1(I)=IPC1(IROFF,MC1(I))
            NCL1(I)=IPC1(ILOFF,MC1(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH1.EQ.IRH1)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NLVALR
            DO J=1,NLVALL
              XCON2(J,I)=0
            END DO
          END DO

C**CONTRACTION SCHEME  '2' (L AND/OR N AND/OR M)

          JPREV=0
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            JNEXT=NONZC2(3,INH2+1)
            ILH2=NONZC2(4,INH2)
            JROFF=NONZC2(1,INH2)
            JLOFF=NONZC2(2,INH2)
            IF(JPREV.NE.IRH2)THEN
C**NEW RHS
              DO J=1,NLVALL
                XK1(J,1)=0
              END DO
            END IF

            DO I=1,NC2
              NCR2(I)=IPC2(JROFF,MC2(I))
              NCL2(I)=IPC2(JLOFF,MC2(I))
            END DO

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX4(NR)=NR1
        MX4(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX4(NR)=NR2
        MX4(NR)=NL2
      END DO
      IR=IRL(NX4(1),NX4(2),NX4(3),NX4(4))
      IL=IRL(MX4(1),MX4(2),MX4(3),MX4(4))
C**TEMPORARY

C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XXX=XA4(I)

            DO J=1,NLVALL
              XK1(J,1)=XK1(J,1)+XXX*CFS2(ILH2,J,KEL,LSML)
            END DO
5002        CONTINUE

            IF(IRH2.NE.JNEXT)THEN

            DO I=1,NLVALR
              DO J=1,NLVALL
                XCON2(J,I)=XCON2(J,I)+XK1(J,1)*CFS2(IRH2,I,KEL,LSMR)
              END DO
            END DO

            END IF
            JPREV=IRH2
          END DO
CC
          DO I=1,NLVALR
            DO J=1,NLVALL
              DO K=1,NKVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6002      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKR,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKL,ILR,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILR,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO

5555  CONTINUE
5556  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM3B.EQ.0)THEN
        IF(IND.NE.0)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM3B=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI4A(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,NN2,MM2,NN3,MM3,
     2NN4,MM4,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP4,ISIZE4,
     3XA4,XRA4,
     3NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,
     4JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     5KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I17,XKAN,MAXQU,MAXPOW,NP4,
     6CP4,JP4,NTOT4,MAX4,INDK,INDL,INDN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,15),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,15),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP4(ISIZE4,4),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I17,JCI4,JCI4),TEMP2(I17,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(I17,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(I17,JCI1,JCI1,MM1),Y0(I17,JCI2,JCI2,MM2)
      DIMENSION Z0(I17,JCI3,JCI3,MM3),W0(I17,JCI4,JCI4,MM4)
      DIMENSION X(17),Y(17),C(17)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA4(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),JP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI4A(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,MM1/MD1,NN2,MM2,
     2MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,IP,ISIZMX,IPC,ICSIZE,
     3IPSIZE,
     3KCONT,XA,ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK,TEMP,XCON,NVAL,
     4ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,
     5W0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,
     6KEL21,NVSYMX,I17,XKAN,MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,MAX4,INDK,
     7INDL,INDN)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI4A(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MH1,MM1,NN2,MH2,MM2,
     2NN3,MH3,MM3,NN4,MH4,MM4,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,
     3ISIZE,IP4,
     3ISIZE4,XA4,XRA4,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,
     4TEMP3,JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VP,VPR,VC,VCR,VR,VRR,
     5J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I17,XKAN,
     6MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,MAX4,INDK,INDL,INDN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,15),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,15),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3),H4(NN4,MH4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP4(ISIZE4,4),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I17,JCI4,JCI4),TEMP2(I17,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(I17,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(I17,JCI1,JCI1,MM1),Y0(I17,JCI2,JCI2,MM2)
      DIMENSION Z0(I17,JCI3,JCI3,MM3),W0(I17,JCI4,JCI4,MM4)
      DIMENSION X(17),Y(17),C(17)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA4(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),JP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM4A.EQ.1)THEN
        WRITE(IOUT,*)'Calculating VCCI4A'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

C** FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1', 'MOD2',
C** 'MOD3' AND 'MOD4'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C**ANALYTIC
      IND=INDN(MOD1)+INDL(MOD2)+INDK(MOD3)+MOD4
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI4
C**TEMPORARY
      IF(ICONDP.NE.0)THEN
        IF(MOLINC.GT.0)GO TO 99999
        GO TO 6666
      END IF
C**TEMPORARY

      IFACTC=INTFAC(NMODE,ICOUPC,4)
      IFACTL=INTFAC(NMODE,ICOUPL,4)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(74)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(74)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.3)READ(64)VR
        IF(ICOUPC.GT.3)READ(84)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.3)READ(64)VRR
        IF(ICOUPC.GT.3)READ(84)VCR
      END IF

      IF(MOLINC.GT.0)GO TO 99999

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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(2)
          X(10)=X(2)
          X(11)=X(2)
          X(12)=X(2)
          X(13)=X(2)
          X(14)=X(2)
          X(15)=X(2)
          X(16)=X(2)
          X(17)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(1)
            Y(5)=Y(3)
            Y(6)=Y(1)
            Y(7)=Y(3)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(3)
            Y(11)=Y(3)
            Y(12)=Y(3)
            Y(13)=Y(3)
            Y(14)=Y(3)
            Y(15)=Y(3)
            Y(16)=Y(3)
            Y(17)=Y(3)
            DO K=1,I17
              X0(K,NL1,NR1,M1)=Y(18-K)*X(18-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(1)
          X(10)=X(2)
          X(11)=X(1)
          X(12)=X(2)
          X(13)=X(1)
          X(14)=X(1)
          X(15)=X(1)
          X(16)=X(1)
          X(17)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(1)
            Y(11)=Y(3)
            Y(12)=Y(1)
            Y(13)=Y(1)
            Y(14)=Y(1)
            Y(15)=Y(1)
            Y(16)=Y(1)
            Y(17)=Y(1)
            DO K=1,I17
              Y0(K,NL2,NR2,M2)=Y(18-K)*X(18-K)
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
          X(4)=H3(NR3,M3,2)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(4)
          X(10)=X(1)
          X(11)=X(1)
          X(12)=X(1)
          X(13)=X(4)
          X(14)=X(1)
          X(15)=X(4)
          X(16)=X(1)
          X(17)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=H3(NL3,M3,2)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(5)
            Y(11)=Y(1)
            Y(12)=Y(1)
            Y(13)=Y(5)
            Y(14)=Y(5)
            Y(15)=Y(1)
            Y(16)=Y(1)
            Y(17)=Y(1)
            DO K=1,I17
              Z0(K,NL3,NR3,M3)=Y(18-K)*X(18-K)
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
          X(6)=H4(NR4,M4,2)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(1)
          X(10)=X(1)
          X(11)=X(6)
          X(12)=X(1)
          X(13)=X(1)
          X(14)=X(6)
          X(15)=X(1)
          X(16)=X(6)
          X(17)=X(1)
          DO NL4=1,JCI4
            Y(1)=H4(NL4,M4,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=H4(NL4,M4,2)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            Y(11)=Y(1)
            Y(12)=Y(7)
            Y(13)=Y(1)
            Y(14)=Y(1)
            Y(15)=Y(7)
            Y(16)=Y(7)
            Y(17)=Y(1)
            DO K=1,I17
              W0(K,NL4,NR4,M4)=Y(18-K)*X(18-K)
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
                    DO K=1,I17
                      TEMP3(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=0.D0
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
C**START 3-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO IRHS4=1,JCI4
            DO IRHS3=1,JCI3
              DO ILHS4=1,JCI4
                DO ILHS3=1,JCI3
                  DO K=1,I17
                    TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=0.D0
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 2-MODE INTEGRATION
CCCC      DO M3=1,MM3/MD3
          DO M3=1,MM3
            DO IRHS4=1,JCI4
              DO ILHS4=1,JCI4
                DO K=1,I17
                  TEMP1(K,ILHS4,IRHS4)=0.D0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
CCCC        DO M4=1,MM4/MD4
            DO M4=1,MM4
              DO I=1,17
                C(I)=0
              END DO
              IF(JCOUPL.GT.0)THEN
                IF(ICOUPC.GT.3)THEN
                  C(1)=VC(M4,M3,M2,M1,5)*IFACTC
                  C(2)=VC(M4,M3,M2,M1,6)*IFACTC
                  C(3)=VC(M4,M3,M2,M1,6)*IFACTC
                  C(4)=VC(M4,M3,M2,M1,7)*IFACTC
                  C(5)=VC(M4,M3,M2,M1,7)*IFACTC
                  C(6)=VC(M4,M3,M2,M1,8)*IFACTC
                  C(7)=VC(M4,M3,M2,M1,8)*IFACTC
                  C(8)=VC(M4,M3,M2,M1,9)*IFACTC
                  C(9)=VC(M4,M3,M2,M1,10)*IFACTC
                  C(10)=VC(M4,M3,M2,M1,10)*IFACTC
                  C(11)=VC(M4,M3,M2,M1,11)*IFACTC
                  C(12)=VC(M4,M3,M2,M1,11)*IFACTC
                  C(13)=VC(M4,M3,M2,M1,12)*IFACTC
                  C(14)=VC(M4,M3,M2,M1,13)*IFACTC
                  C(15)=VC(M4,M3,M2,M1,13)*IFACTC
                  C(16)=VC(M4,M3,M2,M1,14)*IFACTC
                  C(17)=VC(M4,M3,M2,M1,15)*IFACTC
                  IF(J21.GT.1)
     1            C(17)=C(17)+VR(KROT,M4,M3,M2,M1)*IFACTC
                END IF
                IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(17)=C(17)+
     1          VP(M4,M3,M2,M1)*IFACTL
              ELSE
                IF(ICOUPC.GT.3)THEN
                  C(1)=VCR(M4,M3,M2,M1,5)*IFACTC
                  C(2)=VCR(M4,M3,M2,M1,6)*IFACTC
                  C(3)=VCR(M4,M3,M2,M1,6)*IFACTC
                  C(4)=VCR(M4,M3,M2,M1,7)*IFACTC
                  C(5)=VCR(M4,M3,M2,M1,7)*IFACTC
                  C(6)=VCR(M4,M3,M2,M1,8)*IFACTC
                  C(7)=VCR(M4,M3,M2,M1,8)*IFACTC
                  C(8)=VCR(M4,M3,M2,M1,9)*IFACTC
                  C(9)=VCR(M4,M3,M2,M1,10)*IFACTC
                  C(10)=VCR(M4,M3,M2,M1,10)*IFACTC
                  C(11)=VCR(M4,M3,M2,M1,11)*IFACTC
                  C(12)=VCR(M4,M3,M2,M1,11)*IFACTC
                  C(13)=VCR(M4,M3,M2,M1,12)*IFACTC
                  C(14)=VCR(M4,M3,M2,M1,13)*IFACTC
                  C(15)=VCR(M4,M3,M2,M1,13)*IFACTC
                  C(16)=VCR(M4,M3,M2,M1,14)*IFACTC
                  C(17)=VCR(M4,M3,M2,M1,15)*IFACTC
                  IF(J21.GT.1)
     1            C(17)=C(17)+VRR(KROT,M4,M3,M2,M1)*IFACTC
                END IF
                IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(17)=C(17)+
     1          VPR(M4,M3,M2,M1)*IFACTL
              END IF
              DO IRHS4=1,JCI4
                DO ILHS4=1,JCI4
                  DO K=1,I17
                    TEMP1(K,ILHS4,IRHS4)=TEMP1(K,ILHS4,IRHS4)+
     1              W0(K,ILHS4,IRHS4,M4)*C(18-K)
                  END DO
                END DO
              END DO
            END DO
C**END 1-MODE INTEGRATION
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO K=1,I17
                      TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                Z0(K,ILHS3,IRHS3,M3)*TEMP1(K,ILHS4,IRHS4)
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
                      DO K=1,I17
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

C**NSIZE IS NO. UNIQUE INTEGRALS (4-DIM)
        DO IRHS=1,NSIZE
          NR1=IP4(IRHS,1)
          NR2=IP4(IRHS,2)
          NR3=IP4(IRHS,3)
          NR4=IP4(IRHS,4)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP4(ILHS,1)
            NL2=IP4(ILHS,2)
            NL3=IP4(ILHS,3)
            NL4=IP4(ILHS,4)
            DO K=1,I17
              XA4(ILHS+J0)=XA4(ILHS+J0)+
     1        TEMP3(K,NL2,NL3,NL4,NR2,NR3,NR4)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 4-MODE INTEGRATION
      CALL MATOUT(XA4,XRA4,NSIZE,24)
      GO TO 7777
6666  CALL MATIN(XA4,XRA4,NSIZE,24)
7777  CONTINUE

C     IF(MOLINC.GT.0)RETURN

C***********************************ALGORITHM FROM VCI4

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA)
      IOFF=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFF=IOFF+NCSIZE(I)
      END DO
      J0R=IOFF*(IOFF+1)/2
C**TEMPORARY-LAN
      JS=1
C**TEMPORARY-LAN
C     DO 9999 ISM1=1,NVSYM
      DO 9999 ISM1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR1=IPC(IROFF,MODE1)
        NR2=IPC(IROFF,MODE2)
        NR3=IPC(IROFF,MODE3)
        NR4=IPC(IROFF,MODE4)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP4(IR,1).AND.NR2.EQ.IP4(IR,2).AND.
     1       NR3.EQ.IP4(IR,3).AND.NR4.EQ.IP4(IR,4))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1      MODE4.AND.(IPC(IROFF,K).NE.IPC(ILOFF,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPC(ILOFF,MODE1)
          NL2=IPC(ILOFF,MODE2)
          NL3=IPC(ILOFF,MODE3)
          NL4=IPC(ILOFF,MODE4)
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
          XYZ=XA4(I)
          ZYX=0
          IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
            DO I=1,NP4(IND)
              K=JP4(I,IND,1)+1
              L=JP4(I,IND,2)+1
              N=JP4(I,IND,3)+1
              M=JP4(I,IND,4)+1
              ZYX=ZYX+CP4(I,IND)*XKAN(NL1,NR1,K,MOD1)*
     1        XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2        XKAN(NL4,NR4,M,MOD4)
            END DO
C**ANALYTIC
          END IF
3000      CONTINUE
          XK(ILHS,IRHS)=(XYZ+ZYX)*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
C**TEMPORARY-LAN
        J0=IROFF*(IROFF-1)/2+IOFF-J0R
C**TEMPORARY-LAN
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IP(IROFF,KCONT)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IP(IROFF,MCONT).NE.IP(ILOFF,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IP(ILOFF,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
99999 CONTINUE

      IF(ITIM4A.EQ.1)THEN
        CALL TIMIT(3)
        ITIM4A=2
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI4B(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,NN2,MM2,NN3,MM3,
     2NN4,MM4,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,
     3KCONT2,XA,ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK1,XK2,IRL,KTEMP,
     4XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,
     5JCI4,X0,Y0,Z0,W0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,
     6CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I17,XKAN,
     7MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,MAX4,INDK,INDL,INDN,MAXVAL,MAXQ4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,15),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,15),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP4(ISIZE4,4),IRL(MAXQ4,MAXQ4,MAXQ4,MAXQ4)
      DIMENSION TEMP1(I17,JCI4,JCI4),TEMP2(I17,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(I17,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(I17,JCI1,JCI1,MM1),Y0(I17,JCI2,JCI2,MM2)
      DIMENSION Z0(I17,JCI3,JCI3,MM3),W0(I17,JCI4,JCI4,MM4)
      DIMENSION X(17),Y(17),C(17)
C****************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA4(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C****************************************************************
C**ANALYTIC
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),JP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI4B(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,MM1/MD1,NN2,MM2,
     2MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,IP,ISIZMX,IPC1,ICSIZ1,
     3IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP4,ISIZE4,XA4,
     4XRA4,NSIZE,XK1,XK2,IRL,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,
     5TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VP,VPR,VC,VCR,
     6VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,
     7NVSYMX,IND,NONZC1,NONZC2,I17,XKAN,MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,
     8MAX4,INDK,INDL,INDN,MAXVAL,MAXQ4)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI4B(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MH1,MM1,NN2,MH2,MM2,
     2NN3,MH3,MM3,NN4,MH4,MM4,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,
     3IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK1,XK2,
     4IRL,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,
     5JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     6KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I17,
     7XKAN,MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,MAX4,INDK,INDL,INDN,MAXVAL,
     8MAXQ4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,15),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,15),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA4(1)

C**TEMPORARY
      DIMENSION IRL(MAXQ4,MAXQ4,MAXQ4,MAXQ4),NX4(4),MX4(4)
C**TEMPORARY

      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3),H4(NN4,MH4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP4(ISIZE4,4)
      DIMENSION TEMP1(I17,JCI4,JCI4),TEMP2(I17,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(I17,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(I17,JCI1,JCI1,MM1),Y0(I17,JCI2,JCI2,MM2)
      DIMENSION Z0(I17,JCI3,JCI3,MM3),W0(I17,JCI4,JCI4,MM4)
      DIMENSION X(17),Y(17),C(17)
C****************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA4(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C****************************************************************
C**ANALYTIC
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),JP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM4A.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCI4B'
          CALL FLUSH(IOUT)
          CALL TIMIT(1)
        END IF
      END IF

C**************************************************************
C**************************************************************
C**BOTH CONTRACTION SCHEMES HAVE AT LEAST ONE MODE
C**CHECK SCHEME '1'
      I1=KCONT1
      NC1=0
      DO NN=1,ICONT(I1)
        IF(MOD1.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE1
          JC1(NC1)=1
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD2.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE2
          JC1(NC1)=2
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD3.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE3
          JC1(NC1)=3
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD4.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE4
          JC1(NC1)=4
        END IF
      END DO
C**CHECK SCHEME '2'
      I2=KCONT2
      NC2=0
      DO NN=1,ICONT(I2)
        IF(MOD1.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE1
          JC2(NC2)=1
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD2.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE2
          JC2(NC2)=2
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD3.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE3
          JC2(NC2)=3
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD4.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE4
          JC2(NC2)=4
        END IF
      END DO
C**NKMOD1 IS TOTAL NUMBER OF MODES FOR SCHEME INVOLVING 'MOD1' (K)
C**NKMOD2 IS TOTAL NUMBER OF MODES FOR THE OTHER SCHEME (ONE OR ALL
C**OF 'MOD2' (L), 'MOD3' (N) AND 'MOD4' (M))
C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L),
C**'MOD3' (N) OR 'MOD4' (M)
C**IF THEY ARE IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
C**************************************************************
      IF(IND.NE.0)GO TO 7777
C**************************************************************

C**FIND MAX QUANTUM
      MAXQ4=0
      DO IR=1,NSIZE
        NR1=IP4(IR,1)
        NR2=IP4(IR,2)
        NR3=IP4(IR,3)
        NR4=IP4(IR,4)
        IF(NR1.GT.MAXQ4)MAXQ4=NR1
        IF(NR2.GT.MAXQ4)MAXQ4=NR2
        IF(NR3.GT.MAXQ4)MAXQ4=NR3
        IF(NR4.GT.MAXQ4)MAXQ4=NR4
      END DO

C**ANALYTIC
      IND4=INDN(MOD1)+INDL(MOD2)+INDK(MOD3)+MOD4
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI4
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,4)
      IFACTL=INTFAC(NMODE,ICOUPL,4)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(74)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(74)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.3)READ(64)VR
        IF(ICOUPC.GT.3)READ(84)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.3)READ(64)VRR
        IF(ICOUPC.GT.3)READ(84)VCR
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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(2)
          X(10)=X(2)
          X(11)=X(2)
          X(12)=X(2)
          X(13)=X(2)
          X(14)=X(2)
          X(15)=X(2)
          X(16)=X(2)
          X(17)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(1)
            Y(5)=Y(3)
            Y(6)=Y(1)
            Y(7)=Y(3)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(3)
            Y(11)=Y(3)
            Y(12)=Y(3)
            Y(13)=Y(3)
            Y(14)=Y(3)
            Y(15)=Y(3)
            Y(16)=Y(3)
            Y(17)=Y(3)
            DO K=1,I17
              X0(K,NL1,NR1,M1)=Y(18-K)*X(18-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(1)
          X(10)=X(2)
          X(11)=X(1)
          X(12)=X(2)
          X(13)=X(1)
          X(14)=X(1)
          X(15)=X(1)
          X(16)=X(1)
          X(17)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(1)
            Y(11)=Y(3)
            Y(12)=Y(1)
            Y(13)=Y(1)
            Y(14)=Y(1)
            Y(15)=Y(1)
            Y(16)=Y(1)
            Y(17)=Y(1)
            DO K=1,I17
              Y0(K,NL2,NR2,M2)=Y(18-K)*X(18-K)
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
          X(4)=H3(NR3,M3,2)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(4)
          X(10)=X(1)
          X(11)=X(1)
          X(12)=X(1)
          X(13)=X(4)
          X(14)=X(1)
          X(15)=X(4)
          X(16)=X(1)
          X(17)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=H3(NL3,M3,2)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(5)
            Y(11)=Y(1)
            Y(12)=Y(1)
            Y(13)=Y(5)
            Y(14)=Y(5)
            Y(15)=Y(1)
            Y(16)=Y(1)
            Y(17)=Y(1)
            DO K=1,I17
              Z0(K,NL3,NR3,M3)=Y(18-K)*X(18-K)
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
          X(6)=H4(NR4,M4,2)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(1)
          X(10)=X(1)
          X(11)=X(6)
          X(12)=X(1)
          X(13)=X(1)
          X(14)=X(6)
          X(15)=X(1)
          X(16)=X(6)
          X(17)=X(1)
          DO NL4=1,JCI4
            Y(1)=H4(NL4,M4,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=H4(NL4,M4,2)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            Y(11)=Y(1)
            Y(12)=Y(7)
            Y(13)=Y(1)
            Y(14)=Y(1)
            Y(15)=Y(7)
            Y(16)=Y(7)
            Y(17)=Y(1)
            DO K=1,I17
              W0(K,NL4,NR4,M4)=Y(18-K)*X(18-K)
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
                    DO K=1,I17
                      TEMP3(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=0
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
C**START 3-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO IRHS4=1,JCI4
            DO IRHS3=1,JCI3
              DO ILHS4=1,JCI4
                DO ILHS3=1,JCI3
                  DO K=1,I17
                    TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=0
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 2-MODE INTEGRATION
CCCC      DO M3=1,MM3/MD3
          DO M3=1,MM3
            DO IRHS4=1,JCI4
              DO ILHS4=1,JCI4
                DO K=1,I17
                  TEMP1(K,ILHS4,IRHS4)=0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
CCCC        DO M4=1,MM4/MD4
            DO M4=1,MM4
              DO I=1,17
                C(I)=0
              END DO
              IF(JCOUPL.GT.0)THEN
                IF(ICOUPC.GT.3)THEN
                  C(1)=VC(M4,M3,M2,M1,5)*IFACTC
                  C(2)=VC(M4,M3,M2,M1,6)*IFACTC
                  C(3)=VC(M4,M3,M2,M1,6)*IFACTC
                  C(4)=VC(M4,M3,M2,M1,7)*IFACTC
                  C(5)=VC(M4,M3,M2,M1,7)*IFACTC
                  C(6)=VC(M4,M3,M2,M1,8)*IFACTC
                  C(7)=VC(M4,M3,M2,M1,8)*IFACTC
                  C(8)=VC(M4,M3,M2,M1,9)*IFACTC
                  C(9)=VC(M4,M3,M2,M1,10)*IFACTC
                  C(10)=VC(M4,M3,M2,M1,10)*IFACTC
                  C(11)=VC(M4,M3,M2,M1,11)*IFACTC
                  C(12)=VC(M4,M3,M2,M1,11)*IFACTC
                  C(13)=VC(M4,M3,M2,M1,12)*IFACTC
                  C(14)=VC(M4,M3,M2,M1,13)*IFACTC
                  C(15)=VC(M4,M3,M2,M1,13)*IFACTC
                  C(16)=VC(M4,M3,M2,M1,14)*IFACTC
                  C(17)=VC(M4,M3,M2,M1,15)*IFACTC
                  IF(J21.GT.1)
     1            C(17)=C(17)+VR(KROT,M4,M3,M2,M1)*IFACTC
                END IF
                IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(17)=C(17)+
     1          VP(M4,M3,M2,M1)*IFACTL
              ELSE
                IF(ICOUPC.GT.3)THEN
                  C(1)=VCR(M4,M3,M2,M1,5)*IFACTC
                  C(2)=VCR(M4,M3,M2,M1,6)*IFACTC
                  C(3)=VCR(M4,M3,M2,M1,6)*IFACTC
                  C(4)=VCR(M4,M3,M2,M1,7)*IFACTC
                  C(5)=VCR(M4,M3,M2,M1,7)*IFACTC
                  C(6)=VCR(M4,M3,M2,M1,8)*IFACTC
                  C(7)=VCR(M4,M3,M2,M1,8)*IFACTC
                  C(8)=VCR(M4,M3,M2,M1,9)*IFACTC
                  C(9)=VCR(M4,M3,M2,M1,10)*IFACTC
                  C(10)=VCR(M4,M3,M2,M1,10)*IFACTC
                  C(11)=VCR(M4,M3,M2,M1,11)*IFACTC
                  C(12)=VCR(M4,M3,M2,M1,11)*IFACTC
                  C(13)=VCR(M4,M3,M2,M1,12)*IFACTC
                  C(14)=VCR(M4,M3,M2,M1,13)*IFACTC
                  C(15)=VCR(M4,M3,M2,M1,13)*IFACTC
                  C(16)=VCR(M4,M3,M2,M1,14)*IFACTC
                  C(17)=VCR(M4,M3,M2,M1,15)*IFACTC
                  IF(J21.GT.1)
     1            C(17)=C(17)+VRR(KROT,M4,M3,M2,M1)*IFACTC
                END IF
                IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(17)=C(17)+
     1          VPR(M4,M3,M2,M1)*IFACTL
              END IF
              DO IRHS4=1,JCI4
                DO ILHS4=1,JCI4
                  DO K=1,I17
                    TEMP1(K,ILHS4,IRHS4)=TEMP1(K,ILHS4,IRHS4)+
     1              W0(K,ILHS4,IRHS4,M4)*C(18-K)
                  END DO
                END DO
              END DO
            END DO
C**END 1-MODE INTEGRATION
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO K=1,I17
                      TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                Z0(K,ILHS3,IRHS3,M3)*TEMP1(K,ILHS4,IRHS4)
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
                      DO K=1,I17
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

C**NSIZE IS NO. UNIQUE INTEGRALS (4-DIM)
        DO IRHS=1,NSIZE
          NR1=IP4(IRHS,1)
          NR2=IP4(IRHS,2)
          NR3=IP4(IRHS,3)
          NR4=IP4(IRHS,4)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP4(ILHS,1)
            NL2=IP4(ILHS,2)
            NL3=IP4(ILHS,3)
            NL4=IP4(ILHS,4)
            DO K=1,I17
              XA4(ILHS+J0)=XA4(ILHS+J0)+
     1        TEMP3(K,NL2,NL3,NL4,NR2,NR3,NR4)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 4-MODE INTEGRATION
      CALL MATOUT(XA4,XRA4,NSIZE,24)
      GO TO 7777
6666  CALL MATIN(XA4,XRA4,NSIZE,24)
7777  CONTINUE

      IF(IND.NE.0)THEN
        DO IR=1,NSIZE
          NR=0
          DO J=1,NC1
            J1=JC1(J)
            NR1=IP4(IR,J1)
            NR=NR+1
            NX4(NR)=NR1
          END DO
          DO J=1,NC2
            J2=JC2(J)
            NR2=IP4(IR,J2)
            NR=NR+1
            NX4(NR)=NR2
          END DO
          IRL(NX4(1),NX4(2),NX4(3),NX4(4))=IR
        END DO
      END IF

C***********************************ALGORITHM FROM VCI4

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFFR=IOFFR+NCSIZE(I)
      END DO
      J0R=IOFFR*(IOFFR+1)/2
C**TEMPORARY-LAN
      JSR=1
C**TEMPORARY-LAN
C     DO 9999 ISMR1=1,NVSYM
      DO 9999 ISMR1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE

        IF(NONC1.GT.NONC2)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          KCCCL=KCSZL
CC        IF(NONC1.GT.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             KCCCL=KCSZL
C           ELSE
C             KCCCL=IRH1
C           END IF
CC        END IF
          DO ILH1=1,KCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=IROFF
            NONZC1(2,NONZ1)=ILOFF
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
        IF(NONZ1.EQ.0)GO TO 5556
        NONZC1(3,NONZ1+1)=0
        END IF

      END IF

CCCC  IF(NONZ1.EQ.0)GO TO 5556
CCCC  NONZC1(3,NONZ1+1)=0

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE

        IF(NONC1.LE.NONC2)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          LCCCL=LCSZL
CC        IF(NONC1.LE.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             LCCCL=LCSZL
C           ELSE
C             LCCCL=IRH2
C           END IF
CC        END IF
          DO ILH2=1,LCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=JROFF
            NONZC2(2,NONZ2)=JLOFF
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
        IF(NONZ2.EQ.0)GO TO 5555
        NONZC2(3,NONZ2+1)=0
        END IF

      END IF

CCCC  IF(NONZ2.EQ.0)GO TO 5555
CCCC  NONZC2(3,NONZ2+1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM4A.EQ.0)THEN
C       WRITE(IOUT,*)'VCCI4B - ISMR1,ISML1 = ',ISMR1,ISML1
C       WRITE(IOUT,*)'ICONDP = ',ICONDP
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C       WRITE(IOUT,*)'NONC1,NONC2 = ',NONC1,NONC2
C     END IF

      IF(NONC1.LE.NONC2)GO TO 3003
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR/ M)

      DO IRH2=1,LCSZR
        JROFF=IRH2+LCOFFR
        DO I=1,NKVALR
          DO J=1,NKVALL
            DO K=1,NLVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          LCCCL=LCSZL
        ELSE
          LCCCL=IRH2
        END IF
        DO ILH2=1,LCCCL
          JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
          DO K=1,NKMOD2
            IGOT=0
            DO I=1,NC2
              IF(K.EQ.MC2(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 6001
          END DO
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH2.EQ.IRH2)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NKVALR
            DO J=1,NKVALL
              XCON2(J,I)=0
            END DO
          END DO
C
C**CONTRACTION SCHEME  '1' (K)
C
C         DO IRH1=1,KCSZR
C           IROFF=IRH1+KCOFFR
C           DO J=1,NKVALL
C             XK1(J,1)=0
C           END DO
C           DO ILH1=1,KCSZL
C             ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
C             DO K=1,NKMOD1
C               IGOT=0
C               DO I=1,NC1
C                 IF(K.EQ.MC1(I))IGOT=1
C               END DO
C               IF(IGOT.EQ.0.AND.
C    1          (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 5001
C             END DO

          IPREV=0
          DO INH1=1,NONZ1
              IRH1=NONZC1(3,INH1)
              INEXT=NONZC1(3,INH1+1)
              ILH1=NONZC1(4,INH1)
              IROFF=NONZC1(1,INH1)
              ILOFF=NONZC1(2,INH1)
              IF(IPREV.NE.IRH1)THEN
C**NEW RHS
                DO J=1,NKVALL
                  XK1(J,1)=0
                END DO
              END IF

              DO I=1,NC1
                NCR1(I)=IPC1(IROFF,MC1(I))
                NCL1(I)=IPC1(ILOFF,MC1(I))
              END DO

C**FIND RHS INDEX
C             DO IR=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCR1(I).NE.IP4(IR,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCR2(I).NE.IP4(IR,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 1000
C             END DO
1000          CONTINUE
C**FIND LHS INDEX
C             DO IL=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCL1(I).NE.IP4(IL,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCL2(I).NE.IP4(IL,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 2000
C             END DO
2000          CONTINUE

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX4(NR)=NR1
        MX4(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX4(NR)=NR2
        MX4(NR)=NL2
      END DO
      IR=IRL(NX4(1),NX4(2),NX4(3),NX4(4))
      IL=IRL(MX4(1),MX4(2),MX4(3),MX4(4))
C**TEMPORARY

C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
C             XYZ=XA4(I)
C             ZYX=0
C             IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
C               NR1=IP4(MR,1)
C               NL1=IP4(ML,1)
C               NR2=IP4(MR,2)
C               NL2=IP4(ML,2)
C               NR3=IP4(MR,3)
C               NL3=IP4(ML,3)
C               NR4=IP4(MR,4)
C               NL4=IP4(ML,4)
C               DO I=1,NP4(IND4)
C                 K=JP4(I,IND4,1)+1
C                 L=JP4(I,IND4,2)+1
C                 N=JP4(I,IND4,3)+1
C                 M=JP4(I,IND4,4)+1
C                 ZYX=ZYX+CP4(I,IND4)*XKAN(NL1,NR1,K,MOD1)*
C    1            XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
C    2            XKAN(NL4,NR4,M,MOD4)
C               END DO
C**ANALYTIC
C             END IF
C             XXX=XYZ+ZYX
              XXX=XA4(I)

              DO K=1,NKVALL
                XK1(K,1)=XK1(K,1)+XXX*CFS1(ILH1,K,KEL,KSML)
              END DO
5001          CONTINUE

C           END DO
            IF(IRH1.NE.INEXT)THEN

            DO K=1,NKVALR
              DO L=1,NKVALL
                XCON2(L,K)=XCON2(L,K)+XK1(L,1)*CFS1(IRH1,K,KEL,KSMR)
              END DO
            END DO

            END IF
            IPREV=IRH1
          END DO
CC
          DO I=1,NKVALR
            DO J=1,NKVALL
              DO K=1,NLVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6001      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILL,IKR,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKR,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILR,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILR,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO
C*************************************

      GO TO 5555

3003  CONTINUE

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '1' (K)

      DO IRH1=1,KCSZR
        IROFF=IRH1+KCOFFR
        DO I=1,NLVALR
          DO J=1,NLVALL
            DO K=1,NKVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          KCCCL=KCSZL
        ELSE
          KCCCL=IRH1
        END IF
        DO ILH1=1,KCCCL
          ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
          DO K=1,NKMOD1
            IGOT=0
            DO I=1,NC1
              IF(K.EQ.MC1(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 6002
          END DO
          DO I=1,NC1
            NCR1(I)=IPC1(IROFF,MC1(I))
            NCL1(I)=IPC1(ILOFF,MC1(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH1.EQ.IRH1)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NLVALR
            DO J=1,NLVALL
              XCON2(J,I)=0
            END DO
          END DO

C**CONTRACTION SCHEME  '2' (L AND/OR N AND/OR M)

C         DO IRH2=1,LCSZR
C           JROFF=IRH2+LCOFFR
C           DO J=1,NLVALL
C             XK1(J,1)=0
C           END DO
C           DO ILH2=1,LCSZL
C             JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
C             DO K=1,NKMOD2
C               IGOT=0
C               DO I=1,NC2
C                 IF(K.EQ.MC2(I))IGOT=1
C               END DO
C               IF(IGOT.EQ.0.AND.
C    1          (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 5002
C             END DO

          JPREV=0
          DO INH2=1,NONZ2
              IRH2=NONZC2(3,INH2)
              JNEXT=NONZC2(3,INH2+1)
              ILH2=NONZC2(4,INH2)
              JROFF=NONZC2(1,INH2)
              JLOFF=NONZC2(2,INH2)
              IF(JPREV.NE.IRH2)THEN
C**NEW RHS
                DO J=1,NLVALL
                  XK1(J,1)=0
                END DO
              END IF

              DO I=1,NC2
                NCR2(I)=IPC2(JROFF,MC2(I))
                NCL2(I)=IPC2(JLOFF,MC2(I))
              END DO

C**FIND RHS INDEX
C             DO IR=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCR1(I).NE.IP4(IR,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCR2(I).NE.IP4(IR,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 1001
C             END DO
1001          CONTINUE
C**FIND LHS INDEX
C             DO IL=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCL1(I).NE.IP4(IL,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCL2(I).NE.IP4(IL,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 2001
C             END DO
2001          CONTINUE

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX4(NR)=NR1
        MX4(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX4(NR)=NR2
        MX4(NR)=NL2
      END DO
      IR=IRL(NX4(1),NX4(2),NX4(3),NX4(4))
      IL=IRL(MX4(1),MX4(2),MX4(3),MX4(4))
C**TEMPORARY

C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
C             XYZ=XA4(I)
C             ZYX=0
C             IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
C               NR1=IP4(MR,1)
C               NL1=IP4(ML,1)
C               NR2=IP4(MR,2)
C               NL2=IP4(ML,2)
C               NR3=IP4(MR,3)
C               NL3=IP4(ML,3)
C               NR4=IP4(MR,4)
C               NL4=IP4(ML,4)
C               DO I=1,NP4(IND4)
C                 K=JP4(I,IND4,1)+1
C                 L=JP4(I,IND4,2)+1
C                 N=JP4(I,IND4,3)+1
C                 M=JP4(I,IND4,4)+1
C                 ZYX=ZYX+CP4(I,IND4)*XKAN(NL1,NR1,K,MOD1)*
C    1            XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
C    2            XKAN(NL4,NR4,M,MOD4)
C               END DO
C**ANALYTIC
C             END IF
C             XK1(ILH2,IRH2)=XYZ+ZYX
C             XXX=XYZ+ZYX
              XXX=XA4(I)

              DO J=1,NLVALL
                XK1(J,1)=XK1(J,1)+XXX*CFS2(ILH2,J,KEL,LSML)
              END DO
5002          CONTINUE

C           END DO
            IF(IRH2.NE.JNEXT)THEN

            DO I=1,NLVALR
              DO J=1,NLVALL
                XCON2(J,I)=XCON2(J,I)+XK1(J,1)*CFS2(IRH2,I,KEL,LSMR)
              END DO
            END DO

            END IF
            JPREV=IRH2
          END DO
CC
          DO I=1,NLVALR
            DO J=1,NLVALL
              DO K=1,NKVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6002      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKR,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKL,ILR,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILR,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO

5555  CONTINUE
5556  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM4A.EQ.0)THEN
        IF(IND.EQ.1)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM4A=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI4C(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,NN2,MM2,NN3,MM3,
     2NN4,MM4,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP4,ISIZE4,
     3XA4,XRA4,NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,
     4TEMP3,JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VP,VPR,VC,VCR,VR,VRR,J21,
     5KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I17,XKAN,MAXQU,
     6MAXPOW,NP4,CP4,JP4,NTOT4,MAX4,INDK,INDL,INDN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,15),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,15),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE),IRL(ICSIZE)
      DIMENSION IPC(IPSIZE,1),IP4(ISIZE4,4),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I17,JCI4,JCI4),TEMP2(I17,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(I17,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(I17,JCI1,JCI1,MM1),Y0(I17,JCI2,JCI2,MM2)
      DIMENSION Z0(I17,JCI3,JCI3,MM3),W0(I17,JCI4,JCI4,MM4)
      DIMENSION X(17),Y(17),C(17)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA4(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),JP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI4C(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,MM1/MD1,NN2,MM2,
     2MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,IP,ISIZMX,IPC,ICSIZE,
     3IPSIZE,KCONT,XA,ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK,TEMP,IRL,XCON,
     4NVAL,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,
     5W0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,
     6KEL21,NVSYMX,I17,XKAN,MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,MAX4,INDK,
     7INDL,INDN)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI4C(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MH1,MM1,NN2,MH2,MM2,
     2NN3,MH3,MM3,NN4,MH4,MM4,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,
     3ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,
     4IEND,TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VP,VPR,VC,
     5VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,
     6I17,XKAN,MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,MAX4,INDK,INDL,INDN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,15),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,15),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3),H4(NN4,MH4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE),IRL(ICSIZE)
      DIMENSION IPC(IPSIZE,1),IP4(ISIZE4,4),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(I17,JCI4,JCI4),TEMP2(I17,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(I17,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(I17,JCI1,JCI1,MM1),Y0(I17,JCI2,JCI2,MM2)
      DIMENSION Z0(I17,JCI3,JCI3,MM3),W0(I17,JCI4,JCI4,MM4)
      DIMENSION X(17),Y(17),C(17)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA4(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),JP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM4A.EQ.1)THEN
        WRITE(IOUT,*)'Calculating VCCI4C'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

C** FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1', 'MOD2',
C** 'MOD3' AND 'MOD4'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      NLMODE=ICONT(2)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        NLMODE=ICONT(1)
        MCONT=1
      END IF

C**ANALYTIC
      IND=INDN(MOD1)+INDL(MOD2)+INDK(MOD3)+MOD4
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI4
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,4)
      IFACTL=INTFAC(NMODE,ICOUPL,4)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(74)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(74)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.3)READ(64)VR
        IF(ICOUPC.GT.3)READ(84)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.3)READ(64)VRR
        IF(ICOUPC.GT.3)READ(84)VCR
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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(2)
          X(10)=X(2)
          X(11)=X(2)
          X(12)=X(2)
          X(13)=X(2)
          X(14)=X(2)
          X(15)=X(2)
          X(16)=X(2)
          X(17)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(1)
            Y(5)=Y(3)
            Y(6)=Y(1)
            Y(7)=Y(3)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(3)
            Y(11)=Y(3)
            Y(12)=Y(3)
            Y(13)=Y(3)
            Y(14)=Y(3)
            Y(15)=Y(3)
            Y(16)=Y(3)
            Y(17)=Y(3)
            DO K=1,I17
              X0(K,NL1,NR1,M1)=Y(18-K)*X(18-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(1)
          X(10)=X(2)
          X(11)=X(1)
          X(12)=X(2)
          X(13)=X(1)
          X(14)=X(1)
          X(15)=X(1)
          X(16)=X(1)
          X(17)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(1)
            Y(11)=Y(3)
            Y(12)=Y(1)
            Y(13)=Y(1)
            Y(14)=Y(1)
            Y(15)=Y(1)
            Y(16)=Y(1)
            Y(17)=Y(1)
            DO K=1,I17
              Y0(K,NL2,NR2,M2)=Y(18-K)*X(18-K)
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
          X(4)=H3(NR3,M3,2)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(4)
          X(10)=X(1)
          X(11)=X(1)
          X(12)=X(1)
          X(13)=X(4)
          X(14)=X(1)
          X(15)=X(4)
          X(16)=X(1)
          X(17)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=H3(NL3,M3,2)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(5)
            Y(11)=Y(1)
            Y(12)=Y(1)
            Y(13)=Y(5)
            Y(14)=Y(5)
            Y(15)=Y(1)
            Y(16)=Y(1)
            Y(17)=Y(1)
            DO K=1,I17
              Z0(K,NL3,NR3,M3)=Y(18-K)*X(18-K)
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
          X(6)=H4(NR4,M4,2)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(1)
          X(10)=X(1)
          X(11)=X(6)
          X(12)=X(1)
          X(13)=X(1)
          X(14)=X(6)
          X(15)=X(1)
          X(16)=X(6)
          X(17)=X(1)
          DO NL4=1,JCI4
            Y(1)=H4(NL4,M4,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=H4(NL4,M4,2)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            Y(11)=Y(1)
            Y(12)=Y(7)
            Y(13)=Y(1)
            Y(14)=Y(1)
            Y(15)=Y(7)
            Y(16)=Y(7)
            Y(17)=Y(1)
            DO K=1,I17
              W0(K,NL4,NR4,M4)=Y(18-K)*X(18-K)
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
                    DO K=1,I17
                      TEMP3(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=0.D0
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
C**START 3-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO IRHS4=1,JCI4
            DO IRHS3=1,JCI3
              DO ILHS4=1,JCI4
                DO ILHS3=1,JCI3
                  DO K=1,I17
                    TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=0.D0
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 2-MODE INTEGRATION
CCCC      DO M3=1,MM3/MD3
          DO M3=1,MM3
            DO IRHS4=1,JCI4
              DO ILHS4=1,JCI4
                DO K=1,I17
                  TEMP1(K,ILHS4,IRHS4)=0.D0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
CCCC        DO M4=1,MM4/MD4
            DO M4=1,MM4
              DO I=1,17
                C(I)=0
              END DO
              IF(JCOUPL.GT.0)THEN
                IF(ICOUPC.GT.3)THEN
                  C(1)=VC(M4,M3,M2,M1,5)*IFACTC
                  C(2)=VC(M4,M3,M2,M1,6)*IFACTC
                  C(3)=VC(M4,M3,M2,M1,6)*IFACTC
                  C(4)=VC(M4,M3,M2,M1,7)*IFACTC
                  C(5)=VC(M4,M3,M2,M1,7)*IFACTC
                  C(6)=VC(M4,M3,M2,M1,8)*IFACTC
                  C(7)=VC(M4,M3,M2,M1,8)*IFACTC
                  C(8)=VC(M4,M3,M2,M1,9)*IFACTC
                  C(9)=VC(M4,M3,M2,M1,10)*IFACTC
                  C(10)=VC(M4,M3,M2,M1,10)*IFACTC
                  C(11)=VC(M4,M3,M2,M1,11)*IFACTC
                  C(12)=VC(M4,M3,M2,M1,11)*IFACTC
                  C(13)=VC(M4,M3,M2,M1,12)*IFACTC
                  C(14)=VC(M4,M3,M2,M1,13)*IFACTC
                  C(15)=VC(M4,M3,M2,M1,13)*IFACTC
                  C(16)=VC(M4,M3,M2,M1,14)*IFACTC
                  C(17)=VC(M4,M3,M2,M1,15)*IFACTC
                  IF(J21.GT.1)
     1            C(17)=C(17)+VR(KROT,M4,M3,M2,M1)*IFACTC
                END IF
                IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(17)=C(17)+
     1          VP(M4,M3,M2,M1)*IFACTL
              ELSE
                IF(ICOUPC.GT.3)THEN
                  C(1)=VCR(M4,M3,M2,M1,5)*IFACTC
                  C(2)=VCR(M4,M3,M2,M1,6)*IFACTC
                  C(3)=VCR(M4,M3,M2,M1,6)*IFACTC
                  C(4)=VCR(M4,M3,M2,M1,7)*IFACTC
                  C(5)=VCR(M4,M3,M2,M1,7)*IFACTC
                  C(6)=VCR(M4,M3,M2,M1,8)*IFACTC
                  C(7)=VCR(M4,M3,M2,M1,8)*IFACTC
                  C(8)=VCR(M4,M3,M2,M1,9)*IFACTC
                  C(9)=VCR(M4,M3,M2,M1,10)*IFACTC
                  C(10)=VCR(M4,M3,M2,M1,10)*IFACTC
                  C(11)=VCR(M4,M3,M2,M1,11)*IFACTC
                  C(12)=VCR(M4,M3,M2,M1,11)*IFACTC
                  C(13)=VCR(M4,M3,M2,M1,12)*IFACTC
                  C(14)=VCR(M4,M3,M2,M1,13)*IFACTC
                  C(15)=VCR(M4,M3,M2,M1,13)*IFACTC
                  C(16)=VCR(M4,M3,M2,M1,14)*IFACTC
                  C(17)=VCR(M4,M3,M2,M1,15)*IFACTC
                  IF(J21.GT.1)
     1            C(17)=C(17)+VRR(KROT,M4,M3,M2,M1)*IFACTC
                END IF
                IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(17)=C(17)+
     1          VPR(M4,M3,M2,M1)*IFACTL
              END IF
              DO IRHS4=1,JCI4
                DO ILHS4=1,JCI4
                  DO K=1,I17
                    TEMP1(K,ILHS4,IRHS4)=TEMP1(K,ILHS4,IRHS4)+
     1              W0(K,ILHS4,IRHS4,M4)*C(18-K)
                  END DO
                END DO
              END DO
            END DO
C**END 1-MODE INTEGRATION
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO K=1,I17
                      TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                Z0(K,ILHS3,IRHS3,M3)*TEMP1(K,ILHS4,IRHS4)
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
                      DO K=1,I17
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

C**NSIZE IS NO. UNIQUE INTEGRALS (4-DIM)
        DO IRHS=1,NSIZE
          NR1=IP4(IRHS,1)
          NR2=IP4(IRHS,2)
          NR3=IP4(IRHS,3)
          NR4=IP4(IRHS,4)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP4(ILHS,1)
            NL2=IP4(ILHS,2)
            NL3=IP4(ILHS,3)
            NL4=IP4(ILHS,4)
            DO K=1,I17
              XA4(ILHS+J0)=XA4(ILHS+J0)+
     1        TEMP3(K,NL2,NL3,NL4,NR2,NR3,NR4)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 4-MODE INTEGRATION
      CALL MATOUT(XA4,XRA4,NSIZE,24)
      GO TO 7777
6666  CALL MATIN(XA4,XRA4,NSIZE,24)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI4

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        KCOFF=ICOFF1
        LCOFF=ICOFF2
        KCSZ=ICSZ1
        LCSZ=ICSZ2
        ICNTZ=ICCONT(1)
        NKVAL=NCVAL(1,ISM1)
        NLVAL=NCVAL(2,ISM2)
      ELSE
        ISM=ISM2
        KCOFF=ICOFF2
        LCOFF=ICOFF1
        KCSZ=ICSZ2
        LCSZ=ICSZ1
        ICNTZ=ICCONT(2)
        NKVAL=NCVAL(2,ISM2)
        NLVAL=NCVAL(1,ISM1)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

C
      DO IKRL=1,KCSZ
        IOFFRL=IKRL+KCOFF
        NR1=IPC(IOFFRL,MODE1)
        NR2=IPC(IOFFRL,MODE2)
        NR3=IPC(IOFFRL,MODE3)
        NR4=IPC(IOFFRL,MODE4)
        DO IR=1,NSIZE
          IF(NR1.EQ.IP4(IR,1).AND.NR2.EQ.IP4(IR,2).AND.
     1    NR3.EQ.IP4(IR,3).AND.NR4.EQ.IP4(IR,4))GO TO 1000
        END DO
1000    IRL(IKRL)=IR
      END DO
C
C**CASE 1..INTEGRATION MODES IN CONTRACTED SCHEME.
C**      ..OVERLAP OF NON-CONTRACTED BASIS.
      IF(ICNTZ.GT.0)THEN
        DO IRHS=1,KCSZ
          IOFFR=IRHS+KCOFF
          IR=IRL(IRHS)
          DO ILHS=1,IRHS
            IOFFL=ILHS+KCOFF
C**OVERLAP OF REMAINING STATES
            IS=1
            DO K=1,NKMODE
              IF(IS.EQ.0)GO TO 2000
              IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1        MODE4.AND.(IPC(IOFFR,K).NE.IPC(IOFFL,K)))IS=0
            END DO
C**OVERLAP OF REMAINING STATES
            IL=IRL(ILHS)
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA4(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              DO I=1,NP4(IND)
                K=JP4(I,IND,1)+1
                L=JP4(I,IND,2)+1
                N=JP4(I,IND,3)+1
                M=JP4(I,IND,4)+1
                ZYX=ZYX+CP4(I,IND)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)
              END DO
C**ANALYTIC
            END IF
2000        CONTINUE
            XK(ILHS,IRHS)=(XYZ+ZYX)*IS
            XK(IRHS,ILHS)=XK(ILHS,IRHS)
          END DO
        END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',KCSZ,NKVAL,KCSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
        CALL DGEMM('T','N',NKVAL,NKVAL,KCSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

        IRHS=0
        DO IKR=1,NKVAL
          DO ILR=1,LCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            DO IKL=1,IKR
              ILHS=(IKL-1)*LCSZ+ILR
              XYZ=XCON(IKL,IKR)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ
            END DO
          END DO
        END DO

C**CASE 2..INTEGRATION MODES IN NON-CONTRACTED SCHEME.
C**      ..OVERLAP OF CONTRACTED BASIS.
      ELSE
        IRHS=0
        DO ILR=1,NLVAL
          DO IKR=1,KCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            IOFFR=IKR+KCOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            IR=IRL(IKR)
            DO IKL=1,IKR
              ILHS=(ILR-1)*KCSZ+IKL
              IOFFL=IKL+KCOFF
C**OVERLAP OF REMAINING STATES
              IS=1
              DO K=1,NKMODE
                IF(IS.EQ.0)GO TO 3000
                IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.
     1          K.NE.MODE4.AND.(IPC(IOFFR,K).NE.IPC(IOFFL,K)))IS=0
              END DO
C**OVERLAP OF REMAINING STATES
              IL=IRL(IKL)
C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
              XYZ=XA4(I)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ*IS
3000          CONTINUE
            END DO
          END DO
        END DO
      END IF

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
      IF(ITIM4A.EQ.1)THEN
        CALL TIMIT(3)
        ITIM4A=2
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI4D(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,NN2,MM2,NN3,MM3,
     2NN4,MM4,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,
     3KCONT2,XA,ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK1,XK2,TEMP,KTEMP,
     4XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,
     5JCI4,X0,Y0,Z0,W0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,
     6CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I17,XKAN,
     7MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,MAX4,INDK,INDL,INDN,MAXVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,15),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,15),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP4(ISIZE4,4),TEMP(KTEMP,1)
      DIMENSION TEMP1(I17,JCI4,JCI4),TEMP2(I17,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(I17,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(I17,JCI1,JCI1,MM1),Y0(I17,JCI2,JCI2,MM2)
      DIMENSION Z0(I17,JCI3,JCI3,MM3),W0(I17,JCI4,JCI4,MM4)
      DIMENSION X(17),Y(17),C(17)
C****************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA4(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,1)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C****************************************************************
C**ANALYTIC
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),JP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      CALL VDCI4D(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,MM1/MD1,NN2,MM2,
     2MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,IP,ISIZMX,IPC1,ICSIZ1,
     3IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP4,ISIZE4,XA4,
     4XRA4,NSIZE,XK1,XK2,TEMP,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,
     5TEMP1,TEMP2,TEMP3,JCI1,JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VP,VPR,VC,VCR,
     6VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,
     7NVSYMX,IND,NONZC1,NONZC2,I17,XKAN,MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,
     8MAX4,INDK,INDL,INDN,MAXVAL)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI4D(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MH1,MM1,NN2,MH2,MM2,
     2NN3,MH3,MM3,NN4,MH4,MM4,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,
     3IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP4,ISIZE4,XA4,XRA4,NSIZE,XK1,XK2,
     4TEMP,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,JCI1,
     5JCI2,JCI3,JCI4,X0,Y0,Z0,W0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,
     6KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I17,
     7XKAN,MAXQU,MAXPOW,NP4,CP4,JP4,NTOT4,MAX4,INDK,INDL,INDN,MAXVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,15),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,15),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA4(1)
      DIMENSION MODINT(NMODE)
CCCC  DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3),H4(NN4,MH4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP4(ISIZE4,4),TEMP(KTEMP,1)
      DIMENSION TEMP1(I17,JCI4,JCI4),TEMP2(I17,JCI3,JCI4,JCI3,JCI4)
      DIMENSION TEMP3(I17,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION X0(I17,JCI1,JCI1,MM1),Y0(I17,JCI2,JCI2,MM2)
      DIMENSION Z0(I17,JCI3,JCI3,MM3),W0(I17,JCI4,JCI4,MM4)
      DIMENSION X(17),Y(17),C(17)
C****************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA4(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,1)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C****************************************************************
C**ANALYTIC
      DIMENSION NP4(NTOT4),CP4(MAX4,NTOT4),JP4(MAX4,NTOT4,4)
      DIMENSION INDK(1),INDL(1),INDN(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM4A.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCI4D'
          CALL FLUSH(IOUT)
          CALL TIMIT(1)
        END IF
      END IF

C**************************************************************
C**************************************************************
C**BOTH CONTRACTION SCHEMES HAVE AT LEAST ONE MODE
C**CHECK SCHEME '1'
      I1=KCONT1
      NC1=0
      DO NN=1,ICONT(I1)
        IF(MOD1.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE1
          JC1(NC1)=1
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD2.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE2
          JC1(NC1)=2
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD3.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE3
          JC1(NC1)=3
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD4.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE4
          JC1(NC1)=4
        END IF
      END DO
C**CHECK SCHEME '2'
      I2=KCONT2
      NC2=0
      DO NN=1,ICONT(I2)
        IF(MOD1.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE1
          JC2(NC2)=1
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD2.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE2
          JC2(NC2)=2
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD3.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE3
          JC2(NC2)=3
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD4.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE4
          JC2(NC2)=4
        END IF
      END DO
C**NKMOD1 IS TOTAL NUMBER OF MODES FOR SCHEME INVOLVING 'MOD1' (K)
C**NKMOD2 IS TOTAL NUMBER OF MODES FOR THE OTHER SCHEME (ONE OR ALL
C**OF 'MOD2' (L), 'MOD3' (N) AND 'MOD4' (M))
C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L),
C**'MOD3' (N) OR 'MOD4' (M)
C**IF THEY ARE IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
C**************************************************************
      IF(IND.NE.0)GO TO 7777
C**************************************************************

C**ANALYTIC
      IND4=INDN(MOD1)+INDL(MOD2)+INDK(MOD3)+MOD4
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI4
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,4)
      IFACTL=INTFAC(NMODE,ICOUPL,4)
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(74)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(74)VPR
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(J21.GT.1.AND.ICOUPC.GT.3)READ(64)VR
        IF(ICOUPC.GT.3)READ(84)VC
      ELSE
        IF(J21.GT.1.AND.ICOUPC.GT.3)READ(64)VRR
        IF(ICOUPC.GT.3)READ(84)VCR
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
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
CCCC  DO M1=1,MM1/MD1
      DO M1=1,MM1
        DO NR1=1,JCI1
          X(1)=H1(NR1,M1,2)*MD
          X(2)=H1(NR1,M1,1)*MD
          X(3)=X(1)
          X(4)=X(2)
          X(5)=X(1)
          X(6)=X(2)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(2)
          X(10)=X(2)
          X(11)=X(2)
          X(12)=X(2)
          X(13)=X(2)
          X(14)=X(2)
          X(15)=X(2)
          X(16)=X(2)
          X(17)=X(2)
          DO NL1=1,JCI1
            Y(1)=H1(NL1,M1,2)
            Y(2)=Y(1)
            Y(3)=H1(NL1,M1,1)
            Y(4)=Y(1)
            Y(5)=Y(3)
            Y(6)=Y(1)
            Y(7)=Y(3)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(3)
            Y(11)=Y(3)
            Y(12)=Y(3)
            Y(13)=Y(3)
            Y(14)=Y(3)
            Y(15)=Y(3)
            Y(16)=Y(3)
            Y(17)=Y(3)
            DO K=1,I17
              X0(K,NL1,NR1,M1)=Y(18-K)*X(18-K)
            END DO
          END DO
        END DO
      END DO
CCCC  DO M2=1,MM2/MD2
      DO M2=1,MM2
        DO NR2=1,JCI2
          X(1)=H2(NR2,M2,1)
          X(2)=H2(NR2,M2,2)
          X(3)=X(1)
          X(4)=X(1)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(2)
          X(9)=X(1)
          X(10)=X(2)
          X(11)=X(1)
          X(12)=X(2)
          X(13)=X(1)
          X(14)=X(1)
          X(15)=X(1)
          X(16)=X(1)
          X(17)=X(1)
          DO NL2=1,JCI2
            Y(1)=H2(NL2,M2,1)
            Y(2)=Y(1)
            Y(3)=H2(NL2,M2,2)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(3)
            Y(9)=Y(3)
            Y(10)=Y(1)
            Y(11)=Y(3)
            Y(12)=Y(1)
            Y(13)=Y(1)
            Y(14)=Y(1)
            Y(15)=Y(1)
            Y(16)=Y(1)
            Y(17)=Y(1)
            DO K=1,I17
              Y0(K,NL2,NR2,M2)=Y(18-K)*X(18-K)
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
          X(4)=H3(NR3,M3,2)
          X(5)=X(1)
          X(6)=X(1)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(4)
          X(10)=X(1)
          X(11)=X(1)
          X(12)=X(1)
          X(13)=X(4)
          X(14)=X(1)
          X(15)=X(4)
          X(16)=X(1)
          X(17)=X(1)
          DO NL3=1,JCI3
            Y(1)=H3(NL3,M3,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=H3(NL3,M3,2)
            Y(6)=Y(1)
            Y(7)=Y(1)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(5)
            Y(11)=Y(1)
            Y(12)=Y(1)
            Y(13)=Y(5)
            Y(14)=Y(5)
            Y(15)=Y(1)
            Y(16)=Y(1)
            Y(17)=Y(1)
            DO K=1,I17
              Z0(K,NL3,NR3,M3)=Y(18-K)*X(18-K)
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
          X(6)=H4(NR4,M4,2)
          X(7)=X(1)
          X(8)=X(1)
          X(9)=X(1)
          X(10)=X(1)
          X(11)=X(6)
          X(12)=X(1)
          X(13)=X(1)
          X(14)=X(6)
          X(15)=X(1)
          X(16)=X(6)
          X(17)=X(1)
          DO NL4=1,JCI4
            Y(1)=H4(NL4,M4,1)
            Y(2)=Y(1)
            Y(3)=Y(1)
            Y(4)=Y(1)
            Y(5)=Y(1)
            Y(6)=Y(1)
            Y(7)=H4(NL4,M4,2)
            Y(8)=Y(1)
            Y(9)=Y(1)
            Y(10)=Y(1)
            Y(11)=Y(1)
            Y(12)=Y(7)
            Y(13)=Y(1)
            Y(14)=Y(1)
            Y(15)=Y(7)
            Y(16)=Y(7)
            Y(17)=Y(1)
            DO K=1,I17
              W0(K,NL4,NR4,M4)=Y(18-K)*X(18-K)
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
                    DO K=1,I17
                      TEMP3(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=0
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
C**START 3-MODE INTEGRATION
CCCC    DO M2=1,MM2/MD2
        DO M2=1,MM2
          DO IRHS4=1,JCI4
            DO IRHS3=1,JCI3
              DO ILHS4=1,JCI4
                DO ILHS3=1,JCI3
                  DO K=1,I17
                    TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=0
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 2-MODE INTEGRATION
CCCC      DO M3=1,MM3/MD3
          DO M3=1,MM3
            DO IRHS4=1,JCI4
              DO ILHS4=1,JCI4
                DO K=1,I17
                  TEMP1(K,ILHS4,IRHS4)=0
                END DO
              END DO
            END DO
C**START 1-MODE INTEGRATION
CCCC        DO M4=1,MM4/MD4
            DO M4=1,MM4
              DO I=1,17
                C(I)=0
              END DO
              IF(JCOUPL.GT.0)THEN
                IF(ICOUPC.GT.3)THEN
                  C(1)=VC(M4,M3,M2,M1,5)*IFACTC
                  C(2)=VC(M4,M3,M2,M1,6)*IFACTC
                  C(3)=VC(M4,M3,M2,M1,6)*IFACTC
                  C(4)=VC(M4,M3,M2,M1,7)*IFACTC
                  C(5)=VC(M4,M3,M2,M1,7)*IFACTC
                  C(6)=VC(M4,M3,M2,M1,8)*IFACTC
                  C(7)=VC(M4,M3,M2,M1,8)*IFACTC
                  C(8)=VC(M4,M3,M2,M1,9)*IFACTC
                  C(9)=VC(M4,M3,M2,M1,10)*IFACTC
                  C(10)=VC(M4,M3,M2,M1,10)*IFACTC
                  C(11)=VC(M4,M3,M2,M1,11)*IFACTC
                  C(12)=VC(M4,M3,M2,M1,11)*IFACTC
                  C(13)=VC(M4,M3,M2,M1,12)*IFACTC
                  C(14)=VC(M4,M3,M2,M1,13)*IFACTC
                  C(15)=VC(M4,M3,M2,M1,13)*IFACTC
                  C(16)=VC(M4,M3,M2,M1,14)*IFACTC
                  C(17)=VC(M4,M3,M2,M1,15)*IFACTC
                  IF(J21.GT.1)
     1            C(17)=C(17)+VR(KROT,M4,M3,M2,M1)*IFACTC
                END IF
                IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(17)=C(17)+
     1          VP(M4,M3,M2,M1)*IFACTL
              ELSE
                IF(ICOUPC.GT.3)THEN
                  C(1)=VCR(M4,M3,M2,M1,5)*IFACTC
                  C(2)=VCR(M4,M3,M2,M1,6)*IFACTC
                  C(3)=VCR(M4,M3,M2,M1,6)*IFACTC
                  C(4)=VCR(M4,M3,M2,M1,7)*IFACTC
                  C(5)=VCR(M4,M3,M2,M1,7)*IFACTC
                  C(6)=VCR(M4,M3,M2,M1,8)*IFACTC
                  C(7)=VCR(M4,M3,M2,M1,8)*IFACTC
                  C(8)=VCR(M4,M3,M2,M1,9)*IFACTC
                  C(9)=VCR(M4,M3,M2,M1,10)*IFACTC
                  C(10)=VCR(M4,M3,M2,M1,10)*IFACTC
                  C(11)=VCR(M4,M3,M2,M1,11)*IFACTC
                  C(12)=VCR(M4,M3,M2,M1,11)*IFACTC
                  C(13)=VCR(M4,M3,M2,M1,12)*IFACTC
                  C(14)=VCR(M4,M3,M2,M1,13)*IFACTC
                  C(15)=VCR(M4,M3,M2,M1,13)*IFACTC
                  C(16)=VCR(M4,M3,M2,M1,14)*IFACTC
                  C(17)=VCR(M4,M3,M2,M1,15)*IFACTC
                  IF(J21.GT.1)
     1            C(17)=C(17)+VRR(KROT,M4,M3,M2,M1)*IFACTC
                END IF
                IF(IWHICH.GE.0.OR.MOLINC.LE.0)C(17)=C(17)+
     1          VPR(M4,M3,M2,M1)*IFACTL
              END IF
              DO IRHS4=1,JCI4
                DO ILHS4=1,JCI4
                  DO K=1,I17
                    TEMP1(K,ILHS4,IRHS4)=TEMP1(K,ILHS4,IRHS4)+
     1              W0(K,ILHS4,IRHS4,M4)*C(18-K)
                  END DO
                END DO
              END DO
            END DO
C**END 1-MODE INTEGRATION
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO K=1,I17
                      TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                TEMP2(K,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                Z0(K,ILHS3,IRHS3,M3)*TEMP1(K,ILHS4,IRHS4)
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
                      DO K=1,I17
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

C**NSIZE IS NO. UNIQUE INTEGRALS (4-DIM)
        DO IRHS=1,NSIZE
          NR1=IP4(IRHS,1)
          NR2=IP4(IRHS,2)
          NR3=IP4(IRHS,3)
          NR4=IP4(IRHS,4)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP4(ILHS,1)
            NL2=IP4(ILHS,2)
            NL3=IP4(ILHS,3)
            NL4=IP4(ILHS,4)
            DO K=1,I17
              XA4(ILHS+J0)=XA4(ILHS+J0)+
     1        TEMP3(K,NL2,NL3,NL4,NR2,NR3,NR4)*X0(K,NL1,NR1,M1)
            END DO
          END DO
        END DO
      END DO
C**END 4-MODE INTEGRATION
      CALL MATOUT(XA4,XRA4,NSIZE,24)
      GO TO 7777
6666  CALL MATIN(XA4,XRA4,NSIZE,24)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI4

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
      JSR=1
      DO 9999 ISMR1=1,NVSYM
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=IROFF
            NONZC1(2,NONZ1)=ILOFF
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
      END IF
      IF(NONZ1.EQ.0)GO TO 5555
      NONZC1(3,NONZ1+1)=0

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
            DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=JROFF
            NONZC2(2,NONZ2)=JLOFF
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
      END IF
      IF(NONZ2.EQ.0)GO TO 5555
      NONZC2(3,NONZ2+1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM4A.EQ.0)THEN
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C     END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CASE 1..K IS IN CONTRACTED SCHEME, L/N/M IS IN UNCONTRACTED SCHEME
      IF(ICCONT(1).GT.0)THEN

        IF(NONC1.LE.NONC2)GO TO 3003
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        IPREV=0
        DO INH1=1,NONZ1
          IRH1=NONZC1(3,INH1)
C         IF(INH1.NE.NONZ1)THEN
            INEXT=NONZC1(3,INH1+1)
C         ELSE
C           INEXT=0
C         END IF
          ILH1=NONZC1(4,INH1)
          IF(IPREV.NE.IRH1)THEN
            DO INH2=1,NONZ2
              DO IKL=1,NKVALL
                XK2(IKL,INH2)=0
              END DO
            END DO
          END IF
          IROFF=NONZC1(1,INH1)
          ILOFF=NONZC1(2,INH1)
          DO I=1,NC1
            NCR1(I)=IPC1(IROFF,MC1(I))
            NCL1(I)=IPC1(ILOFF,MC1(I))
          END DO

C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR/ M)
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            ILH2=NONZC2(4,INH2)
            JROFF=NONZC2(1,INH2)
            JLOFF=NONZC2(2,INH2)
            DO I=1,NC2
              NCR2(I)=IPC2(JROFF,MC2(I))
              NCL2(I)=IPC2(JLOFF,MC2(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP4(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP4(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1003
            END DO
1003        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP4(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP4(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2003
            END DO
2003        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA4(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP4(MR,1)
              NL1=IP4(ML,1)
              NR2=IP4(MR,2)
              NL2=IP4(ML,2)
              NR3=IP4(MR,3)
              NL3=IP4(ML,3)
              NR4=IP4(MR,4)
              NL4=IP4(ML,4)
              DO I=1,NP4(IND4)
                K=JP4(I,IND4,1)+1
                L=JP4(I,IND4,2)+1
                N=JP4(I,IND4,3)+1
                M=JP4(I,IND4,4)+1
                ZYX=ZYX+CP4(I,IND4)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)
              END DO
C**ANALYTIC
            END IF
            DO IKL=1,NKVALL
              XK2(IKL,INH2)=XK2(IKL,INH2)+(XYZ+ZYX)*
     1        CFS1(ILH1,IKL,KEL,KSML)
            END DO
          END DO
          IF(IRH1.NE.INEXT)THEN
            DO INH2=1,NONZ2
              IRH2=NONZC2(3,INH2)
              ILH2=NONZC2(4,INH2)
              IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
                DO IKR=1,NKVALR
                  CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                  IRHS=(IKR-1)*LCSZR+IRH2
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO IKL=1,NKVALL
                    ILHS=(IKL-1)*LCSZL+ILH2
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(IKL,INH2)
                  END DO
                END DO
C**OFF-DIAGONAL BLOCK - END
              ELSE
C**DIAGONAL BLOCK - START
                DO IKR=1,NKVALR
                  CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                  IRHS=(IKR-1)*LCSZR+IRH2
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO IKL=1,IKR
                    IF(IKL.EQ.IKR.AND.ILH2.GT.IRH2)GO TO 1004
                    ILHS=(IKL-1)*LCSZL+ILH2
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(IKL,INH2)
1004                CONTINUE
                  END DO
                END DO
C**DIAGONAL BLOCK - END
              END IF
            END DO
          END IF
          IPREV=IRH1
        END DO
        GO TO 5555

3003    CONTINUE
        DO IRH1=1,KCSZR
          DO ILH1=1,KCSZL
            XK1(ILH1,IRH1)=0
          END DO
        END DO
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR/ M)

        DO INH2=1,NONZ2
          IRH2=NONZC2(3,INH2)
          ILH2=NONZC2(4,INH2)
          JROFF=NONZC2(1,INH2)
          JLOFF=NONZC2(2,INH2)
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

C**CONTRACTION SCHEME  '1' (K)
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            ILH1=NONZC1(4,INH1)
            IROFF=NONZC1(1,INH1)
            ILOFF=NONZC1(2,INH1)
            DO I=1,NC1
              NCR1(I)=IPC1(IROFF,MC1(I))
              NCL1(I)=IPC1(ILOFF,MC1(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP4(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP4(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1000
            END DO
1000        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP4(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP4(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2000
            END DO
2000        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA4(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP4(MR,1)
              NL1=IP4(ML,1)
              NR2=IP4(MR,2)
              NL2=IP4(ML,2)
              NR3=IP4(MR,3)
              NL3=IP4(ML,3)
              NR4=IP4(MR,4)
              NL4=IP4(ML,4)
              DO I=1,NP4(IND4)
                K=JP4(I,IND4,1)+1
                L=JP4(I,IND4,2)+1
                N=JP4(I,IND4,3)+1
                M=JP4(I,IND4,4)+1
                ZYX=ZYX+CP4(I,IND4)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)
              END DO
C**ANALYTIC
            END IF
            XK1(ILH1,IRH1)=XYZ+ZYX
          END DO

C************************************DGEMM (RHS)
          CALL DGEMM('N','N',KCSZL,NKVALR,KCSZR,1.0D0,XK1(1,1),
     &    KTEMP,
     &    CFS1(1,1,KEL,KSMR),ISIZXX,0.0D0,TEMP(1,1),KTEMP)

C************************************DGEMM (LHS)
          CALL DGEMM('T','N',NKVALL,NKVALR,KCSZL,1.0D0,
     &    CFS1(1,1,KEL,KSML),
     &           ISIZXX,TEMP,KTEMP,0.0D0,XCON2(1,1),MAXVAL)

          IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
            DO IKR=1,NKVALR
              IRHS=(IKR-1)*LCSZR+IRH2
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO IKL=1,NKVALL
                ILHS=(IKL-1)*LCSZL+ILH2
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(IKL,IKR)
              END DO
            END DO
C**OFF-DIAGONAL BLOCK - END
          ELSE
C**DIAGONAL BLOCK - START
            DO IKR=1,NKVALR
              IRHS=(IKR-1)*LCSZR+IRH2
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO IKL=1,IKR
                IF(IKL.EQ.IKR.AND.ILH2.GT.IRH2)GO TO 1002
                ILHS=(IKL-1)*LCSZL+ILH2
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(IKL,IKR)
1002            CONTINUE
              END DO
            END DO
C**DIAGONAL BLOCK - END
          END IF
        END DO

C**CASE 2..L/N/M IS IN CONTRACTED SCHEME, K IS IN UNCONTRACTED SCHEME
      ELSE

        IF(NONC1.GT.NONC2)GO TO 4004
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L AND/OR N)
        IPREV=0
        DO INH2=1,NONZ2
          IRH2=NONZC2(3,INH2)
C         IF(INH2.NE.NONZ2)THEN
            INEXT=NONZC2(3,INH2+1)
C         ELSE
C           INEXT=0
C         END IF
          ILH2=NONZC2(4,INH2)
          IF(IPREV.NE.IRH2)THEN
            DO INH1=1,NONZ1
              DO ILL=1,NLVALL
                XK2(ILL,INH1)=0
              END DO
            END DO
          END IF
          JROFF=NONZC2(1,INH2)
          JLOFF=NONZC2(2,INH2)
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            ILH1=NONZC1(4,INH1)
            IROFF=NONZC1(1,INH1)
            ILOFF=NONZC1(2,INH1)
            DO I=1,NC1
              NCR1(I)=IPC1(IROFF,MC1(I))
              NCL1(I)=IPC1(ILOFF,MC1(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP4(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP4(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1005
            END DO
1005        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP4(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP4(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2005
            END DO
2005        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA4(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP4(MR,1)
              NL1=IP4(ML,1)
              NR2=IP4(MR,2)
              NL2=IP4(ML,2)
              NR3=IP4(MR,3)
              NL3=IP4(ML,3)
              NR4=IP4(MR,4)
              NL4=IP4(ML,4)
              DO I=1,NP4(IND4)
                K=JP4(I,IND4,1)+1
                L=JP4(I,IND4,2)+1
                N=JP4(I,IND4,3)+1
                M=JP4(I,IND4,4)+1
                ZYX=ZYX+CP4(I,IND4)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)
              END DO
C**ANALYTIC
            END IF
            DO ILL=1,NLVALL
              XK2(ILL,INH1)=XK2(ILL,INH1)+(XYZ+ZYX)*
     1        CFS2(ILH2,ILL,KEL,LSML)
            END DO
          END DO
          IF(IRH2.NE.INEXT)THEN
            DO INH1=1,NONZ1
              IRH1=NONZC1(3,INH1)
              ILH1=NONZC1(4,INH1)
              IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
                DO ILR=1,NLVALR
                  CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                  IRHS=(ILR-1)*KCSZR+IRH1
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO ILL=1,NLVALL
                    ILHS=(ILL-1)*KCSZL+ILH1
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(ILL,INH1)
                  END DO
                END DO
C**OFF-DIAGONAL BLOCK - END
              ELSE
C**DIAGONAL BLOCK - START
                DO ILR=1,NLVALR
                  CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                  IRHS=(ILR-1)*KCSZR+IRH1
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO ILL=1,ILR
                    IF(ILL.EQ.ILR.AND.ILH1.GT.IRH1)GO TO 2004
                    ILHS=(ILL-1)*KCSZL+ILH1
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(ILL,INH1)
2004                CONTINUE
                  END DO
                END DO
C**DIAGONAL BLOCK - END
              END IF
            END DO
          END IF
          IPREV=IRH2
        END DO
        GO TO 5555

4004    CONTINUE
        DO IRH2=1,LCSZR
          DO ILH2=1,LCSZL
            XK1(ILH2,IRH2)=0
          END DO
        END DO

C**CONTRACTION SCHEME '1' (K)
        DO INH1=1,NONZ1
          IRH1=NONZC1(3,INH1)
          ILH1=NONZC1(4,INH1)
          JROFF=NONZC1(1,INH1)
          JLOFF=NONZC1(2,INH1)
          DO I=1,NC1
            NCR1(I)=IPC1(JROFF,MC1(I))
            NCL1(I)=IPC1(JLOFF,MC1(I))
          END DO

C**CONTRACTION SCHEME  '2' (L AND/OR N)
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            ILH2=NONZC2(4,INH2)
            IROFF=NONZC2(1,INH2)
            ILOFF=NONZC2(2,INH2)
            DO I=1,NC2
              NCR2(I)=IPC2(IROFF,MC2(I))
              NCL2(I)=IPC2(ILOFF,MC2(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP4(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP4(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1001
            END DO
1001        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP4(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP4(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2001
            END DO
2001        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA4(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP4(MR,1)
              NL1=IP4(ML,1)
              NR2=IP4(MR,2)
              NL2=IP4(ML,2)
              NR3=IP4(MR,3)
              NL3=IP4(ML,3)
              NR4=IP4(MR,4)
              NL4=IP4(ML,4)
              DO I=1,NP4(IND4)
                K=JP4(I,IND4,1)+1
                L=JP4(I,IND4,2)+1
                N=JP4(I,IND4,3)+1
                M=JP4(I,IND4,4)+1
                ZYX=ZYX+CP4(I,IND4)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)
              END DO
C**ANALYTIC
            END IF
            XK1(ILH2,IRH2)=XYZ+ZYX
          END DO

C************************************DGEMM (RHS)
          CALL DGEMM('N','N',LCSZL,NLVALR,LCSZR,1.0D0,XK1(1,1),
     &    KTEMP,CFS2(1,1,KEL,LSMR),ISIZXX,0.0D0,TEMP(1,1),KTEMP)

C************************************DGEMM (LHS)
          CALL DGEMM('T','N',NLVALL,NLVALR,LCSZL,1.0D0,
     &    CFS2(1,1,KEL,LSML),ISIZXX,TEMP,KTEMP,0.0D0,XCON2(1,1),MAXVAL)

          IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
            DO ILR=1,NLVALR
              IRHS=(ILR-1)*KCSZR+IRH1
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO ILL=1,NLVALL
                ILHS=(ILL-1)*KCSZL+ILH1
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(ILL,ILR)
              END DO
            END DO
C**OFF-DIAGONAL BLOCK - END
          ELSE
C**DIAGONAL BLOCK - START
            DO ILR=1,NLVALR
              IRHS=(ILR-1)*KCSZR+IRH1
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO ILL=1,ILR
                IF(ILL.EQ.ILR.AND.ILH1.GT.IRH1)GO TO 2002
                ILHS=(ILL-1)*KCSZL+ILH1
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(ILL,ILR)
2002            CONTINUE
              END DO
            END DO
C**DIAGONAL BLOCK - END
          END IF
        END DO
      END IF
C*************************************

5555  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE
      IF(ITIM4A.EQ.0)THEN
        IF(IND.EQ.1)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM4A=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCV4A(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,ITMODE,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,HTAU,XQTAU,NN1,
     2MM1,NN2,MM2,NN3,MM3,NN4,MM4,NNTAU,MMTAU,IP,ISIZMX,IPC,ICSIZE,
     3IPSIZE,
     3KCONT,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,NSIZE,XK,TEMP,XCON,NVAL,
     4ISTART,IEND,TEMP2,TEMP3,TEMP4,JCI1,JCI2,JCI3,JCI4,JCIM,X0,Y0,Z0,
     5W0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS,ISIZXX,
     6NVALX,KEL21,NVSYMX,I36) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,21),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,21),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      DIMENSION X(36),Y(36),C(36)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC(IPSIZE,1),IP5(ISIZE5,5)
      DIMENSION XA(1),XA5(1),XK(ICSIZE,ICSIZE)
      DIMENSION TEMP(ICSIZE,NVAL),XCON(NVAL,NVAL)
CCCC  DIMENSION TEMP1(I36,JCI1,JCI2,JCI3,JCI4,JCI1,JCI2,JCI3,JCI4)
      DIMENSION TEMP2(I36,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION TEMP3(I36,JCI3,JCI4,JCI3,JCI4),TEMP4(I36,JCI4,JCI4)
      DIMENSION X0(I36,JCI1,JCI1,MM1),Y0(I36,JCI2,JCI2,MM2)
      DIMENSION Z0(I36,JCI3,JCI3,MM3),W0(I36,JCI4,JCI4,MM4)
      DIMENSION T0(I36,JCIM,JCIM,MMTAU)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      IF(N3.EQ.NT.AND.MDT.EQ.2)MD3=1
      IF(N4.EQ.NT.AND.MDT.EQ.2)MD4=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      IF(N1T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N1T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N2T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N3T=ISYMP(N3,NT)
      IF(N3T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N12T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N13T=ISYMP(N13,NT)
      IF(N13T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N23T=ISYMP(N23,NT)
      IF(N23T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      N123T=ISYMP(N123,NT)
      IF(N123T.EQ.N4.AND.MDT.EQ.2)MD4=1
      CALL VDCV4A(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,ITMODE,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,HTAU,XQTAU,NN1,
     2MM1,MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,
     3NNTAU,MMTAU,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP5,
     4ISIZE5,XA5,
     4XRA5,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,TEMP2,TEMP3,TEMP4,JCI1,
     5JCI2,JCI3,JCI4,JCIM,X0,Y0,Z0,W0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,
     6MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,I36)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCV4A(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,ITMODE,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,HTAU,XQTAU,NN1,
     2MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,NN4,MH4,MM4,NNTAU,MMTAU,IP,
     3ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,
     4NSIZE,XK,
     4TEMP,XCON,NVAL,ISTART,IEND,TEMP2,TEMP3,TEMP4,JCI1,JCI2,JCI3,JCI4,
     5JCIM,X0,Y0,Z0,W0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,
     6CFS,ISIZXX,NVALX,KEL21,NVSYMX,I36)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,21),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,21),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      DIMENSION X(36),Y(36),C(36)
      DIMENSION MODINT(NMODE)
C     DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3),H4(NN4,MH4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC(IPSIZE,1),IP5(ISIZE5,5)
      DIMENSION XA(1),XA5(1),XK(ICSIZE,ICSIZE)
      DIMENSION TEMP(ICSIZE,NVAL),XCON(NVAL,NVAL)
CCCC  DIMENSION TEMP1(I36,JCI1,JCI2,JCI3,JCI4,JCI1,JCI2,JCI3,JCI4)
      DIMENSION TEMP2(I36,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION TEMP3(I36,JCI3,JCI4,JCI3,JCI4),TEMP4(I36,JCI4,JCI4)
      DIMENSION X0(I36,JCI1,JCI1,MM1),Y0(I36,JCI2,JCI2,MM2)
      DIMENSION Z0(I36,JCI3,JCI3,MM3),W0(I36,JCI4,JCI4,MM4)
      DIMENSION T0(I36,JCIM,JCIM,MMTAU)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/NCREC/NREC(6),MAXBUF(6)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/MAXLAN/LANMAX
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VCCV4A'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1', 'MOD2',
C**'MOD3', 'MOD4' AND 'TAU'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0CV4
      REWIND 30
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=JNTFAC(NMODE,ICOUPL,4)
      IFACTC=JNTFAC(NMODE,ICOUPC,4)
C**INTRINSIC
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1

      KA=KROT/2
      LMAX=1+MOD(KA,2)
      FACTRC=0.D0
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      IF(N3.EQ.NT.AND.MDT.EQ.2)MD3=1
      IF(N4.EQ.NT.AND.MDT.EQ.2)MD4=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      IF(N1T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N1T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N2T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N3T=ISYMP(N3,NT)
      IF(N3T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N12T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N13T=ISYMP(N13,NT)
      IF(N13T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N23T=ISYMP(N23,NT)
      IF(N23T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      N123T=ISYMP(N123,NT)
      IF(N123T.EQ.N4.AND.MDT.EQ.2)MD4=1
      MD=MD1*MD2*MD3*MD4*MDT

C**FORM INDIVIDUAL INTEGRATION TERMS (START)
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          DO NRR=1,JCIM
            NR=NRR+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,36
              X(I)=X(1)
            END DO
            X(9)=HTAU(NR,MTAU,2,LMAX)*MD
            X(18)=X(9)
            X(25)=X(9)
            X(30)=X(9)
            X(33)=X(9)
            X(35)=X(9)
            DO NLL=1,JCIM
              NL=NLL+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO I=2,36
                Y(I)=Y(1)
              END DO
              Y(10)=HTAU(NL,MTAU,2,LMAX)
              Y(19)=Y(10)
              Y(26)=Y(10)
              Y(31)=Y(10)
              Y(34)=Y(10)
              Y(35)=Y(10)
              DO K=1,I36
                T0(K,NLL,NRR,MTAU)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NRR=1,JCIM
            NR=2*NRR-MOD(NRR,2)+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,21
              X(I)=X(1)
            END DO
            X(5)=HTAU(NR,MTAU,2,LMAX)*MD
            X(10)=X(5)
            X(14)=X(5)
            X(17)=X(5)
            X(19)=X(5)
            X(20)=HTAU(NR,MTAU,3,LMAX)*MD
            X(22)=HTAU(NR,MTAU,1,LMAX)*FACTRC*MD
            DO NLL=1,JCIM
              NL=2*NLL-MOD(NLL,2)+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO K=1,22
                T0(K,NLL,NRR,MTAU)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M1=1,MM1/MD1
      DO M1=1,MM1
        IF(.NOT.LINEAR)THEN
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,36
              X(I)=X(2)
            END DO
            X(11)=X(1)
            X(13)=X(1)
            X(15)=X(1)
            X(17)=X(1)
            X(19)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              Y(2)=H1(NL1,M1,2)
              DO I=3,36
                Y(I)=Y(1)
              END DO
              Y(11)=Y(2)
              Y(12)=Y(2)
              Y(14)=Y(2)
              Y(16)=Y(2)
              Y(18)=Y(2)
              DO K=1,I36
                X0(K,NL1,NR1,M1)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,22
              X(I)=X(2)
            END DO
            X(6)=H1(NR1,M1,3)
            X(7)=X(1)
            X(8)=X(1)
            X(9)=X(1)
            X(10)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              DO K=1,22
                X0(K,NL1,NR1,M1)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M2=1,MM2/MD2
      DO M2=1,MM2
        IF(.NOT.LINEAR)THEN
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            DO I=2,36
              X(I)=X(1)
            END DO
            X(3)=H2(NR2,M2,2)
            X(12)=X(3)
            X(20)=X(3)
            X(22)=X(3)
            X(24)=X(3)
            X(26)=X(3)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO I=2,36
                Y(I)=Y(1)
              END DO
              Y(4)=H2(NL2,M2,2)
              Y(13)=Y(4)
              Y(20)=Y(4)
              Y(21)=Y(4)
              Y(23)=Y(4)
              Y(25)=Y(4)
              DO K=1,I36
                Y0(K,NL2,NR2,M2)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            X(2)=H2(NR2,M2,2)
            DO I=3,22
              X(I)=X(1)
            END DO
            X(7)=X(2)
            X(11)=H2(NR2,M2,3)
            X(12)=X(2)
            X(13)=X(2)
            X(14)=X(2)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO K=1,22
                Y0(K,NL2,NR2,M2)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M3=1,MM3/MD3
      DO M3=1,MM3
        IF(.NOT.LINEAR)THEN
          DO NR3=1,JCI3
            X(1)=H3(NR3,M3,1)
            DO I=2,36
              X(I)=X(1)
            END DO
            X(5)=H3(NR3,M3,2)
            X(14)=X(5)
            X(21)=X(5)
            X(27)=X(5)
            X(29)=X(5)
            X(31)=X(5)
            DO NL3=1,JCI3
              Y(1)=H3(NL3,M3,1)
              DO I=2,36
                Y(I)=Y(1)
              END DO
              Y(6)=H3(NL3,M3,2)
              Y(15)=Y(6)
              Y(22)=Y(6)
              Y(27)=Y(6)
              Y(28)=Y(6)
              Y(30)=Y(6)
              DO K=1,I36
                  Z0(K,NL3,NR3,M3)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR3=1,JCI3
            X(1)=H3(NR3,M3,1)
            DO I=2,22
              X(I)=X(1)
            END DO
            X(3)=H3(NR3,M3,2)
            X(8)=X(3)
            X(12)=X(3)
            X(15)=H3(NR3,M3,3)
            X(16)=X(3)
            X(17)=X(3)
            DO NL3=1,JCI3
              Y(1)=H3(NL3,M3,1)
              DO K=1,22
                  Z0(K,NL3,NR3,M3)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M4=1,MM4/MD4
      DO M4=1,MM4
        IF(.NOT.LINEAR)THEN
          DO NR4=1,JCI4
            X(1)=H4(NR4,M4,1)
            DO I=2,36
              X(I)=X(1)
            END DO
            X(7)=H4(NR4,M4,2)
            X(16)=X(7)
            X(23)=X(7)
            X(28)=X(7)
            X(32)=X(7)
            X(34)=X(7)
            DO NL4=1,JCI4
              Y(1)=H4(NL4,M4,1)
              DO I=2,36
                Y(I)=Y(1)
              END DO
              Y(8)=H4(NL4,M4,2)
              Y(17)=Y(8)
              Y(24)=Y(8)
              Y(29)=Y(8)
              Y(32)=Y(8)
              Y(33)=Y(8)
              DO K=1,I36
                W0(K,NL4,NR4,M4)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR4=1,JCI4
            X(1)=H4(NR4,M4,1)
            DO I=2,22
              X(I)=X(1)
            END DO
            X(4)=H4(NR4,M4,2)
            X(9)=X(4)
            X(13)=X(4)
            X(16)=X(4)
            X(18)=H4(NR4,M4,3)
            X(19)=X(4)
            DO NL4=1,JCI4
              Y(1)=H4(NL4,M4,1)
              DO K=1,22
                W0(K,NL4,NR4,M4)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)

      LSTART=1
      KLAN=LANMAX*(LANMAX+1)/2
      LSIZE=1
      ISKIP=0
5555  CONTINUE
      IF(ISKIP.NE.0)THEN
        REWIND 64
        REWIND 74
        REWIND 84
      END IF
      LEND=LSTART-1
      KXA=0
      DO I=LSTART,NSIZE
        KXA=KXA+LSIZE
        IF(KXA.GT.KLAN)THEN
          KXA=KXA-LSIZE
          GO TO 5504
        END IF
        LSIZE=LSIZE+1
        LEND=LEND+1
      END DO
5504  CONTINUE
      IF(LEND.EQ.0)STOP 'LANMAX TOO SMALL'
C**ZEROISE MATRIX
      CALL DIAGL(KXA,XA5)

C**LOOP ROUND TAU (START 5-MODE INTEGRATION)
      ITAU=INIT-INCTAU
      DO MTAU=1,MMTAU/MDT
        ITAU=ITAU+INCTAU
CCCC    IF(ITAU.GT.362)ITAU=ITAU-360
        IF(ITAU.GT.722)ITAU=ITAU-720
 
C***********************************************************

        IF(JCOUPC.GE.0)THEN
          IF(J21.GT.1.AND.ICOUPC.GE.4)READ(64)VR
          IF(ICOUPC.GE.4)READ(84)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPC.GE.4)READ(64)VRR
          IF(ICOUPC.GE.4)READ(84)VCR
        END IF
        IF(JCOUPL.GT.0)THEN
          READ(74)VP
        ELSE
          READ(74)VPR
        END IF
 
C***********************************************************

C**********************************
CCCC    DO IRHS4=1,JCI4
CCCC      DO IRHS3=1,JCI3
CCCC        DO IRHS2=1,JCI2
CCCC          DO IRHS1=1,JCI1
CCCC            DO ILHS4=1,JCI4
CCCC              DO ILHS3=1,JCI3
CCCC                DO ILHS2=1,JCI2
CCCC                  DO ILHS1=1,JCI1
CCCC                    DO K=1,I36
CCCC      TEMP1(K,ILHS1,ILHS2,ILHS3,ILHS4,IRHS1,IRHS2,IRHS3,IRHS4)=0.D0
CCCC                    END DO
CCCC                  END DO
CCCC                END DO
CCCC              END DO
CCCC            END DO
CCCC          END DO
CCCC        END DO
CCCC      END DO
CCCC    END DO
C**********************************
C**START 4-MODE INTEGRATION
C       DO M1=1,MM1/MD1
        DO M1=1,MM1
          DO IRHS4=1,JCI4
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO K=1,I36
                      TEMP2(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=0.D0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 3-MODE INTEGRATION
C         DO M2=1,MM2/MD2
          DO M2=1,MM2
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO K=1,I36
                      TEMP3(K,ILHS3,ILHS4,IRHS3,IRHS4)=0.D0
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
C           DO M3=1,MM3/MD3
            DO M3=1,MM3
              DO IRHS4=1,JCI4
                DO ILHS4=1,JCI4
                  DO K=1,I36
                    TEMP4(K,ILHS4,IRHS4)=0.D0
                  END DO
                END DO
              END DO
C**START 1-MODE INTEGRATION
C             DO M4=1,MM4/MD4
              DO M4=1,MM4
                DO I=1,36
                  C(I)=0
                END DO
                IF(JCOUPL.GT.0)THEN
                  IF(ICOUPC.GE.4)THEN
                    C(1)=VC(M4,M3,M2,M1,1)*IFACTC
                    C(2)=VC(M4,M3,M2,M1,1)*IFACTC
                    C(3)=VC(M4,M3,M2,M1,2)*IFACTC
                    C(4)=VC(M4,M3,M2,M1,2)*IFACTC
                    C(5)=VC(M4,M3,M2,M1,3)*IFACTC
                    C(6)=VC(M4,M3,M2,M1,3)*IFACTC
                    C(7)=VC(M4,M3,M2,M1,4)*IFACTC
                    C(8)=VC(M4,M3,M2,M1,4)*IFACTC
                    C(9)=VC(M4,M3,M2,M1,5)*IFACTC
                    C(10)=VC(M4,M3,M2,M1,5)*IFACTC
                    C(11)=VC(M4,M3,M2,M1,6)*IFACTC
                    C(12)=VC(M4,M3,M2,M1,7)*IFACTC
                    C(13)=VC(M4,M3,M2,M1,7)*IFACTC
                    C(14)=VC(M4,M3,M2,M1,8)*IFACTC
                    C(15)=VC(M4,M3,M2,M1,8)*IFACTC
                    C(16)=VC(M4,M3,M2,M1,9)*IFACTC
                    C(17)=VC(M4,M3,M2,M1,9)*IFACTC
                    C(18)=VC(M4,M3,M2,M1,10)*IFACTC
                    C(19)=VC(M4,M3,M2,M1,10)*IFACTC
                    C(20)=VC(M4,M3,M2,M1,11)*IFACTC
                    C(21)=VC(M4,M3,M2,M1,12)*IFACTC
                    C(22)=VC(M4,M3,M2,M1,12)*IFACTC
                    C(23)=VC(M4,M3,M2,M1,13)*IFACTC
                    C(24)=VC(M4,M3,M2,M1,13)*IFACTC
                    C(25)=VC(M4,M3,M2,M1,14)*IFACTC
                    C(26)=VC(M4,M3,M2,M1,14)*IFACTC
                    C(27)=VC(M4,M3,M2,M1,15)*IFACTC
                    C(28)=VC(M4,M3,M2,M1,16)*IFACTC
                    C(29)=VC(M4,M3,M2,M1,16)*IFACTC
                    C(30)=VC(M4,M3,M2,M1,17)*IFACTC
                    C(31)=VC(M4,M3,M2,M1,17)*IFACTC
                    C(32)=VC(M4,M3,M2,M1,18)*IFACTC
                    C(33)=VC(M4,M3,M2,M1,19)*IFACTC
                    C(34)=VC(M4,M3,M2,M1,19)*IFACTC
                    C(35)=VC(M4,M3,M2,M1,20)*IFACTC
                    C(36)=VC(M4,M3,M2,M1,21)*IFACTC
                    IF(J21.GT.1)C(36)=C(36)+VR(KROT,M4,M3,M2,M1)*IFACTC
                  END IF
                  C(36)=C(36)+VP(M4,M3,M2,M1)*IFACTL
                ELSE
                  IF(ICOUPC.GE.4)THEN
                    C(1)=VCR(M4,M3,M2,M1,1)*IFACTC
                    C(2)=VCR(M4,M3,M2,M1,1)*IFACTC
                    C(3)=VCR(M4,M3,M2,M1,2)*IFACTC
                    C(4)=VCR(M4,M3,M2,M1,2)*IFACTC
                    C(5)=VCR(M4,M3,M2,M1,3)*IFACTC
                    C(6)=VCR(M4,M3,M2,M1,3)*IFACTC
                    C(7)=VCR(M4,M3,M2,M1,4)*IFACTC
                    C(8)=VCR(M4,M3,M2,M1,4)*IFACTC
                    C(9)=VCR(M4,M3,M2,M1,5)*IFACTC
                    C(10)=VCR(M4,M3,M2,M1,5)*IFACTC
                    C(11)=VCR(M4,M3,M2,M1,6)*IFACTC
                    C(12)=VCR(M4,M3,M2,M1,7)*IFACTC
                    C(13)=VCR(M4,M3,M2,M1,7)*IFACTC
                    C(14)=VCR(M4,M3,M2,M1,8)*IFACTC
                    C(15)=VCR(M4,M3,M2,M1,8)*IFACTC
                    C(16)=VCR(M4,M3,M2,M1,9)*IFACTC
                    C(17)=VCR(M4,M3,M2,M1,9)*IFACTC
                    C(18)=VCR(M4,M3,M2,M1,10)*IFACTC
                    C(19)=VCR(M4,M3,M2,M1,10)*IFACTC
                    C(20)=VCR(M4,M3,M2,M1,11)*IFACTC
                    C(21)=VCR(M4,M3,M2,M1,12)*IFACTC
                    C(22)=VCR(M4,M3,M2,M1,12)*IFACTC
                    C(23)=VCR(M4,M3,M2,M1,13)*IFACTC
                    C(24)=VCR(M4,M3,M2,M1,13)*IFACTC
                    C(25)=VCR(M4,M3,M2,M1,14)*IFACTC
                    C(26)=VCR(M4,M3,M2,M1,14)*IFACTC
                    C(27)=VCR(M4,M3,M2,M1,15)*IFACTC
                    C(28)=VCR(M4,M3,M2,M1,16)*IFACTC
                    C(29)=VCR(M4,M3,M2,M1,16)*IFACTC
                    C(30)=VCR(M4,M3,M2,M1,17)*IFACTC
                    C(31)=VCR(M4,M3,M2,M1,17)*IFACTC
                    C(32)=VCR(M4,M3,M2,M1,18)*IFACTC
                    C(33)=VCR(M4,M3,M2,M1,19)*IFACTC
                    C(34)=VCR(M4,M3,M2,M1,19)*IFACTC
                    C(35)=VCR(M4,M3,M2,M1,20)*IFACTC
                    C(36)=VCR(M4,M3,M2,M1,21)*IFACTC
                    IF(J21.GT.1)C(36)=C(36)+VRR(KROT,M4,M3,M2,M1)*
     1              IFACTC
                  END IF
                  C(36)=C(36)+VPR(M4,M3,M2,M1)*IFACTL
                END IF
                DO IRHS4=1,JCI4
                  DO ILHS4=1,JCI4
                    DO K=1,I36
                      TEMP4(K,ILHS4,IRHS4)=TEMP4(K,ILHS4,IRHS4)+
     1                W0(K,ILHS4,IRHS4,M4)*C(37-K)
                    END DO
                  END DO
                END DO
              END DO
C**END 1-MODE INTEGRATION
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS4=1,JCI4
                    DO ILHS3=1,JCI3
                      DO K=1,I36
                        TEMP3(K,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                  TEMP3(K,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                  Z0(K,ILHS3,IRHS3,M3)*TEMP4(K,ILHS4,IRHS4)
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
                        DO K=1,I36
                          TEMP2(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=
     1                    TEMP2(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)+
     2                    Y0(K,ILHS2,IRHS2,M2)*
     3                    TEMP3(K,ILHS3,ILHS4,IRHS3,IRHS4)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 3-MODE INTEGRATION
C**********************************
CCCC      DO IRHS4=1,JCI4
CCCC        DO IRHS3=1,JCI3
CCCC          DO IRHS2=1,JCI2
CCCC            DO IRHS1=1,JCI1
CCCC              DO ILHS4=1,JCI4
CCCC                DO ILHS3=1,JCI3
CCCC                  DO ILHS2=1,JCI2
CCCC                    DO ILHS1=1,JCI1
CCCC                      DO K=1,I36
CCCC          TEMP1(K,ILHS1,ILHS2,ILHS3,ILHS4,IRHS1,IRHS2,IRHS3,IRHS4)=
CCCC 1        TEMP1(K,ILHS1,ILHS2,ILHS3,ILHS4,IRHS1,IRHS2,IRHS3,IRHS4)+
CCCC 2        X0(K,ILHS1,IRHS1,M1)*
CCCC 3        TEMP2(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)
CCCC                      END DO
CCCC                    END DO
CCCC                  END DO
CCCC                END DO
CCCC              END DO
CCCC            END DO
CCCC          END DO
CCCC        END DO
CCCC      END DO
CCCC    END DO
C**********************************

C**NSIZE IS NO. UNIQUE INTEGRALS (4-DIM)
          J0=0
          DO IRHS=LSTART,LEND
            NR1=IP5(IRHS,1)
            NR2=IP5(IRHS,2)
            NR3=IP5(IRHS,3)
            NR4=IP5(IRHS,4)
            IRTAU=IP5(IRHS,5)
C           J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP5(ILHS,1)
              NL2=IP5(ILHS,2)
              NL3=IP5(ILHS,3)
              NL4=IP5(ILHS,4)
              ILTAU=IP5(ILHS,5)
              DO K=1,I36
CCCC            XA5(ILHS+J0)=XA5(ILHS+J0)+
CCCC 1          TEMP1(K,NL1,NL2,NL3,NL4,NR1,NR2,NR3,NR4)*
CCCC 1          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
                XA5(ILHS+J0)=XA5(ILHS+J0)+
     1          TEMP2(K,NL2,NL3,NL4,NR2,NR3,NR4)*X0(K,NL1,NR1,M1)*
     1          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
              END DO
            END DO
            J0=J0+IRHS
          END DO
C**********************************
        END DO
C**END 4-MODE INTEGRATION
C**********************************
      END DO
C**END TAU LOOP (5-MODE INTEGRATION)

      WRITE(24)LSTART,LEND,KXA
      WRITE(30)LSTART,LEND,KXA
      CALL MTOUT(XA5,XRA5,KXA,24)
      CALL MTOUT(XA5,XRA5,KXA,30)
      IF(LEND.NE.NSIZE)THEN
        LSTART=LEND+1
        ISKIP=1
        GO TO 5555
      END IF
      GO TO 7777
6666  CONTINUE
1111  READ(24)LSTART,LEND,KXA
      WRITE(30)LSTART,LEND,KXA
      CALL MTIN(XA5,XRA5,KXA,24)
      CALL MTOUT(XA5,XRA5,KXA,30)
      IF(LEND.NE.NSIZE)GO TO 1111
7777  REWIND 30
C**INITIAL READ
      READ(30)LSTART,LEND,KXA
      CALL MTIN(XA5,XRA5,KXA,30)

C***********************************ALGORITHM FROM VCV4

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR1=IPC(IROFF,MODE1)
        NR2=IPC(IROFF,MODE2)
        NR3=IPC(IROFF,MODE3)
        NR4=IPC(IROFF,MODE4)
        NRTAU=IPC(IROFF,ITMODE)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP5(IR,1).AND.NR2.EQ.IP5(IR,2).AND.
     1    NR3.EQ.IP5(IR,3).AND.NR4.EQ.IP5(IR,4).AND.NRTAU.EQ.
     2    IP5(IR,5))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1      MODE4.AND.K.NE.ITMODE.AND.(IPC(IROFF,K).NE.IPC(ILOFF,K)))
     2      IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPC(ILOFF,MODE1)
          NL2=IPC(ILOFF,MODE2)
          NL3=IPC(ILOFF,MODE3)
          NL4=IPC(ILOFF,MODE4)
          NLTAU=IPC(ILOFF,ITMODE)
C**FIND LHS INDEX
          DO IL=1,NSIZE
            IF(NL1.EQ.IP5(IL,1).AND.NL2.EQ.IP5(IL,2).AND.
     1      NL3.EQ.IP5(IL,3).AND.NL4.EQ.IP5(IL,4).AND.NLTAU.EQ.
     2      IP5(IL,5))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
CC        I=MR*(MR-1)/2+ML
          IF(MR.GE.LSTART.AND.MR.LE.LEND)THEN
C**ALREADY GOT IT
            GO TO 2222
          END IF
          IF(MR.LT.LSTART)THEN
C**PASSED IT
            REWIND 30
            READ(30)LSTART,LEND,KXA
            CALL MTIN(XA5,XRA5,KXA,30)
          END IF
3333      CONTINUE
          IF(MR.GT.LEND)THEN
C**NOT YET REACHED IT
            READ(30)LSTART,LEND,KXA
            CALL MTIN(XA5,XRA5,KXA,30)
            GO TO 3333
          END IF
2222      CONTINUE
          I=0
          DO K=LSTART,MR-1
            I=I+K
          END DO
          I=I+ML
          XYZ=XA5(I)
C**GET MATRIX ELEMENT
3000      CONTINUE
          XK(ILHS,IRHS)=XYZ*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
        J0=(IROFF)*(IROFF-1)/2+IOFF
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IP(IROFF,KCONT)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IP(IROFF,MCONT).NE.IP(ILOFF,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IP(ILOFF,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM4A.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM4A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCV4B(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,HTAU,XQTAU,NN1,
     2MM1,NN2,MM2,NN3,MM3,NN4,MM4,NNTAU,MMTAU,IP,ISIZMX,IPC1,ICSIZ1,
     3IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP5,ISIZE5,XA5,
     4XRA5,NSIZE,XK1,
     4XK2,IRL,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP2,TEMP3,TEMP4,
     5JCI1,JCI2,JCI3,JCI4,JCIM,X0,Y0,Z0,W0,T0,VP,VPR,VC,VCR,VR,VRR,J21,
     6KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,
     7NONZC1,NONZC2,I36,MAXVAL,MAXQ5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,21),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,21),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP5(ISIZE5,5),IRL(MAXQ5,MAXQ5,MAXQ5,MAXQ5,MAXQ5)
CCCC  DIMENSION TEMP1(I36,JCI1,JCI2,JCI3,JCI4,JCI1,JCI2,JCI3,JCI4)
      DIMENSION TEMP2(I36,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION TEMP3(I36,JCI3,JCI4,JCI3,JCI4),TEMP4(I36,JCI4,JCI4)
      DIMENSION X0(I36,JCI1,JCI1,MM1),Y0(I36,JCI2,JCI2,MM2)
      DIMENSION Z0(I36,JCI3,JCI3,MM3),W0(I36,JCI4,JCI4,MM4)
      DIMENSION T0(I36,JCIM,JCIM,MMTAU)
      DIMENSION X(36),Y(36),C(36)
C**************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA5(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C**************************************************************
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      IF(N3.EQ.NT.AND.MDT.EQ.2)MD3=1
      IF(N4.EQ.NT.AND.MDT.EQ.2)MD4=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      IF(N1T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N1T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N2T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N3T=ISYMP(N3,NT)
      IF(N3T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N12T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N13T=ISYMP(N13,NT)
      IF(N13T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N23T=ISYMP(N23,NT)
      IF(N23T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      N123T=ISYMP(N123,NT)
      IF(N123T.EQ.N4.AND.MDT.EQ.2)MD4=1
      CALL VDCV4B(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,HTAU,XQTAU,NN1,
     2MM1,MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,
     3NNTAU,MMTAU,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,
     4KCONT1,KCONT2,XA,
     4ISIZE,IP5,ISIZE5,XA5,XRA5,NSIZE,XK1,XK2,IRL,KTEMP,XCON2,NVAL1,
     5NVAL2,ISTART,IEND,TEMP2,TEMP3,TEMP4,JCI1,JCI2,JCI3,JCI4,JCIM,X0,
     6Y0,Z0,W0,T0,VP,VPR,VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,
     7CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,I36,MAXVAL,
     8MAXQ5)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCV4B(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MODE1,MODE2,
     1MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,HTAU,XQTAU,NN1,MH1,
     2MM1,NN2,MH2,MM2,NN3,MH3,MM3,NN4,MH4,MM4,NNTAU,MMTAU,IP,ISIZMX,
     3IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP5,
     4ISIZE5,XA5,
     4XRA5,NSIZE,XK1,XK2,IRL,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,
     5TEMP2,TEMP3,TEMP4,JCI1,JCI2,JCI3,JCI4,JCIM,X0,Y0,Z0,W0,T0,VP,VPR,
     6VC,VCR,VR,VRR,J21,KROT,MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,
     7KEL21,NVSYMX,IND,NONZC1,NONZC2,I36,MAXVAL,MAXQ5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,21),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,21),
     1VRR(J21,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD

C**TEMPORARY
      DIMENSION IRL(MAXQ5,MAXQ5,MAXQ5,MAXQ5,MAXQ5),NX5(5),MX5(5)
C**TEMPORARY

      DIMENSION MODINT(NMODE)
C     DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3),H4(NN4,MH4,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION HTAU(NNTAU,MMTAU,3,1),XQTAU(MMTAU)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP5(ISIZE5,5),TEMP(KTEMP,1)
CCCC  DIMENSION TEMP1(I36,JCI1,JCI2,JCI3,JCI4,JCI1,JCI2,JCI3,JCI4)
      DIMENSION TEMP2(I36,JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION TEMP3(I36,JCI3,JCI4,JCI3,JCI4),TEMP4(I36,JCI4,JCI4)
      DIMENSION X0(I36,JCI1,JCI1,MM1),Y0(I36,JCI2,JCI2,MM2)
      DIMENSION Z0(I36,JCI3,JCI3,MM3),W0(I36,JCI4,JCI4,MM4)
      DIMENSION T0(I36,JCIM,JCIM,MMTAU)
      DIMENSION X(36),Y(36),C(36)
C************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA5(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C************************************************************
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MDUM,INIT,INCTAU
CCCC  COMMON/SCOORD/DSTAU(362)
      COMMON/SCOORD/DSTAU(722)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TYPE/LINEAR
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODX
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/NCREC/NREC(6),MAXBUF(6)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
      COMMON/MAXLAN/LANMAX
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM4B.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCV4B'
          CALL TIMIT(1)
          CALL FLUSH(IOUT)
        END IF
      END IF

C**************************************************************
C**************************************************************
C**BOTH CONTRACTION SCHEMES HAVE AT LEAST ONE MODE
C**CHECK SCHEME '1'
      I1=KCONT1
      NC1=0
      DO NN=1,ICONT(I1)
        IF(MOD1.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE1
          JC1(NC1)=1
        END IF
        IF(MOD2.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE2
          JC1(NC1)=2
        END IF
        IF(MOD3.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE3
          JC1(NC1)=3
        END IF
        IF(MOD4.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE4
          JC1(NC1)=4
        END IF
        IF(NNMODE.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE5
          JC1(NC1)=5
        END IF
      END DO
C**CHECK SCHEME '2'
      I2=KCONT2
      NC2=0
      DO NN=1,ICONT(I2)
        IF(MOD1.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE1
          JC2(NC2)=1
        END IF
        IF(MOD2.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE2
          JC2(NC2)=2
        END IF
        IF(MOD3.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE3
          JC2(NC2)=3
        END IF
        IF(MOD4.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE4
          JC2(NC2)=4
        END IF
        IF(NNMODE.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE5
          JC2(NC2)=5
        END IF
      END DO
C**NKMOD1 IS TOTAL NUMBER OF MODES FOR SCHEME INVOLVING 'MOD1' (K)
C**NKMOD2 IS TOTAL NUMBER OF MODES FOR THE OTHER SCHEME (ONE OR ALL
C**OF 'MOD2' (L), 'MOD3' (N) AND 'MOD4' (M))
C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L),
C**'MOD3' (N) OR 'MOD4' (M)
C**IF THEY ARE IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
C**************************************************************
      IF(IND.NE.0)GO TO 7777
C**************************************************************

C**FIND MAX QUANTUM
      MAXQ5=0
      DO IR=1,NSIZE
        NR1=IP5(IR,1)
        NR2=IP5(IR,2)
        NR3=IP5(IR,3)
        NR4=IP5(IR,4)
        NR5=IP5(IR,5)
        IF(NR1.GT.MAXQ5)MAXQ5=NR1
        IF(NR2.GT.MAXQ5)MAXQ5=NR2
        IF(NR3.GT.MAXQ5)MAXQ5=NR3
        IF(NR4.GT.MAXQ5)MAXQ5=NR4
        IF(NR5.GT.MAXQ5)MAXQ5=NR5
      END DO

C***********************************ALGORITHM FROM V0CV4
      REWIND 30
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=JNTFAC(NMODE,ICOUPL,4)
      IFACTC=JNTFAC(NMODE,ICOUPC,4)
C**INTRINSIC
      IF(MOLINC.GT.0)IFACTL=1
      IF(MOLINC.GT.0)IFACTC=1

      KA=KROT/2
      LMAX=1+MOD(KA,2)
      FACTRC=0.D0
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
      IF(N1.EQ.NT.AND.MDT.EQ.2)MD1=1
      IF(N2.EQ.NT.AND.MDT.EQ.2)MD2=1
      IF(N3.EQ.NT.AND.MDT.EQ.2)MD3=1
      IF(N4.EQ.NT.AND.MDT.EQ.2)MD4=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2.AND.MDT.EQ.2)MD2=1
      IF(N1T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N1T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N2T=ISYMP(N2,NT)
      IF(N2T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N2T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N3T=ISYMP(N3,NT)
      IF(N3T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      N12T=ISYMP(N12,NT)
      IF(N12T.EQ.N3.AND.MDT.EQ.2)MD3=1
      IF(N12T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N13T=ISYMP(N13,NT)
      IF(N13T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N23T=ISYMP(N23,NT)
      IF(N23T.EQ.N4.AND.MDT.EQ.2)MD4=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      N123T=ISYMP(N123,NT)
      IF(N123T.EQ.N4.AND.MDT.EQ.2)MD4=1
      MD=MD1*MD2*MD3*MD4*MDT

C**FORM INDIVIDUAL INTEGRATION TERMS (START)
      DO MTAU=1,MMTAU/MDT
        IF(.NOT.LINEAR)THEN
          DO NRR=1,JCIM
            NR=NRR+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,36
              X(I)=X(1)
            END DO
            X(9)=HTAU(NR,MTAU,2,LMAX)*MD
            X(18)=X(9)
            X(25)=X(9)
            X(30)=X(9)
            X(33)=X(9)
            X(35)=X(9)
            DO NLL=1,JCIM
              NL=NLL+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO I=2,36
                Y(I)=Y(1)
              END DO
              Y(10)=HTAU(NL,MTAU,2,LMAX)
              Y(19)=Y(10)
              Y(26)=Y(10)
              Y(31)=Y(10)
              Y(34)=Y(10)
              Y(35)=Y(10)
              DO K=1,I36
                T0(K,NLL,NRR,MTAU)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NRR=1,JCIM
            NR=2*NRR-MOD(NRR,2)+1-MOD(KA,2)
            X(1)=HTAU(NR,MTAU,1,LMAX)*MD
            DO I=2,21
              X(I)=X(1)
            END DO
            X(5)=HTAU(NR,MTAU,2,LMAX)*MD
            X(10)=X(5)
            X(14)=X(5)
            X(17)=X(5)
            X(19)=X(5)
            X(20)=HTAU(NR,MTAU,3,LMAX)*MD
            X(22)=HTAU(NR,MTAU,1,LMAX)*FACTRC*MD
            DO NLL=1,JCIM
              NL=2*NLL-MOD(NLL,2)+1-MOD(KA,2)
              Y(1)=HTAU(NL,MTAU,1,LMAX)
              DO K=1,22
                T0(K,NLL,NRR,MTAU)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M1=1,MM1/MD1
      DO M1=1,MM1
        IF(.NOT.LINEAR)THEN
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,36
              X(I)=X(2)
            END DO
            X(11)=X(1)
            X(13)=X(1)
            X(15)=X(1)
            X(17)=X(1)
            X(19)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              Y(2)=H1(NL1,M1,2)
              DO I=3,36
                Y(I)=Y(1)
              END DO
              Y(11)=Y(2)
              Y(12)=Y(2)
              Y(14)=Y(2)
              Y(16)=Y(2)
              Y(18)=Y(2)
              DO K=1,I36
                X0(K,NL1,NR1,M1)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR1=1,JCI1
            X(1)=H1(NR1,M1,2)
            X(2)=H1(NR1,M1,1)
            DO I=3,22
              X(I)=X(2)
            END DO
            X(6)=H1(NR1,M1,3)
            X(7)=X(1)
            X(8)=X(1)
            X(9)=X(1)
            X(10)=X(1)
            DO NL1=1,JCI1
              Y(1)=H1(NL1,M1,1)
              DO K=1,22
                X0(K,NL1,NR1,M1)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M2=1,MM2/MD2
      DO M2=1,MM2
        IF(.NOT.LINEAR)THEN
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            DO I=2,36
              X(I)=X(1)
            END DO
            X(3)=H2(NR2,M2,2)
            X(12)=X(3)
            X(20)=X(3)
            X(22)=X(3)
            X(24)=X(3)
            X(26)=X(3)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO I=2,36
                Y(I)=Y(1)
              END DO
              Y(4)=H2(NL2,M2,2)
              Y(13)=Y(4)
              Y(20)=Y(4)
              Y(21)=Y(4)
              Y(23)=Y(4)
              Y(25)=Y(4)
              DO K=1,I36
                Y0(K,NL2,NR2,M2)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR2=1,JCI2
            X(1)=H2(NR2,M2,1)
            X(2)=H2(NR2,M2,2)
            DO I=3,22
              X(I)=X(1)
            END DO
            X(7)=X(2)
            X(11)=H2(NR2,M2,3)
            X(12)=X(2)
            X(13)=X(2)
            X(14)=X(2)
            DO NL2=1,JCI2
              Y(1)=H2(NL2,M2,1)
              DO K=1,22
                Y0(K,NL2,NR2,M2)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M3=1,MM3/MD3
      DO M3=1,MM3
        IF(.NOT.LINEAR)THEN
          DO NR3=1,JCI3
            X(1)=H3(NR3,M3,1)
            DO I=2,36
              X(I)=X(1)
            END DO
            X(5)=H3(NR3,M3,2)
            X(14)=X(5)
            X(21)=X(5)
            X(27)=X(5)
            X(29)=X(5)
            X(31)=X(5)
            DO NL3=1,JCI3
              Y(1)=H3(NL3,M3,1)
              DO I=2,36
                Y(I)=Y(1)
              END DO
              Y(6)=H3(NL3,M3,2)
              Y(15)=Y(6)
              Y(22)=Y(6)
              Y(27)=Y(6)
              Y(28)=Y(6)
              Y(30)=Y(6)
              DO K=1,I36
                  Z0(K,NL3,NR3,M3)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR3=1,JCI3
            X(1)=H3(NR3,M3,1)
            DO I=2,22
              X(I)=X(1)
            END DO
            X(3)=H3(NR3,M3,2)
            X(8)=X(3)
            X(12)=X(3)
            X(15)=H3(NR3,M3,3)
            X(16)=X(3)
            X(17)=X(3)
            DO NL3=1,JCI3
              Y(1)=H3(NL3,M3,1)
              DO K=1,22
                  Z0(K,NL3,NR3,M3)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C     DO M4=1,MM4/MD4
      DO M4=1,MM4
        IF(.NOT.LINEAR)THEN
          DO NR4=1,JCI4
            X(1)=H4(NR4,M4,1)
            DO I=2,36
              X(I)=X(1)
            END DO
            X(7)=H4(NR4,M4,2)
            X(16)=X(7)
            X(23)=X(7)
            X(28)=X(7)
            X(32)=X(7)
            X(34)=X(7)
            DO NL4=1,JCI4
              Y(1)=H4(NL4,M4,1)
              DO I=2,36
                Y(I)=Y(1)
              END DO
              Y(8)=H4(NL4,M4,2)
              Y(17)=Y(8)
              Y(24)=Y(8)
              Y(29)=Y(8)
              Y(32)=Y(8)
              Y(33)=Y(8)
              DO K=1,I36
                W0(K,NL4,NR4,M4)=Y(37-K)*X(37-K)
              END DO
            END DO
          END DO
        ELSE
          DO NR4=1,JCI4
            X(1)=H4(NR4,M4,1)
            DO I=2,22
              X(I)=X(1)
            END DO
            X(4)=H4(NR4,M4,2)
            X(9)=X(4)
            X(13)=X(4)
            X(16)=X(4)
            X(18)=H4(NR4,M4,3)
            X(19)=X(4)
            DO NL4=1,JCI4
              Y(1)=H4(NL4,M4,1)
              DO K=1,22
                W0(K,NL4,NR4,M4)=Y(1)*X(K)
              END DO
            END DO
          END DO
        END IF
      END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)

      LSTART=1
      KLAN=LANMAX*(LANMAX+1)/2
      LSIZE=1
      ISKIP=0

555   CONTINUE
      IF(ISKIP.NE.0)THEN
        REWIND 64
        REWIND 74
        REWIND 84
      END IF
      LEND=LSTART-1
      KXA=0
      DO I=LSTART,NSIZE
        KXA=KXA+LSIZE
        IF(KXA.GT.KLAN)THEN
          KXA=KXA-LSIZE
          GO TO 554
        END IF
        LSIZE=LSIZE+1
        LEND=LEND+1
      END DO
554   CONTINUE
      IF(LEND.EQ.0)STOP 'LANMAX TOO SMALL'
C**ZEROISE MATRIX
      CALL DIAGL(KXA,XA5)
      REWIND 30

C**LOOP ROUND TAU (START 5-MODE INTEGRATION)
      ITAU=INIT-INCTAU
      DO MTAU=1,MMTAU/MDT
        ITAU=ITAU+INCTAU
CCCC    IF(ITAU.GT.362)ITAU=ITAU-360
        IF(ITAU.GT.722)ITAU=ITAU-720
 
C***********************************************************

        IF(ICOUPC.GT.0)THEN
          IF(J21.GT.1.AND.ICOUPC.GE.4)READ(64)VR
          IF(ICOUPC.GE.4)READ(84)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPC.GE.4)READ(64)VRR
          IF(ICOUPC.GE.4)READ(84)VCR
        END IF
        IF(JCOUPL.GT.0)THEN
          READ(74)VP
        ELSE
          READ(74)VPR
        END IF
 
C***********************************************************

C**********************************
CCCC    DO IRHS4=1,JCI4
CCCC      DO IRHS3=1,JCI3
CCCC        DO IRHS2=1,JCI2
CCCC          DO IRHS1=1,JCI1
CCCC            DO ILHS4=1,JCI4
CCCC              DO ILHS3=1,JCI3
CCCC                DO ILHS2=1,JCI2
CCCC                  DO ILHS1=1,JCI1
CCCC                    DO K=1,I36
CCCC      TEMP1(K,ILHS1,ILHS2,ILHS3,ILHS4,IRHS1,IRHS2,IRHS3,IRHS4)=0.D0
CCCC                    END DO
CCCC                  END DO
CCCC                END DO
CCCC              END DO
CCCC            END DO
CCCC          END DO
CCCC        END DO
CCCC      END DO
CCCC    END DO
C**********************************
C**START 4-MODE INTEGRATION
C       DO M1=1,MM1/MD1
        DO M1=1,MM1
          DO IRHS4=1,JCI4
            DO IRHS3=1,JCI3
              DO IRHS2=1,JCI2
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO ILHS2=1,JCI2
                      DO K=1,I36
                      TEMP2(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=0.D0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**START 3-MODE INTEGRATION
C         DO M2=1,MM2/MD2
          DO M2=1,MM2
            DO IRHS4=1,JCI4
              DO IRHS3=1,JCI3
                DO ILHS4=1,JCI4
                  DO ILHS3=1,JCI3
                    DO K=1,I36
                      TEMP3(K,ILHS3,ILHS4,IRHS3,IRHS4)=0.D0
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 2-MODE INTEGRATION
C           DO M3=1,MM3/MD3
            DO M3=1,MM3
              DO IRHS4=1,JCI4
                DO ILHS4=1,JCI4
                  DO K=1,I36
                    TEMP4(K,ILHS4,IRHS4)=0.D0
                  END DO
                END DO
              END DO
C**START 1-MODE INTEGRATION
C             DO M4=1,MM4/MD4
              DO M4=1,MM4
                DO I=1,36
                  C(I)=0
                END DO
                IF(JCOUPL.GT.0)THEN
                  IF(ICOUPC.GE.4)THEN
                    C(1)=VC(M4,M3,M2,M1,1)*IFACTC
                    C(2)=VC(M4,M3,M2,M1,1)*IFACTC
                    C(3)=VC(M4,M3,M2,M1,2)*IFACTC
                    C(4)=VC(M4,M3,M2,M1,2)*IFACTC
                    C(5)=VC(M4,M3,M2,M1,3)*IFACTC
                    C(6)=VC(M4,M3,M2,M1,3)*IFACTC
                    C(7)=VC(M4,M3,M2,M1,4)*IFACTC
                    C(8)=VC(M4,M3,M2,M1,4)*IFACTC
                    C(9)=VC(M4,M3,M2,M1,5)*IFACTC
                    C(10)=VC(M4,M3,M2,M1,5)*IFACTC
                    C(11)=VC(M4,M3,M2,M1,6)*IFACTC
                    C(12)=VC(M4,M3,M2,M1,7)*IFACTC
                    C(13)=VC(M4,M3,M2,M1,7)*IFACTC
                    C(14)=VC(M4,M3,M2,M1,8)*IFACTC
                    C(15)=VC(M4,M3,M2,M1,8)*IFACTC
                    C(16)=VC(M4,M3,M2,M1,9)*IFACTC
                    C(17)=VC(M4,M3,M2,M1,9)*IFACTC
                    C(18)=VC(M4,M3,M2,M1,10)*IFACTC
                    C(19)=VC(M4,M3,M2,M1,10)*IFACTC
                    C(20)=VC(M4,M3,M2,M1,11)*IFACTC
                    C(21)=VC(M4,M3,M2,M1,12)*IFACTC
                    C(22)=VC(M4,M3,M2,M1,12)*IFACTC
                    C(23)=VC(M4,M3,M2,M1,13)*IFACTC
                    C(24)=VC(M4,M3,M2,M1,13)*IFACTC
                    C(25)=VC(M4,M3,M2,M1,14)*IFACTC
                    C(26)=VC(M4,M3,M2,M1,14)*IFACTC
                    C(27)=VC(M4,M3,M2,M1,15)*IFACTC
                    C(28)=VC(M4,M3,M2,M1,16)*IFACTC
                    C(29)=VC(M4,M3,M2,M1,16)*IFACTC
                    C(30)=VC(M4,M3,M2,M1,17)*IFACTC
                    C(31)=VC(M4,M3,M2,M1,17)*IFACTC
                    C(32)=VC(M4,M3,M2,M1,18)*IFACTC
                    C(33)=VC(M4,M3,M2,M1,19)*IFACTC
                    C(34)=VC(M4,M3,M2,M1,19)*IFACTC
                    C(35)=VC(M4,M3,M2,M1,20)*IFACTC
                    C(36)=VC(M4,M3,M2,M1,21)*IFACTC
                    IF(J21.GT.1)C(36)=C(36)+VR(KROT,M4,M3,M2,M1)*IFACTC
                  END IF
                 C(36)=C(36)+VP(M4,M3,M2,M1)*IFACTL
                ELSE
                  IF(ICOUPC.GE.4)THEN
                    C(1)=VCR(M4,M3,M2,M1,1)*IFACTC
                    C(2)=VCR(M4,M3,M2,M1,1)*IFACTC
                    C(3)=VCR(M4,M3,M2,M1,2)*IFACTC
                    C(4)=VCR(M4,M3,M2,M1,2)*IFACTC
                    C(5)=VCR(M4,M3,M2,M1,3)*IFACTC
                    C(6)=VCR(M4,M3,M2,M1,3)*IFACTC
                    C(7)=VCR(M4,M3,M2,M1,4)*IFACTC
                    C(8)=VCR(M4,M3,M2,M1,4)*IFACTC
                    C(9)=VCR(M4,M3,M2,M1,5)*IFACTC
                    C(10)=VCR(M4,M3,M2,M1,5)*IFACTC
                    C(11)=VCR(M4,M3,M2,M1,6)*IFACTC
                    C(12)=VCR(M4,M3,M2,M1,7)*IFACTC
                    C(13)=VCR(M4,M3,M2,M1,7)*IFACTC
                    C(14)=VCR(M4,M3,M2,M1,8)*IFACTC
                    C(15)=VCR(M4,M3,M2,M1,8)*IFACTC
                    C(16)=VCR(M4,M3,M2,M1,9)*IFACTC
                    C(17)=VCR(M4,M3,M2,M1,9)*IFACTC
                    C(18)=VCR(M4,M3,M2,M1,10)*IFACTC
                    C(19)=VCR(M4,M3,M2,M1,10)*IFACTC
                    C(20)=VCR(M4,M3,M2,M1,11)*IFACTC
                    C(21)=VCR(M4,M3,M2,M1,12)*IFACTC
                    C(22)=VCR(M4,M3,M2,M1,12)*IFACTC
                    C(23)=VCR(M4,M3,M2,M1,13)*IFACTC
                    C(24)=VCR(M4,M3,M2,M1,13)*IFACTC
                    C(25)=VCR(M4,M3,M2,M1,14)*IFACTC
                    C(26)=VCR(M4,M3,M2,M1,14)*IFACTC
                    C(27)=VCR(M4,M3,M2,M1,15)*IFACTC
                    C(28)=VCR(M4,M3,M2,M1,16)*IFACTC
                    C(29)=VCR(M4,M3,M2,M1,16)*IFACTC
                    C(30)=VCR(M4,M3,M2,M1,17)*IFACTC
                    C(31)=VCR(M4,M3,M2,M1,17)*IFACTC
                    C(32)=VCR(M4,M3,M2,M1,18)*IFACTC
                    C(33)=VCR(M4,M3,M2,M1,19)*IFACTC
                    C(34)=VCR(M4,M3,M2,M1,19)*IFACTC
                    C(35)=VCR(M4,M3,M2,M1,20)*IFACTC
                    C(36)=VCR(M4,M3,M2,M1,21)*IFACTC
                    IF(J21.GT.1)C(36)=C(36)+VRR(KROT,M4,M3,M2,M1)
                  END IF
                  C(36)=C(36)+VPR(M4,M3,M2,M1)*IFACTL
                END IF
                DO IRHS4=1,JCI4
                  DO ILHS4=1,JCI4
                    DO K=1,I36
                      TEMP4(K,ILHS4,IRHS4)=TEMP4(K,ILHS4,IRHS4)+
     1                W0(K,ILHS4,IRHS4,M4)*C(37-K)
                    END DO
                  END DO
                END DO
              END DO
C**END 1-MODE INTEGRATION
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS4=1,JCI4
                    DO ILHS3=1,JCI3
                      DO K=1,I36
                        TEMP3(K,ILHS3,ILHS4,IRHS3,IRHS4)=
     1                  TEMP3(K,ILHS3,ILHS4,IRHS3,IRHS4)+
     2                  Z0(K,ILHS3,IRHS3,M3)*TEMP4(K,ILHS4,IRHS4)
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
                        DO K=1,I36
                          TEMP2(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=
     1                    TEMP2(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)+
     2                    Y0(K,ILHS2,IRHS2,M2)*
     3                    TEMP3(K,ILHS3,ILHS4,IRHS3,IRHS4)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
C**END 3-MODE INTEGRATION
C**********************************
CCCC      DO IRHS4=1,JCI4
CCCC        DO IRHS3=1,JCI3
CCCC          DO IRHS2=1,JCI2
CCCC            DO IRHS1=1,JCI1
CCCC              DO ILHS4=1,JCI4
CCCC                DO ILHS3=1,JCI3
CCCC                  DO ILHS2=1,JCI2
CCCC                    DO ILHS1=1,JCI1
CCCC                      DO K=1,I36
CCCC          TEMP1(K,ILHS1,ILHS2,ILHS3,ILHS4,IRHS1,IRHS2,IRHS3,IRHS4)=
CCCC 1        TEMP1(K,ILHS1,ILHS2,ILHS3,ILHS4,IRHS1,IRHS2,IRHS3,IRHS4)+
CCCC 2        X0(K,ILHS1,IRHS1,M1)*
CCCC 3        TEMP2(K,ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)
CCCC                      END DO
CCCC                    END DO
CCCC                  END DO
CCCC                END DO
CCCC              END DO
CCCC            END DO
CCCC          END DO
CCCC        END DO
CCCC      END DO
CCCC    END DO
C**********************************
 
C**NSIZE IS NO. UNIQUE INTEGRALS (4-DIM)
          J0=0
          DO IRHS=LSTART,LEND
            NR1=IP5(IRHS,1)
            NR2=IP5(IRHS,2)
            NR3=IP5(IRHS,3)
            NR4=IP5(IRHS,4)
            IRTAU=IP5(IRHS,5)
C           J0=IRHS*(IRHS-1)/2
            DO ILHS=1,IRHS
              NL1=IP5(ILHS,1)
              NL2=IP5(ILHS,2)
              NL3=IP5(ILHS,3)
              NL4=IP5(ILHS,4)
              ILTAU=IP5(ILHS,5)
              DO K=1,I36
CCCC            XA5(ILHS+J0)=XA5(ILHS+J0)+
CCCC 1          TEMP1(K,NL1,NL2,NL3,NL4,NR1,NR2,NR3,NR4)*
CCCC 1          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
                XA5(ILHS+J0)=XA5(ILHS+J0)+
     1          TEMP2(K,NL2,NL3,NL4,NR2,NR3,NR4)*X0(K,NL1,NR1,M1)*
     1          T0(K,ILTAU,IRTAU,MTAU)*DSTAU(ITAU)
              END DO
            END DO
            J0=J0+IRHS
          END DO
C**********************************
        END DO
C**END 4-MODE INTEGRATION
C**********************************
      END DO
C**END TAU LOOP (5-MODE INTEGRATION)

      WRITE(24)LSTART,LEND,KXA
      WRITE(30)LSTART,LEND,KXA
      CALL MTOUT(XA5,XRA5,KXA,24)
      CALL MTOUT(XA5,XRA5,KXA,30)
      IF(LEND.NE.NSIZE)THEN
        LSTART=LEND+1
        ISKIP=1
        GO TO 555
      END IF
      GO TO 7777
6666  CONTINUE
1111  READ(24)LSTART,LEND,KXA
      WRITE(30)LSTART,LEND,KXA
      CALL MTIN(XA5,XRA5,KXA,24)
      CALL MTOUT(XA5,XRA5,KXA,30)
      IF(LEND.NE.NSIZE)GO TO 1111
7777  REWIND 30
C**INITIAL READ
      READ(30)LSTART,LEND,KXA
      CALL MTIN(XA5,XRA5,KXA,30)

      IF(IND.NE.0)THEN
        DO IR=1,NSIZE
          NR=0
          DO J=1,NC1
            J1=JC1(J)
            NR1=IP5(IR,J1)
            NR=NR+1
            NX5(NR)=NR1
          END DO
          DO J=1,NC2
            J2=JC2(J)
            NR2=IP5(IR,J2)
            NR=NR+1
            NX5(NR)=NR2
          END DO
          IRL(NX5(1),NX5(2),NX5(3),NX5(4),NX5(5))=IR
        END DO
      END IF

C***********************************ALGORITHM FROM VCV4

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFFR=IOFFR+NCSIZE(I)
      END DO
      J0R=IOFFR*(IOFFR+1)/2
C**TEMPORARY-LAN
      JSR=1
C**TEMPORARY-LAN
C     DO 9999 ISMR1=1,NVSYM
      DO 9999 ISMR1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE

        IF(NONC1.GT.NONC2)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          KCCCL=KCSZL
CC        IF(NONC1.GT.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             KCCCL=KCSZL
C           ELSE
C             KCCCL=IRH1
C           END IF
CC        END IF
          DO ILH1=1,KCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=IROFF
            NONZC1(2,NONZ1)=ILOFF
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
        IF(NONZ1.EQ.0)GO TO 5556
        NONZC1(3,NONZ1+1)=0
        END IF

      END IF

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M AND/OR TAU)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE

        IF(NONC1.LE.NONC2)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M AND/OR TAU)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          LCCCL=LCSZL
CC        IF(NONC1.LE.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             LCCCL=LCSZL
C           ELSE
C             LCCCL=IRH2
C           END IF
CC        END IF
          DO ILH2=1,LCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=JROFF
            NONZC2(2,NONZ2)=JLOFF
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
        IF(NONZ2.EQ.0)GO TO 5555
        NONZC2(3,NONZ2+1)=0
        END IF

      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM4B.EQ.0)THEN
C       WRITE(IOUT,*)'NKVALL,NKVALR = ',NKVALL,NKVALR
C       WRITE(IOUT,*)'NLVALL,NLVALR = ',NLVALL,NLVALR
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C       WRITE(IOUT,*)'NONC1,NONC2 = ',NONC1,NONC2
C     END IF

      IF(NONC1.LE.NONC2)GO TO 3003
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR/ M AND/OR TAU)

      DO IRH2=1,LCSZR
        JROFF=IRH2+LCOFFR
        DO I=1,NKVALR
          DO J=1,NKVALL
            DO K=1,NLVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          LCCCL=LCSZL
        ELSE
          LCCCL=IRH2
        END IF
        DO ILH2=1,LCCCL
          JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
          DO K=1,NKMOD2
            IGOT=0
            DO I=1,NC2
              IF(K.EQ.MC2(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 6001
          END DO
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH2.EQ.IRH2)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NKVALR
            DO J=1,NKVALL
              XCON2(J,I)=0
            END DO
          END DO
C
C**CONTRACTION SCHEME  '1' (K)
C
          IPREV=0
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            INEXT=NONZC1(3,INH1+1)
            ILH1=NONZC1(4,INH1)
            IROFF=NONZC1(1,INH1)
            ILOFF=NONZC1(2,INH1)
            IF(IPREV.NE.IRH1)THEN
C**NEW RHS
              DO J=1,NKVALL
                XK1(J,1)=0
              END DO
            END IF

            DO I=1,NC1
              NCR1(I)=IPC1(IROFF,MC1(I))
              NCL1(I)=IPC1(ILOFF,MC1(I))
            END DO

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX5(NR)=NR1
        MX5(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX5(NR)=NR2
        MX5(NR)=NL2
      END DO
      IR=IRL(NX5(1),NX5(2),NX5(3),NX5(4),NX5(5))
      IL=IRL(MX5(1),MX5(2),MX5(3),MX5(4),MX5(5))
C**TEMPORARY

C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XXX=XA5(I)

            DO K=1,NKVALL
              XK1(K,1)=XK1(K,1)+XXX*CFS1(ILH1,K,KEL,KSML)
            END DO
5001        CONTINUE

            IF(IRH1.NE.INEXT)THEN

            DO K=1,NKVALR
              DO L=1,NKVALL
                XCON2(L,K)=XCON2(L,K)+XK1(L,1)*CFS1(IRH1,K,KEL,KSMR)
              END DO
            END DO

            END IF
            IPREV=IRH1
          END DO
CC
          DO I=1,NKVALR
            DO J=1,NKVALL
              DO K=1,NLVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6001      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILL,IKR,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKR,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILR,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILR,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO
C*************************************

      GO TO 5555

3003  CONTINUE
C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '1' (K)

      DO IRH1=1,KCSZR
        IROFF=IRH1+KCOFFR
        DO I=1,NLVALR
          DO J=1,NLVALL
            DO K=1,NKVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          KCCCL=KCSZL
        ELSE
          KCCCL=IRH1
        END IF
        DO ILH1=1,KCCCL
          ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
          DO K=1,NKMOD1
            IGOT=0
            DO I=1,NC1
              IF(K.EQ.MC1(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 6002
          END DO
          DO I=1,NC1
            NCR1(I)=IPC1(IROFF,MC1(I))
            NCL1(I)=IPC1(ILOFF,MC1(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH1.EQ.IRH1)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NLVALR
            DO J=1,NLVALL
              XCON2(J,I)=0
            END DO
          END DO

C**CONTRACTION SCHEME  '2' (L AND/OR N AND/OR M AND/OR J)

          JPREV=0
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            JNEXT=NONZC2(3,INH2+1)
            ILH2=NONZC2(4,INH2)
            JROFF=NONZC2(1,INH2)
            JLOFF=NONZC2(2,INH2)
            IF(JPREV.NE.IRH2)THEN
C**NEW RHS
              DO J=1,NLVALL
                XK1(J,1)=0
              END DO
            END IF

            DO I=1,NC2
              NCR2(I)=IPC2(JROFF,MC2(I))
              NCL2(I)=IPC2(JLOFF,MC2(I))
            END DO

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX5(NR)=NR1
        MX5(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX5(NR)=NR2
        MX5(NR)=NL2
      END DO
      IR=IRL(NX5(1),NX5(2),NX5(3),NX5(4),NX5(5))
      IL=IRL(MX5(1),MX5(2),MX5(3),MX5(4),MX5(5))
C**TEMPORARY

C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XXX=XA5(I)

            DO J=1,NLVALL
              XK1(J,1)=XK1(J,1)+XXX*CFS2(ILH2,J,KEL,LSML)
            END DO
5002        CONTINUE

            IF(IRH2.NE.JNEXT)THEN

            DO I=1,NLVALR
              DO J=1,NLVALL
                XCON2(J,I)=XCON2(J,I)+XK1(J,1)*CFS2(IRH2,I,KEL,LSMR)
              END DO
            END DO

            END IF
            JPREV=IRH2
          END DO
CC
          DO I=1,NLVALR
            DO J=1,NLVALL
              DO K=1,NKVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6002      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKR,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKL,ILR,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILR,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO

5555  CONTINUE
5556  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM4B.EQ.0)THEN
        IF(IND.NE.0)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM4B=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI5A(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,
     1MODE2,MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,
     2NN1,MM1,NN2,MM2,NN3,MM3,NN4,MM4,NN5,MM5,IP,ISIZMX,IPC,ICSIZE,
     3IPSIZE,KCONT,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,NSIZE,XK,TEMP,XCON,NVAL,
     4ISTART,IEND,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,JCI1,JCI2,JCI3,JCI4,
     5JCI5,X0,Y0,Z0,W0,U0,VP,VPR,MODINT,KEL,NS,
     6CFS,ISIZXX,NVALX,KEL21,NVSYMX,
C**ANALYTIC
     7XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM)
C**ANALYTIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM5,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H4(NN4,MM4,3),H5(NN5,MM5,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQ5(MM5)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP5(ISIZE5,5),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(JCI5,JCI5),TEMP2(JCI4,JCI5,JCI4,JCI5)
      DIMENSION TEMP3(JCI3,JCI4,JCI5,JCI3,JCI4,JCI5)
C     DIMENSION TEMP4(JCI2*JCI3*JCI4*JCI5,JCI2*JCI3*JCI4*JCI5)
C     DIMENSION TEMP5(JCI2,JCI3,JCI4,JCI5)
      DIMENSION X0(JCI1,JCI1,MM1),Y0(JCI2,JCI2,MM2)
      DIMENSION Z0(JCI3,JCI3,MM3),W0(JCI4,JCI4,MM4)
      DIMENSION U0(JCI5,JCI5,MM5)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA5(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP5(NTOT5),CP5(MAX5,NTOT5),JP5(MAX5,NTOT5,5)
      DIMENSION INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
C---TEMP ARRAY INITIALIZE---ACTUALLY THIS CAN BE DONE SOMEWHERE ELSE
C---BEFORE THIS SUB GET CALLED, THEN TEMP4 AND TEMP CAN BE CALLED.
C     I1=0
C     DO I2=1,JCI2
C      DO I3=1,JCI3
C       DO I4=1,JCI4
C        DO I5=1,JCI5
C         I1=I1+1
C         TEMP5(I2,I3,I4,I5)=I1
C        END DO
C       END DO
C      END DO
C     END DO
C---END OF TEMP ARRAY INITIALIZATION
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      MD5=MODINT(MOD5)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(MOD5.EQ.ISYM(I,J))N5=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N1.EQ.N5)MD5=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N2.EQ.N5)MD5=1
      IF(N3.EQ.N4)MD4=1
      IF(N3.EQ.N5)MD5=1
      IF(N4.EQ.N5)MD5=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      IF(N12.EQ.N5)MD5=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      IF(N13.EQ.N5)MD5=1
      N14=ISYMP(N1,N4)
      IF(N14.EQ.N5)MD5=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      IF(N23.EQ.N5)MD5=1
      N24=ISYMP(N2,N4)
      IF(N24.EQ.N5)MD5=1
      N34=ISYMP(N3,N4)
      IF(N34.EQ.N5)MD5=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      IF(N123.EQ.N5)MD5=1
      N124=ISYMP(N12,N4)
      IF(N124.EQ.N5)MD5=1
      N134=ISYMP(N13,N4)
      IF(N134.EQ.N5)MD5=1
      N234=ISYMP(N23,N4)
      IF(N234.EQ.N5)MD5=1
      N1234=ISYMP(N12,N34)
      IF(N1234.EQ.N5)MD5=1
      CALL VDCI5A(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,MODE2,
     1MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,NN1,MM1,
     2MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,NN5,MM5,
     3MM5/MD5,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP5,ISIZE5,
     4XA5,XRA5,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,TEMP3,
     5TEMP4,TEMP5,JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,W0,U0,VP,VPR,
     6MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,
C**ANALYTIC
     7XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM)
C**ANALYTIC
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI5A(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,
     1MODE2,MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,
     2NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,NN4,MH4,MM4,NN5,MH5,MM5,
     3IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,
     4NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,TEMP3,TEMP4,
     5TEMP5,JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,W0,U0,VP,VPR,
     6MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,
C**ANALYTIC
     7XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM)
C**ANALYTIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM5,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION H4(NN4,MH4,3),H5(NN5,MH5,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQ5(MM5)
      DIMENSION IP(ISIZMX,NNMODE)
      DIMENSION IPC(IPSIZE,1),IP5(ISIZE5,5),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(JCI5,JCI5),TEMP2(JCI4,JCI5,JCI4,JCI5)
      DIMENSION TEMP3(JCI3,JCI4,JCI5,JCI3,JCI4,JCI5)
C     DIMENSION TEMP4(JCI2*JCI3*JCI4*JCI5,JCI2*JCI3*JCI4*JCI5)
C     DIMENSION TEMP5(JCI2,JCI3,JCI4,JCI5)
      DIMENSION X0(JCI1,JCI1,MM1),Y0(JCI2,JCI2,MM2)
      DIMENSION Z0(JCI3,JCI3,MM3),W0(JCI4,JCI4,MM4)
      DIMENSION U0(JCI5,JCI5,MM5)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA5(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP5(NTOT5),CP5(MAX5,NTOT5),JP5(MAX5,NTOT5,5)
      DIMENSION INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM5A.EQ.1)THEN
        WRITE(IOUT,*)'Calculating VCCI5A'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

C** FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1', 'MOD2',
C** 'MOD3', MOD4' AND 'MOD5'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C**ANALYTIC
      IND=INDM(MOD1)+INDN(MOD2)+INDL(MOD3)+INDK(MOD4)+MOD5
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI5
C**TEMPORARY
      IF(ICONDP.NE.0)THEN
        IF(MOLINC.GT.0)GO TO 99999
        GO TO 6666
      END IF
C**TEMPORARY

      IFACTL=INTFAC(NMODE,ICOUPL,5)
      IF(MOLINC.GT.0)IFACTL=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(75)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(75)VPR
      END IF

      IF(MOLINC.GT.0)GO TO 99999

C***********************************************************

      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      MD5=MODINT(MOD5)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(MOD5.EQ.ISYM(I,J))N5=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N1.EQ.N5)MD5=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N2.EQ.N5)MD5=1
      IF(N3.EQ.N4)MD4=1
      IF(N3.EQ.N5)MD5=1
      IF(N4.EQ.N5)MD5=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      IF(N12.EQ.N5)MD5=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      IF(N13.EQ.N5)MD5=1
      N14=ISYMP(N1,N4)
      IF(N14.EQ.N5)MD5=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      IF(N23.EQ.N5)MD5=1
      N24=ISYMP(N2,N4)
      IF(N24.EQ.N5)MD5=1
      N34=ISYMP(N3,N4)
      IF(N34.EQ.N5)MD5=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      IF(N123.EQ.N5)MD5=1
      N124=ISYMP(N12,N4)
      IF(N124.EQ.N5)MD5=1
      N134=ISYMP(N13,N4)
      IF(N134.EQ.N5)MD5=1
      N234=ISYMP(N23,N4)
      IF(N234.EQ.N5)MD5=1
      N1234=ISYMP(N12,N34)
      IF(N1234.EQ.N5)MD5=1
      MD=MD1*MD2*MD3*MD4*MD5

C***********************************************************
      IF(IWHICH.GE.0)THEN
C**POTENTIAL
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
        DO M1=1,MM1
          DO NR1=1,JCI1
            X=H1(NR1,M1,1)*IFACTL*MD
            DO NL1=1,JCI1
              Y=H1(NL1,M1,1)
              X0(NL1,NR1,M1)=Y*X
            END DO
          END DO
        END DO
        DO M2=1,MM2
          DO NR2=1,JCI2
            X=H2(NR2,M2,1)
            DO NL2=1,JCI2
              Y=H2(NL2,M2,1)
              Y0(NL2,NR2,M2)=Y*X
            END DO
          END DO
        END DO
        DO M3=1,MM3
          DO NR3=1,JCI3
            X=H3(NR3,M3,1)
            DO NL3=1,JCI3
              Y=H3(NL3,M3,1)
              Z0(NL3,NR3,M3)=Y*X
            END DO
          END DO
        END DO
        DO M4=1,MM4
          DO NR4=1,JCI4
            X=H4(NR4,M4,1)
            DO NL4=1,JCI4
              Y=H4(NL4,M4,1)
              W0(NL4,NR4,M4)=Y*X
            END DO
          END DO
        END DO
        DO M5=1,MM5
          DO NR5=1,JCI5
            X=H5(NR5,M5,1)
            DO NL5=1,JCI5
              Y=H5(NL5,M5,1)
              U0(NL5,NR5,M5)=Y*X
            END DO
          END DO
        END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**START 5-MODE INTEGRATION
        DO M1=1,MM1
C**START 4-MODE INTEGRATION
          DO M2=1,MM2
            DO IRHS5=1,JCI5
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS5=1,JCI5
                    DO ILHS4=1,JCI4
                      DO ILHS3=1,JCI3
                        TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)=0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 3-MODE INTEGRATION
            DO M3=1,MM3
              DO IRHS5=1,JCI5
                DO IRHS4=1,JCI4
                  DO ILHS5=1,JCI5
                    DO ILHS4=1,JCI4
                      TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)=0
                    END DO
                  END DO
                END DO
              END DO
C**START 2-MODE INTEGRATION
              DO M4=1,MM4
                DO IRHS5=1,JCI5
                  DO ILHS5=1,JCI5
                    TEMP1(ILHS5,IRHS5)=0
                  END DO
                END DO
C**START 1-MODE INTEGRATION
                DO M5=1,MM5
                  IF(JCOUPL.GT.0)THEN
                    IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1              C=VP(M5,M4,M3,M2,M1)
                  ELSE
                    IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1              C=VPR(M5,M4,M3,M2,M1)
                  END IF
                  DO IRHS5=1,JCI5
                    DO ILHS5=1,JCI5
                      TEMP1(ILHS5,IRHS5)=TEMP1(ILHS5,IRHS5)+
     1                U0(ILHS5,IRHS5,M5)*C
                    END DO
                  END DO
                END DO
C**END 1-MODE INTEGRATION
                DO IRHS5=1,JCI5
                  DO IRHS4=1,JCI4
                    DO ILHS5=1,JCI5
                      DO ILHS4=1,JCI4
                        TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)=
     1                  TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)+
     2                  W0(ILHS4,IRHS4,M4)*TEMP1(ILHS5,IRHS5)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
C**END 2-MODE INTEGRATION
              DO IRHS5=1,JCI5
                DO IRHS4=1,JCI4
                  DO IRHS3=1,JCI3
                    DO ILHS5=1,JCI5
                      DO ILHS4=1,JCI4
                        DO ILHS3=1,JCI3
                          TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)=
     1                    TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)+
     2                    Z0(ILHS3,IRHS3,M3)*
     3                    TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 3-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (5-DIM)
            DO IRHS=1,NSIZE
              NR1=IP5(IRHS,1)
              NR2=IP5(IRHS,2)
              NR3=IP5(IRHS,3)
              NR4=IP5(IRHS,4)
              NR5=IP5(IRHS,5)
              J0=IRHS*(IRHS-1)/2
              DO ILHS=1,IRHS
                NL1=IP5(ILHS,1)
                NL2=IP5(ILHS,2)
                NL3=IP5(ILHS,3)
                NL4=IP5(ILHS,4)
                NL5=IP5(ILHS,5)
                XA5(ILHS+J0)=XA5(ILHS+J0)+
     1          TEMP3(NL3,NL4,NL5,NR3,NR4,NR5)*X0(NL1,NR1,M1)*
     2          Y0(NL2,NR2,M2)
              END DO
            END DO
          END DO
C**END 4-MODE INTEGRATION
        END DO
C**END 5-MODE INTEGRATION
      END IF
      CALL MATOUT(XA5,XRA5,NSIZE,25)
      GO TO 7777
6666  CALL MATIN(XA5,XRA5,NSIZE,25)
7777  CONTINUE

C     IF(MOLINC.GT.0)RETURN

C***********************************ALGORITHM FROM VCI5

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA)
      IOFF=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFF=IOFF+NCSIZE(I)
      END DO
      J0R=IOFF*(IOFF+1)/2
C**TEMPORARY-LAN
      JS=1
C**TEMPORARY-LAN
C     DO 9999 ISM1=1,NVSYM
      DO 9999 ISM1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        ICOFF=ICOFF1
        ICSZ=ICSZ1
        NKVAL=NCVAL(1,ISM1)
      ELSE
        ISM=ISM2
        ICOFF=ICOFF2
        ICSZ=ICSZ2
        NKVAL=NCVAL(2,ISM2)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZ
        IROFF=IRHS+ICOFF
        NR1=IPC(IROFF,MODE1)
        NR2=IPC(IROFF,MODE2)
        NR3=IPC(IROFF,MODE3)
        NR4=IPC(IROFF,MODE4)
        NR5=IPC(IROFF,MODE5)
C**FIND RHS INDEX
        DO IR=1,NSIZE
          IF(NR1.EQ.IP5(IR,1).AND.NR2.EQ.IP5(IR,2).AND.
     1       NR3.EQ.IP5(IR,3).AND.NR4.EQ.IP5(IR,4).AND.
     2       NR5.EQ.IP5(IR,5))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,IRHS
          ILOFF=ILHS+ICOFF
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1     MODE4.AND.K.NE.MODE5.AND.(IPC(IROFF,K).NE.IPC(ILOFF,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPC(ILOFF,MODE1)
          NL2=IPC(ILOFF,MODE2)
          NL3=IPC(ILOFF,MODE3)
          NL4=IPC(ILOFF,MODE4)
          NL5=IPC(ILOFF,MODE5)
C**FIND LHS INDEX
          DO IL=1,NSIZE
            IF(NL1.EQ.IP5(IL,1).AND.NL2.EQ.IP5(IL,2).AND.
     1         NL3.EQ.IP5(IL,3).AND.NL4.EQ.IP5(IL,4).AND.
     2         NL5.EQ.IP5(IL,5))GO TO 2000
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
          XYZ=XA5(I)
          IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
            XYZ=0
            DO I=1,NP5(IND)
              K=JP5(I,IND,1)+1
              L=JP5(I,IND,2)+1
              N=JP5(I,IND,3)+1
              M=JP5(I,IND,4)+1
              J=JP5(I,IND,5)+1
              XYZ=XYZ+CP5(I,IND)*XKAN(NL1,NR1,K,MOD1)*
     1        XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2        XKAN(NL4,NR4,M,MOD4)*XKAN(NL5,NR5,J,MOD5)
            END DO
C**ANALYTIC
          END IF
3000      CONTINUE
          XK(ILHS,IRHS)=(XYZ+ZYX)*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZ,NKVAL,ICSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVAL,NKVAL,ICSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      DO IRHS=JSTART,JEND
        IROFF=IRHS+IOFF
C**TEMPORARY-LAN
        J0=IROFF*(IROFF-1)/2+IOFF-J0R
C**TEMPORARY-LAN
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IP(IROFF,KCONT)
        DO ILHS=1,IRHS
          ILOFF=ILHS+IOFF
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IP(IROFF,MCONT).NE.IP(ILOFF,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IP(ILOFF,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+J0)=XA(ILHS+J0)+XYZ
4000      CONTINUE
        END DO
      END DO

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
99999 CONTINUE

      IF(ITIM5A.EQ.1)THEN
        CALL TIMIT(3)
        ITIM5A=2
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI5B(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,
     1MODE2,MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,
     2NN1,MM1,NN2,MM2,NN3,MM3,NN4,MM4,NN5,MM5,IP,ISIZMX,IPC1,ICSIZ1,
     3IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP5,ISIZE5,XA5,
     4XRA5,NSIZE,XK1,XK2,IRL,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,
     5TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,
     6W0,U0,VP,VPR,MODINT,KEL,NS,CFS1,CFS2,
     7ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,
C**ANALYTIC
     8XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM,
C**ANALYTIC
     9MAXVAL,MAXQ5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM5,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H4(NN4,MM4,3),H5(NN5,MM5,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP5(ISIZE5,5),IRL(MAXQ5,MAXQ5,MAXQ5,MAXQ5,MAXQ5)
      DIMENSION TEMP1(JCI5,JCI5),TEMP2(JCI4,JCI5,JCI4,JCI5)
      DIMENSION TEMP3(JCI3,JCI4,JCI5,JCI3,JCI4,JCI5)
C     DIMENSION TEMP4(JCI2*JCI3*JCI4*JCI5,JCI2*JCI3*JCI4*JCI5)
C     DIMENSION TEMP5(JCI2,JCI3,JCI4,JCI5)
      DIMENSION X0(JCI1,JCI1,MM1),Y0(JCI2,JCI2,MM2)
      DIMENSION Z0(JCI3,JCI3,MM3),W0(JCI4,JCI4,MM4)
      DIMENSION U0(JCI5,JCI5,MM5)
C******************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA5(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C******************************************************************
C**ANALYTIC
      DIMENSION NP5(NTOT5),CP5(MAX5,NTOT5),JP5(MAX5,NTOT5,5)
      DIMENSION INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
C---TEMP ARRAY INITIALIZE---ACTUALLY THIS CAN BE DONE SOMEWHERE ELSE
C---BEFORE THIS SUB GET CALLED, THEN TEMP4 AND TEMP CAN BE CALLED.
C     I1=0
C     DO I2=1,JCI2
C      DO I3=1,JCI3
C       DO I4=1,JCI4
C        DO I5=1,JCI5
C         I1=I1+1
C         TEMP5(I2,I3,I4,I5)=I1
C        END DO
C       END DO
C      END DO
C     END DO
C---END OF TEMP ARRAY INITIALIZATION
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      MD5=MODINT(MOD5)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(MOD5.EQ.ISYM(I,J))N5=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N1.EQ.N5)MD5=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N2.EQ.N5)MD5=1
      IF(N3.EQ.N4)MD4=1
      IF(N3.EQ.N5)MD5=1
      IF(N4.EQ.N5)MD5=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      IF(N12.EQ.N5)MD5=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      IF(N13.EQ.N5)MD5=1
      N14=ISYMP(N1,N4)
      IF(N14.EQ.N5)MD5=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      IF(N23.EQ.N5)MD5=1
      N24=ISYMP(N2,N4)
      IF(N24.EQ.N5)MD5=1
      N34=ISYMP(N3,N4)
      IF(N34.EQ.N5)MD5=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      IF(N123.EQ.N5)MD5=1
      N124=ISYMP(N12,N4)
      IF(N124.EQ.N5)MD5=1
      N134=ISYMP(N13,N4)
      IF(N134.EQ.N5)MD5=1
      N234=ISYMP(N23,N4)
      IF(N234.EQ.N5)MD5=1
      N1234=ISYMP(N12,N34)
      IF(N1234.EQ.N5)MD5=1
      CALL VDCI5B(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,MODE2,
     1MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,NN1,MM1,
     2MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,NN5,MM5,
     3MM5/MD5,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,
     4KCONT2,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,NSIZE,XK1,XK2,IRL,KTEMP,
     5XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,
     6JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,W0,U0,VP,VPR,
     6MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,
     7NONZC1,NONZC2,
C**ANALYTIC
     8XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM,
C**ANALYTIC
     9MAXVAL,MAXQ5)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI5B(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,
     1MODE2,MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,
     2NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,NN4,MH4,MM4,NN5,MH5,MM5,
     3IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,
     3KCONT2,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,NSIZE,XK1,XK2,IRL,KTEMP,
     4XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,
     5JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,W0,U0,VP,VPR,
     6MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,
     7NONZC1,NONZC2,
C**ANALYTIC
     8XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM,
C**ANALYTIC
     9MAXVAL,MAXQ5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM5,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)

C**TEMPORARY
      DIMENSION IRL(MAXQ5,MAXQ5,MAXQ5,MAXQ5,MAXQ5),NX5(5),MX5(5)
C**TEMPORARY

      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION H4(NN4,MH4,3),H5(NN5,MH5,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQ5(MM5)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP5(ISIZE5,5)
      DIMENSION TEMP1(JCI5,JCI5),TEMP2(JCI4,JCI5,JCI4,JCI5)
      DIMENSION TEMP3(JCI3,JCI4,JCI5,JCI3,JCI4,JCI5)
C     DIMENSION TEMP4(JCI2*JCI3*JCI4*JCI5,JCI2*JCI3*JCI4*JCI5)
C     DIMENSION TEMP5(JCI2,JCI3,JCI4,JCI5)
      DIMENSION X0(JCI1,JCI1,MM1),Y0(JCI2,JCI2,MM2)
      DIMENSION Z0(JCI3,JCI3,MM3),W0(JCI4,JCI4,MM4)
      DIMENSION U0(JCI5,JCI5,MM5)
C******************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA5(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,MAXVAL,MAXVAL,2)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C******************************************************************
C**ANALYTIC
      DIMENSION NP5(NTOT5),CP5(MAX5,NTOT5),JP5(MAX5,NTOT5,5)
      DIMENSION INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM5A.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCI5B'
          CALL FLUSH(IOUT)
          CALL TIMIT(1)
        END IF
      END IF

C**************************************************************
C**************************************************************
C**BOTH CONTRACTION SCHEMES HAVE AT LEAST ONE MODE
C**CHECK SCHEME '1'
      I1=KCONT1
      NC1=0
      DO NN=1,ICONT(I1)
        IF(MOD1.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE1
          JC1(NC1)=1
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD2.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE2
          JC1(NC1)=2
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD3.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE3
          JC1(NC1)=3
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD4.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE4
          JC1(NC1)=4
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD5.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE5
          JC1(NC1)=5
        END IF
      END DO
C**CHECK SCHEME '2'
      I2=KCONT2
      NC2=0
      DO NN=1,ICONT(I2)
        IF(MOD1.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE1
          JC2(NC2)=1
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD2.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE2
          JC2(NC2)=2
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD3.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE3
          JC2(NC2)=3
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD4.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE4
          JC2(NC2)=4
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD5.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE5
          JC2(NC2)=5
        END IF
      END DO
C**NKMOD1 IS TOTAL NUMBER OF MODES FOR SCHEME INVOLVING 'MOD1' (K)
C**NKMOD2 IS TOTAL NUMBER OF MODES FOR THE OTHER SCHEME (ONE OR ALL
C**OF 'MOD2' (L), 'MOD3' (N), 'MOD4' (M) AND 'MOD5' (J))
C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L),
C**'MOD3' (N) OR 'MOD4' (M),'MOD5' (J)
C**IF THEY ARE IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
C**************************************************************
      IF(IND.NE.0)GO TO 7777
C**************************************************************

C**FIND MAX QUANTUM
      MAXQ5=0
      DO IR=1,NSIZE
        NR1=IP5(IR,1)
        NR2=IP5(IR,2)
        NR3=IP5(IR,3)
        NR4=IP5(IR,4)
        NR5=IP5(IR,5)
        IF(NR1.GT.MAXQ5)MAXQ5=NR1
        IF(NR2.GT.MAXQ5)MAXQ5=NR2
        IF(NR3.GT.MAXQ5)MAXQ5=NR3
        IF(NR4.GT.MAXQ5)MAXQ5=NR4
        IF(NR5.GT.MAXQ5)MAXQ5=NR5
      END DO

C**ANALYTIC
      IND5=INDM(MOD1)+INDN(MOD2)+INDL(MOD3)+INDK(MOD4)+MOD5
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI5
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=INTFAC(NMODE,ICOUPL,5)
      IF(MOLINC.GT.0)IFACTL=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(75)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(75)VPR
      END IF

C***********************************************************

      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      MD5=MODINT(MOD5)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(MOD5.EQ.ISYM(I,J))N5=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N1.EQ.N5)MD5=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N2.EQ.N5)MD5=1
      IF(N3.EQ.N4)MD4=1
      IF(N3.EQ.N5)MD5=1
      IF(N4.EQ.N5)MD5=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      IF(N12.EQ.N5)MD5=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      IF(N13.EQ.N5)MD5=1
      N14=ISYMP(N1,N4)
      IF(N14.EQ.N5)MD5=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      IF(N23.EQ.N5)MD5=1
      N24=ISYMP(N2,N4)
      IF(N24.EQ.N5)MD5=1
      N34=ISYMP(N3,N4)
      IF(N34.EQ.N5)MD5=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      IF(N123.EQ.N5)MD5=1
      N124=ISYMP(N12,N4)
      IF(N124.EQ.N5)MD5=1
      N134=ISYMP(N13,N4)
      IF(N134.EQ.N5)MD5=1
      N234=ISYMP(N23,N4)
      IF(N234.EQ.N5)MD5=1
      N1234=ISYMP(N12,N34)
      IF(N1234.EQ.N5)MD5=1
      MD=MD1*MD2*MD3*MD4*MD5

C***********************************************************
      IF(IWHICH.GE.0)THEN
C**POTENTIAL
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
        DO M1=1,MM1
          DO NR1=1,JCI1
            X=H1(NR1,M1,1)*IFACTL*MD
            DO NL1=1,JCI1
              Y=H1(NL1,M1,1)
              X0(NL1,NR1,M1)=Y*X
            END DO
          END DO
        END DO
        DO M2=1,MM2
          DO NR2=1,JCI2
            X=H2(NR2,M2,1)
            DO NL2=1,JCI2
              Y=H2(NL2,M2,1)
              Y0(NL2,NR2,M2)=Y*X
            END DO
          END DO
        END DO
        DO M3=1,MM3
          DO NR3=1,JCI3
            X=H3(NR3,M3,1)
            DO NL3=1,JCI3
              Y=H3(NL3,M3,1)
              Z0(NL3,NR3,M3)=Y*X
            END DO
          END DO
        END DO
        DO M4=1,MM4
          DO NR4=1,JCI4
            X=H4(NR4,M4,1)
            DO NL4=1,JCI4
              Y=H4(NL4,M4,1)
              W0(NL4,NR4,M4)=Y*X
            END DO
          END DO
        END DO
        DO M5=1,MM5
          DO NR5=1,JCI5
            X=H5(NR5,M5,1)
            DO NL5=1,JCI5
              Y=H5(NL5,M5,1)
              U0(NL5,NR5,M5)=Y*X
            END DO
          END DO
        END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**START 5-MODE INTEGRATION
        DO M1=1,MM1
C**START 4-MODE INTEGRATION
          DO M2=1,MM2
            DO IRHS5=1,JCI5
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS5=1,JCI5
                    DO ILHS4=1,JCI4
                      DO ILHS3=1,JCI3
                        TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)=0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 3-MODE INTEGRATION
            DO M3=1,MM3
              DO IRHS5=1,JCI5
                DO IRHS4=1,JCI4
                  DO ILHS5=1,JCI5
                    DO ILHS4=1,JCI4
                      TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)=0
                    END DO
                  END DO
                END DO
              END DO
C**START 2-MODE INTEGRATION
              DO M4=1,MM4
                DO IRHS5=1,JCI5
                  DO ILHS5=1,JCI5
                    TEMP1(ILHS5,IRHS5)=0
                  END DO
                END DO
C**START 1-MODE INTEGRATION
                DO M5=1,MM5
                  IF(JCOUPL.GT.0)THEN
                    IF(IWHICH.GE.0.OR.MOLINC.LE.0)C=VP(M5,M4,M3,M2,M1)
                  ELSE
                    IF(IWHICH.GE.0.OR.MOLINC.LE.0)C=VPR(M5,M4,M3,M2,M1)
                  END IF
                  DO IRHS5=1,JCI5
                    DO ILHS5=1,JCI5
                      TEMP1(ILHS5,IRHS5)=TEMP1(ILHS5,IRHS5)+
     1                U0(ILHS5,IRHS5,M5)*C
                    END DO
                  END DO
                END DO
C**END 1-MODE INTEGRATION
                DO IRHS5=1,JCI5
                  DO IRHS4=1,JCI4
                    DO ILHS5=1,JCI5
                      DO ILHS4=1,JCI4
                        TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)=
     1                  TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)+
     2                  W0(ILHS4,IRHS4,M4)*TEMP1(ILHS5,IRHS5)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
C**END 2-MODE INTEGRATION
              DO IRHS5=1,JCI5
                DO IRHS4=1,JCI4
                  DO IRHS3=1,JCI3
                    DO ILHS5=1,JCI5
                      DO ILHS4=1,JCI4
                        DO ILHS3=1,JCI3
                          TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)=
     1                    TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)+
     2                    Z0(ILHS3,IRHS3,M3)*
     3                    TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 3-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (5-DIM)
            DO IRHS=1,NSIZE
              NR1=IP5(IRHS,1)
              NR2=IP5(IRHS,2)
              NR3=IP5(IRHS,3)
              NR4=IP5(IRHS,4)
              NR5=IP5(IRHS,5)
              J0=IRHS*(IRHS-1)/2
              DO ILHS=1,IRHS
                NL1=IP5(ILHS,1)
                NL2=IP5(ILHS,2)
                NL3=IP5(ILHS,3)
                NL4=IP5(ILHS,4)
                NL5=IP5(ILHS,5)
                XA5(ILHS+J0)=XA5(ILHS+J0)+
     1          TEMP3(NL3,NL4,NL5,NR3,NR4,NR5)*X0(NL1,NR1,M1)*
     2          Y0(NL2,NR2,M2)
              END DO
            END DO
          END DO
C**END 4-MODE INTEGRATION
        END DO
C**END 5-MODE INTEGRATION
      END IF
      CALL MATOUT(XA5,XRA5,NSIZE,25)
      GO TO 7777
6666  CALL MATIN(XA5,XRA5,NSIZE,25)
7777  CONTINUE

      IF(IND.NE.0)THEN
        DO IR=1,NSIZE
          NR=0
          DO J=1,NC1
            J1=JC1(J)
            NR1=IP5(IR,J1)
            NR=NR+1
            NX5(NR)=NR1
          END DO
          DO J=1,NC2
            J2=JC2(J)
            NR2=IP5(IR,J2)
            NR=NR+1
            NX5(NR)=NR2
          END DO
          IRL(NX5(1),NX5(2),NX5(3),NX5(4),NX5(5))=IR
        END DO
      END IF

C***********************************ALGORITHM FROM VCI5

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
C**TEMPORARY-LAN
      DO I=1,ISYMST-1
        IOFFR=IOFFR+NCSIZE(I)
      END DO
      J0R=IOFFR*(IOFFR+1)/2
C**TEMPORARY-LAN
      JSR=1
C**TEMPORARY-LAN
C     DO 9999 ISMR1=1,NVSYM
      DO 9999 ISMR1=ISYMST,ISYMFN
C**TEMPORARY-LAN
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE

        IF(NONC1.GT.NONC2)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          KCCCL=KCSZL
CC        IF(NONC1.GT.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             KCCCL=KCSZL
C           ELSE
C             KCCCL=IRH1
C           END IF
CC        END IF
          DO ILH1=1,KCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=IROFF
            NONZC1(2,NONZ1)=ILOFF
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
        IF(NONZ1.EQ.0)GO TO 5556
        NONZC1(3,NONZ1+1)=0
        END IF

      END IF

CCCC  IF(NONZ1.EQ.0)GO TO 5556
CCCC  NONZC1(3,NONZ1+1)=0

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE

        IF(NONC1.LE.NONC2)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          LCCCL=LCSZL
CC        IF(NONC1.LE.NONC2)THEN
C           IF(ISML1.NE.ISMR1)THEN
C             LCCCL=LCSZL
C           ELSE
C             LCCCL=IRH2
C           END IF
CC        END IF
          DO ILH2=1,LCCCL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=JROFF
            NONZC2(2,NONZ2)=JLOFF
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
        IF(NONZ2.EQ.0)GO TO 5555
        NONZC2(3,NONZ2+1)=0
        END IF

      END IF

CCCC  IF(NONZ2.EQ.0)GO TO 5555
CCCC  NONZC2(3,NONZ2+1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM5A.EQ.0)THEN
C       WRITE(IOUT,*)'VCCI5B - ISMR1,ISML1 = ',ISMR1,ISML1
C       WRITE(IOUT,*)'ICONDP = ',ICONDP
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C     END IF

      IF(NONC1.LE.NONC2)GO TO 3003
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M AND/OR J)

      DO IRH2=1,LCSZR
        JROFF=IRH2+LCOFFR
        DO I=1,NKVALR
          DO J=1,NKVALL
            DO K=1,NLVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          LCCCL=LCSZL
        ELSE
          LCCCL=IRH2
        END IF
        DO ILH2=1,LCCCL
          JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
          DO K=1,NKMOD2
            IGOT=0
            DO I=1,NC2
              IF(K.EQ.MC2(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 6001
          END DO
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH2.EQ.IRH2)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NKVALR
            DO J=1,NKVALL
              XCON2(J,I)=0
            END DO
          END DO
C
C**CONTRACTION SCHEME  '1' (K)
C
C         DO IRH1=1,KCSZR
C           IROFF=IRH1+KCOFFR
C           DO J=1,NKVALL
C             XK1(J,1)=0
C           END DO
C           DO ILH1=1,KCSZL
C             ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
C             DO K=1,NKMOD1
C               IGOT=0
C               DO I=1,NC1
C                 IF(K.EQ.MC1(I))IGOT=1
C               END DO
C               IF(IGOT.EQ.0.AND.
C    1          (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 5001
C             END DO

          IPREV=0
          DO INH1=1,NONZ1
              IRH1=NONZC1(3,INH1)
              INEXT=NONZC1(3,INH1+1)
              ILH1=NONZC1(4,INH1)
              IROFF=NONZC1(1,INH1)
              ILOFF=NONZC1(2,INH1)
              IF(IPREV.NE.IRH1)THEN
C**NEW RHS
                DO J=1,NKVALL
                  XK1(J,1)=0
                END DO
              END IF

              DO I=1,NC1
                NCR1(I)=IPC1(IROFF,MC1(I))
                NCL1(I)=IPC1(ILOFF,MC1(I))
              END DO

C**FIND RHS INDEX
C             DO IR=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCR1(I).NE.IP5(IR,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCR2(I).NE.IP5(IR,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 1000
C             END DO
1000          CONTINUE
C**FIND LHS INDEX
C             DO IL=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCL1(I).NE.IP5(IL,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCL2(I).NE.IP5(IL,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 2000
C             END DO
2000          CONTINUE

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX5(NR)=NR1
        MX5(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX5(NR)=NR2
        MX5(NR)=NL2
      END DO
      IR=IRL(NX5(1),NX5(2),NX5(3),NX5(4),NX5(5))
      IL=IRL(MX5(1),MX5(2),MX5(3),MX5(4),MX5(5))
C**TEMPORARY

C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
C             XYZ=XA5(I)
C             ZYX=0
C             IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
C               NR1=IP5(MR,1)
C               NL1=IP5(ML,1)
C               NR2=IP5(MR,2)
C               NL2=IP5(ML,2)
C               NR3=IP5(MR,3)
C               NL3=IP5(ML,3)
C               NR4=IP5(MR,4)
C               NL4=IP5(ML,4)
C               NR5=IP5(MR,5)
C               NL5=IP5(ML,5)
C               DO I=1,NP5(IND5)
C                 K=JP5(I,IND5,1)+1
C                 L=JP5(I,IND5,2)+1
C                 N=JP5(I,IND5,3)+1
C                 M=JP5(I,IND5,4)+1
C                 J=JP5(I,IND5,5)+1
C                 ZYX=ZYX+CP5(I,IND5)*XKAN(NL1,NR1,K,MOD1)*
C    1            XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
C    2            XKAN(NL4,NR4,M,MOD4)*XKAN(NL5,NR5,J,MOD5)
C               END DO
C**ANALYTIC
C             END IF
C             XXX=XYZ+ZYX
              XXX=XA5(I)

              DO K=1,NKVALL
                XK1(K,1)=XK1(K,1)+XXX*CFS1(ILH1,K,KEL,KSML)
              END DO
5001          CONTINUE

C           END DO
            IF(IRH1.NE.INEXT)THEN

            DO K=1,NKVALR
              DO L=1,NKVALL
                XCON2(L,K)=XCON2(L,K)+XK1(L,1)*CFS1(IRH1,K,KEL,KSMR)
              END DO
            END DO

            END IF
            IPREV=IRH1
          END DO
CC
          DO I=1,NKVALR
            DO J=1,NKVALL
              DO K=1,NLVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS2(ILH2,K,KEL,LSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6001      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN

C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              DO ILR=1,NLVALR
                CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILL,IKR,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKR,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              CCCC=CFS2(IRH2,ILR,KEL,LSMR)
              DO IKR=1,NKVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(ILL,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS2(IRH2,ILL,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(ILR,IKL,IKR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS2(IRH2,ILR,KEL,LSMR)*XK2(ILR,IKL,IKR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO
C*************************************

      GO TO 5555

3003  CONTINUE

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CONTRACTION SCHEME '1' (K)

      DO IRH1=1,KCSZR
        IROFF=IRH1+KCOFFR
        DO I=1,NLVALR
          DO J=1,NLVALL
            DO K=1,NKVALL
              XK2(K,J,I,1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              XK2(K,J,I,2)=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            END DO
          END DO
        END DO
        IF(ISML1.NE.ISMR1)THEN
          KCCCL=KCSZL
        ELSE
          KCCCL=IRH1
        END IF
        DO ILH1=1,KCCCL
          ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
          DO K=1,NKMOD1
            IGOT=0
            DO I=1,NC1
              IF(K.EQ.MC1(I))IGOT=1
            END DO
            IF(IGOT.EQ.0.AND.
     1      (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))GO TO 6002
          END DO
          DO I=1,NC1
            NCR1(I)=IPC1(IROFF,MC1(I))
            NCL1(I)=IPC1(ILOFF,MC1(I))
          END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IS3=1
          IF(ILH1.EQ.IRH1)IS3=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO I=1,NLVALR
            DO J=1,NLVALL
              XCON2(J,I)=0
            END DO
          END DO

C**CONTRACTION SCHEME  '2' (L AND/OR N AND/OR M AND/OR J)

C         DO IRH2=1,LCSZR
C           JROFF=IRH2+LCOFFR
C           DO J=1,NLVALL
C             XK1(J,1)=0
C           END DO
C           DO ILH2=1,LCSZL
C             JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
C             DO K=1,NKMOD2
C               IGOT=0
C               DO I=1,NC2
C                 IF(K.EQ.MC2(I))IGOT=1
C               END DO
C               IF(IGOT.EQ.0.AND.
C    1          (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))GO TO 5002
C             END DO

          JPREV=0
          DO INH2=1,NONZ2
              IRH2=NONZC2(3,INH2)
              JNEXT=NONZC2(3,INH2+1)
              ILH2=NONZC2(4,INH2)
              JROFF=NONZC2(1,INH2)
              JLOFF=NONZC2(2,INH2)
              IF(JPREV.NE.IRH2)THEN
C**NEW RHS
                DO J=1,NLVALL
                  XK1(J,1)=0
                END DO
              END IF

              DO I=1,NC2
                NCR2(I)=IPC2(JROFF,MC2(I))
                NCL2(I)=IPC2(JLOFF,MC2(I))
              END DO

C**FIND RHS INDEX
C             DO IR=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCR1(I).NE.IP5(IR,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCR2(I).NE.IP5(IR,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 1001
C             END DO
1001          CONTINUE
C**FIND LHS INDEX
C             DO IL=1,NSIZE
C               IGOT=1
C               DO I=1,NC1
C                 IF(NCL1(I).NE.IP5(IL,JC1(I)))IGOT=0
C               END DO
C               DO I=1,NC2
C                 IF(NCL2(I).NE.IP5(IL,JC2(I)))IGOT=0
C               END DO
C               IF(IGOT.EQ.1)GO TO 2001
C             END DO
2001          CONTINUE

C**TEMPORARY
      NR=0
      DO J=1,NC1
        NR1=NCR1(J)
        NL1=NCL1(J)
        NR=NR+1
        NX5(NR)=NR1
        MX5(NR)=NL1
      END DO
      DO J=1,NC2
        NR2=NCR2(J)
        NL2=NCL2(J)
        NR=NR+1
        NX5(NR)=NR2
        MX5(NR)=NL2
      END DO
      IR=IRL(NX5(1),NX5(2),NX5(3),NX5(4),NX5(5))
      IL=IRL(MX5(1),MX5(2),MX5(3),MX5(4),MX5(5))
C**TEMPORARY

C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
C             XYZ=XA5(I)
C             ZYX=0
C             IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
C               NR1=IP5(MR,1)
C               NL1=IP5(ML,1)
C               NR2=IP5(MR,2)
C               NL2=IP5(ML,2)
C               NR3=IP5(MR,3)
C               NL3=IP5(ML,3)
C               NR4=IP5(MR,4)
C               NL4=IP5(ML,4)
C               NR5=IP5(MR,5)
C               NL5=IP5(ML,5)
C               DO I=1,NP5(IND5)
C                 K=JP5(I,IND5,1)+1
C                 L=JP5(I,IND5,2)+1
C                 N=JP5(I,IND5,3)+1
C                 M=JP5(I,IND5,4)+1
C                 J=JP5(I,IND5,5)+1
C                 ZYX=ZYX+CP5(I,IND5)*XKAN(NL1,NR1,K,MOD1)*
C    1            XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
C    2            XKAN(NL4,NR4,M,MOD4)*XKAN(NL5,NR5,J,MOD5)
C               END DO
C**ANALYTIC
C             END IF
C             XK1(ILH2,IRH2)=XYZ+ZYX
C             XXX=XYZ+ZYX
              XXX=XA5(I)

              DO J=1,NLVALL
                XK1(J,1)=XK1(J,1)+XXX*CFS2(ILH2,J,KEL,LSML)
              END DO
5002          CONTINUE

C           END DO
            IF(IRH2.NE.JNEXT)THEN

            DO I=1,NLVALR
              DO J=1,NLVALL
                XCON2(J,I)=XCON2(J,I)+XK1(J,1)*CFS2(IRH2,I,KEL,LSMR)
              END DO
            END DO

            END IF
            JPREV=IRH2
          END DO
CC
          DO I=1,NLVALR
            DO J=1,NLVALL
              DO K=1,NKVALL
                XK2(K,J,I,1)=XK2(K,J,I,1)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(J,I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                XK2(K,J,I,2)=XK2(K,J,I,2)+CFS1(ILH1,K,KEL,KSML)*
     1          XCON2(I,J)*IS3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

              END DO
            END DO
          END DO
6002      CONTINUE
        END DO

        IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,NKVALL
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,NLVALL
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
C**OFF-DIAGONAL BLOCK - END
        ELSE
C**DIAGONAL BLOCK - START
          IF(KCONT1.EQ.1)THEN
            IRHS=0
            DO IKR=1,NKVALR
              CCCC=CFS1(IRH1,IKR,KEL,KSMR)
              DO ILR=1,NLVALR
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO IKL=1,IKR-1
                  DO ILL=1,NLVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO ILL=1,ILR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKR,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKR,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          ELSE
            IRHS=0
            DO ILR=1,NLVALR
              DO IKR=1,NKVALR
                CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                IRHS=IRHS+1
                IROFF=IRHS+IOFFR
C**TEMPORARY-LAN
                J0=IROFF*(IROFF-1)/2+IOFFL-J0R
C**TEMPORARY-LAN
                ILHS=0
                DO ILL=1,ILR-1
                  DO IKL=1,NKVALL
                    ILHS=ILHS+1
                    XA(ILHS+J0)=XA(ILHS+J0)+
     1              CCCC*XK2(IKL,ILL,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2             +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILL,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  END DO
                END DO
                DO IKL=1,IKR
                  ILHS=ILHS+1
                  XA(ILHS+J0)=XA(ILHS+J0)+
     1            CCCC*XK2(IKL,ILR,ILR,1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     2           +CFS1(IRH1,IKL,KEL,KSMR)*XK2(IKR,ILR,ILR,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                END DO
              END DO
            END DO
          END IF
C**DIAGONAL BLOCK - END

        END IF
      END DO

5555  CONTINUE
5556  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE
      IF(ITIM5A.EQ.0)THEN
        IF(IND.EQ.1)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM5A=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI5C(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,
     1MODE2,MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,
     2NN1,MM1,NN2,MM2,NN3,MM3,NN4,MM4,NN5,MM5,IP,ISIZMX,IPC,ICSIZE,
     3IPSIZE,KCONT,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,NSIZE,XK,TEMP,IRL,XCON,
     4NVAL,ISTART,IEND,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,JCI1,JCI2,JCI3,
     5JCI4,JCI5,X0,Y0,Z0,W0,U0,VP,VPR,MODINT,KEL,NS,CFS,ISIZXX,NVALX,
     6KEL21,NVSYMX,
C**ANALYTIC
     7XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM)
C**ANALYTIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM5,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H4(NN4,MM4,3),H5(NN5,MM5,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQ5(MM5)
      DIMENSION IP(ISIZMX,NNMODE),IRL(ICSIZE)
      DIMENSION IPC(IPSIZE,1),IP5(ISIZE5,5),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(JCI5,JCI5),TEMP2(JCI4,JCI5,JCI4,JCI5)
      DIMENSION TEMP3(JCI3,JCI4,JCI5,JCI3,JCI4,JCI5)
C     DIMENSION TEMP4(JCI2*JCI3*JCI4*JCI5,JCI2*JCI3*JCI4*JCI5)
C     DIMENSION TEMP5(JCI2,JCI3,JCI4,JCI5)
      DIMENSION X0(JCI1,JCI1,MM1),Y0(JCI2,JCI2,MM2)
      DIMENSION Z0(JCI3,JCI3,MM3),W0(JCI4,JCI4,MM4)
      DIMENSION U0(JCI5,JCI5,MM5)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA5(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP5(NTOT5),CP5(MAX5,NTOT5),JP5(MAX5,NTOT5,5)
      DIMENSION INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
C---TEMP ARRAY INITIALIZE---ACTUALLY THIS CAN BE DONE SOMEWHERE ELSE
C---BEFORE THIS SUB GET CALLED, THEN TEMP4 AND TEMP CAN BE CALLED.
C     I1=0
C     DO I2=1,JCI2
C      DO I3=1,JCI3
C       DO I4=1,JCI4
C        DO I5=1,JCI5
C         I1=I1+1
C         TEMP5(I2,I3,I4,I5)=I1
C        END DO
C       END DO
C      END DO
C     END DO
C---END OF TEMP ARRAY INITIALIZATION
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      MD5=MODINT(MOD5)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(MOD5.EQ.ISYM(I,J))N5=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N1.EQ.N5)MD5=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N2.EQ.N5)MD5=1
      IF(N3.EQ.N4)MD4=1
      IF(N3.EQ.N5)MD5=1
      IF(N4.EQ.N5)MD5=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      IF(N12.EQ.N5)MD5=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      IF(N13.EQ.N5)MD5=1
      N14=ISYMP(N1,N4)
      IF(N14.EQ.N5)MD5=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      IF(N23.EQ.N5)MD5=1
      N24=ISYMP(N2,N4)
      IF(N24.EQ.N5)MD5=1
      N34=ISYMP(N3,N4)
      IF(N34.EQ.N5)MD5=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      IF(N123.EQ.N5)MD5=1
      N124=ISYMP(N12,N4)
      IF(N124.EQ.N5)MD5=1
      N134=ISYMP(N13,N4)
      IF(N134.EQ.N5)MD5=1
      N234=ISYMP(N23,N4)
      IF(N234.EQ.N5)MD5=1
      N1234=ISYMP(N12,N34)
      IF(N1234.EQ.N5)MD5=1
      CALL VDCI5C(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,MODE2,
     1MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,NN1,MM1,
     2MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,NN5,MM5,
     3MM5/MD5,IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP5,ISIZE5,
     4XA5,XRA5,NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,
     5TEMP3,TEMP4,TEMP5,JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,W0,U0,VP,VPR,
     6MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,
C**ANALYTIC
     7XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM)
C**ANALYTIC
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI5C(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,
     1MODE2,MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,
     2NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,NN4,MH4,MM4,NN5,MH5,MM5,
     3IP,ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,
     4NSIZE,XK,TEMP,IRL,XCON,NVAL,ISTART,IEND,TEMP1,TEMP2,TEMP3,TEMP4,
     5TEMP5,JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,W0,U0,VP,VPR,
     6MODINT,KEL,NS,CFS,ISIZXX,NVALX,KEL21,NVSYMX,
C**ANALYTIC
     7XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM)
C**ANALYTIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM5,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION H4(NN4,MH4,3),H5(NN5,MH5,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQ5(MM5)
      DIMENSION IP(ISIZMX,NNMODE),IRL(ICSIZE)
      DIMENSION IPC(IPSIZE,1),IP5(ISIZE5,5),TEMP(ICSIZE,NVAL)
      DIMENSION TEMP1(JCI5,JCI5),TEMP2(JCI4,JCI5,JCI4,JCI5)
      DIMENSION TEMP3(JCI3,JCI4,JCI5,JCI3,JCI4,JCI5)
C     DIMENSION TEMP4(JCI2*JCI3*JCI4*JCI5,JCI2*JCI3*JCI4*JCI5)
C     DIMENSION TEMP5(JCI2,JCI3,JCI4,JCI5)
      DIMENSION X0(JCI1,JCI1,MM1),Y0(JCI2,JCI2,MM2)
      DIMENSION Z0(JCI3,JCI3,MM3),W0(JCI4,JCI4,MM4)
      DIMENSION U0(JCI5,JCI5,MM5)
      DIMENSION XA(1),XK(ICSIZE,ICSIZE),XA5(1),XCON(NVAL,NVAL)
      DIMENSION CFS(ISIZXX,NVALX,KEL21,NVSYMX)
C**ANALYTIC
      DIMENSION NP5(NTOT5),CP5(MAX5,NTOT5),JP5(MAX5,NTOT5,5)
      DIMENSION INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM5A.EQ.1)THEN
        WRITE(IOUT,*)'Calculating VCCI5C'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

C** FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1', 'MOD2',
C** 'MOD3', MOD4' AND 'MOD5'
C**IF THEY ARE IN SCHEME 1.....
      NKMODE=ICONT(1)
      NLMODE=ICONT(2)
      MCONT=2
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        NLMODE=ICONT(1)
        MCONT=1
      END IF

C**ANALYTIC
      IND=INDM(MOD1)+INDN(MOD2)+INDL(MOD3)+INDK(MOD4)+MOD5
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI5
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=INTFAC(NMODE,ICOUPL,5)
      IF(MOLINC.GT.0)IFACTL=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(75)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(75)VPR
      END IF

C***********************************************************

      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      MD5=MODINT(MOD5)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(MOD5.EQ.ISYM(I,J))N5=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N1.EQ.N5)MD5=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N2.EQ.N5)MD5=1
      IF(N3.EQ.N4)MD4=1
      IF(N3.EQ.N5)MD5=1
      IF(N4.EQ.N5)MD5=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      IF(N12.EQ.N5)MD5=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      IF(N13.EQ.N5)MD5=1
      N14=ISYMP(N1,N4)
      IF(N14.EQ.N5)MD5=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      IF(N23.EQ.N5)MD5=1
      N24=ISYMP(N2,N4)
      IF(N24.EQ.N5)MD5=1
      N34=ISYMP(N3,N4)
      IF(N34.EQ.N5)MD5=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      IF(N123.EQ.N5)MD5=1
      N124=ISYMP(N12,N4)
      IF(N124.EQ.N5)MD5=1
      N134=ISYMP(N13,N4)
      IF(N134.EQ.N5)MD5=1
      N234=ISYMP(N23,N4)
      IF(N234.EQ.N5)MD5=1
      N1234=ISYMP(N12,N34)
      IF(N1234.EQ.N5)MD5=1
      MD=MD1*MD2*MD3*MD4*MD5

C***********************************************************
      IF(IWHICH.GE.0)THEN
C**POTENTIAL
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
        DO M1=1,MM1
          DO NR1=1,JCI1
            X=H1(NR1,M1,1)*IFACTL*MD
            DO NL1=1,JCI1
              Y=H1(NL1,M1,1)
              X0(NL1,NR1,M1)=Y*X
            END DO
          END DO
        END DO
        DO M2=1,MM2
          DO NR2=1,JCI2
            X=H2(NR2,M2,1)
            DO NL2=1,JCI2
              Y=H2(NL2,M2,1)
              Y0(NL2,NR2,M2)=Y*X
            END DO
          END DO
        END DO
        DO M3=1,MM3
          DO NR3=1,JCI3
            X=H3(NR3,M3,1)
            DO NL3=1,JCI3
              Y=H3(NL3,M3,1)
              Z0(NL3,NR3,M3)=Y*X
            END DO
          END DO
        END DO
        DO M4=1,MM4
          DO NR4=1,JCI4
            X=H4(NR4,M4,1)
            DO NL4=1,JCI4
              Y=H4(NL4,M4,1)
              W0(NL4,NR4,M4)=Y*X
            END DO
          END DO
        END DO
        DO M5=1,MM5
          DO NR5=1,JCI5
            X=H5(NR5,M5,1)
            DO NL5=1,JCI5
              Y=H5(NL5,M5,1)
              U0(NL5,NR5,M5)=Y*X
            END DO
          END DO
        END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**START 5-MODE INTEGRATION
        DO M1=1,MM1
C**START 4-MODE INTEGRATION
          DO M2=1,MM2
            DO IRHS5=1,JCI5
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS5=1,JCI5
                    DO ILHS4=1,JCI4
                      DO ILHS3=1,JCI3
                        TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)=0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 3-MODE INTEGRATION
            DO M3=1,MM3
              DO IRHS5=1,JCI5
                DO IRHS4=1,JCI4
                  DO ILHS5=1,JCI5
                    DO ILHS4=1,JCI4
                      TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)=0
                    END DO
                  END DO
                END DO
              END DO
C**START 2-MODE INTEGRATION
              DO M4=1,MM4
                DO IRHS5=1,JCI5
                  DO ILHS5=1,JCI5
                    TEMP1(ILHS5,IRHS5)=0
                  END DO
                END DO
C**START 1-MODE INTEGRATION
                DO M5=1,MM5
                  IF(JCOUPL.GT.0)THEN
                    IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1              C=VP(M5,M4,M3,M2,M1)
                  ELSE
                    IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1              C=VPR(M5,M4,M3,M2,M1)
                  END IF
                  DO IRHS5=1,JCI5
                    DO ILHS5=1,JCI5
                      TEMP1(ILHS5,IRHS5)=TEMP1(ILHS5,IRHS5)+
     1                U0(ILHS5,IRHS5,M5)*C
                    END DO
                  END DO
                END DO
C**END 1-MODE INTEGRATION
                DO IRHS5=1,JCI5
                  DO IRHS4=1,JCI4
                    DO ILHS5=1,JCI5
                      DO ILHS4=1,JCI4
                        TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)=
     1                  TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)+
     2                  W0(ILHS4,IRHS4,M4)*TEMP1(ILHS5,IRHS5)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
C**END 2-MODE INTEGRATION
              DO IRHS5=1,JCI5
                DO IRHS4=1,JCI4
                  DO IRHS3=1,JCI3
                    DO ILHS5=1,JCI5
                      DO ILHS4=1,JCI4
                        DO ILHS3=1,JCI3
                          TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)=
     1                    TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)+
     2                    Z0(ILHS3,IRHS3,M3)*
     3                    TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 3-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (5-DIM)
            DO IRHS=1,NSIZE
              NR1=IP5(IRHS,1)
              NR2=IP5(IRHS,2)
              NR3=IP5(IRHS,3)
              NR4=IP5(IRHS,4)
              NR5=IP5(IRHS,5)
              J0=IRHS*(IRHS-1)/2
              DO ILHS=1,IRHS
                NL1=IP5(ILHS,1)
                NL2=IP5(ILHS,2)
                NL3=IP5(ILHS,3)
                NL4=IP5(ILHS,4)
                NL5=IP5(ILHS,5)
                XA5(ILHS+J0)=XA5(ILHS+J0)+
     1          TEMP3(NL3,NL4,NL5,NR3,NR4,NR5)*X0(NL1,NR1,M1)*
     2          Y0(NL2,NR2,M2)
              END DO
            END DO
          END DO
C**END 4-MODE INTEGRATION
        END DO
C**END 5-MODE INTEGRATION
      END IF
      CALL MATOUT(XA5,XRA5,NSIZE,25)
      GO TO 7777
6666  CALL MATIN(XA5,XRA5,NSIZE,25)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI5

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**OFFSET FOR FINAL MATRIX (XA)
      IOFF=0
      JS=1
      DO 9999 ISM1=1,NVSYM
      KSR=ISM1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISM1)
      IF(NS.EQ.1)THEN
        ISM2=ISM1
      END IF
      IF(NS.EQ.2)THEN
        ISM2=ISM1+JS
        JS=-JS
      END IF
      IF(NS.EQ.3)THEN
        ISM2=ISM1+2*JS
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISM2=NVSYM+1-ISM1
        ISM2=5+NVSYM*KSR-ISM1
      END IF
      IF(NS.EQ.5)THEN
        ISM2=ISM1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISM2=ISM1+JS+4*LSR
        JS=-JS
      END IF
      IF(NS.EQ.7)THEN
        ISM2=ISM1+2*JS+4*LSR
        JJJ=MOD(ISM1,2)+1
        JS=JS*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISM2=5+NVSYM*KSR-ISM1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISM2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISM1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISM2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISM=ISM1
        KCOFF=ICOFF1
        LCOFF=ICOFF2
        KCSZ=ICSZ1
        LCSZ=ICSZ2
        ICNTZ=ICCONT(1)
        NKVAL=NCVAL(1,ISM1)
        NLVAL=NCVAL(2,ISM2)
      ELSE
        ISM=ISM2
        KCOFF=ICOFF2
        LCOFF=ICOFF1
        KCSZ=ICSZ2
        LCSZ=ICSZ1
        ICNTZ=ICCONT(2)
        NKVAL=NCVAL(2,ISM2)
        NLVAL=NCVAL(1,ISM1)
      END IF
C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

C
      DO IKRL=1,KCSZ
        IOFFRL=IKRL+KCOFF
        NR1=IPC(IOFFRL,MODE1)
        NR2=IPC(IOFFRL,MODE2)
        NR3=IPC(IOFFRL,MODE3)
        NR4=IPC(IOFFRL,MODE4)
        NR5=IPC(IOFFRL,MODE5)
        DO IR=1,NSIZE
          IF(NR1.EQ.IP5(IR,1).AND.NR2.EQ.IP5(IR,2).AND.
     1    NR3.EQ.IP5(IR,3).AND.NR4.EQ.IP5(IR,4).AND.
     2    NR5.EQ.IP5(IR,5))GO TO 1000
        END DO
1000    IRL(IKRL)=IR
      END DO
C
C**CASE 1..INTEGRATION MODES IN CONTRACTED SCHEME.
C**      ..OVERLAP OF NON-CONTRACTED BASIS.
      IF(ICNTZ.GT.0)THEN
        DO IRHS=1,KCSZ
          IOFFR=IRHS+KCOFF
          IR=IRL(IRHS)
          DO ILHS=1,IRHS
            IOFFL=ILHS+KCOFF
C**OVERLAP OF REMAINING STATES
            IS=1
            DO K=1,NKMODE
              IF(IS.EQ.0)GO TO 2000
              IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1        MODE4.AND.K.NE.MODE5.AND.(IPC(IOFFR,K).NE.
     2        IPC(IOFFL,K)))IS=0
            END DO
C**OVERLAP OF REMAINING STATES
            IL=IRL(ILHS)
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA5(I)
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              XYZ=0
              DO I=1,NP5(IND)
                K=JP5(I,IND,1)+1
                L=JP5(I,IND,2)+1
                N=JP5(I,IND,3)+1
                M=JP5(I,IND,4)+1
                J=JP5(I,IND,5)+1
                XYZ=XYZ+CP5(I,IND)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)*XKAN(NL5,NR5,J,MOD5)
              END DO
C**ANALYTIC
            END IF
2000        CONTINUE
            XK(ILHS,IRHS)=(XYZ+ZYX)*IS
            XK(IRHS,ILHS)=XK(ILHS,IRHS)
          END DO
        END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',KCSZ,NKVAL,KCSZ,1.0D0,XK(1,1),ICSIZE,
     &           CFS(1,1,KEL,ISM),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
        CALL DGEMM('T','N',NKVAL,NKVAL,KCSZ,1.0D0,CFS(1,1,KEL,ISM),
     &           ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

        IRHS=0
        DO IKR=1,NKVAL
          DO ILR=1,LCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            DO IKL=1,IKR
              ILHS=(IKL-1)*LCSZ+ILR
              XYZ=XCON(IKL,IKR)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ
            END DO
          END DO
        END DO

C**CASE 2..INTEGRATION MODES IN NON-CONTRACTED SCHEME.
C**      ..OVERLAP OF CONTRACTED BASIS.
      ELSE
        IRHS=0
        DO ILR=1,NLVAL
          DO IKR=1,KCSZ
            IRHS=IRHS+1
            IROFF=IRHS+IOFF
            IOFFR=IKR+KCOFF
            J0=IROFF*(IROFF-1)/2+IOFF
            IR=IRL(IKR)
            DO IKL=1,IKR
              ILHS=(ILR-1)*KCSZ+IKL
              IOFFL=IKL+KCOFF
C**OVERLAP OF REMAINING STATES
              IS=1
              DO K=1,NKMODE
                IF(IS.EQ.0)GO TO 3000
                IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.
     1          K.NE.MODE4.AND.K.NE.MODE5.AND.(IPC(IOFFR,K).NE.
     2          IPC(IOFFL,K)))IS=0
              END DO
C**OVERLAP OF REMAINING STATES
              IL=IRL(IKL)
C**GET MATRIX ELEMENT
              MR=IR
              ML=IL
              IF(IR.LT.IL)THEN
                MR=IL
                ML=IR
              END IF
              I=MR*(MR-1)/2+ML
              XYZ=XA5(I)
              XA(J0+ILHS)=XA(J0+ILHS)+XYZ*IS
3000          CONTINUE
            END DO
          END DO
        END DO
      END IF

C**UPDATE POINTER FOR XA
      IOFF=IOFF+NCSIZE(ISM1)
9998  CONTINUE
9999  CONTINUE
      IF(ITIM5A.EQ.1)THEN
        CALL TIMIT(3)
        ITIM5A=2
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCCI5D(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,
     1MODE2,MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,
     2NN1,MM1,NN2,MM2,NN3,MM3,NN4,MM4,NN5,MM5,IP,ISIZMX,IPC1,ICSIZ1,
     3IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,KCONT2,XA,ISIZE,IP5,ISIZE5,XA5,
     4XRA5,NSIZE,XK1,XK2,TEMP,KTEMP,XCON2,NVAL1,NVAL2,ISTART,IEND,
     5TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,
     6W0,U0,VP,VPR,MODINT,KEL,NS,CFS1,CFS2,
     7ISIZXX,NVALX,KEL21,NVSYMX,IND,NONZC1,NONZC2,
C**ANALYTIC
     8XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM,
C**ANALYTIC
     9MAXVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM5,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION H4(NN4,MM4,3),H5(NN5,MM5,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP5(ISIZE5,5),TEMP(KTEMP,1)
      DIMENSION TEMP1(JCI5,JCI5),TEMP2(JCI4,JCI5,JCI4,JCI5)
      DIMENSION TEMP3(JCI3,JCI4,JCI5,JCI3,JCI4,JCI5)
C     DIMENSION TEMP4(JCI2*JCI3*JCI4*JCI5,JCI2*JCI3*JCI4*JCI5)
C     DIMENSION TEMP5(JCI2,JCI3,JCI4,JCI5)
      DIMENSION X0(JCI1,JCI1,MM1),Y0(JCI2,JCI2,MM2)
      DIMENSION Z0(JCI3,JCI3,MM3),W0(JCI4,JCI4,MM4)
      DIMENSION U0(JCI5,JCI5,MM5)
C****************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA5(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,1)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C****************************************************************
C**ANALYTIC
      DIMENSION NP5(NTOT5),CP5(MAX5,NTOT5),JP5(MAX5,NTOT5,5)
      DIMENSION INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
C---TEMP ARRAY INITIALIZE---ACTUALLY THIS CAN BE DONE SOMEWHERE ELSE
C---BEFORE THIS SUB GET CALLED, THEN TEMP4 AND TEMP CAN BE CALLED.
C     I1=0
C     DO I2=1,JCI2
C      DO I3=1,JCI3
C       DO I4=1,JCI4
C        DO I5=1,JCI5
C         I1=I1+1
C         TEMP5(I2,I3,I4,I5)=I1
C        END DO
C       END DO
C      END DO
C     END DO
C---END OF TEMP ARRAY INITIALIZATION
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      MD5=MODINT(MOD5)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(MOD5.EQ.ISYM(I,J))N5=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N1.EQ.N5)MD5=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N2.EQ.N5)MD5=1
      IF(N3.EQ.N4)MD4=1
      IF(N3.EQ.N5)MD5=1
      IF(N4.EQ.N5)MD5=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      IF(N12.EQ.N5)MD5=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      IF(N13.EQ.N5)MD5=1
      N14=ISYMP(N1,N4)
      IF(N14.EQ.N5)MD5=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      IF(N23.EQ.N5)MD5=1
      N24=ISYMP(N2,N4)
      IF(N24.EQ.N5)MD5=1
      N34=ISYMP(N3,N4)
      IF(N34.EQ.N5)MD5=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      IF(N123.EQ.N5)MD5=1
      N124=ISYMP(N12,N4)
      IF(N124.EQ.N5)MD5=1
      N134=ISYMP(N13,N4)
      IF(N134.EQ.N5)MD5=1
      N234=ISYMP(N23,N4)
      IF(N234.EQ.N5)MD5=1
      N1234=ISYMP(N12,N34)
      IF(N1234.EQ.N5)MD5=1
      CALL VDCI5D(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,MODE2,
     1MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,NN1,MM1,
     2MM1/MD1,NN2,MM2,MM2/MD2,NN3,MM3,MM3/MD3,NN4,MM4,MM4/MD4,NN5,MM5,
     3MM5/MD5,IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,
     4KCONT2,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,NSIZE,XK1,XK2,TEMP,KTEMP,
     5XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,
     6JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,W0,U0,VP,VPR,
     6MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,
     7NONZC1,NONZC2,
C**ANALYTIC
     8XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM,
C**ANALYTIC
     9MAXVAL)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCI5D(NMODE,NNMODE,MOD1,MOD2,MOD3,MOD4,MOD5,MODE1,
     1MODE2,MODE3,MODE4,MODE5,H1,XQ1,H2,XQ2,H3,XQ3,H4,XQ4,H5,XQ5,
     2NN1,MH1,MM1,NN2,MH2,MM2,NN3,MH3,MM3,NN4,MH4,MM4,NN5,MH5,MM5,
     3IP,ISIZMX,IPC1,ICSIZ1,IPSIZ1,IPC2,ICSIZ2,IPSIZ2,KCONT1,
     3KCONT2,XA,ISIZE,IP5,ISIZE5,XA5,XRA5,NSIZE,XK1,XK2,TEMP,KTEMP,
     4XCON2,NVAL1,NVAL2,ISTART,IEND,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,
     5JCI1,JCI2,JCI3,JCI4,JCI5,X0,Y0,Z0,W0,U0,VP,VPR,
     6MODINT,KEL,NS,CFS1,CFS2,ISIZXX,NVALX,KEL21,NVSYMX,IND,
     7NONZC1,NONZC2,
C**ANALYTIC
     8XKAN,MAXQU,MAXPOW,NP5,CP5,JP5,NTOT5,MAX5,INDK,INDL,INDN,INDM,
C**ANALYTIC
     9MAXVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL LANCZ,LANZA,LANZB
      REAL*8 VP(MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM5,MM4,MM3,MM2,MM1)
      REAL*4 XRA5(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H1(NN1,MH1,3),H2(NN2,MH2,3),H3(NN3,MH3,3)
      DIMENSION H4(NN4,MH4,3),H5(NN5,MH5,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),XQ5(MM5)
      DIMENSION IP(ISIZMX,NNMODE),IPC1(IPSIZ1,1),IPC2(IPSIZ2,1)
      DIMENSION IP5(ISIZE5,5),TEMP(KTEMP,1)
      DIMENSION TEMP1(JCI5,JCI5),TEMP2(JCI4,JCI5,JCI4,JCI5)
      DIMENSION TEMP3(JCI3,JCI4,JCI5,JCI3,JCI4,JCI5)
C     DIMENSION TEMP4(JCI2*JCI3*JCI4*JCI5,JCI2*JCI3*JCI4*JCI5)
C     DIMENSION TEMP5(JCI2,JCI3,JCI4,JCI5)
      DIMENSION X0(JCI1,JCI1,MM1),Y0(JCI2,JCI2,MM2)
      DIMENSION Z0(JCI3,JCI3,MM3),W0(JCI4,JCI4,MM4)
      DIMENSION U0(JCI5,JCI5,MM5)
C****************************************************************
      DIMENSION XA(1),XK1(KTEMP,KTEMP),XA5(1),XCON2(MAXVAL,MAXVAL)
      DIMENSION XK2(MAXVAL,1)
      DIMENSION CFS1(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFS2(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION NONZC1(4,1),NONZC2(4,1)
C****************************************************************
C**ANALYTIC
      DIMENSION NP5(NTOT5),CP5(MAX5,NTOT5),JP5(MAX5,NTOT5,5)
      DIMENSION INDK(1),INDL(1),INDN(1),INDM(1)
      DIMENSION XKAN(MAXQU,MAXQU,MAXPOW,1)
C**ANALYTIC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
C************************
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM5A.EQ.0)THEN
        IF(IND.EQ.0)THEN
          WRITE(IOUT,*)'Calculating VCCI5D'
          CALL FLUSH(IOUT)
          CALL TIMIT(1)
        END IF
      END IF

C**************************************************************
C**************************************************************
C**BOTH CONTRACTION SCHEMES HAVE AT LEAST ONE MODE
C**CHECK SCHEME '1'
      I1=KCONT1
      NC1=0
      DO NN=1,ICONT(I1)
        IF(MOD1.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE1
          JC1(NC1)=1
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD2.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE2
          JC1(NC1)=2
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD3.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE3
          JC1(NC1)=3
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD4.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE4
          JC1(NC1)=4
        END IF
      END DO
      DO NN=1,ICONT(I1)
        IF(MOD5.EQ.JCONT(I1,NN))THEN
          NC1=NC1+1
          MC1(NC1)=MODE5
          JC1(NC1)=5
        END IF
      END DO
C**CHECK SCHEME '2'
      I2=KCONT2
      NC2=0
      DO NN=1,ICONT(I2)
        IF(MOD1.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE1
          JC2(NC2)=1
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD2.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE2
          JC2(NC2)=2
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD3.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE3
          JC2(NC2)=3
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD4.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE4
          JC2(NC2)=4
        END IF
      END DO
      DO NN=1,ICONT(I2)
        IF(MOD5.EQ.JCONT(I2,NN))THEN
          NC2=NC2+1
          MC2(NC2)=MODE5
          JC2(NC2)=5
        END IF
      END DO
C**NKMOD1 IS TOTAL NUMBER OF MODES FOR SCHEME INVOLVING 'MOD1' (K)
C**NKMOD2 IS TOTAL NUMBER OF MODES FOR THE OTHER SCHEME (ONE OR ALL
C**OF 'MOD2' (L), 'MOD3' (N), 'MOD4' (M) AND 'MOD5' (J))
C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1' (K)
C**IF IT IS IN SCHEME 1.....
      NKMOD1=ICONT(1)
C**....IT IS NOT IN SCHEME 2
      IF(KCONT1.EQ.2)THEN
        NKMOD1=ICONT(2)
      END IF
C**SECOND DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD2' (L),
C**'MOD3' (N) OR 'MOD4' (M),'MOD5' (J)
C**IF THEY ARE IN SCHEME 1.....
      NKMOD2=ICONT(1)
C**....THEY ARE NOT IN SCHEME 2
      IF(KCONT2.EQ.2)THEN
        NKMOD2=ICONT(2)
      END IF
C**************************************************************
      IF(IND.NE.0)GO TO 7777
C**************************************************************

C**ANALYTIC
      IND5=INDM(MOD1)+INDN(MOD2)+INDL(MOD3)+INDK(MOD4)+MOD5
C**ANALYTIC
C***********************************ALGORITHM FROM V0CI5
      IF(ICONDP.NE.0)GO TO 6666

      IFACTL=INTFAC(NMODE,ICOUPL,5)
      IF(MOLINC.GT.0)IFACTL=1
      IF(IWHICH.EQ.0)IFATCL=1
      IF(IWHICH.LT.0.AND.MOLINC.LT.0)IFACTL=1

C***********************************************************

      IF(JCOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(75)VP
      ELSE
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)READ(75)VPR
      END IF

C***********************************************************

      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      MD5=MODINT(MOD5)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(MOD5.EQ.ISYM(I,J))N5=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N1.EQ.N5)MD5=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N2.EQ.N5)MD5=1
      IF(N3.EQ.N4)MD4=1
      IF(N3.EQ.N5)MD5=1
      IF(N4.EQ.N5)MD5=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      IF(N12.EQ.N5)MD5=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      IF(N13.EQ.N5)MD5=1
      N14=ISYMP(N1,N4)
      IF(N14.EQ.N5)MD5=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      IF(N23.EQ.N5)MD5=1
      N24=ISYMP(N2,N4)
      IF(N24.EQ.N5)MD5=1
      N34=ISYMP(N3,N4)
      IF(N34.EQ.N5)MD5=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      IF(N123.EQ.N5)MD5=1
      N124=ISYMP(N12,N4)
      IF(N124.EQ.N5)MD5=1
      N134=ISYMP(N13,N4)
      IF(N134.EQ.N5)MD5=1
      N234=ISYMP(N23,N4)
      IF(N234.EQ.N5)MD5=1
      N1234=ISYMP(N12,N34)
      IF(N1234.EQ.N5)MD5=1
      MD=MD1*MD2*MD3*MD4*MD5

C***********************************************************
      IF(IWHICH.GE.0)THEN
C**POTENTIAL
C**FORM INDIVIDUAL INTEGRATION TERMS (START)
        DO M1=1,MM1
          DO NR1=1,JCI1
            X=H1(NR1,M1,1)*IFACTL*MD
            DO NL1=1,JCI1
              Y=H1(NL1,M1,1)
              X0(NL1,NR1,M1)=Y*X
            END DO
          END DO
        END DO
        DO M2=1,MM2
          DO NR2=1,JCI2
            X=H2(NR2,M2,1)
            DO NL2=1,JCI2
              Y=H2(NL2,M2,1)
              Y0(NL2,NR2,M2)=Y*X
            END DO
          END DO
        END DO
        DO M3=1,MM3
          DO NR3=1,JCI3
            X=H3(NR3,M3,1)
            DO NL3=1,JCI3
              Y=H3(NL3,M3,1)
              Z0(NL3,NR3,M3)=Y*X
            END DO
          END DO
        END DO
        DO M4=1,MM4
          DO NR4=1,JCI4
            X=H4(NR4,M4,1)
            DO NL4=1,JCI4
              Y=H4(NL4,M4,1)
              W0(NL4,NR4,M4)=Y*X
            END DO
          END DO
        END DO
        DO M5=1,MM5
          DO NR5=1,JCI5
            X=H5(NR5,M5,1)
            DO NL5=1,JCI5
              Y=H5(NL5,M5,1)
              U0(NL5,NR5,M5)=Y*X
            END DO
          END DO
        END DO
C**FORM INDIVIDUAL INTEGRATION TERMS (END)
C**START 5-MODE INTEGRATION
        DO M1=1,MM1
C**START 4-MODE INTEGRATION
          DO M2=1,MM2
            DO IRHS5=1,JCI5
              DO IRHS4=1,JCI4
                DO IRHS3=1,JCI3
                  DO ILHS5=1,JCI5
                    DO ILHS4=1,JCI4
                      DO ILHS3=1,JCI3
                        TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)=0
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**START 3-MODE INTEGRATION
            DO M3=1,MM3
              DO IRHS5=1,JCI5
                DO IRHS4=1,JCI4
                  DO ILHS5=1,JCI5
                    DO ILHS4=1,JCI4
                      TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)=0
                    END DO
                  END DO
                END DO
              END DO
C**START 2-MODE INTEGRATION
              DO M4=1,MM4
                DO IRHS5=1,JCI5
                  DO ILHS5=1,JCI5
                    TEMP1(ILHS5,IRHS5)=0
                  END DO
                END DO
C**START 1-MODE INTEGRATION
                DO M5=1,MM5
                  IF(JCOUPL.GT.0)THEN
                    IF(IWHICH.GE.0.OR.MOLINC.LE.0)C=VP(M5,M4,M3,M2,M1)
                  ELSE
                    IF(IWHICH.GE.0.OR.MOLINC.LE.0)C=VPR(M5,M4,M3,M2,M1)
                  END IF
                  DO IRHS5=1,JCI5
                    DO ILHS5=1,JCI5
                      TEMP1(ILHS5,IRHS5)=TEMP1(ILHS5,IRHS5)+
     1                U0(ILHS5,IRHS5,M5)*C
                    END DO
                  END DO
                END DO
C**END 1-MODE INTEGRATION
                DO IRHS5=1,JCI5
                  DO IRHS4=1,JCI4
                    DO ILHS5=1,JCI5
                      DO ILHS4=1,JCI4
                        TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)=
     1                  TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)+
     2                  W0(ILHS4,IRHS4,M4)*TEMP1(ILHS5,IRHS5)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
C**END 2-MODE INTEGRATION
              DO IRHS5=1,JCI5
                DO IRHS4=1,JCI4
                  DO IRHS3=1,JCI3
                    DO ILHS5=1,JCI5
                      DO ILHS4=1,JCI4
                        DO ILHS3=1,JCI3
                          TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)=
     1                    TEMP3(ILHS3,ILHS4,ILHS5,IRHS3,IRHS4,IRHS5)+
     2                    Z0(ILHS3,IRHS3,M3)*
     3                    TEMP2(ILHS4,ILHS5,IRHS4,IRHS5)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
C**END 3-MODE INTEGRATION

C**NSIZE IS NO. UNIQUE INTEGRALS (5-DIM)
            DO IRHS=1,NSIZE
              NR1=IP5(IRHS,1)
              NR2=IP5(IRHS,2)
              NR3=IP5(IRHS,3)
              NR4=IP5(IRHS,4)
              NR5=IP5(IRHS,5)
              J0=IRHS*(IRHS-1)/2
              DO ILHS=1,IRHS
                NL1=IP5(ILHS,1)
                NL2=IP5(ILHS,2)
                NL3=IP5(ILHS,3)
                NL4=IP5(ILHS,4)
                NL5=IP5(ILHS,5)
                XA5(ILHS+J0)=XA5(ILHS+J0)+
     1          TEMP3(NL3,NL4,NL5,NR3,NR4,NR5)*X0(NL1,NR1,M1)*
     2          Y0(NL2,NR2,M2)
              END DO
            END DO
          END DO
C**END 4-MODE INTEGRATION
        END DO
C**END 5-MODE INTEGRATION
      END IF
      CALL MATOUT(XA5,XRA5,NSIZE,25)
      GO TO 7777
6666  CALL MATIN(XA5,XRA5,NSIZE,25)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI5

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
      JSR=1
      DO 9999 ISMR1=1,NVSYM
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NS.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NS.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JSR=-JSR
      END IF
      IF(NS.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NS.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NS.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NS.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSMR=ISMR1
        KCOFFR=ICOFF1
        KCSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
        LSMR=ISMR2
        LCOFFR=ICOFF2
        LCSZR=ICSZ2
        NLVALR=NCVAL(2,ISMR2)
      ELSE
        KSMR=ISMR2
        KCOFFR=ICOFF2
        KCSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
        LSMR=ISMR1
        LCOFFR=ICOFF1
        LCSZR=ICSZ1
        NLVALR=NCVAL(1,ISMR1)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NS.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NS.EQ.2)THEN
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NS.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NS.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NS.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NS.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NS.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        KSML=ISML1
        KCOFFL=ICOFF1
        KCSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
        LSML=ISML2
        LCOFFL=ICOFF2
        LCSZL=ICSZ2
        NLVALL=NCVAL(2,ISML2)
      ELSE
        KSML=ISML2
        KCOFFL=ICOFF2
        KCSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
        LSML=ISML1
        LCOFFL=ICOFF1
        LCSZL=ICSZ1
        NLVALL=NCVAL(1,ISML1)
      END IF

C**********************************************************
C**********************************************************
C**FIRST FIND TOTAL OF NON-ZERO TERMS IN K
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 5000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 5000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
5000        CONTINUE
          END DO
        END DO
        IF(NONZ1.GT.NONC1)NONC1=NONZ1
        GO TO 4444
      ELSE
        NONZ1=0
C**CONTRACTION SCHEME  '1' (K)
        DO IRH1=1,KCSZR
          IROFF=IRH1+KCOFFR
          DO ILH1=1,KCSZL
            ILOFF=ILH1+KCOFFL
C**OVERLAP OF REMAINING STATES FOR '1'
            IS1=1
            DO K=1,NKMOD1
              IF(IS1.EQ.0)GO TO 3000
              IGOT=0
              DO I=1,NC1
                IF(K.EQ.MC1(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC1(IROFF,K).NE.IPC1(ILOFF,K)))IS1=0
            END DO
            IF(IS1.EQ.0)GO TO 3000
C**OVERLAP OF REMAINING STATES FOR '1'
            NONZ1=NONZ1+1
            NONZC1(1,NONZ1)=IROFF
            NONZC1(2,NONZ1)=ILOFF
            NONZC1(3,NONZ1)=IRH1
            NONZC1(4,NONZ1)=ILH1
3000        CONTINUE
          END DO
        END DO
C**FIND UNIQUE TERMS
C       NONT1=1
C       NONZT1(1,1)=NONZC1(1,1)
C       NONZT1(2,1)=NONZC1(2,1)
C       DO I=2,NONZ1
C         IROFF=NONZC1(1,I)
C         ILOFF=NONZC1(2,I)
C         DO K=1,NONT1
C           MROFF=NONZT1(1,K)
C           MLOFF=NONZT1(2,K)
C           IF(IROFF.EQ.MROFF.AND.ILOFF.EQ.MLOFF)GO TO 3001
C         END DO
C         NONT1=NONT1+1
C         NONZT1(1,NONT1)=IROFF
C         NONZT1(2,NONT1)=ILOFF
3001      CONTINUE
C       END DO
C       IF(NONZ1.EQ.0)NONT1=0
      END IF
      IF(NONZ1.EQ.0)GO TO 5555
      NONZC1(3,NONZ1+1)=0

4444  CONTINUE
C**********************************************************
C**********************************************************
C**SECOND FIND TOTAL OF NON-ZERO TERMS IN L
C**SAVE INDICES FOR RHS AND LHS
C**********************************************************
C**********************************************************
      IF(IND.EQ.0)THEN
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M AND/OR/J)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 6000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 6000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
6000        CONTINUE
          END DO
        END DO
        IF(NONZ2.GT.NONC2)NONC2=NONZ2
        GO TO 5555
      ELSE
        NONZ2=0
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR M AND/OR/J)
        DO IRH2=1,LCSZR
          JROFF=IRH2+LCOFFR
          DO ILH2=1,LCSZL
            JLOFF=ILH2+LCOFFL
C**OVERLAP OF REMAINING STATES FOR '2'
            IS2=1
            DO K=1,NKMOD2
              IF(IS2.EQ.0)GO TO 4000
              IGOT=0
              DO I=1,NC2
                IF(K.EQ.MC2(I))IGOT=1
              END DO
              IF(IGOT.EQ.0.AND.
     1        (IPC2(JROFF,K).NE.IPC2(JLOFF,K)))IS2=0
            END DO
            IF(IS2.EQ.0)GO TO 4000
C**OVERLAP OF REMAINING STATES FOR '2'
            NONZ2=NONZ2+1
            NONZC2(1,NONZ2)=JROFF
            NONZC2(2,NONZ2)=JLOFF
            NONZC2(3,NONZ2)=IRH2
            NONZC2(4,NONZ2)=ILH2
4000        CONTINUE
          END DO
        END DO
C**FIND UNIQUE TERMS
C       NONT2=1
C       NONZT2(1,1)=NONZC2(1,1)
C       NONZT2(2,1)=NONZC2(2,1)
C       DO I=2,NONZ2
C         JROFF=NONZC2(1,I)
C         JLOFF=NONZC2(2,I)
C         DO K=1,NONT2
C           MROFF=NONZT2(1,K)
C           MLOFF=NONZT2(2,K)
C           IF(JROFF.EQ.MROFF.AND.JLOFF.EQ.MLOFF)GO TO 4001
C         END DO
C         NONT2=NONT2+1
C         NONZT2(1,NONT2)=JROFF
C         NONZT2(2,NONT2)=JLOFF
4001      CONTINUE
C       END DO
C       IF(NONZ2.EQ.0)NONT2=0
      END IF
      IF(NONZ2.EQ.0)GO TO 5555
      NONZC2(3,NONZ2+1)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF(ITIM5A.EQ.0)THEN
C       WRITE(IOUT,*)'KCSZR,KCSZL = ',KCSZR,KCSZL
C       WRITE(IOUT,*)'LCSZR,LCSZL = ',LCSZR,LCSZL
C       WRITE(IOUT,*)'KCSZR*KCSZL = ',KCSZR*KCSZL,
C    1  'LCSZR*LCSZL = ',LCSZR*LCSZL
C       WRITE(IOUT,*)'NONZ1,NONZ2 = ',NONZ1,NONZ2
C     END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C**********************************************************
C**********************************************************
C**THIRD SET UP FINAL MATRIX
C**********************************************************
C**********************************************************

C**CASE 1..K IS IN CONTRACTED SCHEME, L/N/M IS IN UNCONTRACTED SCHEME
      IF(ICCONT(1).GT.0)THEN

        IF(NONC1.LE.NONC2)GO TO 3003
C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
        IPREV=0
        DO INH1=1,NONZ1
          IRH1=NONZC1(3,INH1)
C         IF(INH1.NE.NONZ1)THEN
            INEXT=NONZC1(3,INH1+1)
C         ELSE
C           INEXT=0
C         END IF
          ILH1=NONZC1(4,INH1)
          IF(IPREV.NE.IRH1)THEN
            DO INH2=1,NONZ2
              DO IKL=1,NKVALL
                XK2(IKL,INH2)=0
              END DO
            END DO
          END IF
          IROFF=NONZC1(1,INH1)
          ILOFF=NONZC1(2,INH1)
          DO I=1,NC1
            NCR1(I)=IPC1(IROFF,MC1(I))
            NCL1(I)=IPC1(ILOFF,MC1(I))
          END DO

C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR/ M)
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            ILH2=NONZC2(4,INH2)
            JROFF=NONZC2(1,INH2)
            JLOFF=NONZC2(2,INH2)
            DO I=1,NC2
              NCR2(I)=IPC2(JROFF,MC2(I))
              NCL2(I)=IPC2(JLOFF,MC2(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP5(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP5(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1003
            END DO
1003        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP5(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP5(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2003
            END DO
2003        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA5(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP5(MR,1)
              NL1=IP5(ML,1)
              NR2=IP5(MR,2)
              NL2=IP5(ML,2)
              NR3=IP5(MR,3)
              NL3=IP5(ML,3)
              NR4=IP5(MR,4)
              NL4=IP5(ML,4)
              NR5=IP5(MR,5)
              NL5=IP5(ML,5)
              DO I=1,NP5(IND5)
                K=JP5(I,IND5,1)+1
                L=JP5(I,IND5,2)+1
                N=JP5(I,IND5,3)+1
                M=JP5(I,IND5,4)+1
                J=JP5(I,IND5,5)+1
                ZYX=ZYX+CP5(I,IND5)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)*XKAN(NL5,NR5,J,MOD5)
              END DO
C**ANALYTIC

            END IF
            DO IKL=1,NKVALL
              XK2(IKL,INH2)=XK2(IKL,INH2)+(XYZ+ZYX)*
     1        CFS1(ILH1,IKL,KEL,KSML)
            END DO
          END DO
          IF(IRH1.NE.INEXT)THEN
            DO INH2=1,NONZ2
              IRH2=NONZC2(3,INH2)
              ILH2=NONZC2(4,INH2)
              IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
                DO IKR=1,NKVALR
                  CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                  IRHS=(IKR-1)*LCSZR+IRH2
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO IKL=1,NKVALL
                    ILHS=(IKL-1)*LCSZL+ILH2
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(IKL,INH2)
                  END DO
                END DO
C**OFF-DIAGONAL BLOCK - END
              ELSE
C**DIAGONAL BLOCK - START
                DO IKR=1,NKVALR
                  CCCC=CFS1(IRH1,IKR,KEL,KSMR)
                  IRHS=(IKR-1)*LCSZR+IRH2
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO IKL=1,IKR
                    IF(IKL.EQ.IKR.AND.ILH2.GT.IRH2)GO TO 1004
                    ILHS=(IKL-1)*LCSZL+ILH2
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(IKL,INH2)
1004                CONTINUE
                  END DO
                END DO
C**DIAGONAL BLOCK - END
              END IF
            END DO
          END IF
          IPREV=IRH1
        END DO
        GO TO 5555

3003    CONTINUE
        DO IRH1=1,KCSZR
          DO ILH1=1,KCSZL
            XK1(ILH1,IRH1)=0
          END DO
        END DO
C**CONTRACTION SCHEME '2' (L AND/OR N AND/OR/ M)

        DO INH2=1,NONZ2
          IRH2=NONZC2(3,INH2)
          ILH2=NONZC2(4,INH2)
          JROFF=NONZC2(1,INH2)
          JLOFF=NONZC2(2,INH2)
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

C**CONTRACTION SCHEME  '1' (K)
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            ILH1=NONZC1(4,INH1)
            IROFF=NONZC1(1,INH1)
            ILOFF=NONZC1(2,INH1)
            DO I=1,NC1
              NCR1(I)=IPC1(IROFF,MC1(I))
              NCL1(I)=IPC1(ILOFF,MC1(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP5(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP5(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1000
            END DO
1000        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP5(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP5(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2000
            END DO
2000        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA5(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP5(MR,1)
              NL1=IP5(ML,1)
              NR2=IP5(MR,2)
              NL2=IP5(ML,2)
              NR3=IP5(MR,3)
              NL3=IP5(ML,3)
              NR4=IP5(MR,4)
              NL4=IP5(ML,4)
              NR5=IP5(MR,5)
              NL5=IP5(ML,5)
              DO I=1,NP5(IND5)
                K=JP5(I,IND5,1)+1
                L=JP5(I,IND5,2)+1
                N=JP5(I,IND5,3)+1
                M=JP5(I,IND5,4)+1
                J=JP5(I,IND5,5)+1
                ZYX=ZYX+CP5(I,IND5)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)*XKAN(NL5,NR5,J,MOD5)
              END DO
C**ANALYTIC
            END IF
            XK1(ILH1,IRH1)=XYZ+ZYX
          END DO

C************************************DGEMM (RHS)
          CALL DGEMM('N','N',KCSZL,NKVALR,KCSZR,1.0D0,XK1(1,1),
     &    KTEMP,
     &    CFS1(1,1,KEL,KSMR),ISIZXX,0.0D0,TEMP(1,1),KTEMP)

C************************************DGEMM (LHS)
          CALL DGEMM('T','N',NKVALL,NKVALR,KCSZL,1.0D0,
     &    CFS1(1,1,KEL,KSML),
     &           ISIZXX,TEMP,KTEMP,0.0D0,XCON2(1,1),MAXVAL)

          IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
            DO IKR=1,NKVALR
              IRHS=(IKR-1)*LCSZR+IRH2
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO IKL=1,NKVALL
                ILHS=(IKL-1)*LCSZL+ILH2
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(IKL,IKR)
              END DO
            END DO
C**OFF-DIAGONAL BLOCK - END
          ELSE
C**DIAGONAL BLOCK - START
            DO IKR=1,NKVALR
              IRHS=(IKR-1)*LCSZR+IRH2
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO IKL=1,IKR
                IF(IKL.EQ.IKR.AND.ILH2.GT.IRH2)GO TO 1002
                ILHS=(IKL-1)*LCSZL+ILH2
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(IKL,IKR)
1002            CONTINUE
              END DO
            END DO
C**DIAGONAL BLOCK - END
          END IF
        END DO

C**CASE 2..L/N/M IS IN CONTRACTED SCHEME, K IS IN UNCONTRACTED SCHEME
      ELSE

        IF(NONC1.GT.NONC2)GO TO 4004
C**CONTRACTION SCHEME INVOLVING 'MOD2' (L AND/OR N)
        IPREV=0
        DO INH2=1,NONZ2
          IRH2=NONZC2(3,INH2)
C         IF(INH2.NE.NONZ2)THEN
            INEXT=NONZC2(3,INH2+1)
C         ELSE
C           INEXT=0
C         END IF
          ILH2=NONZC2(4,INH2)
          IF(IPREV.NE.IRH2)THEN
            DO INH1=1,NONZ1
              DO ILL=1,NLVALL
                XK2(ILL,INH1)=0
              END DO
            END DO
          END IF
          JROFF=NONZC2(1,INH2)
          JLOFF=NONZC2(2,INH2)
          DO I=1,NC2
            NCR2(I)=IPC2(JROFF,MC2(I))
            NCL2(I)=IPC2(JLOFF,MC2(I))
          END DO

C**CONTRACTION SCHEME INVOLVING 'MOD1' (K)
          DO INH1=1,NONZ1
            IRH1=NONZC1(3,INH1)
            ILH1=NONZC1(4,INH1)
            IROFF=NONZC1(1,INH1)
            ILOFF=NONZC1(2,INH1)
            DO I=1,NC1
              NCR1(I)=IPC1(IROFF,MC1(I))
              NCL1(I)=IPC1(ILOFF,MC1(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP5(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP5(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1005
            END DO
1005        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP5(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP5(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2005
            END DO
2005        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA5(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP5(MR,1)
              NL1=IP5(ML,1)
              NR2=IP5(MR,2)
              NL2=IP5(ML,2)
              NR3=IP5(MR,3)
              NL3=IP5(ML,3)
              NR4=IP5(MR,4)
              NL4=IP5(ML,4)
              NR5=IP5(MR,5)
              NL5=IP5(ML,5)
              DO I=1,NP5(IND5)
                K=JP5(I,IND5,1)+1
                L=JP5(I,IND5,2)+1
                N=JP5(I,IND5,3)+1
                M=JP5(I,IND5,4)+1
                J=JP5(I,IND5,5)+1
                ZYX=ZYX+CP5(I,IND5)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)*XKAN(NL5,NR5,J,MOD5)
              END DO
C**ANALYTIC
            END IF
            DO ILL=1,NLVALL
              XK2(ILL,INH1)=XK2(ILL,INH1)+(XYZ+ZYX)*
     1        CFS2(ILH2,ILL,KEL,LSML)
            END DO
          END DO
          IF(IRH2.NE.INEXT)THEN
            DO INH1=1,NONZ1
              IRH1=NONZC1(3,INH1)
              ILH1=NONZC1(4,INH1)
              IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
                DO ILR=1,NLVALR
                  CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                  IRHS=(ILR-1)*KCSZR+IRH1
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO ILL=1,NLVALL
                    ILHS=(ILL-1)*KCSZL+ILH1
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(ILL,INH1)
                  END DO
                END DO
C**OFF-DIAGONAL BLOCK - END
              ELSE
C**DIAGONAL BLOCK - START
                DO ILR=1,NLVALR
                  CCCC=CFS2(IRH2,ILR,KEL,LSMR)
                  IRHS=(ILR-1)*KCSZR+IRH1
                  IROFF=IRHS+IOFFR
                  J0=IROFF*(IROFF-1)/2+IOFFL
                  DO ILL=1,ILR
                    IF(ILL.EQ.ILR.AND.ILH1.GT.IRH1)GO TO 2004
                    ILHS=(ILL-1)*KCSZL+ILH1
                    XA(ILHS+J0)=XA(ILHS+J0)+CCCC*XK2(ILL,INH1)
2004                CONTINUE
                  END DO
                END DO
C**DIAGONAL BLOCK - END
              END IF
            END DO
          END IF
          IPREV=IRH2
        END DO
        GO TO 5555

4004    CONTINUE
        DO IRH2=1,LCSZR
          DO ILH2=1,LCSZL
            XK1(ILH2,IRH2)=0
          END DO
        END DO

C**CONTRACTION SCHEME '1' (K)
        DO INH1=1,NONZ1
          IRH1=NONZC1(3,INH1)
          ILH1=NONZC1(4,INH1)
          JROFF=NONZC1(1,INH1)
          JLOFF=NONZC1(2,INH1)
          DO I=1,NC1
            NCR1(I)=IPC1(JROFF,MC1(I))
            NCL1(I)=IPC1(JLOFF,MC1(I))
          END DO

C**CONTRACTION SCHEME  '2' (L AND/OR N)
          DO INH2=1,NONZ2
            IRH2=NONZC2(3,INH2)
            ILH2=NONZC2(4,INH2)
            IROFF=NONZC2(1,INH2)
            ILOFF=NONZC2(2,INH2)
            DO I=1,NC2
              NCR2(I)=IPC2(IROFF,MC2(I))
              NCL2(I)=IPC2(ILOFF,MC2(I))
            END DO
C**FIND RHS INDEX
            DO IR=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCR1(I).NE.IP5(IR,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCR2(I).NE.IP5(IR,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 1001
            END DO
1001        CONTINUE
C**FIND LHS INDEX
            DO IL=1,NSIZE
              IGOT=1
              DO I=1,NC1
                IF(NCL1(I).NE.IP5(IL,JC1(I)))IGOT=0
              END DO
              DO I=1,NC2
                IF(NCL2(I).NE.IP5(IL,JC2(I)))IGOT=0
              END DO
              IF(IGOT.EQ.1)GO TO 2001
            END DO
2001        CONTINUE
C**GET MATRIX ELEMENT
            MR=IR
            ML=IL
            IF(IR.LT.IL)THEN
              MR=IL
              ML=IR
            END IF
            I=MR*(MR-1)/2+ML
            XYZ=XA5(I)
            ZYX=0
            IF(IWHICH.LT.0.AND.MOLINC.GT.0)THEN
C**ANALYTIC
              NR1=IP5(MR,1)
              NL1=IP5(ML,1)
              NR2=IP5(MR,2)
              NL2=IP5(ML,2)
              NR3=IP5(MR,3)
              NL3=IP5(ML,3)
              NR4=IP5(MR,4)
              NL4=IP5(ML,4)
              NR5=IP5(MR,5)
              NL5=IP5(ML,5)
              DO I=1,NP5(IND5)
                K=JP5(I,IND5,1)+1
                L=JP5(I,IND5,2)+1
                N=JP5(I,IND5,3)+1
                M=JP5(I,IND5,4)+1
                J=JP5(I,IND5,5)+1
                ZYX=ZYX+CP5(I,IND5)*XKAN(NL1,NR1,K,MOD1)*
     1          XKAN(NL2,NR2,L,MOD2)*XKAN(NL3,NR3,N,MOD3)*
     2          XKAN(NL4,NR4,M,MOD4)*XKAN(NL5,NR5,J,MOD5)
              END DO
C**ANALYTIC
            END IF
            XK1(ILH2,IRH2)=XYZ+ZYX
          END DO

C************************************DGEMM (RHS)
          CALL DGEMM('N','N',LCSZL,NLVALR,LCSZR,1.0D0,XK1(1,1),
     &    KTEMP,CFS2(1,1,KEL,LSMR),ISIZXX,0.0D0,TEMP(1,1),KTEMP)

C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NLVALL,NLVALR,LCSZL,1.0D0,
     &    CFS2(1,1,KEL,LSML),ISIZXX,TEMP,KTEMP,0.0D0,XCON2(1,1),MAXVAL)

          IF(ISML1.NE.ISMR1)THEN
C**OFF-DIAGONAL BLOCK - START
            DO ILR=1,NLVALR
              IRHS=(ILR-1)*KCSZR+IRH1
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO ILL=1,NLVALL
                ILHS=(ILL-1)*KCSZL+ILH1
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(ILL,ILR)
              END DO
            END DO
C**OFF-DIAGONAL BLOCK - END
          ELSE
C**DIAGONAL BLOCK - START
            DO ILR=1,NLVALR
              IRHS=(ILR-1)*KCSZR+IRH1
              IROFF=IRHS+IOFFR
              J0=IROFF*(IROFF-1)/2+IOFFL
              DO ILL=1,ILR
                IF(ILL.EQ.ILR.AND.ILH1.GT.IRH1)GO TO 2002
                ILHS=(ILL-1)*KCSZL+ILH1
                XA(ILHS+J0)=XA(ILHS+J0)+XCON2(ILL,ILR)
2002            CONTINUE
              END DO
            END DO
C**DIAGONAL BLOCK - END
          END IF
        END DO
      END IF
C*************************************

5555  CONTINUE

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE
      IF(ITIM5A.EQ.0)THEN
        IF(IND.EQ.1)THEN
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          ITIM5A=1
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE RDREC(XA5,XRA5,IRHS,INPC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 XA5(IRHS)
      REAL*4 XRA5(IRHS)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      IF(JCOUPL.GT.0)THEN
        READ(INPC)XA5
      ELSE
        READ(INPC)XRA5
        DO I=1,IRHS
          K=IRHS+1-I
          XA5(K)=XRA5(K)
        END DO
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE WTREC(XA5,XRA5,IRHS,IOUTC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 XA5(IRHS)
      REAL*4 XRA5(IRHS)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      IF(JCOUPL.GT.0)THEN
        WRITE(IOUTC)XA5
      ELSE
        DO I=1,IRHS
          XRA5(I)=XA5(I)
        END DO
        WRITE(IOUTC)XRA5
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VCMI1(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MM,IPL,IPR,
     1ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZEL,ISIZER,IP1,ISIZE1,XA1,
     2XRA1,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,VM,VMR,J21,IABC,MODINT,
     3NSL,NSR,CFSL,CFSR,ISIZXX,NVALX,KEL21,NVSYMX)
C****************************************
C**NOTE:   CALLED ONLY IF ICOUPC > 0
C****************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM,6)
      REAL*4 VMR(MM,6)
      REAL*4 XRA1(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H(NN,MM,3),XQ(MM)
      DIMENSION XA(ISIZEL,ISIZER)
      DIMENSION XA1(1),IP1(ISIZE1,1)
      DIMENSION IPL(ISIZMX,NNMODE),IPR(ISIZMX,NNMODE)
      DIMENSION XK(ICSIZE,ICSIZE),IPC(IPSIZE,1),TEMP(ICSIZE,NVAL)
      DIMENSION XCON(NVAL,NVAL)
      DIMENSION CFSL(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFSR(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      MD=MODINT(MOD1)
      CALL VDCMI1(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MM,MM/MD,IPL,IPR,
     1ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZEL,ISIZER,IP1,ISIZE1,XA1,
     2XRA1,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,VM,VMR,J21,IABC,MODINT,
     3NSL,NSR,CFSL,CFSR,ISIZXX,NVALX,KEL21,NVSYMX)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDCMI1(NMODE,NNMODE,MOD1,MODE1,H,XQ,NN,MH,MM,IPL,IPR,
     1ISIZMX,IPC,ICSIZE,IPSIZE,KCONT,XA,ISIZEL,ISIZER,IP1,ISIZE1,XA1,
     2XRA1,NSIZE,XK,TEMP,XCON,NVAL,ISTART,IEND,VM,VMR,J21,IABC,MODINT,
     3NSL,NSR,CFSL,CFSR,ISIZXX,NVALX,KEL21,NVSYMX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM,6)
      REAL*4 VMR(MM,6)
      REAL*4 XRA1(1)
      DIMENSION MODINT(NMODE)
      DIMENSION H(NN,MM,3),XQ(MM)
      DIMENSION XA(ISIZEL,ISIZER)
      DIMENSION XA1(1),IP1(ISIZE1,1)
      DIMENSION IPL(ISIZMX,NNMODE),IPR(ISIZMX,NNMODE)
      DIMENSION XK(ICSIZE,ICSIZE),IPC(IPSIZE,1),TEMP(ICSIZE,NVAL)
      DIMENSION XCON(NVAL,NVAL)
      DIMENSION CFSL(ISIZXX,NVALX,KEL21,NVSYMX)
      DIMENSION CFSR(ISIZXX,NVALX,KEL21,NVSYMX)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP

      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VCMI1'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
      END IF

C**FIRST DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'MOD1'
C**IF IT IS IN SCHEME 1.....
      NKMODE=ICONT(1)
      MCONT=2
C**....IT IS NOT IN SCHEME 2
      IF(KCONT.EQ.2)THEN
        NKMODE=ICONT(2)
        MCONT=1
      END IF

C***********************************ALGORITHM FROM V0MI1
      IF(ICONDP.NE.0)GO TO 6666

      IFACTC=INTFAC(NMODE,ICOUPC,1)
      IFACTL=INTFAC(NMODE,ICOUPL,1)
      IF(IWHICH.EQ.0)IFATCL=1

C***********************************************************

      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GT.0)READ(91)VM
      ELSE
        IF(ICOUPC.GT.0)READ(91)VMR
      END IF

C***********************************************************

      MD=MODINT(MOD1)

C***********************************************************

C**ONE-MODE COUPLING INTEGRAL FOR 'ACTIVE' SCHEME
CCCC  DO M=1,MM/MD
      DO M=1,MM
        IF(JCOUPL.GT.0)THEN
          TERM=VM(M,IABC)*IFACTC
        ELSE
          TERM=VMR(M,IABC)*IFACTC
        END IF
C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          NR=IP1(IRHS,1)
          X=TERM*H(NR,M,1)*MD
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP1(ILHS,1)
            Y=H(NL,M,1)
            XA1(ILHS+J0)=XA1(ILHS+J0)+Y*X
          END DO
        END DO
      END DO
      CALL MATOUT(XA1,XRA1,NSIZE,21)
      GO TO 7777
6666  CALL MATIN(XA1,XRA1,NSIZE,21)
7777  CONTINUE

C***********************************ALGORITHM FROM VCI1

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**RHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFR=0
      JSR=1
      DO 9999 ISMR1=1,NVSYM
      KSR=ISMR1/5
      LSR=(-1)**(KSR+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISMR1)
      IF(NSR.EQ.1)THEN
        ISMR2=ISMR1
      END IF
      IF(NSR.EQ.2)THEN
        ISMR2=ISMR1+JSR
        JRS=-JRS
      END IF
      IF(NSR.EQ.3)THEN
        ISMR2=ISMR1+2*JSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NSR.EQ.4)THEN
C       ISMR2=NVSYM+1-ISMR1
        ISMR2=5+NVSYM*KSR-ISMR1
      END IF
      IF(NSR.EQ.5)THEN
        ISMR2=ISMR1+4*LSR
      END IF
      IF(NSR.EQ.6)THEN
        ISMR2=ISMR1+JSR+4*LSR
        JSR=-JSR
      END IF
      IF(NSR.EQ.7)THEN
        ISMR2=ISMR1+2*JSR+4*LSR
        JJJ=MOD(ISMR1,2)+1
        JSR=JSR*(-1)**JJJ
      END IF
      IF(NSR.EQ.8)THEN
        ISMR2=5+NVSYM*KSR-ISMR1+4*LSR
      END IF
      IF(ICSZ1.EQ.0)GO TO 9998
      ICSZ2=ISIZC(2,ISMR2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9998
C**OFFSETS FOR XK
      DO JJJ=1,ISMR1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISMR2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT.EQ.1)THEN
        ISMR=ISMR1
        ICOFFR=ICOFF1
        ICSZR=ICSZ1
        NKVALR=NCVAL(1,ISMR1)
      ELSE
        ISMR=ISMR2
        ICOFFR=ICOFF2
        ICSZR=ICSZ2
        NKVALR=NCVAL(2,ISMR2)
      END IF

C******************************************************
C**START SYMMETRY LOOP FOR TOTAL NUMBER OF SPECIES
C******************************************************
C**LHS OFFSET FOR FINAL MATRIX (XA AND IP)
      IOFFL=0
      JSL=1
      DO 9997 ISML1=1,ISMR1
      KSL=ISML1/5
      LSL=(-1)**(KSL+2)
      ICOFF1=0
      ICOFF2=0
      ICSZ1=ISIZC(1,ISML1)
      IF(NSL.EQ.1)THEN
        ISML2=ISML1
      END IF
      IF(NSL.EQ.2)THEN 
        ISML2=ISML1+JSL
        JSL=-JSL
      END IF
      IF(NSL.EQ.3)THEN
        ISML2=ISML1+2*JSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NSL.EQ.4)THEN
C       ISML2=NVSYM+1-ISML1
        ISML2=5+NVSYM*KSL-ISML1
      END IF
      IF(NSL.EQ.5)THEN
        ISML2=ISML1+4*LSL
      END IF
      IF(NSL.EQ.6)THEN
        ISML2=ISML1+JSL+4*LSL
        JSL=-JSL
      END IF
      IF(NSL.EQ.7)THEN
        ISML2=ISML1+2*JSL+4*LSL
        JJJ=MOD(ISML1,2)+1
        JSL=JSL*(-1)**JJJ
      END IF
      IF(NSL.EQ.8)THEN
        ISML2=5+NVSYM*KSL-ISML1+4*LSL
      END IF
      IF(ICSZ1.EQ.0)GO TO 9996
      ICSZ2=ISIZC(2,ISML2)
C**SYMMETRY BLOCK FOR NON-ACTIVE SCHEME MAY NOT EXIST
      IF(ICSZ2.EQ.0)GO TO 9996
C**OFFSETS FOR XK
      DO JJJ=1,ISML1-1
        ICOFF1=ICOFF1+ISIZC(1,JJJ)
      END DO
      DO JJJ=1,ISML2-1
        ICOFF2=ICOFF2+ISIZC(2,JJJ)
      END DO
      IF(KCONT1.EQ.1)THEN
        ISML=ISML1
        ICOFFL=ICOFF1
        ICSZL=ICSZ1
        NKVALL=NCVAL(1,ISML1)
      ELSE
        ISML=ISML2
        ICOFFL=ICOFF2
        ICSZL=ICSZ2
        NKVALL=NCVAL(2,ISML2)
      END IF

C**START AND END OF EACH SYMMETRY BLOCK
      JSTART=1
      JEND=NCSIZE(ISM1)

      DO IRHS=1,ICSZR
        IROFF=IRHS+ICOFFR
        NR=IPC(IROFF,MODE1)
C**FIND RHS INDEX (TRIVIAL CASE)
        DO IR=1,NSIZE
          IF(NR.EQ.IP1(IR,1))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ICSZL
          ILOFF=ILHS+ICOFFL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NKMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.(IPC(IROFF,K).NE.IPC(ILOFF,K)))
     1      IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL=IPC(ILOFF,MODE1)
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
          XYZ=XA1(I)
3000      CONTINUE
          XK(ILHS,IRHS)=XYZ*IS
          XK(IRHS,ILHS)=XK(ILHS,IRHS)
        END DO
      END DO

C************************************************************
C************************************DGEMM (RHS)
        CALL DGEMM('N','N',ICSZL,NKVALR,ICSZR,1.0D0,XK(1,1),ICSIZE,
     &      CFSR(1,1,KEL,ISMR),ISIZXX,0.0D0,TEMP,ICSIZE)
C************************************DGEMM (LHS)
         CALL DGEMM('T','N',NKVALL,NKVALR,ICSZL,1.0D0,
     &      CFSL(1,1,KEL,ISML),ISIZXX,TEMP,ICSIZE,0.0D0,XCON(1,1),NVAL)
C************************************************************

      DO IRHS=1,ISIZER
        IROFF=IRHS+IOFFR
C**RHS INDEX FOR ACTIVE 'MOD1'
        NR=IPR(IROFF,KCONT)
        DO ILHS=1,ISIZEL
          ILOFF=ILHS+IOFFL
C**OVERLAP OF CONTRACTION FUNCTIONS IN 'NON-ACTIVE' SCHEME
          IF(IPR(IRHS,MCONT).NE.IPL(ILHS,MCONT))GO TO 4000
C**LHS INDEX FOR ACTIVE 'MOD1'
          NL=IPL(ILHS,KCONT)
C**GET MATRIX ELEMENT
          XYZ=XCON(NL,NR)
          XA(ILHS+ILOFF,IRHS+IROFF)=XA(ILHS+ILOFF,IRHS+IROFF)+
     1    XYZ*IHERM
4000      CONTINUE
        END DO
      END DO

C**UPDATE LHS POINTER FOR XA AND IP
      IOFFL=IOFFL+NCSIZE(ISML1)
9996  CONTINUE
9997  CONTINUE

C**UPDATE RHS POINTER FOR XA AND IP
      IOFFR=IOFFR+NCSIZE(ISMR1)
9998  CONTINUE
9999  CONTINUE

      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM1A=1
      END IF
      RETURN
      END
