C***********************************************************
C***********************************************************
C**MEMO
C***********************************************************
C***********************************************************
      PARAMETER(MAXSIZ=260000000,MSEG=200)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FILASS/IOUT,INP,MOUTIN
      COMMON/MOLPRO/MOLPRO,MOUT,MINP
      COMMON/CWORK/W(MAXSIZ)
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
      COMMON/TRANSF/LTRAN
      character*4 tpdn
      tpdn='temp'

C**TEMPORARY
      OPEN(7,FILE=''//tpdn//'/molpro7',FORM='FORMATTED',
     1status='UNKNOWN')
      OPEN(20,FILE=''//tpdn//'/lanczos20',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(21,FILE=''//tpdn//'/vibints21',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(22,FILE=''//tpdn//'/vibints22',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(23,FILE=''//tpdn//'/vibints23',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(24,FILE=''//tpdn//'/vibints24',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(25,FILE=''//tpdn//'/vibints25',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(26,FILE=''//tpdn//'/vibints26',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(29,FILE=''//tpdn//'/scratch29',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(30,FILE=''//tpdn//'/scratch30',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(31,FILE=''//tpdn//'/rotints31',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(32,FILE=''//tpdn//'/rotints32',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(33,FILE=''//tpdn//'/rotints33',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(34,FILE=''//tpdn//'/rotints34',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(35,FILE=''//tpdn//'/rotints35',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(36,FILE=''//tpdn//'/rotints36',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(37,FILE=''//tpdn//'/rotints37',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(38,FILE=''//tpdn//'/rotints38',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(39,FILE=''//tpdn//'/rotints39',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(40,FILE=''//tpdn//'/SCFenergy40',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(41,FILE=''//tpdn//'/restart41',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(42,FILE=''//tpdn//'/restart42',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(43,FILE=''//tpdn//'/restart43',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(44,FILE=''//tpdn//'/restart44',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(45,FILE=''//tpdn//'/CIenergy45',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(46,FILE=''//tpdn//'/lanczos46',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(47,FILE=''//tpdn//'/lanczos47',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(48,FILE=''//tpdn//'/lanczos48',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(49,FILE=''//tpdn//'/lanczos49',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(50,FILE=''//tpdn//'/lanczos50',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(51,FILE=''//tpdn//'/lanczos51',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(52,FILE=''//tpdn//'/lanczos52',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(53,FILE=''//tpdn//'/lanczos53',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(54,FILE=''//tpdn//'/lanczos54',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(55,FILE=''//tpdn//'/lanczos55',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(56,FILE=''//tpdn//'/lanczos56',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(57,FILE=''//tpdn//'/lanczos57',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(58,FILE=''//tpdn//'/scratch58',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(59,FILE=''//tpdn//'/scratch59',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(60,FILE=''//tpdn//'/dump60',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(61,FILE=''//tpdn//'/rotation61',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(62,FILE=''//tpdn//'/rotation62',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(63,FILE=''//tpdn//'/rotation63',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(64,FILE=''//tpdn//'/rotation64',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(65,FILE=''//tpdn//'/scratch65',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(70,FILE=''//tpdn//'/scratch70',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(71,FILE=''//tpdn//'/potential71',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(72,FILE=''//tpdn//'/potential72',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(73,FILE=''//tpdn//'/potential73',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(74,FILE=''//tpdn//'/potential74',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(75,FILE=''//tpdn//'/potential75',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(76,FILE=''//tpdn//'/potential76',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(81,FILE=''//tpdn//'/coriolis81',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(82,FILE=''//tpdn//'/coriolis82',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(83,FILE=''//tpdn//'/coriolis83',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(84,FILE=''//tpdn//'/coriolis84',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(91,FILE=''//tpdn//'/coriolis91',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(92,FILE=''//tpdn//'/coriolis92',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(93,FILE=''//tpdn//'/coriolis93',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(94,FILE=''//tpdn//'/coriolis94',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(161,FILE=''//tpdn//'/rotation161',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(162,FILE=''//tpdn//'/rotation162',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(163,FILE=''//tpdn//'/rotation163',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(164,FILE=''//tpdn//'/rotation164',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(171,FILE=''//tpdn//'/potential171',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(172,FILE=''//tpdn//'/potential172',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(173,FILE=''//tpdn//'/potential173',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(174,FILE=''//tpdn//'/potential174',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(181,FILE=''//tpdn//'/coriolis181',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(182,FILE=''//tpdn//'/coriolis182',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(183,FILE=''//tpdn//'/coriolis183',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(184,FILE=''//tpdn//'/coriolis184',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(191,FILE=''//tpdn//'/coriolis191',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(192,FILE=''//tpdn//'/coriolis192',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(193,FILE=''//tpdn//'/coriolis193',FORM='UNFORMATTED',
     1status='UNKNOWN')
      OPEN(194,FILE=''//tpdn//'/coriolis194',FORM='UNFORMATTED',
     1status='UNKNOWN')
C**TEMPORARY
      LTRAN=MAXSIZ
      INP=1
      IOUT=2
      MINP=3
      MOUT=4
      MOUTIN=7
C     CALL TIMIT(2)
      NADD=MSEG
      DO 123 N=1,NADD
      KADD(N)=0
123   LADD(N)=0
      DO 9999 N=1,MAXSIZ
9999  W(N)=0
      NSEG=0
      KFREE=MAXSIZ
      LFREE=1
      KINF=MAXSIZ
      CALL VSCF(W)
      WRITE(IOUT,66) MAXSIZ-KINF,MAXSIZ
 66   FORMAT('USED MEMORY:',I10,3X,'OF',I10)
      CALL TIMIT(4)
      STOP 'END OF MULTIMODE'
      END
C***********************************************************
C***********************************************************
      SUBROUTINE MEMO(MM,LA1,KA1,LA2,KA2,LA3,KA3,LA4,KA4,LA5,KA5)
      PARAMETER(MSEG=200)
      COMMON/FILASS/IOUT,INP
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
100   FORMAT(1X,'TOO MANY ENTRIES TO MEMO',/)
101   FORMAT(/,1X,'NO. ARRAYS MUST BE NON-ZERO',/)
102   FORMAT(/,1X,'ERROR IN MEMO CALLING SEQUENCE',/)
      M=MM
      IF(M.EQ.0)THEN
      WRITE(IOUT,101)
      WRITE(IOUT,102)
      STOP 'NO. ARRAYS ZERO'
      END IF
      IF(IABS(M).GT.5)THEN
      WRITE(IOUT,100)
      WRITE(IOUT,102)
      STOP 'TOO MANY ARRAYS'
      END IF
      MADD=1
 1    IF(IABS(MM).GE.1) CALL MEM1(M,LA1,KA1)
      MADD=2
 2    IF(IABS(MM).GE.2) CALL MEM1(M,LA2,KA2)
      MADD=3
 3    IF(IABS(MM).GE.3) CALL MEM1(M,LA3,KA3)
      MADD=4
 4    IF(IABS(MM).GE.4) CALL MEM1(M,LA4,KA4)
      MADD=5
 5    IF(IABS(MM).GE.5) CALL MEM1(M,LA5,KA5)
      KINF=MIN0(KINF,KFREE)
      RETURN
      END
C***************************************************
C***************************************************
      SUBROUTINE MEM1(M,LA,KA)
C.....RESERVES OR FREES MEMORY OF LENGTH KA AT ADDRESS LA
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(MSEG=200)
      COMMON/FILASS/IOUT,INP
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
      COMMON/CWORK/W(1)
500   FORMAT(/,1X,'ARRAY NUMBER ',I3,' ADDRESS ',I8,' SIZE ',I8,
     1' ALREADY EXISTS',/)
501   FORMAT(/,1X,'ARRAY TO BE ADDED MUST HAVE POSITIVE LENGTH',/)
502   FORMAT(/,1X,'MAXIMUM NO. ARRAYS REACHED.....',
     1'ADJUST MSEG PARAMETER (C.F. ALSO LADD,KADD)',/)
503   FORMAT(/,1X,'MEMORY EXHAUSTED.....',/,
     1         1X,I10,' AVAILABLE',4X,I10,' REQUESTED',4X,
     2         'REQUIRES EXTRA ',I10,/,1X,'ADJUST MAXSIZ PARAMETER',/)

504   FORMAT(/,1X,'ERROR IN MEMO CALLING SEQUENCE',/)
505   FORMAT(/,1X,'ARRAY ADDRESS ',I8,' SIZE ',I8,' NOT FOUND',/)
506   FORMAT(/,1X,'ARRAY TO BE ADDED MUST HAVE ZERO ADDRESS',/)
507   FORMAT(/,1X,'ARRAY TO BE DELETED HAS NON-POSITIVE ADDRESS ',I8,/)
508   FORMAT(/,1X,'ARRAY AT ADDRESS ',I8,' SIZE ',I8,' DOES NOT EXIST',
     1/)
509   FORMAT(/,1X,'ARRAY NUMBER ',I3,' ADDRESS ',I8,' LENGTH ',I8,
     1' HAS WRONG SIZE ',I8,/)
510   FORMAT(/,1X,'ARRAY ADDRESS NAME NOT IN COMMON/CADDR/',/)
511   FORMAT(/,1X,'ERROR IN MEMO SEQUENCE NUMBER ',I2,/)
      IF(M.LT.0) GOTO 100
C....................................CHECK FOR AVAILABLE MEMORY
      IF(LA.NE.0)THEN
        WRITE(IOUT,511)MADD
        WRITE(IOUT,506)
        NUM=0
        DO 1 I=1,NADD
        IF(LADD(I).NE.LA)GO TO 1
        NUM=I
1       CONTINUE
        IF(NUM.EQ.0)THEN
          WRITE(IOUT,505)LA,KA
          WRITE(IOUT,504)
        ELSE
          WRITE(IOUT,500)NUM,LA,KA
        END IF
        STOP 'ARRAY ADDRESS NON-ZERO'
      END IF
C***********
      IF(KA.LE.0)THEN
        WRITE(IOUT,511)MADD
        WRITE(IOUT,501)
        STOP 'ARRAY LENGTH NON-POSITIVE'
      END IF
C***********
      IF(KA.GT.KFREE)THEN
        WRITE(IOUT,511)MADD
        WRITE(IOUT,503)KFREE,KA,KA-KFREE
        STOP 'MEMORY EXHAUSTED'
      END IF
C***********
      IF(NSEG.EQ.MSEG)THEN
        WRITE(IOUT,511)MADD
        WRITE(IOUT,502)
        STOP 'TOO MANY ARRAYS'
      END IF
C***********
      NSEG=NSEG+1
      LA=LFREE
      NUM=0
      DO 2 I=1,NADD
      IF(LADD(I).NE.LA)GO TO 2
      NUM=I
2     CONTINUE
      IF(NUM.EQ.0)THEN
        WRITE(IOUT,511)MADD
        WRITE(IOUT,510)
        WRITE(IOUT,504)
        STOP 'ARRAY NOT FOUND'
      ELSE
        KADD(NUM)=KA
        DO I=1,KA
          W(LA-1+I)=0
        END DO
      END IF
      KFREE=KFREE-KA
      LFREE=LFREE+KA
      RETURN
C....................................FIND SEGMENT TO BE FREED
100   CONTINUE
      IF(LA.LE.0)THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,507)LA
      STOP 'ARRAY ADDRESS NOT POSITIVE'
      END IF
C***********
      NUM=0
      DO 3 I=1,NADD
      IF(LADD(I).NE.LA)GO TO 3
      NUM=I
3     CONTINUE
      IF(NUM.EQ.0)THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,508)LA,KA
      WRITE(IOUT,510)
      WRITE(IOUT,504)
      STOP 'ARRAY NOT FOUND'
      ELSE
C***********
      IF(KA.NE.KADD(NUM))THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,509)NUM,LA,KADD(NUM),KA
      WRITE(IOUT,504)
      STOP 'ARRAY HAS WRONG SIZE'
      END IF
      END IF
C***********
      KTOT=0
      DO 4 I=1,NADD
      IF(LADD(I).GT.LA)KTOT=KTOT+KADD(I)
4     CONTINUE
      IF(KTOT.GT.0)THEN
      DO 5 I=1,KTOT
      W(LA-1+I)=W(LA-1+I+KA)
5     CONTINUE
      END IF
      DO 6 I=1,NADD
      IF(LADD(I).GT.LA)LADD(I)=LADD(I)-KA
6     CONTINUE
      LA=0
      LFREE=LFREE-KA
      KFREE=KFREE+KA
      NSEG=NSEG-1
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE MEMBUG(MM,LA1,KA1,LA2,KA2,LA3,KA3,LA4,KA4,LA5,KA5)
      PARAMETER(MSEG=200)
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
      COMMON/CINCR/INCR
      COMMON/CLBUG/LBUG(MSEG)
      COMMON/CKBUG/KBUG(MSEG)
100   FORMAT(1X,'TOO MANY ENTRIES TO MEMBUG',/)
101   FORMAT(/,1X,'NO. ARRAYS MUST BE NON-ZERO',/)
102   FORMAT(/,1X,'ERROR IN MEMO CALLING SEQUENCE',/)
103   FORMAT(/,1X,'MEMO DEBUG FACILITY (CALLING ROUTINE)',/)
104   FORMAT(/,1X,'MEMO DEBUG FACILITY (CALLED ROUTINE)',/)
105   FORMAT(/,1X,'MEMO DEBUG FACILITY (INITIALISE)',/)
      IF(IPRINT.GT.1)THEN
        IF(MM.EQ.0)WRITE(IOUT,105)
        IF(MM.GT.0)WRITE(IOUT,103)
        IF(MM.LT.0)WRITE(IOUT,104)
      END IF
      M=MM
      IF(M.EQ.0)THEN
        INCR=0
        DO I=1,MSEG
          LBUG(I)=0
          KBUG(I)=0
        END DO
        RETURN
      END IF
      IF(IABS(M).GT.5)THEN
      WRITE(IOUT,100)
      WRITE(IOUT,102)
      STOP 'TOO MANY ARRAYS'
      END IF
      MADD=1
 1    IF(IABS(MM).GE.1) CALL MEM2(M,LA1,KA1)
      MADD=2
 2    IF(IABS(MM).GE.2) CALL MEM2(M,LA2,KA2)
      MADD=3
 3    IF(IABS(MM).GE.3) CALL MEM2(M,LA3,KA3)
      MADD=4
 4    IF(IABS(MM).GE.4) CALL MEM2(M,LA4,KA4)
      MADD=5
 5    IF(IABS(MM).GE.5) CALL MEM2(M,LA5,KA5)
      KINF=MIN0(KINF,KFREE)
      RETURN
      END
C***************************************************
C***************************************************
      SUBROUTINE MEM2(MM,LA,KA)
C.....FINDS ARRAY AT ADDRESS LA, OFFSET KA (CALLING ROUTINE)
C.....FINDS ARRAY AT INDEX LA, SIZE KA (CALLED ROUTINE)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(MSEG=200)
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
      COMMON/CWORK/W(1)
      COMMON/CINCR/INCR
      COMMON/CLBUG/LBUG(MSEG)
      COMMON/CKBUG/KBUG(MSEG)
      COMMON/TRANSF/LTRAN
500   FORMAT(/,1X,'ARRAY NUMBER ',I3,' ADDRESS ',I8,' SIZE ',I8,
     1' ALREADY EXISTS',/)
501   FORMAT(/,1X,'ARRAY TO BE ADDED MUST HAVE POSITIVE LENGTH',/)
502   FORMAT(/,1X,'MAXIMUM NO. ARRAYS REACHED.....',
     1'ADJUST MSEG PARAMETER (C.F. ALSO LADD,KADD)',/)
503   FORMAT(/,1X,'MEMORY EXHAUSTED.....',/,
     1         1X,I10,' AVAILABLE',4X,I10,' REQUESTED',4X,
     2         'REQUIRES EXTRA ',I10,/,1X,'ADJUST MAXSIZ PARAMETER',/)

504   FORMAT(/,1X,'ERROR IN MEMO CALLING SEQUENCE',/)
505   FORMAT(/,1X,'ARRAY ADDRESS ',I8,' SIZE ',I8,' NOT FOUND',/)
506   FORMAT(/,1X,'ARRAY TO BE ADDED MUST HAVE ZERO ADDRESS',/)
C****************************************************************
507   FORMAT(/,1X,'ARRAY HAS NON-POSITIVE ADDRESS ',I8,/)
508   FORMAT(/,1X,'ARRAY AT ADDRESS ',I8,' DOES NOT EXIST',
     1/)
509   FORMAT(1X,'ARRAY NUMBER ',I3,' INDEX ',I3,' OFFSET ',I8,
     1' ACTUAL LENGTH ',I8)
510   FORMAT(1X,'ARRAY NUMBER ',I3,' INDEX ',I3,' REQUIRED LENGTH ',I8,
     1' ACTUAL LENGTH ',I8)
511   FORMAT(/,1X,'ERROR IN MEMO SEQUENCE NUMBER ',I2,/)
512   FORMAT(/,1X,'REQUESTED SIZE FOR INDEX ',I3,' TOO BIG',/)
513   FORMAT(/,1X,'CORE COULD BE OVERWRITTEN ',I3,' INDEX ',I3,/)
C****************************************************************
C**FOR CALLING ROUTINE, LA IS START ADDRESS
      IF(MM.GT.0)GO TO 100
C**FOR CALLED ROUTINE, LA IS INDEX NUMBER
C.............FIND SEGMENT AT INDEX LA, SIZE KA
      IF(IPRINT.GT.1)WRITE(IOUT,510)LBUG(LA),LA,KA,KBUG(LA)
      IF(KA.GT.KBUG(LA))THEN
        WRITE(IOUT,512)LA
        STOP 'ARRAY DIMENSION TOO BIG'
      END IF
      RETURN
C.............FIND SEGMENT AT ADRESS LA, OFFSET KA
100   CONTINUE
      IF(LA.LE.0)THEN
        WRITE(IOUT,511)MADD
        WRITE(IOUT,507)LA
        RETURN
      END IF
C***********
      NUM=0
      DO 3 I=1,NADD
      IF(LADD(I).NE.LA)GO TO 3
      NUM=I
3     CONTINUE
      IF(NUM.EQ.0)THEN
        WRITE(IOUT,511)MADD
        WRITE(IOUT,508)LA
        RETURN
      ELSE
C***********
        INCR=INCR+1
        IF(INCR.GT.MSEG)STOP 'TOO MANY DEBUG ARRAYS'
        IF(LA+KADD(NUM).GT.LTRAN)THEN
          WRITE(IOUT,513)NUM,INCR
          STOP 'DANGER...OVERWRITE'
        END IF
        IF(IPRINT.GT.1)WRITE(IOUT,509)NUM,INCR,KA,KADD(NUM)-KA
        LBUG(INCR)=NUM
        KBUG(INCR)=KADD(NUM)-KA
        RETURN
      END IF
C***********
      WRITE(IOUT,504)
      RETURN
      END
