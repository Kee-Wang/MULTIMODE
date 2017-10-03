C***********************************************************
C***********************************************************
C**MEMO
C***********************************************************
C***********************************************************
C**16 GB
C      PARAMETER(MAXSIZ=2000000000,MSEG=200)
C**8 GB
C      PARAMETER(MAXSIZ=1000000000,MSEG=200)
C**4 GB
C      PARAMETER(MAXSIZ=500000000,MSEG=200)
C**2 GB
      PARAMETER(MAXSIZ=260000000,MSEG=200)
      IMPLICIT REAL*8 (A-H,O-Z)
 
      character*4 tpdn
      character*40 infile,outfile
      logical file_exists

      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7,INP8,
     1INP9,INP10
      COMMON/MOLPRO/MOLPRO,MOUT,MINP
      COMMON/CWORK/W(MAXSIZ)
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
      COMMON/TRANSF/LTRAN

      tpdn='temp'

      call system('mkdir '//tpdn)

      OPEN(3,FILE=''//tpdn//'/molpro3',
     1FORM='FORMATTED',status='UNKNOWN')
      OPEN(4,FILE=''//tpdn//'/molpro4',
     1FORM='FORMATTED',status='UNKNOWN')
      OPEN(7,FILE=''//tpdn//'/molpro7',
     1FORM='FORMATTED',status='UNKNOWN')
      OPEN(10,FILE=''//tpdn//'/restart10',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(11,FILE=''//tpdn//'/restart11',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(20,FILE=''//tpdn//'/lanczos20',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(21,FILE=''//tpdn//'/vibints21',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(22,FILE=''//tpdn//'/vibints22',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(23,FILE=''//tpdn//'/vibints23',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(24,FILE=''//tpdn//'/vibints24',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(25,FILE=''//tpdn//'/vibints25',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(26,FILE=''//tpdn//'/vibints26',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(29,FILE=''//tpdn//'/scratch29',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(30,FILE=''//tpdn//'/scratch30',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(31,FILE=''//tpdn//'/rotints31',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(32,FILE=''//tpdn//'/rotints32',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(33,FILE=''//tpdn//'/rotints33',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(34,FILE=''//tpdn//'/rotints34',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(35,FILE=''//tpdn//'/rotints35',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(36,FILE=''//tpdn//'/rotints36',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(37,FILE=''//tpdn//'/rotints37',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(38,FILE=''//tpdn//'/rotints38',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(39,FILE=''//tpdn//'/rotints39',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(40,FILE=''//tpdn//'/SCFenergy40',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(41,FILE=''//tpdn//'/restart41',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(42,FILE=''//tpdn//'/restart42',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(43,FILE=''//tpdn//'/restart43',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(44,FILE=''//tpdn//'/restart44',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(45,FILE=''//tpdn//'/CIenergy45',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(46,FILE=''//tpdn//'/lanczos46',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(47,FILE=''//tpdn//'/lanczos47',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(48,FILE=''//tpdn//'/lanczos48',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(49,FILE=''//tpdn//'/lanczos49',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(50,FILE=''//tpdn//'/lanczos50',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(51,FILE=''//tpdn//'/lanczos51',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(52,FILE=''//tpdn//'/lanczos52',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(53,FILE=''//tpdn//'/lanczos53',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(54,FILE=''//tpdn//'/lanczos54',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(55,FILE=''//tpdn//'/lanczos55',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(56,FILE=''//tpdn//'/lanczos56',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(57,FILE=''//tpdn//'/lanczos57',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(58,FILE=''//tpdn//'/scratch58',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(59,FILE=''//tpdn//'/scratch59',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(60,FILE=''//tpdn//'/dump60',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(61,FILE=''//tpdn//'/rotation61',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(62,FILE=''//tpdn//'/rotation62',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(63,FILE=''//tpdn//'/rotation63',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(64,FILE=''//tpdn//'/rotation64',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(65,FILE=''//tpdn//'/scratch65',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(70,FILE=''//tpdn//'/scratch70',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(71,FILE=''//tpdn//'/potential71',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(72,FILE=''//tpdn//'/potential72',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(73,FILE=''//tpdn//'/potential73',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(74,FILE=''//tpdn//'/potential74',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(75,FILE=''//tpdn//'/potential75',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(76,FILE=''//tpdn//'/potential76',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(81,FILE=''//tpdn//'/coriolis81',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(82,FILE=''//tpdn//'/coriolis82',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(83,FILE=''//tpdn//'/coriolis83',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(84,FILE=''//tpdn//'/coriolis84',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(85,FILE=''//tpdn//'/scratch85',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(86,FILE=''//tpdn//'/coriolis86',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(91,FILE=''//tpdn//'/coriolis91',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(92,FILE=''//tpdn//'/coriolis92',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(93,FILE=''//tpdn//'/coriolis93',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(94,FILE=''//tpdn//'/coriolis94',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(95,FILE=''//tpdn//'/coriolis95',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(96,FILE=''//tpdn//'/intense96',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(97,FILE=''//tpdn//'/intense97',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(98,FILE=''//tpdn//'/intense98',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(99,FILE=''//tpdn//'/intense99',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(101,FILE=''//tpdn//'/intense101',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(102,FILE=''//tpdn//'/intense102',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(103,FILE=''//tpdn//'/intense103',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(104,FILE=''//tpdn//'/intense104',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(105,FILE=''//tpdn//'/intense105',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(106,FILE=''//tpdn//'/intense106',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(111,FILE=''//tpdn//'/intense111',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(112,FILE=''//tpdn//'/intense112',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(113,FILE=''//tpdn//'/intense113',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(114,FILE=''//tpdn//'/intense114',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(115,FILE=''//tpdn//'/intense115',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(151,FILE=''//tpdn//'/intense151',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(152,FILE=''//tpdn//'/intense152',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(153,FILE=''//tpdn//'/intense153',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(154,FILE=''//tpdn//'/intense154',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(155,FILE=''//tpdn//'/intense155',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(161,FILE=''//tpdn//'/rotation161',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(162,FILE=''//tpdn//'/rotation162',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(163,FILE=''//tpdn//'/rotation163',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(164,FILE=''//tpdn//'/rotation164',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(171,FILE=''//tpdn//'/potential171',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(172,FILE=''//tpdn//'/potential172',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(173,FILE=''//tpdn//'/potential173',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(174,FILE=''//tpdn//'/potential174',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(181,FILE=''//tpdn//'/coriolis181',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(182,FILE=''//tpdn//'/coriolis182',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(183,FILE=''//tpdn//'/coriolis183',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(184,FILE=''//tpdn//'/coriolis184',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(191,FILE=''//tpdn//'/coriolis191',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(192,FILE=''//tpdn//'/coriolis192',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(193,FILE=''//tpdn//'/coriolis193',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(194,FILE=''//tpdn//'/coriolis194',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(195,FILE=''//tpdn//'/intense195',
     1FORM='UNFORMATTED',status='UNKNOWN')
      OPEN(200,FILE=''//tpdn//'/intense200',
     1FORM='UNFORMATTED',status='UNKNOWN')
      LTRAN=MAXSIZ
C input and output file names
       call getarg(1,infile)
       call getarg(2,outfile)

        if(infile.eq.'') stop 'blank input file name '
        if(outfile.eq.'') stop 'blank output file name '

          inquire(file=infile,exist=file_exists)

            if(.not.file_exists) then
               print*,'INPUT FILE ',infile,' DOES NOT EXIST '
               stop
            end if
      INP=1
      IOUT=2
           OPEN(INP,FILE=infile,status='OLD')
           OPEN(IOUT,FILE=outfile,status='unknown')

      MINP=3
      MOUT=4
      MOUTIN=7
      INP4=99
      INP5=100
      INP6=106
      INP7=107
      INP8=108
      INP9=109
      INP10=110
      CALL TIMIT(2)
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
 66   FORMAT(' USED MEMORY:',I10,3X,'OF',I10)
      CALL TIMIT(4)
      CALL TIMIT(5)
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
      END
C***************************************************
C***************************************************
      SUBROUTINE MEMPRT
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MSEG=200)
      COMMON/FILASS/IOUT,INP,MOUTIN
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      WRITE(IOUT,*)
      WRITE(IOUT,*)'MEMO ARRAY NO., PLUS CORRESPONDING ADDRESS & SIZE'
      WRITE(IOUT,*)
      DO I=1,MSEG
        WRITE(IOUT,*)I,LADD(I),KADD(I)
      END DO
      WRITE(IOUT,*)'FREE SPACE = ',MIN0(KFREE,KINF)
      WRITE(IOUT,*)
      CONTINUE
      RETURN
      END

