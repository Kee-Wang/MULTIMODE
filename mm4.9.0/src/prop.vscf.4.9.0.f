C****************************************************************
C****************************************************************
C**PROPERTIES
C****************************************************************
C****************************************************************
      SUBROUTINE RESUME(W,LDUMP)
      PARAMETER (MSEG=200)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FILASS/IOUT,INP
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
**ORIGINAL SETTINGS
C     COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
C    1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
C    2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
C    3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
C    4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
C    5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
C    6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
C    7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
C    8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
C    9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
C    1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
C    2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
C    3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
C    4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
C    5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
C    6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
C    7LDP4A,LDP4B,LDP4C,LDP4D
**ORIGINAL SETTINGS
      COMMON/TRANSF/LTRAN
C**ONLY KEEP THOSE SET UP PRIOR TO CALL
C**************************************************ASSIGN TOTAL STORAGE
C**LXK
      LADD(7)=0
      KADD(7)=0
C**LWK,LEVAL
      DO N=12,13
        LADD(N)=0
        KADD(N)=0
      END DO
C**LESCF TO LQM
      DO 123 N=26,97
      KADD(N)=0
      LADD(N)=0
123   CONTINUE
C**LNDUMP
      KADD(99)=0
      LADD(99)=0
      NADD=MSEG
C**LC0 TO LIP6
      DO 124 N=102,NADD
      KADD(N)=0
      LADD(N)=0
124   CONTINUE
      CALL PROP(W,LDUMP)
      WRITE(IOUT,66) LTRAN-KINF,LTRAN
 66   FORMAT('USED MEMORY:',I10,3X,'OF',I10)
      CALL TIMIT(4)
      STOP 'END OF MULTIMODE'
      END
C****************************************************************
C****************************************************************
      SUBROUTINE PROP(W,LDUMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR,TRIAT
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
C**CURRENT SETTINGS
      COMMON/CADDR/
     &LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LOLDH,LXK,LOLDXQ,LXW,LNBF,
     1LMBF,LSCF,LSX,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LMVB,LH,LXQ,LTEMP,LCONTR,
     3LDUM2(10),
     4LDUM3(10),
     5LISIZE,LW21,LSS,LSSX,LX21,LE21,LIDUMP,LNDUMP,LIP,LVCI,
     6LDUM6(10),
     7LDUM7(10),
     8LXA,LXA1,LXA2,LXA3,LXA4,LEVAL,LIP1,LIP2,LIP3,LIP4,
     9LPD1,LPD2,LPD3,LPD4,LPD5,LPD6,LPD7,LMXBAS,LKDUMP,LNVF,
     1LJP1,LJP2,LJP3,LJP4,LV1,LV2,LV3,LV4,LPD8,LPD9,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LDUM10(2)
C**CURRENT SETTINGS
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/DIPNO/NUMDIP
C*****
      COMMON/TITLE/TITLE
      COMMON/HERM/IHERM
      COMMON/CHECK/MCHECK
      COMMON/SADDLE/JNORM
      COMMON/PATH/ISCFCI
      COMMON/CYCLE/ICYCLE
      COMMON/DUMP/JDUMP,IDUMP,KDUMP,MDUMP,LLDUMP
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/CURVE/ICURV,ICORR
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/REACTL/JREACT
      COMMON/ENTER/IENTER,IENTMX(5),NTOT1,NTOT2,NTOT3,NTOT4,NTOT5
      COMMON/SINCOS/ICS
      COMMON/AVCON/AVCON
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FITTER/MFIT,MFIT1(4),XFIT1,MFIT2(4),XFIT2,MFIT3(4),XFIT3,
     1MFIT4(4),XFIT4
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TRIATO/TRIAT
      COMMON/LANTOL/TOLLAN
      COMMON/CYCLES/NCYCLE
      COMMON/MAXLAN/LANMAX,LLAN20,INP20
      COMMON/GIVEN/LGIV,IGIV
      COMMON/MATSIZ/MATSIZ
      COMMON/TYPE/LINEAR
      COMMON/ECKIND/INDECK
      COMMON/ROTIND/ANGSAV(3),INDROT
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/TOLS/TOL,EPS
      COMMON/ECKCNT/ICNT,INTC
      COMMON/ABINIT/ABINIT
      COMMON/EVL/EVL,CUT
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/NCPOT/NPOT
      COMMON/FACTOR/FACTOR(5),FACTS(5)
      COMMON/VCIMAX/NMAX
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/RPHROT/IROTV
      COMMON/ESTATE/IORDER
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/AXES/MX(3)
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/BASIS/NBAS(5,2),MAXSUM(5,2)
      COMMON/TBASIS/NTBAS(5,2),NTAU(5)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/SYMMP/ISYMP(10,10)
      COMMON/CONTDP/ICONDP
      COMMON/CSAVES/NNMODS,NAMODS,NVMODS,ICOUPS,JREACS,NREACS
      COMMON/NCREC/NREC(4),MAXBUF(4),NUNITR,NUNITW
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CSIZES/ISIZM1,ISIZM2,NVAL1,NVAL2,ICSIZ1,ICSIZ2,
     1IPSIZ1,IPSIZ2
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100)
      COMMON/CONTZ/NONC1,NONC2
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/UNITEX/I75,I76
C************************
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/MAXPT/MBFMAX,MBFMX1,MBFMX2,MBFMX3,MBFMX4,MBFMIN
      COMMON/DISC/IDISC
      COMMON/MODES/NMODE,NATOM
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/DISCSZ/KLC0,KEJK0,KLV1,KLC1,KEJK1,KLV2,KLC2,KEJK2,
     1KLV3,KLC3,KEJK3,KLV4,KLC4,KEJK4
      COMMON/SIZES/KTEMP,ISIZE1,ISIZE2,ISIZE3,ISIZE4,ISIZE5,ISIZE6,
     1ISIZE,JSIZE,ISIZMX
      COMMON/KJSZS/KJP1,KJP2,KJP3,KJP4,KJP5,KJP6
      COMMON/KPSZS/KIP1,KIP2,KIP3,KIP4,KIP5,KIP6,JCC1,JCC2,JCC3,JCC4,
     1JCC5
      COMMON/IPSZS/KPPP1,KPPP2,KPPP3,KPPP4,KPPP5
      COMMON/SIZEJ/JSIZE1(2),JSIZE2(2),JSIZE3(2),JSIZE4(2),JSIZE5(2)
      COMMON/TOTALS/ITOT1(2),ITOT2(2),ITOT3(2),ITOT4(2),ITOT5(2),ITOT
      COMMON/TOTK/KTOT(5,2)
      COMMON/MATRIX/NVAL,NVALR,KSTEP,KSIGN,NVALCF
C**************************************************************
C**************************************************************
      NVAL=NVALV
      NUMDIP=3
      IF(TRIAT)NUMDIP=2
C**FIRST SWEEP OF FILE (60)
      REWIND 60
C**************************************************************
C**************************************************************
C**NUMBER OF MODES (INPUT)
      READ(60)MMODE
      IF(MMODE.NE.NMODE)STOP 'WRONG MOLECULE'
      CALL MEMO(1,LMVB,NMODE,0,0,0,0,0,0,0,0)
      MAXCON=0
      MAXPRM=0
      MAXPTS=0
      DO MODE=1,NMODE
C**NUMBER CONTRACTED AND PRIMITIVE FUNCTIONS (INPUT)
        READ(60)NVV,NBB
        IF(NVV.GT.MAXCON)MAXCON=NVV
        CALL PUTINT(W(LNVF),MODE,NVV)
CC      NVF(MODE)=NVV
        IF(NBB.GT.MAXPRM)MAXPRM=NBB
        CALL PUTINT(W(LNBF),MODE,NBB)
CC      NBF(MODE)=NBB
        CALL MEMO(1,LTEMP,NBB,0,0,0,0,0,0,0,0)
        DO I=1,NVV
C**CONTRACTION COEFFICIENTS FROM PRIMITIVES (INPUT)
          CALL IN60(NBB,W(LTEMP))
CCCC      DO IY=1,NBB
CCCC        CONTR(IY,I,MODE)=W(IY)
CCCC      END DO
        END DO
        CALL MEMO(-1,LTEMP,NBB,0,0,0,0,0,0,0,0)
C**OMEGA AND NUMBER INTEGRATION POINTS (INPUT)
        READ(60)YLAM,MBB
        W(LOMEGA+MODE-1)=1/(YLAM*YLAM)
CC      OMEGA(MODE)=1/(YLAM*YLAM)
        IF(MBB.GT.MAXPTS)MAXPTS=MBB
        CALL PUTINT(W(LMVB),MODE,MBB)
CC      MVB(MODE)=MBB
        CALL MEMO(1,LTEMP,MBB,0,0,0,0,0,0,0,0)
C**INTEGRATION POINTS (INPUT)
        CALL INP60R(MBB,W(LTEMP))
CCCC    READ(60)(XQ(M,MODE),M=1,MBB)
        CALL MEMO(-1,LTEMP,MBB,0,0,0,0,0,0,0,0)
      END DO
      CALL MEMO(5,LCONTR,MAXPRM*MAXCON*NMODE,LXQ,MAXPTS*NMODE,LSCF,
     1MAXCON*MAXCON*NMODE,LH,MAXCON*MAXPTS*3*NMODE,LSX,MAXCON*MAXCON*
     2NMODE)
C  XQ is now the normal mode coordinate at point m mode mode
C**********************************
      DO MODE=1,NMODE
        CALL GETINT(W(LNVF),MODE,NVV)
CC      NVV=NVF(MODE)
        CALL GETINT(W(LMVB),MODE,MBB)
CC      MBB=MVB(MODE)
        DO I=1,NVV
          CALL MEMO(1,LTEMP,NVV,0,0,0,0,0,0,0,0)
C**SCF COEFFICIENTS FROM CONTRACTED FUNCTIONS (INPUT)
c  SCF is the square matrix of coeffs relating the 1d contracted basis
c  to the VSCF/CI basis.  (VSCF for zpe.)  i.e., for Ith vscf/VCI function
c  sum is over IY contracted 1d functions.
          CALL IN60(NVV,W(LTEMP))
          CALL IN3R(W(LSCF),MAXCON,MAXCON,I,MODE,W(LTEMP),NVV)
CC        DO IY=1,NVV
CC          SCF(IY,I,MODE)=W(IY)
CC        END DO
          CALL MEMO(-1,LTEMP,NVV,0,0,0,0,0,0,0,0)
        END DO
        DO K=1,3
          DO M=1,MBB
            CALL MEMO(1,LTEMP,NVV,0,0,0,0,0,0,0,0)
C**VSCF/VCI FUNCTIONS, 1ST AND 2ND DERIVATIVES AT INTEGRATION POINTS
C**WITH SQRT(WEIGHT) INCLUDED (INPUT)
            CALL IN60(NVV,W(LTEMP))
            CALL IN4(W(LH),MAXCON,MAXPTS,M,K,MODE,W(LTEMP),NVV)
CC          DO IY=1,NVV
CC            H(IY,M,K,MODE)=W(IY)
CC          END DO
            CALL MEMO(-1,LTEMP,NVV,0,0,0,0,0,0,0,0)
          END DO
        END DO
      END DO
C**********************************
C**NUMBER VCI SYMMETRIES (INPUT)
      READ(60)NVSYM,ICI,NMAX,(NBAS(I,1),I=1,5),(MAXSUM(I,1),I=1,5)
      CALL MEMO(3,LISIZE,NVSYM,LIDUMP,NVSYM,LNDUMP,NVSYM,0,0,0,0)
      MAXCI=0
      MAXDMP=0
      DO IVSYM=1,NVSYM
C**SIZE OF VCI MATRIX (INPUT)
        READ(60)JSIZE
        NTOT(IVSYM)=JSIZE
        IF(JSIZE.GT.MAXCI)MAXCI=JSIZE
        CALL PUTINT(W(LISIZE),IVSYM,JSIZE)
CC      ISIZE(IVSYM)=JSIZE
        DO I=1,JSIZE
          CALL MEMO(1,LTEMP,NMODE,0,0,0,0,0,0,0,0)
C**VCI BASIS (INPUT)
          CALL INP60I(NMODE,W(LTEMP))
c  IP(I,J,IVSYM): I refers to the eigenvector index, J is the mode
c  and IP is the corrsponding basis function no. for the Jth mode.
CCCC      READ(60)(IP(I,J,IVSYM),J=1,NMODE)
          CALL MEMO(-1,LTEMP,NMODE,0,0,0,0,0,0,0,0)
        END DO
C**NUMBER of  CI FUNCTIONS (INPUT)
        READ(60)JDUMP
        IF(JDUMP.GT.MAXDMP)MAXDMP=JDUMP
        CALL PUTINT(W(LIDUMP),IVSYM,JDUMP)
CC      IDUMP(IVSYM)=JDUMP
        CALL PUTINT(W(LNDUMP),IVSYM,0)
CC      NDUMP(IVSYM)=0
        DO J=1,JDUMP
C**SPECIFIC CI FUNCTION AND ENERGY (INPUT)
          READ(60)K
CCCC      KDUMP(J,IVSYM)=K
          IF(K.NE.0)THEN
            READ(60)ENERGY
CCCC        EVAL(J,IVSYM)=ENERGY
            CALL GETINT(W(LNDUMP),IVSYM,II)
            II=II+1
            CALL PUTINT(W(LNDUMP),IVSYM,II)
CC          NDUMP(IVSYM)=NDUMP(IVSYM)+1
            CALL MEMO(1,LTEMP,JSIZE,0,0,0,0,0,0,0,0)
C**CI COEFFICIENTS FOR Jth EIGENFUNCTION (INPUT)
            CALL INP60R(JSIZE,W(LTEMP))
CCCC        READ(60)(VCI(I,J,IVSYM),I=1,JSIZE)
            CALL MEMO(-1,LTEMP,JSIZE,0,0,0,0,0,0,0,0)
          END IF
        END DO
      END DO
      CALL MEMO(5,LIP,MAXCI*NMODE*NVSYM,LVCI,MAXCI*MAXDMP*NVSYM,
     1LKDUMP,MAXDMP*NVSYM,LXK,MAXDMP*MAXDMP,LEVAL,MAXDMP*NVSYM)
C**************************************************************
C**************************************************************
C**SECOND SWEEP OF FILE (60)
      REWIND 60
C**************************************************************
C**************************************************************
C**NUMBER OF MODES (INPUT)
      READ(60)NMODE
      DO MODE=1,NMODE
C**NUMBER CONTRACTED AND PRIMITIVE FUNCTIONS (INPUT)
        READ(60)NVV,NBB
        CALL MEMO(1,LTEMP,NBB,0,0,0,0,0,0,0,0)
        DO I=1,NVV
C**CONTRACTION COEFFICIENTS FROM PRIMITIVES (INPUT)
c  contract ho primitives to eigs of 1d cuts in normal coords.  Contract
c  from NBB HO functions to NVV contracted 1d functions.  CONTR is the
c  trans matrix from the HOs to the contracted 1d functions.
          CALL IN60(NBB,W(LTEMP))
          CALL IN3R(W(LCONTR),MAXPRM,MAXCON,I,MODE,W(LTEMP),NBB)
CC        DO IY=1,NBB
CC          CONTR(IY,I,MODE)=W(IY)
CC        END DO
        END DO
        CALL MEMO(-1,LTEMP,NBB,0,0,0,0,0,0,0,0)
C**OMEGA AND NUMBER INTEGRATION POINTS (INPUT)
        READ(60)YLAM,MBB
        CALL MEMO(1,LTEMP,MBB,0,0,0,0,0,0,0,0)
C**INTEGRATION POINTS (INPUT)
        CALL INP60R(MBB,W(LTEMP))
        CALL IN2R(W(LXQ),MAXPTS,MODE,W(LTEMP),MBB)
CC      READ(60)(XQ(M,MODE),M=1,MBB)
        CALL MEMO(-1,LTEMP,MBB,0,0,0,0,0,0,0,0)
      END DO
C  XQ is now the normal mode coordinate at point m mode mode
C**********************************
      DO MODE=1,NMODE
        CALL GETINT(W(LNVF),MODE,NVV)
CC      NVV=NVF(MODE)
        CALL GETINT(W(LMVB),MODE,MBB)
CC      MBB=MVB(MODE)
        DO I=1,NVV
          CALL MEMO(1,LTEMP,NVV,0,0,0,0,0,0,0,0)
C**SCF COEFFICIENTS FROM CONTRACTED FUNCTIONS (INPUT)
          CALL IN60(NVV,W(LTEMP))
          CALL MEMO(-1,LTEMP,NVV,0,0,0,0,0,0,0,0)
        END DO
        DO K=1,3
          DO M=1,MBB
            CALL MEMO(1,LTEMP,NVV,0,0,0,0,0,0,0,0)
C**VSCF/VCI FUNCTIONS, 1ST AND 2ND DERIVATIVES AT INTEGRATION POINTS
**WITH SQRT(WEIGHT) INCLUDED (INPUT)
            CALL IN60(NVV,W(LTEMP))
            CALL MEMO(-1,LTEMP,NVV,0,0,0,0,0,0,0,0)
          END DO
        END DO
      END DO
C**********************************
C**NUMBER VCI SYMMETRIES (INPUT)
      READ(60)NVSYM,ICI,NMAX,(NBAS(I,1),I=1,5),(MAXSUM(I,1),I=1,5)
      DO IVSYM=1,NVSYM
C**SIZE OF VCI MATRIX (INPUT)
        READ(60)JSIZE
        DO I=1,JSIZE
          CALL MEMO(1,LTEMP,NMODE,0,0,0,0,0,0,0,0)
C**VCI BASIS (INPUT)
          CALL INP60I(NMODE,W(LTEMP))
c  IP(I,J,IVSYM): I refers to the eigenvector index, J is the mode
c  and IP is the corrsponding basis function no. for the Jth mode.
          CALL IN3I(W(LIP),MAXCI,NMODE,I,IVSYM,W(LTEMP))
CC        READ(60)(IP(I,J,IVSYM),J=1,NMODE)
          CALL MEMO(-1,LTEMP,NMODE,0,0,0,0,0,0,0,0)
        END DO
C**NUMBER of  CI FUNCTIONS (INPUT)
        READ(60)JDUMP
        DO J=1,JDUMP
C**SPECIFIC CI FUNCTION AND ENERGY (INPUT)
          READ(60)K
          CALL IN2I(W(LKDUMP),MAXDMP,J,IVSYM,K)
CC        KDUMP(J,IVSYM)=K
          IF(K.NE.0)THEN
            READ(60)ENERGY
            CALL IN2D(W(LEVAL),MAXDMP,J,IVSYM,ENERGY)
CC          EVAL(J,IVSYM)=ENERGY
            CALL MEMO(1,LTEMP,JSIZE,0,0,0,0,0,0,0,0)
C**CI COEFFICIENTS FOR Jth EIGENFUNCTION (INPUT)
            CALL INP60R(JSIZE,W(LTEMP))
            CALL IN3R(W(LVCI),MAXCI,MAXDMP,J,IVSYM,W(LTEMP),JSIZE)
CC          READ(60)(VCI(I,J,IVSYM),I=1,JSIZE)
            CALL MEMO(-1,LTEMP,JSIZE,0,0,0,0,0,0,0,0)
          END IF
        END DO
      END DO
C**********************************************************************
C**********************************************************************
C**DUMP DIPOLE GRIDS
C**********************************************************************
C**********************************************************************
      IF(LDUMP.EQ.1)
     1CALL GRIDS(W(1),W(LXQ),MAXPTS,W(LMVB),NMODE,NATOM)
C**IF LDUMP<0:
C**GRIDS SET UP WITH (MBF) GAUSS POINTS GIVEN BY (OLDXQ)
C**WILL HAVE TO GENERATE PRIMITIVES, CONTRACTED FNS. (OLDH)
C**AND VSCF FUNCTIONS AT NEW (MBF) GAUSS POINTS.
      LDUMP=IABS(LDUMP)

C**********************************************************************
C**TEMPORARY TEST OVERLAPS
C**********************************************************************
CCCC  CALL PROV(W(LNVF),W(LMVB),NMODE,W(LSX),MAXCON,W(LH),MAXPTS)
CCCC  CALL CIOV(W(LIDUMP),W(LNDUMP),W(LISIZE),NVSYM,W(LXK),W(LKDUMP),
CCCC 1MAXDMP,W(LVCI),W(LIP),MAXCI,NMODE)
CCCC  STOP 'TEST OVERLAP'
C**********************************************************************
C**TEMPORARY TEST OVERLAPS
C**********************************************************************

      IF(LDUMP.NE.1)GO TO 9999
C**********************************************************************
C**********************************************************************
C**GET BASIS SETS FOR BASIC INTEGRALS
C**********************************************************************
C**********************************************************************
4000  CONTINUE
C**RESET NMAX IF UNRESTRICTED
      IF(NMAX.EQ.0)NMAX=-ICI*NMODE
      NBFMIN=100000
      DO K=1,NMODE
        CALL GETINT(W(LNVF),K,IK3)
        IF(IK3.LT.NBFMIN)NBFMIN=IK3
      END DO
C**********************************************************************
C**********************************************************************
C**                                                                 VCI
C**********************************************************************
C**********************************************************************
      IF(ICI.LT.0.AND.NMAX.GE.0)THEN
        JCI=-ICI
        IF(NMAX.GT.0)THEN
C********************************************************************
C**CONVERT NMAX FROM NUMBER OF FUNCTIONS TO ORIGINAL NUMBER OF QUANTA
          NMAX=NMAX-NMODE
C********************************************************************
C**ONLY GET NMAX +1 (MAXIMUM) INTEGRALS
          IF(NMAX+1.LT.JCI)JCI=NMAX+1
C**OR NBFMIN (MAXIMUM) INTEGRALS (WHICHEVER IS SMALLER)
          IF(JCI.GT.NBFMIN)JCI=NBFMIN
C**CONVERT NMAX FROM QUANTA TO NUMBER OF FUNCTIONS
          NMAX=NMAX+NMODE
        END IF
      END IF
      IF(ICI.LT.0)THEN
C********************************************************************
C**RESET MAXJ
        IF(NMAX.LT.0)THEN
          J=-NMAX
          IF(J.GT.5)J=5
          IF(J.LT.0)J=0
          MAXJ=J
        END IF
C**SAVE NTOT
        DO I=1,NVSYM
          NSYM(I)=NTOT(I)
        END DO
C********************************************************************
        NVSYMX=NVSYM
        NVSYM=1
C**NNMAX USED AS INDICATOR AS TO WHICH ALGORITHM
        NNMAX=NMAX
        IF(NNMAX.GE.0)THEN
C********************************************************************
C**NMAX SAVED IN NMAXMX (RE-LOADED LATER)
          NMAXMX=NMAX
C********************************************************************
          DO I=1,5
C**MAX POSSIBLE QUANTA
            NBAS(I,1)=JCI-1
C**MAX POSSIBLE SUM
            MAXSUM(I,1)=I*NBAS(I,1)
            IF(NNMAX.GT.0)MAXSUM(I,1)=MIN0(NNMAX,I*NBAS(I,1))
          END DO
C**INDIVIDUAL QUANTUM
          CALL DEFMAX(W(LMXBAS),NMODE,JCI)
          JCI1=JCI
          JCI2=JCI
          JCI3=JCI
          JCI4=JCI
          JCI5=JCI
        ELSE
C**MAX POSSIBLE QUANTA
          IF(IREACT.NE.0)CALL TMAXIN(W(LMXBAS),NSMODE,MAXJ)
          JCI1=0
          JCI2=0
          JCI3=0
          JCI4=0
          JCI5=0
          DO K=1,NCONT
            ICI1=NBAS(1,K)+1
            ICI2=MAX0(NBAS(1,K)+1,NBAS(2,K)+1)
            ICI3=MAX0(NBAS(1,K)+1,NBAS(2,K)+1,NBAS(3,K)+1)
            ICI4=MAX0(NBAS(1,K)+1,NBAS(2,K)+1,NBAS(3,K)+1,NBAS(4,K)+1)
            ICI5=MAX0(NBAS(1,K)+1,NBAS(2,K)+1,NBAS(3,K)+1,NBAS(4,K)+1,
     1      NBAS(5,K)+1)
            IF(ICI1.GT.JCI1)JCI1=ICI1
            IF(ICI2.GT.JCI2)JCI2=ICI2
            IF(ICI3.GT.JCI3)JCI3=ICI3
            IF(ICI4.GT.JCI4)JCI4=ICI4
            IF(ICI5.GT.JCI5)JCI5=ICI5
          END DO
        END IF
        JCC1=JCI1
        JCC2=JCI2
        JCC3=JCI3
        JCC4=JCI4
        JCC5=JCI5
C*******************************************************
C**GET SIZES OF 1-,2-,3-,4-,5- MODE BASIC BASES - START
C*******************************************************
        NNMODS=NNMODE
        IF(ICOUPL.GT.0.OR.(NNMAX.LT.0.AND.ICI.LT.0))THEN
          NMAX1=MAXSUM(1,1)+1
          NMAX2=MAXSUM(1,2)+1
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE1=JCI1
C**FIRST GET SIZE (ISIZE1 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE1,1,JCI1,0)
          NSIZE1=ISIZE1
          KIP1=ISIZE1
          CALL MEMO(1,LIP1,KIP1,0,0,0,0,0,0,0,0)
        END IF
C
        IF(ICOUPL.GT.1.OR.(NNMAX.LT.0.AND.ICI.LT.-1))THEN
          NMAX1=MAX0(MAXSUM(1,1)+2,MAXSUM(2,1)+2)
          NMAX2=MAX0(MAXSUM(1,2)+2,MAXSUM(2,2)+2)
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE2=JCI2**2
C**FIRST GET SIZE (ISIZE2 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE2,2,JCI2,0)
          NSIZE2=ISIZE2
          KIP2=ISIZE2*2
          CALL MEMO(1,LIP2,KIP2,0,0,0,0,0,0,0,0)
        END IF
C
        IF(ICOUPL.GT.2.OR.(NNMAX.LT.0.AND.ICI.LT.-2))THEN
          NMAX1=MAX0(MAXSUM(1,1)+3,MAXSUM(2,1)+3,MAXSUM(3,1)+3)
          NMAX2=MAX0(MAXSUM(1,2)+3,MAXSUM(2,2)+3,MAXSUM(3,2)+3)
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE3=JCI3**3
C**FIRST GET SIZE (ISIZE3 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE3,3,JCI3,0)
          NSIZE3=ISIZE3
          KIP3=ISIZE3*3
          CALL MEMO(1,LIP3,KIP3,0,0,0,0,0,0,0,0)
        END IF
C  
        IF(ICOUPL.GT.3.OR.(NNMAX.LT.0.AND.ICI.LT.-3))THEN
          NMAX1=MAX0(MAXSUM(1,1)+4,MAXSUM(2,1)+4,MAXSUM(3,1)+4,
     1    MAXSUM(4,1)+4)
          NMAX2=MAX0(MAXSUM(1,2)+4,MAXSUM(2,2)+4,MAXSUM(3,2)+4,
     1    MAXSUM(4,2)+4)
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE4=JCI4**4
C**FIRST GET SIZE (ISIZE4 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE4,4,JCI4,0)
          NSIZE4=ISIZE4
          KIP4=ISIZE4*4
          CALL MEMO(1,LIP4,KIP4,0,0,0,0,0,0,0,0)
        END IF
C
        IF(ICOUPL.GT.4.OR.(NNMAX.LT.0.AND.ICI.LT.-4))THEN
          NMAX1=MAX0(MAXSUM(1,1)+5,MAXSUM(2,1)+5,MAXSUM(3,1)+5,
     1    MAXSUM(4,1)+5,MAXSUM(5,1)+5)
          NMAX2=MAX0(MAXSUM(1,2)+5,MAXSUM(2,2)+5,MAXSUM(3,2)+5,
     1    MAXSUM(4,2)+5,MAXSUM(5,2)+5)
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE5=JCI5**5
C**FIRST GET SIZE (ISIZE5 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE5,5,JCI5,0)
          NSIZE5=ISIZE5
          KIP5=ISIZE5*5
          CALL MEMO(1,LIP5,KIP5,0,0,0,0,0,0,0,0)
        END IF
        NNMODE=NNMODS
C*******************************************************
C**GET SIZES OF 1-,2-,3-,4-,5- MODE BASIC BASES - END
C*******************************************************
        NVSYM=NVSYMX
C**RESTORE NTOT
        DO I=1,NVSYM
          NTOT(I)=NSYM(I)
        END DO
      END IF
C**********************************************************************
C**********************************************************************
C**GET BASIC INTEGRALS FOR VIBRATIONAL BASIS
C**********************************************************************
C**********************************************************************
4500  CONTINUE
      ITIM=-1
      IF(ICOUPL.GE.0)THEN
        REWIND 21
        REWIND 31
        REWIND 41
      END IF
      IF(ICOUPL.GT.1)THEN
        REWIND 22
        REWIND 32
        REWIND 42
      END IF
      IF(ICOUPL.GT.2)THEN
        REWIND 23
        REWIND 33
        REWIND 43
      END IF
      IF(ICOUPL.GT.3)THEN
        REWIND 24
        REWIND 34
        REWIND 44
      END IF
C****************************
C**A1,B2,B1 COMPONENTS OF DIPOLE
C****************************
      DO 3300 IDIP=1,NUMDIP
      ITIM=ITIM+1
      IF(ICOUPL.GT.0)THEN
        REWIND 61
        REWIND 71
        REWIND 81
      END IF
      IF(ICOUPL.GT.1)THEN
        REWIND 62
        REWIND 72
        REWIND 82
      END IF
      IF(ICOUPL.GT.2)THEN
        REWIND 63
        REWIND 73
        REWIND 83
      END IF
      IF(ICOUPL.GT.3)THEN
        REWIND 64
        REWIND 74
        REWIND 84
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NMODE
        IF(K.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
C       CALL MEMO(1,LIP1,KIP1,0,0,0,0,0,0,0,0)
        CALL GETBP1(W(LIP1),ISIZE1,NSMODE,K,W(LMXBAS),NSIZE1,0)
        KXA1=NSIZE1*(NSIZE1+1)/2
        CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGZ(NSIZE1,NSIZE1,W(LXA1),NSIZE1,W(LXA1),ICI)
        CALL GETINT(W(LMVB),K,IK2)
        CALL GETINT(W(LNVF),K,IK3)
        CALL V0DP1(NMODE,1,MAXCON,MAXPTS,W(LH+K2),W(LXQ+K3),W(LXA1),
     1  NSIZE1,IK3,IK2,W(LIP1),ISIZE1,W(LV1),W(LV1),IDIP)
        CALL MATOUT(W(LXA1),W(LXA1),NSIZE1,10*IDIP+11)
        CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
        IF(ICOUPL.EQ.1)GO TO 3301
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
C         CALL MEMO(1,LIP2,KIP2,0,0,0,0,0,0,0,0)
          CALL GETBP2(W(LIP2),ISIZE2,NSMODE,K,L,W(LMXBAS),NSIZE2,0)
          KXA2=NSIZE2*(NSIZE2+1)/2
          JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
          JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
          KTEMP=JCI2*JCI2
          CALL MEMO(2,LXA2,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
C**ZEROISE MATRIX
          CALL DIAGZ(NSIZE2,NSIZE2,W(LXA2),NSIZE2,W(LXA2),ICI)
          CALL GETINT(W(LMVB),L,IL2)
          CALL GETINT(W(LNVF),L,IL3)
          CALL V0DP2(NMODE,1,2,MAXCON,MAXPTS,W(LH+K2),W(LXQ+K3),
     1    W(LH+L2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXA2),NSIZE2,W(LIP2),
     2    ISIZE2,W(LTEMP),JCI1,JCI2,W(LV2),W(LV2),W(LIP1),ISIZE1,
CC   3    IDIP)
     3    NSIZE1,IDIP)
          CALL MATOUT(W(LXA2),W(LXA2),NSIZE2,10*IDIP+12)
          CALL MEMO(-2,LXA2,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
          IF(ICOUPL.EQ.2)GO TO 3302
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
C           CALL MEMO(1,LIP3,KIP3,0,0,0,0,0,0,0,0)
            CALL GETBP3(W(LIP3),ISIZE3,NSMODE,K,L,N,W(LMXBAS),NSIZE3,0)
            KXA3=NSIZE3*(NSIZE3+1)/2
            JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
            JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
            JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
            KTEMP=JCI2*JCI3*JCI2*JCI3
            CALL MEMO(2,LXA3,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE3,NSIZE3,W(LXA3),NSIZE3,W(LXA3),ICI)
            CALL GETINT(W(LMVB),N,IN2)
            CALL GETINT(W(LNVF),N,IN3)
            CALL V0DP3(NMODE,1,2,3,MAXCON,MAXPTS,W(LH+K2),W(LXQ+K3),
     1      W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,
     2      IN2,W(LXA3),NSIZE3,W(LIP3),ISIZE3,W(LTEMP),JCI1,JCI2,JCI3,
CC   3      W(LV3),W(LV3),W(LIP2),ISIZE2,IDIP)
     3      W(LV3),W(LV3),W(LIP2),ISIZE2,NSIZE2,IDIP)
            CALL MATOUT(W(LXA3),W(LXA3),NSIZE3,10*IDIP+13)
            CALL MEMO(-2,LXA3,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
            IF(ICOUPL.EQ.3)GO TO 3303
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.K.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
C             CALL MEMO(1,LIP4,KIP4,0,0,0,0,0,0,0,0)
              CALL GETBP4(W(LIP4),ISIZE4,NSMODE,K,L,N,M,W(LMXBAS),
     2        NSIZE4,0)
              KXA4=NSIZE4*(NSIZE4+1)/2
              JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
              JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
              JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
              JCI4=MAXBFN(W(LMXBAS),NMODE,M,1)
              KTEMP=JCI2*JCI3*JCI4*JCI2*JCI3*JCI4
              CALL MEMO(2,LXA4,KXA4,LTEMP,KTEMP,0,0,0,0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(ISIZE4,ISIZE4,W(LXA4),ISIZE4,W(LXA4),ICI)
              CALL GETINT(W(LMVB),M,IM2)
              CALL GETINT(W(LNVF),M,IM3)
              CALL V0DP4(NMODE,1,2,3,4,MAXCON,MAXPTS,W(LH+K2),
     1        W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2        W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA4),NSIZE4,
     4        W(LIP4),ISIZE4,W(LTEMP),JCI1,JCI2,JCI3,JCI4,W(LV4),W(LV4),
CC   6        W(LIP3),ISIZE3,IDIP)
     6        W(LIP3),ISIZE3,NSIZE3,IDIP)
              CALL MATOUT(W(LXA4),W(LXA4),ISIZE4,10*IDIP+14)
              CALL MEMO(-2,LXA4,KXA4,LTEMP,KTEMP,0,0,0,0,0,0)
              IF(ICOUPL.EQ.4)GO TO 3304
C**5-MODE COUPLING HERE IF NEEDED
3304  CONTINUE
              M2=M2+3*MAXCON*MAXPTS
              M3=M3+MAXPTS
            END DO
3303  CONTINUE
            N2=N2+3*MAXCON*MAXPTS
            N3=N3+MAXPTS
          END DO
3302  CONTINUE
          L2=L2+3*MAXCON*MAXPTS
          L3=L3+MAXPTS
        END DO
3301  CONTINUE
        K2=K2+3*MAXCON*MAXPTS
        K3=K3+MAXPTS
      END DO
3300  CONTINUE
      IF(ICOUPL.GT.0)THEN
        CALL MEMO(-1,LV1,KLV1,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.1)THEN
        CALL MEMO(-1,LV2,KLV2,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.2)THEN
        CALL MEMO(-1,LV3,KLV3,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.3)THEN
        CALL MEMO(-1,LV4,KLV4,0,0,0,0,0,0,0,0)
      END IF
      ITIM=-1
C**********************************************************************
C**********************************************************************
C**                                                SET UP DIPOLE MATRIX
C**********************************************************************
C**********************************************************************
      ISIZMX=MAXCI
C****************************
C**LOOP OVER SYMMETRIES (RHS)
C****************************
      DO 4400 MSR=1,NVSYM
      ITIM=ITIM+1
      ISIZER=NTOT(MSR)
C****************************
C**LOOP OVER SYMMETRIES (LHS)
C****************************
      DO 4400 MSL=1,MSR
      ISIZEL=NTOT(MSL)
      IF(ISIZER.EQ.0.OR.ISIZEL.EQ.0)GO TO 4401
      KXA=ISIZEL*ISIZER
      KTEMP=ISIZEL*MAXDMP
      CALL MEMO(2,LXA,KXA,LTEMP,KTEMP,0,0,0,0,0,0)
      IF(ICOUPL.GE.0)THEN
        REWIND 21
        REWIND 31
        REWIND 41
      END IF
      IF(ICOUPL.GT.1)THEN
        REWIND 22
        REWIND 32
        REWIND 42
      END IF
      IF(ICOUPL.GT.2)THEN
        REWIND 23
        REWIND 33
        REWIND 43
      END IF
      IF(ICOUPL.GT.3)THEN
        REWIND 24
        REWIND 34
        REWIND 44
      END IF
C****************************
C**A1,B2,B1 COMPONENTS OF DIPOLE
C****************************
      DO 9000 IDIP=1,NUMDIP
C*****************************************************************
C*****************************************************************
C**NEEDS SELECTION CRITERION HERE FOR SYMMETRY-ALLOWED DIPOLES
C*****************************************************************
C*****************************************************************
C**TEST IF C2V
      IF(NVSYM.EQ.4.OR.(NVSYM.EQ.2.AND.TRIAT))THEN
C**A1 COMPONENT
        IF(IDIP.EQ.1)THEN
          IF(MSL.NE.MSR)GO TO 9001
        END IF
C**B2 COMPONENT
        IF(IDIP.EQ.2)THEN
          IF(MSL.EQ.2)GO TO 9001
          IF(MSL.NE.MSR-1)GO TO 9001
        END IF
C**B1 COMPONENT
        IF(IDIP.EQ.3)THEN
          IF(MSL.NE.MSR-2)GO TO 9001
        END IF
      END IF
C**ZEROISE MATRIX
      CALL DIAGZ(ISIZEL,ISIZER,W(LXA),IDUM,W(LXA),0)
      CALL VDP0(NMODE,W(LXA),ISIZEL,ISIZER,IDIP,MSL,MSR)
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NMODE
        IF(K.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        CALL GETINT(W(LMVB),K,IK2)
        CALL GETINT(W(LNVF),K,IK3)
        CALL GETBP1(W(LIP1),ISIZE1,NSMODE,K,W(LMXBAS),NSIZE1,0)
C**READ INTO MATRIX
        KXA1=NSIZE1*(NSIZE1+1)/2
        CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
        CALL MATIN(W(LXA1),W(LXA1),NSIZE1,10*IDIP+11)
        CALL VDP1(NMODE,K,MAXCON,MAXPTS,W(LH+K2),W(LXQ+K3),W(LXA),
     1  W(LXA1),NSIZE1,IK3,IK2,W(LIP),ISIZMX,ISIZEL,ISIZER,MSL,MSR,
     2  W(LIP1),ISIZE1)
        CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
C       CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
        IF(ICOUPL.EQ.1)GO TO 701
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          CALL GETINT(W(LMVB),L,IL2)
          CALL GETINT(W(LNVF),L,IL3)
          CALL GETBP2(W(LIP2),ISIZE2,NSMODE,K,L,W(LMXBAS),NSIZE2,0)
C**READ INTO MATRIX
          KXA2=NSIZE2*(NSIZE2+1)/2
          CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
          CALL MATIN(W(LXA2),W(LXA2),NSIZE2,10*IDIP+12)
          CALL VDP2(NMODE,K,L,MAXCON,MAXPTS,W(LH+K2),W(LXQ+K3),
     1    W(LH+L2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXA),W(LXA2),NSIZE2,
     3    W(LIP),ISIZMX,ISIZEL,ISIZER,MSL,MSR,W(LIP2),ISIZE2)
          CALL MEMO(-1,LXA2,KXA2,0,0,0,0,0,0,0,0)
C         CALL MEMO(-1,LIP2,KIP2,0,0,0,0,0,0,0,0)
          IF(ICOUPL.EQ.2)GO TO 702
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            CALL GETINT(W(LMVB),N,IN2)
            CALL GETINT(W(LNVF),N,IN3)
            CALL GETBP3(W(LIP3),ISIZE3,NSMODE,K,L,N,W(LMXBAS),NSIZE3,0)
C**READ INTO MATRIX
            KXA3=NSIZE3*(NSIZE3+1)/2
            CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
            CALL MATIN(W(LXA3),W(LXA3),NSIZE3,10*IDIP+13)
            CALL VDP3(NMODE,K,L,N,MAXCON,MAXPTS,W(LH+K2),W(LXQ+K3),
     1      W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,
     2      IN2,W(LXA),W(LXA3),NSIZE3,W(LIP),ISIZMX,ISIZEL,ISIZER,MSL,
     3      MSR,W(LIP3),ISIZE3)
            CALL MEMO(-1,LXA3,KXA3,0,0,0,0,0,0,0,0)
C           CALL MEMO(-1,LIP3,KIP3,0,0,0,0,0,0,0,0)
            IF(ICOUPL.EQ.3)GO TO 703
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.K.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              CALL GETINT(W(LMVB),M,IM2)
              CALL GETINT(W(LNVF),M,IM3)
              CALL GETBP4(W(LIP4),ISIZE4,NSMODE,K,L,N,M,W(LMXBAS),
     2        NSIZE4,0)
C**READ INTO MATRIX
              KXA4=NSIZE4*(NSIZE4+1)/2
              CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
              CALL MATIN(W(LXA4),W(LXA4),NSIZE4,10*IDIP+14)
              CALL VDP4(NMODE,K,L,N,M,MAXCON,MAXPTS,W(LH+K2),W(LXQ+K3),
     1        W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     2        IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA),W(LXA4),NSIZE4,
     4        W(LIP),ISIZMX,ISIZEL,ISIZER,MSL,MSR,W(LIP4),ISIZE4)
              CALL MEMO(-1,LXA4,KXA4,0,0,0,0,0,0,0,0)
C             CALL MEMO(-1,LIP4,KIP4,0,0,0,0,0,0,0,0)
              IF(ICOUPL.EQ.4)GO TO 704
C**5-MODE COUPLING HERE IF NEEDED
704   CONTINUE
              M2=M2+3*MAXCON*MAXPTS
              M3=M3+MAXPTS
            END DO
703   CONTINUE
            N2=N2+3*MAXCON*MAXPTS
            N3=N3+MAXPTS
          END DO
702   CONTINUE
          L2=L2+3*MAXCON*MAXPTS
          L3=L3+MAXPTS
        END DO
701   CONTINUE
        K2=K2+3*MAXCON*MAXPTS
        K3=K3+MAXPTS
      END DO
C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
      WRITE(IOUT,*)
      CALL DIPTOT(W(LNDUMP),NVSYM,W(LXK),MAXDMP,W(LVCI),MAXCI,W(LXA),
     1ISIZEL,ISIZER,MSL,MSR,W(LTEMP),W(LEVAL),IDIP)
9001  CONTINUE
9000  CONTINUE
      CALL MEMO(-2,LXA,KXA,LTEMP,KTEMP,0,0,0,0,0,0)
4401  CONTINUE
4400  CONTINUE
C**************************************************************
9999  CONTINUE
C**************************************************************
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE IN60(NB,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(NB)
      READ(60)W
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE INP60R(MB,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(MB)
      READ(60)(W(M),M=1,MB)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE INP60I(MB,IW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IW(MB)
      READ(60)(IW(M),M=1,MB)
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE IN2R(XQ,MAXPTS,MODE,W,MBB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XQ(MAXPTS,1),W(1)
      DO M=1,MBB
        XQ(M,MODE)=W(M)
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE IN2D(EVAL,MAXDMP,J,ISYM,ENERGY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EVAL(MAXDMP,1)
      EVAL(J,ISYM)=ENERGY
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE IN2I(KDUMP,MAXDMP,J,ISYM,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION KDUMP(MAXDMP,1)
      KDUMP(J,ISYM)=K
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE IN3R(SCF,MAXCN1,MAXCN2,I,MODE,W,NVV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SCF(MAXCN1,MAXCN2,1),W(1)
      DO IY=1,NVV
        SCF(IY,I,MODE)=W(IY)
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE IN3I(IP,MAXCI,NMODE,I,ISYM,IW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IP(MAXCI,NMODE,1),IW(1)
      DO J=1,NMODE
        IP(I,J,ISYM)=IW(J)
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE IN4(H,MAXCON,MAXPTS,M,K,MODE,W,NVV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(MAXCON,MAXPTS,3,1),W(1)
      DO IY=1,NVV
        H(IY,M,K,MODE)=W(IY)
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE PUTINT(NV,MODE,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NV(1)
      NV(MODE)=I
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETINT(NV,MODE,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NV(1)
      I=NV(MODE)
      RETURN
      END
C****************************************************************
C****************************************************************
C**TEMPORARY TEST OVERLAPS
      SUBROUTINE PROV(NV,MVB,NMODE,SX,MAXCON,H,MAXPTS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NV(NMODE),MVB(NMODE),SX(MAXCON,MAXCON,NMODE)
      DIMENSION H(MAXCON,MAXPTS,3,NMODE)
      COMMON/FILASS/IOUT,INP
C**SCF OVERLAP MATRICES
      DO MODE=1,NMODE
        NVV=NV(MODE)
        DO IX=1,NVV
          DO IY=1,NVV
            SX(IY,IX,MODE)=0
          END DO
        END DO
        MBB=MVB(MODE)
        DO M=1,MBB
          DO IX=1,NVV
            DO IY=1,NVV
              OVER=H(IY,M,1,MODE)*H(IX,M,1,MODE)
              SX(IY,IX,MODE)=SX(IY,IX,MODE)+OVER
            END DO
          END DO
        END DO
        WRITE(IOUT,*)'MODE ',MODE
        DO IX=1,NVV
          WRITE(IOUT,*)(SX(IY,IX,MODE),IY=1,NVV)
          CALL FLUSH(IOUT)
        END DO
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE CIOV(IDUMP,NDUMP,ISIZE,NVSYM,XK,KDUMP,MAXDMP,VCI,
     1IP,MAXCI,NMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IDUMP(NVSYM),NDUMP(NVSYM),ISIZE(NVSYM)
      DIMENSION XK(MAXDMP,MAXDMP),KDUMP(MAXDMP,NVSYM)
      DIMENSION VCI(MAXCI,MAXDMP,NVSYM),IP(MAXCI,NMODE,NVSYM)
      COMMON/FILASS/IOUT,INP
C**TEST CI OVERLAP MATRIX
      DO ISYM=1,NVSYM
        JDUMP=IDUMP(ISYM)
        MAT=NDUMP(ISYM)
        JSIZE=ISIZE(ISYM)
        DO IX=1,MAT
          DO IY=1,MAT
            XK(IY,IX)=0
          END DO
        END DO
C**RHS
        IRHS=0
        DO IRDUMP=1,JDUMP
          KR=KDUMP(IRDUMP,ISYM)
          IF(KR.EQ.0)GO TO 9999
          IRHS=IRHS+1
          DO JRSIZE=1,JSIZE
            COEFFR=VCI(JRSIZE,IRDUMP,ISYM)
C**LHS
            ILHS=0
            DO ILDUMP=1,JDUMP
              KL=KDUMP(ILDUMP,ISYM)
              IF(KL.EQ.0)GO TO 8888
              ILHS=ILHS+1
              DO JLSIZE=1,JSIZE
                COEFFL=VCI(JLSIZE,ILDUMP,ISYM)
                IOVER=1
                DO MODE=1,NMODE
                  IF(IOVER.EQ.0)GO TO 7777
                  IF(IP(JLSIZE,MODE,ISYM).NE.IP(JRSIZE,MODE,ISYM))
     1            IOVER=0
                END DO
7777  CONTINUE
                XK(ILHS,IRHS)=XK(ILHS,IRHS)+IOVER*COEFFL*COEFFR
              END DO
8888  CONTINUE
            END DO
          END DO
9999  CONTINUE
        END DO
        WRITE(IOUT,*)'SYMM ',ISYM
        DO IRHS=1,MAT
          WRITE(IOUT,*)(XK(ILHS,IRHS),ILHS=1,MAT)
        END DO
      END DO
      RETURN
      END
C**TEMPORARY TEST OVERLAPS
C****************************************************************
C****************************************************************
      SUBROUTINE GRIDS(W,XQ,MAXPTS,MVB,NMODE,NATOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(1)
      DIMENSION XQ(MAXPTS,NMODE),MVB(NMODE)
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
C**CURRENT SETTINGS
      COMMON/CADDR/
     &LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LOLDH,LXK,LOLDXQ,LXW,LNBF,
     1LMBF,LSCF,LSX,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LMVB,LH,LXQ,LTEMP,LCONTR,
     3LDUM2(10),
     4LDUM3(10),
     5LISIZE,LW21,LSS,LSSX,LX21,LE21,LIDUMP,LNDUMP,LIP,LVCI,
     6LDUM6(10),
     7LDUM7(10),
     8LXA,LXA1,LXA2,LXA3,LXA4,LEVAL,LIP1,LIP2,LIP3,LIP4,
     9LPD1,LPD2,LPD3,LPD4,LPD5,LPD6,LPD7,LMXBAS,LKDUMP,LNVF,
     1LJP1,LJP2,LJP3,LJP4,LV1,LV2,LV3,LV4,LPD8,LPD9,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LDUM10(2)
C**CURRENT SETTINGS
      COMMON/FILASS/IOUT,INP
      COMMON/TRANSF/LTRAN
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/NCPOT/NPOT
      COMMON/DISC/IDISC
      COMMON/DISCSZ/KLC0,KEJK0,KLV1,KLC1,KEJK1,KLV2,KLC2,KEJK2,
     1KLV3,KLC3,KEJK3,KLV4,KLC4,KEJK4
C**************************************************ASSIGN TOTAL STORAGE
C**********************************************************************
C**********************************************************************
C**                            DUMP DIPOLE DATA TO DISC
C**********************************************************************
C**********************************************************************
C******************************
      CALL DUMDP0(NMODE,NATOM,W(LQQ),W(LX0),NPOT,W(LIPOT),W(LJPOT),
     1W(LCPOT))
      IF(ICOUPL.EQ.0)GO TO 3000
      MBFMX1=MAXPTS
      IF(JCOUPL.GT.0)THEN
        KLV1=MBFMX1
      ELSE
        KLV1=(1+MBFMX1)/2
      END IF
      CALL MEMO(1,LV1,KLV1,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        REWIND 61
        REWIND 71
        REWIND 81
        DO K=1,NMODE
          IK2=MVB(K)
C**WRITE DIPOLE GRIDS TO DISC
          CALL DUMDP1(XQ(1,K),IK2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     1    W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,
     2    W(LV1),W(LV1))
        END DO
      END IF
C******************************
      IF(ICOUPL.EQ.1)GO TO 3000
      MBFMX2=MBFMX1*MAXPTS
      IF(JCOUPL.GT.0)THEN
        KLV2=MBFMX2
      ELSE
        KLV2=(1+MBFMX2)/2
      END IF
      CALL MEMO(1,LV2,KLV2,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        REWIND 62
        REWIND 72
        REWIND 82
        DO K=1,NMODE
          IK2=MVB(K)
          DO L=1,K-1
            IL2=MVB(L)
C**WRITE DIPOLE GRIDS TO DISC
            CALL DUMDP2(XQ(1,K),XQ(1,L),IK2,IL2,NMODE,NATOM,
     1      W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),
     2      W(LJPOT),W(LCPOT),K,L,W(LV2),W(LV2))
          END DO
        END DO
      END IF
C******************************
      IF(ICOUPL.EQ.2)GO TO 3000
      MBFMX3=MBFMX2*MAXPTS
      IF(JCOUPL.GT.0)THEN
        KLV3=MBFMX3
      ELSE
        KLV3=(1+MBFMX3)/2
      END IF
      CALL MEMO(1,LV3,KLV3,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        REWIND 63
        REWIND 73
        REWIND 83
        DO K=1,NMODE
          IK2=MVB(K)
          DO L=1,K-1
            IL2=MVB(L)
            DO N=1,L-1
              IN2=MVB(N)
C**WRITE DIPOLE GRIDS TO DISC
              CALL DUMDP3(XQ(1,K),XQ(1,L),XQ(1,N),IK2,IL2,IN2,
     1        NMODE,NATOM,W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     2        NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,L,N,W(LV3),W(LV3))
            END DO
          END DO
        END DO
      END IF
C******************************
      IF(ICOUPL.EQ.3)GO TO 3000
      MBFMX4=MBFMX3*MAXPTS
      IF(JCOUPL.GT.0)THEN
        KLV4=MBFMX4
      ELSE
        KLV4=(1+MBFMX4)/2
      END IF
      CALL MEMO(1,LV4,KLV4,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        REWIND 64
        REWIND 74
        REWIND 84
        DO K=1,NMODE
          IK2=MVB(K)
          DO L=1,K-1
            IL2=MVB(L)
            DO N=1,L-1
              IN2=MVB(N)
              DO M=1,N-1
                IM2=MVB(M)
C**WRITE DIPOLE GRIDS TO DISC
                CALL DUMDP4(XQ(1,K),XQ(1,L),XQ(1,N),XQ(1,M),
     1          IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     2          W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     3          K,L,N,M,W(LV4),W(LV4))
              END DO
            END DO
          END DO
        END DO
      END IF
C******************************
      IF(ICOUPL.EQ.4)GO TO 3000
C**5-MODE AND HIGHER
3000  CONTINUE
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DUMDP0(NMODE,NATOM,QQ,X0,NPOT,IPOT,JPOT,CPOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(3)
      REAL*4 VR(3)
      DIMENSION QQ(NMODE),X0(NATOM,3)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/DIPNO/NUMDIP
      COMMON/WHICH/IWHICH
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/DIPREF/V,VR
      DO K=1,NMODE
        QQ(K)=0
      END DO
      DO IDIP=1,NUMDIP
        IF(IWHICH.GT.0)THEN
          CALL GETDIP(VDP,NATOM,X0,RR,IDIP)
        ELSE
          IF(IWHICH.EQ.0)THEN
            CALL GETDP0(VDP,NPOT,IPOT,JPOT,CPOT,NMODE,QQ,IDIP)
          ELSE
            CALL GETQDT(VDP,NMODE,QQ,IDIP)
          END IF
        END IF
        IF(JCOUPL.GT.0)THEN
          V(IDIP)=VDP
        ELSE
          VR(IDIP)=VDP
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DUMDP1(XQ,MM,NMODE,NATOM,QQ,RR,XX,X0,XL,XM,NPOT,IPOT,
     1JPOT,CPOT,MODE,V,VR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRIAT
      REAL*8 V(MM)
      REAL*4 VR(MM)
      DIMENSION XQ(MM),RR(NATOM,NATOM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/FILASS/IOUT,INP
      COMMON/WHICH/IWHICH
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/TRIATO/TRIAT
      COMMON/DIPNO/NUMDIP
      DO K=1,NMODE
        QQ(K)=0
      END DO
      DO IDIP=1,NUMDIP
        DO M=1,MM
          QQ(MODE)=XQ(M)
          IF(IWHICH.GT.0)THEN
            DO I=1,NATOM
              DO J=1,3
                XX(I,J)=X0(I,J)+XL(I,MODE,J)*QQ(MODE)/
     1          SQRT(XM(I))
              END DO
            END DO
            CALL GETDIP(VDP,NATOM,XX,RR,IDIP)
          ELSE
            IF(IWHICH.EQ.0)THEN
              CALL GETDP1(VDP,NPOT,IPOT,JPOT,CPOT,NMODE,QQ,IDIP)
            ELSE
              CALL GETQDT(VDP,NMODE,QQ,IDIP)
            END IF
          END IF
          IF(JCOUPL.GT.0)THEN
            V(M)=VDP
          ELSE
            VR(M)=VDP
          END IF
        END DO
        IF(JCOUPL.GT.0)THEN
          WRITE(10*IDIP+51)V
        ELSE
          WRITE(10*IDIP+51)VR
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DUMDP2(XQ1,XQ2,MM1,MM2,NMODE,NATOM,QQ,RR,XX,X0,XL,XM,
     1NPOT,IPOT,JPOT,CPOT,MODE1,MODE2,V,VR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRIAT
      REAL*8 V(MM2,MM1)
      REAL*4 VR(MM2,MM1)
      DIMENSION XQ1(MM1),XQ2(MM2),RR(NATOM,NATOM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/WHICH/IWHICH
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/TRIATO/TRIAT
      COMMON/DIPNO/NUMDIP
      DO K=1,NMODE
        QQ(K)=0
      END DO
      DO IDIP=1,NUMDIP
        DO M1=1,MM1
          QQ(MODE1)=XQ1(M1)
          DO M2=1,MM2
            QQ(MODE2)=XQ2(M2)
            IF(IWHICH.GT.0)THEN
              DO I=1,NATOM
                DO J=1,3
                  XX(I,J)=X0(I,J)+XL(I,MODE1,J)*QQ(MODE1)/
     1            SQRT(XM(I))
                  XX(I,J)=XX(I,J)+XL(I,MODE2,J)*QQ(MODE2)/
     1            SQRT(XM(I))
                END DO
              END DO
              CALL GETDIP(VDP,NATOM,XX,RR,IDIP)
            ELSE
              IF(IWHICH.EQ.0)THEN
                CALL GETDP2(VDP,NPOT,IPOT,JPOT,CPOT,NMODE,QQ,IDIP)
              ELSE
                CALL GETQDT(VDP,NMODE,QQ,IDIP)
              END IF
            END IF
            IF(JCOUPL.GT.0)THEN
              V(M2,M1)=VDP
            ELSE
              VR(M2,M1)=VDP
            END IF
          END DO
        END DO
        IF(JCOUPL.GT.0)THEN
          WRITE(10*IDIP+52)V
        ELSE
          WRITE(10*IDIP+52)VR
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DUMDP3(XQ1,XQ2,XQ3,MM1,MM2,MM3,NMODE,NATOM,QQ,RR,XX,
     1X0,XL,XM,NPOT,IPOT,JPOT,CPOT,MODE1,MODE2,MODE3,V,VR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRIAT
      REAL*8 V(MM3,MM2,MM1)
      REAL*4 VR(MM3,MM2,MM1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),RR(NATOM,NATOM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/WHICH/IWHICH
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/TRIATO/TRIAT
      COMMON/DIPNO/NUMDIP
      DO K=1,NMODE
        QQ(K)=0
      END DO
      DO IDIP=1,NUMDIP
        DO M1=1,MM1
          QQ(MODE1)=XQ1(M1)
          DO M2=1,MM2
            QQ(MODE2)=XQ2(M2)
            DO M3=1,MM3
              QQ(MODE3)=XQ3(M3)
              IF(IWHICH.GT.0)THEN
                DO I=1,NATOM
                  DO J=1,3
                    XX(I,J)=X0(I,J)+XL(I,MODE1,J)*QQ(MODE1)/
     1              SQRT(XM(I))
                    XX(I,J)=XX(I,J)+XL(I,MODE2,J)*QQ(MODE2)/
     1              SQRT(XM(I))
                    XX(I,J)=XX(I,J)+XL(I,MODE3,J)*QQ(MODE3)/
     1              SQRT(XM(I))
                  END DO
                END DO
                CALL GETDIP(VDP,NATOM,XX,RR,IDIP)
              ELSE
                IF(IWHICH.EQ.0)THEN
                  CALL GETDP3(VDP,NPOT,IPOT,JPOT,CPOT,NMODE,QQ,IDIP)
                ELSE
                  CALL GETQDT(VDP,NMODE,QQ,IDIP)
                END IF
              END IF
              IF(JCOUPL.GT.0)THEN
                V(M3,M2,M1)=VDP
              ELSE
                VR(M3,M2,M1)=VDP
              END IF
            END DO
          END DO
        END DO
        IF(JCOUPL.GT.0)THEN
          WRITE(10*IDIP+53)V
        ELSE
          WRITE(10*IDIP+53)VR
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DUMDP4(XQ1,XQ2,XQ3,XQ4,MM1,MM2,MM3,MM4,NMODE,NATOM,QQ,
     1RR,XX,X0,XL,XM,NPOT,IPOT,JPOT,CPOT,MODE1,MODE2,MODE3,MODE4,V,VR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRIAT
      REAL*8 V(MM4,MM3,MM2,MM1)
      REAL*4 VR(MM4,MM3,MM2,MM1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      DIMENSION RR(NATOM,NATOM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/WHICH/IWHICH
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/TRIATO/TRIAT
      COMMON/DIPNO/NUMDIP
      DO K=1,NMODE
        QQ(K)=0
      END DO
      DO IDIP=1,NUMDIP
        DO M1=1,MM1
          QQ(MODE1)=XQ1(M1)
          DO M2=1,MM2
            QQ(MODE2)=XQ2(M2)
            DO M3=1,MM3
              QQ(MODE3)=XQ3(M3)
              DO M4=1,MM4
                QQ(MODE4)=XQ4(M4)
                IF(IWHICH.GT.0)THEN
                  DO I=1,NATOM
                    DO J=1,3
                      XX(I,J)=X0(I,J)+XL(I,MODE1,J)*QQ(MODE1)/
     1                SQRT(XM(I))
                      XX(I,J)=XX(I,J)+XL(I,MODE2,J)*QQ(MODE2)/
     1                SQRT(XM(I))
                      XX(I,J)=XX(I,J)+XL(I,MODE3,J)*QQ(MODE3)/
     1                SQRT(XM(I))
                      XX(I,J)=XX(I,J)+XL(I,MODE4,J)*QQ(MODE4)/
     1                SQRT(XM(I))
                    END DO
                  END DO
                  CALL GETDIP(VDP,NATOM,XX,RR,IDIP)
                ELSE
                  IF(IWHICH.EQ.0)THEN
                    CALL GETDP4(VDP,NPOT,IPOT,JPOT,CPOT,NMODE,
     1              QQ,IDIP)
                  ELSE
                    CALL GETQDT(VDP,NMODE,QQ,IDIP)
                  END IF
                END IF
                IF(JCOUPL.GT.0)THEN
                  V(M4,M3,M2,M1)=VDP
                ELSE
                  VR(M4,M3,M2,M1)=VDP
                END IF
              END DO
            END DO
          END DO
        END DO
        IF(JCOUPL.GT.0)THEN
          WRITE(10*IDIP+54)V
        ELSE
          WRITE(10*IDIP+54)VR
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0DP1(NMODE,MODE,MAXCON,MAXPTS,H,XQ,XA,NSIZE,
     1NN,MM,IP,ISIZE,VP,VPR,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM)
      REAL*4 VPR(MM)
      DIMENSION H(MAXCON,MAXPTS,3),IP(ISIZE,NMODE)
      DIMENSION XA(1),XQ(MM)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0DP1'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.1)IFACT=1
        IF(ICOUPL.EQ.2)IFACT=-(NMODE-2)
        IF(ICOUPL.EQ.3)THEN
          IFACT=(NMODE-3)*(NMODE-1)
          DO I=2,NMODE-2
            IFACT=IFACT-I
          END DO
        END IF
        IF(ICOUPL.EQ.4)THEN
          IFACT=-(NMODE-4)*(NMODE-3)*(NMODE-1)
          DO I=1,NMODE-4
            IFACT=IFACT-I*(NMODE-3-I)
          END DO
          DO I=1,NMODE-4
            IFACT=IFACT+(NMODE-2)*I
          END DO
          DO I=2,NMODE-2
            IFACT=IFACT+(NMODE-4)*I
          END DO
        END IF
      END IF
      IF(JCOUPL.GT.0)THEN
        READ(10*IDIP+51)VP
      ELSE
        READ(10*IDIP+51)VPR
      END IF
      DO M=1,MM
        IF(JCOUPL.GT.0)THEN
          VV=VP(M)*IFACT
        ELSE
          VV=VPR(M)*IFACT
        END IF
C**NSIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,NSIZE
          NR=IP(IRHS,MODE)
          X=VV*H(NR,M,1)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP(ILHS,MODE)
            Y=H(NL,M,1)
            XA(ILHS+J0)=XA(ILHS+J0)+Y*X
          END DO
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
      SUBROUTINE VDP0(NMODE,XA,ISIZEL,ISIZER,IDIP,MSL,MSR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(3)
      REAL*4 VR(3)
      DIMENSION XA(ISIZEL,ISIZER)
      COMMON/WHICH/IWHICH
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/DIPREF/V,VR
      COMMON/FACTOR/FACTOR(5),FACTS(5)
      IF(MSL.NE.MSR)RETURN
C**ONE-MODE CONSTANTS
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.1)IFACT=1
        IF(ICOUPL.EQ.2)IFACT=-(NMODE-2)
        IF(ICOUPL.EQ.3)THEN
          IFACT=(NMODE-3)*(NMODE-1)
          DO I=2,NMODE-2
            IFACT=IFACT-I
          END DO
        END IF
        IF(ICOUPL.EQ.4)THEN
          IFACT=-(NMODE-4)*(NMODE-3)*(NMODE-1)
          DO I=1,NMODE-4
            IFACT=IFACT-I*(NMODE-3-I)
          END DO
          DO I=1,NMODE-4
            IFACT=IFACT+(NMODE-2)*I
          END DO
          DO I=2,NMODE-2
            IFACT=IFACT+(NMODE-4)*I
          END DO
        END IF
      END IF
      FACTS(1)=IFACT
C**TWO-MODE CONSTANTS
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.2)IFACT=1
        IF(ICOUPL.EQ.3)IFACT=-(NMODE-3)
        IF(ICOUPL.EQ.4)THEN
          IFACT=(NMODE-4)*(NMODE-3)
          DO I=1,NMODE-4
            IFACT=IFACT-I
          END DO
        END IF
      END IF
      FACTS(2)=IFACT
C**THREE-MODE CONSTANTS
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.3)IFACT=1
        IF(ICOUPL.EQ.4)IFACT=-(NMODE-4)
      END IF
      FACTS(3)=IFACT
C**FOUR-MODE CONSTANTS
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.4)IFACT=1
      END IF
      FACTS(4)=IFACT
C**(FIVE-MODE CONSTANTS)
      FACT=1
      DO I=1,ICOUPL
        FACT=FACT-FACTOR(I)*FACTS(I)
      END DO
      DO IRHS=1,ISIZER
        IF(JCOUPL.GT.0)XA(IRHS,IRHS)=XA(IRHS,IRHS)+V(IDIP)*FACT
        IF(JCOUPL.LT.0)XA(IRHS,IRHS)=XA(IRHS,IRHS)+VR(IDIP)*FACT
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VDP1(NMODE,MODE,MAXCON,MAXPTS,H,XQ,XA,XK,NSIZE,
     1NN,MM,IP,ISIZMX,ISIZEL,ISIZER,MSL,MSR,IP1,ISIZE1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(MAXCON,MAXPTS,3),XQ(MM),IP(ISIZMX,NMODE,1)
      DIMENSION XA(ISIZEL,1),XK(1),IP1(ISIZE1,1)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VDP1'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR=IP(IRHS,MODE,MSR)
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
            IF(K.NE.MODE.AND.(IP(IRHS,K,MSR).NE.IP(ILHS,K,MSL)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL=IP(ILHS,MODE,MSL)
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
      SUBROUTINE V0DP2(NMODE,MODE1,MODE2,MAXCON,MAXPTS,H1,XQ1,H2,XQ2,
     1NN1,MM1,NN2,MM2,XA,NSIZE,IP,ISIZE,TEMP,JCI1,JCI2,VP,VPR,
CC   2IP1,ISIZE1,IDIP)
     2IP1,ISIZE1,NSIZE1,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM2,MM1)
      REAL*4 VPR(MM2,MM1)
      DIMENSION H1(MAXCON,MAXPTS,3),H2(MAXCON,MAXPTS,3),IP(ISIZE,NMODE)
      DIMENSION IP1(ISIZE1,1)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2)
      DIMENSION TEMP(JCI2,JCI2)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM2A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0DP2'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.2)IFACT=1
        IF(ICOUPL.EQ.3)IFACT=-(NMODE-3)
        IF(ICOUPL.EQ.4)THEN
          IFACT=(NMODE-4)*(NMODE-3)
          DO I=1,NMODE-4
            IFACT=IFACT-I
          END DO
        END IF
      END IF
      IF(JCOUPL.GT.0)THEN
        READ(10*IDIP+52)VP
      ELSE
        READ(10*IDIP+52)VPR
      END IF
      DO M1=1,MM1
        DO IRHS2=1,JCI2
          DO ILHS2=1,JCI2
            TEMP(ILHS2,IRHS2)=0
          END DO
        END DO
        IF(JCOUPL.GT.0)THEN
          DO M2=1,MM2
CC          DO IX1=1,ISIZE1
            DO IX1=1,NSIZE1
              IRHS2=IP1(IX1,1)
              Y1=H2(IRHS2,M2,1)*IFACT
              X1=VP(M2,M1)*Y1
CC            DO IY1=1,ISIZE1
              DO IY1=1,NSIZE1
                ILHS2=IP1(IY1,1)
                Y0=H2(ILHS2,M2,1)
                TEMP(ILHS2,IRHS2)=TEMP(ILHS2,IRHS2)+Y0*X1
              END DO
            END DO
          END DO
        ELSE
          DO M2=1,MM2
CC          DO IX1=1,ISIZE1
            DO IX1=1,NSIZE1
              IRHS2=IP1(IX1,1)
              Y1=H2(IRHS2,M2,1)*IFACT
              X1=VPR(M2,M1)*Y1
CC            DO IY1=1,ISIZE1
              DO IY1=1,NSIZE1
                ILHS2=IP1(IY1,1)
                Y0=H2(ILHS2,M2,1)
                TEMP(ILHS2,IRHS2)=TEMP(ILHS2,IRHS2)+Y0*X1
              END DO
            END DO
          END DO
        END IF
C**NSIZE IS NO. UNIQUE INTEGRALS (2-DIM)
        DO IRHS=1,NSIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          X1=H1(NR1,M1,1)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            Y=H1(NL1,M1,1)
            XA(ILHS+J0)=XA(ILHS+J0)+Y*TEMP(NL2,NR2)*X1
          END DO
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
      SUBROUTINE VDP2(NMODE,MODE1,MODE2,MAXCON,MAXPTS,H1,XQ1,H2,XQ2,
     1NN1,MM1,NN2,MM2,XA,XK,NSIZE,IP,ISIZMX,ISIZEL,ISIZER,MSL,MSR,IP2,
     2ISIZE2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H1(MAXCON,MAXPTS,3),H2(MAXCON,MAXPTS,3)
      DIMENSION IP(ISIZMX,NMODE,1)
      DIMENSION XA(ISIZEL,1),XK(1),IP2(ISIZE2,2)
      DIMENSION XQ1(MM1),XQ2(MM2)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM2A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VDP2'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IP(IRHS,MODE1,MSR)
        NR2=IP(IRHS,MODE2,MSR)
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
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.(IP(IRHS,K,MSR).NE.
     1      IP(ILHS,K,MSL)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IP(ILHS,MODE1,MSL)
          NL2=IP(ILHS,MODE2,MSL)
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
      SUBROUTINE V0DP3(NMODE,MODE1,MODE2,MODE3,MAXCON,MAXPTS,H1,XQ1,H2,
     1XQ2,H3,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,XA,NSIZE,IP,ISIZE,TEMP,JCI1,
CC   2JCI2,JCI3,VP,VPR,IP2,ISIZE2,IDIP)
     2JCI2,JCI3,VP,VPR,IP2,ISIZE2,NSIZE2,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1)
      DIMENSION H1(MAXCON,MAXPTS,3),H2(MAXCON,MAXPTS,3)
      DIMENSION H3(MAXCON,MAXPTS,3)
      DIMENSION IP(ISIZE,NMODE),IP2(ISIZE2,2)
      DIMENSION TEMP(JCI2,JCI3,JCI2,JCI3)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM3A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0DP3'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF

      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.3)IFACT=1
        IF(ICOUPL.EQ.4)IFACT=-(NMODE-4)
      END IF
      IF(JCOUPL.GT.0)THEN
        READ(10*IDIP+53)VP
      ELSE
        READ(10*IDIP+53)VPR
      END IF
      DO M1=1,MM1
        DO IRHS3=1,JCI3
          DO IRHS2=1,JCI2
            DO ILHS3=1,JCI3
              DO ILHS2=1,JCI2
                TEMP(ILHS2,ILHS3,IRHS2,IRHS3)=0
              END DO
            END DO
          END DO
        END DO
        IF(JCOUPL.GT.0)THEN
CC        DO IX2=1,ISIZE2
          DO IX2=1,NSIZE2
            IRHS2=IP2(IX2,1)
            IRHS3=IP2(IX2,2)
            DO M2=1,MM2
              Y1=H2(IRHS2,M2,1)*IFACT
              DO M3=1,MM3
                Z1=H3(IRHS3,M3,1)
                X1=Z1*VP(M3,M2,M1)*Y1
CC              DO IY2=1,ISIZE2
                DO IY2=1,NSIZE2
                  ILHS2=IP2(IY2,1)
                  ILHS3=IP2(IY2,2)
                  Y0=H2(ILHS2,M2,1)
                  Z0=H3(ILHS3,M3,1)
                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3)=
     1            TEMP(ILHS2,ILHS3,IRHS2,IRHS3)+Y0*Z0*X1
                END DO
              END DO
            END DO
          END DO
        ELSE
CC        DO IX2=1,ISIZE2
          DO IX2=1,NSIZE2
            IRHS2=IP2(IX2,1)
            IRHS3=IP2(IX2,2)
            DO M2=1,MM2
              Y1=H2(IRHS2,M2,1)*IFACT
              DO M3=1,MM3
                Z1=H3(IRHS3,M3,1)
                X1=Z1*VPR(M3,M2,M1)*Y1
CC              DO IY2=1,ISIZE2
                DO IY2=1,NSIZE2
                  ILHS2=IP2(IY2,1)
                  ILHS3=IP2(IY2,2)
                  Y0=H2(ILHS2,M2,1)
                  Z0=H3(ILHS3,M3,1)
                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3)=
     1            TEMP(ILHS2,ILHS3,IRHS2,IRHS3)+Y0*Z0*X1
                END DO
              END DO
            END DO
          END DO
        END IF
C**NSIZE IS NO. UNIQUE INTEGRALS (3-DIM)
        DO IRHS=1,NSIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          NR3=IP(IRHS,MODE3)
          X1=H1(NR1,M1,1)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            NL3=IP(ILHS,MODE3)
            Y=H1(NL1,M1,1)
            XA(ILHS+J0)=XA(ILHS+J0)+Y*TEMP(NL2,NL3,NR2,NR3)*X1
          END DO
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
      SUBROUTINE VDP3(NMODE,MODE1,MODE2,MODE3,MAXCON,MAXPTS,H1,XQ1,H2,
     1XQ2,H3,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,XA,XK,NSIZE,IP,ISIZMX,ISIZEL,
     2ISIZER,MSL,MSR,IP3,ISIZE3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H1(MAXCON,MAXPTS,3),H2(MAXCON,MAXPTS,3)
      DIMENSION H3(MAXCON,MAXPTS,3)
      DIMENSION IP(ISIZMX,NMODE,1)
      DIMENSION XA(ISIZEL,1),XK(1),IP3(ISIZE3,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM3A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VDP3'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IP(IRHS,MODE1,MSR)
        NR2=IP(IRHS,MODE2,MSR)
        NR3=IP(IRHS,MODE3,MSR)
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
     1      IP(IRHS,K,MSR).NE.IP(ILHS,K,MSL)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IP(ILHS,MODE1,MSL)
          NL2=IP(ILHS,MODE2,MSL)
          NL3=IP(ILHS,MODE3,MSL)
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
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM3A.EQ.0)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        ITIM3A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0DP4(NMODE,MODE1,MODE2,MODE3,MODE4,MAXCON,MAXPTS,H1,
     1XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,NN2,MM2,NN3,MM3,NN4,MM4,XA,
CC   3NSIZE,IP,ISIZE,TEMP,JCI1,JCI2,JCI3,JCI4,VP,VPR,IP3,ISIZE3,IDIP)
     3NSIZE,IP,ISIZE,TEMP,JCI1,JCI2,JCI3,JCI4,VP,VPR,IP3,ISIZE3,
     4NSIZE3,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1)
      DIMENSION H1(MAXCON,MAXPTS,3),H2(MAXCON,MAXPTS,3)
      DIMENSION H3(MAXCON,MAXPTS,3),H4(MAXCON,MAXPTS,3)
      DIMENSION IP(ISIZE,NMODE),IP3(ISIZE3,3)
      DIMENSION TEMP(JCI2,JCI3,JCI4,JCI2,JCI3,JCI4)
      DIMENSION XA(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0DP4'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.4)IFACT=1
      END IF
      IF(JCOUPL.GT.0)THEN
        READ(10*IDIP+54)VP
      ELSE
        READ(10*IDIP+54)VPR
      END IF
      DO M1=1,MM1
        DO IRHS4=1,JCI4
          DO IRHS3=1,JCI3
            DO IRHS2=1,JCI2
              DO ILHS4=1,JCI4
                DO ILHS3=1,JCI3
                  DO ILHS2=1,JCI2
                    TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=0
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(JCOUPL.GT.0)THEN
CC        DO IX3=1,ISIZE3
          DO IX3=1,NSIZE3
            IRHS2=IP3(IX3,1)
            IRHS3=IP3(IX3,2)
            IRHS4=IP3(IX3,3)
            DO M2=1,MM2
              Y1=H2(IRHS2,M2,1)*IFACT
              DO M3=1,MM3
                Z1=H3(IRHS3,M3,1)
                DO M4=1,MM4
                  W1=H4(IRHS4,M4,1)
                  X1=VP(M4,M3,M2,M1)*Y1*Z1*W1
CC                DO IY3=1,ISIZE3
                  DO IY3=1,NSIZE3
                    ILHS2=IP3(IY3,1)
                    ILHS3=IP3(IY3,2)
                    ILHS4=IP3(IY3,3)
                    Y0=H2(ILHS2,M2,1)
                    Z0=H3(ILHS3,M3,1)
                    W0=H4(ILHS4,M4,1)
                    TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=
     1              TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)+
     2              Y0*Z0*W0*X1
                  END DO
                END DO
              END DO
            END DO
          END DO
        ELSE
CC        DO IX3=1,ISIZE3
          DO IX3=1,NSIZE3
            IRHS2=IP3(IX3,1)
            IRHS3=IP3(IX3,2)
            IRHS4=IP3(IX3,3)
            DO M2=1,MM2
              Y1=H2(IRHS2,M2,1)*IFACT
              DO M3=1,MM3
                Z1=H3(IRHS3,M3,1)
                DO M4=1,MM4
                  W1=H4(IRHS4,M4,1)
                  X1=VPR(M4,M3,M2,M1)*Y1*Z1*W1
CC                DO IY3=1,ISIZE3
                  DO IY3=1,NSIZE3
                    ILHS2=IP3(IY3,1)
                    ILHS3=IP3(IY3,2)
                    ILHS4=IP3(IY3,3)
                    Y0=H2(ILHS2,M2,1)
                    Z0=H3(ILHS3,M3,1)
                    W0=H4(ILHS4,M4,1)
                    TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)=
     1              TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4)+
     2              Y0*Z0*W0*X1
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
C**NSIZE IS NO. UNIQUE INTEGRALS (4-DIM)
        DO IRHS=1,NSIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          NR3=IP(IRHS,MODE3)
          NR4=IP(IRHS,MODE4)
          X1=H1(NR1,M1,1)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            NL3=IP(ILHS,MODE3)
            NL4=IP(ILHS,MODE4)
            Y=H1(NL1,M1,1)
            XA(ILHS+J0)=XA(ILHS+J0)+Y*TEMP(NL2,NL3,NL4,NR2,NR3,NR4)*X1
          END DO
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
      SUBROUTINE VDP4(NMODE,MODE1,MODE2,MODE3,MODE4,MAXCON,MAXPTS,H1,
     1XQ1,H2,XQ2,H3,XQ3,H4,XQ4,NN1,MM1,NN2,MM2,NN3,MM3,NN4,MM4,XA,XK,
     2NSIZE,IP,ISIZMX,ISIZEL,ISIZER,MSL,MSR,IP4,ISIZE4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H1(MAXCON,MAXPTS,3),H2(MAXCON,MAXPTS,3)
      DIMENSION H3(MAXCON,MAXPTS,3),H4(MAXCON,MAXPTS,3)
      DIMENSION IP(ISIZMX,NMODE,1)
      DIMENSION XA(ISIZEL,1),XK(1),IP4(ISIZE4,4)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VDP4'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IP(IRHS,MODE1,MSR)
        NR2=IP(IRHS,MODE2,MSR)
        NR3=IP(IRHS,MODE3,MSR)
        NR4=IP(IRHS,MODE4,MSR)
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
     1      MODE4.AND.(IP(IRHS,K,MSR).NE.IP(ILHS,K,MSL)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IP(ILHS,MODE1,MSL)
          NL2=IP(ILHS,MODE2,MSL)
          NL3=IP(ILHS,MODE3,MSL)
          NL4=IP(ILHS,MODE4,MSL)
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
      SUBROUTINE DIPTOT(NDUMP,NVSYM,XK,MAXDMP,VCI,MAXCI,XA,
     1ISIZEL,ISIZER,MSL,MSR,TEMP,EVAL,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NDUMP(NVSYM),EVAL(MAXDMP,NVSYM)
      DIMENSION XK(MAXDMP,MAXDMP)
      DIMENSION VCI(MAXCI,MAXDMP,NVSYM),XA(ISIZEL,ISIZER)
      DIMENSION TEMP(ISIZEL,1)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
101   FORMAT(//,1X,'LEFT-HAND SYMMETRY: A1',1X,I5,' FUNCTIONS',/)
102   FORMAT(//,1X,'LEFT-HAND SYMMETRY: B2',1X,I5,' FUNCTIONS',/)
103   FORMAT(//,1X,'LEFT-HAND SYMMETRY: B1',1X,I5,' FUNCTIONS',/)
104   FORMAT(//,1X,'LEFT-HAND SYMMETRY: A2',1X,I5,' FUNCTIONS',/)
201   FORMAT(//,1X,'RIGHT-HAND SYMMETRY: A1',1X,I5,' FUNCTIONS',/)
202   FORMAT(//,1X,'RIGHT-HAND SYMMETRY: B2',1X,I5,' FUNCTIONS',/)
203   FORMAT(//,1X,'RIGHT-HAND SYMMETRY: B1',1X,I5,' FUNCTIONS',/)
204   FORMAT(//,1X,'RIGHT-HAND SYMMETRY: A2',1X,I5,' FUNCTIONS',/)
300   FORMAT(5D14.6)
301   FORMAT(//,1X,'A1-COMPONENT OF DIPOLE',/)
302   FORMAT(//,1X,'B2-COMPONENT OF DIPOLE',/)
303   FORMAT(//,1X,'B1-COMPONENT OF DIPOLE',/)
C**RHS
C**CALL MATRIX MULT. ROUTINE MXMA TO SET UP TEMP(IY,I2)
C     CALL MXMA(XA(1,1),1,ISIZEL,VCI(1,1,MSR),1,MAXCI,
C    &      TEMP(1,1),1,ISIZEL,ISIZEL,ISIZER,NDUMP(MSR))
      CALL DGEMM('N','N',ISIZEL,NDUMP(MSR),ISIZER,1.0D0,XA(1,1),
     &       ISIZEL,VCI(1,1,MSR),MAXCI,0.0D0,TEMP,ISIZEL)
C**LHS
C**CALL MXMA TO MULT. TEMP() BY LHS CFS
C     CALL MXMA(VCI(1,1,MSL),MAXCI,1,TEMP(1,1),1,ISIZEL,
C    &      XK(1,1),1,MAXDMP,NVAL,ISIZEL,NVAL)
      CALL DGEMM('T','N',NDUMP(MSL),NDUMP(MSR),ISIZEL,1.0D0,
     & VCI(1,1,MSL),MAXCI,TEMP,ISIZEL,0.0D0,XK(1,1),MAXDMP)
C**WRITE LHS SYMMETRY AND ENERGIES
      IF(MSL.EQ.1)WRITE(IOUT,101)NDUMP(MSL)
      IF(MSL.EQ.2)WRITE(IOUT,102)NDUMP(MSL)
      IF(MSL.EQ.3)WRITE(IOUT,103)NDUMP(MSL)
      IF(MSL.EQ.4)WRITE(IOUT,104)NDUMP(MSL)
      WRITE(IOUT,300)(EVAL(I,MSL)*WAVENM,I=1,NDUMP(MSL))
C**WRITE RHS SYMMETRY AND ENERGIES
      IF(MSR.EQ.1)WRITE(IOUT,201)NDUMP(MSR)
      IF(MSR.EQ.2)WRITE(IOUT,202)NDUMP(MSR)
      IF(MSR.EQ.3)WRITE(IOUT,203)NDUMP(MSR)
      IF(MSR.EQ.4)WRITE(IOUT,204)NDUMP(MSR)
      WRITE(IOUT,300)(EVAL(I,MSR)*WAVENM,I=1,NDUMP(MSR))
C**WRITE DIPOLE COMPONENT
      IF(IDIP.EQ.1)WRITE(IOUT,301)
      IF(IDIP.EQ.2)WRITE(IOUT,302)
      IF(IDIP.EQ.3)WRITE(IOUT,303)
C**WRITE MATRIX
      DO ILHS=1,NDUMP(MSL)
        WRITE(IOUT,300)(XK(ILHS,IRHS),IRHS=1,NDUMP(MSR))
        WRITE(IOUT,*)
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDP0(V,NPOT,IPOT,JPOT,CPOT,NMODE,QQ,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPOT(NPOT,6,3),JPOT(NPOT,6,3),CPOT(NPOT,3),QQ(1)
      V=0
      DO I=1,NPOT
C**MAXIMUM OF 6 MODES COUPLED
        IF(IPOT(I,1,IDIP).EQ.0)THEN
          V=V+CPOT(I,IDIP)
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDP1(V,NPOT,IPOT,JPOT,CPOT,NMODE,QQ,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPOT(NPOT,6,3),JPOT(NPOT,6,3),CPOT(NPOT,3),QQ(1)
      V=0
      DO I=1,NPOT
C**MAXIMUM OF 6 MODES COUPLED
        IF(IPOT(I,2,IDIP).EQ.0)THEN
          K=IPOT(I,1,IDIP)
          L=JPOT(I,1,IDIP)
          TERM=QQ(K)**L
          V=V+CPOT(I,IDIP)*TERM
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDP2(V,NPOT,IPOT,JPOT,CPOT,NMODE,QQ,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPOT(NPOT,6,3),JPOT(NPOT,6,3),CPOT(NPOT,3),QQ(NMODE)
      V=0
      DO I=1,NPOT
        TERM=1
C**MAXIMUM OF 6 MODES COUPLED
        IF(IPOT(I,2,IDIP).NE.0.AND.IPOT(I,3,IDIP).EQ.0)THEN
          DO J=1,2
            K=IPOT(I,J,IDIP)
            L=JPOT(I,J,IDIP)
            TERM=TERM*QQ(K)**L
          END DO
          V=V+CPOT(I,IDIP)*TERM
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDP3(V,NPOT,IPOT,JPOT,CPOT,NMODE,QQ,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPOT(NPOT,6,3),JPOT(NPOT,6,3),CPOT(NPOT,3),QQ(NMODE)
      V=0
      DO I=1,NPOT
        TERM=1
C**MAXIMUM OF 6 MODES COUPLED
        IF(IPOT(I,3,IDIP).NE.0.AND.IPOT(I,4,IDIP).EQ.0)THEN
          DO J=1,3
            K=IPOT(I,J,IDIP)
            L=JPOT(I,J,IDIP)
            TERM=TERM*QQ(K)**L
          END DO
          V=V+CPOT(I,IDIP)*TERM
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDP4(V,NPOT,IPOT,JPOT,CPOT,NMODE,QQ,IDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPOT(NPOT,6,3),JPOT(NPOT,6,3),CPOT(NPOT,3),QQ(NMODE)
      V=0
      DO I=1,NPOT
        TERM=1
C**MAXIMUM OF 6 MODES COUPLED
        IF(IPOT(I,4,IDIP).NE.0)THEN
          DO J=1,4
            K=IPOT(I,J,IDIP)
            L=JPOT(I,J,IDIP)
            TERM=TERM*QQ(K)**L
          END DO
          V=V+CPOT(I,IDIP)*TERM
        END IF
      END DO
      RETURN
      END
