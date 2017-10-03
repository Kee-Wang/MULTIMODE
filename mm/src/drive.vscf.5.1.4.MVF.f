C****************************************************************
C****************************************************************
C**DRIVE
C****************************************************************
C****************************************************************
      SUBROUTINE VSCF(W)
C     PROGRAM VSCF (VERSION 5.1.4)
C     4.9.14
C************************************************************
C                                                           *
C     COPYRIGHT:  S. CARTER T/A 'MULTIMODE'                 *
C                                                           *
C                 39, GROVE HILL                            *
C                 CAVERSHAM                                 *
C                 READING, RG4 8PS                          *
C                 ENGLAND                                   *
C                                                           *
C************************************************************
C
C..VERSION 1
C..1.0 VSCF, in Normal or Internal coordinates (1.10.96)
C..1.1 V-CI, including symmetry (14.10.96)
C..1.2 SCF-CI, including symmetry (1.11.96)
C..1.3 Diagonal Adiabatic Rotations (1.10.97)
C
C..VERSION 2
C..2.0 Davidson/Liu (core) algorithm for VCI (1.12.97)
C..2.1 New VCI integration scheme (1.1.98)
C..2.2 Complete (Whitehead-Handy) Rotations (1.2.98)
C..2.3 Choice of REAL*4 or REAL*8 disc grids (1.3.98)
C..2.4 Selective K-diagonal steps in rotations (1.4.98)
C..2.5 Davidson/Liu (disc) algorithm for VCI (1.6.98)
C
C..VERSION 3
C..3.0 Improved VCI algorithm for MAXBAS,MAXSUM (1.10.98)
C..3.1 Contract primitives to give HEG integration (10.10.98)
C..3.2 HEG integration for SCF virtuals (deleted)
C..3.3 Write SCF and VCI coefficients to disc (1.11.98)
C..3.4 Selective VCI symmetries ...
C..    ... and QL/Givens choice for Davidson/Liu (11.11.98)
C
C..VERSION 4 
C..4.0 Include dump and restart procedure (15.2.99)
C..4.1 Restart for dipole matrix elements (28.2.99)
C..4.2 Include D2h symmetry into CI for J=0 ...
C..    ... and include A2, B1 and B2 integration (15.3.99)
C..4.3 Reaction path hamiltonian for J=0 (1.4.99)
C..4.4 Link to MOLPRO and HERMITE interpolation ...
C..    ... and Least squares fits of potential to 5-mode (1.7.99)
C..4.5 Two-scheme VCI basis (1.10.99)
C..4.6 Two-contraction schemes J=0 (incomplete...do not use)
C..4.7 General RPH for 'analytic' and 'abinitio' potentials (1.10.02)
C..4.8 Coriolis overcounting and left-right differentiation ...
C..    ... and ICOUPC = 0,1,-->,ICOUPL for coriolis (1.11.02)
C..4.9.0
C..    Reaction path hamiltonian for J>0 ...
C..    ... and extension of integration of potential to 6 modes (1.7.03)
C..4.9.1
C..    1- & 2-mode VCI perturbation tolerance in Lanczos matrix (1.4.04)
C..    Least squares fits of 5-mode intrinsic RPH potential (1.5.04)
C..    Extended properties, including RPH (1.11.04)
C..    6-mode CI basis (1.1.05)
C..    Inclusion of sparsity in iteration diagonalisation technique (1.6.05)
C..4.9.2
C..    Completion of two-contraction schemes (1.1.06)
C..    Use of intrinsic potentials and coriolis (1.1.08)
C..    Effective potentials involving 2 modes (1.9.08)
C..    VSCF with primitives (1.9.08)
C
C..VERSION 5
C..5.0 Include P,Q,R-branch intensities for MM-SR (1.1.09)
C..5.1.0
C..    Include P,Q,R-branch intensities for MM-RPH (1.1.10)
C..5.1.1
C..    Full symmetry-adapted version for both MM-SR and MM-RPH (1.7.10)
C..5.1.2
C..    Nuclear spin statistics included at input (1.6.11)
C..    Include Partition Functions (16.11.11)
C..    Allow spectrum to be redrawn in different energy ranges (16.1.12)
C..5.1.3
C..    Include Raman Spectroscopy (1.7.12)
C..    Include integral or half-integral torsion if RPH (1.9.12)
C..    Allow MODINT>2 for torsion (1.1.13)
C..5.1.4
C..    Modification of Q for refinements (4.9.14)
C
C************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 TITLE(80)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SPACE
      CHARACTER*80 VERSN,FILEV,FILER
      LOGICAL LGIV,LINEAR,LANCZ,LANZA,LANZB,TRIAT,ABINIT
      LOGICAL PLANAR(3)
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
     1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
     4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
     5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
     6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
     7LDP4A,LDP4B,LDP4C,LDP4D,LXP,LJP6,LANCNT,LANK,LXJ,LVP5,
     8LC5,LC6,LVP0,LTMPRD,LMVB,LNFC,LTEMP6,LTEMP7,LTEMP8,LTEMP9,
     9LXMODQ,LMVF
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/TITLE/TITLE
      COMMON/HERM/IHERM
      COMMON/CHECK/MCHECK
      COMMON/SADDLE/JNORM
      COMMON/PATH/ISCFCI
      COMMON/CYCLE/ICYCLE
      COMMON/TORS/QTAU
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU,IFLAUD
      COMMON/REACTL/JREACT,IFITRP
      COMMON/NREACT/NREACT
      COMMON/DEGEN/IDEGEN,NDEGEN
      COMMON/NORML/INORM
      COMMON/ENTER/IENTER,IENTMX(5),NTOT1,NTOT2,NTOT3,NTOT4,NTOT5
      COMMON/SINCOS/ICS
      COMMON/CHKMAX/MAXCHK(6)
      COMMON/AVCON/AVCON
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FITTER/MFIT,MFIT1(4),XFIT1,MFIT2(4),XFIT2,MFIT3(4),XFIT3,
     1MFIT4(4),XFIT4,MFIT5(4),XFIT5
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TRIATO/TRIAT
      COMMON/LANTOL/TOLLAN
      COMMON/CYCLES/NCYCLE
      COMMON/ISET/ISET
C**TEMPORARY (DIMENSIONS)
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LAN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TLLCUT,EVLCUT,TOLMAT,LANCYC,LANINC
      COMMON/GIVEN/LGIV,IGIV
      COMMON/MATSIZ/MATSIZ
      COMMON/TYPE/LINEAR,PLANAR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/TOLS/TOL,EPS
      COMMON/ECKCNT/ICNT,INTC
      COMMON/ABINIT/ABINIT
      COMMON/TMPMAX/MAXTMP
      COMMON/EVL/EVL,CUT
      COMMON/ANALYT/MAXQU,MAXPOW
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/NCPOT/NPOT
      COMMON/FACTOR/FACTOR(6),FACTS(6)
      COMMON/VCIMAX/NMAX
      COMMON/MXNMAX/MXNMAX,MAXJ
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/RPHROT/IROTV
      COMMON/ESTATE/IORDER
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/AXES/MX(3),MXDIP(5),MXROT(3),ISYMT
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/TBASIS/NTBAS(6,2),NTAU(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10),NSYMSL(10),MSYMSL(10,100),
     1NSYMSR(10,100),MSYMSR(10,100,10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CSAVES/NNMODS,NAMODS,NVMODS,ICOUPS,ICOUCS,JREACS,NREACS
      COMMON/XSAVES/ICOUPX,ICOUCX,JCOUPS,JCOUCS
      COMMON/NCREC/NREC(6),MAXBUF(6),NUNITR,NUNITW
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/MVAL/MVAL1,MVAL2
      COMMON/CSIZES/ISIZM1,ISIZM2,NVAL1,NVAL2,ICSIZ1,ICSIZ2,
     1IPSIZ1,IPSIZ2
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTY/LLCNT
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/UNITEX/I75,I76,I85,I86
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
      COMMON/KLSZS/KLJP1,KLJP2,KLJP3,KLJP4,KLJP5,KLJP6
      COMMON/KJSZS/KJP1,KJP2,KJP3,KJP4,KJP5,KJP6
      COMMON/KPSZS/KIP1,KIP2,KIP3,KIP4,KIP5,KIP6,JCC1,JCC2,JCC3,JCC4,
     1JCC5,JCC6
      COMMON/IPSZS/KPPP1,KPPP2,KPPP3,KPPP4,KPPP5,KPPP6
      COMMON/SIZEJ/JSIZE1(2),JSIZE2(2),JSIZE3(2),JSIZE4(2),JSIZE5(2),
     1JSIZE6(2)
      COMMON/TOTALS/ITOT1(2),ITOT2(2),ITOT3(2),ITOT4(2),ITOT5(2),
     1ITOT6(2),ITOT
      COMMON/TOTK/KTOT(6,2)
C**RPH ODD K
      COMMON/LMAX/LMAX
C**RPH ODD K
      COMMON/MATRIX/NVALV,NVALR,KSTEP,KSIGN,NVALCF
      DATA SPACE/'  '/
      DATA ISYMP/1,2,3,4,5,6,7,8,0,0,
     1           2,1,4,3,6,5,8,7,0,0,
     2           3,4,1,2,7,8,5,6,0,0,
     3           4,3,2,1,8,7,6,5,0,0,
     4           5,6,7,8,1,2,3,4,0,0,
     5           6,5,8,7,2,1,4,3,0,0,
     6           7,8,5,6,3,4,1,2,0,0,
     7           8,7,6,5,4,3,2,1,0,0,
     8           0,0,0,0,0,0,0,0,0,0,
     9           0,0,0,0,0,0,0,0,0,0/
      DATA VERSN/'MULTIMODE VERSION 5.1.4 (4 September 2014)'/
C**********************************************************************
100   FORMAT(80A1)
199   FORMAT(///,1X,A80,///)
200   FORMAT(//,1X,80A1)
205   FORMAT(/,1X,'NUMBER OF MODES = ',I4,/,
     1         1X,'NUMBER OF SCF STATES = ',I4,/,
     3         1X,'PRINT LEVEL = ',I4,/,
     4         1X,'TOLERANCE FOR SCF = ',D20.12,/,
     5         1X,'POTENTIAL COUPLING OF ',I4,' MODES',/,
     6         1X,'CORIOLIS COUPLING OF ',I4,' MODES',/)
210   FORMAT(/,1X,'IDISC = ',I4,/,
     1         1X,'(WRITE POTENTIAL AND CORIOLIS INFO. TO DISC (0)',/,
     2         1X,'POTENTIAL AND CORIOLIS DISCS ALREADY EXIST (1))',/)
215   FORMAT(/,1X,'SCF CALCULATION ONLY',/)
220   FORMAT(//,1X,'START SCF CALCULATION',/)
225   FORMAT(/,1X,'SCF PLUS CI CALCULATION')
230   FORMAT(/,1X,'NUMBER OF CI ENERGIES REQUIRED = ',I4,/)
235   FORMAT(//,1X,'START CI CALCULATION',/)
240   FORMAT(/,1X,'ZERO POINT ENERGY = ',F10.2,/)
245   FORMAT(//,1X,'OVERLAPS OF ORIGINAL SCF STATE FUNCTIONS',/)
250   FORMAT(/,1X,'OVERLAPS WITH STATE ',I4)
255   FORMAT(//,1X,'TEST OF SCHMIDT ORTHOGONALISATION',/)
260   FORMAT(//,1X,'FINAL VIBRATIONAL (K-DIAGONAL) CI ENERGIES',/)
265   FORMAT(//,1X,'NORMAL COORDINATE POTENTIAL',/)
270   FORMAT(//,1X,'INTERNAL COORDINATE POTENTIAL',/,
     1          1X,'PERMUTATIONS = ',5(F10.3,1X),/)
275   FORMAT(//,1X,'CI MATRIX INVOLVES ',I4,' SCF FUNCTIONS',/)
280   FORMAT(//,1X,'MAXIMUM QUANTA IN VCI BASIS',/,1X,I5,//,
     1          1X,'MAXIMUM SUM QUANTA (0: UNRESTRICTED)',/,1X,I5,/)
285   FORMAT(//,1X,'SIZE OF VIBRATIONAL CI MATRIX IS ',I6,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
290   FORMAT(//,1X,'VIBRATIONAL SYMMETRY ',I2)
295   FORMAT(//,1X,'NUMBER OF ROVIBRATIONAL SYMMETRIES = ',I3,/,
     1          1X,'NUMBER OF VIBRATIONAL SYMMETRIES = ',I3,/,
     2          1X,'NUMBER OF MODE SYMMETRIES = ',I3,/)
300   FORMAT(1X,'MODE SYMMETRY ',I3,'   MODES ',50I3)
305   FORMAT(//,1X,'TOTAL ANGULAR MOMENTUM J = ',I3,/)
310   FORMAT(/,1X,'J = 0 SCF CYCLE',/)
315   FORMAT(/,1X,'J = ',I3,' SCF CYCLE',/)
320   FORMAT('*************************')
325   FORMAT(50(1H*))
330   FORMAT(/,1X,'LANCZOS PERTURBATION WITH ',I5,' 1- & 2- DIM BASIS',
     1/,1X,'TOLERANCE ',F10.6,/,1X,'NUMBER INCREMENTS FOR MATRIX ',I5,/)
350   FORMAT(//,1X,'LINK TO PROGRAM NORMALS',/)
351   FORMAT(/,1X,'POTENTIAL IS A MINIMUM AT REFERENCE STRUCTURE',/)
352   FORMAT(/,1X,'POTENTIAL IS A SADDLE POINT AT REFERENCE STRUCTURE',
     1/)
355   FORMAT(/,1X,'CANNOT USE INORM WITH IWHICH = 0',/)
360   FORMAT(/,1X,'NUMBER OF LANCZOS CYCLES = ',I3)
365   FORMAT(/,1X,'LANCZOS TOLERANCE = ',F10.6)
370   FORMAT(/,1X,'LANCZOS (HALF) MATRIX MAXIMUM ORDER = ',I6,/)
371   FORMAT(/,1X,'LANCZOS BLOCK-CONVERGENCE SIZE = ',I6,/,
     11X,'LANCZOS MATRIX ELEMENT TOLERANCE = ',D10.4,/)
375   FORMAT(//,1X,'ROTATIONAL SYMMETRY ',I2)
380   FORMAT(//,1X,'SIZE OF RO-VIBRATIONAL CI MATRIX IS ',I5,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
385   FORMAT(//,1X,'FINAL RO-VIBRATIONAL CI ENERGIES',/)
390   FORMAT(1X,'BLOCK ',I2,' VIBRATIONAL SYMMETRY ',I2)
395   FORMAT(//,1X,'SIZES OF CI SYMMETRY BLOCKS: ',/,1X,10I10,/)
400   FORMAT(//,1X,'SHOULD NOT OCCUR',/)
405   FORMAT(//,1X,'MINIMIZATION OF POTENTIAL',/)
410   FORMAT(//,1X,'DEGREE OF POLYNOMIAL: ',I3,/)
415   FORMAT(//,1X,'VCI BASIS DETAILS',/)
420   FORMAT(//,1X,'CONTRACTION SCHEME: ',I2,
     1' MAXIMUM QUANTA FOR 1- 2- 3- 4- 5- 6-MODE COUPLED STATES',/,1X,
     26I5,/)
425   FORMAT(1X,'CONTRACTION SCHEME: ',I2,
     1' MAXIMUM SUM QUANTA FOR 1- 2- 3- 4- 5- 6-MODE STATES',/,1X,
     26I5,/)
430   FORMAT(/,1X,'FIND MINIMUM',/)
435   FORMAT(/,1X,'REFERENCE ENERGY (J=0): ',F20.10,' CM-1',/)
440   FORMAT(1X,'MAXIMUM QUANTA FOR INDIVIDUAL MODES',/)
445   FORMAT(1X,'COUPLED INTEGRATION FOR INDIVIDUAL MODES',/)
450   FORMAT(/,1X,'INCREMENT OF K IN ROVIBRATIONAL ANALYSIS: ',I3,/)
455   FORMAT(//,1X,'NUMBER OF REQUIRED SYMMETRIES: ',I2,/,
     1          1X,'SPECIFIC SYMMETRIES: ',8I4,/)
460   FORMAT(//,1X,'CALCULATION TO PERFORM PRELIMINARY CHECKS ONLY')
465   FORMAT(//,1X,'CURVELINEAR COORDINATE FOR MODE ',I3,/)
470   FORMAT(//,1X,'LINK TO MOLPRO WITH LINK INDEX = ',I2,/)
475   FORMAT(1X,'GAUSS 1-DIM INFINITIES: ',/)
480   FORMAT(8F10.5)
485   FORMAT(//,1X,'NUMBER OF CONTRACTION SCHEMES ',I2,/,
     11X,'NUMBER CONTRACTED FUNCTIONS = ',I4,3X,I4,/)
490   FORMAT(1X,'CONTRACTION SCHEME ',I3,'   MODES ',50I3)
495   FORMAT(//,1X,60(1H*),/,1X,'START CONTRACTION SCHEME ',I2,/,1X,
     160(1H*),//)
500   FORMAT(//,1X,60(1H*),/,1X,'COMBINE CONTRACTION SCHEMES 1 AND 2',
     1/,1X,60(1H*),//)
501   FORMAT(//,1X,'CONTRACTION SCHEME ',I1,/,
     11X,'ADJUSTED VALUES NVAL = ',10I6,/)
502   FORMAT(//,1X,'ADJUSTED VALUE NVAL FOR SCHEME: ',I1,' = ',I4,/)
503   FORMAT(//,1X,'MINIMUM ENERGY IS ',D20.12,' CM-1',/)
504   FORMAT(1X,'HEG POINTS USED IN INCREMENTS OF ',I1,/)
505   FORMAT(//,1X,'CONTRACTION SCHEME ',I1,/,
     11X,'ISIZE = ',10I6,/)
506   FORMAT(//,1X,50(1H*),/,
     11X,'HEG INTEGRATION POINTS FOR MODE ',I2,/)
507   FORMAT(//,1X,'CANNOT USE MOLPRO FACILITY FOR REACTION PATH')
508   FORMAT(//,1X,'CANNOT USE DUMP/RESTART FOR ICOUPL > 4')
509   FORMAT(//,1X,'CANNOT USE DUMP/RESTART WITH CONTRACTION SCHEME')
510   FORMAT(//,1X,'CANNOT USE DUMP/RESTART IF REACTION PATH')
511   FORMAT(//,1X,'CANNOT USE FITTED POTENTIAL WITH RPH')
512   FORMAT(//,1X,'ENERGY AT REFERENCE IS ',D20.12,' CM-1',/)
513   FORMAT(1X,'ENERGY AT REFERENCE IS ',D20.12,' HARTREES',/)
514   FORMAT(1X,'POTENTIAL COUPLING OF ',I4,' MODES',/,
     1       1X,'CORIOLIS COUPLING OF ',I4,' MODES',/)
515   FORMAT(1X,'USE INTRINSIC POTENTIALS',/)
576   FORMAT(//////' OUTPUT VIBRATION FILE [.VIB]: ',//,A80)
577   FORMAT(/' OUTPUT ROTATION FILE [.ROT]: ',//,A80,//)

C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                      INITIAL INPUT AND SET-UP STAGE
C**********************************************************************
C**********************************************************************
C**********************************************************************
      WRITE(IOUT,199)VERSN
      QTAU=0
      DO I=1,6
        NREC(I)=-1
        MAXBUF(I)=0
      END DO
      MAXTMP=0
      I61=61
      I62=62
      I63=63
      I64=64
      I71=71
      I72=72
      I73=73
      I74=74
      I75=75
      I76=76
      I81=81
      I82=82
      I83=83
      I84=84
      I85=85
      I86=86
      I91=91
      I92=92
      I93=93
      I94=94
      CHSYM(1)='A1g'
      CHSYM(2)='B2g'
      CHSYM(3)='B1g'
      CHSYM(4)='A2g'
      CHSYM(5)='A1u'
      CHSYM(6)='B2u'
      CHSYM(7)='B1u'
      CHSYM(8)='A2u'
      IDEGEN=0
      NDEGEN=0
      IROTV=0
      VMIN=1.D+20
      INIT=0
      EVLCUT=0
      LINEAR=.FALSE.
      LANCZ=.FALSE.
      LANZA=.FALSE.
      LANZB=.FALSE.
      LGIV=.FALSE.
      TRIAT=.FALSE.
      WAVENM=0.21947463067D+06
      ATTOJ=0.229371D0
      BOHR=0.52917725D0
      ELMASS=1.82288853D+03
      RAD=180/DACOS(-1.D0)
      KCONT=-1
      NCONT=1
      DO K=1,10
        NCSIZE(K)=0
      END DO
      DO K=1,2
        ICONT(K)=0
        ICCONT(K)=0
        NCOUPL(K)=0
        NCOUPC(K)=0
        DO I=1,10
          NCVAL(K,I)=0
          ISIZC(K,I)=0
        END DO
        DO I=1,6
          NTBAS(I,K)=-1
          NBAS(I,K)=-1
          MAXSUM(I,K)=0
        END DO
      END DO
      MAXPOW=0
      DO I=1,5
        IENTMX(I)=0
        NTAU(I)=-1
      END DO
      DO I=1,100
        SYMBOL(I)=SPACE
      END DO
      ISIZE1=0
      ISIZE2=0
      ISIZE3=0
      ISIZE4=0
      ISIZE5=0
      ISIZE6=0
      ICID=0
      IRET=0
      IORDER=0
      ICNT=0
      INTC=0
      ETA=1.D-78
      EPS=1.D-15
      TOL=ETA/EPS
      READ(INP,*)
C****************************************************************************
C**TITLE
C-------
C**
C**         TITLE: Any suitable title...maximum 80 characters
C****************************************************************************
      READ(INP,100)TITLE
      WRITE(IOUT,200)TITLE
      WRITE(IOUT,200)
      READ(INP,*)
C****************************************************************************
C**NATOM,NSTAT,CONV,ICOUPL,ICOUPC,ISCFCI,IWHICH,IDISC,NROTTR,JMAX,MCHECK,INORM
C-----------------------------------------------------------------------------
C**
C**NATOM:   Number of atoms (number of Normal Coordinates NMODE=3N-6, unless
C**         facility NROTTR is used - see below)
C**NSTAT:   Number of SCF states required
C**         NSTAT>0: Input NSTAT specific SCF states depending on NDEF (below)
C**         NSTAT<0: Program generates -NSTAT SCF states in increasing energy
C**CONV:    Convergence threshold for SCF states in cm-1
C**ICOUPL:  Number of modes coupled in potential integration (1,2,3,4,5,6)
C**         ICOUPL>0: use REAL*8 grid
C**         ICOUPL<0: use REAL*4 grid (half the disc storage as REAL*8)
C**         For contraction scheme, always final coupling
C**ICOUPC:  Number of modes coupled in coriolis integration (1,2,3,4)
C**         ICOUPC>0: use REAL*8 grid
C**         ICOUPC<0: use REAL*4 grid (half the disc storage as REAL*8)
C**         Note: ICOUPC can NOT be greater than ICOUPL
C**         For contraction scheme, always final coupling
C**ISCFCI:  Number of final CI energies required (if positive) for VSCF-CI or
C**         VCI (see ICI below)
C**         ISCFCI>0: SCF + CI calculation
C**         ISCFCI=0: SCF only
C**         ISCFCI<0: Terminates after preliminary checks (see MCHECK,INORM,
C**         IREACT, MOLPRO below)
C**IWHICH:  Type of input potential
C**         IWHICH=0: Conventional Normal coordinate force field...Normal mode
C**         analysis and force field are assumed to be supplied by user
C**         IWHICH>0: Internal coordinate (global) potential
C**         IWHICH<0: Normal coordinate potential from fitted potential
C**                   if MOLPRO>0
C**         IWHICH<0: Normal coordinate potential from user-supplied routine
C**                   if MOLPRO=0
C**         IWHICH<0: Normal coordinate potential from interpolated potential
C**                   if MOLPRO<0
C**IDISC:   IDISC=0: Write potential and coriolis grid values to disc
C**         IDISC>0: Assumes grid values already on disc
C**         Rigid-rotor energies written to (61),(62),(63),(64) for
C**         ICOUPC=1,2,3,4
C**         Potential written to units (71),(72),(73),(74),(75),(76) for
C**         ICOUPL=1,2,3,4,5,6
C**         Coriolis (vibration) written to units (81),(82),(83),(84) for
C**         ICOUPC=1,2,3,4
C**         Coriolis (rotation) written to units (91),(92),(93),(94) for
C**         ICOUPC=1,2,3,4
C**NROTTR:  Number of rotational and translational degrees of freedom to be
C**         included in 'vibrational' modes
C**         (number of modes = NMODE=3*NATOM-6+NROTTR)
C**JMAX:    Total angular momentum quantum number
C**         JMAX>0: perform exact Whitehead-Handy if ICI < 0 (below)
C**         JMAX<0: perform adiabatic rotation if ICI < 0 (below)
C**         For VSCF only (ISCFCI = 0) or VSCF-CI (ICI > 0) perform adiabatic
C**         rotation
C**MCHECK:  MCHECK=0: Use input coordinate system as principal axis system
C**         MCHECK>0: Transform equilibrium geometry to principal axis system,
C**         and calculate rotational constants. Program continues with geometry
C**         in principal axis system (terminate if ISCFCI < 0)
C**         MCHECK<0: Finds minimum of potential (terminate if ISCFCI < 0)
C**INORM:   INORM=0: Do not perform Normal mode analysis
C**         INORM>0: Perform Normal mode analysis with internal coordinate
C**         (minimum) potential (IWHICH=1) (terminate if ISCFCI < 0)
C**         INORM<0: Perform Normal mode analysis with internal coordinate
C**         (saddle point) potential (IWHICH=1) (terminate if ISCFCI < 0)
C**
C**         Important note:  Normal modes are ordered in increasing energy.
C****************************************************************************
      READ(INP,*)NATOM,NSTAT,CONV,JCOUPL,JCOUPC,ISCFCI,IWHICH,IDISC,
     1NROTTR,KMAX,MCHECK,JNORM
      IF(JNORM.GT.0)WRITE(IOUT,351)
      IF(JNORM.LT.0)WRITE(IOUT,352)
      INORM=IABS(JNORM)
      JMAX=IABS(KMAX)
      IF(JMAX.GT.49)STOP 'ASSIGNMENTS TO MAXIMUM OF J=49'
      JCOUPS=JCOUPL
      JCOUCS=JCOUPC
      ICOUPL=IABS(JCOUPL)
      ICOUPC=IABS(JCOUPC)
      IF(JCOUPL.GE.0)JCOUPC=ICOUPC
      IF(JCOUPL.LT.0)JCOUPC=-ICOUPC
      IF(ICOUPL.GT.5.AND.IWHICH.LT.0)
     1STOP 'FITTED POTENTIAL MUST HAVE ICOUPL < 6'
      IF(ICOUPC.GT.ICOUPL)STOP 'ICOUPC CANNOT BE GREATER THAN ICOUPL'
C     TRIAT=(NATOM.EQ.3)
      IF(NSTAT.LT.0)THEN
        IORDER=-NSTAT
        NSTAT=IORDER
      END IF
C**ZERO POINT ONLY BY DEFAULT
      IF(NSTAT.EQ.0)IORDER=1
C**'TOTAL' MODES FOR NON-LINEAR MOLECULE IS 3N - 6
C**              FOR LINEAR MOLECULE IT IS 3N - 5

C**===========================================================
C**      IF(NROTTR.LT.0)NROTTR=-1
C**      NMODE=3*NATOM-6+IABS(NROTTR)
      NMODE=3*NATOM-6+NROTTR    ! Chen Qu
C**===========================================================

C**'NON-LINEAR' MODES FOR NON-LINEAR MOLECULES ALWAYS
      IF(NMODE.GT.50)STOP '>50 MODES - YOU MUST BE JOKING!'
C**TOTAL NUMBER MODES
      NONLIN=NMODE
C**'SCF' MODES ALWAYS 3N-6 IF LINEAR, OR TOTAL OTHERWISE\

C**===========================================================
C**      IF(NROTTR.LT.0)THEN
      if (nrottr .eq. -1) then
        NSMODE=3*NATOM-6
      ELSE
        NSMODE=NMODE
      END IF
C**===========================================================

C**KEEP COPY
      NNMODE=NSMODE
C**'VIBRATIONAL' MODES FOR NON-LINEAR MOLECULES ALWAYS 
C**TOTAL NUMBER MODES
      NVMODE=NMODE
C**'ACTIVE' NORMAL-MODES FOR NON-LINEAR AND NON-RPH MOLECULES ALWAYS 
C**TOTAL NUMBER MODES
      NAMODE=NMODE
      READ(INP,*)
C****************************************************************************
C**ICI,NMAX,CUT,EVLJ0,NVALR,KSTEP,IPRINT,MATSIZ,IREACT,MOLPRO,MOLINC
C-------------------------------------------------------------------
C**
C**ICI:     CI basis definition
C**         ICI>0: Number of specific VSCF states in VSCF-CI calculation
C**         ICI<0: -ICI quanta in each mode in VCI calculation (if NMAX.GE.0):
C**NMAX:    NMAX>0: Maximum sum of quanta in VCI basis
C**         NMAX=0: Unrestricted sum of quanta
C**         NMAX<0: Maximum number coupled modes in basis ... maximum -6
C**         (see MAXSUM, MAXBAS below)
C**CUT:     Cutoff energy (cm-1) for CI printout
C**EVLJ0:   Reference energy (cm-1) for J>0 CI calculations.  In general, this
C**         will be the zero point energy for J=0, which is printed at the end
C**         of the output for J=0.  This value will depend on the value of the
C**         user's potential at the minimum.  The user may wish to scale his
C**         potential such that the value at the minimum is zero, in which case
C**         the value printed for J=0 is the true zero point energy in cm-1
C**         If the reference geometry (see X0 below) is not at the potential
C**         minimum, the true zero point energy can only be determined if the
C**         potential at the reference geometry is defined to be zero, in this
C**         case the zero point energy will be the value printed for J=0,
C**         scaled by the energy difference between that at the reference
C**         geometry and that at the true minimum.
C**NVALR:   Number of rovibrational levels required (referred to J=0 stack)
C**KSTEP:   Increment between successive K (0 to 2*J+1)
C**         For KSTEP = 0, only one K-diagonal step is carried out, for Ka = 0
C**IPRINT:  Level of output
C**         IPRINT=3: Print everything
C**         IPRINT=2: Omit rotation matrix element output if J>0
C**         IPRINT=1: Omit SCF output, and K-diagonal output if J>0
C**         IPRINT=0: Omit integration grid and normal coordinate output
C**         IPRINT>0: Print full details of Reaction Path; print ONE (largest)
C**         CI coefficient
C**         IPRINT>0: Print ONE (largest) CI coefficient
C**         IPRINT<0: Print -IPRINT largest CI coefficients
C**MATSIZ:  Generate VCI or VSCF matrix size
C**         MATSIZ=0: Normal run
C**         MATSIZ>0: Calculate current size of VCI matrix and terminate
C**                   Calculate (energy-based) VSCF states and terminate
C**IREACT:  Mode number of floppy mode for reaction path hamiltonian
C**         IREACT=0: Normal run
C**         IREACT=mode: Actual (curvilinear) mode in reaction path
C**         If IREACT>0: terminate after 'Normal mode' analysis if ISCFCI < 0
C**         If IREACT<0: denotes 'abinitio' algorithm - terminate after
C**         'Normal mode' analysis if ISCFCI < 0 (OBSOLETE)
C**         If IREACT<0: use integral torsion if RPH and integrate to 2pi.
C**         (NB If IREACT>0: use half-integral torsion and integrate to 4pi)
C**         
C**MOLPRO:  Link stage with abinitio package (such as MOLPRO)
C**         MOLPRO=0: Normal run
C**         MOLPRO=1: Produce MOLPRO file for 1-dim Gauss grids using input
C**                   harmonic (omega) force field (IWHICH=0)
C**         MOLPRO=2: Fit 1-dim energies from MOLPRO=1 to 'full' or 'even'
C**                   polynomials depending on ISYM(I)
C**         MOLPRO=3: Generate 2-dim HEG grids using input fitted Gauss V(Q)
C**                   force field (IWHICH=0)
C**         MOLPRO=4: Fit 2-dim energies from MOLPRO=3
C**         MOLPRO=5: Generate 3-dim HEG grids using input fitted HEG V(Q1,Q2)
C**                   force field (IWHICH=0)
C**         MOLPRO=6: Fit 3-dim energies from MOLPRO=5
C**         MOLPRO=7: Generate 4-D HEG grids using input fitted HEG V(Q1,Q2,Q3)
C**                   force field (IWHICH=0)
C**         MOLPRO=8: Fit 4-dim energies from MOLPRO=7
C**         MOLPRO=9: Generate 5-D HEG grids using input fitted HEG
C**                   V(Q1,Q2,Q3,Q4) force field (IWHICH=0)
C**         MOLPRO=10: Fit 5-dim energies from MOLPRO=9
C**************************************************************************
C**                   To produce grids from (and fits of) global potential,
C**                   set IWHICH>0 (internal coordinates) for above steps
C**************************************************************************
C**         If MOLPRO>0: Fit potentials (terminate if ISCFCI < 0)
C**         If MOLPRO<0: Interpolate potentials (terminate if ISCFCI < 0)
C**************************************************************************
C**                   To calculate vibrational energies using fitted potential
C**                   in normal coordinates Y = gamma.Q set IWHICH<0, MOLPRO>0
C**                   (see XTANPM below)
C**                   To calculate vibrational energies using interpolated
C**                   potential in normal coordinates Q set IWHICH<0, MOLPRO<0
C**                   To calculate vibrational energies using independently
C**                   generated potential in normal coordinates Y = gamma.Q or
C**                   Q set IWHICH<0, MOLPRO=0 (see USER.DOC)
C**************************************************************************
C**MOLINC:  HEG increment used in fits (1 OR 2 = all or alternate,respectively)
C**         If MOLPRO>0: Fit potentials (terminate if ISCFCI < 0)
C**         If MOLPRO=0: Dump intrinsic potentials instead of full potentials
C****************************************************************************
      READ(INP,*)ICI,NMAX,CUT,EVLJ0,NVALR,KSTEP,JPRINT,MATSIZ,IREACT,
     1MOLPRO,MOLINC
      IF(ICI.EQ.0.AND.ISCFCI.GT.0)STOP 'ICI UNDEFINED FOR CI'
C**CURRENTLY, ONLY SINGLE RPH COORDINATE
      IFLAUD=2
      IF(IREACT.LT.0)THEN
C**INTEGER TAU THROUGHOUT
        IFLAUD=1
        IREACT=-IREACT
      END IF
C**IREACT ALWAYS POSITIVE OR ZERO.
      ABINIT=(IREACT.LT.0)
      IREACT=IABS(IREACT)
C**ENFORCE TRANSFORMATION TO Ir REPRESENTATION IF CALCULATING N. COORDS.
      IF(MCHECK.EQ.0.AND.IREACT.EQ.0.AND.JNORM.NE.0)MCHECK=1
CC    IF(MCHECK.EQ.0.AND.JNORM.NE.0)MCHECK=1
      IF(IREACT.NE.0)NREACT=1
      IF(IREACT.NE.0)THEN
        IF(ICOUPL.GT.5)STOP 'MM-RPH ONLY SUPPORTS ICOUPL=5'
        IF(ICOUPC.GT.4)STOP 'MM-RPH ONLY SUPPORTS ICOUPC=4'
      END IF
      IF(ICOUPL.GT.6)STOP 'MM-SR ONLY SUPPORTS ICOUPL=6'
      IF(ICOUPC.GT.6)STOP 'MM-SR ONLY SUPPORTS ICOUPC=6'
C*********************************************************
      IFITRP=0
      IF(IREACT.NE.0)THEN
        IF(IWHICH.EQ.0)STOP 'ILLEGAL POTENTIAL TYPE WITH RPH'
        IF(IWHICH.GT.0.AND.MOLPRO.NE.0)THEN
          WRITE(IOUT,507)
          STOP 'ILLEGAL USE OF IREACT/MOLPRO'
        END IF
        IF(IWHICH.LT.0)THEN
          IF(MOLPRO.NE.0)THEN
            WRITE(IOUT,511)
            STOP 'ILLEGAL USE OF IREACT/MOLPRO'
          ELSE
            IFITRP=IABS(ICOUPL)
            IF(MOLINC.LT.0)IFITRP=-IFITRP
          END IF
        END IF
      END IF
C*********************************************************
C**MOLINC ALWAYS POSITIVE AFTER INPUT
      MOLINC=IABS(MOLINC)
      IF(MOLINC.GT.2)MOLINC=2
C**INTRINSIC
      IF(MOLINC.LT.1)MOLINC=0
      IF(IWHICH.LT.0)THEN
        IF(MOLPRO.GT.0)THEN
C**USING MOLPRO FITTED POTENTIAL SET MOLPRO=0 (MOLINC > 0)
        ELSE
C**USING USER-SUPPLIED POTENTIAL SET MOLINC=0 (MOLPRO=0)
          IF(MOLPRO.EQ.0.AND.IFITRP.EQ.0)MOLINC=0
C**USING HERMITE INTERPOLATED POTENTIAL SET MOLPRO=0 (MOLINC < 0)
          IF(MOLPRO.NE.0)MOLINC=-MOLINC
        END IF
        MOLPRO=0
      ELSE
        IF(MOLINC.GT.0)WRITE(IOUT,515)
      END IF
      READ(INP,*)
      IF(MOLPRO.NE.0.OR.IWHICH.LT.0)THEN
C****************************************************************************
C**XTANPM(NMODE)
C---------------
C**
C**Omit data for this input if not fitting potentials (MOLPRO = 0)
C**AND not using fitted potentials (IWHICH = 0 or IWHICH > 0)
C**
C**XTANPM(I), I = 1,NMODE
C**         For fitted potentials, the coordinate tanh(gamma.Q) is used.
C**         XTANPM sets the asymptotic value of the Gauss coordinates
C**         QMAX(MODE) such that GAMMA(MODE)*QMAX(MODE) = XTANPM(MODE)
C****************************************************************************
        CALL MEMO(1,LXTANH,NSMODE,0,0,0,0,0,0,0,0)
        READ(INP,*)(W(LXTANH-1+I),I=1,NSMODE)
      END IF
      READ(INP,*)
C****************************************************************************
C**XMODQ(NMODE)
C--------------
C**
C**Scaling factor of HEG points Q for potential
C**
C**XMODQ(I),I=1,NMODE
C**     Scaling factor XMODQ(I) for mode Q(I)
C****************************************************************************
      CALL MEMO(1,LXMODQ,NSMODE,0,0,0,0,0,0,0,0)
      READ(INP,*)(W(LXMODQ-1+I),I=1,NSMODE)
C**RESET IWHICH POSITIVE (IWHICH > 0, MOLINC > 0)
      IF(IFITRP.NE.0)IWHICH=-IWHICH
      IPRINT=IABS(JPRINT)
      IF(KSTEP.LE.0)KSTEP=0
      IF(JMAX.NE.0.AND.ICI.LT.0)WRITE(IOUT,450)KSTEP
      IF(KSTEP.EQ.0)KSTEP=100000
      IF(ICI.GE.0.AND.NMAX.LT.0)NMAX=0
C**TAKE COPY OF NMAX
      NNMAX=NMAX
      WRITE(IOUT,305)JMAX
      WRITE(IOUT,205)NMODE,NSTAT,IPRINT,CONV,ICOUPL,ICOUPC
      WRITE(IOUT,210)IDISC
      IF(INORM.NE.0.AND.IREACT.EQ.0)WRITE(IOUT,350)
      IF(IREACT.GT.0)THEN
        WRITE(IOUT,465)IREACT
        IF(IFLAUD.EQ.1)WRITE(IOUT,*)'USING INTEGER TAU FOR ALL K'
        IF(IFLAUD.EQ.2)WRITE(IOUT,*)'USING HALF-INTEGER TAU FOR ODD K'
      END IF
      IF(IREACT.GT.NSMODE)STOP 'ERROR IN IREACT'
      IF(MCHECK.LT.0)WRITE(IOUT,430)
      IF(INORM.NE.0.AND.IWHICH.EQ.0)THEN
        WRITE(IOUT,355)
        STOP 'ILLEGAL USE OF INORM'
      END IF
C**FOR VCI, RESET ICI TO NUMBER OF FUNCTIONS
      IF(ICI.LT.0)ICI=ICI-1
C**RESET ICI IF TOO BIG (VSCF-CI)
      IF(ICI.GT.0.AND.ICI.GT.IABS(NSTAT))ICI=IABS(NSTAT)
C**RESET ISCFCI IF NOT NORMALS RUN OR MINIMUM CHECK
      IF(ISCFCI.LT.0.AND.INORM.EQ.0.AND.MCHECK.EQ.0.AND.IREACT.EQ.0.
     1AND.MOLPRO.EQ.0.AND.IWHICH.GE.0)ISCFCI=0
C**RESET ICI IF SCF ONLY
      IF(ISCFCI.EQ.0)THEN
        ICI=0
        WRITE(IOUT,215)
      ELSE
C**RESET ISCFCI IF TOO BIG (VSCF-CI)
        IF(ICI.GT.0.AND.ISCFCI.GT.NSTAT)ISCFCI=NSTAT
        IF(ICI.GT.0.AND.ISCFCI.GT.ICI)ISCFCI=ICI
        IF(ISCFCI.LT.0)THEN
          WRITE(IOUT,460)
        ELSE
          WRITE(IOUT,225)
          IF(JMAX.EQ.0)THEN
            WRITE(IOUT,230)ISCFCI
          ELSE
            IF(NVALR.LE.0)NVALR=ISCFCI
C**SCALE NUMBER OF FUNCTIONS CALCULATED IF J>0
            IF(TRIAT)NVALR=(JMAX+1)*NVALR
            IF(.NOT.TRIAT)NVALR=(2*JMAX+1)*NVALR
          END IF
        END IF
      END IF
      IF(IFITRP.NE.0)THEN
        WRITE(IOUT,*)'FITTING ',IABS(IFITRP),'- MODE RPH POTENTIAL'
        WRITE(IOUT,504)MOLINC
        WRITE(IOUT,475)
        WRITE(IOUT,480)(W(LXTANH-1+I),I=1,NSMODE)
      END IF
      IF(MOLPRO.NE.0)THEN
        WRITE(IOUT,470)MOLPRO
        WRITE(IOUT,504)MOLINC
        WRITE(IOUT,475)
        WRITE(IOUT,480)(W(LXTANH-1+I),I=1,NSMODE)
      END IF
      WRITE(IOUT,435)EVLJ0
C**SELECTIVE BASIS SETS (INITIALLY VCI, BUT SCFCI NEEDS NSTAT MODS)
C**J IS NUMBER OF COUPLED MODES
      J=-NMAX
C**MAXIMUM ALLOWED IS 6 (AT THE MOMENT)
      IF(J.GT.6)J=6
      IF(J.LT.0)J=0
      IF(J.EQ.0.AND.ICI.LT.0)THEN
        JCI=-ICI-1
        WRITE(IOUT,415)
        WRITE(IOUT,280)JCI,NNMAX
      END IF
      MAXJ=J
      ICONT(1)=NSMODE
      ICCONT(1)=NSMODE
      DO I=1,NSMODE
        JCONT(1,I)=I
      END DO
C**DEFAULTS
      NVAL1=ISCFCI
      NVAL2=ISCFCI
      NVALCF=ISCFCI
      READ(INP,*)
      IF(MAXJ.GT.0.AND.ISCFCI.GT.0.AND.ICI.LT.0)THEN
C****************************************************************************
C**NCONT,NVAL1,NVAL2
C-------------------
C**
C**Omit data for this input if SCF only (ISCFCI.LE.0)
C**Omit data for this input if VSCF-CI (ICI > 0)
C**Omit data for this input if unrestricted VCI (NMAX = 0)
C**Omit data for this input if restricted VCI (NMAX > 0)
C**
C**NCONT:
C**         If ICI<0: NCONT is number of contraction schemes (Maximum 2)
C**         NCONT=+2:  Complete algorithm for 2 contraction schemes
C**         NCONT=+1:  Original (single-scheme) algorithm
C**NVAL1,NVAL2:
C**         If ICI<0: NVAL1,NVAL2 are maximum number of contracted functions
C**         from schemes 1 and 2 required in complete algorithm
C**         NVAL1:  Number of functions required from contraction scheme 1
C**         NVAL2:  Number of functions required from contraction scheme 2
C**         NVALn>0:  Do all symmetries for contraction scheme 'n'
C**         NVALn<0:  Do specific symmetries for contraction scheme 'n'
C**
C****************************************************************************
        READ(INP,*)KCONT,MVAL1,MVAL2
        KCONT=IABS(KCONT)
        NVAL1=IABS(MVAL1)
        NVAL2=IABS(MVAL2)
        IF(KCONT.GT.2)KCONT=2
C**ISCFCI USED AS ALTERNATIVE IN CASE OF INPUT ERROR
        IF(NVAL1.LE.0)NVAL1=ISCFCI
        IF(NVAL2.LE.0)NVAL2=ISCFCI
      END IF
      READ(INP,*)
      IF(MAXJ.GT.0.AND.ISCFCI.GT.0.AND.ICI.LT.0)THEN
C****************************************************************************
C**ICONT(NCONT),  JCONT(NCONT,ICONT(NCONT)),  NCOUPL(NCONT),NCOUPC(NCONT)
C------------------------------------------------------------------------
C**
C**Omit data for this input if SCF only (ISCFCI = 0)
C**Omit data for this input if preliminary checks (ISCFCI < 0)
C**Omit data for this input if VSCF-CI (ICI > 0)
C**Omit data for this input if unrestricted VCI (NMAX = 0)
C**Omit data for this input if restricted VCI (NMAX > 0)
C**
C**ICONT(I), (JCONT(I,J),J=1,ICONT(I)), NCOUPL(I),NCOUPC(I)
C**
C**         Input one record for each contraction scheme I=1,NCONT
C**ICONT:   Number of modes in contraction scheme NCONT
C**JCONT:   Actual modes in contraction scheme NCONT
C**NCOUPL:  Potential coupling for scheme NCONT
C**NCOUPC:  Coriolis coupling for scheme NCONT
C****************************************************************************
        IF(KCONT.GT.0)THEN
          K=1
          DO I=1,KCONT
            READ(INP,*)ICONT(I),(JCONT(I,J),J=1,IABS(ICONT(I))),
     1      NCOUPL(I),NCOUPC(I)
C**CHECK NCOUPL,NCOUPC
            NUMBER=ICONT(I)
            JREACT=0
            MREACT=0
            DO II=1,NUMBER
              J=JCONT(I,II)
              IF(J.EQ.IREACT)JREACT=IREACT
            END DO
            IF(JREACT.NE.0)MREACT=1
C           NAMODE=NUMBER-NREACT
            NCOUPL(I)=IABS(NCOUPL(I))
            NCOUPC(I)=IABS(NCOUPC(I))
            NCOUPL(I)=MIN0(NCOUPL(I),NAMODE)
            NCOUPC(I)=MIN0(NCOUPC(I),NAMODE)
            IF(JCOUPL.GE.0)THEN
              NCOUPL(I)=IABS(NCOUPL(I))
              NCOUPC(I)=IABS(NCOUPC(I))
C             IF(NCOUPL(I).GT.IABS(ICONT(I)))NCOUPL(I)=IABS(ICONT(I))
C             IF(NCOUPC(I).GT.IABS(ICONT(I)))NCOUPC(I)=IABS(ICONT(I))
            ELSE
              NCOUPL(I)=-IABS(NCOUPL(I))
              NCOUPC(I)=-IABS(NCOUPC(I))
C             IF(-NCOUPL(I).GT.IABS(ICONT(I)))NCOUPL(I)=-IABS(ICONT(I))
C             IF(-NCOUPC(I).GT.IABS(ICONT(I)))NCOUPC(I)=-IABS(ICONT(I))
            END IF
C**END CHECK NCOUPL,NCOUPC
            IF(ICONT(I).LT.0)K=K-1
            ICCONT(I)=ICONT(I)
            ICONT(I)=IABS(ICCONT(I))
          END DO
          IF(K.LT.0)KCONT=-KCONT
        ELSE
          ICONT(1)=NSMODE
          ICCONT(1)=NSMODE
          DO I=1,NSMODE
            JCONT(1,I)=I
          END DO
        END IF
        NCONT=IABS(KCONT)
        IF(NCONT.LT.2)THEN
          NCONT=1
          KCONT=-1
          NCOUPL(1)=JCOUPL
          NCOUPC(1)=JCOUPC
        END IF
        IF(KCONT.EQ.-2.AND.MOLINC.LE.0)THEN
          NCOUPL(1)=JCOUPL 
          NCOUPC(1)=JCOUPC
          NCOUPL(2)=JCOUPL
          NCOUPC(2)=JCOUPC
        END IF
        IF(KCONT.EQ.2)THEN
          IF(ICOUPL.GT.5)THEN
            IF(JCOUPL.GT.0)JCOUPL=5
            IF(JCOUPL.LT.0)JCOUPL=-5
            WRITE(IOUT,*)
            WRITE(IOUT,*)'RECOMBINATION ICOUPL CHANGED TO 5'
            WRITE(IOUT,*)
          END IF
          ICOUPL=IABS(JCOUPL)
          IF(ICOUPC.GT.4)THEN
            IF(JCOUPC.GT.0)JCOUPC=4
            IF(JCOUPC.LT.0)JCOUPC=-4
            WRITE(IOUT,*)
            WRITE(IOUT,*)'RECOMBINATION ICOUPC CHANGED TO 4'
            WRITE(IOUT,*)
          END IF
          ICOUPC=IABS(JCOUPC)
        END IF
        WRITE(IOUT,485)NCONT,NVAL1,NVAL2
        DO I=1,NCONT
          WRITE(IOUT,490)I,(JCONT(I,J),J=1,ICONT(I))
          WRITE(IOUT,514)NCOUPL(I),NCOUPC(I)
        END DO
C**TEMPORARY
C**SAVE ORIGINAL VALUES
        ICOUPS=ICOUPL
        ICOUCS=ICOUPC
C**SET MAXIMA
        ICOUPX=MAX0(ICOUPL,IABS(NCOUPL(1)),IABS(NCOUPL(2)))
        ICOUCX=MAX0(ICOUPC,IABS(NCOUPC(1)),IABS(NCOUPC(2)))
        ICOUPL=ICOUPX
        ICOUPC=ICOUCX
C**TEMPORARY
      END IF
      CALL MEMO(1,LMXBAS,NMODE*6*2,0,0,0,0,0,0,0,0)
      READ(INP,*)
      IF(MAXJ.GT.0.AND.ISCFCI.GT.0.AND.ICI.LT.0)THEN
        WRITE(IOUT,415)
C****************************************************************************
C**MAXSUM(-NMAX,NCONT)
C---------------------
C**
C**Omit data for this input if SCF only (ISCFCI = 0)
C**Omit data for this input if preliminary checks (ISCFCI < 0)
C**Omit data for this input if VSCF-CI (ICI > 0)
C**Omit data for this input if unrestricted VCI (NMAX = 0)
C**Omit data for this input if restricted VCI (NMAX > 0)
C**
C**(MAXSUM(I,K),I=1,-NMAX)
C**
C**         Input one record for each contraction scheme K=1,NCONT
C**MAXSUM:  Maximum sum quanta in M-mode basis ... restricted only
C**         Input single record of -NMAX integers (maximum 6)
C**         for each contraction scheme K
C**
C**         MAXSUM(M,-NMAX) are defined as follows:
C**         MAXSUM(1,K) is maximum sum quanta in one-mode basis sets
C**         MAXSUM(2,K) is maximum sum quanta in two-mode basis sets
C**         MAXSUM(3,K) is maximum sum quanta in three-mode basis sets
C**         MAXSUM(4,K) is maximum sum quanta in four-mode basis sets
C**         MAXSUM(5,K) is maximum sum quanta in five-mode basis sets
C**         MAXSUM(6,K) is maximum sum quanta in six-mode basis sets
C****************************************************************************
        DO K=1,NCONT
          READ(INP,*)(MAXSUM(I,K),I=1,MAXJ)
        END DO
      END IF
      READ(INP,*)
      IF(MAXJ.GT.0.AND.ISCFCI.GT.0.AND.ICI.LT.0)THEN
C****************************************************************************
C**MAXBAS(NMODE,-NMAX)
C---------------------
C**
C**Omit data for this input if SCF only (ISCFCI = 0)
C**Omit data for this input if preliminary checks (ISCFCI < 0)
C**Omit data for this input if VSCF-CI (ICI > 0)
C**Omit data for this input if unrestricted VCI (NMAX = 0)
C**Omit data for this input if restricted VCI (NMAX > 0)
C**
C**(MAXBAS(M,I),M=1,NMODE)
C**
C**         Input one record for each I-mode basis I=1,-NMAX
C**MAXBAS:  Defines selective quanta in VCI basis
C**         Input -NMAX records (maximum 6) of NMODE integers
C**
C**         MAXBAS(M,1) is maximum quanta for mode M in one-mode basis sets
C**         MAXBAS(M,2) is maximum quanta for mode M in two-mode basis sets
C**         MAXBAS(M,3) is maximum quanta for mode M in three-mode basis sets
C**         MAXBAS(M,4) is maximum quanta for mode M in four-mode basis sets
C**         MAXBAS(M,5) is maximum quanta for mode M in five-mode basis sets
C**         MAXBAS(M,6) is maximum quanta for mode M in six-mode basis sets
C****************************************************************************
        IF(ICI.LT.0)WRITE(IOUT,440)
        CALL MAXIN(W(LMXBAS),NMODE,MAXJ,ICI)
        IF(ICI.LT.0)THEN
C**NBAS:    Maximum no. quanta in MAXJ-Mode basis
          DO K=1,NCONT
            WRITE(IOUT,420)K,(NBAS(I,K),I=1,MAXJ)
            WRITE(IOUT,425)K,(MAXSUM(I,K),I=1,MAXJ)
          END DO
        END IF
      END IF
C**KEEP COPY
      IF(NMAX.LT.0.AND.ICI.LT.0)ICI=NMAX
      CALL MEMO(1,LMODNT,NMODE,0,0,0,0,0,0,0,0)
      READ(INP,*)
C****************************************************************************
C**MODINT(NMODE)
C---------------
C**
C**(MODINT(M),M=1,NMODE)
C**
C**MODINT:    Integration factor (1 or 2) for NMODE modes
C**MODINT=1:  Mode is totally symmetric (A1g)
C**MODINT=2:  Mode is antisymmetric (B2g,B1g,A2g,A1u,B2u,B1u,A2u)
C**MODINT>2:  Torsional mode is antisymmetric and integration is to 4pi
C****************************************************************************
      CALL INMDIN(W(LMODNT),NMODE,0)
      CALL MEMO(1,LMEFF,2*NMODE,0,0,0,0,0,0,0,0)
      READ(INP,*)
C****************************************************************************
C**MEFF(NMODE)
C-------------
C**
C**(MEFF(M),M=1,NMODE)
C**
C**MEFF:      Corresponding mode forming effective 1-dim potential for
C**           current mode.
C**MEFF(M)=N  An effective 1-dim potential for mode M will be generated in
C**           order to define the 1-dim basis.  The effective potential will
C**           be with respect to mode N
C****************************************************************************
      CALL INMDIN(W(LMEFF),NMODE,1)
      CALL INMDIN(W(LMEFF+NMODE),NMODE,-1)
      READ(INP,*)
C****************************************************************************
C**NCYCLE,TOLLAN,LANMAX,IGIV,LAN20,LANCUT,TOLCUT,LANCYC,TOLMAT,LANINC
C--------------------------------------------------------------------
C**
C**NCYCLE:  Number of Lanczos cycles (Lanczos ignored if NCYCLE=0)
C**         NCYCLE>0: (use sparse disc usage always - slower if PT)
C**         NCYCLE<0: (use dense disc if 1 & 2-dim PT - faster if PT)
C**TOLLAN:  Tolerance for Lanczos eigenvalues
C**LANMAX:  Maximum order of half-matrix used to build Lanczos matrix
C**LGIV:    LGIV=0: Use QL algorithm in Lanczos
C**         LGIV>0: Use GIVENS in Lanczos
C**LAN20:   Size below which uses single disc file fort.20
C**         If size exceeds LAN20, disc output spread over 5 discs to
C**         avoid Linux 2GB problem
C**         LAN20>0: create all Lanczos files and continue
C**         LAN20<0: restart Lanczos diagonalisation using previous files
C**
C**         NB Restart for single VCI symmetry block only.....
C**         .....job with LAN20 > 0 can be `killed' once message:
C**         `Calculating Lanczos' is reached.
C**LANCUT:  LANCUT>0: number of 1 & 2 dim functions to keep
C**         LANCUT<0: stop after set-up (single symmetry only)
C**TOLCUT:  Tolerance for perturbation test
C**         TOLCUT>0: 1 & 2 dim perturbaton basis
C**         TOLCUT<0: 1,2 & 3 dim perturbation basis
C**LANCYC:  Size of inner matrix (cf ISCFCI)
C**TOLMAT:  Tolerance below which matrix element = 0 exactly
C**LANINC:  No. of increments for final matrix C' if PT 
C****************************************************************************
      LANCUT=0
      LLNCUT=0
      TOLCUT=0
      TLLCUT=0
      READ(INP,*)NCYCLE,TOLLAN,LANMAX,IGIV,LLAN20,LLNCUT,TLLCUT,LANCYC,
     1TOLMAT,LANINC
      LANZA=(NCYCLE.GT.0)
      LANZB=(NCYCLE.LT.0)
      NCYCLE=IABS(NCYCLE)
      IF(ICI.LT.0)WRITE(IOUT,360)NCYCLE
      LANCZ=(NCYCLE.GT.0)
      IF(LANCZ)THEN
        LAN20=IABS(LLAN20)
        IF(ICI.LT.0)THEN
          WRITE(IOUT,365)TOLLAN
          WRITE(IOUT,370)LANMAX
          IF(LANCYC.GT.ISCFCI)LANCYC=ISCFCI
          IF(LANCYC.LE.0)LANCYC=ISCFCI
          WRITE(IOUT,371)LANCYC,TOLMAT
        END IF
        LANCUT=IABS(LLNCUT)
        TOLCUT=DABS(TLLCUT)
        IF(LANINC.LT.1)LANINC=1
        IF(TLLCUT.EQ.0)THEN
          TOLCUT=0
          LANCUT=0
          LLNCUT=0
        END IF
        IF(LANCUT.GT.0)WRITE(IOUT,330)LANCUT,TOLCUT,LANINC
      ELSE
        TLLCUT=0
        LLNCUT=0
        TOLCUT=0
        LANCUT=0
      END IF
      READ(INP,*)
C****************************************************************************
C**NRSYM,NWSYM,ISYMPG
C--------------------
C**
C**NRSYM:   Number of rovibrational symmetry species (1 or 2 or 4 or 8)
C**NWSYM:   Number of mode symmetry species
C**         For no symmetry NRSYM=NVSYM=NWSYM=1 (C1)...can be used for 
C**         molecules of any symmetry and size
C**         For triatomic Cs molecules NRSYM=2 (A',A"), NVSYM=NWSYM=1 (A')
C**         For four-atom or larger Cs molecules NRSYM=NVSYM=NWSYM=2 (A',A")
C**         For four-atom or larger C2 molecules NRSYM=NVSYM=NWSYM=2 (A,B)
C**         For triatomic C2v molecules NRSYM=4 (A1,B2,B1,A2), 
C**         NVSYM=NWSYM=2 (A1,B2)
C**         For four-atom C2v molecules NRSYM=NVSYM=4 (A1,B2,B1,A2), 
C**         NWSYM=3 (A1,B2,B1) if out-of-plane
C**         NWSYM=3 (A1,B2,A2) if torsion
C**         For larger C2v molecules NRSYM=NVSYM=4 (A1,B2,B1,A2), 
C**         NWSYM=4 (A1,B2,B1,A2)
C**         For larger molecules, where reduction of symmetry to C2v is 
C**         possible NRSYM=NVSYM=4 (A1,B2,B1,A2), NWSYM=4 (A1,B2,B1,A2)
C**         Where there is a centre of symmetry, C2v is converted to pseudo D2h
C**         NRSYM=NVSYM=8 (A1g,B2g,B1g,A2g; A1u,B2u,B1u,A2u), NWSYM=1 -> 8
C**ISYMPG:  Point Group indicator (default setting ISYMPG=0), including:
C**         C1, Cs or C2, C2v or D2 or C2h, D2h
C****************************************************************************
      READ(INP,*)NRSYM,NWSYM,ISYMPG
      NVSYM=NRSYM
      IF(TRIAT.AND.NRSYM.GT.1)NVSYM=NRSYM/2
      WRITE(IOUT,295)NRSYM,NVSYM,NWSYM
      READ(INP,*)
C****************************************************************************
C**NSYM(NWSYM),ISYM(NWSYM,NSYM(NWSYM))    (order of NWSYM:  A1,B2,B1,A2)
C-----------------------------------------------------------------------
C**
C**NSYM(I), (ISYM(I,J),J=1,NSYM(I))
C**
C**         Input one record for each mode symmetry I=1,NWSYM
C**         For triatomic C2v molecules, NWSYM=1 corresponds to A1,
C**         NWSYM=2 corresponds to B2
C**         For four-atom C2V molecules, NWSYM=1 corresponds to A1,
C**         NWSYM=2 corresponds to B2, NWSYM=3 corresponds to B1
C**         etc.
C**         For D2h molecules, NWSYM=1 -> 4 are 'g'; NWSYM=5 -> 8 are 'u'
C**         The user must know the symmetry of the Normal modes which are
C**         ordered either in increasing energy (see INORM=1),
C**         or according to the user's convention (see INORM=0)
C**
C**NSYM:    Number of modes of current (NWSYM) symmetry
C**ISYM:    NSYM specific modes of current (NWSYM) symmetry
C**
C**         For example, in the case of H2O, there are 2 mode vibrations
C**         of A1 symmetry (modes 1 and 2, say) corresponding to NWSYM=1.
C**         There is a single mode vibration of B2 symmetry (mode 3, say)
C**         corresponding to NWSYM=2.
C**         Hence for NWSYM=1 (totally symmetric), the input required is:
C**         NSYM=2,  ISYM(I)=1 2
C**         For I=2, the input required is:
C**         NSYM=1,  ISYM(I)=3
C****************************************************************************
      DO I=1,NWSYM
        READ(INP,*)NSYM(I),(ISYM(I,J),J=1,NSYM(I))
        DO J=1,NSYM(I)
          IF(ISYM(I,J).GT.NVMODE.OR.ISYM(I,J).LE.0)STOP 'ERROR IN ISYM'
        END DO
        WRITE(IOUT,300)I,(ISYM(I,J),J=1,NSYM(I))
      END DO
      READ(INP,*)
C****************************************************************************
C**MVSYM,MWSYM(MVSYM)     (order of MWSYM:  A1,B2,B1,A2)
C-------------------------------------------------------
C**
C**MVSYM,(MWSYM(I),I=1,MVSYM)
C**
C**MVSYM:   Single integer corresponding to the number of actual vibrational
C**         (JMAX=0) or rovibrational (JMAX>0) symmetry species required
C**MWSYM:   MVSYM integers corresponding to the Specific vibrational (JMAX=0)
C**         or rovibrational (JMAX>0) symmetry species required
C**
C**         For J = 0 get vibrational energies of specified symmetry.
C**         For J > 0 all vibrational symmetries are included in K-blocks, but
C**         final rovibrational calculations are done for symmetries specified.
C**
C**Note:    For a triatomic molecule, the number of possible rovibrational
C**         species is twice the number of possible vibrational species
C****************************************************************************
      READ(INP,*)MVSYM,(MWSYM(I),I=1,MVSYM)
      WRITE(IOUT,455)MVSYM,(MWSYM(I),I=1,MVSYM)
      IF(JMAX.EQ.0.AND.MVSYM.GT.NVSYM)MVSYM=NVSYM
      IF(JMAX.GT.0.AND.MVSYM.GT.NRSYM.AND..NOT.TRIAT)MVSYM=NRSYM
      IF(JMAX.GT.0.AND.MVSYM.GT.2*NVSYM.AND.TRIAT)MVSYM=2*NVSYM
      DO I=1,MVSYM
        IF(JMAX.EQ.0.AND.MWSYM(I).GT.NVSYM.OR.MWSYM(I).LT.0)MWSYM(I)=0
        IF(JMAX.GT.0.AND..NOT.TRIAT.AND.(MWSYM(I).GT.NRSYM.OR.
     1  MWSYM(I).LT.0))MWSYM(I)=0
        IF(JMAX.GT.0.AND.TRIAT.AND.(MWSYM(I).GT.2*NVSYM.OR.
     1  MWSYM(I).LT.0))MWSYM(I)=0
      END DO
      READ(INP,*)
C****************************************************************************
C**MXDIP(3)     (Irreducible representations of dipole components)
C-----------------------------------------------------------------
C**
C**(MXDIP(I),I=1,3)
C**
C**MXDIP(1): I.R. of molecule-fixed X component of dipole
C**MXDIP(2): I.R. of molecule-fixed Y component of dipole
C**MXDIP(3): I.R. of molecule-fixed Z component of dipole
C**
C**Note:     I.R. are given according to standard MM convention
C**          C1:  A=1; C2: A=1,B=2; Cs: A'=1,A"=2; C2v: A1=1,B2=2,B1=3,A2=4
C**          D2:  A=1,B3=2,B2=3,B1=4; C2h: Ag=1,Bg=2,Au=3,Bu=4
C**          D2h: A1g=1,B2g=1,B1g=3,A2g=4,A1u=5,B2u=6,B1u=7,A2u=8
C****************************************************************************
      READ(INP,*)(MXDIP(I),I=1,3)
      READ(INP,*)
C****************************************************************************
C**MX(3)     (Relationships between symmetry axes and principle axes)
C--------------------------------------------------------------
C**
C**(MX(I),I=1,3)
C**
C**MX(1): Symmetry axis (x,y,z) corresponding to principal X-axis
C**MX(2): Symmetry axis (x,y,z) corresponding to principal Y-axis
C**MX(3): Symmetry axis (x,y,z) corresponding to principal Z-axis
C**
C**Note:     Throughout, we will refer to principle axes as X,Y,Z
C**                      we will refer to symmetry axes as x,y,z
C****************************************************************************
      READ(INP,*)(MX(I),I=1,3)
      READ(INP,*)
C****************************************************************************
C**MXROT(3)     (Irreducible representations of Rx, Ry, Rz)
C**             (Some Group Tables refer to Jx, Jy, Jz)
C--------------------------------------------------------------
C**
C**(MXROT(I),I=1,3)
C**
C**MXROT(1): I.R. of rotation about symmetry x-axis (Rx)
C**MXROT(2): I.R. of rotation about symmetry y-axis (Ry)
C**MXROT(3): I.R. of rotation about symmetry z-axis (Rz)
C**
C**Note:     I.R. are given according to standard MM convention
C**          C1:  A=1; C2: A=1,B=2; Cs: A'=1,A"=2; C2v: A1=1,B2=2,B1=3,A2=4
C**          D2:  A=1,B3=2,B2=3,B1=4
C**          C2h: A1g=1,B2g=1,B1g=3,A2g=4,A1u=5,B2u=6,B1u=7,A2u=8
C****************************************************************************
      READ(INP,*)(MXROT(I),I=1,3)
      READ(INP,*)
C****************************************************************************
C**
C**ISYMT
C**
C**ISYMT:   I.R. of sin(M.tau/2).   For use only with MM-RPH.
C**         ISYMT=1 if MM-SR
C****************************************************************************
      READ(INP,*)ISYMT
      IF(IREACT.EQ.0)ISYMT=1
C**NOW DETERMINE I.R. OF FIRST 4 ROTATIONAL FUNCTIONS
      CALL SORTIT
      DO I=1,MVSYM
        JDUMP(I)=0
      END DO
      KDUMP=0
      READ(INP,*)
      IF(ICI.LT.0.AND.ISCFCI.GT.0)THEN
C****************************************************************************
C**LDUMP,MDUMP,IDUMP
C-------------------
C**
C**Omit data for this input if SCF only (ISCFCI = 0)
C**Omit data for this input if preliminary checks (ISCFCI < 0)
C**Omit data for this input if VSCF-CI (ICI > 0)
C**
C**LDUMP:  Restart parameter for property fuunctions
C**        Currently the following seven properties are supported:
C**        LDUMP=1: input dipole functions and evaluate J=0 matrix elements
C**        LDUMP=2: expectation values of Q
C**        LDUMP=3: expectation values of Q^2
C**        LDUMP=4: wavefunction plots
C**        LDUMP=5: user-defined property
C**        LDUMP=6: F-C factors (assumes 2 electronic states)
C**        LDUMP=7: P,Q,R-branch intensities
C**        LDUMP=8: Raman Spectrum
C**        LDUMP=9: Partition Functions
C**        LDUMP<0: Redraw spectrum using previous files...
C**        ...LDUMP=-7 (IR), LDUMP=-8 (Raman)
C**MDUMP:  Number of (central) grid points used to generate Eckart frame
C**        MDUMP>0: program continues with VCI
C**        MDUMP<0: program halts after generating required points
C**        (obsolete)
C**IDUMP:  Number of J=0 VCI functions written to disc (ICI < 0)
C**        LDUMP=0: IDUMP consecutive functions for each symmetry
C**OBSOLETE
C********* LDUMP<0: NDUMP(IABS(IDUMP)) specific functions for each symmetry
C**OBSOLETE
C****************************************************************************
        READ(INP,*)LDUMP,KDUMP,(JDUMP(I),I=1,MVSYM)
      END IF
      IF((LDUMP.EQ.1.OR.IABS(LDUMP).EQ.7).AND.JMAX.EQ.0)THEN
        ISET=0
        DO I=1,MVSYM
          IF(JDUMP(I).NE.0)ISET=1
        END DO
        IF(ISET.EQ.0)THEN
          DO I=1,MVSYM
            JDUMP(I)=ISCFCI
          END DO
        END IF
      END IF
C**SCALE NUMBER FUNCTIONS DUMPED IF J>0
      DO I=1,MVSYM
        IF(NATOM.EQ.3)JDUMP(I)=(JMAX+1)*JDUMP(I)
        IF(.NOT.TRIAT)JDUMP(I)=(2*JMAX+1)*JDUMP(I)
      END DO
      IDUMP=0
C**GET MAXIMUM NO. DUMPED FUNCTIONS
      DO I=1,MVSYM
        IF(JDUMP(I).GT.IDUMP)IDUMP=JDUMP(I)
      END DO
      MDUMP=IABS(KDUMP)
C     IF((IDUMP.NE.0.OR.LDUMP.GT.0).AND.ICOUPL.GT.4)THEN
C       WRITE(IOUT,508)
C       STOP 'ILLEGAL USE OF DUMP/RESTART'
C     END IF
      IF((IDUMP.NE.0.OR.LDUMP.GT.0).AND.KCONT.GT.0)THEN
        WRITE(IOUT,509)
        STOP 'ILLEGAL USE OF DUMP/RESTART'
      END IF
C     IF((IDUMP.NE.0.OR.LDUMP.GT.0).AND.IREACT.NE.0)THEN
C       WRITE(IOUT,510)
C       STOP 'ILLEGAL USE OF DUMP/RESTART'
C     END IF
C**CAN'T BE GREATER THAN MAXIMUM ASKED FOR
      DO I=1,MVSYM
        IF(KMAX.LE.0)THEN
C**J=0
          IF(JDUMP(I).GT.ISCFCI)JDUMP(I)=ISCFCI
        ELSE
C**J>0
          IF(JDUMP(I).GT.NVALR)JDUMP(I)=NVALR
        END IF
      END DO
C**OBSOLETE
CC    IF(LDUMP.LT.0)THEN
CC      READ(INP,*)
C****************************************************************************
C**NDUMP(JDUMP(MVSYM),MVSYM)
C---------------------
C**
C**(NDUMP(I,MVSYM),I=1,JDUMP(MVSYM))
C**
C**        Input only if LDUMP < 0
C**        Input MVSYM records of JDUMP(MVSYM) function numbers
C**NDUMP:  JDUMP(MVSYM) specific functions for symmetry species MVSYM to 
C**        dump to disc
C**        Function numbers are those given in the J=0 output on fort.2
C****************************************************************************
CC      CALL MEMO(1,LNDUMP,MVSYM*IDUMP,0,0,0,0,0,0,0,0)
C**ANY OUT-OF-RANGE FUNCTION WILL BE GIVEN VALUE OF 'NEGATIVE'
CC      CALL INDUMP(W(LNDUMP),MVSYM,IDUMP,ISCFCI,JDUMP)
CC    END IF
C**OBSOLETE
      IF(IABS(LDUMP).EQ.7.AND.JMAX.EQ.0)LDUMP=1
      IF(IDUMP.NE.0.AND.(LDUMP.NE.0.AND.LDUMP.NE.1))LDUMP=0
C**READ REMAINING STANDARD INPUT
      KNBF=NMODE
      KMBF=NMODE
      KXM=NATOM
      KX0=NATOM*3
      KXL=NATOM*NMODE*3
CCCC  IF(IREACT.GT.0)KXL=NATOM*NMODE*3*362*2
      IF(IREACT.GT.0)KXL=NATOM*NMODE*3*722*2
      CALL MEMO(5,LNBF,KNBF,LMBF,KMBF,LXM,KXM,LX0,KX0,LXL,KXL)
      CALL MEMO(4,LYL,KXL,LY0,KX0,LNVF,KNBF,LMVF,KMBF,0,0)
C**MAXIMUM OF 6 MODES COUPLED!!
C**RESET NPOT TO NMODE IF INTERNAL FORCE FIELD
      IF(IWHICH.NE.0)NPOT=NMODE
      KOMEGA=NMODE
CCCC  IF(IREACT.GT.0)KOMEGA=NMODE*362
      IF(IREACT.GT.0)KOMEGA=NMODE*722
      KISTAT=NSTAT*NMODE
      CALL MEMO(2,LOMEGA,KOMEGA,LISTAT,KISTAT,0,0,0,0,0,0)
      CALL MEMO(4,LJSTAT,KISTAT,LESTAT,NSTAT,LKSTAT,NMODE,LWSTAT,NSTAT,
     10,0)
      KXX=NATOM*3
      KRR=NATOM*NATOM
      IF(IREACT.GT.0)THEN
CCCC    KXX=NATOM*3*362*2
CCCC    KRR=NATOM*NATOM*364*9
        KXX=NATOM*3*722*2
        KRR=NATOM*NATOM*724*9
      END IF
      CALL MEMO(3,LXX,KXX,LRR,KRR,LTEMP5,NSMODE,0,0,0,0)
C**RESET ICI IF TOO BIG (VIRTUALS)

C**=======================================================================
      if (nrottr .eq. -1 .and. .not. triat) ireact=-nsmode
C**      IF(NROTTR.LT.0.AND..NOT.TRIAT)IREACT=-NSMODE
      CALL INPUT(NMODE,W(LOMEGA),W(LXX)
     1,W(LNBF),W(LMBF),NSTAT,W(LISTAT),NATOM,W(LXM),W(LX0),W(LXL),
     2W(LRR),ICI,W(LMXBAS),MAXJ,W(LJSTAT),W(LKSTAT),W(LESTAT),
     3W(LWSTAT),W(LNVF),W(LMVF),W(LTEMP5),INORM,MCHECK,NSMODE,NVMODE,0,
     4W)
C**=======================================================================

C**'ACTIVE' NORMAL-MODES FOR RPH MOLECULES 
C**(3N - 6 - 1)
      IF(IREACT.GT.0)NAMODE=NMODE-NREACT

C     IF(LINEAR)THEN
C**BATH
C**NOT ALWAYS TREATED AS LINEAR IF NROTTR > 0
C       IF(NROTTR.GT.0)LINBND=0
C**BATH
C**'ACTIVE' NORMAL-MODES FOR LINEAR MOLECULES 
C**(3N - 5 MINUS NUMBER OF LINEAR BENDS)
C       NAMODE=NMODE-LINBND
C**'VIBRATIONAL' NORMAL-MODES FOR LINEAR MOLECULE 
C**(3N - 5 MINUS NUMBER OF LINEAR BENDS)
C       NVMODE=NMODE-LINBND
C     END IF

      IF(ICI.GT.0)WRITE(IOUT,275)ICI
      CALL FLUSH(IOUT)

C***************************************
      IF((IWHICH.LT.0.AND.MOLINC.NE.0).OR.MOLPRO.NE.0.OR.IFITRP.NE.0)
     1  THEN
        IF(ICOUPL.GT.1.OR.IABS(MOLPRO).GT.4)
     1  CALL MEMO(1,LINDK,NAMODE,0,0,0,0,0,0,0,0)
        IF(ICOUPL.GT.2.OR.IABS(MOLPRO).GT.6)
     1  CALL MEMO(1,LINDL,NAMODE,0,0,0,0,0,0,0,0)
        IF(ICOUPL.GT.3.OR.IABS(MOLPRO).GT.8)
     1  CALL MEMO(1,LINDN,NAMODE,0,0,0,0,0,0,0,0)
        IF(ICOUPL.GT.4.OR.IABS(MOLPRO).GT.10)
     1  CALL MEMO(1,LINDM,NAMODE,0,0,0,0,0,0,0,0)
        NTOT1=NSMODE
        CALL GETIND(W(LINDK),W(LINDL),W(LINDN),W(LINDM),NTOT1,
     1  NTOT2,NTOT3,NTOT4,NTOT5)
        IF(ICOUPL.GT.0.OR.IABS(MOLPRO).GT.2)
     1  CALL MEMO(1,LNP1,NTOT1,0,0,0,0,0,0,0,0)
        IF(ICOUPL.GT.1.OR.IABS(MOLPRO).GT.4)
     1  CALL MEMO(1,LNP2,NTOT2,0,0,0,0,0,0,0,0)
        IF(ICOUPL.GT.2.OR.IABS(MOLPRO).GT.6)
     1  CALL MEMO(1,LNP3,NTOT3,0,0,0,0,0,0,0,0)
        IF(ICOUPL.GT.3.OR.IABS(MOLPRO).GT.8)
     1  CALL MEMO(1,LNP4,NTOT4,0,0,0,0,0,0,0,0)
        IF(ICOUPL.GT.4.OR.IABS(MOLPRO).GT.10)
     1  CALL MEMO(1,LNP5,NTOT5,0,0,0,0,0,0,0,0)
      END IF
      READ(INP,*)
C****************************************************************************
C**MAXCHK(ICOUPL)
C---------------
C**
C**(MAXCHK(M),M=1,ICOUPL)
C**
C**MAXCHK:    Inspection of M-mode grids for spurious maxima
C**           MAXCHK(M)=0:  Do nothing
C**           MAXCHK(M)=+N:  Report maxima found for mode N at M-mode
C**           MAXCHK(M)=-N:  Attempt correction for mode IABS(N) at M-mode
C**
C****************************************************************************
      READ(INP,*)(MAXCHK(I),I=1,ICOUPL)
      CALL FLUSH(IOUT)
      READ(INP,*)
C**********************************************************************
C**USER-DEFINED POTENTIAL
C**At this point, control passes to the user to input any data he requires in
C**the routine USERIN (see Manual)
C**********************************************************************
      WRITE(IOUT,*)
      WRITE(IOUT,*)'***********************START USER PARAMETERS'
      CALL USERIN
      IF(INORM.EQ.0)THEN
        CALL GETPOT(V0,NATOM,W(LX0),W(LRR))
        WRITE(IOUT,*)
        WRITE(IOUT,*)'V AT EQU ',V0*WAVENM
        WRITE(*,*)'V AT EQU ',V0*WAVENM
        WRITE(IOUT,*)
      END IF

C**HERMITE INTERPOLATION
      IF((IWHICH.LT.0.AND.MOLINC.LT.0).OR.MOLPRO.LT.0)THEN
        IENTER=0
2222    CONTINUE
        MMX1=IENTMX(1)
        MAX2=IENTMX(2)
        MAX3=IENTMX(3)
        MAX4=IENTMX(4)
        MAX5=IENTMX(5)
        CALL HERMIN(NSMODE,W(LXTANH),W(LNP1),W(LCP1),W(LVP1),W(LDP1),
     1  NTOT1,MMX1,W(LNP2),W(LCP2),W(LVP2),W(LDP2A),W(LDP2B),NTOT2,
     2  MAX2,W(LNP3),W(LCP3),W(LVP3),W(LDP3A),W(LDP3B),W(LDP3C),NTOT3,
     3  MAX3,W(LNP4),W(LCP4),W(LVP4),W(LDP4A),W(LDP4B),W(LDP4C),
     4  W(LDP4D),NTOT4,MAX4,W(LINDK),W(LINDL),W(LINDN),W(LINDM))
        IF(IENTER.EQ.1)
     1  CALL MEMO(3,LCP1,IENTMX(1)*NTOT1,LVP1,IENTMX(1)*NTOT1,
     2  LDP1,IENTMX(1)*NTOT1,0,0,0,0)
        IF(IENTER.EQ.2)
     1  CALL MEMO(4,LCP2,IENTMX(2)*NTOT1,
     2  LVP2,IENTMX(2)*IENTMX(2)*NTOT2,LDP2A,IENTMX(2)*IENTMX(2)*NTOT2,
     3  LDP2B,IENTMX(2)*IENTMX(2)*NTOT2,0,0)
        IF(IENTER.EQ.4)
     1  CALL MEMO(5,LCP3,IENTMX(3)*NTOT1,
     2  LVP3,IENTMX(3)*IENTMX(3)*IENTMX(3)*NTOT3,
     3  LDP3A,IENTMX(3)*IENTMX(3)*IENTMX(3)*NTOT3,
     4  LDP3B,IENTMX(3)*IENTMX(3)*IENTMX(3)*NTOT3,
     5  LDP3C,IENTMX(3)*IENTMX(3)*IENTMX(3)*NTOT3)
        IF(IENTER.EQ.6)
     1  CALL MEMO(5,LCP4,IENTMX(4)*NTOT1,
     2  LVP4,IENTMX(4)*IENTMX(4)*IENTMX(4)*IENTMX(4)*NTOT4,
     3  LDP4A,IENTMX(4)*IENTMX(4)*IENTMX(4)*IENTMX(4)*NTOT4,
     4  LDP4B,IENTMX(4)*IENTMX(4)*IENTMX(4)*IENTMX(4)*NTOT4,
     5  LDP4C,IENTMX(4)*IENTMX(4)*IENTMX(4)*IENTMX(4)*NTOT4)
        IF(IENTER.EQ.6)
     1  CALL MEMO(1,LDP4D,IENTMX(4)*IENTMX(4)*IENTMX(4)*IENTMX(4)*
     2  NTOT4,0,0,0,0,0,0,0,0)
        IF(IENTER.NE.0)GO TO 2222
      END IF
C**MOLPRO FITS
      IF((IWHICH.LT.0.AND.MOLINC.GT.0).OR.MOLPRO.GT.0)THEN
        IENTER=0
3333    CONTINUE
        MMX1=IENTMX(1)
        MAX2=IENTMX(2)
        MAX3=IENTMX(3)
        MAX4=IENTMX(4)
        MAX5=IENTMX(5)
        CALL FITSIN(NSMODE,W(LXTANH),W(LNP1),W(LCP1),W(LMP1),NTOT1,
     1  MMX1,W(LNP2),W(LCP2),W(LMP2),NTOT2,MAX2,W(LNP3),
     2  W(LCP3),W(LMP3),NTOT3,MAX3,W(LNP4),W(LCP4),W(LMP4),NTOT4,
     3  MAX4,W(LNP5),W(LCP5),W(LMP5),NTOT5,MAX5,
     3  W(LINDK),W(LINDL),W(LINDN),W(LINDM))
        IF(IENTER.EQ.1)
     1  CALL MEMO(2,LCP1,IENTMX(1)*NTOT1,LMP1,
     2  IENTMX(1)*NTOT1,0,0,0,0,0,0)
        IF(IENTER.EQ.3)
     1  CALL MEMO(2,LCP2,IENTMX(2)*NTOT2,LMP2,
     2  IENTMX(2)*NTOT2*2,0,0,0,0,0,0)
        IF(IENTER.EQ.5)
     1  CALL MEMO(2,LCP3,IENTMX(3)*NTOT3,LMP3,
     2  IENTMX(3)*NTOT3*3,0,0,0,0,0,0)
        IF(IENTER.EQ.7)
     1  CALL MEMO(2,LCP4,IENTMX(4)*NTOT4,LMP4,
     2  IENTMX(4)*NTOT4*4,0,0,0,0,0,0)
        IF(IENTER.EQ.9)
     1  CALL MEMO(2,LCP5,IENTMX(5)*NTOT5,
     2  LMP5,IENTMX(5)*NTOT5*5,0,0,0,0,0,0)
        IF(IENTER.NE.0)GO TO 3333
        MAXPOW=MAXPOW+1
      END IF
      WRITE(IOUT,*)
      WRITE(IOUT,*)'***********************END USER PARAMETERS'
      CALL FLUSH(IOUT)
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                                  END OF INPUT STAGE    
C**********************************************************************
C**********************************************************************
C**********************************************************************
      KXZ=NMODE*NMODE*3
      IF(IREACT.GT.0)THEN
CCCC    KXZ=NMODE*NMODE*3*362
        KXZ=NMODE*NMODE*3*722
      END IF
      KYZ=NMODE*NMODE*3
      CALL MEMO(2,LYZ,KYZ,LXZ,KXZ,0,0,0,0,0,0)
      IF(INORM.NE.0)THEN
        CALL GETNOR(NATOM,NMODE,W(LXM),W(LX0),W(LOMEGA),
     1  W(LXL),W(LXX),W(LRR),W(LIPOT),W(LJPOT),W(LCPOT),NPOT)
      END IF

C**TEMPORARY
CCCC  CALL NEWMIN(W(LX0),W(LXL),W(LXM),NATOM,NMODE)
C**TEMPORARY


      IF(LINEAR)THEN
C**BATH
C**NOT ALWAYS TREATED AS LINEAR IF NROTTR > 0
        IF(NROTTR.GT.0)LINBND=0
C**BATH
C**'ACTIVE' NORMAL-MODES FOR LINEAR MOLECULES 
C**(3N - 5 MINUS NUMBER OF LINEAR BENDS)
        NAMODE=NMODE-LINBND
C**'VIBRATIONAL' NORMAL-MODES FOR LINEAR MOLECULE 
C**(3N - 5 MINUS NUMBER OF LINEAR BENDS)
        NVMODE=NMODE-LINBND
      END IF

      CALL FLUSH(IOUT)
      IF(LINEAR)THEN
C**MOVE ALL LINEAR PARAMETERS TO TOP
        CALL MEMO(1,LTEMP,NATOM*3,0,0,0,0,0,0,0,0)
C**BATH
C**NOT ALWAYS TREATED AS LINEAR IF NROTTR > 0
        IF(NROTTR.LT.0)
     1  CALL BLIP(NATOM,NMODE,W(LNBF),W(LMBF),W(LNVF),W(LOMEGA),W(LXL),
     2  W(LMODNT),W(LMEFF),W(LTEMP),W(LMXBAS),NSMODE,MAXJ,ICI,1)
C**BATH
        CALL MEMO(-1,LTEMP,NATOM*3,0,0,0,0,0,0,0,0)
      END IF
      DO ICPL=1,ICOUPL
        FACTOR(ICPL)=1.D0
        DO I=2,NAMODE
          FACTOR(ICPL)=FACTOR(ICPL)*I
        END DO
        DO I=1,ICPL
          FACTOR(ICPL)=FACTOR(ICPL)/I
        END DO
        IDENOM=NAMODE-ICPL
        DO I=1,IDENOM
          FACTOR(ICPL)=FACTOR(ICPL)/I
        END DO
        IF(IREACT.GT.0)FACTOR(ICPL)=FACTOR(ICPL)+1
      END DO
      IF(IWHICH.EQ.0)THEN
        WRITE(IOUT,265)
      ELSE
        WRITE(IOUT,270)(FACTOR(ICPL),ICPL=1,ICOUPL)
      END IF
      CALL FLUSH(IOUT)
      IF(INORM.NE.0)THEN
        IF(ISCFCI.LT.0.AND.MCHECK.EQ.0.AND.IREACT.EQ.0.AND.MOLPRO.EQ.0)
     1  THEN 
          WRITE(6,*) 'NORMALS RUN ONLY'
          RETURN
        END IF
      END IF
C**GET MEMORY FOR INERTIA TENSOR CALCULATIONS
      KB=NMODE*NMODE
      KAA=NMODE*3*3
      KAB=NMODE*3
      KBB=NMODE
      IF(IREACT.GT.0)THEN
CCCC    KB=NMODE*NMODE*363
CCCC    KAA=NMODE*3*3*362
CCCC    KAB=NMODE*3*362
CCCC    KBB=NMODE*363
        KB=NMODE*NMODE*723
        KAA=NMODE*3*3*722
        KAB=NMODE*3*722
        KBB=NMODE*723
      END IF
      KQQ=NMODE
      CALL MEMO(5,LAB,KAB,LB,KB,LAA,KAA,LBB,KBB,LQQ,KQQ)
      JREACT=IREACT
      IF(IREACT.GT.0)THEN
C**********************************************************************
C**********************************************************************
C**THIS IS THE MAIN CURVELINEAR ALGORITHM
CCCC    CALL REACT(NATOM,NMODE,W(LXM),W(LX0),W(LOMEGA),
CCCC 1  W(LXL),W(LXL+NATOM*NMODE*3*362),W(LXX),W(LXX+NATOM*3*362),
CCCC 2  W(LRR),W(LYZ),W(LRR+NATOM*NATOM*9),W(LRR+NATOM*NATOM*9*363),
CCCC 3  W(LXZ),W(LB+NMODE*NMODE),W(LAA),W(LAB),W(LBB+NMODE),IREACT)
        CALL REACT(NATOM,NMODE,W(LXM),W(LX0),W(LOMEGA),
     1  W(LXL),W(LXL+NATOM*NMODE*3*722),W(LXX),W(LXX+NATOM*3*722),
     2  W(LRR),W(LYZ),W(LRR+NATOM*NATOM*9),W(LRR+NATOM*NATOM*9*723),
     3  W(LXZ),W(LB+NMODE*NMODE),W(LAA),W(LAB),W(LBB+NMODE),IREACT)
        IF(ISCFCI.LT.0)THEN
          WRITE(6,*)'REACTION PATH'
          RETURN
        END IF
        IRET=0
        CALL MEMO(1,LTEMP,NATOM*3,0,0,0,0,0,0,0,0)
C**REARRANGE SO THAT CURVELINEAR COORDINATE IS MODE 'NMODE'
        CALL FLIP(NATOM,NMODE,W(LNBF),W(LMBF),W(LNVF),W(LMVF),
     1  W(LOMEGA),W(LXL),W(LMODNT),W(LTEMP),W(LMXBAS),W(LMEFF),
     2  W(LXMODQ),MAXJ,ICI,IREACT)
        CALL MEMO(-1,LTEMP,NATOM*3,0,0,0,0,0,0,0,0)
C**AT THIS POINT, CURVELINEAR TAU HAS BEEN MOVED TO MODE 3N-6
C**IREACT IS NOW MEANINGLESS.....SET IT TO 3N-6 (NMODE)
        JREACT=NMODE
C**********************************************************************
C**********************************************************************
      END IF
C**IREACT NOW POSITIVE OR ZERO
      IREACT=IABS(JREACT)
C*****************************************************HEG
C*****************************************************HEG
C**SCF CYCLE
      ICYCLE=0
      KH=0
      KXQ=0
      KXW=0
      KXV=0
      KTEMP=0
      DO K=1,NSMODE
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        KH=KH+I1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS, AND RPH (J>0)
        IF(K.GT.NONLIN)KH=KH+I1
        IF(JREACT.GT.0.AND.JMAX.GT.0.AND.K.EQ.NSMODE)KH=KH+I1
        KXQ=KXQ+I2
        KXW=KXW+I2
        KXV=KXV+I2
        CALL INTARR(W(LNBF),W(LNBF),K,I1,I2,I3)
        IF(I1.GT.KTEMP)KTEMP=I1
      END DO
      CALL MEMO(5,LH,KH,LXQ,KXQ,LXW,KXW,LXV,KXV,LQM,NMODE)
      CALL MEMO(2,LXP,KH,LXJ,KXW,0,0,0,0,0,0)
      DO K=1,NMODE
        W(LQM-1+K)=0
      END DO
C**TEMPORARY-LIN
      CALL MEMO(2,LTEMP,KTEMP*2,LSUP4,KTEMP,0,0,0,0,0,0)
C**TEMPORARY-LIN
C**DUMP VIBRATIONAL INFO. TO INP4
      IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0.AND.LDUMP.EQ.0)THEN
C**OUTPUT/INPUT VIBRATION FILE(S) [.VIB]
        READ(INP,*)
        READ(INP,'(A)')FILEV
        WRITE(IOUT,576)FILEV
        OPEN(UNIT=INP4,FILE=FILEV,FORM='UNFORMATTED',
     1  STATUS='UNKNOWN')
        IF(JMAX.GT.0)THEN
C**OUTPUT/INPUT ROTATION FILE(S) [.ROT]
          READ(INP,*)
          READ(INP,'(A)')FILER
          WRITE(IOUT,577)FILER
          OPEN(UNIT=INP5,FILE=FILER,FORM='UNFORMATTED',
     1    STATUS='UNKNOWN')
        END IF
C
        REWIND INP4
C**NEEDS INP4,INP6 IF F-C
C       REWIND INP6
C**RE-POSITION VIBRATIONAL DISC IF REQUIRED
        IF(JMAX.GT.0)CALL DUMMVS(W,W(LMXBAS),NMODE,JMAX)
C
        REWIND INP5
C**NEEDS INP5,INP7 IF F-C
C       REWIND INP7
C**RE-POSITION ROTATIONAL DISC IF REQUIRED
        IF(JMAX.GT.1)CALL DUMMRS(W,JMAX,NMODE)
C
        WRITE(INP4)NMODE
      END IF
      CALL MEMO(2,LWK,NSMODE*4,LWKC1,NSMODE,0,0,0,0,0,0)
C*********************************************************************
C*********************************************************************
      K1=0
      K2=0
      K3=0
      DO K=1,NSMODE
        KX1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        CALL INTARR(W(LNVF),W(LMBF),K,IX1,IX2,IX3)
        KXK=I2
        KWRK=I2*(4*I2+7)
        CALL MEMO(4,LXK,KXK*KXK,LEVAL,KXK,LYK,KXK*KXK,LWRK,KWRK,0,0)
C**NO TORSION FOR TRIATOMIC (NSMODE=NAMODE=3)
        IF(K.EQ.IREACT.AND.NSMODE.NE.NAMODE)THEN
C**GET TORSION INTEGRATION POINTS AND WEIGHTS AND TAU FUNCTIONS
C**GET I3 (ODD) FUNCTIONS C0,S1,C1,S2,C2,........
          CALL TORSPT(I2,W(LH+K2),W(LXQ+K3),W(LXW+K3),W(LXP+K1),
     1    I3/2,0,W(LNVF),NATOM,W(LXX),W(LXM),W(LRR))
CCCC      MMTAU=I2
CCCC      NNTAU=I3
        ELSE
          IF(K.LE.NONLIN)THEN
C**GET HERMITE INTEGRATION POINTS AND WEIGHTS AND H.O. FUNCTIONS
            CALL HERMPT(I3,IX3,I2,W(LH+K2),W(LXQ+K3),W(LXW+K3),
     1      W(LXV+K3),W(LXP+K1),W(LOMEGA+KX1),W(LXK),K,NMODE,NATOM,
     2      W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),
     3      W(LJPOT),W(LCPOT),W(LTEMP),MSYM,W(LXTANH),W(LWK),W(LWKC1),
     4      W(LNP1),W(LCP1),W(LMP1),W(LVP1),W(LDP1),IENTMX(1))
          ELSE
C**GET LAGUERRE INTEGRATION POINTS AND WEIGHTS AND MORSE-LIKE FUNCTIONS
            CALL MORSPT(I3,I2,W(LH+K2),W(LXQ+K3),W(LXW+K3),W(LXV+K3),
     1      W(LOMEGA+KX1),W(LXK),K,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     2      W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     3      W(LTEMP),MSYM,W(LXTANH),W(LWK),W(LWKC1),W(LNP1),W(LCP1),
     4      W(LMP1),W(LVP1),W(LDP1),IENTMX(1),0)
          END IF

          CALL FLUSH(IOUT)
C*********************************************************************
C*********************************************************************
C**MOLPRO LINK INDEX = 2 (FIT 1-DIM GAUSS)
C
C**FIT ABINITIO POTENTIAL
          IF(MOLPRO.EQ.2)THEN
            CALL MEMO(1,LWKC2,I2*I2,0,0,0,0,0,0,0,0)
            CALL GETQ1(W(LXK),W(LEVAL),W(LYK),I2,W(LXQ+K3),
     1      W(LTEMP),I2,K,MSYM,W(LXTANH),W(LWKC2))
            CALL MEMO(-1,LWKC2,I2*I2,0,0,0,0,0,0,0,0)
          END IF
C**INTERPOLATE ABINITIO POTENTIAL
          IF(MOLPRO.EQ.-2)THEN
            CALL GETH1(W(LXQ+K3),I2,W(LTEMP),K,MSYM,W(LXTANH))
          END IF
C
C*********************************************************************
C*********************************************************************
        END IF
        K1=K1+I2*IX3
        K2=K2+I1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
        IF(K.GT.NONLIN)K2=K2+I1
        K3=K3+I2
        CALL MEMO(-4,LXK,KXK*KXK,LEVAL,KXK,LYK,KXK*KXK,LWRK,KWRK,0,0)
      END DO
C*********************************************************************
      WRITE(IOUT,503)VMIN*WAVENM
      IF(IWHICH.GT.0.AND..NOT.ABINIT)THEN
        CALL GETPOT(VPT,NATOM,W(LX0),W(LRR))
        WRITE(IOUT,512)VPT*WAVENM
        WRITE(IOUT,513)VPT
      END IF
      CALL FLUSH(IOUT)
      IF((IABS(MOLPRO).EQ.1.OR.IABS(MOLPRO).EQ.2).AND.IWHICH.EQ.0.
     1AND.ISCFCI.LT.0)THEN
        WRITE(6,*) 'MOLPRO ONLY'
        RETURN
      END IF

C*********************************************************************
      K2=0
      K3=0
      KX1=0
      KX2=0
      KX3=0
      DO K=1,NSMODE
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        KXK=I2
        KWRK=I2*(4*I2+7)
        CALL MEMO(5,LXK,KXK*KXK,LEVAL,KXK,LYK,KXK*KXK,LWRK,KWRK,
     1  LSK,KXK*KXK)
C**NO TORSION FOR TRIATOMIC (NSMODE=NAMODE=3)
        IF(K.EQ.IREACT.AND.NSMODE.NE.NAMODE)THEN
C**GET TORSION INTEGRATION POINTS AND WEIGHTS AND TAU FUNCTIONS
C**GET I3 (EVEN) FUNCTIONS (ZERO, C0,   S1,   C1,........) EVEN K
C**GET I3 (EVEN) FUNCTIONS (S/2,  C/2, 3S/2, 3C/2,........) ODD K
          CALL TORSPT(I2,W(LH+KX2),W(LXQ+KX3),W(LXW+KX3),W(LXP+KX1),
     1    I3/2,1,W(LNVF),NATOM,W(LXX),W(LXM),W(LRR))
          MMTAU=I2
          NNTAU=I3
        ELSE
          IF(IPRINT.GT.0)WRITE(IOUT,506)K
C**SKIP CONTRACTION?
          NOCON=INTFUN(W(LTEMP5),K)
          IF(NOCON.EQ.0)STOP 'ERROR IN NVF'
          IF(NOCON.GT.0)THEN
C**
            MZ2=0
            CALL VCHECK(W(LNBF),W(LMBF),W(LMEFF),NSMODE,K,KZ1,IZ1,MZ1,
     1      K3)
C**TEMPORARY
      WRITE(IOUT,*)'KZ1,IZ1,MZ1 = ',KZ1,IZ1,MZ1
      CALL FLUSH(IOUT)
C**TEMPORARY
            IF(MZ1.NE.0)THEN
              CALL VCHECK(W(LNBF),W(LMBF),W(LMEFF+NSMODE),NSMODE,K,KZ2,
     1        IZ2,MZ2,K3)
C**TEMPORARY
      WRITE(IOUT,*)'KZ2,IZ2,MZ2 = ',KZ2,IZ2,MZ2
      CALL FLUSH(IOUT)
C**TEMPORARY
              IF(MZ2.EQ.0)THEN
                IF(MZ1.NE.IREACT)
     1          CALL VEFF1(W(LXQ+K3),I2,NMODE,NATOM,W(LQQ),W(LRR),
     2          W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),
     3          W(LCPOT),K,W(LXW+K3),W(LXW+K3),W(LMODNT),W(LXTANH),
     4          W(LXQ+KZ1),IZ1,MZ1,W(LNP1),W(LCP1),W(LMP1),W(LVP1),
     5          W(LDP1),IENTMX(1))
                IF(MZ1.EQ.IREACT)
     1          CALL VEFFV1(W(LXQ+K3),I2,NMODE,NATOM,W(LQQ),W(LRR),
     2          W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),
     3          W(LCPOT),K,W(LXW+K3),W(LXW+K3),W(LMODNT),W(LXTANH),
     4          W(LXQ+KZ1),IZ1)
              ELSE
                CALL VEFF2(W(LXQ+K3),I2,NMODE,NATOM,W(LQQ),W(LRR),
     1          W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),
     2          W(LCPOT),K,W(LXW+K3),W(LXW+K3),W(LMODNT),W(LXTANH),
     3          W(LXQ+KZ1),IZ1,MZ1,W(LXQ+KZ2),IZ2,MZ2,W(LNP1),W(LCP1),
     4          W(LMP1),W(LVP1),W(LDP1),IENTMX(1))
              END IF
            ELSE
              CALL DUMPT1(W(LXQ+K3),I2,NMODE,NATOM,W(LQQ),W(LRR),
     1        W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),
     2        W(LCPOT),K,W(LXW+K3),W(LXW+K3),-1,W(LMODNT),W(LXTANH),
     3        W(LNP1),W(LCP1),W(LMP1),W(LVP1),W(LDP1),IENTMX(1),
     4        W(LXMODQ))
            END IF

C**TEMPORARY-LIN
C**FORM HEG POINTS
            CALL HEG0(I3,I2,K,NVB,MVB,W(LNVF),W(LMVF),NSMODE)
C**CONTRACT CURRENT BASIS
            LINLP=0
            IF(LINEAR.AND.K.GT.NONLIN)LINLP=JMAX
            DO L=0,LINLP
              IF(LINEAR.AND.K.GT.NONLIN)THEN
                WRITE(IOUT,*)
                WRITE(IOUT,*)'CONTRACTION OF LINEAR ANGLE BEND : L = ',
     1          L,' FOR MODE ',K
                WRITE(IOUT,*)'NB = ',I3,' MB = ',I2
              END IF
              LL=L
              IF(LL.GT.0)
C**GET LAGUERRE INTEGRATION POINTS AND WEIGHTS AND MORSE-LIKE FUNCTIONS
     1        CALL MORSPT(I3,I2,W(LH+K2),W(LXQ+K3),W(LXW+K3),W(LXV+K3),
     2        W(LOMEGA+K1),W(LXK),K,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     3        W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     4        W(LTEMP),MSYM,W(LXTANH),W(LWK),W(LWKC1),W(LNP1),W(LCP1),
     5        W(LMP1),W(LVP1),W(LDP1),IENTMX(1),LL)
              CALL CONTRA(I3,I2,W(LH+K2),W(LXQ+K3),W(LXW+K3),
     1        W(LOMEGA+K1),W(LSUP4),W(LXK),W(LSK),W(LYK),K,W(LNVF),
     2        NMODE,W(LEVAL),W(LMODNT))
              IDEGEN=0
              NDEGEN=0
C**USE NVB CONTRACTED PRIMITIVES AT MVB 'HEG' POINTS
              CALL HEG1(LL,I3,I2,W(LH+K2),W(LXQ+K3),W(LXJ+K3),
     1        W(LOMEGA+K1),W(LXK),K,W(LTEMP),W(LH+KX2),NVB,MVB,
     2        W(LXP+KX1),W(LXV+KX3),W(LXQ+KX3),W(LMBF),W(LNBF),NMODE,
     3        W(LYK),W(LSUP4),NATOM,W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),
     4        W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LXTANH),W(LWK),
     5        W(LWKC1),W(LNP1),W(LCP1),W(LMP1),W(LVP1),W(LDP1),
     6        IENTMX(1))
              IF(MZ1.NE.0.AND.L.EQ.LINLP)THEN
                WRITE(IOUT,*)
                WRITE(IOUT,*)'MODE ',K,' EFFECTIVE TO MODE ',MZ1,MZ2
              END IF
            END DO
C**TEMPORARY-LIN 
C**NBF, MBF MODIFIED 
C** 
          ELSE 
            CALL HMOVE(W(LH+K2),W(LH+KX2),I3,I3,I2) 
            CALL QMOVE(W(LXQ+K3),W(LXQ+KX3),I2) 
          END IF 
C** 
        END IF
        K2=K2+I1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
        IF(K.GT.NONLIN)K2=K2+I1
        K3=K3+I2
        CALL INTARR(W(LNBF),W(LMBF),K,IX1,IX2,IX3)
        KX1=KX1+IX2*IX3
        KX2=KX2+IX1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
        IF(K.GT.NONLIN)KX2=KX2+IX1
        KX3=KX3+IX2
C       CALL MEMO(-1,LXK,KXK*KXK,0,0,0,0,0,0,0,0)
        CALL MEMO(-5,LXK,KXK*KXK,LEVAL,KXK,LYK,KXK*KXK,LWRK,KWRK,
     1  LSK,KXK*KXK)
      END DO
      CALL MEMO(-1,LXJ,KXW,0,0,0,0,0,0,0,0)
C**************************************************

      WRITE(IOUT,503)VMIN*WAVENM
      IF(IREACT.NE.0)THEN
C**GET POINTERS FOR TAU FUNCTIONS
        K2TAU=0
        K3TAU=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          K2TAU=K2TAU+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(K.GT.NONLIN)K2TAU=K2TAU+IK1
          K3TAU=K3TAU+IK2
        END DO
      END IF
      IF(IREACT.GT.0)CALL MEMO(-1,LWKC1,NSMODE,0,0,0,0,0,0,0,0)
C*********************************************************************
C*********************************************************************
C**MOLPRO LINK INDEX = 3 AND INDEX = 4 (CREATE AND FIT 2-DIM HEG)
C
      IF(IABS(MOLPRO).EQ.3.OR.IABS(MOLPRO).EQ.4)THEN
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            CALL MEMO(2,LV2,3*(2+IK2)*(2+IL2),LWRK,IK2+IL2,0,0,0,0,0,0)
            CALL MOLDP2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1      W(LQQ),W(LXX),W(LX0),W(LXL),W(LXM),K,L,W(LV2),MSYM1,MSYM2,
     2      W(LXTANH),W(LWK),W(LWRK),W(LWRK+IK2),W(LNP1),W(LCP1),
     3      W(LMP1),W(LVP1),W(LDP1),NTOT1,IENTMX(1))
            CALL MEMO(-1,LWRK,IK2+IL2,0,0,0,0,0,0,0,0)
            IF(MOLPRO.EQ.4)THEN
              KWRK=MFIT*(4*MFIT+7)
C             KXK=IK2*IL2
              KXK=MFIT
              CALL MEMO(5,LXK,KXK*KXK,LEVAL,4*KXK,LYK,KXK,LWRK,KWRK,
     1        LWKC2,MFIT*MFIT)
              CALL GETQ2(W(LXK),W(LEVAL),W(LYK),IK2,W(LXQ+K3),IL2,
     1        W(LXQ+L3),W(LV2),K,L,MSYM1,MSYM2,NMODE,W(LQQ),W(LXTANH),
     2        W(LWRK),MFIT,W(LWKC1),W(LWKC2),W(LNP1),W(LCP1),W(LMP1),
     3        NTOT1,IENTMX(1))
              CALL MEMO(-5,LXK,KXK*KXK,LEVAL,4*KXK,LYK,KXK,LWRK,KWRK,
     1        LWKC2,MFIT*MFIT)
            END IF
            CALL MEMO(-1,LV2,3*(2+IK2)*(2+IL2),0,0,0,0,0,0,0,0)
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C*********************************************************************
C*********************************************************************
C**MOLPRO LINK INDEX = 5 AND INDEX = 6 (CREATE AND FIT 3-DIM HEG)
C
      IF(IABS(MOLPRO).EQ.5.OR.IABS(MOLPRO).EQ.6)THEN
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N3=0
            DO N=1,L-1
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
              CALL MEMO(2,LV2,4*(2+IK2)*(2+IL2)*(2+IN2),
     1        LWRK,IK2+IL2+IN2,0,0,0,0,0,0)
              CALL MOLDP3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,IN2,
     1        NMODE,NATOM,W(LQQ),W(LXX),W(LX0),W(LXL),W(LXM),K,L,N,
     2        W(LV2),MSYM1,MSYM2,MSYM3,W(LXTANH),W(LWK),W(LWRK),
     3        W(LWRK+IK2),W(LWRK+IK2+IL2),W(LNP1),W(LCP1),W(LMP1),
     4        W(LVP1),W(LDP1),NTOT1,IENTMX(1),W(LNP2),W(LCP2),W(LMP2),
     5        W(LVP2),W(LDP2A),W(LDP2B),NTOT2,IENTMX(2),
     5        W(LINDK))
              CALL MEMO(-1,LWRK,IK2+IL2+IN2,0,0,0,0,0,0,0,0)
              IF(MOLPRO.EQ.6)THEN
                KWRK=MFIT*(4*MFIT+7)
C               KXK=IK2*IL2*IN2
                KXK=MFIT
                CALL MEMO(5,LXK,KXK*KXK,LEVAL,5*KXK,LYK,KXK,LWRK,KWRK,
     1          LWKC2,MFIT*MFIT)
                CALL GETQ3(W(LXK),W(LEVAL),W(LYK),IK2,W(LXQ+K3),IL2,
     1          W(LXQ+L3),IN2,W(LXQ+N3),W(LV2),K,L,N,MSYM1,MSYM2,MSYM3,
     2          NMODE,W(LQQ),W(LXTANH),W(LWRK),MFIT,W(LWKC1),W(LWKC2),
     3          W(LNP1),W(LCP1),W(LMP1),NTOT1,IENTMX(1),W(LNP2),
     4          W(LCP2),W(LMP2),NTOT2,IENTMX(2),W(LINDK))
                CALL MEMO(-5,LXK,KXK*KXK,LEVAL,5*KXK,LYK,KXK,LWRK,KWRK,
     1          LWKC2,MFIT*MFIT)
              END IF
              CALL MEMO(-1,LV2,4*(2+IK2)*(2+IL2)*(2+IN2),
     1        0,0,0,0,0,0,0,0)
              N3=N3+IN2
            END DO
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C*********************************************************************
C*********************************************************************
C**MOLPRO LINK INDEX = 7 AND INDEX = 8 (CREATE AND FIT 4-DIM HEG)
C
      INDFIT=0
      IF(IFITRP.EQ.4)REWIND 24
      IF(IABS(MOLPRO).EQ.7.OR.IABS(MOLPRO).EQ.8.OR.IABS(IFITRP).EQ.4)
     1THEN
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N3=0
            DO N=1,L-1
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
              M3=0
              DO M=1,N-1
                CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
                IF(IREACT.LE.0)THEN
                  CALL MEMO(2,LV2,5*(2+IK2)*(2+IL2)*(2+IN2)*(2+IM2),
     1            LWRK,IK2+IL2+IN2+IM2,0,0,0,0,0,0)
                  CALL MOLDP4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),W(LXQ+M3),
     1            IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),W(LXX),W(LX0),
     2            W(LXL),W(LXM),K,L,N,M,W(LV2),MSYM1,MSYM2,MSYM3,MSYM4,
     3            W(LXTANH),W(LWK),W(LWRK),W(LWRK+IK2),W(LWRK+IK2+IL2),
     4            W(LWRK+IK2+IL2+IN2),W(LNP1),W(LCP1),W(LMP1),
     5            W(LVP1),W(LDP1),NTOT1,IENTMX(1),W(LNP2),W(LCP2),
     6            W(LMP2),W(LVP2),W(LDP2A),W(LDP2B),NTOT2,IENTMX(2),
     7            W(LNP3),W(LCP3),W(LMP3),W(LVP3),W(LDP3A),W(LDP3B),
     8            W(LDP3C),NTOT3,IENTMX(3),
     8            W(LINDK),W(LINDL))
                  CALL MEMO(-1,LWRK,IK2+IL2+IN2+IM2,0,0,0,0,0,0,0,0)
                  IF(MOLPRO.EQ.8)THEN
                    KWRK=MFIT*(4*MFIT+7)
C                   KXK=IK2*IL2*IN2*IM2
                    KXK=MFIT
                    CALL MEMO(5,LXK,KXK*KXK,LEVAL,6*KXK,LYK,KXK,
     1              LWRK,KWRK,LWKC2,MFIT*MFIT)
                    CALL GETQ4(W(LXK),W(LEVAL),W(LYK),IK2,W(LXQ+K3),
     1              IL2,W(LXQ+L3),IN2,W(LXQ+N3),IM2,W(LXQ+M3),W(LV2),K,
     2              L,N,M,MSYM1,MSYM2,MSYM3,MSYM4,NMODE,W(LQQ),
     3              W(LXTANH),W(LWRK),MFIT,W(LWKC1),W(LWKC2),W(LNP1),
     4              W(LCP1),W(LMP1),NTOT1,IENTMX(1),W(LNP2),W(LCP2),
     5              W(LMP2),NTOT2,IENTMX(2),W(LNP3),W(LCP3),W(LMP3),
     6              NTOT3,IENTMX(3),W(LINDK),W(LINDL))
                    CALL MEMO(-5,LXK,KXK*KXK,LEVAL,6*KXK,LYK,KXK,
     1              LWRK,KWRK,LWKC2,MFIT*MFIT)
                  END IF
                  CALL MEMO(-1,LV2,5*(2+IK2)*(2+IL2)*(2+IN2)*(2+IM2),
     1            0,0,0,0,0,0,0,0)
                ELSE
C**INITIAL CALL (IND=0) SETS UP MFIT,NFIT AND MH1,MH2,MH3,MH4
                  IND=0
                  IF(IFITRP.LT.0)THEN
                    MDIM=IK2*IL2*IN2*IM2
                    CALL MEMO(1,LYK,MDIM,0,0,0,0,0,0,0,0)
                  END IF
                  CALL MOLVP4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),W(LXQ+M3),
     1            IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),W(LXX),W(LXL),
     2            W(LXM),W(LWRK),K,L,N,M,W(LV2),W(LXK),MOLINC,
     3            W(LXTANH),W(LEVAL),W(LYK),W(LWKC1),MDIM,W(LWKC2),
     4            W(LMODNT),MH1,MH2,MH3,MH4,IND,W(LNP4),NTOT4,MAXPOW,
     5            INDFIT)
                  IF(IFITRP.LT.0)CALL MEMO(-1,LYK,MDIM,0,0,0,0,0,0,0,0)
                  IF(IFITRP.GT.0)THEN
                    KWRK=MFIT*(4*MFIT+7)
                    KXK=MFIT
                    CALL MEMO(2,LV2,MFIT,LWRK,NATOM*NATOM,0,0,0,0,0,0)
                    CALL MEMO(5,LXK,KXK*KXK,LEVAL,3*KXK,LYK,KXK,
     1              LWKC1,KWRK,LWKC2,MFIT*MFIT)
                    IND=1
                    CALL MOLVP4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1              W(LXQ+M3),IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),
     2              W(LXX),W(LXL),W(LXM),W(LWRK),K,L,N,M,W(LV2),W(LXK),
     3              MOLINC,W(LXTANH),W(LEVAL),W(LYK),W(LWKC1),MFIT,
     4              W(LWKC2),W(LMODNT),MH1,MH2,MH3,MH4,IND,W(LNP4),
     5              NTOT4,MAXPOW,INDFIT)
                    CALL MEMO(-5,LXK,KXK*KXK,LEVAL,3*KXK,LYK,KXK,
     1              LWKC1,KWRK,LWKC2,MFIT*MFIT)
                    CALL MEMO(-2,LV2,MFIT,LWRK,NATOM*NATOM,
     1              0,0,0,0,0,0)
                  END IF
                END IF
                M3=M3+IM2
              END DO
              N3=N3+IN2
            END DO
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
      IF(IFITRP.EQ.-4)THEN
        MAXPOW=MAXPOW+1
        REWIND 24
      END IF
C*********************************************************************
C*********************************************************************
C**MOLPRO LINK INDEX = 9 AND INDEX = 10 (CREATE AND FIT 5-DIM HEG)
C
      IF(IABS(MOLPRO).EQ.9.OR.IABS(MOLPRO).EQ.10)THEN
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N3=0
            DO N=1,L-1
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
              M3=0
              DO M=1,N-1
                CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
                  J3=0
                  DO J=1,M-1
                  CALL INTARR(W(LNBF),W(LMBF),J,IJ1,IJ2,IJ3)
                  CALL MEMO(2,LV2,2*IK2*IL2*IN2*IM2*IJ2,
     1            LWRK,IK2+IL2+IN2+IM2+IJ2,0,0,0,0,0,0)
                  CALL MOLDP5(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),W(LXQ+M3),
     1            W(LXQ+J3),IK2,IL2,IN2,IM2,IJ2,NMODE,NATOM,W(LQQ),
     2            W(LXX),W(LX0),W(LXL),W(LXM),K,L,N,M,J,W(LV2),MSYM1,
     3            MSYM2,MSYM3,MSYM4,MSYM5,W(LXTANH),W(LWK),W(LWRK),
     4            W(LWRK+IK2),W(LWRK+IK2+IL2),W(LWRK+IK2+IL2+IN2),
     5            W(LWRK+IK2+IL2+IN2+IM2),W(LNP1),W(LCP1),W(LMP1),
     6            NTOT1,IENTMX(1),W(LNP2),W(LCP2),W(LMP2),NTOT2,
     7            IENTMX(2),W(LNP3),W(LCP3),W(LMP3),NTOT3,IENTMX(3),
     8            W(LNP4),W(LCP4),W(LMP4),NTOT4,IENTMX(4),W(LINDK),
     9            W(LINDL),W(LINDN))
                  CALL MEMO(-1,LWRK,IK2+IL2+IN2+IM2+IJ2,
     1            0,0,0,0,0,0,0,0)
                  IF(MOLPRO.EQ.10)THEN
                    KWRK=MFIT*(4*MFIT+7)
C                   KXK=IK2*IL2*IN2*IM2*IJ2
                    KXK=MFIT
                    CALL MEMO(5,LXK,KXK*KXK,LEVAL,7*KXK,LYK,KXK,
     1              LWRK,KWRK,LWKC2,MFIT*MFIT)
                    CALL GETQ5(W(LXK),W(LEVAL),W(LYK),IK2,W(LXQ+K3),
     1              IL2,W(LXQ+L3),IN2,W(LXQ+N3),IM2,W(LXQ+M3),IJ2,
     2              W(LXQ+J3),W(LV2),K,L,N,M,J,MSYM1,MSYM2,MSYM3,MSYM4,
     3              MSYM5,NMODE,W(LQQ),W(LXTANH),W(LWRK),MFIT,W(LWKC1),
     4              W(LWKC2),W(LNP1),W(LCP1),W(LMP1),NTOT1,IENTMX(1),
     5              W(LNP2),W(LCP2),W(LMP2),NTOT2,IENTMX(2),W(LNP3),
     6              W(LCP3),W(LMP3),NTOT3,IENTMX(3),W(LNP4),W(LCP4),
     7              W(LMP4),NTOT4,IENTMX(4),W(LINDK),W(LINDL),W(LINDN))
                    CALL MEMO(-5,LXK,KXK*KXK,LEVAL,7*KXK,LYK,KXK,
     1              LWRK,KWRK,LWKC2,MFIT*MFIT)
                  END IF
                  CALL MEMO(-1,LV2,2*IK2*IL2*IN2*IM2*IJ2,
     1            0,0,0,0,0,0,0,0)
                  J3=J3+IJ2
                END DO
                M3=M3+IM2
              END DO
              N3=N3+IN2
            END DO
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C
      IF(IREACT.LE.0)CALL MEMO(-1,LWKC1,NSMODE,0,0,0,0,0,0,0,0)
      CALL MEMO(-1,LWK,NSMODE*4,0,0,0,0,0,0,0,0)
C***********************
      CALL MEMO(-2,LTEMP,KTEMP*2,LSUP4,KTEMP,0,0,0,0,0,0)
      IF(MOLPRO.NE.0.AND.ISCFCI.LT.0)THEN
        WRITE(6,*) 'MOLPRO ONLY'
        RETURN
      END IF
      IF(IWHICH.LT.0.AND.ISCFCI.LT.0)THEN
        WRITE(6,*) 'HEG ONLY'
        RETURN
      END IF
      IF(IFITRP.GT.0)THEN
        WRITE(6,*) 'RPH FIT ONLY'
        RETURN
      END IF
C*****************************************************HEG
C*****************************************************HEG
      MBFMIN=100000
      DO K=1,NSMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(IK2.LT.MBFMIN)MBFMIN=IK2
      END DO
C**FIND MINIMUM OF POTENTIAL (INCLUDE LINEAR MOLECULES)
      IF(MCHECK.LT.0.AND.JREACT.LE.0)THEN
        WRITE(IOUT,405)
C**START WITH CUBIC FIT
        N=4
1111    CONTINUE
        IF(N.LE.MBFMIN)WRITE(IOUT,410)N
        INDEX=1
C**4*N IS BIG ENOUGH FOR E04FBF FOR N=3
        CALL MEMO(4,LXK,N*N,LEVAL,N,LSUP4,N,LWRK,4*N,0,0)
        K3=0
C**ONLY LOOK AT 'ACTIVE' NORMAL MODES
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
          MODD=K
          CALL FITQ1(W(LXK),W(LEVAL),W(LSUP4),W(LWRK),N,W(LXQ+K3),
     1    W(LXV+K3),I2,W(LQM),NMODE,K,INDEX,W(LXTANH))
          K3=K3+I2
        END DO
        CALL MEMO(-4,LXK,N*N,LEVAL,N,LSUP4,N,LWRK,4*N,0,0)
        IF(INDEX.EQ.0)THEN
          CALL GETMIN(W(LQM),NMODE,NATOM,W(LXX),W(LX0),W(LXM),W(LXL))
          K3=0
C**ONLY LOOK AT 'ACTIVE' NORMAL MODES
          DO K=1,NAMODE
            CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
            CALL GETV(I2,W(LXQ+K3),W(LXV+K3),K,NMODE,NATOM,W(LQQ),
     1      W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),
     2      W(LJPOT),W(LCPOT),W(LNP1),W(LCP1),W(LMP1),W(LVP1),W(LDP1),
     3      IENTMX(1))
            K3=K3+I2
          END DO
          N=N+2
          GO TO 1111
        END IF
        CALL CHECKM(W(LXM),W(LX0),W(LXX),W(LRR),NATOM)
      END IF
C**RELATE PRINCIPAL AXES TO INPUT AXES
C     CALL PAXES(W(LXM),W(LX0),W(LXX),W(LRR),NATOM,1)
      IF(MCHECK.LT.0.AND.JREACT.LE.0)THEN
        WRITE(IOUT,350)
        CALL FLUSH(IOUT)
        CALL GETNOR(NATOM,NMODE,W(LXM),W(LX0),W(LOMEGA),
     1  W(LXL),W(LXX),W(LRR),W(LIPOT),W(LJPOT),W(LCPOT),NPOT)
        IF(LINEAR)THEN
C**MOVE LINEAR VECTORS TO TOP
          CALL MEMO(1,LTEMP,NATOM*3,0,0,0,0,0,0,0,0)
          CALL BLIP(NATOM,NMODE,W(LNBF),W(LMBF),W(LNVF),W(LOMEGA),
     1    W(LXL),W(LMODNT),W(LMEFF),W(LTEMP),W(LMXBAS),NSMODE,MAXJ,ICI,
     2    0)
          CALL MEMO(-1,LTEMP,NATOM*3,0,0,0,0,0,0,0,0)
        END IF
        IF(ISCFCI.LT.0)THEN
          WRITE(6,*) 'MINIMUM POTENTIAL RUN ONLY'
          RETURN
        END IF
      END IF
C**RELATE PRINCIPAL AXES TO INPUT AXES
      CALL PAXES(W(LXM),W(LX0),W(LXX),W(LRR),W(LXL),NMODE,NATOM,1)
      CALL INPUT(NMODE,W(LOMEGA),W(LXX)
     1,W(LNBF),W(LMBF),NSTAT,W(LISTAT),NATOM,W(LXM),W(LX0),W(LXL),
     2W(LRR),ICI,W(LMXBAS),MAXJ,W(LJSTAT),W(LKSTAT),W(LESTAT),
     3W(LWSTAT),W(LNVF),W(LMVF),W(LTEMP5),INORM,MCHECK,NSMODE,NVMODE,1,
     4W)
C**GENERATE ECKART GRID
C     IF(MDUMP.NE.0)THEN
C       IRET=1
C       CALL ECKART(NATOM,NMODE,W(LX0),W(LXL),W(LXM),W(LXX),W(LYZ),
C    1  W(LRR),W(LMBF),W(LXQ),W(LQQ),NAMODE)
C       IF(KDUMP.LT.0)THEN
C         WRITE(6,*) 'ECKART'
C         RETURN
C       END IF
C       IRET=0
C     END IF
C**CALCULATE ZETA CONSTANTS AND EQUILIBRIUM BOND LENGTHS
      KR0=NATOM*NATOM
      CALL MEMO(1,LR0,KR0,0,0,0,0,0,0,0,0)
      CALL RZETA(NATOM,NMODE,W(LXL),W(LXZ),W(LX0),W(LR0),W(LXM),NMODE)
      CALL FLUSH(IOUT)
      KESCF=NSTAT
      CALL MEMO(1,LESCF,KESCF,0,0,0,0,0,0,0,0)
      ITIM=-1
      J21=2*JMAX+1
      CALL MEMO(2,LSSX,J21*9*J21,LSS,J21*9*J21,0,0,0,0,0,0)
      CALL MEMO(3,LW21,J21,LX21,J21*J21,LE21,J21,0,0,0,0)
C**GET ROTATIONAL MATRIX ELEMENTS FOR K-DIAGONAL
      CALL ROTEL(J21-1,W(LSS),W(LSSX),J21,1,1,1)
      IF(MATSIZ.NE.0.AND.ICI.LT.0.AND.LDUMP.EQ.0)GO TO 4000
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                          RESTART STAGE (PROPERTIES)
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**READ CI COEFFS AND CALCULATE PROPERTY
      IF(LDUMP.NE.0)THEN
        NNMODS=NNMODE
        NAMODS=NAMODE
        NVMODS=NVMODE
        LLCNT=1
        IF(KCONT.GT.0)LLCNT=NCONT
        IF(LLCNT.EQ.1)ICONT(1)=NSMODE
        WRITE(IOUT,*)
      LLDUMP=IABS(LDUMP)
      IF(LLDUMP.EQ.1)WRITE(IOUT,*)'RESTART FOR DIPOLE MATRIX ELEMENTS'
      IF(LLDUMP.EQ.2)WRITE(IOUT,*)'RESTART FOR EXPECTATION VALUE OF Q'
      IF(LLDUMP.EQ.3)WRITE(IOUT,*)'RESTART FOR EXPECTATION VALUE OF Q^2'
      IF(LLDUMP.EQ.4)WRITE(IOUT,*)'RESTART FOR WAVEFUNCTION PLOTS'
      IF(LLDUMP.EQ.5)WRITE(IOUT,*)'RESTART FOR USER-DEFINED PROPERTY'
      IF(LLDUMP.EQ.6)WRITE(IOUT,*)'RESTART FOR FRANCK-CONDON FACTORS'
      IF(LLDUMP.EQ.7)WRITE(IOUT,*)'RESTART FOR P,Q,R-BRANCH INTENSITIES'
      IF(LLDUMP.EQ.8)WRITE(IOUT,*)'RESTART FOR RAMAN INTENSITIES'
      IF(LLDUMP.EQ.9)WRITE(IOUT,*)'RESTART FOR PARTITION FUNCTIONS'
        WRITE(IOUT,*)
        CALL RESUME(W,LDUMP)
      END IF

C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                                 END OF SET-UP STAGE    
C**********************************************************************
C**********************************************************************
C**********************************************************************
      WRITE(IOUT,325)
      WRITE(IOUT,325)
      WRITE(IOUT,220)
      CALL FLUSH(IOUT)
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                      DUMP POTENTIAL AND CORIOLIS DATA TO DISC (SCF)
C**********************************************************************
C**********************************************************************
C**********************************************************************
      ICOUPL=ICOUPS
      ICOUPC=ICOUCS
      CALL DISCDP(W,KCONT,0)
      ICOUPL=ICOUPX
      ICOUPC=ICOUCX
      CALL MEMO(-1,LXMODQ,NSMODE,0,0,0,0,0,0,0,0)
      IF(JREACT.GT.0)THEN
C**FIRST ORTHONORMALISE TAU FUNCTIONS (COS THEN SIN)
        LMAX=1
        IF(JMAX.GT.0)LMAX=2
        DO IND=1,LMAX
          IF(IND.EQ.1)WRITE(IOUT,*)'INTEGRAL TORSION FOR EVEN Ka'
          IF(IND.EQ.2)WRITE(IOUT,*)'HALF-INTEGRAL TORSION FOR ODD Ka'
          DO ICS=1,2
            IF(ICS.EQ.1)WRITE(IOUT,*)'COS FUNCTIONS'
            IF(ICS.EQ.2)WRITE(IOUT,*)'SIN FUNCTIONS'
            REWIND 71
            REWIND 81
            CALL INTARR(W(LNBF),W(LMBF),NSMODE,IK1TAU,IK2TAU,IK3TAU)
            KXK=IK3TAU/2-(2-IND)*(ICS-1)
            KWK=IK1TAU
            CALL MEMO(3,LYK,IK3TAU*IK3TAU,LXK,KXK*KXK,LOV,KXK*KXK,
     1      0,0,0,0)
C**ZEROISE HAMILTONIAN AND OVERLAP MATRICES
            CALL MTZERO(W(LXK),KXK)
            CALL MTZERO(W(LOV),KXK)
            CALL TAUM0(W(LH+K2TAU),W(LXQ+K3TAU),IK3TAU,IK2TAU,W(LXK),
     1      W(LOV),KXK,NMODE,W(LC0),W(LC0),W(LMODNT),MDTAU,IND)
            KXA=KXK*KXK
            CALL MEMO(4,LWK,KWK,LSUP4,KXK,LEVAL,KXK,LXA,KXA,0,0)
            WRITE(IOUT,*)'Calculating DIAG'
            CALL TIMIT(1)
            CALL FLUSH(IOUT)
            IRET=1
            IROTV=1
            CALL DIAG(W(LXK),W(LXK),KXK,KXK,1,W(LSUP4),W(LEVAL),W(LWK),
     1      KXK,KXK,W(LISTAT),NSTAT,NMODE,W(LOV),W(LXA),W(LXK),IDUM,
     2      IDUM,IDUM)
            IRET=0
            IROTV=0
            CALL TIMIT(3)
            CALL FLUSH(IOUT)
            CALL MEMO(-3,LSUP4,KXK,LEVAL,KXK,LXA,KXA,0,0,0,0)
C**LINEAR COMBINATIONS SHOULD BE ORTHONORMAL
            CALL REFORM(W(LH+K2TAU),IK3TAU,W(LXK),KXK,IK3TAU,IK2TAU,
     1      W(LWK),W(LYK),W(LOV),0,IND,W(LEVAL),W(LMODNT))
            CALL MEMO(-4,LXK,KXK*KXK,LWK,KWK,LYK,IK3TAU*IK3TAU,
     1      LOV,KXK*KXK,0,0)
          END DO
        END DO
C**CHANGE NBF TO NVF (LIKE HEG CONTRACTION)
        CALL INTARR(W(LNVF),W(LMBF),NSMODE,KK1TAU,KK2TAU,KK3TAU)
        CALL MEMO(1,LTEMP,KK3TAU*IK2TAU*6,0,0,0,0,0,0,0,0)
        CALL TAUMID(W(LH+K2TAU),IK3TAU,W(LTEMP),KK3TAU,IK2TAU)
C       CALL TAUBAS(W(LH+K2TAU),IK3TAU,W(LNBF),W(LH+K2TAU),KK3TAU,
C    1  W(LNVF),NSMODE,IK2TAU)
        CALL TAUBAS(W(LTEMP),KK3TAU,W(LNBF),W(LH+K2TAU),KK3TAU,
     1  W(LNVF),NSMODE,IK2TAU)
        CALL MEMO(-1,LTEMP,KK3TAU*IK2TAU*6,0,0,0,0,0,0,0,0)
      END IF

C**ANALYTIC
      IF((IWHICH.LT.0.AND.MOLINC.GT.0).OR.IFITRP.NE.0)THEN
        MAXQU=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          IF(IK3.GT.MAXQU)MAXQU=IK3
        END DO
        CALL MEMO(1,LXKAN,MAXQU*MAXQU*NAMODE*MAXPOW,0,0,0,0,0,0,0,0)
        K2=0
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          CALL GETANI(W(LXKAN),MAXQU,MAXPOW,K,W(LH+K2),W(LXQ+K3),IK3,
     1    IK2,W(LXTANH))
          K2=K2+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(K.GT.NONLIN)K2=K2+IK1
          K3=K3+IK2
        END DO
      END IF
C**ANALYTIC

C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                      START OF VIBRATIONAL SCF STAGE
C**********************************************************************
C**********************************************************************
C**********************************************************************
C     ICOUPS=ICOUPL
C     ICOUCS=ICOUPC
      ICOUPL=ICOUPS
      ICOUPC=ICOUCS
      ICPSS=ICOUPS
      ICCSS=ICOUCS
      ICPXSV=ICOUPX
      ICCXSV=ICOUCX
      NCL1SV=NCOUPL(1)
      NCL2SV=NCOUPL(2)
      NCC1SV=NCOUPC(1)
      NCC2SV=NCOUPC(2)
      IF(IREACT.EQ.0.AND.ICOUPX.GT.5)ICOUPX=5
      IF(IREACT.NE.0.AND.ICOUPX.GT.4)ICOUPX=4
      IF(ICOUCX.GT.4)ICOUCX=4
      IF(IREACT.EQ.0.AND.ICOUPS.GT.5)ICOUPL=5
      IF(IREACT.NE.0.AND.ICOUPS.GT.4)ICOUPL=4
      IF(ICOUCS.GT.4)ICOUPC=4
      IF(IREACT.EQ.0.AND.ICOUPS.GT.5)ICOUPS=5
      IF(IREACT.NE.0.AND.ICOUPS.GT.4)ICOUPS=4
      IF(ICOUCS.GT.4)ICOUCS=4
      IF(IREACT.EQ.0.AND.NCOUPL(1).GT.5)NCOUPL(1)=5
      IF(IREACT.NE.0.AND.NCOUPL(1).GT.4)NCOUPL(1)=4
      IF(IREACT.EQ.0.AND.NCOUPL(2).GT.5)NCOUPL(2)=5
      IF(IREACT.NE.0.AND.NCOUPL(2).GT.4)NCOUPL(2)=4
      IF(NCOUPC(1).GT.4)NCOUPC(1)=4
      IF(NCOUPC(2).GT.4)NCOUPC(2)=4
      CALL VRSCF(W,KCONT,NSTAT,CONV,K2TAU,K3TAU,MDTAU)
      IDEGEN=0
      NDEGEN=0
      NCOUPL(1)=NCL1SV
      NCOUPL(2)=NCL2SV
      NCOUPC(1)=NCC1SV
      NCOUPC(2)=NCC2SV
      ICOUPX=ICPXSV
      ICOUCX=ICCXSV
      ICOUPL=ICOUPX
      ICOUPC=ICOUCX
      ICOUPS=ICPSS
      ICOUCS=ICCSS
C**********************************************************************
C**********************************************************************
C**                                                  RETURN IF SCF ONLY
      IF(ISCFCI.EQ.0)RETURN
      CALL MEMO(-1,LTEMP5,NSMODE)
C**********************************************************************
C**MAKE COPY (59) -> (INP4) AND WRITE VIRTUAL FUNCTIONS + 1ST, 2ND DERIVS.
      IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0)THEN
        REWIND 59
        K2=0
        DO MODE=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),MODE,IMODE1,IMODE2,IMODE3)
          CALL MEMO(1,LWK,IMODE3,0,0,0,0,0,0,0,0)
CCCC      CALL RECYCL(W(LWK),IMODE3,59,60)
          INNCR=0
          LHDIM=IMODE3
          IF(IREACT.GT.0.AND.MODE.EQ.NSMODE)THEN
            INNCR=1
C           LHDIM=IK3TAU
          END IF
C**RPH ODD K
          CALL OUTVCI(IMODE3,IMODE2,W(LH+K2),LHDIM,W(LWK),INNCR,LMAX)
C**RPH ODD K
          CALL MEMO(-1,LWK,IMODE3,0,0,0,0,0,0,0,0)
          K2=K2+IMODE1
        END DO
      END IF
      CALL FLUSH(IOUT)
      IF(ISCFCI.EQ.0)THEN
        WRITE(IOUT,215)
        RETURN
      ELSE
        WRITE(IOUT,325)
        WRITE(IOUT,325)
        WRITE(IOUT,235)
      END IF
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                        END OF VIBRATIONAL SCF STAGE
C**********************************************************************
C**********************************************************************
C**********************************************************************
4000  CONTINUE
C**RESET NMAX IF UNRESTRICTED
      IF(NNMAX.EQ.0)NMAX=-ICI*NVMODE
      NBFMIN=100000
      DO K=1,NVMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(IK3.LT.NBFMIN)NBFMIN=IK3
      END DO
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                                  START OF VCI STAGE
C**********************************************************************
C**********************************************************************
C**********************************************************************
      DO I=1,NVSYM
        NTOTL1(I)=0
        NTOTL2(I)=0
        NTOTL3(I)=0
      END DO
      IF(ICI.LT.0.AND.NMAX.GE.0)THEN
        JCI=-ICI
        IF(NMAX.GT.0)THEN
C**ONLY GET NMAX +1 (MAXIMUM) INTEGRALS
          IF(NMAX+1.LT.JCI)JCI=NMAX+1
C**OR NBFMIN (MAXIMUM) INTEGRALS (WHICHEVER IS SMALLER)
          IF(JCI.GT.NBFMIN)JCI=NBFMIN
C**CONVERT NMAX FROM QUANTA TO NUMBER OF FUNCTIONS
          NMAX=NMAX+NVMODE
        END IF
      END IF
      IF(ICI.LT.0)THEN
        NVSYMX=NVSYM
        NVSYM=1
C**NNMAX USED AS INDICATOR AS TO WHICH ALGORITHM
        NNMAX=NMAX
        IF(NNMAX.GE.0)THEN
          DO I=1,6
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
          JCI6=JCI
        ELSE
C**MAX POSSIBLE QUANTA
          IF(IREACT.NE.0)CALL TMAXIN(W(LMXBAS),NMODE,MAXJ)
          JCI1=0
          JCI2=0
          JCI3=0
          JCI4=0
          JCI5=0
          JCI6=0
          DO K=1,NCONT
            ICI1=NBAS(1,K)+1
            ICI2=MAX0(NBAS(1,K)+1,NBAS(2,K)+1)
            ICI3=MAX0(NBAS(1,K)+1,NBAS(2,K)+1,NBAS(3,K)+1)
            ICI4=MAX0(NBAS(1,K)+1,NBAS(2,K)+1,NBAS(3,K)+1,NBAS(4,K)+1)
            ICI5=MAX0(NBAS(1,K)+1,NBAS(2,K)+1,NBAS(3,K)+1,NBAS(4,K)+1,
     1      NBAS(5,K)+1)
            ICI6=MAX0(NBAS(1,K)+1,NBAS(2,K)+1,NBAS(3,K)+1,NBAS(4,K)+1,
     1      NBAS(5,K)+1,NBAS(6,K)+1)
            IF(ICI1.GT.JCI1)JCI1=ICI1
            IF(ICI2.GT.JCI2)JCI2=ICI2
            IF(ICI3.GT.JCI3)JCI3=ICI3
            IF(ICI4.GT.JCI4)JCI4=ICI4
            IF(ICI5.GT.JCI5)JCI5=ICI5
            IF(ICI6.GT.JCI6)JCI6=ICI6
          END DO
        END IF
        JCC1=JCI1
        JCC2=JCI2
        JCC3=JCI3
        JCC4=JCI4
        JCC5=JCI5
        JCC6=JCI6
C*******************************************************
C**GET SIZES OF 1-,2-,3-,4-,5-,6- MODE BASIC BASES - START
C*******************************************************
        IF(ICOUPL.GT.0.OR.(NNMAX.LT.0.AND.ICI.LT.0))THEN
C         NMAX1=MAXSUM(1,1)+1
C         NMAX2=MAXSUM(1,2)+1
          NMAX1=0
          NMAX2=0
          DO I=1,MAXJ
            NMAX1=MAX0(NMAX1,MAXSUM(I,1)+1)
            NMAX2=MAX0(NMAX2,MAXSUM(I,2)+1)
          END DO
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE1=JCI1
C**FIRST GET 1-D SIZE FOR 'NORMAL' (ISIZE1 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE1,1,JCI1,0)
          NSIZE1=ISIZE1
          KIP1=ISIZE1
        END IF
C
        IF((ICOUPL.GT.1.OR.(NNMAX.LT.0.AND.ICI.LT.-1)).OR.
     1    (ICOUPL.EQ.1.AND.IREACT.NE.0))THEN
C         NMAX1=MAX0(MAXSUM(1,1)+2,MAXSUM(2,1)+2)
C         NMAX2=MAX0(MAXSUM(1,2)+2,MAXSUM(2,2)+2)
          NMAX1=0
          NMAX2=0
          DO I=1,MAXJ
            NMAX1=MAX0(NMAX1,MAXSUM(I,1)+2)
            NMAX2=MAX0(NMAX2,MAXSUM(I,2)+2)
          END DO
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE2=JCI2**2
C**FIRST GET 2-D SIZE FOR 'NORMAL' (ISIZE2 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE2,2,JCI2,0)
          NSIZE2=ISIZE2
          IF(LINEAR.AND..NOT.TRIAT)THEN
            NMAX1=MAXSUM(1,1)+1
            NMAX2=MAXSUM(1,2)+1
            NMAX=MAX0(NMAX1,NMAX2)
            LSIZE2=ISIZE1
C**FIRST GET 1-D SIZE FOR 'LINEAR' (LSIZE2 MODIFIED IN GETSZ FOR NMAX)
            CALL GETSZ(LSIZE2,1,JCI1,1)
            KLJP2=LSIZE2*2
          END IF
          KIP2=ISIZE2*2
        END IF
C
        IF((ICOUPL.GT.2.OR.(NNMAX.LT.0.AND.ICI.LT.-2)).OR.
     1    (ICOUPL.EQ.2.AND.IREACT.NE.0))THEN
C         NMAX1=MAX0(MAXSUM(1,1)+3,MAXSUM(2,1)+3,MAXSUM(3,1)+3)
C         NMAX2=MAX0(MAXSUM(1,2)+3,MAXSUM(2,2)+3,MAXSUM(3,2)+3)
          NMAX1=0
          NMAX2=0
          DO I=1,MAXJ
            NMAX1=MAX0(NMAX1,MAXSUM(I,1)+3)
            NMAX2=MAX0(NMAX2,MAXSUM(I,2)+3)
          END DO
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE3=JCI3**3
C**FIRST GET 3-D SIZE FOR 'NORMAL' (ISIZE3 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE3,3,JCI3,0)
          NSIZE3=ISIZE3
          IF(LINEAR.AND..NOT.TRIAT)THEN
            NMAX1=MAX0(MAXSUM(1,1)+2,MAXSUM(2,1)+2)
            NMAX2=MAX0(MAXSUM(1,2)+2,MAXSUM(2,2)+2)
            NMAX=MAX0(NMAX1,NMAX2)
            LSIZE3=ISIZE2
C**FIRST GET 2-D SIZE FOR 'LINEAR' (LSIZE3 MODIFIED IN GETSZ FOR NMAX)
            CALL GETSZ(LSIZE3,2,JCI2,1)
            KLJP3=LSIZE3*3
          END IF
          KIP3=ISIZE3*3
        END IF
C
        IF((ICOUPL.GT.3.OR.(NNMAX.LT.0.AND.ICI.LT.-3)).OR.
     1    (ICOUPL.EQ.3.AND.IREACT.NE.0))THEN
C         NMAX1=MAX0(MAXSUM(1,1)+4,MAXSUM(2,1)+4,MAXSUM(3,1)+4,
C    1    MAXSUM(4,1)+4)
C         NMAX2=MAX0(MAXSUM(1,2)+4,MAXSUM(2,2)+4,MAXSUM(3,2)+4,
C    1    MAXSUM(4,2)+4)
          NMAX1=0
          NMAX2=0
          DO I=1,MAXJ
            NMAX1=MAX0(NMAX1,MAXSUM(I,1)+4)
            NMAX2=MAX0(NMAX2,MAXSUM(I,2)+4)
          END DO
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE4=JCI4**4
C**FIRST GET 4-D SIZE FOR 'NORMAL' (ISIZE4 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE4,4,JCI4,0)
          NSIZE4=ISIZE4
          IF(LINEAR.AND..NOT.TRIAT)THEN
            NMAX1=MAX0(MAXSUM(1,1)+3,MAXSUM(2,1)+3,MAXSUM(3,1)+3)
            NMAX2=MAX0(MAXSUM(1,2)+3,MAXSUM(2,2)+3,MAXSUM(3,2)+3)
            NMAX=MAX0(NMAX1,NMAX2)
            LSIZE4=ISIZE3 
C**FIRST GET 3-D SIZE FOR 'LINEAR' (LSIZE4 MODIFIED IN GETSZ FOR NMAX)
            CALL GETSZ(LSIZE4,3,JCI3,1)
            KLJP4=LSIZE4*4
          END IF
          KIP4=ISIZE4*4
        END IF
C
        IF((ICOUPL.GT.4.OR.(NNMAX.LT.0.AND.ICI.LT.-4)).OR.
     1    (ICOUPL.EQ.4.AND.IREACT.NE.0))THEN
C         NMAX1=MAX0(MAXSUM(1,1)+5,MAXSUM(2,1)+5,MAXSUM(3,1)+5,
C    1    MAXSUM(4,1)+5,MAXSUM(5,1)+5)
C         NMAX2=MAX0(MAXSUM(1,2)+5,MAXSUM(2,2)+5,MAXSUM(3,2)+5,
C    1    MAXSUM(4,2)+5,MAXSUM(5,2)+5)
          NMAX1=0
          NMAX2=0
          DO I=1,MAXJ
            NMAX1=MAX0(NMAX1,MAXSUM(I,1)+5)
            NMAX2=MAX0(NMAX2,MAXSUM(I,2)+5)
          END DO
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE5=JCI5**5
C**FIRST GET 5-D SIZE FOR 'NORMAL' (ISIZE5 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE5,5,JCI5,0)
          NSIZE5=ISIZE5
          IF(LINEAR.AND..NOT.TRIAT)THEN
            NMAX1=MAX0(MAXSUM(1,1)+4,MAXSUM(2,1)+4,MAXSUM(3,1)+4,
     1      MAXSUM(4,1)+4)
            NMAX2=MAX0(MAXSUM(1,2)+4,MAXSUM(2,2)+4,MAXSUM(3,2)+4,
     1      MAXSUM(4,2)+4)
            NMAX=MAX0(NMAX1,NMAX2)
            LSIZE5=ISIZE4
C**FIRST GET 4-D SIZE FOR 'LINEAR' (LSIZE5 MODIFIED IN GETSZ FOR NMAX)
            CALL GETSZ(LSIZE5,4,JCI4,1)
            KLJP5=LSIZE5*5
          END IF
          KIP5=ISIZE5*5
        END IF
C
        IF((ICOUPL.GT.5.OR.(NNMAX.LT.0.AND.ICI.LT.-5)).OR.
     1    (ICOUPL.EQ.5.AND.IREACT.NE.0))THEN
C         NMAX1=MAX0(MAXSUM(1,1)+6,MAXSUM(2,1)+6,MAXSUM(3,1)+6,
C    1    MAXSUM(4,1)+6,MAXSUM(5,1)+6,MAXSUM(6,1)+6)
C         NMAX2=MAX0(MAXSUM(1,2)+6,MAXSUM(2,2)+6,MAXSUM(3,2)+6,
C    1    MAXSUM(4,2)+6,MAXSUM(5,2)+6,MAXSUM(6,2)+6)
          NMAX1=0
          NMAX2=0
          DO I=1,MAXJ
            NMAX1=MAX0(NMAX1,MAXSUM(I,1)+6)
            NMAX2=MAX0(NMAX2,MAXSUM(I,2)+6)
          END DO
          NMAX=MAX0(NMAX1,NMAX2)
          ISIZE6=JCI6**6
C**FIRST GET 6-D SIZE FOR 'NORMAL' (ISIZE6 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE6,6,JCI6,0)
          NSIZE6=ISIZE6
          IF(LINEAR.AND..NOT.TRIAT)THEN
            NMAX1=MAX0(MAXSUM(1,1)+5,MAXSUM(2,1)+5,MAXSUM(3,1)+5,
     1      MAXSUM(4,1)+5,MAXSUM(5,1)+5)
            NMAX2=MAX0(MAXSUM(1,2)+5,MAXSUM(2,2)+5,MAXSUM(3,2)+5,
     1      MAXSUM(4,2)+5,MAXSUM(5,2)+5)
            NMAX=MAX0(NMAX1,NMAX2)
            LSIZE6=ISIZE5
C**FIRST GET 5-D SIZE FOR 'LINEAR' (LSIZE6 MODIFIED IN GETSZ FOR NMAX)
            CALL GETSZ(LSIZE6,5,JCI5,1)
            KLJP6=LSIZE6*6
          END IF
          KIP6=ISIZE6*6
        END IF
C*******************************************************
C**GET SIZES OF 1-,2-,3-,4-,5-,6- MODE BASIC BASES - END
C*******************************************************
        NMAX=NNMAX
        IF(NNMAX.GE.0)THEN
C******************************
C**GET FULL CI BASIS HERE IN IP
C******************************
          NVSYM=NVSYMX
          ISIZE=(JCI)**NVMODE
          JSIZE=ISIZE
C**FIRST GET SIZE (ISIZE MODIFIED IN GETSZ FOR NMAX)
          IF(NMAX.GT.0)CALL GETSZ(ISIZE,NVMODE,JCI,0)
          KIP=ISIZE*NVMODE
          KJP=KIP
          KSIZE=ISIZE
          CALL MEMO(2,LIP,KIP,LJP,KJP,0,0,0,0,0,0)
C**TEST IF LINEAR (IF SO, ISIZE MODIFIED IN GETIP FOR TORSION)
          IF(LINEAR.AND..NOT.TRIAT)THEN
            CALL GETIP(W(LIP),W(LJP),ISIZE,ISIZE,JSIZE,NVMODE,JCI,ICI,
     1      -1)
            CALL MEMO(-2,LIP,KIP,LJP,KJP,0,0,0,0,0,0)
C**ISIZE CONTAINS CORRECT (LINEAR) TOTAL
C**NOW NEED NSMODE QUANTA TO INCLUDE TORSION
            MAXSIZ=MAX0(ISIZE,KSIZE)
            ISIZE=MAXSIZ
            KIP=ISIZE*NSMODE
            KJP=MAXSIZ*NSMODE
            CALL MEMO(2,LIP,KIP,LJP,KJP,0,0,0,0,0,0)
            CALL GETIP(W(LIP),W(LJP),ISIZE,MAXSIZ,JSIZE,NVMODE,JCI,ICI,
     1      1)
C**RESET JSIZE
            JSIZE=MAXSIZ
          ELSE
            CALL GETIP(W(LIP),W(LJP),ISIZE,ISIZE,JSIZE,NVMODE,JCI,ICI,
     1      0)
C**RESET JSIZE
            JSIZE=ISIZE
          END IF
          WRITE(IOUT,395)(NTOT(K),K=1,NVSYM)
          IF(MATSIZ.NE.0)THEN
            WRITE(6,*) 'VCI MATRIX SIZE (1)'
            RETURN
          END IF
        ELSE
C************************************
C**GET SELECTIVE CI BASIS SIZES FIRST
C************************************
          ITOT=0
          DO K=1,2
            ITOT1(K)=1
            ITOT2(K)=0
            ITOT3(K)=0
            ITOT4(K)=0
            ITOT5(K)=0
            ITOT6(K)=0
            JSIZE1(K)=0
            JSIZE2(K)=0
            JSIZE3(K)=0
            JSIZE4(K)=0
            JSIZE5(K)=0
            JSIZE6(K)=0
          END DO
          KJP1=0
          KJP2=0
          KJP3=0
          KJP4=0
          KJP5=0
          KJP6=0
C**LOOP ROUND CONTRACTION SCHEMES
          DO K=1,NCONT
            NCMODE=ICONT(K)
            JCI1=NBAS(1,K)+1
            NMAX=MAXSUM(1,K)+1
            KSIZE1=JCI1
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
            CALL GETSZ(KSIZE1,1,JCI1,0)
C**SUBTRACT ZERO POINT BASIS
            KSIZE1=KSIZE1-1
C**ADD 1 TO ONE-MODE MATRIX FOR ZERO POINT BASIS
            JSIZE1(K)=KSIZE1*NCMODE+1
            KJP1=KJP1+JSIZE1(K)*NCMODE
            IF(KJP1.GT.0.AND.K.EQ.NCONT)
     1      CALL MEMO(1,LJP1,KJP1,0,0,0,0,0,0,0,0)
C************************************
            IF(ICI.LT.-1)THEN
              JCI2=NBAS(2,K)+1
              NMAX=MAXSUM(2,K)+1
              KSIZE1=JCI2
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
              CALL GETSZ(KSIZE1,1,JCI2,0)
              KSIZE1=KSIZE1-1
              NMAX=MAXSUM(2,K)+2
              KSIZE2=JCI2*JCI2
C**FIRST GET SIZE (KSIZE2 MODIFIED IN GETSZ FOR NMAX)
              CALL GETSZ(KSIZE2,2,JCI2,0)
              KSIZE2=KSIZE2-2*KSIZE1-1
              JSIZE2(K)=KSIZE2*NCMODE*(NCMODE-1)/2
              KJP2=KJP2+JSIZE2(K)*NCMODE
              IF(KJP2.GT.0.AND.K.EQ.NCONT)
     1        CALL MEMO(1,LJP2,KJP2,0,0,0,0,0,0,0,0)
C************************************
              IF(ICI.LT.-2)THEN
                JCI3=NBAS(3,K)+1
                NMAX=MAXSUM(3,K)+1
                KSIZE1=JCI3
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
                CALL GETSZ(KSIZE1,1,JCI3,0)
                KSIZE1=KSIZE1-1
                NMAX=MAXSUM(3,K)+2
                KSIZE2=JCI3*JCI3
C**FIRST GET SIZE (KSIZE2 MODIFIED IN GETSZ FOR NMAX)
                CALL GETSZ(KSIZE2,2,JCI3,0)
                KSIZE2=KSIZE2-2*KSIZE1-1
                NMAX=MAXSUM(3,K)+3
                KSIZE3=JCI3*JCI3*JCI3
C**FIRST GET SIZE (KSIZE3 MODIFIED IN GETSZ FOR NMAX)
                CALL GETSZ(KSIZE3,3,JCI3,0)
                KSIZE3=KSIZE3-3*KSIZE2-3*KSIZE1-1
                JSIZE3(K)=KSIZE3*NCMODE*(NCMODE-1)*(NCMODE-2)/(3*2)
                KJP3=KJP3+JSIZE3(K)*NCMODE
                IF(KJP3.GT.0.AND.K.EQ.NCONT)
     1          CALL MEMO(1,LJP3,KJP3,0,0,0,0,0,0,0,0)
C************************************
                IF(ICI.LT.-3)THEN
                  JCI4=NBAS(4,K)+1
                  NMAX=MAXSUM(4,K)+1
                  KSIZE1=JCI4
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(KSIZE1,1,JCI4,0)
                  KSIZE1=KSIZE1-1
                  NMAX=MAXSUM(4,K)+2
                  KSIZE2=JCI4*JCI4
C**FIRST GET SIZE (KSIZE2 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(KSIZE2,2,JCI4,0)
                  KSIZE2=KSIZE2-2*KSIZE1-1
                  NMAX=MAXSUM(4,K)+3
                  KSIZE3=JCI4*JCI4*JCI4
C**FIRST GET SIZE (KSIZE3 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(KSIZE3,3,JCI4,0)
                  KSIZE3=KSIZE3-3*KSIZE2-3*KSIZE1-1
                  NMAX=MAXSUM(4,K)+4
                  KSIZE4=JCI4*JCI4*JCI4*JCI4
C**FIRST GET SIZE (KSIZE4 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(KSIZE4,4,JCI4,0)
                  KSIZE4=KSIZE4-4*KSIZE3-6*KSIZE2-4*KSIZE1-1
                  JSIZE4(K)=KSIZE4*NCMODE*(NCMODE-1)*(NCMODE-2)*
     1            (NCMODE-3)/(4*3*2)
                  KJP4=KJP4+JSIZE4(K)*NCMODE
                  IF(KJP4.GT.0.AND.K.EQ.NCONT)
     1            CALL MEMO(1,LJP4,KJP4,0,0,0,0,0,0,0,0)
C************************************
                  IF(ICI.LT.-4)THEN
                    JCI5=NBAS(5,K)+1
                    NMAX=MAXSUM(5,K)+1
                    KSIZE1=JCI5
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
                    CALL GETSZ(KSIZE1,1,JCI5,0)
                    KSIZE1=KSIZE1-1
                    NMAX=MAXSUM(5,K)+2
                    KSIZE2=JCI5*JCI5
C**FIRST GET SIZE (KSIZE2 MODIFIED IN GETSZ FOR NMAX)
                    CALL GETSZ(KSIZE2,2,JCI5,0)
                    KSIZE2=KSIZE2-2*KSIZE1-1
                    NMAX=MAXSUM(5,K)+3
                    KSIZE3=JCI5*JCI5*JCI5
C**FIRST GET SIZE (KSIZE3 MODIFIED IN GETSZ FOR NMAX)
                    CALL GETSZ(KSIZE3,3,JCI5,0)
                    KSIZE3=KSIZE3-3*KSIZE2-3*KSIZE1-1
                    NMAX=MAXSUM(5,K)+4
                    KSIZE4=JCI5*JCI5*JCI5*JCI5
C**FIRST GET SIZE (KSIZE4 MODIFIED IN GETSZ FOR NMAX)
                    CALL GETSZ(KSIZE4,4,JCI5,0)
                    KSIZE4=KSIZE4-4*KSIZE3-6*KSIZE2-4*KSIZE1-1
                    NMAX=MAXSUM(5,K)+5
                    KSIZE5=JCI5*JCI5*JCI5*JCI5*JCI5
C**FIRST GET SIZE (KSIZE5 MODIFIED IN GETSZ FOR NMAX)
                    CALL GETSZ(KSIZE5,5,JCI5,0)
                    KSIZE5=KSIZE5-5*KSIZE4-10*KSIZE3-10*KSIZE2-
     1              5*KSIZE1-1
                    JSIZE5(K)=KSIZE5*NCMODE*(NCMODE-1)*(NCMODE-2)*
     1              (NCMODE-3)*(NCMODE-4)/(5*4*3*2)
                    KJP5=KJP5+JSIZE5(K)*NCMODE
                    IF(KJP5.GT.0.AND.K.EQ.NCONT)
     1              CALL MEMO(1,LJP5,KJP5,0,0,0,0,0,0,0,0)
C************************************
                    IF(ICI.LT.-5)THEN
                      JCI6=NBAS(6,K)+1
                      NMAX=MAXSUM(6,K)+1
                      KSIZE1=JCI6
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
                      CALL GETSZ(KSIZE1,1,JCI6,0)
                      KSIZE1=KSIZE1-1
                      NMAX=MAXSUM(6,K)+2
                      KSIZE2=JCI6*JCI6
C**FIRST GET SIZE (KSIZE2 MODIFIED IN GETSZ FOR NMAX)
                      CALL GETSZ(KSIZE2,2,JCI6,0)
                      KSIZE2=KSIZE2-2*KSIZE1-1
                      NMAX=MAXSUM(6,K)+3
                      KSIZE3=JCI6*JCI6*JCI6
C**FIRST GET SIZE (KSIZE3 MODIFIED IN GETSZ FOR NMAX)
                      CALL GETSZ(KSIZE3,3,JCI6,0)
                      KSIZE3=KSIZE3-3*KSIZE2-3*KSIZE1-1
                      NMAX=MAXSUM(6,K)+4
                      KSIZE4=JCI6*JCI6*JCI6*JCI6
C**FIRST GET SIZE (KSIZE4 MODIFIED IN GETSZ FOR NMAX)
                      CALL GETSZ(KSIZE4,4,JCI6,0)
                      KSIZE4=KSIZE4-4*KSIZE3-6*KSIZE2-4*KSIZE1-1
                      NMAX=MAXSUM(6,K)+5
                      KSIZE5=JCI6*JCI6*JCI6*JCI6*JCI6
C**FIRST GET SIZE (KSIZE5 MODIFIED IN GETSZ FOR NMAX)
                      CALL GETSZ(KSIZE5,5,JCI6,0)
                      KSIZE5=KSIZE5-5*KSIZE4-10*KSIZE3-10*KSIZE2-
     1                5*KSIZE1-1
                      NMAX=MAXSUM(6,K)+6
                      KSIZE6=JCI6*JCI6*JCI6*JCI6*JCI6*JCI6
C**FIRST GET SIZE (KSIZE6 MODIFIED IN GETSZ FOR NMAX)
                      CALL GETSZ(KSIZE6,6,JCI6,0)
                      KSIZE6=KSIZE6-6*KSIZE5-15*KSIZE4-20*KSIZE3-
     1                15*KSIZE2-6*KSIZE1-1
                      JSIZE6(K)=KSIZE6*NCMODE*(NCMODE-1)*(NCMODE-2)*
     1                (NCMODE-3)*(NCMODE-4)*(NCMODE-5)/(6*5*4*3*2)
                      KJP6=KJP6+JSIZE6(K)*NCMODE
                      IF(KJP6.GT.0.AND.K.EQ.NCONT)
     1                CALL MEMO(1,LJP6,KJP6,0,0,0,0,0,0,0,0)
C************************************
                    END IF
                  END IF
                END IF
              END IF
            END IF
C**END LOOP ROUND CONTRACTION SCHEMES
          END DO
C**GET SELECTIVE CI BASES 
          K1=0
          K2=0
          K3=0
          K4=0
          K5=0
          K6=0
C**LOOP ROUND CONTRACTION SCHEMES
          DO KBAS=1,NCONT
            NCMODE=ICONT(KBAS)
            DO K=1,NCMODE
              KK=JCONT(KBAS,K)
              JCI1=NBAS(1,KBAS)+1
              NMAX=MAXSUM(1,KBAS)+1
              IF(JSIZE1(KBAS).GT.0)
     1        CALL GETJP1(W(LJP1+K1),JSIZE1(KBAS),NMODE,NVMODE,K,JCI1,
     2        NMAX,ITOT1(KBAS),W(LMXBAS),NCMODE,KK)
              IF(ICI.GT.-2)GO TO 6601
              DO L=1,K-1
                LL=JCONT(KBAS,L)
                JCI2=NBAS(2,KBAS)+1
                NMAX=MAXSUM(2,KBAS)+2
                IF(JSIZE2(KBAS).GT.0)
     1          CALL GETJP2(W(LJP2+K2),JSIZE2(KBAS),NMODE,NVMODE,
     2          NVMODE,K,L,JCI2,JCI2,NMAX,ITOT2(KBAS),W(LMXBAS),NCMODE,
     3          KK,LL)
                IF(ICI.GT.-3)GO TO 6602
                DO N=1,L-1
                  NN=JCONT(KBAS,N)
                  JCI3=NBAS(3,KBAS)+1
                  NMAX=MAXSUM(3,KBAS)+3
                  IF(JSIZE3(KBAS).GT.0)
     1            CALL GETJP3(W(LJP3+K3),JSIZE3(KBAS),NMODE,NVMODE,
     2            NVMODE,NVMODE,K,L,N,JCI3,JCI3,JCI3,NMAX,ITOT3(KBAS),
     3            W(LMXBAS),NCMODE,KK,LL,NN)
                  IF(ICI.GT.-4)GO TO 6603
                  DO M=1,N-1
                    MM=JCONT(KBAS,M)
                    JCI4=NBAS(4,KBAS)+1
                    NMAX=MAXSUM(4,KBAS)+4
                    IF(JSIZE4(KBAS).GT.0)
     1              CALL GETJP4(W(LJP4+K4),JSIZE4(KBAS),NMODE,NVMODE,
     2              NVMODE,NVMODE,NVMODE,K,L,N,M,JCI4,JCI4,JCI4,JCI4,
     3              NMAX,ITOT4(KBAS),W(LMXBAS),NCMODE,KK,LL,NN,MM)
                    IF(ICI.GT.-5)GO TO 6604
                    DO I=1,M-1
                      II=JCONT(KBAS,I)
                      JCI5=NBAS(5,KBAS)+1
                      NMAX=MAXSUM(5,KBAS)+5
                      IF(JSIZE5(KBAS).GT.0)
     1                CALL GETJP5(W(LJP5+K5),JSIZE5(KBAS),NMODE,NVMODE,
     2                NVMODE,NVMODE,NVMODE,NVMODE,K,L,N,M,I,JCI5,JCI5,
     3                JCI5,JCI5,JCI5,NMAX,ITOT5(KBAS),W(LMXBAS),
     4                NCMODE,KK,LL,NN,MM,II)
                      IF(ICI.GT.-6)GO TO 6605
                      DO J=1,I-1
                        JJ=JCONT(KBAS,J)
                        JCI6=NBAS(6,KBAS)+1
                        NMAX=MAXSUM(6,KBAS)+6
                        IF(JSIZE6(KBAS).GT.0)
     1                  CALL GETJP6(W(LJP6+K6),JSIZE6(KBAS),NMODE,
     2                  NVMODE,NVMODE,NVMODE,NVMODE,NVMODE,NVMODE,K,L,
     3                  N,M,I,J,JCI6,JCI6,JCI6,JCI6,JCI6,JCI6,NMAX,
     4                  ITOT6(KBAS),W(LMXBAS),NCMODE,KK,LL,NN,MM,II,
     5                  JJ)
                      END DO
6605  CONTINUE
                    END DO
6604  CONTINUE
                  END DO
6603  CONTINUE
                END DO
6602  CONTINUE
              END DO
6601  CONTINUE
            END DO
            K1=K1+JSIZE1(1)*NCMODE
            K2=K2+JSIZE2(1)*NCMODE
            K3=K3+JSIZE3(1)*NCMODE
            K4=K4+JSIZE4(1)*NCMODE
            K5=K5+JSIZE5(1)*NCMODE
            K6=K6+JSIZE6(1)*NCMODE
C**END LOOP ROUND CONTRACTION SCHEMES
          END DO
        END IF
      END IF
C*****************
C**SKIP IF VSCF-CI
C*****************
      IF(ICI.GE.0)GO TO 4500
      KEL21=J21
      IF(JMAX.NE.0)THEN
        KEL21=0
        DO KROT=1,J21,KSTEP
          KEL21=KEL21+1
        END DO
      END IF
C**********************************************************************
C**********************************************************************
C**                                      LOOP ROUND CONTRACTION SCHEMES
C**********************************************************************
C**********************************************************************
      LCOUNT=1
      IF(KCONT.GT.0)LCOUNT=NCONT
C**GET SPACE FOR FUTURE ASSIGNMENTS (MAX NVAL1,NVAL2 PER SCHEME)
      KSCFCI=MAX0(NVAL1,NVAL2)
      KCASS=KSCFCI*KEL21*NVSYMX*LCOUNT
      CALL MEMO(1,LCASS,KCASS,0,0,0,0,0,0,0,0)
C**SAVE IMPORTANT QUANTITIES
      NNMODS=NNMODE
      NAMODS=NAMODE
      NVMODS=NVMODE
C     ICOUPS=ICOUPL
C     ICOUCS=ICOUPC
      JCOUPS=JCOUPL
      JCOUCS=JCOUPC
      JREACS=JREACT
      NREACS=NREACT
C**LOOP ROUND CONTRACTION SCHEMES (NON-LINEAR ONLY AT PRESENT)
      DO 9999 LCONT=1,LCOUNT
      WRITE(IOUT,495)LCONT
C*********************************************
C*********************************************
C**    PERFORM VCI FOR INDIVIDUAL CONTRACTIONS
C*********************************************
C*********************************************
      IF(LCOUNT.EQ.2)THEN
        JCOUPL=NCOUPL(LCONT)
        JCOUPC=NCOUPC(LCONT)
        ICOUPL=IABS(JCOUPL)
        ICOUPC=IABS(JCOUPC)
        IF(JCOUPL.GE.0)JCOUPC=ICOUPC
        IF(JCOUPL.LT.0)JCOUPC=-ICOUPC
      END IF
      IF(KCONT.LT.0)THEN
        ICOUPL=ICOUPS
        ICOUPC=ICOUCS
      END IF
      IF(ICOUPL.GT.5.AND.IWHICH.LT.0)
     1STOP 'FITTED POTENTIAL MUST HAVE ICOUPL < 6'
      IF(ICOUPC.GT.ICOUPL)STOP 'ICOUPC CANNOT BE GREATER THAN ICOUPL'
      IF(ICOUPC.GT.6)STOP 'ICOUPC CANNOT EXCEED 6'
      CALL VCI(W,KCONT,LCONT,NNMAX,KLP,KSCFCI,NVSYMX,EVLJ0,K2TAU,K3TAU,
     1W(LMXBAS),NMODE)
      IF(KCONT.LT.0)THEN
        ICOUPL=ICOUPX
        ICOUPC=ICOUCX
      END IF
9999  CONTINUE
      IF(MATSIZ.NE.0)STOP 'VCI MATRIX SIZE(2)'
      IF(MVAL1.LT.0.OR.MVAL2.LT.0)THEN
        WRITE(IOUT,225)
        RETURN
      END IF
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                                    END OF VCI STAGE
C**********************************************************************
C**********************************************************************
C**********************************************************************
      CALL MEMO(-1,LXW,KXW,0,0,0,0,0,0,0,0)
C**RESTORE IMPORTANT QUANTITIES
      NAMODE=NAMODS
      NNMODE=NNMODS
      NVMODE=NVMODS
      ICOUPL=ICOUPS
      ICOUPC=ICOUCS
      JCOUPL=JCOUPS
      JCOUPC=JCOUCS
      JREACT=JREACS
      NREACT=NREACS
      ISIZXX=ISIZMX
      LCOUNT=-LCOUNT+1
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                        COMBINE TWO CONTRACTION SCHEMES IF REQUIRED
C**********************************************************************
C**********************************************************************
C**********************************************************************
      IF(LCOUNT.LT.0)CALL COMBIN(W,KLP,KSCFCI,EVLJ0)
      IF(JTHIS.EQ.0)GO TO 4444
      IF(KMAX.LT.0)GO TO 4444
      KIP=ISIZMX*NNMODE
      CALL MEMO(2,LIPL,KIP,LIPR,KIP,0,0,0,0,0,0)
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                                     VCI - ROTATIONS
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**RESTORE NUMBER OF CONTRACTION SCHEMES FOR ROTATIONS
      LCOUNT=1
      IF(KCONT.GT.0)LCOUNT=NCONT
      IF(ICOUPL.GT.4)ICOUPL=4
      IF(ICOUPC.GT.4)ICOUPC=4
C     IF(ICOUPL.GT.3)ICOUPL=3
C     IF(ICOUPC.GT.3)ICOUPC=3
      CALL VIBROT(W)
4444  CONTINUE
      CALL FLUSH(IOUT)
C**ANALYTIC
      IF((IWHICH.LT.0.AND.MOLINC.GT.0).OR.IFITRP.NE.0)
     1CALL MEMO(-1,LXKAN,MAXQU*MAXQU*NAMODE*MAXPOW,0,0,0,0,0,0,0,0)
C**ANALYTIC
      GO TO 5000
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                                             VSCF-CI
C**********************************************************************
C**********************************************************************
C**********************************************************************
4500  CONTINUE
      CALL VSCFCI(W,NSTAT,EVLJ0)
5000  CONTINUE
      IF(JTHIS.NE.0)WRITE(IOUT,240)EVL
      WRITE(IOUT,225)
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMMVS(W,MAXBAS,NMODE,NUMJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRIAT
      DIMENSION W(1),MAXBAS(NMODE,6,2)
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
     1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
     4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
     5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
     6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
     7LDP4A,LDP4B,LDP4C,LDP4D,LXP,LJP6,LANCNT,LANK,LXJ,LVP5,
     8LC5,LC6,LVP0,LTMPRD,LMVB,LNFC,LTEMP6,LTEMP7,LTEMP8,LTEMP9
C**TEMPORARY (DIMENSIONS)
C**RPH ODD K
      COMMON/REACTL/JREACT,IFITRP
C**RPH ODD K
      COMMON/TRIATO/TRIAT
      COMMON/FSYMM/NFTOT(10,2),NRTOT(10,2,2)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MXSYM,MWSYM(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/FCONTS/NCONT(2)
      COMMON/FBASIS/NBAS(6,2,2),MAXSUM(6,2,2)
      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7
      DO JJ=1,NUMJ
        MAXJJ=JJ-1
C**************************************************************
C**************************************************************
        CALL MEMO(2,LNFC,NMODE*2,LMVB,NMODE,0,0,0,0,0,0)
        READ(INP4)MMODE
        IF(MMODE.NE.NMODE)STOP 'WRONG MOLECULE'
        DO MODE=1,NMODE
C**NUMBER CONTRACTED AND PRIMITIVE FUNCTIONS (INPUT)
          READ(INP4)NVV,NBB
          CALL PUTINT(W(LNFC),MODE,NVV)
C**OMEGA AND NUMBER INTEGRATION POINTS (INPUT)
          READ(INP4)YLAM,MBB
          CALL PUTINT(W(LMVB),MODE,MBB)
          CALL MEMO(1,LTMPRD,MBB,0,0,0,0,0,0,0,0)
C**INTEGRATION POINTS (INPUT)
          CALL INP60R(MBB,W(LTMPRD),INP4)
          CALL MEMO(-1,LTMPRD,MBB,0,0,0,0,0,0,0,0)
        END DO
C**********************************
        DO MODE=1,NMODE
          CALL GETINT(W(LNFC),MODE,NVV)
CC        NVV=NFC(MODE)
          CALL GETINT(W(LMVB),MODE,MBB)
CC        MBB=MVB(MODE)
C**RPH ODD K
          INDMX=1
          IF(JREACT.GT.0.AND.MODE.EQ.NMODE.AND.MAXJJ.GT.0)INDMX=2
          DO IND=1,INDMX
            DO K=1,3
              DO M=1,MBB
                CALL MEMO(1,LTMPRD,NVV,0,0,0,0,0,0,0,0)
C**VSCF/VCI FUNCTIONS, 1ST AND 2ND DERIVATIVES AT INTEGRATION POINTS
**WITH SQRT(WEIGHT) INCLUDED (INPUT)
                CALL IN60(NVV,W(LTMPRD),INP4)
                CALL MEMO(-1,LTMPRD,NVV,0,0,0,0,0,0,0,0)
              END DO
            END DO
          END DO
C**RPH ODD K
        END DO
C**********************************
        READ(INP4)JTHIS,KSTEP,MVSYM,
     1  (MWSYM(I),I=1,MVSYM),ICI,NMAX,NCONT(1),
     2  ((NBAS(I,J,1),I=1,6),J=1,2),
     3  ((MAXSUM(I,J,1),I=1,6),J=1,2),
     4  ((MAXBAS(I,J,1),I=1,NMODE),J=1,IABS(NMAX))
        J21=2*JTHIS+1
C**TEMPORARY
      WRITE(IOUT,*)'JTHIS,MVSYM = ',JTHIS,MVSYM
      CALL FLUSH(IOUT)
        MMSS=MVSYM
        IF(JTHIS.NE.0)MMSS=NVSYM
        DO IYSYM=1,MMSS
C       DO IYSYM=1,MVSYM
C**TEMPORARY
          IVSYM=MWSYM(IYSYM)
C**SIZE OF VCI MATRIX (INPUT)
          READ(INP4)JSIZE
          NFTOT(IVSYM,1)=JSIZE
          IF(JSIZE.EQ.0)GO TO 1112
          DO I=1,JSIZE
            CALL MEMO(1,LTMPRD,NMODE,0,0,0,0,0,0,0,0)
            CALL INP60I(NMODE,W(LTMPRD),INP4)
            CALL MEMO(-1,LTMPRD,NMODE,0,0,0,0,0,0,0,0)
          END DO
C**LOOP OVER Ka
          DO KROT=1,J21,KSTEP
C**NUMBER OF CI FUNCTIONS (INPUT)
            READ(INP4)JDUMP
            DO J=1,JDUMP
C**SPECIFIC CI FUNCTION AND ENERGY (INPUT)
              READ(INP4)K
              IF(K.GE.0)THEN
                READ(INP4)ENERGY
                CALL MEMO(1,LTMPRD,2+NMODE,0,0,0,0,0,0,0,0)
C**ASSIGNMENTS FOR Jth EIGENFUNCTION (INPUT)
                CALL INP60A(NMODE,1,1,W(LTMPRD),1,1,1,INP4)
                CALL MEMO(-1,LTMPRD,2+NMODE,0,0,0,0,0,0,0,0)
                CALL MEMO(1,LTMPRD,JSIZE,0,0,0,0,0,0,0,0)
C**CI COEFFICIENTS FOR Jth EIGENFUNCTION (INPUT)
                CALL INP60R(JSIZE,W(LTMPRD),INP4)
                CALL MEMO(-1,LTMPRD,JSIZE,0,0,0,0,0,0,0,0)
              END IF
            END DO
          END DO
1112      CONTINUE
        END DO
        CALL MEMO(-2,LNFC,NMODE*2,LMVB,NMODE,0,0,0,0,0,0)
C
        WRITE(IOUT,*)'READING VIBRATIONAL DUMMY DATA FOR J = ',JTHIS
        IF(JTHIS.NE.MAXJJ)THEN
          WRITE(IOUT,*)'OUTPUT FOR J= ',MAXJJ,' IS MISSING'
          STOP 'WRONG J'
        END IF
        CALL FLUSH(IOUT)
C
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMMRS(W,NUMJ,NMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(1)
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
     1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
     4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
     5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
     6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
     7LDP4A,LDP4B,LDP4C,LDP4D,LXP,LJP6,LANCNT,LANK,LXJ,LVP5,
     8LC5,LC6,LVP0,LTMPRD,LMVB,LNFC,LTEMP6,LTEMP7,LTEMP8,LTEMP9
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7
      DO JJ=1,NUMJ
        MAXJJ=JJ-1
        IF(MAXJJ.EQ.0)GO TO 9999
        DO IYSYM=1,MVSYM
C**SIZE OF VCI MATRIX
          READ(INP5)JSIZE
          IF(JSIZE.EQ.0)GO TO 1113
C**SIZE OF EACH Ka BLOCK
C**NUMBER OF CI FUNCTIONS
          READ(INP5)NVAL,JDUMM,JDUMP
          J21=JSIZE/NVAL
          DO J=1,JDUMP
C**SPECIFIC CI FUNCTION AND ENERGY (INPUT)
            READ(INP5)K
            IF(K.NE.0)THEN
              READ(INP5)ENERGY
              CALL MEMO(1,LTMPRD,2+NMODE,0,0,0,0,0,0,0,0)
C**ASSIGNMENTS FOR Jth EIGENFUNCTION (INPUT)
              CALL INP60A(NMODE,1,1,W(LTMPRD),1,1,1,INP5)
              CALL MEMO(-1,LTMPRD,2+NMODE,0,0,0,0,0,0,0,0)
              DO IK=1,J21
                CALL MEMO(1,LTMPRD,NVAL,0,0,0,0,0,0,0,0)
                CALL INP60R(NVAL,W(LTMPRD),INP5)
                CALL MEMO(-1,LTMPRD,NVAL,0,0,0,0,0,0,0,0)
              END DO
            END IF
          END DO
1113      CONTINUE
        END DO
C
        WRITE(IOUT,*)'READING ROVIBRATIONAL DUMMY DATA FOR J = ',MAXJJ
        CALL FLUSH(IOUT)
C
9999    CONTINUE
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************

C**TEMPORARY
      SUBROUTINE NEWMIN(X0,XL,XM,NATOM,NMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FILASS/IOUT
      DIMENSION X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      WRITE(IOUT,*)'OLD GEOMETRY'
      DO I=1,NATOM
        WRITE(IOUT,*)(X0(I,J),J=1,3)
      END DO
      WRITE(IOUT,*)'NUMBER ATOMS = ',NATOM
      WRITE(IOUT,*)'NUMBER MODES = ',NMODE
      WRITE(IOUT,*)'MASSES'
      WRITE(IOUT,*)(XM(I),I=1,NATOM)
      WRITE(IOUT,*)'VECTORS'
      DO K=1,NMODE
        WRITE(IOUT,*)'MODE = ',K
        DO I=1,NATOM
          WRITE(IOUT,*)(XL(I,K,J),J=1,3)
        END DO
      END DO
      WRITE(IOUT,*)
      WRITE(IOUT,*)'NEW GEOMETRY'
      DO I=1,NATOM
        DO J=1,3
C         X0(I,J)=X0(I,J)+XL(I,1,J)*56.9D0/SQRT(XM(I))
C         X0(I,J)=X0(I,J)-XL(I,2,J)*0.0006D0/SQRT(XM(I))
C         X0(I,J)=X0(I,J)+XL(I,3,J)*0.0002D0/SQRT(XM(I))
C         X0(I,J)=X0(I,J)+XL(I,4,J)*8.9848D0/SQRT(XM(I))
C         X0(I,J)=X0(I,J)+XL(I,5,J)*3.5225D0/SQRT(XM(I))
C         X0(I,J)=X0(I,J)+XL(I,6,J)*3.7399D0/SQRT(XM(I))
C         X0(I,J)=X0(I,J)-XL(I,7,J)*0.8294D0/SQRT(XM(I))
C         X0(I,J)=X0(I,J)-XL(I,8,J)*0.3777D0/SQRT(XM(I))
          X0(I,J)=X0(I,J)+XL(I,9,J)*10.6682D0/SQRT(XM(I))
CCCC      X0(I,J)=X0(I,J)+XL(I,9,J)*19.6682D0/SQRT(XM(I))
        END DO
        WRITE(IOUT,*)(X0(I,J),J=1,3)
      END DO
      RETURN
      END
C**TEMPORARY

C**************************************************************
C*************************************************TEST ROUTINE
      SUBROUTINE HERMTN(X,M,N1,N2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(M,M)
      COMMON/FILASS/IOUT
      WRITE(IOUT,8888)
8888  FORMAT(1X,'HERMITIAN TEST',/)
      DO 1 I=1,N1
      WRITE(IOUT,9999)(X(I,J),J=1,I)
      WRITE(IOUT,9999)(X(J,I),J=1,I)
1     CONTINUE
9999  FORMAT(1X,5D14.5)
      TOL=1.D-10
      WRITE(IOUT,7777)TOL
7777  FORMAT(//,1X,'TOLERANCE TEST OF ',D14.5,/)
      DO 2 I=1,N2
      DO 2 J=1,I
      IF(I.EQ.J)GO TO 2
      IF(DABS(X(I,J)-X(J,I)).GT.TOL)THEN
      WRITE(IOUT,5555)I,J,X(I,J),X(J,I)
5555  FORMAT(1X,2I5,2D14.5)
      END IF
2     CONTINUE
      RETURN
      END
C*************************************************TEST ROUTINE
      SUBROUTINE HMOVE(H,HX,NH,NHX,MH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(NH,MH,3),HX(NHX,MH,3)
      DO K=1,3
        DO M=1,MH
          DO N=1,NHX
            HX(N,M,K)=H(N,M,K)
          END DO
        END DO
      END DO
      RETURN
      END
C*************************************************TEST ROUTINE
      SUBROUTINE QMOVE(XQ,XQX,MM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XQ(MM),XQX(MM)
      DO M=1,MM
        XQX(M)=XQ(M)
      END DO
      RETURN
      END
C*************************************************TEST ROUTINE
      SUBROUTINE NMOVE(NBF,K,NVF)
      DIMENSION NBF(1)
      NBF(K)=NVF
      RETURN
      END
C*************************************************TEST ROUTINE
      INTEGER FUNCTION INTFUN(NOCON,K)
      DIMENSION NOCON(K)
      INTFUN=NOCON(K)
      RETURN
      END
C*************************************************TEST ROUTINE
C**************************************************************
      SUBROUTINE VIBROT(W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 TITLE(80)
      LOGICAL LGIV,LINEAR,LANCZ,LANZA,LANZB,TRIAT
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
     1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
     4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
     5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
     6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
     7LDP4A,LDP4B,LDP4C,LDP4D,LXP,LJP6,LANCNT,LANK,LXJ,LVP5,
     8LC5,LC6,LVP0,LTMPRD,LMVB,LNFC,LTEMP6,LTEMP7,LTEMP8,LTEMP9
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/TITLE/TITLE
      COMMON/HERM/IHERM
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/PATH/ISCFCI
      COMMON/DEGEN/IDEGEN,NDEGEN
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU,IFLAUD
      COMMON/REACTL/JREACT,IFITRP
      COMMON/NREACT/NREACT
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/CYCLES/NCYCLE
      COMMON/TRIATO/TRIAT
      COMMON/LANTOL/TOLLAN
C**TEMPORARY (DIMENSIONS)
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LAN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TOLCUT,EVLCUT,TOLMAT,LANCYC
      COMMON/GIVEN/LGIV,IGIV
      COMMON/TYPE/LINEAR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/TOLS/TOL,EPS
      COMMON/EVL/EVL,CUT
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR(6)
      COMMON/VCIMAX/NMAX
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/ESTATE/IORDER
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/AXES/MX(3),MXDIP(5),MXROT(3),ISYMT
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/TBASIS/NTBAS(6,2),NTAU(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10),NSYMSL(10),MSYMSL(10,100),
     1NSYMSR(10,100),MSYMSR(10,100,10)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CSAVES/NNMODS,NAMODS,NVMODS,ICOUPS,ICOUCS,JREACS,NREACS
      COMMON/NCREC/NREC(6),MAXBUF(6),NUNITR,NUNITW
      COMMON/CSIZES/ISIZM1,ISIZM2,NVAL1,NVAL2,ICSIZ1,ICSIZ2,
     1IPSIZ1,IPSIZ2
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/MAXPT/MBFMAX,MBFMX1,MBFMX2,MBFMX3,MBFMX4,MBFMIN
      COMMON/DISC/IDISC
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/UNITEX/I75,I76,I85,I86
      COMMON/MODES/NMODE,NATOM
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/SIZES/KTEMP,ISIZE1,ISIZE2,ISIZE3,ISIZE4,ISIZE5,ISIZE6,
     1ISIZE,JSIZE,ISIZMX
      COMMON/KLSZS/KLJP1,KLJP2,KLJP3,KLJP4,KLJP5,KLJP6
      COMMON/KJSZS/KJP1,KJP2,KJP3,KJP4,KJP5,KJP6
      COMMON/KPSZS/KIP1,KIP2,KIP3,KIP4,KIP5,KIP6,JCC1,JCC2,JCC3,JCC4,
     1JCC5,JCC6
      COMMON/IPSZS/KPPP1,KPPP2,KPPP3,KPPP4,KPPP5,KPPP6
      COMMON/MATRIX/NVAL,NVALR,KSTEP,KSIGN,NVALCF
C**********************************************************************
100   FORMAT(80A1)
200   FORMAT(//,1X,80A1)
205   FORMAT(/,1X,'NUMBER OF MODES = ',I4,/,
     1         1X,'NUMBER OF SCF STATES = ',I4,/,
     2         1X,'NUMBER POTENTIAL TERMS = ',I4,/,
     3         1X,'PRINT LEVEL = ',I4,/,
     4         1X,'TOLERANCE FOR SCF = ',D20.12,/,
     5         1X,'COUPLING OF ',I4,' MODES',/)
210   FORMAT(/,1X,'IDISC = ',I4,/,
     1         1X,'(WRITE POTENTIAL AND CORIOLIS INFO. TO DISC (0)',/,
     2         1X,'POTENTIAL AND CORIOLIS DISCS ALREADY EXIST (1))',/)
215   FORMAT(/,1X,'SCF CALCULATION ONLY',/)
220   FORMAT(//,1X,'START SCF CALCULATION',/)
225   FORMAT(//,1X,'SCF PLUS CI CALCULATION')
230   FORMAT(/,1X,'NUMBER OF CI ENERGIES REQUIRED = ',I4,/)
235   FORMAT(//,1X,'START CI CALCULATION',/)
240   FORMAT(/,1X,'ZERO POINT ENERGY = ',F10.2,/)
245   FORMAT(//,1X,'OVERLAPS OF ORIGINAL SCF STATE FUNCTIONS',/)
250   FORMAT(/,1X,'OVERLAPS WITH STATE ',I4)
255   FORMAT(//,1X,'TEST OF SCHMIDT ORTHOGONALISATION',/)
260   FORMAT(//,1X,'FINAL VIBRATIONAL (K-DIAGONAL) CI ENERGIES',/)
265   FORMAT(//,1X,'NORMAL COORDINATE POTENTIAL',/)
270   FORMAT(//,1X,'INTERNAL COORDINATE POTENTIAL',/)
275   FORMAT(//,1X,'CI MATRIX INVOLVES ',I4,' SCF FUNCTIONS',/)
280   FORMAT(//,1X,'CI MATRIX INVOLVES ',I4,' VIRTUAL GS FUNCTIONS',/,
     1          1X,'WITH SUM OF QUANTA < ',I4,/)
285   FORMAT(//,1X,'SIZE OF VIBRATIONAL CI MATRIX IS ',I6,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
290   FORMAT(//,1X,'VIBRATIONAL SYMMETRY ',I2)
295   FORMAT(//,1X,'NUMBER OF VIBRATIONAL SYMMETRIES = ',I3,/,
     1          1X,'NUMBER OF MODE SYMMETRIES = ',I3,/)
300   FORMAT(1X,'MODE SYMMETRY ',I3,'   MODES ',20I3)
305   FORMAT(//,1X,'TOTAL ANGULAR MOMENTUM J = ',I3,/)
310   FORMAT(/,1X,'J = 0 SCF CYCLE',/)
315   FORMAT(/,1X,'J = ',I3,' SCF CYCLE',/)
320   FORMAT('*************************')
325   FORMAT(50(1H*))
350   FORMAT(/,1X,'LINK TO PROGRAM NORMALS',/)
360   FORMAT(/,1X,'NUMBER OF LANCZOS CYCLES = ',I3)
365   FORMAT(/,1X,'LANCZOS TOLERANCE = ',F10.6)
370   FORMAT(/,1X,'LANCZOS (HALF) MATRIX MAXIMUM ORDER = ',I6,/)
375   FORMAT(//,1X,'SEQUENCE ',I2,' ROTATIONAL SYMMETRY ',I2,/)
380   FORMAT(//,1X,'SIZE OF RO-VIBRATIONAL CI MATRIX IS ',I5,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
385   FORMAT(//,1X,'FINAL RO-VIBRATIONAL CI ENERGIES',/)
390   FORMAT(1X,'BLOCK ',I2,' SEQUENCE ',I2,' VIBRATIONAL SYMMETRY ',I2)
395   FORMAT(//,1X,'SIZES OF CI SYMMETRY BLOCKS: ',/,1X,10I10,/)
400   FORMAT(//,1X,'SHOULD NOT OCCUR',/)
405   FORMAT(//,1X,'MINIMIZATION OF POTENTIAL',/)
410   FORMAT(//,1X,'DEGREE OF POLYNOMIAL: ',I3,/)
415   FORMAT(//,1X,60(1H*),/,1X,'START CONTRACTION SCHEME ',I2,/,1X,
     160(1H*),//)
420   FORMAT(//,1X,'CONTRACTION SCHEME ',I1,/,
     11X,'ADJUSTED VALUES NVAL = ',10I6,/)
425   FORMAT(//,1X,'CONTRACTION SCHEME ',I1,/,
     11X,'ISIZE = ',10I6,/)
500   FORMAT(//,1X,60(1H*),/,1X,'COMBINE CONTRACTION SCHEMES 1 AND 2',
     1/,1X,60(1H*),//)
C**********************************************************************
C**********************************************************************
C**                                                     VCI - ROTATIONS
C**********************************************************************
C**********************************************************************
C**TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY
C**TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY
CC    CALL ONEMOD(W(LMODNT),NMODE)
C**TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY
C**TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY
C***************************************************
C***************************************************
C**INITIALLY NO LANCZOS
        LANCZ=.FALSE.
C***************************************************
C***************************************************
C**JREACT DIFFERS BETWEEN CONTRACTION SCHEMES
      IF(JREACT.GT.0)THEN
C**GET POINTERS FOR TAU FUNCTIONS
        K2TAU=0
        K3TAU=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          K2TAU=K2TAU+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(K.GT.NONLIN)K2TAU=K2TAU+IK1
          K3TAU=K3TAU+IK2
        END DO
        CALL INTARR(W(LNBF),W(LMBF),NSMODE,IK1TAU,IK2TAU,IK3TAU)
        IF(JCOUPC.GE.0)THEN
          KVM0=12
        ELSE
          KVM0=(1+12)/2
        END IF
        CALL MEMO(1,LVM0,KVM0,0,0,0,0,0,0,0,0)
        IF(IDISC.EQ.0)THEN
C**GET TORSION-ONLY GRID DATA
          CALL DUMVR0(NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),
     1    W(LB),W(LB+NMODE*NMODE),W(LAA),W(LBB),W(LBB+NMODE),W(LXX),
     2    W(LXX+NATOM*3*722),W(LX0),W(LXL),W(LXM),W(LC0),W(LC0),
     3    W(LVM0),W(LVM0),1,W(LVP0),1,W(LMODNT))
        END IF
      END IF

      IF(LDUMP.EQ.0.AND.(JREACT.LE.0.OR.NSMODE.EQ.NAMODE))
     1CALL DUMCR0(NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),
     2W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),1,1)
      IF(ICOUPC.EQ.0)GO TO 3000
      IF(JCOUPC.GE.0)THEN
        KVM1=MBFMX1*7
        IF(JREACT.NE.0)KVM1=MBFMX1*15
      ELSE
        KVM1=(1+MBFMX1*7)/2
        IF(JREACT.NE.0)KVM1=(1+MBFMX1*15)/2
      END IF
      CALL MEMO(1,LVM1,KVM1,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**WRITE CORIOLIS GRIDS TO DISC (91)
          IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
            CALL DUMCR1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),
     1      W(LB),W(LAA),W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,W(LC1),
     2      W(LC1),W(LVM1),W(LVM1),1,1,W(LMODNT))
          ELSE
            CALL DUMVR1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LXZ),
     1      W(LAB),W(LB),W(LB+NMODE*NMODE),W(LAA),W(LBB),
     2      W(LBB+NMODE),W(LXX),W(LXX+NATOM*3*722),W(LX0),W(LXL),
     3      W(LXM),K,W(LC1),W(LC1),W(LVM1),W(LVM1),1,W(LXQ+K3TAU),
     4      W(LVP0),W(LVP1),KVP,1,W(LMODNT))
          END IF
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPC.EQ.1)GO TO 3000
      IF(JCOUPC.GE.0)THEN
        KVM2=MBFMX2*13
        IF(JREACT.NE.0)KVM2=MBFMX2*18
      ELSE
        KVM2=(1+MBFMX2*13)/2
        IF(JREACT.NE.0)KVM2=(1+MBFMX2*18)/2
      END IF
      CALL MEMO(1,LVM2,KVM2,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**WRITE CORIOLIS GRIDS TO DISC (92)
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1      (NSMODE.EQ.NAMODE))THEN
              CALL DUMCR2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1        W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LXX),W(LX0),
     2        W(LXL),W(LXM),K,L,W(LC2),W(LC2),W(LVM2),W(LVM2),1,1,
     3        W(LMODNT),W(LNBF),W(LMBF),W(LVP1),W(LV1),W(LV1),1)
            ELSE
              CALL DUMVR2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1        W(LQQ),W(LXZ),W(LAB),W(LB),W(LB+NMODE*NMODE),W(LAA),
     2        W(LBB),W(LBB+NMODE),W(LXX),W(LXX+NATOM*3*722),W(LX0),
     3        W(LXL),W(LXM),K,L,W(LC2),W(LC2),W(LVM2),W(LVM2),1,
     4        W(LXQ+K3TAU),W(LVP0),W(LVP1),W(LVP2),KVP,1,W(LMODNT),
     5        W(LNBF),W(LMBF),W(LC1),W(LC1),IK2TAU)
            END IF
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPC.EQ.2)GO TO 3000
      IF(JCOUPC.GE.0)THEN
        KVM3=MBFMX3*16
        IF(JREACT.NE.0)KVM3=MBFMX3*21
      ELSE
        KVM3=(1+MBFMX3*16)/2
        IF(JREACT.NE.0)KVM3=(1+MBFMX3*21)/2
      END IF
      CALL MEMO(1,LVM3,KVM3,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N3=0
            DO N=1,L-1
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**WRITE CORIOLIS GRIDS TO DISC (93)
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
                CALL DUMCR3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,IN2,
     1          NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),
     2          W(LXX),W(LX0),W(LXL),W(LXM),K,L,N,W(LC3),W(LC3),
     3          W(LVM3),(LVM3),1,1,W(LMODNT),W(LNBF),W(LMBF),W(LVP1),
     4          W(LVP2),W(LV1),W(LV1),W(LV2),W(LV2),1)
              ELSE
                CALL DUMVR3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,
     1          IN2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     2          W(LB+NMODE*NMODE),W(LAA),W(LBB),W(LBB+NMODE),W(LXX),
     3          W(LXX+NATOM*3*722),W(LX0),W(LXL),W(LXM),K,L,N,W(LC3),
     4          W(LC3),W(LVM3),W(LVM3),1,W(LXQ+K3TAU),W(LVP0),
     5          W(LVP1),W(LVP2),W(LVP3),KVP,1,W(LMODNT),W(LNBF),
     6          W(LMBF),W(LC1),W(LC1),W(LC2),W(LC2),IK2TAU)
              END IF
              N3=N3+IN2
            END DO
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPC.EQ.3)GO TO 3000
      IF(JCOUPC.GE.0)THEN
        KVM4=MBFMX4*19
        IF(JREACT.NE.0)KVM4=MBFMX4*24
      ELSE
        KVM4=(1+MBFMX4*19)/2
        IF(JREACT.NE.0)KVM4=(1+MBFMX4*24)/2
      END IF
      CALL MEMO(1,LVM4,KVM4,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N3=0
            DO N=1,L-1
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
              M3=0
              DO M=1,N-1
                CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**WRITE CORIOLIS GRIDS TO DISC (94)
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2          THEN
                  CALL DUMCR4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),W(LXQ+M3),
     1            IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),
     2            W(LB),W(LAA),W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,L,
     3            N,M,W(LC4),W(LC4),W(LVM4),W(LVM4),1,1,W(LMODNT),
     4            W(LNBF),W(LMBF),W(LVP1),W(LVP2),W(LVP3),W(LV1),
     5            W(LV1),W(LV2),W(LV2),W(LV3),W(LV3),1)
                ELSE
                  CALL DUMVR4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1            W(LXQ+M3),IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),
     2            W(LXZ),W(LAB),W(LB),W(LB+NMODE*NMODE),W(LAA),
     3            W(LBB),W(LBB+NMODE),W(LXX),W(LXX+NATOM*3*722),
     4            W(LX0),W(LXL),W(LXM),K,L,N,M,W(LC4),W(LC4),W(LVM4),
     5            W(LVM4),1,W(LXQ+K3TAU),W(LMODNT))
                END IF
                M3=M3+IM2
              END DO
              N3=N3+IN2
            END DO
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPC.EQ.4)GO TO 3000
C**5-MODE AND HIGHER
3000  CONTINUE

C**SAVE IMPORTANT QUANTITIES
      NNMODS=NNMODE
      NAMODS=NAMODE
      NVMODS=NVMODE
C     ICOUPS=ICOUPL ??
C     ICOUCS=ICOUPC ??
C     JCOUPS=JCOUPL ??
C     JCOUCS=JCOUPC ??
      JREACS=JREACT
      NREACS=NREACT
C**LOOP ROUND CONTRACTION SCHEMES (NON-LINEAR ONLY AT PRESENT)
      DO 5000 LCONT=1,LCOUNT
      WRITE(IOUT,415)LCONT
C*********************************************
C*********************************************
C**    PERFORM VMI FOR INDIVIDUAL CONTRACTIONS
C*********************************************
C*********************************************
      JCOUPL=NCOUPL(LCONT)
      JCOUPC=NCOUPC(LCONT)
      ICOUPL=IABS(JCOUPL)
      ICOUPC=IABS(JCOUPC)

C**RESET IMPORTANT QUANTITIES FOR EACH CONTRACTION IF REQUIRED
      IF(LCOUNT.GT.1)THEN
        NUMBER=ICONT(LCONT)
        JREACT=0
        NREACT=0
        DO I=1,NUMBER
          J=JCONT(LCONT,I)
          IF(J.EQ.IREACT)JREACT=IREACT
        END DO
        IF(JREACT.NE.0)NREACT=1
        NAMODE=NUMBER-NREACT
C**NUMBER COULD BE SINGLE (RPH) MODE.......DO NOTHING
        IF(NAMODE.EQ.0)GO TO 5005
        NVMODE=NUMBER
C**NNMODE LIKE NSMODE.......TOTAL MODES IN VCI
        NNMODE=NUMBER
        ICOUPL=MIN0(ICOUPS,NAMODE)
        ICOUPC=MIN0(ICOUCS,NAMODE)
        DO ICPL=1,ICOUPL
          FACTOR(ICPL)=1.D0
          DO I=2,NAMODE
            FACTOR(ICPL)=FACTOR(ICPL)*I
          END DO
          DO I=1,ICPL
            FACTOR(ICPL)=FACTOR(ICPL)/I
          END DO
          IDENOM=NAMODE-ICPL
          DO I=1,IDENOM
            FACTOR(ICPL)=FACTOR(ICPL)/I
          END DO
          IF(JREACT.GT.0)FACTOR(ICPL)=FACTOR(ICPL)+1
        END DO
C**********************************************************************
C**********************************************************************
C**                                    DUMP CORIOLIS DATA TO DISC (VCI)
C**                                 FOR EACH CONTRACTION SCHEME IN TURN
C**********************************************************************
C**********************************************************************
        IF(JREACS.GT.0.AND.JREACT.LE.0)THEN
          I91=191
          I92=192
          I93=193
          I94=194
C????????????????????????????????????
CCCC      CALL DISCDP(W,KCONT,LCONT)
C????????????????????????????????????
        ELSE
          I91=91
          I92=92
          I93=93
          I94=94
        END IF
C**********************************************************************
C**********************************************************************
      END IF

C     REWIND 31
C     REWIND 32
C     REWIND 33
C     REWIND 34
C     REWIND 35
C     REWIND 36
C     REWIND 37
C     REWIND 38
C     REWIND 39

      REWIND 58
      REWIND 59

      ITIM=-1
      DO 9999 IABC=1,9
C
C**FIND RELEVANT SYMMETRIES FOR THIS IABC
      KSTART=1
      KEND=2*JTHIS+1
      IF(.NOT.TRIAT)CALL GETSYM(IABC,KSTART,KEND)
C
C**K=0/K=0
      REWIND 33
C**K=1/K=0
      REWIND 34
C**K=0/K=1
      REWIND 35
C**K=1/K=1
      REWIND 36

C     IF(ICOUPL.GT.1)CALL MEMO(1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
C**********************************************************************
C**   LOOP ROUND Ka = 0,1 IF RPH FOR INTEGER OR HALF-INTEGER TORSION
C**********************************************************************
C     ITIM=-1
      LMAXR=1
      LMAXL=1
      IF(JREACT.GT.0)THEN
        LMAXR=IFLAUD
        LMAXL=IFLAUD
C       LMAXR=2
C       LMAXL=2
      END IF
C**TEMPORARY-LIN
      IF(LINEAR)THEN
        LMAXR=JMAX+1
        LMAXL=JMAX+1
      END IF
      DO 3330 KROTR=1,LMAXR
      REWIND 70
C**READ INTO ELEMENT (2) ... RHS
      K2=0
      DO K=1,NAMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(K.GT.NONLIN)THEN
          DO IDUM=1,KROTR
            CALL IN70(W(LH+K2+IK1),IK3,IK2)
          END DO
        END IF
        K2=K2+IK1
        IF(K.GT.NONLIN)K2=K2+IK1
      END DO
      REWIND 70
      DO 3300 KROTL=1,LMAXL
C**READ INTO ELEMENT (1) ... LHS
      K2=0
      DO K=1,NAMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(K.GT.NONLIN)THEN
          CALL IN70(W(LH+K2),IK3,IK2)
        END IF
        K2=K2+IK1
        IF(K.GT.NONLIN)K2=K2+IK1
      END DO
C**TEMPORARY-LIN
      IF(ICOUPC.GE.0)REWIND 21
      IF(ICOUPC.GT.1)REWIND 22
      IF(ICOUPC.GT.2)REWIND 23
      IF(ICOUPC.GT.3)REWIND 24
      IND=MOD(KROTR-1,2)*LMAXR+KROTL+32
C**********************************************************************
C**   GET BASIC INTEGRALS FOR FILES  (31)-(39) IF J>0
C**********************************************************************
C     DO 4444 IABC=1,9
      ITIM=ITIM+1
      IF(ICOUPC.GE.0)REWIND 91
      IF(ICOUPC.GT.1)REWIND 92
      IF(ICOUPC.GT.2)REWIND 93
      IF(ICOUPC.GT.3)REWIND 94
C**BOTH LINEAR AND RPH
      IF(IREACT.GT.0)THEN
        CALL INTARR(W(LNBF),W(LMBF),NSMODE,IK1TAU,IK2TAU,IK3TAU)
C**GET TORSION-ONLY INTEGRALS IF RPH
        IF(JREACT.GT.0)THEN
C*****************************   LIKE V0MI1 (NMODE)
C**
CC        IF(IABC.LT.7)THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
            CALL MEMO(1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
            CALL GETBP1(W(LIP1),ISIZE1,NMODE,NNMODE,W(LMXBAS),NSIZE1,
     1      0)
            KXA1=NSIZE1*(NSIZE1+1)/2
            CALL MEMO(1,LXA0,KXA1,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE1,NSIZE1,W(LXA0),NSIZE1,W(LXA0),ICI)
            CALL V0MV0(NAMODE,1,W(LH+K2TAU),W(LXQ+K3TAU),W(LXA0),
     1      NSIZE1,IK3TAU,IK2TAU,W(LIP1),ISIZE1,W(LVM0),W(LVM0),J21,
     2      KROTL,KROTR,IABC,W(LMODNT))
            CALL MATOUT(W(LXA0),W(LXA0),NSIZE1,21)
            CALL MEMO(-1,LXA0,KXA1,0,0,0,0,0,0,0,0)
            CALL MEMO(-1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
CC        END IF
C**
C*****************************   LIKE V0MI1 (NMODE)
        END IF
      END IF
     
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K=0
      K2=0
      K3=0
      KRK=1
      DO KK=1,NAMODE
        IF(KK.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        IF(LCOUNT.GT.1)THEN
          IF(KRK.LE.NAMODE)THEN
            KNEXT=JCONT(LCONT,KRK)
            IF(KK.EQ.KNEXT)KRK=KRK+1
          ELSE
            KNEXT=0
          END IF
        ELSE
          KNEXT=KK
        END IF
        CALL INTARR(W(LNBF),W(LMBF),KK,IK1,IK2,IK3)
C**NEXT K
        K=KNEXT
C**TEMPORARY-LIN
        MK2=0
        IF(K.GT.NONLIN)MK2=IK1
C**TEMPORARY-LIN
        K1=K-1
C**DO WE WANT THIS KK?
        IF(KK.NE.KNEXT)THEN
C**DUMMY READ IF MISSING MODE KK
          IF((JREACT.LE.0.AND.KK.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     1    THEN
            IF(ICOUPC.GT.0.AND.IREACT.EQ.0)
     1      CALL DUMRM1(IK2,W(LVM1),W(LVM1))
          ELSE
CCCC        IF(ICOUPC.GT.0)CALL DUMRM1(IK2,IK2TAU,W(LVM1),W(LVM1))
          END IF
          GO TO 771
        END IF
        IF(ICOUPC.EQ.0)GO TO 7701
C**SKIP IF UNWANTED KK
        IF(KNEXT.EQ.0.OR.KK.NE.KNEXT)GO TO 771
C**CORIOLIS
        IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
          IF(IABC.LT.7)THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
            CALL MEMO(1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
            CALL GETBP1(W(LIP1),ISIZE1,NMODE,K,W(LMXBAS),NSIZE1,0)
            KXA1=NSIZE1*(NSIZE1+1)/2
            CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE1,NSIZE1,W(LXA1),NSIZE1,W(LXA1),ICI)
C**TEMPORARY-LIN
            CALL V0MI1(NAMODE,1,K,W(LH+K2),W(LH+K2+MK2),W(LXQ+K3),
     1      W(LXA1),NSIZE1,IK3,IK2,W(LIP1),ISIZE1,W(LVM1),W(LVM1),J21,
     2      KROTL,KROTR,IABC,W(LMODNT))
C**TEMPORARY-LIN
            CALL MATOUT(W(LXA1),W(LXA1),NSIZE1,21)
            CALL MEMO(-2,LXA1,KXA1,LIP1,KPPP1,0,0,0,0,0,0)
          END IF
        ELSE
C*****************************   LIKE V0MI2 (K+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
          IF(K.GT.NONLIN)THEN
            KPPP2=KLJP2
            IND2=1
          ELSE
            KPPP2=KIP2
            IND2=0
          END IF
          CALL MEMO(1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
          IF(K.GT.NONLIN)THEN
            CALL GETBP1(W(LIP2),LSIZE2,NMODE,K,W(LMXBAS),NSIZE2,IND2)
            MSIZE2=NSIZE2
          ELSE
            CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,NSMODE,W(LMXBAS),
     1      NSIZE2,IND2)
            MSIZE2=ISIZE2
          END IF
          KXA2=NSIZE2*(NSIZE2+1)/2
          JCI=MAXBFN(W(LMXBAS),NMODE,K,1)
          IF(JREACT.GT.0)THEN
            JCIM=MAX0(NTAU(1)+1,NTAU(2)+1)
          ELSE
            JCIM1=2*(NBAS(1,1)+1)-1
            JCIM2=2*(NBAS(1,2)+1)-1
            JCIM=MAX0(JCIM1,JCIM2)
          END IF
          KTEMP=JCI*JCI*5
          CALL MEMO(2,LXA1,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
          KXK0=5*JCI*JCI*IK2
          KXL0=5*JCIM*JCIM*IK2TAU
          CALL MEMO(2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
C**ZEROISE MATRIX
          CALL DIAGZ(NSIZE2,NSIZE2,W(LXA1),NSIZE2,W(LXA1),ICI)
          CALL V0MV1(NAMODE,1,2,K,W(LH+K2),W(LXQ+K3),W(LH+K2TAU),
     1    W(LXQ+K3TAU),IK3,IK2,IK3TAU,IK2TAU,W(LXA1),NSIZE2,W(LIP2),
     2    MSIZE2,W(LTEMP),JCI,JCIM,W(LXK0),W(LXL0),W(LVM1),
     3    W(LVM1),J21,KROTL,KROTR,IABC,W(LMODNT))
          CALL MATOUT(W(LXA1),W(LXA1),NSIZE2,21)
          CALL MEMO(-2,LXA1,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
          CALL MEMO(-2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0MI2 (K+NMODE)
        END IF
771   CONTINUE
        IF(ICOUPC.EQ.1)GO TO 7701
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L=0
        L2=0
        L3=0
        LRL=1
        DO LL=1,KK-1
          IF(LL.EQ.1.AND.KK.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          IF(LCOUNT.GT.1)THEN
            IF(LRL.LE.NAMODE)THEN
              LNEXT=JCONT(LCONT,LRL)
              IF(LL.EQ.LNEXT)LRL=LRL+1
            ELSE
              LNEXT=0
            END IF
          ELSE
            LNEXT=LL
          END IF
          CALL INTARR(W(LNBF),W(LMBF),LL,IL1,IL2,IL3)
C**NEXT L
          L=LNEXT
C**TEMPORARY-LIN
          ML2=0
          IF(L.GT.NONLIN)ML2=IL1
C**TEMPORARY-LIN
C**DO WE WANT THIS KK AND THIS LL?
          IF(KK.NE.KNEXT.OR.LL.NE.LNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1      (NSMODE.EQ.NAMODE))THEN
              IF(IREACT.EQ.0)CALL DUMRM2(IK2,IL2,W(LVM2),W(LVM2))
            ELSE
CCCC          CALL DUMRM2(IK2,IL2,IK2TAU,W(LVM2),W(LVM2))
            END IF
            GO TO 772
          END IF
C**CORIOLIS
          IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1    (NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
            CALL MEMO(1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
            CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,L,W(LMXBAS),NSIZE2,0)
            KXA2=NSIZE2*(NSIZE2+1)/2
            JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
            JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
            KTEMP=JCI2*JCI2*5
            CALL MEMO(2,LXA2,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
C???????????????????????????
C           KTEMP1=JCI2*JCI2*4
C           CALL MEMO(1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
C???????????????????????????
            KXK0=5*JCI1*JCI1*IK2
            KXL0=5*JCI2*JCI2*IL2
            CALL MEMO(2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE2,NSIZE2,W(LXA2),NSIZE2,W(LXA2),ICI)
C**TEMPORARY-LIN
            CALL V0MI2(NAMODE,1,2,K,L,W(LH+K2),W(LH+K2+MK2),W(LXQ+K3),
     1      W(LH+L2),W(LH+L2+ML2),
     1      W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXA2),NSIZE2,W(LIP2),ISIZE2,
     2      W(LTEMP),W(LTEMP1),JCI1,JCI2,W(LXK0),W(LXL0),W(LVM2),
     3      W(LVM2),J21,KROTL,KROTR,IABC,W(LMODNT))
C**TEMPORARY-LIN
            CALL MATOUT(W(LXA2),W(LXA2),NSIZE2,22)
            CALL MEMO(-3,LXA2,KXA2,LIP2,KPPP2,LTEMP,KTEMP,0,0,0,0)
            CALL MEMO(-2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
C???????????????????????????
C           CALL MEMO(-1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
C???????????????????????????
          ELSE
C*****************************   LIKE V0MI3 (K+L+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
            IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
              KPPP3=KLJP3
              IND3=1
            ELSE
              KPPP3=KIP3
              IND3=0
            END IF
            CALL MEMO(1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
            IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
              CALL GETBP2(W(LIP3),LSIZE3,NMODE,K,L,W(LMXBAS),NSIZE3,
     1        IND3)
              MSIZE3=NSIZE3
            ELSE
              CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,NSMODE,W(LMXBAS),
     1        NSIZE3,IND3)
              MSIZE3=ISIZE3
            END IF
            KXA3=NSIZE3*(NSIZE3+1)/2
            JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
            JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
            IF(JREACT.GT.0)THEN
              JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1)
            ELSE
              JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1)-1
              JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1)-1
              JCIM=MAX0(JCIM1,JCIM2)
            END IF
            KTEMP=JCI1*JCI2*JCI1*JCI2*7
            CALL MEMO(2,LXA2,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
            KTEMP1=JCI2*JCI2*7
            CALL MEMO(1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
            KXK0=7*JCI1*JCI1*IK2
            KXL0=7*JCI2*JCI2*IL2
            KXN0=7*JCIM*JCIM*IK2TAU
            CALL MEMO(3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE3,NSIZE3,W(LXA2),NSIZE3,W(LXA2),ICI)
            CALL V0MV2(NAMODE,1,2,3,K,L,W(LH+K2),W(LXQ+K3),W(LH+L2),
     1      W(LXQ+L3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IK3TAU,
     2      IK2TAU,W(LXA2),NSIZE3,W(LIP3),MSIZE3,W(LTEMP),W(LTEMP1),
     3      JCI1,JCI2,JCIM,W(LXK0),W(LXL0),W(LXN0),W(LVM2),
     4      W(LVM2),J21,KROTL,KROTR,IABC,W(LMODNT))
            CALL MATOUT(W(LXA2),W(LXA2),NSIZE3,22)
            CALL MEMO(-2,LXA2,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
            CALL MEMO(-1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
            CALL MEMO(-3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
            CALL MEMO(-1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0MI3 (K+L+NMODE)
          END IF
772    CONTINUE
          IF(ICOUPC.EQ.2)GO TO 7702
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N=0
          N2=0
          N3=0
          NRN=1
          DO NN=1,LL-1
            IF(NN.EQ.1.AND.LL.EQ.2.AND.KK.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            IF(LCOUNT.GT.1)THEN
              IF(NRN.LE.NAMODE)THEN
                NNEXT=JCONT(LCONT,NRN)
                IF(NN.EQ.NNEXT)NRN=NRN+1
              ELSE
                NNEXT=0
              END IF
            ELSE
              NNEXT=NN
            END IF
            CALL INTARR(W(LNBF),W(LMBF),NN,IN1,IN2,IN3)
C**NEXT N
            N=NNEXT
C**TEMPORARY-LIN
            MN2=0
            IF(N.GT.NONLIN)MN2=IN1
C**TEMPORARY-LIN
C**DO WE WANT THIS KK AND THIS LL AND THIS NN?
            IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
                IF(IREACT.EQ.0)CALL DUMRM3(IK2,IL2,IN2,W(LVM3),W(LVM3))
              ELSE
CCCC            CALL DUMRM3(IK2,IL2,IN2,IK2TAU,W(LVM3),W(LVM3))
              END IF
              GO TO 773
            END IF
C**CORIOLIS
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1      N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
              CALL MEMO(1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
              CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,N,W(LMXBAS),NSIZE3,
     1        0)
              KXA3=NSIZE3*(NSIZE3+1)/2
              JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
              JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
              JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
              KTEMP=JCI2*JCI3*JCI2*JCI3*7
              CALL MEMO(2,LXA3,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
              KTEMP1=JCI3*JCI3*7
              CALL MEMO(1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
C???????????????????????????????????????
C             KTEMP2=JCI2*JCI3*JCI2*JCI3*3
C             CALL MEMO(1,LTEMP2,KTEMP2,0,0,0,0,0,0,0,0)
C???????????????????????????????????????
              KXK0=7*JCI1*JCI1*IK2
              KXL0=7*JCI2*JCI2*IL2
              KXN0=7*JCI3*JCI3*IN2
              CALL MEMO(3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(NSIZE3,NSIZE3,W(LXA3),NSIZE3,W(LXA3),ICI)
C**TEMPORARY-LIN
              CALL V0MI3(NAMODE,1,2,3,K,L,N,W(LH+K2),W(LH+K2+MK2),
     1        W(LXQ+K3),W(LH+L2),W(LH+L2+ML2),W(LXQ+L3),W(LH+N2),
     2        W(LH+N2+MN2),W(LXQ+N3),IK3,IK2,IL3,IL2,
     2        IN3,IN2,W(LXA3),NSIZE3,W(LIP3),ISIZE3,W(LTEMP),W(LTEMP1),
     3        W(LTEMP2),JCI1,JCI2,JCI3,W(LXK0),W(LXL0),W(LXN0),W(LVM3),
     4        W(LVM3),J21,KROTL,KROTR,IABC,W(LMODNT))
C**TEMPORARY-LIN
              CALL MATOUT(W(LXA3),W(LXA3),NSIZE3,23)
              CALL MEMO(-3,LXA3,KXA3,LIP3,KPPP3,LTEMP,KTEMP,0,0,0,0)
C???????????????????????????????????????
C             CALL MEMO(-1,LTEMP2,KTEMP2,0,0,0,0,0,0,0,0)
C???????????????????????????????????????
              CALL MEMO(-1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
              CALL MEMO(-3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
            ELSE
C*****************************   LIKE V0MI4 (K+L+N+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
              IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                KPPP4=KLJP4
                IND4=1
              ELSE
                KPPP4=KIP4
                IND4=0
              END IF
              CALL MEMO(1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
              IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                CALL GETBP3(W(LIP4),LSIZE4,NMODE,K,L,N,W(LMXBAS),
     1          NSIZE4,IND4)
                MSIZE4=NSIZE4
              ELSE
                CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,NSMODE,
     1          W(LMXBAS),NSIZE4,IND4)
                MSIZE4=ISIZE4
              END IF

!              KXA4=NSIZE4*(NSIZE4+1)/2
              if (mod(nsize4,2) .eq. 0) then
                 kxa4 = nsize4 / 2
                 kxa4 = kxa4 * (nsize4 + 1)
              else
                 kxa4 = (nsize4 + 1) / 2
                 kxa4 = kxa4 * nsize4
              end if

              JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
              JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
              JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
              IF(JREACT.GT.0)THEN
                JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1,NTAU(4)+1)
              ELSE
                JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1,NBAS(3,1)+1)-1
                JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1,NBAS(3,2)+1)-1
                JCI=MAX0(JCIM1,JCIM2)
              END IF
              KTEMP=JCI1*JCI2*JCI3*JCI1*JCI2*JCI3*9
              CALL MEMO(2,LXA3,KXA4,LTEMP,KTEMP,0,0,0,0,0,0)
              KTEMP1=JCI3*JCI3*9
              KTEMP2=JCI2*JCI3*JCI2*JCI3*9
              CALL MEMO(2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
              KXK0=9*JCI1*JCI1*IK2
              KXL0=9*JCI2*JCI2*IL2
              KXN0=9*JCI3*JCI3*IN2
              KXM0=9*JCIM*JCIM*IK2TAU
              CALL MEMO(4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(NSIZE4,NSIZE4,W(LXA3),NSIZE4,W(LXA3),ICI)
              CALL V0MV3(NAMODE,1,2,3,4,K,L,N,W(LH+K2),W(LXQ+K3),
     1        W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+K2TAU),
     2        W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,IK3TAU,IK2TAU,
     3        W(LXA3),NSIZE4,W(LIP4),MSIZE4,W(LTEMP),W(LTEMP1),
     4        W(LTEMP2),JCI1,JCI2,JCI3,JCIM,W(LXK0),W(LXL0),W(LXN0),
     5        W(LXM0),W(LVM3),W(LVM3),J21,KROTL,KROTR,IABC,W(LMODNT))
              CALL MATOUT(W(LXA3),W(LXA3),NSIZE4,23)
              CALL MEMO(-2,LXA3,KXA4,LTEMP,KTEMP,0,0,0,0,0,0)
              CALL MEMO(-2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
              CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,0,0)
              CALL MEMO(-1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0MI4 (K+L+N+NMODE)
            END IF
773   CONTINUE
            IF(ICOUPC.EQ.3)GO TO 7703
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M=0
            M2=0
            M3=0
            MRM=1
            DO MM=1,NN-1
              IF(MM.EQ.1.AND.NN.EQ.2.AND.LL.EQ.3.AND.KK.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              IF(LCOUNT.GT.1)THEN
                IF(MRM.LE.NAMODE)THEN
                  MNEXT=JCONT(LCONT,MM)
                  IF(MM.EQ.MNEXT)MRM=MRM+1
                ELSE
                  MNEXT=0
                END IF
              ELSE
                MNEXT=MM
              END IF
              CALL INTARR(W(LNBF),W(LMBF),MM,IM1,IM2,IM3)
C**NEXT M
              M=MNEXT
C**TEMPORARY-LIN
              MM2=0
              IF(M.GT.NONLIN)MM2=IM1
C**TEMPORARY-LIN
C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM?
              IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
     1        MM.NE.MNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2          THEN
                  IF(IREACT.EQ.0)CALL DUMRM4(IK2,IL2,IN2,IM2,W(LVM4),
     1            W(LVM4))
                ELSE
CCCC              CALL DUMRM4(IK2,IL2,IN2,IM2,IK2TAU,W(LVM4),W(LVM4))
                END IF
                GO TO 774
              END IF
C**CORIOLIS
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2        THEN
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
                CALL MEMO(1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
                CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,M,W(LMXBAS),
     2          NSIZE4,0)

!                KXA4=NSIZE4*(NSIZE4+1)/2
                if (mod(nsize4,2) .eq. 0) then
                   kxa4 = nsize4 / 2
                   kxa4 = kxa4 * (nsize4 + 1)
                else
                   kxa4 = (nsize4 + 1) / 2
                   kxa4 = kxa4 * nsize4
                end if
 
                JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
                JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
                JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
                JCI4=MAXBFN(W(LMXBAS),NMODE,M,1)
                KTEMP=JCI2*JCI3*JCI4*JCI2*JCI3*JCI4*9
                CALL MEMO(2,LXA4,KXA4,LTEMP,KTEMP,0,0,0,0,0,0)
                KTEMP1=JCI4*JCI4*9
                KTEMP2=JCI3*JCI4*JCI3*JCI4*9
                CALL MEMO(2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
C???????????????????????????????????????????????????
C               KTEMP3=JCI2*JCI3*JCI4*JCI2*JCI3*JCI4*4
C               CALL MEMO(3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,KTEMP3,
C    1          0,0,0,0)
C???????????????????????????????????????????????????
                KXK0=9*JCI1*JCI1*IK2
                KXL0=9*JCI2*JCI2*IL2
                KXN0=9*JCI3*JCI3*IN2
                KXM0=9*JCI4*JCI4*IM2
                CALL MEMO(4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,0,
     1          0)
C**ZEROISE MATRIX
                CALL DIAGZ(NSIZE4,NSIZE4,W(LXA4),NSIZE4,W(LXA4),ICI)
C**TEMPORARY-LIN
                CALL V0MI4(NAMODE,1,2,3,4,K,L,N,M,W(LH+K2),
     1          W(LH+K2+MK2),W(LXQ+K3),W(LH+L2),W(LH+L2+ML2),W(LXQ+L3),
     2          W(LH+N2),W(LH+N2+MN2),W(LXQ+N3),W(LH+M2),W(LH+M2+MM2),
     2          W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA4),
     3          NSIZE4,W(LIP4),ISIZE4,W(LTEMP),W(LTEMP1),W(LTEMP2),
     4          W(LTEMP3),JCI1,JCI2,JCI3,JCI4,W(LXK0),W(LXL0),W(LXN0),
     5          W(LXM0),W(LVM4),W(LVM4),J21,IABC,W(LMODNT))
C**TEMPORARY-LIN
                CALL MATOUT(W(LXA4),W(LXA4),NSIZE4,24)
                CALL MEMO(-3,LXA4,KXA4,LIP4,KPPP4,LTEMP,KTEMP,0,0,0,0)
                CALL MEMO(-2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
C???????????????????????????????????????????????????
C               CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,KTEMP3,
C    1          0,0,0,0)
C???????????????????????????????????????????????????
                CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          0,0)
              ELSE
C*****************************   LIKE V0MI5 (K+L+N+M+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR FIVE MODES
                IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1          M.GT.NONLIN)THEN
                  KPPP5=KLJP5
                  IND5=1
                ELSE
                  KPPP5=KIP5
                  IND5=0
                END IF
                CALL MEMO(1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
                IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1          M.GT.NONLIN)THEN
                  CALL GETBP4(W(LIP5),LSIZE5,NMODE,K,L,N,M,W(LMXBAS),
     1            NSIZE5,IND5)
                  MSIZE5=NSIZE5
                ELSE
                  CALL GETBP5(W(LIP5),ISIZE5,NMODE,K,L,N,M,NSMODE,
     1            W(LMXBAS),NSIZE5,IND5)
                  MSIZE5=ISIZE5
                END IF
                KXA5=NSIZE5*(NSIZE5+1)/2
                JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
                JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
                JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
                JCI4=MAXBFN(W(LMXBAS),NMODE,M,1)
                IF(JREACT.GT.0)THEN
                  JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1,NTAU(4)+1,
     1            NTAU(5)+1)
                ELSE
                  JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1,NBAS(3,1)+1,
     1            NBAS(4,1)+1)-1
                  JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1,NBAS(3,2)+1,
     1            NBAS(4,2)+1)-1
                  JCI=MAX0(JCIM1,JCIM2)
                END IF
CCCC            KTEMP=JCI1*JCI2*JCI3*JCI4*JCI1*JCI2*JCI3*JCI4*2
                CALL MEMO(1,LXA4,KXA5,0,0,0,0,0,0,0,0)
                KTEMP1=JCI4*JCI4*11
                KTEMP2=JCI3*JCI4*JCI3*JCI4*11
                KTEMP3=JCI2*JCI3*JCI4*JCI2*JCI3*JCI4*11
                CALL MEMO(3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,KTEMP3,
     1          0,0,0,0)
                KXK0=11*JCI1*JCI1*IK2
                KXL0=11*JCI2*JCI2*IL2
                KXN0=11*JCI3*JCI3*IN2
                KXM0=11*JCI4*JCI4*IM2
                KXP0=11*JCIM*JCIM*IK2TAU
                CALL MEMO(5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          LXP0,KXP0)
C**ZEROISE MATRIX
                CALL DIAGZ(NSIZE5,NSIZE5,W(LXA4),NSIZE5,W(LXA4),ICI)
                CALL V0MV4(NAMODE,1,2,3,4,5,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1          W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2          W(LXQ+M3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,
     3          IN2,IM3,IM2,IK3TAU,IK2TAU,W(LXA4),NSIZE5,W(LIP5),
     4          MSIZE5,W(LTEMP1),W(LTEMP2),W(LTEMP3),JCI1,JCI2,JCI3,
     5          JCI4,JCIM,W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXP0),
     6          W(LVM4),W(LVM4),J21,KROTL,KROTR,IABC,W(LMODNT))
                CALL MATOUT(W(LXA4),W(LXA4),NSIZE5,24)
                CALL MEMO(-1,LXA4,KXA5,0,0,0,0,0,0,0,0)
                CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,KTEMP3,
     1          0,0,0,0)
                CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          LXP0,KXP0)
                CALL MEMO(-1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0MI5 (K+L+N+M+NMODE)
              END IF
774   CONTINUE
              IF(ICOUPC.EQ.4)GO TO 7704
C**5-MODE COUPLING HERE IF NEEDED
7704  CONTINUE
              M2=M2+IM1
              M3=M3+IM2
            END DO
7703  CONTINUE
            N2=N2+IN1
            N3=N3+IN2
          END DO
7702  CONTINUE
          L2=L2+IL1
          L3=L3+IL2
        END DO
7701  CONTINUE
        K2=K2+IK1
        K3=K3+IK2
      END DO
C4444  CONTINUE
C     IF(ICOUPL.GT.1)CALL MEMO(-1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
C???????????????????????????????????????????????????????????
CCCC  IF(ICOUPL.GT.0)CALL MEMO(-1,LVM1,KVM1,0,0,0,0,0,0,0,0)
CCCC  IF(ICOUPL.GT.1)CALL MEMO(-1,LVM2,KVM2,0,0,0,0,0,0,0,0)
CCCC  IF(ICOUPL.GT.2)CALL MEMO(-1,LVM3,KVM3,0,0,0,0,0,0,0,0)
CCCC  IF(ICOUPL.GT.3)CALL MEMO(-1,LVM4,KVM4,0,0,0,0,0,0,0,0)
C???????????????????????????????????????????????????????????
      IF(KROTL.EQ.1.AND.KROTR.EQ.1)ITIM=-1
C**********************************************************************
C**                                          LOOP OVER SYMMETRIES (LHS)
C**********************************************************************
      DO 5500 NSL=1,NVSYM
      ISIZEL=NTOT(NSL)
      IF(ISIZEL.EQ.0)GO TO 5501
      CALL PUTJP(W(LJP),JSIZE,W(LIPL),ISIZMX,NMODE,NSL)
C**********************************************************************
C**                                          LOOP OVER SYMMETRIES (RHS)
C**********************************************************************
      DO 5502 NSR=NSL,NVSYM
      ISIZER=NTOT(NSR)
      IF(ISIZER.EQ.0)GO TO 5503
      CALL PUTJP(W(LJP),JSIZE,W(LIPR),ISIZMX,NMODE,NSR)
C**
C**SEARCH FOR MATCH FOR SYMMETRIES NSL,NSR IN RO-VIB BLOCKS
      IGOT=0
      DO NROT=1,NRSYM
        NL=NSYMSL(NROT)
C**NO LHS SYMMETRIES
        IF(NL.EQ.0)GO TO 5503
        DO NLS=1,NL
          NSSL=MSYMSL(NROT,NLS)
          NR=NSYMSR(NROT,NLS)
C**LHS SYMMETRY WRONG OR NO CORRESPONDING RHS SYMMETRIES
          IF(NSSL.NE.NSL.OR.NR.EQ.0)GO TO 5504
          DO NRS=1,NR
            NSSR=MSYMSR(NROT,NLS,NR)
C**FOUND A MATCH FOR THIS NSL,NSR
            IF(NSR.EQ.NSSR)IGOT=1
          END DO
C**FOUND MATCH - PROCESS IT
          IF(IGOT.NE.0)GO TO 5505
5504      CONTINUE
        END DO
      END DO
C**IF NO MATCH, DON'T NEED TO DO THIS NSL,NSR
      IF(IGOT.EQ.0)GO TO 5503
5505  CONTINUE
C**
      IF(ICOUPC.GE.0)THEN
        REWIND 21
      END IF
      IF(ICOUPC.GT.1)THEN
        REWIND 22
      END IF
      IF(ICOUPC.GT.2)THEN
        REWIND 23
      END IF
      IF(ICOUPC.GT.3)THEN
        REWIND 24
      END IF
C**********************************************************************
C**   SET UP FILES  (31)-(39) FOR COMMON INTEGRALS
C**********************************************************************
C     DO 9999 IABC=1,9
      ITIM=ITIM+1
      KXA=ISIZEL*ISIZER
CC    KXA=ISIZEL
      CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
      CALL DIAGZ(ISIZEL,ISIZER,W(LXA),IDUM,W(LXA),0)
CC    CALL DIAGZ(ISIZEL,1,W(LXA),IDUM,W(LXA),0)
C**WRITE ZERO MATRIX TO (32)
CC    REWIND 32
CC    DO I=1,ISIZER
CC      CALL OUTXA(W(LXA),ISIZEL)
CC    END DO
      IF(IABC.LT.7)THEN
        IHERM=1
      ELSE
        IHERM=-1
      END IF
C**BOTH LINEAR AND RPH
      IF(IREACT.GT.0)THEN
        CALL INTARR(W(LNBF),W(LMBF),NNMODE,IK1TAU,IK2TAU,IK3TAU)
C**GET TORSION-ONLY INTEGRALS IF RPH
        IF(JREACT.GT.0)THEN
C****************************************  LIKE VMI1 (NMODE)
C**
C         IF(IABC.LT.7)THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
            CALL MEMO(1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
            CALL GETBP1(W(LIP1),ISIZE1,NMODE,NNMODE,W(LMXBAS),NSIZE1,
     1      0)
C**READ INTO MATRIX
            KXA1=NSIZE1*(NSIZE1+1)/2
            CALL MEMO(1,LXA0,KXA1,0,0,0,0,0,0,0,0)
            CALL MATIN(W(LXA0),W(LXA0),NSIZE1,21)
            CALL VMV0(NNMODE,W(LXA),W(LXA0),NSIZE1,W(LIPL),W(LIPR),
     1      ISIZMX,ISIZEL,ISIZER,W(LIP1),ISIZE1)
            CALL MEMO(-1,LXA0,KXA1,0,0,0,0,0,0,0,0)
            CALL MEMO(-1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
C         END IF
C**
C****************************************  LIKE VMI1 (NMODE)
        END IF
      END IF
      IF(JREACT.LE.0.OR.NSMODE.EQ.NAMODE)THEN
        IF(IABC.LT.7.AND.NSL.EQ.NSR)
     1  CALL VMI0(W(LXA),ISIZEL,ISIZER,IABC)
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K=0
      K2=0
      K3=0
      DO KK=1,NAMODE
        IF(KK.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        IF(LCOUNT.GT.1)THEN
          KNEXT=JCONT(LCONT,KK)
        ELSE
          KNEXT=KK
        END IF
        IF(KNEXT-K.GT.1)THEN
C**UPDATE K2,K3 FOR MISSING MODES
          DO KADJ=K+1,KNEXT-1
            CALL INTARR(W(LNBF),W(LMBF),KADJ,IK1,IK2,IK3)
            K2=K2+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
            IF(KADJ.GT.NONLIN)K2=K2+IK1
            K3=K3+IK2
          END DO
        END IF
C**NEXT K
        K=KNEXT
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(ICOUPC.EQ.0)GO TO 801
C**CORIOLIS
        IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
          IF(IABC.LT.7)THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
            CALL MEMO(1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
            CALL GETBP1(W(LIP1),ISIZE1,NMODE,K,W(LMXBAS),NSIZE1,0)
C**READ INTO MATRIX
            CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
            CALL MATIN(W(LXA1),W(LXA1),NSIZE1,21)
            CALL VMI1(NNMODE,K,W(LXA),W(LXA1),NSIZE1,W(LIPL),W(LIPR),
     1      ISIZMX,ISIZEL,ISIZER,W(LIP1),ISIZE1)
            CALL MEMO(-2,LXA1,KXA1,LIP1,KPPP1,0,0,0,0,0,0)
          END IF
        ELSE
C****************************************  LIKE VMI2 (K+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
          IF(K.GT.NONLIN)THEN
            KPPP2=KLJP2
            IND2=1
          ELSE
            KPPP2=KIP2
            IND2=0
          END IF
          CALL MEMO(1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
          IF(K.GT.NONLIN)THEN
            CALL GETBP1(W(LIP2),LSIZE2,NMODE,K,W(LMXBAS),NSIZE2,IND2)
            MSIZE2=NSIZE2
          ELSE
            CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,NSMODE,W(LMXBAS),
     1      NSIZE2,IND2)
            MSIZE2=ISIZE2
          END IF
C**READ INTO MATRIX
          KXA2=NSIZE2*(NSIZE2+1)/2
          CALL MEMO(1,LXA1,KXA2,0,0,0,0,0,0,0,0)
          CALL MATIN(W(LXA1),W(LXA1),NSIZE2,21)
          CALL VMV1(NNMODE,KK,W(LXA),W(LXA1),NSIZE2,W(LIPL),W(LIPR),
     1    ISIZMX,ISIZEL,ISIZER,W(LIP2),MSIZE2)
          CALL MEMO(-1,LXA1,KXA2,0,0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VMI2 (K+NMODE)
        END IF
        IF(ICOUPC.EQ.1)GO TO 801
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L=0
        L2=0
        L3=0
        DO LL=1,KK-1
          IF(LL.EQ.1.AND.KK.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          IF(LCOUNT.GT.1)THEN
            LNEXT=JCONT(LCONT,LL)
          ELSE
            LNEXT=LL
          END IF
          IF(LNEXT-L.GT.1)THEN
C**UPDATE L3 FOR MISSING MODES
            DO LADJ=L+1,LNEXT-1
              CALL INTARR(W(LNBF),W(LMBF),LADJ,IL1,IL2,IL3)
              L2=L2+IL1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
              IF(LADJ.GT.NONLIN)L2=L2+IL1
              L3=L3+IL2
            END DO
          END IF
C**NEXT L
          L=LNEXT
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**CORIOLIS
          IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1    (NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
            CALL MEMO(1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
            CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,L,W(LMXBAS),NSIZE2,0)
C**READ INTO MATRIX
            CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
            CALL MATIN(W(LXA2),W(LXA2),NSIZE2,22)
            CALL VMI2(NNMODE,K,L,W(LXA),W(LXA2),NSIZE2,W(LIPL),W(LIPR),
     1      ISIZMX,ISIZEL,ISIZER,W(LIP2),ISIZE2)
            CALL MEMO(-2,LXA2,KXA2,LIP2,KPPP2,0,0,0,0,0,0)
          ELSE
C****************************************  LIKE VMI3 (K+L+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
            IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
              KPPP3=KLJP3
              IND3=1
            ELSE
              KPPP3=KIP3
              IND3=0
            END IF
            CALL MEMO(1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
            IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
              CALL GETBP2(W(LIP3),LSIZE3,NMODE,K,L,W(LMXBAS),NSIZE3,
     1        IND3)
              MSIZE3=NSIZE3
            ELSE
              CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,NSMODE,W(LMXBAS),
     1        NSIZE3,IND3)
              MSIZE3=ISIZE3
            END IF
C**READ INTO MATRIX
            KXA3=NSIZE3*(NSIZE3+1)/2
            CALL MEMO(1,LXA2,KXA3,0,0,0,0,0,0,0,0)
            CALL MATIN(W(LXA2),W(LXA2),NSIZE3,22)
            CALL VMV2(NNMODE,KK,LL,W(LXA),W(LXA2),NSIZE3,W(LIPL),
     1      W(LIPR),ISIZMX,ISIZEL,ISIZER,W(LIP3),MSIZE3)
            CALL MEMO(-1,LXA2,KXA3,0,0,0,0,0,0,0,0)
            CALL MEMO(-1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VMI3 (K+L+NMODE)
          END IF
          IF(ICOUPC.EQ.2)GO TO 802
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N=0
          N2=0
          N3=0
          DO NN=1,LL-1
            IF(NN.EQ.1.AND.LL.EQ.2.AND.KK.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            IF(LCOUNT.GT.1)THEN
              NNEXT=JCONT(LCONT,NN)
            ELSE
              NNEXT=NN
            END IF
            IF(NNEXT-N.GT.1)THEN
C**UPDATE N3 FOR MISSING MODES
              DO NADJ=N+1,NNEXT-1
                CALL INTARR(W(LNBF),W(LMBF),NADJ,IN1,IN2,IN3)
                N2=N2+IN1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
                IF(NADJ.GT.NONLIN)N2=N2+IN1
                N3=N3+IN2
              END DO
            END IF
C**NEXT N
            N=NNEXT
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**CORIOLIS
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1      N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
              CALL MEMO(1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
              CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,N,W(LMXBAS),NSIZE3,
     2        0)
C**READ INTO MATRIX
              CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
              CALL MATIN(W(LXA3),W(LXA3),NSIZE3,23)
              CALL VMI3(NNMODE,K,L,N,W(LXA),W(LXA3),NSIZE3,W(LIPL),
     1        W(LIPR),ISIZMX,ISIZEL,ISIZER,W(LIP3),ISIZE3)
              CALL MEMO(-2,LXA3,KXA3,LIP3,KPPP3,0,0,0,0,0,0)
            ELSE
C****************************************  LIKE VMI4 (K+L+N+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
              IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                KPPP4=KLJP4
                IND4=1
              ELSE
                KPPP4=KIP4
                IND4=0
              END IF
              CALL MEMO(1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
              IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                CALL GETBP3(W(LIP4),LSIZE4,NMODE,K,L,N,W(LMXBAS),
     1          NSIZE4,IND4)
                MSIZE4=NSIZE4
              ELSE
                CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,NSMODE,
     1          W(LMXBAS),NSIZE4,IND4)
                MSIZE4=ISIZE4
              END IF
C**READ INTO MATRIX

!              KXA4=NSIZE4*(NSIZE4+1)/2
              if (mod(nsize4,2) .eq. 0) then
                 kxa4 = nsize4 / 2
                 kxa4 = kxa4 * (nsize4 + 1)
              else
                 kxa4 = (nsize4 + 1) / 2
                 kxa4 = kxa4 * nsize4
              end if

              CALL MEMO(1,LXA3,KXA4,0,0,0,0,0,0,0,0)
              CALL MATIN(W(LXA3),W(LXA3),NSIZE4,23)
              CALL VMV3(NNMODE,KK,LL,NN,W(LXA),W(LXA3),NSIZE4,W(LIPL),
     1        W(LIPR),ISIZMX,ISIZEL,ISIZER,W(LIP4),MSIZE4)
              CALL MEMO(-1,LXA3,KXA4,0,0,0,0,0,0,0,0)
              CALL MEMO(-1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VMI4 (K+L+N+NMODE)
            END IF
            IF(ICOUPC.EQ.3)GO TO 803
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M=0
            M2=0
            M3=0
            DO MM=1,NN-1
              IF(MM.EQ.1.AND.NN.EQ.2.AND.LL.EQ.3.AND.KK.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              IF(LCOUNT.GT.1)THEN
                MNEXT=JCONT(LCONT,MM)
              ELSE
                MNEXT=MM
              END IF
              IF(MNEXT-M.GT.1)THEN
C**UPDATE M3 FOR MISSING MODES
                DO MADJ=M+1,MNEXT-1
                  CALL INTARR(W(LNBF),W(LMBF),MADJ,IM1,IM2,IM3)
                  M2=M2+IM1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
                  IF(MADJ.GT.NONLIN)M2=M2+IM1
                  M3=M3+IM2
                END DO
              END IF 
C**NEXT M
              M=MNEXT
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**CORIOLIS
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
                CALL MEMO(1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
                CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,M,W(LMXBAS),
     2          NSIZE4,0)
C**READ INTO MATRIX
                CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                CALL MATIN(W(LXA4),W(LXA4),NSIZE4,24)
                CALL VMI4(NNMODE,K,L,N,M,W(LXA),W(LXA4),NSIZE4,W(LIPL),
     1          W(LIPR),ISIZMX,ISIZEL,ISIZER,W(LIP4),ISIZE4)
                CALL MEMO(-2,LXA4,KXA4,LIP4,KPPP4,0,0,0,0,0,0)
              ELSE
C****************************************  LIKE VMI5 (K+L+N+M+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR FIVE MODES
                IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1          M.GT.NONLIN)THEN
                  KPPP5=KLJP5
                  IND5=1
                ELSE
                  KPPP5=KIP5
                  IND5=0
                END IF
                CALL MEMO(1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
                IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1          M.GT.NONLIN)THEN
                  CALL GETBP4(W(LIP5),LSIZE5,NMODE,K,L,N,M,W(LMXBAS),
     1            NSIZE5,IND5)
                  MSIZE5=NSIZE5
                ELSE
                  CALL GETBP5(W(LIP5),ISIZE5,NMODE,K,L,N,M,NSMODE,
     1            W(LMXBAS),NSIZE5,IND5)
                  MSIZE5=ISIZE5
                END IF
C**READ INTO MATRIX
                KXA5=NSIZE5*(NSIZE5+1)/2
                CALL MEMO(1,LXA4,KXA5,0,0,0,0,0,0,0,0)
                CALL MATIN(W(LXA4),W(LXA4),NSIZE5,24)
                CALL VMV4(NNMODE,KK,LL,NN,MM,W(LXA),W(LXA4),NSIZE5,
     1          W(LIPL),W(LIPR),ISIZMX,ISIZEL,ISIZER,W(LIP5),MSIZE5)
                CALL MEMO(-1,LXA4,KXA5,0,0,0,0,0,0,0,0)
                CALL MEMO(-1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VMI5 (K+L+N+M+NMODE)
              END IF
              IF(ICOUPC.EQ.4)GO TO 804
C**5-MODE COUPLING HERE IF NEEDED
804   CONTINUE
              M2=M2+IM1
              M3=M3+IM2
            END DO
803   CONTINUE
            N2=N2+IN1
            N3=N3+IN2
          END DO
802   CONTINUE
          L2=L2+IL1
          L3=L3+IL2
        END DO
801   CONTINUE
        K2=K2+IK1
        K3=K3+IK2
      END DO
C     CALL DUMVM(W(LXA),ISIZEL,ISIZER,33)
      CALL DUMVM(W(LXA),ISIZEL,ISIZER,IND)
      CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
C9999  CONTINUE
5503  CONTINUE
5502  CONTINUE
5501  CONTINUE
5500  CONTINUE
3300  CONTINUE
3330  CONTINUE
C5005  CONTINUE
C5000  CONTINUE

      ITIM=0

C**RESTORE IMPORTANT QUANTITIES
      NAMODE=NAMODS
      NNMODE=NNMODS
      NVMODE=NVMODS
C     ICOUPL=ICOUPS ??
C     ICOUPC=ICOUCS ??
C     JCOUPL=JCOUPS ??
C     JCOUPC=JCOUCS ??
      JREACT=JREACS
      NREACT=NREACS
C*********************************************
C*********************************************
C**COMBINE TWO CONTRACTION SCHEMES IF REQUIRED
C*********************************************
C*********************************************
      ISIZXX=ISIZMX
      LCOUNT=-LCOUNT+1
      IF(LCOUNT.LT.0)THEN
        WRITE(IOUT,500)
C**********************************************************************
C**********************************************************************
C**                               WRITE DISCS AGAIN FOR COMPLETE SYSTEM
C**********************************************************************
C**********************************************************************
        I91=91
        I92=92
        I93=93
        I94=94
C??????????????????????????????
CCCC    CALL DISCDP(W,KCONT,0)
C??????????????????????????????
        DO I=1,2
          WRITE(IOUT,420)I,(NCVAL(I,J),J=1,NVSYM)
        END DO
        DO I=1,2
          WRITE(IOUT,425)I,(ISIZC(I,J),J=1,NVSYM)
        END DO
        ISIZXX=MAX0(ISIZM1,ISIZM2)
        NVALCF=MAX0(NVAL1,NVAL2)
CC      CALL MEMO(1,LASSIG,ISIZXX*KEL21*NVSYM*3,0,0,0,0,0,0,0,0)
CC      CALL MEMO(2,LCFS,ISIZXX*NVALCF*KEL21*NVSYM*2,LEVCI,
CC   1  NVALCF*KEL21*NVSYM*2,0,0,0,0,0,0)
C**RELOAD COEFFICIENTS
        REWIND 30
        CALL CFENER(W(LCFS),W(LEVCI),ISIZXX,KEL21,NVSYM,NVALCF,1)
        CALL CFENER(W(LCFS),W(LEVCI),ISIZXX,KEL21,NVSYM,NVALCF,2)
        DO ICPL=1,ICOUPL
          FACTOR(ICPL)=1.D0
          DO I=2,NAMODE
            FACTOR(ICPL)=FACTOR(ICPL)*I
          END DO
          DO I=1,ICPL
            FACTOR(ICPL)=FACTOR(ICPL)/I
          END DO 
          IDENOM=NAMODE-ICPL
          DO I=1,IDENOM
            FACTOR(ICPL)=FACTOR(ICPL)/I
          END DO
          IF(JREACT.GT.0)FACTOR(ICPL)=FACTOR(ICPL)+1
        END DO
C**********************************************************************
C***************************************************
C***************************************************
C**INITIALLY NO LANCZOS
        LANCZ=.FALSE.
C***************************************************
C***************************************************
C**BUILD MATRIX FROM TWO EARLIER CONTRACTION SCHEMES
C**FOR EACH OF SCHEME 1 SYMMETRIES, THERE ARE NVAL1 STORED FUNCTIONS
C**FOR EACH OF SCHEME 2 SYMMETRIES, THERE ARE NVAL2 STORED FUNCTIONS
        ITIM=-1
C**********************************************************************
C**   SET UP FILES  (31)-(39) FOR COMMON INTEGRALS
C**********************************************************************
C       REWIND 31
C       REWIND 32
C       REWIND 33
C       REWIND 34
C       REWIND 35
C       REWIND 36
C       REWIND 37
C       REWIND 38
C       REWIND 39
C**********************************************************************
C**                                          LOOP OVER SYMMETRIES (LHS)
C**********************************************************************
        ICONDP=0
        DO 4400 NSL=1,NVSYM

C**NVAL1,NVAL2 ARE SIZES OF TWO EARLIER CONTRACTION SCHEMES
C**JUST IN CASE SINGLE CONTRACTION SCHEME REQUESTED (BAD USE OF 'MM')
        IF(NVAL2.EQ.0)NVAL2=1
        JSIZEL=NVAL1*NVAL2*NVSYM
        KJPL=JSIZEL*2
C**GET MAXIMUM (2-DIM) BASIS
        CALL MEMO(1,LJPL,KJPL,0,0,0,0,0,0,0,0)

C**DETERMINE TOTAL SIZE OF SYMMETRY-ADAPTED MATRIX
        ISIZEL=0
C**JUST IN CASE SINGLE CONTRACTION SCHEME REQUESTED (BAD USE OF 'MM')
        IF(ISIZC(2,1).EQ.0)ISIZEL=NCVAL(1,1)
        IF(NSL.EQ.1)THEN
          DO I=1,NVSYM
            IF(ISIZC(1,I)*ISIZC(2,I).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPL),JSIZEL,NCVAL(1,I),NCVAL(2,I),ISIZEL)
            END IF
          END DO
        END IF
        IF(NSL.EQ.2)THEN
          K=1
          DO I=1,NVSYM
            IF(ISIZC(1,I)*ISIZC(2,I+K).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPL),JSIZEL,NCVAL(1,I),NCVAL(2,I+K),
     1        ISIZEL)
            END IF
            K=-K
          END DO
        END IF
        IF(NSL.EQ.3)THEN
          K=1
          DO I=1,NVSYM
            IF(ISIZC(1,I)*ISIZC(2,I+2*K).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+2*K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPL),JSIZEL,NCVAL(1,I),NCVAL(2,I+2*K),
     1        ISIZEL)
            END IF
            KK=MOD(I,2)+1
            K=K*(-1)**KK
          END DO
        END IF
        IF(NSL.EQ.4)THEN
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
C           K=NVSYM+1-I
            K=5+NVSYM*KSR-I
            IF(ISIZC(1,I)*ISIZC(2,K).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPL),JSIZEL,NCVAL(1,I),NCVAL(2,K),ISIZEL)
            END IF
          END DO
        END IF
        IF(NSL.EQ.5)THEN
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            K=I+4*LSR
            IF(ISIZC(1,I)*ISIZC(2,K).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPL),JSIZEL,NCVAL(1,I),NCVAL(2,K),ISIZEL)
            END IF
          END DO
        END IF
        IF(NSL.EQ.6)THEN
          K=1
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            IF(ISIZC(1,I)*ISIZC(2,I+K+4*LSR).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+K+4*LSR)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPL),JSIZEL,NCVAL(1,I),NCVAL(2,I+K+4*LSR),
     1        ISIZEL)
            END IF
            K=-K
          END DO
        END IF
        IF(NSL.EQ.7)THEN
          K=1
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            IF(ISIZC(1,I)*ISIZC(2,I+2*K+4*LSR).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+2*K+4*LSR)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPL),JSIZEL,NCVAL(1,I),
     1        NCVAL(2,I+2*K+4*LSR),ISIZEL)
            END IF
            KK=MOD(I,2)+1
            K=K*(-1)**KK
          END DO
        END IF
        IF(NSL.EQ.8)THEN
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            K=5+NVSYM*KSR-I
            IF(ISIZC(1,I)*ISIZC(2,K+4*LSR).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,K+4*LSR)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPL),JSIZEL,NCVAL(1,I),NCVAL(2,K+4*LSR),
     1        ISIZEL)
            END IF
          END DO
        END IF
        IF(ISIZEL.EQ.0)GO TO 4401
CC      CALL PUTJP(W(LJPL),JSIZE,W(LIPL),ISIZMX,NMODE,NSL)
C**********************************************************************
C**                                          LOOP OVER SYMMETRIES (RHS)
C**********************************************************************
        DO 4402 NSR=NSL,NVSYM

C**NVAL1,NVAL2 ARE SIZES OF TWO EARLIER CONTRACTION SCHEMES
C**JUST IN CASE SINGLE CONTRACTION SCHEME REQUESTED (BAD USE OF 'MM')
        IF(NVAL2.EQ.0)NVAL2=1
        JSIZER=NVAL1*NVAL2*NVSYM
        KJPR=JSIZER*2
C**GET MAXIMUM (2-DIM) BASIS
        CALL MEMO(1,LJPR,KJPR,0,0,0,0,0,0,0,0)

C**DETERMINE TOTAL SIZE OF SYMMETRY-ADAPTED MATRIX
        ISIZER=0
C**JUST IN CASE SINGLE CONTRACTION SCHEME REQUESTED (BAD USE OF 'MM')
        IF(ISIZC(2,1).EQ.0)ISIZE=NCVAL(1,1)
        IF(NSR.EQ.1)THEN
          DO I=1,NVSYM
            IF(ISIZC(1,I)*ISIZC(2,I).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPR),JSIZER,NCVAL(1,I),NCVAL(2,I),ISIZER)
            END IF
          END DO
        END IF
        IF(NSR.EQ.2)THEN
          K=1
          DO I=1,NVSYM
            IF(ISIZC(1,I)*ISIZC(2,I+K).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPR),JSIZER,NCVAL(1,I),NCVAL(2,I+K),
     1        ISIZER)
            END IF
            K=-K
          END DO
        END IF
        IF(NSR.EQ.3)THEN
          K=1
          DO I=1,NVSYM
            IF(ISIZC(1,I)*ISIZC(2,I+2*K).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+2*K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPR),JSIZER,NCVAL(1,I),NCVAL(2,I+2*K),
     1        ISIZER)
            END IF
            KK=MOD(I,2)+1
            K=K*(-1)**KK
          END DO
        END IF
        IF(NSR.EQ.4)THEN
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
C           K=NVSYM+1-I
            K=5+NVSYM*KSR-I
            IF(ISIZC(1,I)*ISIZC(2,K).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPR),JSIZER,NCVAL(1,I),NCVAL(2,K),ISIZER)
            END IF
          END DO
        END IF
        IF(NSR.EQ.5)THEN
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            K=I+4*LSR
            IF(ISIZC(1,I)*ISIZC(2,K).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPR),JSIZER,NCVAL(1,I),NCVAL(2,K),ISIZER)
            END IF
          END DO
        END IF
        IF(NSR.EQ.6)THEN
          K=1
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            IF(ISIZC(1,I)*ISIZC(2,I+K+4*LSR).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+K+4*LSR)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPR),JSIZER,NCVAL(1,I),NCVAL(2,I+K+4*LSR),
     1        ISIZER)
            END IF
            K=-K
          END DO
        END IF
        IF(NSR.EQ.7)THEN
          K=1
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            IF(ISIZC(1,I)*ISIZC(2,I+2*K+4*LSR).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+2*K+4*LSR)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPR),JSIZER,NCVAL(1,I),
     1        NCVAL(2,I+2*K+4*LSR),ISIZER)
            END IF
            KK=MOD(I,2)+1
            K=K*(-1)**KK
          END DO
        END IF
        IF(NSR.EQ.8)THEN
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            K=5+NVSYM*KSR-I
            IF(ISIZC(1,I)*ISIZC(2,K+4*LSR).NE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,K+4*LSR)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJPR),JSIZER,NCVAL(1,I),NCVAL(2,K+4*LSR),
     1        ISIZER)
            END IF
          END DO
        END IF
        IF(ISIZER.EQ.0)GO TO 4403
CC      CALL PUTJP(W(LJPR),JSIZER,W(LIPR),ISIZMX,NMODE,NSR)
        IF(ICOUPC.GE.0)THEN
          REWIND 21
        END IF
        IF(ICOUPC.GT.1)THEN
          REWIND 22
        END IF
        IF(ICOUPC.GT.2)THEN
          REWIND 23
        END IF
        IF(ICOUPC.GT.3)THEN
          REWIND 24
        END IF
C       DO 9898 IABC=1,9
        ITIM=ITIM+1
        IF(ICOUPC.GT.0)REWIND 91
        IF(ICOUPC.GT.1)REWIND 92
        IF(ICOUPC.GT.2)REWIND 93
        IF(ICOUPC.GT.3)REWIND 94
        KXA=ISIZEL*ISIZER
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGZ(ISIZEL,ISIZER,W(LXA),IDUM,W(LXA),0)
        IF(IABC.LT.7)THEN
          IHERM=1
        ELSE
          IHERM=-1
        END IF
        CALL VMI0(W(LXA),ISIZEL,ISIZER,IABC)
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
        K2=0
        K3=0
        DO K=1,NAMODE
          IF(K.EQ.1.AND.ITIM.EQ.0)THEN
            ITIM1A=0
            ITIM1B=0
          END IF
          K1=K-1
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'K'
          DO KK=1,ICONT(1)
            IF(K.EQ.JCONT(1,KK))THEN
              KKCONT=1
              KMODX=KK
              KONBAS=0
              KKXK=ICSIZ1
              KKIP=IPSIZ1
              KVAL=NVAL1
            END IF
          END DO
          DO KK=1,ICONT(2)
            IF(K.EQ.JCONT(2,KK))THEN
              KKCONT=2
              KMODX=KK
              KONBAS=KLP
              KKXK=ICSIZ2
              KKIP=IPSIZ2
              KVAL=NVAL2
            END IF
          END DO
          IF(ICOUPC.EQ.0)GO TO 5551
C**MATRIX FOR CENTRAL (DGEMM) INTEGRALS FOR THIS SCHEME
          CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
C******************************************************************
C**CORIOLIS AND POTENTIAL
          IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
            IF(IABC.LT.7)THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
              NREC(1)=NREC(1)+1
              CALL MEMO(1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
              CALL GETBP1(W(LIP1),ISIZE1,NMODE,K,W(LMXBAS),NSIZE1,0)
              KXA1=NSIZE1*(NSIZE1+1)/2
              CALL MEMO(3,LXA1,KXA1,LTEMP,KKXK*KVAL,LWRK,
     1        KVAL*KVAL,0,0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(NSIZE1,NSIZE1,W(LXA1),NSIZE1,W(LXA1),ICI)
              CALL VCMI1(NAMODE,NNMODE,K,KMODX,W(LH+K2),W(LXQ+K3),
     1        IK3,IK2,W(LJPL),W(JIPR),ISIZMX,W(LKP+KONBAS),KKXK,KKIP,
     2        KKCONT,W(LXA),ISIZEL,ISIZER,W(LIP1),ISIZE1,W(LXA1),
     3        W(LXA1),NSIZE1,W(LXK),W(LTEMP),W(LWRK),KVAL,ISTART,IEND,
     4        W(LVM1),W(LVM1),J21,IABC,W(LMODNT),NSL,NSR,
     5        W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
CCC  6        W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     7        ISIZXX,NVALCF,KEL21,NVSYM)
              CALL MEMO(-3,LXA1,KXA1,LTEMP,KKXK*KVAL,LWRK,
     1        KVAL*KVAL,0,0,0,0)
              CALL MEMO(-1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
            END IF
          END IF
          CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
          IF(ICOUPC.EQ.1)GO TO 5551

5551   CONTINUE
          K2=K2+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(K.GT.NONLIN)K2=K2+IK1
          K3=K3+IK2
        END DO
C       CALL DUMVM(W(LXA),ISIZEL,ISIZER,30+IABC)
        CALL DUMVM(W(LXA),ISIZEL,ISIZER,33)
        CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
C9898    CONTINUE
        ICONDP=1
4403    CONTINUE
        CALL MEMO(-1,LJPR,KJPR,0,0,0,0,0,0,0,0)
4402    CONTINUE
4401    CONTINUE
        CALL MEMO(-1,LJPL,KJPL,0,0,0,0,0,0,0,0)
4400    CONTINUE
      END IF
C*********************************************
C*********************************************
C**COMBINE TWO CONTRACTION SCHEMES IF REQUIRED
C*********************************************
C*********************************************

C**********************************************************************
C**                                 SET UP ROVIBRATIONAL MATRIX FOR K>0
C**********************************************************************
C     IF(NVSYM.EQ.1)THEN
C**H2CO C1
C       NRSYM=1
C**H2O
C       IF(TRIAT)NRSYM=2
C     END IF
C     IF(NVSYM.EQ.2)THEN
C**H2CO Cs
C       NRSYM=2
C**H2O
C       IF(TRIAT)NRSYM=4
C     END IF
C     IF(NVSYM.EQ.4)THEN
C**H2CO C2v
C**H2O2 C2v IF JREACT>0
C       NRSYM=4
C     END IF
C**NRSYM IS TOTAL NUMBER ROTATIONAL SYMMETRIES
C     NRSYM=NVSYM

C**INDEX IS A COUNT OF THE ACTUAL ROTATIONAL SYMMETRIES REQUIRED
      INDEX=1
      DO 2222 NROT=1,NRSYM
C**********************************************************************
C**                                                LOOP OVER SYMMETRIES
C**********************************************************************
C**SKIP IF FINISHED
      IF(INDEX.GT.MVSYM)GO TO 2222
C**SKIP UNTIL FIND NEXT SYMMETRY THAT IS REQUIRED
      IF(NROT.NE.MWSYM(INDEX))GO TO 2222
      INDEX=INDEX+1
C     JS=MOD(JTHIS,2)+1
C     WRITE(IOUT,375)NROT,ISYMP(NROT,ISYMD(1,JS))
      IF(TRIAT)THEN
C**START B2 (IF EXISTS)
        NS=2
        IF(NROT.LE.NVSYM)THEN
          ISIZE=(JTHIS+1)*NVAL
          K21=JTHIS+1
C**ROTATION ELEMENTS 1,3,5,...
          KOFF=1
C**START A1 (A')
          IF(NROT.EQ.1)NS=1
        ELSE
          ISIZE=JTHIS*NVAL
          K21=JTHIS
C**ROTATION ELEMENTS 2,4,6,...
          KOFF=2
C**START A1 (A')
          IF(NROT.EQ.NVSYM+1)NS=1
        END IF
      ELSE
C**NS0 WILL BE H2CO SYMMETRY THROUGHOUT
        NS0=NROT
        NS=NROT
        K21=2*JTHIS+1
        ISIZE=K21*NVAL
      END IF
C**ISTART, IEND ARE START AND END COLUMNS
      ISTART=1
      KSTART=1
      IF(LANCZ)THEN
        CALL MEMO(3,LWK,ISIZE,LYK,ISIZE,LZK,ISIZE,0,0,0,0)
        CALL MEMO(2,LANCNT,ISIZE,LANK,ISIZE,0,0,0,0,0,0)
        KLAN=LANMAX*(LANMAX+1)/2
        LSIZE=1
        KSIZE=1
        IEND=ISIZE
      ELSE
        KXA=ISIZE*(ISIZE+1)/2
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
        IF(IABC.EQ.1)THEN
C**ZEROISE MATRIX FIRST TIME ROUND
          CALL DIAGZ(ISIZE,ISIZE,W(LXA),ISIZE,W(LXA),ICI)
        ELSE
C**LOAD PARTIAL MATRIX FOR NROT FROM (59)
C         CALL ROTMAT(W(LXA),KXA,0)
        END IF
        KEND=K21
      END IF
C
C**MODIFY ISIZE FOR 'MISSING' BLOCKS (IF ANY)
      DO KROT=KSTART,KEND
        KROTT=KROT
8       IF(KROTT.GT.4)THEN
          KROTT=KROTT-4
          GO TO 8
        END IF
C**NROT IS ACTUAL (EVEN-J) SYMMETRY
        NS=ISYMP(ISYMD(KROTT,1),NROT)
C**NS IS ACTUAL VIBRATIONAL SYMMETRY
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KA=KROT/2
        IF(MOD(KA,2).NE.0)NS=ISYMP(NS,ISYMT)
C**NS IS STORED SYMMETRY
        IF(NTOT(NS).EQ.0)ISIZE=ISIZE-NVAL
      END DO
C     WRITE(IOUT,380)ISIZE,CUT
C**WRITE DUMP FILE
C     WRITE(INP5)ISIZE
C**LOAD PARTIAL MATRIX FOR NROT FROM (59)
      IF(ISIZE.NE.0.AND.IABC.NE.1)CALL ROTMAT(W(LXA),KXA,0)
C
      ISKIP=0
7777  CONTINUE
      IF(LANCZ)THEN
C***********************************************
C**REPOSITIONING IF TOO BIG MISSING (ISKIP.NE.0)
C***********************************************
        IEND=ISTART-1
        KEND=KSTART-1
        KXA=0
        KCOL=0
        DO I=ISTART,ISIZE
          KXA=KXA+LSIZE
C**CAN'T GET IT ALL IN
          IF(KXA.GT.KLAN)THEN
C**SUBTRACT TOTAL SIZE 'K' BLOCK SO FAR
            KXA=KXA-KSIZE
C**SET LSIZE TO VALUE FOR START OF 'K' BLOCK
            LSIZE=LSIZE-KCOL
C**STARTING VALUE OF KSIZE FOR NEXT 'K' BLOCK (SIZE FIRST COLUMN)
            KSIZE=LSIZE
C**RESET IEND TO LAST COLUMN PREVIOUS BLOCK
            IEND=IEND-KCOL
            GO TO 7776
          END IF
C**UPDATE FOR NEXT COLUMN
          LSIZE=LSIZE+1
C**KSIZE: TOTAL SIZE 'K' BLOCK INCLUDING NEXT COLUMN
          KSIZE=KSIZE+LSIZE
C**IEND: CURRENT END COLUMN
          IEND=IEND+1
C**KCOL: NUMBER COLUMNS IN 'K' BLOCK SO FAR
          KCOL=KCOL+1
          IF(MOD(KCOL,NVAL).EQ.0)THEN
C**STARTING VALUES FOR NEXT BLOCK
            KCOL=0
            KSIZE=LSIZE
C**KEND: CURRENT KROT(LHS) FOR COMPLETED BLOCK
            KEND=KEND+1
          END IF
        END DO
7776    CONTINUE
        IF(IEND.EQ.0)STOP 'LANMAX TOO SMALL'
        WRITE(IOUT,*)'ISTART,IEND,KSTART,KEND,KXA',
     1  ISTART,IEND,KSTART,KEND,KXA
        CALL FLUSH(IOUT)
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGL(KXA,W(LXA))
      END IF
      CALL MEMO(2,LTEMP,ISIZMX*NVALCF,LXK,ISIZMX*ISIZMX,0,0,0,0,0,0)
C     CALL MEMO(3,LTEMP,ISIZMX*NVALCF,LXK,ISIZMX,LZK,NVALCF*NVALCF,
C    10,0,0,0)
C**********************************************************************
C**                                                       LOOP ROUND Ka
C**********************************************************************
      IOFF=0
      DO 8888 KROT=KSTART,KEND
      KROTL=KROT
      KROTTL=KROTL
1     IF(KROTTL.GT.4)THEN
        KROTTL=KROTTL-4
        GO TO 1
      END IF
C**NROT IS ACTUAL (EVEN-J) SYMMETRY
      NSL=ISYMP(ISYMD(KROTTL,1),NROT)
C**NSL IS ACTUAL SYMMETRY
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
      KAL=KROTL/2
      IF(MOD(KAL,2).NE.0)NSL=ISYMP(NSL,ISYMT)
C**NSL IS STORED SYMMETRY
      IF(TRIAT)NSL=NS
C**WRITE ACTUAL SYMMETRY
C     IF(MOD(KAL,2).NE.0)THEN
C       WRITE(IOUT,390)KROT,NSL,ISYMP(NSL,ISYMT)
C     ELSE
C       WRITE(IOUT,390)KROT,NSL,NSL
C     END IF
C**'MISSING' LHS BLOCK
      IF(NTOT(NSL).EQ.0)GO TO 8888
      IF(TRIAT)THEN
        JROT=2*(KROT-1)+KOFF
      ELSE
        JROT=KROT
      END IF
C**SET UP DIAGONAL ELEMENTS
      IF(IABC.EQ.1)THEN
      CALL DIAGEL(ISIZE,KROT,W(LXA),W(LEVCI),NVALCF,KEL21,JROT,NSL,
     1ISTART,KSTART,IOFF)
      END IF
      IOFF=IOFF+NVAL
C**ISIZEL IS SIZE LHS CI BASIS
C**ISIZER IS SIZE RHS CI BASIS
C**NSL IS SYMMETRY LHS
C**NSR IS SYMMETRY RHS
C**KROTL DENOTES LHS BLOCK OFFSET
C**KROTR DENOTES RHS BLOCK OFFSET
C**JROTL IS LHS ROTATIONAL FUNCTION
C**JROTR IS RHS ROTATIONAL FUNCTION
C**IOFF IS OFFSET FOR EACH NON-ZERO BLOCK OF LENGTH NVAL
C***************************************************************
C**START WITH DIAGONAL BLOCKS .... K/K (K/K+2 IF K=1)
C***************************************************************
      ISIZEL=NTOT(NSL)
      IF(TRIAT)THEN
        JROTL=2*(KROTL-1)+KOFF
      ELSE
        JROTL=KROTL
      END IF
      JROTR=JROTL
      KSIGN=1
      IF(KSTEP.NE.1)THEN
        NSR=NSL
        ISIZER=ISIZEL
        KROTR=KROTL
        IELL=1
        IELX=(JROTL-1)/KSTEP
        IELL=IELL+IELX*KSTEP
        IELR=IELL
        KSIGN=-1
C**NEGATE OLD DIAGONAL ELEMENTS FOR Jx**2 + Jy**2 + Jz**2 (1, 2 AND 3)
C**Jx**2
        IF(IABC.EQ.1)
     &  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,IELL,IELR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,1)
C**Jy**2
        IF(IABC.EQ.2)
     &  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,IELL,IELR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,2)
C**Jz**2
        IF(IABC.EQ.3)
     &  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,IELL,IELR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,3)
C**SET UP NEW DIAGONAL ELEMENTS FOR Jx**2 + Jy**2 + Jz**2 (1, 2 AND 3)
        KSIGN=1
C**Jx**2
        IF(IABC.EQ.1)
     &  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,1)
C**Jy**2
        IF(IABC.EQ.2)
     &  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,2)
C**Jz**2
        IF(IABC.EQ.3)
     &  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,3)
      END IF
C***************************************************************
C**OFF-DIAGONAL BLOCKS .... K/K+1
C***************************************************************
      IF(TRIAT)THEN
        NSR=NSL+1
        IF(NSR.GT.NVSYM)NSR=1
        KINCR=1
        KROTR=KROTL+KINCR
        JROTR=2*(KROTR-1)+KOFF
      ELSE
        KINCR=3
        KROTR=KROTL+2
        KROTTR=KROTR
2       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 2
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        JROTR=JROTL+2
      END IF
      ISIZER=NTOT(NSR)
      IF(MOD(KROT,2).EQ.0.AND.KROT+KINCR.LE.K21)THEN
C**TEMPORARY
      IF(IABC.EQ.4.OR.IABC.EQ.8)THEN
      WRITE(IOUT,*)'IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR'
      WRITE(IOUT,*)IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR
      END IF
C**TEMPORARY
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jx (7)
        IF(ISIZER.NE.0.AND.IABC.EQ.8)
     1  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     4  JROTL,JROTR,NSL,NSR,ISTART,KSTART,8)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JyJz]+ (5)
        IF(ISIZER.NE.0.AND.IABC.EQ.4)
     1  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     4  JROTL,JROTR,NSL,NSR,ISTART,KSTART,4)
        IF(.NOT.TRIAT)THEN
          KROTR=KROTL+KINCR
          KROTTR=KROTR
3         IF(KROTTR.GT.4)THEN
            KROTTR=KROTTR-4
            GO TO 3
          END IF
          NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
          KAR=KROTTR/2
          IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
          JROTR=JROTL+KINCR
          ISIZER=NTOT(NSR)
C**TEMPORARY
      IF(IABC.EQ.5.OR.IABC.EQ.7)THEN
      WRITE(IOUT,*)'IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR'
      WRITE(IOUT,*)IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR
      END IF
C**TEMPORARY
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jy (8)
          IF(ISIZER.NE.0.AND.IABC.EQ.7)
     1    CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2    (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3    ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     4    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,7)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJz]+ (4)
          IF(ISIZER.NE.0.AND.IABC.EQ.5)
     1    CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2    (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3    ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     4    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,5)
        END IF
      END IF
      IF(.NOT.TRIAT)THEN
        KINCR=2
        KROTR=KROTL+KINCR
        KROTTR=KROTR
4       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 4
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        JROTR=JROTL+KINCR
        ISIZER=NTOT(NSR)
      END IF
      IF(MOD(KROT,2).EQ.1.AND.KROT+KINCR.LE.K21)THEN
C**TEMPORARY
      IF(IABC.EQ.4.OR.IABC.EQ.8)THEN
      WRITE(IOUT,*)'IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR'
      WRITE(IOUT,*)IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR
      END IF
C**TEMPORARY
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jx (7)
        IF(ISIZER.NE.0.AND.IABC.EQ.8)
     1  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     4  JROTL,JROTR,NSL,NSR,ISTART,KSTART,8)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JyJz]+ (5)
        IF(ISIZER.NE.0.AND.IABC.EQ.4)
     1  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     4  JROTL,JROTR,NSL,NSR,ISTART,KSTART,4)
        IF(.NOT.TRIAT)THEN
          KROTR=KROTL+1
          KROTTR=KROTR
5         IF(KROTTR.GT.4)THEN
            KROTTR=KROTTR-4
            GO TO 5
          END IF
          NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
          KAR=KROTTR/2
          IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
          JROTR=JROTL+1
          ISIZER=NTOT(NSR)
C**TEMPORARY
      IF(IABC.EQ.5.OR.IABC.EQ.7)THEN
      WRITE(IOUT,*)'IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR'
      WRITE(IOUT,*)IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR
      END IF
C**TEMPORARY
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jy (8)
          IF(ISIZER.NE.0.AND.IABC.EQ.7)
     1    CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2    (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3    ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     4    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,7)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJz]+ (4)
          IF(ISIZER.NE.0.AND.IABC.EQ.5)
     1    CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2    (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3    ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     4    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,5)
        END IF
      END IF
C***************************************************************
C**OFF-DIAGONAL BLOCKS .... K/K+2
C***************************************************************
      NSR=NSL
      IF(TRIAT)THEN
        KINCR=2
        KROTR=KROTL+KINCR
        JROTR=2*(KROTR-1)+KOFF
      ELSE
        KINCR=4
        KROTR=KROTL+KINCR
        JROTR=JROTL+KINCR
      END IF
      ISIZER=NTOT(NSR)
      IF(KROT+KINCR.LE.K21)THEN
C**TEMPORARY
      IF(IABC.EQ.1.OR.IABC.EQ.2)THEN
      WRITE(IOUT,*)'IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR'
      WRITE(IOUT,*)IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR
      END IF
C**TEMPORARY
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jx**2 + Jy**2 (1 AND 2)
C**Jx**2
        IF(ISIZER.NE.0.AND.IABC.EQ.1)
     1  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     4  JROTL,JROTR,NSL,NSR,ISTART,KSTART,1)
C**Jy**2
        IF(ISIZER.NE.0.AND.IABC.EQ.2)
     1  CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2  (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3  ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     4  JROTL,JROTR,NSL,NSR,ISTART,KSTART,2)
      END IF
C***************************************************************
C**OFF-DIAGONAL BLOCKS .... Jz AND [JxJy]+
C***************************************************************
      IF(.NOT.TRIAT)THEN
        KINCR=1
        KROTR=KROTL+KINCR
        KROTTR=KROTR
6       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 6
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        JROTR=JROTL+KINCR
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).EQ.0)THEN
C**TEMPORARY
      IF(IABC.EQ.6.OR.IABC.EQ.9)THEN
      WRITE(IOUT,*)'IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR'
      WRITE(IOUT,*)IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR
      END IF
C**TEMPORARY
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jz (9)
          IF(ISIZER.NE.0.AND.IABC.EQ.9)
     1    CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2    (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3    ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     4    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,9)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJy]+ (6)
          IF(ISIZER.NE.0.AND.IABC.EQ.6)
     1    CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2    (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3    ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     4    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,6)
          KINCR=5
          KROTR=KROTL+KINCR
          JROTR=JROTL+KINCR
          IF(KROT+KINCR.LE.K21)THEN
C**TEMPORARY
      IF(IABC.EQ.6)THEN
      WRITE(IOUT,*)'IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR'
      WRITE(IOUT,*)IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR
      END IF
C**TEMPORARY
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJy]+ (6)
            IF(ISIZER.NE.0.AND.IABC.EQ.6)
     1      CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2      (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,
     3      W(LTEMP),ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,
     4      KROTR,JROTL,JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,6)
          END IF
        END IF
        KINCR=4
        KROTR=KROTL+3
        KROTTR=KROTR
7       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 7
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        JROTR=JROTL+3
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).NE.0.AND.KROT+KINCR.LE.K21)THEN
C**TEMPORARY
      IF(IABC.EQ.6)THEN
      WRITE(IOUT,*)'IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR'
      WRITE(IOUT,*)IABC,JROTL,JROTR,KROTL,KROTR,NSL,NSR
      END IF
C**TEMPORARY
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJy]+ (6)
          IF(ISIZER.NE.0.AND.IABC.EQ.6)
     1    CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVALCF*KEL21*
     2    (NSL-1)),W(LCFS+ISIZMX*NVALCF*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     3    ISIZEL,ISIZER,NVALCF,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     4    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,6)
        END IF
      END IF
C***************************************************************
C**UPDATE VIBRATIONAL SYMMETRY FOR NEXT TIME IF REQUIRED
C***************************************************************
      IF(TRIAT)THEN
        NS=NS+1
        IF(NS.GT.NVSYM)NS=1
      END IF
8888  CONTINUE
      CALL MEMO(-2,LTEMP,ISIZMX*NVALCF,LXK,ISIZMX*ISIZMX,0,0,0,0,0,0)
C     CALL MEMO(-3,LTEMP,ISIZMX*NVALCF,LXK,ISIZMX,LZK,NVALCF*NVALCF,
C    10,0,0,0)
      IF(ISIZE.EQ.0)GO TO 2223
      ISKIP=0
      IF(LANCZ)THEN
        CALL LANOUT(W(LXA),W(LWK),W(LWK),W(LYK),W(LZK),W(LZK),ISIZE,
     1  ISTART,IEND,W(LANCNT),W(LANK))
        CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
        IF(IEND.NE.ISIZE)THEN
          ISTART=IEND+1
          KSTART=KEND+1
          ISKIP=1
        END IF
      END IF
      IF(ISKIP.NE.0)GO TO 7777
      IF(LANCZ)CALL MEMO(-1,LZK,ISIZE,0,0,0,0,0,0,0,0)
C**DUMP IABC PART FOR THIS NROT MATRIX TO DISC 58
      CALL ROTMAT(W(LXA),KXA,1)
2223  CONTINUE
      IF(.NOT.LANCZ)CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
2222  CONTINUE
      REWIND 58
      REWIND 59
C**THERE ARE NROT (PARTIAL) MATRICES ON DISC 58
C**INDEX IS A COUNT OF THE ACTUAL ROTATIONAL SYMMETRIES REQUIRED
      INDEX=1
      DO 2224 NROT=1,NRSYM
C**********************************************************************
C**                                                LOOP OVER SYMMETRIES
C**********************************************************************
C**SKIP IF FINISHED
      IF(INDEX.GT.MVSYM)GO TO 2224
C**SKIP UNTIL FIND NEXT SYMMETRY THAT IS REQUIRED
      IF(NROT.NE.MWSYM(INDEX))GO TO 2224
      INDEX=INDEX+1
      ISIZE=K21*NVAL
      KXA=ISIZE*(ISIZE+1)/2
      CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**MODIFY ISIZE FOR 'MISSING' BLOCKS (IF ANY)
      DO KROT=KSTART,KEND
        KROTT=KROT
9       IF(KROTT.GT.4)THEN
          KROTT=KROTT-4
          GO TO 9
        END IF
C**NROT IS ACTUAL (EVEN-J) SYMMETRY
        NS=ISYMP(ISYMD(KROTT,1),NROT)
C**NS IS ACTUAL VIBRATIONAL SYMMETRY
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KA=KROT/2
        IF(MOD(KA,2).NE.0)NS=ISYMP(NS,ISYMT)
C**NS IS STORED SYMMETRY
C**TEMPORARY
C     WRITE(IOUT,*)'THIS K,SYMM = ',KROT,NS
C     CALL FLUSH(IOUT)
C**TEMPORARY
        IF(NTOT(NS).EQ.0)ISIZE=ISIZE-NVAL
      END DO
      IF(ISIZE.EQ.0)GO TO 2225
C**RELOAD MATRIX FOR THIS NROT FROM (58), WRITE IT TO (59)
      CALL ROTMAT(W(LXA),KXA,2)
2225  CONTINUE
      IF(.NOT.LANCZ)CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
2224  CONTINUE
      REWIND 58
      REWIND 59
9999  CONTINUE
5005  CONTINUE
5000  CONTINUE
C**THERE ARE NOW NROT (COMPLETE) MATRICES ON DISC 59
C**TEMPORARY
      WRITE(IOUT,*)
      WRITE(IOUT,*)'NVAL = ',NVAL
      WRITE(IOUT,*)'SIZES VIBRATIONAL BLOCKS'
      DO NROT=1,NVSYM
        WRITE(IOUT,*)NROT,NTOT(NROT)
      END DO
      WRITE(IOUT,*)
C**TEMPORARY

C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
C**INDEX IS A COUNT OF THE ACTUAL ROTATIONAL SYMMETRIES REQUIRED
      INDEX=1
      DO 2226 NROT=1,NRSYM
C**********************************************************************
C**                                                LOOP OVER SYMMETRIES
C**********************************************************************
C**SKIP IF FINISHED
      IF(INDEX.GT.MVSYM)GO TO 2226
C**SKIP UNTIL FIND NEXT SYMMETRY THAT IS REQUIRED
      IF(NROT.NE.MWSYM(INDEX))GO TO 2226
      INDEX=INDEX+1
      JS=MOD(JTHIS,2)+1
      WRITE(IOUT,375)NROT,ISYMP(NROT,ISYMD(1,JS))
C**NS0 WILL BE H2CO SYMMETRY THROUGHOUT
      NS0=NROT
      NS=NROT
      K21=2*JTHIS+1
      ISIZE=K21*NVAL
      KXA=ISIZE*(ISIZE+1)/2
      CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
      CALL DIAGZ(ISIZE,ISIZE,W(LXA),ISIZE,W(LXA),ICI)
      KEND=K21
C**MODIFY ISIZE FOR 'MISSING' BLOCKS (IF ANY)
      DO KROT=KSTART,KEND
        KROTT=KROT
808     IF(KROTT.GT.4)THEN
          KROTT=KROTT-4
          GO TO 808
        END IF
C**NROT IS ACTUAL (EVEN-J) SYMMETRY
        NS=ISYMP(ISYMD(KROTT,1),NROT)
C**NS IS ACTUAL VIBRATIONAL SYMMETRY
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KA=KROT/2
        IF(MOD(KA,2).NE.0)NS=ISYMP(NS,ISYMT)
C**NS IS STORED SYMMETRY
        IF(NTOT(NS).EQ.0)ISIZE=ISIZE-NVAL
C**WRITE ACTUAL SYMMETRY
        IF(MOD(KA,2).NE.0)THEN
          WRITE(IOUT,390)KROT,NS,ISYMP(NS,ISYMT)
        ELSE
          WRITE(IOUT,390)KROT,NS,NS
        END IF
      END DO
      WRITE(IOUT,380)ISIZE,CUT
C**WRITE DUMP FILE
      IF(ISCFCI.GT.0.AND.IDUMP.NE.0)WRITE(INP5)ISIZE
      IF(ISIZE.EQ.0)GO TO 2227
C**RELOAD MATRIX FOR THIS NROT FROM (59)
      CALL ROTMAT(W(LXA),KXA,0)
      REWIND 60
      WRITE(IOUT,385)
      IF(NVALR.GT.ISIZE)NVALR=ISIZE
      NVECR=NVALR
      KXK=ISIZE*NVALR
      KSUP4=5*ISIZE
      KWK=ISIZE
      KEVAL=NVALR
      CALL TIMIT(1)
      ICID=1
      IF(LANCZ)THEN
C       LGIV=(IGIV.NE.0)
C       KXKL=LANCYC*NCYCLE*LANCYC*NCYCLE
C       IF(LGIV)KXKL=LANCYC*NCYCLE*LANCYC
C       CALL MEMO(3,LXK,ISIZE*NVALR,LVK,ISIZE*LANCYC,
C    1  LZK,ISIZE*LANCYC,0,0,0,0)
C       CALL MEMO(4,LSK,ISIZE*LANCYC,LEN,NVALR,LEL,NVALR,LNV,NVALR,
C    1  0,0)
C       CALL MEMO(4,LWRK,LANCYC*NCYCLE,
C    1  LSUP4,5*LANCYC*NCYCLE,LEVAL,LANCYC*NCYCLE,LXKL,KXKL,0,0)
C       KXA=LANCYC*NCYCLE*(LANCYC*NCYCLE+1)/2
C       IF(LGIV)CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C       CALL MEMO(1,LXKLC,KXA,0,0,0,0,0,0,0,0)
C       WRITE(IOUT,*)'Calculating LANCZOS'
C       IF(LGIV)THEN
C         CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
C    1    ISIZE,NVALR,W(LJP),JSIZE,NMODE,W(LKP),W(LKP+KLP),W(LASSIG),
C    2    ISIZMX,KEL21,NROT,
C    2    W(LEN),W(LEL),W(LNV),W(LXA),W(LXKLC),W(LXKL),W(LWRK),
C    3    W(LSUP4),W(LEVAL),NCYCLE,NVSYM,W(LANCNT),W(LANK),LANCYC)
C       ELSE
C         CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
C    1    ISIZE,NVALR,W(LJP),JSIZE,NMODE,W(LKP),W(LKP+KLP),W(LASSIG),
C    2    ISIZMX,KEL21,NROT,
C    2    W(LEN),W(LEL),W(LNV),W(LXKL),W(LXKLC),W(LXKL),W(LWRK),
C    3    W(LSUP4),W(LEVAL),NCYCLE,NVSYM,W(LANCNT),W(LANK),LANCYC)
C       END IF
C       CALL TIMIT(3)
C       CALL FLUSH(IOUT)
C**DUMP ROTATIONAL INFO. TO (INP5)
C       IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
C         IDDDP=JDUMP(NROT)
C         IF(NVALR.LT.JDUMP(NROT))IDDDP=NVALR
C         CALL OUTCIR(ISIZE,W(LXK),W(LXK),W(LNDUMP),NROT,LDUMP,IDDDP,
C    1    W(LEL),NMODE)
C       END IF
C       CALL MEMO(-5,LXK,ISIZE*NVALR,LVK,ISIZE*LANCYC,
C    1  LZK,ISIZE*LANCYC,LWK,ISIZE,LYK,ISIZE)
C       CALL MEMO(-5,LSK,ISIZE*LANCYC,LEN,NVALR,LEL,NVALR,LNV,NVALR,
C    1  LXKL,KXKL)
C       CALL MEMO(-4,LXKLC,KXA,LWRK,LANCYC*NCYCLE,
C    1  LSUP4,5*LANCYC*NCYCLE,LEVAL,LANCYC*NCYCLE,0,0)
C       IF(LGIV)CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
C       CALL MEMO(-2,LANCNT,ISIZE,LANK,ISIZE,0,0,0,0,0,0)
C       LGIV=.FALSE.
      ELSE
        LGIV=.TRUE.
        CALL MEMO(4,LWK,KWK,LEVAL,KEVAL,LSUP4,KSUP4,LXK,KXK,0,0)
        WRITE(IOUT,*)'Calculating DIAG'
        CALL FLUSH(IOUT)
        CALL DIAG(W(LXA),W(LXK),ISIZE,ISIZE,0,W(LSUP4),W(LEVAL),W(LWK),
     1  NVALR,NVECR,W(LJP),JSIZE,NMODE,W(LXK),W(LXK),W(LASSIG),ISIZMX,
     2  KEL21,NROT)
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
C**DUMP CI COEFFICIENTS TO (INP5)
        IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
          IDDDP=JDUMP(NROT)
          IF(NVALR.LT.JDUMP(NS))IDDDP=NVALR
          CALL OUTCIR(ISIZE,W(LXK),W(LXK),W(LNDUMP),NROT,LDUMP,IDDDP,
     1    W(LWK),NMODE)
        END IF
        CALL MEMO(-4,LXK,KXK,LWK,KWK,LEVAL,KEVAL,LSUP4,KSUP4,0,0)
        LGIV=.FALSE.
      END IF
2227  CONTINUE
      IF(.NOT.LANCZ)CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
2226  CONTINUE
      ICID=0
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETSYM(IABC,KSTART,KEND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10),NSYMSL(10),MSYMSL(10,100),
     1NSYMSR(10,100),MSYMSR(10,100,10)
      COMMON/AXES/MX(3),MXDIP(5),MXROT(3),ISYMT
      COMMON/MATRIX/NVAL,NVALR,KSTEP,KSIGN,NVALCF
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/FILASS/IOUT
C**************************************************************
      K21=2*JTHIS+1
      INDEX=1
      DO 2222 NROT=1,NRSYM
      NSYMSL(NROT)=0
      DO I=1,100
        NSYMSR(NROT,I)=0
      END DO
C**SKIP IF FINISHED
      IF(INDEX.GT.MVSYM)GO TO 2222
C**SKIP UNTIL FIND NEXT SYMMETRY THAT IS REQUIRED
      IF(NROT.NE.MWSYM(INDEX))GO TO 2222
      INDEX=INDEX+1
      DO 8888 KROT=KSTART,KEND
      KROTL=KROT
      KROTTL=KROTL
1     IF(KROTTL.GT.4)THEN
        KROTTL=KROTTL-4
        GO TO 1
      END IF
C**NROT IS ACTUAL (EVEN-J) SYMMETRY
      NSL=ISYMP(ISYMD(KROTTL,1),NROT)
C**NSL IS ACTUAL SYMMETRY
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
      KAL=KROTL/2
      IF(MOD(KAL,2).NE.0)NSL=ISYMP(NSL,ISYMT)
C**NSL IS STORED SYMMETRY
      IF(NTOT(NSL).EQ.0)GO TO 8888
      IF(KSTEP.EQ.1.AND.IABC.EQ.3)GO TO 8888
      NSYMSL(NROT)=NSYMSL(NROT)+1
      MSYMSL(NROT,NSYMSL(NROT))=NSL
C**IABC=1,2,3
      IF(KSTEP.NE.1.AND.IABC.LE.3)THEN
        NSR=NSL
        NSYMSR(NROT,NSYMSL(NROT))=NSYMSR(NROT,NSYMSL(NROT))+1
        MSYMSR(NROT,NSYMSL(NROT),NSYMSR(NROT,NSYMSL(NROT)))=NSR
      END IF
      KINCR=3
C**IABC=4,8
      IF(IABC.EQ.4.OR.IABC.EQ.8)THEN
        KROTR=KROTL+2
        KROTTR=KROTR
2       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 2
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).EQ.0.AND.KROT+KINCR.LE.K21.AND.ISIZER.NE.0)THEN
        NSYMSR(NROT,NSYMSL(NROT))=NSYMSR(NROT,NSYMSL(NROT))+1
        MSYMSR(NROT,NSYMSL(NROT),NSYMSR(NROT,NSYMSL(NROT)))=NSR
        END IF
      END IF
C**IABC=5,7
      IF(IABC.EQ.5.OR.IABC.EQ.7)THEN
        KROTR=KROTL+KINCR
        KROTTR=KROTR
3       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 3
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).EQ.0.AND.KROT+KINCR.LE.K21.AND.ISIZER.NE.0)THEN
        NSYMSR(NROT,NSYMSL(NROT))=NSYMSR(NROT,NSYMSL(NROT))+1
        MSYMSR(NROT,NSYMSL(NROT),NSYMSR(NROT,NSYMSL(NROT)))=NSR
        END IF
      END IF
      KINCR=2
C**IABC=4,8
      IF(IABC.EQ.4.OR.IABC.EQ.8)THEN
        KROTR=KROTL+KINCR
        KROTTR=KROTR
4       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 4
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).EQ.1.AND.KROT+KINCR.LE.K21.AND.ISIZER.NE.0)THEN
        NSYMSR(NROT,NSYMSL(NROT))=NSYMSR(NROT,NSYMSL(NROT))+1
        MSYMSR(NROT,NSYMSL(NROT),NSYMSR(NROT,NSYMSL(NROT)))=NSR
        END IF
      END IF
C**IABC=5,7
      IF(IABC.EQ.5.OR.IABC.EQ.7)THEN
        KROTR=KROTL+1
        KROTTR=KROTR
5       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 5
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).EQ.1.AND.KROT+KINCR.LE.K21.AND.ISIZER.NE.0)THEN
        NSYMSR(NROT,NSYMSL(NROT))=NSYMSR(NROT,NSYMSL(NROT))+1
        MSYMSR(NROT,NSYMSL(NROT),NSYMSR(NROT,NSYMSL(NROT)))=NSR
        END IF
      END IF
      KINCR=4
C**IABC=1,2
      IF(IABC.EQ.1.OR.IABC.EQ.2)THEN
        KROTR=KROTL+KINCR
        NSR=NSL
        ISIZER=NTOT(NSR)
        IF(KROT+KINCR.LE.K21.AND.ISIZER.NE.0)THEN
        NSYMSR(NROT,NSYMSL(NROT))=NSYMSR(NROT,NSYMSL(NROT))+1
        MSYMSR(NROT,NSYMSL(NROT),NSYMSR(NROT,NSYMSL(NROT)))=NSR
        END IF
      END IF
      KINCR=1
C**IABC=6,9
      IF(IABC.EQ.6.OR.IABC.EQ.9)THEN
        KROTR=KROTL+KINCR
        KROTTR=KROTR
6       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 6
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).EQ.0.AND.ISIZER.NE.0)THEN
        NSYMSR(NROT,NSYMSL(NROT))=NSYMSR(NROT,NSYMSL(NROT))+1
        MSYMSR(NROT,NSYMSL(NROT),NSYMSR(NROT,NSYMSL(NROT)))=NSR
        END IF
      END IF
      KINCR=5
C**IABC=6
      IF(IABC.EQ.6)THEN
        KROTR=KROTL+KINCR
        JROTR=JROTL+KINCR
        IF(MOD(KROT,2).EQ.0.AND.ISIZER.NE.0.AND.KROT+KINCR.LE.K21)THEN
        NSYMSR(NROT,NSYMSL(NROT))=NSYMSR(NROT,NSYMSL(NROT))+1
        MSYMSR(NROT,NSYMSL(NROT),NSYMSR(NROT,NSYMSL(NROT)))=NSR
        END IF
      END IF
      KINCR=4
C**IABC=6
      IF(IABC.EQ.6)THEN
        KROTR=KROTL+3
        KROTTR=KROTR
7       IF(KROTTR.GT.4)THEN
          KROTTR=KROTTR-4
          GO TO 7
        END IF
        NSR=ISYMP(ISYMD(KROTTR,1),NROT)
C**ODD KA HAVE DIFFERENT TORSIONAL SYMMETRIES IF RPH
        KAR=KROTTR/2
        IF(MOD(KAR,2).NE.0)NSR=ISYMP(NSR,ISYMT)
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).NE.0.AND.KROT+KINCR.LE.K21.AND.ISIZER.NE.0)THEN
        NSYMSR(NROT,NSYMSL(NROT))=NSYMSR(NROT,NSYMSL(NROT))+1
        MSYMSR(NROT,NSYMSL(NROT),NSYMSR(NROT,NSYMSL(NROT)))=NSR
        END IF
      END IF
8888  CONTINUE
2222  CONTINUE
      WRITE(IOUT,*)'IABC = ',IABC
      DO 1111 NROT=1,NRSYM
      WRITE(IOUT,*)'LEFT-HAND SIDE'
      N=NSYMSL(NROT)
      WRITE(IOUT,*)N
      WRITE(IOUT,*)'RIGHT-HAND SIDE'
      DO I=1,N
        M=NSYMSR(NROT,I)
        WRITE(IOUT,*)I,': ',M
        WRITE(IOUT,*)'LEFT: ',MSYMSL(NROT,I),' RIGHT: ',
     1  (MSYMSR(NROT,I,K),K=1,M)
      END DO
1111  CONTINUE
      RETURN
      END
C**************************************************************
C**************************************************************
C**TEMPORARY
      SUBROUTINE PRXA(XA,ISIZE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(1)
      COMMON/FILASS/IOUT
      K=0
      DO J=1,ISIZE
        WRITE(IOUT,*)(XA(K+I),I=1,J)
        K=K+J
      END DO
      RETURN
      END
C**TEMPORARY
C**************************************************************
      SUBROUTINE ONEMOD(MODINT,NMODE)
      DIMENSION MODINT(NMODE)
      DO K=1,NMODE
        MODINT(K)=1
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETIND(INDK,INDL,INDN,INDM,NMODE,NTOT2,NTOT3,NTOT4,
     1NTOT5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      DIMENSION INDK(NMODE),INDL(NMODE),INDN(NMODE),INDM(NMODE)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM
      COMMON/FILASS/IOUT
      DO I=1,NMODE
        IF(ICOUPL.GT.1.OR.MOLPRO.GT.4)
     1  INDK(I)=(I-2)*(I-1)/2
        IF(ICOUPL.GT.2.OR.MOLPRO.GT.6)
     1  INDL(I)=(I-3)*(I-2)*(I-1)/6
        IF(ICOUPL.GT.3.OR.MOLPRO.GT.8)
     1  INDN(I)=(I-4)*(I-3)*(I-2)*(I-1)/24
        IF(ICOUPL.GT.4.OR.MOLPRO.GT.10)
     1  INDM(I)=(I-5)*(I-4)*(I-3)*(I-2)*(I-1)/120
      END DO
      IF(ICOUPL.GT.1.OR.MOLPRO.GT.4)NTOT2=INDK(NMODE)+NMODE-1
      IF(ICOUPL.GT.2.OR.MOLPRO.GT.6)NTOT3=INDL(NMODE)+INDK(NMODE)
      IF(ICOUPL.GT.3.OR.MOLPRO.GT.8)NTOT4=INDN(NMODE)+INDL(NMODE)
      IF(ICOUPL.GT.4.OR.MOLPRO.GT.10)NTOT5=INDM(NMODE)+INDN(NMODE)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE FITQ1(XK,EVAL,SOL,WRK,N,XQ,XV,M,QM,NMODE,MODE,INDEX,
     1XTANPM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XK(N,N),EVAL(N,1),SOL(N,1),WRK(4*N),XQ(M),XV(M)
      DIMENSION QM(NMODE),XTANPM(NMODE)
      COMMON/XACC/XACC
      COMMON/NFIT/NFIT
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/FILASS/IOUT,INP
      EXTERNAL FSFUN1,FPRINT
100   FORMAT(/,1X,'DGESV IFAIL = ',I3)
101   FORMAT(D20.10,5X,I3)
102   FORMAT(5X,'SHOULD NOT OCCUR')
103   FORMAT(1X,'SUMSQ = ',D20.10,/)
104   FORMAT(/,1X,'E04FBF IFAIL = ',I3)
105   FORMAT(/,7X,'COEFFICIENT',7X,'POWER',/)
106   FORMAT(/,1X,'MODE: ',I3,/)
107   FORMAT(/,1X,'VALUE OF Q AT MINIMUM: ',D20.10,/)
108   FORMAT(/,1X,'FITTING OF POTENTIAL FOR MODE: ',I3)
109   FORMAT(//,4X,'OBSERVED',7X,'CALCULATED',/)
110   FORMAT(1X,F12.8,5X,F12.8)
111   FORMAT(1X,'SYMMETRY ',A2)
112   FORMAT(/,1X,'RE-FIT POTENTIAL FOR MODE: ',I3)
113   FORMAT(//,4X,'LIMIT(-)',7X,'LIMIT(+)',7X,'XTANH PARAM = ',F12.8,
     1/)
      IF(M.LT.4)STOP 'TOO FEW POINTS'
      IF(M.LT.N)RETURN
      WRITE(IOUT,108)MODE
      ICYCL=0
      INDEX=0
      ISTART=(M-N)/2
      NM=N
      MN=1
      IF(MODE.GT.NONLIN)MN=2
      DO I=1,NM
        EVAL(I,1)=XV(ISTART+I)
        DO J=1,NM
          XK(I,J)=XQ(ISTART+I)**(MN*(J-1))
        END DO
      END DO
      IFAIL=1
      CALL DGESV(NM,1,XK,N,SOL,EVAL,N,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
        STOP 'ERROR IN DGESV'
      END IF
      DO I=1,NM
        SOL(I,1)=EVAL(I,1)
      END DO
      WRITE(IOUT,105)
      DO J=1,NM
        WRITE(IOUT,101)SOL(J,1),MN*(J-1)
      END DO
      CALL FLUSH(IOUT)
C**************************************************************
      IFAIL=1
      MFIT=1
      NFIT=N
      XACC=1.D-15
      FTOL=XACC
      XTOL=1.D-10
C**TEMPORARY
C     STEP=10.D0
      STEP=1.D0
C**TEMPORARY
      MAXCAL=1000
      KPRINT=1000
      IWXY=2*MFIT*(MFIT+MFIT)+2*MFIT+5*MFIT
      IF(IWXY.GT.3*N)THEN
        WRITE(IOUT,102)
        STOP
      END IF
C**ANSWER NEAR ZERO, SO INITIAL GUESS SET TO ZERO
      QM(MODE)=0
      CALL E04FBF(MFIT,MFIT,QM(MODE),SOL,SUMSQ,FTOL,XTOL,STEP,WRK,3*N,
     1FSFUN1,FPRINT,KPRINT,MAXCAL,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,104)IFAIL
        STOP 'ERROR IN E04FBF'
      END IF
      WRITE(IOUT,107)QM(MODE)
C     WRITE(IOUT,103)SUMSQ
C**************************************************************
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE FSFUN1(M,N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      DIMENSION X(N),F(M)
      COMMON/NFIT/NFIT
      IF(MODD.LE.NONLIN)THEN
        F(1)=F(2)
        DO I=3,NFIT
          K=I-1
          J=K-1
          F(1)=F(1)+F(I)*K*X(1)**J
        END DO
      ELSE
        F(1)=0
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE FPRINT(M,N,X,F,S,NCALL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**MONITOR PROGRESS OF E04FBF IF REQUIRED
      DIMENSION X(N),F(M)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETMIN(QM,NMODE,NATOM,XX,X0,XM,XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QM(NMODE),XX(NATOM,3),X0(NATOM,3),XM(NATOM)
      DIMENSION XL(NATOM,NMODE,3)
      COMMON/FILASS/IOUT,INP
100   FORMAT(50(1H#))
      WRITE(IOUT,100)
C**CARTESIAN COORDINATES
      DO I=1,NATOM
        DO J=1,3
          XX(I,J)=X0(I,J)
          DO K=1,NMODE
            XX(I,J)=XX(I,J)+XL(I,K,J)*QM(K)/SQRT(XM(I))
          END DO
        END DO
      END DO
      DO I=1,NATOM
        DO J=1,3
          X0(I,J)=XX(I,J)
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETV(M,XQ,XV,MODE,NMODE,NATOM,QQ,RR,XX,X0,XL,XM,
     1NPOT,IPOT,JPOT,CPOT,NP1,CP1,IP1,VP1,DP1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      LOGICAL ABINIT
      DIMENSION NP1(1),CP1(MMX1,1),IP1(MMX1,1),VP1(MMX1,1),DP1(MMX1,1)
      DIMENSION DUMMY(NMODE)
      DIMENSION XQ(M),XV(M)
      DIMENSION RR(NATOM,NATOM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/ABINIT/ABINIT
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      DO J=1,M
C**ENERGY THIS POINT
        DO K=1,NMODE
          QQ(K)=0
        END DO
        QQ(MODE)=XQ(J)
        IF(IWHICH.GT.0)THEN
          DO I=1,NATOM
            DO K=1,3
              XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQ(MODE)/
     1        SQRT(XM(I))
            END DO
          END DO
          IF(ABINIT)THEN
            CALL GETAPT(VPT,NATOM,XX,RR,QQ,NMODE,0)
          ELSE
            CALL GETPOT(VPT,NATOM,XX,RR)
          END IF
        ELSE
          IF(IWHICH.EQ.0)THEN
            CALL GETPT1(VPT,NPOT,IPOT,JPOT,CPOT,0,QQ)
          ELSE
            IF(MOLINC.GT.0)
     1      CALL GETPQT(VPT,NMODE,QQ,XTANPM,NP1,CP1,IP1,NMODE,MMX1)
            IF(MOLINC.EQ.0)
     1      CALL GETQPT(VPT,NMODE,QQ,XTANPM)
            IF(MOLINC.LT.0)
     1      CALL GETQP1(VPT,DUMMY,NMODE,QQ,NP1,CP1,VP1,DP1,NMODE,MMX1)
          END IF
        END IF
        XV(J)=VPT
C**ENERGY THIS POINT
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETNOR(NATOM,NMODE,XM,X0,OMEGA,XL,XX,RR,
     1IPOT,JPOT,CPOT,NPOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      DIMENSION XM(NATOM),X0(NATOM,3),XL(NATOM,NMODE,3),XX(NATOM,3)
      DIMENSION RR(NATOM,NATOM),OMEGA(NMODE)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
220   FORMAT(/,1X,'MODE',I4,'  OMEGA = ',F20.12,/)
245   FORMAT(//,1X,'DISPLACEMENTS OF NORMAL MODES, L(alpha i,k)')
250   FORMAT(/,1X,'MODE(k) = ',I4)
255   FORMAT(1X,'ATOM(i) ',I2,':  x =',F12.6,'  y =',F12.6,'  z =',
     1F12.6)
      CALL NORMAL(NATOM,NMODE,XM,X0,OMEGA,XL,XX,RR)
C*************************************************************
C*************************************************************
C**CALL USER-SUPPLIED ROUTINE TO GET 'CLEAN' VECTORS IF E OR F
      CALL ROTATE(NATOM,X0,OMEGA,NMODE,XL,WAVENM,NMODE,1)
      CALL ROTATE(NATOM,X0,OMEGA,NMODE,XL,WAVENM,NMODE,2)
      LINBND=0
      IF(LINEAR)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'NORMALS FOR LINEAR MOLECULE'
C**COUNT NUMBER OF LINEAR ANGLE BENDS IF LINEAR
        DO J=1,NMODE
          IF(J.LT.NMODE)THEN
            IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.1.D-3)
     1      LINBND=LINBND+1
          END IF
        END DO
        WRITE(IOUT,*)'NO. LINEAR ANGLE BENDS = ',LINBND
      END IF
C*************************************************************
C*************************************************************
      IF(IPRINT.GT.0)WRITE(IOUT,245)
      DO J=1,NMODE
        IF(IPRINT.GT.0)WRITE(IOUT,250)J
        DO I=1,NATOM
          IF(IPRINT.GT.0)WRITE(IOUT,255)I,(XL(I,J,K),K=1,3)
        END DO
      END DO
      DO K=1,NMODE
        DO I=1,6
          IPOT(K,I)=0
          JPOT(K,I)=0
        END DO
        IPOT(K,1)=K
        JPOT(K,1)=2
        CPOT(K)=OMEGA(K)*OMEGA(K)/2
        IF(IPRINT.GT.0)WRITE(IOUT,220)IPOT(K,1),OMEGA(K)*WAVENM
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE BLIP(NATOM,NMODE,NBF,MBF,NVF,OMEGA,XL,MODINT,MEFF,
     1XLS,MAXBAS,NSMODE,MAXJ,ICI,IREACT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NBF(NMODE),MBF(NMODE),NVF(NMODE),MODINT(NMODE)
      DIMENSION XLS(NATOM,3),MAXBAS(NMODE,MAXJ),MAXS(4)
      DIMENSION MEFF(NMODE,2)
C**EQUILIBRIUM VALUES ONLY
      DIMENSION OMEGA(NMODE),XL(NATOM,NMODE,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/CHECK/MCHECK
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/SUPPL/JSYM(10,100),MCONT(2,100)
C************************
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
220   FORMAT(/,1X,'MODE',I4,'  OMEGA = ',F20.12,/)
245   FORMAT(//,1X,'DISPLACEMENTS OF NORMAL MODES, L(alpha i,k)')
250   FORMAT(/,1X,'MODE(k) = ',I4)
255   FORMAT(1X,'ATOM(i) ',I2,':  x =',F12.6,'  y =',F12.6,'  z =',
     1F12.6)
      IF(IREACT.NE.0)THEN
C**SAVE ORIGINAL SYMMETRIES FIRST TIME IN
        DO I=1,NWSYM
          DO J=1,NSYM(I)
            JSYM(I,J)=ISYM(I,J)
          END DO
        END DO
C**SAVE ORIGINAL CONTRACTIONS FIRST TIME IN
        DO I=1,NCONT
          DO J=1,ICONT(I)
            MCONT(I,J)=JCONT(I,J)
          END DO
        END DO
      ELSE
C**RE-LOAD ORIGINAL SYMMETRIES SECOND TIME IN
        DO I=1,NWSYM
          DO J=1,NSYM(I)
            ISYM(I,J)=JSYM(I,J)
          END DO
        END DO
C**RE-LOAD ORIGINAL CONTRACTIONS SECOND TIME IN
        DO I=1,NCONT
          DO J=1,ICONT(I)
            JCONT(I,J)=MCONT(I,J)
          END DO
        END DO
      END IF
C**NUMBER OF NON-LINEAR MODES = 3N-5 - 2*LINBND
      NONLIN=NMODE-2*LINBND
C**REPEAT OPERATION UNTIL ALL LINEAR BENDS MOVED
1000  CONTINUE
C**IND SET TO ZERO UNTIL SOMETHING MOVED
      IND=0
      DO N=1,NMODE-1
        IF(IND.NE.0)GO TO 1001
C**SEE IF ALL BEEN MOVED ALREADY
        IF(N.LE.NONLIN)THEN
C**IF NOT, FIND LINEAR PAIR
          IF(DABS(OMEGA(N)-OMEGA(N+1))*WAVENM.LT.1.D-3)THEN
C**SKIP AFTER MOVING 2 SETS FROM POSITION N
            IND=1
            DO IDUM=1,2
C**SAVE CURRENT LINEAR PARAMETERS
C             IF(IDUM.EQ.1)THEN
                NBFS=NBF(N)
                MBFS=MBF(N)
                NVFS=NVF(N)
                MEFFS1=MEFF(N,1)
                MEFFS2=MEFF(N,2)
                MODNTS=MODINT(N)
                IF(ICI.LT.0)THEN
                  DO I=1,MAXJ
                    MAXS(I)=MAXBAS(N,I)
                  END DO
                END IF
C             END IF
C************************************
              OMEGAS=OMEGA(N)
              DO I=1,NATOM
                DO K=1,3
                  XLS(I,K)=XL(I,N,K)
C**TEMPORARY-LIN
      IF(IDUM.EQ.1)THEN
        IF(I.EQ.1.AND.DABS(XLS(I,K)).GT.1.D-3)THEN
          IF(XLS(I,K).GT.0)ISG=1
          IF(XLS(I,K).LT.0)ISG=-1
        END IF
      ELSE
        IF(I.EQ.1.AND.DABS(XLS(I,K)).GT.1.D-3)THEN
          IF(XLS(I,K).GT.0)JSG=1
          IF(XLS(I,K).LT.0)JSG=-1
        END IF
        IF(ISG.EQ.JSG)THEN
          IJSG=1
        ELSE
          IJSG=-1
        END IF
      END IF
C**TEMPORARY-LIN
                END DO
              END DO
C**TEMPORARY-LIN
      IF(IDUM.EQ.2)THEN
        DO I=1,NATOM
          DO K=1,3
            XLS(I,K)=IJSG*XLS(I,K)
          END DO
        END DO
      END IF
C**TEMPORARY-LIN
C**MOVE REMAINDER DOWN
              DO N0=N+1,NMODE
C**IF IREACT = 0, EVERYTHING EXCEPT OMEGA AND XL ALREADY MOVED
C               IF(IREACT.NE.0.AND.IDUM.EQ.1.AND.N0.LE.NAMODE)THEN
                  NBF(N0-1)=NBF(N0)
                  MBF(N0-1)=MBF(N0)
                  NVF(N0-1)=NVF(N0)
                  MEFF(N0-1,1)=MEFF(N0,1)
                  MEFF(N0-1,2)=MEFF(N0,2)
                  MODINT(N0-1)=MODINT(N0)
                  IF(ICI.LT.0)THEN
                    DO I=1,MAXJ
                      MAXBAS(N0-1,I)=MAXBAS(N0,I)
                    END DO
                  END IF
C               END IF
C************************************
                OMEGA(N0-1)=OMEGA(N0)
                DO I=1,NATOM
                  DO K=1,3
                    XL(I,N0-1,K)=XL(I,N0,K)
                  END DO
                END DO
              END DO
C**MOVE LINEAR PARAMETERS TO TOP
C             IF(IREACT.NE.0.AND.IDUM.EQ.1)THEN
C               NBF(NAMODE)=NBFS
C               MBF(NAMODE)=MBFS
C               NVF(NAMODE)=NVFS
C               MEFF(NAMODE,1)=MEFFS1
C               MEFF(NAMODE,2)=MEFFS2
C               MODINT(NAMODE)=MODNTS
C               IF(ICI.LT.0)THEN
C                 DO I=1,MAXJ
C                   MAXBAS(NAMODE,I)=MAXS(I)
C                 END DO
C               END IF
C             END IF
              NBF(NMODE)=NBFS
              MBF(NMODE)=MBFS
              NVF(NMODE)=NVFS
              MEFF(NMODE,1)=MEFFS1
              MEFF(NMODE,2)=MEFFS2
              MODINT(NMODE)=MODNTS
              IF(ICI.LT.0)THEN
                DO I=1,MAXJ
                  MAXBAS(NMODE,I)=MAXS(I)
                END DO
              END IF
C************************************
              OMEGA(NMODE)=OMEGAS
              DO I=1,NATOM
                DO K=1,3
                  XL(I,NMODE,K)=XLS(I,K)
                END DO
              END DO
C**SYMMETRY DEFINED BY EQUILIBRIUM VECTORS (NMODE OF THEM)
              IR=0
              JR=0
              DO I=1,NWSYM
                DO J=1,NSYM(I)
                  IF(N.EQ.ISYM(I,J))THEN
C**FOUND LINEAR SYMMETRY (ALWAYS N)
                    IR=I
                    JR=J
                  END IF
                END DO
              END DO
C**MOVE REMAIINDER DOWN ONE PLACE
              DO N0=N+1,NMODE
                DO I=1,NWSYM
                  DO J=1,NSYM(I)
C**MOVE IT INTO TOP (NMODE) POSITION
                    IF(N0.EQ.ISYM(I,J))THEN
                      ISYM(I,J)=ISYM(I,J)-1
                    END IF
                  END DO
                END DO
              END DO
              IF(IR.NE.0.AND.JR.NE.0)ISYM(IR,JR)=NMODE
C**CONTRACTIONS DEFINED BY EQUILIBRIUM VECTORS (NMODE OF THEM)
              IR=0
              JR=0
              DO I=1,NCONT
                DO J=1,ICONT(I)
                  IF(N.EQ.JCONT(I,J))THEN
C**FOUND LINEAR SYMMETRY (ALWAYS N)
                    IR=I
                    JR=J
                  END IF
                END DO
              END DO
C**MOVE REMAIINDER DOWN ONE PLACE
              DO N0=N+1,NMODE
                DO I=1,NCONT
                  DO J=1,ICONT(I)
C**MOVE IT INTO TOP (NMODE) POSITION
                    IF(N0.EQ.JCONT(I,J))THEN
                      JCONT(I,J)=JCONT(I,J)-1
                    END IF
                  END DO
                END DO
              END DO
              IF(IR.NE.0.AND.JR.NE.0)JCONT(IR,JR)=NMODE
            END DO
          END IF
        END IF
1001    CONTINUE
      END DO
C**REPEAT FOR NEXT LINEAR BEND
      IF(IND.NE.0)GO TO 1000
C**************************************************************
C**NOW SORT INTO X,Y PAIRS
      DO N=NONLIN+1,NMODE-1,2
C**INTERCHANGE IF FIRST NOT 'X' (FIRST HALF)
        IF(XL(1,N,2).NE.0)THEN
C**SAVE SECOND SET LINEAR PARAMETERS
          OMEGAS=OMEGA(N+1)
          DO I=1,NATOM
            DO K=1,3
              XLS(I,K)=XL(I,N+1,K)
            END DO
          END DO
C**FIND SECOND SYMMETRY
          IRS=0
          JRS=0
          DO I=1,NWSYM
            DO J=1,NSYM(I)
              IF(N+1.EQ.ISYM(I,J))THEN
                IRS=I
                JRS=J
              END IF
            END DO
          END DO
C**FIND FIRST SYMMETRY
          IR=0
          JR=0
          DO I=1,NWSYM
            DO J=1,NSYM(I)
              IF(N.EQ.ISYM(I,J))THEN
                IR=I
                JR=J
              END IF
            END DO
          END DO
C**IF SECOND EXISTS, SET IT TO FIRST
          IF(IRS.NE.0.AND.JRS.NE.0)ISYM(IRS,JRS)=N
C**IF FIRST EXISTS, SET IT TO SECOND
          IF(IR.NE.0.AND.JR.NE.0)ISYM(IR,JR)=N+1
C**FIND SECOND CONTRACTION
          IRS=0
          JRS=0
          DO I=1,NCONT
            DO J=1,ICONT(I)
              IF(N+1.EQ.JCONT(I,J))THEN
                IRS=I
                JRS=J
              END IF
            END DO
          END DO
C**FIND FIRST CONTRACTION
          IR=0
          JR=0
          DO I=1,NCONT
            DO J=1,ICONT(I)
              IF(N.EQ.JCONT(I,J))THEN
                IR=I
                JR=J
              END IF
            END DO
          END DO
C**IF SECOND EXISTS, SET IT TO FIRST
          IF(IRS.NE.0.AND.JRS.NE.0)JCONT(IRS,JRS)=N
C**IF FIRST EXISTS, SET IT TO SECOND
          IF(IR.NE.0.AND.JR.NE.0)JCONT(IR,JR)=N+1
C**MOVE FIRST SET INTO SECOND SET
          OMEGA(N+1)=OMEGA(N)
          DO I=1,NATOM
            DO K=1,3
              XL(I,N+1,K)=XL(I,N,K)
            END DO
          END DO
C**MOVE SECOND SET INTO FIRST SET
          OMEGA(N)=OMEGAS
          DO I=1,NATOM
            DO K=1,3
              XL(I,N,K)=XLS(I,K)
            END DO
          END DO
        END IF
      END DO
C**************************************************************
C**NOW PUT ALL X TOGETHER
      DO IDUM=1,LINBND-1
        N0=NONLIN+1+IDUM
        N=NONLIN+1+2*IDUM
C**SAVE X SET LINEAR PARAMETERS
        OMEGAS=OMEGA(N)
        DO I=1,NATOM
          DO K=1,3
            XLS(I,K)=XL(I,N,K)
          END DO
        END DO
C**FIND X SYMMETRY
        IRS=0
        JRS=0
        DO I=1,NWSYM
          DO J=1,NSYM(I)
            IF(N.EQ.ISYM(I,J))THEN
              IRS=I
              JRS=J
            END IF
          END DO
        END DO
C**FIND Y SYMMETRY
        IR=0
        JR=0
        DO I=1,NWSYM
          DO J=1,NSYM(I)
            IF(N0.EQ.ISYM(I,J))THEN
              IR=I
              JR=J
            END IF
          END DO
        END DO
C**IF X EXISTS, SET IT TO Y
        IF(IRS.NE.0.AND.JRS.NE.0)ISYM(IRS,JRS)=N0
C**IF Y EXISTS, SET IT TO X
        IF(IR.NE.0.AND.JR.NE.0)ISYM(IR,JR)=N
C**FIND X CONTRACTION
        IRS=0
        JRS=0
        DO I=1,NCONT
          DO J=1,ICONT(I)
            IF(N.EQ.JCONT(I,J))THEN
              IRS=I
              JRS=J
            END IF
          END DO
        END DO
C**FIND Y CONTRACTION
        IR=0
        JR=0
        DO I=1,NCONT
          DO J=1,ICONT(I)
            IF(N0.EQ.JCONT(I,J))THEN
              IR=I
              JR=J
            END IF
          END DO
        END DO
C**IF X EXISTS, SET IT TO Y
        IF(IRS.NE.0.AND.JRS.NE.0)JCONT(IRS,JRS)=N0
C**IF Y EXISTS, SET IT TO X
        IF(IR.NE.0.AND.JR.NE.0)JCONT(IR,JR)=N
C**MOVE Y SET INTO X SET
        OMEGA(N)=OMEGA(N0)
        DO I=1,NATOM
          DO K=1,3
            XL(I,N,K)=XL(I,N0,K)
          END DO
        END DO
C**MOVE X SET INTO Y SET
        OMEGA(N0)=OMEGAS
        DO I=1,NATOM
          DO K=1,3
            XL(I,N0,K)=XLS(I,K)
          END DO
        END DO
      END DO
C**TEMPORARY
C     GO TO 9999
C**TEMPORARY
C**************************************************************
C**NOW ORDER Y ACCORDING TO X; IDUM POINTS AT X
      DO IDUM=NONLIN+1,NONLIN+LINBND
C**FIND LINEAR PAIR; N0 POINTS AT MATCHING Y
        DO N0=NONLIN+LINBND+1,NONLIN+LINBND+LINBND
          IF(DABS(OMEGA(IDUM)-OMEGA(N0))*WAVENM.LT.1.D-3)THEN
C**SAVE Y SET LINEAR PARAMETERS
            OMEGAS=OMEGA(N0)
            DO I=1,NATOM
              DO K=1,3
                XLS(I,K)=XL(I,N0,K)
              END DO
            END DO
C**FIND CORRECT Y POSITION
            N=LINBND+IDUM
C**INTERCHANGE Y POSITIONS
            OMEGA(N0)=OMEGA(N)
            DO I=1,NATOM
              DO K=1,3
                XL(I,N0,K)=XL(I,N,K)
              END DO
            END DO
            OMEGA(N)=OMEGAS
            DO I=1,NATOM
              DO K=1,3
                XL(I,N,K)=XLS(I,K)
              END DO
            END DO
          END IF
        END DO
      END DO
C**TEMPORARY
9999  CONTINUE
C**TEMPORARY
      WRITE(IOUT,*)'***********************************************'
      WRITE(IOUT,*)'*********NORMAL MODE ORDERING CHANGED**********'
      WRITE(IOUT,*)'***********************************************'
      DO I=1,NMODE
        WRITE(IOUT,*)'MODE ',I
        WRITE(IOUT,*)'NBF ',NBF(I)
        WRITE(IOUT,*)'MBF ',MBF(I)
        WRITE(IOUT,*)'NVF ',NVF(I)
        WRITE(IOUT,*)'MEFF ',MEFF(I,1)
        WRITE(IOUT,*)'     ',MEFF(I,2)
        WRITE(IOUT,*)'MODINT ',MODINT(I)
        WRITE(IOUT,*)'OMEGA ',OMEGA(I)*WAVENM
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

C     WRITE(IOUT,245)
C     DO J=1,NMODE
C       WRITE(IOUT,250)J
C       DO I=1,NATOM
C         WRITE(IOUT,255)I,(XL(I,J,K),K=1,3)
C       END DO
C     END DO
C     DO K=1,NMODE
C       WRITE(IOUT,220)K,OMEGA(K)*WAVENM
C     END DO
C     DO I=1,NAMODE
C       WRITE(IOUT,*)NBF(I),MBF(I),NVF(I)
C     END DO
C     DO I=1,NWSYM
C       WRITE(IOUT,*)'MODE SYMMETRY ',I,' NO. MODES ',NSYM(I)
C       WRITE(IOUT,*)'MODES ',(ISYM(I,J),J=1,NSYM(I))
C     END DO
C     DO I=1,NCONT
C       WRITE(IOUT,*)'CONTRACTION SCHEME ',I,' NO. MODES ',ICONT(I)
C       WRITE(IOUT,*)'MODES ',(JCONT(I,J),J=1,ICONT(I))
C     END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE INPUT(NMODE,OMEGA,XX,NBF,MBF,NSTAT
     1,ISTAT,NATOM,XM,X0,XL,RR,ICI,MAXBAS,MAXJ,JSTAT,KSTAT,ESTAT,WSTAT,
     2NVF,MVF,NOCON,INORM,MCHECK,NSMODE,NVMODE,INDEX,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      DIMENSION JSTAT(NSTAT,NMODE),KSTAT(NMODE),ESTAT(NSTAT)
      DIMENSION WSTAT(NSTAT),RR(NATOM,NATOM)
      DIMENSION NVF(NSMODE),MVF(NSMODE),MAXBAS(NMODE,MAXJ),NOCON(NSMODE)
      DIMENSION NBF(NSMODE),MBF(NSMODE),OMEGA(NMODE),ISTAT(NSTAT,NMODE)
      DIMENSION XM(NATOM),X0(NATOM,3),XL(NATOM,NMODE,3),XX(NATOM,3)
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CADDR/LIPOT,LJPOT,LCPOT
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM
      COMMON/TYPE/LINEAR
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/NCPOT/NPOT
      COMMON/WHICH/IWHICH
      COMMON/PATH/ISCFCI
      COMMON/ESTATE/IORDER
      COMMON/FILASS/IOUT,INP
      COMMON/MATSIZ/MATSIZ
      COMMON/MOMI/XK(3,3),XMU(3,3)
      DIMENSION YJ(3)
100   FORMAT(100A2)
200   FORMAT(/,1X,'FOR MODE ',I4,/,
     1         1X,'NUMBER PRIMITIVE FUNCTIONS = ',I4,/,
     2         1X,'NUMBER GAUSS POINTS = ',I4,/,
C*****************************************************HEG
C*****************************************************HEG
     3         1X,'NUMBER CONTRACTED FUNCTIONS = ',I4,/,
     4         1X,'NUMBER HEG POINTS = ',I4,/)
C*****************************************************HEG
C*****************************************************HEG
205   FORMAT(1X,'NO. OF GAUSS POINTS CHANGED FROM ',I4,' TO ',I4)
2055  FORMAT(1X,'NO. OF HEG POINTS CHANGED FROM ',I4,' TO ',I4)
210   FORMAT(/,1X,'POTENTIAL TERMS: IPOT, JPOT, CPOT',/)
215   FORMAT(1X,6I2,1X,6I2,1X,D20.12)
220   FORMAT(/,1X,'MODE',I4,'  OMEGA = ',F7.2,/)
225   FORMAT(//,1X,'NORMAL COORDINATE ANALYSIS',/)
230   FORMAT(1X,'ATOM ',I2,':  MASS =',F12.7)
235   FORMAT(//,1X,'EQUILIBRIUM GEOMETRY',/)
240   FORMAT(1X,'ATOM ',I2,':  X0 =',F12.7,'  Y0 =',F12.7,'  Z0 =',
     1F12.7)
245   FORMAT(//,1X,'DISPLACEMENTS OF NORMAL MODES, L(alpha i,k)')
250   FORMAT(/,1X,'MODE(k) = ',I4)
255   FORMAT(1X,'ATOM(i) ',I2,':  x =',F12.6,'  y =',F12.6,'  z =',
     1F12.6)
260   FORMAT(1X,'NO. OF CONTRACTED FUNCTIONS CHANGED FROM ',I4,' TO ',
     1I4)
265   FORMAT(1X,'ATOM ',I2,':  SYMBOL =',A2)
270   FORMAT(//,1X,'NUMBER POTENTIAL TERMS = ',I4,/)
      IF(INDEX.NE.0)GO TO 9999
C**NO. PRIMITIVE QUANTA AND ASSOCIATED INTEGRATION POINTS
      READ(INP,*)
C****************************************************************************
C**NBF(NMODE),MBF(NMODE),NVF(NMODE),MVF(NMODE)
C--------------------------------------------
C**
C**            Input NMODE records of three integers
C**NBF(K):     Number of quanta in harmonic-oscillator basis for mode K
C**MBF(K):     Number of Gauss Hermite integration points for mode K
C**NVF(K):     Number of quanta of contracted numerical functions for mode K
C**MVF(K):     Number of HEG integration points for mode K
C**Note:       For VCI (ICI < 0) NVF(K) will be overwritten by -ICI or by the
C**            maximum quantum within MAXBAS, depending on the setting of NMAX
C****************************************************************************
      DO K=1,NMODE
        READ(INP,*)NBF(K),MBF(K),NVF(K),MVF(K)
C**SKIP CONTRACTION?
        NOCON(K)=NVF(K)
        NVF(K)=IABS(NVF(K))
C**NUMBER OF FUNCTIONS ONE GREATER THAN NUMBER OF QUANTA
        NBF(K)=NBF(K)+1
        NVF(K)=NVF(K)+1
C**MAKE SURE 'EVEN' NUMBER QUANTA FOR TAU
CC      IF(K.EQ.IABS(IREACT))THEN
CC        IF(MOD(NBF(K),2).EQ.0)NBF(K)=NBF(K)-1
CC        IF(MOD(NVF(K),2).EQ.0)NVF(K)=NVF(K)-1
CC      END IF
C**MAKE SURE 'ODD' NUMBER QUANTA FOR TAU
        IF(K.EQ.IABS(IREACT))THEN
          IF(MOD(NBF(K),2).NE.0)NBF(K)=NBF(K)+1
          IF(MOD(NVF(K),2).NE.0)NVF(K)=NVF(K)+1
        END IF
        IF(ICI.LT.0.AND.ISCFCI.GT.0)THEN
          IF(MAXJ.EQ.0)THEN
            IF(-ICI.GT.NBF(K))ICI=-NBF(K)
            NVF(K)=-ICI
          ELSE
            CALL HEG2(MAXBAS,NMODE,MAXJ,NVF,K)
          END IF
        END IF
        IF(NVF(K).LE.0)NVF(K)=NBF(K)
        IF(NVF(K).GT.NBF(K))THEN
          WRITE(IOUT,260)NBF(K),NVF(K)
          NBF(K)=NVF(K)
        END IF
C**KINETIC ENERGY INTEGRALS (AT LEAST) EXACT IF MBF.GE.NBF
        IF(MBF(K).LT.NBF(K))THEN
          WRITE(IOUT,205)MBF(K),NBF(K)
          MBF(K)=NBF(K)
        END IF
C**KINETIC ENERGY INTEGRALS (AT LEAST) EXACT IF MVF.GE.NVF
        IF(MVF(K).LT.NVF(K))THEN
          WRITE(IOUT,2055)MVF(K),NVF(K)
          MVF(K)=NVF(K)
        END IF
CC      IF(K.NE.IABS(IREACT))THEN
CC        IF(NVF(K).LT.NBF(K))THEN
CC          INCR=MBF(K)-NBF(K)
CC          IF(NBF(K).LT.NVF(K)+INCR)NBF(K)=NVF(K)+INCR
CC          IF(MBF(K).LT.NBF(K)+INCR)MBF(K)=NBF(K)+INCR
CC        END IF
CC      END IF
C**MAKE SURE 'ODD' NUMBER QUANTA FOR TAU
        IF(K.EQ.IABS(IREACT))THEN
          IF(MOD(NVF(K),2).NE.0)NVF(K)=NVF(K)+1
        END IF
C       IF(MOLPRO.EQ.0)THEN
C**MAKE NO, POINTS EVEN FOR STANDARD RUN
          IF(MOD(MBF(K),2).NE.0)MBF(K)=MBF(K)+1
          IF(MOD(MVF(K),2).NE.0)MVF(K)=MVF(K)+1
C       ELSE
C**MAKE NO, POINTS ODD FOR MOLPRO RUN
C         IF(MOD(MBF(K),2).EQ.0)MBF(K)=MBF(K)+1
C       END IF
        WRITE(IOUT,200)K,NBF(K),MBF(K),NVF(K),MVF(K)
C**SKIP CONTRACTION?
        IF(NOCON(K).GE.0)NOCON(K)=NVF(K)
        IF(NOCON(K).LT.0)NOCON(K)=-NVF(K)
      END DO
C**********************************************************
C**********************************************************
      READ(INP,*)
C****************************************************************************
C**SYMBOL(NATOM)
C--------------
C**
C**(SYMBOL(I),I=1,NATOM)
C**
C**SYMBOL:  Nuclear Atomic Symbol (H for hydrogen, etc.) in A2 FORMAT
C****************************************************************************
      READ(INP,100)(SYMBOL(I),I=1,NATOM)
      DO I=1,NATOM
        WRITE(IOUT,265)I,SYMBOL(I)
      END DO
C**********************************************************
C**********************************************************
      READ(INP,*)
      IF(IORDER.EQ.0)THEN
C****************************************************************************
C**NDEF
C------
C**
C**         SCF state definition
C**
        READ(INP,*)NDEF
C**
C**Read data for this input ONLY if NSTAT > 0
C**         NSTAT +ve and NDEF -ve:
C**         Input NSTAT records ISTAT(NSTAT,NMODE) corresponding to NSTAT
C**         specific SCF state definitions (see Jelski et al in
C**         J. Comput. Chem., 17, 1645 (1996)).
C**
C**         NSTAT +ve and NDEF=0 or NDEF +ve:
C**         Input NSTAT records N1 N2 where N1 are the number of quanta for
C**         mode N2
C**         All remaining modes will take on NDEF quanta (see Jelski et al)
C**
C**         NSTAT -ve:
C**         Input NOTHING. In this case NDEF is a dummy parameter.
C**         It is recommended that NSTAT is set -ve for all VCI calculations
C**         (see ICI -ve above)
C**
C**         NB. For linear molecules, ....nb1,nb2,...,nl1,nl2,...
C****************************************************************************
C***************************************
        READ(INP,*)
C**ORIGINAL JELSKI
        IF(NDEF.LT.0)THEN
          DO I=1,NSTAT
            READ(INP,*)(ISTAT(I,K),K=1,NMODE)
          END DO
        ELSE
          DO K=1,NMODE
            ISTAT(1,K)=NDEF
          END DO
          DO I=2,NSTAT
            READ(INP,*)N1,N2
            DO K=1,NMODE
              IF(K.NE.N2)THEN
                ISTAT(I,K)=NDEF
              ELSE
                ISTAT(I,K)=N1
              END IF
            END DO
          END DO
        END IF
      END IF
C**********************************************************
C**********************************************************
      READ(INP,*)
C****************************************************************************
C**NPOT
C-----
C**
C**Read data for this input ONLY for Normal Coordinate Potential (IWHICH=0)
C**NPOT:    Number of user-supplied Normal coordinate force field terms
C**         IPOT,JPOT,CPOT (below).
C****************************************************************************
      IF(IWHICH.EQ.0)READ(INP,*)NPOT
      WRITE(IOUT,270)NPOT
      KIPOT=NPOT*6
      KJPOT=NPOT*6
      KCPOT=NPOT
      IF(NPOT.NE.0)
     1CALL MEMO(3,LIPOT,KIPOT,LJPOT,KJPOT,LCPOT,KCPOT,0,0,0,0)
C**********************************************************
C**********************************************************
      READ(INP,*)
C****************************************************************************
C**IPOT(NPOT),JPOT(NPOT),CPOT(NPOT)
C----------------------------------
C**
C**Read data for this input EITHER for Normal Coordinate Potential (IWHICH=0)
C**OR for definition of NMODE values of user-supplied values of omega (INORM=0)
C**         Normal coordinate force field, defined by:
C**         Jelski et al in J. Comput. Chem., 17, 1645 (1996).
C**         A maximum of 6 coupled modes is allowed.
C**IPOT:    (first 6 integers): defines the Normal mode(s) involved.
C**JPOT:    (second 6 integers): defines the number of quanta in the
C**         corresponding mode(s)
C**CPOT:    defines the potential contribution F(ijkl..) in hartrees, except
C**         for the harmonic terms (see below)
C**
C**These parameters are input under the following conditions:
C**IWHICH=0: Input NPOT terms IPOT,JPOT,CPOT where the first NMODE terms are the
C**          NMODE harmonic force constants F(ii)
C**IWHICH=1: Input NMODE terms corresponding to the NMODE harmonic force
C**          constants F(ii) if INORM=0
C**          Input NOTHING if INORM=1
C**
C**Note:     F(ii) is omega(i) in cm-1
C****************************************************************************
      IF(IWHICH.EQ.0.OR.INORM.EQ.0)THEN
        WRITE(IOUT,210)
        CALL INPOT(W(LIPOT),W(LJPOT),W(LCPOT),NPOT,OMEGA,NMODE)
      END IF
C**********************************************************
C**********************************************************
C**READ NORMAL COORDINATE DETAILS
      WRITE(IOUT,225)
      READ(INP,*)
C****************************************************************************
C**XM(NATOM)
C-----------
C**
C**(XM(I),I=1,NATOM)
C**XM:      Nuclear masses (U)
C****************************************************************************
      READ(INP,*)(XM(I),I=1,NATOM)
      DO I=1,NATOM
        WRITE(IOUT,230)I,XM(I)
        XM(I)=XM(I)*ELMASS
      END DO
C**EQUILIBRIUM CARTESIAN COORDINATES
      READ(INP,*)
C****************************************************************************
C**X0(NATOM,3)
C-------------
C**
C**(X0(NATOM,I),I=1,3)
C**X0:      Input NATOM records of reference structure x,y,z coordinates
C****************************************************************************
      DO I=1,NATOM
        READ(INP,*)(X0(I,J),J=1,3)
      END DO
      WRITE(IOUT,235)
      DO I=1,NATOM
        WRITE(IOUT,240)I,(X0(I,J),J=1,3)
      END DO
C**NORMAL COORDINATE DISPLACEMENTS
      READ(INP,*)
C****************************************************************************
C**XL(NATOM,NMODE,3)
C-------------------
C**
C**Omit data for this input if evaluating Normal Coordinates (INORM=1)
C**
C**XL(NATOM,NMODE,I),   I=1,3
C**         Mass-weighted Normal Coordinate displacement
C**XL:      Input NMODE groups of NATOM records of x,y,z mass-weighted Normal
C**         Coordinate displacements (INORM=0)
C****************************************************************************
C     IF(INORM.EQ.0.AND.IREACT.EQ.0)THEN
      IF(INORM.EQ.0)THEN
        WRITE(IOUT,245)
        DO J=1,NMODE
          DO I=1,NATOM
            READ(INP,*)(XL(I,J,K),K=1,3)
          END DO
          WRITE(IOUT,250)J
          DO I=1,NATOM
            WRITE(IOUT,255)I,(XL(I,J,K),K=1,3)
          END DO
        END DO
        LINBND=0
        IF(LINEAR)THEN
C**COUNT NUMBER OF LINEAR ANGLE BENDS IF LINEAR
          DO J=1,NMODE
            IF(J.LT.NMODE)THEN
              IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.1.D-3)
     1        LINBND=LINBND+1
            END IF
          END DO
        END IF
C**NUMBER OF NON-LINEAR MODES = 3N-5 - 2*LINBND
        NONLIN=NMODE-2*LINBND
      END IF
C**TEST FOR LINEARITY
      CALL PAXES(XM,X0,XX,RR,XL,NMODE,NATOM,0)
      WRITE(IOUT,*)
      IF(LINEAR)WRITE(IOUT,*)'LINEAR MOLECULE'
      IF(.NOT.LINEAR)WRITE(IOUT,*)'BENT MOLECULE'
C**CHECK C.OF M. AND PRINCIPAL AXES
CC    IF(MCHECK.GT.0.AND.IREACT.EQ.0)THEN
      IF(MCHECK.GT.0)THEN
        CALL CHECKM(XM,X0,XX,RR,NATOM)
CC      IF(ISCFCI.LT.0)STOP 'C.M. CHECK'
        IF(ISCFCI.LT.0.AND.IREACT.EQ.0)STOP 'C.M. CHECK'
      END IF
CC    IF(INORM.EQ.0.AND.MCHECK.GT.0.AND.IREACT.EQ.0)THEN
      IF(INORM.EQ.0.AND.MCHECK.GT.0)THEN
C**TRANSFORM N. COORD. VECTORS IF TRANSFORMED TO PA SYSTEM
C**XMU ACTS AS TRANSFORMATION MATRIX
        DO M=1,NMODE
          DO I=1,NATOM
            DO J=1,3
              YJ(J)=0
              DO K=1,3
                YJ(J)=YJ(J)+XMU(K,J)*XL(I,M,K)
              END DO
            END DO
            DO K=1,3
              XL(I,M,K)=YJ(K)
            END DO
          END DO
        END DO
        WRITE(IOUT,*)
        WRITE(IOUT,*)'TRANSFORMED N. COORD. VECTORS'
        WRITE(IOUT,*)
        WRITE(IOUT,245)
        DO J=1,NMODE
          WRITE(IOUT,250)J
          DO I=1,NATOM
            WRITE(IOUT,255)I,(XL(I,J,K),K=1,3)
          END DO
        END DO
      END IF
      RETURN
C**********************************************************
C**********************************************************
9999  CONTINUE
      IF(IORDER.NE.0)THEN
C**SET UP NSTAT STATES BASED ON ENERGY
        DO K=1,NSMODE
          DO I=1,NSTAT
            ISTAT(I,K)=0
          END DO
        END DO
C**START STACK OFF WITH PURE MODE NV1 (LOWEST ENERGY)
        IF(.NOT.LINEAR)THEN
C**NON-LINEAR
          NMULT=1
          IF(IREACT.EQ.0)NV1=1
          IF(IREACT.GT.0)NV1=NSMODE
        ELSE
C**LINEAR
          NV1=NONLIN+1
          NMULT=2
        END IF
        ESTAT(1)=0
        DO I=2,NSTAT
          J=I-1
          ESTAT(I)=OMEGA(NV1)*J*NMULT
          ISTAT(I,NV1)=J
        END DO
C**ADD REMAINING MODES
        KBK=0
        DO KVK=2,NVMODE
          NV1=NV1+1
          IF(NV1.LE.NVMODE)THEN
            K=NV1
          ELSE
            KBK=KBK+1
            K=KBK
          END IF
          NMULT=1
          IF(K.GT.NONLIN)NMULT=2
C**SAVE ISTAT AND ESTAT IN JSTAT AND WSTAT
          DO I=1,NSTAT
            WSTAT(I)=ESTAT(I)
            DO L=1,NSMODE
              JSTAT(I,L)=ISTAT(I,L)
            END DO
          END DO
C**FOR MODE K, ADD N1 QUANTA TO ORIGINAL STATE I1 IN ISTAT
          MORES=0
          I1=0
1000      CONTINUE
          I1=I1+1
          MOREQ=0
          N1=0
2000      CONTINUE
          N1=N1+1
C**KSTAT CONTAINS NEW STATE, ETEMP NEW ENERGY
          DO L=1,NSMODE
            KSTAT(L)=ISTAT(I1,L)
          END DO
          KSTAT(K)=ISTAT(I1,K)+N1
          ETEMP=0
          DO L=1,NVMODE
            ETEMP=ETEMP+KSTAT(L)*OMEGA(L)*NMULT
          END DO
C**SEE IF NEW STATE REQUIRED
          DO I=1,NSTAT
            IF(WSTAT(I).GT.ETEMP)THEN
C**NEED IT...INSERT IT AND MOVE THE REST UP
C**FIRST MOVE THEM
              I2=NSTAT-I
              DO J=1,I2
                WSTAT(NSTAT+1-J)=WSTAT(NSTAT-J)
                DO L=1,NMODE
                  JSTAT(NSTAT+1-J,L)=JSTAT(NSTAT-J,L)
                END DO
              END DO
C**THEN INSERT
              WSTAT(I)=ETEMP
              DO L=1,NSMODE
                JSTAT(I,L)=KSTAT(L)
              END DO
              GO TO 3000
            ELSE
C**MOREQ SET IF CAN'T INSERT NEW STATE
              IF(I.EQ.NSTAT)THEN
                MOREQ=1
                GO TO 3000
              END IF
            END IF
          END DO
3000      CONTINUE
C**TEST IF POSSIBLY ONE MORE QUANTUM
          IF(MOREQ.EQ.0)GO TO 2000
C**TEST IF POSSIBLY ONE MORE STATE
          IF(I1.EQ.NSTAT)MORES=1
          IF(MORES.EQ.0)GO TO 1000
C**MOVE NEW STATES BACK TO ISTAT AND ESTAT
          DO I=1,NSTAT
            ESTAT(I)=WSTAT(I)
            DO L=1,NSMODE
              ISTAT(I,L)=JSTAT(I,L)
            END DO
          END DO
        END DO
C**IF LINEAR GENERATE NSMODE QUANTA FROM (N) FOR LINEAR BENDS
        IF(LINEAR)THEN
C**CURRENTLY J=0 ONLY
          I=0
          DO IVI=1,NSTAT
            I=I+1
            IF(I.GT.NSTAT)GO TO 9000
            MINV=1000
            DO K=NONLIN+1,NONLIN+LINBND
              IF(ISTAT(I,K).LT.MINV)MINV=ISTAT(I,K)
            END DO
C**MINV IS MAXIMUM VALUE OF L
C**IF ODD, MINIMUM VALUE OF L IS 1
C**IF EVEN, MINIMUM VALUE OF L IS 0
            IF(MOD(MINV,2).EQ.0)THEN
              MINL=0
            ELSE
              MINL=1
            END IF
8000  CONTINUE
            ISTAT(I,NSMODE)=4*MINL
            I=I+1
            IF(I.GT.NSTAT)GO TO 9000
            MINL=MINL+2
            IF(MINL.LE.MINV)THEN
C**NEXT HIGHER L
              IF(I.LT.NSTAT)THEN
C**MOVE REMAINDER UP ONE PLACE
                DO JVJ=I+1,NSTAT
                  J=NSTAT-JVJ+I+1
                  DO K=1,NONLIN+LINBND
                    ISTAT(J,K)=ISTAT(J-1,K)
                  END DO
                END DO
              END IF
C**DUPLICATE
              DO K=1,NONLIN+LINBND
                ISTAT(I,K)=ISTAT(I-1,K)
              END DO
              GO TO 8000
            END IF
          END DO
9000  CONTINUE
        END IF
      ELSE
C**THROW OUT FORBIDDEN STATES IF LINEAR
        NXSTAT=0
        DO I=1,NSTAT
          ISUM=0
          DO K=NONLIN+1,NVMODE
            ISUM=ISUM+ISTAT(I,K)
          END DO
          IF(MOD(ISUM,2).EQ.0)THEN
            NXSTAT=NXSTAT+1
            DO K=1,NSMODE
              ISTAT(NXSTAT,K)=ISTAT(I,K)
            END DO
          END IF
        END DO
        NSTAT=NXSTAT
C**IF LINEAR GENERATE NSMODE QUANTA FROM (N,L) FOR LINEAR BENDS
        IF(LINEAR)THEN
C**CURRENTLY J=0 ONLY
          DO I=1,NSTAT
            DO K=NONLIN+1,NONLIN+LINBND
              ISTAT(I,K)=ISTAT(I,K)-ISTAT(I,K+LINBND)
              ISTAT(I,K)=ISTAT(I,K)/2
            END DO
            ISUM=0
            DO K=NONLIN+LINBND+1,NMODE
C**ERROR TRAPS REQUIRED
              ISUM=ISUM+ISTAT(I,K)
            END DO
C**NOT FOR LINEAR TRIATOMIC
            IF(NSMODE.GT.NONLIN+LINBND)THEN
              ISTAT(I,NSMODE)=2*ISUM
            END IF
          END DO
        END IF
      END IF
      WRITE(IOUT,*)
      IF(IORDER.NE.0)THEN
        WRITE(IOUT,*)NSTAT,' ENERGIES AND COMPUTED SCF STATES'
      ELSE
        WRITE(IOUT,*)NSTAT,' COMPUTED SCF STATES'
      END IF
      WRITE(IOUT,*)
      DO I=1,NSTAT
        IF(IORDER.NE.0)
     1  WRITE(IOUT,*)'STATE: ',I,'    ENERGY: ',ESTAT(I)*WAVENM
        WRITE(IOUT,*)(ISTAT(I,K),K=1,NSMODE)
      END DO
      IF(MATSIZ.NE.0.AND.ISCFCI.EQ.0)STOP 'VSCF MATRIX SIZE'
C**********************************************************
C**********************************************************
C**RESET QUANTA TO FUNCTION POINTERS
      DO I=1,NSTAT
        DO K=1,NSMODE
          ISTAT(I,K)=ISTAT(I,K)+1
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE INPOT(IPOT,JPOT,CPOT,NPOT,OMEGA,NMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT),OMEGA(NMODE)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
215   FORMAT(1X,6I2,1X,6I2,1X,D20.12)
220   FORMAT(/,1X,'MODE',I4,'  OMEGA = ',F7.2,/)
      DO I=1,NPOT
C**MAXIMUM OF 6 MODES COUPLED
        READ(INP,*)(IPOT(I,J),J=1,6),(JPOT(I,J),J=1,6),CPOT(I)
        WRITE(IOUT,215)(IPOT(I,J),J=1,6),(JPOT(I,J),J=1,6),
     1  CPOT(I)
        IF(I.LE.NMODE)THEN
CCCC      IF(IWHICH.EQ.0)OMEGA(I)=SQRT(2*CPOT(I))
CCCC      IF(IWHICH.NE.0)OMEGA(I)=CPOT(I)/WAVENM
          OMEGA(I)=CPOT(I)/WAVENM
CCCC      IF(IWHICH.EQ.0)CPOT(I)=OMEGA(I)*OMEGA(I)/2
          CPOT(I)=OMEGA(I)*OMEGA(I)/2
          WRITE(IOUT,220)IPOT(I,1),OMEGA(I)*WAVENM
        END IF
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE CHECKM(XM,X0,XX,RR,NATOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINEAR
      DIMENSION XM(NATOM),X0(NATOM,3),XX(NATOM,3),RR(NATOM,NATOM)
      DIMENSION XCM(3),WR(3),E(3)
      DIMENSION MROT(3)
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/REACTL/JREACT,IFITRP
      COMMON/FILASS/IOUT,INP
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/TYPE/LINEAR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/MOMI/XK(3,3),XMU(3,3)
100   FORMAT(//,1X,'CENTRE OF MASS CONDITIONS: X,Y,Z ',//,1X,3D20.12,/)
150   FORMAT(//,1X,'PRINCIPAL AXIS CONDITIONS',/)
200   FORMAT(3D20.12)
250   FORMAT(//,1X,'VECTORS',/)
300   FORMAT(//,1X,'C.OF.M COORDINATES: X,Y,Z ',/)
350   FORMAT(//,1X,'ORIGINAL COORDINATES: X,Y,Z ',/)
400   FORMAT(//,1X,'PRINCIPAL AXIS COORDINATES: X,Y,Z ',/)
450   FORMAT(//,1X,'ROTATIONAL CONSTANTS',/)
500   FORMAT(//,1X,'MOMENT OF INERTIA TENSOR',/)
550   FORMAT(//,1X,'CHECK PRINCIPAL AXES COORDINATES')
600   FORMAT(//,1X,'ENERGY (Eh) AT MINIMUM ',F12.6,/)
650   FORMAT(//,1X,'ROTATED PRINCIPAL AXIS COORDINATES: X,Y,Z ',/)
      IF(IRET.EQ.0.OR.JPRINT.GT.0)WRITE(IOUT,550)
C**RETURN WITH ROTATIONAL CONSTANTS
C**ORIGINAL COORDINATES
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,350)
        DO I=1,NATOM
          WRITE(IOUT,200)(X0(I,J),J=1,3)
        END DO
      END IF
      CALL ROTC(X0,XM,NATOM,A,B,C,BBAR,1)
      IF(IRET.EQ.0.OR.JPRINT.GT.0)WRITE(IOUT,450)
      IF(IRET.EQ.0.OR.JPRINT.GT.0)WRITE(IOUT,200)A*WAVENM,B*WAVENM,
     1C*WAVENM
      INDRTG=0
8888  CONTINUE
C**FIND CENTRE OF MASS
      DO I=1,3
        XCM(I)=0
      END DO
      XMASS=0
      DO I=1,NATOM
        XMASS=XMASS+XM(I)
      END DO
      DO I=1,3
        DO J=1,NATOM
          XCM(I)=XCM(I)+XM(J)*X0(J,I)
        END DO
        XCM(I)=XCM(I)/XMASS
      END DO
      DO I=1,3
        DO J=1,NATOM
          XX(J,I)=X0(J,I)-XCM(I)
          IF(DABS(XX(J,I)).LT.1.D-10)XX(J,I)=0.D0
        END DO
      END DO
C**COORDINATES RELATIVE TO C. OF M.
      IF(IRET.EQ.0.OR.JPRINT.GT.0)WRITE(IOUT,300)
      DO I=1,NATOM
        IF(IRET.EQ.0.OR.JPRINT.GT.0)WRITE(IOUT,200)(XX(I,J),J=1,3)
      END DO
      IF(INDRTG.NE.0)GO TO 7777
C**********************************************
      IF(IRET.GT.0)RETURN
C**********************************************
      CALL ROTC(XX,XM,NATOM,A,B,C,BBAR,1)
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,450)
        WRITE(IOUT,200)A*WAVENM,B*WAVENM,C*WAVENM
      END IF
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
            X=X+XX(I,J)*XX(I,J) ! chen
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
C**MOMENT OF INERTIA TENSOR
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,500)
        DO I=1,3
          WRITE(IOUT,200)(XK(J,I),J=1,3)
        END DO
      END IF
C**DIAGONALISE MOMENT OF INERTIA MATRIX
      IF(IRET.EQ.0)IRET=1
      CALL DIAG(XK,XK,3,3,-1,WR,E,WR,3,3,XK,3,3,XK,XK,XK,IDUM,IDUM,
     1IDUM)
C**E ARE IN DESCENDING ORDER
      IF(IRET.GT.0)IRET=0
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
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,450)
        WRITE(IOUT,200)A*WAVENM,B*WAVENM,C*WAVENM
C**WRITE VECTORS
        WRITE(IOUT,250)
        DO I=1,3
          WRITE(IOUT,200)(XK(J,I),J=1,3)
        END DO
      END IF
C**TEMPORARY
C**INTERCHANGE X & Y
C     IF(JPRINT.GT.0.AND.JREACT.NE.0)
C    1WRITE(IOUT,*)'INTERCHANGING X & Y IN CHECKM'
C     DO J=1,3
C       XSAVE=XK(J,1)
C       XK(J,1)=XK(J,2)
C       XK(J,2)=XSAVE
C     END DO
C**TEMPORARY
C**MOVE TO PRINCIPAL AXES
      DO I=1,NATOM
        DO J=1,3
          X0(I,J)=0
          DO K=1,3
            X0(I,J)=X0(I,J)+XK(K,J)*XX(I,K)
          END DO
        END DO
      END DO
C**KEEP COPY
      DO J=1,3
        DO K=1,3
          XMU(K,J)=XK(K,J)
        END DO
      END DO
CCCC  IF(INDRTG.NE.0)GO TO 7777
C**ORIGINAL (PRINCIPAL AXIS) COORDINATES
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,400)
        DO I=1,NATOM
          WRITE(IOUT,200)(X0(I,J),J=1,3)
        END DO
      END IF
C*******************************************
C*******************************************
C**CALL USER-SUPPLIED ROUTINE TO GET 'CLEAN'
C**GEOMETRY IF SYMMETRIC OR SPHERICAL TOP
      CALL RTGEOM(NATOM,X0,RR,E,WAVENM,INDRTG)
      IF(INDRTG.NE.0)GO TO 8888
      GO TO 9999
C*******************************************
C*******************************************
7777  CONTINUE
C**ROTATED (PRINCIPAL AXIS) COORDINATES
      IF(INDRTG.NE.0)THEN
        IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
          WRITE(IOUT,650)
          DO I=1,NATOM
            WRITE(IOUT,200)(X0(I,J),J=1,3)
          END DO
        END IF
      END IF
C****************************************************
C**CHECK RESULTS
C****************************************************
9999  CONTINUE
      DO I=1,3
        XCM(I)=0
      END DO
      DO I=1,3
        DO J=1,NATOM
          XCM(I)=XCM(I)+XM(J)*X0(J,I)
        END DO
      END DO
      IF(IRET.EQ.0.OR.JPRINT.GT.0)WRITE(IOUT,100)(XCM(I),I=1,3)
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
            X=X+X0(I,J)*X0(I,J)  ! chen
          END DO
          XK(IX,IX)=XK(IX,IX)+XM(I)*X
        END DO
      END DO
      DO IX=1,3
        DO IY=1,IX  ! chen
          DO I=1,NATOM
            XK(IY,IX)=XK(IY,IX)-XM(I)*X0(I,IY)*X0(I,IX)
          END DO
        END DO
      END DO
      DO IX=1,3
        DO IY=1,IX-1  ! chen
          XK(IX,IY)=XK(IY,IX)
        END DO
      END DO
C**MOMENT OF INERTIA TENSOR
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,150)
        DO I=1,3
          WRITE(IOUT,200)(XK(J,I),J=1,3)
        END DO
      END IF
      IF(XK(1,1).GT.1.D-35)THEN
        A=1/(2*XK(1,1))
      ELSE
        A=1.D+35
      END IF
      IF(XK(2,2).GT.1.D-35)THEN
        B=1/(2*XK(2,2))
      ELSE
        B=1.D+35
      END IF
      IF(XK(3,3).GT.1.D-35)THEN
        C=1/(2*XK(3,3))
      ELSE
        C=1.D+35
      END IF
C**SEARCH FOR LARGEST - REQUIRE 1 TO BE SMALLEST, 3 TO BE LARGEST
C**CORRESPONDING TO X,Y,Z
C**NOTE: XK(1,1) WOULD BE LARGEST, XK(3,3) WOULD BE SMALLEST
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'RE-ORDER ROTATIONAL CONSTANTS FOR X,Y,Z ASCENDING'
      END IF
      MROT(1)=1
      MROT(2)=2
      MROT(3)=3
      DO I=1,2
        DO J=I+1,3
          IF(XK(J,J).GT.XK(I,I))THEN
            SAVE=XK(I,I)
            XK(I,I)=XK(J,J)
            XK(J,J)=SAVE
            MSAVE=MROT(I)
            MROT(I)=MROT(J)
            MROT(J)=MSAVE
C**REORDER GEOMETRY
            DO K=1,NATOM
              SAVE=X0(K,I)
              X0(K,I)=X0(K,J)
              X0(K,J)=SAVE
            END DO
          END IF
        END DO
      END DO
C**RECALCULATE A,B,C
      IF(XK(1,1).GT.1.D-35)THEN
        A=1/(2*XK(1,1))
      ELSE
        A=1.D+35
      END IF
      IF(XK(2,2).GT.1.D-35)THEN
        B=1/(2*XK(2,2))
      ELSE
        B=1.D+35
      END IF
      IF(XK(3,3).GT.1.D-35)THEN
        C=1/(2*XK(3,3))
      ELSE
        C=1.D+35
      END IF
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,450)
        WRITE(IOUT,200)A*WAVENM,B*WAVENM,C*WAVENM
      END IF
C**NOW INTERCHANGE X & Y
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'INTERCHANGING X & Y TO Ir REPRESENTATION'
      END IF
      XSAVE=XK(1,1)
      XK(1,1)=XK(2,2)
      XK(2,2)=XSAVE
      DO I=1,NATOM
        XSAVE=X0(I,1)
        X0(I,1)=X0(I,2)
        X0(I,2)=XSAVE
      END DO
C**RECALCULATE A,B,C
      IF(XK(1,1).GT.1.D-35)THEN
        A=1/(2*XK(1,1))
      ELSE
        A=1.D+35
      END IF
      IF(XK(2,2).GT.1.D-35)THEN
        B=1/(2*XK(2,2))
      ELSE
        B=1.D+35
      END IF
      IF(XK(3,3).GT.1.D-35)THEN
        C=1/(2*XK(3,3))
      ELSE
        C=1.D+35
      END IF
      IF(IRET.EQ.0.OR.JPRINT.GT.0)THEN
        WRITE(IOUT,450)
        WRITE(IOUT,200)A*WAVENM,B*WAVENM,C*WAVENM
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE SORTIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/AXES/MX(3),MXDIP(5),MXROT(3)
C**SYMMETRIES OF ROTATIONAL FUNCTIONS D(0),iD(1+),D(1-),iD(2-)
C**FOR EVEN J
C**D(0)~1 THEREFORE ALWAYS TOTALLY SYMMETRIC (D(0:0) ~ 1)
      ISYMD(1,1)=1
C**iD(1+)~JX
      ISYMD(2,1)=MXROT(MX(1))
C**D(1-)~JY
      ISYMD(3,1)=MXROT(MX(2))
C**iD(2-)~JZ
      ISYMD(4,1)=MXROT(MX(3))
C**FOR ODD J
C**D(0) LIKE JZ
      ISYMD(1,2)=MXROT(MX(3))
C**iD(1+)~JY
      ISYMD(2,2)=MXROT(MX(2))
C**D(1-)~JX
      ISYMD(3,2)=MXROT(MX(1))
C**iD(2-) TOTALLY SYMMETRIC
      ISYMD(4,2)=1
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE PAXES(XM,X0,XX,RR,XL,NMODE,NATOM,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**TEMPORARY...LINEAR
      LOGICAL LINEAR
C**TEMPORARY...LINEAR
      LOGICAL PLANAR(3)
      DIMENSION XM(NATOM),X0(NATOM,3),XX(NATOM,3),RR(NATOM,NATOM)
      DIMENSION XL(NATOM,NMODE,3)
      DIMENSION XK(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/CHECK/MCHECK
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/AXES/MX(3),MXDIP(5),MXROT(3),ISYMT
C**TEMPORARY...LINEAR
      COMMON/TYPE/LINEAR,PLANAR
C**TEMPORARY...LINEAR
      COMMON/FILASS/IOUT,INP
200   FORMAT(3D20.12)
250   FORMAT(//,1X,'PAXES - MOMENT OF INERTIA TENSOR',/)
300   FORMAT(/,1X,'PRINCIPAL X,Y,Z AXES CORRESPOND TO SYMMETRY ',
     1'x,y,z AXES',3I3,/)
301   FORMAT(/,1X,'PRINCIPAL X,Y,Z AXES CORRESPOND TO SYMMETRIES ',3I3,
     1/)
302   FORMAT(/,1X,'I.R. OF ROTATIONAL EVEN-J FUNCTIONS D(0),iD(1+),'
     1'D(1-),iD(2-) ',4I3,/)
303   FORMAT(/,1X,'CHECK OF SPACE-FIXED Z-AXIS COMPONENTS (EVEN J)',
     13I3,/)
304   FORMAT(/,1X,'SYMMETRY OF SIN(M.TAU/2)',I3,/)
305   FORMAT(/,1X,'I.R. OF ROTATIONAL ODD-J FUNCTIONS D(0),iD(1+),'
     1'D(1-),iD(2-) ',4I3,/)
306   FORMAT(/,1X,'CHECK OF SPACE-FIXED Z-AXIS COMPONENTS (ODD J)',
     13I3,/)
350   FORMAT(/,1X,'LINEAR MOLECULE WAS ORIGINALLY NON-LINEAR')
355   FORMAT(/,1X,'NON-LINEAR MOLECULE WAS ORIGINALLY LINEAR')
      IF(IND.NE.0)THEN
        IF(LINEAR)ILIN=1
        IF(.NOT.LINEAR)ILIN=0
      END IF
C**SEARCH FOR LINEARITY ALONG 'Z'
      LINEAR=.TRUE.
      DO I=1,NATOM
C**LOOK AT 'X' AND 'Y'
        DO J=1,2
          IF(X0(I,J).NE.0)LINEAR=.FALSE.
        END DO
      END DO
      IF(IND.EQ.0)RETURN
C
      IF(LINEAR.AND.ILIN.EQ.0)THEN
        WRITE(IOUT,350)
        STOP 'LINEAR'
      END IF
      IF(.NOT.LINEAR.AND.ILIN.NE.0)THEN
        WRITE(IOUT,355)
        STOP 'NOT LINEAR'
      END IF
C     IF(MCHECK.EQ.0)CALL CHECKM(XM,X0,XX,RR,NATOM)
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
            X=X+X0(I,J)*X0(I,J)  ! chen
          END DO
          XK(IX,IX)=XK(IX,IX)+XM(I)*X
        END DO
      END DO
      do ix=1,3  ! chen
         do i=1,natom
            xk(ix,ix)=xk(ix,ix)-xm(i)*x0(i,ix)*x0(i,ix)
         end do
      end do
C      DO IX=1,3
C        DO IY=1,IX
C          XK(IX,IY)=XK(IY,IX)
C        END DO
C      END DO
C**MOMENT OF INERTIA TENSOR
      WRITE(IOUT,250)
      DO I=1,3
        WRITE(IOUT,200)(XK(J,I),J=1,3)
      END DO
C
C**DETERMINE SYMMETRY PROPERTIES OF PRINCIPAL AXES
C**P.A. AXES ARE IN THE ORDER: X=1(C), Y=2(B), Z=3(A)
C**SYMMETRY SPECIES ARE IN THE ORDER: A1=1, B2=2, B1=3, A2=4 (C2v)
C**                                   A=1, B=2 (C2)
C**                                   A'=1, A"=2 (Cs)
C**                                   A=1 (C1)
C**DIPOLES INPUT IN ORDER X=1, Y=2, Z=3
C**MXDIP(1)=SYMMETRY OF X (A1=1, B2=2, B1=3)
C**MXDIP(2)=SYMMETRY OF Y (       ''       )
C**MXDIP(3)=SYMMETRY OF Z (       ''       )
C
C**TEST FOR PLANARITY
      DO K=1,3
        PLANAR(K)=.TRUE.
      END DO
      DO K=1,3
        DO I=1,NATOM
          IF(XM(I).NE.0)PLANAR(K)=.FALSE.
        END DO
      END DO
C
C**SET DEFAULTS TO C1 SYMMETRY
C     MXDIP(1)=1
C     MXDIP(2)=1
C     MXDIP(3)=1
C     IF(NVSYM.NE.4)GO TO 50
C**C2V: FIND SYMMETRY AXIS
C**FIRST FIND 2 EQUIVALENT ATOMS
C     DO I=1,NATOM-1
C       DO J=I+1,NATOM
C         IF(XM(I).NE XM(J))GO TO 2
C**FOUND ATOM I EQUIVALENT TO ATOM J
C         IF(X0(I,1).EQ.X0(J,1))THEN
C**X-COORDS. EQUAL - CHECK IF SYMMETRY AXIS
C**FIRST SEE IF PLANAR IN X
C           IF(PLANAR(1))THEN
C             IF((X0(I,2).EQ.-X0(J,2)).AND.(X0(I,3).EQ.-X0(J,3)))THEN
C               MXDIP(1)=1
C**IN PLANAR CASE, FOUND TRUE C2 AXIS
C               GO TO 10
C             END IF
C           ELSE
C             IF((X0(I,2).EQ.-X0(J,2)).AND.(X0(I,3).EQ.-X0(J,3)))THEN
C               MXDIP(1)=1
C**IN NON-PLANAR CASE, MIGHT NOT BE TRUE C2 AXIS
C               GO TO 1
C             END IF
C           END IF
C         END IF
C         IF(X0(I,2).EQ.X0(J,2))THEN
C**Y-COORDS. EQUAL - CHECK IF SYMMETRY AXIS
C**FIRST SEE IF PLANAR IN Y
C           IF(PLANAR(2))THEN
C             IF((X0(I,1).EQ.-X0(J,1)).AND.(X0(I,3).EQ.-X0(J,3)))THEN
C               MXDIP(1)=2
C**IN PLANAR CASE, FOUND TRUE C2 AXIS
C               GO TO 10
C             END IF
C           ELSE
C             IF((X0(I,1).EQ.-X0(J,1)).AND.(X0(I,3).EQ.-X0(J,3)))THEN
C               MXDIP(1)=2
C**IN NON-PLANAR CASE, MIGHT NOT BE TRUE C2 AXIS
C               GO TO 1
C             END IF
C           END IF
C         END IF
C         IF(X0(I,3).EQ.X0(J,3))THEN
C**Z-COORDS. EQUAL - CHECK IF SYMMETRY AXIS
C**FIRST SEE IF PLANAR IN Z
C           IF(PLANAR(3))THEN
C             IF((X0(I,2).EQ.-X0(J,2)).AND.(X0(I,1).EQ.-X0(J,1)))THEN
C               MXDIP(1)=3
C**IN PLANAR CASE, FOUND TRUE C2 AXIS
C               GO TO 10
C             END IF
C           ELSE
C             IF((X0(I,2).EQ.-X0(J,2)).AND.(X0(I,1).EQ.-X0(J,1)))THEN
C               MXDIP(1)=3
C**IN NON-PLANAR CASE, MIGHT NOT BE TRUE C2 AXIS
C               GO TO 1
C             END IF
C           END IF
C         END IF
C**DON'T APPEAR TO HAVE FOUND C2 AXIS YET - KEEP TRYING
C         GO TO 2
C**FOUND POSSIBLE C2 SYMMETRY AXIS - CHECK IT OUT
1         CONTINUE
C**C2V: FIND VERTICAL PLANE OF SYMMETRY
C     GO TO 100
50    CONTINUE
C     IF(NVSYM.NE.2)GO TO 100
C**C2: FIND 2 EQUIVALENT ATOMS
C**CS: FIND PLANE OF SYMMETRY
100   CONTINUE
      WRITE(IOUT,301)MXDIP(1),MXDIP(2),MXDIP(3)
      WRITE(IOUT,300)MX(1),MX(2),MX(3)
      WRITE(IOUT,302)(ISYMD(I,1),I=1,4)
      WRITE(IOUT,303)ISYMP(MXDIP(1),ISYMD(3,1)),ISYMP(MXDIP(2),
     1ISYMD(2,1)),ISYMP(MXDIP(3),ISYMD(1,1))
      WRITE(IOUT,305)(ISYMD(I,2),I=1,4)
      WRITE(IOUT,306)ISYMP(MXDIP(1),ISYMD(3,2)),ISYMP(MXDIP(2),
     1ISYMD(2,2)),ISYMP(MXDIP(3),ISYMD(1,2))
      WRITE(IOUT,304)ISYMT
      IF(LDUMP.EQ.8)THEN
        MXDIP(1)=ISYMD(1,1)
        MXDIP(2)=ISYMD(2,1)
        MXDIP(3)=ISYMD(3,1)
        MXDIP(4)=ISYMD(4,1)
        MXDIP(5)=ISYMD(1,1)
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE HERMPT(NV,NVF,M,H,XQ,XJ,XV,XP,OMEGA,XQJ,MODE,NMODE,
     1NATOM,QQ,RR,XX,X0,XL,XM,NPOT,IPOT,JPOT,CPOT,ENERGY,MSYM,XTANPM,
     2WK,VINF,NP1,CP1,IP1,VP1,DP1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ABINIT
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      DIMENSION NP1(1),CP1(MMX1,1),IP1(MMX1,1),VP1(MMX1,1),DP1(MMX1,1)
      DIMENSION DUMMY(NMODE)
      DIMENSION H(NV,M,3),XQ(M),XJ(M),XP(NVF,M),XV(M),ENERGY(M,2)
      DIMENSION XQJ(M,M),XTANPM(NMODE),WK(2,NMODE,2)
      DIMENSION RR(NATOM,NATOM),QQ(NMODE),VINF(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/FILASS/IOUT,INP
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/CHKMAX/MAXCHK(6)
      COMMON/VMIN/VMIN
      COMMON/ABINIT/ABINIT
      COMMON/REACTL/JREACT,IFITRP
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/WHICH/IWHICH
      COMMON/PRINT/IPRINT,JPRINT
C**************************************************************
201   FORMAT(//,1X,50(1H*),/,
     11X,'INTEGRATION POINTS AND WEIGHTS FOR MODE ',I2,/)
202   FORMAT(13X,'Q(A0)',11X,'WEIGHT(EXP(-Q**2).D(Q))',3X,
     1'POTENTIAL(CM-1)',/)
203   FORMAT(2X,F20.12,2X,D25.12,2X,F15.2)
204   FORMAT(I4)
205   FORMAT(1X,'MODE: ',I2,' SYMMETRY: ',A3,' POINT: ',I3)
206   FORMAT(A2,1X,3F20.10)
207   FORMAT(41X,F12.8)
208   FORMAT(2F20.10)
C*****************************************************************
      PIQ=1.33133536380038D0
      YLAM=1/SQRT(OMEGA)
      DO 10 J=1,M
      DO 10 I=1,M
10    XQJ(I,J)=0.D0
      M1=M-1
      DO 20 J=1,M1
      I=J+1
      XQJ(J,I)=SQRT(DFLOAT(J)/2.D0)
      XQJ(I,J)=XQJ(J,I)
20    CONTINUE
30    CALL DIAG(XQJ,XQJ,M,M,-1,XJ,XQ,XJ,M,M,XQJ,M,M,XQJ,XQJ,XQJ,IDUM,
     1IDUM,IDUM)
      IF(IPRINT.GT.0)WRITE(IOUT,201)MODE
      IF (IPRINT.GT.0) WRITE(IOUT,202)
      DO I=1,NWSYM
        DO K=1,NSYM(I)
          IF(ISYM(I,K).EQ.MODE)MSYM=I
        END DO
      END DO
      IF(MSYM.GT.4)MSYM=MSYM-4
C***********************
      K=0
C******************************************************
      DO 60 J=1,M
      XJ(J)=PIQ*PIQ*XQJ(1,J)*XQJ(1,J)
      Q=XQ(J)
      XQ(J)=Q*YLAM
C.....................................................??????
CC    IF((MOLPRO.NE.0.OR.IWHICH.LT.0).AND.J.EQ.1)
      IF((IABS(MOLPRO).EQ.1.OR.IABS(MOLPRO).EQ.2.OR.IFITRP.NE.0)
     1.AND.J.EQ.1)XTANPM(MODE)=XTANPM(MODE)/DABS(XQ(1))
C.....................................................??????
C**ENERGY THIS POINT
      DO K=1,NMODE
        QQ(K)=0
      END DO
      QQ(MODE)=XQ(J)
C*********************************************************************
C*********************************************************************
C**MOLPRO LINK INDEX = 1 AND INDEX = 2
C
      IF(IABS(MOLPRO).EQ.1.OR.IABS(MOLPRO).EQ.2)THEN
C**ONLY DUMP IF UNIQUE SYMMETRY 
        ISKIP=0
        IF(J.LE.M/2.OR.(MSYM.EQ.1))THEN
          IF(IABS(MOLPRO).EQ.1)THEN
            IF(IWHICH.GT.0)THEN
C**FITTING GLOBAL POTENTIAL
              DO I=1,NATOM
                DO K=1,3
                  XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQ(MODE)/
     1            SQRT(XM(I))
                END DO
              END DO
              CALL GETPOT(VVV,NATOM,XX,RR)
              IF(J.EQ.1.OR.J.LT.M)THEN
                QQQ=QQ(MODE)+1.D-6
                DO I=1,NATOM
                  DO K=1,3
                    XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQQ/
     1              SQRT(XM(I))
                  END DO
                END DO
                CALL GETPOT(VVP,NATOM,XX,RR)
              END IF
              IF(J.EQ.1)DDD=(VVP-VVV)*1.D6
              IF(J.EQ.M.OR.J.GT.1)THEN
                QQQ=QQ(MODE)-1.D-6
                DO I=1,NATOM
                  DO K=1,3
                    XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQQ/
     1              SQRT(XM(I))
                  END DO
                END DO
                CALL GETPOT(VVM,NATOM,XX,RR)
              END IF
              IF(J.EQ.M)DDD=(VVV-VVM)*1.D6
              IF(J.GT.1.AND.J.LT.M)DDD=(VVP-VVM)*1.D6/2

              WRITE(MOUT,204)NATOM
              WRITE(MOUT,205)MODE,CHSYM(MSYM),J
              IF(MOLPRO.GT.0)WRITE(MOUT,208)VVV
              IF(MOLPRO.LT.0)WRITE(MOUT,208)VVV,DDD
            ELSE
C**FITTING ABINITIO POTENTIAL
              WRITE(MOUT,204)NATOM
              WRITE(MOUT,205)MODE,CHSYM(MSYM),J
              WRITE(MOUT,*)'**********************************'
            END IF
          END IF

          IF(MOLPRO.EQ.2)THEN
C**FIT ABINITIO POTENTIAL
            READ(MINP,204)NATOM
            READ(MINP,205)MODE,CHSYM(MSYM),JDUM
            READ(MINP,*)ENERGY(J,1)
          END IF
          IF(MOLPRO.EQ.-2)THEN
C**INTERPOLATE ABINITIO POTENTIAL
            READ(MINP,204)NATOM
            READ(MINP,205)MODE,CHSYM(MSYM),JDUM
            READ(MINP,*)ENERGY(J,1),ENERGY(J,2)
          END IF
        ELSE
C**SKIP THIS POINT ON SYMMETRY GROUNDS
          ISKIP=1
        END IF
      END IF
C***********************
      DO I=1,NATOM
        DO K=1,3
          XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQ(MODE)/
     1    SQRT(XM(I))
        END DO
C**MOLPRO LINK INDEX = 1 AND INDEX = 2
C
        IF(IABS(MOLPRO).EQ.1.AND.ISKIP.EQ.0)THEN
          WRITE(MOUT,206)SYMBOL(I),(XX(I,K)*BOHR,K=1,3)
        END IF
        IF(IABS(MOLPRO).EQ.2.AND.ISKIP.EQ.0)THEN
          READ(MINP,206)SYMBAD,(XX(I,K),K=1,3)
          DO K=1,3
            XX(I,K)=XX(I,K)/BOHR
          END DO
        END IF
C
      END DO
C***********************
      IF(IWHICH.GT.0)THEN
        IF(ABINIT)THEN
          CALL GETAPT(VPT,NATOM,XX,RR,QQ,NMODE,0)
        ELSE
          CALL GETPOT(VPT,NATOM,XX,RR)
        END IF
      ELSE
        IF(IWHICH.EQ.0)THEN
          CALL GETPT1(VPT,NPOT,IPOT,JPOT,CPOT,0,QQ)
        ELSE
          IF(MOLINC.GT.0)
     1    CALL GETPQT(VPT,NMODE,QQ,XTANPM,NP1,CP1,IP1,NMODE,MMX1)
          IF(MOLINC.EQ.0)
     1    CALL GETQPT(VPT,NMODE,QQ,XTANPM)
          IF(MOLINC.LT.0)
     1    CALL GETQP1(VPT,DUMMY,NMODE,QQ,NP1,CP1,VP1,DP1,NMODE,MMX1)
        END IF
      END IF
      XV(J)=VPT
      IF (IPRINT.GT.0) WRITE(IOUT,203)XQ(J),XJ(J),VPT*WAVENM
C**ENERGY THIS POINT
      IF(-MAXCHK(1).NE.MODE.AND.VPT.LT.VMIN)VMIN=VPT
      IF(J.EQ.1)THEN
        VMAX=VPT
        WK(1,MODE,2)=VPT
      END IF
      IF(J.EQ.M)THEN
        WK(2,MODE,2)=VPT
      END IF
      IF(J.EQ.M.AND.VPT.LT.VMAX)VMAX=VPT
60    CONTINUE
C******************************************************************
      MODCHK=IABS(MAXCHK(1))
      IF(MODE.EQ.MODCHK)THEN
        JCPSV=JCOUPL
        JCOUPL=1
        MX=M/2
        IF(MAXCHK(1).GT.0)THEN
          WRITE(IOUT,*)'START GRID INSPECTION FOR MODE ',MODCHK
        ELSE
          WRITE(IOUT,*)'START HOLE SEARCH FOR MODE ',MODCHK
        END IF
CC      IF(MAXCHK(1).LT.0)THEN
C
          VLAST=-1.0D+20
C**MIN=0 UNTIL FIRST MINIMUM (+1) or MAXIMUM (-1) REACHED
          MIN=0
          GRAD=0
          VDIFF=0
C**FIRST GO TO NEGATIVE Q
          DO J=MX,1,-1
            VDP=VFUN1(XV,XV,M,J,VMIN,WAVENM)
C**PREVIOUS GRADIENT
            GRADPR=GRAD
            IF(J.NE.MX)THEN
C**GRAD -VE IF V RISING, +VE IF V FALLING
              GRAD=WAVENM*(VDP-VLAST)/(XQ(J)-XQ(J+1))
            END IF
C**TEST FOR V RISING OR FALLING BEFORE FIRST MINIMUM
            IF(VDP.GT.VLAST.OR.MIN.EQ.0)THEN
              VDIFF=VDP-VLAST
              CALL SET71(XV,XV,M,VDP,J)
              IF(VDP.LT.VMIN)THEN
                WRITE(IOUT,*)'SUSPECT V FOR MODE ',MODE
                WRITE(IOUT,*)'M1 ',J
                WRITE(IOUT,*)'Q1 ',XQ(J)
                WRITE(IOUT,*)'POTENTIAL ',(VDP)*WAVENM
                WRITE(IOUT,*)'*****************************'
                WRITE(IOUT,*)
C**SET VMIN IF NOT REACHED FIRST MINIMUM
                IF(MIN.EQ.0)VMIN=VDP
              END IF
C**TEST FOR MINIMUM
              IF(GRADPR.GE.0.AND.GRAD.LT.0)THEN
                WRITE(IOUT,*)'MINIMUM FOUND FOR MODE ',MODE
                WRITE(IOUT,*)'M1 ',J+1
                WRITE(IOUT,*)'Q1 ',XQ(J+1)
                WRITE(IOUT,*)'POTENTIAL ',(VLAST)*WAVENM
                WRITE(IOUT,*)'*****************************'
                WRITE(IOUT,*)
                MIN=1
              END IF
              VLAST=VDP
C**TEST FOR V FALLING AFTER FIRST MINIMUM
            ELSE
              IF(MAXCHK(1).GT.0)THEN
                CALL SET71(XV,XV,M,VDP,J)
              ELSE
                CALL SET71(XV,XV,M,VLAST+VDIFF,J)
              END IF
              IF(GRADPR.LT.0.AND.GRAD.GE.0)THEN
                WRITE(IOUT,*)'MAXIMUM FOUND FOR MODE ',MODE
                WRITE(IOUT,*)'M1 ',J+1
                WRITE(IOUT,*)'Q1 ',XQ(J+1)
                WRITE(IOUT,*)'POTENTIAL ',(VLAST)*WAVENM
                WRITE(IOUT,*)'*****************************'
                WRITE(IOUT,*)
                MIN=-1
              END IF
              IF(VLAST+VDIFF.LT.VMIN)THEN
                WRITE(IOUT,*)'SUSPECT V FOR MODE ',MODE
                WRITE(IOUT,*)'M1 ',J+1
                WRITE(IOUT,*)'Q1 ',XQ(J+1)
                WRITE(IOUT,*)'POTENTIAL ',(VLAST)*WAVENM
                WRITE(IOUT,*)'*****************************'
                WRITE(IOUT,*)
              END IF
              IF(MAXCHK(1).LT.0)VLAST=VLAST+VDIFF
            END IF
          END DO
C
          VLAST=-1.0D+20
C**MIN=0 UNTIL FIRST MINIMUM (+1) or MAXIMUM (-1) REACHED
          MIN=0
          GRAD=0
          VDIFF=0
C**NOW GO TO POSITIVE Q
          DO J=MX+1,2*MX
            VDP=VFUN1(XV,XV,M,J,VMIN,WAVENM)
C**PREVIOUS GRADIENT
            GRADPR=GRAD
            IF(J.NE.MX+1)THEN
C**GRAD +VE IF V RISING, -VE IF V FALLING
              GRAD=WAVENM*(VDP-VLAST)/(XQ(J)-XQ(J-1))
            END IF
C**TEST FOR V RISING OR FALLING BEFORE FIRST MINIMUM
            IF(VDP.GT.VLAST.OR.MIN.EQ.0)THEN
              VDIFF=VDP-VLAST
              CALL SET71(XV,XV,M,VDP,J)
              IF(VDP.LT.VMIN)THEN
                WRITE(IOUT,*)'SUSPECT V FOR MODE ',MODE
                WRITE(IOUT,*)'M1 ',J
                WRITE(IOUT,*)'Q1 ',XQ(J)
                WRITE(IOUT,*)'POTENTIAL ',(VDP)*WAVENM
                WRITE(IOUT,*)'*****************************'
                WRITE(IOUT,*)
C**SET VMIN IF NOT REACHED FIRST MINIMUM
                IF(MIN.EQ.0)VMIN=VDP
              END IF
C**TEST FOR MINIMUM
              IF(GRADPR.LT.0.AND.GRAD.GE.0)THEN
                WRITE(IOUT,*)'MINIMUM FOUND FOR MODE ',MODE
                WRITE(IOUT,*)'M1 ',J-1
                WRITE(IOUT,*)'Q1 ',XQ(J-1)
                WRITE(IOUT,*)'POTENTIAL ',(VLAST)*WAVENM
                WRITE(IOUT,*)'*****************************'
                WRITE(IOUT,*)
                MIN=1
              END IF
              VLAST=VDP
C**TEST FOR V FALLING AFTER FIRST MINIMUM
            ELSE
              IF(MAXCHK(1).GT.0)THEN
                CALL SET71(XV,XV,M,VDP,J)
              ELSE
                CALL SET71(XV,XV,M,VLAST+VDIFF,J)
              END IF
              IF(GRADPR.GE.0.AND.GRAD.LT.0)THEN
                WRITE(IOUT,*)'MAXIMUM FOUND FOR MODE ',MODE
                WRITE(IOUT,*)'M1 ',J-1
                WRITE(IOUT,*)'Q1 ',XQ(J-1)
                WRITE(IOUT,*)'POTENTIAL ',(VLAST)*WAVENM
                WRITE(IOUT,*)'*****************************'
                WRITE(IOUT,*)
                MIN=-1
              END IF
              IF(VLAST+VDIFF.LT.VMIN)THEN
                WRITE(IOUT,*)'SUSPECT V FOR MODE ',MODE
                WRITE(IOUT,*)'M1 ',J-1
                WRITE(IOUT,*)'Q1 ',XQ(J-1)
                WRITE(IOUT,*)'POTENTIAL ',(VLAST)*WAVENM
                WRITE(IOUT,*)'*****************************'
                WRITE(IOUT,*)
              END IF
              IF(MAXCHK(1).LT.0)VLAST=VLAST+VDIFF
            END IF
          END DO
C
CC      END IF
        DO J=1,M
          VDP=VFUN1(XV,XV,M,J,VMIN,WAVENM)-VMIN
          IF(VDP.LT.0.AND.MAXCHK(1).LT.0)THEN
            WRITE(IOUT,*)'STILL HAVE SUSPECT V FOR MODE ',
     1      MODE
            WRITE(IOUT,*)'M1 ',J
            WRITE(IOUT,*)'Q1 ',XQ(J)
            WRITE(IOUT,*)'POTENTIAL ',(VDP+VMIN)*WAVENM
            WRITE(IOUT,*)'*****************************'
            WRITE(IOUT,*)
          END IF
        END DO
        IF(MAXCHK(1).GT.0)THEN
          WRITE(IOUT,*)' ORIGINAL POTENTIAL'
          WRITE(IOUT,201)MODE
          WRITE(IOUT,202)
          DO J=1,M
            WRITE(IOUT,203)XQ(J),XJ(J),VFUN1(XV,XV,M,J,VMIN,WAVENM)
     1      *WAVENM
          END DO
          WRITE(IOUT,*)'END GRID INSPECTION FOR MODE ',MODCHK
        ELSE
          WRITE(IOUT,*)' AMMEMDED POTENTIAL'
          WRITE(IOUT,201)MODE
          WRITE(IOUT,202)
          DO J=1,M
            WRITE(IOUT,203)XQ(J),XJ(J),VFUN1(XV,XV,M,J,VMIN,WAVENM)
     1      *WAVENM
          END DO
          WRITE(IOUT,*)'END HOLE SEARCH FOR MODE ',MODCHK
        END IF
        WRITE(IOUT,*)
        JCOUPL=JCPSV
      END IF

      VINF(MODE)=VMIN-VMAX/5
      WK(1,MODE,1)=XQ(1)
      WK(2,MODE,1)=XQ(M)
C******************************************************
      DO 100 J=1,M
      Q=(XQ(J))/YLAM
      H(1,J,1)=1.D0/PIQ
      H(2,J,1)=2.D0*Q/(PIQ*SQRT(2.D0))
      IF(NV.GE.3)THEN
      DO 70 I=3,NV
      I1=I-1
      I2=I-2
70    H(I,J,1)=2*(Q*H(I1,J,1)-I2*H(I2,J,1)/SQRT(2.D0*I2))/SQRT(2.D0*I1)
      END IF
      H(1,J,2)=0.D0
      H(2,J,2)=2.D0/(PIQ*SQRT(2.D0))
      H(1,J,3)=0.D0
      H(2,J,3)=0.D0
      IF(NV.GE.3)THEN
      DO 80 I=3,NV
      I1=I-1
      I2=I-2
      H(I,J,2)=2.D0*I1*H(I1,J,1)/SQRT(2.D0*I1)
      H(I,J,3)=4.D0*I1*I2*H(I2,J,1)/SQRT(4.D0*I1*I2)
80    CONTINUE
      END IF
      DO 90 I=1,NV
      H(I,J,3)=H(I,J,3)-2*Q*H(I,J,2)+(Q*Q-1.D0)
     1*H(I,J,1)
90    H(I,J,2)=H(I,J,2)-Q*H(I,J,1)
100   CONTINUE
C******************************************************
      DO 110 J=1,M
      XJS=SQRT(XJ(J))
      Q=XQ(J)/YLAM
      DO 110 I=1,NV
      H(I,J,1)=H(I,J,1)*XJS
      H(I,J,2)=H(I,J,2)*XJS/YLAM
      H(I,J,3)=H(I,J,3)*XJS/(YLAM*YLAM)
110   CONTINUE
C**WEIGHT
      DO J=1,M
        XJS=SQRT(XJ(J))
        Q=XQ(J)/YLAM
        DO I=1,NVF
          XP(I,J)=XJS/EXP(-Q*Q/2)
        END DO
      END DO
C******************************************************
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE MORSPT(NBF1,M1,H1,XQ1,AN,XV,OMEGA,XK,MODE,NMODE,
     1NATOM,QQ,RR,XX,X0,XL,XM,NPOT,IPOT,JPOT,CPOT,ENERGY,MSYM,XTANPM,
     2WK,VINF,NP1,CP1,IP1,VP1,DP1,MMX1,LPREV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ABINIT
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      DIMENSION NP1(1),CP1(MMX1,1),IP1(MMX1,1),VP1(MMX1,1),DP1(MMX1,1)
      DIMENSION DUMMY(NMODE)
      DIMENSION H1(NBF1,M1,3,2),XQ1(M1),AN(M1),ENERGY(M1)
      DIMENSION XK(M1,M1),XV(M1),XTANPM(NMODE),WK(2,NMODE)
      DIMENSION RR(NATOM,NATOM),QQ(NMODE),VINF(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/FILASS/IOUT,INP
      COMMON/VMIN/VMIN
      COMMON/ABINIT/ABINIT
      COMMON/WHICH/IWHICH
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C**************************************************************
500   FORMAT(//,1X,50(1H*),/,
     11X,'LINEAR-BEND INTEGRATION POINTS AND WEIGHTS FOR MODE ',I2,/)
501   FORMAT(8X,'Q(A0)',10X,'MORSE WEIGHT(Q**L.EXP(-Q).D(Q))',3X,
     1'POTENTIAL(CM-1)',/)
502   FORMAT(2X,F20.12,2X,D25.12,2X,F15.2)
503   FORMAT(5X,'S =',F10.3,' IS TOO BIG - MORSE FAILED',/)
504   FORMAT(//,1X,'2*SMAX= ',F7.2,' QMAX= ',F7.2,/)
C**************************************************************
      YLAM=1/SQRT(OMEGA)
      IF(LPREV.NE.0)GO TO 1000
C**SINGLE SET OF POINTS FOR L=0
      S2MAX=0
      SMAX=0
      IF(IPRINT.GT.0)THEN
        WRITE(IOUT,500)MODE
        WRITE(IOUT,501)
        CALL FLUSH(IOUT)
      END IF
      DO I=1,M1
        DO J=1,M1
          XK(J,I)=0.D0
        END DO
      END DO
      M=M1-1
      DO J=1,M
        I=J+1
        XK(J,I)=-SQRT(J*(J+S2MAX))
        XK(I,J)=XK(J,I)
      END DO
      DO J=1,M1
        XK(J,J)=2*J-1+S2MAX
      END DO
      CALL DIAG(XK,XK,M1,M1,-1,AN,XQ1,AN,M1,M1,
     1XK,M1,M1,XK,XK,XK,IDUM,IDUM,IDUM)
C**TWO SETS OF FUNCTIONS FOR L=0 AND L=1
      DO IJK=1,2
        S2=IJK-1
        S=S2/2
        DO J=1,M1
          IF(IJK.EQ.1)THEN
            Q=DABS(XQ1(J))
            QBAR=SQRT(Q)
            XQ1(J)=QBAR*YLAM
            XJS=XK(1,J)
            XJS2=XJS*XJS
            XJS=SQRT(XJS2)
C**TEMPORARY
C           XJS=XJS/SQRT(2.D0)
C**TEMPORARY
            XK(1,J)=XJS
          ELSE
            QBAR=XQ1(J)/YLAM
            Q=QBAR*QBAR
            XJS=XK(1,J)
          END IF
C**ENERGY THIS POINT
          DO K=1,NMODE
            QQ(K)=0
          END DO
          QQ(MODE)=XQ1(J)
          DO I=1,NATOM
            DO K=1,3
              XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQ(MODE)/
     1        SQRT(XM(I))
            END DO
          END DO
          IF(IWHICH.GT.0)THEN
            IF(ABINIT)THEN
              CALL GETAPT(VPT,NATOM,XX,RR,QQ,NMODE,0)
            ELSE
              CALL GETPOT(VPT,NATOM,XX,RR)
            END IF
          ELSE
            IF(IWHICH.EQ.0)THEN
              CALL GETPT1(VPT,NPOT,IPOT,JPOT,CPOT,0,QQ)
            ELSE
              IF(MOLINC.GT.0)
     1        CALL GETPQT(VPT,NMODE,QQ,XTANPM,NP1,CP1,IP1,NMODE,MMX1)
              IF(MOLINC.EQ.0)
     1        CALL GETQPT(VPT,NMODE,QQ,XTANPM)
              IF(MOLINC.LT.0)
     1        CALL GETQP1(VPT,DUMMY,NMODE,QQ,NP1,CP1,VP1,DP1,NMODE,
     2        MMX1)
            END IF
          END IF
          XV(J)=VPT
          IF(IJK.EQ.1)THEN
            IF (IPRINT.GT.0) 
     1      WRITE(IOUT,502)XQ1(J),XJS2,VPT*WAVENM
C**ENERGY THIS POINT
            IF(VPT.LT.VMIN)VMIN=VPT
            IF(J.EQ.1)VMAX=VPT
            IF(J.EQ.M.AND.VPT.LT.VMAX)VMAX=VPT
          END IF

C**START MORSE
C         DO II=1,NBF1
CC**NBF GOES FROM NBF1 TO 1 (INDEX FOR ARRAY)
C           NBF=NBF1+1-II
CC**N GOES FROM NBF1-1 TO 0 (FOR MORSE N IS QUANTUM...STARTS 0)
C           N=NBF-1
CC**FOR MORSE N STARTS 0 (N + 2S = 2S)
CC**FOR LINEAR BEND N STARTS S (N + S = 2S)
CC**N + 2S (N=0) IS THE SAME AS N + S (N=S)
C           A=N+S2
C           AL=A+1.D0
C           CS=S2+1
C           BOT=1.D0
C           H1(1,J,1,IJK)=1
C           H1(1,J,2,IJK)=0
C           H1(1,J,3,IJK)=0
C           IF(NBF.GE.2)THEN
C             H1(2,J,1,IJK)=S2+1-Q
C             H1(2,J,3,IJK)=0
C           END IF
C           DO I=1,NBF
C             BOT=BOT*I
C             IF(I.GE.3)THEN
C               I1=I-1
C               I2=I-2
C               H1(I,J,1,IJK)=((2*I1+S2-1-Q)*H1(I1,J,1,IJK)-(I1+S2-1)*
C     1         H1(I2,J,1,IJK))/I1
C             END IF
C           END DO
C           AN(NBF)=SQRT(BOT)*SQRT(2.D0)/(SQRT(NBF*1.D0))
C1000       CS=CS-1.D0
C           IF(CS.LE.1.0)GO TO 2000
C           AL=AL-1.D0
C           AN(NBF)=(AN(NBF)/SQRT(AL))*SQRT(CS)
C           GO TO 1000
C2000       GAMMA=FACTL(CS,IFAIL)
C           AN(NBF)=AN(NBF)*SQRT(GAMMA)
C           IFAIL=0
C           DUMAL=AL
C           DUMM=FACTL(DUMAL,IFAIL)
C           IF(IFAIL.NE.0)THEN
C             WRITE(IOUT,503)S2/2
C             STOP 'ERROR IN MORSPT'
C           END IF
C3000       AL=AL-1.D0
C           IF(AL.LE.1.0)GO TO 4000
C           AN(NBF)=AN(NBF)/SQRT(AL)
C           GO TO 3000
C4000       GAMMA=FACTL(AL,IFAIL)
C           AN(NBF)=AN(NBF)/SQRT(GAMMA)
C           IF(NBF.GE.2)THEN
C             H1(NBF,J,2,IJK)=(N*H1(NBF,J,1,IJK)-(N+S2)*H1(N,J,1,IJK))
C     1       /Q
C           END IF
C           IF(NBF.GE.3)THEN
C              H1(NBF,J,3,IJK)=-((S2+1-Q)*H1(NBF,J,2,IJK)+N*
C     1        H1(NBF,J,1,IJK))/Q
C           END IF
C           H1(NBF,J,1,IJK)=H1(NBF,J,1,IJK)*AN(NBF)
C           H1(NBF,J,2,IJK)=H1(NBF,J,2,IJK)*AN(NBF)
C           H1(NBF,J,3,IJK)=H1(NBF,J,3,IJK)*AN(NBF)
C         END DO
C         FACT=S-SMAX
C         EXPO=Q**FACT
C         DO NBF=1,NBF1
C           H1(NBF,J,3,IJK)=(H1(NBF,J,3,IJK)+H1(NBF,J,2,IJK)*(2*S/Q-1)+
C     1     H1(NBF,J,1,IJK)*(S*(S-1)/(Q*Q)-S/Q+0.25D0))*EXPO
C           H1(NBF,J,2,IJK)=(H1(NBF,J,2,IJK)+H1(NBF,J,1,IJK)*
C     1     (S/Q-0.5D0))*EXPO
C           H1(NBF,J,1,IJK)=H1(NBF,J,1,IJK)*EXPO
C         END DO
C         DO I=1,NBF1
C           H1(I,J,3,IJK)=4*QBAR*(Q*H1(I,J,3,IJK)+H1(I,J,2,IJK))
C           H1(I,J,2,IJK)=2*Q*H1(I,J,2,IJK)
C         END DO
C         DO I=1,NBF1
C           H1(I,J,2,IJK)=H1(I,J,2,IJK)/QBAR
C           H1(I,J,3,IJK)=H1(I,J,3,IJK)/QBAR
C           H1(I,J,1,IJK)=H1(I,J,1,IJK)*XJS
C           H1(I,J,2,IJK)=H1(I,J,2,IJK)*XJS/YLAM
C           H1(I,J,3,IJK)=H1(I,J,3,IJK)*XJS/(YLAM*YLAM)
C         END DO
C**END MORSE

C**START LAGUERRE
C         DO 12 I=1,M2
          X=QBAR*QBAR
          IF(IJK.EQ.1)THEN
C         H1(1,J,1,IJK)=SQRT(2.D0)
C         H1(2,J,1,IJK)=SQRT(2.D0)*(1-X)
          H1(1,J,1,IJK)=1
          H1(2,J,1,IJK)=(1-X)
          H1(1,J,2,IJK)=0
          H1(1,J,3,IJK)=0
          END IF
C         I1=NV2+1
C         I2=NV2+2
          IF(IJK.EQ.2)THEN
C         H1(1,J,1,IJK)=SQRT(2.D0)
C         H1(2,J,1,IJK)=SQRT(2.D0)*(2-X)
          H1(1,J,1,IJK)=1
          H1(2,J,1,IJK)=(2-X)
          H1(1,J,2,IJK)=0
          H1(1,J,3,IJK)=0
          END IF
C         I1=NV2*2+1
C         I2=NV2*2+2
C         F(I1,I,1)=SQRT(2.)
C         F(I2,I,1)=SQRT(2.)*(3.-X)
C         F(I1,I,2)=0.
C         F(I1,I,3)=0.
C         JK=JJ+1
C         DO 17 K=1,JK
C         IA=K-1
          IA=IJK-1
          DO 17 I=3,NBF1
          I1=I-1
          I2=I-2
C         IJ=(K-1)*NV2+J
C         IJ1=IJ-1
C         IJ2=IJ-2
17        H1(I,J,1,IJK)=((2*I2+1+IA-X)*H1(I1,J,1,IJK)-
     1    (I2+IA)*H1(I2,J,1,IJK))/I1
C         DO 18 K=1,JK
C         IA=K-1
          IA=IJK-1
          DO 18 I=2,NBF1
          I1=I-1
C         IJ=(K-1)*NV2+J
C         IJ1=IJ-1
          H1(I,J,2,IJK)=(I1*H1(I,J,1,IJK)-(I1+IA)*H1(I1,J,1,IJK))/X
18        H1(I,J,3,IJK)=((I1-1)*H1(I,J,2,IJK)-(I1+IA)*H1(I1,J,2,IJK))/X
C         DO 19 K=1,JK
          DO 19 I=1,NBF1
C         IJ=(K-1)*NV2+J
          H1(I,J,3,IJK)=2*(H1(I,J,2,IJK)+2*X*H1(I,J,3,IJK))
19        H1(I,J,2,IJK)=H1(I,J,2,IJK)*2*QBAR
C         DO 21 K=1,JK
C         AJ=K-1
          AJ=IJK-1
          DO 21 I=1,NBF1
C         L=(K-1)*NV2+J
          H1(I,J,3,IJK)=(AJ+AJ+1)*(H1(I,J,2,IJK)-QBAR*H1(I,J,1,IJK))+
     1    QBAR*(H1(I,J,3,IJK)-2*QBAR*H1(I,J,2,IJK)+(X-1)*H1(I,J,1,IJK))
          H1(I,J,2,IJK)=(AJ-X)*H1(I,J,1,IJK)+QBAR*H1(I,J,2,IJK)
CCCC      F(L,I,4)=F(L,I,1)
          H1(I,J,1,IJK)=QBAR*H1(I,J,1,IJK)
          DO K=1,3
            H1(I,J,K,IJK)=H1(I,J,K,IJK)*QBAR**(IA-1)
          END DO
21        CONTINUE
          T=1
          B=1
C         DO 23 K=2,JK
          IF(IJK.EQ.2)THEN
C         L=K-1
          L=1
          B=B*L
          Z=T/B
          DO 24 I=1,NBF1
          N=I-1
C         IJ=(K-1)*NV2+J
C         DO 25 M=1,4
          DO 25 M=1,3
25        H1(I,J,M,IJK)=H1(I,J,M,IJK)*SQRT(Z)
          AAT=N+1
          AAB=L+N+1
24        Z=Z*AAT/AAB
          END IF
23        CONTINUE
          DO I=1,NBF1
C           H1(I,J,2,IJK)=H1(I,J,2,IJK)/QBAR
C           H1(I,J,3,IJK)=H1(I,J,3,IJK)/QBAR
            H1(I,J,1,IJK)=H1(I,J,1,IJK)*XJS
            H1(I,J,2,IJK)=H1(I,J,2,IJK)*XJS/YLAM
            H1(I,J,3,IJK)=H1(I,J,3,IJK)*XJS/(YLAM*YLAM)
          END DO
12        CONTINUE
C**END LAGUERRE
        END DO
        IF(IJK.EQ.1)THEN
          VINF(MODE)=VMIN-VMAX/5
          WK(1,MODE)=XQ1(1)
          WK(2,MODE)=XQ1(M)
        END IF
      END DO

C**VERY TEMPORARY
1000  CONTINUE
C**ASSUME PSI(L) IN STORAGE POSITION 1
C**GENERATE PSI(L+1) IN STORAGE POSITION 2
      L=LPREV
      DO J=1,M1
        Q=XQ1(J)
        QBAR=Q/YLAM
        Q2=Q*Q
        DO I=1,NBF1
C**CONSTRUCT Q.L(n,0).exp(-Q^2)/2.exp(il.delta) etc.
C**CONVERT TO DERIVATIVES WRT QBAR
C         H1(I,J,1,2)=H1(I,J,1,1)
C         H1(I,J,2,2)=H1(I,J,2,1)*YLAM
C         H1(I,J,3,2)=H1(I,J,3,1)*YLAM*YLAM
C**FORM SECOND DERIVATIVE
C         H1(I,J,3,2)=H1(I,J,3,2)-H1(I,J,2,2)/QBAR
C         H1(I,J,3,2)=QBAR*H1(I,J,3,2)+2*H1(I,J,2,2)
C**FORM FIRST DERIVATIVE
C         H1(I,J,2,2)=QBAR*H1(I,J,2,2)+H1(I,J,1,2)
C**FORM FUNCTION
C         H1(I,J,1,2)=QBAR*H1(I,J,1,2)
C**FORM SECOND DERIVATIVE COMBINATION
C         H1(I,J,3,2)=H1(I,J,3,2)+H1(I,J,2,2)/QBAR-H1(I,J,1,2)/
C    1    (QBAR*QBAR)
C**CONVERT TO DERIVATIVES WRT Q
C         H1(I,J,2,2)=H1(I,J,2,2)/YLAM
C         H1(I,J,3,2)=H1(I,J,3,2)/(YLAM*YLAM)
          H1(I,J,1,2)=QBAR*H1(I,J,1,1)
          H1(I,J,2,2)=QBAR*H1(I,J,2,1)+H1(I,J,1,1)/YLAM
          H1(I,J,3,2)=QBAR*H1(I,J,3,1)+2*H1(I,J,2,1)/YLAM-
     1    2*L*H1(I,J,1,1)/(QBAR*YLAM*YLAM)
        END DO
      END DO
      DO I=1,NBF1
C**SCHMIDT
        IF(I.GT.1)THEN
          DO J=I-1,1,-1
            X=0
            DO M=1,M1
              X=X+H1(I,M,1,2)*H1(J,M,1,2)
            END DO
            DO M=1,M1
              DO K=1,3
                H1(I,M,K,2)=H1(I,M,K,2)-X*H1(J,M,K,2)
              END DO
            END DO
          END DO
        END IF
C**NORMALISE
        X=0
        DO M=1,M1
          X=X+H1(I,M,1,2)*H1(I,M,1,2)
        END DO
        A=1/SQRT(X)
        DO M=1,M1
          DO K=1,3
            H1(I,M,K,2)=H1(I,M,K,2)*A
          END DO
        END DO
      END DO
C**CHECK ORTHONORMALITY
C     DO I=1,NBF1
C       DO J=1,NBF1
C         XK(J,I)=0
C         DO M=1,M1
C           XK(J,I)=XK(J,I)+H1(J,M,1,2)*H1(I,M,1,2)
C         END DO
C       END DO
C       WRITE(IOUT,*)(XK(J,I),J=1,I)
C     END DO
C**VERY TEMPORARY

      RETURN
      END
C***********************************************************
C***********************************************************
      DOUBLE PRECISION FUNCTION FACTL(C,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(8)
      COMMON/FILASS/IOUT,INP
      DATA B/-0.577191652D0,0.988205891D0,-0.897056937D0,0.918206857D0,
     1-0.756704078D0,0.482199394D0,-0.193527818D0,0.035868343D0/
C**************************************************************
500   FORMAT(5X,'S NEGATIVE - MORSE FAILED',/)
C**************************************************************
      BIG=1.D+75
      CORIG=C
      IF(C.LT.0.D0)THEN
        WRITE(IOUT,500)
        STOP 'ERROR IN FACTL - C NEGATIVE'
      END IF
      FACTL=1.D0
1000  IF(C.LE.1.0D0)GO TO 2000
      BIGX=BIG/C
      IF(FACTL.GT.BIGX)GO TO 3000
      FACTL=FACTL*C
      C=C-1.0
      GO TO 1000
2000  FACT=1.D0
      X=1.D0
      DO I=1,8
        X=X*C
        FACT=FACT+B(I)*X
      END DO
      FACTL=FACTL*FACT
      RETURN
3000  IFAIL=1
      WRITE(IOUT,*)'FACTL - C',CORIG
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE HEG0(NB,M,MODE,NVB,MVB,NVF,MVF,NMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      DIMENSION NVF(NMODE),MVF(NMODE)
      COMMON/PATH/ISCFCI
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM
      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7
      COMMON/PRINT/IPRINT,JPRINT
201   FORMAT(//,1X,50(1H*),/,
     11X,'HEG INTEGRATION POINTS FOR MODE ',I2,/)
      IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0)THEN
        WRITE(INP4)NVF(MODE),NB
      END IF
C**DETERMINE NEW NUMBER OF POINTS FROM NB,M,NVF
CCCC  MVB=NVF(MODE)+M-NB
      MVB=MVF(MODE)
C     IF(MOLPRO.EQ.0)THEN
C**MAKE IT EVEN FOR STANDARD RUN
        IF(MOD(MVB,2).NE.0)MVB=MVB+1
C     ELSE
C**MAKE IT ODD FOR MOLPRO RUN
C       IF(MOD(MVB,2).EQ.0)MVB=MVB+1
C     END IF
C**CAN'T BE BIGGER THAT NUMBER OF PRIMITIVES
C**IF SO, RETAIN M GAUSS POINTS
      IF(MVB.GT.NB)MVB=M
C**THERE ARE M ORIGINAL GAUSS POINTS
C**THERE ARE NB ORIGINAL PRIMITIVES
C**THERE ARE MVB CONTRACTED FUNCTIONS (AND HEG POINTS)
C**WE USE NVB CONTRACTED FUNCTIONS
      NVB=NVF(MODE)
C     IF(IPRINT.GT.0)WRITE(IOUT,201)MODE
      IF(IPRINT.GT.0)THEN
        WRITE(IOUT,*)' NO. GAUSS POINTS',M
        WRITE(IOUT,*)' NO. PRIMITIVES',NB
        WRITE(IOUT,*)' NO. HEG POINTS',MVB
        WRITE(IOUT,*)' NO. CONTRACTED FUNCTIONS USED',NVB
        WRITE(IOUT,*)
        IF(MVB.EQ.M)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)' *********RETAIN GAUSS INTEGRATION*********'
        END IF
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE HEG1(LL,NB,M,H,XQ,XJ,OMEGA,XQJ,MODE,HX,HY,NVB,
     1MVB,XP,XV,XQX,MBF,NBF,NMODE,YK,AN,NATOM,QQ,RR,XX,X0,XL,XM,NPOT,
     2IPOT,JPOT,CPOT,XTANPM,WK,VINF,NP1,CP1,IP1,VP1,DP1,MMX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ABINIT,LINEAR
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SYMBAD
      DIMENSION NP1(1),CP1(MMX1,1),IP1(MMX1,1),VP1(MMX1,1),DP1(MMX1,1)
      DIMENSION DUMMY(NMODE),XP(NVB,MVB)
      DIMENSION H(NB,M,3,1),XQ(M),XJ(M),XV(M),MBF(NMODE),NBF(NMODE)
      DIMENSION HX(NB,M,3),HY(NVB,MVB,3,1)
      DIMENSION XQX(NB),XTANPM(NMODE)
      DIMENSION XQJ(NB,NB),YK(NB,NB),AN(NB),WK(2,NMODE,2)
      DIMENSION RR(NATOM,NATOM),QQ(NMODE),VINF(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/REACTL/JREACT,IFITRP
      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7
      COMMON/VMIN/VMIN
      COMMON/WHICH/IWHICH
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/PATH/ISCFCI
      COMMON/ABINIT/ABINIT
      COMMON/TYPE/LINEAR
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/CYCLE/ICYCLE
      COMMON/CIDIAG/ICID,ICI,JCI
C**TEMPORARY (DIMENSIONS)
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
C**************************************************************
200   FORMAT(1X,'IFAIL IN F04AAF = ',I3)
201   FORMAT(//,1X,'HEG INTEGRATION POINTS FOR MODE ',I2,/)
202   FORMAT(13X,'Q(A0)',11X,'POTENTIAL(CM-1)',/)
203   FORMAT(2X,F20.12,2X,F15.2)
204   FORMAT(5X,'S =',F10.3,' IS TOO BIG - MORSE FAILED',/)
205   FORMAT(//,1X,'ENERGY REFERENCE (V0) FOR X = V/(V-V0) = ',
     1F12.4,' CM-1',/)
206   FORMAT(13X,'Q(A0)',7X,'TANH(GAMMA.Q)',4X,'POTENTIAL(CM-1)',4X,/)
C    1'   V/(V-V0)',/)
207   FORMAT(2X,F20.12,2X,F10.4,8X,F15.2)
C     8X,F10.4)
208   FORMAT(1X,'GAUSS EXTREMES',/,
     11X,'Q(1) = ',F20.12,' Q(M) = ',F20.12)
209   FORMAT(1X,'TANH(GAMMA.Q(1)) = ',F8.2,' TANH(GAMMA.Q(M)) = ',F8.2)
210   FORMAT(1X,'V(1) = ',F20.2,' V(M) = ',F20.2,/)
C*****************************************************************
      PIQ=1.33133536380038D0
      YLAM=1/SQRT(OMEGA)
      IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0.AND.LDUMP.EQ.0)
     1WRITE(INP4)YLAM,MVB
C**TEMPORARY-LIN
      IF(MVB.NE.M.AND.MODE.LE.NONLIN)THEN
        WRITE(IOUT,*)'FORM HEG CONTRACTED FUNCTIONS FOR L = ',LL
        IF(LL.EQ.0)THEN
C**TEMPORARY-LIN
C**GET HEG POINTS IF SUFFICIENT PRIMITIVES
          DO J=1,MVB
            DO I=1,MVB
              XQJ(I,J)=0
            END DO
          END DO
C**GET EXPECTATION VALUE OF Q FOR MODE 'MODE' WITH GAUSS POINTS
          DO MM=1,M
            Q=XQ(MM)/YLAM
            DO IX=1,MVB
              X=H(IX,MM,1,1)*Q
              DO IY=1,MVB
                XQJ(IY,IX)=XQJ(IY,IX)+H(IY,MM,1,1)*X
              END DO
            END DO
          END DO
          CALL DIAG(XQJ,XQJ,NB,MVB,-1,XJ,XQX,XJ,MVB,MVB,XQJ,MVB,MVB,
     1    XQJ,XQJ,XQJ,IDUM,IDUM,IDUM)
C**TEMPORARY-LIN
        END IF
CC      IF(IPRINT.GT.0)WRITE(IOUT,201)MODE
        IF(IPRINT.GT.0)THEN
          IF(MOLPRO.EQ.0.AND.IWHICH.GE.0)THEN
            WRITE(IOUT,202)
          END IF
        END IF
C**TEMPORARY-LIN

CC      IF(MODE.LE.NONLIN)THEN
C**NON-LINEAR MODES
          DO J=1,MVB
            Q=XQX(J)
            XQX(J)=Q*YLAM
C**ENERGY THIS POINT
            DO K=1,NMODE
              QQ(K)=0
            END DO
            QQ(MODE)=XQX(J)
            IF(IWHICH.GT.0)THEN
              DO I=1,NATOM
                DO K=1,3
                  XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQ(MODE)/
     1            SQRT(XM(I))
                END DO
              END DO
              IF(ABINIT)THEN
                CALL GETAPT(VPT,NATOM,XX,RR,QQ,NMODE,0)
              ELSE
                CALL GETPOT(VPT,NATOM,XX,RR)
              END IF
            ELSE
              IF(IWHICH.EQ.0)THEN
                CALL GETPT1(VPT,NPOT,IPOT,JPOT,CPOT,0,QQ)
              ELSE
                IF(MOLINC.GT.0)
     1          CALL GETPQT(VPT,NMODE,QQ,XTANPM,NP1,CP1,IP1,NMODE,MMX1)
                IF(MOLINC.EQ.0)
     1          CALL GETQPT(VPT,NMODE,QQ,XTANPM)
                IF(MOLINC.LT.0)
     1          CALL GETQP1(VPT,DUMMY,NMODE,QQ,NP1,CP1,VP1,DP1,NMODE,
     2          MMX1)
              END IF
            END IF
            XV(J)=VPT
            IF(IPRINT.GT.0)THEN
              IF(MOLPRO.EQ.0.AND.IWHICH.GE.0)THEN
C               WRITE(IOUT,203)XQX(J),(VPT-VMIN)*WAVENM
                WRITE(IOUT,203)XQX(J),VPT*WAVENM
              END IF
            END IF
            IF(VPT.LT.VMIN)VMIN=VPT
C**ENERGY THIS POINT
            IF(J.EQ.1)VMAX=VPT
            IF(J.EQ.MVB.AND.VPT.LT.VMAX)VMAX=VPT
          END DO
          IF(MOLPRO.NE.0.OR.IWHICH.LT.0.OR.IFITRP.NE.0)THEN
            VINF(MODE)=VINF(MODE)+VMAX/10
            IF(IPRINT.GT.0)THEN
C             WRITE(IOUT,205)VINF(MODE)*WAVENM
              WRITE(IOUT,208)WK(1,MODE,1),WK(2,MODE,1)
              WRITE(IOUT,209)XTANH(XTANPM(MODE)*WK(1,MODE,1)),
     1                       XTANH(XTANPM(MODE)*WK(2,MODE,1))
              WRITE(IOUT,210)WAVENM*WK(1,MODE,2),WAVENM*WK(2,MODE,2)
              WRITE(IOUT,206)
            END IF
            DO J=1,MVB
C**ENERGY THIS POINT
              DO K=1,NMODE
                QQ(K)=0
              END DO
              QQ(MODE)=XQX(J)
              IF(IWHICH.GT.0)THEN
                DO I=1,NATOM
                  DO K=1,3
                    XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQ(MODE)/
     1              SQRT(XM(I))
                  END DO
                END DO
                IF(ABINIT)THEN
                  CALL GETAPT(VPT,NATOM,XX,RR,QQ,NMODE,0)
                ELSE
                  CALL GETPOT(VPT,NATOM,XX,RR)
                END IF
              ELSE
                IF(IWHICH.EQ.0)THEN
                  CALL GETPT1(VPT,NPOT,IPOT,JPOT,CPOT,0,QQ)
                ELSE
                  IF(MOLINC.GT.0)
     1            CALL GETPQT(VPT,NMODE,QQ,XTANPM,NP1,CP1,IP1,NMODE,
     2            MMX1)
                  IF(MOLINC.EQ.0)
     1            CALL GETQPT(VPT,NMODE,QQ,XTANPM)
                  IF(MOLINC.LT.0)
     1            CALL GETQP1(VPT,DUMMY,NMODE,QQ,NP1,CP1,VP1,DP1,NMODE,
     2            MMX1)
                END IF
              END IF
              XV(J)=VPT
              IF(IPRINT.GT.0)THEN
                IF(MOLPRO.NE.0.OR.IWHICH.LT.0.OR.IFITRP.NE.0)THEN
                  WRITE(IOUT,207)XQX(J),XTANH(XTANPM(MODE)*XQX(J)),
     1            (VPT-VMIN)*WAVENM
C    2            ,(VPT-VMIN)/(VPT-VMIN-VINF(MODE))
                END IF
              END IF
C**ENERGY THIS POINT
            END DO
          END IF
C**GET PRIMITIVES (NB PRIMITIVES) AT HEG POINTS (MVB POINTS)
          DO 100 J=1,MVB
          Q=(XQX(J))/YLAM
          HX(1,J,1)=1.D0/PIQ
          HX(2,J,1)=2.D0*Q/(PIQ*SQRT(2.D0))
          IF(NB.GE.3)THEN
          DO 70 I=3,NB
          I1=I-1
          I2=I-2
70        HX(I,J,1)=2*(Q*HX(I1,J,1)-I2*HX(I2,J,1)/SQRT(2.D0*I2))/
     1    SQRT(2.D0*I1)
          END IF
          HX(1,J,2)=0.D0
          HX(2,J,2)=2.D0/(PIQ*SQRT(2.D0))
          HX(1,J,3)=0.D0
          HX(2,J,3)=0.D0
          IF(NB.GE.3)THEN
          DO 80 I=3,NB
          I1=I-1
          I2=I-2
          HX(I,J,2)=2.D0*I1*HX(I1,J,1)/SQRT(2.D0*I1)
          HX(I,J,3)=4.D0*I1*I2*HX(I2,J,1)/SQRT(4.D0*I1*I2)
80        CONTINUE
          END IF
          DO 90 I=1,NB
          HX(I,J,3)=HX(I,J,3)-2*Q*HX(I,J,2)+(Q*Q-1.D0)
     1    *HX(I,J,1)
90        HX(I,J,2)=HX(I,J,2)-Q*HX(I,J,1)
100       CONTINUE
          DO 110 J=1,MVB
          DO 110 I=1,NB
          HX(I,J,2)=HX(I,J,2)/YLAM
          HX(I,J,3)=HX(I,J,3)/(YLAM*YLAM)
110       CONTINUE
C**FORM NEW CONTRACTED FUNCTIONS AT HEG PTS IN OLD STORAGE LOCATIONS
          DO K=1,3
            DO J=1,MVB
              DO I=1,NVB
                HY(I,J,K,1)=0
                DO L=1,NB
                  HY(I,J,K,1)=HY(I,J,K,1)+YK(L,I)*HX(L,J,K)
                END DO
              END DO
            END DO
          END DO
C**WEIGHTS FOR FUNCTION I (FUNCTION*SQRT(WT) IN XQJ)
          DO I=1,NVB
            DO J=1,MVB
C**XJ IS SQRT(WEIGHT)
              XJ(J)=XQJ(I,J)/HY(I,J,1,1)
              Q=(XQX(J))/YLAM
              XP(I,J)=XJ(J)/EXP(-Q*Q/2)
            END DO
            DO K=1,3
              DO J=1,MVB
                HY(I,J,K,1)=HY(I,J,K,1)*XJ(J)
              END DO
            END DO
          END DO
C**TEMPORARY-LIN
CC      ELSE
C**LINEAR MODES
CC      END IF
      ELSE
C**RETAIN GAUSS POINTS AND GAUSS CONTRACTED FUNCTIONS
        WRITE(IOUT,*)
        WRITE(IOUT,*)' *********RETAIN GAUSS INTEGRATION*********'
C       IJK=1+LL
        IJK=1
        DO K=1,3
          DO J=1,M
            DO I=1,NVB
              IF(MODE.GT.NONLIN)HX(I,J,K)=H(I,J,K,IJK)
              IF(MODE.LE.NONLIN)HY(I,J,K,IJK)=H(I,J,K,IJK)
            END DO
          END DO
        END DO
        IF(.NOT.LINEAR.OR.LL.EQ.JMAX)THEN
          DO J=1,M
            XQX(J)=XQ(J)
          END DO
        END IF
C**TEMPORARY-LIN
      END IF

      IF(MODE.GT.NONLIN)THEN
        DO K=1,3
          DO J=1,M
            WRITE(70)(HX(I,J,K),I=1,NVB)
          END DO
          DO J=1,M
            DO I=1,NB
              H(I,J,K,1)=H(I,J,K,2)
            END DO
          END DO
        END DO
      END IF
C**TEMPORARY-LIN
      IF(LL.EQ.0)THEN
C**TEMPORARY-LIN
        IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0.AND.LDUMP.EQ.0)THEN
          WRITE(INP4)(XQX(J),J=1,MVB)
        END IF
C**TEMPORARY-LIN
      END IF
C**TEMPORARY-LIN
C**TEMPORARY (TEST ORTHOGONALITY)
CCCC  IDUM=1
CCCC  IF(MODE.GT.NONLIN)IDUM=2
CCCC  DO IJK=1,IDUM
CCCC    DO J=1,NVB
CCCC      DO I=1,NVB
CCCC        XQJ(I,J)=0
CCCC      END DO
CCCC    END DO
CCCC    DO MM=1,MVB
CCCC      DO IX=1,NVB
CCCC        X=HY(IX,MM,1,IJK)
CCCC        DO IY=1,NVB
CCCC          XQJ(IY,IX)=XQJ(IY,IX)+HY(IY,MM,1,IJK)*X
CCCC        END DO
CCCC      END DO
CCCC    END DO
CCCC    DO IX=1,NVB
CCCC      WRITE(IOUT,*)(XQJ(IY,IX),IY=1,NVB)
CCCC    END DO
CCCC  END DO
C**TEMPORARY (TEST ORTHOGONALITY)
C**CHANGE NUMBER OF POINTS ARRAY
      IF(MVB.NE.M.AND.MODE.LE.NONLIN)THEN
        MBF(MODE)=MVB
      ELSE
        MBF(MODE)=M
      END IF
C**CHANGE NUMBER OF FUNCTIONS ARRAY
      NBF(MODE)=NVB
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE HEG2(MAXBAS,NMODE,J,NVF,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NMODE,J),NVF(NMODE)
      NVF(MODE)=0
      DO K=1,J
        IF(MAXBAS(MODE,K).GT.NVF(MODE))NVF(MODE)=MAXBAS(MODE,K)
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRV0(MMTAU,VC,VCR,VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP,VC(3),VR(J21)
      REAL*4 VPR,VCR(3),VRR(J21)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      MDT=MODINT(NSMODE)
      DO MTAU=1,MMTAU/MDT
        IF(JCOUPL.GE.0)THEN
          READ(71)VP
        ELSE
          READ(71)VPR
        END IF
        IF(JCOUPC.GE.0)THEN
          IF(J21.GT.1.AND.ICOUPC.GE.0)READ(61)VR
          IF(ICOUPC.GE.0)READ(81)VC
        ELSE
          IF(J21.GT.1.AND.ICOUPC.GE.0)READ(61)VRR
          IF(ICOUPC.GE.0)READ(81)VCR
        END IF
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRD1(MOD1,MM,VP,VPR,VC,VCR,VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM),VC(MM),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM),VRR(J21,MM)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FILASS/IOUT,INP
      MD=MODINT(MOD1)
      CALL IN6781(VR,VRR,VP,VPR,VC,VCR,MM/MD,J21)

C***********************************************************

C     IF(JCOUPL.GT.0)THEN
C       IF(J21.GT.1.AND.ICOUPL.EQ.1)READ(61)VR
C       READ(71)VP
C       READ(81)VC
C     ELSE
C       IF(J21.GT.1.AND.ICOUPL.EQ.1)READ(61)VRR
C       READ(71)VPR
C       READ(81)VCR
C     END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRM1(MM,VM,VMR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM,6)
      REAL*4 VMR(MM,6)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FILASS/IOUT,INP

C***********************************************************

      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.1)READ(91)VM
      ELSE
        IF(ICOUPC.GE.1)READ(91)VMR
      END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRV1(MOD1,MM,MMTAU,VP,VPR,VC,VCR,VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM),VC(MM,6),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM,6),VRR(J21,MM)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      MDT=MODINT(NSMODE)
      MD1=MODINT(MOD1)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(NSMODE.EQ.ISYM(I,J))NT=I
        END DO
      END DO
      IF(N1.EQ.NT)MD1=1
      CALL DDMRV1(MOD1,MM/MD1,MMTAU,VP,VPR,VC,VCR,VR,VRR,J21,MODINT)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DDMRV1(MOD1,MM,MMTAU,VP,VPR,VC,VCR,VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      REAL*8 VP(MM),VC(MM,6),VR(J21,MM)
      REAL*4 VPR(MM),VCR(MM,6),VRR(J21,MM)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/XSAVES/ICOUPX,ICOUCX,JCOUPS,JCOUCS
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/FILASS/IOUT
      ICCCCC=ICOUPC
      IF(MOLINC.EQ.0)THEN
        ICCCCC=ICOUCX
      END IF
      MDT=MODINT(NSMODE)
      DO MTAU=1,MMTAU/MDT
        IF(JCOUPL.GE.0)THEN
          READ(71)VP
        ELSE
          READ(71)VPR
        END IF
        IF(JCOUPC.GE.0)THEN
          IF(J21.GT.1.AND.ICCCCC.GE.1)READ(61)VR
          IF(ICCCCC.GT.0)THEN
            READ(81)VC
          END IF
        ELSE
          IF(J21.GT.1.AND.ICCCCC.GE.1)READ(61)VRR
          IF(ICCCCC.GT.0)READ(81)VCR
        END IF
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRD2(MOD1,MOD2,MM1,MM2,VP,VPR,VC,VCR,VR,VRR,J21,
     1MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,6),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,6),VRR(J21,MM2,MM1)
      DIMENSION MODINT(1)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
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
      CALL IN6782(VR,VRR,VP,VPR,VC,VCR,MM1/MD1,MM2/MD2,J21)

C***********************************************************

C     IF(JCOUPL.GT.0)THEN
C       IF(J21.GT.1.AND.ICOUPL.EQ.2)READ(62)VR
C       READ(72)VP
C       READ(82)VC
C     ELSE
C       IF(J21.GT.1.AND.ICOUPL.EQ.2)READ(62)VRR
C       READ(72)VPR
C       READ(82)VCR
C     END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRM2(MM1,MM2,VM,VMR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM2,MM1,12)
      REAL*4 VMR(MM2,MM1,12)
C**TEMPORARY (DIMENSIONS)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FILASS/IOUT,INP

C***********************************************************

      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.2)READ(92)VM
      ELSE
        IF(ICOUPC.GE.2)READ(92)VMR
      END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRV2(MOD1,MOD2,MM1,MM2,MMTAU,VP,VPR,VC,VCR,VR,VRR,
     1J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,10),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,10),VRR(J21,MM2,MM1)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
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
      IF(N1.EQ.NT)MD1=1
      IF(N2.EQ.NT)MD2=1
      N1T=ISYMP(N1,NT)
      IF(N1T.EQ.N2)MD2=1
      CALL DDMRV2(MOD1,MOD2,MM1/MD1,MM2/MD2,MMTAU,VP,VPR,VC,VCR,VR,VRR,
     1J21,MODINT)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DDMRV2(MOD1,MOD2,MM1,MM2,MMTAU,VP,VPR,VC,VCR,VR,VRR,
     1J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      REAL*8 VP(MM2,MM1),VC(MM2,MM1,10),VR(J21,MM2,MM1)
      REAL*4 VPR(MM2,MM1),VCR(MM2,MM1,10),VRR(J21,MM2,MM1)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/XSAVES/ICOUPX,ICOUCX,JCOUPS,JCOUCS
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/FILASS/IOUT
      ICCCCC=ICOUPC
      IF(MOLINC.EQ.0)THEN
        ICCCCC=ICOUCX
      END IF
      MDT=MODINT(NSMODE)
      DO MTAU=1,MMTAU/MDT
        IF(JCOUPL.GT.0)THEN
          READ(72)VP
        ELSE
          READ(72)VPR
        END IF
        IF(JCOUPC.GE.0)THEN
          IF(J21.GT.1.AND.ICCCCC.GE.2)READ(62)VR
          IF(ICCCCC.GE.2)THEN
            READ(82)VC
          END IF
        ELSE
          IF(J21.GT.1.AND.ICCCCC.GE.2)READ(62)VRR
          IF(ICCCCC.GE.2)READ(82)VCR
        END IF
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRD3(MOD1,MOD2,MOD3,MM1,MM2,MM3,VP,VPR,VC,VCR,VR,
     1VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,10),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,10),VRR(J21,MM3,MM2,MM1)
      DIMENSION MODINT(1)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
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
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N2.EQ.N3)MD3=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      CALL IN6783(VR,VRR,VP,VPR,VC,VCR,MM1/MD1,MM2/MD2,MM3/MD3,J21)

C***********************************************************

C     IF(JCOUPL.GT.0)THEN
C       IF(J21.GT.1.AND.ICOUPL.EQ.3)READ(63)VR
C       READ(73)VP
C       READ(83)VC
C     ELSE
C       IF(J21.GT.1.AND.ICOUPL.EQ.3)READ(63)VRR
C       READ(73)VPR
C       READ(83)VCR
C     END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRM3(MM1,MM2,MM3,VM,VMR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM3,MM2,MM1,15)
      REAL*4 VMR(MM3,MM2,MM1,15)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FILASS/IOUT,INP

C***********************************************************

      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.3)READ(93)VM
      ELSE
        IF(ICOUPC.GE.3)READ(93)VMR
      END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRV3(MOD1,MOD2,MOD3,MM1,MM2,MM3,MMTAU,VP,VPR,VC,VCR,
     1VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,15),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,15),VRR(J21,MM3,MM2,MM1)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
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
      IF(N1.EQ.NT)MD1=1
      IF(N2.EQ.NT)MD2=1
      IF(N3.EQ.NT)MD3=1
      CALL DDMRV3(MOD1,MOD2,MOD3,MM1/MD1,MM2/MD2,MM3/MD3,MMTAU,VP,VPR,
     1VC,VCR,VR,VRR,J21,MODINT)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DDMRV3(MOD1,MOD2,MOD3,MM1,MM2,MM3,MMTAU,VP,VPR,VC,VCR,
     1VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      REAL*8 VP(MM3,MM2,MM1),VC(MM3,MM2,MM1,15),VR(J21,MM3,MM2,MM1)
      REAL*4 VPR(MM3,MM2,MM1),VCR(MM3,MM2,MM1,15),VRR(J21,MM3,MM2,MM1)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/XSAVES/ICOUPX,ICOUCX,JCOUPS,JCOUCS
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      ICCCCC=ICOUPC
      IF(MOLINC.EQ.0)THEN
        ICCCCC=ICOUCX
      END IF
      MDT=MODINT(NSMODE)
      DO MTAU=1,MMTAU/MDT
        IF(JCOUPL.GT.0)THEN
          READ(73)VP
        ELSE
          READ(73)VPR
        END IF
        IF(JCOUPC.GT.0)THEN
          IF(J21.GT.1.AND.ICCCCC.GE.3)READ(63)VR
          IF(ICCCCC.GE.3)READ(83)VC
        ELSE
          IF(J21.GT.1.AND.ICCCCC.GE.3)READ(63)VRR
          IF(ICCCCC.GE.3)READ(83)VCR
        END IF
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRD4(MOD1,MOD2,MOD3,MOD4,MM1,MM2,MM3,MM4,VP,VPR,VC,
     1VCR,VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,15),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,15),
     1VRR(J21,MM4,MM3,MM2,MM1)
      DIMENSION MODINT(1)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
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
      CALL IN6784(VR,VRR,VP,VPR,VC,VCR,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,
     1J21)

C***********************************************************

C     IF(JCOUPL.GT.0)THEN
C       IF(J21.GT.1.AND.ICOUPL.EQ.4)READ(64)VR
C       READ(74)VP
C       READ(84)VC
C     ELSE
C       IF(J21.GT.1.AND.ICOUPL.EQ.4)READ(64)VRR
C       READ(74)VPR
C       READ(84)VCR
C     END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRM4(MM1,MM2,MM3,MM4,VM,VMR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM4,MM3,MM2,MM1,18)
      REAL*4 VMR(MM4,MM3,MM2,MM1,18)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FILASS/IOUT

C***********************************************************

      IF(JCOUPC.GE.0)THEN
        IF(ICOUPC.GE.4)READ(94)VM
      ELSE
        IF(ICOUPC.GE.4)READ(94)VMR
      END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRV4(MOD1,MOD2,MOD3,MOD4,MM1,MM2,MM3,MM4,MMTAU,VP,
     1VPR,VC,VCR,VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,21),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,21),
     1VRR(J21,MM4,MM3,MM2,MM1)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
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
      IF(N1.EQ.NT)MD1=1
      IF(N2.EQ.NT)MD2=1
      IF(N3.EQ.NT)MD3=1
      IF(N4.EQ.NT)MD4=1
      CALL DDMRV4(MOD1,MOD2,MOD3,MOD4,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,
     1MMTAU,VP,VPR,VC,VCR,VR,VRR,J21,MODINT)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DDMRV4(MOD1,MOD2,MOD3,MOD4,MM1,MM2,MM3,MM4,MMTAU,VP,
     1VPR,VC,VCR,VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      REAL*8 VP(MM4,MM3,MM2,MM1),VC(MM4,MM3,MM2,MM1,21),
     1VR(J21,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM4,MM3,MM2,MM1),VCR(MM4,MM3,MM2,MM1,21),
     1VRR(J21,MM4,MM3,MM2,MM1)
      DIMENSION MODINT(1)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/XSAVES/ICOUPX,ICOUCX,JCOUPS,JCOUCS
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      ICCCCC=ICOUPC
      IF(MOLINC.EQ.0)THEN
        ICCCCC=ICOUCX
      END IF
      MDT=MODINT(NSMODE)
      DO MTAU=1,MMTAU/MDT
        IF(JCOUPL.GT.0)THEN
          READ(74)VP
        ELSE
          READ(74)VPR
        END IF
        IF(JCOUPC.GT.0)THEN
          IF(J21.GT.1..AND.ICCCCC.GE.4)READ(64)VR
          IF(ICCCCC.GE.4)READ(84)VC
        ELSE
          IF(J21.GT.1..AND.ICCCCC.GE.4)READ(64)VRR
          IF(ICCCCC.GE.4)READ(84)VCR
        END IF
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRD5(MOD1,MOD2,MOD3,MOD4,MOD5,MM1,MM2,MM3,MM4,MM5,
     1VP,VPR,VC,VCR,VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM5,MM4,MM3,MM2,MM1),VC(MM5,MM4,MM3,MM2,16),
     1VR(J21,MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM5,MM4,MM3,MM2,MM1),VCR(MM5,MM4,MM3,MM2,16),
     1VRR(J21,MM5,MM4,MM3,MM2,MM1)
      DIMENSION MODINT(1)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/FILASS/IOUT
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
      CALL IN6785(VP,VPR,VC,VCR,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,
     1MM5/MD5)

C***********************************************************

C     IF(JCOUPL.GT.0)THEN
C       IF(J21.GT.1.AND.ICOUPL.EQ.5)READ(65)VR
C       READ(75)VP
C     ELSE
C       IF(J21.GT.1.AND.ICOUPL.EQ.5)READ(65)VRR
C       READ(75)VPR
C     END IF
C     IF(JCOUPC.GT.0)THEN
C       IF(ICOUPC.GT.4)READ(85)VC
C     ELSE
C       IF(ICOUPC.GT.4)READ(85)VCR
C     END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DUMRD6(MOD1,MOD2,MOD3,MOD4,MOD5,MOD6,MM1,MM2,MM3,MM4,
     1MM5,MM6,VP,VPR,VC,VCR,VR,VRR,J21,MODINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VP(MM6,MM5,MM4,MM3,MM2,MM1),
     1VC(MM6,MM5,MM4,MM3,22),VR(J21,MM6,MM5,MM4,MM3,MM2,MM1)
      REAL*4 VPR(MM6,MM5,MM4,MM3,MM2,MM1),
     1VCR(MM6,MM5,MM4,MM3,22),VRR(J21,MM6,MM5,MM4,MM3,MM2,MM1)
      DIMENSION MODINT(1)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      MD1=MODINT(MOD1)
      MD2=MODINT(MOD2)
      MD3=MODINT(MOD3)
      MD4=MODINT(MOD4)
      MD5=MODINT(MOD5)
      MD6=MODINT(MOD6)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MOD1.EQ.ISYM(I,J))N1=I
          IF(MOD2.EQ.ISYM(I,J))N2=I
          IF(MOD3.EQ.ISYM(I,J))N3=I
          IF(MOD4.EQ.ISYM(I,J))N4=I
          IF(MOD5.EQ.ISYM(I,J))N5=I
          IF(MOD6.EQ.ISYM(I,J))N6=I
        END DO
      END DO
      IF(N1.EQ.N2)MD2=1
      IF(N1.EQ.N3)MD3=1
      IF(N1.EQ.N4)MD4=1
      IF(N1.EQ.N5)MD5=1
      IF(N1.EQ.N6)MD6=1
      IF(N2.EQ.N3)MD3=1
      IF(N2.EQ.N4)MD4=1
      IF(N2.EQ.N5)MD5=1
      IF(N2.EQ.N6)MD6=1
      IF(N3.EQ.N4)MD4=1
      IF(N3.EQ.N5)MD5=1
      IF(N3.EQ.N6)MD6=1
      IF(N4.EQ.N5)MD5=1
      IF(N4.EQ.N6)MD6=1
      IF(N5.EQ.N6)MD6=1
      N12=ISYMP(N1,N2)
      IF(N12.EQ.N3)MD3=1
      IF(N12.EQ.N4)MD4=1
      IF(N12.EQ.N5)MD5=1
      IF(N12.EQ.N6)MD6=1
      N13=ISYMP(N1,N3)
      IF(N13.EQ.N4)MD4=1
      IF(N13.EQ.N5)MD5=1
      IF(N13.EQ.N6)MD6=1
      N14=ISYMP(N1,N4)
      IF(N14.EQ.N5)MD5=1
      IF(N14.EQ.N6)MD6=1
      N15=ISYMP(N1,N5)
      IF(N15.EQ.N6)MD6=1
      N23=ISYMP(N2,N3)
      IF(N23.EQ.N4)MD4=1
      IF(N23.EQ.N5)MD5=1
      IF(N23.EQ.N6)MD6=1
      N24=ISYMP(N2,N4)
      IF(N24.EQ.N5)MD5=1
      IF(N24.EQ.N6)MD6=1
      N25=ISYMP(N2,N5)
      IF(N25.EQ.N6)MD6=1
      N34=ISYMP(N3,N4)
      IF(N34.EQ.N5)MD5=1
      IF(N34.EQ.N6)MD6=1
      N35=ISYMP(N3,N5)
      IF(N35.EQ.N6)MD6=1
      N123=ISYMP(N12,N3)
      IF(N123.EQ.N4)MD4=1
      IF(N123.EQ.N5)MD5=1
      IF(N123.EQ.N6)MD6=1
      N124=ISYMP(N12,N4)
      IF(N124.EQ.N5)MD5=1
      IF(N124.EQ.N6)MD6=1
      N125=ISYMP(N12,N5)
      IF(N125.EQ.N6)MD6=1
      N134=ISYMP(N13,N4)
      IF(N134.EQ.N5)MD5=1
      IF(N134.EQ.N6)MD6=1
      N135=ISYMP(N13,N5)
      IF(N135.EQ.N6)MD6=1
      N145=ISYMP(N14,N5)
      IF(N145.EQ.N6)MD6=1
      N234=ISYMP(N23,N4)
      IF(N234.EQ.N5)MD5=1
      IF(N234.EQ.N6)MD6=1
      N235=ISYMP(N23,N5)
      IF(N235.EQ.N6)MD6=1
      N245=ISYMP(N24,N5)
      IF(N245.EQ.N6)MD6=1
      N345=ISYMP(N34,N5)
      IF(N345.EQ.N6)MD6=1
      N1234=ISYMP(N12,N34)
      IF(N1234.EQ.N5)MD5=1
      IF(N1234.EQ.N6)MD6=1
      N1235=ISYMP(N12,N35)
      IF(N1235.EQ.N6)MD6=1
      N1245=ISYMP(N12,N45)
      IF(N1245.EQ.N6)MD6=1
      N1345=ISYMP(N13,N45)
      IF(N1345.EQ.N6)MD6=1
      N2345=ISYMP(N23,N45)
      IF(N2345.EQ.N6)MD6=1
      N12345=ISYMP(N123,N45)
      IF(N12345.EQ.N6)MD6=1
      CALL IN6786(VP,VPR,VC,VCR,MM1/MD1,MM2/MD2,MM3/MD3,MM4/MD4,MM5/MD5,
     1MM6/MD6)

C***********************************************************

C     IF(JCOUPL.GT.0)THEN
C       IF(J21.GT.1.AND.ICOUPL.EQ.6)READ(66)VR
C       READ(76)VP
C     ELSE
C       IF(J21.GT.1.AND.ICOUPL.EQ.6)READ(66)VRR
C       READ(76)VPR
C     END IF
C     IF(JCOUPC.GT.0)THEN
C       IF(ICOUPC.GT.5)READ(86)VC
C     ELSE
C       IF(ICOUPC.GT.5)READ(86)VCR
C     END IF

C***********************************************************

      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DISCDP(W,KCONT,LCONT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100),SPACE
      LOGICAL LINEAR
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
     1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
     4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
     5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
     6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
     7LDP4A,LDP4B,LDP4C,LDP4D,LXP,LJP6,LANCNT,LANK,LXJ,LVP5,
     8LC5,LC6,LVP0,LTMPRD,LMVB,LNFC,LTEMP6,LTEMP7,LTEMP8,LTEMP9,
     9LXMODQ
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/FILASS/IOUT,INP
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/DISC/IDISC
      COMMON/MODES/NMODE,NATOM
      COMMON/NCPOT/NPOT
      COMMON/WHICH/IWHICH
      COMMON/ENTER/IENTER,IENTMX(5),NTOT1,NTOT2,NTOT3,NTOT4,NTOT5
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/REACTL/JREACT,IFITRP
C**TEMPORARY (DIMENSIONS)
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/CSAVES/NNMODS,NAMODS,NVMODS,ICOUPS,ICOUCS,JREACS,NREACS
      COMMON/XSAVES/ICOUPX,ICOUCX,JCOUPS,JCOUCS
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/MAXPT/MBFMAX,MBFMX1,MBFMX2,MBFMX3,MBFMX4,MBFMIN
      COMMON/DISCSZ/KLC0,KEJK0,KLV1,KLC1,KEJK1,KLV2,KLC2,KEJK2,
     1KLV3,KLC3,KEJK3,KLV4,KLC4,KEJK4
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
C**TEMPORARY
      COMMON/SIZEJ/JSIZE1(2),JSIZE2(2),JSIZE3(2),JSIZE4(2),JSIZE5(2),
     1JSIZE6(2)
      COMMON/TOTALS/ITOT1(2),ITOT2(2),ITOT3(2),ITOT4(2),ITOT5(2),
     1ITOT6(2),ITOT
C**TEMPORARY
C************************
      IF(MOLINC.EQ.0)THEN
        ICOUPL=ICOUPX
        ICOUPC=ICOUCX
      END IF
C**TEMPORARY
      WRITE(IOUT,*)'DISCDP - KCONT,LCONT,JREACT = ',KCONT,LCONT,JREACT
      WRITE(IOUT,*)'ICOUPL,ICOUPC,NCOUPL(1),NCOUPC(1),NCOUPL(2),
     1NCOUPC(2)'
      WRITE(IOUT,*)ICOUPL,ICOUPC,NCOUPL(1),NCOUPC(1),NCOUPL(2),
     1NCOUPC(2)
      WRITE(IOUT,*)'JCOUPL,JCOUPC,ICOUPX,ICOUCX,ICOUPS,ICOUCS'
      WRITE(IOUT,*)JCOUPL,JCOUPC,ICOUPX,ICOUCX,ICOUPS,ICOUCS
      WRITE(IOUT,*)'NMODE,NAMODE,NSMODE,NVMODE,NNMODE = '
      WRITE(IOUT,*)NMODE,NAMODE,NSMODE,NVMODE,NNMODE
      CALL FLUSH(IOUT)
C**TEMPORARY
C**JREACT DIFFERS BETWEEN CONTRACTION SCHEMES
      KEJK0=J21
      IF(LCONT.EQ.0)
     1CALL MEMO(1,LEJK0,KEJK0,0,0,0,0,0,0,0,0)
      IF(JREACT.GT.0)THEN
        CALL INTARR(W(LNBF),W(LMBF),NSMODE,IK1TAU,IK2TAU,IK3TAU)
        K3TAU=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          K3TAU=K3TAU+IK2
        END DO
        KLC0=3
        IF(LCONT.EQ.0)
     1  CALL MEMO(1,LC0,KLC0,0,0,0,0,0,0,0,0)
        IF(IDISC.EQ.0)THEN
C**GET TORSION-ONLY GRID DATA
          CALL MEMO(1,LVP0,IK2TAU,0,0,0,0,0,0,0,0)
          CALL DUMVT0(NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     1    W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     2    W(LXQ+K3TAU),W(LVP0),1,W(LMODNT),W(LXMODQ))
          CALL MEMO(-1,LVP0,IK2TAU,0,0,0,0,0,0,0,0)
          CALL MEMO(1,LVP0,3*IK2TAU,0,0,0,0,0,0,0,0)
          CALL DUMVR0(NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),
     1    W(LB),W(LB+NMODE*NMODE),W(LAA),W(LBB),W(LBB+NMODE),W(LXX),
     2    W(LXX+NATOM*3*722),W(LX0),W(LXL),W(LXM),W(LC0),W(LC0),
     3    W(LVM0),W(LVM0),0,W(LVP0),1,W(LMODNT))
          CALL MEMO(-1,LVP0,3*IK2TAU,0,0,0,0,0,0,0,0)
C         IF(JMAX.GT.0.AND.ICOUPL.EQ.0)
          IF(JMAX.GT.0)
     1    CALL GETMV0(NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),
     2    W(LB),W(LB+NMODE*NMODE),W(LAA),W(LBB),W(LBB+NMODE),W(LXX),
     3    W(LXX+NATOM*3*722),W(LX0),W(LXL),W(LXM),W(LEJK0),W(LEJK0),
     4    W(LSS),J21,W(LX21),W(LW21),W(LE21),W(LMODNT))
        END IF
      END IF
C******************************
      IF(LDUMP.EQ.0.AND.(JREACT.LE.0.OR.NSMODE.EQ.NAMODE))THEN
        CALL DUMPT0(NMODE,NATOM,W(LQQ),W(LRR),W(LX0),W(LXMODQ))
        CALL DUMCR0(NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),
     1  W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),0,1)
      END IF
      IF(LDUMP.GT.0.AND.(JREACT.LE.0.OR.NSMODE.EQ.NAMODE))
     1CALL DUMDP0(NMODE,NATOM,W(LQQ),W(LX0),NPOT,W(LIPOT),W(LJPOT),
     2W(LCPOT),K,LDUMP)
C**GET EQUILIBRIUM ROTATIONAL CONSTANTS
      IF(.NOT.LINEAR.AND.LDUMP.EQ.0.AND.JMAX.GT.0.AND.(JREACT.LE.0.OR.
     1NSMODE.EQ.NAMODE))
     2CALL GETMI0(NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),
     3W(LXX),W(LX0),W(LXL),W(LXM),W(LEJK0),W(LEJK0),W(LSS),J21,W(LX21),
     4W(LW21),W(LE21))
      IF(ICOUPL.EQ.0)GO TO 9999
C**INTRINSIC
      IF(KCONT.LT.0)THEN
        JCOUPS=JCOUPL
        JCOUCS=JCOUPC
C       ICOUPX=ICOUPL
C       ICOUCX=ICOUPC
      END IF
C******************************
      MBFMX1=0
C**NAMODE DIFFERS BETWEEN CONTRACTION SCHEMES
      DO KK=1,NAMODE
        IF(LCONT.NE.0)THEN
          K=JCONT(LCONT,KK)
        ELSE
          K=KK
        END IF
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(IK2.GT.MBFMX1)MBFMX1=IK2
      END DO
      MBFMAX=MBFMX1
      IF(JCOUPL.GT.0)THEN
        KLV1=MBFMX1
      ELSE
        KLV1=(1+MBFMX1)/2
      END IF
C**INTRINSIC
      KLVR1=KLV1
      IF(JREACT.NE.0.AND.MOLINC.GT.0.AND.ICOUPX.GE.2)KLVR1=KLV1*IK2TAU
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.1)THEN
          KLC1=MBFMX1
          IF(JREACT.NE.0)KLC1=MBFMX1*6
        END IF
      ELSE
        IF(ICOUCX.GE.1)THEN
          KLC1=(1+MBFMX1)/2
          IF(JREACT.NE.0)KLC1=(1+MBFMX1*6)/2
        END IF
      END IF
C**INTRINSIC
      KLCR1=KLC1
      IF(JREACT.NE.0.AND.MOLINC.GT.0.AND.ICOUCX.GE.2)KLCR1=KLC1*IK2TAU
      IF(LCONT.EQ.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1  CALL MEMO(1,LV1,KLVR1,0,0,0,0,0,0,0,0)
        IF(ICOUCX.GE.1)CALL MEMO(1,LC1,KLCR1,0,0,0,0,0,0,0,0)
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.1)KEJK1=MBFMX1*J21
      ELSE
        IF(ICOUCX.GE.1)KEJK1=(1+MBFMX1*J21)/2
      END IF
      IF(LCONT.EQ.0)THEN
      IF(ICOUCX.GE.1)CALL MEMO(1,LEJK1,KEJK1,0,0,0,0,0,0,0,0)
      END IF
      IF(IDISC.EQ.0)THEN
        K=0
        K3=0
        DO KK=1,NAMODE
          IF(LCONT.NE.0)THEN
            KNEXT=JCONT(LCONT,KK)
          ELSE
            KNEXT=KK
          END IF
          IF(KNEXT-K.GT.1)THEN
C**UPDATE K3 FOR MISSING MODES
            DO KADJ=K+1,KNEXT-1
              CALL INTARR(W(LNBF),W(LMBF),KADJ,IK1,IK2,IK3)
              K3=K3+IK2
            END DO
          END IF
C**NEXT K
          K=KNEXT

          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**INTRINSIC
          IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
            CALL GROUP1(K,ICONT(1),JCONT,1,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(1)
              JCOUPC=NCOUPC(1)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 3031
            END IF
            CALL GROUP1(K,ICONT(2),JCONT,2,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(2)
              JCOUPC=NCOUPC(2)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 3031
            END IF
            JCOUPL=JCOUPS
            JCOUPC=JCOUCS
            ICOUPL=ICOUPS
            ICOUPC=ICOUCS
3031        CONTINUE
          END IF

          IF(ICOUPL.LT.1)GO TO 3301
          IF(LDUMP.GT.0)THEN
C**WRITE DIPOLE GRIDS TO DISC
            CALL DUMDP1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     1      W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,
     2      W(LV1),W(LV1),LDUMP)
          ELSE
C**WRITE POTENTIAL GRIDS TO DISC
C**DUMPT1 CALLED FOR LINEAR TRIATOMICS (NSMODE=NAMODE=3), 
C**NON-LINEAR WITHOUT RPH (JREACT=0), 
C**LINEAR STRETCH-ONLY (NONLIN) MODES
            IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
              IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1        CALL DUMPT1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LRR),
     2        W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),
     3        W(LCPOT),K,W(LV1),W(LV1),1,W(LMODNT),W(LXTANH),W(LNP1),
     4        W(LCP1),W(LMP1),W(LVP1),W(LDP1),IENTMX(1),W(LXMODQ))
            ELSE
C**DUMVT1 CALLED FOR NON-LINEAR RPH (NSMODE=NAMODE+1, JREACT>0),
C**LINEAR BEND-ONLY MODES (K>NONLIN) FOR TETRAATOMICS (NSMODE>NAMODE)
C**INTRINSIC
              KVP=MBFMAX
              CALL MEMO(2,LVP0,IK2TAU,LVP1,KVP*IK2TAU,0,0,0,0,0,0)
              CALL DUMVT1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LRR),
     1        W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),
     2        W(LCPOT),K,W(LV1),W(LV1),W(LXQ+K3TAU),W(LVP0),W(LVP1),
     3        KVP,1,W(LMODNT),W(LXMODQ))
              CALL MEMO(-2,LVP0,IK2TAU,LVP1,KVP*IK2TAU,0,0,0,0,0,0)
            END IF
            IF(ICOUPC.LT.1)GO TO 3301
C**WRITE CORIOLIS GRIDS TO DISC
            IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
              CALL DUMCR1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LXZ),
     1        W(LAB),W(LB),W(LAA),W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,
     2        W(LC1),W(LC1),W(LVM1),W(LVM1),0,1,W(LMODNT))
            ELSE
              IF(MOLPRO.EQ.0.AND.MOLINC.GT.0)
     1        CALL MEMO(2,LVP0,3*IK2TAU,LVP1,6*KVP*IK2TAU,0,0,0,0,0,0)
              CALL DUMVR1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LXZ),
     1        W(LAB),W(LB),W(LB+NMODE*NMODE),W(LAA),W(LBB),
     2        W(LBB+NMODE),W(LXX),W(LXX+NATOM*3*722),W(LX0),W(LXL),
     3        W(LXM),K,W(LC1),W(LC1),W(LVM1),W(LVM1),0,W(LXQ+K3TAU),
     4        W(LVP0),W(LVP1),KVP,1,W(LMODNT))
              IF(MOLPRO.EQ.0.AND.MOLINC.GT.0)
     1        CALL MEMO(-2,LVP0,3*IK2TAU,LVP1,6*KVP*IK2TAU,0,0,0,0,0,0)
            END IF
          END IF
          IF(ICOUPC.LT.1)GO TO 3301
C         IF(JMAX.GT.0.AND.ICOUPL.EQ.1)THEN
          IF(JMAX.GT.0)THEN
            IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
              CALL GETMI1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LXZ),
     1        W(LAB),W(LB),W(LAA),W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,
     2        W(LEJK1),W(LEJK1),W(LSS),J21,W(LX21),W(LW21),W(LE21),
     3        W(LMODNT))
            ELSE
              CALL GETMV1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LXZ),
     1        W(LAB),W(LB),W(LB+NMODE*NMODE),W(LAA),W(LBB),
     2        W(LBB+NMODE),W(LXX),W(LXX+NATOM*3*722),W(LX0),W(LXL),
     3        W(LXM),K,W(LEJK1),W(LEJK1),W(LSS),J21,W(LX21),W(LW21),
     4        W(LE21),W(LMODNT))
            END IF
          END IF
3301  CONTINUE
          K3=K3+IK2
        END DO
      END IF
C******************************
C**INTRINSIC
C     IF(ICOUPL.EQ.1)GO TO 9999
      IF(ICOUPX.EQ.1)GO TO 9999
      MBFMX2=0
      DO KK=1,NAMODE
        IF(LCONT.NE.0)THEN
          K=JCONT(LCONT,KK)
        ELSE
          K=KK
        END IF
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        DO LL=1,KK-1
          IF(LCONT.NE.0)THEN
            L=JCONT(LCONT,LL)
          ELSE
            L=LL
          END IF
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
          IF(IK2*IL2.GT.MBFMX2)MBFMX2=IK2*IL2
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        KLV2=MBFMX2
      ELSE
        KLV2=(1+MBFMX2)/2
      END IF
C**INTRINSIC
      IF(MOLINC.GT.0.AND.ICOUCX.GE.3)KLV2=6*KLV2
      KLVR2=KLV2
      IF(JREACT.NE.0.AND.MOLINC.GT.0.AND.ICOUPX.GE.3)KLVR2=KLV2*IK2TAU
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.2)THEN
          KLC2=MBFMX2*6
          IF(IREACT.NE.0)KLC2=MBFMX2*10
        END IF
      ELSE
        IF(ICOUCX.GE.2)THEN
          KLC2=(1+MBFMX2*6)/2
          IF(IREACT.NE.0)KLC2=(1+MBFMX2*10)/2
        END IF
      END IF
C**INTRINSIC
      KLCR2=KLC2
      IF(JREACT.NE.0.AND.MOLINC.GT.0.AND.ICOUCX.GE.3)KLCR2=KLC2*IK2TAU
      IF(LCONT.EQ.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1  CALL MEMO(1,LV2,KLVR2,0,0,0,0,0,0,0,0)
        IF(ICOUCX.GE.2)CALL MEMO(1,LC2,KLCR2,0,0,0,0,0,0,0,0)
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.2)KEJK2=MBFMX2*J21
      ELSE
        IF(ICOUCX.GE.2)KEJK2=(1+MBFMX2*J21)/2
      END IF
      IF(LCONT.EQ.0)THEN
      IF(ICOUCX.GE.2)CALL MEMO(1,LEJK2,KEJK2,0,0,0,0,0,0,0,0)
      END IF
      IF(IDISC.EQ.0)THEN
        K=0
        K3=0
        DO KK=1,NAMODE
          IF(LCONT.NE.0)THEN
            KNEXT=JCONT(LCONT,KK)
          ELSE
            KNEXT=KK
          END IF
          IF(KNEXT-K.GT.1)THEN
C**UPDATE K3 FOR MISSING MODES
            DO KADJ=K+1,KNEXT-1
              CALL INTARR(W(LNBF),W(LMBF),KADJ,IK1,IK2,IK3)
              K3=K3+IK2
            END DO
          END IF
C**NEXT K
          K=KNEXT
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L=0
          L3=0
          DO LL=1,KK-1
            IF(LCONT.NE.0)THEN
              LNEXT=JCONT(LCONT,LL)
            ELSE
              LNEXT=LL
            END IF
            IF(LNEXT-L.GT.1)THEN
C**UPDATE L3 FOR MISSING MODES
              DO LADJ=L+1,LNEXT-1
                CALL INTARR(W(LNBF),W(LMBF),LADJ,IL1,IL2,IL3)
                L3=L3+IL2
              END DO
            END IF
C**NEXT L
            L=LNEXT

            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**INTRINSIC
            IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
              CALL GROUP2(K,L,ICONT(1),JCONT,1,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(1)
                JCOUPC=NCOUPC(1)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 3032
              END IF
              CALL GROUP2(K,L,ICONT(2),JCONT,2,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(2)
                JCOUPC=NCOUPC(2)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 3032
              END IF
              JCOUPL=JCOUPS
              JCOUPC=JCOUCS
              ICOUPL=ICOUPS
              ICOUPC=ICOUCS
3032          CONTINUE
            END IF

            IF(ICOUPL.LT.2)GO TO 3302
            IF(LDUMP.GT.0)THEN
C**WRITE DIPOLE GRIDS TO DISC
              CALL DUMDP2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1        W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),
     2        W(LJPOT),W(LCPOT),K,L,W(LV2),W(LV2),LDUMP)
            ELSE
C**WRITE POTENTIAL GRIDS TO DISC
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1        (NSMODE.EQ.NAMODE))THEN
C**INTRINSIC
                KVP=MBFMAX
                CALL MEMO(1,LVP1,KVP*2,0,0,0,0,0,0,0,0)
                IENTMX(1)=KVP
                IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1          CALL DUMPT2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     2          W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),NPOT,
     3          W(LIPOT),W(LJPOT),W(LCPOT),K,L,W(LV2),W(LV2),1,
     4          W(LMODNT),W(LXTANH),W(LNP1),W(LCP1),W(LVP1),W(LDP1),
     5          NTOT1,IENTMX(1),W(LNP2),W(LCP2),W(LVP2),W(LDP2A),
     6          W(LDP2B),NTOT2,IENTMX(2),W(LINDK),W(LNBF),W(LMBF),
     7          W(LV1),W(LV1),W(LXMODQ))
                CALL MEMO(-1,LVP1,KVP*2,0,0,0,0,0,0,0,0)
              ELSE
C**INTRINSIC
                KVP=MBFMAX
                CALL MEMO(3,LVP0,IK2TAU,LVP1,2*KVP*IK2TAU,
     1          LVP2,KVP*KVP*IK2TAU,0,0,0,0)
                CALL DUMVT2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1          W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),NPOT,
     2          W(LIPOT),W(LJPOT),W(LCPOT),K,L,W(LV2),W(LV2),
     3          W(LXQ+K3TAU),W(LVP0),W(LVP1),W(LVP2),KVP,1,
     4          W(LMODNT),W(LNBF),W(LMBF),W(LV1),W(LV1),IK2TAU,
     5          W(LXMODQ))
                CALL MEMO(-3,LVP0,IK2TAU,LVP1,2*KVP*IK2TAU,
     1          LVP2,KVP*KVP*IK2TAU,0,0,0,0)
              END IF
              IF(ICOUPC.LT.2)GO TO 3302
C**WRITE CORIOLIS GRIDS TO DISC
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1        (NSMODE.EQ.NAMODE))THEN
C**INTRINSIC
                KVP=MBFMAX
                CALL MEMO(1,LVP1,KVP*2,0,0,0,0,0,0,0,0)
                CALL DUMCR2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1          W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LXX),W(LX0),
     2          W(LXL),W(LXM),K,L,W(LC2),W(LC2),W(LVM2),W(LVM2),0,1,
C    3          W(LMODNT),W(LNBF),W(LMBF),W(LVP1),W(LV1),W(LV1),KVP)
     3          W(LMODNT),W(LNBF),W(LMBF),W(LVP1),W(LC1),W(LC1),KVP)
                CALL MEMO(-1,LVP1,KVP*2,0,0,0,0,0,0,0,0)
              ELSE
                IF(MOLPRO.EQ.0.AND.MOLINC.GT.0)
     1          CALL MEMO(3,LVP0,3*IK2TAU,LVP1,6*KVP*IK2TAU*2,
     2          LVP2,10*KVP*KVP*IK2TAU,0,0,0,0)
                CALL DUMVR2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1          W(LQQ),W(LXZ),W(LAB),W(LB),W(LB+NMODE*NMODE),W(LAA),
     2          W(LBB),W(LBB+NMODE),W(LXX),W(LXX+NATOM*3*722),W(LX0),
     3          W(LXL),W(LXM),K,L,W(LC2),W(LC2),W(LVM2),W(LVM2),0,
     4          W(LXQ+K3TAU),W(LVP0),W(LVP1),W(LVP2),KVP,1,W(LMODNT),
     5          W(LNBF),W(LMBF),W(LC1),W(LC1),IK2TAU)
                IF(MOLPRO.EQ.0.AND.MOLINC.GT.0)
     1          CALL MEMO(-3,LVP0,3*IK2TAU,LVP1,6*KVP*IK2TAU*2,
     2          LVP2,10*KVP*KVP*IK2TAU,0,0,0,0)
              END IF
            END IF
            IF(ICOUPC.LT.2)GO TO 3302
C           IF(JMAX.GT.0.AND.ICOUPL.EQ.2)THEN
            IF(JMAX.GT.0)THEN
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1        (NSMODE.EQ.NAMODE))THEN
                CALL GETMI2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1          W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LXX),W(LX0),
     2          W(LXL),W(LXM),K,L,W(LEJK2),W(LEJK2),W(LSS),J21,W(LX21),
     3          W(LW21),W(LE21),W(LMODNT))
              ELSE
                CALL GETMV2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1          W(LQQ),W(LXZ),W(LAB),W(LB),W(LB+NMODE*NMODE),W(LAA),
     2          W(LBB),W(LBB+NMODE),W(LXX),W(LXX+NATOM*3*722),W(LX0),
     3          W(LXL),W(LXM),K,L,W(LEJK2),W(LEJK2),W(LSS),J21,W(LX21),
     4          W(LW21),W(LE21),W(LMODNT))
              END IF
            END IF
3302  CONTINUE
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C******************************
C**INTRINSIC
C     IF(ICOUPL.EQ.2)GO TO 9999
      IF(ICOUPX.EQ.2)GO TO 9999
      MBFMX3=0
      DO KK=1,NAMODE
        IF(LCONT.NE.0)THEN
          K=JCONT(LCONT,KK)
        ELSE
          K=KK
        END IF
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        DO LL=1,KK-1
          IF(LCONT.NE.0)THEN
            L=JCONT(LCONT,LL)
          ELSE
            L=LL
          END IF
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
          DO NN=1,LL-1
            IF(LCONT.NE.0)THEN
              N=JCONT(LCONT,NN)
            ELSE
              N=NN
            END IF
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
            IF(IK2*IL2*IN2.GT.MBFMX3)MBFMX3=IK2*IL2*IN2
          END DO
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        KLV3=MBFMX3
      ELSE
        KLV3=(1+MBFMX3)/2
      END IF
C**INTRINSIC
      IF(MOLINC.GT.0.AND.ICOUCX.GE.4)KLV3=10*KLV3
      KLVR3=KLV3
      IF(JREACT.NE.0.AND.MOLINC.GT.0.AND.ICOUPX.GE.4)KLVR3=KLV3*IK2TAU
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.3)THEN
          KLC3=MBFMX3*10
          IF(IREACT.NE.0)KLC3=MBFMX3*15
        END IF
      ELSE
        IF(ICOUCX.GE.3)THEN
          KLC3=(1+MBFMX3*10)/2
          IF(IREACT.NE.0)KLC3=(1+MBFMX3*15)/2
        END IF
      END IF
C**INTRINSIC
      KLCR3=KLC3
      IF(JREACT.NE.0.AND.MOLINC.GT.0.AND.ICOUCX.GE.4)KLCR3=KLC2*IK2TAU
      IF(LCONT.EQ.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1  CALL MEMO(1,LV3,KLVR3,0,0,0,0,0,0,0,0)
        IF(ICOUCX.GE.3)CALL MEMO(1,LC3,KLCR3,0,0,0,0,0,0,0,0)
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.3)KEJK3=MBFMX3*J21
      ELSE
        IF(ICOUCX.GE.3)KEJK3=(1+MBFMX3*J21)/2
      END IF
      IF(LCONT.EQ.0)THEN
        IF(ICOUCX.GE.3)CALL MEMO(1,LEJK3,KEJK3,0,0,0,0,0,0,0,0)
      END IF
      IF(IDISC.EQ.0)THEN
        K=0
        K3=0
        DO KK=1,NAMODE
          IF(LCONT.NE.0)THEN
            KNEXT=JCONT(LCONT,KK)
          ELSE
            KNEXT=KK
          END IF
          IF(KNEXT-K.GT.1)THEN
C**UPDATE K3 FOR MISSING MODES
            DO KADJ=K+1,KNEXT-1
              CALL INTARR(W(LNBF),W(LMBF),KADJ,IK1,IK2,IK3)
              K3=K3+IK2
            END DO
          END IF
C**NEXT K
          K=KNEXT
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L=0
          L3=0
          DO LL=1,KK-1
            IF(LCONT.NE.0)THEN
              LNEXT=JCONT(LCONT,LL)
            ELSE
              LNEXT=LL
            END IF
            IF(LNEXT-L.GT.1)THEN
C**UPDATE L3 FOR MISSING MODES
              DO LADJ=L+1,LNEXT-1
                CALL INTARR(W(LNBF),W(LMBF),LADJ,IL1,IL2,IL3)
                L3=L3+IL2
              END DO
            END IF
C**NEXT L
            L=LNEXT
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N=0
            N3=0
            DO NN=1,LL-1
              IF(LCONT.NE.0)THEN
                NNEXT=JCONT(LCONT,NN)
              ELSE
                NNEXT=NN
              END IF
              IF(NNEXT-N.GT.1)THEN
C**UPDATE N3 FOR MISSING MODES
                DO NADJ=N+1,NNEXT-1
                  CALL INTARR(W(LNBF),W(LMBF),NADJ,IN1,IN2,IN3)
                  N3=N3+IN2
                END DO
              END IF
C**NEXT N
              N=NNEXT
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**INTRINSIC
              IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
                CALL GROUP3(K,L,N,ICONT(1),JCONT,1,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(1)
                  JCOUPC=NCOUPC(1)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 3033
                END IF
                CALL GROUP3(K,L,N,ICONT(2),JCONT,2,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(2)
                  JCOUPC=NCOUPC(2)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 3033
                END IF
                JCOUPL=JCOUPS
                JCOUPC=JCOUCS
                ICOUPL=ICOUPS
                ICOUPC=ICOUCS
3033            CONTINUE
              END IF
              IF(ICOUPL.LT.3)GO TO 3303
              IF(LDUMP.GT.0)THEN
C**WRITE DIPOLE GRIDS TO DISC
                CALL DUMDP3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,IN2,
     1          NMODE,NATOM,W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     2          NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,L,N,W(LV3),W(LV3),
     3          LDUMP)
              ELSE
C**WRITE POTENTIAL GRIDS TO DISC
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**INTRINSIC
                  KVP=MBFMAX
                  CALL MEMO(2,LVP1,KVP*3,LVP2,KVP*KVP*3,
     1            0,0,0,0,0,0)
                  IENTMX(1)=KVP
                  IENTMX(2)=KVP
                  IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1            CALL DUMPT3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,
     2            IN2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),
     3            W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,L,N,W(LV3),
     4            W(LV3),1,W(LMODNT),W(LXTANH),W(LNP1),W(LCP1),W(LVP1),
     5            W(LDP1),NTOT1,IENTMX(1),W(LNP2),W(LCP2),W(LVP2),
     6            W(LDP2A),W(LDP2B),NTOT2,IENTMX(2),W(LNP3),W(LCP3),
     7            W(LVP3),W(LDP3A),W(LDP3B),W(LDP3C),NTOT3,IENTMX(3),
     8            W(LINDK),W(LINDL),W(LNBF),W(LMBF),W(LV1),W(LV1),
     9            W(LV2),W(LV2),W(LXMODQ))
                  CALL MEMO(-2,LVP1,KVP*3,LVP2,KVP*KVP*3,
     1            0,0,0,0,0,0)
                ELSE
C**INTRINSIC
                  KVP=MBFMAX
                  CALL MEMO(4,LVP0,IK2TAU,LVP1,3*KVP*IK2TAU,
     1            LVP2,3*KVP*KVP*IK2TAU,LVP3,KVP*KVP*KVP*IK2TAU,0,0)
                  CALL DUMVT3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,
     1            IN2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),
     2            W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,L,N,W(LV3),
     3            W(LV3),W(LXQ+K3TAU),W(LVP0),W(LVP1),W(LVP2),W(LVP3),
     4            KVP,1,W(LMODNT),W(LNBF),W(LMBF),W(LV1),W(LV1),W(LV2),
     5            W(LV2),IK2TAU,W(LXMODQ))
                  CALL MEMO(-4,LVP0,IK2TAU,LVP1,3*KVP*IK2TAU,
     1            LVP2,3*KVP*KVP*IK2TAU,LVP3,KVP*KVP*KVP*IK2TAU,0,0)
                END IF
                IF(ICOUPC.LT.3)GO TO 3303
C**WRITE CORIOLIS GRIDS TO DISC
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**INTRINSIC
                  CALL MEMO(2,LVP1,KVP*3,LVP2,KVP*KVP*6*3,
     1            0,0,0,0,0,0)
                  CALL DUMCR3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,
     1            IN2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),
     2            W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,L,N,W(LC3),
     3            W(LC3),W(LVM3),W(LVM3),0,1,W(LMODNT),W(LNBF),W(LMBF),
C    4            W(LVP1),W(LVP2),W(LV1),W(LV1),W(LV2),W(LV2),KVP)
     4            W(LVP1),W(LVP2),W(LC1),W(LC1),W(LC2),W(LC2),KVP)
                  CALL MEMO(-2,LVP1,KVP*3,LVP2,KVP*KVP*6*3,
     1            0,0,0,0,0,0)
                ELSE
                  IF(MOLPRO.EQ.0.AND.MOLINC.GT.0)
     1            CALL MEMO(4,LVP0,3*IK2TAU,LVP1,6*KVP*IK2TAU*3,
     2            LVP2,10*KVP*KVP*IK2TAU*3,LVP3,15*KVP*KVP*KVP*IK2TAU,
     3            0,0)
                  CALL DUMVR3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,
     1            IN2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     2            W(LB+NMODE*NMODE),W(LAA),W(LBB),W(LBB+NMODE),W(LXX),
     3            W(LXX+NATOM*3*722),W(LX0),W(LXL),W(LXM),K,L,N,W(LC3),
     4            W(LC3),W(LVM3),W(LVM3),0,W(LXQ+K3TAU),W(LVP0),
     5            W(LVP1),W(LVP2),W(LVP3),KVP,1,W(LMODNT),W(LNBF),
     6            W(LMBF),W(LC1),W(LC1),W(LC2),W(LC2),IK2TAU)
                  IF(MOLPRO.EQ.0.AND.MOLINC.GT.0)
     1            CALL MEMO(-4,LVP0,3*IK2TAU,LVP1,6*KVP*IK2TAU*3,
     2            LVP2,10*KVP*KVP*IK2TAU*3,LVP3,15*KVP*KVP*KVP*IK2TAU,
     2            0,0)
                END IF
              END IF
              IF(ICOUPC.LT.3)GO TO 3303
C             IF(JMAX.GT.0.AND.ICOUPL.EQ.3)THEN
              IF(JMAX.GT.0)THEN
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
                  CALL GETMI3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,
     1            IN2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),
     2            W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,L,N,W(LEJK3),
     3            W(LEJK3),W(LSS),J21,W(LX21),W(LW21),W(LE21),
     4            W(LMODNT))
                ELSE
                  CALL GETMV3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,
     1            IN2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     2            W(LB+NMODE*NMODE),W(LAA),W(LBB),W(LBB+NMODE),W(LXX),
     3            W(LXX+NATOM*3*722),W(LX0),W(LXL),W(LXM),K,L,N,
     4            W(LEJK3),W(LEJK3),W(LSS),J21,W(LX21),W(LW21),
     5            W(LE21),W(LMODNT))
                END IF
              END IF
3303  CONTINUE
              N3=N3+IN2
            END DO
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C******************************
C**INTRINSIC
C     IF(ICOUPL.EQ.3)GO TO 9999
      IF(ICOUPX.EQ.3)GO TO 9999
      IF(ICOUCX.LT.4.AND.IFITRP.GT.3)GO TO 9999
      MBFMX4=0
      DO KK=1,NAMODE
        IF(LCONT.NE.0)THEN
          K=JCONT(LCONT,KK)
        ELSE
          K=KK
        END IF
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        DO LL=1,KK-1
          IF(LCONT.NE.0)THEN
            L=JCONT(LCONT,LL)
          ELSE
            L=LL
          END IF
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
          DO NN=1,LL-1
            IF(LCONT.NE.0)THEN
              N=JCONT(LCONT,NN)
            ELSE
              N=NN
            END IF
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
            DO MM=1,NN-1
              IF(LCONT.NE.0)THEN
                M=JCONT(LCONT,MM)
              ELSE
                M=MM
              END IF
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
              IF(IK2*IL2*IN2*IM2.GT.MBFMX4)MBFMX4=IK2*IL2*IN2*IM2
            END DO
          END DO
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        KLV4=MBFMX4
      ELSE
        KLV4=(1+MBFMX4)/2
      END IF
C**INTRINSIC
      IF(MOLINC.GT.0.AND.ICOUCX.GE.5)KLV4=15*KLV4
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.4)THEN
          KLC4=MBFMX4*15
          IF(IREACT.NE.0)KLC4=MBFMX4*21
        END IF
      ELSE
        IF(ICOUCX.GE.4)THEN
          KLC4=(1+MBFMX4*15)/2
          IF(IREACT.NE.0)KLC4=(1+MBFMX4*21)/2
        END IF
      END IF
      IF(LCONT.EQ.0)THEN
        IF((IWHICH.GE.0.OR.MOLINC.LE.0).AND.IABS(IFITRP).LT.4)
     1  CALL MEMO(1,LV4,KLV4,0,0,0,0,0,0,0,0)
        IF(ICOUCX.GE.4)CALL MEMO(1,LC4,KLC4,0,0,0,0,0,0,0,0)
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.4)KEJK4=MBFMX4*J21
      ELSE
        IF(ICOUCX.GE.4)KEJK4=(1+MBFMX4*J21)/2
      END IF
      IF(LCONT.EQ.0)THEN
        IF(ICOUCX.GE.4)CALL MEMO(1,LEJK4,KEJK4,0,0,0,0,0,0,0,0)
      END IF
      IF(IDISC.EQ.0)THEN
        K=0
        K3=0
        DO KK=1,NAMODE
          IF(LCONT.NE.0)THEN
            KNEXT=JCONT(LCONT,KK)
          ELSE
            KNEXT=KK
          END IF
          IF(KNEXT-K.GT.1)THEN
C**UPDATE K3 FOR MISSING MODES
            DO KADJ=K+1,KNEXT-1
              CALL INTARR(W(LNBF),W(LMBF),KADJ,IK1,IK2,IK3)
              K3=K3+IK2
            END DO
          END IF
C**NEXT K
          K=KNEXT
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L=0
          L3=0
          DO LL=1,KK-1
            IF(LCONT.NE.0)THEN
              LNEXT=JCONT(LCONT,LL)
            ELSE
              LNEXT=LL
            END IF
            IF(LNEXT-L.GT.1)THEN
C**UPDATE L3 FOR MISSING MODES
              DO LADJ=L+1,LNEXT-1
                CALL INTARR(W(LNBF),W(LMBF),LADJ,IL1,IL2,IL3)
                L3=L3+IL2
              END DO
            END IF
C**NEXT L
            L=LNEXT
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N=0
            N3=0
            DO NN=1,LL-1
              IF(LCONT.NE.0)THEN
                NNEXT=JCONT(LCONT,NN)
              ELSE
                NNEXT=NN
              END IF
              IF(NNEXT-N.GT.1)THEN
C**UPDATE N3 FOR MISSING MODES
                DO NADJ=N+1,NNEXT-1
                  CALL INTARR(W(LNBF),W(LMBF),NADJ,IN1,IN2,IN3)
                  N3=N3+IN2
                END DO
              END IF
C**NEXT N
              N=NNEXT
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
              M=0
              M3=0
              DO MM=1,NN-1
                IF(LCONT.NE.0)THEN
                  MNEXT=JCONT(LCONT,MM)
                ELSE
                  MNEXT=MM
                END IF
                IF(MNEXT-M.GT.1)THEN
C**UPDATE M3 FOR MISSING MODES
                  DO MADJ=M+1,MNEXT-1
                    CALL INTARR(W(LNBF),W(LMBF),MADJ,IM1,IM2,IM3)
                    M3=M3+IM2
                  END DO
                END IF
C**NEXT M
                M=MNEXT
                CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**INTRINSIC
                IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
                  CALL GROUP4(K,L,N,M,ICONT(1),JCONT,1,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(1)
                    JCOUPC=NCOUPC(1)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 3034
                  END IF
                  CALL GROUP4(K,L,N,M,ICONT(2),JCONT,2,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(2)
                    JCOUPC=NCOUPC(2)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 3034
                  END IF
                  JCOUPL=JCOUPS
                  JCOUPC=JCOUCS
                  ICOUPL=ICOUPS
                  ICOUPC=ICOUCS
3034              CONTINUE
                END IF
                IF(ICOUPL.LT.4)GO TO 3304
                IF(LDUMP.GT.0)THEN
C**WRITE DIPOLE GRIDS TO DISC
                  CALL DUMDP4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),W(LXQ+M3),
     1            IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     2            W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     3            K,L,N,M,W(LV4),W(LV4),LDUMP)
                ELSE
C**WRITE POTENTIAL GRIDS TO DISC
                  IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1            N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2            THEN
C**INTRINSIC
                    KVP=MBFMAX
                    CALL MEMO(3,LVP1,KVP*4,LVP2,KVP*KVP*6,
     1              LVP3,KVP*KVP*KVP*4,0,0,0,0)
                    IENTMX(1)=KVP
                    IENTMX(2)=KVP
                    IENTMX(3)=KVP
                    IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1              CALL DUMPT4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     2              W(LXQ+M3),IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),
     3              W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),
     4              W(LJPOT),W(LCPOT),K,L,N,M,W(LV4),W(LV4),1,
     5              W(LMODNT),W(LXTANH),W(LNP1),W(LCP1),W(LVP1),
     5              W(LDP1),NTOT1,IENTMX(1),W(LNP2),W(LCP2),W(LVP2),
     6              W(LDP2A),W(LDP2B),NTOT2,IENTMX(2),W(LNP3),W(LCP3),
     7              W(LVP3),W(LDP3A),W(LDP3B),W(LDP3C),NTOT3,IENTMX(3),
     8              W(LNP4),W(LCP4),W(LVP4),W(LDP4A),W(LDP4B),W(LDP4C),
     9              W(LDP4D),NTOT4,IENTMX(4),W(LINDK),W(LINDL),
     1              W(LINDN),W(LNBF),W(LMBF),W(LV1),W(LV1),W(LV2),
     2              W(LV2),W(LV3),W(LV3),W(LXMODQ))
                    CALL MEMO(-3,LVP1,KVP*4,LVP2,KVP*KVP*6,
     1              LVP3,KVP*KVP*KVP*4,0,0,0,0)
                  ELSE
                    IF(IABS(IFITRP).LT.4)THEN
C**INTRINSIC
                      KVP=MBFMAX
                      CALL MEMO(5,LVP0,IK2TAU,LVP1,4*KVP*IK2TAU,
     1                LVP2,6*KVP*KVP*IK2TAU,LVP3,4*KVP*KVP*KVP*IK2TAU,
     2                LVP4,KVP*KVP*KVP*KVP*IK2TAU)
                      CALL DUMVT4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1                W(LXQ+M3),IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),
     2                W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),
     3                W(LJPOT),W(LCPOT),K,L,N,M,W(LV4),W(LV4),
     4                W(LXQ+K3TAU),W(LVP0),W(LVP1),W(LVP2),W(LVP3),
     5                W(LVP4),KVP,1,W(LMODNT),W(LNBF),W(LMBF),W(LV1),
     6                W(LV1),W(LV2),W(LV2),W(LV3),W(LV3),IK2TAU,
     7                W(LXMODQ))
                      CALL MEMO(-5,LVP0,IK2TAU,LVP1,4*KVP*IK2TAU,
     1                LVP2,6*KVP*KVP*IK2TAU,LVP3,4*KVP*KVP*KVP*IK2TAU,
     2                LVP4,KVP*KVP*KVP*KVP*IK2TAU)
                    END IF
                  END IF
                  IF(ICOUPC.LT.4)GO TO 3304
C**WRITE CORIOLIS GRIDS TO DISC
                  IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1            N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2            THEN
C**INTRINSIC
                    CALL MEMO(3,LVP1,KVP*4,LVP2,KVP*KVP*6*6,
     1              LVP3,KVP*KVP*KVP*10*4,0,0,0,0)
                    CALL DUMCR4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1              W(LXQ+M3),IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),
     2              W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LXX),W(LX0),
     3              W(LXL),W(LXM),K,L,N,M,W(LC4),W(LC4),W(LVM4),
     4              W(LVM4),0,1,W(LMODNT),W(LNBF),W(LMBF),W(LVP1),
C    5              W(LVP2),W(LVP3),W(LV1),W(LV1),W(LV2),W(LV2),
C    6              W(LV3),W(LV3),KVP)
     5              W(LVP2),W(LVP3),W(LC1),W(LC1),W(LC2),W(LC2),
     6              W(LC3),W(LC3),KVP)
                    CALL MEMO(-3,LVP1,KVP*4,LVP2,KVP*KVP*6*6,
     1              LVP3,KVP*KVP*KVP*10*4,0,0,0,0)
                  ELSE
                    IF(MOLPRO.EQ.0.AND.MOLINC.GT.0)
     1              CALL MEMO(5,LVP0,3*IK2TAU,LVP1,6*KVP*IK2TAU*4,
     2              LVP2,10*KVP*KVP*IK2TAU*6,LVP3,
     3              15*KVP*KVP*KVP*IK2TAU*4,LVP4,
     4              21*KVP*KVP*KVP*KVP*IK2TAU)
                    CALL DUMVR4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1              W(LXQ+M3),IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),
     2              W(LXZ),W(LAB),W(LB),W(LB+NMODE*NMODE),W(LAA),
     3              W(LBB),W(LBB+NMODE),W(LXX),W(LXX+NATOM*3*722),
     4              W(LX0),W(LXL),W(LXM),K,L,N,M,W(LC4),W(LC4),W(LVM4),
     5              W(LVM4),0,W(LXQ+K3TAU),W(LVP0),W(LVP1),W(LVP2),
     6              W(LVP3),W(LVP4),KVP,1,W(LMODNT),W(LNBF),W(LMBF),
     7              W(LC1),W(LC1),W(LC2),W(LC2),W(LC3),W(LC3),IK2TAU)
                    IF(MOLPRO.EQ.0.AND.MOLINC.GT.0)
     1              CALL MEMO(-5,LVP0,3*IK2TAU,LVP1,6*KVP*IK2TAU*4,
     2              LVP2,10*KVP*KVP*IK2TAU*6,LVP3,
     3              15*KVP*KVP*KVP*IK2TAU*4,LVP4,
     4              21*KVP*KVP*KVP*KVP*IK2TAU)
                  END IF
                END IF
                IF(ICOUPC.LT.4)GO TO 3304
C               IF(JMAX.GT.0.AND.ICOUPL.EQ.4)THEN
                IF(JMAX.GT.0)THEN
                  IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1            N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2            THEN
                    CALL GETMI4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1              W(LXQ+M3),IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),
     2              W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LXX),W(LX0),
     3              W(LXL),W(LXM),K,L,N,M,W(LEJK4),W(LEJK4),W(LSS),J21,
     4              W(LX21),W(LW21),W(LE21),W(LMODNT))
                  ELSE
                    CALL GETMV4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1              W(LXQ+M3),IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),
     2              W(LXZ),W(LAB),W(LB),W(LB+NMODE*NMODE),W(LAA),
     3              W(LBB),W(LBB+NMODE),W(LXX),W(LXX+NATOM*3*722),
     4              W(LX0),W(LXL),W(LXM),K,L,N,M,W(LEJK4),W(LEJK4),
     5              W(LSS),J21,W(LX21),W(LW21),W(LE21),W(LMODNT))
                  END IF
                END IF
3304  CONTINUE
                M3=M3+IM2
              END DO
              N3=N3+IN2
            END DO
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C******************************
C**INTRINSIC
C     IF(ICOUPL.EQ.4)GO TO 9999
      IF(ICOUPX.EQ.4)GO TO 9999
      MBFMX5=0
      DO KK=1,NAMODE
        IF(LCONT.NE.0)THEN
          K=JCONT(LCONT,KK)
        ELSE
          K=KK
        END IF
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        DO LL=1,KK-1
          IF(LCONT.NE.0)THEN
            L=JCONT(LCONT,LL)
          ELSE
            L=LL
          END IF
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
          DO NN=1,LL-1
            IF(LCONT.NE.0)THEN
              N=JCONT(LCONT,NN)
            ELSE
              N=NN
            END IF
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
            DO MM=1,NN-1
              IF(LCONT.NE.0)THEN
                M=JCONT(LCONT,MM)
              ELSE
                M=MM
              END IF
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
              DO JJ=1,MM-1
                IF(LCONT.NE.0)THEN
                  J=JCONT(LCONT,JJ)
                ELSE
                  J=JJ
                END IF
                CALL INTARR(W(LNBF),W(LMBF),J,IJ1,IJ2,IJ3)
                IF(IK2*IL2*IN2*IM2*IJ2.GT.MBFMX5)
     1          MBFMX5=IK2*IL2*IN2*IM2*IJ2
              END DO
            END DO
          END DO
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        KLV5=MBFMX5
      ELSE
        KLV5=(1+MBFMX5)/2
      END IF
C**INTRINSIC
      IF(MOLINC.GT.0.AND.ICOUCX.GE.6)KLV5=16*KLV5
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.5)THEN
          KLC5=MBFMX4*16
CC        IF(IREACT.NE.0)KLC5=MBFMX4*21
        END IF
      ELSE
        IF(ICOUCX.GE.5)THEN
          KLC5=(1+MBFMX4*16)/2
CC        IF(IREACT.NE.0)KLC5=(1+MBFMX4*21)/2
        END IF
      END IF
      IF(LCONT.EQ.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1  CALL MEMO(1,LV5,KLV5,0,0,0,0,0,0,0,0)
        IF(ICOUCX.GE.5)CALL MEMO(1,LC5,KLC5,0,0,0,0,0,0,0,0)
      END IF
CC    IF(JCOUPC.GE.0)THEN
CC      IF(ICOUCX.GE.5)KEJK5=MBFMX5*J21
CC    ELSE
CC      IF(ICOUCX.GE.5)KEJK5=(1+MBFMX5*J21)/2
CC    END IF
CC    IF(LCONT.EQ.0)THEN
CC      IF(ICOUCX.GE.5)CALL MEMO(1,LEJK5,KEJK5,0,0,0,0,0,0,0,0)
CC    END IF
      IF(IDISC.EQ.0)THEN
        K=0
        K3=0
        DO KK=1,NAMODE
          IF(LCONT.NE.0)THEN
            KNEXT=JCONT(LCONT,KK)
          ELSE
            KNEXT=KK
          END IF
          IF(KNEXT-K.GT.1)THEN
C**UPDATE K3 FOR MISSING MODES
            DO KADJ=K+1,KNEXT-1
              CALL INTARR(W(LNBF),W(LMBF),KADJ,IK1,IK2,IK3)
              K3=K3+IK2
            END DO
          END IF
C**NEXT K
          K=KNEXT
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L=0
          L3=0
          DO LL=1,KK-1
            IF(LCONT.NE.0)THEN
              LNEXT=JCONT(LCONT,LL)
            ELSE
              LNEXT=LL
            END IF
            IF(LNEXT-L.GT.1)THEN
C**UPDATE L3 FOR MISSING MODES
              DO LADJ=L+1,LNEXT-1
                CALL INTARR(W(LNBF),W(LMBF),LADJ,IL1,IL2,IL3)
                L3=L3+IL2
              END DO
            END IF
C**NEXT L
            L=LNEXT
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N=0
            N3=0
            DO NN=1,LL-1
              IF(LCONT.NE.0)THEN
                NNEXT=JCONT(LCONT,NN)
              ELSE
                NNEXT=NN
              END IF
              IF(NNEXT-N.GT.1)THEN
C**UPDATE N3 FOR MISSING MODES
                DO NADJ=N+1,NNEXT-1
                  CALL INTARR(W(LNBF),W(LMBF),NADJ,IN1,IN2,IN3)
                  N3=N3+IN2
                END DO
              END IF
C**NEXT N
              N=NNEXT
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
              M=0
              M3=0
              DO MM=1,NN-1
                IF(LCONT.NE.0)THEN
                  MNEXT=JCONT(LCONT,MM)
                ELSE
                  MNEXT=MM
                END IF
                IF(MNEXT-M.GT.1)THEN
C**UPDATE M3 FOR MISSING MODES
                  DO MADJ=M+1,MNEXT-1
                    CALL INTARR(W(LNBF),W(LMBF),MADJ,IM1,IM2,IM3)
                    M3=M3+IM2
                  END DO
                END IF
C**NEXT M
                M=MNEXT
                CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
                J=0
                J3=0
                DO JJ=1,MM-1
                  IF(LCONT.NE.0)THEN
                    JNEXT=JCONT(LCONT,JJ)
                  ELSE
                    JNEXT=JJ
                  END IF
                  IF(JNEXT-J.GT.1)THEN
C**UPDATE J3 FOR MISSING MODES
                    DO JADJ=J+1,JNEXT-1
                      CALL INTARR(W(LNBF),W(LMBF),JADJ,IJ1,IJ2,IJ3)
                      J3=J3+IJ2
                    END DO
                  END IF
C**NEXT J
                  J=JNEXT
                  CALL INTARR(W(LNBF),W(LMBF),J,IJ1,IJ2,IJ3)
C**INTRINSIC
                  IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
                    CALL GROUP5(K,L,N,M,J,ICONT(1),JCONT,1,IND)
                    IF(IND.NE.0)THEN
                      JCOUPL=NCOUPL(1)
                      JCOUPC=NCOUPC(1)
                      ICOUPL=IABS(JCOUPL)
                      ICOUPC=IABS(JCOUPC)
                      GO TO 3035
                    END IF
                    CALL GROUP5(K,L,N,M,J,ICONT(2),JCONT,2,IND)
                    IF(IND.NE.0)THEN
                      JCOUPL=NCOUPL(2)
                      JCOUPC=NCOUPC(2)
                      ICOUPL=IABS(JCOUPL)
                      ICOUPC=IABS(JCOUPC)
                      GO TO 3035
                    END IF
                    JCOUPL=JCOUPS
                    JCOUPC=JCOUCS
                    ICOUPL=ICOUPS
                    ICOUPC=ICOUCS
3035                CONTINUE
                  END IF
                  IF(ICOUPL.LT.5)GO TO 3305
                  IF(LDUMP.EQ.0)THEN
C**WRITE POTENTIAL GRIDS TO DISC
                    IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.
     1              AND.N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN)
     2              .OR.(NSMODE.EQ.NAMODE))THEN
C**INTRINSIC
                      KVP=MBFMAX
                      IF(MOLINC.GT.0)
     1                CALL MEMO(4,LVP1,KVP*5,LVP2,KVP*KVP*10,
     2                LVP3,KVP*KVP*KVP*10,LVP4,KVP*KVP*KVP*KVP*5,0,0)
                      IENTMX(1)=KVP
                      IENTMX(2)=KVP
                      IENTMX(3)=KVP
                      IENTMX(4)=KVP
                      IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1                CALL DUMPT5(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     2                W(LXQ+M3),W(LXQ+J3),IK2,IL2,IN2,IM2,IJ2,NMODE,
     3                NATOM,W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     4                NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,L,N,M,J,W(LV5),
     5                W(LV5),1,W(LMODNT),W(LXTANH),W(LVP1),IENTMX(1),
     6                W(LVP2),IENTMX(2),W(LVP3),IENTMX(3),W(LVP4),
     7                IENTMX(4),W(LNBF),W(LMBF),W(LV1),W(LV1),W(LV2),
     8                W(LV2),W(LV3),W(LV3),W(LV4),W(LV4),W(LXMODQ))
                      IF(MOLINC.GT.0)
     1                CALL MEMO(-4,LVP1,KVP*5,LVP2,KVP*KVP*10,
     2                LVP3,KVP*KVP*KVP*10,LVP4,KVP*KVP*KVP*KVP*5,0,0)
                    ELSE
                      CALL DUMVT5(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1                W(LXQ+M3),W(LXQ+J3),IK2,IL2,IN2,IM2,IJ2,NMODE,
     2                NATOM,W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     3                K,L,N,M,J,W(LV5),W(LV5),W(LMODNT),IK2TAU,
     4                W(LXMODQ))
                    END IF
                    IF(ICOUPC.LT.5)GO TO 3305
C**WRITE CORIOLIS GRIDS TO DISC
                    IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1              N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN).OR.
     2              (NSMODE.EQ.NAMODE))THEN
C**INTRINSIC
                      CALL MEMO(4,LVP1,KVP*5,LVP2,KVP*KVP*10*6,
     1                LVP3,KVP*KVP*KVP*10*10,LVP4,KVP*KVP*KVP*KVP*15*5,
     2                0,0)
                      CALL DUMCR5(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1                W(LXQ+M3),W(LXQ+J3),IK2,IL2,IN2,IM2,IJ2,NMODE,
     2                NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),
     3                W(LXX),W(LX0),W(LXL),W(LXM),K,L,N,M,J,W(LC5),
     4                W(LC5),W(LVM5),W(LVM5),0,1,W(LMODNT),W(LNBF),
     5                W(LMBF),W(LVP1),W(LVP2),W(LVP3),W(LVP4),W(LV1),
C    6                W(LV1),W(LV2),W(LV2),W(LV3),W(LV3),W(LV4),W(LV4),
     6                W(LC1),W(LC2),W(LC2),W(LC3),W(LC3),W(LC4),W(LC4),
     7                KVP)
                      CALL MEMO(-4,LVP1,KVP*5,LVP2,KVP*KVP*10*6,
     1                LVP3,KVP*KVP*KVP*10*10,LVP4,KVP*KVP*KVP*KVP*15*5,
     2                0,0)
                    END IF
                  END IF
3305  CONTINUE
                  J3=J3+IJ2
                END DO
                M3=M3+IM2
              END DO
              N3=N3+IN2
            END DO
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
C******************************
C**INTRINSIC
C     IF(ICOUPL.EQ.5)GO TO 9999
      IF(ICOUPX.EQ.5)GO TO 9999
      MBFMX6=0
      DO KK=1,NAMODE
        IF(LCONT.NE.0)THEN
          K=JCONT(LCONT,KK)
        ELSE
          K=KK
        END IF
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        DO LL=1,KK-1
          IF(LCONT.NE.0)THEN
            L=JCONT(LCONT,LL)
          ELSE
            L=LL
          END IF
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
          DO NN=1,LL-1
            IF(LCONT.NE.0)THEN
              N=JCONT(LCONT,NN)
            ELSE
              N=NN
            END IF
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
            DO MM=1,NN-1
              IF(LCONT.NE.0)THEN
                M=JCONT(LCONT,MM)
              ELSE
                M=MM
              END IF
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
              DO JJ=1,MM-1
                IF(LCONT.NE.0)THEN
                  J=JCONT(LCONT,JJ)
                ELSE
                  J=JJ
                END IF
                CALL INTARR(W(LNBF),W(LMBF),J,IJ1,IJ2,IJ3)
                DO II=1,JJ-1
                  IF(LCONT.NE.0)THEN
                    I=JCONT(LCONT,II)
                  ELSE
                    I=II
                  END IF
                  CALL INTARR(W(LNBF),W(LMBF),I,II1,II2,II3)
                  IF(IK2*IL2*IN2*IM2*IJ2*II2.GT.MBFMX6)
     1            MBFMX6=IK2*IL2*IN2*IM2*IJ2*II2
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        KLV6=MBFMX6
      ELSE
        KLV6=(1+MBFMX6)/2
      END IF
      IF(JCOUPC.GE.0)THEN
        IF(ICOUCX.GE.6)THEN
          KLC6=MBFMX4*22
CC        IF(IREACT.NE.0)KLC6=MBFMX4*21
        END IF
      ELSE
        IF(ICOUCX.GE.6)THEN
          KLC6=(1+MBFMX4*22)/2
CC        IF(IREACT.NE.0)KLC6=(1+MBFMX4*21)/2
        END IF
      END IF
      IF(LCONT.EQ.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1  CALL MEMO(1,LV6,KLV6,0,0,0,0,0,0,0,0)
        IF(ICOUCX.GE.6)CALL MEMO(1,LC6,KLC6,0,0,0,0,0,0,0,0)
      END IF
      IF(IDISC.EQ.0)THEN
        K=0
        K3=0
        DO KK=1,NAMODE
          IF(LCONT.NE.0)THEN
            KNEXT=JCONT(LCONT,KK)
          ELSE
            KNEXT=KK
          END IF
          IF(KNEXT-K.GT.1)THEN
C**UPDATE K3 FOR MISSING MODES
            DO KADJ=K+1,KNEXT-1
              CALL INTARR(W(LNBF),W(LMBF),KADJ,IK1,IK2,IK3)
              K3=K3+IK2
            END DO
          END IF
C**NEXT K
          K=KNEXT
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L=0
          L3=0
          DO LL=1,KK-1
            IF(LCONT.NE.0)THEN
              LNEXT=JCONT(LCONT,LL)
            ELSE
              LNEXT=LL
            END IF
            IF(LNEXT-L.GT.1)THEN
C**UPDATE L3 FOR MISSING MODES
              DO LADJ=L+1,LNEXT-1
                CALL INTARR(W(LNBF),W(LMBF),LADJ,IL1,IL2,IL3)
                L3=L3+IL2
              END DO
            END IF
C**NEXT L
            L=LNEXT
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N=0
            N3=0
            DO NN=1,LL-1
              IF(LCONT.NE.0)THEN
                NNEXT=JCONT(LCONT,NN)
              ELSE
                NNEXT=NN
              END IF
              IF(NNEXT-N.GT.1)THEN
C**UPDATE N3 FOR MISSING MODES
                DO NADJ=N+1,NNEXT-1
                  CALL INTARR(W(LNBF),W(LMBF),NADJ,IN1,IN2,IN3)
                  N3=N3+IN2
                END DO
              END IF
C**NEXT N
              N=NNEXT
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
              M=0
              M3=0
              DO MM=1,NN-1
                IF(LCONT.NE.0)THEN
                  MNEXT=JCONT(LCONT,MM)
                ELSE
                  MNEXT=MM
                END IF
                IF(MNEXT-M.GT.1)THEN
C**UPDATE M3 FOR MISSING MODES
                  DO MADJ=M+1,MNEXT-1
                    CALL INTARR(W(LNBF),W(LMBF),MADJ,IM1,IM2,IM3)
                    M3=M3+IM2
                  END DO
                END IF
C**NEXT M
                M=MNEXT
                CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
                J=0
                J3=0
                DO JJ=1,MM-1
                  IF(LCONT.NE.0)THEN
                    JNEXT=JCONT(LCONT,JJ)
                  ELSE
                    JNEXT=JJ
                  END IF
                  IF(JNEXT-J.GT.1)THEN
C**UPDATE J3 FOR MISSING MODES
                    DO JADJ=J+1,JNEXT-1
                      CALL INTARR(W(LNBF),W(LMBF),JADJ,IJ1,IJ2,IJ3)
                      J3=J3+IJ2
                    END DO
                  END IF
C**NEXT J
                  J=JNEXT
                  CALL INTARR(W(LNBF),W(LMBF),J,IJ1,IJ2,IJ3)
                  I=0
                  I3=0
                  DO II=1,JJ-1
                    IF(LCONT.NE.0)THEN
                      INEXT=JCONT(LCONT,II)
                    ELSE
                      INEXT=II
                    END IF
                    IF(INEXT-I.GT.1)THEN
C**UPDATE I3 FOR MISSING MODES
                      DO IADJ=I+1,INEXT-1
                        CALL INTARR(W(LNBF),W(LMBF),IADJ,II1,II2,II3)
                        I3=I3+II2
                      END DO
                    END IF
C**NEXT I
                    I=INEXT
                    CALL INTARR(W(LNBF),W(LMBF),I,II1,II2,II3)
C**INTRINSIC
                    IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
                      CALL GROUP6(K,L,N,M,J,I,ICONT(1),JCONT,1,IND)
                      IF(IND.NE.0)THEN
                        JCOUPL=NCOUPL(1)
                        JCOUPC=NCOUPC(1)
                        ICOUPL=IABS(JCOUPL)
                        ICOUPC=IABS(JCOUPC)
                        GO TO 3036
                      END IF
                      CALL GROUP6(K,L,N,M,J,I,ICONT(2),JCONT,2,IND)
                      IF(IND.NE.0)THEN
                        JCOUPL=NCOUPL(2)
                        JCOUPC=NCOUPC(2)
                        ICOUPL=IABS(JCOUPL)
                        ICOUPC=IABS(JCOUPC)
                        GO TO 3036
                      END IF
                      JCOUPL=JCOUPS
                      JCOUPC=JCOUCS
                      ICOUPL=ICOUPS
                      ICOUPC=ICOUCS
3036                  CONTINUE
                    END IF
                    IF(ICOUPL.LT.6)GO TO 3306
                    IF(LDUMP.EQ.0)THEN
                      IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.
     1                AND.N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN.
     2                AND.I.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**INTRINSIC
                        KVP=MBFMAX
                        IF(MOLINC.GT.0)
     1                  CALL MEMO(5,LVP1,KVP*6,LVP2,KVP*KVP*15,
     2                  LVP3,KVP*KVP*KVP*20,LVP4,KVP*KVP*KVP*KVP*15,
     3                  LVP5,KVP*KVP*KVP*KVP*KVP*6)
                        IENTMX(1)=KVP
                        IENTMX(2)=KVP
                        IENTMX(3)=KVP
                        IENTMX(4)=KVP
                        IENTMX(5)=KVP
                        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1                  CALL DUMPT6(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     2                  W(LXQ+M3),W(LXQ+J3),W(LXQ+I3),IK2,IL2,IN2,IM2,
     3                  IJ2,II2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     4                  W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),
     5                  W(LCPOT),K,L,N,M,J,I,W(LV6),W(LV6),1,W(LMODNT),
     6                  W(LXTANH),W(LVP1),IENTMX(1),W(LVP2),IENTMX(2),
     7                  W(LVP3),IENTMX(3),W(LVP4),IENTMX(4),W(LVP5),
     8                  IENTMX(5),W(LNBF),W(LMBF),W(LV1),W(LV1),W(LV2),
     9                  W(LV2),W(LV3),W(LV3),W(LV4),W(LV4),W(LV5),
     1                  W(LV5),W(LXMODQ))
                        IF(MOLINC.GT.0)
     1                  CALL MEMO(-5,LVP1,KVP*6,LVP2,KVP*KVP*15,
     2                  LVP3,KVP*KVP*KVP*20,LVP4,KVP*KVP*KVP*KVP*15,
     3                  LVP5,KVP*KVP*KVP*KVP*KVP*6)
                      END IF
                      IF(ICOUPC.LT.6)GO TO 3306
C**WRITE CORIOLIS GRIDS TO DISC
                      IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.
     1                AND.N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN.
     2                AND.I.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**INTRINSIC
                        CALL MEMO(5,LVP1,KVP*6,LVP2,KVP*KVP*15*6,
     1                  LVP3,KVP*KVP*KVP*10*20,LVP4,
     2                  KVP*KVP*KVP*KVP*15*15,LVP5,
     3                  KVP*KVP*KVP*KVP*KVP*16*6)
                        REWIND 181
                        REWIND 182
                        REWIND 183
                        REWIND 184
                        REWIND 185
                        REWIND 186
                        CALL DUMCR6(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),
     1                  W(LXQ+M3),W(LXQ+J3),W(LXQ+I3),IK2,IL2,IN2,IM2,
     2                  IJ2,II2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     3                  W(LAA),W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,L,
     4                  N,M,J,I,W(LC6),W(LC6),W(LVM6),W(LVM6),0,1,
     5                  W(LMODNT),W(LNBF),W(LMBF),W(LVP1),W(LVP2),
C    6                  W(LVP3),W(LVP4),W(LVP5),W(LV1),W(LV1),W(LV2),
C    7                  W(LV2),W(LV3),W(LV3),W(LV4),W(LV4),W(LV5),
C    8                  W(LV5),KVP)
     6                  W(LVP3),W(LVP4),W(LVP5),W(LC1),W(LC1),W(LC2),
     7                  W(LC2),W(LC3),W(LC3),W(LC4),W(LC4),W(LC5),
     8                  W(LC5),KVP)
                        CALL MEMO(-5,LVP1,KVP*6,LVP2,KVP*KVP*15*6,
     1                  LVP3,KVP*KVP*KVP*10*20,LVP4,
     2                  KVP*KVP*KVP*KVP*15*15,LVP5,
     3                  KVP*KVP*KVP*KVP*KVP*16*6)
                      END IF
                    END IF
3306  CONTINUE
                    I3=I3+II2
                  END DO
                  J3=J3+IJ2
                END DO
                M3=M3+IM2
              END DO
              N3=N3+IN2
            END DO
            L3=L3+IL2
          END DO
          K3=K3+IK2
        END DO
      END IF
9999  CONTINUE
C     ICOUPL=ICOUPS
C     ICOUPC=ICOUCS
C**INTRINSIC
      IF(LCONT.EQ.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)THEN
          IF(ICOUPX.GT.1)THEN
            CALL MEMO(-1,LV1,KLVR1,0,0,0,0,0,0,0,0)
            CALL MEMO(1,LV1,KLV1,0,0,0,0,0,0,0,0)
          END IF
          IF(ICOUCX.GT.1)THEN
            CALL MEMO(-1,LC1,KLCR1,0,0,0,0,0,0,0,0)
            CALL MEMO(1,LC1,KLC1,0,0,0,0,0,0,0,0)
          END IF
          IF(ICOUPX.GT.2)THEN
            CALL MEMO(-1,LV2,KLVR2,0,0,0,0,0,0,0,0)
            CALL MEMO(1,LV2,KLV2,0,0,0,0,0,0,0,0)
          END IF
          IF(ICOUCX.GT.2)THEN
            CALL MEMO(-1,LC2,KLCR2,0,0,0,0,0,0,0,0)
            CALL MEMO(1,LC2,KLC2,0,0,0,0,0,0,0,0)
          END IF
          IF(ICOUPX.GT.3)THEN
            CALL MEMO(-1,LV3,KLVR3,0,0,0,0,0,0,0,0)
            CALL MEMO(1,LV3,KLV3,0,0,0,0,0,0,0,0)
          END IF
          IF(ICOUCX.GT.3)THEN
            CALL MEMO(-1,LC3,KLCR3,0,0,0,0,0,0,0,0)
            CALL MEMO(1,LC3,KLC3,0,0,0,0,0,0,0,0)
          END IF
        END IF
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE VRSCF(W,KCONT,NSTAT,CONV,K2TAU,K3TAU,MDTAU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      LOGICAL LINEAR
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
     1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
     4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
     5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
     6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
     7LDP4A,LDP4B,LDP4C,LDP4D,LXP,LJP6,LANCNT,LANK,LXJ,LVP5,
     8LC5,LC6,LVP0,LTMPRD,LMVB,LNFC,LTEMP6,LTEMP7,LTEMP8,LTEMP9
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/TYPE/LINEAR
      COMMON/PATH/ISCFCI
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/CYCLE/ICYCLE
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM5A,ITIM5B,ITIM6A,ITIM6B,ITIM
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/DEGEN/IDEGEN,NDEGEN
      COMMON/WHICH/IWHICH
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/MODES/NMODE,NATOM
      COMMON/VIBMOD/NSMODE,NVMODE,NNMODE
      COMMON/NORMOD/NAMODE,LINBND,NONLIN,MODD
      COMMON/ANALYT/MAXQU,MAXPOW
      COMMON/CSAVES/NNMODS,NAMODS,NVMODS,ICOUPS,ICOUCS,JREACS,NREACS
      COMMON/XSAVES/ICOUPX,ICOUCX,JCOUPS,JCOUCS
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/ENTER/IENTER,IENTMX(5),NTOT1,NTOT2,NTOT3,NTOT4,NTOT5
      COMMON/EVL/EVL,CUT
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU,IFLAUD
      COMMON/REACTL/JREACT,IFITRP
C**RPH ODD K
      COMMON/LMAX/LMAX
C**RPH ODD K
C**********************************************************************
215   FORMAT(/,1X,'SCF CALCULATION ONLY',/)
310   FORMAT(/,1X,'J = 0 SCF CYCLE',/)
315   FORMAT(/,1X,'J = ',I3,' SCF CYCLE',/)
320   FORMAT('*************************')
C**********************************************************************
      ITIM1A=0
      ITIM1B=0
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                                      VSCF PROCEDURE
C**********************************************************************
C**********************************************************************
C**********************************************************************
      IF(MOLINC.EQ.0)THEN
        ICOUPL=ICOUPX
        ICOUPC=ICOUCX
      END IF
C**TEMPORARY
      WRITE(IOUT,*)'VRSCF - KCONT = ',KCONT
      WRITE(IOUT,*)'ICOUPL,ICOUPC,NCOUPL(1),NCOUPC(1),NCOUPL(2),
     1NCOUPC(2)'
      WRITE(IOUT,*)ICOUPL,ICOUPC,NCOUPL(1),NCOUPC(1),NCOUPL(2),
     1NCOUPC(2)
      WRITE(IOUT,*)'JCOUPL,JCOUPC,ICOUPX,ICOUCX,ICOUPS,ICOUCS'
      WRITE(IOUT,*)JCOUPL,JCOUPC,ICOUPX,ICOUCX,ICOUPS,ICOUCS
      CALL FLUSH(IOUT)
C**TEMPORARY
CC    IF(KCONT.LT.0)ICCCS=ICOUPC
      ICCCS=ICOUPC
      DO 1500 ISTATE=1,NSTAT
C**********************************************************************
C**                                       LOOP OVER VSCF STATES IN TURN
C**********************************************************************
      ITIM=ITIM+1
C**DO SCF CYCLE FOR J=0 FIRST TIME
      J21=1
      JTHIS=0
      WRITE(IOUT,320)
C**RPH ODD K
      WRITE(IOUT,315)JTHIS
C**RPH ODD K
      WRITE(IOUT,320)
      CALL FLUSH(IOUT)
8000  CONTINUE
C**TEMPORARY-LIN
      REWIND 70
C**RPH ODD K
      LMAX=0
C**RPH ODD K
C**TEMPORARY-LIN
      IF(J21.EQ.1)E0=1.D+20
      KCORIG=JTHIS
      DO 7000 KROT=1,J21
C**********************************************************************
C**                                      LOOP OVER ROTATIONAL FUNCTIONS
C**********************************************************************
      count=0.0
C**KA,KC FOR VSCF
      KA=KROT/2
      KC=KCORIG-KA
C**RPH ODD K
      IF(KA.GT.0.AND.MOD(KA,2).EQ.MOD(KROT,2))KC=KC+1
C**RPH ODD K
C**DON'T DO BOTH KA>0 (BOTH IDENTICAL AT K-DIAGONAL EXCEPT Ka=1)
      IF(MOD(KROT,2).EQ.0)GO TO 7000
C**ONLY DO KA=0 IF INTEGER TAU THROUGHOUT RPH
C     IF(KROT.GT.1.AND.IFLAUD.EQ.1)GO TO 7000
C**CURRENTLY ONLY DO Ka=0 (INTEGER) AND Ka=1 (HALF-INTEGER)
C**TEMPORARY-LIN
C**RPH ODD K
      IF(KROT.GT.3)GO TO 7000
      LMAX=LMAX+1
C**RPH ODD K
      K2=0
      DO K=1,NAMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(K.GT.NONLIN)THEN
C**NEED TO REWRITE (70) AFTER Ka VSCF OPTIMISATION
          CALL IN70(W(LH+K2),IK3,IK2)
        END IF
        K2=K2+IK1
        IF(K.GT.NONLIN)K2=K2+IK1
      END DO
C**TEMPORARY-LIN
      ICYCLE=1
555   CONTINUE
C**********************************************************************
C**                               RETURN TO LABEL 555 UNTIL CONVERGENCE
C**********************************************************************
      ITIM=ITIM+1
C**FIRST STATE IS ZERO POINT LEVEL...DISPLAY ZERO POINT ENERGY
C**ENERGIES OF REMAINING STATES ARE RELATIVE TO ZERO POINT
C**OF PURE VIBRATIONAL LEVEL
      IF(ISTATE.EQ.1.AND.J21.EQ.1)EVL=0
      MNLK2=0
C     MNLK3=0
      REWIND 58
      REWIND 59
      DO 1000 MODE=1,NVMODE
C**********************************************************************
C**                 LOOP OVER ALL MODES....FORM NEW FUNCTIONS THIS MODE
C**********************************************************************
      CALL INTARR(W(LNBF),W(LMBF),MODE,IMODE1,IMODE2,IMODE3)
C**POINT TO L(0) OR L(1) IF LINEAR
      MLNK2=0
C     IF(MODE.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
CC      CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,MLNK2,IMODE1)
C       MLNK2=IMODE1
C     END IF
      KKA=MLNK2/IMODE1
      KXK=IMODE3
      INCT=MOD(IFLAUD,2)*MOD(KA,2)
C**EVEN Ka HAVE MISSING SIN(0.TAU)
      IF(JREACT.GT.0.AND.MODE.EQ.NSMODE)KXK=IMODE3-1+MOD(KA,2)-INCT
      CALL MEMO(3,LXK,KXK*KXK,LYK,KXK*KXK,LOV,KXK*KXK,0,0,0,0)
      CALL MTZERO(W(LXK),KXK)
      IF(ICOUPX.GE.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 71
      END IF
      IF(ICOUCX.GE.0)THEN
        IF(J21.GT.1)REWIND 61
        REWIND 81
      END IF
      IF(ICOUPX.GT.1)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 72
      END IF
      IF(ICOUCX.GT.1)THEN
        IF(J21.GT.1)REWIND 62
        REWIND 82
      END IF
      IF(ICOUPX.GT.2)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 73
      END IF
      IF(ICOUCX.GT.2)THEN
        IF(J21.GT.1)REWIND 63
        REWIND 83
      END IF
      IF(ICOUPX.GT.3)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 74
      END IF
      IF(ICOUCX.GT.3)THEN
        IF(J21.GT.1)REWIND 64
        REWIND 84
      END IF
      IF(ICOUPX.GT.4)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 75
      END IF
      IF(ICOUCX.GT.4)REWIND 85
      IF(ICOUPX.GT.5)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 76
      END IF
      IF(ICOUCX.GT.5)REWIND 86
      IF(IREACT.GT.0)THEN
        CALL INTARR(W(LNBF),W(LMBF),NSMODE,IK1TAU,IK2TAU,IK3TAU)
C**GET TORSION-ONLY INTEGRALS
        IF(JREACT.GT.0)THEN
          IF(MODE.EQ.NSMODE)THEN
            CALL TMIS0(W(LH+K2TAU),W(LXQ+K3TAU),IK3TAU,IK2TAU,W(LXK),
     1      KXK,NAMODE,W(LC0),W(LC0),W(LEJK0),
     2      W(LEJK0),J21,KROT,W(LMODNT))
          ELSE
            CALL TVAT0(W(LISTAT),NSTAT,NAMODE,ISTATE,
     1      W(LH+K2TAU),W(LXQ+K3TAU),IK3TAU,IK2TAU,
     2      W(LXK),KXK,W(LC0),W(LC0),W(LEJK0),W(LEJK0),
     3      J21,KROT,W(LMODNT))
          END IF
        END IF
      END IF
      IF(JREACT.LE.0.OR.NSMODE.EQ.NAMODE)
     1CALL THAT0(NAMODE,W(LXK),KXK,W(LEJK0),W(LEJK0),J21,KROT)
      IF(JREACT.GT.0)THEN
        IND=0
        CALL MEMO(1,LXP0,IK3TAU*IK3TAU*IK2TAU,0,0,0,0,0,0,0,0)
        IF(IFITRP.NE.0)
     1  CALL MEMO(3,LCP4,IENTMX(4)*MDTAU*NTOT4,LMP4,IENTMX(4)*4*NTOT4,
     2  LVP4,IENTMX(4)*NTOT4,0,0,0,0)
      END IF
C**INTRINSIC
      IF(KCONT.LT.0)THEN
        JCOUPS=JCOUPL
        JCOUCS=JCOUPC
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NAMODE
        IF(MODE.EQ.1.AND.K.EQ.1.AND.ITIM.EQ.0)THEN
          ITIM1A=0
          ITIM1B=0
        END IF
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**POINT TO L(0) OR L(1) IF LINEAR
        MK2=0
C       IF(K.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
CC        CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,MK2,IK1)
C         MK2=IK1
C       END IF
C******************************************************************
C**POSSIBLY LINEAR (JREACT<0)
        IF(JREACT.LE.0.AND.ICOUPX.EQ.0)THEN
C**'KINETIC ENERGY'
          IF(K.EQ.MODE)THEN
            CALL THISKE(W(LH+K2+MK2),W(LXQ+K3),IK3,IK2,K,W(LXK),
     1      W(LMODNT),W(LOMEGA+K1))
          ELSE
            CALL THATKE(W(LISTAT),NSTAT,NMODE,ISTATE,K,W(LH+K2+MK2),
     1      W(LXQ+K3),IK3,IK2,W(LXK),KXK,W(LMODNT),W(LOMEGA+K1))
          END IF
        END IF
C**INTRINSIC
        IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
          CALL GROUP1(K,ICONT(1),JCONT,1,IND)
          IF(IND.NE.0)THEN
            JCOUPL=NCOUPL(1)
            JCOUPC=NCOUPC(1)
            ICOUPL=IABS(JCOUPL)
            ICOUPC=IABS(JCOUPC)
            GO TO 3031
          END IF
          CALL GROUP1(K,ICONT(2),JCONT,2,IND)
          IF(IND.NE.0)THEN
            JCOUPL=NCOUPL(2)
            JCOUPC=NCOUPC(2)
            ICOUPL=IABS(JCOUPL)
            ICOUPC=IABS(JCOUPC)
            GO TO 3031
          END IF
          JCOUPL=JCOUPS
          JCOUPC=JCOUCS
          ICOUPL=ICOUPS
          ICOUPC=ICOUCS
3031      CONTINUE
CC        IF(KCONT.LT.0.AND.MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
          IF(MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
     1    ICOUPC=ICCCS
        END IF
        IF(ICOUPL.LT.1)GO TO 5501
C******************************************************************
C**CORIOLIS AND POTENTIAL
        IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
          IF(K.EQ.MODE)THEN
            CALL THIS1(W(LH+K2+MK2),W(LXQ+K3),IK3,IK2,W(LXK),
     1      NAMODE,K,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),W(LEJK1),J21,
     2      KROT,W(LMODNT),W(LOMEGA+K1),W(LXKAN),MAXQU,MAXPOW,W(LNP1),
     3      W(LCP1),W(LMP1),NTOT1,IENTMX(1))
          ELSE
            CALL THAT1(W(LISTAT),NSTAT,NAMODE,ISTATE,K,W(LH+K2+MK2),
     1      W(LXQ+K3),IK3,IK2,W(LXK),KXK,W(LV1),W(LV1),W(LC1),W(LC1),
     2      W(LEJK1),W(LEJK1),J21,KROT,W(LMODNT),W(LOMEGA+K1),W(LXKAN),
     3      MAXQU,MAXPOW,W(LNP1),W(LCP1),W(LMP1),NTOT1,IENTMX(1))
          END IF
        ELSE
C**IF REACT, ALWAYS INTEGRATE OVER TAU
          CALL MEMO(1,LXK0,IK2*7,0,0,0,0,0,0,0,0)
          IF(K.EQ.MODE)THEN
            CALL TVIS1(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1      W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IK3TAU,IK2TAU,W(LXK0),
     2      W(LXK),NAMODE,ISTATE,K,W(LV1),W(LV1),W(LC1),W(LC1),
     3      W(LEJK1),W(LEJK1),J21,KROT,W(LMODNT),W(LOMEGA+K1))
          ELSE
            IF(MODE.EQ.NSMODE)THEN
              CALL TMIS1(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1        W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IK3TAU,IK2TAU,W(LXK0),
     2        W(LXK),KXK,NAMODE,ISTATE,K,W(LV1),W(LV1),W(LC1),W(LC1),
     3        W(LEJK1),W(LEJK1),J21,KROT,W(LMODNT))
            ELSE
              CALL TVAT1(W(LISTAT),NSTAT,NAMODE,ISTATE,K,W(LH+K2+MK2),
     1        W(LXQ+K3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IK3TAU,IK2TAU,
     2        W(LXK0),W(LXK),KXK,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
     3        W(LEJK1),J21,KROT,W(LMODNT),W(LOMEGA+K1))
            END IF
          END IF
          CALL MEMO(-1,LXK0,IK2*7,0,0,0,0,0,0,0,0)
        END IF
C**INTRINSIC
5501    IF(ICOUPX.EQ.1)GO TO 501
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(MODE.EQ.1.AND.L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)THEN
            ITIM2A=0
            ITIM2B=0
          END IF
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**POINT TO L(0) OR L(1) IF LINEAR
          ML2=0
C         IF(L.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
CC          CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,ML2,IL1)
C           ML2=IL1
C         END IF
C**INTRINSIC
          IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
            CALL GROUP2(K,L,ICONT(1),JCONT,1,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(1)
              JCOUPC=NCOUPC(1)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 3032
            END IF
            CALL GROUP2(K,L,ICONT(2),JCONT,2,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(2)
              JCOUPC=NCOUPC(2)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 3032
            END IF
            JCOUPL=JCOUPS
            JCOUPC=JCOUCS
            ICOUPL=ICOUPS
            ICOUPC=ICOUCS
3032        CONTINUE
CC          IF(KCONT.LT.0.AND.MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
            IF(MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
     1      ICOUPC=ICCCS
          END IF
          IF(ICOUPL.LT.2)GO TO 5502
C**CORIOLIS AND POTENTIAL
          IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1    (NSMODE.EQ.NAMODE))THEN
            CALL MEMO(2,LXK0,IK2*4,LXL0,IL2*4,0,0,0,0,0,0)
            IF(L.EQ.MODE)THEN
              CALL THIS2B(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1        W(LH+L2+ML2),
     1        W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXK0),W(LXL0),W(LXK),NAMODE,
     2        ISTATE,K,L,W(LV2),W(LV2),W(LC2),W(LC2),W(LEJK2),W(LEJK2),
     3        J21,KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,W(LNP2),W(LCP2),
     4        W(LMP2),NTOT2,IENTMX(2),W(LINDK))
            ELSE
              IF(K.EQ.MODE)THEN
                CALL THIS2A(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1          W(LH+L2+ML2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXK0),W(LXL0),
     2          W(LXK),NAMODE,ISTATE,K,L,W(LV2),W(LV2),W(LC2),W(LC2),
     3          W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),W(LXKAN),MAXQU,
     4          MAXPOW,W(LNP2),W(LCP2),W(LMP2),NTOT2,IENTMX(2),
     5          W(LINDK))
              ELSE
                CALL THAT2(W(LISTAT),NSTAT,NAMODE,ISTATE,K,L,
     1          W(LH+K2+MK2),W(LXQ+K3),W(LH+L2+ML2),W(LXQ+L3),IK3,
     2          IK2,IL3,IL2,W(LXK0),
     2          W(LXL0),W(LXK),KXK,W(LV2),W(LV2),W(LC2),W(LC2),
     3          W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),W(LXKAN),MAXQU,
     4          MAXPOW,W(LNP2),W(LCP2),W(LMP2),NTOT2,IENTMX(2),
     5          W(LINDK))
              END IF
            END IF
            CALL MEMO(-2,LXK0,IK2*4,LXL0,IL2*4,0,0,0,0,0,0)
          ELSE
C**IF REACT, ALWAYS INTEGRATE OVER TAU
            CALL MEMO(2,LXK0,IK2*11,LXL0,IL2*11,0,0,0,0,0,0)
            IF(L.EQ.MODE)THEN
              CALL TVIS2B(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1        W(LH+L2+ML2),W(LXQ+L3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,
     2        IL3,IL2,IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXK),NAMODE,
     3        ISTATE,K,L,W(LV2),W(LV2),W(LC2),W(LC2),W(LEJK2),W(LEJK2),
     4        J21,KROT,W(LMODNT))
            ELSE
              IF(K.EQ.MODE)THEN
                CALL TVIS2A(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1          W(LH+L2+ML2),W(LXQ+L3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,
     2          IK2,IL3,IL2,IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXK),
     3          NAMODE,ISTATE,K,L,W(LV2),W(LV2),W(LC2),W(LC2),W(LEJK2),
     4          W(LEJK2),J21,KROT,W(LMODNT))
              ELSE
                IF(MODE.EQ.NSMODE)THEN
                  CALL TMIS2(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1            W(LH+L2+ML2),W(LXQ+L3),W(LH+K2TAU),W(LXQ+K3TAU),
     2            IK3,IK2,
     2            IL3,IL2,IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXK),KXK,
     3            NAMODE,ISTATE,K,L,W(LV2),W(LV2),W(LC2),W(LC2),
     4            W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT))
                ELSE
                  CALL TVAT2(W(LISTAT),NSTAT,NAMODE,ISTATE,K,L,
     1            W(LH+K2+MK2),W(LXQ+K3),W(LH+L2+ML2),W(LXQ+L3),
     2            W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IK3TAU,
     3            IK2TAU,W(LXK0),W(LXL0),W(LXK),KXK,W(LV2),W(LV2),
     4            W(LC2),W(LC2),W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT))
                END IF
              END IF
            END IF
            CALL MEMO(-2,LXK0,IK2*11,LXL0,IL2*11,0,0,0,0,0,0)
          END IF
C**INTRINSIC
5502      IF(ICOUPX.EQ.2)GO TO 502
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(MODE.EQ.1.AND.N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.
     1      ITIM.EQ.0)THEN
              ITIM3A=0
              ITIM3B=0
            END IF
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**POINT TO L(0) OR L(1) IF LINEAR
            MN2=0
C           IF(N.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
CC            CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,MN2,IN1)
C             MN2=IN1
C           END IF
C**INTRINSIC
            IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
              CALL GROUP3(K,L,N,ICONT(1),JCONT,1,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(1)
                JCOUPC=NCOUPC(1)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 3033
              END IF
              CALL GROUP3(K,L,N,ICONT(2),JCONT,2,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(2)
                JCOUPC=NCOUPC(2)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 3033
              END IF
              JCOUPL=JCOUPS
              JCOUPC=JCOUCS
              ICOUPL=ICOUPS
              ICOUPC=ICOUCS
3033          CONTINUE
CC            IF(KCONT.LT.0.AND.MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
              IF(MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
     1        ICOUPC=ICCCS
            END IF
            IF(ICOUPL.LT.3)GO TO 5503
C**CORIOLIS AND POTENTIAL
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1      N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
              CALL MEMO(3,LXK0,IK2*7,LXL0,IL2*7,LXN0,IN2*7,0,0,0,0)
              IF(N.EQ.MODE)THEN
                CALL THIS3C(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1          W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),IK3,IK2,
     2          IL3,IL2,
     2          IN3,IN2,W(LXK0),W(LXL0),W(LXN0),W(LXK),NAMODE,ISTATE,K,
     3          L,N,W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,
     4          KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,W(LNP3),W(LCP3),
     5          W(LMP3),NTOT3,IENTMX(3),W(LINDK),W(LINDL))
              ELSE
                IF(L.EQ.MODE)THEN
                  CALL THIS3B(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1            W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),IK3,
     2            IK2,IL3,
     2            IL2,IN3,IN2,W(LXK0),W(LXL0),W(LXN0),W(LXK),NAMODE,
     3            ISTATE,K,L,N,W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),
     4            W(LEJK3),J21,KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,
     5            W(LNP3),W(LCP3),W(LMP3),NTOT3,IENTMX(3),W(LINDK),
     6            W(LINDL))
                ELSE
                  IF(K.EQ.MODE)
     1            CALL THIS3A(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     2            W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),IK3,
     2            IK2,IL3,
     3            IL2,IN3,IN2,W(LXK0),W(LXL0),W(LXN0),W(LXK),NAMODE,
     4            ISTATE,K,L,N,W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),
     5            W(LEJK3),J21,KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,
     6            W(LNP3),W(LCP3),W(LMP3),NTOT3,IENTMX(3),W(LINDK),
     7            W(LINDL))
                  IF(K.NE.MODE)
     1            CALL THAT3(W(LISTAT),NSTAT,NAMODE,ISTATE,K,L,N,
     2            W(LH+K2+MK2),W(LXQ+K3),W(LH+L2+ML2),W(LXQ+L3),
     3            W(LH+N2+MN2),
     3            W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,IN2,W(LXK0),W(LXL0),
     4            W(LXN0),W(LXK),KXK,W(LV3),W(LV3),W(LC3),W(LC3),
     5            W(LEJK3),W(LEJK3),J21,KROT,W(LMODNT),W(LXKAN),MAXQU,
     6            MAXPOW,W(LNP3),W(LCP3),W(LMP3),NTOT3,IENTMX(3),
     7            W(LINDK),W(LINDL))
                END IF
              END IF
              CALL MEMO(-3,LXK0,IK2*7,LXL0,IL2*7,LXN0,IN2*7,0,0,0,0)
            ELSE
C**IF REACT, ALWAYS INTEGRATE OVER TAU
              CALL MEMO(3,LXK0,IK2*16,LXL0,IL2*16,LXN0,IN2*16,0,0,0,0)
              IF(N.EQ.MODE)THEN
                CALL TVIS3C(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1          W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),
     2          W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,
     3          IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXN0),W(LXK),NAMODE,
     4          ISTATE,K,L,N,W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),
     5          W(LEJK3),J21,KROT,W(LMODNT))
              ELSE
                IF(L.EQ.MODE)THEN
                  CALL TVIS3B(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1            W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),
     2            W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,
     3            IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXN0),W(LXK),NAMODE,
     4            ISTATE,K,L,N,W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),
     5            W(LEJK3),J21,KROT,W(LMODNT))
                ELSE
                  IF(K.EQ.MODE)
     1            CALL TVIS3A(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     2            W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),
     3            W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,
     4            IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXN0),W(LXK),NAMODE,
     5            ISTATE,K,L,N,W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),
     6            W(LEJK3),J21,KROT,W(LMODNT))
                  IF(K.NE.MODE)THEN
                    IF(MODE.EQ.NSMODE)
     1              CALL TMIS3(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     2              W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),
     2              W(LH+K2TAU),
     3              W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,IK3TAU,IK2TAU,
     4              W(LXK0),W(LXL0),W(LXN0),W(LXK),KXK,NAMODE,ISTATE,K,
     5              L,N,W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),
     6              J21,KROT,W(LMODNT))
                    IF(MODE.NE.NSMODE)
     1              CALL TVAT3(W(LISTAT),NSTAT,NAMODE,ISTATE,K,L,N,
     2              W(LH+K2+MK2),W(LXQ+K3),W(LH+L2+ML2),W(LXQ+L3),
     3              W(LH+N2+MN2),W(LXQ+N3),W(LH+K2TAU),W(LXQ+K3TAU),
     4              IK3,IK2,IL3,IL2,IN3,IN2,IK3TAU,IK2TAU,W(LXK0),
     5              W(LXL0),W(LXN0),W(LXK),KXK,W(LV3),W(LV3),W(LC3),
     6              W(LC3),W(LEJK3),W(LEJK3),J21,KROT,W(LMODNT))
                  END IF
                END IF
              END IF
              CALL MEMO(-3,LXK0,IK2*16,LXL0,IL2*16,LXN0,IN2*16,0,0,0,0)
            END IF
C**INTRINSIC
5503        IF(ICOUPX.EQ.3)GO TO 503
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(MODE.EQ.1.AND.M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.
     1        K.EQ.4.AND.ITIM.EQ.0)THEN
                ITIM4A=0
                ITIM4B=0
              END IF
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**POINT TO L(0) OR L(1) IF LINEAR
              MM2=0
C             IF(M.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
CC              CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,MM2,IM1)
C               MM2=IM1
C             END IF
C**INTRINSIC
              IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
                CALL GROUP4(K,L,N,M,ICONT(1),JCONT,1,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(1)
                  JCOUPC=NCOUPC(1)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 3034
                END IF
                CALL GROUP4(K,L,N,M,ICONT(2),JCONT,2,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(2)
                  JCOUPC=NCOUPC(2)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 3034
                END IF
                JCOUPL=JCOUPS
                JCOUPC=JCOUCS
                ICOUPL=ICOUPS
                ICOUPC=ICOUCS
3034            CONTINUE
CC              IF(KCONT.LT.0.AND.MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
                IF(MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
     1          ICOUPC=ICCCS
              END IF
              IF(ICOUPL.LT.4)GO TO 5504
C**CORIOLIS AND POTENTIAL
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2        THEN
                CALL MEMO(4,LXK0,IK2*11,LXL0,IL2*11,LXN0,IN2*11,LXM0,
     1          IM2*11,0,0)
                IF(M.EQ.MODE)THEN
                  CALL THIS4D(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1            W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2            W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXK0),
     3            W(LXL0),W(LXN0),W(LXM0),W(LXK),NAMODE,ISTATE,K,L,N,M,
     4            W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,
     5            KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,W(LNP4),W(LCP4),
     6            W(LMP4),NTOT4,IENTMX(4),W(LINDK),W(LINDL),W(LINDN))
                ELSE
                  IF(N.EQ.MODE)THEN
                    CALL THIS4C(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1              W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2              W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXK0),
     3              W(LXL0),W(LXN0),W(LXM0),W(LXK),NAMODE,ISTATE,K,L,N,
     4              M,W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
     5              J21,KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,W(LNP4),
     6              W(LCP4),W(LMP4),NTOT4,IENTMX(4),W(LINDK),W(LINDL),
     7              W(LINDN))
                  ELSE
                    IF(L.EQ.MODE)THEN
                      CALL THIS4B(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1                W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2                W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,
     3                W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXK),NAMODE,
     4                ISTATE,K,L,N,M,W(LV4),W(LV4),W(LC4),W(LC4),
     5                W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT),W(LXKAN),
     6                MAXQU,MAXPOW,W(LNP4),W(LCP4),W(LMP4),NTOT4,
     7                IENTMX(4),W(LINDK),W(LINDL),W(LINDN))
                    ELSE
                      IF(K.EQ.MODE)
     1                CALL THIS4A(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     2                W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     3                W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,
     4                W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXK),NAMODE,
     5                ISTATE,K,L,N,M,W(LV4),W(LV4),W(LC4),W(LC4),
     6                W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT),W(LXKAN),
     7                MAXQU,MAXPOW,W(LNP4),W(LCP4),W(LMP4),NTOT4,
     8                IENTMX(4),W(LINDK),W(LINDL),W(LINDN))
                      IF(K.NE.MODE)
     1                CALL THAT4(W(LISTAT),NSTAT,NAMODE,ISTATE,K,L,N,M,
     2                W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     3                W(LXQ+N3),W(LH+M2),W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,
     4                IN2,IM3,IM2,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     5                W(LXK),KXK,W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),
     6                W(LEJK4),J21,KROT,W(LMODNT),W(LXKAN),MAXQU,
     7                MAXPOW,W(LNP4),W(LCP4),W(LMP4),NTOT4,IENTMX(4),
     8                W(LINDK),W(LINDL),W(LINDN))
                    END IF
                  END IF
                END IF
                CALL MEMO(-4,LXK0,IK2*11,LXL0,IL2*11,LXN0,IN2*11,LXM0,
     1          IM2*11,0,0)
              ELSE
C**IF REACT, ALWAYS INTEGRATE OVER TAU
                IF(IABS(IFITRP).LT.4)THEN
                  CALL MEMO(4,LXK0,IK2*22,LXL0,IL2*22,LXN0,IN2*22,LXM0,
     1            IM2*22,0,0)
                  IF(M.EQ.MODE)THEN
                    CALL TVIS4D(W(LISTAT),NSTAT,W(LH+K2+MK2),W(LXQ+K3),
     1              W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),
     2              W(LH+M2+MM2),W(LXQ+M3),W(LH+K2TAU),W(LXQ+K3TAU),
     3              IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,IK3TAU,IK2TAU,
     4              W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXK),NAMODE,
     5              ISTATE,K,L,N,M,W(LV4),W(LV4),W(LC4),W(LC4),
     6              W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT))
                  ELSE
                    IF(N.EQ.MODE)THEN
                      CALL TVIS4C(W(LISTAT),NSTAT,W(LH+K2+MK2),
     1                W(LXQ+K3),W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),
     2                W(LXQ+N3),W(LH+M2+MM2),W(LXQ+M3),W(LH+K2TAU),
     3                W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,
     4                IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     5                W(LXK),NAMODE,ISTATE,K,L,N,M,W(LV4),W(LV4),
     6                W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,KROT,
     7                W(LMODNT))
                    ELSE
                      IF(L.EQ.MODE)THEN
                        CALL TVIS4B(W(LISTAT),NSTAT,W(LH+K2+MK2),
     1                  W(LXQ+K3),W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),
     2                  W(LXQ+N3),W(LH+M2+MM2),W(LXQ+M3),W(LH+K2TAU),
     3                  W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,
     4                  IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     5                  W(LXK),NAMODE,ISTATE,K,L,N,M,W(LV4),W(LV4),
     6                  W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,KROT,
     7                  W(LMODNT))
                      ELSE
                        IF(K.EQ.MODE)
     1                  CALL TVIS4A(W(LISTAT),NSTAT,W(LH+K2+MK2),
     2                  W(LXQ+K3),W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),
     3                  W(LXQ+N3),W(LH+M2+MM2),W(LXQ+M3),W(LH+K2TAU),
     4                  W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,
     5                  IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6                  W(LXK),NAMODE,ISTATE,K,L,N,M,W(LV4),W(LV4),
     7                  W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,KROT,
     8                  W(LMODNT))
                        IF(K.NE.MODE)THEN
                          IF(MODE.EQ.NSMODE)
     1                    CALL TMIS4(W(LISTAT),NSTAT,W(LH+K2),
     2                    W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     3                    W(LXQ+N3),W(LH+M2),W(LXQ+M3),W(LH+K2TAU),
     4                    W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,
     5                    IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXN0),
     6                    W(LXM0),W(LXK),KXK,NAMODE,ISTATE,K,L,N,M,
     7                    W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),
     8                    W(LEJK4),J21,KROT,W(LMODNT))
                          IF(MODE.NE.NSMODE)
     1                    CALL TVAT4(W(LISTAT),NSTAT,NAMODE,ISTATE,K,L,
     2                    N,M,W(LH+K2+MK2),W(LXQ+K3),W(LH+L2+ML2),
     3                    W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),
     4                    W(LH+M2+MM2),W(LXQ+M3),W(LH+K2TAU),
     5                    W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,
     6                    IK3TAU,IK2TAU,W(LXK0),W(LXL0),W(LXN0),
     7                    W(LXM0),W(LXK),KXK,W(LV4),W(LV4),W(LC4),
     8                    W(LC4),W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT))
                        END IF
                      END IF
                    END IF
                  END IF
                  CALL MEMO(-4,LXK0,IK2*22,LXL0,IL2*22,LXN0,IN2*22,
     1            LXM0,IM2*22,0,0)
                ELSE
                  IF(M.EQ.MODE)THEN
                    CALL TVVS4D(W(LISTAT),NSTAT,IM3,W(LH+K2TAU),
     1              W(LXQ+K3TAU),IK3TAU,IK2TAU,W(LXK),ISTATE,K,L,N,M,
     2              KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,NAMODE,
     3              W(LNP4),W(LCP4),W(LMP4),W(LVP4),NTOT4,IENTMX(4),
     4              W(LINDK),W(LINDL),W(LINDN),IND,MDTAU)
                  ELSE
                    IF(N.EQ.MODE)THEN
                      CALL TVVS4C(W(LISTAT),NSTAT,IN3,W(LH+K2TAU),
     1                W(LXQ+K3TAU),IK3TAU,IK2TAU,W(LXK),ISTATE,K,L,N,M,
     2                KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,NAMODE,
     3                W(LNP4),W(LCP4),W(LMP4),W(LVP4),NTOT4,IENTMX(4),
     4                W(LINDK),W(LINDL),W(LINDN),IND,MDTAU)
                    ELSE
                      IF(L.EQ.MODE)THEN
                        CALL TVVS4B(W(LISTAT),NSTAT,IL3,W(LH+K2TAU),
     1                  W(LXQ+K3TAU),IK3TAU,IK2TAU,W(LXK),ISTATE,K,L,N,
     2                  M,KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,NAMODE,
     3                  W(LNP4),W(LCP4),W(LMP4),W(LVP4),NTOT4,
     4                  IENTMX(4),W(LINDK),W(LINDL),W(LINDN),IND,MDTAU)
                      ELSE
                        IF(K.EQ.MODE)
     1                  CALL TVVS4A(W(LISTAT),NSTAT,IK3,W(LH+K2TAU),
     2                  W(LXQ+K3TAU),IK3TAU,IK2TAU,W(LXK),ISTATE,K,L,N,
     3                  M,KROT,W(LMODNT),W(LXKAN),MAXQU,MAXPOW,NAMODE,
     4                  W(LNP4),W(LCP4),W(LMP4),W(LVP4),NTOT4,
     5                  IENTMX(4),W(LINDK),W(LINDL),W(LINDN),IND,MDTAU)
                        IF(K.NE.MODE)THEN
                          IF(MODE.EQ.NSMODE)THEN
                            CALL MEMO(1,LTEMP,IENTMX(4)*IK3TAU*IK3TAU,
     1                      0,0,0,0,0,0,0,0)
                            CALL TMVS4(W(LISTAT),NSTAT,W(LH+K2TAU),
     1                      W(LXQ+K3TAU),IK3TAU,IK2TAU,W(LXK),KXK,
     2                      ISTATE,K,L,N,M,KROT,W(LMODNT),W(LXKAN),
     3                      MAXQU,MAXPOW,NAMODE,W(LNP4),W(LCP4),
     4                      W(LMP4),W(LTEMP),NTOT4,IENTMX(4),W(LINDK),
     5                      W(LINDL),W(LINDN),IND,MDTAU,W(LXP0))
                            CALL MEMO(-1,LTEMP,IENTMX(4)*IK3TAU*IK3TAU,
     1                      0,0,0,0,0,0,0,0)
                          END IF
                          IF(MODE.NE.NSMODE)
     1                    CALL TVVT4(W(LISTAT),NSTAT,ISTATE,K,L,N,M,
     2                    W(LH+K2TAU),W(LXQ+K3TAU),IK3TAU,IK2TAU,
     3                    W(LXK),KXK,KROT,W(LMODNT),W(LXKAN),MAXQU,
     4                    MAXPOW,NAMODE,W(LNP4),W(LCP4),W(LMP4),
     5                    W(LVP4),NTOT4,IENTMX(4),W(LINDK),W(LINDL),
     6                    W(LINDN),IND,MDTAU)
                        END IF
                      END IF
                    END IF
                  END IF
                END IF
              END IF
C**INTRINSIC
5504          IF(ICOUPX.EQ.4)GO TO 504
C**INTEGRATE OVER FIVE NORMAL COORDINATES
              J2=0
              J3=0
              DO J=1,M-1
                IF(MODE.EQ.1.AND.J.EQ.1.AND.M.EQ.2.AND.N.EQ.3.AND.
     1          L.EQ.4.AND.K.EQ.5.AND.ITIM.EQ.0)THEN
                  ITIM5A=0
                  ITIM5B=0
                END IF
                CALL INTARR(W(LNBF),W(LMBF),J,IJ1,IJ2,IJ3)
C**POINT TO L(0) OR L(1) IF LINEAR
                MJ2=0
C               IF(J.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
CC                CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,MJ2,IJ1)
C                 MJ2=IJ1
C               END IF
C**INTRINSIC
                IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
                  CALL GROUP5(K,L,N,M,J,ICONT(1),JCONT,1,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(1)
                    JCOUPC=NCOUPC(1)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 3035
                  END IF
                  CALL GROUP5(K,L,N,M,J,ICONT(2),JCONT,2,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(2)
                    JCOUPC=NCOUPC(2)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 3035
                  END IF
                  JCOUPL=JCOUPS
                  JCOUPC=JCOUCS
                  ICOUPL=ICOUPS
                  ICOUPC=ICOUCS
3035              CONTINUE
CC                IF(KCONT.LT.0.AND.MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
                  IF(MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
     1            ICOUPC=ICCCS
                END IF
                IF(ICOUPL.LT.5)GO TO 5505
C**CORIOLIS AND POTENTIAL
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN).OR.
     2          (NSMODE.EQ.NAMODE))THEN
                  CALL MEMO(5,LXK0,IK2,LXL0,IL2,LXN0,IN2,LXM0,
     1            IM2,LXJ0,IJ2)
                  IF(J.EQ.MODE)THEN
                    CALL THIS5E(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1              W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2              W(LXQ+M3),W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,
     3              IN2,IM3,IM2,IJ3,IJ2,W(LXK0),W(LXL0),W(LXN0),
     4              W(LXM0),W(LXJ0),W(LXK),NAMODE,ISTATE,K,L,N,M,J,
     5              W(LV5),W(LV5),W(LMODNT),W(LXKAN),MAXQU,MAXPOW,
     6              W(LNP5),W(LCP5),W(LMP5),NTOT5,MAX5,W(LINDK),
     7              W(LINDL),W(LINDN),W(LINDM))
                  ELSE
                    IF(M.EQ.MODE)THEN
                      CALL THIS5D(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1                W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2                W(LXQ+M3),W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,
     3                IN2,IM3,IM2,IJ3,IJ2,W(LXK0),W(LXL0),W(LXN0),
     4                W(LXM0),W(LXJ0),W(LXK),NAMODE,ISTATE,K,L,N,M,J,
     5                W(LV5),W(LV5),W(LMODNT),W(LXKAN),MAXQU,MAXPOW,
     6                W(LNP5),W(LCP5),W(LMP5),NTOT5,MAX5,W(LINDK),
     7                W(LINDL),W(LINDN),W(LINDM))
                    ELSE
                      IF(N.EQ.MODE)THEN
                        CALL THIS5C(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1                  W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2                  W(LXQ+M3),W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,
     3                  IN3,IN2,IM3,IM2,IJ3,IJ2,W(LXK0),W(LXL0),
     4                  W(LXN0),W(LXM0),W(LXJ0),W(LXK),NAMODE,ISTATE,K,
     5                  L,N,M,J,W(LV5),W(LV5),W(LMODNT),W(LXKAN),MAXQU,
     7                  MAXPOW,W(LNP5),W(LCP5),W(LMP5),NTOT5,MAX5,
     8                  W(LINDK),W(LINDL),W(LINDN),W(LINDM))
                      ELSE
                        IF(L.EQ.MODE)THEN
                          CALL THIS5B(W(LISTAT),NSTAT,W(LH+K2),
     1                    W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2                    W(LXQ+N3),W(LH+M2),W(LXQ+M3),W(LH+J2),
     3                    W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,
     4                    IJ3,IJ2,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     5                    W(LXJ0),W(LXK),NAMODE,ISTATE,K,L,N,M,J,
     6                    W(LV5),W(LV5),W(LMODNT),W(LXKAN),MAXQU,
     8                    MAXPOW,W(LNP5),W(LCP5),W(LMP5),NTOT5,MAX5,
     9                    W(LINDK),W(LINDL),W(LINDN),W(LINDM))
                        ELSE
                          IF(K.EQ.MODE)
     1                    CALL THIS5A(W(LISTAT),NSTAT,W(LH+K2),
     2                    W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     3                    W(LXQ+N3),W(LH+M2),W(LXQ+M3),W(LH+J2),
     4                    W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,
     5                    IJ3,IJ2,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6                    W(LXJ0),W(LXK),NAMODE,ISTATE,K,L,N,M,J,
     7                    W(LV5),W(LV5),W(LMODNT),W(LXKAN),MAXQU,
     9                    MAXPOW,W(LNP5),W(LCP5),W(LMP5),NTOT5,MAX5,
     1                    W(LINDK),W(LINDL),W(LINDN),W(LINDM))
                          IF(K.NE.MODE)
     1                    CALL THAT5(W(LISTAT),NSTAT,NAMODE,ISTATE,K,L,
     2                    N,M,J,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),
     3                    W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     4                    W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,IN2,
     5                    IM3,IM2,IJ3,IJ2,W(LXK0),W(LXL0),W(LXN0),
     6                    W(LXM0),W(LXJ0),W(LXK),KXK,W(LV5),W(LV5),
     8                    W(LMODNT),W(LXKAN),MAXQU,MAXPOW,W(LNP5),
     9                    W(LCP5),W(LMP5),NTOT5,MAX5,W(LINDK),W(LINDL),
     1                    W(LINDN),W(LINDM))
                        END IF
                      END IF
                    END IF
                  END IF
                  CALL MEMO(-5,LXK0,IK2,LXL0,IL2,LXN0,IN2,LXM0,
     1            IM2,LXJ0,IJ2)
                ELSE
C**REACTION PATH
                END IF
5505            IF(ICOUPX.EQ.5)GO TO 505
C**INTEGRATE OVER SIX NORMAL COORDINATES
                I2=0
                I3=0
                DO I=1,J-1
                  IF(MODE.EQ.1.AND.I.EQ.1.AND.J.EQ.2.AND.M.EQ.3.AND.
     1            N.EQ.4.AND.L.EQ.5.AND.K.EQ.6.AND.ITIM.EQ.0)THEN
                    ITIM6A=0
                    ITIM6B=0
                  END IF
                  CALL INTARR(W(LNBF),W(LMBF),I,II1,II2,II3)
C**POINT TO L(0) OR L(1) IF LINEAR
                  MI2=0
C                 IF(I.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
CC                  CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,MI2,II1)
C                   MI2=II1
C                 END IF
C**INTRINSIC
                  IF(KCONT.LT.0.OR.MOLINC.GT.0)THEN
                    CALL GROUP6(K,L,N,M,J,I,ICONT(1),JCONT,1,IND)
                    IF(IND.NE.0)THEN
                      JCOUPL=NCOUPL(1)
                      JCOUPC=NCOUPC(1)
                      ICOUPL=IABS(JCOUPL)
                      ICOUPC=IABS(JCOUPC)
                      GO TO 3036
                    END IF
                    CALL GROUP6(K,L,N,M,J,I,ICONT(2),JCONT,2,IND)
                    IF(IND.NE.0)THEN
                      JCOUPL=NCOUPL(2)
                      JCOUPC=NCOUPC(2)
                      ICOUPL=IABS(JCOUPL)
                      ICOUPC=IABS(JCOUPC)
                      GO TO 3036
                    END IF
                    JCOUPL=JCOUPS
                    JCOUPC=JCOUCS
                    ICOUPL=ICOUPS
                    ICOUPC=ICOUCS
3036                CONTINUE
CC                  IF(KCONT.LT.0.AND.MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
                    IF(MOLINC.EQ.0.AND.ICOUPC.GT.ICCCS)
     1              ICOUPC=ICCCS
                  END IF
                  IF(ICOUPL.LT.6)GO TO 5506
C**CORIOLIS AND POTENTIAL
                  IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1            N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN.AND.
     2            I.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
                    CALL MEMO(5,LXK0,IK2,LXL0,IL2,LXN0,IN2,LXM0,
     1              IM2,LXJ0,IJ2)
                    CALL MEMO(1,LXI0,II2,0,0,0,0,0,0,0,0)
                    IF(I.EQ.MODE)THEN
                      CALL THIS6F(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1                W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2                W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+I2),W(LXQ+I3),
     3                IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,II3,II2,
     4                W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),W(LXI0),
     5                W(LXK),NAMODE,ISTATE,K,L,N,M,J,I,W(LV6),W(LV6),
     6                W(LMODNT))
                    ELSE
                      IF(J.EQ.MODE)THEN
                        CALL THIS6E(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1                  W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2                  W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+I2),
     3                  W(LXQ+I3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,IJ3,
     4                  IJ2,II3,II2,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     5                  W(LXJ0),W(LXI0),W(LXK),NAMODE,ISTATE,K,L,N,M,J,
     6                  I,W(LV6),W(LV6),W(LMODNT))
                      ELSE
                        IF(M.EQ.MODE)THEN
                          CALL THIS6D(W(LISTAT),NSTAT,W(LH+K2),
     1                    W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2                    W(LXQ+N3),W(LH+M2),W(LXQ+M3),W(LH+J2),
     3                    W(LXQ+J3),W(LH+I2),W(LXQ+I3),IK3,IK2,IL3,IL2,
     4                    IN3,IN2,IM3,IM2,IJ3,IJ2,II3,II2,W(LXK0),
     5                    W(LXL0),W(LXN0),W(LXM0),W(LXJ0),W(LXI0),
     6                    W(LXK),NAMODE,ISTATE,K,L,N,M,J,I,W(LV6),
     7                    W(LV6),W(LMODNT))
                        ELSE
                          IF(N.EQ.MODE)THEN
                            CALL THIS6C(W(LISTAT),NSTAT,W(LH+K2),
     1                      W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2                      W(LXQ+N3),W(LH+M2),W(LXQ+M3),W(LH+J2),
     3                      W(LXQ+J3),W(LH+I2),W(LXQ+I3),IK3,IK2,IL3,
     4                      IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,II3,II2,
     5                      W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),
     6                      W(LXI0),W(LXK),NAMODE,ISTATE,K,L,N,M,J,I,
     7                      W(LV6),W(LV6),W(LMODNT))
                          ELSE
                            IF(L.EQ.MODE)THEN
                              CALL THIS6B(W(LISTAT),NSTAT,W(LH+K2),
     1                        W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2                        W(LXQ+N3),W(LH+M2),W(LXQ+M3),W(LH+J2),
     3                        W(LXQ+J3),W(LH+I2),W(LXQ+I3),IK3,IK2,IL3,
     4                        IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,II3,II2,
     5                        W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),
     6                        W(LXI0),W(LXK),NAMODE,ISTATE,K,L,N,M,J,I,
     7                        W(LV6),W(LV6),W(LMODNT))
                            ELSE
                              IF(K.EQ.MODE)
     1                        CALL THIS6A(W(LISTAT),NSTAT,W(LH+K2),
     2                        W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     3                        W(LXQ+N3),W(LH+M2),W(LXQ+M3),W(LH+J2),
     4                        W(LXQ+J3),W(LH+I2),W(LXQ+I3),IK3,IK2,IL3,
     5                        IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,II3,II2,
     6                        W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),
     7                        W(LXI0),W(LXK),NAMODE,ISTATE,K,L,N,M,J,I,
     8                        W(LV6),W(LV6),W(LMODNT))
                              IF(K.NE.MODE)
     1                        CALL THAT6(W(LISTAT),NSTAT,NAMODE,ISTATE,
     2                        K,L,N,M,J,I,W(LH+K2),W(LXQ+K3),W(LH+L2),
     3                        W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     4                        W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+I2),
     5                        W(LXQ+I3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,
     6                        IM2,IJ3,IJ2,II3,II2,W(LXK0),W(LXL0),
     7                        W(LXN0),W(LXM0),W(LXJ0),W(LXI0),W(LXK),
     8                        KXK,W(LV6),W(LV6),W(LMODNT))
                            END IF
                          END IF
                        END IF
                      END IF
                    END IF
                    CALL MEMO(-5,LXK0,IK2,LXL0,IL2,LXN0,IN2,LXM0,
     1              IM2,LXJ0,IJ2)
                    CALL MEMO(-1,LXI0,II2,0,0,0,0,0,0,0,0)
                  ELSE
C**REACTION PATH
                  END IF
5506  CONTINUE
                  I2=I2+II1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
                  IF(I.GT.NONLIN)I2=I2+II1
                  I3=I3+II2
                END DO
505   CONTINUE
                J2=J2+IJ1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
                IF(J.GT.NONLIN)J2=J2+IJ1
                J3=J3+IJ2
              END DO
504   CONTINUE
              M2=M2+IM1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
              IF(M.GT.NONLIN)M2=M2+IM1
              M3=M3+IM2
            END DO
503   CONTINUE
            N2=N2+IN1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
            IF(N.GT.NONLIN)N2=N2+IN1
            N3=N3+IN2
          END DO
502   CONTINUE
          L2=L2+IL1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(L.GT.NONLIN)L2=L2+IL1
          L3=L3+IL2
        END DO
501   CONTINUE
        K2=K2+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
        IF(K.GT.NONLIN)K2=K2+IK1
        K3=K3+IK2
      END DO
500   CONTINUE
CC    IF(KCONT.LT.0)ICOUPC=ICCCS
      ICOUPC=ICCCS
      IF(JREACT.GT.0)THEN
        IND=1
        CALL MEMO(-1,LXP0,IK3TAU*IK3TAU*IK2TAU,0,0,0,0,0,0,0,0)
        IF(IFITRP.NE.0)
     1  CALL MEMO(-3,LCP4,IENTMX(4)*MDTAU*NTOT4,LMP4,IENTMX(4)*4*NTOT4,
     2  LVP4,IENTMX(4)*NTOT4,0,0,0,0)
      END IF
      IF(IABS(IFITRP).EQ.4)REWIND 24
C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
C**DIAGONALISE FOR THIS MODE..NO NEED FOR ENERGIES UNTIL FINAL MODE
      KWK=KXK
      KEVAL=KXK
      CALL MEMO(3,LWK,KWK,LEVAL,KEVAL,LTEMP,KEVAL,0,0,0,0)
      IF(ITIM.EQ.1.AND.ITIM1B.EQ.2)THEN
        WRITE(IOUT,*)'Calculating DIAG'
        CALL TIMIT(1)
      END IF
C     CALL INTARR(W(LMODNT),W(LMODNT),MODE,IDUM,IDUM,MD)
      DO I=1,NWSYM
        DO J=1,NSYM(I)
          IF(MODE.EQ.ISYM(I,J))NMSYM=I
        END DO
      END DO
      IF(NMSYM.EQ.1)THEN
        MD=1
      ELSE
        MD=2
      END IF
      DO ILOOP=1,MD
        IDEGEN=ILOOP
        CALL REDEFH(W(LXK),W(LYK),KXK,MD)
        IF(MODE.NE.NVMODE)THEN
          CALL DIAG(W(LYK),W(LYK),KXK,NDEGEN,-1,W(LWK),W(LTEMP),W(LWK),
     1    KXK,KXK,W(LISTAT),NSTAT,NSMODE,W(LYK),W(LYK),W(LYK),IDUM,
     2    IDUM,IDUM)
        ELSE
          CALL DIAG(W(LYK),W(LYK),KXK,NDEGEN,0,W(LWK),W(LTEMP),W(LWK),
     1    KXK,KXK,W(LISTAT),NSTAT,NSMODE,W(LYK),W(LYK),W(LYK),IDUM,
     2    IDUM,IDUM)
C**FIND ENERGY FOR THIS STATE....CORRESPONDS TO QUANTUM OF FINAL MODE
          CALL REDEFE(W(LTEMP),W(LEVAL),KXK,MD)
          IF(ILOOP.EQ.MD)
     1    CALL ENERGY(W(LISTAT),NSTAT,NSMODE,ISTATE,MODE,W(LEVAL),E,
     2    W(LESCF+ISTATE-1),W(LWK))
        END IF
        IF(ITIM.EQ.1.AND.ITIM1B.EQ.2)CALL TIMIT(3)
        CALL FLUSH(IOUT)
C**FORM NEW FUNCTIONS THIS MODE IF REQUIRED
        IF(.NOT.LINEAR)THEN
          IF(KROT.EQ.1.OR.(MODE.EQ.NSMODE.AND.JREACT.GT.0))
     1    CALL REFORM(W(LH+MNLK2+MLNK2),IMODE3,W(LYK),KXK,IMODE3,
     2    IMODE2,W(LWK),W(LXK),W(LOV),MODE,1+KA,W(LEVAL),W(LMODNT))
        ELSE
          IF(KROT.EQ.1.OR.MOD(KA,2).NE.0.OR.
     1    (MODE.EQ.NSMODE.AND.JREACT.GT.0))
     1    CALL REFORM(W(LH+MNLK2),IMODE3,W(LYK),KXK,IMODE3,IMODE2,
     2    W(LWK),W(LXK),W(LOV),MODE,1+KKA,W(LEVAL),W(LMODNT))
        END IF
      END DO
      CALL MEMO(-2,LEVAL,KEVAL,LTEMP,KEVAL,0,0,0,0,0,0)
      MNLK2=MNLK2+IMODE1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
      IF(MODE.GT.NONLIN)MNLK2=MNLK2+IMODE1
C     MNLK3=MNLK3+IMODE2
      CALL MEMO(-4,LXK,KXK*KXK,LWK,KWK,LYK,KXK*KXK,LOV,KXK*KXK,0,0)

C**ANALYTIC
      IF((IWHICH.LT.0.AND.MOLINC.GT.0).OR.IFITRP.NE.0)THEN
        K2=0
        K3=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          IF(K.EQ.MODE)
     1    CALL GETANI(W(LXKAN),MAXQU,MAXPOW,K,W(LH+K2),W(LXQ+K3),IK3,
     2    IK2,W(LXTANH))
          K2=K2+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(K.GT.NONLIN)K2=K2+IK1
          K3=K3+IK2
        END DO
C**ANALYTIC
      END IF

1000  CONTINUE
C**********************************************************************
C**                                           TEST FOR SELF-CONSISTENCY
C**********************************************************************
      IF(ISTATE.EQ.1.AND.J21.EQ.1)THEN
        CALL PRSCF(W(LISTAT),NSTAT,NMODE,ISTATE,E0,EVL,JTHIS,KA,KC,0)
        E0E=E0-EVL
      ELSE
        CALL PRSCF(W(LISTAT),NSTAT,NMODE,ISTATE,E0,E,JTHIS,KA,KC,0)
        E0E=E0-E
      END IF
c  change convergence criterion to converge abs diff in cm-1
c  to conv, eg. conv = 5e-2 cm-2.  5/12/97
      if(count.gt.15.0)then
        write(IOUT,*)ISTATE,'   did not converge'
        go to 999
      end if
      count=count + 1.0
      IF(dabs(e0e).gt.conv)then
        IF(ISTATE.EQ.1.AND.J21.EQ.1)E0=EVL
        IF(ISTATE.NE.1.OR.J21.NE.1)E0=E
        CALL FLUSH(IOUT)
C**MAKE COPY (59) -> (58)
C       IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0)THEN
C         REWIND 58
C         REWIND 59
C         DO MODE=1,NSMODE
C           CALL INTARR(W(LNBF),W(LMBF),MODE,IMODE1,IMODE2,IMODE3)
C           KXXK=IMODE3
C           IF(IREACT.GT.0.AND.MODE.EQ.NSMODE)KXXK=KXK
C           CALL MEMO(1,LWK,KXXK,0,0,0,0,0,0,0,0)
C           CALL RECYCL(W(LWK),KXXK,59,58)
C           CALL MEMO(-1,LWK,KXXK,0,0,0,0,0,0,0,0)
C         END DO
C       END IF
        ICYCLE=ICYCLE+1
        GO TO 555
      END IF
C**********************************************************************
C**                                                 SCF STATE CONVERGED
C**********************************************************************
999   continue
      IF(ISTATE.EQ.1.AND.J21.EQ.1)THEN
        CALL PRSCF(W(LISTAT),NSTAT,NSMODE,ISTATE,E0,EVL,JTHIS,KA,KC,40)
      ELSE
        CALL PRSCF(W(LISTAT),NSTAT,NSMODE,ISTATE,E0,E,JTHIS,KA,KC,40)
      END IF
C**VIBRATIONAL GROUND STATE K=0 SCF FUNCTION ONLY IF VIRTUALS
C     IF(ICI.LT.0.AND.J21.NE.1)GO TO 1600
      IF(J21.NE.1.AND.KROT.NE.J21)THEN
        WRITE(IOUT,320)
        WRITE(IOUT,320)
      END IF

C**CONTRACT BASIS IF NVF<0
      K2=0
      KX2=0
      IMOVE=0
      DO K=1,NVMODE
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        NOCON=INTFUN(W(LTEMP5),K)
C**NBF MODIFIED
        IF(NOCON.LT.0)THEN
          CALL NMOVE(W(LNBF),K,-NOCON)
        END IF
        CALL INTARR(W(LNBF),W(LMBF),K,IX1,IX2,IX3)
        IF(NOCON.LT.0.OR.IMOVE.NE.0)THEN
          CALL HMOVE(W(LH+K2),W(LH+KX2),I3,IX3,I2)
          IMOVE=1
        END IF
        K2=K2+I1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
        IF(K.GT.NONLIN)K2=K2+I1
        KX2=KX2+IX1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
        IF(K.GT.NONLIN)KX2=KX2+IX1
      END DO
7000  CONTINUE
C**********************************************************************
C**                                          TIDY UP FOR THIS SCF STATE
C**********************************************************************
      IF(J21.EQ.1)THEN
        IF(JMAX.NE.0)THEN
          J21=2*JMAX+1
          JTHIS=JMAX
C**RPH ODD K
          WRITE(IOUT,320)
          WRITE(IOUT,315)JMAX
          WRITE(IOUT,320)
          GO TO 8000
C**RPH ODD K
        END IF
      END IF
      IF(ISCFCI.NE.0)THEN
        K2=0
        DO K=1,NSMODE
          CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
C**USE SCF FUNCTIONS
          IF(ICI.GT.0)
     1    CALL OUT29(W(LH+K2),I3,I2,W(LISTAT),NSTAT,NSMODE,
     2    ISTATE,K)
C**USE VIRTUAL FUNCTIONS FOR VIBRATIONAL GROUND STATE IN CORE
          K2=K2+I1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(K.GT.NONLIN)K2=K2+I1
        END DO
      END IF
C**FINISHED THIS STATE
C     CALL TIMIT(3)
      ITIM=0
C**VIBRATIONAL GROUND STATE SCF FUNCTION ONLY IF VIRTUALS
      IF(ICI.LT.0)GO TO 1600
C**ONLY NEED FIRST ICI SCF STATES IF CI
      IF(ISTATE.EQ.ICI)GO TO 1600
1500  CONTINUE
      WRITE(IOUT,215)
      RETURN
C**********************************************************************
C**                                                  RETURN IF SCF ONLY
C**********************************************************************
1600  CONTINUE
      RETURN
      END
C**************************************************************
C**TEMPORARY
      SUBROUTINE PRH(H,N,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(N,M,3)
      COMMON/FILASS/IOUT
      DO K=1,3
        WRITE(IOUT,*)'K = ',K
        DO J=1,M
          WRITE(IOUT,*)'J = ',J
          WRITE(IOUT,*)(H(I,J,K),I=1,N)
        END DO
      END DO
      RETURN
      END
C**TEMPORARY
C**************************************************************
      SUBROUTINE VCI(W,KCONT,LCONT,NNMAX,KLP,KSCFCI,NVSYMX,EVLJ0,K2TAU,
     1K3TAU,MAXBAS,KMXBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      LOGICAL LGIV,LINEAR,LANCZ,LANZA,LANZB,TRIAT,ABINIT,LOGTMP
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1),MAXBAS(KMXBAS,6)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
     1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
     4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
     5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
     6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
     7LDP4A,LDP4B,LDP4C,LDP4D,LXP,LJP6,LANCNT,LANK,LXJ,LVP5,
     8LC5,LC6,LVP0,LTMPRD,LMVB,LNFC,LTEMP6,LTEMP7,LTEMP8,LTEMP9
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/HERM/IHERM
      COMMON/CHECK/MCHECK
      COMMON/SADDLE/JNORM
      COMMON/PATH/ISCFCI
      COMMON/CYCLE/ICYCLE
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/REACTL/JREACT,IFITRP
      COMMON/ANALYT/MAXQU,MAXPOW
      COMMON/DEGEN/IDEGEN,NDEGEN
      COMMON/ENTER/IENTER,IENTMX(5),NTOT1,NTOT2,NTOT3,NTOT4,NTOT5
      COMMON/SINCOS/ICS
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FITTER/MFIT,MFIT1(4),XFIT1,MFIT2(4),XFIT2,MFIT3(4),XFIT3,
     1MFIT4(4),XFIT4,MFIT5(4),XFIT5
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TRIATO/TRIAT
      COMMON/LANTOL/TOLLAN
      COMMON/CYCLES/NCYCLE
C**TEMPORARY (DIMENSIONS)
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LAN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TLLCUT,EVLCUT,TOLMAT,LANCCC,LANINC
      COMMON/GIVEN/LGIV,IGIV
      COMMON/MATSIZ/MATSIZ
      COMMON/TYPE/LINEAR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP,MOUTIN,INP4,INP5,INP6,INP7
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/TOLS/TOL,EPS
      COMMON/ECKCNT/ICNT,INTC
      COMMON/ABINIT/ABINIT
      COMMON/EVL/EVL,CUT
      COMMON/COUPLE/ICOUPL,JCOUPL,ICOUPC,JCOUPC
      COMMON/WHICH/IWHICH
      COMMON/NCPOT/NPOT
      COMMON/FACTOR/FACTOR(6),FACTS(6)
      COMMON/VCIMAX/NMAX
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/RPHROT/IROTV
      COMMON/ESTATE/IORDER
      COMMON/JKAKC/JTHIS,KA,KC,KKROT
      COMMON/AXES/MX(3),MXDIP(5),MXROT(3),ISYMT
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/TBASIS/NTBAS(6,2),NTAU(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/TSYMM/INTSM(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CSAVES/NNMODS,NAMODS,NVMODS,ICOUPS,ICOUCS,JREACS,NREACS
      COMMON/XSAVES/ICOUPX,ICOUCX,JCOUPS,JCOUCS
      COMMON/NCREC/NREC(6),MAXBUF(6),NUNITR,NUNITW
      COMMON/CSIZES/ISIZM1,ISIZM2,NVAL1,NVAL2,ICSIZ1,ICSIZ2,
     1IPSIZ1,IPSIZ2
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/MVAL/MVAL1,MVAL2
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/UNITEX/I75,I76,I85,I86
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
      COMMON/KLSZS/KLJP1,KLJP2,KLJP3,KLJP4,KLJP5,KLJP6
      COMMON/KJSZS/KJP1,KJP2,KJP3,KJP4,KJP5,KJP6
      COMMON/KPSZS/KIP1,KIP2,KIP3,KIP4,KIP5,KIP6,JCC1,JCC2,JCC3,JCC4,
     1JCC5,JCC6
      COMMON/IPSZS/KPPP1,KPPP2,KPPP3,KPPP4,KPPP5,KPPP6
      COMMON/SIZEJ/JSIZE1(2),JSIZE2(2),JSIZE3(2),JSIZE4(2),JSIZE5(2),
     1JSIZE6(2)
      COMMON/TOTALS/ITOT1(2),ITOT2(2),ITOT3(2),ITOT4(2),ITOT5(2),
     1ITOT6(2),ITOT
      COMMON/TOTK/KTOT(6,2)
      COMMON/MATRIX/NVALV,NVALR,KSTEP,KSIGN,NVALCF
C**************************************************************
240   FORMAT(/,1X,'ZERO POINT ENERGY = ',F10.2,/)
260   FORMAT(//,1X,'FINAL VIBRATIONAL (K-DIAGONAL) CI ENERGIES',/)
285   FORMAT(//,1X,'SIZE OF VIBRATIONAL CI MATRIX IS ',I6,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
290   FORMAT(//,1X,'SEQUENCE ',I2,' VIBRATIONAL SYMMETRY ',I2)
395   FORMAT(//,1X,'SIZES OF CI SYMMETRY BLOCKS: ',/,1X,10I10,/)
400   FORMAT(//,1X,'SHOULD NOT OCCUR',/)
502   FORMAT(//,1X,'ADJUSTED VALUE NVAL FOR SCHEME: ',I1,' = ',I4,/)
C**************************************************************
C**TEMPORARY
      IF(LCOUNT.EQ.2)THEN
        LOGTMP=LANCZ
C       LANCZ=.FALSE.
      END IF
C**TEMPORARY
      IF((MATSIZ.NE.0.AND.ICI.LT.0).OR.
     1(LANCZ.AND.LLAN20.LT.0))GO TO 3331
C**TEMPORARY
      CALL MEMPRT
C**TEMPORARY

C**RESET IMPORTANT QUANTITIES FOR EACH CONTRACTION IF REQUIRED
      IF(LCOUNT.GT.1)THEN
        NUMBER=ICONT(LCONT)
        JREACT=0
        NREACT=0
        DO I=1,NUMBER
          J=JCONT(LCONT,I)
          IF(J.EQ.IREACT)JREACT=IREACT
        END DO
        IF(JREACT.NE.0)NREACT=1
        NAMODE=NUMBER-NREACT
        NVMODE=NUMBER
C**NNMODE LIKE NSMODE.......TOTAL MODES IN VCI
        NNMODE=NUMBER

CC      ICOUPL=MIN0(ICOUPS,NAMODE)
CC      ICOUPC=MIN0(ICOUCS,NAMODE)
        ICOUPL=MIN0(ICOUPL,NAMODE)
        ICOUPC=MIN0(ICOUPC,NAMODE)
        NCOUPL(LCONT)=MIN0(NCOUPL(LCONT),NAMODE)
        NCOUPC(LCONT)=MIN0(NCOUPC(LCONT),NAMODE)

        DO ICPL=1,ICOUPL
          FACTOR(ICPL)=1.D0
          DO I=2,NAMODE
            FACTOR(ICPL)=FACTOR(ICPL)*I
          END DO
          DO I=1,ICPL
            FACTOR(ICPL)=FACTOR(ICPL)/I
          END DO
          IDENOM=NAMODE-ICPL
          DO I=1,IDENOM
            FACTOR(ICPL)=FACTOR(ICPL)/I
          END DO
          IF(JREACT.GT.0)FACTOR(ICPL)=FACTOR(ICPL)+1
        END DO
C**********************************************************************
C**********************************************************************
C**                      DUMP POTENTIAL AND CORIOLIS DATA TO DISC (VCI)
C**                                 FOR EACH CONTRACTION SCHEME IN TURN
C**********************************************************************
C**********************************************************************
        IF(JREACS.GT.0.AND.JREACT.LE.0)THEN
          I61=161
          I62=162
          I63=163
          I64=164
          I71=171
          I72=172
          I73=173
          I74=174
          I75=175
          I76=176
          I81=181
          I82=182
          I83=183
          I84=184
          I85=185
          I86=186
          I91=191
          I92=192
          I93=193
          I94=194
          CALL DISCDP(W,KCONT,LCONT)
        ELSE
          I61=61
          I62=62
          I63=63
          I64=64
          I71=71
          I72=72
          I73=73
          I74=74
          I75=75
          I76=76
          I81=81
          I82=82
          I83=83
          I84=84
          I85=85
          I86=86
          I91=91
          I92=92
          I93=93
          I94=94
        END IF
C**********************************************************************
C**********************************************************************
      END IF
C**SAVE INPUT VALUES
      ICSSSL=ICOUPL
      ICSSSC=ICOUPC
C**TEMPORARY
      WRITE(IOUT,*)'VCI - KCONT,LCONT,JREACT = ',KCONT,LCONT,JREACT
C     ICPXSV=ICOUPX
C     ICCXSV=ICOUCX
      IF(KCONT.EQ.2)THEN
C       ICOUPX=ICOUPL
C       ICOUCX=ICOUPC
      END IF
      WRITE(IOUT,*)'ICOUPL,ICOUPC,NCOUPL(1),NCOUPC(1),NCOUPL(2),
     1NCOUPC(2)'
      WRITE(IOUT,*)ICOUPL,ICOUPC,NCOUPL(1),NCOUPC(1),NCOUPL(2),
     1NCOUPC(2)
C**TEMPORARY
C     ICOUPX=MAX0(ICOUPL,IABS(NCOUPL(1)),IABS(NCOUPL(2)))
C     ICOUCX=MAX0(ICOUPC,IABS(NCOUPC(1)),IABS(NCOUPC(2)))
C**TEMPORARY
      WRITE(IOUT,*)'JCOUPL,JCOUPC,ICOUPX,ICOUCX,ICOUPS,ICOUCS'
      WRITE(IOUT,*)JCOUPL,JCOUPC,ICOUPX,ICOUCX,ICOUPS,ICOUCS
      WRITE(IOUT,*)'NMODE,NAMODE,NSMODE,NVMODE,NNMODE = ',
     1NMODE,NAMODE,NSMODE,NVMODE,NNMODE
      WRITE(IOUT,*)'NAMODS,NVMODS,NNMODS = ',NAMODS,NVMODS,NNMODS
      CALL FLUSH(IOUT)
C**TEMPORARY

      ITIM=-1
C**********************************************************************
C**                                  GET BASIC INTEGRALS FOR K-DIAGONAL
C**********************************************************************
      IF(ICOUPX.GE.0)REWIND 21
      IF(ICOUPX.GT.1)REWIND 22
      IF(ICOUPX.GT.2)REWIND 23
      IF(ICOUPX.GT.3)REWIND 24
      IF(ICOUPX.GT.4)REWIND 25
      IF(ICOUPX.GT.5)REWIND 26
C**TEMPORARY-LIN
      REWIND 70
C**TEMPORARY-LIN
      DO 3330 KROT=1,J21,KSTEP
      KKROT=KROT
      KA=KROT/2
      KC=KCORIG-KA
C**RPH ODD K
      IF(KA.GT.0.AND.MOD(KA,2).EQ.MOD(KROT,2))KC=KC+1
C**RPH ODD K
C**TEMPORARY-LIN
      IF(KROT.EQ.1.OR.MOD(KROT,2).EQ.0)THEN
        K2=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          IF(K.GT.NONLIN)THEN
C**NEED TO READ ACCORDING TO KSTEP
            CALL IN70(W(LH+K2),IK3,IK2)
          END IF
          K2=K2+IK1
          IF(K.GT.NONLIN)K2=K2+IK1
        END DO
      END IF
C**TEMPORARY-LIN
C**********************************************************************
C**                                                       LOOP ROUND Ka
C**********************************************************************
      ITIM=ITIM+1
CCCC  IF(JREACS.LE.0.OR.(JREACS.GT.0.AND.JREACT.GT.0))THEN
        IF(ICOUPX.GE.0)THEN
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND I71
        END IF
        IF(ICOUCX.GE.0)THEN
          IF(J21.GT.1)REWIND I61
          REWIND I81
        END IF
        IF(ICOUPX.GT.1)THEN
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND I72
        END IF
        IF(ICOUCX.GT.1)THEN
          IF(J21.GT.1)REWIND I62
          REWIND I82
        END IF
        IF(ICOUPX.GT.2)THEN
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND I73
        END IF
        IF(ICOUCX.GT.2)THEN
          IF(J21.GT.1)REWIND I63
          REWIND I83
        END IF
        IF(ICOUPX.GT.3)THEN
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND I74
        END IF
        IF(ICOUCX.GT.3)THEN
          IF(J21.GT.1)REWIND I64
          REWIND I84
        END IF
        IF(ICOUPX.GT.4)THEN
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND I75
        END IF
        IF(ICOUCX.GT.4)REWIND I85
        IF(ICOUPX.GT.5)THEN
          IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND I76
        END IF
        IF(ICOUCX.GT.5)REWIND I86
CCCC  END IF
C**BOTH LINEAR AND RPH
      IF(IREACT.GT.0)THEN
        CALL INTARR(W(LNBF),W(LMBF),NSMODE,IK1TAU,IK2TAU,IK3TAU)
C**DUMMY READ IF RPH BUT NOT THIS TIME
CCCC    IF(JREACS.GT.0.AND.JREACT.LE.0)THEN
CCCC      CALL DUMRV0(IK2TAU,W(LC0),W(LC0),W(LEJK0),W(LEJK0),J21,
CCCC 1    W(LMODNT))
CCCC    END IF
C**GET TORSION-ONLY INTEGRALS IF RPH
        IF(JREACT.GT.0)THEN
C*****************************   LIKE V0CI1 (NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
          CALL MEMO(1,LIP1,KIP1,0,0,0,0,0,0,0,0)
          CALL GETBP1(W(LIP1),ISIZE1,NMODE,NSMODE,W(LMXBAS),NSIZE1,0)
          KXA1=NSIZE1*(NSIZE1+1)/2
          CALL MEMO(1,LXA0,KXA1,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
          CALL DIAGZ(NSIZE1,NSIZE1,W(LXA0),NSIZE1,W(LXA0),ICI)
          CALL V0CV0(NAMODE,1,W(LH+K2TAU),W(LXQ+K3TAU),W(LXA0),NSIZE1,
     1    IK3TAU,IK2TAU,W(LIP1),ISIZE1,W(LC0),W(LC0),W(LEJK0),W(LEJK0),
     2    J21,KROT,W(LMODNT))
          CALL MATOUT(W(LXA0),W(LXA0),NSIZE1,21)
          CALL MEMO(-1,LXA0,KXA1,0,0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0CI1 (NMODE)
        END IF
      END IF
C**NUMBER COULD BE SINGLE (RPH) MODE.......DO NOTHING MOE
      IF(NAMODE.EQ.0)GO TO 3332
C**INTRINSIC
      IF(KCONT.LT.0)THEN
        JCOUPS=JCOUPL
        JCOUCS=JCOUPC
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K=0
      K2=0
      K3=0
      KRK=1
      DO KK=1,NAMODS
        IF(KK.EQ.1.AND.ITIM.EQ.0)THEN
          ITIM1A=0
          ITIM1B=0
        END IF
        IF(LCOUNT.GT.1)THEN
          IF(KRK.LE.NAMODE)THEN
            KNEXT=JCONT(LCONT,KRK)
            IF(KK.EQ.KNEXT)KRK=KRK+1
          ELSE
            KNEXT=0
          END IF
        ELSE
          KNEXT=KK
        END IF
        CALL INTARR(W(LNBF),W(LMBF),KK,IK1,IK2,IK3)
C**POINT TO L(0) OR L(1) IF LINEAR
        MK2=0
C       IF(KK.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
C         CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,MK2,IK1)
C         MK2=IK1
C       END IF
C**FINISHED WITH KK IF KNEXT=0
CCCC    IF(KNEXT.EQ.0)GO TO 330
C**NEXT K
        K=KNEXT
        K1=K-1

C**DO WE WANT THIS KK?
C       IF(KK.NE.KNEXT)THEN
C**DUMMY READ IF MISSING MODE KK
C         IF((JREACT.LE.0.AND.KK.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
C    1    THEN
C           IF(ICOUPX.GT.0.AND.IREACT.EQ.0)
C    1      CALL DUMRD1(KK,IK2,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
C    2      W(LEJK1),J21,W(LMODNT))
C         ELSE
C           IF(ICOUPX.GT.0)CALL DUMRV1(KK,IK2,IK2TAU,W(LV1),W(LV1),
C    1      W(LC1),W(LC1),W(LEJK1),W(LEJK1),J21,W(LMODNT))
C         END IF
CC        GO TO 330
C         GO TO 331
C       END IF


C******************************************************************
        KPPP1=KIP1
C**POSSIBLY LINEAR (JREACT<0)
        IF(JREACT.LE.0.AND.ICOUPX.EQ.0)THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
          CALL MEMO(1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
          CALL GETBP1(W(LIP1),ISIZE1,NMODE,K,W(LMXBAS),NSIZE1,0)
          KXA1=NSIZE1*(NSIZE1+1)/2
C**ZEROISE MATRIX
          CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL DIAGZ(NSIZE1,NSIZE1,W(LXA1),NSIZE1,W(LXA1),ICI)
C**'KINETIC ENERGY'
          CALL V0CIKE(NAMODE,1,K,W(LH+K2+MK2),W(LXQ+K3),W(LXA1),NSIZE1,
     1    IK3,IK2,W(LOMEGA+K1),W(LMODNT),W(LIP1),ISIZE1)
          CALL MATOUT(W(LXA1),W(LXA1),NSIZE1,21)
          CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
          GO TO 3301
        END IF
C330   CONTINUE
C******************************************************************


C**INTRINSIC
C       IF(KCONT.LT.0)THEN
          CALL GROUP1(KK,ICONT(1),JCONT,1,IND)
          IF(IND.NE.0)THEN
            JCOUPL=NCOUPL(1)
            JCOUPC=NCOUPC(1)
            ICOUPL=IABS(JCOUPL)
            ICOUPC=IABS(JCOUPC)
            GO TO 3031
          END IF
          CALL GROUP1(KK,ICONT(2),JCONT,2,IND)
          IF(IND.NE.0)THEN
            JCOUPL=NCOUPL(2)
            JCOUPC=NCOUPC(2)
            ICOUPL=IABS(JCOUPL)
            ICOUPC=IABS(JCOUPC)
            GO TO 3031
          END IF
          JCOUPL=JCOUPS
          JCOUPC=JCOUCS
          ICOUPL=ICOUPS
          ICOUPC=ICOUCS
3031      CONTINUE
C       END IF

        IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.1)GO TO 331


C**DO WE WANT THIS KK?
        IF(KK.NE.KNEXT.OR.(KK.EQ.KNEXT.AND.ICOUPL.LT.1))THEN
C**DUMMY READ IF MISSING MODE KK
          IF((JREACT.LE.0.AND.KK.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     1    THEN
            IF(ICOUPX.GT.0.AND.IREACT.EQ.0)
     1      CALL DUMRD1(KK,IK2,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
     2      W(LEJK1),J21,W(LMODNT))
          ELSE
            IF(ICOUPX.GT.0)CALL DUMRV1(KK,IK2,IK2TAU,W(LV1),W(LV1),
     1      W(LC1),W(LC1),W(LEJK1),W(LEJK1),J21,W(LMODNT))
          END IF
          GO TO 331
        END IF


C**SKIP IF UNWANTED KK
        IF(KNEXT.EQ.0.OR.KK.NE.KNEXT)GO TO 331
C******************************************************************
        KPPP1=KIP1
C**CORIOLIS AND POTENTIAL
        IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
          CALL MEMO(1,LIP1,KPPP1,0,0,0,0,0,0,0,0)
          CALL GETBP1(W(LIP1),ISIZE1,NMODE,K,W(LMXBAS),NSIZE1,0)
          KXA1=NSIZE1*(NSIZE1+1)/2
          CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
          CALL DIAGZ(NSIZE1,NSIZE1,W(LXA1),NSIZE1,W(LXA1),ICI)
C         IF(ICOUPC.GE.1)
C    1    CALL V0CI1(NAMODE,1,K,W(LH+K2+MK2),W(LXQ+K3),W(LXA1),NSIZE1,
C    2    IK3,IK2,W(LIP1),ISIZE1,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
C    3    W(LEJK1),J21,KROT,W(LMODNT),W(LOMEGA+K1),1)
          CALL V0CI1(NAMODE,1,K,W(LH+K2+MK2),W(LXQ+K3),W(LXA1),NSIZE1,
     1    IK3,IK2,W(LIP1),ISIZE1,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
     4    W(LEJK1),J21,KROT,W(LMODNT),W(LOMEGA+K1),0)
          CALL MATOUT(W(LXA1),W(LXA1),NSIZE1,21)
          CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
        ELSE
C*****************************   LIKE V0CI2 (K+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
          IF(K.GT.NONLIN)THEN
            KPPP2=KLJP2
            IND2=1
          ELSE
            KPPP2=KIP2
            IND2=0
          END IF
          CALL MEMO(1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
          IF(K.GT.NONLIN)THEN
            CALL GETBP1(W(LIP2),LSIZE2,NMODE,K,W(LMXBAS),NSIZE2,IND2)
            MSIZE2=NSIZE2
          ELSE
            CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,NSMODE,W(LMXBAS),
     1      NSIZE2,IND2)
            MSIZE2=ISIZE2
          END IF
          KXA2=NSIZE2*(NSIZE2+1)/2
          JCI=MAXBFN(W(LMXBAS),NMODE,K,1)
          IF(JREACT.GT.0)THEN
            JCIM=MAX0(NTAU(1)+1,NTAU(2)+1)
          ELSE
            JCIM1=2*(NBAS(1,1)+1)-1
            JCIM2=2*(NBAS(1,2)+1)-1
            JCIM=MAX0(JCIM1,JCIM2)
          END IF
C**BEWARE LINEAR
          I9=1
          IF(ICOUPC.GT.0)I9=9
C**BEWARE LINEAR
          KTEMP=JCI*JCI*I9
          KXK0=I9*JCI*JCI*IK2
          KXL0=I9*JCIM*JCIM*IK2TAU
          CALL MEMO(2,LXA1,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
          CALL MEMO(2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
C**ZEROISE MATRIX
          CALL DIAGZ(NSIZE2,NSIZE2,W(LXA1),NSIZE2,W(LXA1),ICI)
          CALL V0CV1(NAMODE,1,2,K,W(LH+K2+MK2),W(LXQ+K3),W(LH+K2TAU),
     1    W(LXQ+K3TAU),W(LXA1),NSIZE2,IK3,IK2,IK3TAU,IK2TAU,W(LIP2),
     2    MSIZE2,W(LTEMP),JCI,JCIM,W(LXK0),W(LXL0),W(LV1),W(LV1),
     3    W(LC1),W(LC1),W(LEJK1),W(LEJK1),J21,KROT,W(LMODNT),
     4    W(LOMEGA+K1),I9)
          CALL MATOUT(W(LXA1),W(LXA1),NSIZE2,21)
          CALL MEMO(-2,LXA1,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
          CALL MEMO(-2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0CI2 (K+NMODE)
        END IF
331   CONTINUE
C**INTRINSIC
        IF(ICOUPX.EQ.1)GO TO 3301
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L=0
        L2=0
        L3=0
        LRL=1
        DO LL=1,KK-1
          IF(LL.EQ.1.AND.KK.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          IF(LCOUNT.GT.1)THEN
            IF(LRL.LE.NAMODE)THEN
              LNEXT=JCONT(LCONT,LRL)
              IF(LL.EQ.LNEXT)LRL=LRL+1
            ELSE
              LNEXT=0
            END IF
          ELSE
            LNEXT=LL
          END IF
          CALL INTARR(W(LNBF),W(LMBF),LL,IL1,IL2,IL3)
C**POINT TO L(0) OR L(1) IF LINEAR
          ML2=0
C         IF(LL.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
C           CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,ML2,IL1)
C           ML2=IL1
C         END IF
C**FINISHED WITH LL IF LNEXT=0
CCCC      IF(LNEXT.EQ.0)GO TO 332
C**NEXT L
          L=LNEXT

C**DO WE WANT THIS KK AND THIS LL?
C         IF(KK.NE.KNEXT.OR.LL.NE.LNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL
C           IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
C    1      (NSMODE.EQ.NAMODE))THEN
C             IF(ICOUPX.GT.1.AND.IREACT.EQ.0)
C    1        CALL DUMRD2(KK,LL,IK2,IL2,W(LV2),W(LV2),
C    2        W(LC2),W(LC2),W(LEJK2),W(LEJK2),J21,W(LMODNT))
C           ELSE
C             IF(ICOUPX.GT.1)
C    1        CALL DUMRV2(KK,LL,IK2,IL2,IK2TAU,W(LV2),W(LV2),W(LC2),
C    2        W(LC2),W(LEJK2),W(LEJK2),J21,W(LMODNT))
C           END IF
C           GO TO 332
C         END IF


C**INTRINSIC
C         IF(KCONT.LT.0)THEN
            CALL GROUP2(KK,LL,ICONT(1),JCONT,1,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(1)
              JCOUPC=NCOUPC(1)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 3032
            END IF
            CALL GROUP2(KK,LL,ICONT(2),JCONT,2,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(2)
              JCOUPC=NCOUPC(2)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 3032
            END IF
            JCOUPL=JCOUPS
            JCOUPC=JCOUCS
            ICOUPL=ICOUPS
            ICOUPC=ICOUCS
3032        CONTINUE
C         END IF

          IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.2)GO TO 332


C**DO WE WANT THIS KK AND THIS LL?
          IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.(KK.EQ.KNEXT.AND.
     1    LL.EQ.LNEXT.AND.ICOUPL.LT.2))THEN
C**DUMMY READ IF MISSING MODES KK AND LL
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1      (NSMODE.EQ.NAMODE))THEN
              IF(ICOUPX.GT.1.AND.IREACT.EQ.0)
     1        CALL DUMRD2(KK,LL,IK2,IL2,W(LV2),W(LV2),
     1        W(LC2),W(LC2),W(LEJK2),W(LEJK2),J21,W(LMODNT))
            ELSE
              IF(ICOUPX.GT.1)
     1        CALL DUMRV2(KK,LL,IK2,IL2,IK2TAU,W(LV2),W(LV2),W(LC2),
     2        W(LC2),W(LEJK2),W(LEJK2),J21,W(LMODNT))
            END IF
            GO TO 332
          END IF


C**SKIP IF UNWANTED LL
          IF(LNEXT.EQ.0.OR.LL.NE.LNEXT)GO TO 332
C******************************************************************
C**CORIOLIS AND POTENTIAL
          IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1    (NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
            KPPP2=KIP2
            CALL MEMO(1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
            CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,L,W(LMXBAS),NSIZE2,0)
            KXA2=NSIZE2*(NSIZE2+1)/2
            CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
            JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
            JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
            I5=1
            IF(ICOUPC.GT.1)I5=5
            KTEMP=JCI2*JCI2*I5
            KXK0=I5*JCI1*JCI1*IK2
            KXL0=I5*JCI2*JCI2*IL2
            CALL MEMO(1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
            CALL MEMO(2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE2,NSIZE2,W(LXA2),NSIZE2,W(LXA2),ICI)
C           IF(ICOUPC.GE.2)
C    1      CALL V0CI2(NAMODE,1,2,K,L,W(LH+K2+MK2),W(LXQ+K3),
C    2      W(LH+L2+ML2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXA2),NSIZE2,
C    3      W(LIP2),ISIZE2,W(LTEMP),JCI1,JCI2,W(LXK0),W(LXL0),W(LV2),
C    4      W(LV2),W(LC2),W(LC2),W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),
C    5      I5,1)
            CALL V0CI2(NAMODE,1,2,K,L,W(LH+K2+MK2),W(LXQ+K3),
     1      W(LH+L2+ML2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXA2),NSIZE2,
     2      W(LIP2),ISIZE2,W(LTEMP),JCI1,JCI2,W(LXK0),W(LXL0),W(LV2),
     3      W(LV2),W(LC2),W(LC2),W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),
     4      I5,0)
            CALL MATOUT(W(LXA2),W(LXA2),NSIZE2,22)
            CALL MEMO(-1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
            CALL MEMO(-2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
            CALL MEMO(-1,LXA2,KXA2,0,0,0,0,0,0,0,0)
            CALL MEMO(-1,LIP2,KIP2,0,0,0,0,0,0,0,0)
          ELSE
C*****************************   LIKE V0CI3 (K+L+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
            IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
              KPPP3=KLJP3
              IND3=1
            ELSE
              KPPP3=KIP3
              IND3=0
            END IF
            CALL MEMO(1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
            IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
              CALL GETBP2(W(LIP3),LSIZE3,NMODE,K,L,W(LMXBAS),NSIZE3,
     1        IND3)
              MSIZE3=NSIZE3
            ELSE
              CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,NSMODE,W(LMXBAS),
     1        NSIZE3,IND3)
              MSIZE3=ISIZE3
            END IF
            KXA3=NSIZE3*(NSIZE3+1)/2
            JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
            JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
            IF(JREACT.GT.0)THEN
              JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1)
            ELSE
              JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1)-1
              JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1)-1
              JCIM=MAX0(JCIM1,JCIM2)
            END IF
C**BEWARE LINEAR
            I16=1
            IF(ICOUPC.GT.1)I16=16
C**BEWARE LINEAR
            KTEMP=JCI1*JCI2*JCI1*JCI2*I16
            KTEMP1=JCI2*JCI2*I16
            KXK0=I16*JCI1*JCI1*IK2
            KXL0=I16*JCI2*JCI2*IL2
            KXN0=I16*JCIM*JCIM*IK2TAU
            CALL MEMO(2,LXA2,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
            CALL MEMO(1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
            CALL MEMO(3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE3,NSIZE3,W(LXA2),NSIZE3,W(LXA2),ICI)
            CALL V0CV2(NAMODE,1,2,3,K,L,W(LH+K2+MK2),W(LXQ+K3),
     1      W(LH+L2+ML2),
     1      W(LXQ+L3),W(LH+K2TAU),W(LXQ+K3TAU),W(LXA2),NSIZE3,IK3,IK2,
     2      IL3,IL2,IK3TAU,IK2TAU,W(LIP3),MSIZE3,W(LTEMP),W(LTEMP1),
     3      JCI1,JCI2,JCIM,W(LXK0),W(LXL0),W(LXN0),W(LV2),W(LV2),
     4      W(LC2),W(LC2),W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),I16)
            CALL MATOUT(W(LXA2),W(LXA2),NSIZE3,22)
            CALL MEMO(-2,LXA2,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
            CALL MEMO(-1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
            CALL MEMO(-3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
            CALL MEMO(-1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0CI3 (K+L+NMODE)
          END IF
332   CONTINUE
C**INTRINSIC
          IF(ICOUPX.EQ.2)GO TO 3302
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N=0
          N2=0
          N3=0
          NRN=1
          DO NN=1,LL-1
            IF(NN.EQ.1.AND.LL.EQ.2.AND.KK.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            IF(LCOUNT.GT.1)THEN
              IF(NRN.LE.NAMODE)THEN
                NNEXT=JCONT(LCONT,NRN)
                IF(NN.EQ.NNEXT)NRN=NRN+1
              ELSE
                NNEXT=0
              END IF
            ELSE
              NNEXT=NN
            END IF
            CALL INTARR(W(LNBF),W(LMBF),NN,IN1,IN2,IN3)
C**POINT TO L(0) OR L(1) IF LINEAR
            MN2=0
C           IF(NN.GT.NONLIN.AND.MOD(KA,2).NE.0)THEN
C             CALL SELECT(W(LISTAT),NSTAT,NAMODE,ISTATE,MN2,IN1)
C             MN2=IN1
C           END IF
C**FINISHED WITH NN IF NNEXT=0
CCCC        IF(NNEXT.EQ.0)GO TO 333
C**NEXT N
            N=NNEXT


C**DO WE WANT THIS KK AND THIS LL AND THIS NN?
C           IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN
C             IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
C    1        N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C               IF(IREACT.EQ.0)CALL DUMRD3(KK,LL,NN,IK2,IL2,IN2,W(LV3),
C    1          W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,W(LMODNT))
C             ELSE
CC              CALL DUMRV3(KK,LL,NN,IK2,IL2,IN2,IK2TAU,W(LV3),W(LV3),
CC   1          W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,W(LMODNT))
C             END IF
C             GO TO 333
C           END IF


C**INTRINSIC
C           IF(KCONT.LT.0)THEN
              CALL GROUP3(KK,LL,NN,ICONT(1),JCONT,1,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(1)
                JCOUPC=NCOUPC(1)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 3033
              END IF
              CALL GROUP3(KK,LL,NN,ICONT(2),JCONT,2,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(2)
                JCOUPC=NCOUPC(2)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 3033
              END IF
              JCOUPL=JCOUPS
              JCOUPC=JCOUCS
              ICOUPL=ICOUPS
              ICOUPC=ICOUCS
3033          CONTINUE
C           END IF

            IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.3)GO TO 333


C**DO WE WANT THIS KK AND THIS LL AND THIS NN?
            IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
     1      (KK.EQ.KNEXT.AND.LL.EQ.LNEXT.AND.NN.EQ.NNEXT.AND.
     2      ICOUPL.LT.3))THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
                IF(ICOUPX.GT.2.AND.IREACT.EQ.0)
     1          CALL DUMRD3(KK,LL,NN,IK2,IL2,IN2,W(LV3),
     2          W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,W(LMODNT))
              ELSE
C               CALL DUMRV3(KK,LL,NN,IK2,IL2,IN2,IK2TAU,W(LV3),W(LV3),
C    1          W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,W(LMODNT))
              END IF
              GO TO 333
            END IF


C**SKIP IF UNWANTED NN
            IF(NNEXT.EQ.0.OR.NN.NE.NNEXT)GO TO 333
C******************************************************************
C**CORIOLIS AND POTENTIAL
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1      N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
              KPPP3=KIP3
              CALL MEMO(1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
              CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,N,W(LMXBAS),NSIZE3,
     1        0)
              KXA3=NSIZE3*(NSIZE3+1)/2
              CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
              JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
              JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
              JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
              I10=1
              IF(ICOUPC.GT.2)I10=10
              KTEMP=JCI2*JCI3*JCI2*JCI3*I10
              KTEMP1=JCI3*JCI3*I10
              KXK0=I10*JCI1*JCI1*IK2
              KXL0=I10*JCI2*JCI2*IL2
              KXN0=I10*JCI3*JCI3*IN2
              CALL MEMO(1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
              CALL MEMO(1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
              CALL MEMO(3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(NSIZE3,NSIZE3,W(LXA3),NSIZE3,W(LXA3),ICI)
C             IF(ICOUPC.GE.3)
C    1        CALL V0CI3(NAMODE,1,2,3,K,L,N,W(LH+K2+MK2),W(LXQ+K3),
C    2        W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),IK3,IK2,
C    3        IL3,IL2,IN3,IN2,W(LXA3),NSIZE3,W(LIP3),ISIZE3,W(LTEMP),
C    4        W(LTEMP1),JCI1,JCI2,JCI3,W(LXK0),W(LXL0),W(LXN0),W(LV3),
C    5        W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,KROT,
C    6        W(LMODNT),I10,1)
              CALL V0CI3(NAMODE,1,2,3,K,L,N,W(LH+K2+MK2),W(LXQ+K3),
     1        W(LH+L2+ML2),W(LXQ+L3),W(LH+N2+MN2),W(LXQ+N3),IK3,IK2,
     2        IL3,IL2,IN3,IN2,W(LXA3),NSIZE3,W(LIP3),ISIZE3,W(LTEMP),
     3        W(LTEMP1),JCI1,JCI2,JCI3,W(LXK0),W(LXL0),W(LXN0),W(LV3),
     4        W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,KROT,
     5        W(LMODNT),I10,0)
              CALL MATOUT(W(LXA3),W(LXA3),NSIZE3,23)
              CALL MEMO(-1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
              CALL MEMO(-1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
              CALL MEMO(-3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
              CALL MEMO(-1,LXA3,KXA3,0,0,0,0,0,0,0,0)
              CALL MEMO(-1,LIP3,KIP3,0,0,0,0,0,0,0,0)
            ELSE
C*****************************   LIKE V0CI4 (K+L+N+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
              IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                KPPP4=KLJP4
                IND4=1
              ELSE
                KPPP4=KIP4
                IND4=0
              END IF
              CALL MEMO(1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
              IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                CALL GETBP3(W(LIP4),LSIZE4,NMODE,K,L,N,W(LMXBAS),
     1          NSIZE4,IND4)
                MSIZE4=NSIZE4
              ELSE
                CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,NSMODE,
     1          W(LMXBAS),NSIZE4,IND4)
                MSIZE4=ISIZE4
              END IF

!              KXA4=NSIZE4*(NSIZE4+1)/2
              if (mod(nsize4,2) .eq. 0) then
                 kxa4 = nsize4 / 2
                 kxa4 = kxa4 * (nsize4 + 1)
              else
                 kxa4 = (nsize4 + 1) / 2
                 kxa4 = kxa4 * nsize4
              end if

              JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
              JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
              JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
              IF(JREACT.GT.0)THEN
                JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1,NTAU(4)+1)
              ELSE
                JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1,NBAS(3,1)+1)-1
                JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1,NBAS(3,2)+1)-1
                JCIM=MAX0(JCIM1,JCIM2)
              END IF
C**BEWARE LINEAR
              I25=1
              IF(ICOUPC.GT.2)I25=25
C**BEWARE LINEAR
              KTEMP1=JCI3*JCI3*I25
              KXK0=I25*JCI1*JCI1*IK2
              KXL0=I25*JCI2*JCI2*IL2
              KXN0=I25*JCI3*JCI3*IN2
              KXM0=I25*JCIM*JCIM*IK2TAU
              CALL MEMO(3,LXA3,KXA4,LTEMP1,KTEMP1,LTEMP5,NSIZE4,0,0,0,0)
              CALL MEMO(4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(NSIZE4,NSIZE4,W(LXA3),NSIZE4,W(LXA3),ICI)
C**GET ICOUNT2
              CALL V0CV3(NAMODE,1,2,3,4,K,L,N,W(LH+K2),W(LXQ+K3),
     1        W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+K2TAU),
     2        W(LXQ+K3TAU),W(LXA3),NSIZE4,IK3,IK2,IL3,IL2,IN3,IN2,
     3        IK3TAU,IK2TAU,W(LIP4),MSIZE4,W(LTEMP1),W(LTEMP2),
     4        W(LTEMP3),JCI1,JCI2,JCI3,JCIM,W(LXK0),W(LXL0),W(LXN0),
     5        W(LXM0),W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),
     6        J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP4),
     7        W(LTEMP5),-1,I25)
              KTEMP=ICOUNT2*2
              KTEMP2=ICOUNT2*ICOUNT2*I25
              CALL MEMO(2,LTEMP,KTEMP,LTEMP2,KTEMP2,0,0,0,0,0,0)
C**GET ICOUNT3
              CALL V0CV3(NAMODE,1,2,3,4,K,L,N,W(LH+K2),W(LXQ+K3),
     1        W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+K2TAU),
     2        W(LXQ+K3TAU),W(LXA3),NSIZE4,IK3,IK2,IL3,IL2,IN3,IN2,
     3        IK3TAU,IK2TAU,W(LIP4),MSIZE4,W(LTEMP1),W(LTEMP2),
     4        W(LTEMP3),JCI1,JCI2,JCI3,JCIM,W(LXK0),W(LXL0),W(LXN0),
     5        W(LXM0),W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),
     6        J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP4),
     7        W(LTEMP5),1,I25)
              KTEMP4=ICOUNT3*2
              KTEMP3=ICOUNT3*ICOUNT3*I25
              CALL MEMO(2,LTEMP3,KTEMP3,LTEMP4,KTEMP4,0,0,0,0,0,0)
C**GET XA3
              CALL V0CV3(NAMODE,1,2,3,4,K,L,N,W(LH+K2),W(LXQ+K3),
     1        W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+K2TAU),
     2        W(LXQ+K3TAU),W(LXA3),NSIZE4,IK3,IK2,IL3,IL2,IN3,IN2,
     3        IK3TAU,IK2TAU,W(LIP4),MSIZE4,W(LTEMP1),W(LTEMP2),
     4        W(LTEMP3),JCI1,JCI2,JCI3,JCIM,W(LXK0),W(LXL0),W(LXN0),
     5        W(LXM0),W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),
     6        J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP4),
     7        W(LTEMP5),0,I25)
              CALL MATOUT(W(LXA3),W(LXA3),NSIZE4,23)
              CALL MEMO(-2,LXA3,KXA4,LTEMP3,KTEMP3,0,0,0,0,0,0)
              CALL MEMO(-2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
              CALL MEMO(-3,LTEMP,KTEMP,LTEMP4,KTEMP4,LTEMP5,NSIZE4,
     1        0,0,0,0)
              CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,0,0)
              CALL MEMO(-1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0CI4 (K+L+N+NMODE)
            END IF
333   CONTINUE
C**INTRINSIC
            IF(ICOUPX.EQ.3)GO TO 3303
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M=0
            M2=0
            M3=0
            MRM=1
            DO MM=1,NN-1
              IF(MM.EQ.1.AND.NN.EQ.2.AND.LL.EQ.3.AND.KK.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              IF(LCOUNT.GT.1)THEN
                IF(MRM.LE.NAMODE)THEN
                  MNEXT=JCONT(LCONT,MRM)
                  IF(MM.EQ.MNEXT)MRM=MRM+1
                ELSE
                  MNEXT=0
                END IF
              ELSE
                MNEXT=MM
              END IF
              CALL INTARR(W(LNBF),W(LMBF),MM,IM1,IM2,IM3)
C**FINISHED WITH MM IF MNEXT=0
CCCC          IF(MNEXT.EQ.0)GO TO 334
C**NEXT M
              M=MNEXT

C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM?
C             IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
C    1        MM.NE.MNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM
C               IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
C    1          N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
C    2          THEN
C                 IF(IREACT.EQ.0)CALL DUMRD4(KK,LL,NN,MM,IK2,IL2,IN2,
C    1            IM2,W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
C    2            J21,W(LMODNT))
C               ELSE
CC                CALL DUMRV4(KK,LL,NN,MM,IK2,IL2,IN2,IM2,IK2TAU,
CC   1            W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,
CC   2            W(LMODNT))
C               END IF
C               GO TO 334
C             END IF

C**INTRINSIC
C             IF(KCONT.LT.0)THEN
                CALL GROUP4(KK,LL,NN,MM,ICONT(1),JCONT,1,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(1)
                  JCOUPC=NCOUPC(1)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 3034
                END IF
                CALL GROUP4(KK,LL,NN,MM,ICONT(2),JCONT,2,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(2)
                  JCOUPC=NCOUPC(2)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 3034
                END IF
                JCOUPL=JCOUPS
                JCOUPC=JCOUCS
                ICOUPL=ICOUPS
                ICOUPC=ICOUCS
3034            CONTINUE
C             END IF

              IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.4)GO TO 334

C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM?
              IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
     1        MM.NE.MNEXT.OR.(KK.EQ.KNEXT.AND.LL.EQ.LNEXT.AND.
     2        NN.EQ.NNEXT.AND.MM.EQ.MNEXT.AND.ICOUPL.LT.4))THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2          THEN
                  IF(ICOUPX.GT.3.AND.IREACT.EQ.0)CALL DUMRD4(KK,LL,NN,
     1            MM,IK2,IL2,IN2,
     2            IM2,W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
     3            J21,W(LMODNT))
                ELSE
C                 CALL DUMRV4(KK,LL,NN,MM,IK2,IL2,IN2,IM2,IK2TAU,
C    1            W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,
C    2            W(LMODNT))
                END IF
                GO TO 334
              END IF

C**SKIP IF UNWANTED MM
              IF(MNEXT.EQ.0.OR.MM.NE.MNEXT)GO TO 334
C******************************************************************
C**CORIOLIS AND POTENTIAL
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2        THEN
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
                KPPP4=KIP4
                CALL MEMO(1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
                CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,M,W(LMXBAS),
     2          NSIZE4,0)

!                KXA4=NSIZE4*(NSIZE4+1)/2
                if (mod(nsize4,2) .eq. 0) then
                   kxa4 = nsize4 / 2
                   kxa4 = kxa4 * (nsize4 + 1)
                else
                   kxa4 = (nsize4 + 1) / 2
                   kxa4 = kxa4 * nsize4
                end if

                JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
                JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
                JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
                JCI4=MAXBFN(W(LMXBAS),NMODE,M,1)
                I17=1
                IF(ICOUPC.GT.3)I17=17
                KTEMP1=JCI4*JCI4*I17
                KXK0=I17*JCI1*JCI1*IK2
                KXL0=I17*JCI2*JCI2*IL2
                KXN0=I17*JCI3*JCI3*IN2
                KXM0=I17*JCI4*JCI4*IM2

                CALL MEMO(3,LXA4,KXA4,LTEMP1,KTEMP1,LTEMP5,NSIZE4,
     1          0,0,0,0)

                CALL MEMO(4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          0,0)

C**ZEROISE MATRIX
                CALL DIAGZ(NSIZE4,NSIZE4,W(LXA4),NSIZE4,W(LXA4),ICI)
C**GET ICOUNT2
                CALL V0CI4(NAMODE,1,2,3,4,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1          W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2          W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA4),
     3          NSIZE4,W(LIP4),ISIZE4,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     4          JCI1,JCI2,JCI3,JCI4,W(LXK0),W(LXL0),W(LXN0),
     5          W(LXM0),W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
     6          J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP4),
     7          W(LTEMP5),-1,I17)
                KTEMP=ICOUNT2*2
                KTEMP2=ICOUNT2*ICOUNT2*I17
                CALL MEMO(2,LTEMP,KTEMP,LTEMP2,KTEMP2,0,0,0,0,0,0)
C**GET ICOUNT3
                CALL V0CI4(NAMODE,1,2,3,4,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1          W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2          W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA4),
     3          NSIZE4,W(LIP4),ISIZE4,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     4          JCI1,JCI2,JCI3,JCI4,W(LXK0),W(LXL0),W(LXN0),
     5          W(LXM0),W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
     6          J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP4),
     7          W(LTEMP5),1,I17)
                KTEMP4=ICOUNT3*2
                KTEMP3=ICOUNT3*ICOUNT3*I17
                CALL MEMO(2,LTEMP3,KTEMP3,LTEMP4,KTEMP4,0,0,0,0,0,0)
C**GET XA4
                CALL V0CI4(NAMODE,1,2,3,4,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1          W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2          W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA4),
     3          NSIZE4,W(LIP4),ISIZE4,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     4          JCI1,JCI2,JCI3,JCI4,W(LXK0),W(LXL0),W(LXN0),
     5          W(LXM0),W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
     6          J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP4),
     7          W(LTEMP5),0,I17)
                CALL MATOUT(W(LXA4),W(LXA4),NSIZE4,24)
                CALL MEMO(-1,LTEMP3,KTEMP3,0,0,0,0,0,0,0,0)
                CALL MEMO(-2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
                CALL MEMO(-3,LTEMP,KTEMP,LTEMP4,KTEMP4,LTEMP5,NSIZE4,
     1          0,0,0,0)
                CALL MEMO(-1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          0,0)
                CALL MEMO(-1,LIP4,KIP4,0,0,0,0,0,0,0,0)
              ELSE
C*****************************   LIKE V0CI5 (K+L+N+M+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR FIVE MODES
                IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1            M.GT.NONLIN)THEN
                  KPPP5=KLJP5
                  IND5=1
                ELSE
                  KPPP5=KIP5
                  IND5=0
                END IF
                CALL MEMO(1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
                IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1            M.GT.NONLIN)THEN
                  CALL GETBP4(W(LIP5),LSIZE5,NMODE,K,L,N,M,W(LMXBAS),
     1            NSIZE5,IND5)
                  MSIZE5=NSIZE5
                ELSE
                  CALL GETBP5(W(LIP5),ISIZE5,NMODE,K,L,N,M,NSMODE,
     1            W(LMXBAS),NSIZE5,IND5)
                  MSIZE5=ISIZE5
                END IF
                KXA5=NSIZE5*(NSIZE5+1)/2
                JCI1=MAXBFN(W(LMXBAS),NSMODE,K,1)
                JCI2=MAXBFN(W(LMXBAS),NSMODE,L,1)
                JCI3=MAXBFN(W(LMXBAS),NSMODE,N,1)
                JCI4=MAXBFN(W(LMXBAS),NSMODE,M,1)
                IF(JREACT.GT.0)THEN
                  JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1,NTAU(4)+1,
     1            NTAU(5)+1)
                ELSE
                  JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1,NBAS(3,1)+1,
     1            NBAS(4,1)+1)-1
                  JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1,NBAS(3,2)+1,
     1            NBAS(4,2)+1)-1
                  JCIM=MAX0(JCIM1,JCIM2)
                END IF
C**BEWARE LINEAR
                I36=1
                IF(ICOUPC.GT.3)I36=36
C**BEWARE LINEAR
                KTEMP1=JCI4*JCI4*I36
                KXK0=I36*JCI1*JCI1*IK2
                KXL0=I36*JCI2*JCI2*IL2
                KXN0=I36*JCI3*JCI3*IN2
                KXM0=I36*JCI4*JCI4*IM2
                KXP0=I36*JCIM*JCIM*IK2TAU
                CALL MEMO(3,LXA4,KXA5,LTEMP1,KTEMP1,LTEMP7,NSIZE5,
     1          0,0,0,0)
                CALL MEMO(5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          LXP0,KXP0)
C**ZEROISE MATRIX
                CALL DIAGZ(NSIZE5,NSIZE5,W(LXA4),NSIZE5,W(LXA4),ICI)
C**GET ICOUNT2
                CALL V0CV4(NAMODE,1,2,3,4,5,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1          W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2          W(LXQ+M3),W(LH+K2TAU),W(LXQ+K3TAU),W(LXA4),NSIZE5,IK3,
     3          IK2,IL3,IL2,IN3,IN2,IM3,IM2,IK3TAU,IK2TAU,W(LIP5),
     4          MSIZE5,W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),JCI1,
     5          JCI2,JCI3,JCI4,JCIM,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6          W(LXP0),W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
     7          J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP5),
     8          ICOUNT4,W(LTEMP6),W(LTEMP7),-2,I36)
                KTEMP=ICOUNT2*2
                KTEMP2=ICOUNT2*ICOUNT2*I36
                CALL MEMO(2,LTEMP,KTEMP,LTEMP2,KTEMP2,0,0,0,0,0,0)
C**GET ICOUNT3
                CALL V0CV4(NAMODE,1,2,3,4,5,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1          W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2          W(LXQ+M3),W(LH+K2TAU),W(LXQ+K3TAU),W(LXA4),NSIZE5,IK3,
     3          IK2,IL3,IL2,IN3,IN2,IM3,IM2,IK3TAU,IK2TAU,W(LIP5),
     4          MSIZE5,W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),JCI1,
     5          JCI2,JCI3,JCI4,JCIM,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6          W(LXP0),W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),  
     7          J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP5),
     8          ICOUNT4,W(LTEMP6),W(LTEMP7),-1,I36)
                KTEMP5=ICOUNT3*2
                KTEMP3=ICOUNT3*ICOUNT3*I36
                CALL MEMO(2,LTEMP3,KTEMP3,LTEMP5,KTEMP5,0,0,0,0,0,0)
C**GET ICOUNT4
                CALL V0CV4(NAMODE,1,2,3,4,5,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1          W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2          W(LXQ+M3),W(LH+K2TAU),W(LXQ+K3TAU),W(LXA4),NSIZE5,IK3,
     3          IK2,IL3,IL2,IN3,IN2,IM3,IM2,IK3TAU,IK2TAU,W(LIP5),
     4          MSIZE5,W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),JCI1,
     5          JCI2,JCI3,JCI4,JCIM,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6          W(LXP0),W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
     7          J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP5),
     8          ICOUNT4,W(LTEMP6),W(LTEMP7),1,I36)
                KTEMP6=ICOUNT4*2
                KTEMP4=ICOUNT4*ICOUNT4*I36
                CALL MEMO(2,LTEMP4,KTEMP4,LTEMP6,KTEMP6,0,0,0,0,0,0)
C**GET XA4
                CALL V0CV4(NAMODE,1,2,3,4,5,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1          W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2          W(LXQ+M3),W(LH+K2TAU),W(LXQ+K3TAU),W(LXA4),NSIZE5,IK3,
     3          IK2,IL3,IL2,IN3,IN2,IM3,IM2,IK3TAU,IK2TAU,W(LIP5),
     4          MSIZE5,W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),JCI1,
     5          JCI2,JCI3,JCI4,JCIM,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6          W(LXP0),W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
     7          J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),ICOUNT3,W(LTEMP5),
     8          ICOUNT4,W(LTEMP6),W(LTEMP7),0,I36)
                CALL MATOUT(W(LXA4),W(LXA4),NSIZE5,24)
                CALL MEMO(-3,LXA4,KXA5,LTEMP1,KTEMP1,LTEMP7,NSIZE5,
     1          0,0,0,0)
                CALL MEMO(-2,LTEMP,KTEMP,LTEMP2,KTEMP2,0,0,0,0,0,0)
                CALL MEMO(-2,LTEMP3,KTEMP3,LTEMP5,KTEMP5,0,0,0,0,0,0)
                CALL MEMO(-2,LTEMP4,KTEMP4,LTEMP6,KTEMP6,0,0,0,0,0,0)
                CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          LXP0,KXP0)
                CALL MEMO(-1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0CI5 (K+L+N+M+NMODE)
              END IF
334   CONTINUE
C**INTRINSIC
              IF(ICOUPX.EQ.4)GO TO 3304
C**INTEGRATE OVER FIVE NORMAL COORDINATES
              J=0
              J2=0
              J3=0
              JRJ=1
              DO JJ=1,MM-1
                IF(JJ.EQ.1.AND.MM.EQ.2.AND.NN.EQ.3.AND.LL.EQ.4.AND.
     1          KK.EQ.5.AND.ITIM.EQ.0)ITIM5A=0
                IF(LCOUNT.GT.1)THEN
                  IF(JRJ.LE.NAMODE)THEN
                    JNEXT=JCONT(LCONT,JRJ)
                    IF(JJ.EQ.JNEXT)JRJ=JRJ+1
                  ELSE
                    JNEXT=0
                  END IF
                ELSE
                  JNEXT=JJ
                END IF
                CALL INTARR(W(LNBF),W(LMBF),JJ,IJ1,IJ2,IJ3)
C**FINISHED WITH JJ IF JNEXT=0
CCCC          IF(JNEXT.EQ.0)GO TO 335
C**NEXT J
                J=JNEXT

C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM AND THIS JJ?
C               IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
C    1          MM.NE.MNEXT.OR.JJ.NE.JNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM AND JJ
C                 IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
C    1            N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN).OR.
C    2            (NSMODE.EQ.NAMODE))THEN
C                   IF(IREACT.EQ.0)CALL DUMRD5(KK,LL,NN,MM,JJ,IK2,IL2,
C    1              IN2,IM2,IJ2,W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),
C    2              W(LEJK5),J21,W(LMODNT))
C                 ELSE
CC                  CALL DUMRV5(KK,LL,NN,MM,JJ,IK2,IL2,IN2,IM2,IJ2,
CC   1              IK2TAU,
CC   1              W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),W(LEJK5),J21,
CC   2              W(LMODNT))
C                 END IF
C                 GO TO 335
C               END IF

C**INTRINSIC
C               IF(KCONT.LT.0)THEN
                  CALL GROUP5(KK,LL,NN,MM,JJ,ICONT(1),JCONT,1,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(1)
                    JCOUPC=NCOUPC(1)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 3035
                  END IF
                  CALL GROUP5(KK,LL,NN,MM,JJ,ICONT(2),JCONT,2,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(2)
                    JCOUPC=NCOUPC(2)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 3035
                  END IF
                  JCOUPL=JCOUPS
                  JCOUPC=JCOUCS
                  ICOUPL=ICOUPS
                  ICOUPC=ICOUCS
3035              CONTINUE
C               END IF

                IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.5)GO TO 335

C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM AND THIS JJ?
                IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
     1          MM.NE.MNEXT.OR.JJ.NE.JNEXT.OR.(KK.EQ.KNEXT.AND.
     2          LL.EQ.LNEXT.AND.NN.EQ.NNEXT.AND.MM.EQ.MNEXT.AND.
     3          JJ.EQ.JNEXT.AND.ICOUPL.LT.5))THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM AND JJ
                  IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1            N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN).OR.
     2            (NSMODE.EQ.NAMODE))THEN
                    IF(ICOUPX.GT.4.AND.IREACT.EQ.0)CALL DUMRD5(KK,LL,
     1              NN,MM,JJ,IK2,IL2,
     2              IN2,IM2,IJ2,W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),
     3              W(LEJK5),J21,W(LMODNT))
                  ELSE
CC                  CALL DUMRV5(KK,LL,NN,MM,JJ,IK2,IL2,IN2,IM2,IJ2,
CC   1              IK2TAU,
CC   1              W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),W(LEJK5),J21,
CC   2              W(LMODNT))
                  END IF
                  GO TO 335
                END IF

C**SKIP IF UNWANTED JJ
                IF(JNEXT.EQ.0.OR.JJ.NE.JNEXT)GO TO 335
C******************************************************************
C**CORIOLIS AND POTENTIAL
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN).OR.
     2          (NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR FIVE MODES
                  KPPP5=KIP5
                  CALL MEMO(1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
                  CALL GETBP5(W(LIP5),ISIZE5,NMODE,K,L,N,M,J,
     1            W(LMXBAS),NSIZE5,0)
                  KXA5=NSIZE5*(NSIZE5+1)/2
                  JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
                  JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
                  JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
                  JCI4=MAXBFN(W(LMXBAS),NMODE,M,1)
                  JCI5=MAXBFN(W(LMXBAS),NMODE,J,1)
                  I26=1
                  IF(ICOUPC.GT.4)I26=26
                  KTEMP1=JCI5*JCI5*I26
                  KXK0=JCI1*JCI1*IK2*I26
                  KXL0=JCI2*JCI2*IL2*I26
                  KXN0=JCI3*JCI3*IN2*I26
                  KXM0=JCI4*JCI4*IM2*I26
                  KXJ0=JCI5*JCI5*IJ2*I26
                  CALL MEMO(3,LXA5,KXA5,LTEMP1,KTEMP1,LTEMP7,NSIZE5,
     1            0,0,0,0)
                  CALL MEMO(5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,
     1            LXM0,KXM0,LXJ0,KXJ0)
C**ZEROISE MATRIX
                  CALL DIAGZ(NSIZE5,NSIZE5,W(LXA5),NSIZE5,W(LXA5),ICI)
C**GET ICOUNT2
                  CALL V0CI5(NAMODE,1,2,3,4,5,K,L,N,M,J,W(LH+K2),
     1            W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2            W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,
     3            IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,W(LXA5),NSIZE5,W(LIP5),
     4            ISIZE5,W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),JCI1,
     5            JCI2,JCI3,JCI4,JCI5,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6            W(LXJ0),W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),
     7            W(LEJK5),J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),
     8            ICOUNT3,W(LTEMP5),ICOUNT4,W(LTEMP6),W(LTEMP7),-2,I26)
                  KTEMP=ICOUNT2*2
                  KTEMP2=ICOUNT2*ICOUNT2*I26
                  CALL MEMO(2,LTEMP,KTEMP,LTEMP2,KTEMP2,0,0,0,0,0,0)
C**GET ICOUNT3
                  CALL V0CI5(NAMODE,1,2,3,4,5,K,L,N,M,J,W(LH+K2),
     1            W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2            W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,
     3            IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,W(LXA5),NSIZE5,W(LIP5),
     4            ISIZE5,W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),JCI1,
     5            JCI2,JCI3,JCI4,JCI5,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6            W(LXJ0),W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),
     7            W(LEJK5),J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),
     8            ICOUNT3,W(LTEMP5),ICOUNT4,W(LTEMP6),W(LTEMP7),-1,I26)
                  KTEMP5=ICOUNT3*2
                  KTEMP3=ICOUNT3*ICOUNT3*I26
                  CALL MEMO(2,LTEMP3,KTEMP3,LTEMP5,KTEMP5,0,0,0,0,0,0)
C**GET ICOUNT4
                  CALL V0CI5(NAMODE,1,2,3,4,5,K,L,N,M,J,W(LH+K2),
     1            W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2            W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,
     3            IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,W(LXA5),NSIZE5,W(LIP5),
     4            ISIZE5,W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),JCI1,
     5            JCI2,JCI3,JCI4,JCI5,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6            W(LXJ0),W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),
     7            W(LEJK5),J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),
     8            ICOUNT3,W(LTEMP5),ICOUNT4,W(LTEMP6),W(LTEMP7),1,I26)
                  KTEMP6=ICOUNT4*2
                  KTEMP4=ICOUNT4*ICOUNT4*I26
                  CALL MEMO(2,LTEMP4,KTEMP4,LTEMP6,KTEMP6,0,0,0,0,0,0)
C**GET XA5
                  CALL V0CI5(NAMODE,1,2,3,4,5,K,L,N,M,J,W(LH+K2),
     1            W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2            W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,
     3            IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,W(LXA5),NSIZE5,W(LIP5),
     4            ISIZE5,W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),JCI1,
     5            JCI2,JCI3,JCI4,JCI5,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     6            W(LXJ0),W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),
     7            W(LEJK5),J21,KROT,W(LMODNT),ICOUNT2,W(LTEMP),
     8            ICOUNT3,W(LTEMP5),ICOUNT4,W(LTEMP6),W(LTEMP7),0,I26)
                  CALL MATOUT(W(LXA5),W(LXA5),NSIZE5,25)
                  CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,
     1            LXM0,KXM0,LXJ0,KXJ0)
                  CALL MEMO(-3,LTEMP,KTEMP,LTEMP1,KTEMP1,LTEMP2,
     1            KTEMP2,0,0,0,0)
                  CALL MEMO(-2,LTEMP3,KTEMP3,LTEMP5,KTEMP5,0,0,0,0,0,0)
                  CALL MEMO(-2,LTEMP4,KTEMP4,LTEMP6,KTEMP6,0,0,0,0,0,0)
                  CALL MEMO(-2,LXA5,KXA5,LTEMP7,NSIZE5,0,0,0,0,0,0)
                  CALL MEMO(-1,LIP5,KIP5,0,0,0,0,0,0,0,0)
                ELSE
C*****************************   LIKE V0CI6 (K+L+N+M+J+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR SIX MODES
                  IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1              M.GT.NONLIN.OR.J.GT.NONLIN)THEN
                    KPPP6=KLJP6
                    IND6=1
                  ELSE
                    KPPP6=KIP6
                    IND6=0
                  END IF
                  CALL MEMO(1,LIP6,KPPP6,0,0,0,0,0,0,0,0)
                  IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1              M.GT.NONLIN.OR.J.GT.NONLIN)THEN
                    CALL GETBP5(W(LIP6),LSIZE6,NMODE,K,L,N,M,J,
     1              W(LMXBAS),NSIZE6,IND6)
                    MSIZE6=NSIZE6
                  ELSE
                    CALL GETBP6(W(LIP6),ISIZE6,NMODE,K,L,N,M,J,NSMODE,
     1              W(LMXBAS),NSIZE6,IND6)
                    MSIZE6=ISIZE6
                  END IF
C                 KXA6=NSIZE6*(NSIZE6+1)/2
                  IF(MOD(NSIZE6,2).EQ.0)THEN
                    KXA6=NSIZE6/2
                    KXA6=KXA6*(NSIZE6+1)
                  ELSE
                    KXA6=(NSIZE6+1)/2
                    KXA6=KXA6*NSIZE6
                  END IF
                  JCI1=MAXBFN(W(LMXBAS),NSMODE,K,1)
                  JCI2=MAXBFN(W(LMXBAS),NSMODE,L,1)
                  JCI3=MAXBFN(W(LMXBAS),NSMODE,N,1)
                  JCI4=MAXBFN(W(LMXBAS),NSMODE,M,1)
                  JCI5=MAXBFN(W(LMXBAS),NSMODE,J,1)
                  IF(JREACT.GT.0)THEN
                    JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1,NTAU(4)+1,
     1              NTAU(5)+1,NTAU(6)+1)
                  ELSE
                    JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1,NBAS(3,1)+1,
     1              NBAS(4,1)+1,NBAS(5,1)+1)-1
                    JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1,NBAS(3,2)+1,
     1              NBAS(4,2)+1,NBAS(5,2)+1)-1
                    JCIM=MAX0(JCIM1,JCIM2)
                  END IF
C**BEWARE LINEAR
C                 I36=1
C                 IF(ICOUPC.GT.4)I36=36
C**BEWARE LINEAR
                  KTEMP1=JCI5*JCI5
                  KXK0=JCI1*JCI1*IK2
                  KXL0=JCI2*JCI2*IL2
                  KXN0=JCI3*JCI3*IN2
                  KXM0=JCI4*JCI4*IM2
                  KXJ0=JCI5*JCI5*IJ2
                  KXP0=JCIM*JCIM*IK2TAU
                  CALL MEMO(3,LXA5,KXA6,LTEMP1,KTEMP1,LTEMP9,NSIZE6,
     1            0,0,0,0)
                  CALL MEMO(5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1            LXJ0,KXJ0)
                  CALL MEMO(1,LXP0,KXP0,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                  CALL DIAGZ(NSIZE6,NSIZE6,W(LXA5),NSIZE6,W(LXA5),ICI)
C**GET ICOUNT2
                  CALL V0CV5(NAMODE,1,2,3,4,5,6,K,L,N,M,J,W(LH+K2),
     1            W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2            W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+K2TAU),
     3            W(LXQ+K3TAU),W(LXA5),NSIZE6,IK3,IK2,IL3,IL2,IN3,IN2,
     4            IM3,IM2,IJ3,IJ2,IK3TAU,IK2TAU,W(LIP6),MSIZE6,
     5            W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),
     6            JCI1,JCI2,JCI3,JCI4,JCI5,JCIM,W(LXK0),W(LXL0),
     7            W(LXN0),W(LXM0),W(LXJ0),W(LXP0),W(LV5),W(LV5),W(LC5),
     8            W(LC5),W(LEJK5),W(LEJK5),J21,KROT,W(LMODNT),ICOUNT2,
     9            W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1            ICOUNT5,W(LTEMP8),W(LTEMP9),-3)
                  KTEMP=ICOUNT2*2
                  KTEMP2=ICOUNT2*ICOUNT2
                  CALL MEMO(2,LTEMP,KTEMP,LTEMP2,KTEMP2,0,0,0,0,0,0)
C**GET ICOUNT3
                  CALL V0CV5(NAMODE,1,2,3,4,5,6,K,L,N,M,J,W(LH+K2),
     1            W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2            W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+K2TAU),
     3            W(LXQ+K3TAU),W(LXA5),NSIZE6,IK3,IK2,IL3,IL2,IN3,IN2,
     4            IM3,IM2,IJ3,IJ2,IK3TAU,IK2TAU,W(LIP6),MSIZE6,
     5            W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),
     6            JCI1,JCI2,JCI3,JCI4,JCI5,JCIM,W(LXK0),W(LXL0),
     7            W(LXN0),W(LXM0),W(LXJ0),W(LXP0),W(LV5),W(LV5),W(LC5),
     8            W(LC5),W(LEJK5),W(LEJK5),J21,KROT,W(LMODNT),ICOUNT2,
     9            W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1            ICOUNT5,W(LTEMP8),W(LTEMP9),-2)
                  KTEMP6=ICOUNT3*2
                  KTEMP3=ICOUNT3*ICOUNT3
                  CALL MEMO(2,LTEMP3,KTEMP3,LTEMP6,KTEMP6,0,0,0,0,0,0)
C**GET ICOUNT4
                  CALL V0CV5(NAMODE,1,2,3,4,5,6,K,L,N,M,J,W(LH+K2),
     1            W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2            W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+K2TAU),
     3            W(LXQ+K3TAU),W(LXA5),NSIZE6,IK3,IK2,IL3,IL2,IN3,IN2,
     4            IM3,IM2,IJ3,IJ2,IK3TAU,IK2TAU,W(LIP6),MSIZE6,
     5            W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),
     6            JCI1,JCI2,JCI3,JCI4,JCI5,JCIM,W(LXK0),W(LXL0),
     7            W(LXN0),W(LXM0),W(LXJ0),W(LXP0),W(LV5),W(LV5),W(LC5),
     8            W(LC5),W(LEJK5),W(LEJK5),J21,KROT,W(LMODNT),ICOUNT2,
     9            W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1            ICOUNT5,W(LTEMP8),W(LTEMP9),-1)
                  KTEMP7=ICOUNT4*2
                  KTEMP4=ICOUNT4*ICOUNT4
                  CALL MEMO(2,LTEMP4,KTEMP4,LTEMP7,KTEMP7,0,0,0,0,0,0)
C**GET ICOUNT5
                  CALL V0CV5(NAMODE,1,2,3,4,5,6,K,L,N,M,J,W(LH+K2),
     1            W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2            W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+K2TAU),
     3            W(LXQ+K3TAU),W(LXA5),NSIZE6,IK3,IK2,IL3,IL2,IN3,IN2,
     4            IM3,IM2,IJ3,IJ2,IK3TAU,IK2TAU,W(LIP6),MSIZE6,
     5            W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),
     6            JCI1,JCI2,JCI3,JCI4,JCI5,JCIM,W(LXK0),W(LXL0),
     7            W(LXN0),W(LXM0),W(LXJ0),W(LXP0),W(LV5),W(LV5),W(LC5),
     8            W(LC5),W(LEJK5),W(LEJK5),J21,KROT,W(LMODNT),ICOUNT2,
     9            W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1            ICOUNT5,W(LTEMP8),W(LTEMP9),1)
                  KTEMP8=ICOUNT5*2
                  KTEMP5=ICOUNT5*ICOUNT5
                  CALL MEMO(2,LTEMP5,KTEMP5,LTEMP8,KTEMP8,0,0,0,0,0,0)
C**GET XA5
                  CALL V0CV5(NAMODE,1,2,3,4,5,6,K,L,N,M,J,W(LH+K2),
     1            W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2            W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+K2TAU),
     3            W(LXQ+K3TAU),W(LXA5),NSIZE6,IK3,IK2,IL3,IL2,IN3,IN2,
     4            IM3,IM2,IJ3,IJ2,IK3TAU,IK2TAU,W(LIP6),MSIZE6,
     5            W(LTEMP1),W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),
     6            JCI1,JCI2,JCI3,JCI4,JCI5,JCIM,W(LXK0),W(LXL0),
     7            W(LXN0),W(LXM0),W(LXJ0),W(LXP0),W(LV5),W(LV5),W(LC5),
     8            W(LC5),W(LEJK5),W(LEJK5),J21,KROT,W(LMODNT),ICOUNT2,
     9            W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1            ICOUNT5,W(LTEMP8),W(LTEMP9),0)
                  CALL MATOUT(W(LXA5),W(LXA5),NSIZE6,25)
                  CALL MEMO(-3,LXA5,KXA6,LTEMP1,KTEMP1,LTEMP9,NSIZE6,
     1            0,0,0,0)
                  CALL MEMO(-2,LTEMP,KTEMP,LTEMP2,KTEMP2,0,0,0,0,0,0)
                  CALL MEMO(-2,LTEMP3,KTEMP3,LTEMP6,KTEMP6,0,0,0,0,0,0)
                  CALL MEMO(-2,LTEMP4,KTEMP4,LTEMP7,KTEMP7,0,0,0,0,0,0)
                  CALL MEMO(-2,LTEMP5,KTEMP5,LTEMP8,KTEMP8,0,0,0,0,0,0)
                  CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1            LXJ0,KXJ0)
                  CALL MEMO(-1,LXP0,KXP0,0,0,0,0,0,0,0,0)
                  CALL MEMO(-1,LIP6,KPPP6,0,0,0,0,0,0,0,0)
C**
C*****************************   LIKE V0CI6 (K+L+N+M+J+NMODE)

                END IF
335   CONTINUE
C**INTRINSIC
                IF(ICOUPX.EQ.5)GO TO 3305
C**INTEGRATE OVER SIX NORMAL COORDINATES
                I=0
                I2=0
                I3=0
                IRI=1
                DO II=1,JJ-1
                  IF(II.EQ.1.AND.JJ.EQ.2.AND.MM.EQ.3.AND.NN.EQ.4.AND.
     1            LL.EQ.5.AND.KK.EQ.6.AND.ITIM.EQ.0)ITIM6A=0
                  IF(LCOUNT.GT.1)THEN
                    IF(IRI.LE.NAMODE)THEN
                      INEXT=JCONT(LCONT,IRI)
                      IF(II.EQ.INEXT)IRI=IRI+1
                    ELSE
                      INEXT=0
                    END IF
                  ELSE
                    INEXT=II
                  END IF
                  CALL INTARR(W(LNBF),W(LMBF),II,II1,II2,II3)
C**NEXT I
                  I=INEXT

C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM AND THIS JJ
C**AND THIS II?
C                 IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
C    1            MM.NE.MNEXT.OR.JJ.NE.JNEXT.OR.II.NE.INEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM AND JJ AND II
C                   IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.
C    1              AND.N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN.
C    2              AND.I.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C                     IF(IREACT.EQ.0)CALL DUMRD6(KK,LL,NN,MM,JJ,II,IK2,
C    1                IL2,IN2,IM2,IJ2,II2,W(LV6),W(LV6),W(LC6),W(LC6),
C    2                W(LEJK6),W(LEJK6),J21,W(LMODNT))
C                   END IF
C                   GO TO 336
C                 END IF

C**INTRINSIC
C                 IF(KCONT.LT.0)THEN
                    CALL GROUP6(K,L,N,M,J,I,ICONT(1),JCONT,1,IND)
                    IF(IND.NE.0)THEN
                      JCOUPL=NCOUPL(1)
                      JCOUPC=NCOUPC(1)
                      ICOUPL=IABS(JCOUPL)
                      ICOUPC=IABS(JCOUPC)
                      GO TO 3036
                    END IF
                    CALL GROUP6(K,L,N,M,J,I,ICONT(2),JCONT,2,IND)
                    IF(IND.NE.0)THEN
                      JCOUPL=NCOUPL(2)
                      JCOUPC=NCOUPC(2)
                      ICOUPL=IABS(JCOUPL)
                      ICOUPC=IABS(JCOUPC)
                      GO TO 3036
                    END IF
                    JCOUPL=JCOUPS
                    JCOUPC=JCOUCS
                    ICOUPL=ICOUPS
                    ICOUPC=ICOUCS
3036                CONTINUE
C                 END IF
                  IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.6)
     1            GO TO 336


C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM AND THIS JJ
C**AND THIS II?
                  IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
     1            MM.NE.MNEXT.OR.JJ.NE.JNEXT.OR.II.NE.INEXT.OR.
     2            (KK.EQ.KNEXT.AND.LL.EQ.LNEXT.AND.NN.EQ.NNEXT.AND.
     3            MM.EQ.MNEXT.AND.JJ.EQ.JNEXT.AND.II.EQ.INEXT.AND.
     4            ICOUPL.LT.6))THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM AND JJ AND II
                    IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.
     1              AND.N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN.
     2              AND.I.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
                      IF(ICOUPX.GT.5.AND.IREACT.EQ.0)
     1                CALL DUMRD6(KK,LL,NN,MM,JJ,II,IK2,
     2                IL2,IN2,IM2,IJ2,II2,W(LV6),W(LV6),W(LC6),W(LC6),
     3                W(LEJK6),W(LEJK6),J21,W(LMODNT))
                    END IF
                    GO TO 336
                  END IF

C**SKIP IF UNWANTED II
                  IF(INEXT.EQ.0.OR.II.NE.INEXT)GO TO 336
C******************************************************************
C**CORIOLIS AND POTENTIAL
                  IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1            N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN.AND.
     2            I.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR SIX MODES
                    KPPP6=KIP6
                    CALL MEMO(1,LIP6,KPPP6,0,0,0,0,0,0,0,0)
                    CALL GETBP6(W(LIP6),ISIZE6,NMODE,K,L,N,M,J,I,
     1              W(LMXBAS),NSIZE6,0)
                    KXA6=NSIZE6*(NSIZE6+1)/2
                    JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
                    JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
                    JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
                    JCI4=MAXBFN(W(LMXBAS),NMODE,M,1)
                    JCI5=MAXBFN(W(LMXBAS),NMODE,J,1)
                    JCI6=MAXBFN(W(LMXBAS),NMODE,I,1)
                    I37=1
                    IF(ICOUPC.GT.5)I37=37
                    KTEMP1=JCI6*JCI6*I37
                    KXK0=JCI1*JCI1*IK2*I37
                    KXL0=JCI2*JCI2*IL2*I37
                    KXN0=JCI3*JCI3*IN2*I37
                    KXM0=JCI4*JCI4*IM2*I37
                    KXJ0=JCI5*JCI5*IJ2*I37
                    KXI0=JCI6*JCI6*II2*I37
                    CALL MEMO(3,LXA6,KXA6,LTEMP1,KTEMP1,LTEMP9,
     1              NSIZE6,0,0,0,0)
                    CALL MEMO(5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,
     1              LXM0,KXM0,LXJ0,KXJ0)
                    CALL MEMO(1,LXI0,KXI0,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                    CALL DIAGZ(NSIZE6,NSIZE6,W(LXA6),NSIZE6,W(LXA6),
     1              ICI)
C**GET ICOUNT2
                    CALL V0CI6(NAMODE,1,2,3,4,5,6,K,L,N,M,J,I,W(LH+K2),
     1              W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2              W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+I2),
     3              W(LXQ+I3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,
     4              II3,II2,W(LXA6),NSIZE6,W(LIP6),ISIZE6,W(LTEMP1),
     5              W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),JCI1,JCI2,
     6              JCI3,JCI4,JCI5,JCI6,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     7              W(LXJ0),W(LXI0),W(LV6),W(LV6),W(LC6),W(LC6),
     8              W(LEJK6),W(LEJK6),J21,KROT,W(LMODNT),ICOUNT2,
     9              W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1              ICOUNT5,W(LTEMP8),W(LTEMP9),-3,I37)
                    KTEMP=ICOUNT2*2
                    KTEMP2=ICOUNT2*ICOUNT2*I37
                    CALL MEMO(2,LTEMP,KTEMP,LTEMP2,KTEMP2,0,0,0,0,0,0)
C**GET ICOUNT3
                    CALL V0CI6(NAMODE,1,2,3,4,5,6,K,L,N,M,J,I,W(LH+K2),
     1              W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2              W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+I2),
     3              W(LXQ+I3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,
     4              II3,II2,W(LXA6),NSIZE6,W(LIP6),ISIZE6,W(LTEMP1),
     5              W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),JCI1,JCI2,
     6              JCI3,JCI4,JCI5,JCI6,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     7              W(LXJ0),W(LXI0),W(LV6),W(LV6),W(LC6),W(LC6),
     8              W(LEJK6),W(LEJK6),J21,KROT,W(LMODNT),ICOUNT2,
     9              W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1              ICOUNT5,W(LTEMP8),W(LTEMP9),-2,I37)
                    KTEMP6=ICOUNT3*2
                    KTEMP3=ICOUNT3*ICOUNT3*I37
                    CALL MEMO(2,LTEMP3,KTEMP3,LTEMP6,KTEMP6,0,0,0,0,
     1              0,0)
C**GET ICOUNT4
                    CALL V0CI6(NAMODE,1,2,3,4,5,6,K,L,N,M,J,I,W(LH+K2),
     1              W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2              W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+I2),
     3              W(LXQ+I3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,
     4              II3,II2,W(LXA6),NSIZE6,W(LIP6),ISIZE6,W(LTEMP1),
     5              W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),JCI1,JCI2,
     6              JCI3,JCI4,JCI5,JCI6,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     7              W(LXJ0),W(LXI0),W(LV6),W(LV6),W(LC6),W(LC6),
     8              W(LEJK6),W(LEJK6),J21,KROT,W(LMODNT),ICOUNT2,
     9              W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1              ICOUNT5,W(LTEMP8),W(LTEMP9),-1,I37)
                    KTEMP7=ICOUNT4*2
                    KTEMP4=ICOUNT4*ICOUNT4*I37
                    CALL MEMO(2,LTEMP4,KTEMP4,LTEMP7,KTEMP7,0,0,0,0,
     1              0,0)
C**GET ICOUNT5
                    CALL V0CI6(NAMODE,1,2,3,4,5,6,K,L,N,M,J,I,W(LH+K2),
     1              W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2              W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+I2),
     3              W(LXQ+I3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,
     4              II3,II2,W(LXA6),NSIZE6,W(LIP6),ISIZE6,W(LTEMP1),
     5              W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),JCI1,JCI2,
     6              JCI3,JCI4,JCI5,JCI6,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     7              W(LXJ0),W(LXI0),W(LV6),W(LV6),W(LC6),W(LC6),
     8              W(LEJK6),W(LEJK6),J21,KROT,W(LMODNT),ICOUNT2,
     9              W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1              ICOUNT5,W(LTEMP8),W(LTEMP9),1,I37)
                    KTEMP8=ICOUNT5*2
                    KTEMP5=ICOUNT5*ICOUNT5*I37
                    CALL MEMO(2,LTEMP5,KTEMP5,LTEMP8,KTEMP8,0,0,0,0,
     1              0,0)
C**GET XA6
                    CALL V0CI6(NAMODE,1,2,3,4,5,6,K,L,N,M,J,I,W(LH+K2),
     1              W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     2              W(LH+M2),W(LXQ+M3),W(LH+J2),W(LXQ+J3),W(LH+I2),
     3              W(LXQ+I3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,IJ3,IJ2,
     4              II3,II2,W(LXA6),NSIZE6,W(LIP6),ISIZE6,W(LTEMP1),
     5              W(LTEMP2),W(LTEMP3),W(LTEMP4),W(LTEMP5),JCI1,JCI2,
     6              JCI3,JCI4,JCI5,JCI6,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     7              W(LXJ0),W(LXI0),W(LV6),W(LV6),W(LC6),W(LC6),
     8              W(LEJK6),W(LEJK6),J21,KROT,W(LMODNT),ICOUNT2,
     9              W(LTEMP),ICOUNT3,W(LTEMP6),ICOUNT4,W(LTEMP7),
     1              ICOUNT5,W(LTEMP8),W(LTEMP9),0,I37)
                    CALL MATOUT(W(LXA6),W(LXA6),NSIZE6,26)
                    CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,
     1              LXM0,KXM0,LXJ0,KXJ0)
                    CALL MEMO(-1,LXI0,KXI0,0,0,0,0,0,0,0,0)
                    CALL MEMO(-3,LTEMP,KTEMP,LTEMP1,KTEMP1,LTEMP2,
     1              KTEMP2,0,0,0,0)
                    CALL MEMO(-2,LTEMP3,KTEMP3,LTEMP6,KTEMP6,0,0,0,0,
     1              0,0)
                    CALL MEMO(-2,LTEMP4,KTEMP4,LTEMP7,KTEMP7,0,0,0,0,
     1              0,0)
                    CALL MEMO(-2,LTEMP5,KTEMP5,LTEMP8,KTEMP8,0,0,0,0,
     1              0,0)
                    CALL MEMO(-2,LXA6,KXA6,LTEMP9,NSIZE6,0,0,0,0,0,0)
                    CALL MEMO(-1,LIP6,KIP6,0,0,0,0,0,0,0,0)
                  END IF
336   CONTINUE
3306  CONTINUE
                  I2=I2+II1
                  IF(I.GT.NONLIN)I2=I2+II1
                  I3=I3+II2
                END DO
3305  CONTINUE
                J2=J2+IJ1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
                IF(J.GT.NONLIN)J2=J2+IJ1
                J3=J3+IJ2
              END DO
3304  CONTINUE
              M2=M2+IM1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
              IF(M.GT.NONLIN)M2=M2+IM1
              M3=M3+IM2
            END DO
3303  CONTINUE
            N2=N2+IN1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
            IF(N.GT.NONLIN)N2=N2+IN1
            N3=N3+IN2
          END DO
3302  CONTINUE
          L2=L2+IL1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(L.GT.NONLIN)L2=L2+IL1
          L3=L3+IL2
        END DO
3301  CONTINUE
        K2=K2+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
        IF(K.GT.NONLIN)K2=K2+IK1
        K3=K3+IK2
      END DO
C**INTRINSIC
C     IF(KCONT.LT.0)THEN
C       JCOUPL=JCOUPS
C       JCOUPC=JCOUCS
C       ICOUPL=ICOUPS
C       ICOUPC=ICOUCS
C     END IF
3332  CONTINUE
3330  CONTINUE
      ITIM=-1
3331  CONTINUE

C****************************************************************
C****************************************************************

C**********************************************************************
C**                               PUT COMPLETE SELECTIVE BASIS TOGETHER
C**********************************************************************
      IF(NNMAX.LT.0.AND.KCONT.LT.0)THEN
C**PUT FULL SELECTIVE BASIS TOGETHER
        NVSYM=NVSYMX
        ISIZE=1
C**LOOP ROUND CONTRACTION SCHEMES
        DO K=1,NCONT
          ISIZE=ISIZE*(ITOT1(K)+ITOT2(K)+ITOT3(K)+ITOT4(K)+ITOT5(K)+
     1    ITOT6(K))
          KTOT(6,K)=ITOT6(K)
          KTOT(5,K)=ITOT5(K)
          KTOT(4,K)=ITOT4(K)
          KTOT(3,K)=ITOT3(K)
          KTOT(2,K)=ITOT2(K)
          KTOT(1,K)=ITOT1(K)
          IF(K.EQ.1)THEN
            KIP=ISIZE*NNMODE
            KJP=KIP
            JSIZE=ISIZE
          END IF
          IF(K.EQ.NCONT.AND.NCONT.GT.1)THEN
C**GET INITIAL SIZE OVER BOTH CONTRACTION BASES (ALLOW NMAX = -6)
C**ALLOW ONLY FOR (-NMAX)-MODE BASIS HERE...THE REST WILL FOLLOW
C**KEEP 6-MODE WITH ZPE
            DO I1=1,6
              K1=KTOT(I1,1)
              IF(I1.EQ.1)K1=K1-1
              DO I2=7-I1,6
                K2=KTOT(I2,2)
                IF(I2.EQ.1)K2=K2-1
                ISIZE=ISIZE-K1*K2
              END DO
            END DO
            IF(ICI.GT.-6)THEN
C**REMOVE 6-MODE WITH ZPE
              ISIZE=ISIZE-KTOT(6,1)-KTOT(6,2)
C**KEEP 5-MODE WITH ZPE
              DO I1=1,5
                K1=KTOT(I1,1)
                IF(I1.EQ.1)K1=K1-1
                I2=6-I1
                K2=KTOT(I2,2)
                IF(I2.EQ.1)K2=K2-1
                ISIZE=ISIZE-K1*K2
              END DO
              IF(ICI.GT.-5)THEN
C**REMOVE 5-MODE WITH ZPE
                ISIZE=ISIZE-KTOT(5,1)-KTOT(5,2)
C**KEEP 4-MODE WITH ZPE
                DO I1=1,4
                  K1=KTOT(I1,1)
                  IF(I1.EQ.1)K1=K1-1
                  I2=5-I1
                  K2=KTOT(I2,2)
                  IF(I2.EQ.1)K2=K2-1
                  ISIZE=ISIZE-K1*K2
                END DO
                IF(ICI.GT.-4)THEN
C**REMOVE 4-MODE WITH ZPE
                  ISIZE=ISIZE-KTOT(4,1)-KTOT(4,2)
C**KEEP 3-MODE WITH ZPE
                  DO I1=1,3
                    K1=KTOT(I1,1)
                    IF(I1.EQ.1)K1=K1-1
                    I2=4-I1
                    K2=KTOT(I2,2)
                    IF(I2.EQ.1)K2=K2-1
                    ISIZE=ISIZE-K1*K2
                  END DO
                  IF(ICI.GT.-3)THEN
C**REMOVE 3-MODE WITH ZPE
                    ISIZE=ISIZE-KTOT(3,1)-KTOT(3,2)
C**KEEP 2-MODE WITH ZPE
                    DO I1=1,2
                      K1=KTOT(I1,1)
                      IF(I1.EQ.1)K1=K1-1
                      I2=3-I1
                      K2=KTOT(I2,2)
                      IF(I2.EQ.1)K2=K2-1
                      ISIZE=ISIZE-K1*K2
                    END DO
                    IF(ICI.GT.-2)THEN
C**REMOVE 2-MODE WITH ZPE
                      ISIZE=ISIZE-KTOT(2,1)-KTOT(2,2)
C**KEEP 1-MODE WITH ZPE
                      DO I1=1,1
                        K1=KTOT(I1,1)
                        IF(I1.EQ.1)K1=K1-1
                        I2=2-I1
                        K2=KTOT(I2,2)
                        IF(I2.EQ.1)K2=K2-1
                        ISIZE=ISIZE-K1*K2
                      END DO
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
C**END LOOP ROUND CONTRACTION SCHEMES
        END DO
        CALL MEMO(1,LIP,KIP,0,0,0,0,0,0,0,0)
C**LOOP TWICE........FIRST TIME, FILL IN SCHEME (1) ONLY AND MODIFY
C**                  ITOT ACCORDING TO NMAX
C**          ........SECOND TIME, GET ARRAYS FOR BOTH SCHEMES AND FILL 
C**                  THEM IN
        DO KBAS=1,NCONT
C**PUT ZERO QUANTA EVERYWHERE
          CALL GETJP0(W(LIP),JSIZE,NNMODE)
C**SET UP CONTRACTION SCHEME(1) BASIS INITIALLY
          NCMODE=ICONT(1)
          IF(ICI.LT.0.AND.KJP1.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP1),JSIZE1(1),NCMODE,ITOT,
     2    ITOT1(1),1)
          IF(ICI.LT.-1.AND.KJP2.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP2),JSIZE2(1),NCMODE,ITOT,
     2    ITOT2(1),1)
          IF(ICI.LT.-2.AND.KJP3.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP3),JSIZE3(1),NCMODE,ITOT,
     2    ITOT3(1),1)
          IF(ICI.LT.-3.AND.KJP4.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP4),JSIZE4(1),NCMODE,ITOT,
     2    ITOT4(1),1)
          IF(ICI.LT.-4.AND.KJP5.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP5),JSIZE5(1),NCMODE,ITOT,
     2    ITOT5(1),1)
          IF(ICI.LT.-5.AND.KJP6.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP6),JSIZE6(1),NCMODE,ITOT,
     2    ITOT6(1),1)
C**ADD SECOND CONTRACTION SCHEME IF EXISTS
          IF(NCONT.GT.1)THEN
C**PUT IN 1-MODE FOR SCHEME (2)
            IF(ICI.LT.0)
     1      CALL GETCP(W(LIP),JSIZE,W(LJP1),W(LJP2),W(LJP3),W(LJP4),
     2      W(LJP5),W(LJP6),JSIZE1(1),JSIZE2(1),JSIZE3(1),JSIZE4(1),
     3      JSIZE5(1),JSIZE6(1),W(LJP1+JSIZE1(1)*NCMODE),JSIZE1(2),
     4      ITOT,ITOT1(1),ITOT2(1),ITOT3(1),ITOT4(1),ITOT5(1),ITOT6(1),
     5      ITOT1(2)-1,ICI,1,0,KBAS-1,W(LMXBAS),NMODE)
C**PUT IN 2-MODE FOR SCHEME (2)
            IF(ICI.LT.-1)
     1      CALL GETCP(W(LIP),JSIZE,W(LJP1),W(LJP2),W(LJP3),W(LJP4),
     2      W(LJP5),W(LJP6),JSIZE1(1),JSIZE2(1),JSIZE3(1),JSIZE4(1),
     3      JSIZE5(1),JSIZE6(1),W(LJP2+JSIZE2(1)*NCMODE),JSIZE2(2),
     4      ITOT,ITOT1(1),ITOT2(1),ITOT3(1),ITOT4(1),ITOT5(1),ITOT6(1),
     5      ITOT2(2),ICI,0,1,KBAS-1,W(LMXBAS),NMODE)
C**PUT IN 3-MODE FOR SCHEME (2)
            IF(ICI.LT.-2)
     1      CALL GETCP(W(LIP),JSIZE,W(LJP1),W(LJP2),W(LJP3),W(LJP4),
     2      W(LJP5),W(LJP6),JSIZE1(1),JSIZE2(1),JSIZE3(1),JSIZE4(1),
     3      JSIZE5(1),JSIZE6(1),W(LJP3+JSIZE3(1)*NCMODE),JSIZE3(2),
     4      ITOT,ITOT1(1),ITOT2(1),ITOT3(1),ITOT4(1),ITOT5(1),ITOT6(1),
     5      ITOT3(2),ICI,0,2,KBAS-1,W(LMXBAS),NMODE)
C**PUT IN 4-MODE FOR SCHEME (2)
            IF(ICI.LT.-3)
     1      CALL GETCP(W(LIP),JSIZE,W(LJP1),W(LJP2),W(LJP3),W(LJP4),
     2      W(LJP5),W(LJP6),JSIZE1(1),JSIZE2(1),JSIZE3(1),JSIZE4(1),
     3      JSIZE5(1),JSIZE6(1),W(LJP4+JSIZE4(1)*NCMODE),JSIZE4(2),
     4      ITOT,ITOT1(1),ITOT2(1),ITOT3(1),ITOT4(1),ITOT5(1),ITOT6(1),
     5      ITOT4(2),ICI,0,3,KBAS-1,W(LMXBAS),NMODE)
C**PUT IN 5-MODE FOR SCHEME (2)
            IF(ICI.LT.-4)
     1      CALL GETCP(W(LIP),JSIZE,W(LJP1),W(LJP2),W(LJP3),W(LJP4),
     2      W(LJP5),W(LJP6),JSIZE1(1),JSIZE2(1),JSIZE3(1),JSIZE4(1),
     3      JSIZE5(1),JSIZE6(1),W(LJP5+JSIZE5(1)*NCMODE),JSIZE5(2),
     4      ITOT,ITOT1(1),ITOT2(1),ITOT3(1),ITOT4(1),ITOT5(1),ITOT6(1),
     5      ITOT5(2),ICI,0,4,KBAS-1,W(LMXBAS),NMODE)
C**PUT IN 6-MODE FOR SCHEME (2)
            IF(ICI.LT.-5)
     1      CALL GETCP(W(LIP),JSIZE,W(LJP1),W(LJP2),W(LJP3),W(LJP4),
     2      W(LJP5),W(LJP6),JSIZE1(1),JSIZE2(1),JSIZE3(1),JSIZE4(1),
     3      JSIZE5(1),JSIZE6(1),W(LJP6+JSIZE6(1)*NCMODE),JSIZE6(2),
     4      ITOT,ITOT1(1),ITOT2(1),ITOT3(1),ITOT4(1),ITOT5(1),ITOT6(1),
     5      ITOT6(2),ICI,0,5,KBAS-1,W(LMXBAS),NMODE)
          END IF
          IF(KBAS.LT.NCONT)THEN
            CALL MEMO(-1,LIP,KIP,0,0,0,0,0,0,0,0)
            ISIZE=ITOT
            JSIZE=ISIZE
            KIP=ISIZE*NNMODE
            KJP=KIP
            CALL MEMO(1,LIP,KIP,0,0,0,0,0,0,0,0)
            ITOT=0
          END IF
        END DO
C**END OF KBAS LOOP
        IF(ICI.LT.0.AND.KJP1.GT.0)
     1  CALL MEMO(-1,LJP1,KJP1,0,0,0,0,0,0,0,0)
        IF(ICI.LT.-1.AND.KJP2.GT.0)
     1  CALL MEMO(-1,LJP2,KJP2,0,0,0,0,0,0,0,0)
        IF(ICI.LT.-2.AND.KJP3.GT.0)
     1  CALL MEMO(-1,LJP3,KJP3,0,0,0,0,0,0,0,0)
        IF(ICI.LT.-3.AND.KJP4.GT.0)
     1  CALL MEMO(-1,LJP4,KJP4,0,0,0,0,0,0,0,0)
        IF(ICI.LT.-4.AND.KJP5.GT.0)
     1  CALL MEMO(-1,LJP5,KJP5,0,0,0,0,0,0,0,0)
        IF(ICI.LT.-5.AND.KJP6.GT.0)
     1  CALL MEMO(-1,LJP6,KJP6,0,0,0,0,0,0,0,0)
        CALL MEMO(1,LJP,KJP,0,0,0,0,0,0,0,0)
        CALL GETIP(W(LIP),W(LJP),ISIZE,ISIZE,JSIZE,NNMODE,0,0,0)
        WRITE(IOUT,395)(NTOT(K),K=1,NVSYM)
        IF(ITOT.EQ.0)THEN
          WRITE(IOUT,400)
          STOP 'SHOULD NOT OCCUR'
        END IF
      ELSE
        IF(NNMAX.LT.0)THEN
CCCC      IF(LCONT.GT.1)CALL MEMO(-2,LIP,KIP,LJP,KJP,0,0,0,0,0,0)
C**COPY INDIVIDUAL CONTRACTION BASES ACROSS
          NCMODE=0
          IF(LCONT.EQ.2)NCMODE=ICONT(1)
          NVSYM=NVSYMX
          ISIZE=ITOT1(LCONT)+ITOT2(LCONT)+ITOT3(LCONT)+ITOT4(LCONT)+
     1    ITOT5(LCONT)+ITOT6(LCONT)
          JSIZE=ISIZE
          KIP=ISIZE*NNMODE
          KJP=KIP
          CALL MEMO(1,LIP,KIP,0,0,0,0,0,0,0,0)
C**PUT ZERO QUANTA EVERYWHERE
          CALL GETJP0(W(LIP),JSIZE,NNMODE)
          ITOT=0
          IF(ICI.LT.0.AND.KJP1.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP1+JSIZE1(1)*NCMODE),
     2    JSIZE1(LCONT),NNMODE,ITOT,ITOT1(LCONT),0)
          IF(ICI.LT.-1.AND.KJP2.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP2+JSIZE2(1)*NCMODE),
     2    JSIZE2(LCONT),NNMODE,ITOT,ITOT2(LCONT),0)
          IF(ICI.LT.-2.AND.KJP3.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP3+JSIZE3(1)*NCMODE),
     2    JSIZE3(LCONT),NNMODE,ITOT,ITOT3(LCONT),0)
          IF(ICI.LT.-3.AND.KJP4.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP4+JSIZE4(1)*NCMODE),
     2    JSIZE4(LCONT),NNMODE,ITOT,ITOT4(LCONT),0)
          IF(ICI.LT.-4.AND.KJP5.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP5+JSIZE5(1)*NCMODE),
     2    JSIZE5(LCONT),NNMODE,ITOT,ITOT5(LCONT),0)
          IF(ICI.LT.-5.AND.KJP6.GT.0)
     1    CALL GETJP(W(LIP),JSIZE,W(LJP6+JSIZE6(1)*NCMODE),
     2    JSIZE6(LCONT),NNMODE,ITOT,ITOT6(LCONT),0)
          CALL MEMO(1,LJP,KJP,0,0,0,0,0,0,0,0)
          CALL GETIP(W(LIP),W(LJP),ISIZE,ISIZE,JSIZE,NNMODE,0,0,LCONT)
          IF(LCONT.EQ.1)THEN
            KLP=KJP
            CALL MEMO(1,LLP,KLP,0,0,0,0,0,0,0,0)
            DO I=1,KLP
              W(LLP-1+I)=W(LJP-1+I)
            END DO
C**JUST IN CASE SINGLE CONTRACTION SCHEME REQUESTED (BAD USE OF 'MM')
            IF(LCOUNT.EQ.2.AND.ICONT(2).EQ.0)THEN
              KKP=KJP
              CALL MEMO(1,LKP,KLP+KKP,0,0,0,0,0,0,0,0)
C**SAVE BASIS FOR FIRST CONTRACTION SCHEME
              DO I=1,KLP
                W(LKP-1+I)=W(LLP-1+I)
              END DO
            END IF
          ELSE
            KKP=KJP
            CALL MEMO(1,LKP,KLP+KKP,0,0,0,0,0,0,0,0)
C**SAVE TOTAL BASIS FOR FIRST CONTRACTION SCHEME
            DO I=1,KLP
              W(LKP-1+I)=W(LLP-1+I)
            END DO
C**SAVE TOTAL BASIS FOR SECOND CONTRACTION SCHEME
            DO I=1,KKP
              W(LKP-1+I+KLP)=W(LJP-1+I)
            END DO
            CALL MEMO(-1,LLP,KLP,0,0,0,0,0,0,0,0)
          END IF
          WRITE(IOUT,395)(NTOT(K),K=1,NVSYM)
          IF(ITOT.EQ.0)THEN
            WRITE(IOUT,400)
            STOP 'SHOULD NOT OCCUR'
          END IF
C**ICSIZE IS TOTAL SIZE OF BASIS FOR CONTRACTION SCHEMES
          MAXSIZ=0
          DO I=1,NVSYM
            IF(NTOT(I).GT.MAXSIZ)MAXSIZ=NTOT(I)
          END DO
          IF(LCONT.EQ.1)ICSIZ1=MAXSIZ
          IF(LCONT.EQ.2)ICSIZ2=MAXSIZ
          IF(LCONT.EQ.1)IPSIZ1=ISIZE
          IF(LCONT.EQ.2)IPSIZ2=ISIZE
        END IF
      END IF

C****************************************************************
C****************************************************************

C**********************************************************************
C**                                            SET UP K-DIAGONAL MATRIX
C**********************************************************************
C     IF(MATSIZ.NE.0.AND.LLNCUT.EQ.0)RETURN
      IF(MATSIZ.NE.0)RETURN
      CALL MEMO(-1,LIP,KIP,0,0,0,0,0,0,0,0)
C**MAXIMUM NTOT(NS)
      ISIZMX=0
      DO NS=1,NVSYM
        IF(NTOT(NS).GT.ISIZMX)ISIZMX=NTOT(NS)
      END DO
      CALL MEMO(1,LASSIG,ISIZMX*KEL21*NVSYM*3,0,0,0,0,0,0,0,0)
C**MINIMUM NTOT(NS)
      ISIZMN=200000
      DO NS=1,NVSYM
        IF(NTOT(NS).LT.ISIZMN.AND.NTOT(NS).NE.0)ISIZMN=NTOT(NS)
      END DO
C**FOR J>0 OR CONTRACTIONS, NEED TO SAVE COEFFICIENTS, DIMENSION NVALCF
      IF(LCOUNT.EQ.2.OR.JTHIS.NE.0)THEN
        IF(LCONT.EQ.1)THEN
          IF(ISIZMX.LT.NVAL1)NVAL1=ISIZMX
          NVALCF=NVAL1
        END IF
        IF(LCONT.EQ.2)THEN
          IF(ISIZMX.LT.NVAL2)NVAL2=ISIZMX
          NVALCF=NVAL2
        END IF
      ELSE
        NVAL2=0
      END IF
C**MAXIMUM OF ISIZM1,ISIZM2 USED IN COMBINATION OF 2 CONTRACTIONS
      IF(LCONT.EQ.1)ISIZM1=ISIZMX
      ISIZM2=0
      IF(LCONT.EQ.2)ISIZM2=ISIZMX
      KIP=ISIZMX*NNMODE
      CALL MEMO(1,LIP,KIP,0,0,0,0,0,0,0,0)
      IF(JTHIS.NE.0.OR.LCOUNT.GT.1)THEN
        CALL MEMO(2,LCFS,ISIZMX*NVALCF*KEL21*NVSYM,LEVCI,
     1  NVALCF*KEL21*NVSYM,0,0,0,0,0,0)
      END IF
C**SELECTIVE SYMMETRIES IF J=0
      MMSS=MVSYM
      DO I=1,MVSYM
        INTSM(I)=MWSYM(I)
      END DO
C**MUST DO ALL 'VIBRATIONAL' SYMMETRIES IF J>0 OR FURTHER CONTRACTIONS
C**(DEPENDING ON SIGN OF INPUT NVAL1,NVAL2)
      IF(JTHIS.NE.0.OR.(LCOUNT.GT.1.AND.MVAL1.GE.0.AND.MVAL2.GE.0))THEN
        MMSS=NVSYM
        DO I=1,NVSYM
          INTSM(I)=I
        END DO
      END IF
      IF(ISCFCI.GT.0.AND.IDUMP.NE.0)WRITE(INP4)JMAX,KSTEP,
C    1MMSS,(INTSM(I),I=1,MMSS),
     1MVSYM,(INTSM(I),I=1,MVSYM),
     2ICI,NNMAX,NCONT,
     3((NBAS(I,J),I=1,6),J=1,2),((MAXSUM(I,J),I=1,6),J=1,2),
     4((MAXBAS(I,J),I=1,NSMODE),J=1,IABS(NNMAX))
C**********************************************************************
C**                                                LOOP OVER SYMMETRIES
C**********************************************************************
C**ROTATION MATRIX IS K TIMES MINIMUM VALUE OF NVAL (NO CONTRACTIONS)
      IF(LCONT.EQ.1)NVALV=NVAL1
      IF(LCONT.EQ.2)NVALV=NVAL2
      IF(JTHIS.NE.0.AND.LCOUNT.EQ.1)THEN
        IF(NVALV.GT.ISIZMN)NVALV=ISIZMN
        NVAL1=NVALV
        WRITE(IOUT,*)
        WRITE(IOUT,*)'ROTATION USES NVALV = ',NVALV
      END IF
      DO 4400 MS=1,MMSS
      IF(JTHIS.EQ.0.AND.(LCOUNT.NE.2.OR.MVAL1.LT.0.OR.MVAL2.LT.0))THEN
        NS=MWSYM(MS)
      ELSE
        NS=MS
      END IF
CCCC  WRITE(IOUT,290)NS
      IF(NS.EQ.0)GO TO 4401
      ITIM=ITIM+1
      ISIZE=NTOT(NS)
      ISIZEX=ISIZE
      IF(LCOUNT.GT.1)ISIZC(LCONT,NS)=ISIZE
C     WRITE(IOUT,285)ISIZE,CUT
      CALL FLUSH(IOUT)
      IF(ISIZE.NE.0)THEN
C**MAKE SURE NVAL WITHIN RANGE OF MATRIX SIZE
        IF(LCONT.EQ.1)NCVAL(LCONT,NS)=NVAL1
        IF(LCONT.EQ.2)NCVAL(LCONT,NS)=NVAL2
        IF(NCVAL(LCONT,NS).GT.ISIZE)NCVAL(LCONT,NS)=ISIZE
        NVAL=NCVAL(LCONT,NS)
C**TEMPORARY
C       IF(NVAL.EQ.ISIZE)THEN
C         NVAL=NVAL-1
C         NCVAL(LCONT,NS)=NVAL
C       END IF
C**TEMPORARY
C**NVALV IS NUMBER DUMPED FOR EACH CONTRACTION SCHEME
C       IF(LCOUNT.GT.1)NVALV=NVAL
        NVALV=NVAL
C       WRITE(IOUT,502)LCONT,NVAL
C**SORT INTO SYMMETRY BASES
        CALL PUTJP(W(LJP),JSIZE,W(LIP),ISIZMX,NNMODE,NS)
      END IF
C**DUMP BASIS TO (INP4)
      IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
        CALL OUTBAS(ISIZE,W(LIP),ISIZMX,NNMODE)
      END IF
      IF(ISIZE.EQ.0)GO TO 4401

C**************************************************
      IF(.NOT.LANCZ)THEN
        IF(ICOUPX.GE.0)THEN
          REWIND 21
        END IF
        IF(ICOUPX.GT.1)THEN
          REWIND 22
        END IF
        IF(ICOUPX.GT.2)THEN
          REWIND 23
        END IF
        IF(ICOUPX.GT.3)THEN
          REWIND 24
        END IF
        IF(ICOUPX.GT.4)THEN
          REWIND 25
        END IF
        IF(ICOUPX.GT.5)THEN
          REWIND 26
        END IF
      END IF
C**************************************************

      KCORIG=JTHIS
      KEL=0
C**********************************************************************
C**                                                       LOOP ROUND Ka
C**********************************************************************
      DO 9000 KROT=1,J21,KSTEP
      KKROT=KROT
      KEL=KEL+1
C**KA,KC FOR VCI
      KA=KROT/2
      KC=KCORIG-KA
      IF(KMAX.GT.0)THEN
        IF(KA.GT.0.AND.MOD(KA,2).EQ.MOD(KROT,2))KC=KC+1
C       IF(MOD(KA,2).NE.0)THEN
C         KC=KC+MOD(KROT+1,2)
C       ELSE
CC        KC=KC+MOD(KROT+1,2)
C         IF(KA.NE.0)KC=KC+MOD(KROT+1,2)
C       END IF
      ELSE
        IF(MOD(KROT,2).EQ.0)KC=KC+1
      END IF
      KROTT=KROT
1     IF(KROTT.GT.4)THEN
        KROTT=KROTT-4
        GO TO 1
      END IF
      JS=MOD(JTHIS,2)+1
      IF(MOD(KA,2).EQ.0)THEN
        WRITE(IOUT,290)NS,ISYMP(NS,ISYMD(KROTT,JS))
      ELSE
        WRITE(IOUT,290)NS,ISYMP(ISYMP(NS,ISYMT),ISYMD(KROTT,JS))
      END IF
C**ISTART, IEND ARE START AND END COLUMNS
      ISTART=1
C**SET DEFAULT MARKER (LANIND) FOR VCIn AND VCVn
      LANIND=0
      IF(LANCZ)THEN
        LAN12D=0
        LANCUT=IABS(LLNCUT)
        IF(LANCUT.NE.0)THEN
          LANIND=1
          NVALSV=NVAL
C**DEFAULT 1 & 2 MODE BASIS
          LAN12D=NTOTL1(NS)+NTOTL2(NS)
C**SPECIAL CASE MAY NEED 3 MODE BASIS
          IF(TLLCUT.LT.0)LAN12D=LAN12D+NTOTL3(NS)
          ISIZE=LAN12D
        END IF
        CALL MEMO(2,LWK,ISIZEX,LYK,ISIZEX,0,0,0,0,0,0)
        CALL MEMO(2,LANCNT,ISIZEX,LANK,ISIZEX,0,0,0,0,0,0)
        IF(LLAN20.LT.0)GO TO 700
        CALL MEMO(1,LZK,ISIZEX,0,0,0,0,0,0,0,0)
        KLAN=LANMAX*(LANMAX+1)/2
        LSIZE=1
      ELSE
C**TEMPORARY - HERMITIAN
        KXA=ISIZE*(ISIZE+1)/2
C       KXA=ISIZE*ISIZE
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGZ(ISIZE,ISIZE,W(LXA),ISIZE,W(LXA),ICI)
C**TEMPORARY - HERMITIAN
        IEND=ISIZE
      END IF
      ISKIP=0
6666  CONTINUE
      IF(LANCZ)THEN
        IF(ICOUPX.GE.0)THEN
          REWIND 21
        END IF
        IF(ICOUPX.GT.1)THEN
          REWIND 22
        END IF
        IF(ICOUPX.GT.2)THEN
          REWIND 23
        END IF
        IF(ICOUPX.GT.3)THEN
          REWIND 24
        END IF
        IF(ICOUPX.GT.4)THEN
          REWIND 25
        END IF
        IF(ICOUPX.GT.5)THEN
          REWIND 26
        END IF
        IF(ISKIP.NE.0)THEN
C**RE-POSITION DISCS
          DO KKK=1,KROT-KSTEP,KSTEP
C************************************************************
C*******************************************TO BE COMPLETED??
C************************************************************
            IF(IREACT.NE.0)THEN
              CALL MEMO(1,LXA0,KXA1,0,0,0,0,0,0,0,0)
              CALL MATIN(W(LXA0),W(LXA0),NSIZE1,21)
              CALL MEMO(-1,LXA0,KXA1,0,0,0,0,0,0,0,0)
            END IF
            DO K=1,NAMODE
              IF(ICOUPX.GT.0)THEN
                CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
                CALL MATIN(W(LXA1),W(LXA1),NSIZE1,21)
                CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
              END IF
              DO L=1,K-1
                IF(ICOUPX.GT.1)THEN
                  CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
                  CALL MATIN(W(LXA2),W(LXA2),NSIZE2,22)
                  CALL MEMO(-1,LXA2,KXA2,0,0,0,0,0,0,0,0)
                END IF
                DO N=1,L-1
                  IF(ICOUPX.GT.2.AND.IREACT.EQ.0)THEN
                    CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
                    CALL MATIN(W(LXA3),W(LXA3),NSIZE3,23)
                    CALL MEMO(-1,LXA3,KXA3,0,0,0,0,0,0,0,0)
                  END IF
                  DO M=1,N-1
                    IF(ICOUPX.GT.3.AND.IREACT.EQ.0)THEN
                      CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                      CALL MATIN(W(LXA4),W(LXA4),NSIZE4,24)
                      CALL MEMO(-1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                    END IF
                    DO J=1,M-1
                      IF(ICOUPX.GT.4)THEN
                        CALL MEMO(1,LXA5,KXA5,0,0,0,0,0,0,0,0)
                        CALL MATIN(W(LXA5),W(LXA5),NSIZE5,25)
                        CALL MEMO(-1,LXA5,KXA5,0,0,0,0,0,0,0,0)
                      END IF
                      DO I=1,J-1
                        IF(ICOUPX.GT.5)THEN
                          CALL MEMO(1,LXA6,KXA6,0,0,0,0,0,0,0,0)
                          CALL MATIN(W(LXA6),W(LXA6),NSIZE6,26)
                          CALL MEMO(-1,LXA6,KXA6,0,0,0,0,0,0,0,0)
                        END IF
C************************************************************
C*******************************************TO BE COMPLETED??
C************************************************************
4006  CONTINUE
                      END DO
4005  CONTINUE
                    END DO
4004  CONTINUE
                  END DO
4003  CONTINUE
                END DO
4002  CONTINUE
              END DO
4001  CONTINUE
            END DO
          END DO
        END IF
        IF(LANIND.NE.2)THEN
          LASIZE=ISIZE
          IEND=ISTART-1
          KXA=0
          DO I=ISTART,LASIZE
            KXA=KXA+LSIZE
            IF(KXA.GT.KLAN)THEN
              KXA=KXA-LSIZE
              GO TO 6665
            END IF
            LSIZE=LSIZE+1
            IEND=IEND+1
          END DO
6665      CONTINUE
          IF(IEND.EQ.0)STOP 'LANMAX TOO SMALL'
          WRITE(IOUT,*)'ISTART,IEND,KXA',ISTART,IEND,KXA
          CALL FLUSH(IOUT)
          CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
          CALL DIAGL(KXA,W(LXA))
        END IF
      END IF
C**BOTH LINEAR AND RPH
      IF(IREACT.GT.0)THEN
        CALL INTARR(W(LNBF),W(LMBF),NSMODE,IK1TAU,IK2TAU,IK3TAU)
C**GET TORSION-ONLY INTEGRALS IF RPH
        IF(JREACT.GT.0)THEN
C****************************************  LIKE VCI1 (NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
          CALL MEMO(1,LIP1,KIP1,0,0,0,0,0,0,0,0)
          CALL GETBP1(W(LIP1),ISIZE1,NMODE,NSMODE,W(LMXBAS),NSIZE1,0)
C**READ INTO MATRIX
          KXA1=NSIZE1*(NSIZE1+1)/2
          CALL MEMO(1,LXA0,KXA1,0,0,0,0,0,0,0,0)
          CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
          CALL MATIN(W(LXA0),W(LXA0),NSIZE1,21)
          CALL VCV0(NNMODE,W(LXA),W(LXA0),
     1    NSIZE1,W(LIP),ISIZMX,ISIZE,W(LIP1),ISIZE1,ISTART,
     2    IEND,W(LWRK),W(LOV),LANIND,W(LXA),LAN12D,W(LYK),
     3    W(LTEMP),ICOUNT,ISIZEX)
          CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
          CALL MEMO(-1,LXA0,KXA1,0,0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VCI1 (NMODE)
        END IF
      END IF
C**COULD BE SINGLE (RPH) MODE - DO NOTHING MORE
      IF(NAMODE.EQ.0)GO TO 777
      IF(ITIM.EQ.0)ITIM1A=0
C**RESTORE INPUT VALUES
      ICOUPL=ICSSSL
      ICOUPC=ICSSSC
      IF(JREACT.LE.0.OR.NSMODE.EQ.NAMODE)THEN
        CALL VCI0(NAMODE,W(LXA),ISIZE,ISTART,IEND,LANIND,
     1  W(LXA),LAN12D,W(LYK),W(LTEMP),ICOUNT,ISIZEX,W(LEJK0),W(LEJK0),
     2  J21,KROT)
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K=0
      K2=0
      K3=0
      KRK=1
      DO KK=1,NAMODS
C       IF(KK.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        IF(LCOUNT.GT.1)THEN
          IF(KRK.LE.NAMODE)THEN
            KNEXT=JCONT(LCONT,KRK)
            IF(KK.EQ.KNEXT)KRK=KRK+1
          ELSE
            KNEXT=0
          END IF
          KKK=KRK-1
        ELSE
          KNEXT=KK
          KKK=KK
        END IF
C       IF(KNEXT-K.GT.1)THEN
C**UPDATE K2,K3 FOR MISSING MODES
C         DO KADJ=K+1,KNEXT-1
C           CALL INTARR(W(LNBF),W(LMBF),KADJ,IK1,IK2,IK3)
C           K2=K2+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
C           IF(KADJ.GT.NONLIN)K2=K2+IK1
C           K3=K3+IK2
C         END DO
C       END IF
        CALL INTARR(W(LNBF),W(LMBF),KK,IK1,IK2,IK3)
C**NEXT K
        K=KNEXT
        K1=K-1


C**DO WE WANT THIS KK?
        IF(KK.NE.KNEXT)THEN
C**DUMMY READ IF MISSING MODE KK
          IF((JREACT.LE.0.AND.KK.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     1    THEN
CC          IF(ICOUPL.GT.0.AND.IREACT.EQ.0)
CC   1      CALL DUMRD1(KK,IK2,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
CC   2      W(LEJK1),J21,W(LMODNT))
CC        ELSE
CC          IF(ICOUPL.GT.0)CALL DUMRV1(KK,IK2,IK2TAU,W(LV1),W(LV1),
CC   1      W(LC1),W(LC1),W(LEJK1),W(LEJK1),J21,W(LMODNT))
          END IF
C         GO TO 70
          GO TO 71
        END IF


C******************************************************************
C**POSSIBLY LINEAR (JREACT<0)
        IF(JREACT.LE.0.AND.ICOUPX.EQ.0)THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
          CALL MEMO(1,LIP1,KIP1,0,0,0,0,0,0,0,0)
          CALL GETBP1(W(LIP1),ISIZE1,NMODE,K,W(LMXBAS),NSIZE1,0)
          KXA1=NSIZE1*(NSIZE1+1)/2
C**READ INTO MATRIX
          CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
          CALL MATIN(W(LXA1),W(LXA1),NSIZE1,21)
C         CALL VCIKE(NNMODE,KK,W(LH+K2),W(LXQ+K3),W(LXA),W(LXA1),
          CALL VCIKE(NNMODE,KKK,W(LH+K2),W(LXQ+K3),W(LXA),W(LXA1),
     1    NSIZE1,IK3,IK2,W(LOMEGA+K1),W(LIP),ISIZMX,ISIZE,W(LIP1),
     2    ISIZE1,ISTART,IEND,W(LWRK),W(LOV),LANIND,
     3    W(LXA),LAN12D,W(LYK),W(LTEMP),ICOUNT,ISIZEX)
          CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
          CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
          GO TO 701
        END IF
C70    CONTINUE
C******************************************************************


C**INTRINSIC
C       IF(KCONT.LT.0)THEN
          CALL GROUP1(KK,ICONT(1),JCONT,1,IND)
          IF(IND.NE.0)THEN
            JCOUPL=NCOUPL(1)
            JCOUPC=NCOUPC(1)
            ICOUPL=IABS(JCOUPL)
            ICOUPC=IABS(JCOUPC)
            GO TO 7071
          END IF
          CALL GROUP1(KK,ICONT(2),JCONT,2,IND)
          IF(IND.NE.0)THEN
            JCOUPL=NCOUPL(2)
            JCOUPC=NCOUPC(2)
            ICOUPL=IABS(JCOUPL)
            ICOUPC=IABS(JCOUPC)
            GO TO 7071
          END IF
          JCOUPL=JCOUPS
          JCOUPC=JCOUCS
          ICOUPL=ICOUPS
          ICOUPC=ICOUCS
7071      CONTINUE
C       END IF
        IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.1)GO TO 71
CCCC    IF(ICOUPL.LT.1)GO TO 701


C**SKIP IF UNWANTED KK
        IF(KNEXT.EQ.0.OR.KK.NE.KNEXT)GO TO 71
C******************************************************************
C**CORIOLIS AND POTENTIAL
        IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
          CALL MEMO(1,LIP1,KIP1,0,0,0,0,0,0,0,0)
          CALL GETBP1(W(LIP1),ISIZE1,NMODE,K,W(LMXBAS),NSIZE1,0)
C**READ INTO MATRIX
          KXA1=NSIZE1*(NSIZE1+1)/2
          CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
          CALL MATIN(W(LXA1),W(LXA1),NSIZE1,21)
C         CALL VCI1(NNMODE,KK,W(LXA),W(LXA1),NSIZE1,W(LIP),ISIZMX,
          CALL VCI1(NNMODE,KKK,W(LXA),W(LXA1),NSIZE1,W(LIP),ISIZMX,
     1    ISIZE,W(LIP1),ISIZE1,ISTART,IEND,W(LWRK),W(LOV),LANIND,
     2    W(LXA),LAN12D,W(LYK),W(LTEMP),ICOUNT,ISIZEX,W(LXKAN),MAXQU,
     2    MAXPOW,J21,K,W(LNP1),W(LCP1),W(LMP1),NTOT1,IENTMX(1))
          CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
          CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
        ELSE
C****************************************  LIKE VCI2 (K+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
          IF(K.GT.NONLIN)THEN
            KPPP2=KLJP2
            IND2=1
          ELSE
            KPPP2=KIP2
            IND2=0
          END IF
          CALL MEMO(1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
          IF(K.GT.NONLIN)THEN
            CALL GETBP1(W(LIP2),LSIZE2,NMODE,K,W(LMXBAS),NSIZE2,IND2)
            MSIZE2=NSIZE2
          ELSE
            CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,NSMODE,W(LMXBAS),
     1      NSIZE2,IND2)
            MSIZE2=ISIZE2
          END IF
C**READ INTO MATRIX
          KXA2=NSIZE2*(NSIZE2+1)/2
          CALL MEMO(1,LXA1,KXA2,0,0,0,0,0,0,0,0)
          CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
          CALL MATIN(W(LXA1),W(LXA1),NSIZE2,21)
          CALL VCV1(NNMODE,KKK,W(LXA),W(LXA1),
     1    NSIZE2,W(LIP),ISIZMX,ISIZE,W(LIP2),MSIZE2,ISTART,
     2    IEND,W(LWRK),W(LOV),LANIND,W(LXA),LAN12D,W(LYK),
     3    W(LTEMP),ICOUNT,ISIZEX)
          CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
          CALL MEMO(-1,LXA1,KXA2,0,0,0,0,0,0,0,0)
          CALL MEMO(-1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VCI2 (K+NMODE)
        END IF
71    CONTINUE
C**INTRINSIC
        IF(ICOUPX.EQ.1)GO TO 701
C       IF(ICOUPL.EQ.1)GO TO 701
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L=0
        L2=0
        L3=0
        LRL=1
        DO LL=1,KK-1
          IF(LL.EQ.1.AND.KK.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          IF(LCOUNT.GT.1)THEN
            IF(LRL.LE.NAMODE)THEN
              LNEXT=JCONT(LCONT,LRL)
              IF(LL.EQ.LNEXT)LRL=LRL+1
            ELSE
              LNEXT=0
            END IF
            LLL=LRL-1
          ELSE
            LNEXT=LL
            LLL=LL
          END IF
C         IF(LNEXT-L.GT.1)THEN
C**UPDATE L3 FOR MISSING MODES
C           DO LADJ=L+1,LNEXT-1
C             CALL INTARR(W(LNBF),W(LMBF),LADJ,IL1,IL2,IL3)
C             L2=L2+IL1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
C             IF(LADJ.GT.NONLIN)L2=L2+IL1
C             L3=L3+IL2
C           END DO
C         END IF
          CALL INTARR(W(LNBF),W(LMBF),LL,IL1,IL2,IL3)
C**NEXT L
          L=LNEXT


C**DO WE WANT THIS KK AND THIS LL?
          IF(KK.NE.KNEXT.OR.LL.NE.LNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1      (NSMODE.EQ.NAMODE))THEN
CC            IF(COUPX.GT.1.AND.IREACT.EQ.0)
CC   1        CALL DUMRD2(KK,LL,IK2,IL2,W(LV2),W(LV2),
CC   2        W(LC2),W(LC2),W(LEJK2),W(LEJK2),J21,W(LMODNT))
CC          ELSE
CC            IF(ICOUPX.GT.1)
CC   1        CALL DUMRV2(KK,LL,IK2,IL2,IK2TAU,W(LV2),W(LV2),W(LC2),
CC   2        W(LC2),W(LEJK2),W(LEJK2),J21,W(LMODNT))
            END IF
            GO TO 72
          END IF


C**INTRINSIC
C         IF(KCONT.LT.0)THEN
            CALL GROUP2(KK,LL,ICONT(1),JCONT,1,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(1)
              JCOUPC=NCOUPC(1)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 7072
            END IF
            CALL GROUP2(KK,LL,ICONT(2),JCONT,2,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(2)
              JCOUPC=NCOUPC(2)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 7072
            END IF
            JCOUPL=JCOUPS
            JCOUPC=JCOUCS
            ICOUPL=ICOUPS
            ICOUPC=ICOUCS
7072        CONTINUE
C         END IF
          IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.2)GO TO 72
CCCC      IF(ICOUPL.LT.2)GO TO 72

C**SKIP IF UNWANTED LL
          IF(LNEXT.EQ.0.OR.LL.NE.LNEXT)GO TO 72
C******************************************************************
C**CORIOLIS AND POTENTIAL
          IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1    (NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
            CALL MEMO(1,LIP2,KIP2,0,0,0,0,0,0,0,0)
            CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,L,W(LMXBAS),NSIZE2,0)
C**READ INTO MATRIX
            KXA2=NSIZE2*(NSIZE2+1)/2
            CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
            CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
            CALL MATIN(W(LXA2),W(LXA2),NSIZE2,22)
C           CALL VCI2(NNMODE,KK,LL,W(LXA),W(LXA2),NSIZE2,W(LIP),ISIZMX,
            CALL VCI2(NNMODE,KKK,LLL,W(LXA),W(LXA2),NSIZE2,W(LIP),
     1      ISIZMX,ISIZE,W(LIP2),ISIZE2,ISTART,IEND,W(LWRK),W(LOV),
     2      LANIND,W(LXA),LAN12D,W(LYK),W(LTEMP),ICOUNT,ISIZEX,
     3      W(LXKAN),MAXQU,MAXPOW,J21,K,L,W(LNP2),W(LCP2),W(LMP2),
     4      NTOT2,IENTMX(2),W(LINDK))
            CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
            CALL MEMO(-1,LXA2,KXA2,0,0,0,0,0,0,0,0)
            CALL MEMO(-1,LIP2,KIP2,0,0,0,0,0,0,0,0)
          ELSE
C****************************************  LIKE VCI3 (K+L+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
            IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
              KPPP3=KLJP3
              IND3=1
            ELSE
              KPPP3=KIP3
              IND3=0
            END IF
            CALL MEMO(1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
            IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
              CALL GETBP2(W(LIP3),LSIZE3,NMODE,K,L,W(LMXBAS),NSIZE3,
     1        IND3)
              MSIZE3=NSIZE3
            ELSE
              CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,NSMODE,W(LMXBAS),
     1        NSIZE3,IND3)
              MSIZE3=ISIZE3
            END IF
C**READ INTO MATRIX
            KXA3=NSIZE3*(NSIZE3+1)/2
            CALL MEMO(1,LXA2,KXA3,0,0,0,0,0,0,0,0)
            CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
            CALL MATIN(W(LXA2),W(LXA2),NSIZE3,22)
            CALL VCV2(NNMODE,KKK,LLL,W(LXA),W(LXA2),
     1      NSIZE3,W(LIP),ISIZMX,ISIZE,W(LIP3),MSIZE3,ISTART,
     2      IEND,W(LWRK),W(LOV),LANIND,W(LXA),LAN12D,W(LYK),
     3      W(LTEMP),ICOUNT,ISIZEX)
            CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
            CALL MEMO(-1,LXA2,KXA3,0,0,0,0,0,0,0,0)
            CALL MEMO(-1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VCI3 (K+L+NMODE)
          END IF
72    CONTINUE
C**INTRINSIC
          IF(ICOUPX.EQ.2)GO TO 702
C         IF(ICOUPL.EQ.2)GO TO 702
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N=0
          N2=0
          N3=0
          NRN=1
          DO NN=1,LL-1
            IF(NN.EQ.1.AND.LL.EQ.2.AND.KK.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            IF(LCOUNT.GT.1)THEN
              IF(NRN.LE.NAMODE)THEN
                NNEXT=JCONT(LCONT,NRN)
                IF(NN.EQ.NNEXT)NRN=NRN+1
              ELSE
                NNEXT=0
              END IF
              NNN=NRN-1
            ELSE
              NNEXT=NN
              NNN=NN
            END IF
C           IF(NNEXT-N.GT.1)THEN
C**UPDATE N3 FOR MISSING MODES
C             DO NADJ=N+1,NNEXT-1
C               CALL INTARR(W(LNBF),W(LMBF),NADJ,IN1,IN2,IN3)
C**DUMMY READ IF MISSING MODES KK AND LL AND NN
C               IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
C    1          N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
CC                IF(IREACT.EQ.0)CALL DUMRD3(KK,LL,NN,IK2,IL2,IN2,
CC   1            W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,
CC   2            W(LMODNT))
C               ELSE
C                 CALL DUMRV3(KK,LL,NADJ,IK2,IL2,IN2,IK2TAU,W(LV3),
C    1            W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,W(LMODNT))
C               END IF
C               N2=N2+IN1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
C               IF(NADJ.GT.NONLIN)N2=N2+IN1
C               N3=N3+IN2
C             END DO
C           END IF
            CALL INTARR(W(LNBF),W(LMBF),NN,IN1,IN2,IN3)
C**NEXT N
            N=NNEXT


C**DO WE WANT THIS KK AND THIS LL AND THIS NN?
C           IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN
C             IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
C    1        N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
CC              IF(IREACT.EQ.0)CALL DUMRD3(KK,LL,NN,IK2,IL2,IN2,W(LV3),
CC   1          W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,W(LMODNT))
C             ELSE
C               CALL DUMRV3(KK,LL,NN,IK2,IL2,IN2,IK2TAU,W(LV3),W(LV3),
C    1          W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,W(LMODNT))
C             END IF
C             GO TO 73
C           END IF


C**INTRINSIC
C           IF(KCONT.LT.0)THEN
              CALL GROUP3(KK,LL,NN,ICONT(1),JCONT,1,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(1)
                JCOUPC=NCOUPC(1)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 7073
              END IF
              CALL GROUP3(KK,LL,NN,ICONT(2),JCONT,2,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(2)
                JCOUPC=NCOUPC(2)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 7073
              END IF
              JCOUPL=JCOUPS
              JCOUPC=JCOUCS
              ICOUPL=ICOUPS
              ICOUPC=ICOUCS
7073          CONTINUE
C           END IF
            IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.3)GO TO 73


C**DO WE WANT THIS KK AND THIS LL AND THIS NN?
            IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
     1      (KK.EQ.KNEXT.AND.LL.EQ.LNEXT.AND.NN.EQ.NNEXT.AND.
     2      ICOUPL.LT.3))THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
CC              IF(IREACT.EQ.0)CALL DUMRD3(KK,LL,NN,IK2,IL2,IN2,W(LV3),
CC   1          W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,W(LMODNT))
              ELSE
                IF(ICOUPX.GT.2)
     1          CALL DUMRV3(KK,LL,NN,IK2,IL2,IN2,IK2TAU,W(LV3),W(LV3),
     2          W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,W(LMODNT))
              END IF
              GO TO 73
            END IF


C**SKIP IF UNWANTED NN
            IF(NNEXT.EQ.0.OR.NN.NE.NNEXT)GO TO 73
C******************************************************************
C**CORIOLIS AND POTENTIAL
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1      N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
              CALL MEMO(1,LIP3,KIP3,0,0,0,0,0,0,0,0)
              CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,N,W(LMXBAS),NSIZE3,
     2        0)
C**READ INTO MATRIX
              KXA3=NSIZE3*(NSIZE3+1)/2
              CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
              CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
              CALL MATIN(W(LXA3),W(LXA3),NSIZE3,23)
C             CALL VCI3(NNMODE,KK,LL,NN,W(LXA),W(LXA3),NSIZE3,W(LIP),
              CALL VCI3(NNMODE,KKK,LLL,NNN,W(LXA),W(LXA3),NSIZE3,
     1        W(LIP),ISIZMX,ISIZE,W(LIP3),ISIZE3,ISTART,IEND,W(LWRK),
     2        W(LOV),LANIND,W(LXA),LAN12D,W(LYK),W(LTEMP),ICOUNT,
     3        ISIZEX,W(LXKAN),MAXQU,MAXPOW,J21,K,L,N,W(LNP3),W(LCP3),
     4        W(LMP3),NTOT3,IENTMX(3),W(LINDK),W(LINDL))
              CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
              CALL MEMO(-1,LXA3,KXA3,0,0,0,0,0,0,0,0)
              CALL MEMO(-1,LIP3,KIP3,0,0,0,0,0,0,0,0)
            ELSE
C****************************************  LIKE VCI4 (K+L+N+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
              IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                KPPP4=KLJP4
                IND4=1
              ELSE
                KPPP4=KIP4
                IND4=0
              END IF
              CALL MEMO(1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
              IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                CALL GETBP3(W(LIP4),LSIZE4,NMODE,K,L,N,W(LMXBAS),
     1          NSIZE4,IND4)
                MSIZE4=NSIZE4
              ELSE
                CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,NSMODE,
     1          W(LMXBAS),NSIZE4,IND4)
                MSIZE4=ISIZE4
              END IF
C**READ INTO MATRIX

!              KXA4=NSIZE4*(NSIZE4+1)/2
              if (mod(nsize4,2) .eq. 0) then
                 kxa4 = nsize4 / 2
                 kxa4 = kxa4 * (nsize4 + 1)
              else
                 kxa4 = (nsize4 + 1) / 2
                 kxa4 = kxa4 * nsize4
              end if

              CALL MEMO(1,LXA3,KXA4,0,0,0,0,0,0,0,0)
              CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
              CALL MATIN(W(LXA3),W(LXA3),NSIZE4,23)
              CALL VCV3(NNMODE,KKK,LLL,NNN,W(LXA),W(LXA3),
     1        NSIZE4,W(LIP),ISIZMX,ISIZE,W(LIP4),MSIZE4,ISTART,
     2        IEND,W(LWRK),W(LOV),LANIND,W(LXA),LAN12D,W(LYK),
     3        W(LTEMP),ICOUNT,ISIZEX)
              CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
              CALL MEMO(-1,LXA3,KXA4,0,0,0,0,0,0,0,0)
              CALL MEMO(-1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VCI4 (K+L+N+NMODE)
            END IF
73    CONTINUE
C**INTRINSIC
            IF(ICOUPX.EQ.3)GO TO 703
C           IF(ICOUPL.EQ.3)GO TO 703
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M=0
            M2=0
            M3=0
            MRM=1
            DO MM=1,NN-1
              IF(MM.EQ.1.AND.NN.EQ.2.AND.LL.EQ.3.AND.KK.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              IF(LCOUNT.GT.1)THEN
                IF(MRM.LE.NAMODE)THEN
                  MNEXT=JCONT(LCONT,MRM)
                  IF(MM.EQ.MNEXT)MRM=MRM+1
                ELSE
                  MNEXT=0
                END IF
                MMM=MRM-1
              ELSE
                MNEXT=MM
                MMM=MM
              END IF
C             IF(MNEXT-M.GT.1)THEN
C**UPDATE M3 FOR MISSING MODES
C               DO MADJ=M+1,MNEXT-1
C                 CALL INTARR(W(LNBF),W(LMBF),MADJ,IM1,IM2,IM3)
C                 M2=M2+IM1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
C                 IF(MADJ.GT.NONLIN)M2=M2+IM1
C                 M3=M3+IM2
C               END DO
C             END IF
              CALL INTARR(W(LNBF),W(LMBF),MM,IM1,IM2,IM3)
C**NEXT M
              M=MNEXT

C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM?
C             IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
C    1        MM.NE.MNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM
C               IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
C    1          N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
C    2          THEN
CC                IF(IREACT.EQ.0)CALL DUMRD4(KK,LL,NN,MM,IK2,IL2,IN2,
CC   1            IM2,W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
CC   2            J21,W(LMODNT))
C               ELSE
C                 CALL DUMRV4(KK,LL,NN,MM,IK2,IL2,IN2,IM2,IK2TAU,
C    1            W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,
C    2            W(LMODNT))
C               END IF
C               GO TO 74
C             END IF

C**INTRINSIC
C             IF(KCONT.LT.0)THEN
                CALL GROUP4(KK,LL,NN,MM,ICONT(1),JCONT,1,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(1)
                  JCOUPC=NCOUPC(1)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 7074
                END IF
                CALL GROUP4(KK,LL,NN,MM,ICONT(2),JCONT,2,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(2)
                  JCOUPC=NCOUPC(2)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 7074
                END IF
                JCOUPL=JCOUPS
                JCOUPC=JCOUCS
                ICOUPL=ICOUPS
                ICOUPC=ICOUCS
7074            CONTINUE
C             END IF
              IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.4)GO TO 74

C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM?
              IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
     1        MM.NE.MNEXT.OR.(KK.EQ.KNEXT.AND.LL.EQ.LNEXT.AND.
     2        NN.EQ.NNEXT.AND.MM.EQ.MNEXT.AND.ICOUPL.LT.4))THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2          THEN
C                 IF(IREACT.EQ.0)CALL DUMRD4(KK,LL,NN,MM,IK2,IL2,IN2,
C    1            IM2,W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),
C    2            J21,W(LMODNT))
                ELSE
                  IF(ICOUPX.GT.3)
     1            CALL DUMRV4(KK,LL,NN,MM,IK2,IL2,IN2,IM2,IK2TAU,
     2            W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,
     3            W(LMODNT))
                END IF
                GO TO 74
              END IF

C**SKIP IF UNWANTED MM
              IF(MNEXT.EQ.0.OR.MM.NE.MNEXT)GO TO 74
C******************************************************************
C**CORIOLIS AND POTENTIAL
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
                CALL MEMO(1,LIP4,KIP4,0,0,0,0,0,0,0,0)
                CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,M,W(LMXBAS),
     2          NSIZE4,0)
C**READ INTO MATRIX

!                KXA4=NSIZE4*(NSIZE4+1)/2
                if (mod(nsize4,2) .eq. 0) then
                   kxa4 = nsize4 / 2
                   kxa4 = kxa4 * (nsize4 + 1)
                else
                   kxa4 = (nsize4 + 1) / 2
                   kxa4 = kxa4 * nsize4
                end if

                CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                CALL MATIN(W(LXA4),W(LXA4),NSIZE4,24)
C               CALL VCI4(NNMODE,KK,LL,NN,MM,W(LXA),W(LXA4),NSIZE4,
                CALL VCI4(NNMODE,KKK,LLL,NNN,MMM,W(LXA),W(LXA4),NSIZE4, 
     1          W(LIP),ISIZMX,ISIZE,W(LIP4),ISIZE4,ISTART,IEND,W(LWRK),
     2          W(LOV),LANIND,W(LXA),LAN12D,W(LYK),W(LTEMP),ICOUNT,
     3          ISIZEX,W(LXKAN),MAXQU,MAXPOW,J21,K,L,N,M,W(LNP4),
     4          W(LCP4),W(LMP4),NTOT4,IENTMX(4),W(LINDK),W(LINDL),
     5          W(LINDN))
                CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                CALL MEMO(-1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                CALL MEMO(-1,LIP4,KIP4,0,0,0,0,0,0,0,0)
              ELSE
C****************************************  LIKE VCI5 (K+L+N+M+NMODE)
C**
                IF(IABS(IFITRP).LT.4)THEN
C**GET BASIC INTEGRAL BASIS FOR FIVE MODES
                  IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1              M.GT.NONLIN)THEN
                    KPPP5=KLJP5
                    IND5=1
                  ELSE
                    KPPP5=KIP5
                    IND5=0
                  END IF
                  CALL MEMO(1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
                  IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1              M.GT.NONLIN)THEN
                    CALL GETBP4(W(LIP5),LSIZE5,NSMODE,K,L,N,M,W(LMXBAS),
     1              NSIZE5,IND5)
                    MSIZE5=NSIZE5
                  ELSE
                    CALL GETBP5(W(LIP5),ISIZE5,NSMODE,K,L,N,M,NSMODE,
     1              W(LMXBAS),NSIZE5,IND5)
                    MSIZE5=ISIZE5
                  END IF
C**READ INTO MATRIX
                  KXA5=NSIZE5*(NSIZE5+1)/2
                  CALL MEMO(1,LXA4,KXA5,0,0,0,0,0,0,0,0)
                  CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                  CALL MATIN(W(LXA4),W(LXA4),NSIZE5,24)
                  CALL VCV4(NNMODE,KKK,LLL,NNN,MMM,W(LXA),W(LXA4),
     1            NSIZE5,W(LIP),ISIZMX,ISIZE,W(LIP5),MSIZE5,ISTART,
     2            IEND,W(LWRK),W(LOV),LANIND,W(LXA),LAN12D,W(LYK),
     3            W(LTEMP),ICOUNT,ISIZEX)
                  CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                  CALL MEMO(-1,LXA4,KXA5,0,0,0,0,0,0,0,0)
                  CALL MEMO(-1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
                ELSE
                  JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1,NTAU(4)+1,
     1            NTAU(5)+1)
                  KXP0=I36*JCIM*JCIM*IK2TAU
                  CALL MEMO(1,LXP0,KXP0,0,0,0,0,0,0,0,0)
C                 CALL VVV4(NNMODE,KK,LL,NN,MM,W(LXA),W(LIP),ISIZMX,
                  CALL VVV4(NNMODE,KKK,LLL,NNN,MMM,W(LXA),
     1            W(LIP),ISIZMX,
     1            ISIZE,ISTART,IEND,W(LWRK),ICOUNT,NL1MAX,NL2MAX,
     2            NL3MAX,NL4MAX,NR1MAX,NR2MAX,NR3MAX,NR4MAX,
     3            W(LH+K2TAU),W(LXQ+K3TAU),IK3TAU,IK2TAU,MDTAU,JCIM,
     4            W(LXP0),W(LMODNT),KROT,I36,0,W(LXKAN),MAXQU,MAXPOW,
     5            NAMODE,W(LCP4),W(LMP4),W(LNP4),IENTMX(4))
C**ON EXIT...
C**ICOUNT ARE NON-ZERO INTEGRALS
C**MDTAU ARE NUMBER RPH QUADRATURE POINTS SET ON FIRST LOOP WITH IND=0
                  CALL MEMO(3,LCP4,IENTMX(4)*MDTAU,LMP4,IENTMX(4)*4,
     1            LVP4,IENTMX(4)*JCIM*JCIM,0,0,0,0)
                  CALL MEMO(1,LWRK,3*ICOUNT,0,0,0,0,0,0,0,0)
C                 CALL VVV4(NNMODE,KK,LL,NN,MM,W(LXA),W(LIP),ISIZMX,
                  CALL VVV4(NNMODE,KKK,LLL,NNN,MMM,W(LXA),
     1            W(LIP),ISIZMX,
     1            ISIZE,ISTART,IEND,W(LWRK),ICOUNT,NL1MAX,NL2MAX,
     2            NL3MAX,NL4MAX,NR1MAX,NR2MAX,NR3MAX,NR4MAX,
     3            W(LH+K2TAU),W(LXQ+K3TAU),IK3TAU,IK2TAU,MDTAU,JCIM,
     4            W(LXP0),W(LMODNT),KROT,I36,1,W(LXKAN),MAXQU,MAXPOW,
     5            NAMODE,W(LCP4),W(LMP4),W(LVP4),IENTMX(4))
                  CALL MEMO(-1,LWRK,3*ICOUNT,0,0,0,0,0,0,0,0)
                  CALL MEMO(-3,LCP4,IENTMX(4)*MDTAU,LMP4,IENTMX(4)*4,
     1            LVP4,IENTMX(4)*JCIM*JCIM,0,0,0,0)
                  CALL MEMO(-1,LXP0,KXP0,0,0,0,0,0,0,0,0)
                END IF
C**
C****************************************  LIKE VCI5 (K+L+N+M+NMODE)
              END IF
74    CONTINUE
C**INTRINSIC
              IF(ICOUPX.EQ.4)GO TO 704
C             IF(ICOUPL.EQ.4)GO TO 704
C**INTEGRATE OVER FIVE NORMAL COORDINATES
              J=0
              J2=0
              J3=0
              JRJ=1
              DO JJ=1,MM-1
                IF(JJ.EQ.1.AND.MM.EQ.2.AND.NN.EQ.3.AND.LL.EQ.4.AND.
     1          KK.EQ.5.AND.ITIM.EQ.0)ITIM5A=0
                IF(LCOUNT.GT.1)THEN
                  IF(JRJ.LE.NAMODE)THEN
                    JNEXT=JCONT(LCONT,JRJ)
                    IF(JJ.EQ.JNEXT)JRJ=JRJ+1
                  ELSE
                    JNEXT=0
                  END IF
                  JJJ=JRJ-1
                ELSE
                  JNEXT=JJ
                  JJJ=JJ
                END IF
C               IF(JNEXT-J.GT.1)THEN
C**UPDATE J3 FOR MISSING MODES
C                 DO JADJ=J+1,JNEXT-1
C                   CALL INTARR(W(LNBF),W(LMBF),JADJ,IJ1,IJ2,IJ3)
C                   J2=J2+IJ1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
C                   IF(JADJ.GT.NONLIN)J2=J2+IJ1
C                   J3=J3+IJ2
C                 END DO
C               END IF
                CALL INTARR(W(LNBF),W(LMBF),JJ,IJ1,IJ2,IJ3)
C**NEXT J
                J=JNEXT

C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM AND THIS JJ?
C               IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
C    1          MM.NE.MNEXT.OR.JJ.NE.JNEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM AND JJ
C                 IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
C    1            N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN).OR.
C    2            (NSMODE.EQ.NAMODE))THEN
C                   IF(IREACT.EQ.0)CALL DUMRD5(KK,LL,NN,MM,JJ,IK2,IL2,
C    1              IN2,IM2,IJ2,W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),
C    2              W(LEJK5),J21,W(LMODNT))
C                 END IF
C                 GO TO 75
C               END IF

C**INTRINSIC
C               IF(KCONT.LT.0)THEN
                  CALL GROUP5(KK,LL,NN,MM,JJ,ICONT(1),JCONT,1,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(1)
                    JCOUPC=NCOUPC(1)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 7075
                  END IF
                  CALL GROUP5(KK,LL,NN,MM,JJ,ICONT(2),JCONT,2,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(2)
                    JCOUPC=NCOUPC(2)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 7075
                  END IF
                  JCOUPL=JCOUPS
                  JCOUPC=JCOUCS
                  ICOUPL=ICOUPS
                  ICOUPC=ICOUCS
7075              CONTINUE
C               END IF
                IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.5)GO TO 75


C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM AND THIS JJ?
                IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
     1          MM.NE.MNEXT.OR.JJ.NE.JNEXT.OR.(KK.EQ.KNEXT.AND.
     2          LL.EQ.LNEXT.AND.NN.EQ.NNEXT.AND.MM.EQ.MNEXT.AND.
     3          JJ.EQ.JNEXT.AND.ICOUPL.LT.5))THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM AND JJ
                  IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1            N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN).OR.
     2            (NSMODE.EQ.NAMODE))THEN
                    IF(ICOUPX.GT.4.AND.IREACT.EQ.0)CALL DUMRD5(KK,LL,
     1              NN,MM,JJ,IK2,IL2,
     2              IN2,IM2,IJ2,W(LV5),W(LV5),W(LC5),W(LC5),W(LEJK5),
     3              W(LEJK5),J21,W(LMODNT))
                  END IF
                  GO TO 75
                END IF


C**SKIP IF UNWANTED JJ
                IF(JNEXT.EQ.0.OR.JJ.NE.JNEXT)GO TO 75
C******************************************************************
C**CORIOLIS AND POTENTIAL
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN).OR.
     2          (NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR FIVE MODES
                  CALL MEMO(1,LIP5,KIP5,0,0,0,0,0,0,0,0)
                  CALL GETBP5(W(LIP5),ISIZE5,NMODE,K,L,N,M,J,
     1            W(LMXBAS),NSIZE5,0)
C**READ INTO MATRIX
                  KXA5=NSIZE5*(NSIZE5+1)/2
                  CALL MEMO(1,LXA5,KXA5,0,0,0,0,0,0,0,0)
                  CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                  CALL MATIN(W(LXA5),W(LXA5),NSIZE5,25)
C                 CALL VCI5(NNMODE,KK,LL,NN,MM,JJ,W(LXA),W(LXA5),
                  CALL VCI5(NNMODE,KKK,LLL,NNN,MMM,JJJ,W(LXA),W(LXA5),
     1            NSIZE5,W(LIP),ISIZMX,ISIZE,W(LIP5),ISIZE5,ISTART,
     2            IEND,W(LWRK),W(LOV),LANIND,W(LXA),LAN12D,W(LYK),
     3            W(LTEMP),ICOUNT,ISIZEX,W(LXKAN),MAXQU,MAXPOW,J21,K,L,
     4            N,M,J,W(LNP5),W(LCP5),W(LMP5),NTOT5,IENTMX(5),
     5            W(LINDK),W(LINDL),W(LINDN),W(LINDM))
                  CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                  CALL MEMO(-1,LXA5,KXA5,0,0,0,0,0,0,0,0)
                  CALL MEMO(-1,LIP5,KIP5,0,0,0,0,0,0,0,0)
                ELSE
C****************************************  LIKE VCI6 (K+L+N+M+J+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR SIX MODES
                  IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1              M.GT.NONLIN.OR.J.GT.NONLIN)THEN
                    KPPP6=KLJP6
                    IND6=1
                  ELSE
                    KPPP6=KIP6
                    IND6=0
                  END IF
                  CALL MEMO(1,LIP6,KPPP6,0,0,0,0,0,0,0,0)
                  IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1              M.GT.NONLIN.OR.J.GT.NONLIN)THEN
                    CALL GETBP5(W(LIP6),LSIZE6,NSMODE,K,L,N,M,J,
     1              W(LMXBAS),NSIZE6,IND6)
                    MSIZE6=NSIZE6
                  ELSE
                    CALL GETBP6(W(LIP6),ISIZE6,NSMODE,K,L,N,M,J,NSMODE,
     1              W(LMXBAS),NSIZE6,IND6)
                    MSIZE6=ISIZE6
                  END IF
C**READ INTO MATRIX
                  KXA6=NSIZE6*(NSIZE6+1)/2
                  CALL MEMO(1,LXA5,KXA6,0,0,0,0,0,0,0,0)
                  CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                  CALL MATIN(W(LXA5),W(LXA5),NSIZE6,25)
                  CALL VCV5(NNMODE,KKK,LLL,NNN,MMM,JJJ,W(LXA),W(LXA5),
     1            NSIZE6,W(LIP),ISIZMX,ISIZE,W(LIP6),MSIZE6,ISTART,
     2            IEND,W(LWRK),W(LOV),LANIND,W(LXA),LAN12D,W(LYK),
     3            W(LTEMP),ICOUNT,ISIZEX)
                  CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                  CALL MEMO(-1,LXA5,KXA6,0,0,0,0,0,0,0,0)
                  CALL MEMO(-1,LIP6,KPPP6,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VCI6 (K+L+N+M+J+NMODE)
                END IF
75    CONTINUE
C**INTRINSIC
                IF(ICOUPX.EQ.5)GO TO 705
C               IF(ICOUPL.EQ.5)GO TO 705
C**INTEGRATE OVER SIX NORMAL COORDINATES
                I=0
                I2=0
                I3=0
                IRI=1
                DO II=1,JJ-1
                  IF(II.EQ.1.AND.JJ.EQ.2.AND.MM.EQ.3.AND.NN.EQ.4.AND.
     1            LL.EQ.5.AND.KK.EQ.6.AND.ITIM.EQ.0)ITIM6A=0
                  IF(LCOUNT.GT.1)THEN
                    IF(IRI.LE.NAMODE)THEN
                      INEXT=JCONT(LCONT,IRI)
                      IF(II.EQ.INEXT)IRI=IRI+1
                    ELSE
                      INEXT=0
                    END IF
                    III=IRI-1
                  ELSE
                    INEXT=II
                    III=II
                  END IF
C                 IF(INEXT-I.GT.1)THEN
C**UPDATE I3 FOR MISSING MODES
C                   DO IADJ=I+1,INEXT-1
C                     CALL INTARR(W(LNBF),W(LMBF),IADJ,II1,II2,II3)
C                     I2=I2+II1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
C                     IF(IADJ.GT.NONLIN)I2=I2+II1
C                     I3=I3+II2
C                   END DO
C                 END IF
                  CALL INTARR(W(LNBF),W(LMBF),II,II1,II2,II3)
C**NEXT I
                  I=INEXT

C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM AND THIS JJ
C**AND THIS II?
C                 IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
C    1            MM.NE.MNEXT.OR.JJ.NE.JNEXT.OR.II.NE.INEXT)THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM AND JJ AND II
C                   IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.
C    1              AND.N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN.
C    2              AND.I.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C                     IF(IREACT.EQ.0)CALL DUMRD6(KK,LL,NN,MM,JJ,II,IK2,
C    1                IL2,IN2,IM2,IJ2,II2,W(LV6),W(LV6),W(LC6),W(LC6),
C    2                W(LEJK6),W(LEJK6),J21,W(LMODNT))
C                   END IF
C                   GO TO 76
C                 END IF

C**INTRINSIC
C                 IF(KCONT.LT.0)THEN
                    CALL GROUP6(K,L,N,M,J,I,ICONT(1),JCONT,1,IND)
                    IF(IND.NE.0)THEN
                      JCOUPL=NCOUPL(1)
                      JCOUPC=NCOUPC(1)
                      ICOUPL=IABS(JCOUPL)
                      ICOUPC=IABS(JCOUPC)
                      GO TO 7076
                    END IF
                    CALL GROUP6(K,L,N,M,J,I,ICONT(2),JCONT,2,IND)
                    IF(IND.NE.0)THEN
                      JCOUPL=NCOUPL(2)
                      JCOUPC=NCOUPC(2)
                      ICOUPL=IABS(JCOUPL)
                      ICOUPC=IABS(JCOUPC)
                      GO TO 7076
                    END IF
                    JCOUPL=JCOUPS
                    JCOUPC=JCOUCS
                    ICOUPL=ICOUPS
                    ICOUPC=ICOUCS
7076                CONTINUE
C                 END IF
                  IF((KCONT.LT.0.OR.MOLINC.GT.0).AND.ICOUPL.LT.6)
     1            GO TO 76


C**DO WE WANT THIS KK AND THIS LL AND THIS NN AND THIS MM AND THIS JJ
C**AND THIS II?
                  IF(KK.NE.KNEXT.OR.LL.NE.LNEXT.OR.NN.NE.NNEXT.OR.
     1            MM.NE.MNEXT.OR.JJ.NE.JNEXT.OR.II.NE.INEXT.OR.
     2            (KK.EQ.KNEXT.AND.LL.EQ.LNEXT.AND.NN.EQ.NNEXT.AND.
     3            MM.EQ.MNEXT.AND.JJ.EQ.JNEXT.AND.II.EQ.INEXT.AND.
     4            ICOUPL.LT.6))THEN
C**DUMMY READ IF MISSING MODES KK AND LL AND NN AND MM AND JJ AND II
                    IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.
     1              AND.N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN.
     2              AND.I.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
                      IF(ICOUPX.GT.5.AND.IREACT.EQ.0)
     1                CALL DUMRD6(KK,LL,NN,MM,JJ,II,IK2,
     2                IL2,IN2,IM2,IJ2,II2,W(LV6),W(LV6),W(LC6),W(LC6),
     3                W(LEJK6),W(LEJK6),J21,W(LMODNT))
                    END IF
                    GO TO 76
                  END IF


C**SKIP IF UNWANTED II
                  IF(INEXT.EQ.0.OR.II.NE.INEXT)GO TO 76
C******************************************************************
C**CORIOLIS AND POTENTIAL
                  IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1            N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN.AND.
     2            I.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR SIX MODES
                    CALL MEMO(1,LIP6,KIP6,0,0,0,0,0,0,0,0)
                    CALL GETBP6(W(LIP6),ISIZE6,NMODE,K,L,N,M,J,I,
     1              W(LMXBAS),NSIZE6,0)
C**READ INTO MATRIX
                    KXA6=NSIZE6*(NSIZE6+1)/2
                    CALL MEMO(1,LXA6,KXA6,0,0,0,0,0,0,0,0)
                    CALL MEMO(2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                    CALL MATIN(W(LXA6),W(LXA6),NSIZE6,26)
C                   CALL VCI6(NNMODE,KK,LL,NN,MM,JJ,II,W(LXA),
                    CALL VCI6(NNMODE,KKK,LLL,NNN,MMM,JJJ,III,W(LXA),
     1              W(LXA6),NSIZE6,W(LIP),ISIZMX,ISIZE,W(LIP6),ISIZE6,
     2              ISTART,IEND,W(LWRK),W(LOV),LANIND,W(LXA),LAN12D,
     3              W(LYK),W(LTEMP),ICOUNT,ISIZEX)
                    CALL MEMO(-2,LWRK,ISIZEX,LOV,ISIZEX,0,0,0,0,0,0)
                    CALL MEMO(-1,LXA6,KXA6,0,0,0,0,0,0,0,0)
                    CALL MEMO(-1,LIP6,KIP6,0,0,0,0,0,0,0,0)
                  END IF
76    CONTINUE
706   CONTINUE
                  I2=I2+II1
                  I3=I3+II2
                END DO
705   CONTINUE
                J2=J2+IJ1
                J3=J3+IJ2
              END DO
704   CONTINUE
              M2=M2+IM1
              M3=M3+IM2
            END DO
703   CONTINUE
            N2=N2+IN1
            N3=N3+IN2
          END DO
702   CONTINUE
          L2=L2+IL1
          L3=L3+IL2
        END DO
701   CONTINUE
        K2=K2+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
        IF(K.GT.NONLIN)K2=K2+IK1
        K3=K3+IK2
      END DO
777   CONTINUE
      ISKIP=0
      IF(LANCZ.AND.LANIND.NE.2)THEN
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        WRITE(IOUT,*)'Calculating LANOUT'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
        CALL LANOUT(W(LXA),W(LWK),W(LWK),W(LYK),W(LZK),W(LZK),ISIZE,
     1  ISTART,IEND,W(LANCNT),W(LANK))
        IF(LANIND.EQ.0.OR.LANIND.EQ.3)
     1  CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        IF(IEND.NE.LASIZE)THEN
          ISTART=IEND+1
          ISKIP=1
        END IF
      END IF
      IF(ISKIP.NE.0)GO TO 6666

C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
700   CONTINUE
C**TEMPORARY - HERMITIAN
C     CALL HERMTN(W(LXA),ISIZE,ISIZE,ISIZE)
C     STOP 'HERMITIAN'
C**TEMPORARY - HERMITIAN

      IF(LANIND.EQ.0)THEN
        KXK=ISIZE*NVAL
        KWK=ISIZE
        KSUP4=5*ISIZE
        KEVAL=NVAL
      END IF

      IF(MS.EQ.1.AND.LCOUNT.EQ.1)THEN
        EVL=EVLJ0
      ELSE
        IF(MS.EQ.1)EVL=0
      END IF

      IF(LANCZ)THEN

C**************************************************************
C**************************************************************
C**TEST TO SEE IF LANCUT SET...GET SIZE OF 1 & 2 MODE VCI BASIS
        IF(LANIND.EQ.1)THEN
C**REWIND 51 (TO BE RE-WRITTEN IN FULL LATER)
          REWIND 51
C*******************************************
          WRITE(IOUT,*)
          WRITE(IOUT,*)
          IF(TLLCUT.GT.0)
     1    WRITE(IOUT,*)'SYMMETRY = ',NS,' SIZE 1 & 2 DIM = ',LAN12D
          IF(TLLCUT.LT.0)
     1    WRITE(IOUT,*)'SYMMETRY = ',NS,' SIZE 1,2 & 3 DIM = ',LAN12D
          WRITE(IOUT,*)
          WRITE(IOUT,*)'*******************************************'
          IF(TLLCUT.GT.0)
     1    WRITE(IOUT,*)'1 & 2 DIMENSIONAL DIAGONALISATION'
          IF(TLLCUT.LT.0)
     1    WRITE(IOUT,*)'1,2 & 3 DIMENSIONAL DIAGONALISATION'
          WRITE(IOUT,*)'*******************************************'
          WRITE(IOUT,*)
          CALL FLUSH(IOUT)
C**DIAGONALISE 1 & 2 MODE PART OF FULL MATRIX, SIZE LAN12D...
C**...SAVE LANCUT EIGENVALUES & EIGENVECTORS
C**
          KXK=LAN12D
          LGIV=.TRUE.
          CALL MEMO(4,LEVAL,KXK,LSUP4,5*KXK,LWRK,
     1    KXK,LXK,KXK*LANCUT)
          ICID=-1
          IF(LANCUT.LE.KXK)THEN
            NVAL=LANCUT
            NVEC=LANCUT
          ELSE
            NVAL=KXK
            NVEC=KXK
          END IF
          IF(MS.EQ.1)EVL=0
          WRITE(IOUT,*)'Calculating DIAG'
          CALL TIMIT(1)
          CALL FLUSH(IOUT)
          CALL DIAG(W(LXA),W(LXK),KXK,KXK,0,W(LSUP4),W(LEVAL),W(LWRK),
     1    NVAL,NVEC,W(LIP),ISIZMX,NNMODE,W(LXK),W(LXK),W(LASSIG),
     2    ISIZMX,KEL21,NS)
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          IF(EVLCUT.EQ.0)EVLCUT=EVL
          CALL MEMO(-3,LSUP4,5*KXK,LXA,KXA,LWRK,KXK,0,0,0,0,0)
          ISKIP=0
          ISIZE=(ISIZEX-LAN12D)/LANINC
          ISIZEZ=ISIZE*LANINC
          CALL MEMO(1,LXA,LAN12D*ISIZE,0,0,0,0,0,0,0,0)
          CALL DIAGZ(LAN12D,ISIZE,W(LXA),0,W(LXA),0)
          ISTART=LAN12D+1
          IEND=LAN12D+ISIZE
          LANIND=2
          REWIND 52
          ICOUNT=0
          ICOUNTX=0
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          GO TO 6666
        END IF

        IF(LANIND.EQ.2)THEN
          CALL MEMO(1,LOV,LANCUT,0,0,0,0,0,0,0,0)
          IF(ICOUNTX.EQ.0)CALL MEMO(1,LTEMP1,ISIZEX,0,0,0,0,0,0,0,0)
          CALL TIMIT(3)
          CALL FLUSH(IOUT)
          CALL LANPT(W(LXK),W(LEVAL),LAN12D,W(LTEMP1),W(LXA),
     1    W(LYK),ISIZEX,W(LOV),W(LTEMP),ICOUNT,ICOUNTX,W(LJP),JSIZE,
     2    W(LIP),ISIZMX,NS,NNMODE,LANCUT,ISTART,IEND,0)
          IF(ICOUNT.EQ.ICOUNTX)GO TO 6663
          CALL MEMO(1,LTEMP,LAN12D*(ICOUNT-ICOUNTX),0,0,0,0,0,0,0,0)
          CALL LANPT(W(LXK),W(LEVAL),LAN12D,W(LTEMP1),W(LXA),
     1    W(LYK),ISIZEX,W(LOV),W(LTEMP),ICOUNT,ICOUNTX,W(LJP),JSIZE,
     2    W(LIP),ISIZMX,NS,NNMODE,LANCUT,ISTART,IEND,1)
          IF(LANINC.GT.1)
     1    CALL MEMO(-1,LTEMP,LAN12D*(ICOUNT-ICOUNTX),0,0,0,0,0,0,0,0)
          ICOUNTX=ICOUNT

6663      CONTINUE
          CALL MEMO(-1,LOV,LANCUT,0,0,0,0,0,0,0,0)
          IF(IEND.EQ.ISIZEX)GO TO 6664
          ISTART=IEND+1
          IEND=IEND+ISIZE
          IF(IEND.EQ.ISIZEZ+LAN12D)THEN
            CALL MEMO(-1,LXA,LAN12D*ISIZE,0,0,0,0,0,0,0,0)
            ISIZE=ISIZEX-IEND+ISIZE
            CALL MEMO(1,LXA,LAN12D*ISIZE,0,0,0,0,0,0,0,0)
            IEND=ISIZEX
          END IF
          CALL DIAGZ(LAN12D,ISIZE,W(LXA),0,W(LXA),0)
          GO TO 6666

6664      CONTINUE
          IF(LANINC.GT.1.AND.ICOUNT.GT.0)
     1    CALL MEMO(1,LTEMP,LAN12D*ICOUNT,0,0,0,0,0,0,0,0)
          IF(ICOUNT.GT.0)
     1    CALL LANPT(W(LXK),W(LEVAL),LAN12D,W(LTEMP1),W(LXA),
     2    W(LYK),ISIZEX,W(LOV),W(LTEMP),ICOUNT,ICOUNTX,W(LJP),JSIZE,
     3    W(LIP),ISIZMX,NS,NNMODE,LANCUT,ISTART,IEND,2)
          CALL MEMO(-1,LXA,LAN12D*ISIZE,0,0,0,0,0,0,0,0)
          CALL MEMO(-2,LXK,KXK*LANCUT,LEVAL,KXK,0,0,0,0,0,0)
          ISKIP=0
          ISTART=LAN12D+1
          ISIZE=ICOUNT+LAN12D
          LANIND=3
          IF(ICOUNT.GT.0)GO TO 6666
        END IF

        IF(LANCZ)THEN
          CALL MEMO(-1,LZK,ISIZEX,0,0,0,0,0,0,0,0)
          IF(LANIND.EQ.3)THEN
            IF(ICOUNT.GT.0)CALL MEMO(-1,LTEMP,LAN12D*ICOUNT,
     1      0,0,0,0,0,0,0,0)
            CALL MEMO(-1,LTEMP1,ISIZEX,0,0,0,0,0,0,0,0)
          END IF
          LAN12D=0

          IF(MS.EQ.1.AND.LCOUNT.EQ.1)THEN
            EVL=EVLJ0
          ELSE
            IF(MS.EQ.1)EVL=0
          END IF

          WRITE(IOUT,*)
          WRITE(IOUT,*)
          WRITE(IOUT,*)'SYMMETRY = ',NS,' SIZE FULL DIM = ',ISIZE
          WRITE(IOUT,*)
          WRITE(IOUT,*)'*******************************************'
          WRITE(IOUT,*)'FULL-DIMENSIONAL DIAGONALISATION'
          WRITE(IOUT,*)'*******************************************'
          WRITE(IOUT,*)
          IF(LANIND.EQ.3)NVAL=NVALSV
          NVEC=NVAL
        END IF
C**************************************************************
C**************************************************************
        WRITE(IOUT,260)
        CALL FLUSH(IOUT)
        IF(MATSIZ.NE.0)GO TO 9000
C**************************************************************
C**************************************************************
        REWIND 60
        LANCYC=LANCCC
        IF(LANCCC.GT.NVAL)LANCYC=NVAL
C**TEMPORARY
        IF(LCOUNT.GT.1)LANCYC=NVAL
C**TEMPORARY
C
        LGIV=(IGIV.NE.0)
        KXKL=LANCYC*NCYCLE*LANCYC*NCYCLE
        IF(LGIV)KXKL=LANCYC*NCYCLE*LANCYC
        CALL MEMO(3,LXK,ISIZE*NVAL,LVK,ISIZE*LANCYC,
     1  LZK,ISIZE*LANCYC,0,0,0,0)
        CALL MEMO(4,LSK,ISIZE*LANCYC,LEN,NVAL,LEL,NVAL,LNV,NVAL,0,0)
        CALL MEMO(4,LWRK,LANCYC*NCYCLE,LSUP4,5*LANCYC*NCYCLE,LEVAL,
     1  LANCYC*NCYCLE,LXKL,KXKL,0,0)
        KXA=LANCYC*NCYCLE*(LANCYC*NCYCLE+1)/2
        IF(LGIV)CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
        CALL MEMO(1,LXKLC,KXA,0,0,0,0,0,0,0,0)
        WRITE(IOUT,*)'Calculating LANCZOS'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
        IF(LGIV)THEN
          CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1    ISIZE,NVAL,W(LIP),ISIZMX,NNMODE,W(LKP),W(LKP+KLP),W(LASSIG),
     2    ISIZMX,KEL21,NS,
     2    W(LEN),W(LEL),W(LNV),W(LXA),W(LXKLC),W(LXKL),W(LWRK),
     3    W(LSUP4),W(LEVAL),NCYCLE,NVSYM,W(LANCNT),W(LANK),LANCYC)
        ELSE
          CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1    ISIZE,NVAL,W(LIP),ISIZMX,NNMODE,W(LKP),W(LKP+KLP),W(LASSIG),
     2    ISIZMX,KEL21,NS,
     2    W(LEN),W(LEL),W(LNV),W(LXKL),W(LXKLC),W(LXKL),W(LWRK),
     3    W(LSUP4),W(LEVAL),NCYCLE,NVSYM,W(LANCNT),W(LANK),LANCYC)
        END IF
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
C       EVLJ0=EVL
        IF(JTHIS.NE.0.OR.LCOUNT.GT.1)THEN
C**SAVE VCI ENERGIES AND COEFFICIENTS FOR Ka>0
          CALL ENERCF(W(LCFS),ISIZMX,W(LEVCI),ISIZE,NVALCF,KEL21,KEL,
     1    NS,W(LXK),W(LXK),W(LEL))
        END IF
C**DUMP CI COEFFICIENTS TO (INP4)
        IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
          IF(JTHIS.EQ.0)THEN
            IDDDP=JDUMP(NS)
            IF(NVAL.LT.JDUMP(NS))IDDDP=NVAL
            CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,LDUMP,IDDDP,
     1      W(LEL),NMODE)
          ELSE
            CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,NVAL,NVAL,
     1      W(LEL),NMODE)
          END IF
        END IF
        CALL MEMO(-5,LXK,ISIZE*NVAL,LVK,ISIZE*LANCYC,
     1  LZK,ISIZE*LANCYC,LWK,ISIZEX,LYK,ISIZEX)
        CALL MEMO(-5,LSK,ISIZE*LANCYC,LEN,NVAL,LEL,NVAL,LNV,NVAL,
     1  LXKL,KXKL)
        CALL MEMO(-4,LXKLC,KXA,LWRK,LANCYC*NCYCLE,
     1  LSUP4,5*LANCYC*NCYCLE,LEVAL,LANCYC*NCYCLE,0,0)
        IF(LGIV)CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LANCNT,ISIZEX,LANK,ISIZEX,0,0,0,0,0,0)
        LGIV=.FALSE.
      ELSE
        LGIV=.TRUE.
        CALL MEMO(4,LXK,KXK,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL,0,0)
        WRITE(IOUT,*)'Calculating DIAG'
        WRITE(IOUT,285)ISIZE,CUT
        WRITE(IOUT,502)LCONT,NVAL
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
        CALL DIAG(W(LXA),W(LXK),ISIZE,ISIZE,1,W(LSUP4),W(LEVAL),W(LWK),
     1  NVAL,NVEC,W(LIP),ISIZMX,NNMODE,W(LXK),W(LXK),W(LASSIG),ISIZMX,
     2  KEL21,NS)
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        IF(JTHIS.NE.0.OR.LCOUNT.GT.1)THEN
C**SAVE VCI ENERGIES AND COEFFICIENTS FOR Ka>0
          CALL ENERCF(W(LCFS),ISIZMX,W(LEVCI),ISIZE,NVALCF,KEL21,KEL,
     1    NS,W(LXK),W(LXK),W(LWK))
        END IF
C**DUMP CI COEFFICIENTS TO (INP4)
        IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
          IF(JTHIS.EQ.0)THEN
            IDDDP=JDUMP(NS)
            IF(NVAL.LT.JDUMP(NS))IDDDP=NVAL
            CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,LDUMP,IDDDP,
     1      W(LWK),NMODE)
          ELSE
            CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,NVAL,NVAL,
     1      W(LWK),NMODE)
          END IF
        END IF
        CALL MEMO(-4,LXK,KXK,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL,0,0)
        LGIV=.FALSE.
        CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
      END IF
C**SAVE PRIMARY ASSIGNMENT
      CALL SVASS(W(LASSIG),ISIZMX,W(LCASS),KSCFCI,KEL21,NVSYM,LCONT,NS,
     1NVAL,KA,KC)
9000  CONTINUE
4401  CONTINUE
4400  CONTINUE
C**DELETE ANY LEFT-OVER ARRAYS FROM FIRST CONTRACTION SCHEME
      IF(LCOUNT.GT.1)THEN
        CALL MEMO(-2,LIP,KIP,LJP,KJP,0,0,0,0,0,0)
        CALL MEMO(-1,LASSIG,ISIZMX*KEL21*NVSYM*3,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LCFS,ISIZMX*NVALCF*KEL21*NVSYM,LEVCI,
     1  NVALCF*KEL21*NVSYM,0,0,0,0,0,0)
      END IF
C     ICOUPX=ICPXSV
C     ICOUCX=ICCXSV
      WRITE(IOUT,240)EVL
C**TEMPORARY
      IF(ICOUNT.EQ.2)LANCZ=LOGTMP
C**TEMPORARY
      RETURN
      END
C**************************************************************

C**TEMPORARY
      SUBROUTINE PRMATZ(LAN12D,ISIZE,XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XL(LAN12D,1)
      COMMON/FILASS/IOUT
      WRITE(IOUT,*)'PRMATZ'
      DO J=1,10
      WRITE(IOUT,*)'COLUMN ',J
      WRITE(IOUT,*)(XL(I,J),I=1,LAN12D)
      END DO
      RETURN
      END
C**TEMPORARY

      SUBROUTINE GROUP1(K,ICONT,JCONT,NCONT,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JCONT(2,100)
      IND=0
      DO I=1,ICONT
        IF(K.EQ.JCONT(NCONT,I))IND=1
      END DO
      RETURN
      END
C**************************************************************
      SUBROUTINE GROUP2(K,L,ICONT,JCONT,NCONT,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JCONT(2,100)
      IND1=0
      IND2=0
      IND=0
      DO I=1,ICONT
        IF(K.EQ.JCONT(NCONT,I))IND1=1
        IF(L.EQ.JCONT(NCONT,I))IND2=1
      END DO
      IF(IND1.NE.0.AND.IND2.NE.0)IND=1
      RETURN
      END
C**************************************************************
      SUBROUTINE GROUP3(K,L,N,ICONT,JCONT,NCONT,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JCONT(2,100)
      IND1=0
      IND2=0
      IND3=0
      IND=0
      DO I=1,ICONT
        IF(K.EQ.JCONT(NCONT,I))IND1=1
        IF(L.EQ.JCONT(NCONT,I))IND2=1
        IF(N.EQ.JCONT(NCONT,I))IND3=1
      END DO
      IF(IND1.NE.0.AND.IND2.NE.0.AND.IND3.NE.0)IND=1
      RETURN
      END
C**************************************************************
      SUBROUTINE GROUP4(K,L,N,M,ICONT,JCONT,NCONT,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JCONT(2,100)
      IND1=0
      IND2=0
      IND3=0
      IND4=0
      IND=0
      DO I=1,ICONT
        IF(K.EQ.JCONT(NCONT,I))IND1=1
        IF(L.EQ.JCONT(NCONT,I))IND2=1
        IF(N.EQ.JCONT(NCONT,I))IND3=1
        IF(M.EQ.JCONT(NCONT,I))IND4=1
      END DO
      IF(IND1.NE.0.AND.IND2.NE.0.AND.IND3.NE.0.AND.IND4.NE.0)IND=1
      RETURN
      END
C**************************************************************
      SUBROUTINE GROUP5(K,L,N,M,J,ICONT,JCONT,NCONT,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JCONT(2,100)
      IND1=0
      IND2=0
      IND3=0
      IND4=0
      IND5=0
      IND=0
      DO I=1,ICONT
        IF(K.EQ.JCONT(NCONT,I))IND1=1
        IF(L.EQ.JCONT(NCONT,I))IND2=1
        IF(N.EQ.JCONT(NCONT,I))IND3=1
        IF(M.EQ.JCONT(NCONT,I))IND4=1
        IF(J.EQ.JCONT(NCONT,I))IND5=1
      END DO
      IF(IND1.NE.0.AND.IND2.NE.0.AND.IND3.NE.0.AND.IND4.NE.0.
     1AND.IND5.NE.0)IND=1
      RETURN
      END
C**************************************************************
      SUBROUTINE GROUP6(K,L,N,M,J,I,ICONT,JCONT,NCONT,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JCONT(2,100)
      IND1=0
      IND2=0
      IND3=0
      IND4=0
      IND5=0
      IND6=0
      IND=0
      DO II=1,ICONT
        IF(K.EQ.JCONT(NCONT,II))IND1=1
        IF(L.EQ.JCONT(NCONT,II))IND2=1
        IF(N.EQ.JCONT(NCONT,II))IND3=1
        IF(M.EQ.JCONT(NCONT,II))IND4=1
        IF(J.EQ.JCONT(NCONT,II))IND5=1
        IF(I.EQ.JCONT(NCONT,II))IND6=1
      END DO
      IF(IND1.NE.0.AND.IND2.NE.0.AND.IND3.NE.0.AND.IND4.NE.0.
     1AND.IND5.NE.0.AND.IND6.NE.0)IND=1
      RETURN
      END
C**************************************************************
      SUBROUTINE VSCFCI(W,NSTAT,EVLJ0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      LOGICAL LGIV,LINEAR,LANCZ,LANZA,LANZB,TRIAT,ABINIT
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
     1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
     4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
     5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
     6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
     7LDP4A,LDP4B,LDP4C,LDP4D,LXP,LJP6,LANCNT,LANK,LXJ,LVP5,
     8LC5,LC6,LVP0,LTMPRD,LMVB,LNFC,LTEMP6,LTEMP7,LTEMP8,LTEMP9
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/HERM/IHERM
      COMMON/CHECK/MCHECK
      COMMON/SADDLE/JNORM
      COMMON/PATH/ISCFCI
      COMMON/CYCLE/ICYCLE
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/DEGEN/IDEGEN,NDEGEN
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/REACTL/JREACT,IFITRP
      COMMON/SINCOS/ICS
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FITTER/MFIT,MFIT1(4),XFIT1,MFIT2(4),XFIT2,MFIT3(4),XFIT3,
     1MFIT4(4),XFIT4,MFIT5(4),XFIT5
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TRIATO/TRIAT
      COMMON/LANTOL/TOLLAN
      COMMON/CYCLES/NCYCLE
C**TEMPORARY (DIMENSIONS)
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LAN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TOLCUT,EVLCUT,TOLMAT,LANCYC
      COMMON/GIVEN/LGIV,IGIV
      COMMON/MATSIZ/MATSIZ
      COMMON/TYPE/LINEAR
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
      COMMON/FACTOR/FACTOR(6),FACTS(6)
      COMMON/VCIMAX/NMAX
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/RPHROT/IROTV
      COMMON/ESTATE/IORDER
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/AXES/MX(3),MXDIP(5),MXROT(3)
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/TBASIS/NTBAS(6,2),NTAU(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CSAVES/NNMODS,NAMODS,NVMODS,ICOUPS,ICOUCS,JREACS,NREACS
      COMMON/NCREC/NREC(6),MAXBUF(6),NUNITR,NUNITW
      COMMON/CSIZES/ISIZM1,ISIZM2,NVAL1,NVAL2,ICSIZ1,ICSIZ2,
     1IPSIZ1,IPSIZ2
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/UNITEX/I75,I76,I85,I86
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
      COMMON/KLSZS/KLJP1,KLJP2,KLJP3,KLJP4,KLJP5,KLJP6
      COMMON/KJSZS/KJP1,KJP2,KJP3,KJP4,KJP5,KJP6
      COMMON/KPSZS/KIP1,KIP2,KIP3,KIP4,KIP5,KIP6,JCC1,JCC2,JCC3,JCC4,
     1JCC5,JCC6
      COMMON/IPSZS/KPPP1,KPPP2,KPPP3,KPPP4,KPPP5,KPPP6
      COMMON/SIZEJ/JSIZE1(2),JSIZE2(2),JSIZE3(2),JSIZE4(2),JSIZE5(2),
     1JSIZE6(2)
      COMMON/TOTALS/ITOT1(2),ITOT2(2),ITOT3(2),ITOT4(2),ITOT5(2),
     1ITOT6(2),ITOT
      COMMON/TOTK/KTOT(6,2)
      COMMON/MATRIX/NVALV,NVALR,KSTEP,KSIGN,NVALCF
C**************************************************************
260   FORMAT(//,1X,'FINAL VIBRATIONAL (K-DIAGONAL) CI ENERGIES',/)
285   FORMAT(//,1X,'SIZE OF VIBRATIONAL CI MATRIX IS ',I6,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
290   FORMAT(//,1X,'VIBRATIONAL SYMMETRY ',I2)
395   FORMAT(//,1X,'SIZES OF CI SYMMETRY BLOCKS: ',/,1X,10I10,/)
C**************************************************************
C**************************************************************
C**MODULE FOR VSCF-CI
C**************************************************************
C**************************************************************
      IF(ICI.GT.0.AND.ISCFCI.GT.0)THEN
        ISIZE=ICI
        JSIZE=ISIZE
        KIP=ISIZE*NVMODE
        KJP=KIP
        CALL MEMO(2,LIP,KIP,LJP,KJP,0,0,0,0,0,0)
C**GET CI BASIS HERE IN IP
        CALL MOVEIP(W(LISTAT),NSTAT,W(LIP),ISIZE,NVMODE)
        CALL GETIP(W(LIP),W(LJP),ISIZE,ISIZE,JSIZE,NVMODE,ICI,ICI,0)
        WRITE(IOUT,395)(NTOT(K),K=1,NVSYM)
        CALL FLUSH(IOUT)
        IF(MATSIZ.NE.0)THEN
          WRITE(6,*) 'VSCF-CI MATRIX SIZE'
          RETURN
        END IF
        CALL MEMO(-2,LIP,KIP,LISTAT,NSTAT*NMODE,0,0,0,0,0,0)
C**RESET NSTAT IF VSCF-CI
        NSTAT=0
        DO MS=1,MVSYM
          NS=MWSYM(MS)
          IF(NS.NE.0)NSTAT=NSTAT+NTOT(NS)
        END DO
        CALL MEMO(1,LISTAT,NSTAT*NVMODE,0,0,0,0,0,0,0,0)
C**ONLY DO THOSE VSCF STATES WANTED IN VSCF-CI
C**BUG BUG BUG BUG
C       DO MS=1,MVSYM
C         NS=MWSYM(MS)
C         IF(NS.NE.0.AND.NTOT(NS).NE.0)
C    1    CALL PUTJP(W(LJP),JSIZE,W(LISTAT),NSTAT,NVMODE,NS)
C       END DO
        CALL MOVEIP(W(LJP),JSIZE,W(LISTAT),NSTAT,NVMODE)
C**BUG BUG BUG BUG
      END IF
C**GET FUNCTIONS SPACE FOR TWO SCF STATES
      KHLR=0
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        KHLR=KHLR+I2*3
      END DO
      CALL MEMO(2,LHL,KHLR,LHR,KHLR,0,0,0,0,0,0)
      NTOTMX=0
      DO NS=1,NVSYM
        IF(NTOT(NS).GT.NTOTMX)NTOTMX=NTOT(NS)
      END DO
      ISIZMX=NTOTMX
      CALL MEMO(1,LASSIG,ISIZMX*KEL21*NVSYM*3,0,0,0,0,0,0,0,0)
C     CALL MEMO(-1,LIP,KIP,0,0,0,0,0,0,0,0)
      KS=NVMODE
      CALL MEMO(1,LS,KS,0,0,0,0,0,0,0,0)
C**LOOP OVER SYMMETRIES
      ITIM=-1
      DO 6000 MS=1,MVSYM
C**********************************************************************
C**                                               LOOP ROUND SYMMETRIES
C**********************************************************************
      NS=MWSYM(MS)
      WRITE(IOUT,290)NS
      IF(NS.EQ.0)GO TO 6005
      ISIZE=NTOT(NS)
      WRITE(IOUT,285)ISIZE,CUT
      IF(ISIZE.EQ.0)GO TO 6005
      KIP=ISIZE*NSMODE
      KOV=ISIZE*ISIZE
      CALL MEMO(2,LIP,KIP,LOV,KOV,0,0,0,0,0,0)
      CALL PUTJP(W(LJP),JSIZE,W(LIP),ISIZE,NSMODE,NS)
      KCORIG=JTHIS
      DO 9500 KROT=1,J21
C**********************************************************************
C**                                                       LOOP ROUND Ka
C**********************************************************************
C**KA,KC FOR VSCF-CI
      KA=KROT/2
      KC=KCORIG-KA
C     IF(MOD(KA,2).NE.0)THEN
C       KC=KC+MOD(KROT,2)
C     ELSE
C       KC=KC+MOD(KROT+1,2)
C     END IF
C**RPH ODD K
      IF(KA.GT.0.AND.MOD(KA,2).EQ.MOD(KROT,2))KC=KC+1
C**RPH ODD K
C**ZEROISE OVERLAP MATRIX
      CALL DIAGZ(ISIZE,ISIZE,W(LOV),ISIZE,W(LOV),ICI)
C**SQUARE MATRIX NEEDED
      KXK=ISIZE*ISIZE
      CALL MEMO(1,LXK,KXK,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
      CALL DIAGZ(ISIZE,ISIZE,W(LXK),ISIZE,W(LXK),ICI)
C********************************************
C**LOOP OVER LHS VSCF STATES TO BE CALCULATED
C********************************************
      DO 2500 ISTATE=1,ISIZE
      REWIND 29
      INDEX=0
      DO IDUMMY=1,ICI
        IF(INDEX.EQ.0)THEN
C**READ LHS FUNCTIONS FOR THIS STATE
          K2=0
          DO K=1,NSMODE
            CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
            CALL IN29(W(LHL+K2),I2)
            K2=K2+I2*3
          END DO
          CALL COMPARE(W(LISTAT),NSTAT,W(LIP),ISIZE,NSMODE,IDUMMY,
     1    ISTATE,INDEX)
        END IF
      END DO
C********************************************
C**LOOP OVER RHS VSCF STATES TO BE CALCULATED
C********************************************
      DO 2000 JSTATE=ISTATE,ISIZE
      IF(ISTATE.NE.JSTATE)THEN
        REWIND 29
        JNDEX=0
        DO JDUMMY=1,ICI
          IF(JNDEX.EQ.0)THEN
C**READ RHS FUNCTIONS FOR THIS STATE
            K2=0
            DO K=1,NSMODE
              CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
              CALL IN29(W(LHR+K2),I2)
              K2=K2+I2*3
            END DO
            CALL COMPARE(W(LISTAT),NSTAT,W(LIP),ISIZE,NSMODE,JDUMMY,
     1      JSTATE,JNDEX)
          END IF
        END DO
      ELSE
C**COPY LHS FUNCTION
        DO K=1,KHLR
          W(LHR-1+K)=W(LHL-1+K)
        END DO
      END IF
C**************************************
C**GET OVERLAPS FOR THIS PAIR OF STATES
C**************************************
      K2=0
      DO K=1,NSMODE
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        IF(K.NE.IREACT)THEN
          CALL OVERLP(W(LHL+K2),W(LHR+K2),I2,W(LS-1+K))
        ELSE
          CALL OVERVP(W(LHL+K2),W(LHR+K2),I2,W(LS-1+K))
        END IF
        K2=K2+I2*3
      END DO
C**SET UP OVERLAP MATRIX
      CALL SCFOV(W(LS),NVMODE,W(LOV),ISIZE,ISTATE,JSTATE)
      ITIM=ITIM+1
      IF(ICOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 71
      END IF
      IF(ICOUPC.GT.0)THEN
        IF(J21.GT.1)REWIND 61
        REWIND 81
      END IF
      IF(ICOUPL.GT.1)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 72
      END IF
      IF(ICOUPC.GT.1)THEN
        IF(J21.GT.1)REWIND 62
        REWIND 82
      END IF
      IF(ICOUPL.GT.2)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 73
      END IF
      IF(ICOUPC.GT.2)THEN
        IF(J21.GT.1)REWIND 63
        REWIND 83
      END IF
      IF(ICOUPL.GT.3)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 74
      END IF
      IF(ICOUPC.GT.3)THEN
        IF(J21.GT.1)REWIND 64
        REWIND 84
      END IF
      IF(ICOUPL.GT.4)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 75
      END IF
      IF(IREACT.NE.0)THEN
C**GET POINTERS FOR TAU FUNCTIONS
        K2TAU=0
        K3TAU=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          K2TAU=K2TAU+IK2*3
          K3TAU=K3TAU+IK2
        END DO
        CALL INTARR(W(LNBF),W(LMBF),NSMODE,IK1TAU,IK2TAU,IK3TAU)
        CALL CV0(NAMODE,W(LHL+K2TAU),W(LHR+K2TAU),W(LXQ+K3TAU),W(LS),
     1  W(LXK),ISIZE,IK2TAU,ISTATE,JSTATE,W(LC0),W(LC0),W(LEJK0),
     2  W(LEJK0),J21,KROT,W(LMODNT))
      END IF
      IF(IREACT.EQ.0)
     1CALL CI0(NMODE,W(LS),W(LXK),ISIZE,ISTATE,JSTATE,W(LEJK0),
     2W(LEJK0),J21,KROT)
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NAMODE
        IF(K.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(IREACT.EQ.0)THEN
C**'KINETIC ENERGY'
          CALL CIKE(NMODE,K,W(LHL+K2),W(LHR+K2),W(LXQ+K3),W(LS),W(LXK),
     1    ISIZE,IK2,W(LOMEGA+K1),ISTATE,JSTATE)
        END IF
        IF(ICOUPL.EQ.0)GO TO 601
C**CORIOLIS AND POTENTIAL
        IF(IREACT.EQ.0)THEN
          CALL CI1(NMODE,K,W(LHL+K2),W(LHR+K2),W(LXQ+K3),W(LS),W(LXK),
     1    ISIZE,IK2,ISTATE,JSTATE,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
     2    W(LEJK1),J21,KROT,W(LMODNT))
        ELSE
C**IF REACT, ALWAYS INTEGRATE OVER TAU
          CALL CV1(NAMODE,K,W(LHL+K2),W(LHR+K2),W(LXQ+K3),W(LHL+K2TAU),
     1    W(LHR+K2TAU),W(LXQ+K3TAU),W(LS),W(LXK),ISIZE,IK2,IK2TAU,
     2    ISTATE,JSTATE,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),W(LEJK1),
     3    J21,KROT,W(LMODNT))
        END IF
        IF(ICOUPL.EQ.1)GO TO 601
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**CORIOLIS AND POTENTIAL
          IF(IREACT.EQ.0)THEN
            CALL CI2(NMODE,K,L,W(LHL+K2),W(LHR+K2),W(LXQ+K3),W(LHL+L2),
     1      W(LHR+L2),W(LXQ+L3),IK2,IL2,W(LS),W(LXK),ISIZE,NATOM,
     2      W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LRR),W(LXX),
     3      W(LX0),W(LXL),W(LXM),ISTATE,JSTATE,NPOT,W(LIPOT),W(LJPOT),
     4      W(LCPOT),W(LV2),W(LV2),W(LC2),W(LC2),W(LEJK2),W(LEJK2),J21,
     5      KROT,W(LMODNT))
          ELSE
C**IF REACT, ALWAYS INTEGRATE OVER TAU
            CALL CV2(NAMODE,K,L,W(LHL+K2),W(LHR+K2),W(LXQ+K3),
     1      W(LHL+L2),W(LHR+L2),W(LXQ+L3),W(LHL+K2TAU),W(LHR+K2TAU),
     2      W(LXQ+K3TAU),IK2,IL2,IK2TAU,W(LS),W(LXK),ISIZE,ISTATE,
     3      JSTATE,W(LV2),W(LV2),W(LC2),W(LC2),W(LEJK2),W(LEJK2),J21,
     5      KROT,W(LMODNT))
          END IF
          IF(ICOUPL.EQ.2)GO TO 602
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**CORIOLIS AND POTENTIAL
            IF(IREACT.EQ.0)THEN
              CALL CI3(NMODE,K,L,N,W(LHL+K2),W(LHR+K2),W(LXQ+K3),
     1        W(LHL+L2),W(LHR+L2),W(LXQ+L3),W(LHL+N2),W(LHR+N2),
     2        W(LXQ+N3),IK2,IL2,IN2,W(LS),W(LXK),ISIZE,NATOM,W(LQQ),
     3        W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LRR),W(LXX),
     4        W(LX0),W(LXL),W(LXM),ISTATE,JSTATE,
     5        NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV3),W(LV3),W(LC3),
     6        W(LC3),W(LEJK3),W(LEJK3),J21,KROT,W(LMODNT))
            ELSE
C**IF REACT, ALWAYS INTEGRATE OVER TAU
              CALL CV3(NAMODE,K,L,N,W(LHL+K2),W(LHR+K2),W(LXQ+K3),
     1        W(LHL+L2),W(LHR+L2),W(LXQ+L3),W(LHL+N2),W(LHR+N2),
     2        W(LXQ+N3),W(LHL+K2TAU),W(LHR+K2TAU),W(LXQ+K3TAU),IK2,IL2,
     3        IN2,IK2TAU,W(LS),W(LXK),ISIZE,ISTATE,JSTATE,W(LV3),
     4        W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,KROT,
     5        W(LMODNT))
            END IF
            IF(ICOUPL.EQ.3)GO TO 603
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.K.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**CORIOLIS AND POTENTIAL
              IF(IREACT.EQ.0)THEN
                CALL CI4(NMODE,K,L,N,M,W(LHL+K2),W(LHR+K2),W(LXQ+K3),
     1          W(LHL+L2),W(LHR+L2),W(LXQ+L3),W(LHL+N2),W(LHR+N2),
     2          W(LXQ+N3),W(LHL+M2),W(LHR+M2),W(LXQ+M3),IK2,IL2,IN2,
     3          IM2,W(LS),W(LXK),ISIZE,NATOM,W(LQQ),
     4          W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LRR),W(LXX),
     5          W(LX0),W(LXL),W(LXM),ISTATE,JSTATE,
     6          NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV4),W(LV4),W(LC4),
     7          W(LC4),W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT))
              ELSE
C**IF REACT, ALWAYS INTEGRATE OVER TAU
                CALL CV4(NAMODE,K,L,N,M,W(LHL+K2),W(LHR+K2),W(LXQ+K3),
     1          W(LHL+L2),W(LHR+L2),W(LXQ+L3),W(LHL+N2),W(LHR+N2),
     2          W(LXQ+N3),W(LHL+M2),W(LHR+M2),W(LXQ+M3),W(LHL+K2TAU),
     3          W(LHR+K2TAU),W(LXQ+K3TAU),IK2,IL2,IN2,IM2,IK2TAU,W(LS),
     4          W(LXK),ISIZE,ISTATE,JSTATE,W(LV4),W(LV4),W(LC4),W(LC4),
     5          W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT))
              END IF
              IF(ICOUPL.EQ.4)GO TO 604
C**5-MODE COUPLING HERE IF NEEDED
604   CONTINUE
              M2=M2+IM2*3
              M3=M3+IM2
            END DO
603   CONTINUE
            N2=N2+IN2*3
            N3=N3+IN2
          END DO
602   CONTINUE
          L2=L2+IL2*3
          L3=L3+IL2
        END DO
601   CONTINUE
        K2=K2+IK2*3
        K3=K3+IK2
      END DO
2000  CONTINUE
2500  CONTINUE
C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
      NVAL=ISCFCI
      NVEC=ISCFCI
      IF(ISIZE.LT.ISCFCI)THEN
        NVAL=ISIZE
        NVEC=ISIZE
      END IF
      WRITE(IOUT,260)
      KXA=ISIZE*ISIZE
      KWK=ISIZE
      KSUP4=ISIZE
      KEVAL=ISIZE
CCCC  IF(NS.EQ.1.AND.J21.EQ.1)EVL=0
CCCC  IF(NS.EQ.1.AND.J21.GT.1)EVL=EVLJ0
      IF(MS.EQ.1)EVL=EVLJ0
      CALL MEMO(4,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL,LXA,KXA,0,0)
      WRITE(IOUT,*)'Calculating DIAG'
      CALL TIMIT(1)
      CALL FLUSH(IOUT)
      CALL DIAG(W(LXK),W(LXK),ISIZE,ISIZE,1,W(LSUP4),W(LEVAL),W(LWK),
     1NVAL,NVEC,W(LIP),ISIZE,NMODE,W(LOV),W(LXA),W(LASSIG),ISIZMX,
     2J21,NS)
      CALL TIMIT(3)
      CALL FLUSH(IOUT)
      CALL MEMO(-5,LXA,KXA,LXK,KXK,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL)
9500  CONTINUE
      CALL MEMO(-2,LIP,KIP,LOV,KOV,0,0,0,0,0,0)
6005  CONTINUE
6000  CONTINUE
      IF(JREACT.NE.0)THEN
        CALL MEMO(-2,LEJK0,KEJK0,LC0,KLC0,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1  CALL MEMO(-1,LV1,KLV1,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPC.GT.0)THEN
        CALL MEMO(-1,LEJK1,KEJK1,0,0,0,0,0,0,0,0)
        CALL MEMO(-1,LC1,KLC1,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.1)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1  CALL MEMO(-1,LV2,KLV2,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPC.GT.1)THEN
        CALL MEMO(-1,LEJK2,KEJK2,0,0,0,0,0,0,0,0)
        CALL MEMO(-1,LC2,KLC2,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.2)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1  CALL MEMO(-1,LV3,KLV3,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPC.GT.2)THEN
        CALL MEMO(-1,LEJK3,KEJK3,0,0,0,0,0,0,0,0)
        CALL MEMO(-1,LC3,KLC3,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.3)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)
     1  CALL MEMO(-1,LV4,KLV4,0,0,0,0,0,0,0,0)
      END IF
      IF(ICOUPC.GT.3)THEN
        CALL MEMO(-1,LEJK4,KEJK4,0,0,0,0,0,0,0,0)
        CALL MEMO(-1,LC4,KLC4,0,0,0,0,0,0,0,0)
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE COMBIN(W,KLP,KSCFCI,EVLJ0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CHSYM(8)
      CHARACTER*2 SYMBOL(100)
      LOGICAL LGIV,LINEAR,LANCZ,LANZA,LANZB,TRIAT,ABINIT
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF,
     1LMODNT,LC0,LVM0,LEJK0,LXA0,LTEMP1,LTEMP2,LTEMP3,LXP0,LXTANH,
     2LMEFF,LIP5,LJP5,LKP,LLP,LTEMP4,LCASS,LWKC1,LWKC2,LWKC3,
     3LJPL,LJPR,LXJ0,LXI0,LXA5,LXA6,LNP1,LCP1,LMP1,LNP2,
     4LCP2,LMP2,LINDK,LNP3,LCP3,LMP3,LINDL,LNP4,LCP4,LMP4,
     5LINDN,LNP5,LCP5,LMP5,LINDM,LTEMP5,LXKAN,LV5,LV6,LIP6,
     6LVP1,LDP1,LVP2,LDP2A,LDP2B,LVP3,LDP3A,LDP3B,LDP3C,LVP4,
     7LDP4A,LDP4B,LDP4C,LDP4D,LXP,LJP6,LANCNT,LANK,LXJ,LVP5,
     8LC5,LC6,LVP0,LTMPRD,LMVB,LNFC,LTEMP6,LTEMP7,LTEMP8,LTEMP9
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/HERM/IHERM
      COMMON/CHECK/MCHECK
      COMMON/SADDLE/JNORM
      COMMON/PATH/ISCFCI
      COMMON/CYCLE/ICYCLE
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/DEGEN/IDEGEN,NDEGEN
      COMMON/REACTN/IREACT,MMTAU,INIT,INCTAU
      COMMON/REACTL/JREACT,IFITRP
      COMMON/ENTER/IENTER,IENTMX(5),NTOT1,NTOT2,NTOT3,NTOT4,NTOT5
      COMMON/SINCOS/ICS
      COMMON/MOLPRO/MOLPRO,MOUT,MINP,SYMBOL,CHSYM,MOLINC
      COMMON/FITTER/MFIT,MFIT1(4),XFIT1,MFIT2(4),XFIT2,MFIT3(4),XFIT3,
     1MFIT4(4),XFIT4,MFIT5(4),XFIT5
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TRIATO/TRIAT
      COMMON/LANTOL/TOLLAN
      COMMON/CYCLES/NCYCLE
C**TEMPORARY (DIMENSIONS)
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/MAXLAN/LANMAX,LLAN20,INP20,LLNCUT,LAN12D,NTOTL1(10),
     1NTOTL2(10),NTOTL3(10),LANDUM,TOLCUT,EVLCUT,TOLMAT,LANCYC
      COMMON/GIVEN/LGIV,IGIV
      COMMON/MATSIZ/MATSIZ
      COMMON/TYPE/LINEAR
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
      COMMON/FACTOR/FACTOR(6),FACTS(6)
      COMMON/VCIMAX/NMAX
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/RPHROT/IROTV
      COMMON/ESTATE/IORDER
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/AXES/MX(3),MXDIP(5),MXROT(3)
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/TBASIS/NTBAS(6,2),NTAU(6)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NRSYM,NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/SYMMP/ISYMP(10,10),ISYMPG,ISYMD(4,2)
      COMMON/CONTDP/ICONDP,ISYMST,ISYMFN
      COMMON/CSAVES/NNMODS,NAMODS,NVMODS,ICOUPS,ICOUCS,JREACS,NREACS
      COMMON/XSAVES/ICOUPX,ICOUCX,JCOUPS,JCOUCS
      COMMON/NCREC/NREC(6),MAXBUF(6),NUNITR,NUNITW
      COMMON/CSIZES/ISIZM1,ISIZM2,NVAL1,NVAL2,ICSIZ1,ICSIZ2,
     1IPSIZ1,IPSIZ2
      COMMON/CVAL/NCVAL(2,10),NCSIZE(10)
      COMMON/CONTX/LCOUNT,ISIZC(2,10)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
      COMMON/CONTT/ICCONT(2)
      COMMON/CONTC/NONC1,NONC2
      COMMON/CONSCH/NC1,NC2,NCL1(5),NCR1(5),NCL2(5),NCR2(5),MC1(5),
     1MC2(5),JC1(5),JC2(5)
      COMMON/UNITNO/I61,I62,I63,I64,I71,I72,I73,I74,I81,I82,I83,I84,
     1I91,I92,I93,I94
      COMMON/UNITEX/I75,I76,I85,I86
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
      COMMON/KLSZS/KLJP1,KLJP2,KLJP3,KLJP4,KLJP5,KLJP6
      COMMON/KJSZS/KJP1,KJP2,KJP3,KJP4,KJP5,KJP6
      COMMON/KPSZS/KIP1,KIP2,KIP3,KIP4,KIP5,KIP6,JCC1,JCC2,JCC3,JCC4,
     1JCC5,JCC6
      COMMON/IPSZS/KPPP1,KPPP2,KPPP3,KPPP4,KPPP5,KPPP6
      COMMON/SIZEJ/JSIZE1(2),JSIZE2(2),JSIZE3(2),JSIZE4(2),JSIZE5(2),
     1JSIZE6(2)
      COMMON/TOTALS/ITOT1(2),ITOT2(2),ITOT3(2),ITOT4(2),ITOT5(2),
     1ITOT6(2),ITOT
      COMMON/TOTK/KTOT(6,2)
      COMMON/MATRIX/NVALV,NVALR,KSTEP,KSIGN,NVALCF
C**************************************************************
240   FORMAT(/,1X,'ZERO POINT ENERGY = ',F10.2,/)
260   FORMAT(//,1X,'FINAL VIBRATIONAL (K-DIAGONAL) CI ENERGIES',/)
285   FORMAT(//,1X,'SIZE OF VIBRATIONAL CI MATRIX IS ',I6,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
290   FORMAT(//,1X,'VIBRATIONAL SYMMETRY ',I2)
500   FORMAT(//,1X,60(1H*),/,1X,'COMBINE CONTRACTION SCHEMES 1 AND 2',
     1/,1X,60(1H*),//)
501   FORMAT(//,1X,'CONTRACTION SCHEME ',I1,/,
     11X,'ADJUSTED VALUES NVAL = ',10I6,/)
505   FORMAT(//,1X,'CONTRACTION SCHEME ',I1,/,
     11X,'ISIZE = ',10I6,/)
514   FORMAT(1X,'POTENTIAL COUPLING OF ',I4,' MODES',/,
     1       1X,'CORIOLIS COUPLING OF ',I4,' MODES',/)
C*********************************************
C*********************************************
C**            COMBINE TWO CONTRACTION SCHEMES
C*********************************************
C*********************************************
      WRITE(IOUT,500)
C     WRITE(IOUT,514)ICOUPL,ICOUPC
C**TEMPORARY
      WRITE(IOUT,*)'COMBIN'
      WRITE(IOUT,*)'ICOUPL,ICOUPC,NCOUPL(1),NCOUPC(1),NCOUPL(2),
     1NCOUPC(2)'
      WRITE(IOUT,*)ICOUPL,ICOUPC,NCOUPL(1),NCOUPC(1),NCOUPL(2),
     1NCOUPC(2)
C**TEMPORARY
C     ICOUPX=MAX0(ICOUPL,IABS(NCOUPL(1)),IABS(NCOUPL(2)))
C     ICOUCX=MAX0(ICOUPC,IABS(NCOUPC(1)),IABS(NCOUPC(2)))
C**TEMPORARY
      WRITE(IOUT,*)'JCOUPL,JCOUPC,ICOUPX,ICOUCX,ICOUPS,ICOUCS'
      WRITE(IOUT,*)JCOUPL,JCOUPC,ICOUPX,ICOUCX,ICOUPS,ICOUCS
      JCOUPS=JCOUPL
      JCOUCS=JCOUPC
      CALL FLUSH(IOUT)
C**TEMPORARY

C**TEMPORARY
      CALL MEMPRT
C**TEMPORARY
      IF(IREACT.NE.0)THEN
C**GET POINTERS FOR TAU FUNCTIONS
        K2TAU=0
        K3TAU=0
        DO K=1,NAMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          K2TAU=K2TAU+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(K.GT.NONLIN)K2TAU=K2TAU+IK1
          K3TAU=K3TAU+IK2
        END DO
      END IF
C**********************************************************************
C**********************************************************************
C**                               WRITE DISCS AGAIN FOR COMPLETE SYSTEM
C**********************************************************************
C**********************************************************************
CCCC  IF(ICOUPL.GT.0)THEN
CCCC    IF(J21.GT.1)REWIND 61
CCCC    REWIND 71
CCCC    REWIND 81
CCCC  END IF
CCCC  IF(ICOUPL.GT.1)THEN
CCCC    IF(J21.GT.1)REWIND 62
CCCC    REWIND 72
CCCC    REWIND 82
CCCC  END IF
CCCC  IF(ICOUPL.GT.2)THEN
CCCC    IF(J21.GT.1)REWIND 63
CCCC    REWIND 73
CCCC    REWIND 83
CCCC  END IF
CCCC  IF(ICOUPL.GT.3)THEN
CCCC    IF(J21.GT.1)REWIND 64
CCCC    REWIND 74
CCCC    REWIND 84
CCCC  END IF
      I61=61
      I62=62
      I63=63
      I64=64
      I71=71
      I72=72
      I73=73
      I74=74
      I75=75
      I76=76
      I81=81
      I82=82
      I83=83
      I84=84
      I85=85
      I86=86
      I91=91
      I92=92
      I93=93
      I94=94
CCCC  CALL DISCDP(W,KCONT,0)
C**********************************************************************
C**********************************************************************
      DO I=1,2
        WRITE(IOUT,501)I,(NCVAL(I,J),J=1,NVSYM)
      END DO
      DO I=1,2
        WRITE(IOUT,505)I,(ISIZC(I,J),J=1,NVSYM)
      END DO
      ISIZXX=MAX0(ISIZM1,ISIZM2)
      NVALCF=MAX0(NVAL1,NVAL2)
      CALL MEMO(1,LASSIG,ISIZXX*KEL21*NVSYM*3,0,0,0,0,0,0,0,0)
      CALL MEMO(2,LCFS,ISIZXX*NVALCF*KEL21*NVSYM*2,LEVCI,
     1NVALCF*KEL21*NVSYM*2,0,0,0,0,0,0)
C**RELOAD COEFFICIENTS
      REWIND 30
      CALL CFENER(W(LCFS),W(LEVCI),ISIZXX,KEL21,NVSYM,NVALCF,1)
      CALL CFENER(W(LCFS),W(LEVCI),ISIZXX,KEL21,NVSYM,NVALCF,2)
      DO ICPL=1,ICOUPX
        FACTOR(ICPL)=1.D0
        DO I=2,NAMODE
          FACTOR(ICPL)=FACTOR(ICPL)*I
        END DO
        DO I=1,ICPL
          FACTOR(ICPL)=FACTOR(ICPL)/I
        END DO
        IDENOM=NAMODE-ICPL
        DO I=1,IDENOM
          FACTOR(ICPL)=FACTOR(ICPL)/I
        END DO
        IF(JREACT.GT.0)FACTOR(ICPL)=FACTOR(ICPL)+1
      END DO
C***************************************************
C***************************************************
C**INITIALLY NO LANCZOS
C***************************************************
C***************************************************
C**BUILD MATRIX FROM TWO EARLIER CONTRACTION SCHEMES
C**FOR EACH OF SCHEME 1 SYMMETRIES, THERE ARE NVAL1 STORED FUNCTIONS
C**FOR EACH OF SCHEME 2 SYMMETRIES, THERE ARE NVAL2 STORED FUNCTIONS
C**SELECTIVE SYMMETRIES IF J=0
      MMSS=MVSYM
C**MUST DO ALL 'VIBRATIONAL' SYMMETRIES IF J>0
      IF(JTHIS.NE.0)MMSS=NVSYM
C**********************************************************************
C**                                                LOOP OVER SYMMETRIES
C**********************************************************************
      NVALV=ISCFCI
      WRITE(IOUT,*)
      WRITE(IOUT,*)'ROTATION USES NVALV = ',NVALV
      ITIM=-1
      ICONDP=0
      DO 5500 MS=1,MMSS
      IF(JTHIS.EQ.0)THEN
        NS=MWSYM(MS)
      ELSE
        NS=MS
      END IF
      WRITE(IOUT,290)NS
      IF(NS.EQ.0)GO TO 5501
      ITIM=ITIM+1
C**NVAL1,NVAL2 ARE SIZES OF TWO EARLIER CONTRACTION SCHEMES
C**JUST IN CASE SINGLE CONTRACTION SCHEME REQUESTED (BAD USE OF 'MM')
      IF(NVAL2.EQ.0)NVAL2=1
      IF(ICCONT(1).GT.0.AND.ICCONT(2).GT.0)JSIZE=NVAL1*NVAL2*NVSYM
      IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
        JSIZE=0
        IF(NS.EQ.1)THEN
          DO I=1,NVSYM
            IF(ICCONT(1).LT.0)THEN
              JSIZE=JSIZE+ISIZC(1,I)*NCVAL(2,I)
            ELSE
              JSIZE=JSIZE+NCVAL(1,I)*ISIZC(2,I)
            END IF
          END DO
        END IF
        IF(NS.EQ.2)THEN
          K=1
          DO I=1,NVSYM
            IF(ICCONT(1).LT.0)THEN
              JSIZE=JSIZE+ISIZC(1,I)*NCVAL(2,I+K)
            ELSE
              JSIZE=JSIZE+NCVAL(1,I)*ISIZC(2,I+K)
            END IF
            K=-K
          END DO
        END IF
        IF(NS.EQ.3)THEN
          K=1
          DO I=1,NVSYM
            IF(ICCONT(1).LT.0)THEN
              JSIZE=JSIZE+ISIZC(1,I)*NCVAL(2,I+2*K)
            ELSE
              JSIZE=JSIZE+NCVAL(1,I)*ISIZC(2,I+2*K)
            END IF
            KK=MOD(I,2)+1
            K=K*(-1)**KK
          END DO
        END IF
        IF(NS.EQ.4)THEN
          DO I=1,NVSYM
            KSR=I/5
            K=5+NVSYM*KSR-I
            IF(ICCONT(1).LT.0)THEN
              JSIZE=JSIZE+ISIZC(1,I)*NCVAL(2,K)
            ELSE
              JSIZE=JSIZE+NCVAL(1,I)*ISIZC(2,K)
            END IF
          END DO
        END IF
        IF(NS.EQ.5)THEN
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            K=I+4*LSR
            IF(ICCONT(1).LT.0)THEN
              JSIZE=JSIZE+ISIZC(1,I)*NCVAL(2,K)
            ELSE
              JSIZE=JSIZE+NCVAL(1,I)*ISIZC(2,K)
            END IF
          END DO
        END IF
        IF(NS.EQ.6)THEN
          K=1
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            IF(ICCONT(1).LT.0)THEN
              JSIZE=JSIZE+ISIZC(1,I)*NCVAL(2,I+K+4*LSR)
            ELSE
              JSIZE=JSIZE+NCVAL(1,I)*ISIZC(2,I+K+4*LSR)
            END IF
            K=-K
          END DO
        END IF
        IF(NS.EQ.7)THEN
          K=1
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            IF(ICCONT(1).LT.0)THEN
              JSIZE=JSIZE+ISIZC(1,I)*NCVAL(2,I+2*K+4*LSR)
            ELSE
              JSIZE=JSIZE+NCVAL(1,I)*ISIZC(2,I+2*K+4*LSR)
            END IF
            KK=MOD(I,2)+1
            K=K*(-1)**KK
          END DO
        END IF
        IF(NS.EQ.8)THEN
          DO I=1,NVSYM
            KSR=I/5
            LSR=(-1)**(KSR+2)
            K=5+NVSYM*KSR-I
            IF(ICCONT(1).LT.0)THEN
              JSIZE=JSIZE+ISIZC(1,I)*NCVAL(2,K+4*LSR)
            ELSE
              JSIZE=JSIZE+NCVAL(1,I)*ISIZC(2,K+4*LSR)
            END IF
          END DO
        END IF
      END IF
      KJP=JSIZE*2
C**GET MAXIMUM (2-DIM) BASIS
      CALL MEMO(1,LJP,KJP,0,0,0,0,0,0,0,0)

C**DETERMINE TOTAL SIZE OF SYMMETRY-ADAPTED MATRIX
      ISIZE=0
C**JUST IN CASE SINGLE CONTRACTION SCHEME REQUESTED (BAD USE OF 'MM')
      IF(ISIZC(2,1).EQ.0)ISIZE=NCVAL(1,1)
      IF(NS.EQ.1)THEN
        DO I=1,NVSYM
          IF(ISIZC(1,I)*ISIZC(2,I).NE.0)THEN
            IF(ICCONT(1).GE.0.AND.ICCONT(2).GE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),NCVAL(2,I),ISIZE)
            ELSE
              IF(ICCONT(1).LT.0)THEN
                NCSIZE(I)=ISIZC(1,I)*NCVAL(2,I)
                CALL GETCP0(W(LJP),JSIZE,ISIZC(1,I),NCVAL(2,I),ISIZE)
              ELSE
                NCSIZE(I)=NCVAL(1,I)*ISIZC(2,I)
                CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),ISIZC(2,I),ISIZE)
              END IF
            END IF
          END IF
        END DO
      END IF
      IF(NS.EQ.2)THEN
        K=1
        DO I=1,NVSYM
          IF(ISIZC(1,I)*ISIZC(2,I+K).NE.0)THEN
            IF(ICCONT(1).GE.0.AND.ICCONT(2).GE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),NCVAL(2,I+K),ISIZE)
            ELSE
              IF(ICCONT(1).LT.0)THEN
                NCSIZE(I)=ISIZC(1,I)*NCVAL(2,I+K)
                CALL GETCP0(W(LJP),JSIZE,ISIZC(1,I),NCVAL(2,I+K),ISIZE)
              ELSE
                NCSIZE(I)=NCVAL(1,I)*ISIZC(2,I+K)
                CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),ISIZC(2,I+K),ISIZE)
              END IF
            END IF
          END IF
          K=-K
        END DO
      END IF
      IF(NS.EQ.3)THEN
        K=1
        DO I=1,NVSYM
          IF(ISIZC(1,I)*ISIZC(2,I+2*K).NE.0)THEN
            IF(ICCONT(1).GE.0.AND.ICCONT(2).GE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+2*K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),NCVAL(2,I+2*K),ISIZE)
            ELSE
              IF(ICCONT(1).LT.0)THEN
                NCSIZE(I)=ISIZC(1,I)*NCVAL(2,I+2*K)
                CALL GETCP0(W(LJP),JSIZE,ISIZC(1,I),NCVAL(2,I+2*K),
     1          ISIZE)
              ELSE
                NCSIZE(I)=NCVAL(1,I)*ISIZC(2,I+2*K)
                CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),ISIZC(2,I+2*K),
     1          ISIZE)
              END IF
            END IF
          END IF
          KK=MOD(I,2)+1
          K=K*(-1)**KK
        END DO
      END IF
      IF(NS.EQ.4)THEN
        DO I=1,NVSYM
          KSR=I/5
          LSR=(-1)**(KSR+2)
C         K=NVSYM+1-I 
          K=5+NVSYM*KSR-I
          IF(ISIZC(1,I)*ISIZC(2,K).NE.0)THEN
            IF(ICCONT(1).GE.0.AND.ICCONT(2).GE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),NCVAL(2,K),ISIZE)
            ELSE
              IF(ICCONT(1).LT.0)THEN
                NCSIZE(I)=ISIZC(1,I)*NCVAL(2,K)
                CALL GETCP0(W(LJP),JSIZE,ISIZC(1,I),NCVAL(2,K),ISIZE)
              ELSE
                NCSIZE(I)=NCVAL(1,I)*ISIZC(2,K)
                CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),ISIZC(2,K),ISIZE)
              END IF
            END IF
          END IF
        END DO
      END IF
      IF(NS.EQ.5)THEN
        DO I=1,NVSYM
          KSR=I/5
          LSR=(-1)**(KSR+2)
          K=I+4*LSR
          IF(ISIZC(1,I)*ISIZC(2,K).NE.0)THEN
            IF(ICCONT(1).GE.0.AND.ICCONT(2).GE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,K)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),NCVAL(2,K),ISIZE)
            ELSE
              IF(ICCONT(1).LT.0)THEN
                NCSIZE(I)=ISIZC(1,I)*NCVAL(2,K)
                CALL GETCP0(W(LJP),JSIZE,ISIZC(1,I),NCVAL(2,K),ISIZE)
              ELSE
                NCSIZE(I)=NCVAL(1,I)*ISIZC(2,K)
                CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),ISIZC(2,K),ISIZE)
              END IF
            END IF
          END IF
        END DO
      END IF
      IF(NS.EQ.6)THEN
        K=1
        DO I=1,NVSYM
          KSR=I/5
          LSR=(-1)**(KSR+2)
          IF(ISIZC(1,I)*ISIZC(2,I+K+4*LSR).NE.0)THEN
            IF(ICCONT(1).GE.0.AND.ICCONT(2).GE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+K+4*LSR)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),NCVAL(2,I+K+4*LSR),
     1        ISIZE)
            ELSE
              IF(ICCONT(1).LT.0)THEN
                NCSIZE(I)=ISIZC(1,I)*NCVAL(2,I+K+4*LSR)
                CALL GETCP0(W(LJP),JSIZE,ISIZC(1,I),NCVAL(2,I+K+4*LSR),
     1          ISIZE)
              ELSE
                NCSIZE(I)=NCVAL(1,I)*ISIZC(2,I+K+4*LSR)
                CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),ISIZC(2,I+K+4*LSR),
     1          ISIZE)
              END IF
            END IF
          END IF
          K=-K
        END DO
      END IF
      IF(NS.EQ.7)THEN
        K=1
        DO I=1,NVSYM
          KSR=I/5
          LSR=(-1)**(KSR+2)
          IF(ISIZC(1,I)*ISIZC(2,I+2*K+4*LSR).NE.0)THEN
            IF(ICCONT(1).GE.0.AND.ICCONT(2).GE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,I+2*K+4*LSR)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),NCVAL(2,I+2*K+4*LSR),
     1        ISIZE)
            ELSE
              IF(ICCONT(1).LT.0)THEN
                NCSIZE(I)=ISIZC(1,I)*NCVAL(2,I+2*K+4*LSR)
                CALL GETCP0(W(LJP),JSIZE,ISIZC(1,I),
     1          NCVAL(2,I+2*K+4*LSR),ISIZE)
              ELSE
                NCSIZE(I)=NCVAL(1,I)*ISIZC(2,I+2*K+4*LSR)
                CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),
     1          ISIZC(2,I+2*K+4*LSR),ISIZE)
              END IF
            END IF
          END IF
          KK=MOD(I,2)+1
          K=K*(-1)**KK
        END DO
      END IF
      IF(NS.EQ.8)THEN
        DO I=1,NVSYM
          KSR=I/5
          LSR=(-1)**(KSR+2)
          K=5+NVSYM*KSR-I
          IF(ISIZC(1,I)*ISIZC(2,K+4*LSR).NE.0)THEN
            IF(ICCONT(1).GE.0.AND.ICCONT(2).GE.0)THEN
              NCSIZE(I)=NCVAL(1,I)*NCVAL(2,K+4*LSR)
C**UPDATE ISIZE AND GET BASIS FOR THIS BLOCK
              CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),NCVAL(2,K+4*LSR),
     1        ISIZE)
            ELSE
              IF(ICCONT(1).LT.0)THEN
                NCSIZE(I)=ISIZC(1,I)*NCVAL(2,K+4*LSR)
                CALL GETCP0(W(LJP),JSIZE,ISIZC(1,I),NCVAL(2,K+4*LSR),
     1          ISIZE)
              ELSE
                NCSIZE(I)=NCVAL(1,I)*ISIZC(2,K+4*LSR)
                CALL GETCP0(W(LJP),JSIZE,NCVAL(1,I),ISIZC(2,K+4*LSR),
     1          ISIZE)
              END IF
            END IF
          END IF
        END DO
      END IF

C**ISIZMX IS TOTAL SIZE BASIS
      ISIZMX=ISIZE
      WRITE(IOUT,285)ISIZE,CUT
      CALL FLUSH(IOUT)
      IF(ISIZE.EQ.0)GO TO 5502
C**GET ACTUAL (2-DIM) BASIS
      KIP=ISIZE*2
      CALL MEMO(1,LIP,KIP,0,0,0,0,0,0,0,0)
C**SORT INTO SYMMETRIES (COPY INTO IP) ??
      CALL GETCPS(W(LJP),JSIZE,W(LIP),ISIZE)
      CALL MEMO(-1,LJP,KJP,0,0,0,0,0,0,0,0)
      IF(ICOUPX.GE.0)REWIND 21
      IF(ICOUPX.GT.1)REWIND 22
      IF(ICOUPX.GT.2)REWIND 23
      IF(ICOUPX.GT.3)REWIND 24
      IF(ICOUPX.GT.4)REWIND 25
      IF(ICOUPX.GT.5)REWIND 26
      KCORIG=JTHIS
      KEL=0
C**********************************************************************
C**                                                       LOOP ROUND Ka
C**********************************************************************
      DO 8000 KROT=1,J21,KSTEP
      KEL=KEL+1
C**KA,KC FOR VCI
      KA=KROT/2
      KC=KCORIG-KA
      IF(KMAX.GT.0)THEN
        IF(MOD(KA,2).NE.0)THEN
          KC=KC+MOD(KROT,2)
        ELSE
C         KC=KC+MOD(KROT+1,2)
          IF(KA.NE.0)KC=KC+MOD(KROT+1,2)
        END IF
      ELSE
        IF(MOD(KROT,2).EQ.0)KC=KC+1
      END IF

C*********************************************************
C*********************************************************
C**ISTART, IEND ARE START AND END COLUMNS
      ISTART=1
C**TEMPORARY-LAN
      ISYMST=1
C**TEMPORARY-LAN
      IF(LANCZ)THEN
C**TEMPORARY-LAN
        ISYMFN=1
        LSIZE0=0
C**TEMPORARY-LAN
        CALL MEMO(3,LWK,ISIZE,LYK,ISIZE,LZK,ISIZE,0,0,0,0)
        CALL MEMO(2,LANCNT,ISIZE,LANK,ISIZE,0,0,0,0,0,0)
        KLAN=LANMAX*(LANMAX+1)/2
        LSIZE=1
      ELSE
C**TEMPORARY-LAN
        ISYMFN=NVSYM
C**TEMPORARY-LAN
        KXA=ISIZE*(ISIZE+1)/2
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGZ(ISIZE,ISIZE,W(LXA),ISIZE,W(LXA),ICI)
C**IEND IS SIZE OF SINGLE UNIT
        IEND=ISIZE
      END IF
      ISKIP=0
5555  CONTINUE
      IF(LANCZ)THEN
        IF(ISKIP.NE.0)THEN
          IF(ICOUPX.GE.0)REWIND 21
          IF(ICOUPX.GT.1)REWIND 22
          IF(ICOUPX.GT.2)REWIND 23
          IF(ICOUPX.GT.3)REWIND 24
          IF(ICOUPX.GT.4)REWIND 25
          IF(ICOUPX.GT.5)REWIND 26
C**RE-POSITION DISCS
          DO KKK=1,KROT-KSTEP,KSTEP
C**********************************************************
C*******************************************TO BE COMPLETED
C**********************************************************
            IF(IREACT.NE.0)THEN
              CALL MEMO(1,LXA0,KXA1,0,0,0,0,0,0,0,0)
              CALL MATIN(W(LXA0),W(LXA0),NSIZE1,21)
              CALL MEMO(-1,LXA0,KXA1,0,0,0,0,0,0,0,0)
            END IF
            DO K=1,NAMODE
              IF(ICOUPX.GT.0)THEN
                CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
                CALL MATIN(W(LXA1),W(LXA1),NSIZE1,21)
                CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
              END IF
              DO L=1,K-1
                IF(ICOUPX.GT.1)THEN
                  CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
                  CALL MATIN(W(LXA2),W(LXA2),NSIZE2,22)
                  CALL MEMO(-1,LXA2,KXA2,0,0,0,0,0,0,0,0)
                END IF
                DO N=1,L-1
                  IF(ICOUPX.GT.2)THEN
                    CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
                    CALL MATIN(W(LXA3),W(LXA3),NSIZE3,23)
                    CALL MEMO(-1,LXA3,KXA3,0,0,0,0,0,0,0,0)
                  END IF
                  DO M=1,N-1
                    IF(ICOUPX.GT.3.AND.ICOUPX.NE.4)THEN
                      CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                      CALL MATIN(W(LXA4),W(LXA4),NSIZE4,24)
                      CALL MEMO(-1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                    END IF
                    DO J=1,M-1
                      IF(ICOUPX.GT.4)THEN
                        CALL MEMO(1,LXA5,KXA5,0,0,0,0,0,0,0,0)
                        CALL MATIN(W(LXA5),W(LXA5),NSIZE5,25)
                        CALL MEMO(-1,LXA5,KXA5,0,0,0,0,0,0,0,0)
                      END IF
                      DO I=1,J-1
                        IF(ICOUPX.GT.5)THEN
                          CALL MEMO(1,LXA6,KXA6,0,0,0,0,0,0,0,0)
                          CALL MATIN(W(LXA6),W(LXA6),NSIZE6,26)
                          CALL MEMO(-1,LXA6,KXA6,0,0,0,0,0,0,0,0)
                        END IF
C**********************************************************
C*******************************************TO BE COMPLETED
C**********************************************************
                      END DO
5005  CONTINUE
                    END DO
5004  CONTINUE
                  END DO
5003  CONTINUE
                END DO
5002  CONTINUE
              END DO
5001  CONTINUE
            END DO
          END DO
        END IF
        IEND=ISTART-1
C**TEMPORARY-LAN
        ISYMST=ISYMFN
        ISYMFN=ISYMST-1
        KXA0=0
        KXA=0
        DO JSYM=ISYMST,NVSYM
          IF(NCSIZE(JSYM).NE.0)THEN
            LSIZE=1
            DO I=1,NCSIZE(JSYM)
              KXA=KXA+LSIZE0+LSIZE
              IF(KXA.GT.KLAN)THEN
                KXA=KXA0
                IEND=LSIZE0
                GO TO 5504
              END IF
              LSIZE=LSIZE+1
              IEND=IEND+1
            END DO
            LSIZE0=LSIZE0+NCSIZE(JSYM)
            KXA0=KXA
          END IF
          ISYMFN=JSYM
        END DO
C**TEMPORARY-LAN
C       KXA=0
C       DO I=ISTART,ISIZE
C         KXA=KXA+LSIZE
C         IF(KXA.GT.KLAN)THEN
C           KXA=KXA-LSIZE
C           GO TO 5504
C         END IF
C         LSIZE=LSIZE+1
C         IEND=IEND+1
C       END DO
5504    CONTINUE
C**TEMPORARY-LAN
        WRITE(IOUT,*)'START & END SYMMETRIES = ',ISYMST,ISYMFN
C**TEMPORARY-LAN
        WRITE(IOUT,*)'ISTART,IEND,KXA',ISTART,IEND,KXA
        IF(IEND.LT.ISTART)STOP 'LANMAX TOO SMALL'
        CALL FLUSH(IOUT)
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGL(KXA,W(LXA))
      END IF
C*********************************************************
C*********************************************************

      IF(ICOUPX.GE.0)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 71
      END IF
      IF(ICOUCX.GE.0)THEN
        IF(J21.GT.1)REWIND 61
        REWIND 81
      END IF
      IF(ICOUPX.GT.1)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 72
      END IF
      IF(ICOUCX.GT.1)THEN
        IF(J21.GT.1)REWIND 62
        REWIND 82
      END IF
      IF(ICOUPX.GT.2)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 73
      END IF
      IF(ICOUCX.GT.2)THEN
        IF(J21.GT.1)REWIND 63
        REWIND 83
      END IF
      IF(ICOUPX.GT.3)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 74
      END IF
      IF(ICOUCX.GT.3)THEN
        IF(J21.GT.1)REWIND 64
        REWIND 84
      END IF
      IF(ICOUPX.GT.4)THEN
        IF(IWHICH.GE.0.OR.MOLINC.LE.0)REWIND 75
      END IF
C**BOTH LINEAR AND RPH
      IF(IREACT.GT.0)THEN
        CALL INTARR(W(LNBF),W(LMBF),NNMODE,IK1TAU,IK2TAU,IK3TAU)
C**GET TORSION-ONLY INTEGRALS IF RPH
        IF(JREACT.GT.0)THEN
C**DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'TAU'
          DO KK=1,ICONT(1)
            IF(NSMODE.EQ.JCONT(1,KK))THEN
              ITCONT=1
              ITMODE=KK
              ITBAS=0
              KTXK=ICSIZ1
              KTIP=IPSIZ1
              ITVAL=NVAL1
            END IF
          END DO
          DO KK=1,ICONT(2)
            IF(NSMODE.EQ.JCONT(2,KK))THEN
              ITCONT=2
              ITMODE=KK
              ITBAS=KLP
              KTXK=ICSIZ2
              KTIP=IPSIZ2
              ITVAL=NVAL2
            END IF
          END DO
C**MATRIX FOR CENTRAL (DGEMM) INTEGRALS FOR THIS SCHEME (MAYBE)
          CALL MEMO(1,LXW,KTXK*KTXK,0,0,0,0,0,0,0,0)
C****************************************  LIKE VCCI1 (NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
          CALL MEMO(1,LIP1,KIP1,0,0,0,0,0,0,0,0)
          CALL GETBP1(W(LIP1),ISIZE1,NMODE,NNMODE,W(LMXBAS),NSIZE1,
     1    0)
          KXA1=NSIZE1*(NSIZE1+1)/2
          CALL MEMO(3,LXA0,KXA1,LTEMP,KTXK*ITVAL,LWRK,ITVAL*ITVAL,
     1    0,0,0,0)
C**ZEROISE MATRIX
          CALL DIAGZ(NSIZE1,NSIZE1,W(LXA0),NSIZE1,W(LXA0),ICI)
          CALL VCCV0(NAMODE,NNMODE,ITMODE,W(LH+K2TAU),W(LXQ+K3TAU),
     1    IK3TAU,IK2TAU,W(LIP),ISIZMX,W(LKP+ITBAS),KTXK,KTIP,ITCONT,
     2    W(LXA),ISIZE,W(LIP1),ISIZE1,W(LXA0),W(LXA0),NSIZE1,W(LXW),
     3    W(LTEMP),W(LWRK),ITVAL,ISTART,IEND,W(LC0),W(LC0),W(LEJK0),
     4    W(LEJK0),J21,KROT,W(LMODNT),KEL,NS,
     5    W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(ITCONT-1)),ISIZXX,NVALCF,
     6    KEL21,NVSYM)
          CALL MEMO(-3,LXA0,KXA1,LTEMP,KTXK*ITVAL,LWRK,ITVAL*ITVAL,
     1    0,0,0,0)
          CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
          CALL MEMO(-1,LXW,KTXK*KTXK,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VCCI1 (NMODE)
        END IF
      END IF
      IF(JREACT.LE.0.OR.NSMODE.EQ.NAMODE)THEN
C**TEMPORARY-LAN
C       CALL VCI0(NAMODE,W(LXA),ISIZE,1,ISIZE,W(LEJK0),W(LEJK0),J21,
C    1  KROT)
        CALL VCI0(NAMODE,W(LXA),ISIZE,ISTART,IEND,W(LEJK0),W(LEJK0),
     2  J21,KROT)
C**TEMPORARY-LAN
        IF(MOLINC.GT.0.AND.ICCONT(1).GE.0.AND.ICCONT(2).GE.0)
     1  CALL VCCI0(W(LIP),ISIZMX,W(LXA),KEL,NS,W(LEVCI),NVALCF,KEL21,
     1  NVSYM)
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NAMODE
        IF(K.EQ.1.AND.ITIM.EQ.0)THEN
          ITIM1A=0
          ITIM1B=0
        END IF
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'K'
        DO KK=1,ICONT(1)
          IF(K.EQ.JCONT(1,KK))THEN
            KKCONT=1
            KMODX=KK
            KONBAS=0
            KKXK=ICSIZ1
            KKIP=IPSIZ1
            KVAL=NVAL1
          END IF
        END DO
        DO KK=1,ICONT(2)
          IF(K.EQ.JCONT(2,KK))THEN
            KKCONT=2
            KMODX=KK
            KONBAS=KLP
            KKXK=ICSIZ2
            KKIP=IPSIZ2
            KVAL=NVAL2
          END IF
        END DO
C**MATRIX FOR CENTRAL (DGEMM) INTEGRALS FOR THIS SCHEME
C******************************************************************
C**POSSIBLY LINEAR (JREACT<0)
        IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
          IF(ICOUPL.EQ.0)THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
            NREC(1)=NREC(1)+1
            CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
            CALL MEMO(1,LIP1,KIP1,0,0,0,0,0,0,0,0)
            CALL GETBP1(W(LIP1),ISIZE1,NMODE,K,W(LMXBAS),NSIZE1,0)
            KXA1=NSIZE1*(NSIZE1+1)/2
            CALL MEMO(3,LXA1,KXA1,LTEMP,KKXK*KVAL,LWRK,KVAL*KVAL,
     1      0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE1,NSIZE1,W(LXA1),NSIZE1,W(LXA1),ICI)
            IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
              CALL MEMO(1,LOV,KKXK,0,0,0,0,0,0,0,0)
              CALL VCCIKC(NAMODE,NNMODE,K,KMODX,W(LH+K2),W(LXQ+K3),IK3,
     1        IK2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),
     2        ISIZE,W(LIP1),ISIZE1,W(LXA1),W(LXA1),NSIZE1,W(LXK),
     3        W(LTEMP),W(LOV),W(LWRK),KVAL,ISTART,IEND,W(LMODNT),
     4        W(LOMEGA+K1),KEL,NS,
     5        W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),ISIZXX,
     6        NVALCF,KEL21,NVSYM)
              CALL MEMO(-1,LOV,KKXK,0,0,0,0,0,0,0,0)
            ELSE
              IF(MOLINC.LE.0)
     1        CALL VCCIKE(NAMODE,NNMODE,K,KMODX,W(LH+K2),W(LXQ+K3),IK3,
     1        IK2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),
     2        ISIZE,W(LIP1),ISIZE1,W(LXA1),W(LXA1),NSIZE1,W(LXK),
     3        W(LTEMP),W(LWRK),KVAL,ISTART,IEND,W(LMODNT),W(LOMEGA+K1),
     4        KEL,NS,W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     5        ISIZXX,NVALCF,KEL21,NVSYM)
            END IF
            CALL MEMO(-3,LXA1,KXA1,LTEMP,KKXK*KVAL,LWRK,KVAL*KVAL,
     1      0,0,0,0)
            CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
            CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
          END IF
        END IF
C**INTRINSIC
        IF(MOLINC.GT.0)THEN
          CALL GROUP1(K,ICONT(1),JCONT,1,IND)
          IF(IND.NE.0)THEN
            JCOUPL=NCOUPL(1)
            JCOUPC=NCOUPC(1)
            ICOUPL=IABS(JCOUPL)
            ICOUPC=IABS(JCOUPC)
            GO TO 3031
          END IF
          CALL GROUP1(K,ICONT(2),JCONT,2,IND)
          IF(IND.NE.0)THEN
            JCOUPL=NCOUPL(2)
            JCOUPC=NCOUPC(2)
            ICOUPL=IABS(JCOUPL)
            ICOUPC=IABS(JCOUPC)
            GO TO 3031
          END IF
          JCOUPL=JCOUPS
          JCOUPC=JCOUCS
          ICOUPL=ICOUPS
          ICOUPC=ICOUCS
3031      CONTINUE
        END IF
        IF(ICOUPL.LT.1)GO TO 331
C**MATRIX FOR CENTRAL (DGEMM) INTEGRALS FOR THIS SCHEME
CCCC    CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
C******************************************************************
C**CORIOLIS AND POTENTIAL
        IF((JREACT.LE.0.AND.K.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR SINGLE MODE
          NREC(1)=NREC(1)+1
          CALL MEMO(1,LIP1,KIP1,0,0,0,0,0,0,0,0)
          CALL GETBP1(W(LIP1),ISIZE1,NMODE,K,W(LMXBAS),NSIZE1,0)
          KXA1=NSIZE1*(NSIZE1+1)/2
          CALL MEMO(3,LXA1,KXA1,LTEMP,KKXK*KVAL,LWRK,
     1    KVAL*KVAL,0,0,0,0)
          CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
          CALL DIAGZ(NSIZE1,NSIZE1,W(LXA1),NSIZE1,W(LXA1),ICI)
          IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
            CALL MEMO(1,LOV,KKXK,0,0,0,0,0,0,0,0)
            CALL VCCI1C(NAMODE,NNMODE,K,KMODX,W(LH+K2),W(LXQ+K3),IK3,
     1      IK2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),
     2      ISIZE,W(LIP1),ISIZE1,W(LXA1),W(LXA1),NSIZE1,W(LXK),
     3      W(LTEMP),W(LOV),W(LWRK),KVAL,ISTART,IEND,W(LV1),W(LV1),
     4      W(LC1),W(LC1),W(LEJK1),W(LEJK1),J21,KROT,W(LMODNT),
     5      W(LOMEGA+K1),KEL,NS,
     6      W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),ISIZXX,NVALCF,
     7      KEL21,NVSYM,W(LXKAN),MAXQU,MAXPOW,W(LNP1),W(LCP1),W(LMP1),
     8      NTOT1,IENTMX(1))
            CALL MEMO(-1,LOV,KKXK,0,0,0,0,0,0,0,0)
          ELSE
            CALL VCCI1(NAMODE,NNMODE,K,KMODX,W(LH+K2),W(LXQ+K3),IK3,
     1      IK2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),
     2      ISIZE,W(LIP1),ISIZE1,W(LXA1),W(LXA1),NSIZE1,W(LXK),
     3      W(LTEMP),W(LWRK),KVAL,ISTART,IEND,W(LV1),W(LV1),W(LC1),
     4      W(LC1),W(LEJK1),W(LEJK1),J21,KROT,W(LMODNT),W(LOMEGA+K1),
     5      KEL,NS,W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),ISIZXX,
     6      NVALCF,KEL21,NVSYM,W(LXKAN),MAXQU,MAXPOW,W(LNP1),W(LCP1),
     7      W(LMP1),NTOT1,IENTMX(1))
          END IF
          CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
          CALL MEMO(-3,LXA1,KXA1,LTEMP,KKXK*KVAL,LWRK,
     1    KVAL*KVAL,0,0,0,0)
          CALL MEMO(-1,LIP1,KIP1,0,0,0,0,0,0,0,0)
        ELSE
C****************************************  LIKE VCCI2 (K+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
          NREC(1)=NREC(1)+1
          MAXXK=MAX0(KKXK,KTXK)
          MAXVAL=MAX0(KVAL,ITVAL)
          JCI=MAXBFN(W(LMXBAS),NMODE,K,1)
          IF(JREACT.GT.0)THEN
            JCIM=MAX0(NTAU(1)+1,NTAU(2)+1)
          ELSE
            JCIM1=2*(NBAS(1,1)+1)-1
            JCIM2=2*(NBAS(1,2)+1)-1
            JCIM=MAX0(JCIM1,JCIM2)
          END IF
          I9=1
          IF(ICOUPC.GT.0)I9=9
          KTEMP1=JCI*JCI*I9
          KXK0=I9*JCI*JCI*IK2
          KXL0=I9*JCIM*JCIM*IK2TAU
          CALL MEMO(1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
          CALL MEMO(2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
          IF(KKCONT.EQ.ITCONT)THEN
            CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
            CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1      0,0,0,0,0,0)
            IF(K.GT.NONLIN)THEN
              KPPP2=KLJP2
              IND2=1
            ELSE
              KPPP2=KIP2
              IND2=0
            END IF
            CALL MEMO(1,LIP2,KPPP2,0,0,0,0,0,0,0,0)
            IF(K.GT.NONLIN)THEN
              CALL GETBP1(W(LIP2),LSIZE2,NMODE,K,W(LMXBAS),NSIZE2,
     1        IND2)
              MSIZE2=NSIZE2
            ELSE
              CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,NNMODE,W(LMXBAS),
     1        NSIZE2,IND2)
              MSIZE2=ISIZE2
            END IF
            KXA2=NSIZE2*(NSIZE2+1)/2
            CALL MEMO(1,LXA1,KXA2,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE2,NSIZE2,W(LXA1),NSIZE2,W(LXA1),ICI)
            CALL VCCV1A(NAMODE,NNMODE,K,KMODX,ITMODE,W(LH+K2),
     1      W(LXQ+K3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IK3TAU,IK2TAU,
     2      W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),
     3      ISIZE,
     3      W(LIP2),MSIZE2,W(LXA1),W(LXA1),NSIZE2,W(LXK),W(LTEMP),
     4      W(LWRK),KVAL,ISTART,IEND,W(LTEMP1),JCI,JCIM,W(LXK0),
     5      W(LXL0),W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),W(LEJK1),
     6      J21,KROT,W(LMODNT),W(LOMEGA+K1),KEL,NS,
     7      W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(ITCONT-1)),ISIZXX,
     8      NVALCF,KEL21,NVSYM,I9)
            CALL MEMO(-2,LIP2,KPPP2,LXA1,KXA2,0,0,0,0,0,0)
            CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
            CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1      0,0,0,0,0,0)
          ELSE
C**NOT DONE FOR LINEAR
            ISIZC2=JCC2**2
            NMAX1=MAX0(MAXSUM(1,1)+2,MAXSUM(2,1)+2)
            NMAX2=MAX0(MAXSUM(1,2)+2,MAXSUM(2,2)+2)
            NMAX=NMAX1+NMAX2
C**FIRST GET 2-D SIZE FOR 'NORMAL' (ISIZC2 MODIFIED IN GETSZ FOR NMAX)
            CALL GETSZ(ISIZC2,2,JCC2,0)
            NSIZE2=ISIZC2
            MSIZE2=NSIZE2
            KCP2=ISIZC2*2
            CALL MEMO(1,LIP2,KCP2,0,0,0,0,0,0,0,0)
            CALL GETCP2(W(LIP2),ISIZC2,NMODE,K,NNMODE,W(LMXBAS),
     1      NSIZE2)
            KXA2=NSIZE2*(NSIZE2+1)/2
            CALL MEMO(1,LXA1,KXA2,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
            CALL DIAGZ(NSIZE2,NSIZE2,W(LXA1),NSIZE2,W(LXA1),ICI)
            NONC1=0
            NONC2=0
            CALL VCCV1B(NAMODE,NNMODE,K,KMODX,ITMODE,W(LH+K2),
     1      W(LXQ+K3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IK3TAU,IK2TAU,
     2      W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,KKIP,W(LKP+ITBAS),KTXK,
     3      KTIP,KKCONT,ITCONT,W(LXA),ISIZE,W(LIP2),MSIZE2,W(LXA1),
     4      W(LXA1),NSIZE2,W(LXK),W(LXW),W(LTEMP),MAXXK,
     5      W(LWRK),KVAL,ITVAL,ISTART,IEND,W(LTEMP1),JCI,JCIM,
     6      W(LXK0),W(LXL0),W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
     7      W(LEJK1),J21,KROT,W(LMODNT),W(LOMEGA+K1),KEL,NS,
     8      W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     9      W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(ITCONT-1)),
     1      ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I9,MAXVAL,MAXQ2)
C***************
            CALL MEMO(1,LXK,MAXQ2*MAXQ2,0,0,0,0,0,0,0,0)
            CALL MEMO(1,LXW,2*MAXVAL*MAXVAL*MAXVAL,0,0,0,0,0,0,0,0)
            CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),0,0,0,0,0,0)
            CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1      0,0,0,0,0,0)
            CALL VCCV1B(NAMODE,NNMODE,K,KMODX,ITMODE,W(LH+K2),
     1      W(LXQ+K3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IK3TAU,IK2TAU,
     2      W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,KKIP,W(LKP+ITBAS),KTXK,
     3      KTIP,KKCONT,ITCONT,W(LXA),ISIZE,W(LIP2),MSIZE2,W(LXA1),
     4      W(LXA1),NSIZE2,W(LXK),W(LXW),W(LTEMP),MAXXK,
     5      W(LWRK),KVAL,ITVAL,ISTART,IEND,W(LTEMP1),JCI,JCIM,
     6      W(LXK0),W(LXL0),W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
     7      W(LEJK1),J21,KROT,W(LMODNT),W(LOMEGA+K1),KEL,NS,
     8      W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     9      W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(ITCONT-1)),
     1      ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I9,MAXVAL,MAXQ2)
            CALL MEMO(-1,LXK,MAXQ2*MAXQ2,0,0,0,0,0,0,0,0)
            CALL MEMO(-1,LXW,2*MAXVAL*MAXVAL*MAXVAL,0,0,0,0,0,0,0,0)
            CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),0,0,0,0,0,0)
            CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1      0,0,0,0,0,0)
            CALL MEMO(-2,LIP2,KCP2,LXA1,KXA2,0,0,0,0,0,0)
          END IF
          CALL MEMO(-2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
          CALL MEMO(-1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
C**
C****************************************  LIKE VCCI2 (K+NMODE)
        END IF
CCCC    CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
331   CONTINUE
C**INTRINSIC
        IF(ICOUPX.LE.1)GO TO 5551
C**MATRIX FOR CENTRAL (DGEMM) INTEGRALS FOR THIS SCHEME
CCCC    CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)THEN
            ITIM2A=0
            ITIM2B=0
          END IF
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'L'
          DO LL=1,ICONT(1)
            IF(L.EQ.JCONT(1,LL))THEN
              LLCONT=1
              LMODX=LL
              LONBAS=0
              KLXK=ICSIZ1
              KLIP=IPSIZ1
              LVAL=NVAL1
            END IF
          END DO
          DO LL=1,ICONT(2)
            IF(L.EQ.JCONT(2,LL))THEN
              LLCONT=2
              LMODX=LL
              LONBAS=KLP
              KLXK=ICSIZ2
              KLIP=IPSIZ2
              LVAL=NVAL2
            END IF
          END DO
C**INTRINSIC
          IF(MOLINC.GT.0)THEN
            CALL GROUP2(K,L,ICONT(1),JCONT,1,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(1)
              JCOUPC=NCOUPC(1)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 3032
            END IF
            CALL GROUP2(K,L,ICONT(2),JCONT,2,IND)
            IF(IND.NE.0)THEN
              JCOUPL=NCOUPL(2)
              JCOUPC=NCOUPC(2)
              ICOUPL=IABS(JCOUPL)
              ICOUPC=IABS(JCOUPC)
              GO TO 3032
            END IF
            JCOUPL=JCOUPS
            JCOUPC=JCOUCS
            ICOUPL=ICOUPS
            ICOUPC=ICOUCS
3032        CONTINUE
          END IF
          IF(ICOUPL.LT.2)GO TO 332
C**CORIOLIS AND POTENTIAL
          IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN).OR.
     1    (NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR TWO MODES
            NREC(2)=NREC(2)+1
            MAXXK=MAX0(KKXK,KLXK)
            MAXVAL=MAX0(KVAL,LVAL)
C           CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
C    1      0,0,0,0,0,0)
            JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
            JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
            I5=1
            IF(ICOUPC.GT.1)I5=5
            KTEMP1=JCI2*JCI2*I5
            KXK0=I5*JCI1*JCI1*IK2
            KXL0=I5*JCI2*JCI2*IL2
            CALL MEMO(1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
            CALL MEMO(2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
            IF(KKCONT.EQ.LLCONT)THEN
              CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
              CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1        0,0,0,0,0,0)
              CALL MEMO(1,LIP2,KIP2,0,0,0,0,0,0,0,0)
              CALL GETBP2(W(LIP2),ISIZE2,NMODE,K,L,W(LMXBAS),NSIZE2,
     1        0)
              KXA2=NSIZE2*(NSIZE2+1)/2
              CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(NSIZE2,NSIZE2,W(LXA2),NSIZE2,W(LXA2),ICI)
              IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
                CALL MEMO(1,LOV,KKXK,0,0,0,0,0,0,0,0)
                CALL VCCI2C(NAMODE,NNMODE,K,L,KMODX,LMODX,W(LH+K2),
     1          W(LXQ+K3),W(LH+L2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LIP),
     2          ISIZMX,W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),ISIZE,
     3          W(LIP2),ISIZE2,W(LXA2),W(LXA2),NSIZE2,W(LXK),W(LTEMP),
     4          W(LOV),W(LWRK),KVAL,ISTART,IEND,W(LTEMP1),JCI1,JCI2,
     5          W(LXK0),W(LXL0),W(LV2),W(LV2),W(LC2),W(LC2),W(LEJK2),
     6          W(LEJK2),J21,KROT,W(LMODNT),KEL,NS,
     7          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),ISIZXX,
     8          NVALCF,KEL21,NVSYM,I5,W(LXKAN),MAXQU,MAXPOW,W(LNP2),
     9          W(LCP2),W(LMP2),NTOT2,IENTMX(2),W(LINDK))
                CALL MEMO(-1,LOV,KKXK,0,0,0,0,0,0,0,0)
              ELSE
                CALL VCCI2A(NAMODE,NNMODE,K,L,KMODX,LMODX,W(LH+K2),
     1          W(LXQ+K3),W(LH+L2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LIP),
     2          ISIZMX,W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),ISIZE,
     3          W(LIP2),ISIZE2,W(LXA2),W(LXA2),NSIZE2,W(LXK),W(LTEMP),
     4          W(LWRK),KVAL,ISTART,IEND,W(LTEMP1),JCI1,JCI2,
     5          W(LXK0),W(LXL0),W(LV2),W(LV2),W(LC2),W(LC2),W(LEJK2),
     6          W(LEJK2),J21,KROT,W(LMODNT),KEL,NS,
     7          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),ISIZXX,
     8          NVALCF,KEL21,NVSYM,I5,W(LXKAN),MAXQU,MAXPOW,W(LNP2),
     9          W(LCP2),W(LMP2),NTOT2,IENTMX(2),W(LINDK))
              END IF
              CALL MEMO(-2,LIP2,KIP2,LXA2,KXA2,0,0,0,0,0,0)
              CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
              CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1        0,0,0,0,0,0)
            ELSE
              ISIZC2=JCC2**2
              NMAX1=MAX0(MAXSUM(1,1)+2,MAXSUM(2,1)+2)
              NMAX2=MAX0(MAXSUM(1,2)+2,MAXSUM(2,2)+2)
              NMAX=NMAX1+NMAX2
C**FIRST GET 2-D SIZE FOR 'NORMAL' (ISIZC2 MODIFIED IN GETSZ FOR NMAX)
              CALL GETSZ(ISIZC2,2,JCC2,0)
              NSIZE2=ISIZC2
              KCP2=ISIZC2*2
              CALL MEMO(1,LIP2,KCP2,0,0,0,0,0,0,0,0)
              CALL GETCP2(W(LIP2),ISIZC2,NMODE,K,L,W(LMXBAS),NSIZE2)
              KXA2=NSIZE2*(NSIZE2+1)/2
              CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(NSIZE2,NSIZE2,W(LXA2),NSIZE2,W(LXA2),ICI)
              NONC1=0
              NONC2=0
              IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
                CALL VCCI2D(NAMODE,NNMODE,K,L,KMODX,LMODX,W(LH+K2),
     1          W(LXQ+K3),W(LH+L2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LIP),
     2          ISIZMX,W(LKP+KONBAS),KKXK,KKIP,W(LKP+LONBAS),KLXK,KLIP,
     3          KKCONT,LLCONT,W(LXA),ISIZE,W(LIP2),ISIZC2,W(LXA2),
     4          W(LXA2),NSIZE2,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),
     5          KVAL,LVAL,ISTART,IEND,W(LTEMP1),
     6          JCI1,JCI2,W(LXK0),W(LXL0),W(LV2),W(LV2),W(LC2),W(LC2),
     7          W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),KEL,NS,
     8          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     9          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(LLCONT-1)),
     1          ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I5,W(LXKAN),
     2          MAXQU,MAXPOW,W(LNP2),W(LCP2),W(LMP2),NTOT2,IENTMX(2),
     3          W(LINDK),MAXVAL)
C***************
CCCC            CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                CALL MEMO(1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                NNC12=MAX0(NONC1,NONC2)
                CALL MEMO(1,LXW,NNC12*MAXVAL,0,0,0,0,0,0,0,0)
                CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1          0,0,0,0,0,0)
                CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1          0,0,0,0,0,0)
                CALL VCCI2D(NAMODE,NNMODE,K,L,KMODX,LMODX,W(LH+K2),
     1          W(LXQ+K3),W(LH+L2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LIP),
     2          ISIZMX,W(LKP+KONBAS),KKXK,KKIP,W(LKP+LONBAS),KLXK,KLIP,
     3          KKCONT,LLCONT,W(LXA),ISIZE,W(LIP2),ISIZC2,W(LXA2),
     4          W(LXA2),NSIZE2,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),
     5          KVAL,LVAL,ISTART,IEND,W(LTEMP1),
     6          JCI1,JCI2,W(LXK0),W(LXL0),W(LV2),W(LV2),W(LC2),W(LC2),
     7          W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),KEL,NS,
     8          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     9          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(LLCONT-1)),
     1          ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I5,W(LXKAN),
     2          MAXQU,MAXPOW,W(LNP2),W(LCP2),W(LMP2),NTOT2,IENTMX(2),
     3          W(LINDK),MAXVAL)
                CALL MEMO(-1,LXW,NNC12*MAXVAL,0,0,0,0,0,0,0,0)
                CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1          0,0,0,0,0,0)
                CALL MEMO(-1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
CCCC            CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1          0,0,0,0,0,0)
              ELSE
                CALL VCCI2B(NAMODE,NNMODE,K,L,KMODX,LMODX,W(LH+K2),
     1          W(LXQ+K3),W(LH+L2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LIP),
     2          ISIZMX,W(LKP+KONBAS),KKXK,KKIP,W(LKP+LONBAS),KLXK,KLIP,
     3          KKCONT,LLCONT,W(LXA),ISIZE,W(LIP2),ISIZC2,W(LXA2),
     4          W(LXA2),NSIZE2,W(LXK),W(LXW),W(LTEMP),
     5          MAXXK,W(LWRK),KVAL,LVAL,ISTART,IEND,W(LTEMP1),
     6          JCI1,JCI2,W(LXK0),W(LXL0),W(LV2),W(LV2),W(LC2),W(LC2),
     7          W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),KEL,NS,
     8          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     9          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(LLCONT-1)),
     1          ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I5,W(LXKAN),
     2          MAXQU,MAXPOW,W(LNP2),W(LCP2),W(LMP2),NTOT2,IENTMX(2),
     3          W(LINDK),MAXVAL,MAXQ2)
C***************
CCCC            CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                CALL MEMO(1,LXK,MAXQ2*MAXQ2,0,0,0,0,0,0,0,0)
                CALL MEMO(1,LXW,2*MAXVAL*MAXVAL*MAXVAL,0,0,0,0,0,0,0,0)
                CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1          0,0,0,0,0,0)
                CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1          0,0,0,0,0,0)
                CALL VCCI2B(NAMODE,NNMODE,K,L,KMODX,LMODX,W(LH+K2),
     1          W(LXQ+K3),W(LH+L2),W(LXQ+L3),IK3,IK2,IL3,IL2,W(LIP),
     2          ISIZMX,W(LKP+KONBAS),KKXK,KKIP,W(LKP+LONBAS),KLXK,KLIP,
     3          KKCONT,LLCONT,W(LXA),ISIZE,W(LIP2),ISIZC2,W(LXA2),
     4          W(LXA2),NSIZE2,W(LXK),W(LXW),W(LTEMP),
     5          MAXXK,W(LWRK),KVAL,LVAL,ISTART,IEND,W(LTEMP1),
     6          JCI1,JCI2,W(LXK0),W(LXL0),W(LV2),W(LV2),W(LC2),W(LC2),
     7          W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),KEL,NS,
     8          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     9          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(LLCONT-1)),
     1          ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I5,W(LXKAN),
     2          MAXQU,MAXPOW,W(LNP2),W(LCP2),W(LMP2),NTOT2,IENTMX(2),
     3          W(LINDK),MAXVAL,MAXQ2)
                CALL MEMO(-1,LXK,MAXQ2*MAXQ2,0,0,0,0,0,0,0,0)
CCCC            CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                CALL MEMO(-1,LXW,2*MAXVAL*MAXVAL*MAXVAL,0,0,0,0,0,0,0,0)
                CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1          0,0,0,0,0,0)
                CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1          0,0,0,0,0,0)
              END IF
              CALL MEMO(-2,LIP2,KCP2,LXA2,KXA2,0,0,0,0,0,0)
            END IF
            CALL MEMO(-2,LXK0,KXK0,LXL0,KXL0,0,0,0,0,0,0)
            CALL MEMO(-1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
C           CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,
C    1      MAXVAL*MAXVAL,0,0,0,0,0,0)
          ELSE
C*****************************   LIKE VCCI3 (K+L+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
            NREC(2)=NREC(2)+1
            MAXXK=MAX0(KKXK,KLXK,KTXK)
            MAXVAL=MAX0(KVAL,LVAL,ITVAL)
            JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
            JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
            IF(JREACT.GT.0)THEN
              JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1)
            ELSE
              JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1)-1
              JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1)-1
              JCIM=MAX0(JCIM1,JCIM2)
            END IF
            I16=1
            IF(ICOUPC.GT.1)I16=16
            KTEMP1=JCI1*JCI2*JCI1*JCI2*I16
            KTEMP2=JCI2*JCI2*I16
            KXK0=I16*JCI1*JCI1*IK2
            KXL0=I16*JCI2*JCI2*IL2
            KXN0=I16*JCIM*JCIM*IK2TAU
            CALL MEMO(2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
            CALL MEMO(3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
            IF(KKCONT.EQ.LLCONT.AND.KKCONT.EQ.ITCONT)THEN
              CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
              CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1        0,0,0,0,0,0)
              IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
                KPPP3=KLJP3
                IND3=1
              ELSE
                KPPP3=KIP3
                IND3=0
              END IF
              CALL MEMO(1,LIP3,KPPP3,0,0,0,0,0,0,0,0)
              IF(K.GT.NONLIN.OR.L.GT.NONLIN)THEN
                CALL GETBP2(W(LIP3),LSIZE3,NMODE,K,L,W(LMXBAS),
     1          NSIZE3,IND3)
                MSIZE3=NSIZE3
              ELSE
                CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,NNMODE,
     1          W(LMXBAS),NSIZE3,IND3)
                MSIZE3=ISIZE3
              END IF
              KXA3=NSIZE3*(NSIZE3+1)/2
              CALL MEMO(1,LXA2,KXA3,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(NSIZE3,NSIZE3,W(LXA2),NSIZE3,W(LXA2),ICI)
              CALL VCCV2A(NAMODE,NNMODE,K,L,KMODX,LMODX,ITMODE,
     1        W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+K2TAU),
     2        W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IK3TAU,IK2TAU,W(LIP),
     3        ISIZMX,W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),ISIZE,
     4        W(LIP3),
     4        MSIZE3,W(LXA2),W(LXA2),NSIZE3,W(LXK),W(LTEMP),W(LWRK),
     5        KVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),JCI1,JCI2,JCIM,
     6        W(LXK0),W(LXL0),W(LXN0),W(LV2),W(LV2),W(LC2),W(LC2),
     7        W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),KEL,NS,
     8        W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(ITCONT-1)),ISIZXX,
     9        NVALCF,KEL21,NVSYM,I16)
              CALL MEMO(-2,LIP3,KPPP3,LXA2,KXA3,0,0,0,0,0,0)
              CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
              CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1        0,0,0,0,0,0)
            ELSE
C**NOT DONE FOR LINEAR
              ISIZC3=JCC3**3
              NMAX1=MAX0(MAXSUM(1,1)+3,MAXSUM(2,1)+3,MAXSUM(3,1)+3)
              NMAX2=MAX0(MAXSUM(1,2)+3,MAXSUM(2,2)+3,MAXSUM(3,2)+3)
              NMAX=NMAX1+NMAX2
C**FIRST GET 3-D SIZE FOR 'NORMAL' (ISIZE3 MODIFIED IN GETSZ FOR NMAX)
              CALL GETSZ(ISIZC3,3,JCC3,0)
              NSIZE3=ISIZC3
              MSIZE3=NSIZE3
              KCP3=ISIZC3*3
              CALL MEMO(1,LIP3,KCP3,0,0,0,0,0,0,0,0)
              CALL GETCP3(W(LIP3),ISIZC3,NMODE,K,L,NNMODE,W(LMXBAS),
     1        NSIZE3)
              KXA3=NSIZE3*(NSIZE3+1)/2
              CALL MEMO(1,LXA2,KXA3,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
              CALL DIAGZ(NSIZE3,NSIZE3,W(LXA2),NSIZE3,W(LXA2),ICI)
              NONC1=0
              NONC2=0
              IF(KKCONT.NE.LLCONT)THEN
                KKKBAS=LONBAS
                KKKCON=LLCONT
                KKKXK=KLXK
                KKKIP=KLIP
                KKKVAL=LVAL
              ELSE
                KKKBAS=ITBAS
                KKKCON=ITCONT
                KKKXK=KTXK
                KKKIP=KTIP
                KKKVAL=ITVAL
              END IF
              CALL VCCV2B(NAMODE,NNMODE,K,L,KMODX,LMODX,ITMODE,
     1        W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+K2TAU),
     2        W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IK3TAU,IK2TAU,W(LIP),
     3        ISIZMX,W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,
     4        KKKIP,KKCONT,
     4        KKKCON,W(LXA),ISIZE,W(LIP3),MSIZE3,W(LXA2),W(LXA2),
     5        NSIZE3,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),
     6        KVAL,KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),JCI1,JCI2,
     7        JCIM,W(LXK0),W(LXL0),W(LXN0),W(LV2),W(LV2),W(LC2),
     8        W(LC2),W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),KEL,NS,
     9        W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     1        W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     2        ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I16,MAXVAL,MAXQ3)
C***************
              CALL MEMO(1,LXK,MAXQ3*MAXQ3*MAXQ3,0,0,0,0,0,0,0,0)
              CALL MEMO(1,LXW,2*MAXVAL*MAXVAL*MAXVAL,0,0,0,0,0,0,0,0)
              CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),0,0,0,0,0,0)
              CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1        0,0,0,0,0,0)
              CALL VCCV2B(NAMODE,NNMODE,K,L,KMODX,LMODX,ITMODE,
     1        W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+K2TAU),
     2        W(LXQ+K3TAU),IK3,IK2,IL3,IL2,IK3TAU,IK2TAU,W(LIP),
     3        ISIZMX,W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,
     4        KKKIP,KKCONT,
     4        KKKCON,W(LXA),ISIZE,W(LIP3),MSIZE3,W(LXA2),W(LXA2),
     5        NSIZE3,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),
     6        KVAL,KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),JCI1,JCI2,
     7        JCIM,W(LXK0),W(LXL0),W(LXN0),W(LV2),W(LV2),W(LC2),
     8        W(LC2),W(LEJK2),W(LEJK2),J21,KROT,W(LMODNT),KEL,NS,
     9        W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     1        W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     2        ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I16,MAXVAL,
     3        MAXQ3)
              CALL MEMO(-1,LXK,MAXQ3*MAXQ3*MAXQ3,0,0,0,0,0,0,0,0)
              CALL MEMO(-1,LXW,2*MAXVAL*MAXVAL*MAXVAL,0,0,0,0,0,0,0,0)
              CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1        0,0,0,0,0,0)
              CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1        0,0,0,0,0,0)
              CALL MEMO(-2,LIP3,KCP3,LXA2,KXA3,0,0,0,0,0,0)
            END IF
            CALL MEMO(-2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
            CALL MEMO(-3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
C**
C*****************************   LIKE VCCI3 (K+L+NMODE)
          END IF
332   CONTINUE
C**INTRINSIC
          IF(ICOUPX.EQ.2)GO TO 5552
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.ITIM.EQ.0)THEN
              ITIM3A=0
              ITIM3B=0
            END IF
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'N'
            DO NN=1,ICONT(1)
              IF(N.EQ.JCONT(1,NN))THEN
                NNCONT=1
                NMODX=NN
                NONBAS=0
                KNXK=ICSIZ1
                KNIP=IPSIZ1
                NVAL=NVAL1
              END IF
            END DO
            DO NN=1,ICONT(2)
              IF(N.EQ.JCONT(2,NN))THEN
                NNCONT=2
                NMODX=NN
                NONBAS=KLP
                KNXK=ICSIZ2
                KNIP=IPSIZ2
                NVAL=NVAL2
              END IF
            END DO
C**INTRINSIC
            IF(MOLINC.GT.0)THEN
              CALL GROUP3(K,L,N,ICONT(1),JCONT,1,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(1)
                JCOUPC=NCOUPC(1)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 3033
              END IF
              CALL GROUP3(K,L,N,ICONT(2),JCONT,2,IND)
              IF(IND.NE.0)THEN
                JCOUPL=NCOUPL(2)
                JCOUPC=NCOUPC(2)
                ICOUPL=IABS(JCOUPL)
                ICOUPC=IABS(JCOUPC)
                GO TO 3033
              END IF
              JCOUPL=JCOUPS
              JCOUPC=JCOUCS
              ICOUPL=ICOUPS
              ICOUPC=ICOUCS
3033          CONTINUE
            END IF
            IF(ICOUPL.LT.3)GO TO 333
C**CORIOLIS AND POTENTIAL
            IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1      N.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR THREE MODES
              NREC(3)=NREC(3)+1
              MAXXK=MAX0(KKXK,KLXK,KNXK)
              MAXVAL=MAX0(KVAL,LVAL,NVAL)
C             CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
C    1        0,0,0,0,0,0)
              JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
              JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
              JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
              I10=1
              IF(ICOUPC.GT.2)I10=10
              KTEMP1=JCI3*JCI3*I10
              KTEMP2=JCI2*JCI3*JCI2*JCI3*I10
              KXK0=I10*JCI1*JCI1*IK2
              KXL0=I10*JCI2*JCI2*IL2
              KXN0=I10*JCI3*JCI3*IN2
              CALL MEMO(2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
              CALL MEMO(3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
              IF(KKCONT.EQ.LLCONT.AND.KKCONT.EQ.NNCONT)THEN
                CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1          0,0,0,0,0,0)
                CALL MEMO(1,LIP3,KIP3,0,0,0,0,0,0,0,0)
                CALL GETBP3(W(LIP3),ISIZE3,NMODE,K,L,N,W(LMXBAS),
     1          NSIZE3,0)
                KXA3=NSIZE3*(NSIZE3+1)/2
                CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                CALL DIAGZ(NSIZE3,NSIZE3,W(LXA3),NSIZE3,W(LXA3),ICI)
                IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
                  CALL MEMO(1,LOV,KKXK,0,0,0,0,0,0,0,0)
                  CALL VCCI3C(NAMODE,NNMODE,K,L,N,KMODX,LMODX,
     1            NMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2            W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,IN2,W(LIP),ISIZMX,
     3            W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),ISIZE,W(LIP3),
     4            ISIZE3,W(LXA3),W(LXA3),NSIZE3,W(LXK),W(LTEMP),
     5            W(LOV),W(LWRK),KVAL,ISTART,IEND,W(LTEMP1),
     6            W(LTEMP2),JCI1,JCI2,JCI3,W(LXK0),W(LXL0),W(LXN0),
     7            W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,
     8            KROT,W(LMODNT),KEL,NS,
     9            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),ISIZXX,
     1            NVALCF,KEL21,NVSYM,I10,W(LXKAN),MAXQU,MAXPOW,W(LNP3),
     2            W(LCP3),W(LMP3),NTOT3,IENTMX(3),W(LINDK),W(LINDL))
                  CALL MEMO(-1,LOV,KKXK,0,0,0,0,0,0,0,0)
                ELSE
                  CALL VCCI3A(NAMODE,NNMODE,K,L,N,KMODX,LMODX,
     1            NMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2            W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,IN2,W(LIP),ISIZMX,
     3            W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),ISIZE,W(LIP3),
     4            ISIZE3,W(LXA3),W(LXA3),NSIZE3,W(LXK),W(LTEMP),
     5            W(LWRK),KVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),
     6            JCI1,JCI2,JCI3,W(LXK0),W(LXL0),W(LXN0),
     7            W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,
     8            KROT,W(LMODNT),KEL,NS,
     9            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),ISIZXX,
     1            NVALCF,KEL21,NVSYM,I10,W(LXKAN),MAXQU,MAXPOW,W(LNP3),
     2            W(LCP3),W(LMP3),NTOT3,IENTMX(3),W(LINDK),W(LINDL))
                END IF
                CALL MEMO(-2,LIP3,KIP3,LXA3,KXA3,0,0,0,0,0,0)
                CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1          0,0,0,0,0,0)
              ELSE
                ISIZC3=JCC3**3
                NMAX1=MAX0(MAXSUM(1,1)+3,MAXSUM(2,1)+3,MAXSUM(3,1)+3)
                NMAX2=MAX0(MAXSUM(1,2)+3,MAXSUM(2,2)+3,MAXSUM(3,2)+3)
                NMAX=NMAX1+NMAX2
C**FIRST GET 3-D SIZE FOR 'NORMAL' (ISIZE3 MODIFIED IN GETSZ FOR NMAX)
                CALL GETSZ(ISIZC3,3,JCC3,0)
                NSIZE3=ISIZC3
                KCP3=ISIZC3*3
                CALL MEMO(1,LIP3,KCP3,0,0,0,0,0,0,0,0)
                CALL GETCP3(W(LIP3),ISIZC3,NMODE,K,L,N,W(LMXBAS),
     1          NSIZE3)
                KXA3=NSIZE3*(NSIZE3+1)/2
                CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                CALL DIAGZ(NSIZE3,NSIZE3,W(LXA3),NSIZE3,W(LXA3),ICI)
                NONC1=0
                NONC2=0
                IF(KKCONT.NE.LLCONT)THEN
                  KKKBAS=LONBAS
                  KKKCON=LLCONT
                  KKKXK=KLXK
                  KKKIP=KLIP
                  KKKVAL=LVAL
                ELSE
                  KKKBAS=NONBAS
                  KKKCON=NNCONT
                  KKKXK=KNXK
                  KKKIP=KNIP
                  KKKVAL=NVAL
                END IF
                IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
                  CALL VCCI3D(NAMODE,NNMODE,K,L,N,KMODX,LMODX,NMODX,
     1            W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2            W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,IN2,W(LIP),ISIZMX,
     3            W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,
     4            KKCONT,KKKCON,W(LXA),ISIZE,W(LIP3),ISIZC3,W(LXA3),
     5            W(LXA3),NSIZE3,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),
     6            KVAL,KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),
     7            JCI1,JCI2,JCI3,W(LXK0),W(LXL0),W(LXN0),W(LV3),W(LV3),
     8            W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,KROT,W(LMODNT),
     9            KEL,NS,W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     2            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     3            ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I10,W(LXKAN),
     4            MAXQU,MAXPOW,W(LNP3),W(LCP3),W(LMP3),NTOT3,IENTMX(3),
     5            W(LINDK),W(LINDL),MAXVAL)
C***************
CCCC              CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                  NNC12=MAX0(NONC1,NONC2)
                  CALL MEMO(1,LXW,NNC12*MAXVAL,0,0,0,0,0,0,0,0)
                  CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1            0,0,0,0,0,0)
                  CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1            0,0,0,0,0,0)
                  CALL VCCI3D(NAMODE,NNMODE,K,L,N,KMODX,LMODX,NMODX,
     1            W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2            W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,IN2,W(LIP),ISIZMX,
     3            W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,
     4            KKCONT,KKKCON,W(LXA),ISIZE,W(LIP3),ISIZC3,W(LXA3),
     5            W(LXA3),NSIZE3,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),
     6            KVAL,KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),
     7            JCI1,JCI2,JCI3,W(LXK0),W(LXL0),W(LXN0),W(LV3),W(LV3),
     8            W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,KROT,W(LMODNT),
     9            KEL,NS,W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     2            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     3            ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I10,
     4            W(LXKAN),MAXQU,MAXPOW,W(LNP3),W(LCP3),W(LMP3),NTOT3,
     5            IENTMX(3),W(LINDK),W(LINDL),MAXVAL)
                  CALL MEMO(-1,LXW,NNC12*MAXVAL,0,0,0,0,0,0,0,0)
                  CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1            0,0,0,0,0,0)
                  CALL MEMO(-1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
CCCC              CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1            0,0,0,0,0,0)
                ELSE
                  CALL VCCI3B(NAMODE,NNMODE,K,L,N,KMODX,LMODX,NMODX,
     1            W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2            W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,IN2,W(LIP),ISIZMX,
     3            W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,
     4            KKCONT,KKKCON,
     4            W(LXA),ISIZE,W(LIP3),ISIZC3,W(LXA3),W(LXA3),NSIZE3,
     5            W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,KKKVAL,
     7            ISTART,IEND,W(LTEMP1),W(LTEMP2),JCI1,JCI2,
     8            JCI3,W(LXK0),W(LXL0),W(LXN0),W(LV3),W(LV3),W(LC3),
     9            W(LC3),W(LEJK3),W(LEJK3),J21,KROT,W(LMODNT),KEL,NS,
     1            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     2            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     3            ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I10,W(LXKAN),
     4            MAXQU,MAXPOW,W(LNP3),W(LCP3),W(LMP3),NTOT3,IENTMX(3),
     5            W(LINDK),W(LINDL),MAXVAL,MAXQ3)
C***************
CCCC              CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(1,LXK,MAXQ3*MAXQ3*MAXQ3,0,0,0,0,0,0,0,0)
                  CALL MEMO(1,LXW,2*MAXVAL*MAXVAL*MAXVAL,
     1            0,0,0,0,0,0,0,0)
                  CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1            0,0,0,0,0,0)
                  CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,
     1            MAXVAL*MAXVAL,0,0,0,0,0,0)
                  CALL VCCI3B(NAMODE,NNMODE,K,L,N,KMODX,LMODX,NMODX,
     1            W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     2            W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,IN2,W(LIP),ISIZMX,
     3            W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,
     4            KKCONT,KKKCON,
     4            W(LXA),ISIZE,W(LIP3),ISIZC3,W(LXA3),W(LXA3),NSIZE3,
     5            W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,KKKVAL,
     7            ISTART,IEND,W(LTEMP1),W(LTEMP2),JCI1,JCI2,
     8            JCI3,W(LXK0),W(LXL0),W(LXN0),W(LV3),W(LV3),W(LC3),
     9            W(LC3),W(LEJK3),W(LEJK3),J21,KROT,W(LMODNT),KEL,NS,
     1            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     2            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     3            ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I10,
     4            W(LXKAN),MAXQU,MAXPOW,W(LNP3),W(LCP3),W(LMP3),NTOT3,
     5            IENTMX(3),W(LINDK),W(LINDL),MAXVAL,MAXQ3)
                  CALL MEMO(-1,LXK,MAXQ3*MAXQ3*MAXQ3,0,0,0,0,0,0,0,0)
CCCC              CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(-1,LXW,2*MAXVAL*MAXVAL*MAXVAL,
     1            0,0,0,0,0,0,0,0)
                  CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1            0,0,0,0,0,0)
                  CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,
     1            MAXVAL*MAXVAL,0,0,0,0,0,0)
                END IF
                CALL MEMO(-2,LIP3,KCP3,LXA3,KXA3,0,0,0,0,0,0)
              END IF
C             CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,
C    1        MAXVAL*MAXVAL,0,0,0,0,0,0)
              CALL MEMO(-2,LTEMP1,KTEMP1,LTEMP2,KTEMP2,0,0,0,0,0,0)
              CALL MEMO(-3,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,0,0,0,0)
            ELSE
C*****************************   LIKE VCCI4 (K+L+N+NMODE)
C**
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
              NREC(3)=NREC(3)+1
              MAXXK=MAX0(KKXK,KLXK,KNXK,KTXK)
              MAXVAL=MAX0(KVAL,LVAL,NVAL,ITVAL)
              JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
              JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
              JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
              IF(JREACT.GT.0)THEN
                JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1,NTAU(4)+1)
              ELSE
                JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1,NBAS(3,1)+1)-1
                JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1,NBAS(3,2)+1)-1
                JCI=MAX0(JCIM1,JCIM2)
              END IF
              I25=1
              IF(ICOUPC.GT.2)I25=25
              KTEMP1=JCI1*JCI2*JCI3*JCI1*JCI2*JCI3*I25
              KTEMP2=JCI2*JCI3*JCI2*JCI3*I25
              KTEMP3=JCI3*JCI3*I25
              KXK0=I25*JCI1*JCI1*IK2
              KXL0=I25*JCI2*JCI2*IL2
              KXN0=I25*JCI3*JCI3*IN2
              KXM0=I25*JCIM*JCIM*IK2TAU
              CALL MEMO(3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,KTEMP3,
     1        0,0,0,0)
              CALL MEMO(4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1        0,0)
              IF(KKCONT.EQ.LLCONT.AND.KKCONT.EQ.NNCONT.AND.KKCONT.EQ.
     1        ITCONT)THEN
                CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1          0,0,0,0,0,0)
                IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                  KPPP4=KLJP4
                  IND4=1
                ELSE
                  KPPP4=KIP4
                  IND4=0
                END IF
                CALL MEMO(1,LIP4,KPPP4,0,0,0,0,0,0,0,0)
                IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN)THEN
                  CALL GETBP3(W(LIP4),LSIZE4,NMODE,K,L,N,W(LMXBAS),
     1            NSIZE4,IND4)
                  MSIZE4=NSIZE4
                ELSE
                  CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,NNMODE,
     1            W(LMXBAS),NSIZE4,IND4)
                  MSIZE4=ISIZE4
                END IF

!                KXA4=NSIZE4*(NSIZE4+1)/2
                if (mod(nsize4,2) .eq. 0) then
                   kxa4 = nsize4 / 2
                   kxa4 = kxa4 * (nsize4 + 1)
                else
                   kxa4 = (nsize4 + 1) / 2
                   kxa4 = kxa4 * nsize4
                end if

                CALL MEMO(1,LXA3,KXA4,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                CALL DIAGZ(NSIZE4,NSIZE4,W(LXA3),NSIZE4,W(LXA3),ICI)
                CALL VCCV3A(NAMODE,NNMODE,K,L,N,KMODX,LMODX,
     1          NMODX,ITMODE,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),
     2          W(LH+N2),W(LXQ+N3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,
     3          IL3,IL2,IN3,IN2,IK3TAU,IK2TAU,W(LIP),ISIZMX,
     4          W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),ISIZE,W(LIP4),
     5          MSIZE4,W(LXA3),W(LXA3),NSIZE4,W(LXK),W(LTEMP),
     6          W(LWRK),KVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),
     7          W(LTEMP3),JCI1,JCI2,JCI3,JCIM,W(LXK0),W(LXL0),
     8          W(LXN0),W(LXM0),W(LV3),W(LV3),W(LC3),W(LC3),W(LEJK3),
     9          W(LEJK3),J21,KROT,W(LMODNT),KEL,NS,
     1          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(ITCONT-1)),ISIZXX,
     2          NVALCF,KEL21,NVSYM,I25)
                CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1          0,0,0,0,0,0)
                CALL MEMO(-2,LIP4,KPPP4,LXA3,KXA4,0,0,
     1          0,0,0,0)
                CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,KTEMP3,
     1          0,0,0,0)
                CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          0,0)
              ELSE
C**NOT DONE FOR LINEAR
                ISIZC4=JCC4**4
                NMAX1=MAX0(MAXSUM(1,1)+4,MAXSUM(2,1)+4,MAXSUM(3,1)+4,
     1          MAXSUM(4,1)+4)
                NMAX2=MAX0(MAXSUM(1,2)+4,MAXSUM(2,2)+4,MAXSUM(3,2)+4,
     1          MAXSUM(4,2)+4)
                NMAX=NMAX1+NMAX2
C**FIRST GET 4-D SIZE FOR 'NORMAL' (ISIZE4 MODIFIED IN GETSZ FOR NMAX)
                CALL GETSZ(ISIZC4,4,JCC4,0)
                NSIZE4=ISIZC4
                MSIZE4=NSIZE4
                KCP4=ISIZC4*4
                CALL MEMO(1,LIP4,KCP4,0,0,0,0,0,0,0,0)
                CALL GETCP4(W(LIP4),ISIZC4,NMODE,K,L,N,NNMODE,
     1          W(LMXBAS),NSIZE4)

!                KXA4=NSIZE4*(NSIZE4+1)/2
                if (mod(nsize4,2) .eq. 0) then
                   kxa4 = nsize4 / 2
                   kxa4 = kxa4 * (nsize4 + 1)
                else
                   kxa4 = (nsize4 + 1) / 2
                   kxa4 = kxa4 * nsize4
                end if

                CALL MEMO(1,LXA3,KXA4,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                CALL DIAGZ(NSIZE4,NSIZE4,W(LXA3),NSIZE4,W(LXA3),ICI)
                NONC1=0
                NONC2=0
                IF(KKCONT.NE.LLCONT)THEN
                  KKKBAS=LONBAS
                  KKKCON=LLCONT
                  KKKXK=KLXK
                  KKKIP=KLIP
                  KKKVAL=LVAL
                ELSE
                  IF(KKCONT.NE.NNCONT)THEN
                    KKKBAS=NONBAS
                    KKKCON=NNCONT
                    KKKXK=KNXK
                    KKKIP=KNIP
                    KKKVAL=NVAL
                  ELSE
                    KKKBAS=ITBAS
                    KKKCON=ITCONT
                    KKKXK=KTXK
                    KKKIP=KTIP
                    KKKVAL=ITVAL
                  END IF
                END IF
                CALL VCCV3B(NAMODE,NNMODE,K,L,N,KMODX,LMODX,
     1          NMODX,ITMODE,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),
     2          W(LH+N2),W(LXQ+N3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,
     3          IL3,IL2,IN3,IN2,IK3TAU,IK2TAU,W(LIP),ISIZMX,
     4          W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,
     5          KKCONT,KKKCON,
     5          W(LXA),ISIZE,W(LIP4),MSIZE4,W(LXA3),W(LXA3),NSIZE4,
     6          W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,KKKVAL,
     7          ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),JCI1,JCI2,
     8          JCI3,JCIM,W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LV3),
     9          W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,KROT,
     1          W(LMODNT),KEL,NS,
     2          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     3          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     4          ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I25,MAXVAL,MAXQ4)
                CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,KTEMP3,
     1          0,0,0,0)
                CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          0,0)
C***************
                CALL MEMO(1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                CALL MEMO(1,LXW,2*MAXVAL*MAXVAL*MAXVAL,0,0,0,0,0,0,0,0)
                CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1          0,0,0,0,0,0)
                CALL MEMO(2,LTEMP,MAXQ4*MAXQ4*MAXQ4*MAXQ4,
     1          LWRK,MAXVAL*MAXVAL,0,0,0,0,0,0)
                CALL VCCV3B(NAMODE,NNMODE,K,L,N,KMODX,LMODX,
     1          NMODX,ITMODE,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),
     2          W(LH+N2),W(LXQ+N3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,
     3          IL3,IL2,IN3,IN2,IK3TAU,IK2TAU,W(LIP),ISIZMX,
     4          W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,
     5          KKCONT,KKKCON,
     5          W(LXA),ISIZE,W(LIP4),MSIZE4,W(LXA3),W(LXA3),NSIZE4,
     6          W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,KKKVAL,
     7          ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),JCI1,JCI2,
     8          JCI3,JCIM,W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LV3),
     9          W(LV3),W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,KROT,
     1          W(LMODNT),KEL,NS,
     2          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     3          W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     4          ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I25,MAXVAL,
     5          MAXQ4)
                CALL MEMO(-1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                CALL MEMO(-1,LXW,2*MAXVAL*MAXVAL*MAXVAL,
     1          0,0,0,0,0,0,0,0)
                CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1          0,0,0,0,0,0)
                CALL MEMO(-2,LTEMP,MAXQ4*MAXQ4*MAXQ4*MAXQ4,
     1          LWRK,MAXVAL*MAXVAL,0,0,0,0,0,0)
                CALL MEMO(-2,LIP4,KCP4,LXA3,KXA4,0,0,0,0,0,0)
              END IF
C**
C*****************************   LIKE VCCI4 (K+L+N+NMODE)
            END IF
333   CONTINUE
C**INTRINSIC
            IF(ICOUPX.EQ.3)GO TO 5553
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.K.EQ.4.AND.
     1        ITIM.EQ.0)THEN
                ITIM4A=0
                ITIM4B=0
              END IF
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'N'
              DO MM=1,ICONT(1)
                IF(M.EQ.JCONT(1,MM))THEN
                  MMCONT=1
                  MMODX=MM
                  MONBAS=0
                  KMXK=ICSIZ1
                  KMIP=IPSIZ1
                  MVAL=NVAL1
                END IF
              END DO
              DO MM=1,ICONT(2)
                IF(M.EQ.JCONT(2,MM))THEN
                  MMCONT=2
                  MMODX=MM
                  MONBAS=KLP
                  KMXK=ICSIZ2
                  KMIP=IPSIZ2
                  MVAL=NVAL2
                END IF
              END DO
C**INTRINSIC
              IF(MOLINC.GT.0)THEN
                CALL GROUP4(K,L,N,M,ICONT(1),JCONT,1,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(1)
                  JCOUPC=NCOUPC(1)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 3034
                END IF
                CALL GROUP4(K,L,N,M,ICONT(2),JCONT,2,IND)
                IF(IND.NE.0)THEN
                  JCOUPL=NCOUPL(2)
                  JCOUPC=NCOUPC(2)
                  ICOUPL=IABS(JCOUPL)
                  ICOUPC=IABS(JCOUPC)
                  GO TO 3034
                END IF
                JCOUPL=JCOUPS
                JCOUPC=JCOUCS
                ICOUPL=ICOUPS
                ICOUPC=ICOUCS
3034            CONTINUE
              END IF
              IF(ICOUPL.LT.4)GO TO 334
C**CORIOLIS AND POTENTIAL
              IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1        N.LE.NONLIN.AND.M.LE.NONLIN).OR.(NSMODE.EQ.NAMODE))
     2        THEN
C**GET BASIC INTEGRAL BASIS FOR FOUR MODES
                NREC(4)=NREC(4)+1
                MAXXK=MAX0(KKXK,KLXK,KNXK,KMXK)
                MAXVAL=MAX0(KVAL,LVAL,NVAL,MVAL)
C               CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
C    1          0,0,0,0,0,0)
                JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
                JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
                JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
                JCI4=MAXBFN(W(LMXBAS),NMODE,M,1)
                I17=1
                IF(ICOUPC.GT.3)I17=17
                KTEMP1=JCI4*JCI4*I17
                KTEMP2=JCI3*JCI4*JCI3*JCI4*I17
                KTEMP3=JCI2*JCI3*JCI4*JCI2*JCI3*JCI4*I17
                KXK0=I17*JCI1*JCI1*IK2
                KXL0=I17*JCI2*JCI2*IL2
                KXN0=I17*JCI3*JCI3*IN2
                KXM0=I17*JCI4*JCI4*IM2
                CALL MEMO(3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
     1          KTEMP3,0,0,0,0)
                CALL MEMO(4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          0,0)
                IF(KKCONT.EQ.LLCONT.AND.KKCONT.EQ.NNCONT.AND.KKCONT.
     1          EQ.MMCONT)THEN
                  CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1            0,0,0,0,0,0)
                  CALL MEMO(1,LIP4,KIP4,0,0,0,0,0,0,0,0)
                  CALL GETBP4(W(LIP4),ISIZE4,NMODE,K,L,N,M,
     1            W(LMXBAS),NSIZE4,0)

!                  KXA4=NSIZE4*(NSIZE4+1)/2
                  if (mod(nsize4,2) .eq. 0) then
                     kxa4 = nsize4 / 2
                     kxa4 = kxa4 * (nsize4 + 1)
                  else
                     kxa4 = (nsize4 + 1) / 2
                     kxa4 = kxa4 * nsize4
                  end if

                  CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                  CALL DIAGZ(NSIZE4,NSIZE4,W(LXA4),NSIZE4,W(LXA4),
     1            ICI)
                  IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
                    CALL MEMO(1,LOV,KKXK,0,0,0,0,0,0,0,0)
                    CALL VCCI4C(NAMODE,NNMODE,K,L,N,M,KMODX,
     1              LMODX,NMODX,MMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),
     2              W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     3              IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LIP),ISIZMX,
     4              W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),ISIZE,
     5              W(LIP4),ISIZE4,W(LXA4),W(LXA4),NSIZE4,W(LXK),
     6              W(LTEMP),W(LOV),W(LWRK),KVAL,ISTART,IEND,
     7              W(LTEMP1),W(LTEMP2),W(LTEMP3),JCI1,JCI2,JCI3,JCI4,
     8              W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LV4),W(LV4),
     9              W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT),
     1              KEL,NS,
     1              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     2              ISIZXX,NVALCF,KEL21,NVSYM,I17,W(LXKAN),MAXQU,
     3              MAXPOW,W(LNP4),W(LCP4),W(LMP4),NTOT4,IENTMX(4),
     4              W(LINDK),W(LINDL),W(LINDN))
                    CALL MEMO(-1,LOV,KKXK,0,0,0,0,0,0,0,0)
                  ELSE
                    CALL VCCI4A(NAMODE,NNMODE,K,L,N,M,KMODX,
     1              LMODX,NMODX,MMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),
     2              W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     3              IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LIP),ISIZMX,
     4              W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),ISIZE,
     5              W(LIP4),ISIZE4,W(LXA4),W(LXA4),NSIZE4,W(LXK),
     6              W(LTEMP),W(LWRK),KVAL,ISTART,IEND,
     7              W(LTEMP1),W(LTEMP2),W(LTEMP3),JCI1,JCI2,JCI3,JCI4,
     8              W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LV4),W(LV4),
     9              W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT),
     1              KEL,NS,
     1              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     2              ISIZXX,NVALCF,KEL21,NVSYM,I17,W(LXKAN),MAXQU,
     3              MAXPOW,W(LNP4),W(LCP4),W(LMP4),NTOT4,IENTMX(4),
     4              W(LINDK),W(LINDL),W(LINDN))
                  END IF
                  CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1            0,0,0,0,0,0)
                  CALL MEMO(-2,LIP4,KIP4,LXA4,KXA4,0,0,0,0,0,
     1            0)
                  CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
     1            KTEMP3,0,0,0,0)
                  CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1            0,0)
                ELSE
                  ISIZC4=JCC4**4
                  NMAX1=MAX0(MAXSUM(1,1)+4,MAXSUM(2,1)+4,
     1            MAXSUM(3,1)+4,MAXSUM(4,1)+4)
                  NMAX2=MAX0(MAXSUM(1,2)+4,MAXSUM(2,2)+4,
     1            MAXSUM(3,2)+4,MAXSUM(4,2)+4)
                  NMAX=NMAX1+NMAX2
C**FIRST GET 4-D SIZE FOR 'NORMAL' (ISIZE4 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(ISIZC4,4,JCC4,0)
                  NSIZE4=ISIZC4
                  KCP4=ISIZC4*4
                  CALL MEMO(1,LIP4,KCP4,0,0,0,0,0,0,0,0)
                  CALL GETCP4(W(LIP4),ISIZC4,NMODE,K,L,N,M,
     1            W(LMXBAS),NSIZE4)

!                  KXA4=NSIZE4*(NSIZE4+1)/2
                  if (mod(nsize4,2) .eq. 0) then
                     kxa4 = nsize4 / 2
                     kxa4 = kxa4 * (nsize4 + 1)
                  else
                     kxa4 = (nsize4 + 1) / 2
                     kxa4 = kxa4 * nsize4
                  end if

                  CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                  CALL DIAGZ(NSIZE4,NSIZE4,W(LXA4),NSIZE4,W(LXA4),
     1            ICI)
                  NONC1=0
                  NONC2=0
                  IF(KKCONT.NE.LLCONT)THEN
                    KKKBAS=LONBAS
                    KKKCON=LLCONT
                    KKKXK=KLXK
                    KKKIP=KLIP
                    KKKVAL=LVAL
                  ELSE
                    IF(KKCONT.NE.NNCONT)THEN
                      KKKBAS=NONBAS
                      KKKCON=NNCONT
                      KKKXK=KNXK
                      KKKIP=KNIP
                      KKKVAL=NVAL
                    ELSE
                      KKKBAS=MONBAS
                      KKKCON=MMCONT
                      KKKXK=KMXK
                      KKKIP=KMIP
                      KKKVAL=MVAL
                    END IF
                  END IF
                  IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
                    CALL VCCI4D(NAMODE,NNMODE,K,L,N,M,KMODX,LMODX,
     1              NMODX,MMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),
     2              W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),IK3,IK2,IL3,
     3              IL2,IN3,IN2,IM3,IM2,W(LIP),ISIZMX,W(LKP+KONBAS),
     4              KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,KKCONT,KKKCON,
     5              W(LXA),ISIZE,W(LIP4),ISIZC4,W(LXA4),W(LXA4),NSIZE4,
     6              W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,
     7              KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     8              JCI1,JCI2,JCI3,JCI4,W(LXK0),W(LXL0),
     9              W(LXN0),W(LXM0),W(LV4),W(LV4),W(LC4),W(LC4),
     1              W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT),KEL,NS,
     2              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     3              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     4              ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I17,W(LXKAN),
     5              MAXQU,MAXPOW,W(LNP4),W(LCP4),W(LMP4),NTOT4,
     6              IENTMX(4),W(LINDK),W(LINDL),W(LINDN),MAXVAL)
                    CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
     1              KTEMP3,0,0,0,0)
                    CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,
     1              KXM0,0,0)
C***************
CCCC                CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                    CALL MEMO(1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                    NNC12=MAX0(NONC1,NONC2)
                    CALL MEMO(1,LXW,NNC12*MAXVAL,
     1              0,0,0,0,0,0,0,0)
                    CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1              0,0,0,0,0,0)
                    CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1              0,0,0,0,0,0)
                    CALL VCCI4D(NAMODE,NNMODE,K,L,N,M,KMODX,LMODX,
     1              NMODX,MMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),
     2              W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),IK3,IK2,IL3,
     3              IL2,IN3,IN2,IM3,IM2,W(LIP),ISIZMX,W(LKP+KONBAS),
     4              KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,KKCONT,KKKCON,
     5              W(LXA),ISIZE,W(LIP4),ISIZC4,W(LXA4),W(LXA4),NSIZE4,
     6              W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,
     7              KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     8              JCI1,JCI2,JCI3,JCI4,W(LXK0),W(LXL0),
     9              W(LXN0),W(LXM0),W(LV4),W(LV4),W(LC4),W(LC4),
     1              W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT),KEL,NS,
     2              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     3              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     4              ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I17,
     5              W(LXKAN),MAXQU,MAXPOW,W(LNP4),W(LCP4),W(LMP4),
     6              NTOT4,IENTMX(4),W(LINDK),W(LINDL),W(LINDN),MAXVAL)
                    CALL MEMO(-1,LXW,NNC12*MAXVAL,
     1              0,0,0,0,0,0,0,0)
                    CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1              0,0,0,0,0,0)
                    CALL MEMO(-1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
CCCC                CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                    CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1              0,0,0,0,0,0)
                  ELSE
                    CALL VCCI4B(NAMODE,NNMODE,K,L,N,M,KMODX,LMODX,
     1              NMODX,MMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),
     2              W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),IK3,IK2,IL3,
     3              IL2,IN3,IN2,IM3,IM2,W(LIP),ISIZMX,W(LKP+KONBAS),
     4              KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,KKCONT,KKKCON,
     5              W(LXA),ISIZE,W(LIP4),ISIZC4,W(LXA4),W(LXA4),NSIZE4,
     6              W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,
     7              KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     8              JCI1,JCI2,JCI3,JCI4,W(LXK0),W(LXL0),
     9              W(LXN0),W(LXM0),W(LV4),W(LV4),W(LC4),W(LC4),
     1              W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT),KEL,NS,
     2              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     3              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     4              ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I17,W(LXKAN),
     5              MAXQU,MAXPOW,W(LNP4),W(LCP4),W(LMP4),NTOT4,
     6              IENTMX(4),W(LINDK),W(LINDL),W(LINDN),MAXVAL,MAXQ4)
                    CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
     1              KTEMP3,0,0,0,0)
                    CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,
     1              KXM0,0,0)
C***************
CCCC                CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                    CALL MEMO(1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                    CALL MEMO(1,LXW,2*MAXVAL*MAXVAL*MAXVAL,
     1              0,0,0,0,0,0,0,0)
                    CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1              0,0,0,0,0,0)
                    CALL MEMO(2,LTEMP,MAXQ4*MAXQ4*MAXQ4*MAXQ4,
     1              LWRK,MAXVAL*MAXVAL,0,0,0,0,0,0)
                    CALL VCCI4B(NAMODE,NNMODE,K,L,N,M,KMODX,LMODX,
     1              NMODX,MMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),
     2              W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),IK3,IK2,IL3,
     3              IL2,IN3,IN2,IM3,IM2,W(LIP),ISIZMX,W(LKP+KONBAS),
     4              KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,KKCONT,KKKCON,
     5              W(LXA),ISIZE,W(LIP4),ISIZC4,W(LXA4),W(LXA4),NSIZE4,
     6              W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,
     7              KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     8              JCI1,JCI2,JCI3,JCI4,W(LXK0),W(LXL0),
     9              W(LXN0),W(LXM0),W(LV4),W(LV4),W(LC4),W(LC4),
     1              W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT),KEL,NS,
     2              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     3              W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     4              ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I17,
     5              W(LXKAN),MAXQU,MAXPOW,W(LNP4),W(LCP4),W(LMP4),
     6              NTOT4,IENTMX(4),W(LINDK),W(LINDL),W(LINDN),MAXVAL,
     7              MAXQ4)
                    CALL MEMO(-1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
CCCC                CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                    CALL MEMO(-1,LXW,2*MAXVAL*MAXVAL*MAXVAL,
     1              0,0,0,0,0,0,0,0)
                    CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1              0,0,0,0,0,0)
                    CALL MEMO(-2,LTEMP,MAXQ4*MAXQ4*MAXQ4*MAXQ4,
     1              LWRK,MAXVAL*MAXVAL,0,0,0,0,0,0)
                  END IF
                  CALL MEMO(-2,LIP4,KCP4,LXA4,KXA4,0,0,0,0,0,
     1            0)
                END IF
C               CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,
C    1          MAXVAL*MAXVAL,0,0,0,0,0,0)
C               CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
C    1          KTEMP3,0,0,0,0)
C               CALL MEMO(-4,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
C    1          0,0)
              ELSE
C*****************************   LIKE VCCI5 (K+L+N+M+NMODE)
                NREC(4)=NREC(4)+1
C**
C**GET BASIC INTEGRAL BASIS FOR FIVE MODES
                MAXXK=MAX0(KKXK,KLXK,KNXK,KMXK,KTXK)
                MAXVAL=MAX0(KVAL,LVAL,NVAL,MVAL,ITVAL)
                JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
                JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
                JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
                JCI4=MAXBFN(W(LMXBAS),NMODE,M,1)
                IF(JREACT.GT.0)THEN
                  JCIM=MAX0(NTAU(1)+1,NTAU(2)+1,NTAU(3)+1,NTAU(4)+1,
     1            NTAU(5)+1)
                ELSE
                  JCIM1=2*MAX0(NBAS(1,1)+1,NBAS(2,1)+1,NBAS(3,1)+1,
     1            NBAS(4,1)+1)-1
                  JCIM2=2*MAX0(NBAS(1,2)+1,NBAS(2,2)+1,NBAS(3,2)+1,
     1            NBAS(4,2)+1)-1
                  JCI=MAX0(JCIM1,JCIM2)
                END IF
                I36=1
                IF(ICOUPC.GT.3)I36=36
                KTEMP1=JCI1*JCI2*JCI3*JCI4*JCI1*JCI2*JCI3*JCI4*I36
CC              CALL MEMO(1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
                KTEMP2=JCI2*JCI3*JCI4*JCI2*JCI3*JCI4*I36
                KTEMP3=JCI3*JCI4*JCI3*JCI4*I36
                KTEMP4=JCI4*JCI4*I36
                KXK0=I36*JCI1*JCI1*IK2
                KXL0=I36*JCI2*JCI2*IL2
                KXN0=I36*JCI3*JCI3*IN2
                KXM0=I36*JCI4*JCI4*IM2
                KXP0=I36*JCIM*JCIM*IK2TAU
                CALL MEMO(3,LTEMP2,KTEMP2,LTEMP3,KTEMP3,LTEMP4,
     1          KTEMP4,0,0,0,0)
                CALL MEMO(5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          LXP0,KXP0)
                IF(KKCONT.EQ.LLCONT.AND.KKCONT.EQ.NNCONT.AND.KKCONT.
     1          EQ.MMCONT.AND.KKCONT.EQ.ITCONT)THEN
                  CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1            0,0,0,0,0,0)
                  IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1            M.GT.NONLIN)THEN
                    KPPP5=KLJP5
                    IND5=1
                  ELSE
                    KPPP5=KIP5
                    IND5=0
                  END IF
                  CALL MEMO(1,LIP5,KPPP5,0,0,0,0,0,0,0,0)
                  IF(K.GT.NONLIN.OR.L.GT.NONLIN.OR.N.GT.NONLIN.OR.
     1            M.GT.NONLIN)THEN
                    CALL GETBP4(W(LIP5),LSIZE5,NMODE,K,L,N,M,
     1              W(LMXBAS),NSIZE5,IND5)
                    MSIZE5=NSIZE5
                  ELSE
                    CALL GETBP5(W(LIP5),ISIZE5,NMODE,K,L,N,M,NNMODE,
     1              W(LMXBAS),NSIZE5,IND5)
                    MSIZE5=ISIZE5
                  END IF
CC                KXA5=NSIZE5*(NSIZE5+1)/2
                  KXA5=LANMAX*(LANMAX+1)/2
                  CALL MEMO(1,LXA4,KXA5,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
CC                CALL DIAGZ(NSIZE5,NSIZE5,W(LXA4),NSIZE5,W(LXA4),
CC   1            ICI)
                  CALL VCCV4A(NAMODE,NNMODE,K,L,N,M,KMODX,
     1            LMODX,NMODX,MMODX,ITMODE,W(LH+K2),W(LXQ+K3),
     2            W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     3            W(LXQ+M3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IL3,IL2,
     4            IN3,IN2,IM3,IM2,IK3TAU,IK2TAU,W(LIP),ISIZMX,
     4            W(LKP+KONBAS),KKXK,KKIP,KKCONT,W(LXA),ISIZE,
     5            W(LIP5),
     5            MSIZE5,W(LXA4),W(LXA4),NSIZE5,W(LXK),W(LTEMP),
     6            W(LWRK),KVAL,ISTART,IEND,W(LTEMP2),W(LTEMP3),
     7            W(LTEMP4),JCI1,JCI2,JCI3,JCI4,JCIM,W(LXK0),W(LXL0),
     8            W(LXN0),W(LXM0),W(LXP0),W(LV4),W(LV4),W(LC4),
     9            W(LC4),W(LEJK4),W(LEJK4),J21,KROT,W(LMODNT),KEL,NS,
     1            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     2            ISIZXX,NVALCF,KEL21,NVSYM,I36)
                  CALL MEMO(-2,LIP5,KPPP5,LXA4,KXA5,0,0,
     1            0,0,0,0)
                  CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1            0,0,0,0,0,0)
                ELSE
C**NOT DONE FOR LINEAR
                  ISIZC5=JCC5**5
                  NMAX1=MAX0(MAXSUM(1,1)+5,MAXSUM(2,1)+5,
     1            MAXSUM(3,1)+5,MAXSUM(4,1)+5,MAXSUM(5,1)+5)
                  NMAX2=MAX0(MAXSUM(1,2)+5,MAXSUM(2,2)+5,
     1            MAXSUM(3,2)+5,MAXSUM(4,2)+5,MAXSUM(5,2)+5)
                  NMAX=NMAX1+NMAX2
C**FIRST GET 4-D SIZE FOR 'NORMAL' (ISIZE4 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(ISIZC5,5,JCC5,0)
                  NSIZE5=ISIZC5
                  MSIZE5=NSIZE5
                  KCP5=ISIZC5*5
                  CALL MEMO(1,LIP5,KCP5,0,0,0,0,0,0,0,0)
                  CALL GETCP5(W(LIP5),ISIZC5,NMODE,K,L,N,M,NNMODE,
     1            W(LMXBAS),NSIZE5)
CC                KXA5=NSIZE5*(NSIZE5+1)/2
                  KXA5=LANMAX*(LANMAX+1)/2
                  CALL MEMO(1,LXA4,KXA5,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
CC                CALL DIAGZ(NSIZE5,NSIZE5,W(LXA4),NSIZE5,W(LXA4),
CC   1            ICI)
                  NONC1=0
                  NONC2=0
                  IF(KKCONT.NE.LLCONT)THEN
                    KKKBAS=LONBAS
                    KKKCON=LLCONT
                    KKKXK=KLXK
                    KKKIP=KLIP
                    KKKVAL=LVAL
                  ELSE
                    IF(KKCONT.NE.NNCONT)THEN
                      KKKBAS=NONBAS
                      KKKCON=NNCONT
                      KKKXK=KNXK
                      KKKIP=KNIP
                      KKKVAL=NVAL
                    ELSE
                      IF(KKCONT.NE.MMCONT)THEN
                        KKKBAS=MONBAS
                        KKKCON=MMCONT
                        KKKXK=KMXK
                        KKKIP=KMIP
                        KKKVAL=MVAL
                      ELSE
                        KKKBAS=ITBAS
                        KKKCON=ITCONT
                        KKKXK=KTXK
                        KKKIP=KTIP
                        KKKVAL=ITVAL
                      END IF
                    END IF
                  END IF
                  CALL VCCV4B(NAMODE,NNMODE,K,L,N,M,KMODX,
     1            LMODX,NMODX,MMODX,ITMODE,W(LH+K2),W(LXQ+K3),
     2            W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     3            W(LXQ+M3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IL3,IL2,
     4            IN3,IN2,IM3,IM2,IK3TAU,IK2TAU,W(LIP),ISIZMX,
     5            W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,
     6            KKCONT,
     6            KKKCON,W(LXA),ISIZE,W(LIP5),MSIZE5,W(LXA4),W(LXA4),
     7            NSIZE5,W(LXK),W(LXW),W(LTEMP),MAXXK,
     8            W(LWRK),KVAL,KKKVAL,ISTART,IEND,W(LTEMP2),
     9            W(LTEMP3),W(LTEMP4),JCI1,JCI2,JCI3,JCI4,JCIM,
     1            W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXP0),W(LV4),
     2            W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,KROT,
     3            W(LMODNT),KEL,NS,
     4            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     5            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     6            ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,I36,MAXVAL,
     7            MAXQ5)
C***************
                  CALL MEMO(1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(1,LXW,2*MAXVAL*MAXVAL*MAXVAL,
     1            0,0,0,0,0,0,0,0)
                  CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1            0,0,0,0,0,0)
                  CALL MEMO(2,LTEMP,MAXQ5*MAXQ5*MAXQ5*MAXQ5*MAXQ5,
     1            LWRK,MAXVAL*MAXVAL,0,0,0,0,0,0)
                  CALL VCCV4B(NAMODE,NNMODE,K,L,N,M,KMODX,
     1            LMODX,NMODX,MMODX,ITMODE,W(LH+K2),W(LXQ+K3),
     2            W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     3            W(LXQ+M3),W(LH+K2TAU),W(LXQ+K3TAU),IK3,IK2,IL3,IL2,
     4            IN3,IN2,IM3,IM2,IK3TAU,IK2TAU,W(LIP),ISIZMX,
     5            W(LKP+KONBAS),KKXK,KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,
     6            KKCONT,
     6            KKKCON,W(LXA),ISIZE,W(LIP5),MSIZE5,W(LXA4),
     7            W(LXA4),NSIZE5,W(LXK),W(LXW),W(LTEMP),MAXXK,
     8            W(LWRK),KVAL,KKKVAL,ISTART,IEND,W(LTEMP2),
     9            W(LTEMP3),W(LTEMP4),JCI1,JCI2,JCI3,JCI4,JCIM,
     1            W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXP0),W(LV4),
     2            W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,KROT,
     3            W(LMODNT),KEL,NS,
     4            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     5            W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     6            ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),I36,MAXVAL,
     7            MAXQ5)
                  CALL MEMO(-1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                  CALL MEMO(-1,LXW,2*MAXVAL*MAXVAL*MAXVAL,
     1            0,0,0,0,0,0,0,0)
                  CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1            0,0,0,0,0,0)
                  CALL MEMO(-2,LIP5,KCP5,LXA4,KXA5,0,0,
     1            0,0,0,0)
                  CALL MEMO(-2,LTEMP,MAXQ5*MAXQ5*MAXQ5*MAXQ5*MAXQ5,
     1            LWRK,MAXVAL*MAXVAL,0,0,0,0,0,0)
                END IF
CC              CALL MEMO(-1,LTEMP1,KTEMP1,0,0,0,0,0,0,0,0)
                CALL MEMO(-3,LTEMP2,KTEMP2,LTEMP3,KTEMP3,LTEMP4,
     1          KTEMP4,0,0,0,0)
                CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,KXM0,
     1          LXP0,KXP0)
C**
C*****************************   LIKE VCCI5 (K+L+N+M+NMODE)
              END IF
334   CONTINUE
C**INTRINSIC
              IF(ICOUPX.EQ.4)GO TO 5554
C**INTEGRATE OVER FIVE NORMAL COORDINATES
              J2=0
              J3=0
              DO J=1,M-1
                IF(J.EQ.1.AND.M.EQ.2.AND.N.EQ.3.AND.L.EQ.4.AND.
     1          K.EQ.5.AND.ITIM.EQ.0)THEN
                  ITIM5A=0
                  ITIM5B=0
                END IF
                CALL INTARR(W(LNBF),W(LMBF),J,IJ1,IJ2,IJ3)
C**DETERMINE WHICH CONTRACTION SCHEME INVOLVES 'N'
                DO JJ=1,ICONT(1)
                  IF(J.EQ.JCONT(1,JJ))THEN
                    JJCONT=1
                    JMODX=JJ
                    JONBAS=0
                    KJXK=ICSIZ1
                    KJIP=IPSIZ1
                    JVAL=NVAL1
                  END IF
                END DO
                DO JJ=1,ICONT(2)
                  IF(J.EQ.JCONT(2,JJ))THEN
                    JJCONT=2
                    JMODX=JJ
                    JONBAS=KLP
                    KJXK=ICSIZ2
                    KJIP=IPSIZ2
                    JVAL=NVAL2
                  END IF
                END DO
C**INTRINSIC
                IF(MOLINC.GT.0)THEN
                  CALL GROUP5(K,L,N,M,J,ICONT(1),JCONT,1,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(1)
                    JCOUPC=NCOUPC(1)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 3035
                  END IF
                  CALL GROUP5(K,L,N,M,J,ICONT(2),JCONT,2,IND)
                  IF(IND.NE.0)THEN
                    JCOUPL=NCOUPL(2)
                    JCOUPC=NCOUPC(2)
                    ICOUPL=IABS(JCOUPL)
                    ICOUPC=IABS(JCOUPC)
                    GO TO 3035
                  END IF
                  JCOUPL=JCOUPS
                  JCOUPC=JCOUCS
                  ICOUPL=ICOUPS
                  ICOUPC=ICOUCS
3035              CONTINUE
                END IF
                IF(ICOUPL.LT.5)GO TO 335
C**CORIOLIS AND POTENTIAL
                IF((JREACT.LE.0.AND.K.LE.NONLIN.AND.L.LE.NONLIN.AND.
     1          N.LE.NONLIN.AND.M.LE.NONLIN.AND.J.LE.NONLIN).OR.
     2          (NSMODE.EQ.NAMODE))THEN
C**GET BASIC INTEGRAL BASIS FOR FIVE MODES
                  NREC(5)=NREC(5)+1
                  MAXXK=MAX0(KKXK,KLXK,KNXK,KMXK,KJXK)
                  MAXVAL=MAX0(KVAL,LVAL,NVAL,MVAL,JVAL)
C                 CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
C    1            0,0,0,0,0,0)
                  JCI1=MAXBFN(W(LMXBAS),NMODE,K,1)
                  JCI2=MAXBFN(W(LMXBAS),NMODE,L,1)
                  JCI3=MAXBFN(W(LMXBAS),NMODE,N,1)
                  JCI4=MAXBFN(W(LMXBAS),NMODE,M,1)
                  JCI5=MAXBFN(W(LMXBAS),NMODE,J,1)
                  KTEMP1=JCI5*JCI5
                  KTEMP2=JCI4*JCI5*JCI4*JCI5
                  KTEMP3=JCI3*JCI4*JCI5*JCI3*JCI4*JCI5
C                 KTEMP4=JCI2*JCI3*JCI4*JCI5*JCI2*JCI3*JCI4*JCI5
C                 KTEMP5=JCI2*JCI3*JCI4*JCI5
C                 CALL MEMO(5,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
C    1            KTEMP3,LTEMP4,KTEMP4,LTEMP5,KTEMP5)
                  CALL MEMO(3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
     1            KTEMP3,0,0,0,0)
                  KXK0=JCI1*JCI1*IK2
                  KXL0=JCI2*JCI2*IL2
                  KXN0=JCI3*JCI3*IN2
                  KXM0=JCI4*JCI4*IM2
                  KXJ0=JCI5*JCI5*IJ2
                  CALL MEMO(5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,
     1            KXM0,LXJ0,KXJ0)
                  IF(KKCONT.EQ.LLCONT.AND.KKCONT.EQ.NNCONT.AND.
     1            KKCONT.EQ.MMCONT.AND.KKCONT.EQ.JJCONT)THEN
                    CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                    CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,MAXVAL*MAXVAL,
     1              0,0,0,0,0,0)
                    CALL MEMO(1,LIP5,KIP5,0,0,0,0,0,0,0,0)
                    CALL GETBP5(W(LIP5),ISIZE5,NMODE,K,L,N,M,J,
     1              W(LMXBAS),NSIZE5,0)
                    KXA5=NSIZE5*(NSIZE5+1)/2
                    CALL MEMO(1,LXA5,KXA5,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                    CALL DIAGZ(NSIZE5,NSIZE5,W(LXA5),NSIZE5,W(LXA5),
     1              ICI)
                    IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
                      CALL MEMO(1,LOV,KKXK,0,0,0,0,0,0,0,0)
                      CALL VCCI5C(NAMODE,NNMODE,K,L,N,M,J,KMODX,LMODX,
     1                NMODX,MMODX,JMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),
     2                W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     3                W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,
     4                IM2,IJ3,IJ2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,
     5                KKIP,KKCONT,W(LXA),ISIZE,W(LIP5),ISIZE5,W(LXA5),
     6                W(LXA5),NSIZE5,W(LXK),W(LTEMP),W(LOV),W(LWRK),
     7                KVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     8                W(LTEMP4),W(LTEMP5),JCI1,JCI2,JCI3,JCI4,JCI5,
     9                W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),W(LV4),
     1                W(LV4),W(LMODNT),KEL,NS,
     3                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     4                ISIZXX,NVALCF,KEL21,NVSYM,
     5                W(LXKAN),MAXQU,MAXPOW,W(LNP5),W(LCP5),W(LMP5),
     6                NTOT5,MAX5,W(LINDK),W(LINDL),W(LINDN),W(LINDM))
                      CALL MEMO(-1,LOV,KKXK,0,0,0,0,0,0,0,0)
                    ELSE
                      CALL VCCI5A(NAMODE,NNMODE,K,L,N,M,J,KMODX,LMODX,
     1                NMODX,MMODX,JMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),
     2                W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     3                W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,
     4                IM2,IJ3,IJ2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,
     5                KKIP,KKCONT,W(LXA),ISIZE,W(LIP5),ISIZE5,W(LXA5),
     6                W(LXA5),NSIZE5,W(LXK),W(LTEMP),W(LWRK),KVAL,
     7                ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     8                W(LTEMP4),W(LTEMP5),JCI1,JCI2,JCI3,JCI4,JCI5,
     9                W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),W(LV4),
     1                W(LV4),W(LMODNT),KEL,NS,
     3                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     4                ISIZXX,NVALCF,KEL21,NVSYM,
     5                W(LXKAN),MAXQU,MAXPOW,W(LNP5),W(LCP5),W(LMP5),
     6                NTOT5,MAX5,W(LINDK),W(LINDL),W(LINDN),W(LINDM))
                    END IF
                    CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                    CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,
     1              MAXVAL*MAXVAL,0,0,0,0,0,0)
                    CALL MEMO(-2,LIP5,KIP5,LXA5,KXA5,0,0,0,0,0,0)
                    CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
     1              KTEMP3,0,0,0,0)
                    CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,
     1              KXM0,LXJ0,KXJ0)
                  ELSE
                    ISIZC5=JCC5**5
                    NMAX1=MAX0(MAXSUM(1,1)+5,MAXSUM(2,1)+5,
     1              MAXSUM(3,1)+5,MAXSUM(4,1)+5,MAXSUM(5,1)+5)
                    NMAX2=MAX0(MAXSUM(1,2)+5,MAXSUM(2,2)+5,
     1              MAXSUM(3,2)+5,MAXSUM(4,2)+5,MAXSUM(5,2)+5)
                    NMAX=NMAX1+NMAX2
C**FIRST GET 5-D SIZE FOR 'NORMAL' (ISIZE5 MODIFIED IN GETSZ FOR NMAX)
                    CALL GETSZ(ISIZC5,5,JCC5,0)
                    NSIZE5=ISIZC5
                    KCP5=ISIZC5*5
                    CALL MEMO(1,LIP5,KCP5,0,0,0,0,0,0,0,0)
                    CALL GETCP5(W(LIP5),ISIZC5,NMODE,K,L,N,M,J,
     1              W(LMXBAS),NSIZE5)
                    KXA5=NSIZE5*(NSIZE5+1)/2
                    CALL MEMO(1,LXA5,KXA5,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
                    CALL DIAGZ(NSIZE5,NSIZE5,W(LXA5),NSIZE5,W(LXA5),
     1              ICI)
                    NONC1=0
                    NONC2=0
                    IF(KKCONT.NE.LLCONT)THEN
                      KKKBAS=LONBAS
                      KKKCON=LLCONT
                      KKKXK=KLXK
                      KKKIP=KLIP
                      KKKVAL=LVAL
                    ELSE
                      IF(KKCONT.NE.NNCONT)THEN
                        KKKBAS=NONBAS
                        KKKCON=NNCONT
                        KKKXK=KNXK
                        KKKIP=KNIP
                        KKKVAL=NVAL
                      ELSE
                        IF(KKCONT.NE.MMCONT)THEN
                          KKKBAS=MONBAS
                          KKKCON=MMCONT
                          KKKXK=KMXK
                          KKKIP=KMIP
                          KKKVAL=MVAL
                        ELSE
                          KKKBAS=JONBAS
                          KKKCON=JJCONT
                          KKKXK=KJXK
                          KKKIP=KJIP
                          KKKVAL=JVAL
                        END IF
                      END IF
                    END IF
                    IF(ICCONT(1).LT.0.OR.ICCONT(2).LT.0)THEN
                      CALL VCCI5D(NAMODE,NNMODE,K,L,N,M,J,KMODX,LMODX,
     1                NMODX,MMODX,JMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),
     2                W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     3                W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,
     4                IM2,IJ3,IJ2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,
     5                KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,KKCONT,KKKCON,
     6                W(LXA),ISIZE,W(LIP5),ISIZC5,W(LXA5),W(LXA5),
     7                NSIZE5,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,
     8                KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     9                W(LTEMP4),W(LTEMP5),JCI1,JCI2,JCI3,JCI4,JCI5,
     1                W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),W(LV4),
     2                W(LV4),W(LMODNT),KEL,NS,
     4                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     5                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     6                ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,
     7                W(LXKAN),MAXQU,MAXPOW,W(LNP5),W(LCP5),W(LMP5),
     8                NTOT5,MAX5,W(LINDK),W(LINDL),W(LINDN),W(LINDM),
     9                MAXVAL)
                      CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
     1                KTEMP3,0,0,0,0)
                      CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,
     1                KXM0,LXJ0,KXJ0)
C***************
CCCC                  CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                      CALL MEMO(1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                      NNC12=MAX0(NONC1,NONC2)
                      CALL MEMO(1,LXW,NNC12*MAXVAL,
     1                0,0,0,0,0,0,0,0)
                      CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1                0,0,0,0,0,0)
                      CALL MEMO(2,LTEMP,MAXXK*MAXVAL,LWRK,
     1                MAXVAL*MAXVAL,0,0,0,0,0,0)
                      CALL VCCI5D(NAMODE,NNMODE,K,L,N,M,J,KMODX,LMODX,
     1                NMODX,MMODX,JMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),
     2                W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     3                W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,
     4                IM2,IJ3,IJ2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,
     5                KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,KKCONT,KKKCON,
     6                W(LXA),ISIZE,W(LIP5),ISIZC5,W(LXA5),W(LXA5),
     7                NSIZE5,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,
     8                KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     9                W(LTEMP4),W(LTEMP5),JCI1,JCI2,JCI3,JCI4,JCI5,
     1                W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),W(LV4),
     2                W(LV4),W(LMODNT),KEL,NS,
     4                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     5                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     6                ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),
     7                W(LXKAN),MAXQU,MAXPOW,W(LNP5),W(LCP5),W(LMP5),
     8                NTOT5,MAX5,W(LINDK),W(LINDL),W(LINDN),W(LINDM),
     9                MAXVAL)
                      CALL MEMO(-1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
CCCC                  CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                      CALL MEMO(-1,LXW,NNC12*MAXVAL,
     1                0,0,0,0,0,0,0,0)
                      CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1                0,0,0,0,0,0)
                      CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,
     1                MAXVAL*MAXVAL,0,0,0,0,0,0)
                    ELSE
                      CALL VCCI5B(NAMODE,NNMODE,K,L,N,M,J,KMODX,LMODX,
     1                NMODX,MMODX,JMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),
     2                W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     3                W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,
     4                IM2,IJ3,IJ2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,
     5                KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,KKCONT,KKKCON,
     6                W(LXA),ISIZE,W(LIP5),ISIZC5,W(LXA5),W(LXA5),
     7                NSIZE5,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,
     8                KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     9                W(LTEMP4),W(LTEMP5),JCI1,JCI2,JCI3,JCI4,JCI5,
     1                W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),W(LV4),
     2                W(LV4),W(LMODNT),KEL,NS,
     4                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     5                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     6                ISIZXX,NVALCF,KEL21,NVSYM,0,IDUM,IDUM,
     7                W(LXKAN),MAXQU,MAXPOW,W(LNP5),W(LCP5),W(LMP5),
     8                NTOT5,MAX5,W(LINDK),W(LINDL),W(LINDN),W(LINDM),
     9                MAXVAL,MAXQ5)
                      CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
     1                KTEMP3,0,0,0,0)
                      CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,
     1                KXM0,LXJ0,KXJ0)
C***************
CCCC                  CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                      CALL MEMO(1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
                      CALL MEMO(1,LXW,2*MAXVAL*MAXVAL*MAXVAL,
     1                0,0,0,0,0,0,0,0)
                      CALL MEMO(2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1                0,0,0,0,0,0)
                      CALL MEMO(2,LTEMP,MAXQ5*MAXQ5*MAXQ5*MAXQ5*MAXQ5,
     1                LWRK,MAXVAL*MAXVAL,0,0,0,0,0,0)
                      CALL VCCI5B(NAMODE,NNMODE,K,L,N,M,J,KMODX,LMODX,
     1                NMODX,MMODX,JMODX,W(LH+K2),W(LXQ+K3),W(LH+L2),
     2                W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     3                W(LH+J2),W(LXQ+J3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,
     4                IM2,IJ3,IJ2,W(LIP),ISIZMX,W(LKP+KONBAS),KKXK,
     5                KKIP,W(LKP+KKKBAS),KKKXK,KKKIP,KKCONT,KKKCON,
     6                W(LXA),ISIZE,W(LIP5),ISIZC5,W(LXA5),W(LXA5),
     7                NSIZE5,W(LXK),W(LXW),W(LTEMP),MAXXK,W(LWRK),KVAL,
     8                KKKVAL,ISTART,IEND,W(LTEMP1),W(LTEMP2),W(LTEMP3),
     9                W(LTEMP4),W(LTEMP5),JCI1,JCI2,JCI3,JCI4,JCI5,
     1                W(LXK0),W(LXL0),W(LXN0),W(LXM0),W(LXJ0),W(LV4),
     2                W(LV4),W(LMODNT),KEL,NS,
     4                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKCONT-1)),
     5                W(LCFS+ISIZXX*NVALCF*KEL21*NVSYM*(KKKCON-1)),
     6                ISIZXX,NVALCF,KEL21,NVSYM,1,W(LHR),W(LHL),
     7                W(LXKAN),MAXQU,MAXPOW,W(LNP5),W(LCP5),W(LMP5),
     8                NTOT5,MAX5,W(LINDK),W(LINDL),W(LINDN),W(LINDM),
     9                MAXVAL,MAXQ5)
                      CALL MEMO(-1,LXK,MAXXK*MAXXK,0,0,0,0,0,0,0,0)
CCCC                  CALL MEMO(1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
                      CALL MEMO(-1,LXW,2*MAXVAL*MAXVAL*MAXVAL,
     1                0,0,0,0,0,0,0,0)
                      CALL MEMO(-2,LHR,4*(NONC1+1),LHL,4*(NONC2+1),
     1                0,0,0,0,0,0)
                      CALL MEMO(-2,LTEMP,MAXQ5*MAXQ5*MAXQ5*MAXQ5*MAXQ5,
     1                LWRK,MAXVAL*MAXVAL,0,0,0,0,0,0)
                    END IF
                    CALL MEMO(-2,LIP5,KCP5,LXA5,KXA5,0,0,0,0,0,0)
                  END IF
C                 CALL MEMO(-2,LTEMP,MAXXK*MAXVAL,LWRK,
C    1            MAXVAL*MAXVAL,0,0,0,0,0,0)
C                 CALL MEMO(-5,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
C    1            KTEMP3,LTEMP4,KTEMP4,LTEMP5,KTEMP5)
C                 CALL MEMO(-3,LTEMP1,KTEMP1,LTEMP2,KTEMP2,LTEMP3,
C    1            KTEMP3,0,0,0,0)
C                 CALL MEMO(-5,LXK0,KXK0,LXL0,KXL0,LXN0,KXN0,LXM0,
C    1            KXM0,LXJ0,KXJ0)
                ELSE
C**NO REACTION PATH
                END IF
335   CONTINUE
                J2=J2+IJ1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
                IF(J.GT.NONLIN)J2=J2+IJ1
                J3=J3+IJ2
              END DO
5554   CONTINUE
              M2=M2+IM1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
              IF(M.GT.NONLIN)M2=M2+IM1
              M3=M3+IM2
            END DO
5553   CONTINUE
            N2=N2+IN1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
            IF(N.GT.NONLIN)N2=N2+IN1
            N3=N3+IN2
          END DO
5552   CONTINUE
          L2=L2+IL1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
          IF(L.GT.NONLIN)L2=L2+IL1
          L3=L3+IL2
        END DO
CCCC    CALL MEMO(-1,LXK,KKXK*KKXK,0,0,0,0,0,0,0,0)
5551   CONTINUE
        K2=K2+IK1
C**RESERVE TWICE STORAGE FOR LINEAR FUNCTIONS
        IF(K.GT.NONLIN)K2=K2+IK1
        K3=K3+IK2
      END DO
      DO I=1,4
        NREC(I)=-1
      END DO
      ISKIP=0
      IF(LANCZ)THEN
        WRITE(IOUT,*)'Calculating LANOUT'

C**TEMPORARY
C     CALL PRXA(W(LXA),ISIZE)
C**TEMPORARY

        CALL TIMIT(1)
        CALL FLUSH(IOUT)
        CALL LANOUT(W(LXA),W(LWK),W(LWK),W(LYK),W(LZK),W(LZK),ISIZE,
     1  ISTART,IEND,W(LANCNT),W(LANK))
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
        IF(IEND.NE.ISIZE)THEN
          ISTART=IEND+1
          ISKIP=1
C**TEMPORARY-LAN
          ISYMFN=ISYMFN+1
          ICONDP=1
C**TEMPORARY-LAN
        END IF
      END IF
      IF(ISKIP.NE.0)GO TO 5555
      IF(LANCZ)CALL MEMO(-1,LZK,ISIZE,0,0,0,0,0,0,0,0)
C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
      REWIND 60
      WRITE(IOUT,260)
      CALL FLUSH(IOUT)
      NVAL=ISCFCI
      NVEC=ISCFCI
      IF(ISIZE.LT.ISCFCI)THEN
        NVAL=ISIZE
        NVEC=ISIZE
      END IF
      KXK=ISIZE*NVAL
      KWK=ISIZE
      KSUP4=5*ISIZE
      KEVAL=NVAL
      IF(MS.EQ.1)EVL=EVLJ0
      CALL TIMIT(1)
      IF(LANCZ)THEN
        IF(ISIZE.LT.LANCYC)THEN
          LANCYC=ISIZE
        END IF
        LGIV=(IGIV.NE.0)
        KXKL=LANCYC*NCYCLE*LANCYC*NCYCLE
        IF(LGIV)KXKL=LANCYC*NCYCLE*LANCYC
        CALL MEMO(3,LXK,ISIZE*NVAL,LVK,ISIZE*LANCYC,
     1  LZK,ISIZE*LANCYC,0,0,0,0)
        CALL MEMO(4,LSK,ISIZE*LANCYC,LEN,NVAL,LEL,NVAL,LNV,NVAL,
     1  0,0)
        CALL MEMO(4,LWRK,LANCYC*NCYCLE,LSUP4,5*LANCYC*NCYCLE,
     1  LEVAL,LANCYC*NCYCLE,LXKL,KXKL,0,0)
        KXA=LANCYC*NCYCLE*(LANCYC*NCYCLE+1)/2
        IF(LGIV)CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
        CALL MEMO(1,LXKLC,KXA,0,0,0,0,0,0,0,0)
        WRITE(IOUT,*)'Calculating LANCZOS'
        CALL TIMIT(1)
        CALL FLUSH(IOUT)
        IF(LGIV)THEN
          CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1    ISIZE,NVAL,W(LIP),ISIZMX,NNMODE,W(LKP),W(LKP+KLP),W(LCASS),
     2    KSCFCI,KEL21,
     2    NS,W(LEN),W(LEL),W(LNV),W(LXA),W(LXKLC),W(LXKL),W(LWRK),
     3    W(LSUP4),W(LEVAL),NCYCLE,NVSYM,W(LANCNT),W(LANK),LANCYC)
        ELSE
          CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1    ISIZE,NVAL,W(LIP),ISIZMX,NNMODE,W(LKP),W(LKP+KLP),W(LCASS),
     2    KSCFCI,KEL21,
     2    NS,W(LEN),W(LEL),W(LNV),W(LXKL),W(LXKLC),W(LXKL),W(LWRK),
     3    W(LSUP4),W(LEVAL),NCYCLE,NVSYM,W(LANCNT),W(LANK),LANCYC)
        END IF
        CALL TIMIT(3)
        CALL FLUSH(IOUT)
        IF(JTHIS.NE.0)THEN
C**SAVE VCI ENERGIES AND COEFFICIENTS FOR Ka>0
C         CALL ENERCF(W(LCFS),ISIZXX,W(LEVCI),ISIZE,NVALCF,KEL21,
C    1    KEL,NS,W(LXK),W(LXK),W(LEL))
        END IF
C**DUMP CI COEFFICIENTS TO (INP4)
        IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
          IF(JTHIS.EQ.0)THEN
            IDDDP=JDUMP(NS)
            IF(NVAL.LT.JDUMP(NS))IDDDP=NVAL
            CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,LDUMP,IDDDP,
     1      W(LEL),NMODE)
          ELSE
            CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,NVAL,NVAL,
     1      W(LEL),NMODE)
          END IF
        END IF
        CALL MEMO(-5,LXK,ISIZE*NVAL,LVK,ISIZE*LANCYC,
     1  LZK,ISIZE*LANCYC,LWK,ISIZE,LYK,ISIZE)
        CALL MEMO(-5,LSK,ISIZE*LANCYC,LEN,NVAL,LEL,NVAL,LNV,NVAL,
     1  LXKL,KXKL)
        CALL MEMO(-4,LXKLC,KXA,LWRK,LANCYC*NCYCLE,
     1  LSUP4,5*LANCYC*NCYCLE,LEVAL,LANCYC*NCYCLE,0,0)
        IF(LGIV)CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LANCNT,ISIZE,LANK,ISIZE,0,0,0,0,0,0)
        LGIV=.FALSE.
      ELSE
        LGIV=.TRUE.
        CALL MEMO(4,LXK,KXK,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL,0,0)
        WRITE(IOUT,*)'Calculating DIAG'

C**TEMPORARY
C     CALL PRXA(W(LXA),ISIZE)
C**TEMPORARY

        CALL FLUSH(IOUT)
        CALL DIAG(W(LXA),W(LXK),ISIZE,ISIZE,NVSYM,W(LSUP4),W(LEVAL),
     1  W(LWK),NVAL,NVEC,W(LIP),ISIZMX,NNMODE,W(LKP),W(LKP+KLP),
     2  W(LCASS),KSCFCI,KEL21,NS)
        IF(JTHIS.NE.0)THEN
C**SAVE VCI ENERGIES AND COEFFICIENTS FOR Ka>0
C         CALL ENERCF(W(LCFS),ISIZXX,W(LEVCI),ISIZE,NVALCF,KEL21,KEL,
C    1    NS,W(LXK),W(LXK),W(LWK))
        END IF
C**DUMP CI COEFFICIENTS TO (INP4)
        IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
          IF(JTHIS.EQ.0)THEN
            IDDDP=JDUMP(NS)
            IF(NVAL.LT.JDUMP(NS))IDDDP=NVAL
            CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,LDUMP,IDDDP,
     1      W(LWK),NMODE)
          ELSE
            CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,NVAL,NVAL,
     1      W(LWK),NMODE)
          END IF
        END IF
        CALL MEMO(-4,LXK,KXK,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL,0,0)
        LGIV=.FALSE.
        CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
      END IF
      CALL TIMIT(3)
      CALL FLUSH(IOUT)
8000  CONTINUE
      ICONDP=1
      CALL MEMO(-1,LIP,KIP,0,0,0,0,0,0,0,0)
      GO TO 5501
5502  CONTINUE
      CALL MEMO(-1,LJP,KJP,0,0,0,0,0,0,0,0)
5501  CONTINUE
5500  CONTINUE
      WRITE(IOUT,240)EVL
      CALL FLUSH(IOUT)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE PRJP(ISTAT,NSTAT,NSMODE,ICI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ISTAT(NSTAT,NSMODE)
      COMMON/FILASS/IOUT,INP
      DO I=1,ICI
      WRITE(IOUT,*)(ISTAT(I,J),J=1,NSMODE)
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE MAXIN(MAXBAS,NMODE,J,ICI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NMODE,6)
C**TEMPORARY (DIMENSIONS)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
C************************
      COMMON/BASIS/NBAS(6,2),MAXSUM(6,2)
      COMMON/TBASIS/NTBAS(6,2),NTAU(6)
      COMMON/FILASS/IOUT,INP
105   FORMAT(1X,I1,'-MODE STATES: ',/,1X,50I3)
C**SET DEFAULTS OF ZERO
      DO K=1,6
        DO I=1,NMODE
          MAXBAS(I,K)=1
        END DO
      END DO
      DO K=1,J
        READ(INP,*)(MAXBAS(I,K),I=1,NMODE)
        IF(ICI.LT.0)WRITE(IOUT,105)K,(MAXBAS(I,K),I=1,NMODE)
        IF(IABS(ICI).NE.100)THEN
          DO I=1,NMODE
            IF(MAXBAS(I,K).LT.0)MAXBAS(I,K)=0
            DO L=1,NCONT
              NUM=ICONT(L)
              DO M=1,NUM
                IF(I.EQ.JCONT(L,M))THEN
                  IF(MAXBAS(I,K).GT.NBAS(K,L))NBAS(K,L)=MAXBAS(I,K)
                END IF
              END DO
            END DO
C**CONVERT TO NO. FUNCTIONS
            MAXBAS(I,K)=MAXBAS(I,K)+1
          END DO
        END IF
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE TMAXIN(MAXBAS,NMODE,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NMODE,6)
C**TEMPORARY (DIMENSIONS)
      COMMON/CONTS/NCONT,ICONT(2),JCONT(2,100),NCOUPL(2),NCOUPC(2)
C************************
      COMMON/TBASIS/NTBAS(6,2),NTAU(6)
      COMMON/FILASS/IOUT,INP
      DO K=1,J
C**DON'T INCLUDE TORSION
        DO I=1,NMODE-1
          DO L=1,NCONT
            NUM=ICONT(L)
            DO M=1,NUM
              IF(I.EQ.JCONT(L,M))THEN
C**MAXBAS EXPRESSED AS FUNCTION, NTBAS EXPRESSED AS QUANTA
                IF(MAXBAS(I,K)-1.GT.NTBAS(K,L))NTBAS(K,L)=MAXBAS(I,K)-1
              END IF
            END DO
          END DO
        END DO
        NTAU(K)=MAXBAS(NMODE,K)-1
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE DEFMAX(MAXBAS,NSMODE,JCI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NSMODE,6)
      DO I=1,6
        DO J=1,NSMODE
          MAXBAS(J,I)=JCI
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      INTEGER FUNCTION MAXBFN(MAXBAS,NMODE,K,J)
      DIMENSION MAXBAS(NMODE,J)
      MAXBFN=MAXBAS(K,J)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE INMDIN(MODINT,NMODE,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MODINT(NMODE)
      COMMON/FILASS/IOUT,INP
100   FORMAT(//,1X,'MODE-SYMMETRY-INTEGRATION PARAMETERS',/,
     11X,39I2,/)
101   FORMAT(//,1X,'EFFECTIVE POTENTIAL PARAMETERS',/,
     11X,39I2,/)
102   FORMAT(1X,39I2,/)
      READ(INP,*)(MODINT(M),M=1,NMODE)
      IF(IND.GT.0)THEN
        WRITE(IOUT,101)(MODINT(M),M=1,NMODE)
      END IF
      IF(IND.LT.0)THEN
        WRITE(IOUT,102)(MODINT(M),M=1,NMODE)
      END IF
      IF(IND.NE.0)RETURN
      DO M=1,NMODE
        IF(MODINT(M).LT.1)MODINT(M)=1
C       IF(MODINT(M).GT.2)MODINT(M)=1
      END DO
      WRITE(IOUT,100)(MODINT(M),M=1,NMODE)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE PRMAT(XA,ISIZE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(1)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      K=0
      DO I=1,ISIZE
        WRITE(IOUT,*)(XA(J+K)*WAVENM,J=1,I)
        K=K+I
      END DO
      RETURN
      END
