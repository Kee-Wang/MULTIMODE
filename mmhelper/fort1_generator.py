
def generate_fort1(filename, p):
    """

    :param filename: The name of output file. By default is `fort.1`
    :param p: An nested array of parameters. p[i] is the (i-1)th line of parameters in orginal `fort.1` input.
    :return: None

    TO-DO:
    1. Input error control
    """
    "Conversion Table"
    TITLE = p[0]
    [NATOM,NSTAT,CONV,ICOUPL,ICOUPC,ISCFCI,IWHICH,IDISC,NROTTR,JMAX,MCHECK,INORM] = p[1]
    [ICI, NMAX, CUT, EVLJ0, NVALR, KSTEP, IPRINT, MATSIZ, IREACT, MOLPRO, MOLINC] = p[2]
    XTANPM = p[3]
    XMODQ = p[4]
    [NCONT, NVAL1, NVAL2] = p[5]
    [ICONT, JCONT] = p[6]
    MAXSUM = p[7]
    MAXBAS = p[8]
    MODINT = p[9]
    MEFF = p[10]
    [NSYM, ISYM] = p[13]
    MXDIP = p[15]
    MX = p[16]
    MXROT = p[17]
    [LDUMP, MDUMP, IDUMP] = p[19]
    [NBF, MBF, NVF, MVF] = p[20]
    COORD = p[26]
    MAXCHK = p[28]
    OVF = p[29]


    "Generating fort.1 file"

    f=open(filename,'w')
    f.write('C**TITLE \n')
    f.write('{:s}\n'.format(p[0]))
    f.write('C**NATOM,NSTAT,CONV,ICOUPL,ICOUPC,ISCFCI,IWHICH,IDISC,NROTTR,JMAX,MCHECK,INORM \n')
    f.write('   {:^6d}{:^5d}{:^5s}{:^7d}{:^7d}{:^7d}{:^7d}{:^6d}{:^8d}{:^5d}{:^7d}{:^5d}\n'.format(*p[1]))
    f.write('C**ICI NMAX  CUT     EVLJ0   NVALR KSTEP IPRINT MATSIZ IREACT MOLPRO MOLINC \n')
    f.write('  {:^5d}{:<6d}{:<8s}{:^5.2f}{:^8d}{:^7d}{:^7d}{:^7d}{:^7d}{:^8d}{:^6d}\n'.format(*p[2]))
    f.write('C**XTANPM (DATA IGNORED IF MOLPRO.EQ.0 AND IWHICH.GE.0)\n')
    if MOLPRO == 0 and IWHICH >= 0:
        pass
    else:
        f.write(p[3]) #Need to work on the format of this!
    f.write('C**XMODQ \n')
    f.write('   {:^5.1f}{:^5.1f}{:^5.1f}\n'.format(*p[4]))
    f.write('C**NCONT NVAL1 NVAL2 (DATA IGNORED IF NMAX.GE.0 OR ISCFCI.LE.0) \n')
    f.write('   {:^5d} {:^5d} {:^5d} \n'.format(*p[5]))
    f.write('C**ICONT JCONT (DATA IGNORED IF NMAX.GE.0 OR ISCFCI.LE.0) \n')
    f.write('   {:<5d} {:<30s} \n'.format(*p[6])) #Here JCONT used string
    f.write('C**MAXSUM (DATA IGNORED IF NMAX.GE.0 OR ISCFCI.LE.0) \n')
    if NMAX >= 0 or ISCFCI <= 0:
        pass
    else:
        f.write(('{:5d}' * len(MAXSUM) + '\n').format(*MAXSUM))
    f.write('C**MAXBAS (DATA IGNORED IF NMAX.GE.0 OR ISCFCI.LE.0) \n')
    if NMAX >= 0 or ISCFCI <= 0:
        pass
    else:
        for nbasis in MAXBAS: #Iterate MAXBAS
            f.write(('{:5d}' * len(nbasis) + '\n').format(*nbasis))
    f.write('C**MODINT\n')
    f.write(('{:5d}' * len(MODINT) + '\n').format(*MODINT))
    f.write('C**MEFF\n')
    for line in MEFF:
        f.write(('{:5d}' * len(line) + '\n').format(*line))
    f.write('C**NCYCLE,TOLLAN,LANMAX,LGIV,LAN20,LANCUT,TOLCUT,LANCYC,TOLMAT, LANINC (LANCZOS)\n')
    f.write('   {:^7d}{:^7s}{:^7d}{:^5d}{:^5d}{:^7d}{:^7s}{:^7d}{:^8s}{:^8d}\n'.format(*p[11]))
    f.write('C**NRSYM NWSYM ISYMPG\n')
    f.write('  {:^6d}{:^6d}{:^6d}\n'.format(*p[12]))
    f.write('C**NSYM ISYM (ORDER A1,B2,B1,A2) \n')
    f.write(('  {:^6d}'+'{:3d}' * len(ISYM) + '\n').format(NSYM, *ISYM))
    f.write('C**MVSYM   MWSYM \n')
    f.write('   {:^5d}   {:^5d}\n'.format(*p[14]))
    f.write('C**MXDIP (I.R. of mu(X), mu(Y), mu(Z))\n')
    f.write(('{:5d}' * len(MXDIP) + '\n').format(*MXDIP))
    f.write('C**MX (symmetry axes x,y,z corresponding to principle axes X,Y,Z)\n')
    f.write(('{:5d}' * len(MX) + '\n').format(*MX))
    f.write('C**MXROT (I.R. of Rx, Ry, Rz) \n')
    f.write(('{:5d}'*len(MXROT)+'\n').format(*MXROT))
    f.write('C**ISYMT (I.R. of Sin(M.tau/2)\n')
    f.write('{:5d}\n'.format(p[18]))
    f.write('C**LDUMP MDUMP IDUMP (DATA IGNORED IF ICI.GE.0)\n')
    if ICI >= 0:
        pass
    else:
        f.write((('  '+'{:^5d}'*2+'  '+'{:>3d}'*len(IDUMP)+'\n').format(LDUMP, MDUMP,*IDUMP)))
    f.write('C**NBF,MBF,NVF,MVF \n')
    for i in range(len(NBF)):
        f.write(('{:5d}'*4+'\n').format(NBF[i],MBF[i],NVF[i],MVF[i]))
    f.write('C**NUCLEAR SYMBOLS\n')
    f.write((' '*4+'{:4s}'*3+'\n').format(*p[21]))
    f.write('C**SCF STATE DEFINITIONS (SEE ORIGINAL DEFINITION BY JELSKI) \n')
    f.write('C**NPOT (ONLY IF IWHICH=0) \n')
    if IWHICH is 0:
        pass
    else:
        pass
    f.write('C**IPOT,JPOT,CPOT (ONLY IF IWHICH=0 OR INORM=0)\n')
    if IWHICH is 0 or INORM is 0:
        pass
    else:
        pass
    f.write('C**NUCLEAR MASSES (XM(I),I=1,NATOM)\n')
    f.write(('{:13.8f}'*3 + '\n').format(*NMASS))
    f.write('C**EQUILIBRIUM CARTESIAN COORDINATES ((X0(I,J),J=x,y,z),I=1,NATOM)\n')
    f.write(COORD)
    f.write('\n')
    f.write('C**NORMAL COORDINATE DISPLACEMENTS (((XL(I,J,K),K=x,y,z),I=1,NATOM),J=1,NMODE)\n')
    f.write('C**MAXCHK \n')
    f.write(('{:5d}' * len(MAXCHK) + '\n').format(*MAXCHK))
    f.write('C**OUTPUT VIBRATION FILE\n')
    f.write('{:<30s}'.format(p[29]))









    f.close
    return



TITLE = 'water example'
NATOM = 3
NSTAT = -1
CONV = '1.D-2'
ICOUPL = 3
ICOUPC = 3
ISCFCI = 200
IWHICH = 1
IDISC = 0
NROTTR = 0
JMAX = 0
MCHECK = 1
INORM = 1
ICI=-1
NMAX=-3
CUT = '1.0D+4'
EVLJ0 = 0.0
NVALR = 20
KSTEP =0
IPRINT = -3
MATSIZ = 0
IREACT = 0
MOLPRO = 0
MOLINC = 0
XTANPM = None
XMODQ = [1.0, 1.0, 1.0]
NCONT = 1
NVAL1 = 0
NVAL2 = 0
ICONT = -3
JCONT = '1 2 3   3 3'
MAXSUM = [12, 12, 12]
MAXBAS = [[11, 12, 12],[12,12,12],[13,12,12]]
MODINT = [1, 1, 1]
MEFF = [[0,0,0],[0,0,0]]
NCYCLE=0
TOLLAN = '1.D-3'
LANMAX = 2000
LGIV = 1
LAN20 = 500
LANCUT = 0
TOLCUT = '0.D0'
LANCYC = 200
TOLMAT = '1.D-10'
LANINC = 4
NRSYM = 1
NWSYM = 1
ISYMPG = 0
NSYM = 3
ISYM = [1,2,3]
MVSYM = 1
MWSYM = 1
MXDIP = [1,1,1]
MX = [1,2,3]
MXROT = [1,1,1]
ISYMT = 1
LDUMP = 0
MDUMP = 0
IDUMP = [0, 0]
NBF = [12, 12, 12]
MBF = [18, 18, 18]
NVF = [12, 12, 12]
MVF = [18, 18, 18]
SYMB = ['H','H','O']
SSD = None
NPOT = None
IPOT = None
JPOT = None
CPOT = None
NMASS = [1.007825, 15.9949, 1.007825]
COORD = """0.143065634079D+01  0.984836585032D+00  0.000000000000D+00
  -0.614302129101D-07 -0.124107359988D+00  0.000000000000D+00
  -0.143065634079D+01  0.984836585032D+00  0.000000000000D+00"""
NCD = None
MAXCHK = [0, 0, 0]
OVF = 'water.vib'




inputs = list()
inputs.append(TITLE)
inputs.append([NATOM,NSTAT,CONV,ICOUPL,ICOUPC,ISCFCI,IWHICH,IDISC,NROTTR,JMAX,MCHECK,INORM])
inputs.append([ICI, NMAX,CUT, EVLJ0, NVALR, KSTEP, IPRINT, MATSIZ, IREACT, MOLPRO, MOLINC])
inputs.append(XTANPM)
inputs.append(XMODQ)
inputs.append([NCONT,NVAL1,NVAL2])
inputs.append([ICONT, JCONT])
inputs.append(MAXSUM)
inputs.append(MAXBAS)
inputs.append(MODINT)
inputs.append(MEFF)
inputs.append([NCYCLE,TOLLAN,LANMAX,LGIV,LAN20,LANCUT,TOLCUT,LANCYC,TOLMAT, LANINC])
inputs.append([NRSYM,NWSYM,ISYMPG])
inputs.append([NSYM, ISYM])
inputs.append([MVSYM, MWSYM])
inputs.append(MXDIP)
inputs.append(MX)
inputs.append(MXROT)
inputs.append(ISYMT)
inputs.append([LDUMP, MDUMP, IDUMP])
inputs.append([NBF,MBF,NVF,MVF])
inputs.append(SYMB)
inputs.append(SSD)
inputs.append(NPOT)
inputs.append([IPOT, JPOT, CPOT])
inputs.append(NMASS)
inputs.append(COORD)
inputs.append(NCD)
inputs.append(MAXCHK)
inputs.append(OVF)

generate_fort1('test',inputs)