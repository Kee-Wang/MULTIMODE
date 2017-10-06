
try:
    # for Python2
    from Tkinter import *
except ImportError:
    # for Python3
    from tkinter import *

class MMHelper():
    def __init__(self, root):
        self.TITLE = 'water example'
        self.NATOM = 3
        self.NSTAT = -1
        self.CONV = '1.D-2'
        self.ICOUPL = 3
        self.ICOUPC = 3
        self.ISCFCI = 200
        self.IWHICH = 1
        self.IDISC = 0
        self.NROTTR = 0
        self.JMAX = 0
        self.MCHECK = 1
        self.INORM = 1
        self.ICI=-1
        self.NMAX=-3
        self.CUT = '1.0D+4'
        self.EVLJ0 = 0.0
        self.NVALR = 20
        self.KSTEP =0
        self.IPRINT = -3
        self.MATSIZ = 0
        self.IREACT = 0
        self.MOLPRO = 0
        self.MOLINC = 0
        self.XTANPM = None
        self.XMODQ = [1.0, 1.0, 1.0]
        self.NCONT = 1
        self.NVAL1 = 0
        self.NVAL2 = 0
        self.ICONT = -3
        self.JCONT = '1 2 3   3 3'
        self.MAXSUM = [12, 12, 12]
        self.MAXBAS = [[11, 12, 12],[12,12,12],[13,12,12]]
        self.MODINT = [1, 1, 1]
        self.MEFF = [[0,0,0],[0,0,0]]
        self.NCYCLE=0
        self.TOLLAN = '1.D-3'
        self.LANMAX = 2000
        self.LGIV = 1
        self.LAN20 = 500
        self.LANCUT = 0
        self.TOLCUT = '0.D0'
        self.LANCYC = 200
        self.TOLMAT = '1.D-10'
        self.LANINC = 4
        self.NRSYM = 1
        self.NWSYM = 1
        self.ISYMPG = 0
        self.NSYM = 3
        self.ISYM = [1,2,3]
        self.MVSYM = 1
        self.MWSYM = 1
        self.MXDIP = [1,1,1]
        self.MX = [1,2,3]
        self.MXROT = [1,1,1]
        self.ISYMT = 1
        self.LDUMP = 0
        self.MDUMP = 0
        self.IDUMP = [0, 0]
        self.NBF = [12, 12, 12]
        self.MBF = [18, 18, 18]
        self.NVF = [12, 12, 12]
        self.MVF = [18, 18, 18]
        self.SYMB = ['H','H','O']
        self.SSD = None
        self.NPOT = None
        self.IPOT = None
        self.JPOT = None
        self.CPOT = None
        self.NMASS = [1.007825, 15.9949, 1.007825]
        self.COORD = """0.143065634079D+01  0.984836585032D+00  0.000000000000D+00
          -0.614302129101D-07 -0.124107359988D+00  0.000000000000D+00
          -0.143065634079D+01  0.984836585032D+00  0.000000000000D+00"""
        self.NCD = None
        self.MAXCHK = [0, 0, 0]
        self.OVF = 'water.vib'
        self.LINEAR_check_button = IntVar()
        self.LINEAR = 0
        self.remember = 0
        self.build_window()


    def build_window(self):
        def longframe(text,default_value = '',info=['Quick Info','Sorry, no info, please see manual.']):
            def infobox(event, info):
                import tkinter
                from tkinter import messagebox
                tkinter.messagebox.showinfo(*info)

            fm = Frame(root)
            fm.pack()
            label = Label(fm, text=text)
            label.pack(side=LEFT)
            entry = Entry(fm)
            entry.pack(side=LEFT)
            entry.insert(0, default_value)

            label.bind('<Button-1>', lambda event: infobox(event,info))
            return entry



        bottomFrame = Frame(root)
        bottomFrame.pack(side=BOTTOM)

        Label(root,text="MULTIMODE's Little Helper ",font=("Helvetica", 25)).pack()
        Label(root,text="for MM5.1.4  --ver. 0.01 \n").pack()
        Label(root,text='(Click keywords before entry for quick explanation)').pack()
        #Label(root,text="This program generates `fort.1` file according to your need. \n All entries are required.\nClick the parameter words for explanations and example.").pack()

        Label(root, text="1. Molecule System Information", bg="red", fg="white",font=("Helvetica", 16)).pack()


        self.TITLE_entry = longframe('Name of project: ',"My MULTMODE Project",['Info','The name used inside `fort.1` file, not the filename.'])
        self.NATOM_entry = longframe('Number of atoms: ',3)
        self.SYMB_entry = longframe('Nuclear Symbols: ','H H',['Example', 'Tips: \n1. Separate with spaces.\n2. Order of symbols has to be consistent with your potential.\nFor example, water could be: \nO H H'])
        self.NMASS_entry = longframe('Nuclear Mass (in a.u.): ','11')


        Label(root, text='Equilibrium Cartesian Coordinates (in Bohr): ').pack()


        def retrieve_input():
            inputValue=COORD.get("1.0","end-1c")
            print(inputValue)

        text_frame = Frame(root, highlightbackground="green", highlightcolor="green", highlightthickness=1, width=100, height=100, bd= 0)
        text_frame.pack()
        self.COORD_text = Text(text_frame, height=5, width=20)
        self.COORD_text.pack(side=TOP)
        self.COORD_text.insert(END, "Copy-Paste to here.\nFormat:\nx1 y1 z1\nx2 y2 z2 \nx3 y3 z3")


        Checkbutton(root, text="Ths is linear molecule", variable = self.LINEAR_check_button).pack()
        #LINEAR.pack()


        Label(root, text="2. VSCF/VCI configurations", bg="red", fg="white",font=("Helvetica", 16)).pack()

        self.ICOUPL_entry = longframe('Number of modes coupled\n in potential integration (ICOUPL): ', 4, ['Manual: ICOUPL', 'The number of modes truncated in VSCF calculation. The higher ICOUPL the more computationally expensive. The typical convergence value is ICOUPL = 4'])

        self.NMAX_entry = longframe('Number of mode coupling\n for VCI basis(NMAX): ',4)
        self.MAXSUM_entry = longframe('MAXSUM',12)
        self.MAXBAS_entry = longframe('MAXBAS for each mode',12)
        self.NBF_entry = longframe('NBF',12)
        self.NROTTR_entry = longframe('NROTTR', '12',['Quick Info',"Number of rotational and translational degrees of freedom to beincluded in 'vibrational' modes \nnumber of modes = NMODE=3*NATOM-6+NROTTR"])




        def prt(even):
            import tkinter
            from tkinter import messagebox
            tkinter.messagebox.showinfo('Executable Missing','Did not find executable, please load `mm.x` executable file.')




        Label(root, text="3. Utilities", bg="red", fg="white",font=("Helvetica", 16)).pack()
        remember = Checkbutton(root, text="Remember Settings")
        remember.pack()


        Button(root,text='Load existing `fort.1`').pack()
        save = Button(root,text='Generate `fort.1`')
        save.pack()
        save.bind('<Button-1>', lambda event: self.update(event))


        #param = [NATOM, NSYMB, NMASS, COORD, LINEAR, ICOUPL, NMAX, MAXSUM, MAXBAS, NBF,NROTTR]  #save.bind('<Button-1>', lambda event: self.update(event,param))

        Button(root,text='Get VCI Matrix Size only').pack()
        b = Button(bottomFrame,text='run `./mm.x fort.1`')
        b.pack(side = LEFT)
        b.bind('<Button-1>', lambda event: prt(event))
        Button(bottomFrame,text='Load `mm.x` executable file').pack(side = LEFT)

    def update(self, even):

        import tkinter
        from tkinter import messagebox
        #[NATOM, NSYMB, NMASS, COORD, LINEAR, ICOUPL, NMAX, MAXSUM, MAXBAS, NBF, NROTTR] = param
        # for p in param:
        #     try:
        #         if len(str(p.get())) is 0: #For Entry element, use str to bypass int
        #             tkinter.messagebox.showinfo('Warning', 'Detected empty entry, operation aborted.')
        #             return
        #         else:
        #             continue
        #     except:
        #         if  p.get("1.0", "end-1c") is 0: #For Text element
        #             tkinter.messagebox.showinfo('Warning','Detected empty entry, operation aborted.')
        #             return
        #         else:
        #             continue

        self.TITLE = self.TITLE_entry.get()
        self.NATOM = int(self.NATOM_entry.get())
        self.SYMB = self.SYMB_entry.get()

        self.NMASS = self.NMASS_entry.get()

        self.COORD = self.COORD_text.get("1.0", "end-1c")
        self.ICOUPL = int(self.ICOUPL_entry.get())
        self.NMAX = int(self.NMAX_entry.get())
        self.MAXSUM = int(self.MAXSUM_entry.get())
        self.MAXBAS = int(self.MAXBAS_entry.get())
        self.LINEAR = int(self.LINEAR_check_button.get())

        #self.NBF_temp = int(self.NBF_entry.get())
        Nnode = 3 * self.NATOM - 6 + self.LINEAR
        print('button',self.LINEAR_check_button.get())
        NBF = list()
        MBF = list()
        NVF = list()
        MVF = list()

        for node in range(Nnode):
            NBF.append(int(self.NBF_entry.get()))
            MBF.append(self.MBF[0])
            NVF.append(self.NVF[0])
            MVF.append(self.MVF[0])
        self.NBF = NBF
        self.MBF = MBF
        self.NVF = NVF
        self.MVF = MVF

        MAXBAS = list()
        MAXSUM = list()

        print('Nnode',Nnode)
        for mode in range(self.NMAX):
            MAXSUM.append(int(self.MAXSUM))
            MAXBAS_temp = list()
            for node in range(Nnode):
                MAXBAS_temp.append(int(self.MAXBAS))
            MAXBAS.append(MAXBAS_temp)

        self.MAXSUM = MAXSUM
        self.MAXBAS = MAXBAS
        #print(self.MAXSUM, self.MAXBAS)
        self.NMAX = -4



        #print('No. Nnode =',Nnode)



        self.NROTTR = int(self.NROTTR_entry.get())

        inputs = list()
        inputs.append(self.TITLE)
        inputs.append([self.NATOM, self.NSTAT, self.CONV, self.ICOUPL, self.ICOUPC, self.ISCFCI, self.IWHICH, self.IDISC, self.NROTTR, self.JMAX, self.MCHECK, self.INORM])
        inputs.append([self.ICI, self.NMAX, self.CUT, self.EVLJ0,self.NVALR, self.KSTEP, self.IPRINT, self.MATSIZ, self.IREACT, self.MOLPRO, self.MOLINC])
        inputs.append(self.XTANPM)
        inputs.append(self.XMODQ)
        inputs.append([self.NCONT, self.NVAL1, self.NVAL2])
        inputs.append([self.ICONT, self.JCONT])
        inputs.append(self.MAXSUM)
        inputs.append(self.MAXBAS)
        inputs.append(self.MODINT)
        inputs.append(self.MEFF)
        inputs.append([self.NCYCLE, self.TOLLAN, self.LANMAX, self.LGIV, self.LAN20, self.LANCUT, self.TOLCUT, self.LANCYC, self.TOLMAT, self.LANINC])
        inputs.append([self.NRSYM, self.NWSYM, self.ISYMPG])
        inputs.append([self.NSYM, self.ISYM])
        inputs.append([self.MVSYM, self.MWSYM])
        inputs.append(self.MXDIP)
        inputs.append(self.MX)
        inputs.append(self.MXROT)
        inputs.append(self.ISYMT)
        inputs.append([self.LDUMP, self.MDUMP, self.IDUMP])
        inputs.append([self.NBF, self.MBF, self.NVF, self.MVF])
        inputs.append(self.SYMB)
        inputs.append(self.SSD)
        inputs.append(self.NPOT)
        inputs.append([self.IPOT, self.JPOT, self.CPOT])
        inputs.append(self.NMASS)
        inputs.append(self.COORD)
        inputs.append(self.NCD)
        inputs.append(self.MAXCHK)
        inputs.append(self.OVF)



        self.generate_fort1(inputs)
        tkinter.messagebox.showinfo('Info','Congratulations, `fort.1` is generated successfully.')



    def generate_fort1(self,p):
        """

        :param p: An nested array of parameters. p[i] is the (i-1)th line of parameters in orginal `fort.1` input.
        :return: None

        TO-DO:
        1. Input error control
        """
        "Conversion Table"
        TITLE = p[0]
        [NATOM, NSTAT, CONV, ICOUPL, ICOUPC, ISCFCI, IWHICH, IDISC, NROTTR, JMAX, MCHECK, INORM] = p[1]
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
        NMASS = p[25]
        COORD = p[26]
        MAXCHK = p[28]
        OVF = p[29]

        "Generating fort.1 file"

        f = open('fort.1', 'w')
        f.write('C**TITLE \n')
        f.write('{:s}\n'.format(p[0]))
        f.write('C**NATOM,NSTAT,CONV,ICOUPL,ICOUPC,ISCFCI,IWHICH,IDISC,NROTTR,JMAX,MCHECK,INORM \n')
        print(p[1])
        f.write('   {:^6d}{:^5d}{:^5s}{:^7d}{:^7d}{:^7d}{:^7d}{:^6d}{:^8d}{:^5d}{:^7d}{:^5d}\n'.format(*p[1]))
        f.write('C**ICI NMAX  CUT     EVLJ0   NVALR KSTEP IPRINT MATSIZ IREACT MOLPRO MOLINC \n')
        f.write('  {:^5d}{:<6d}{:<8s}{:^5.2f}{:^8d}{:^7d}{:^7d}{:^7d}{:^7d}{:^8d}{:^6d}\n'.format(*p[2]))
        f.write('C**XTANPM (DATA IGNORED IF MOLPRO.EQ.0 AND IWHICH.GE.0)\n')
        if MOLPRO == 0 and IWHICH >= 0:
            pass
        else:
            f.write(p[3])  # Need to work on the format of this!
        f.write('C**XMODQ \n')
        f.write('   {:^5.1f}{:^5.1f}{:^5.1f}\n'.format(*p[4]))
        f.write('C**NCONT NVAL1 NVAL2 (DATA IGNORED IF NMAX.GE.0 OR ISCFCI.LE.0) \n')
        f.write('   {:^5d} {:^5d} {:^5d} \n'.format(*p[5]))
        f.write('C**ICONT JCONT (DATA IGNORED IF NMAX.GE.0 OR ISCFCI.LE.0) \n')
        f.write('   {:<5d} {:<30s} \n'.format(*p[6]))  # Here JCONT used string
        f.write('C**MAXSUM (DATA IGNORED IF NMAX.GE.0 OR ISCFCI.LE.0) \n')
        if NMAX >= 0 or ISCFCI <= 0:
            pass
        else:
            f.write(('{:5d}' * len(MAXSUM) + '\n').format(*MAXSUM))
        f.write('C**MAXBAS (DATA IGNORED IF NMAX.GE.0 OR ISCFCI.LE.0) \n')
        if NMAX >= 0 or ISCFCI <= 0:
            pass
        else:
            for nbasis in MAXBAS:  # Iterate MAXBAS
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
        f.write(('  {:^6d}' + '{:3d}' * len(ISYM) + '\n').format(NSYM, *ISYM))
        f.write('C**MVSYM   MWSYM \n')
        f.write('   {:^5d}   {:^5d}\n'.format(*p[14]))
        f.write('C**MXDIP (I.R. of mu(X), mu(Y), mu(Z))\n')
        f.write(('{:5d}' * len(MXDIP) + '\n').format(*MXDIP))
        f.write('C**MX (symmetry axes x,y,z corresponding to principle axes X,Y,Z)\n')
        f.write(('{:5d}' * len(MX) + '\n').format(*MX))
        f.write('C**MXROT (I.R. of Rx, Ry, Rz) \n')
        f.write(('{:5d}' * len(MXROT) + '\n').format(*MXROT))
        f.write('C**ISYMT (I.R. of Sin(M.tau/2)\n')
        f.write('{:5d}\n'.format(p[18]))
        f.write('C**LDUMP MDUMP IDUMP (DATA IGNORED IF ICI.GE.0)\n')
        if ICI >= 0:
            pass
        else:
            f.write((('  ' + '{:^5d}' * 2 + '  ' + '{:>3d}' * len(IDUMP) + '\n').format(LDUMP, MDUMP, *IDUMP)))
        f.write('C**NBF,MBF,NVF,MVF \n')
        for i in range(len(NBF)):
            f.write(('{:5d}' * 4 + '\n').format(NBF[i], MBF[i], NVF[i], MVF[i]))
        f.write('C**NUCLEAR SYMBOLS\n')
        f.write((' ' * 4 + '{:4s}' * 3 + '\n').format(*p[21]))
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
        #f.write(('{:13.8f}' * 3 + '\n').format(*NMASS))
        f.write(('{:s}'  + '\n').format(*NMASS))
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



if __name__ == '__main__':
    root = Tk()

    app = MMHelper(root)

    root.mainloop()




    #
    # def get(event):
    #     return entry_natom.get()
    #
    # var_NATOM = StringVar(root)
    #
    # label_info = Label(fm, text='Molecule System Information').pack(side=LEFT)
    #
    # label_natom = Label(fm, text='Number of atoms: ').pack(side=LEFT)
    # entry_natom = Entry(fm, textvariable=var_NATOM)
    # entry_natom.pack(side=LEFT)
    # entry_natom.bind('<Return>', get)
    # NATOM = get
    # print(NATOM,'tesss')
    #
    # line2=Frame(root)
    # line2.pack()


    # natom_entry = Entry(line2,textvariable=sv)
    # natom_entry.pack(side=LEFT)
    # natom_entry.delete(0, END)
    # natom_entry.insert(0, "Enter value here")
    # NATOM = natom_entry.get()
    # bottomFrame = Frame(root)
    # bottomFrame.pack(side=BOTTOM)
    # def printName(event, NATOM):
    #     import tkinter
    #     from tkinter import messagebox
    #     print('test')
    #     tkinter.messagebox.showinfo('Window Title', NATOM)
    # button1 = Button(bottomFrame, text = 'Get VCI Matrix Size')
    # button2 = Button(bottomFrame, text = 'Print NATOM')
    # button2.bind('<Button-1>',lambda event:printName(event, var_NATOM.get()))
    # button2.pack()
    # button1.pack()







    #topFrame = Frame(root) #
    #topFrame.pack() #If you want to show this frame, you have to pack it in the window. Default is top.
    #label_1 = Label(root,text='Name')
    #theLabel = Label(root,text="MULTIMODE's Little Helper").pack(side=TOP)
    #Entry(text="Entry").pack(side=TOP)


    #bottomFrame = Frame(root)
    #bottomFrame.pack(side=BOTTOM) # Pack this frame on the bottom.
    #button1 = Button(bottomFrame,text='Get VCI Matrix Size')
    #button1.pack(side=BOTTOM)
    #theLabel.pack() #label pops out really fast and close really fast.






     #This means have window continuous on the screen until
