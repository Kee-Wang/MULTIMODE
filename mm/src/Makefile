FC := gfortran
FFLAGS := -fdefault-real-8 -fdefault-double-8 -O

OBJ := \
 user.H2O.o\
 blas.vscf.3.4.o\
 contr.vscf.5.1.1.o\
 curve.vscf.5.1.4.MVF.o\
 diag.vscf.5.1.3.o\
 drive.vscf.5.1.4.MVF.o\
 eckart.vscf.5.1.1.o\
 lanczos.vscf.5.1.4.MVF.o\
 lap.vscf.3.4.o\
 lib.vscf.3.4.o\
 matinv.vscf.2.3.o\
 memo.vscf.5.1.3.o\
 molpro.vscf.5.0.o\
 normal.vscf.4.9.0.o\
 projec.vscf.4.9.1.o\
 prop.vscf.5.1.3.o\
 psplot.o\
 rot.vscf.5.1.4.MVF.o\
 spectrum.o\
 timit.vscf.4.9.1.o\
 user.react.eckart.5.1.3.o\
 user.vscf.4.9.2.o\
 vib.vscf.5.1.4.MVF.o

%.o : %.f
	$(FC) -c $< $(FFLAGS)

%.o : %.f90
	$(FC) -c $< $(FFLAGS)

mm.x : $(OBJ) 
	$(FC) $^ -o $@ $(FFLAGS)

clean :
	-rm -f $(OBJ) mm.x 
