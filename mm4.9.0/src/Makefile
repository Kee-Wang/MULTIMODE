#FC := ifort
#FFLAGS := -r8 -O 
FC := gfortran
FFLAGS := -O -fdefault-real-8 -fdefault-double-8

OBJ := \
user.H2O.o \
memo.vscf.4.9.EMORY.o \
user.vscf.4.9.0.o \
react.vscf.eckart.o \
eckart.vscf.4.9.0.o \
contr.vscf.4.9.0.o \
lanczos.vscf.4.9.0.o \
molpro.vscf.4.9.0.o \
vib.vscf.4.9.0.o \
curve.vscf.4.9.0.o \
normal.vscf.4.9.0.o \
diag.vscf.4.9.0.o \
projec.vscf.4.9.0.o \
rot.vscf.4.9.0.o \
drive.vscf.4.9.0.o \
prop.vscf.4.9.0.o \
blas.vscf.3.4.o \
lib.vscf.3.4.o \
lap.vscf.3.4.o \
matinv.vscf.2.3.o \
timit.vscf.2.3.o

%.o : %.f
	$(FC) -c $(FFLAGS) $<

%.o : %.f90
	$(FC) -c $(FFLAGS) $<

.PHONY : clean

all : h2o.x

h2o.x : $(OBJ) 
	$(FC) $^ -o $@ $(FFLAGS)

clean :
	-rm -f $(OBJ) 
