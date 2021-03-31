FC=gfortran
FFLAGS=-I.

SDIR= .

#OBJ = BPP_DDphi BPP_DDL BPP_DDE BPP_seq BPP_bent_arm

all : BPP_DDphi BPP_DDL BPP_DDE BPP_seq BPP_bent_arm

BPP_DDphi: BPP_DDphi.f90
	$(FC) -o $@ $< $(FFLAGS)

BPP_DDL: BPP_DDL.f90
	$(FC) -o $@ $< $(FFLAGS)

BPP_DDE: BPP_DDE.f90
	$(FC) -o $@ $< $(FFLAGS)

BPP_seq: BPP_seq.f90
	$(FC) -o $@ $< $(FFLAGS)

BPP_bent_arm: BPP_bent_arm.f90
	$(FC) -o $@ $< $(FFLAGS)

.PHONY : clean

clean:
	@rm -f *.o *.mod BPP_*
