CC = icpc
FC = ifort -extend_source -r8
CFLAGS = -O3 -xSSE4.2 -fp-model precise -prec_div
FFLAGS = -nofor-main -lifcore -extend_source -r8 -O3 -xSSE4.2 -fp-model precise -prec_div

HEADERS = ChemBurner.H ChemDriver.H ChemIntegrate.H ChemGlobals.H chemF.H ChemExternal.H

OBJ = ChemBurner.o ChemDriver.o ChemIntegrate.o Jacobian.o chem.o chemF.o


chem: $(OBJ)
	$(CC) $(OBJ) -lifcore -o $@

.SUFFIXES: .F90

.F90.o:
	$(FC) -c $(FFLAGS)  $<

.cpp.o: 
	$(CC) -c $(CFLAGS) $<

$(OBJ):  $(HEADERS)

clean:
	rm -f $(RESULT) *.o *.mod




