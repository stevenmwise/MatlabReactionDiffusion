FOR = gfortran

FFTWLIB = /usr/local/

FLAGS = -IMODF -JMODF -O3
#FLAGS = -IMODF -JMODF -O3 -ffast-math -funroll-loops
#FLAGS = -IMODF -JMODF -fbounds-check

EXE = xnlcg

LIB = ../../LIBRARY

OBJ = OF/global.o OF/problemdef.o OF/utilities.o \
      OF/poissonfft.o OF/problemncg.o OF/searchdirection.o \
      OF/nonlinearcg.o OF/nlcgdriver.o

$(EXE): $(OBJ)
	$(FOR) $(FLAGS) $(OBJ) -o $(EXE) -L$(FFTWLIB)/lib -I$(FFTWLIB)/include -lfftw3 -lm

OF/global.o: $(LIB)/global.f90
	$(FOR) $(FLAGS) -c $(LIB)/global.f90 -o OF/global.o

OF/problemdef.o: problemdef.f90
	$(FOR) $(FLAGS) -c problemdef.f90 -o OF/problemdef.o

OF/utilities.o: $(LIB)/utilities.f90
	$(FOR) $(FLAGS) -c $(LIB)/utilities.f90 -o OF/utilities.o

OF/poissonfft.o: $(LIB)/poissonfft.f90
	$(FOR) $(FLAGS) -c $(LIB)/poissonfft.f90 -o OF/poissonfft.o -L$(FFTWLIB)/lib -I$(FFTWLIB)/include -lfftw3 -lm

OF/problemncg.o: problemncg.f90
	$(FOR) $(FLAGS) -c problemncg.f90 -o OF/problemncg.o

OF/searchdirection.o: searchdirection.f90
	$(FOR) $(FLAGS) -c searchdirection.f90 -o OF/searchdirection.o -L$(FFTWLIB)/lib -I$(FFTWLIB)/include -lfftw3 -lm

OF/nonlinearcg.o: $(LIB)/nonlinearcg.f90
	$(FOR) $(FLAGS) -c $(LIB)/nonlinearcg.f90 -o OF/nonlinearcg.o

OF/nlcgdriver.o: nlcgdriver.f90
	$(FOR) $(FLAGS) -c nlcgdriver.f90 -o OF/nlcgdriver.o

run:
	@./$(EXE)

reset:
	rm $(EXE) MODF/* OF/*

remove:
	rm OUT/*.dat output.dat 
