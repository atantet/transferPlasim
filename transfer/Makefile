OPT=/opt/local/
INCDIR=-I$(PHD)/dev/include -I$(HOME)/local/include -I$(OPT)/include/ -I$(HOME)/local/arpack++/include
LIBDIR=-L$(HOME)/local/lib -L$(OPT)/lib/
LIBSSTD=-lgsl -lgslcblas -lm -larpack -llapack -lsuperlu_4.3 -lcblas -lgfortran
LIBSOMP=$(LIBSSTD) -lgomp

CC=g++
CFLAGSSTD=-Wall -ansi -pedantic-errors -O3
CFLAGSOMP=$(CFLAGSSTD) -fopenmp
LDFLAGSSTD=-undefined dynamic_lookup
LDFLAGSGDB=$(LDFLAGSSTD) -gdb

LDFLAGS=$(LDFLAGSSTD)
#LDFLAGS=$(LDFLAGSGDB)

EXE=spectrumPlasim
CFLAGS=$(CFLAGSSTD)
LIBS=$(LIBSSTD)

#EXE=transferPlasim
#CFLAGS=$(CFLAGSOMP)
#LIBS=$(LIBSOMP)

all:$(EXE)

$(EXE).o:$(EXE).cpp
	$(CC) -c $(CFLAGS) $(INCDIR) $(EXE).cpp

$(EXE):$(EXE).o
	$(CC) $(LDFLAGS) -o $(EXE) $(LIBDIR) $(EXE).o $(LIBS)

clean:
	rm -f *.o $(EXE)
