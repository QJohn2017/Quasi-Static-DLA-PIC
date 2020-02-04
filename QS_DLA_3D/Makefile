PROGRAM.O = Fields_Particles.o Pulse.o Test_Particle.o Magnectic_Fields.o Wake_Fields.o Electron_Pusher.o OneZstep.o OneTstep.o Initialization.o New_Wake_VW.o

FC = gfortran
LD = gfortran

FCFLAGS = -O3  
LDFLAGS = -O3  
INC = -I$ /usr/local/include
LIB = -L$ /usr/local/lib -lnetcdff

EXE = QSDLA

all : program

program: $(EXE)

$(EXE): $(PROGRAM.O)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIB)
%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $@ $< $(INC)
clean:
	rm -r $(PROGRAM.O) $(EXE) *.mod