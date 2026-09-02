# Compiler
FC = ifort

# Compiler flags
FFLAGS = -qopenmp -heap-arrays -mcmodel=large

# Include directories
INCLUDES = -I$(MKLROOT)/include

# Libraries
LIBS = -Wl,--start-group \
       $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
       $(MKLROOT)/lib/intel64/libmkl_sequential.a \
       $(MKLROOT)/lib/intel64/libmkl_core.a \
       -Wl,--end-group \
       -lpthread -lm

# Source files
SRC = \
    LapackRoutines.f90 \
    Setup.f90 \
    TightBinding.f90 \
    Geometry.f90 \
    HartreeFock.f90 \
    Main.f90

# Executable
TARGET = f.out

all:
	$(FC) -o $(TARGET) $(FFLAGS) $(INCLUDES) $(SRC) $(LIBS)

clean:
	rm -f *.o *.mod $(TARGET)

run: all
	mkdir -p output dataFock
	bash -c "source setMemory.sh && ./f.out"
