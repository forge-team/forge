# Run executables from the FORGE root so dataFock/ is the producer's directory.
FC = gfortran
FFLAGS = -O2 -fopenmp -ffree-line-length-none -fallow-argument-mismatch
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LIBS = -framework Accelerate
else
LIBS = -llapack -lblas
endif

BUILD = .build
PROGRAMS = Compute_Spectral Compute_SpectralFull Compute_SpectralOptimize Compute_DrudeFull Compute_KramersKronig compute_localDos
UPSTREAM_PROGRAMS = compute_plotBands Main_OrderParamsFull Main_OrderParamsLowEn
COMMON = LapackRoutines.f90 ../Setup.f90 ../TightBinding.f90 ../Geometry.f90 ../HartreeFock.f90 Hamiltonian.f90 PostProcessingInput.f90
OBJECTS = LapackRoutines.o Setup.o TightBinding.o Geometry.o HartreeFock.o Hamiltonian.o PostProcessingInput.o
UPSTREAM_COMMON = ../LapackRoutines.f90 ../Setup.f90 ../TightBinding.f90 ../Geometry.f90 ../HartreeFock.f90 ../OrderParameter.f90 PostProcessingInput.f90
UPSTREAM_OBJECTS = LapackRoutines.o Setup.o TightBinding.o Geometry.o HartreeFock.o OrderParameter.o PostProcessingInput.o

.PHONY: all check
all: $(addprefix $(BUILD)/,$(PROGRAMS) $(UPSTREAM_PROGRAMS))

$(BUILD)/modules.stamp: $(COMMON) Makefile
	mkdir -p $(BUILD)
	cd $(BUILD) && $(FC) $(FFLAGS) -c $(addprefix ../,$(COMMON))
	touch $@

$(BUILD)/%: %.f90 $(BUILD)/modules.stamp
	cd $(BUILD) && $(FC) $(FFLAGS) -o $* $(addprefix ../,$<) $(OBJECTS) $(LIBS)

# Keep dev's full-matrix LAPACK wrappers separate from the response kernels'
# upper-Hermitian-triangle wrappers; both define module lapack_routines.
$(BUILD)/upstream/modules.stamp: $(UPSTREAM_COMMON) Makefile
	mkdir -p $(BUILD)/upstream
	cd $(BUILD)/upstream && $(FC) $(FFLAGS) -c $(addprefix ../../,$(UPSTREAM_COMMON))
	touch $@

$(addprefix $(BUILD)/,$(UPSTREAM_PROGRAMS)): $(BUILD)/%: %.f90 $(BUILD)/upstream/modules.stamp
	cd $(BUILD)/upstream && $(FC) $(FFLAGS) -o ../$* ../../$< $(UPSTREAM_OBJECTS) $(LIBS)

check:
	python3 tests/check_interfaces.py
	python3 tests/check_screening.py
	python3 tests/check_localdos.py
	python3 tests/check_moved_drivers.py
