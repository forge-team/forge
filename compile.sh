ifort -qopenmp -o f.out LapackRoutines.f90 Setup.f90 TightBinding.f90 Geometry.f90 HartreeFock.f90 Main.f90 -qmkl -lpthread -lm
