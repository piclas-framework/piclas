#!/bin/bash
# ACHTUNG: hdf5 muss vorher ohne MPI kompiliert werden! 
echo 'Compiling extractParticleData with GNU ...'
INC="../share/INTEL-SINGLE/hdf5-1.8.7/hdf5/include"
LIB='../share/INTEL-SINGLE/hdf5-1.8.7/hdf5/lib'
LIB2='../lib/libflexi.a'

# ACHTUNG: hdf5 muss vorher ohne MPI kompiliert werden! 
ifort -fpp -assume bscc -r8 -i4 -traceback -warn unused  -O2 -vec-report0  -I $INC  -DMAXWELL3D -DRK -DGTS -DPARTICLES -c extractParticleData.f90
ifort -r8 -i4 -traceback -assume bscc  -O2  -L $LIB  extractParticleData.o -lhdf5_fortran -lhdf5 -lz  -L $LIB2  -o extractParticleData

rm *.o *.mod
