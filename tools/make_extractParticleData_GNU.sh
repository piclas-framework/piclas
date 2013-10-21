#!/bin/bash
# ACHTUNG: hdf5 muss vorher ohne MPI kompiliert werden! 
echo 'Compiling extractParticleData with GNU ...'
INC="../share/GNU-SINGLE/hdf5-1.8.7/hdf5/include"
LIB='../share/GNU-SINGLE/hdf5-1.8.7/hdf5/lib'
#LIB2='../share/GNU-SINGLE/RECIPES'

gfortran -xf95-cpp-input -fdefault-real-8 -fdefault-double-8 -fbackslash -I $INC -c extractParticleData.f90
gfortran  extractParticleData.o -L $LIB -o extractParticleData -lhdf5_fortran -lhdf5 -lz

echo 'Compilation done.'
echo 'Cleaning ...'
rm *.o
echo '... DONE!'
