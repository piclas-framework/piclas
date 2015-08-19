#!/bin/bash
echo ' -----------------------------------------------'
echo ' -      Compiling CreateElectronicDatabase     -'
echo ' -----------------------------------------------'
gfortran -fdefault-real-8 -fdefault-double-8 -I ../../share/GNU-SINGLE/hdf5-1.8.14/hdf5/include -c create_electronic_database.f90
echo ' linking....'
gfortran -fdefault-real-8 -fdefault-double-8 create_electronic_database.o  -L ../../share/GNU-SINGLE/hdf5-1.8.14/hdf5/lib -o CreateElectronicDatabase -lhdf5_fortran -lhdf5 -lz -ldl
rm create_electronic_database.o
echo ' -----------------------------------------------'
echo ' - CreateElectronicDatabase executable created -'
echo ' -----------------------------------------------'
