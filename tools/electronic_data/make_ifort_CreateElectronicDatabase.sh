#!/bin/bash
echo ' -----------------------------------------------'
echo ' -      Compiling CreateElectronicDatabase     -'
echo ' -----------------------------------------------'
ifort -r8 -I ../../share/hdf5-1.8.7/hdf5/include/ -c create_electronic_database.f90
ifort -r8 create_electronic_database.o -lhdf5_fortran -lhdf5 -lz -L ../../share/hdf5-1.8.7/hdf5/lib/ -o CreateElectronicDatabase
rm create_electronic_database.o
echo ' -----------------------------------------------'
echo ' - CreateElectronicDatabase executable created -'
echo ' -----------------------------------------------'
