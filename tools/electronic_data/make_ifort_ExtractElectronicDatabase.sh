#!/bin/bash
ifort -r8 -I ../../share/hdf5-1.8.7/hdf5/include/ -c -c extract_electronic_database.f90
ifort -r8 extract_electronic_database.o -lhdf5_fortran -lhdf5 -lz -L ../../share/hdf5-1.8.7/hdf5/lib/  -o ExtractElectronicDatabase
rm extract_electronic_database.o
echo ' ------------------------------------------------'
echo ' - ExtractElectronicDatabase executable created -'
echo ' ------------------------------------------------'
