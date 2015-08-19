  #!/bin/bash
  ifort -r8 -I ../../share/INTEL-SINGLE/hdf5-1.8.14/include/ -c -c extract_electronic_database.f90
  ifort -r8 extract_electronic_database.o -lhdf5_fortran -lhdf5 -lz -L ../../share/INTEL-SINGLE/hdf5-1.8.14/lib/  -o ExtractElectronicDatabase # -ldl
rm extract_electronic_database.o
echo ' ------------------------------------------------'
echo ' - ExtractElectronicDatabase executable created -'
echo ' ------------------------------------------------'
