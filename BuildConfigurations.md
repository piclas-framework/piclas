# Build Configurations

## Validated Compiler Combinations

*  Master

| User                | System        | Compiler        | HDF5  | MPI             | CMake  | Makefile | Notes       |
| ------------------- |:-------------:| :--------------:|:-----:|:---------------:|:------:|:-------: |:-----------:|
| Philip Ortwein      | Laptop        | gnu4.9.2        |1.8.16 |openmpi-1.8.4    |3.4.3   |          |             |
|                     | giganto       | intel15.0       |1.8.16 |openmpi-1.8.2    |2.8.12.2|          | no autolist |
|                     | hazelhen      | intel15.0.4.223 |1.8.14 |cray-mpich-7.3.1 |3.4.2   |          | manual tecio|
|                     |               | cray-libsci13.3 |cray   |                 |        |          |             |
| Stephen Copplestone | Laptop        | gnu4.8          |1.8.16 |openmpi-1.6.5    | yes    |          |             |