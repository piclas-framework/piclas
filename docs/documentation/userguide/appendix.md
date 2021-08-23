# Appendix

## Tested compiler combinations

| Dev | Version (Date) |   System    | Compiler |  HDF5  |      MPI      |  CMake   |                  Notes                   |
| --- | :------------: | :---------: | :------: | :----: | :-----------: | :------: | :--------------------------------------: |
| PN  | 1.4.0 (Nov 19) | boltzplatz  | gnu7.4.0 | 1.10.5 | openmpi-3.1.3 | 3.15.3-d |                                          |
| SC  | 1.4.0 (Nov 19) | boltzreggie | gnu9.2.0 | 1.10.5 | openmpi-4.0.1 | 3.15.3-d | Does not work for more than 3 processors |
|     | 1.4.0 (Nov 19) | boltzreggie | gnu9.2.0 | 1.10.5 | openmpi-4.0.2 | 3.15.3-d |                                          |

### Previously tested combinations

These combinations have been tested over a year ago. Thus, their applicability to the current version cannot be guaranteed.

| User |  System  |    Compiler     |  HDF5  |       MPI        |  CMake   |                       Notes                       |
| ---- | :------: | :-------------: | :----: | :--------------: | :------: | :-----------------------------------------------: |
| PO   |  Laptop  |    gnu4.9.2     | 1.8.16 |  openmpi-1.8.4   |  3.4.3   | gnu-sanitizer not working with DSMC, memory leak. |
|      | giganto  |    intel15.0    | 1.8.16 |  openmpi-1.8.2   | 2.8.12.2 |                    no autolist                    |
|      | hazelhen | intel15.0.4.223 | 1.8.14 | cray-mpich-7.3.1 |  3.4.2   |                   manual tecio                    |
|      |          | cray-libsci13.3 |  cray  |                  |          |                                                   |
|      |  Laptop  |    gnu5.2.0     | 1.8.16 |  openmpi-1.10.1  |  3.4.3   | gnu-sanitizer not working with DSMC, memory leak. |
|      |  Laptop  |    gnu7.3.+     | patch1 |  openmpi-3.0.0   |  3.10.+  |       Requires HDF_ROOT instead of HDF5_DIR       |
| SC   |  Laptop  |    gnu4.8.4     | 1.8.16 |  openmpi-1.6.5   |  3.2.2   |                                                   |
|      |  Laptop  |    gnu5.4.0     | 1.8.18 |  openmpi-1.8.8   |  3.5.1   |                                                   |
|      | hazelhen | intel15.0.4.223 | 1.8.14 | cray-mpich-7.3.1 |  3.4.2   |   set tecio path by hand (copy from old PICLas)   |
| WR   |  Laptop  |    gnu5.2.0     | 1.8.16 |  openmpi-1.10.0  |  3.4.3   |  linking only works with gnu5.2.0    --> solved   |
|      |  Laptop  |    gnu4.8.4     | 1.8.16 |  openmpi-1.10.0  |  3.4.3   |                                                   |
|      |  Laptop  |    gnu4.8.4     | 1.8.14 |  openmpi-1.6.5   |  3.4.3   |                                                   |
|      |  Laptop  |    gnu4.8.4     | 1.8.16 |  openmpi-1.6.5   |  3.4.3   |                                                   |
|      |  Laptop  |    gnu5.2.0     | 1.8.16 |  openmpi-1.6.5   |  3.4.3   |                                                   |
|      |  Laptop  |   intel15.0.4   | 1.8.16 |  openmpi-1.10.0  |  3.4.3   |                                                   |

