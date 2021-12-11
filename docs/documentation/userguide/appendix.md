# Appendix

## Tested compiler combinations

The following list summarizes all **tested combinations** of the required libraries (HDF5, OpenMPI, CMake etc.)

| Dev |  PICLas Version  |    System   |  Compiler |  HDF5  |      MPI      |   CMake  |
| --- |  :------------:  | :---------: |  :------: | :----: | :-----------: | :------: |
|  SC | 2.3.0 (Nov 2021) |      PC     | gcc11.2.0 | 1.12.1 | openmpi-4.1.1 |  3.21.3  |
|  SC | 2.2.0 (Nov 2021) |      PC     | gcc10.1.0 | 1.10.5 | openmpi-4.0.2 |  3.17.0  |
|  AM | 2.1.0 (Nov 2021) |      PC     | gcc9.3.0  | 1.10.6 | openmpi-3.1.6 |  3.17.0  |
|  SC | 2.0.0 (Nov 2021) |  boltzhawk  |  gcc9.3.0 | 1.10.5 | openmpi-3.1.6 |  3.17.0  |
|  SC | 2.0.0 (Nov 2021) | boltzreggie |  gcc9.2.0 | 1.10.5 | openmpi-4.0.2 |  3.15.0  |
|  SC | 2.0.0 (Nov 2021) |  boltzplatz |  gcc9.2.0 | 1.10.5 | openmpi-3.1.6 |  3.17.0  |
|  SC | 2.0.0 (Nov 2021) |     hawk    |  gcc9.2.0 | 1.10.5 |    mpt2.23    |  3.16.4  |
|  SC | 2.0.0 (Nov 2021) |     fh1     | intel18.1 |  1.10  |    impi2018   |   3.17   |
|  SC | 2.0.0 (Nov 2021) |     fh2     | intel19.1 |  1.10  |    impi2019   |   3.17   |
|  PN |  1.4.0 (Nov 19)  |  boltzplatz |  gnu7.4.0 | 1.10.5 | openmpi-3.1.3 | 3.15.3-d |
|  SC |  1.4.0 (Nov 19)  | boltzreggie |  gnu9.2.0 | 1.10.5 | openmpi-4.0.2 | 3.15.3-d |

Combinations that can cause problems are listed in the following table

| Dev | PICLas Version |    System   | Compiler |  HDF5  |      MPI      |   CMake  |                                            Notes                                            |
| --- | :------------: | :---------: | :------: | :----: | :-----------: | :------: |          :-----------------------------------------------------------------------:          |
|  SC | 1.4.0 (Nov 19) | boltzreggie | gnu9.2.0 | 1.10.5 | openmpi-4.0.1 | 3.15.3-d | Does not work for more than 3 processors probably due to a problem with the OpenMPI version |
|     |                |             |          |        |               |          |                                                                                             |
