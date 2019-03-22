# Welcome to PICLas

## PICLas Documentation

The documentation is a work in progress, where everybody can contribute. A PDF file can be generated simply by executing *make* in the *doc/userguide/* folder.

1. [Introduction](doc/userguide/000_userguide.md)
2. [Installation](doc/userguide/010_installation.md)
3. [Workflow](doc/userguide/020_workflow.md)
4. [Features & Models](doc/userguide/030_features_models.md)
5. [Tools](doc/userguide/040_tools.md)
6. [Tutorials](doc/userguide/050_tutorials.md)
7. [Unit Tests](doc/userguide/070_unittest.md)
8. [Development Guide](doc/userguide/080_develop_guide.md)
9. [Parameters](doc/userguide/099_parameter.md)

## Recommended Configuration

For an installation guide please refer to the documentation. The tools are known to work with the following software versions:

* OpenMPI > 1.8
* HDF5 > 1.8.13
* CMake > 3.0.0

A list of successful build combinations can be found in [BuildConfigurations.md](BuildConfigurations.md).

### Installation on Cluster using CMake

For installation instruction on forhlr1 see [Install-on-FORHLR1.md](INSTALLFORHLR1.md).

For installation instruction on hazelhen see [Install-on-HAZELHEN.md](INSTALLHAZELHEN.md).
This instruction contains the information on using the cray-tools, too.

## Style Guide

For information about the coding style rules see [STYLEGUIDE.md](STYLEGUIDE.md).

## Particle Tracking

For information about the different particle tracking routines see [Tracking.md](TRACKING.md).

## Regressioncheck

For information about the regressioncheck see [REGGIE.md](REGGIE.md).

## Compiler Flag

For information about the compiler-flag see [COMPILER.md](COMPILER.md).

## Arbitrary Particle Distribution 

For an easy guide to arbitrary particle distributions during an initial emission see [PDFINIT.md](PDFINIT.md).


## Used libraries

PICLas uses several external libraries as well as auxilliary functions from open source projects, including:

* [HDF5](https://www.hdfgroup.org/)
* [MPI](http://www.mcs.anl.gov/research/projects/mpi/)
* [LAPACK](http://www.netlib.org/lapack/)
* [cmake](https://www.cmake.org)
* CTAGS

## I need help or further information:

* [Numerics Research Group](https://nrg.iag.uni-stuttgart.de/)
* [IRS - Numerische Modellierung und Simulation](https://www.irs.uni-stuttgart.de/forschung/numerische_modellierung_und_simulation/index.html)
* [HOPR (High Order Preprocessor)](https://hopr-project.org)


