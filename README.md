<img src="doc/logo.png" width="582" height="287">

# Welcome to PICLas

The code is licensed under the GNU General Public License V3.
The licence can be found in [LICENCE](LICENCE.md) and the list of contributors in [CONTRIBUTORS](CONTRIBUTORS.md).
Among other things, PICLas is a scientific project.
If you use PICLas for publications or presentations in science, please support the project by citing our publications given in [references](REFERENCE.md).

## PICLas Documentation

The documentation is a work in progress. For an installation guide please refer to the documentation.

1. [Introduction](doc/userguide/000_userguide.md)
2. [Installation](doc/userguide/010_installation.md)
3. [Workflow](doc/userguide/020_workflow.md)
4. [Features & Models](doc/userguide/030_features_models.md)
5. [Visualization & Output](doc/userguide/040_visu_output.md)
6. [Tools](doc/userguide/050_tools.md)
7. [Tutorials](doc/userguide/060_tutorials.md)
8. [Development Guide](doc/userguide/080_cluster_guide.md)

A PDF file can be generated simply by executing *make* in the *doc/userguide/* folder. This requires `pandoc`, `pandoc-citeproc` and `LaTeX`.
More information about pandoc and its installation can be found [here](https://pandoc.org/installing.html).

An automatically generated documentation using doxygen can be produced in *doc/doxygen/* by executing the `builddoxy.sh` script.
This requires the `doxygen` package to be installed and is useful to get an overview over the code structure with call graphs.
The resulting documentation can be viewed by opening the `piclas/doc/doxygen/doxygen/html/index.html` file in a browser.

The tools are known to work with the following software versions:

* OpenMPI > 1.8
* HDF5 > 1.8.13
* CMake > 3.0.0

## Regression Testing

For information about the regression testing see [REGGIE](REGGIE.md).

## Used libraries

PICLas uses several external libraries as well as auxiliary functions from open source projects, including:

* [HDF5](https://www.hdfgroup.org/)
* [MPI](http://www.mcs.anl.gov/research/projects/mpi/)
* [LAPACK](http://www.netlib.org/lapack/)
* [cmake](https://www.cmake.org)

## I need help or further information:

* Academic
  * [IAG - Numerics Research Group](https://www.iag.uni-stuttgart.de/en/working-groups/numerical-methods/)
  * [IRS - Numerical Modelling and Simulation Group](https://www.irs.uni-stuttgart.de/forschung/raumtransporttechnologie/numerische_modellierung_und_simulation/)
  * [HOPR (High Order Preprocessor)](https://www.hopr-project.org/index.php/Home)
* Simulation service & training: [boltzplatz - numerical plasma dynamics GmbH](https://boltzplatz.eu)
