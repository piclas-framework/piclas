# Workflow

In this chapter, the complete process of setting up a simulation with **PICLas** is detailed.

(sec:compiler-options)=
## Compiler options
This section describes the main configuration options which can be set when building **PICLas** using CMake.
Some options are dependent on others being enabled (or disabled), so the available ones may change.

The first set of options describe general CMake behaviour:

* ``CMAKE_BUILD_TYPE``: This statically specifies what build type (configuration) will be built in this build tree. Possible values are
    * Release: "Normal" execution.
    * Profile: Performance profiling using gprof.
    * Debug: Debug compiler for detailed error messages during code development.
    * SANI: Sanitizer compiler for even more detailed error messages during code development.

* ``CMAKE_HOSTNAME``: This will display the host name of the machine you are compiling on.

* ``CMAKE_INSTALL_PREFIX``: If “make install” is invoked or INSTALL is built, this directory is prepended onto all install directories
This variable defaults to /usr/local on UNIX.

For some external libraries and programs that **PICLas** uses, the following options apply:

* ``CTAGS_PATH``: This variable specifies the Ctags install directory, an optional program used to jump between tags in the source file.

* ``LIBS_BUILD_HDF5``: This will be set to ON if no pre-built HDF5 installation was found on your machine. In this case a HDF5 version
will be build and used instead. For a detailed description of the installation of HDF5, please refer to Section {ref}`sec:hdf5-installation`.

* ``HDF5_DIR``: If you want to use a pre-built HDF5 library that has been build using the CMake system, this directory should contain
the CMake configuration file for HDF5 (optional).

* ``PICLAS_BUILD_POSTI``: Enables the compilation of additional tools and activates the following options:
  * ``POSTI_BUILD_SUPERB``: Enables the compilation of **superB**, which is allows the computation of magnetic fields based on an
  input of coils and permanent magnets, see Section {ref}`sec:superB`
  * ``POSTI_BUILD_VISU``: Enables the compilation of the post-processing tool **piclas2vtk**, which enables the conversion of
  output files into the VTK format
  * ``POSTI_USE_PARAVIEW``: Enables the compilation of the ParaView plugin, which enables the direct read-in of output files within ParaView

(sec:solver-settings)=
## Solver settings

Before setting up a simulation, the code must be compiled with the desired parameters. The most important compiler options to be set are:

* ``PICLAS_TIMEDISCMETHOD``: Module selection
    * DSMC: Direct Simulation Monte Carlo
    * RK4: Time integration method Runge-Kutta 4th order in time
* ``PICLAS_EQNSYSNAME``: Equation system to be solved
    * maxwell:
    * poisson:
* ``PICLAS_POLYNOMIAL_DEGREE``: Defines the polynomial degree of the solution. The order of convergence follows as $N+1$. Each grid
cell contains $(N+1)^3$ collocation points to represent the solution.
* ``PICLAS_NODETYPE``: The nodal collocation points used during the simulation
    * GAUSS:
    * GAUSS-LOBATTO:
* ``PICLAS_INTKIND8``: Enables simulations with particle numbers above 2 147 483 647
* ``PICLAS_READIN_CONSTANTS``: Enables user-defined natural constants for the speed of light *c0*, permittivity *eps* and
    permeability *mu* of vacuum,
    which must then be supplied in the parameter file. The default if *OFF* and the values for the speed of light c0=299792458.0 [m/s], permittivity
    eps=8.8541878176e-12 [F/m] and permeability mu=1.2566370614e-6 [H/m] are hard-coded.

The options EQNSYSNAME, POLYNOMIAL_DEGREE and NODETYPE can be ignored for a DSMC simulation. For parallel computation the following
flags should be configured:

* ``LIBS_USE_MPI``: Enabling parallel computation. For a detailed description of the installation of MPI, please refer to refer to
                    Section {ref}`sec:installing-mpi`.
* ``PICLAS_LOADBALANCE``: Enable timer-based load-balancing by automatic determination of workload weights for each simulation
                          element.

All other options are set in the parameter file.

## Setup of parameter file(s)

The settings of the simulation are controlled through parameter files, which are given as arguments to the binary. In the case of
PIC simulations the input of a single
parameter file (e.g. *parameter.ini*) is sufficient, while the DSMC method requires the input of a species parameter file (e.g.
*DSMC.ini*). The most recent list of parameters can be found by invoking the help in the console:

    piclas --help

General parameters such the name of project (used for filenames), the mesh file (as produced by HOPR), end time of the simulation (in seconds) and the time step, at which the particle data is written out (in seconds), are:

    ProjectName    = TestCase
    MeshFile       = test_mesh.h5
    TEnd           = 1e-3
    Analyze_dt     = 1e-4
    ManualTimeStep = 1e-4 (over-rides the automatic time step calculation in the Maxwell solver)

Generally following types are used:

~~~~~~~
INTEGER = 1
REAL    = 1.23456
REAL    = 1.23E12
LOGICAL = T         ! True
LOGICAL = F         ! False
STRING  = PICLAS
VECTOR  = (/1.0,2.0,3.0/)
~~~~~~~

The concept of the parameter file is described as followed:

* Each single line is saved and examined for specific variable names
* The examination is case-insensitive
* Comments can be set with symbol "!" in front of the text
* Numbers can also be set by using "pi"
~~~~~~~
    vector = (/1,2Pi,3Pi/)
~~~~~~~
* The order of defined variables is with one exception irrelevant, except for the special case when redefining boundaries.
However, it is preferable to group similar variables together.

The options and underlying models are discussed in Chapter {ref}`features-and-models/index:Features & Models`, while the available 
output options are given in Chapter {ref}`visu_output:Visualization & Output`.
Due to the sheer number of parameters available, it is advisable to build upon an existing parameter file from one of the tutorials
in Chapter {ref}`tutorials/index:Tutorials`.

## Simulation

After the mesh generation, compilation of the binary and setup of the parameter files, the code can be executed by

    piclas parameter.ini [DSMC.ini]

The simulation may be restarted from an existing state file

    piclas parameter.ini [DSMC.ini] [restart_file.h5]

A state file is generated at the end of the simulation and also at every time step defined by `Analyze_dt`. **Note:** When
restarting from an earlier time (or zero), all later state files possibly contained in your directory are deleted!

After a successful simulation, state files will be written out in the HDF5 format preceded by the project name, file type (e.g.
State, DSMCState, DSMCSurfState) and the time stamp:

    TestCase_State_001.5000000000000000.h5
    TestCase_DSMCState_001.5000000000000000.h5

The format and floating point length of the time stamp *001.5000000000000000* can be adjusted with the parameter

    TimeStampLength = 21

where the floating format with length of *F21.14* is used as default value.

### Parallel execution
The simulation code is specifically designed for (massively) parallel execution using the MPI library. For parallel runs, the code
must be compiled with `PICLAS_MPI=ON`. Parallel execution is then controlled using `mpirun`

    mpirun -np [no. processors] piclas parameter.ini [DSMC.ini] [restart_file.h5]

The grid elements are organized along a space-filling curved, which gives a unique one-dimensional element list. In a parallel run,
the mesh is simply divided into parts along the space filling curve. Thus, domain decomposition is done *fully automatic* and is
not limited by e.g. an integer factor between the number of cores and elements. The only limitation is that the number of cores
may not exceed the number of elements.

## Post-processing

**PICLas** comes with a tool for visualization. The piclas2vtk tool converts the HDF5 files generated by **PICLas** to the binary
VTK format, readable by many visualization tools like ParaView and VisIt. The tool is executed by

~~~~~~~
piclas2vtk [posti.ini] output.h5
~~~~~~~

Multiple HDF5 files can be passed to the piclas2vtk tool at once. The (optional) runtime parameters to be set in `posti.ini` are
given in Chapter {ref}`visu_output:Visualization & Output`.
