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
    * Sanitize: Sanitizer compiler for even more detailed error messages during code development.
    * Nitro: Fast compiler option `-Ofast` for even more speed but at the cost of accuracy.

* ``CMAKE_HOSTNAME``: This will display the host name of the machine you are compiling on.

* ``CMAKE_INSTALL_PREFIX``: If “make install” is invoked or INSTALL is built, this directory is prepended onto all install directories
This variable defaults to /usr/local on UNIX.

For some external libraries and programs that **PICLas** uses, the following options apply:

* ``CTAGS_PATH``: This variable specifies the Ctags install directory, an optional program used to jump between tags in the source file.

* `LIBS_BUILD_HOPR`: Enable the compilation of the mesh pre-processor HOPR during the PICLas compilation. The executable `hopr` will be placed in the build/bin/ folder next to the other executables. For more details, on the utilization of HOPR, see {ref}`sec:mesh-generation`.

* `LIBS_DOWNLOAD_HOPR`: Enable downloading the mesh pre-processor HOPR during the PICLas compilation from GitHub. The executable `hopr` will be linked in the build/bin/ folder next to the other executables. For more details, on the utilization of HOPR, see {ref}`sec:mesh-generation`.

* ``LIBS_BUILD_HDF5``: This will be set to ON if no pre-built HDF5 installation was found on your machine. In this case a HDF5 version
will be built and used instead. For a detailed description of the installation of HDF5, please refer to Section {ref}`sec:hdf5-installation`.

* ``HDF5_DIR``: If you want to use a pre-built HDF5 library that has been built using the CMake system, this directory should contain
the CMake configuration file for HDF5 (optional).

* ``PICLAS_MEASURE_MPI_WAIT``: Measure the time that is spent in MPI_WAIT() and output info to MPIW8Time.csv and MPIW8TimeProc.csv

* ``PICLAS_BUILD_POSTI``: Enables the compilation of additional tools and activates the following options:
  * ``POSTI_BUILD_SUPERB``: Enables the compilation of **superB**, which is allows the computation of magnetic fields based on an
  input of coils and permanent magnets, see Section {ref}`sec:superB`
  * ``POSTI_BUILD_VISU``: Enables the compilation of the post-processing tool **piclas2vtk**, which enables the conversion of
  output files into the VTK format
  * ``POSTI_USE_PARAVIEW``: Enables the compilation of the ParaView plugin, which enables the direct read-in of output files within ParaView

* ``PICLAS_SHARED_MEMORY``: Split type for creating new communicators based on colors and keys (requires MPI 3 or higher).
  Options with the prefix OMPI_ are specific to Open MPI.
  * ``MPI_COMM_TYPE_SHARED``: creates one shared memory domain per physical node (default)
  * ``OMPI_COMM_TYPE_CORE``:  creates one shared memory domain per MPI thread
  * ``PICLAS_COMM_TYPE_NODE``: creates one shared memory domain per X numbers of MPI threads defined by ``PICLAS_SHARED_MEMORY_CORES``
    * ``PICLAS_SHARED_MEMORY_CORES``: Number of MPI threads per virtual node (default is 2). Assumes that all MPI threads run on the
      same physical node.

Some settings are not shown in the graphical user interface, but can be changed via command line

* ``PICLAS_INSTRUCTION``: Processor instruction settings (mainly depending on the hardware on which the compilation process is
  performed or the target hardware where piclas will be executed). This variable is set automatically depending on the machine where
  piclas is compiled. CMake prints the value of this parameter during configuration

      -- Compiling Nitro/Release/Profile with [GNU] (v12.2.0) fortran compiler using PICLAS_INSTRUCTION [-march=native] instructions.

  When compiling piclas on one machine and executing the code on a different one, the instruction setting should be set to
  `generic`. This can be accomplished by running

      cmake -DPICLAS_INSTRUCTION=-mtune=generic

  To reset the instruction settings, run cmake again but with

      -DPICLAS_INSTRUCTION=

  which resorts to using the automatic determination depending on the detected machine.

(sec:solver-settings)=
## Solver settings

Before setting up a simulation, the code must be compiled with the desired parameters. The most important compiler options to be set are:

* ``PICLAS_TIMEDISCMETHOD``: Time integration method
    * Leapfrog: 2nd order when only electric fields are relevant (poisson solver)
    * Boris-Leapfrog: 2nd order for electric and magnetic fields (poisson solver)
    * Higuera-Cary: 2nd order for electric and magnetic fields (poisson solver)
    * RK3: Runge-Kutta 3rd order in time
    * RK4: Runge-Kutta 4th order in time
    * RK14: Low storage Runge-Kutta 4, 14 stages version - Niegemann et al 2012
    * DSMC: Direct Simulation Monte Carlo, Section {ref}`sec:DSMC`
    * FP-Flow: Fokker-Planck-based collision operator, Section {ref}`sec:FP-Flow`
    * BGK-Flow: Bhatnagar-Gross-Krook collision operator, Section {ref}`sec:BGK-Flow`
* ``PICLAS_EQNSYSNAME``: Equation system to be solved
    * maxwell:
    * poisson:
* ``PICLAS_POLYNOMIAL_DEGREE``: Defines the polynomial degree of the solution. The order of convergence follows as $N+1$. Each grid
cell contains $(N+1)^3$ collocation points to represent the solution.
* ``PICLAS_NODETYPE``: The nodal collocation points used during the simulation
    * GAUSS: Legendre-Gauss distributed nodes
    * GAUSS-LOBATTO: Legendre-Gauss-Lobatto distributed nodes
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
* The order of defined variables is irrelevant, except for the special case when redefining boundaries.
However, it is preferable to group similar variables together.

The options and underlying models are discussed in Chapter {ref}`userguide/features-and-models/index:Features & Models`, while the available
output options are given in Chapter {ref}`userguide/visu_output:Visualization & Output`.
Due to the sheer number of parameters available, it is advisable to build upon an existing parameter file from one of the tutorials
in Chapter {ref}`userguide/tutorials/index:Tutorials`.

## Simulation

After the mesh generation, compilation of the binary and setup of the parameter files, the code can be executed by

    piclas parameter.ini [DSMC.ini]

The simulation may be restarted from an existing state file

    piclas parameter.ini [DSMC.ini] [restart_file.h5]

The state file , e.g., TestCase_State_000.5000000000000000.h5, contains all the required information to continue the simulation from
this point in time

    piclas parameter.ini DSMC.ini TestCase_State_000.5000000000000000.h5

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

### Profile-guided optimization (PGO)

To further increase performance for production runs, profile-guided optimization can be utilized with the GNU compiler. This requires the execution of a representative simulation run with PICLas compiled using profiling instrumentation. For this purpose, the code has to be configured and compiled using the following additional settings and the `Profile` build type:

    -DPICLAS_PERFORMANCE=ON -DUSE_PGO=ON -DCMAKE_BUILD_TYPE=Profile

A short representative simulation has to be performed, where additional files with the profiling information will be stored. Note that the test run should be relatively short as the code will be substantially slower than the regular `Release` build type. Afterwards, the code can be configured and compiled again for the production runs, using the `Release` build type:

    -DPICLAS_PERFORMANCE=ON -DUSE_PGO=ON -DCMAKE_BUILD_TYPE=Release

Warnings regarding missing profiling files (`-Wmissing-profile`) can be ignored, if they concern modules not relevant for the current simulation method (e.g. `bgk_colloperator.f90` will be missing profile information if only a DSMC simulation has been performed).

## Post-processing

**PICLas** comes with a tool for visualization. The piclas2vtk tool converts the HDF5 files generated by **PICLas** to the binary
VTK format, readable by many visualization tools like ParaView and VisIt. The tool is executed by

~~~~~~~
piclas2vtk [posti.ini] output.h5
~~~~~~~

Multiple HDF5 files can be passed to the piclas2vtk tool at once. The (optional) runtime parameters to be set in `posti.ini` are
given in Chapter {ref}`userguide/visu_output:Visualization & Output`.
