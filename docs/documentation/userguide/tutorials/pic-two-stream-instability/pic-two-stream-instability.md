(sec:tutorial-pic-two-stream-instability)=
# Two Stream Instability (PIC, Poisson's Equation)

This tutorial considers a one-dimensional (1D) setup of the so-called two-stream instability (TSI) using PICLas.
The TSI is a key phenomenon in plasma physics, occurring when two counter-flowing and interpenetrating streams of
charged particles, such as electrons and ions, travel with different (opposite) velocities.
This relative motion induces perturbations that grow exponentially, leading to the amplification of plasma waves and potentially
resulting in turbulent behaviour.

In PICLas, this process can be simulated with the electrostatioc PIC-Poisson solver.
The setup simulates two electron streams moving with opposite velocities $\nu_{1}$=$-\nu_{2}$ of identical densities $n_{1}$=$n_{2}$,
interacting with each other, leading to an instability that generates a fluctuating electric field.
The particle phase-space plot at the end of the simulation results in an eye-like structure.

Begin the tutorial, by copying the *pic-poisson-TSI* directory from the tutorial folder in the top level
directory of the piclas repository to a separate location

    cp -r $PICLAS_PATH/tutorials/pic-poisson-TSI .
    cd pic-poisson-TSI

where the variable `PICLAS_PATH` contains the path to the location of the piclas repository.

## Mesh Generation with HOPR (pre-processing)

In a first step, the mesh file in HDF5 format is created.
The mesh files used by **piclas** are created by supplying an input file *hopr.ini* with the required information for a mesh that
has either been created by an external mesh generator or directly from block-structured information in the *hopr.ini* file itself.
Here, a block-structured grid is created directly from the information in the hopr.ini file.
To create the *.h5* mesh file, run

    hopr hopr.ini

This creates the mesh file *two_stream_instability_mesh.h5* in HDF5 format.
The size of the simulation domain is set to [$4\pi\times0.03\times0.03$] m$^{3}$ and is defined by the single block information
in the line, where each node of the hexahedral element is defined

    Corner = (/0.,0.,0.,,12.566370614359,0.,0.,,12.566370614359,0.03,0.,,0.,0.03,0.,,0.,0.,0.03,,12.566370614359,0.,0.03,,12.566370614359,0.03,0.03,,0.,0.03,0.03/)

The number of mesh elements for the block in each direction can be adjusted by changing the line

    nElems = (/801,1,1/) ! number of elements in each direction (x,y,z)

Each side of the block has to be assigned a boundary index, which corresponds to the boundaries defined in the next steps

    BCIndex = (/5,3,2,4,1,6/)

The field boundaries can directly be defined in the hopr.ini file (contrary to the particle boundary conditions, which are defined
in the parameter.ini).
Periodic boundaries always have to be defined in the hopr.ini.

    !=============================================================================== !
    ! BOUNDARY CONDITIONS
    !=============================================================================== !
    BoundaryName = BC_periodicx+ ! Periodic (+vv1)
    BoundaryType = (/1,0,0,1/)   ! Periodic (+vv1)
    BoundaryName = BC_periodicx- ! Periodic (-vv1)
    BoundaryType = (/1,0,0,-1/)  ! Periodic (-vv1)
    BoundaryName = BC_periodicy+ ! Periodic (+vv2)
    BoundaryType = (/1,0,0,2/)   ! Periodic (+vv2)
    BoundaryName = BC_periodicy- ! Periodic (-vv2)
    BoundaryType = (/1,0,0,-2/)  ! Periodic (-vv2)
    BoundaryName = BC_periodicz+ ! Periodic (+vv3)
    BoundaryType = (/1,0,0,3/)   ! Periodic (+vv3)
    BoundaryName = BC_periodicz- ! Periodic (-vv3)
    BoundaryType = (/1,0,0,-3/)  ! Periodic (-vv3)

    VV = (/12.566370614359 , 0.  , 0./)    ! Displacement vector 1 (x-direction)
    VV = (/0.     , 0.03 , 0./)            ! Displacement vector 2 (y-direction)
    VV = (/0.     , 0.  , 0.03/)           ! Displacement vector 3 (z-direction)

In this case, a fully periodic setup is chosen by defining periodic boundaries on all six sides of the block, reflecting each
positive and negative Cartesian coordinate. In $x$-direction,

    BoundaryName = BC_periodicx+ ! Periodic (+vv1)
    BoundaryType = (/1,0,0,1/)   ! Periodic (+vv1)
    BoundaryName = BC_periodicx- ! Periodic (-vv1)
    BoundaryType = (/1,0,0,-1/)  ! Periodic (-vv1)

where for each of the six boundaries, a name `BoundaryName` and a type `BoundaryType` must be defined (in this order).
The boundary name can be chosen by the user and may be used again in the parameter.ini.
The first "1" in `BoundaryType` corresponds to the type "periodic" and the last entry, here, either "1" or "-1" corresponds to the
first periodic vector that is defined via `VV=(/12.566370614359 , 0.  , 0./)` that handles periodicity in the x-direction and gives the
orientation on the boundary for the vector. Note that each periodic boundary must have one positive and one negative corresponding
boundary for the same periodic vector.


## PIC Simulation with PICLas

Install **piclas** by compiling the source code as described in Chapter {ref}`userguide/installation:Installation`, specifically
described under Section {ref}`userguide/installation:Compiling the code`.
Always build the code in a separate directory located in the piclas top level directory.
For this PIC tutorial, e.g., create a directory *build_poisson_Leapfrog* by running

    cd $PICLAS_PATH

where the variable `PICLAS_PATH` contains the path to the location of the piclas repository.
If the piclas repository is located in the home directory, runlasma Wave

    cd /home/$(whoami)/piclas

and the create the build directory, in which the compilation process will take place

    mkdir build_poisson_Leapfrog

and the directory structure, which can be viewed via

    ls -l

should look like this

     build_poisson_RK3
     cmake
     CMakeListsLib.txt
     CMakeListsMachine.txt
     CMakeLists.txt
     CONTRIBUTORS.md
     docs
     LICENCE.md
     README.md
     REFERENCE.md
     REGGIE.md
     regressioncheck
     share
     SpeciesDatabase.h5
     src
     tools
     tutorials
     unitTests

Always compile the code within the *build* directory, hence, navigate to the *build_poisson_Leapfrog* directory before running cmake

    cd build_poisson_Leapfrog

For this specific tutorial, make sure to set the correct compile flags

    PICLAS_EQNSYSNAME     = poisson
    LIBS_USE_PETSC        = ON
    PICLAS_TIMEDISCMETHOD = Leapfrog


using the ccmake (gui for cmake) or simply run the following command from inside the *build* directory

    cmake ../ -DLIBS_USE_PETSC=0N -DPICLAS_EQNSYSNAME=poisson -DPICLAS_TIMEDISCMETHOD=Leapfrog

to configure the build process and run

    make -j

afterwards to compile the executable.
For this setup, the Poisson solver is used with the Leapfrog time discretization method.
An overview over the available solver and discretization options is given in Section {ref}`sec:solver-settings`.
To run the simulation, the *piclas* binary is used and to visualise the resulting *.h5* output files, the post-processing
too *piclas2vtk* is required to convert them into the standard *.vtu* format.

To avoid having to use the absolute file path of the executables, an alias or a symbolic link may be created.
To create a symbolic link within the directory of the tutorial, run

    ln -s $PICLAS_PATH/build_poisson_Leapfrog/bin/piclas
    ln -s $PICLAS_PATH/build_poisson_Leapfrog/bin/piclas2vtk

where the variable `PICLAS_PATH` contains the path to the location of the piclas repository.
If the piclas repository is located in the home directory, the two commands

    ln -s /home/$(whoami)/piclas/build_poisson_Leapfrog/bin/piclas
    ln -s /home/$(whoami)/piclas/build_poisson_Leapfrog/bin/piclas2vtk

can directly be executed without needing to modify them.

Always check where piclas is located before running the commands to ensure that the correct solver settings have been compiled into
the executable.

### Overview

The simulation setup is defined in *parameter.ini*.
The process of the TSI corresponds to a transfer of initial transversal velocities $\nu_{tr,x}(t=0)=\nu_{1}(t=0)$, into thermal
velocities $\nu_{th,x}$.
The transversal velocity in x-direction is computed as

$$\nu_{tr,x}={\frac{1}{N_{Parts}}\sum_{i=1}^{N_{Parts}}\nu_{x,i}}~,$$

which in term is used to compute the thermal velocity in x-direction as

$$\nu_{th,x}={\frac{1}{N_{Part}-1}\sum_{i=1}^{N_{Parts}}(\nu_{x,i}-\nu_{tr,x})²}~.$$

The initial velocity $\nu_{tr,x}$ is computed as

$$\nu_{tr,x}={}\frac{3}{\sqrt{2}}\omega_{p}~$$

For a specific electron number density, the (cold) plasma frequency of the system is

$$\omega_{p}=\nu_{e}=\sqrt{\frac{e^{2}n_{e}}{\varepsilon_{0}m_{e}}}~,$$

which is the frequency with which the charge density of the electrons oscillates, where
$\varepsilon_{0}$ is the permittivity of vacuum, $e$ is the elementary charge, $n_{e}$ and $m_{e}$
are the electron density and mass, respectively.
For a standard explicit PIC method, the plasma frequency yields the smallest time step that has to be resolved numerically. The
Debye length

$$\lambda_{D}=\sqrt{\frac{\varepsilon_{0}k_{B}T_{e}}{e^{2}n_{e}}}~,$$

where $\varepsilon_{0}$ is the permittivity of vacuum, $k_{B}$ is the Boltzmann constant, $e$ is the
elementary charge and $T_{e}$ and $n_{e}$ are the electron temperature and density, respectively, gives a spatial resolution
constraint.
The restriction on the spatial resolution is the number of elements (and polynomial degree $N$) that are required to resolve
the physical properties of the PIC simulation. If the temporal and/or spatial constraints are violated, the simulation cannot not yield
physical results and might result in a termination of the simulation.


The physical parameters for this test case are summarized in {numref}`tab:pic-two-stream-instability`.

```{table} Physical properties
---
name: tab:pic-two-stream-instability
---
|             Property            |            Value           |
| ------------------------------- |  :-----------------------: |
|      electron mass $m_{e}$      |   $\pu{9.1093826E-31 kg}$  |
|     electron charge $q_{e}$     | $\pm\pu{1.60217653E-19 C}$ |
|         ion mass $m_{i}$        |  $\pu{1.672621637E-27 kg}$ |
|        ion charge $q_{i}$       | $\pm\pu{3.20435306E-19 C}$ |
```

### Numerical setup

The general numerical parameters (defined in *parameter.ini*) are set via

    ! =============================================================================== !
    ! DISCRETIZATION
    ! =============================================================================== !
    N = 2 ! Polynomial degree of the DG method (field solver)

    ! =============================================================================== !
    ! MESH
    ! =============================================================================== !
    MeshFile = two_stream_instability_mesh.h5 ! Relative path to the mesh .h5 file

    ! =============================================================================== !
    ! General
    ! =============================================================================== !
    ProjectName       = two_stream_instability ! Project name that is used for naming state files
    ColoredOutput     = F           ! Turn ANSI terminal colors ON/OFF
    doPrintStatusLine = T           ! Output live of ETA
    TrackingMethod    = TriaTracking

where, among others, the polynomial degree $N$, the path to the mesh file `MeshFile`, project name and the option to print the ETA
to the terminal output in each time step.

The temporal parameters of the simulation are controlled via

    ! =============================================================================== !
    ! CALCULATION
    ! =============================================================================== !
    ManualTimeStep  = 1e-10 ! Fixed pre-defined time step only when using the Poisson solver. Maxwell solver calculates dt that considers the CFL criterion
    tend            = 1e-6  ! Final simulation time
    Analyze_dt      = 1E-8  ! decrease analyze_dt for better resolution of fish eye
    IterDisplayStep = 1000  ! Number of iterations between terminal output showing the current time step iteration
    TimeStampLength = 13    ! Reduces the length of the timestamps in filenames for better postprocessing

where the time step for the field and particle solver is set via `ManualTimeStep`, the final simulation time `tend`, the time
between restart/checkpoint file output `Analyze_dt` (also the output time for specific analysis functions) and the number of time
step iterations `IterDisplayStep` between information output regarding the current status of the simulation that is written to std.out.
The format of the restart/checkpoint files, which are created every time the simulation time reaches a multiple of
`Analyze_dt`, is set via `TimeStampLength = 13`, the names of the output files are shortened for better postprocessing.
This is done, because some visualisation tools might otherwise incorrectly group the files when opening a complete set of output files.
The remaining parameters are selected for the field and particle solver as well as run-time analysis.

#### Boundary conditions

As there are no walls present in the setup, all boundaries are set as periodic boundary conditions for the field as well as the
particle solver. The particle boundary conditions are set by the following

    ! =============================================================================== !
    ! PARTICLE Boundary Conditions
    ! =============================================================================== !
    Part-nBounds              = 6             ! Number of particle boundaries
    Part-Boundary1-SourceName = BC_periodicx+ ! Name of 1st particle BC
    Part-Boundary1-Condition  = periodic      ! Type of 1st particle BC
    Part-Boundary2-SourceName = BC_periodicx- ! ...
    Part-Boundary2-Condition  = periodic      ! ...
    Part-Boundary3-SourceName = BC_periodicy+ ! ...
    Part-Boundary3-Condition  = periodic      ! ...
    Part-Boundary4-SourceName = BC_periodicy- ! ...
    Part-Boundary4-Condition  = periodic      ! ...
    Part-Boundary5-SourceName = BC_periodicz+ ! ...
    Part-Boundary5-Condition  = periodic      ! ...
    Part-Boundary6-SourceName = BC_periodicz- ! ...
    Part-Boundary6-Condition  = periodic      ! ...

    Part-nPeriodicVectors = 3 ! Number of periodic boundary (particle and field) vectors

    Part-FIBGMdeltas = (/12.56637061435917295 , 0.03 , 0.03/) ! Cartesian background mesh (bounding box around the complete simulation domain)
    Part-FactorFIBGM = (/801., 1. , 1./)   ! Division factor that is applied t the "Part-FIBGMdeltas" values to define the dx, dy and dz distances of the Cartesian background mesh

where, the number of boundaries `Part-nBounds` (6 for a 3D cuboid domain) is followed by the names of
the boundaries (given by the *hopr.ini* file) and the type `periodic`.
Furthermore, the periodic vectors must be supplied and the size of the Cartesian background mesh `Part-FIBGMdeltas`, which can be
accompanied by a division factor (that is the number of background cells) in each direction given by `Part-FactorFIBGM`.
Here, the size and number of cells of the background mesh correspond to the actual mesh.

#### Field solver

The settings for the field solver (HDGSEM) are given by

    ! =============================================================================== !
    ! Field Solver: HDGSEM
    ! =============================================================================== !
    IniExactFunc = 0    ! Initial field condition. 0: zero solution vector
    PrecondType  = 10   ! Direct numerical method

For this tutorial the PETSc solver is used where a multitude of different numerical methods to solve the resulting system of linear
equations is given by the implemented PETSc library are available.
The direct solver of PETSc is used by setting `PrecondType=10`.
The parameter `IniExactFunc` sets the initial solution of the field solver and a value of zero means that the electric potential is
initialised with zero everywhere.

The PIC parameters for interpolation (of electric fields to the particle positions) and deposition (mapping of charge properties
from particle locations to the grid) are selected via

    ! =============================================================================== !
    ! PIC: Interpolation/Deposition
    ! =============================================================================== !
    PIC-DoInterpolation  = T                   ! Activate Lorentz forces acting on charged particles
    PIC-DoDeposition     = T                   ! Activate particle charge deposition
    PIC-Deposition-Type  = cell_volweight_mean ! Linear deposition method

Electro(-magnetic) forces that act on charged particles and accelerate these are activated by setting `PIC-DoInterpolation=T`.
The deposition of charges to the grid, which is required to consider source terms in the field solver, is activated via `PIC-DoDeposition=T`.
For charge and current deposition, `PIC-Deposition-Type = cell_volweight_mean` is selected.
Note that current deposition is only performed as post-processing output when using the Poisson solver in PICLas, as this solver
only required the charge deposition as input source terms.
The different available deposition types are described in more detail in Section {ref}`sec:PIC-deposition`.

#### Particle solver

The numerical scheme for tracking the movement of all particles throughout the simulation domain can be switched by

    ! =============================================================================== !
    ! Particle Solver
    ! =============================================================================== !
    TrackingMethod = refmapping ! [INT/STR] Define Method that is used for tracking of particles:
                                ! refmapping (1): reference mapping of particle position with (bi-)linear and
                                ! bezier (curved) description of sides.
                                ! tracing (2): tracing of particle path with (bi-)linear and bezier (curved)
                                ! description of sides.
                                ! triatracking (3): tracing of particle path with triangle-aproximation of
                                ! (bi-)linear sides.

The number of particle species `Part-nSpecies` that are used in the simulation (created initially or during the simulation time
through chemical reactions) defines the number of subsequent parameters that will be defined in the parameter input file

    ! =============================================================================== !
    ! PARTICLE Emission
    ! =============================================================================== !
    Part-nSpecies             = 3    ! Number of particle species

The creation (sometimes labelled emission or initialisation) of particles at the beginning or during the simulation
is controlled via the following parameters. Here, the parameters for the electrons are shown, however, the parameters for the ions
are set analogously and included in the supplied parameter.ini.
For each species, the mass (`Part-SpeciesX-MassIC`), charge (`Part-SpeciesX-ChargeIC`) and weighting factor (`Part-SpeciesX-MacroParticleFactor`)
have to be defined.

    ! -------------------------------------
    ! Electrons 1
    ! -------------------------------------
    Part-Species1-ChargeIC            = -1.60217653E-19 ! Electric charge of species #1
    Part-Species1-MassIC              = 9.1093826E-31   ! Rest mass of species #1
    Part-Species1-MacroParticleFactor = 2.4e5           ! Weighting factor for species #1
    Part-Species1-nInits              = 1               ! Number of initialization/emission regions for species #1

The number of initialization sets is defined by `Part-Species1-nInits`, where each initialization set is accompanied
by a block of parameters that starts from `Part-Species1-Init1-SpaceIC` up to `Part-Species1-Init1-VeloVecIC` and are preceded by the
corresponding `-InitX` counter. In this example we have a single initialization set per species definition.
The `Part-Species1-Init1-SpaceIC =  sin_deviation` flag defines the type of the initialization set, here, the distribution the particles
equidistantly on a line and sinusoidally dislocates them (representing an initial stage of a plasma wave in 1D).
Each type of the initialization set might have a different set of parameters and an overview is given in Section
{ref}`sec:particle-initialization-and-emission`.


    Part-Species1-Init1-ParticleNumber        = 100000           ! Number of simulation particles for species #1 and initialization #1
    Part-Species1-Init1-maxParticleNumber-x   = 100000           ! Number of simulation particles in x-direction for species #1 and initialization #1
    Part-Species1-Init1-SpaceIC               = sin_deviation    ! Sinusoidal distribution is space
    Part-Species1-Init1-BasePointIC           = (/0.0 , 0.015 , 0.015/)
    Part-Species1-Init1-velocityDistribution  = constant         ! Constant velocity distribution
    Part-Species1-Init1-maxParticleNumber-y   = 1                ! Number of particles in y
    Part-Species1-Init1-maxParticleNumber-z   = 1                ! Number of particles in z
    Part-Species1-Init1-MWTemperatureIC       = 0.72429730341e23
    Part-Species1-Init1-Amplitude             = 0.004            ! Specific factor for the sinusoidal distribution is space
    Part-Species1-Init1-WaveNumber            = 0.5              ! Specific factor for the sinusoidal distribution is space
    Part-Species1-Init1-VeloIC                = 1.06e8           ! Velocity magnitude [m/s]
    Part-Species1-Init1-VeloVecIC             = (/1.,0.001,0.001/)  ! Normalized velocity vector

The extent of dislocation is controlled by `Part-Species1-Init1-Amplitude`, which is only set for the electron species as the ion
species is not dislocated (they remain equidistantly distributed).
The parameter `Part-Species1-Init1-WaveNumber` sets the number of sine wave repetitions in the $x$-direction of the domain.
In case of the `SpaceIC=sin_deviation`, the number of simulation particles must be equal to the multiplied values given in
`Part-Species1-Init1-maxParticleNumber-x/y/z` as this emission type allows distributing the particles not only in one, but in all
three Cartesian coordinates, which is not required for this 1D example.

In the same manner, a second stream of electron is defined with  only the  difference of `Part-Species2-Init1-VeloVecIC` = (/-1.,0.001,0.001/) and `Part-Species2-Init1-Amplitude` = -0.004

The ions at rest in the background are defined as well providing a quasi-neutral state.

    ! ------------------------------------------------------------------------------- !
    ! Ions 3
    ! ------------------------------------------------------------------------------- !
    Part-Species3-ChargeIC                     = 3.204e-19
    Part-Species3-MassIC                       = 6.69e-27
    Part-Species3-MacroParticleFactor          = 2.4e5            ! Weighting factor for species #3
    Part-Species3-nInits                       = 1
    !Part-Species3-Init1-maxParticleNumber-x   = 100000           ! Number of simulation particles in x-direction for species 3
    Part-Species3-Init1-maxParticleNumber-y    = 1                ! Number of particles in y
    Part-Species3-Init1-maxParticleNumber-z    = 1                ! Number of particles in z
    Part-Species3-Init1-SpaceIC                = sin_deviation
    Part-Species3-Init1-BasePointIC            = (/0.0 , 0.015 ,0.015/)
    Part-Species3-Init1-velocityDistribution   = constant
    Part-Species3-Init1-VeloIC                 = 0.0
    Part-Species3-Init1-VeloVecIC              = (/0.,0.,0./)
    Part-Species3-Init1-Amplitude              = 0.0             ! Specific factor for the sinusoidal distribution is space
    Part-Species3-Init1-WaveNumber             = 0.              ! Specific factor for the sinusoidal distribution is space

#### Analysis setup

For run-time analysis, piclas offers a multitude of possible settings to activate.

    ! =============================================================================== !
    ! Analysis
    ! =============================================================================== !
    CalcCharge               = T ! writes rel/abs charge error to PartAnalyze.csv
    CalcPotentialEnergy      = T ! writes the potential field energy to FieldAnalyze.csv
    CalcKineticEnergy        = T ! writes the kinetic energy of all particle species to PartAnalyze.csv
    PIC-OutputSource         = T ! writes the deposited charge (RHS of Poisson's equation to XXX_State_000.0000XXX.h5)
    CalcPICTimeStep          = T ! writes the PIC time step restriction to XXX_State_000.0000XXX.h5 (rule of thumb)
    CalcPointsPerDebyeLength = T ! writes the PIC grid step restriction to XXX_State_000.0000XXX.h5 (rule of thumb)
    CalcTotalEnergy          = T ! writes the total energy of the system to PartAnalyze.csv (field and particle)
    CalcInternalEnergy       = T
    CalcTemp                 = T
    CalcVelos                = T
    Part-LorentzType         = 3
    CalcTimeAverage          = T
    CalcPlasmaParameter      =T
    VarNameAvg               = ChargeDensity-Spec01
    VarNameAvg               = ChargeDensity-Spec02
    VarNameAvg               = ChargeDensity-Spec03

The function of each parameter is given in the code comments. Information regarding every parameter can be obtained from
running the command

    piclas --help "CalcCharge"

where each parameter is simply supplied to the *help* module of **piclas**. This help module can also output the complete set of
parameters via `piclas --help` or a subset of them by supplying a section, e.g., `piclas --help "HDG"` for the HDGSEM solver.

### Running the code

The command

    ./piclas parameter.ini | tee std.out

executes the code and dumps all output into the file *std.out*.
To reduce the computation time, the simulation can be run using the Message Passing Interface (MPI) on multiple cores, in this case 8

    mpirun -np 8 piclas parameter.ini | tee std.out

If the run is successful, the contents of the working folder should look like this

    4.0K drwxrwxr-x  4.0K Jun 28 13:07 ./
    4.0K drwxrwxr-x  4.0K Jun 25 23:56 ../
    8.0K -rw-rw-r--  5.8K Jun 28 12:51 ElemTimeStatistics.csv
    120K -rw-rw-r--  113K Jun 28 12:51 FieldAnalyze.csv
    4.0K -rw-rw-r--  2.1K Jun 26 16:49 hopr.ini
    8.0K -rw-rw-r--  5.0K Jun 28 13:07 parameter.ini
    156K -rw-rw-r--  151K Jun 28 12:51 PartAnalyze.csv
     32K -rw-rw-r--   32K Jun 26 16:43 plasma_wave_mesh.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:44 two_stream_instability_State_.000000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:45 two_stream_instability_State_000.000000004.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:45 two_stream_instability_State_000.000000008.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:46 two_stream_instability_State_000.000000012.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:47 two_stream_instability_State_000.000000016.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:48 two_stream_instability_State_000.000000020.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:49 two_stream_instability_State_000.000000024.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:50 two_stream_instability_State_000.000000028.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:50 two_stream_instability_State_000.000000032.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 two_stream_instability_State_000.000000036.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 two_stream_instability_State_000.000000040.h5
    72K -rw-rw-r--   71K Jun 28 12:51 std.out

Multiple additional files should appear, which are named  *Projectname_State_Timestamp.h5* and contain the solution vector of the
equation system variables at each interpolation nodes at the given time, which corresponds to multiples of `Analyze_dt`.
If something goes wrong during the simulation, an output message should be displayed that is also written to the log file *std.out*
indicating what caused the crash.

A successful simulation should display the following lines at the end of the output:

    --------------------------------------------------------------------------------------------
     Sys date   :    03.07.2021 14:34:26
     PID: CALCULATION TIME PER TSTEP/DOF: [ 4.03465E-05 sec ]
     EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]: [ 4.16179E-07 sec/h ]
     Timestep   :    1.0000000E-10
    #Timesteps :    1.0000000E+04
     WRITE STATE TO HDF5 FILE [two_stream_instability_State_000.000001000.h5] ... DONE [ 0.06 sec ] [ 0:00:00:00 ]
    #Particles :    3.0000000E+05    Average particles per proc :    3.7500000E+04    Min :    3.1914000E+04    Max :    4.0600000E+04
     WRITE TIME AVERAGED STATE AND FLUCTUATIONS TO HDF5 FILE... DONE [ 0.01 sec ] [ 0:00:00:00 ]
    --------------------------------------------------------------------------------------------
    #Particles :    3.0000000E+05 (peak)         Average (peak) :    3.7500000E+04    Min :    2.1633000E+04    Max :    1.1188400E+05
    ============================================================================================
     PICLAS FINISHED! [ 1081.47 sec ] [ 0:00:18:01 ]
    ============================================================================================

## Visualization (post-processing)

To visualise the solution, the *State*-files (*landau_damping_State....h5*) must be converted into a format suitable for **ParaView**, **VisIt** or any other
visualisation tool.
The conversion if done using the program **piclas2vtk**, which is found at the same location as the **piclas** executable.
The parameters for **piclas2vtk** are also defined in the **parameter.ini** file under

    ! =============================================================================== !
    ! piclas2vtk
    ! =============================================================================== !
    NVisu         = 6! Polynomial degree used for the visualization when the .h5 file is converted to .vtu/.vtk format. Should be at least N+1
    VisuParticles = T  ! Activate the conversion of particles from .h5 to .vtu/.vtk format. Particles will be displayed as a point cloud with properties, such as velocity, species ID, etc.

where `NVisu` is the polynomial visualization degree on which the field solution is interpolated.
Depending on the used polynomial degree `N` and subsequently the degree of visualization `NVisu`, which should always be higher than
`N`, the resulting electric potential $\Phi$ and its derivative, the electric field strength **E** might show signs of oscillations.
Additionally, the flag `VisuParticles=T` activates the output of particle position, velocity and species to the *vtk*-files.
Running the command

    ./piclas2vtk parameter.ini plasma_wave_State_000.000000*

to generate the corresponding *vtk*-files, which can then be loaded into the visualisation tool.
The thermal and drift velocity plots can be viewed, for example with **ParaView**, by opening `PartAnalyze.csv` in the
**Line Chart View**. The graphs are depicted in {numref}`fig:instability-results`.

```{figure} results/velo.png
---
name: fig:instability-results
width: 700px
---

Transfer of transversal to thermal velocity over time.
```
To visualise the probability distribution function (PDF) in ParaView, a script is located in the tools folder of the piclas
repository under *piclas/tools/paraview/pdf/511*.
This is done by selecting **Tools &rarr; Manage Plugins &rarr; Load New &rarr; pdf_paraview511.xml**.
After loading the script, ensure that the *Auto Load* feature in the list of plugins is selected for the PDF to avoid loading the script
every time ParaView is opened.
Then load the particle file, for example, *two_stream_instability_visuPart_000.000000000.vtu* and apply the
`Probability Distribution Function 5.11 (piclas)` filter to this file.
This can be done by searching the feature directly in ParaView with the shortcut keys `CTRL+SPACE`.
Also adjust the spatial extent by setting `xMax`to 12.566, that is $4\pi$, which is the length of the domain in this tutorial.

An animation can be created from the collection of output files and is shown in {numref}`fig:PDF-gif`.

```{figure} results/compressed.gif
---
name: fig:PDF-gif
width: 700px
---

Temporal evolution of the PDF of the two-stream instability setup.