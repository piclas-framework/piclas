# Landau Damping (PIC, Poisson's Equation)

Landau damping is a fundamental concept in plasma physics, where the amplitude of electrostatic waves decreases over time without energy dissipation through collisions. This effect, first described by Lev Landau in 1946, occurs due to the resonant interaction between the wave and particles in the plasma whose velocities are close to the wave's phase velocity.

In detail, Landau damping happens because some particles in the plasma move at velocities that resonate with the phase velocity of the electrostatic wave. These particles can either gain or lose energy from the wave depending on their initial velocities. If a greater number of particles gain energy (typically slower particles that are accelerated by the wave) than lose energy (faster particles that are decelerated by the wave), the wave loses energy, causing its amplitude to decrease over time. This process is purely kinetic and does not rely on collisions between particles.

The phenomenon is often analyzed using the Vlasov-Poisson equations, which describe the behavior of the plasma in terms of particle velocity distributions and electric fields. Advanced mathematical techniques, such as Laplace transforms and contour integration, are used to solve these equations and derive the damping rate {cite}`Finn2023` {cite}`Dawson1961`. An understanding of the phenomenon using simple derivations "which does not require the application of methods of complex analysis" is shown and discussed in {cite}`Shalchi2021`.


The general information needed to setup the simulation is given in the previous tutorial in {ref}`sec:tutorial-pic-poisson-plasma-wave`.

To begin with the tutorial, copy the *pic-poisson-landau-damping* directory from the tutorial folder in the top level
directory to a separate location

    cp -r $PICLAS_PATH/tutorials/pic-poisson-landau-damping .
    cd pic-poisson-landau-damping

where the variable `PICLAS_PATH` contains the path to the location of the piclas repository.

## Mesh Generation with HOPR (pre-processing)

First, the mesh file has to be created by running

    hopr hopr.ini

within the *pic-poisson-landau-damping* directory.

This creates the mesh file *landau_damping_mesh.h5* in HDF5 format.
The size of the simulation domain is [$4\pi\times1\times1$] m$^{3}$ and is defined by the single block information
in the line, where each node of the hexahedral element is defined

    Corner = (/0.,0.,0.,,12.5663706144,0.,0.,,12.5663706144, ... /)

The number of mesh elements for the block in each direction can be adjusted by changing the line

    nElems = (/30,1,1/) ! number of elements in each direction (x,y,z)

Each side of the block has to be assigned a boundary index, which corresponds to the boundaries defined in the next steps

    BCIndex = (/5,3,2,4,1,6/)


Periodic boundaries have to be defined in the hopr.ini via

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

    VV = (/12.5663706144 , 0.  , 0./)   ! Displacement vector 1 (x-direction)
    VV = (/0.            , 1   , 0./)   ! Displacement vector 2 (y-direction)
    VV = (/       0.     , 0.  , 1 /)   ! Displacement vector 3 (z-direction)


## PIC Simulation with PICLas

Install the required **piclas** executable by compiling the source code as described in Chapter {ref}`userguide/installation:Installation`, specifically
described under Section {ref}`userguide/installation:Compiling the code`.
Always build the code in a separate directory located in the piclas top level directory.
For this PIC tutorial, e.g., create a directory *build_poisson_Leapfrog* in the piclas repository by running

    cd $PICLAS_PATH

where the variable `$PICLAS_PATH` contains the path to the location of the piclas repository.
If the piclas repository is located in the home directory, simply run

    cd /home/$(whoami)/piclas

and the create the build directory, in which the compilation process will take place

    mkdir build_poisson_Leapfrog

and the directory structure, which can be viewed via

    ls -l

should look like this

     build_poisson_Leapfrog
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

    PICLAS_EQNSYSNAME       = poisson
    LIBS_USE_PETSC          = ON
    PICLAS_READIN_CONSTANTS = ON
    PICLAS_TIMEDISCMETHOD   = Leapfrog

which are forwarded to cmake by running the following command from inside the *build_poisson_Leapfrog* directory

    cmake ../ -DPICLAS_READIN_CONSTANTS=ON -DLIBS_USE_PETSC=0N -DPICLAS_EQNSYSNAME=poisson -DPICLAS_TIMEDISCMETHOD=Leapfrog

to configure the build process and run

    make -j

afterwards to compile the executable.
For this setup, the Poisson solver is used with the Leapfrog time discretization method.
An overview over the available solver and discretization options is given in Section {ref}`sec:solver-settings`.
To run the simulation, the *piclas* binary is used and to visualise the resulting *.h5* output files, the post-processing
too *piclas2vtk* is required to convert them into the standard *.vtu* format.

The compile flag `PICLAS_READIN_CONSTANT=ON` enables user-defined natural constants for the speed of light $c$, permittivity $\varepsilon_0$ and
permeability $\mu$ of vacuum, which must then be supplied in the parameter input file in this test case.
The physical constants used for defining the species properties (mass and charge) in this tutorial are also normalized.

To avoid having to use the absolute file path of the executables, an alias or a symbolic link may be created.
To create a symbolic link within the directory of the tutorial, run

    ln -s $PICLAS_PATH/build_poisson_Leapfrog/bin/piclas
    ln -s $PICLAS_PATH/build_poisson_Leapfrog/bin/piclas2vtk

where the variable `PICLAS_PATH` contains the path to the location of the piclas repository.
If the piclas repository is located in the home directory, the two commands

    ln -s /home/$(whoami)/piclas/build_poisson_Leapfrog/bin/piclas
    ln -s /home/$(whoami)/piclas/build_poisson_Leapfrog/bin/piclas2vtk

can directly be executed without needing to modify them.

Note that for a given electron number density, the plasma frequency of the system is given by

$$\omega_{p}=\omega_{e}=\sqrt{\frac{e^{2}n_{e}}{\varepsilon_{0}m_{e}}}~,$$

which is the frequency with which the charge density of the electrons oscillates, where
$\varepsilon_{0}$ is the permittivity of vacuum, $e$ is the elementary charge, $n_{e}$ and $m_{e}$
are the electron density and mass at rest, respectively.
For an explicit PIC method, the plasma frequency yields the smallest time step that has to be resolved numerically.
The Debye length

$$\lambda_{D}=\sqrt{\frac{\varepsilon_{0}k_{B}T_{e}}{e^{2}n_{e}}}~,$$

where $\varepsilon_{0}$ is the permittivity of vacuum, $k_{B}$ is the Boltzmann constant, $e$ is the
elementary charge and $T_{e}$ and $n_{e}$ are the electron temperature and density, respectively, sets the spatial resolution
constraint.
The restriction on the spatial resolution is therefore given by the number of elements (and polynomial degree $N$) that are required to resolve
the Debye length within the PIC simulation. If the temporal and spatial constraints are violated, the simulation cannot remain
stable over time.


### Numerical Setup
The input parameters for the simulation setup are defined in *parameter.ini* and the general numerical parameters are the following

    ! =============================================================================== !
    ! DISCRETIZATION
    ! =============================================================================== !
    N             = 4                       ! Polynomial degree of the DG method (field solver)

    ! =============================================================================== !
    ! MESH
    ! =============================================================================== !
    MeshFile      = landau_damping_mesh.h5  ! Relative path to the mesh .h5 file

    ! =============================================================================== !
    ! General
    ! =============================================================================== !
    ProjectName       = landau_damping      ! Project name that is used for naming state files
    doPrintStatusLine = T                   ! Output live of ETA

    c0  = 1.e8
    eps = 1
    mu  = 1.e-16

where the polynomial degree $N$ (results in a spatial order of convergence of $N+1$), the path to the mesh file
`MeshFile`, project name and the option to print the ETA to the terminal output in each time step.

As mentioned earlier, the compile flag `PICLAS_READIN_CONSTANT=ON` enables user-defined natural constants for the speed of light $c$,
permittivity $\varepsilon_0$ and permeability $\mu$ of vacuum, which are set in the parameter input file via `c0`, `eps`, and `mu`.

The temporal parameters of the simulation are controlled via

    ! =============================================================================== !
    ! CALCULATION
    ! =============================================================================== !
    ManualTimeStep  =  0.1 ! Fixed pre-defined time step only when using the Poisson solver.
    tend            = 50.0 ! Final simulation time
    Analyze_dt      =  1.0 ! Simulation time between analysis output to .h5
    IterDisplayStep = 10   ! Number of iterations between terminal output showing the current time step iteration
    TimeStampLength =  5   ! Reduces the length of the timestamps in filenames for better post-processing

where the time step for the field and particle solver is set via `ManualTimeStep`, the final simulation time `tend`, the time
between restart/checkpoint file output `Analyze_dt` (also the output time for specific analysis functions) and the number of time
step iterations `IterDisplayStep` between information output regarding the current status of the simulation that is written to std.out.
The format of the restart/checkpoint files, which are created every time the simulation time reaches a multiple of
`Analyze_dt`, is set via `TimeStampLength = 5`, the names of the output files are shortened for better postprocessing.
The remaining parameters are selected for the field and particle solver as well as run-time analysis.
This is done, because some visualisation tools might otherwise incorrectly group the files when opening a complete set of output
files.

#### Boundary conditions

As there are no walls present in the setup, all boundaries are set as periodic boundary conditions for the field as well as the
particle solver. The particle boundary conditions are set by the following lines

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

    Part-FIBGMdeltas = (/12.5663706144 , 1. , 1./) ! Cartesian background mesh (bounding box around the complete simulation domain)
    Part-FactorFIBGM = (/30     , 1   , 1/)   ! Division factor that is applied t the "Part-FIBGMdeltas" values to define the dx, dy and dz distances of the Cartesian background mesh

where, the number of boundaries `Part-nBounds` (6 in 3D cuboid) is followed by the names of
the boundaries (given by the hopr.ini file) and the type `periodic`. Furthermore, the periodic vectors must be supplied and the size
of the Cartesian background mesh `Part-FIBGMdeltas`, which can be accompanied by a division factor (that is the number of background cells)
in each direction given by `Part-FactorFIBGM`. Here, the size and number of cells of the background mesh correspond to the actual mesh.

#### Field solver

The settings for the field solver (HDGSEM) are given by

    ! =============================================================================== !
    ! Field Solver: HDGSEM
    ! =============================================================================== !
    epsCG        = 1e-12 ! Stopping criterion (residual) of iterative CG solver (default that is used for the HDGSEM solver)
    maxIterCG    = 1000  ! Maximum number of iterations
    IniExactFunc = 0     ! Initial field condition. 0: zero solution vector

where `epsCG` sets the abort residual and `maxIterCG` sets the maximum number of iterations within the CG solver or the PETSc
solver, depending on which of the two is used.
In this tutorial, the PETSc solver is used as the executable is compiled with `LIBS_USE_PETSC=ON`.
The parameter `IniExactFunc` sets the initial solution of the field solver and a value of zero means that the electric potential is
initialised with zero everywhere.

The PIC parameters for interpolation (of electric fields to the particle positions) and deposition (mapping of charge properties
from particle locations to the grid) are selected via

    ! =============================================================================== !
    ! PIC: Interpolation/Deposition
    ! =============================================================================== !
    PIC-DoInterpolation             = T              ! Activate Lorentz forces acting on charged particles
    PIC-DoDeposition                = T              ! Activate charge deposition to the grid
    PIC-Deposition-Type             = shape_function ! Particle-field coupling method. shape_function_adaptive determines the cut-off radius of the shape function automatically
    PIC-shapefunction-dimension     = 1              ! Sets the shape function 1D (the default is 3D)
    PIC-shapefunction-direction     = 1              ! Sets the axial direction of the 1D shape function (1:x, 2:y, 3:z)
    PIC-shapefunction-alpha         = 10             ! Sets the shape function exponent, which effectively scales the waist diameter of the shape function
    PIC-shapefunction-radius        = 0.5            ! Radius of influence for the shape function deposition method
    PIC-shapefunction-3D-deposition = F              ! Deposit the charge over volume (3D) is true or over a line (1D) or area (2D) if set false

Electro(-magnetic) forces that act on charged particles and accelerate these are activated by setting `PIC-DoInterpolation=T`.
The deposition of charges to the grid, which is required to consider source terms in the field solver, is activated via `PIC-DoDeposition=T`.
For the deposition method, the shape function is selected via `PIC-Deposition-Type=shape_function`, which has multiple additional
parameters to set the dimensionality and form of the shape function.
The dimension `PIC-shapefunction-dimension=1` and direction `PIC-shapefunction-direction=1` sets a 1D shape function in
x-direction and are specific to the one-dimensional setup that is simulated here.
The form of the shape function is adjusted by setting the radius of influence via `PIC-shapefunction-radius = 0.5` and the waist
radius via `PIC-shapefunction-alpha = 10`.
The parameter `PIC-shapefunction-3D-deposition` decides whether the charge is deposited in a volume or on a line/area, depending on
the dimensionality that is set for the shape function.
The different available deposition types are described in more detail in Section {ref}`sec:PIC-deposition` and piclas can display a
help section for the deposition methods by running

    ./piclas --help "PIC Deposition"

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
    Part-nSpecies = 2 ! Number of particle species

The creation (sometimes labelled emission or initialisation) of particles at the beginning or during the simulation
is controlled via the following parameters. Here, the parameters for the electrons are shown, however, the parameters for the ions
are set analogously and included in the supplied parameter.ini.
For each species, the mass (`Part-SpeciesX-MassIC`), charge (`Part-SpeciesX-ChargeIC`) and macro-particle weighting factor
(`Part-SpeciesX-MacroParticleFactor`) have to be defined.

    ! -------------------------------------
    ! Electrons 1
    ! -------------------------------------
    Part-Species1-ChargeIC            = -1.                ! Electric charge of species #1
    Part-Species1-MassIC              = 1.                 ! Rest mass of species #1
    Part-Species1-MacroParticleFactor = 3.14159265358e-4   ! Weighting factor for species #1
    Part-Species1-nInits              = 1                  ! Number of initialisation/emission regions for species #1

The number of initialisation sets is defined by `Part-Species1-nInits`, where each initialisation set is accompanied
by a block of parameters that starts from `Part-Species1-Init1-SpaceIC` up to `Part-Species1-Init1-VeloVecIC` and are preceded by the
corresponding `-InitX` counter. In this example we have a single initialisation set per species definition.
The `Part-Species1-Init1-SpaceIC =  cos_distribution` flag defines the type of the initialisation set, here, the distribution the particles
equidistantly on a line and dislocates them in a cosine pattern, representing an initial stage of a plasma wave in 1D.
Each type of the initialisation set might have a different set of parameters and an overview is given in Section
{ref}`sec:particle-initialization-and-emission`.

    Part-Species1-Init1-SpaceIC               = cos_distribution          ! Cosine distribution is space
    Part-Species1-Init1-ParticleNumber        = 40000                     ! Number of simulation particles for species #1 and initialisation #1
    Part-Species1-Init1-maxParticleNumber-x   = 40000                     ! Number of simulation particles in x-direction for species #1 and initialisation #1
    Part-Species1-Init1-velocityDistribution  = maxwell_distribution_1D   ! Constant velocity distribution
    Part-Species1-Init1-MWTemperatureIC       = 0.72429730341e23          ! Translational temprature
    Part-Species1-Init1-maxParticleNumber-y   = 1                         ! Number of particles in y
    Part-Species1-Init1-maxParticleNumber-z   = 1                         ! Number of particles in z
    Part-Species1-Init1-Amplitude             = 0.05                      ! Specific factor for the sinusoidal distribution is space
    Part-Species1-Init1-WaveNumber            = 0.5                       ! Specific factor for the sinusoidal distribution is space
    Part-Species1-Init1-VeloIC                = 0.                        ! Velocity magnitude [m/s]
    Part-Species1-Init1-VeloVecIC             = (/1.,0.,0./)              ! Normalized velocity vector


To calculate the number of simulation particles of, for example electrons, defined by `Part-Species1-Init1-ParticleNumber`, the given
number density , the selected weighting factor $w_{e}$ and the volume of the
complete domain ($V=4\pi\cdot1\cdot1\pu{m^{3}}$) are utilized.

$$ N_{e,sim} = \frac{n_{e} V}{w_{e}} $$

In this case, however, the number of particles are pre-defined and the weighting factor is derived from the above equation.
The extent of dislocation is controlled by `Part-Species1-Init1-Amplitude`, which is only set for the electron species as the ion
species is not dislocated (they remain equidistantly distributed).
The parameter `Part-Species1-Init1-WaveNumber` sets the number of cosine wave repetitions in the `x`-direction of the domain.
In case of the `SpaceIC=sin_deviation`, the number of simulation particles must be equal to the multiplied values given in
`Part-Species1-Init1-maxParticleNumber-x/y/z` as this emission type allows distributing the particles not only in one, but in all
three Cartesian coordinates, which is not required for this 1D example.


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


The function of each parameter is given in the code comments.
Information regarding every parameter can be obtained from running the command

    piclas --help "CalcCharge"

where each parameter is simply supplied to the *help* module of **piclas**.
This help module can also output the complete set of parameters via

    piclas --help

or a subset of them by supplying a section via 

    piclas --help "HDG"

for the HDGSEM solver.

### Running the code

The command

    ./piclas parameter.ini | tee std.out

executes the code and dumps all output into the file *std.out*.
For faster execution, piclas can be run in parallel via

    mpirun -np 8 ./piclas parameter.ini | tee std.out

which executes piclas by utilising 8 MPI processes (in general do not surpass the actual number of processes on the executing
system).

If the run is successful, the contents of the working folder should look like this

    8.0K -rw-rw-r--  5.8K Jun 28 12:51 ElemTimeStatistics.csv
    120K -rw-rw-r--  113K Jun 28 12:51 FieldAnalyze.csv
    4.0K -rw-rw-r--  2.1K Jun 26 16:49 hopr.ini
    8.0K -rw-rw-r--  5.0K Jun 28 13:07 parameter.ini
    156K -rw-rw-r--  151K Jun 28 12:51 PartAnalyze.csv
     32K -rw-rw-r--   32K Jun 26 16:43 landau_damping_mesh.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:44 landau_damping_State_000.000000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:45 landau_damping_State_000.000000004.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:45 landau_damping_State_000.000000008.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:46 landau_damping_State.000000012.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:47 landau_damping_State.000000016.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:48 landau_damping_State.000000020.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:49 landau_damping_State.000000024.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:50 landau_damping_State.000000028.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:50 landau_damping_State.000000032.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 landau_damping_State.000000036.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 landau_damping_State.000000040.h5
     72K -rw-rw-r--   71K Jun 28 12:51 std.out

Multiple additional files should appear, which are named  *Projectname_State_Timestamp.h5* and contain the solution vector of the
equation system variables at each interpolation nodes at the given time, which corresponds to multiples of `Analyze_dt`.
If something goes wrong during the simulation, an output message should be displayed that is also written to the log file *std.out*
indicating what caused the crash.

A successful simulation should display the following lines at the end of the output:

    --------------------------------------------------------------------------------------------
     Sys date  :    03.07.2021 14:34:26
     PID: CALCULATION TIME PER TSTEP/DOF: [ 1.21226E-04 sec ]
     EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]: [ 7.50886E+02 sec/h ]
     Timestep  :    1.0000000E-01
    #Timesteps :    5.0000000E+02
     WRITE STATE TO HDF5 FILE [landau_damping_State_050.0.h5] ... DONE [ 0.02 sec ] [ 0:00:00:00 ]
    #Particles :    4.0100000E+04    Average particles per proc :    1.0025000E+04    Min :    9.3390000E+03    Max :    1.0702000E+04
    --------------------------------------------------------------------------------------------
    #Particles :    4.0100000E+04 (peak)         Average (peak) :    1.0025000E+04    Min :    9.1730000E+03    Max :    1.1010000E+04
    ============================================================================================
     PICLAS FINISHED! [           60.42 sec ] [     0:00:01:00]
    ============================================================================================


## Visualisation (post-processing)

To visualise the solution, the *State*-files (*landau_damping_State....h5*) must be converted into a format suitable for **ParaView**, **VisIt** or any other
visualisation tool.
The conversion if done using the program **piclas2vtk**, which is found at the same location as the **piclas** executable.

The parameters for **piclas2vtk** are also defined in the **parameter.ini** file under

    ! =============================================================================== !
    ! piclas2vtk
    ! =============================================================================== !
    NVisu         = 10 ! Polynomial degree used for the visualisation when the .h5 file is converted to .vtu/.vtk format. Should be at least N+1
    VisuParticles = T  ! Activate the conversion of particles from .h5 to .vtu/.vtk format. Particles will be displayed as a point cloud with properties, such as velocity, species ID, etc.

where `NVisu` is the polynomial visualisation degree on which the field solution is interpolated.
Depending on the used polynomial degree `N`, the degree of visualisation `NVisu` should always be higher than
`N` because the PIC simulation is always subject to noise that is influenced by the discretization (number of elements and
polynomial degree as well as number of particles) and is visible in the solution results of the simulation.
Additionally, the flag `VisuParticles` activates the output of particle position, velocity and species index to the *vtk*-files.
Runining the command

    ./piclas2vtk parameter.ini landau_damping_State_000.0*

generates the corresponding *vtk*-files, which can then be loaded into the visualisation tool.
The analysis of the numerical results contained in *PartAnalyze.csv* and *FieldAnalyze.csv* using, for example, **ParaView** to display the data, shows details on
* gain in the kinetic energy of particles,
* decay in the electric field energy,
* total energy error converging over time.

The error in total energy, which is calculated from the sum of kinetic and potential electric energy (labelled `008-E-kin+pot` in *PartAnalyze.csv*),
is depicted in {numref}`fig:plasma-wave-res`, showing that the error converges to approximately 0.5 percent.

```{figure} results/land.png
---
name: fig:plasma-wave-res
width: 600px
---

Total energy error over time.
```

The electric field energy, labelled `"002-E-El"` in *FieldAnalyze.csv* is shown in {numref}`fig:plasma-wave-visual`.
Additionally, the analytical solution of the envelope of the electric field energy with a damping rate of
$\gamma = -0.153$ as derived in {cite}`Canosa1973` is displayed.
```{figure} results/damp.png
---
name: fig:plasma-wave-visual
width: 700px
---

Electric field energy over time in semi-log scale.
```
The relation between the electric field energy and the damping rate can be expressed via

$$ E_{el} = a \cdot e^{2\gamma t} $$

where the coefficient `a=0.0173` is found by fitting the envelope to the numerical results.
In ParaView, a *calculator* filter with the expression `0.0173*exp(2*(-0.153)*"001-time")` can be applied to the temporal data in
*FieldAnalyze.csv* to create the envelope plot for a direct comparison with the numerical result.
