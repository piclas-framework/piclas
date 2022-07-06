# Plasma Wave (PIC, Poisson's Equation)

The setup considers a 1D plasma oscillation, which is a common and simple electrostatic PIC benchmark {cite}`Birdsall1991`,
{cite}`Hockney1988`,{cite}`Jacobs2006b`.
In PICLas it can be simulated either with the full Maxwell solver (DGSEM) or with the Poisson solver (HDGSEM), where the latter is
chosen for this this tutorial. In this setup, electrons oscillate around the almost immobile ions, which creates a fluctuating
electric field.

Before beginning with the tutorial, copy the `pic-poisson-plasma-wave` directory from the tutorial folder in the top level
directory to a separate location

    cp -r $PICLAS_PATH/tutorials/pic-poisson-plasma-wave .
    cd pic-poisson-plasma-wave

## Mesh Generation with HOPR (pre-processing)

Before the actual simulation is conducted, a mesh file in the correct HDF5 format has to be supplied.
The mesh files used by **piclas** are created by supplying an input file *hopr.ini* with the required information for a mesh that
has either been created by an external mesh generator or directly from block-structured information in the hopr.ini file itself.
Here, a block-structured grid is created directly from the information in the hopr.ini file.
To create the .h5 mesh file, simply run

    hopr hopr.ini

This creates the mesh file *plasma_wave_mesh.h5* in HDF5 format and is depicted in {numref}`fig:plasma-wave-mesh`.
Alternatively, if you do not want to run **hopr** yourself, you can also use the provided mesh.

The size of the simulation domain is set to [$2\pi\times0.2\times0.2$] m$^{3}$ and is defined by the single block information
in the line, where each node of the hexahedral element is defined

    Corner         =   (/0.,0.,0.,,6.2831,0.,0.,,6.2831, ... /)

The number of mesh elements for the block in each direction can be adjusted by changing the line

    nElems         = (/60,1,1/)                ! number of elements in each direction (x,y,z)

Each side of the block has to be assigned a boundary index, which corresponds to the boundaries defined in the next steps

    BCIndex        = (/5,3,2,4,1,6/)

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

    VV = (/6.2831 , 0.  , 0./)   ! Displacement vector 1 (x-direction)
    VV = (/0.     , 0.2 , 0./)   ! Displacement vector 2 (y-direction)
    VV = (/0.     , 0.  , 0.2/)  ! Displacement vector 3 (z-direction)

In this case a fully periodic setup is chosen by defining periodic boundaries on all six sides of the block, reflecting each
positive and negative Cartesian coordinate. In x-direction,

    BoundaryName = BC_periodicx+ ! Periodic (+vv1)
    BoundaryType = (/1,0,0,1/)   ! Periodic (+vv1)
    BoundaryName = BC_periodicx- ! Periodic (-vv1)
    BoundaryType = (/1,0,0,-1/)  ! Periodic (-vv1)

where for each of the six boundaries, a name `BoundaryName` and a type `BoundaryType` must be defined (in this order).
The boundary name can be chosen by the user and will be used again in the parameter.ini.
The first "1" in `BoundaryType` corresponds to the type "periodic" and the last entry, here, either "1" or "-1" corresponds to the
first periodic vector that is defined via `VV=(/6.2831 , 0.  , 0./)` that handles periodicity in the x-direction and gives the
orientation on the boundary for the vector. Note that each periodic boundary must have one positive and one negative corresponding
boundary for the same periodic vector.

```{figure} mesh/tut-pic-pw-mesh.jpg
---
name: fig:plasma-wave-mesh
---

Mesh with $60\times1\times1$ elements and a size of [$2\pi\times0.2\times0.2$] m$^{3}$.
```

## PIC Simulation with PICLas

Install **piclas** by compiling the source code as described in Chapter {ref}`userguide/installation:Installation` and make sure to set
the correct compile flags

    PICLAS_EQNSYSNAME     = poisson
    PICLAS_TIMEDISCMETHOD = RK3

or simply run the following command from inside the *build* directory

    cmake ../ -DPICLAS_EQNSYSNAME=poisson -DPICLAS_TIMEDISCMETHOD=RK3

to configure the build process and run `make` afterwards to build the executable. For this setup, we have chosen the Poisson solver
and selected the three-stage, third-order low-storage Runge-Kutta time discretization method. An overview over the available solver
and discretization options is given in Section {ref}`sec:solver-settings`. To run the simulation and analyse the results, the *piclas* and *piclas2vtk* executables have to be run. To avoid having to use the entire file path, you can either set aliases for both, copy them to your local tutorial directory or create a link to the files via.

    ln -s $PICLAS_PATH/build/bin/piclas
    ln -s $PICLAS_PATH/build/bin/piclas2vtk

The simulation setup is defined in *parameter.ini*. For a specific electron number density, the plasma frequency of the system is
given by

$$\omega_{p}=\omega_{e}=\sqrt{\frac{e^{2}n_{e}}{\varepsilon_{0}m_{e}}}~,$$

which is the frequency with which the charge density of the electrons oscillates, where
$\varepsilon_{0}$ is the permittivity of vacuum, $e$ is the elementary charge, $n_{e}$ and $m_{e}$
are the electron density and mass, respectively.
For the standard PIC method, the plasma frequency yields the smallest time step that has to be resolved numerically. The
Debye length

$$\lambda_{D}=\sqrt{\frac{\varepsilon_{0}k_{B}T_{e}}{e^{2}n_{e}}}~,$$

where $\varepsilon_{0}$ is the permittivity of vacuum, $k_{B}$ is the Boltzmann constant, $e$ is the
elementary charge and $T_{e}$ and $n_{e}$ are the electron temperature and density, respectively, gives a spatial resolution
constraint. In this test case, however, the electron temperature is not the defining factor for the spatial resolution because of the 1D nature
of the setup. Therefore, the resolution that is required is dictated by the gradient of the electric potential solution, i.e., the
electric field, which accelerates the charged particles and must be adequately resolved.
The restriction on the spatial resolution is simply the number of elements (and polynomial degree $N$) that are required to resolve
the physical properties of the PIC simulation. If the temporal and spatial constraints are violated, the simulation will not yield
physical results and might even result in a termination of the simulation.

The physical parameters for this test case are summarized in {numref}`tab:pic_poisson_plasma_wave_phys`.

```{table} Physical properties
---
name: tab:pic_poisson_plasma_wave_phys
---
|             Property            |            Value           |
| ------------------------------- |  :-----------------------: |
| electron number density $n_{e}$ |     $\pu{8e11 m^{-3}}$     |
|      electron mass $m_{e}$      |   $\pu{9.1093826E-31 kg}$  |
|    ion number density $n_{i}$   |     $\pu{8e11 m^{-3}}$     |
|         ion mass $m_{i}$        |  $\pu{1.672621637E-27 kg}$ |
|  electron/ion charge $q_{i,e}$  | $\pm\pu{1.60217653E-19 C}$ |
```

### General numerical setup

The general numerical parameters are selected by the following

    ! =============================================================================== !
    ! DISCRETIZATION
    ! =============================================================================== !
    N             = 5  ! Polynomial degree

    ! =============================================================================== !
    ! MESH
    ! =============================================================================== !
    MeshFile      = plasma_wave_mesh.h5
    useCurveds    = F

    ! =============================================================================== !
    ! General
    ! =============================================================================== !
    ProjectName      = plasma_wave
    Logging          = F
    WriteErrorFiles  = F
    TrackingMethod   = refmapping

where, among others, the polynomial degree $N$, the path to the mesh file `MeshFile`, project name and particle tracking method
`TrackingMethod` are chosen.

The temporal parameters of the simulation are controlled via

    ! =============================================================================== !
    ! CALCULATION
    ! =============================================================================== !
    ManualTimeStep  = 5e-10
    tend            = 40e-9
    Analyze_dt      = 4e-9
    IterDisplayStep = 50

where the time step for the field and particle solver is set via `ManualTimeStep`, the final simulation time `tend`, the time
between restart/checkpoint file output `Analyze_dt` (also the output time for specific analysis functions) and the number of time
step iterations `IterDisplayStep` between information output regarding the current status of the simulation that is written to std.out.
The remaining parameters are selected for the field and particle solver as well as run-time analysis.

### Boundary conditions

As there are no walls present in the setup, all boundaries are set as periodic boundary conditions for the field as well as the
particle solver. The particle boundary conditions are set by the following lines

    ! =============================================================================== !
    ! PARTICLES Boundary Conditions
    ! =============================================================================== !
    Part-nBounds              = 6
    Part-Boundary1-SourceName = BC_periodicx+
    Part-Boundary1-Condition  = periodic
    Part-Boundary2-SourceName = BC_periodicx-
    Part-Boundary2-Condition  = periodic
    Part-Boundary3-SourceName = BC_periodicy+
    Part-Boundary3-Condition  = periodic
    Part-Boundary4-SourceName = BC_periodicy-
    Part-Boundary4-Condition  = periodic
    Part-Boundary5-SourceName = BC_periodicz+
    Part-Boundary5-Condition  = periodic
    Part-Boundary6-SourceName = BC_periodicz-
    Part-Boundary6-Condition  = periodic

    Part-nPeriodicVectors = 3
    Part-PeriodicVector1  = (/6.2831,0.,0./)
    Part-PeriodicVector2  = (/0.,0.2,0./)
    Part-PeriodicVector3  = (/0.,0.,0.2/)

    Part-FIBGMdeltas = (/6.2831 , 0.2 , 0.2/)
    Part-FactorFIBGM = (/60     , 1   , 1/)

where, the number of boundaries `Part-nBounds` (6 in 3D cuboid) is followed by the names of
the boundaries (given by the hopr.ini file) and the type `periodic`. Furthermore, the periodic vectors must be supplied and the size
of the Cartesian background mesh `Part-FIBGMdeltas`, which can be accompanied by a division factor (i.e. number of background cells)
in each direction given by `Part-FactorFIBGM`. Here, the size and number of cells of the background mesh correspond to the actual mesh.

### Field solver

The settings for the field solver (HDGSEM) are given by

    ! =============================================================================== !
    ! HDGSEM
    ! =============================================================================== !
    epsCG                 = 1e-6
    maxIterCG             = 1000
    IniExactFunc          = 0

where `epsCG` sets the abort residual of the CG solver, `maxIterCG` sets the maximum number of iterations within the CG solver and
`IniExactFunc` set the initial solution of the field solver (here 0 says that nothing is selected).

The PIC parameters for interpolation (of electric fields to the particle positions) and deposition (mapping of charge properties
from particle locations to the grid) are selected via

    ! =============================================================================== !
    ! PIC: Interpolation/Deposition
    ! =============================================================================== !
    PIC-DoInterpolation       = T
    PIC-Interpolation-Type    = particle_position

    PIC-Deposition-Type         = shape_function_adaptive
    PIC-shapefunction-dimension = 1
    PIC-shapefunction-direction = 1
    PIC-shapefunction-alpha     = 4

where the interpolation type `PIC-Interpolation-Type = particle_position` is required and currently the only option. For deposition,
a polynomial shape function with the exponent `PIC-shapefunction-alpha` of the type `PIC-Deposition-Type = shape_function_adaptive` is
selected. The dimension `PIC-shapefunction-dimension`, here 1D and direction `PIC-shapefunction-direction`, are selected specifically
for the one-dimensional setup that is simulated here. The different available deposition types are described in more detail in
Section {ref}`sec:PIC-deposition`.

### Particle solver

For the treatment of particles, the maximum number of particles `Part-maxParticleNumber` that each processor can hold has to be supplied and
the number of particle species `Part-nSpecies` that are used in the simulation (created initially or during the simulation time
through chemical reactions).

    ! =============================================================================== !
    ! PARTICLE Emission
    ! =============================================================================== !
    Part-maxParticleNumber    = 4000
    Part-nSpecies             = 2

The inserting (sometimes labelled emission or initialization) of particles at the beginning or during the course of the simulation
is controlled via the following parameters. Here, only
the parameters for the electrons are shown, however, the parameters for the ions are set analogously and included in the supplied parameter.ini.
For each species, the mass (`Part-SpeciesX-MassIC`), charge (`Part-SpeciesX-ChargeIC`) and weighting factor (`Part-SpeciesX-MacroParticleFactor`)
have to be defined.

    ! -------------------------------------
    ! Electrons - Species 1
    ! -------------------------------------
    Part-Species1-MacroParticleFactor = 5e8
    Part-Species1-ChargeIC            = -1.60217653E-19
    Part-Species1-MassIC              = 9.1093826E-31

The number of initialization sets is defined by `Part-Species1-nInits`, where each initialization set is accompanied
by a block of parameters that starts from `Part-Species1-Init1-SpaceIC` up to `Part-Species1-Init1-VeloVecIC` and are preceded by the
corresponding `-InitX` counter. In this example we have a single initialization set per species definition.
The `Part-Species1-Init1-SpaceIC =  sin_deviation` flag defines the type of the initialization set, here, the distribution the particles
equidistantly on a line and sinusoidally dislocates them (representing an initial stage of a plasma wave in 1D).
Each type of the initialization set might have a different set of parameters and an overview is given in Section
{ref}`sec:particle-initialization-and-emission`.

    Part-Species1-nInits=1

    Part-Species1-Init1-SpaceIC               = sin_deviation
    Part-Species1-Init1-velocityDistribution  = constant
    Part-Species1-Init1-ParticleNumber        = 400
    Part-Species1-Init1-maxParticleNumber-x   = 400
    Part-Species1-Init1-maxParticleNumber-y   = 1
    Part-Species1-Init1-maxParticleNumber-z   = 1
    Part-Species1-Init1-Amplitude             = 0.01
    Part-Species1-Init1-WaveNumber            = 2.
    Part-Species1-Init1-VeloIC                = 0.
    Part-Species1-Init1-VeloVecIC             = (/1.,0.,0./)

To calculate the number of simulation particles of, e.g. electrons, defined by `Part-Species1-Init1-ParticleNumber`, the given
number density $n_{e}$ in {numref}`tab:pic_poisson_plasma_wave_phys`, the selected weighting factor $w_{e}$ and the volume of the 
complete domain ($V=2\pi\cdot0.2\cdot0.2\pu{m^{3}}$) are utilized.

$$ N_{e,sim} = \frac{n_{e} V}{w_{e}} $$

In this case, however, the number of particles are pre-defined and the weighting factor is derived from the above equation.
The extent of dislocation is controlled by `Part-Species1-Init1-Amplitude`, which is only set for the electron species as the ion
species is not dislocated (they remain equidistantly distributed).
The parameter `Part-Species1-Init1-WaveNumber` sets the number of sine wave repetitions in the `x`-direction of the domain.
In case of the `SpaceIC=sin_deviation`, the number of simulation particles must be equal to the multiplied values given in
`Part-Species1-Init1-maxParticleNumber-x/y/z` as this emission type allows distributing the particles not only in one, but in all
three Cartesian coordinates, which is not required for this 1D example.

### Analysis setup

Finally, some parameters for run-time analysis are chosen by setting them `T` (true). Further, with `TimeStampLength = 13`, the names of the output files are shortened for better postprocessing. If this is not done, e.g. Paraview does not sort the files correctly and will display faulty behaviour over time.

    ! =============================================================================== !
    ! Analysis
    ! =============================================================================== !
    TimeStampLength         = 13 ! Reduces the length of the timestamps in filenames for better postprocessing
    CalcCharge               = T ! writes rel/abs charge error to PartAnalyze.csv
    CalcPotentialEnergy      = T ! writes the potential field energy to FieldAnalyze.csv
    CalcKineticEnergy        = T ! writes the kinetic energy of all particle species to PartAnalyze.csv
    PIC-OutputSource         = T ! writes the deposited charge (RHS of Poisson's equation to XXX_State_000.0000XXX.h5)
    CalcPICTimeStep          = T ! writes the PIC time step restriction to XXX_State_000.0000XXX.h5 (rule of thumb)
    CalcPointsPerDebyeLength = T ! writes the PIC grid step restriction to XXX_State_000.0000XXX.h5 (rule of thumb)
    CalcTotalEnergy          = T ! writes the total energy of the system to PartAnalyze.csv (field and particle)

The function of each parameter is given in the code comments. Information regarding every parameter can be obtained from
running the command

    piclas --help "CalcCharge"

where each parameter is simply supplied to the *help* module of **piclas**. This help module can also output the complete set of
parameters via `piclas --help` or a subset of them by supplying a section, e.g., `piclas --help "HDG"` for the HDGSEM solver.

### Running the code

The command

    ./piclas parameter.ini | tee std.out

executes the code and dumps all output into the file *std.out*.
To reduce the computation time, the simulation can be run using the Message Passing Interface (MPI) on multiple cores, in this case 4
	
    mpirun -np 4 piclas parameter.ini | tee std.out

If the run has completed successfully, which should take only a brief moment, the contents of the working folder should look like

    4.0K drwxrwxr-x  4.0K Jun 28 13:07 ./
    4.0K drwxrwxr-x  4.0K Jun 25 23:56 ../
    8.0K -rw-rw-r--  5.8K Jun 28 12:51 ElemTimeStatistics.csv
    120K -rw-rw-r--  113K Jun 28 12:51 FieldAnalyze.csv
    4.0K -rw-rw-r--  2.1K Jun 26 16:49 hopr.ini
    8.0K -rw-rw-r--  5.0K Jun 28 13:07 parameter.ini
    156K -rw-rw-r--  151K Jun 28 12:51 PartAnalyze.csv
     32K -rw-rw-r--   32K Jun 26 16:43 plasma_wave_mesh.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:44 plasma_wave_State_000.000000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:45 plasma_wave_State_000.000000004.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:45 plasma_wave_State_000.000000008.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:46 plasma_wave_State_000.000000012.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:47 plasma_wave_State_000.000000016.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:48 plasma_wave_State_000.000000020.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:49 plasma_wave_State_000.000000024.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:50 plasma_wave_State_000.000000028.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:50 plasma_wave_State_000.000000032.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 plasma_wave_State_000.000000036.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 plasma_wave_State_000.000000040.h5
     72K -rw-rw-r--   71K Jun 28 12:51 std.out

Multiple additional files have been created, which are are named  **Projectname_State_Timestamp.h5**.
They contain the solution vector of the equation system variables at each interpolation nodes at the given time, which corresponds
to multiples of **Analyze_dt**. If these files are not present, something went wrong during the execution of **piclas**.
In that case, check the `std.out` file for an error message.

After a successful completion, the last lines in this file should look as shown below:

    --------------------------------------------------------------------------------------------
    Sys date  :    03.07.2021 14:34:26
    PID: CALCULATION TIME PER TSTEP/DOF: [ 5.85952E-05 sec ]
    EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]: [ 2.38587E-06 sec/h ]
    Timestep  :    5.0000000E-10
    #Timesteps :    8.0000000E+01
    WRITE STATE TO HDF5 FILE [plasma_wave_State_000.000000040.h5] ...DONE  [.008s]
    #Particles :    8.0000000E+02
    --------------------------------------------------------------------------------------------
    ============================================================================================
    PICLAS FINISHED! [           60.42 sec ] [     0:00:01:00]
    ============================================================================================

## Visualization (post-processing)

To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**, **VisIt** or any other
visualisation tool for which the program **piclas2vtk** is used.

The parameters for **piclas2vtk** are stored in the **parameter.ini** file under

    ! =============================================================================== !
    ! piclas2vtk
    ! =============================================================================== !
    NVisu         = 10
    VisuParticles = T

where `NVisu` is the polynomial visualization degree on which the field solution is interpolated.
Depending on the used polynomial degree `N` and subsequently the degree of visualization `NVisu`, which should always be higher than
`N`, the resulting electric potential $\Phi$ and its derivative the electric field strength **E** might show signs of oscillations.
This is because the PIC simulation is always subject to noise that is influenced by the discretization (number of elements and 
polynomial degree as well as number of particles) and is visible in the solution as this is a snapshot of the current simulation.

Additionally, the flag `VisuParticles` activates the output of particle position, velocity and species to the *vtk*-files.

Run the command

    ./piclas2vtk parameter.ini plasma_wave_State_000.000000*

to generate the corresponding *vtk*-files, which can then be loaded into the visualisation tool.

The electric potential field can be viewed, e.g., by opening `plasma_wave_Solution_000.000000040.vtu` and plotting the field
`Phi` along the x-axis, which should look like the following


```{figure} results/tut-pic-pw-results.jpg
---
name: fig:plasma-wave-results
---

Resulting electric potential and field.
```
















