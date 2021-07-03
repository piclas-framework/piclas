## PIC Plasma Wave (Poisson's Equation)
\label{sec:tut_freestream}

The setup considers a 1D plasma oscillation, which is a common and simple electrostatic PIC benchmark [@Birdsall1991], [@Hockney1988],[@Jacobs2006b].
In PICLas it can be simulated either with the full Maxwell solver (DGSEM) or with the Poisson solver (HDGSEM), where the latter is
chosen for this this tutorial.
In this setup, electrons oscillate around the almost immobile ions, which creates a fluctuating electric field


Copy the `pic-poisson-plasma-wave` directory from the tutorial folder in the top level directory to a separate location

        cp -r $PICLAS_PATH/tutorials/pic-poisson-plasma-wave .
        cd pic-poisson-plasma-wave




### Mesh Generation with HOPR (pre-processing)

Before the actual simulation is conducted, a mesh file in the correct HDF5 format has to be supplied.
The mesh files used by **piclas** are created by supplying an input file *hopr.ini* with the required information for a mesh that
has either been created by an external mesh generator or directly from block-structured information in the hopr.ini file itself.
Here, a block-structured grid is created directly from the information in the hopr.ini file.
To create the .h5 mesh file, simply run

    hopr hopr.ini

This creates the mesh file *plasma_wave_mesh.h5* in HDF5 format.
The size of the simulation domain is set to [$2\pi\times0.2\times0.2$] m$^{3}$ and is defined by the single block information in the line

    Corner         =   (/0.,0.,0.,,6.2831,0.,0.,,6.2831, ...........

The number of mesh elements in each direction can be adjusted by changing the line

    nElems         = (/60,1,1/)                ! number of elements in each directio

Alternatively, if you do not want to run **hopr** yourself, you can also use the provided mesh.


Contrary to the particle boundary conditions other than periodic, the field boundaries can directly be defined in the hopr.ini file

    !=============================================================================== !
    ! BOUNDARY CONDITIONS
    !=============================================================================== !
    nUserDefinedBoundaries=6
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

    nVV=3                     ! Number of periodic displacement vectors for each periodic BC
    VV=(/6.2831 , 0.  , 0./)  ! Displacement vector 1 (x-direction)
    VV=(/0.     , 0.2 , 0./)  ! Displacement vector 2 (y-direction)
    VV=(/0.     , 0.  , 0.2/) ! Displacement vector 3 (z-direction)

and in this case a fully periodic setup is chosen by defining periodic boundaries on all six sides of the block, reflecting each
positive and negative Cartesian coordinate. In x-direction,

    BoundaryName = BC_periodicx+ ! Periodic (+vv1)
    BoundaryType = (/1,0,0,1/)   ! Periodic (+vv1)
    BoundaryName = BC_periodicx- ! Periodic (-vv1)
    BoundaryType = (/1,0,0,-1/)  ! Periodic (-vv1)

where for each of the six boundaries, a name `BoundaryNamethe` and a type `BoundaryType` must be defined (in this order)
The first "1" in `BoundaryType` corresponds to the type "periodic" and the last entry, here, either "1" or "-1" corresponds to the
first periodic vector that is defined via `VV=(/6.2831 , 0.  , 0./)` that handles periodicity in the x-direction and gives the
orientation on the boundary for the vector.
Note that each periodic boundary must have one positive and one negative corresponding boundary for the same periodic vector.

![Mesh with $60\times1\times1$ elements and a size of [$2\pi\times0.2\times0.2$] m$^{3}$.\label{fig:pic_poisson_plasma_wave}](tutorials/pic-poisson-plasma-wave/mesh/pic.pdf)





### PIC Simulation with PICLas

Install **piclas** by compiling the source code as described in Chapter \ref{chap:installation} and make sure to set the correct compile
flags

    PICLAS_EQNSYSNAME = poisson
    PICLAS_TIMEDISCMETHOD = RK3

or simply run the following command from inside the *build* directory

    cmake ../ -DPICLAS_EQNSYSNAME=poisson -DPICLAS_TIMEDISCMETHOD=RK3

to configure the build process and run `make` afterwards to build the executable.

The simulation setup is defined in *parameter.ini*. For a specific electron number density, the plasma frequency of the system is
given by

$$\omega_{p}=\omega_{e}=\sqrt{\frac{e^{2}n_{e}}{\varepsilon_{0}m_{e}}}~,$$

which is the frequency with which the charge density of the electrons oscillates, where
$\varepsilon_{0}$ is the permittivity of vacuum, $e$ is the elementary charge, $n_{e}$ and $m_{e}$
are the electron density and mass, respectively.
For the standard PIC method, the plasma frequency yields the smallest time step that has to be resolved numerically, whereas the
Debye length

$$\lambda_{D}=\sqrt{\frac{\varepsilon_{0}k_{B}T_{e}}{e^{2}n_{e}}}~,$$

where $\varepsilon_{0}$ is the permittivity of vacuum, $k_{B}$ is the Boltzmann constant, $e$ is the
elementary charge and $T_{e}$ and $n_{e}$ are the electron temperature and density, respectively, gives a spatial resolution
constraint.
In this test case, however, the electron temperature is not the defining factor for the spatial resolution because of the 1D nature
of the setup. Therefore, the resolution that is required is dictated by the gradient of the electric potential solution, i.e., the
electric field, which accelerates the charged particles, which must be adequately resolved.
The restriction on the spatial resolution is simply the number of elements (and polynomial degree $N$) that are required to resolve
the physical properties of the PIC simulation. If the temporal and spatial constraints are violated, the simulation will not yield
physical results and might even result in a termination of the simulation.

For this test case, the following physical parameters are chosen.

Table: Physical properties \label{tab:pic_poisson_plasma_wave_pyhs}

|             Property            |                Value                |
| ------------------------------- |      :------------------------:     |
| electron number density $n_{e}$ |        \SI{8e11}{\metre^{-3}}       |
|      electron mass $m_{e}$      |    \SI{9.1093826E-31}{\kilogram}    |
|    ion number density $n_{i}$   |        \SI{8e11}{\metre^{-3}}       |
|         ion mass $m_{i}$        |   \SI{1.672621637E-27}{\kilogram}   |
|  electron/ion charge $q_{i,e}$  | $\pm$\SI{1.60217653E-19}{\coulomb}  |


#### Numerical Setup

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
step iterations `IterDisplayStep` between information output regarding the current status of the simulation that is written to STD.out.

The remaining numerical parameters are selected for the field and particle solver and run-time analysis.
As there are no walls present in the setup, all boundaries are set as periodic boundary conditions for the field as well as the
particle solver. The particle boundary conditions are set by the following lines

    ! =============================================================================== !
    ! PARTICLES Boundary Conditions
    ! =============================================================================== !
    Part-maxParticleNumber    = 4000
    Part-nSpecies             = 2

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

where, the maximum number of particles for each processor `Part-maxParticleNumber` has to be supplied, the number of particle
species `Part-nSpecies` that are present and number number of boundaries `Part-nBounds` (a in 3D cuboid), followed by the names of
the boundaries (given by the hopr.ini file) and the type `periodic`. Furthermore, the periodic vectors must be supplied and the size
of the Cartesian background mesh `Part-FIBGMdeltas`, which can be accompanied by a division factor in each direction given by
`Part-FactorFIBGM`.

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

    PIC-shapefunction-dimension = 1
    PIC-shapefunction-direction = 1
    PIC-shapefunction-alpha     = 4

    PIC-Deposition-Type = shape_function_adaptive

where the interpolation type `PIC-Interpolation-Type` is required. For deposition, a polynomial shape function with exponent
`PIC-shapefunction-alpha` of type
`PIC-Deposition-Type` is used. The dimension `PIC-shapefunction-dimension`, here 1D and direction `PIC-shapefunction-direction`, are
selected specifically for the one-dimensional setup that is simulated here.

The inserting (sometimes labelled emission or initialization) of particles at the beginning or during the course of the simulation
is controlled via

    ! =============================================================================== !
    ! PARTICLE Emission
    ! =============================================================================== !

    Part-Species1-MacroParticleFactor = 5e8
    Part-Species2-MacroParticleFactor = 5e8
    Part-Species1-Init1-ParticleNumber      = 400
    Part-Species1-Init1-maxParticleNumber-x = 400
    Part-Species2-Init1-ParticleNumber      = 400
    Part-Species2-Init1-maxParticleNumber-x = 400

    ! -------------------------------------
    ! Electrons 1
    ! -------------------------------------
    Part-Species1-ChargeIC            = -1.60217653E-19
    Part-Species1-MassIC              = 9.1093826E-31

    Part-Species1-nInits=1

    Part-Species1-Init1-SpaceIC               = sin_deviation
    Part-Species1-Init1-velocityDistribution  = constant
    Part-Species1-Init1-maxParticleNumber-y   = 1
    Part-Species1-Init1-maxParticleNumber-z   = 1
    Part-Species1-Init1-Amplitude             = 0.01
    Part-Species1-Init1-WaveNumber            = 2.
    Part-Species1-Init1-VeloIC                = 0.
    Part-Species1-Init1-VeloVecIC             = (/1.,0.,0./)

    ! -------------------------------------
    ! Ions 2
    ! -------------------------------------
    Part-Species2-ChargeIC            = 1.60217653E-19
    Part-Species2-MassIC              = 1.672621637E-27

    Part-Species2-nInits=1

    Part-Species2-Init1-SpaceIC               = sin_deviation
    Part-Species2-Init1-velocityDistribution  = constant
    Part-Species2-Init1-maxParticleNumber-y   = 1
    Part-Species2-Init1-maxParticleNumber-z   = 1
    Part-Species2-Init1-Amplitude             = 0.0
    Part-Species2-Init1-WaveNumber            = 0.
    Part-Species2-Init1-VeloIC                = 0.0
    Part-Species2-Init1-VeloVecIC             = (/0.,0.,0./)
    ! -------------------------------------

where, for each particle species, the number `Part-SpeciesX-nInits` controls how many initialization blocks are to be used.
Each block is accompanied by a set of parameters that start from `Part-SpeciesX-Init1-SpaceIC` up to `Part-SpeciesX-Init1-VeloVecIC`
, which, in the above example, distribute the particles equidistantly on a line and sinusoidally dislocates them (representing
an initial stage of a plasma wave in 1D).

The extent of dislocation is controlled by `Part-SpeciesX-Init1-Amplitude`, which is only set for the electron species as the ion
species is no dislocated (they remain equidistantly distributed).
The parameter `Part-SpeciesX-Init1-WaveNumber` set the number of sin wave repetitions in the `x`-direction of the domain.

The number of simulation particles, given by`Part-SpeciesX-Init1-ParticleNumber` and weighted by
`Part-SpeciesX-MacroParticleFactor`, the multiplication of which gives the number of real physical particles and together with the volume
of the complete domain, the density of each species.

In case of the `SpaceIC=sin\_deviation`, the number of simulation particles must be equal to the multiplied values given in
`Part-SpeciesX-Init1-maxParticleNumber-x/y/z` as this emission type allows distributing the particles not only in one, but in all
three Cartesian coordinates, which is not required for this 1D example.

Furthermore, the masses `Part-SpeciesX-MassIC`and charges `Part-SpeciesX-ChargeIC` of each species are required.

Finally, some parameters for run-time analysis are chosen by setting them `T` (true).

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

where the function of the parameters is given in the code comments. Information regarding every parameter can be obtained from
running the command

    ./piclas --help "CalcCharge"

where each parameter can simply be supplied to the *help* module of **piclas**. This help module can also output the complete set of
parameters via `./piclas --help` or a subset of them by supplying a section, e.g., `./piclas --help "HDG"` for the HDGSEM solver.

#### Running the code

The command

~~~~~~~
piclas parameter.ini > std.out
~~~~~~~

runs the code and dumps all output into the file *std.out*.
If the run has completed successfully, which should take only a brief moment, the contents of the working folder should look like

    4.0K drwxrwxr-x  4.0K Jun 28 13:07 ./
    4.0K drwxrwxr-x  4.0K Jun 25 23:56 ../
    8.0K -rw-rw-r--  5.8K Jun 28 12:51 ElemTimeStatistics.csv
    120K -rw-rw-r--  113K Jun 28 12:51 FieldAnalyze.csv
    4.0K -rw-rw-r--  2.1K Jun 26 16:49 hopr.ini
    8.0K -rw-rw-r--  5.0K Jun 28 13:07 parameter.ini
    156K -rw-rw-r--  151K Jun 28 12:51 PartAnalyze.csv
     32K -rw-rw-r--   32K Jun 26 16:43 plasma_wave_mesh.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:44 plasma_wave_State_000.00000000000000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:45 plasma_wave_State_000.00000000400000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:45 plasma_wave_State_000.00000000800000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:46 plasma_wave_State_000.00000001200000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:47 plasma_wave_State_000.00000001600000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:48 plasma_wave_State_000.00000002000000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:49 plasma_wave_State_000.00000002400000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:50 plasma_wave_State_000.00000002800000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:50 plasma_wave_State_000.00000003200000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 plasma_wave_State_000.00000003600000000.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 plasma_wave_State_000.00000004000000000.h5
     72K -rw-rw-r--   71K Jun 28 12:51 std.out

Multiple additional files have been created, which are are named  **Projectname_State_Timestamp.h5**.
They contain the solution vector of the equation system variables at each interpolation node at the given time, which corresponds
to multiplies of **Analyze_dt**. If these files are not present, something went wrong during the execution of **piclas**.
In that case, check the _std.out_ file for an error message.

After a successful completion, the last lines in this files should look like in figure \ref{fig:freestream_stdout}

<!--![The _std.out_ file after a successful run\label{fig:freestream_stdout}](tutorials/00_freestream/freestream_stdout.png)-->



### Visualization (post-processing)

To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**, **VisIt** or any other
visualisation tool for which the program **piclas2vtk** is used.

The parameters for **piclas2vtk** are stored in the **parameter.ini** file under

    ! =============================================================================== !
    ! piclas2vtk
    ! =============================================================================== !
    NVisu         = 10
    VisuParticles = T

where `NVisu` is the polynomial visualization degree on which the field solution is interpolated. The flag `VisuParticles` activates
the output of particle position, velocity and species to the *vtk*-files.

Run the command

~~~~~~~
piclas2vtk parameter.ini  plasma_wave_State_000.000000*
~~~~~~~
to generate the corresponding *vtk*-files, which can then be loaded into the visualisation tool.

The electric potential field can be viewed, e.g., by opening `plasma_wave_Solution_000.00000040000000000.vtu` and plotting the field
`Phi`, which should look like the following


![Electric potential and field strength for the plasma wave test case.](tutorials/pic-poisson-plasma-wave/results/pic.pdf)
















