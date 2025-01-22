# Negative Streamer Discharge (PIC, Poisson's Equation)

## Drift-Diffusion Electron Fluid

Streamers, as seen in lightning and sprite events, are also observed in industrial applications such as lighting, purification of dirty gases and water, plasma jets and bullets for disinfection, and plasma-assisted combustion. Using this setup, a 1D streamer simulation will be carried out using PICLas.

Before beginning with the tutorial, copy the `pic-poisson-streamer` from the tutorial folder in the top level directory to a separate location.

    cp -r $PICLAS_PATH/tutorials/pic-poisson-streamer .
    cd pic-poisson-streamer

If the piclas repository is located in the home directory,

    cp -r /home/$(whoami)/piclas/tutorials/pic-poisson-streamer .
    cd pic-poisson-streamer

For the first part of the tutorial, switch into the `drift-diffusion` directory via

    cd drift-diffusion

before continuing with the mesh building section. Then, the species database needs to be copied into that directory of a symbolic link needs to be created

    ln -s $PICLAS_PATH/SpeciesDatabase.h5

### Mesh Generation with HOPR (pre-processing)-1

Before the actual simulation is conducted, a mesh file in the correct HDF5 format has to be supplied.
The mesh files used by piclas are created by supplying an input file hopr.ini with the required
information for a mesh that has either been created by an external mesh generator or directly from
block-structured information in the hopr.ini file itself.
Here, a block-structured grid is created directly from the information in the hopr.ini file.

To create the .h5 mesh file, simply run

    hopr hopr.ini

This creates the mesh file *streamer_small_1500_mesh.h5* in HDF5 format.
Alternatively, if you do not want to run **hopr** yourself, you can also use the provided mesh. The only difference from the mesh created in previous tutorial, is the size of the simulation domain in this tutorial which is set to [$4\pi\times1\times1$] m$^{3}$ and is defined by the single block information

in the line, where each node of the hexahedral element is defined. Here, l1 and l2 define the dimensions of the hexahedral elements. While l1 represents the length in the x-direction, l2 indicates the length in the y-direction.

    Corner     = (/0.,0.,0.,,l1,0.,0.,,l1,l2,0.,,0.,l2,0.,,0.,0.,l2,,l1,0.,l2,,l1,l2,l2,,0.,l2,l2/)

The number of mesh elements for the block in each direction can be adjusted by changing the line. ni defines the number of elements to be created in the x-direction. This allows for adjusting the number of elements within a block.

    nElems     = (/ni,1,1/)                  ! number of elements in each direction

Each side of the block has to be assigned a boundary index, which corresponds to the boundaries defined in the next steps

    BCIndex    = (/5,3,2,4,1,6/)             ! Indices of UserDefinedBoundaries


The field boundaries can directly be defined in the hopr.ini file (contrary to the particle boundary conditions, which are defined in the parameter.ini).
Periodic boundaries should be defined in the hopr.ini.

    !=============================================================================== !
    ! BOUNDARY CONDITIONS
    !=============================================================================== !
    BoundaryName = BC_Xleft
    BoundaryType = (/3,0,0,0/)
    BoundaryName = BC_Xright
    BoundaryType = (/3,0,0,0/)
    BoundaryName = BC_periodicy+ ! Periodic (+vv1)
    BoundaryType = (/1,0,0,1/)   ! Periodic (+vv1)
    BoundaryName = BC_periodicy- ! Periodic (-vv1)
    BoundaryType = (/1,0,0,-1/)  ! Periodic (-vv1)
    BoundaryName = BC_periodicz+ ! Periodic (+vv2)
    BoundaryType = (/1,0,0,2/)   ! Periodic (+vv2)
    BoundaryName = BC_periodicz- ! Periodic (-vv2)
    BoundaryType = (/1,0,0,-2/)  ! Periodic (-vv2)

    VV=(/0. , l2 , 0./)    ! Displacement vector 1 (y-direction)
    VV=(/0. , 0. , l2/)   ! Displacement vector 2 (z-direction)


### Simulation with PICLas-1

Install piclas by compiling the source code as described in Chapter Installation, specifically described under Section Compiling the code. Always build the code in a separate directory located in the piclas top level directory. For this PIC tutorial, e.g., create a directory build_drift_diffusion by running

    cd $PICLAS_PATH

where the variable PICLAS_PATH contains the path to the location of the piclas repository. If the piclas repository is located in the home directory,

    cd /home/$(whoami)/piclas

and the create the build directory, in which the compilation process will take place. This command should be run in the piclas folder.

    mkdir build_drift_diffusion

Always compile the code within the build directory, hence, navigate to the build_drift_diffusion directory before running cmake

    cd build_drift_diffusion

For this specific tutorial, make sure to set the correct compile flags

    PICLAS_EQNSYSNAME       = drift_diffusion
    LIBS_USE_PETSC          = ON
    PICLAS_READIN_CONSTANTS = ON
    PICLAS_TIMEDISCMETHOD   = Explicit-FV


using the ccmake (gui for cmake) or simply run the following command from inside the build directory

    cmake .. -DPICLAS_EQNSYSNAME=drift_diffusion -DPICLAS_TIMEDISCMETHOD=Explicit-FV -DLIBS_USE_PETSC=ON

to configure the build process and run

    make -j


afterwards to compile the executable. For this setup, we have chosen the drift_diffusion solver and selected the Explicit-FV time discretization method. An overview over the available solver and discretization options is given in Section Solver settings. To run the simulation and analyse the results, the piclas and piclas2vtk executables have to be run. To avoid having to use the entire file path, you can either set aliases for both, copy them to your local tutorial directory or create a link to the files via

    ln -s $PICLAS_PATH/build_drift_diffusion/bin/piclas
    ln -s $PICLAS_PATH/build_drift_diffusion/bin/piclas2vtk

where the variable PICLAS_PATH contains the path to the location of the piclas repository. If the piclas repository is located in the home directory, the two commands

    ln -s /home/$(whoami)/piclas/build_drift_diffusion/bin/piclas
    ln -s /home/$(whoami)/piclas/build_drift_diffusion/bin/piclas2vtk


#### General numerical setup-1

The general numerical parameters (defined in the parameter.ini) are selected by the following

    ! =============================================================================== !
    ! DISCRETIZATION
    ! =============================================================================== !
    N             = 1  ! Polynomial degree of the DG method (field solver)

    !IniRefState   = 1

    IniExactFunc-FV  = 3  !streamer
    IniRefState-FV   = 1
    ! =============================================================================== !
    ! MESH
    ! =============================================================================== !
    MeshFile      = ./streamer_small_1500_mesh.h5 ! Relative path to the mesh .h5 file

    ! =============================================================================== !
    ! General
    ! =============================================================================== !
    ProjectName       = streamer_N2
    doPrintStatusLine = T                   ! Output live of ETA




where, among others, the polynomial degree $N$, the path to the mesh file `MeshFile`, project name and the option to print the ETA
to the terminal output in each time step.

The temporal parameters of the simulation are controlled via

    ! =============================================================================== !
    ! CALCULATION
    ! =============================================================================== !
    ManualTimeStep  = 1.e-13
    tend            = 7.e-10
    Analyze_dt      = 1.e-10
    IterDisplayStep = 100

where the time step for the field and particle solver is set via `ManualTimeStep`, the final simulation time `tend`, the time
between restart/checkpoint file output `Analyze_dt` (also the output time for specific analysis functions) and the number of time
step iterations `IterDisplayStep` between information output regarding the current status of the simulation that is written to std.out.
The remaining parameters are selected for the field and particle solver as well as run-time analysis.


#### Boundary conditions-1

As there are no walls present in the setup, all boundaries are set as periodic boundary conditions for the field as well as the
particle solver. The particle boundary conditions are set by the following lines

    ! =============================================================================== !
    ! PARTICLE Boundary Conditions
    ! =============================================================================== !
    Part-nBounds              = 6             ! Number of particle boundaries
    Part-Boundary1-SourceName = BC_Xleft      ! Name of 1st particle BC
    Part-Boundary1-Condition  = symmetric     ! Type of 1st particle BC
    Part-Boundary2-SourceName = BC_Xright     ! ...
    Part-Boundary2-Condition  = symmetric     ! ...
    Part-Boundary3-SourceName = BC_periodicy+ ! ...
    Part-Boundary3-Condition  = periodic      ! ...
    Part-Boundary4-SourceName = BC_periodicy- ! ...
    Part-Boundary4-Condition  = periodic      ! ...
    Part-Boundary5-SourceName = BC_periodicz+ ! ...
    Part-Boundary5-Condition  = periodic      ! ...
    Part-Boundary6-SourceName = BC_periodicz- ! ...
    Part-Boundary6-Condition  = periodic      ! ...

    Part-nPeriodicVectors = 2 ! Number of periodic boundary (particle and field) vector


#### Field solver-1

The settings for the field solver (HDGSEM) are given by

    ! =============================================================================== !
    ! Field Solver: HDGSEM
    ! =============================================================================== !
    epsCG                 = 1e-6   ! Stopping criterion (residual) of iterative CG solver (default that is used for the HDGSEM solver)
    maxIterCG             = 10000  ! Maximum number of iterations
    IniExactFunc          = 0      ! Initial field condition. 0: zero solution vector

where `epsCG` sets the abort residual of the CG solver, `maxIterCG` sets the maximum number of iterations within the CG solver and
`IniExactFunc` set the initial solution of the field solver (here 0 says that nothing is selected).

The numerical scheme for tracking the movement of all particles throughout the simulation domain can be switched by

    ! =============================================================================== !
    ! Particle Solver
    ! =============================================================================== !
    TrackingMethod    = triatracking ! Particle tracking method

The PIC parameters for interpolation (of electric fields to the particle positions) and deposition (mapping of charge properties
from particle locations to the grid) are selected via

    ! =============================================================================== !
    ! PIC: Interpolation/Deposition
    ! =============================================================================== !
    PIC-DoInterpolation            = T                        ! Activate Lorentz forces acting on charged particles
    PIC-Interpolation-Type         = particle_position        ! Field interpolation method for Lorentz force calculation
    PIC-Deposition-Type            = cell_volweight_mean
    PIC-shapefunction-radius       = 0.0001
    PIC-shapefunction-dimension    = 1                        ! Shape function specific dimensional setting


where the interpolation type `PIC-Interpolation-Type = particle_position` is currently the only option for specifying how
electro(-magnetic) fields are interpolated to the position of the charged particles.
The dimension `PIC-shapefunction-dimension`, here 1D and direction `PIC-shapefunction-direction`, are selected specifically
for the one-dimensional setup that is simulated here. The different available deposition types are described in more detail in Section {ref}`sec:PIC-deposition`.


#### Particle solver-1

For the treatment of particles, the maximum number of particles `Part-maxParticleNumber` that each processor can hold has to be supplied and
the number of particle species `Part-nSpecies` that are used in the simulation (created initially or during the simulation time
through chemical reactions).

    ! =============================================================================== !
    ! PARTICLE Emission
    ! =============================================================================== !
    Particles-Species-Database = SpeciesDatabase.h5
    Part-nSpecies             = 2        ! Number of particle species


    ! -------------------------------------
    ! Background Gas - N2
    ! -------------------------------------
    Part-Species1-SpeciesName = N2
    Part-Species1-MacroParticleFactor       = 1600.
    BGGas-DriftDiff-Database = Phelps


    Part-Species1-nInits                      = 1

    Part-Species1-Init1-velocityDistribution  = maxwell_lpn
    Part-Species1-Init1-SpaceIC               = background
    Part-Species1-Init1-PartDensity           = 2.506e25
    Part-Species1-Init1-VeloIC                = 0.
    Part-Species1-Init1-VeloVecIC             = (/1.,0.,0./)
    Part-Species1-Init1-MWTemperatureIC       = 298.
    Part-Species1-Init1-TempVib               = 298.
    Part-Species1-Init1-TempRot               = 298.
    Part-Species1-Init1-TempElec              = 298.

    ! -------------------------------------
    ! Ions - N2+
    ! -------------------------------------
    Part-Species2-SpeciesName = N2Ion1
    Part-Species2-MacroParticleFactor = 1600. ! Weighting factor for species #2
    Part-Species2-nInits              = 0               ! Number of initialization/emission regions for species #1

    Particles-HaloEpsVelo = 3.E8


#### Field Boundaries-1

BoundaryName specifies the name of the boundary. Two boundaries are defined here: BC_Xleft (left x boundary) and BC_Xright (right x boundary). BoundaryType defines the type of boundary condition. Each BoundaryType is expressed as a pair of numbers, such as ((/12,0/), (/4,0/)). '12' represents a Neumann boundary condition, which defines a situation where the electric field (E) remains constant. '4' represents a Dirichlet boundary condition, indicating a situation where the potential is zero (for example, a grounded surface). '11' defines a special type of Neumann boundary condition. In this case, it represents conditions where the electric field (E) has a constant value, but unlike number 12, '11' includes a more specific condition.

    ! =============================================================================== !
    ! Field Boundaries
    ! =============================================================================== !
    BoundaryName = BC_Xleft
    BoundaryType = (/12,0/)                ! 12: Neumann
    BoundaryType-FV = (/3,0/) ! Neumann for eFluid


    BoundaryName = BC_Xright
    BoundaryType = (/4,0/) ! 4: Dirichlet with zero potential
    BoundaryType-FV = (/3,0/) ! Neumann for eFluid


#### Analysis setup-1

Finally, some parameters for run-time analysis are chosen by setting them `T` (true). Further, with `TimeStampLength = 18`, the names of the output files are shortened for better postprocessing. If this is not done, e.g. Paraview does not sort the files correctly and will display faulty behaviour over time.

    ! =============================================================================== !
    ! Analysis
    ! =============================================================================== !
    TimeStampLength          = 18 ! Reduces the length of the timestamps in filenames for better postprocessing
    CalcCharge               = T  ! writes rel/abs charge error to PartAnalyze.csv
    CalcPotentialEnergy      = T  ! writes the potential field energy to FieldAnalyze.csv
    CalcKineticEnergy        = T  ! writes the kinetic energy of all particle species to PartAnalyze.csv
    PIC-OutputSource         = T  ! writes the deposited charge (RHS of Poissons equation to XXX_State_000.0000XXX.h5)
    CalcPICTimeStep          = T  ! writes the PIC time step restriction to XXX_State_000.0000XXX.h5 (rule of thumb)
    CalcPointsPerDebyeLength = T  ! writes the PIC grid step restriction to XXX_State_000.0000XXX.h5 (rule of thumb)
    CalcTotalEnergy          = T  ! writes the total energy of the system to PartAnalyze.csv (field and particle)

    CalcElectronTemperature = T
    CalcDebyeLength = T
    CalcPlasmaFrequency = T




The function of each parameter is given in the code comments. Information regarding every parameter can be obtained from running the command

    piclas --help "CalcCharge"

where each parameter is simply supplied to the *help* module of **piclas**. This help module can also output the complete set of
parameters via `piclas --help` or a subset of them by supplying a section, e.g., `piclas --help "HDG"` for the HDGSEM solver.


#### Running the code-1

The command

    ./piclas parameter.ini | tee std.out

executes the code and dumps all output into the file *std.out*.

If the run has completed successfully, which should take only a brief moment, the contents of the working folder should look like



    8.0K -rw-rw-r--  5.8K *** 28 12:51 ElemTimeStatistics.csv
    120K -rw-rw-r--  113K *** 28 12:51 FieldAnalyze.csv
    4.0K -rw-rw-r--  2.1K *** 26 16:49 hopr.ini
    8.0K -rw-rw-r--  5.0K *** 28 13:07 parameter_streamer.ini
    156K -rw-rw-r--  151K *** 28 12:51 PartAnalyze.csv
     32K -rw-rw-r--   32K *** 26 16:43 streamer_N2_mesh.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:44 streamer_N2_DSMCState_000.00000000001000.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:45 streamer_N2_DSMCState_000.00000000002000.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:45 streamer_N2_DSMCState_000.00000000003000.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:46 streamer_N2_DSMCState_000.00000000004000.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:47 streamer_N2_DSMCState_000.00000000005000.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:48 streamer_N2_DSMCState_000.0000000000****.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:49 streamer_N2_DSMCState_000.0000000000****.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:50 streamer_N2_DSMCState_000.0000000000****.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:50 streamer_N2_DSMCState_000.0000000000****.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:51 streamer_N2_State_000.00000000099000.h5
    1.6M -rw-rw-r--  1.6M *** 28 12:51 streamer_N2_State_000.00000000100000.h5
     72K -rw-rw-r--   71K *** 28 12:51 std.out

Multiple additional files have been created, which are are named  **Projectname_State_Timestamp.h5**.
They contain the solution vector of the equation system variables at each interpolation nodes at the given time, which corresponds
to multiples of **Analyze_dt**. If these files are not present, something went wrong during the execution of **piclas**.
In that case, check the `std.out` file for an error message.

After a successful completion, the last lines in this file should look as shown below:

    Timestep  :    1.0000000E-13
    #Timesteps :    7.0000000E+03
    WRITE STATE TO HDF5 FILE [streamer_N2_State_000.00000000070000.h5] ... DONE [ 0.09 sec ] [ 0:00:00:00 ]
    #Particles :    1.5270500E+05    Average particles per proc :    1.9088125E+04    Min :    0.0000000E+00    Max :    3.6427000E+04
    ------------------------------------------------------------------------------------------------------------------------------------
    #Particles :    1.5270500E+05 (peak)         Average (peak) :    1.9088125E+04    Min :    0.0000000E+00    Max :    3.6432000E+04
    ====================================================================================================================================
    PICLAS FINISHED! [ 1003.58 sec ] [ 0:00:16:43 ]
    ====================================================================================================================================



### Visualization (post-processing)-1

To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**, **VisIt** or any other
visualisation tool for which the program **piclas2vtk** is used.

The parameters for **piclas2vtk** are stored in the **parameter.ini** file under

    ! =============================================================================== !
    ! piclas2vtk
    ! =============================================================================== !
    NVisu         = 1 ! Polynomial degree used for the visualization when the .h5 file is converted to .vtu/.vtk format. Should be at least N+1
    VisuParticles = F  ! Activate the conversion of particles from .h5 to .vtu/.vtk format. Particles will be displayed as a point cloud with properties, such as velocity, species ID, etc.

    Part-NumberOfRandomSeeds = 2
    Particles-RandomSeed1    = 1
    Particles-RandomSeed2    = 2


where `NVisu` is the polynomial visualization degree on which the field solution is interpolated.
Depending on the used polynomial degree `N`, the degree of visualization `NVisu` should always be higher than
`N` because the PIC simulation is always subject to noise that is influenced by the discretization (number of elements and
polynomial degree as well as number of particles) and is visible in the solution results of the current simulation.

Additionally, the flag `VisuParticles` activates the output of particle position, velocity and species to the *vtk*-files.

Run the command

    ./piclas2vtk parameter.ini streamer_N2_State_000.00000000*

to generate the corresponding *vtk*-files, which can then be loaded into the visualisation tool.

The electron and ion densities (streamer_N2_Solution_ElemData_000.00000000070000.vtu) and electric field in x-direction (streamer_N2_Solution_000.00000000070000.vtu) at t=0.7 ns (end of simulation) are shown in the plots.


```{figure} results/fluid_density.png
---
name: fig:streamer-fluid
width: 700px
---

Electron and ion density in a planar front in nitrogen gas.
```

```{figure} results/fluid_electric.png
---
name: fig:streamer-elec
width: 700px
---


Electric Field across the domain.
```


## Particle-in-cell/Monte Carlo Collision (PIC/MCC)

In this tutorial, we build upon the foundational concepts introduced in the fluid model tutorial for negative streamer simulations. While the fluid model provides a macroscopic perspective, the MCC model simulates the microscopic particle dynamics through direct Monte Carlo collision techniques. This approach allows us to capture electron collision and ionization events with high fidelity, particularly critical in simulating the formation and propagation of negative streamers in nitrogen gas. Here, we will guide you through setting up the MCC model in piclas for a detailed and accurate streamer analysis.

This tutorial demonstrates how to simulate a negative streamer discharge using the Particle-in-Cell (PIC) with Monte Carlo Collision (MCC) model in PICLas.

A streamer discharge is a form of plasma created when an electric field in a gas accelerates free electrons, leading to further ionization. We will focus on simulating the evolution of a negative streamer propagating in a 1D domain.

Negative streamers are associated with the rapid propagation of electrons under the influence of strong electric fields in gases like nitrogen or air. When these electrons ionize the gas, they create an ionization front which can grow rapidly if the electric field is strong enough.

Before beginning the tutorial, copy the `pic-poisson-streamer` directory from the tutorial folder in the top level
directory to a separate location

    cp -r $PICLAS_PATH/tutorials/pic-poisson-streamer .
    cd pic-poisson-streamer

If the piclas repository is located in the home directory,

    cp -r /home/$(whoami)/piclas/tutorials/pic-poisson-streamer .
    cd pic-poisson-streamer

For the second part of the tutorial, switch into the `pic-mcc` directory via

    cd pic-mcc

before continuing with the mesh building section. Then, the species database needs to be copied into that directory of a symbolic link needs to be created

    ln -s $PICLAS_PATH/SpeciesDatabase.h5

In addition, an electronic database is needed. This can be linked or copied from one of the regressionchecks, e.g.

    cp $PICLAS_PATH/regressioncheck/NIG_Reservoir/MCC_BGG_MultiSpec_XSec_TCE_QK_Chem/DSMCSpecies_electronic_state_full_Data.h5 .


### Mesh Generation with HOPR (pre-processing)-2

Before the actual simulation is conducted, a mesh file in the correct HDF5 format has to be supplied.
The mesh files used by **piclas** are created by supplying an input file *hopr.ini* with the required information for a mesh that
has either been created by an external mesh generator or directly from block-structured information in the hopr.ini file itself.
Here, a block-structured grid is created directly from the information in the hopr.ini file.
To create the .h5 mesh file, simply run

    hopr hopr.ini

This creates the mesh file *pic_mc_streamer_mesh.h5* in HDF5 format.
Alternatively, if you do not want to run **hopr** yourself, you can also use the provided mesh.

The size of the simulation domain is set to [$0.012\times0.00005\times0.00005$] m$^{3}$ and is defined by the single block information
in the line, where each node of the hexahedral element is defined

    Corner         =   (/0.,0.,0.,,0.012,0.,0.,,0.012,0.00005,0.,,0.,0.00005,0.,,0.,0.,0.00005,,0.012,0.,0.00005,,0.012,0.00005,0.00005,,0.,0.00005,0.00005/)

The number of mesh elements for the block in each direction can be adjusted by changing the line

    nElems         = (/4800,1,1/)                ! number of elements in each direction (x,y,z)

Each side of the block has to be assigned a boundary index, which corresponds to the boundaries defined in the next steps

    BCIndex        = (/5,3,2,4,1,6/)

The field boundaries can directly be defined in the hopr.ini file (contrary to the particle boundary conditions, which are defined
in the parameter.ini).
Periodic boundaries always have to be defined in the hopr.ini.

    !=============================================================================== !
    ! BOUNDARY CONDITIONS
    !=============================================================================== !
    BoundaryName = BC_Xleft
    BoundaryType = (/3,0,0,0/)
    BoundaryName = BC_Xright
    BoundaryType = (/3,0,0,0/)
    BoundaryName = BC_periodicy+ ! Periodic (+vv1)
    BoundaryType = (/1,0,0,1/)   ! Periodic (+vv1)
    BoundaryName = BC_periodicy- ! Periodic (-vv1)
    BoundaryType = (/1,0,0,-1/)  ! Periodic (-vv1)
    BoundaryName = BC_periodicz+ ! Periodic (+vv2)
    BoundaryType = (/1,0,0,2/)   ! Periodic (+vv2)
    BoundaryName = BC_periodicz- ! Periodic (-vv2)
    BoundaryType = (/1,0,0,-2/)  ! Periodic (-vv2)

    VV = (/0.     , 0.00005 , 0./)            ! Displacement vector 1 (y-direction)
    VV = (/0.     , 0.  , 0.00005/)           ! Displacement vector 2 (z-direction)


### Simulation with PICLas-2

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

    cmake ../ -DPICLAS_EQNSYSNAME=poisson -DPICLAS_TIMEDISCMETHOD=Leapfrog -DLIBS_USE_PETSC=ON

to configure the build process and run

    make -j

afterwards to compile the executable. For this setup, we have chosen the Poisson solver
and selected the Leapfrog time discretization method. An overview over the available solver
and discretization options is given in Section {ref}`sec:solver-settings`.
To run the simulation and analyse the results, the *piclas* and *piclas2vtk* executables have to be run.
To avoid having to use the entire file path, you can either set aliases for both, copy them to your local tutorial directory or
create a link to the files via

    ln -s $PICLAS_PATH/build_poisson_Leapfrog/bin/piclas
    ln -s $PICLAS_PATH/build_poisson_Leapfrog/bin/piclas2vtk

where the variable `PICLAS_PATH` contains the path to the location of the piclas repository.
If the piclas repository is located in the home directory, the two commands

    ln -s /home/$(whoami)/piclas/build_poisson_Leapfrog/bin/piclas
    ln -s /home/$(whoami)/piclas/build_poisson_Leapfrog/bin/piclas2vtk

can be executed.

Please check where piclas is located before running the commands.

The simulation setup is defined in *parameter.ini*.

The dominant reaction in a negative streamer discharge is the **electron impact ionization** of neutral gas molecules. When a free electron, accelerated by the electric field, collides with a neutral molecule (such as nitrogen in our case), it ionizes the molecule by knocking out an electron. This process generates a positive ion and an additional free electron, contributing to the electron avalanche that sustains and propagates the streamer.

For nitrogen $N_2$, the ionization reaction is:

$$
e^- + N_2 \rightarrow 2e^- + N_2^+
$$


In this reaction:
- Electrons collides with a nitrogen molecule.
- The collision results in the creation of a positive ion $N_2^{+}$ and two free electrons.
- The newly created free electrons continue to ionize other molecules, which leads to the avalanche effect that is characteristic of streamer propagation.

This ionization process is the key mechanism driving the rapid growth of the streamer and the formation of plasma within the discharge.




<!-- For a specific electron number density, the plasma frequency of the system is

$$\omega_{p}=\nu_{e}=\sqrt{\frac{e^{2}n_{e}}{\varepsilon_{0}m_{e}}}~,$$

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

 -->


The physical parameters for this test case are presnt in `DSMC.ini`. These read-in constants are also present in the databse called the `SpeciesDatabase.h5`provided at the top level of the PICLas directory. In our simulation setup we will simply be using the `DSMC.ini` which contains physical properties of all the species in the setup.


#### General numerical setup-2

The general numerical parameters (defined in the parameter.ini) are selected by the following

    ! =============================================================================== !
    ! DISCRETIZATION
    ! =============================================================================== !
    N             = 2  ! Polynomial degree of the DG method (field solver)

    ! =============================================================================== !
    ! MESH
    ! =============================================================================== !
    MeshFile      = streamer_mesh.h5 ! Relative path to the mesh .h5 file

    ! =============================================================================== !
    ! General
    ! =============================================================================== !
    ProjectName       = streamer_N2 ! Project name that is used for naming state files
    ColoredOutput     = F           ! Turn ANSI terminal colors ON/OFF
    doPrintStatusLine = T           ! Output live of ETA
    TrackingMethod    = TriaTracking

where, among others, the polynomial degree $N$, the path to the mesh file `MeshFile`, project name and the option to print the ETA
to the terminal output in each time step.

The temporal parameters of the simulation are controlled via

    ! =============================================================================== !
    ! CALCULATION
    ! =============================================================================== !
    ManualTimeStep                   = 1e-13  ! Fixed pre-defined time step only when using the Poisson solver. Maxwell solver calculates dt that considers the CFL criterion
    tend                             = 7e-10  ! Final simulation time
    Analyze_dt                       = 5E-12  ! decrease analyze_dt for better resolution of fish eye
    IterDisplayStep                  = 500    ! Number of iterations between terminal output showing the current time step iteration
    Part-WriteMacroValues            = T      ! Set [T] to activate ITERATION DEPENDANT h5 output of macroscopic values sampled (Can not be enabled together with Part-TimeFracForSampling)
    Part-IterationForMacroVal        = 100    ! Set number of iterations used for sampling                                      (Can not be enabled together with Part-TimeFracForSampling)

where the time step for the field and particle solver is set via `ManualTimeStep`, the final simulation time `tend`, the time
between restart/checkpoint file output `Analyze_dt` (also the output time for specific analysis functions) and the number of time
step iterations `IterDisplayStep` between information output regarding the current status of the simulation that is written to std.out.
The remaining parameters are selected for the field and particle solver as well as run-time analysis.


#### Boundary conditions-2

As there are no walls present in the setup, all boundaries are set as periodic boundary conditions for the field as well as the
particle solver. The particle boundary conditions are set by the following lines

    ! =============================================================================== !
    ! PARTICLE Boundary Conditions
    ! =============================================================================== !
    Part-nBounds              = 6             ! Number of particle boundaries
    Part-Boundary1-SourceName = BC_Xleft      ! Name of 1st particle BC
    Part-Boundary1-Condition  = symmetric     ! Type of 1st particle BC
    Part-Boundary2-SourceName = BC_Xright     ! ...
    Part-Boundary2-Condition  = symmetric     ! ...
    Part-Boundary3-SourceName = BC_periodicy+ ! ...
    Part-Boundary3-Condition  = periodic      ! ...
    Part-Boundary4-SourceName = BC_periodicy- ! ...
    Part-Boundary4-Condition  = periodic      ! ...
    Part-Boundary5-SourceName = BC_periodicz+ ! ...
    Part-Boundary5-Condition  = periodic      ! ...
    Part-Boundary6-SourceName = BC_periodicz- ! ...
    Part-Boundary6-Condition  = periodic      ! ...

    Part-nPeriodicVectors = 2 ! Number of periodic boundary (particle and field) vectors

    Part-FIBGMdeltas = (/0.0012 , 0.0001, 0.0001/)  ! Cartesian background mesh (bounding box around the complete simulation domain)
    Part-FactorFIBGM = (/1500     , 1   , 1/)       ! Division factor that is applied t the "Part-FIBGMdeltas" values to define the dx, dy and dz distances of the Cartesian background mesh


#### Field solver-2

The settings for the field solver (HDGSEM) are given by

    ! =============================================================================== !
    ! Field Solver: HDGSEM
    ! =============================================================================== !
    epsCG                 = 1e-6   ! Stopping criterion (residual) of iterative CG solver (default that is used for the HDGSEM solver)
    maxIterCG             = 10000  ! Maximum number of iterations
    IniExactFunc          = 0      ! Initial field condition. 0: zero solution vector
    PrecondType           = 2

For this tutroial we use a PETSc Solver where a multitude of different numerical methods to solve the resulting system of linear equations is given by the implemented PETSc library and `epsCG` sets the abort residual of the CG solver, `maxIterCG` sets the maximum number of iterations within the CG solver, the `PrecondType` is set to 2, `IniExactFunc` sets the initial solution of the field solver (here 0 says that nothing is selected).

The PIC parameters for interpolation (of electric fields to the particle positions) and deposition (mapping of charge properties
from particle locations to the grid) are selected via

    ! =============================================================================== !
    ! PIC: Interpolation/Deposition
    ! =============================================================================== !
    PIC-DoInterpolation      = T                   ! Activate Lorentz forces acting on charged particles
    PIC-Interpolation-Type   = particle_position   ! Field interpolation method for Lorentz force calculation

    PIC-Deposition-Type      = cell_volweight_mean ! Linear deposition method



where the interpolation type `PIC-Interpolation-Type = particle_position` is currently the only option for specifying how
electro(-magnetic) fields are interpolated to the position of the charged particles.
For charge and current deposition,
`PIC-Deposition-Type = cell_volweight_mean` is selected.
The different available deposition types are described in more detail in Section {ref}`sec:PIC-deposition`.


#### Particle solver-2


The numerical scheme for tracking the movement of all particles throughout the simulation domain can be switched by

    ! =============================================================================== !
    ! Particle Solver
    ! =============================================================================== !
    TrackingMethod    = triatracking ! Particle tracking method

For the treatment of particles, the number of particle species `Part-nSpecies` that are used in the simulation (created initially or during the simulation time
through chemical reactions).

    ! =============================================================================== !
    ! PARTICLE Emission
    ! =============================================================================== !
    Part-nSpecies             = 3    ! Number of particle species

The inserting (sometimes labelled emission or initialization) of particles at the beginning or during the course of the simulation
is controlled via the following parameter. Here, only
the parameters for the electrons are shown, however, the parameters for the backgroung gas and ions are included in the supplied parameter.ini.
For each species a weighting factor (`Part-SpeciesX-MacroParticleFactor`) have to be defined.

    ! -------------------------------------
    ! Electrons 1
    ! -------------------------------------
    Part-Species3-Init1-PartDensity   = 2e20
    Part-Species3-MacroParticleFactor = 1               ! Weighting factor for species #3
    Part-Species3-nInits              = 1               ! Number of initialization/emission regions for species #3

The number of initialization sets is defined by `Part-Species[$]-nInits`, where each initialization set is accompanied
by a block of parameters that is preceded by the corresponding `-Init[$]` counter. In this example we have a single initialization set. The `Part-Species[$]-Init[$]-SpaceIC = cuboid` flag defines the type of the initialization set. Here, the particles are placed in a cuboid which is spanned by its base plane (`Part-Species[$]-Init[$]-BasePointIC`, `Part-Species[$]-Init[$]-BaseVector1IC`, `Part-Species[$]-Init[$]-BaseVector2IC`), a normal (`Part-Species[$]-Init[$]-NormalIC`) and its height (`Part-Species[$]-Init[$]-CuboidHeightIC`). Each type of the initialization set might have a different set of parameters and an overview is given in Section
{ref}`sec:particle-initialization-and-emission`.

    Part-Species3-Init1-SpaceIC               = cuboid   !(cartesian)
    Part-Species3-Init1-BasePointIC           =(/0.0009 , 0. , 0./)
    Part-Species3-Init1-BaseVector1IC         =(/0, 0.00005 , 0./)
    Part-Species3-Init1-BaseVector2IC         =(/0, 0. , 0.00005/)
    Part-Species3-Init1-NormalIC              =(/1, 0. , 0./)
    Part-Species3-Init1-CuboidHeightIC        = 1e-7
    Part-Species3-Init1-velocityDistribution  = maxwell_lpn !constant( unrealistic temp)
    Part-Species3-Init1-MWTemperatureIC       = 298

    ! -------------------------------------
    ! Ions - N2+
    ! -------------------------------------
    Part-Species2-MacroParticleFactor = 1.              ! Weighting factor for species #2
    Part-Species2-nInits              = 0               ! Number of initialization/emission regions for species #1

In order to reduce the simulation time we use the flag

    Part-vMPF                           = T
    Part-Species2-vMPFMergeThreshold    = 300
    Part-Species3-vMPFMergeThreshold    = 300


The merge routine randomly deletes particles until the desired number of particles is reached and the weighting factor is adopted accordingly.


Since we are dealing with ionization reaction in this setup, DSMC has to be enabled (`UseDSMC = T`). `Particles-DSMC-CollisMode` is an important parameterwhich has to beset to 3 so that chemistry is enabled.

    !===============================================================================
    ! DSMC
    !===============================================================================
    UseDSMC                           = T
    Particles-DSMC-CollisMode         = 3



After the chemistry is enabled, the model for a reaction has to be chosen. In our setup, it would be the cross section collision model. We only have one ionization reaction due to collision occuring in our setup  $$ e^- + N_2 \rightarrow 2e^- + N_2^+ $$ we will then define the reaction by

    DSMC-NumOfReactions                       = 1
    DSMC-Reaction1-ReactionModel              = XSec
    DSMC-Reaction1-Reactants                  = (/1,3,0/)
    DSMC-Reaction1-Products                   = (/2,3,3,0/)


For more information regarding different chemistry models refer to the section {ref}`sec:DSMC-chemistry`.


As we are dealing with a particle model of the streamer case for modelling of particle collisions with the Particle-in-Cell method,
often the Monte Carlo Collision (MCC) algorithm is utilized.
Here, experimentally measured or ab-initio calculated cross-sections are typically utilized to determine the collision probability,
based on the cross-section [m$^2$], the timestep [s], the particle velocity [m/s] and the target gas number density [m$^{-3}$].
To activate the MCC procedure, the collision cross-sections have to be supplied via read-in from a database.


    Particles-CollXSec-Database               = XSec_Database_N2_Ionized.h5
    Particles-CollXSec-NullCollision          = T

In PICLas, the null collision method reduces the computational effort as not every particle has to be checked for a collision.
Cross-section data can be retrieved from the LXCat database and converted with a Python script provided in the tools folder: piclas/tools/crosssection_database. Refer to the section {ref}`sssec:tools-maintain-database-xsec-collision` for more information.

For more information regarding cross-section based collision probability refer to the section {ref}`sec:background-gas`.


#### Field Boundaries-2

To accelarate the electrons towards the left boundary of the domain, i.e x=0, external electric field is provided.



    ! =============================================================================== !
    ! Field Boundaries
    ! =============================================================================== !
    BoundaryName = BC_Xleft
    BoundaryType = (/12,0/)                ! Neumann with fixed E

    BoundaryName = BC_Xright
    BoundaryType = (/4,0/)                 ! 4: Dirichlet with zero potential


#### Analysis setup-2

Finally, some parameters for run-time analysis are chosen by setting them `T` (true). Further, with `TimeStampLength = 18`, the names of the output files are shortened for better postprocessing. If this is not done, e.g. Paraview does not sort the files correctly and will display faulty behaviour over time.

    ! =============================================================================== !
    ! Analysis
    ! =============================================================================== !
    TimeStampLength          = 18 ! Reduces the length of the timestamps in filenames for better postprocessing
    CalcCharge               = T ! writes rel/abs charge error to PartAnalyze.csv
    CalcPotentialEnergy      = T ! writes the potential field energy to FieldAnalyze.csv
    CalcKineticEnergy        = T ! writes the kinetic energy of all particle species to PartAnalyze.csv
    PIC-OutputSource         = T ! writes the deposited charge (RHS of Poisson's equation to XXX_State_000.0000XXX.h5)
    CalcPICTimeStep          = T ! writes the PIC time step restriction to XXX_State_000.0000XXX.h5 (rule of thumb)
    CalcPointsPerDebyeLength = T ! writes the PIC grid step restriction to XXX_State_000.0000XXX.h5 (rule of thumb)
    CalcTotalEnergy          = T ! writes the total energy of the system to PartAnalyze.csv (field and particle)
    CalcInternalEnergy       = T
    CalcElectronTemperature  = T
    CalcDebyeLength          = T
    CalcPlasmaFrequency      = T
    CalcNumDens              = T
    CalcNumSpec              = T

The function of each parameter is given in the code comments. Information regarding every parameter can be obtained from
running the command

    piclas --help "CalcCharge"

where each parameter is simply supplied to the *help* module of **piclas**. This help module can also output the complete set of
parameters via `piclas --help` or a subset of them by supplying a section, e.g., `piclas --help "HDG"` for the HDGSEM solver.


#### Running the code-2

The command

    ./piclas parameter.ini DSMC.ini | tee std.out

executes the code and dumps all output into the file *std.out*.
To reduce the computation time, the simulation can be run using the Message Passing Interface (MPI) on multiple cores, in this case 8

    mpirun -np 8 piclas parameter.ini | tee std.out

If the run has completed successfully, which should take only a brief moment, the contents of the working folder should look like

    4.0K drwxrwxr-x  4.0K Jun 28 13:07 ./
    4.0K drwxrwxr-x  4.0K Jun 25 23:56 ../
    8.0K -rw-rw-r--  5.8K Jun 28 12:51 ElemTimeStatistics.csv
    120K -rw-rw-r--  113K Jun 28 12:51 FieldAnalyze.csv
    4.0K -rw-rw-r--  2.1K Jun 26 16:49 hopr.ini
    8.0K -rw-rw-r--  5.0K Jun 28 13:07 parameter.ini
    156K -rw-rw-r--  151K Jun 28 12:51 PartAnalyze.csv
     32K -rw-rw-r--   32K Jun 26 16:43 streamer_mesh.h5
    1.6M -rw-rw-r--  1.6M Jun 28 12:44 streamer_N2_State_000.00000000000000.h5
    ......
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 streamer_N2_State_000.00000000100000.h5

    1.6M -rw-rw-r--  1.6M Jun 28 12:51 streamer_N2_DSMCState_000.00000000001000.h5
    .......
    1.6M -rw-rw-r--  1.6M Jun 28 12:51 streamer_N2_DSMCState_000.00000000100000.h5

    72K -rw-rw-r--   71K Jun 28 12:51 std.out

Multiple additional files have been created, which are are named  **Projectname_State_Timestamp.h5**.
They contain the solution vector of the equation system variables at each interpolation nodes at the given time, which corresponds
to multiples of **Analyze_dt**. If these files are not present, something went wrong during the execution of **piclas**.
In that case, check the `std.out` file for an error message.

After a successful completion, the last lines in this file should look as shown below:

    ------------------------------------------------------------------------------------------------------------------------------------
    Sys date  :    03.10.2024 09:23:48
    PID: CALCULATION TIME PER TSTEP/DOF: [ 8.60746E-05 sec ]
    EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]: [ 1.52911E-11 sec/h ]
    Timestep  :    1.0000000E-13
    #Timesteps :    1.0000000E+04
    WRITE STATE TO HDF5 FILE [streamer_N2_State_000.00000000100000.h5] ... DONE [ 0.15 sec ] [ 0:00:00:00 ]
    #Particles :    1.2438110E+06    Average particles per proc :    1.5547638E+05    Min :    8.4251000E+04    Max :    5.1103400E+05
    ------------------------------------------------------------------------------------------------------------------------------------
    #Particles :    1.2461060E+06 (peak)         Average (peak) :    1.5576325E+05    Min :    0.0000000E+00    Max :    1.0267060E+06
    ====================================================================================================================================
    PICLAS FINISHED! [ 26482.49 sec ] [ 0:07:21:22 ]
    ====================================================================================================================================


### Visualization (post-processing)-2

To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**, **VisIt** or any other
visualisation tool for which the program **piclas2vtk** is used.

The parameters for **piclas2vtk** are stored in the **parameter.ini** file under

    ! =============================================================================== !
    ! piclas2vtk
    ! =============================================================================== !
    NVisu         = 3! Polynomial degree used for the visualization when the .h5 file is converted to .vtu/.vtk format. Should be at least N+1
    VisuParticles = T  ! Activate the conversion of particles from .h5 to .vtu/.vtk format. Particles will be displayed as a point cloud with properties, such as velocity, species ID, etc.

where `NVisu` is the polynomial visualization degree on which the field solution is interpolated.
Depending on the used polynomial degree `N` and subsequently the degree of visualization `NVisu`, which should always be higher than
`N`.
Additionally, the flag `VisuParticles` activates the output of particle position, velocity and species to the *vtk*-files.

Run the command

    ./piclas2vtk parameter.ini streamer_N2_State_000.00000000*

to generate the corresponding *vtk*-files, which can then be loaded into the visualisation tool.

The electron and ion densities and electric field at t=0.7 ns is shown in the plots.


```{figure} results/stream.png
---
name: fig:streamer-density
width: 700px
---

Electron and ion density in a planar front in nitrogen gas.
```

```{figure} results/electric.png
---
name: fig:streamer-mcc
width: 700px
---


Electric Field across the domain.
```


### Cross-section Database


In this simulation setup, we have integrated a cross-section database that focuses solely on the effective and ionization cross sections from the LXCat database. These two categories of cross-sections are important because they dictate the formation rates and types of species that result from electron collision reactions with nitrogen gas. The effective cross-section represents the probability of any electron-induced reaction occurring, while the ionization cross-section specifically addresses the formation of electrons and N₂+ ions.

For our application in PICLas, we particularly need cross-section data for electron interactions with nitrogen gas molecules. This data was extracted from the LXCat database and converted into an accessible HDF5 format using a Python script located in the tools folder: piclas/tools/crosssection_database.

The resulting HDF5 file, `XSec_Database_N2_Ionized.h5`, referenced in the parameter.ini file, contains energy-dependent cross-section values. These values indicate the specific energies at which effective and ionization collisions occur, providing critical data for predicting species formation due to these interactions.

```{figure} results/lxcat.png
---
name: fig:database
width: 700px
---

Cross-section data for electron interactions with nitrogen (N₂) relevant to negative streamer development, depicting effective and ionization cross-sections across electron energy levels. These values are essential for modeling electron impact ionization and excitation processes in nitrogen, influencing the formation and propagation of negative streamers.
```


## Comparison of Electron Fluid and PIC-MCC

```{figure} results/comparison_elec.png
---
name: fig:elec
width: 700px
---

Electron Density
```


```{figure} results/comparison_ion.png
---
name: fig:densities
width: 700px
---

Ion Density
```

```{figure} results/compare_elec.png
---
name: fig:elecfields
width: 700px
---

Electric Field
```

