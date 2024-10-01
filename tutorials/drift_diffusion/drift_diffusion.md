# Drift Diffusion Model Tutorial

Streamers, as seen in lightning and sprite events, are also observed in industrial applications such as lighting, purification of dirty gases and water, plasma jets and bullets for disinfection, and plasma-assisted combustion. Using this setup, a 1D streamer simulation will be carried out using PICLas.

Before beginning with the tutorial, copy the “drift _diffusion directory” from the tutorial folder in the top level directory to a separate location.

    cp -r $PICLAS_PATH/tutorials/drift_diffusion .
    cd drift_diffusion

## Mesh Generation with HOPR (pre-processing)

Before the actual simulation is conducted, a mesh file in the correct HDF5 format has to be supplied. The mesh files used by piclas are created by supplying an input file hopr.ini with the required information for a mesh that has either been created by an external mesh generator or directly from block-structured information in the hopr.ini file itself. Here, a block-structured grid is created directly from the information in the hopr.ini file. 

To create the .h5 mesh file, simply run

    hopr hopr_streamer.ini

This creates the mesh file *streamer_small_1500_mesh.h5* in HDF5 format.
Alternatively, if you do not want to run **hopr** yourself, you can also use the provided mesh. The only difference from the mesh created in previous tutorial, is the size of the simulation domain in this tutorial which is set to [$4\pi\times1\times1$] m$^{3}$ and is defined by the single block information

in the line, where each node of the hexahedral element is defined

    Corner     = (/0.,0.,0.,,l1,0.,0.,,l1,l2,0.,,0.,l2,0.,,0.,0.,l2,,l1,0.,l2,,l1,l2,l2,,0.,l2,l2/)

The number of mesh elements for the block in each direction can be adjusted by changing the line

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
    VV=(/0. , 0. , l2/)   ! Displacement vector 2 (z-direction
     



## Simulation with PICLas

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
    PICLAS_PETSC            = ON
    PICLAS_READIN_CONSTANTS = ON
    PICLAS_TIMEDISCMETHOD   = Explicit-FV


using the ccmake (gui for cmake) or simply run the following command from inside the build directory

    cmake .. -DPICLAS_EQNSYSNAME=drift_diffusion -DPICLAS_TIMEDISCMETHOD=Explicit-FV -D PICLAS_PETSC=ON

to configure the build process and run

    make -j


afterwards to compile the executable. For this setup, we have chosen the drift_diffusion solver and selected the Explicit-FV time discretization method. An overview over the available solver and discretization options is given in Section Solver settings. To run the simulation and analyse the results, the piclas and piclas2vtk executables have to be run. To avoid having to use the entire file path, you can either set aliases for both, copy them to your local tutorial directory or create a link to the files via

    ln -s $PICLAS_PATH/build_drift_diffusion/bin/piclas
    ln -s $PICLAS_PATH/build_drift_diffusion/bin/piclas2vtk

where the variable PICLAS_PATH contains the path to the location of the piclas repository. If the piclas repository is located in the home directory, the two commands

    ln -s /home/$(whoami)/piclas/build_drift_diffusion/bin/piclas
    ln -s /home/$(whoami)/piclas/build_drift_diffusion/bin/piclas2vtk


### General numerical setup

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
    tend            = 1.e-9 
    Analyze_dt      = 1.e-11  
    IterDisplayStep = 100

where the time step for the field and particle solver is set via `ManualTimeStep`, the final simulation time `tend`, the time
between restart/checkpoint file output `Analyze_dt` (also the output time for specific analysis functions) and the number of time
step iterations `IterDisplayStep` between information output regarding the current status of the simulation that is written to std.out.
The remaining parameters are selected for the field and particle solver as well as run-time analysis.

### Boundary conditions

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

 

### Field solver

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

### Particle solver

For the treatment of particles, the maximum number of particles `Part-maxParticleNumber` that each processor can hold has to be supplied and
the number of particle species `Part-nSpecies` that are used in the simulation (created initially or during the simulation time
through chemical reactions).

    ! =============================================================================== !
    ! PARTICLE Emission
    ! =============================================================================== !
    Particles-Species-Database = ../../piclas/SpeciesDatabase.h5
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

    Particles-HaloEpsVelo = 1.E12

    

### Analysis setup

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

### Running the code

The command

    ./piclas parameter_streamer.ini | tee std.out

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

    --------------------------------------------------------------------------------------------
    Sys date  :    03.07.2021 14:34:26
    PID: CALCULATION TIME PER TSTEP/DOF: [ 1.21226E-04 sec ]
    EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]: [ 7.50886E+02 sec/h ]
    Timestep  :    1.0000000E-01
    #Timesteps :   5.0000000E+02
     WRITE STATE TO HDF5 FILE [streamer_N2_State_000.00000000100000.h5] ... DONE [ 0.02 sec ] [ 0:00:00:00 ]

    #Particles :    4.0100000E+04    Average particles per proc :    1.0025000E+04    Min :    9.3390000E+03    Max :    1.0702000E+04

    #Particles :    4.0100000E+04 (peak)         Average (peak) :    1.0025000E+04    Min :    9.1730000E+03    Max :    1.1010000E+04

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

    ./piclas2vtk parameter_streamer.ini streamer_N2_DSMCState_000.00000000**

to generate the corresponding *vtk*-files, which can then be loaded into the visualisation tool.












