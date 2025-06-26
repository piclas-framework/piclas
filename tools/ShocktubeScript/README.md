# Introduction:

The ShockTube script, "ShockedPICLas," can be used to simulate one-dimensional (1D) standing shock waves using PICLas as they occur in shock tubes. The script automatically builds meshes and parameter files for simulations, runs the simulations, and evaluates them.

This is done in an iterative process with two stages. In the first stage, the shock speed is iterated. The incoming flow hits a standing wall. A shock then forms in front of the wall and moves toward the inlet boundary. The inflow speed is varied to match the desired shock speed. When transitioning to the second stage, the post-shock condition is extracted from the simulation. In the second stage, the reference view changes — the script now tries to produce a standing shock with the post-shock conditions found in the first stage. The flow velocity continues to be adjusted until the shock remains stationary within the computational domain. Then, it samples the result until it has low noise.

# Usage of the Script

To use the script, one must know how to conduct simulations with PICLas and which parameters are necessary. The script forwards most of the simulation parameters to the PICLas simulations unchanged, but some are adjusted based on the simulation case and enabled autoadjusts of simulation parameters.

## Setup of a Simulation

| Parameter                        | Data Type | Default | Description                                                                                                   |
|----------------------------------|-----------|---------|---------------------------------------------------------------------------------------------------------------|
| Shock-Target-Velo                | float     |         | Intended Velocity of the Shock                                                                                |
| Shock-Target-dev                 | float     |         | Allowed deviation of the target velocity                                                                      |
| Shock-Simulation-dt              | float     |         | dt of the simulation, if autoadjust then dt of the first simulation                                           |
| Shock-Simulation-tend            | float     |         | tend of the simulation, if autoadjust then tend of the first simulation                                       |
| Shock-Simulation-MPF             | float     |         | MPF of the simulation, if autoadjust then MPF of the first simulation                                         |
| Shock-Simulation-MaxPartNum      | int       | -1      | MaxPartNum of the simulation, if autoadjust then MaxPartNum of the first simulation, deprecated disabled by -1|
| Shock-Simulation-Median          | float     | 0.01    | allowed relative median in Integrated File (smoothness condition of overall solution)                         |
| Shock-Mesh-dxElem                | float     |         | x extend of a Element/Cell                                                                                    |
| Shock-Mesh-Length                | float     |         | length of the simulation, if autoadjust then length of the first simulation                                   |
| Shock-Mesh-HalfWidth             | float     | 0.05    | half y and z extend of the mesh                                                                               |
| Shock-Start-Velo                 | float     | 0.0     | Initial velocity of the flow field. If 0 then "Shock-Target-Velo" is used                                     |
| Shock-ResponseRatio              | float     | 1.0     | = delta Inflowvelocity/delta Shockspeed, only for continuation of simulations                                 |
| Shock-KeepAllSimulationData      | bool      | False   | If true: for every simulation a new folder is created. If false: last simulation data will be deleted         |
| Shock-UseMPI                     | bool      | False   | If piclas should run with mpi                                                                                 |
| Shock-nMPICores                  | int       | 1       | Only if "Shock-UseMPI": number of used cores                                                                  |
| Shock-PiclasPath                 | str       |         | path to the piclas executeable or cl command                                                                  |
| Shock-HoprPath                   | str       |         | path to the hopr executeable or cl command                                                                    |
| Shock-DSMCSpecies                | str       |         | Name of the DSMCSpecies file. has not to be set if DSMCSpecies file is not needed                             |
| Shock-AdditionalCopyFiles        | str       |         | Only if "KeepAllSimulationData": additional neccessary files like electronic data base                        |
| Shock-Simulation-ElemsPerWrite   | float     | 1       | Average travel distance of the shock in Elems between piclas macrowrites                                      |
| Shock-Simulation-ElectronSpecies | int       | -1      | Species number of electrons =-1 if there is no electron species                                               |
| Shock-Autoadjust-Mesh            | bool      | False   | If true: mesh length will be adjusted by last simulation results, might be unstable                           |
| Shock-Autoadjust-dt              | bool      | False   | If true: dt will be adjusted by last simulation results.                                                      |
| Shock-Autoadjust-tend            | bool      | False   | If true: tend will be adjusted by last simulation results.                                                    |
| Shock-Autoadjust-maxPartNum      | bool      | False   | If true: maxpartnum will be adjusted by last simulation results, deprecated                                   |
| Shock-Autoadjust-MPF             | bool      | False   | If true: MPF will be adjusted by last simulation results.                                                     |
| Shock-PartperElemmin             | float     | 10      | Only if "Shock-Autoadjust-MPF": lowest Simpartnum allowed in Elems                                            |
| Shock-MCSoverMFPmax              | float     | 0.2     | Only if "Shock-Autoadjust-MPF": highest MCSoverMFP allowed in Elems                                           |
| Shock-VeloTempEps                | float     | 5       | Only if "Shock-Autoadjust-dt": maximum considered multiple thermal velocity of the mean thermal velocity      |
| Shock-MaxCollProbMax             | float     | 0.9     | Only if "Shock-Autoadjust-dt": highest MaxCollProbMax allowed in Elems                                        |

The end time and length should be within a reasonable range. For example, the shock should not touch the inflow boundary, and there should be enough time for a shock to develop.

All other PICLas simulation parameters must be set in this file and will be passed on to the simulation. The inflow composition is defined by setting it as init1. The script will set the corresponding velocity and surface flux based on the iterated velocities. As in a normal PICLas simulation, all species parameters must be defined in the DSMCSpecies.ini. It has to be tested to see if it works with the Species Database.


## Example: EAST-50-29

This is the example case. It was a test in the EAST shock tube facility [[1]](#1) which is also an ESA reference case [[2]](#2). The shock tube operates with air and reaches a shock speed of 10,300 m/s, meaning a significant amount of ionization occurs. During the test, the velocity of the shock front was measured, and spectra of the emitted light were recorded at different wavelengths. PICLas was used to simulate the flow field, the PICLas radiation solver generated the resulting spectrum. The results are published in [[3]](#3). The necessary files to conduct the simulation are appended. This section only explains the SHocktube.ini file, which is an extended parameter.ini file of a PICLas simulation.

The shocktube script related parameters are set to:
```
! =============================================================================== !
! Shock Simulation
! =============================================================================== !
Shock-Target-Velo=10300                                                 ! Desired velocity of the shock
Shock-Target-dev=10.0                                                   ! Desired acuracy of the traget velocity
Shock-Simulation-dt=1.230767e-09                                        ! Time step of the simulation, start time step for autoadjust
Shock-Simulation-tend=4.236113e-04                                      ! End time of the simulation, start time step for autoadjust
Shock-Simulation-MPF=3.005187e+16                                       ! Particle Weight for the simulation, start time step for autoadjust (vmpf is not supported)
Shock-Simulation-ElemsPerWrite=10                                       ! How many elements shall the shock front propagate per output (noise,memory)
Shock-Simulation-ElectronSpecies=11                                     ! Species 11 are the electrons
Shock-Mesh-dxElem=0.0005                                                ! x-Resolution of the Mesh
Shock-Mesh-Length=0.3                                                   ! Length of the mesh
Shock-Start-Velo=9587.299262                                            ! The start velocity is set to the the final value or a fast convergence, typically it is not known and has to be iterated by the script
Shock-KeepAllSimulationData=true                                        ! The data of the iterations are not deleted
Shock-UseMPI=true                                                       ! MPI is enabled
Shock-nMPICores=10                                                      ! 10 Cores are used for the simulation
Shock-PiclasPath=~/master.dev/build/bin/piclas                          ! Path to piclas executable
Shock-HoprPath=hopr                                                     ! Path to hopr executable
Shock-DSMCSpecies=DSMCSpecies_Earth_11Spec.ini                          ! DSMC.ini file
Shock-AdditionalCopyFiles=DSMCSpecies_electronic_state_full_Data.h5     ! The electronic level database is needed to calculate the backward reaction rates and electrical exictations
Shock-Autoadjust-Mesh=false                                             ! The mesh will not be adjusted due to instabillities
Shock-Autoadjust-dt=true                                                ! The simulation time step will be adjusted, based on the maximum collision probability and travel distance of the particles
Shock-Autoadjust-tend=true                                              ! The end time of the simulation will be adjusted based on the location of the shock front or the smoothness of the solution
Shock-Autoadjust-MPF=true                                               ! The particle weight will be adjusted, based on the MCS over MPF and minimal particle number
Shock-partperelemmin=10                                                 ! Minimum 10 particles per computational element
```

The initialization of particles is set to the following parameters:

```
! =============================================================================== !
! Piclas-parameters copied to the shocktube simulations
! Init1 are the pre-shock conditions
! =============================================================================== !

Part-nSpecies=11
! =============================================================================== !
! Species1 - N
! =============================================================================== !
Part-Species1-MassIC=2.32600E-26        ! N Molecular Mass
! =============================================================================== !
! Species2 - O
! =============================================================================== !
Part-Species2-MassIC=2.65700E-26         ! O Molecular Mass
! =============================================================================== !
! Species3 - N2
! =============================================================================== !
Part-Species3-MassIC=4.65200E-26         ! N2 Molecular Mass

Part-Species3-nInits=1

Part-Species3-Init1-velocityDistribution=maxwell
Part-Species3-Init1-VeloVecIC=(/1,0.,0/)
Part-Species3-Init1-MWTemperatureIC=300
Part-Species3-Init1-TempVib=300
Part-Species3-Init1-TempRot=300
Part-Species3-Init1-TempElec=300
Part-Species3-Init1-PartDensity=5.15217391304348E+020
! =============================================================================== !
! Species4 - O2
! =============================================================================== !
Part-Species4-MassIC=5.31400E-26        ! O2 Molecular Mass

Part-Species4-nInits=1

Part-Species4-Init1-velocityDistribution=maxwell
Part-Species4-Init1-VeloVecIC=(/1,0.,0/)
Part-Species4-Init1-MWTemperatureIC=300
Part-Species4-Init1-TempVib=300
Part-Species4-Init1-TempRot=300
Part-Species4-Init1-TempElec=300
Part-Species4-Init1-PartDensity=1.3695652173913E+020
! =============================================================================== !
! Species5 - NOexitations
! =============================================================================== !
Part-Species5-MassIC=4.98300E-26          ! NO Molecular Mass
! =============================================================================== !
! Species6 - N+
! =============================================================================== !
Part-Species6-MassIC=2.3259089061644E-26        ! N Molecular Mass
Part-Species6-ChargeIC=1.60217653E-19
! =============================================================================== !
! Species7 - O+
! =============================================================================== !
Part-Species7-MassIC=2.6569089061644E-26         ! O Molecular Mass
Part-Species7-ChargeIC=1.60217653E-19
! =============================================================================== !
! Species8 - N2+
! =============================================================================== !
Part-Species8-MassIC=4.6519089061644E-26          ! N2 Molecular Mass
Part-Species8-ChargeIC=1.60217653E-19
! =============================================================================== !
! Species9 - O2+
! =============================================================================== !
Part-Species9-MassIC=5.3137267184932E-26        ! O2 Molecular Mass
Part-Species9-ChargeIC=1.60217653E-19
! =============================================================================== !
! Species10 - NO+
! =============================================================================== !
Part-Species10-MassIC=4.9828178123288E-26          ! NO Molecular Mass
Part-Species10-ChargeIC=1.60217653E-19
! =============================================================================== !
! Species11 - e
! =============================================================================== !
Part-Species11-MassIC=9.10938356E-31        ! e Mass
Part-Species11-ChargeIC=-1.60217653E-19
```
The mass and non-zero charge is defined for all species. The inflow composition is defined as init1; so the inflow consists of only nitrogen and oxygen molecules. All other species will be formed in and behind the shock front. The script will set the corresponding velocity and surface flux.

Other simulation parameters are forwarded directly to the PICLas simulation and are set to:
```
! =============================================================================== !
! MESH
! =============================================================================== !
TrackingMethod=triatracking                                             ! used tracking method, fast
! =============================================================================== !
! OUTPUT / VISUALIZATION
! =============================================================================== !
ProjectName   = EAST-50-29
IterDisplayStep  = 100exitations
Logging       = F
! =============================================================================== !
! CALCULATION
! =============================================================================== !
CFLscale   = 0.2  ! Scaling of theoretical CFL number

! =============================================================================== !
! DSMC
! =============================================================================== !
UseDSMC=true                                                            ! DSMC collisions are necessary to form a shock
Particles-DSMC-CollisMode=3 !(1:elast coll, 2: elast + rela, 3:chem)    ! chemical reactions re desired
Part-NumberOfRandomSeeds=2
Particles-RandomSeed1=1
Particles-RandomSeed2=2
Particles-ModelForVibrationEnergy=0 !(0:SHO, 1:TSHO)                    ! Simple harmonic ocillator as vibration model
Particles-DSMC-UseNearestNeighbour=true                                 ! nearest neighbor search of the second collision partner
Particles-DSMC-UseOctree=true                                           ! octree is used for fast pairing
Particles-MPIWeight=1000                                                ! for MPI simulations of the second stage
Particles-HaloEpsVelo=100000000                                         ! reduces particle loss to consider the whole computational domain
Particles-DSMC-CalcQualityFactors=true                                  ! needed for autoadjusts
Particles-DSMC-BackwardReacRate = true                                  ! automatic calculation of backward reaction rates
Particles-DSMC-PartitionMaxTemp = 120000.                               ! above maximum expected temperature for chemical reaction
Particles-DSMC-PartitionInterval= 20.                                   ! resolution of chemical reaction database
Particles-DSMC-AmbipolarDiffusion = T                                   ! electrons are pinned to heavy particles, no calculation of electrical field necessary, much faster simulation
Particles-DSMC-ElectronicModel  = 2                                     ! Every particle carries an distribution of electrical exictations, better resolution of exitation states for radiation solver
Particles-DSMCElectronicDatabase = DSMCSpecies_electronic_state_full_Data.h5 ! Electronic level database
EpsMergeElectronicState = 1E-2                                          ! pinning of electrical states for faster simulations
```

## Conducting the simulation

The script works with Python 2 and 3, but has only been tested with Python 3. The script is executed with:
```
python3 Shocktube.py [-h] [--Equicon] [--restart RESTART] InputFile

positional arguments:
  InputFile             File containing input parameters

options:
  -h, --help            show this help message and exit
  --Equicon             Calculate only equilibrium condition
  --restart RESTART, -r RESTART
                        Restart from the given simulation (directory)
```

The code creates a directory for the simulation and copies or creates all the necessary files in the folder. Then, it builds the mesh according to the parameters and starts the simulation. Once the simulation is finished, it extracts the shock speed. If the shock speed is outside the desired range around the target speed, the inflow speed is updated with a second-order interpolation or extrapolation, and the next simulation is started. If autoadjust is enabled, the simulation parameters are also changed for optimal simulation, as defined by the autoadjust parameters.

If the shock speed meets the criteria, the script extracts the post-shock condition for the right boundary for the standing iterations. There, the inflow velocity is fine-tuned to create a perfectly standing shock after some initial oscillations. The criterion is that the shock moves less than half an element width over the whole simulation duration. All outputs after the initial oscillations are combined into a single output. The result is checked for smoothness, and the simulation time is extended if the criterion is not met.

# Implementation

Some notes about the implementation of core elements of the script

## Determination of the shock speed


This follows the Higdon algorithm closely [[4]](#4). The speed of the pressure jump is extracted since it combines all species data and is the fastest parameter to reach equilibrium in the post-shock condition. During the initial development of the shock, some non-linearities appear; these must be excluded. This is accomplished through a series of second-order fits of the shock position; each fit excludes an additional writeout, beginning with the earliest one. When the sign of the quadratic term changes, it is assumed that the noise in the remaining signal is greater than the initial non-linearities. The velocity is the slope of the linear fit of the remaining writeouts.

## Extraction of the Equilibrium condition

First, the post-shock region is determined by either the element with the middle number density or the element with the highest translational temperature. The criterion resulting in the shortest region is used. Then, a fourth-order polynomial or an exponential function is fitted to this region for each parameter. The equilibrium value of the parameter (e.g., the number density of species 1) is the middle value of the turning points of the polynomial or the offset of the exponential function. The initial phase of shock generation produces disturbances in the post-shock region near the wall. These disturbances must be filtered by the fit function.

# References
<a id="1">[1]</a>
Aaron M. Brandis and Brett A. Cruden, “Benchmark shock tube experiments of radiative heating relevant to earth re-entry,” in 55th AIAA Aerospace Sciences Meeting. jan 2017, American Institute of Aeronautics and Astronautics.

<a id="2">[2]</a>
Brett Cruden, “Test-case for radiative heating in lunar return-similar shock tube tests,” in 8th International Workshop on Radiation of High Temperature Gases in Atmospheric Entry. European Space Agency, Mar. 2019.

<a id="3">[3]</a>
Raphael Tietz, Julian Beyer, Marcel Pfeiffer, and Stefanos Fasoulas, "SHOCK TUBE SIMULATIONS WITH THE PIC-DSMC CODE PICLAS," in 2nd International Conference on Flight Vehicles, Aerothermodynamics and Re-entry, Heilbronn, Germany, Jun. 2022

<a id="4">[4]</a>
Kyle J. Higdon, Monte Carlo sensitivity analyses of DSMC parameters for ionizing hypersonic flows, phdthesis, University of Texas, Aug. 2018., pages 18ff