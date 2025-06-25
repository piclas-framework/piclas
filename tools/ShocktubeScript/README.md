# Introduction:

The ShockTube script, "ShockedPICLas," can be used to simulate one-dimensional (1D) standing shock waves using PICLas as they occur in shock tubes. The script automatically builds meshes and parameter files for simulations, runs the simulations, and evaluates them.

This is done in an iterative process with two stages. In the first stage, the shock speed is iterated. The incoming flow hits a standing wall. A shock then forms in front of the wall and moves toward the inlet boundary. The inflow speed is varied to match the desired shock speed. When transitioning to the second stage, the post-shock condition is extracted from the simulation. In the second stage, the reference view changes — the script now tries to produce a standing shock with the post-shock conditions found in the first stage. The flow velocity continues to be adjusted until the shock remains stationary within the computational domain. Then, it samples the result until it has low noise.

# Usage of the Script

To use the script, one must know how to conduct simulations with PICLas and which parameters are necessary. The script forwards most of the simulation parameters to the PICLas simulations unchanged, but some are adjusted based on the simulation  case and enabled autoadjusts of simulation parameters.

## Example: EAST-50-29

This is the example case. It was a test in the EAST shock tube facility [[1]](#1) which is also an ESA reference case [[2]](#2). The results are published in [[3]](#3). The shock tube operates with air and reaches a shock speed of 10,300 m/s, meaning a significant amount of ionization occurs. During the test, the velocity of the shock front was measured, and spectra of the emitted light were recorded at different wavelengths.

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
| Shock-Autoadjust-Mesh            | bool      | False  [[3]](#3) | If true: mesh length will be adjusted by last simulation results, might be unstable                           |
| Shock-Autoadjust-dt              | bool      | False   | If true: dt will be adjusted by last simulation results.                                                      |
| Shock-Autoadjust-tend            | bool      | False   | If true: tend will be adjusted by last simulation results.                                                    |
| Shock-Autoadjust-maxPartNum      | bool      | False   | If true: maxpartnum will be adjusted by last simulation results, deprecated                                   |
| Shock-Autoadjust-MPF             | bool      | False   | If true: MPF will be adjusted by last simulation results.                                                     |
| Shock-PartperElemmin             | float     | 10      | Only if "Shock-Autoadjust-MPF": lowest Simpartnum allowed in Elems                                            |
| Shock-MCSoverMFPmax              | float     | 0.2     | Only if "Shock-Autoadjust-MPF": highest MCSoverMFP allowed in Elems                                           |
| Shock-VeloTempEps                | float     | 5       | Only if "Shock-Autoadjust-dt": maximum considered multiple thermal velocity of the mean thermal velocity      |
| Shock-MaxCollProbMax             | float     | 0.9     | Only if "Shock-Autoadjust-dt": highest MaxCollProbMax allowed in Elems                                        |

The end time and length should be within a reasonable range. For example, the shock should not touch the inflow boundary, and there should be enough time for a shock to develop.

All other PICLas simulation parameters must be set in this file and will be passed on to the simulation. The inflow composition is defined by setting it as init1, as shown in the example file:

```
! =============================================================================== !
! Species3 - N2
! =============================================================================== !
Part-Species3-nInits=1
Part-Species3-Init1-velocityDistribution=maxwell
Part-Species3-Init1-VeloVecIC=(/1,0.,0/)
Part-Species3-Init1-MWTemperatureIC=300
Part-Species3-Init1-TempVib=300
Part-Species3-Init1-TempRot=300
Part-Species3-Init1-TempElec=300
Part-Species3-Init1-PartDensity=5.15217391304348E+020
! =============================================================================== !
! Species4 - O2Shock-PartperElemmin             | float     | 10      | Only if "Shock-Autoadjust-MPF": lowest Simpartnum allowed in Elems                                            |
| Shock-MCSoverMFPmax              | float     | 0.2     | Only if "Shock-Autoadjust-MPF": highest MCSoverMFP allowed in Elems                                           |
| Shock-VeloTempEps                | float     | 5       | Only if "Shock-Autoadjust-dt": maximum considered multiple thermal velocity of the mean thermal velocity      |
| Shock-MaxCollProbMax
! =============================================================================== !
Part-Species4-nInits=1
Part-Species4-Init1-velocityDistribution=maxwell
Part-Species4-Init1-VeloVecIC=(/1,0.,0/)
Part-Species4-Init1-MWTemperatureIC=300
Part-Species4-Init1-TempVib=300
Part-Species4-Init1-TempRot=300
Part-Species4-Init1-TempElec=300
Part-Species4-Init1-PartDensity=1.3695652173913E+020
```
The script will set the corresponding velocity and surface flux. As in a normal PICLas simulation, all species parameters must be defined in the DSMCSpecies.ini. It has to be tested to see if it works with the Species Database.

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


This follows the Higdon algorithm closely [[4]](#4). The speed of the pressure jump is extracted since it combines all species data and is the fastest parameter to reach equilibrium in the post-shock condition. During the initial development of the shock, some nonlinearities appear; these must be excluded. This is accomplished through a series of second-order fits of the shock position; each fit excludes an additional writeout, beginning with the earliest one. When the sign of the quadratic term changes, it is assumed that the noise in the remaining signal is greater than the initial nonlinearities. The velocity is the slope of the linear fit of the remaining writeouts.

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