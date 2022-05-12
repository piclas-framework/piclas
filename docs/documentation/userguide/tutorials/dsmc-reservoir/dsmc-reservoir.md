(sec:tutorial-dsmc-reservoir)=
# Adiabatic Box/Reservoir (DSMC, Relaxation/Chemistry)
An essential feature of DSMC simulations is their ability to treat thermal and chemical non-equilibrium in a physically correct manner. A simple example for such a use case are adiabatic box/reservoir simulations, which start in a non-equilibrium state. The subsequent simulation should lead to the correct relaxation towards equilibrium. Hence, this tutorial provides an example case at thermal non-equilibrium with disabled chemistry and another based on chemical non-equilibrium with chemistry enabled.

Before beginning with the tutorial, copy the `dsmc-reservoir` directory from the tutorial folder in the top level directory to a separate location

    cp -r $PICLAS_PATH/tutorials/dsmc-reservoir .
    cd dsmc-reservoir

## Mesh Generation with HOPR (pre-processing)

Before the actual simulation is conducted, a mesh file in the HDF5 format has to be supplied. The mesh files used by **piclas** are created by supplying an input file *hopr.ini* with the required information for a mesh that has either been created by an external mesh generator or directly from block-structured information in the hopr.ini file itself. Here, a block-structured grid is created directly from the information in the hopr.ini file. To create the .h5 mesh file, simply run

    hopr hopr.ini

This creates the mesh file *dsmc_reservoir_mesh.h5* in HDF5 format, which is depicted in {numref}`fig:dsmc-reservoir-mesh-corners`. Alternatively, if you do not want to run **hopr** yourself, you can also use the provided mesh. After this step, the mesh needs to be copied in the simulation directories

    cp dsmc_reservoir_mesh.h5 ./chemistry-off/
    cp dsmc_reservoir_mesh.h5 ./chemistry-on/

The size of the simulation domain is set to [$\pu{4.64e-6}\times\pu{4.64e-6}\times\pu{4.64e-6}$] m$^{3}$ and is defined by the single block information, where each node of the hexahedral block is defined

    Corner         =   (/0.0,0.0,0.0,, 1.0,0.0,0.0,, ... /)

```{figure} mesh/dsmc-reservoir-mesh-corners.svg
---
name: fig:dsmc-reservoir-mesh-corners
width: 50%
---
Order of the corners to define the used mesh. The first node is placed at the origin.
```

Afterwards this element is scaled via

    postScaleMesh  = T
    meshScale      = 4.64E-6 

The number of mesh elements for the block in each direction can be adjusted by changing the line

    nElems         = (/2,2,1/)

However, for an adiabatic box/reservoir simulation, a single element would be sufficient/optimal. Each side of the block has to be assigned a boundary index, which corresponds to the boundaries defined in the next steps. Due to the fact, that all boundaries in this example should behave similar, only one index is needed.

    BCIndex        = (/1,1,1,1,1,1/)

The field boundaries can directly be defined in the hopr.ini file (contrary to the particle boundary conditions, which are defined
in the parameter.ini). For particle-based methods without electromagnetic fields, the boundary type is set in the parameter.ini.

    !=============================================================================== !
    ! BOUNDARY CONDITIONS
    !=============================================================================== !
      BoundaryName = BC_wall
      BoundaryType = (/4,0,0,0/)

For more information about hopr, visit [https://github.com/hopr-framework/hopr](https://github.com/hopr-framework/hopr).

## Simulation: Chemistry disabled

Install **piclas** by compiling the source code as described in Chapter {ref}`userguide/installation:Installation` and make sure to set the correct compile flags (ie. chose the correct simulation method)

    PICLAS_TIMEDISCMETHOD = RESERVOIR

or simply run the following command from inside the *build* directory

    cmake ../ -DPICLAS_TIMEDISCMETHOD=RESERVOIR

to configure the build process and run `make` afterwards to build the executable. For this setup, the reservoir method, which is based on the DSMC method, is needed to allow for reservoir specific settings. It is recommended to either utilize a separate build folder (e.g. build_DSMC/) or to delete the contents of the folder beforehand to avoid conflicts between different compiler options (e.g. the setting `PICLAS_EQNSYSNAME = poisson` from the plasma wave tutorial is in conflict with the DSMC method). An overview over the available solver and discretization options is given in Section {ref}`sec:solver-settings`. The physical parameters for this test case are summarized in {numref}`tab:dsmc_chem_off_phys`.

```{table} Physical properties at the simulation start
---
name: tab:dsmc_chem_off_phys
---
|                   Property                   |        Value        |
| -------------------------------------------- | :-----------------: |
| Species                                      | $\text{CO}_2$       |
| Molecule mass $m_{\text{CO}_2}$              | $\pu{7.306e-26 kg}$ |
| Number density $n_{\text{CO}_2}$             | $\pu{1e23 m^{-3}}$  |
| Translational temperature $T_{\text{trans}}$ | $\pu{10000 K}$      |
| Rotational temperature $T_{\text{rot}}$      | $\pu{7500 K}$       |
| Vibrational temperature $T_{\text{vib}}$     | $\pu{5000 K}$       |
```

### General numerical setup

The general numerical parameters are selected by the following

    ! =============================================================================== !
    ! MESH
    ! =============================================================================== !
    MeshFile   = dsmc_reservoir_mesh.h5
    ! =============================================================================== !
    ! PROJECT
    ! =============================================================================== !
    ProjectName     = dsmc_reservoir_chemisty_off
    TrackingMethod  = triatracking

where, the path to the mesh file `MeshFile`, project name and particle tracking method `TrackingMethod` are chosen. The temporal parameters of the simulation are controlled via

    ! =============================================================================== !
    ! CALCULATION
    ! =============================================================================== !
    ! Time
    TEnd                  = 2E-6
    ManualTimeStep        = 1.0E-8
    Analyze_dt            = 5E-6
    IterDisplayStep       = 100
    Particles-HaloEpsVelo = 5000

where the final simulation time `tend` [s], the time step for the field and particle solver is set via `ManualTimeStep` [s]. The time between restart/checkpoint file output is defined via `Analyze_dt` (which is also the output time for specific analysis functions in the field solver context). The number of time step iterations `IterDisplayStep` defines the interval between information output regarding the current status of the simulation, which is written to std.out. The `Particles-HaloEpsVelo` [m/s] determines the size of the halo region for MPI communication and should not be smaller than the fastest particles in the simulation.

(sec:tutorial-dsmc-analysis-setup)=
### Analysis setup 

For this case our focus is on the run-time analysis to investigate the transient behavior of the reservoir. The first parameter `Part-AnalyzeStep` allows to perform the output every N$^\text{th}$ iteration to reduce the size of the output file and to increase the computational speed. Different parameters for run-time analysis can be enabled, in this case the number of particles per species (`CalcNumSpec`) and the temperature output (`CalcTemp`). It is also recommended to enable `Particles-DSMC-CalcQualityFactors`, which provides outputs to evaluate the quality of the simulation results such as the mean and maximum collision probabilities. The parameter `TimeStampLength = 13` reduces the output filename length. It can be needed for postprocessing, as e.g. ParaView sometimes does not sort the files correctly if the timestamps are too long. The displayed time solution would then be faulty.

    ! =============================================================================== !
    ! Particle Analysis
    ! =============================================================================== !
    
    Part-AnalyzeStep = 1
    CalcNumSpec = T
    CalcTemp    = T
    Particles-DSMC-CalcQualityFactors = T
    TimeStampLength  = 13

All available options with a short description can be displayed using the help of PICLas:

    piclas --help 'Particle Analyze'

### Boundary conditions

The boundary conditions are set by the following lines

    ! =============================================================================== !
    ! Boundaries
    ! =============================================================================== !
    Part-nBounds              = 1
    Part-Boundary1-SourceName = BC_wall
    Part-Boundary1-Condition  = reflective

where, the number of boundaries `Part-nBounds` is followed by the names of the boundaries (given by the hopr.ini file) and the type `reflective`.

(sec:tutorial-dsmc-particle-solver)=
### Particle solver

For the treatment of particles, the maximum number of particles `Part-maxParticleNumber` that each processor can hold has to be supplied and the number of particle species `Part-nSpecies` that are used in the simulation (created initially or during the simulation time through chemical reactions).

    ! =============================================================================== !
    ! PARTICLES
    ! =============================================================================== !
    Part-maxParticleNumber = 420000
    Part-nSpecies          = 1
    Part-FIBGMdeltas       = (/4.64E-6,4.64E-6,4.64E-6/)

The inserting (sometimes labelled emission or initialization) of particles at the beginning or during the course of the simulation is controlled via the following parameters. For each species, the mass (`Part-Species[$]-MassIC`), charge (`Part-Species[$]-ChargeIC`) and weighting factor $w$ (`Part-Species[$]-MacroParticleFactor`) have to be defined.

    ! =============================================================================== !
    ! Species1 - CO2
    ! =============================================================================== !
    Part-Species1-MassIC              = 7.306E-26
    Part-Species1-ChargeIC            = 0
    Part-Species1-MacroParticleFactor = 5E2

The number of initialization sets is defined by `Part-Species[$]-nInits`, where each initialization set is accompanied
by a block of parameters that is preceded by the corresponding `-Init[$]` counter. In this example we have a single initialization set. The `Part-Species[$]-Init[$]-SpaceIC = cuboid` flag defines the type of the initialization set. Here, the particles are placed in a cuboid which is spanned by its base plane (`Part-Species[$]-Init[$]-BasePointIC`, `Part-Species[$]-Init[$]-BaseVector1IC`, `Part-Species[$]-Init[$]-BaseVector2IC`), a normal (`Part-Species[$]-Init[$]-NormalIC`) and its height (`Part-Species[$]-Init[$]-CuboidHeightIC`). Each type of the initialization set might have a different set of parameters and an overview is given in Section
{ref}`sec:particle-initialization-and-emission`. Here, simulation particles are inserted at a translational temperature (`MWTemperatureIC`) of $\pu{10000 K}$ using a Maxwellian velocity distribution (`velocityDistribution`), a vibrational temperature (`TempVib`) of $\pu{5000 K}$, a rotational temperature (`TempRot`) of $\pu{7500 K}$ at a zero flow velocity, which is defined through a velocity vector (`VeloVecIC`, will be normalized internally) and its magnitude (`VeloIC`).

    Part-Species1-nInits = 1

    Part-Species1-Init1-SpaceIC              = cuboid
    Part-Species1-Init1-velocityDistribution = maxwell_lpn
    Part-Species1-Init1-MWTemperatureIC      = 10000
    Part-Species1-Init1-TempVib              = 5000
    Part-Species1-Init1-TempRot              = 7500
    Part-Species1-Init1-PartDensity          = 1E23
    Part-Species1-Init1-BasePointIC          = (/0.,0.,0./)
    Part-Species1-Init1-BaseVector1IC        = (/4.64E-6,0.,0./)
    Part-Species1-Init1-BaseVector2IC        = (/0.,4.64E-6,0./)
    Part-Species1-Init1-NormalIC             = (/0.,0.,1./)
    Part-Species1-Init1-CuboidHeightIC       = 4.64E-6
    Part-Species1-Init1-VeloIC               = 0
    Part-Species1-Init1-VeloVecIC            = (/0.,0.,1./)

To calculate the number of simulation particles defined by `Part-Species[$]-Init[$]-PartDensity`, the selected weighting factor $w_{\text{CO}_2}$ and the volume of the complete domain $V=(\pu{4.64e-6 m})^3$ are utilized. This value can be used to chose the maximum particle number per processor accordingly.

$$ N_{\text{CO}_2,\text{sim}} = \frac{n_{\text{CO}_2} V}{w_{\text{CO}_2}} $$

(sec:tutorial-dsmc-dsmc-setup)=
### DSMC setup

Finally, DSMC has to be enabled (`UseDSMC = T`) and the particle movement is disabled via `Particles-DSMCReservoirSim = T` to reduce the computational effort. Keep in mind that the latter needs a compiled 
version of piclas using `PICLAS_TIMEDISCMETHOD = RESERVOIR`. Besides these settings `Particles-DSMC-CollisMode` is an important parameter. If set to 1, only elastic collisions 
are performed, if set to 2 relaxation processes are allowed and if set to 3 chemistry is enabled. Additionally, constant values for the rotational (`Particles-DSMC-RotRelaxProb`) and vibrational (`Particles-DSMC-VibRelaxProb`) relaxation probabilities are defined.

    ! =============================================================================== !
    ! DSMC
    ! =============================================================================== !
    UseDSMC                     = T
    Particles-DSMCReservoirSim  = T
    Particles-DSMC-CollisMode   = 2
    Particles-DSMC-RotRelaxProb = 0.2
    Particles-DSMC-VibRelaxProb = 0.02 

Besides the data given in the **parameter.ini**, a proper DSMC simulation needs additional species information, which is defined in the **DSMC.ini**. The species numeration needs to be the same in both files.

    ! =============================================================================== !
    ! Species1, CO2
    ! =============================================================================== !
    Part-Species1-SpeciesName    = CO2
    Part-Species1-InteractionID  = 2
    Part-Species1-PolyatomicMol  = T
    Part-Species1-NumOfAtoms     = 3
    Part-Species1-LinearMolec    = T
    Part-Species1-Tref           = 273
    Part-Species1-dref           = 5.10E-10
    Part-Species1-omega          = 0.24
    Part-Species1-CharaTempVib1  = 959.66
    Part-Species1-CharaTempVib2  = 959.66
    Part-Species1-CharaTempVib3  = 1918.6
    Part-Species1-CharaTempVib4  = 3382
    Part-Species1-CharaTempRot   = 0.6
    Part-Species1-Ediss_eV       = 5.45

The first block from `Part-Species[$]-InteractionID` to `Part-Species[$]-LinearMolec` declares the structure of the species. Available species types set by `Part-Species[$]-InteractionID` are listed in Section {ref}`sec:DSMC-species`. The second one from `Part-Species[$]-Tref` to `Part-Species[$]-omega` are the basis of the collision model utilized, in this case the Variable Hard Sphere (VHS) model. It is important to keep in mind that the $\omega$ in this file differs from the $\omega$ used by Bird. $\omega = \omega_\text{bird1994} - 0.5$. The last block from `Part-Species[$]-CharaTempVib1` to `Part-Species[$]-CharaTempRot` defines the vibrational excitation modes and sets the characteristic rotational temperature (utilized for the partition function). Finally, `Part-Species1-Ediss_eV` defines the dissociation energy of the molecule in [eV].

### Running the code

The command

    ./piclas parameter.ini DSMC.ini | tee std.out

executes the code and dumps all output into the file *std.out*.
If the run has completed successfully, which should take only a brief moment, the contents of the working folder should look like

    drwxrwxr-x 4,0K Okt 21 07:59 ./
    drwxrwxr-x 4,0K Dez  4 11:09 ./
    drwxrwxr-x 4,0K Dez  3 01:02 ../
    -rw-rw-r-- 1,4K Dez  4 11:05 DSMC.ini
    -rw-rw-r-- 8,7K Dez  4 11:09 dsmc_reservoir_chemisty_off_DSMCState_000.000005000.h5
    -rw-rw-r-- 1,8M Dez  4 11:09 dsmc_reservoir_chemisty_off_State_000.000000000.h5
    -rw-rw-r-- 1,8M Dez  4 11:09 dsmc_reservoir_chemisty_off_State_000.000005000.h5
    -rw-rw-r-- 6,9K Nov  1 05:05 dsmc_reservoir_mesh.h5
    -rw-rw-r--  931 Dez  4 11:09 ElemTimeStatistics.csv
    -rw-rw-r-- 7,5K Dez  4 10:49 parameter.ini
    -rw-rw-r-- 587K Dez  4 11:09 PartAnalyze.csv
    -rw-rw-r--  58K Dez  4 11:09 std.out
    -rw-rw-r--  15K Dez  4 11:09 userblock.tmp


Multiple additional files have been created, which are are named  **Projectname_State_Timestamp.h5**.
They contain the particle information (position, velocity, energy) and the solution vector of the equation system variables (if a field solver had been utilized) at each interpolation node at the given time, which corresponds
to multiples of **Analyze_dt**. If these files are not present, something went wrong during the execution of **piclas**.
In that case, check the `std.out` file for an error message and feel free to open an [issue](https://github.com/piclas-framework/piclas/issues).

After a successful completion, the last lines in this file should look as shown below:

    ------------------------------------------------------------------------
     Sys date  :    04.12.2021 11:05:43
     PID: CALCULATION TIME PER TSTEP/DOF: [ 1.77323E-04 sec ]
     EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]: [ 1.19627E-03 sec/h ]
     Timestep  :    2.0000000E-09
    #Timesteps :    2.5000000E+03
     WRITE STATE TO HDF5 FILE [dsmc_reservoir_chemisty_off_State_000.000005000.h5] ...DONE  [.005s]
    #Particles :    1.9979000E+04
    ------------------------------------------------------------------------
    ========================================================================
     PICLAS FINISHED! [            3.82 sec ] [     0:00:00:03]
    ========================================================================

## Visualization, chemistry disabled (post-processing)

In addition to `std.out`, which contains information about the simulation process, three other important file types were created. The simulation results are stored in `Projectname_State_Timestamp.h5`, `Projectname_DSMCState_Timestamp.h5` and `PartAnalyze.csv`. The latter stores the data generated during each the time step (or every `Part-AnalyzeStep`) and can be analyzed with e.g. [gnuplot](http://www.gnuplot.info/), any spreadsheet tool or [ParaView](https://www.paraview.org/). It allows the visualization of the transient development of the selected parameters (see Section {ref}`sec:tutorial-dsmc-analysis-setup`). As an example, the process of the examined relaxation, i.e. the convergence of translational, rotational and vibrational temperature, is shown in {numref}`fig:dsmc-reservoir-temperature-relaxation`.

```{figure} results/dsmc-reservoir-temperature-relaxation.svg
---
name: fig:dsmc-reservoir-temperature-relaxation
width: 50%
---

Temperature relaxation process towards thermal equilibrium.
```

## Simulation: Chemistry enabled

In the next step, we add more species and chemical reactions to the simulation, focusing on the different input parameters. Besides changes at simulation and output generation times the most important modification concerning the **parameter.ini** and **DSMC.ini** is related to chemistry. The definition of more than one species is needed to capture the products of the chemical reactions. The input of new species is similar to the first species as shown in in Section {ref}`sec:tutorial-dsmc-particle-solver` and {ref}`sec:tutorial-dsmc-dsmc-setup`. The values of the general physical properties are listed in {numref}`tab:dsmc_chem_on_phys`. As can be seen the number density was slightly reduced and the species are initialized in thermal equilibrium.

```{table} Physical properties at the simulation start
---
name: tab:dsmc_chem_on_phys
---
|                                         Property                                            |               Value                |
| ------------------------------------------------------------------------------------------- | :--------------------------------: |
| Species                                                                                     | $\text{CO}_2, \text{CO}, \text{O}$ |
| Molecule mass $m_{\text{CO}_2}$                                                             | $\pu{7.306E-26 kg}$                |
| Molecule mass $m_{\text{CO}}$                                                               | $\pu{4.65100E-26 kg}$              |
| Molecule mass $m_{\text{O}}$                                                                | $\pu{2.65700E-26 kg}$              |
| Number density $n_{\text{CO}_2}=n_{\text{CO}}=n_{\text{O}}$                                 | $\pu{1e22 m^{-3}}$                 |
| Translational temperature $T_{\text{trans, CO}_2}=T_{\text{trans, CO}}=T_{\text{trans, O}}$ | $\pu{10000 K}$                     |
| Rotational temperature $T_{\text{rot, CO}_2}=T_{\text{rot, CO}}$                            | $\pu{10000 K}$                     |
| Vibrational temperature $T_{\text{vib, CO}_2}=T_{\text{vib, CO}}$                           | $\pu{10000 K}$                     |
```

### Reactions

The parameter `Particles-DSMC-CollisMode = 3` mentioned in {ref}`sec:tutorial-dsmc-dsmc-setup` must be changed to enable the use of the chemistry module. `Particles-DSMC-BackwardReacRate = T` activates the calculation of the backward reaction to every considered reaction. The reaction rates are added to the output by setting `CalcReacRates`.

    Particles-DSMC-CollisMode        = 3
    Particles-DSMC-BackwardReacRate  = T
    CalcReacRates = T

While the electronic model is disabled by `Particles-DSMC-ElectronicModel = 0`, a `Particles-DSMCElectronicDatabase` has to be provided for the calculation of the partition functions for the equilibrium constant required for the calculation of the backward reaction rates. It contains the species-specific electronic energy levels, for more information see Secion {ref}`sec:DSMC-electronic-relaxation`.

    Particles-DSMC-ElectronicModel   = 0
    Particles-DSMCElectronicDatabase = DSMCSpecies_electronic_state_full_Data.h5

Copy one of the databases, which are used for the regression testing, to your current folder, e.g. from the path below

    cp $PICLAS_PATH/regressioncheck/WEK_Reservoir/CHEM_EQUI_diss_CH4/DSMCSpecies_electronic_state_full_Data.h5 ./

The definition of the reactions needs to be given in the **DSMC.ini** file. Additionally, the heat of formation of each species has to be provided for the calculation of the heat of reaction [K].

    Part-Species1-HeatOfFormation_K = -47328.35
    Part-Species2-HeatOfFormation_K = -13293.70
    Part-Species3-HeatOfFormation_K =  29969.45

As with the species, the number of reactions must first be defined by `DSMC-NumOfReactions`. Next, a model has to be chosen for each reaction via `DSMC-Reaction[$]-ReactionModel`. In this example the TCE model is used. Thus the connected parameters of the Arrhenius equation `DSMC-Reaction[$]-Arrhenius-Prefactor`, `DSMC-Reaction[$]-Arrhenius-Powerfactor` and `DSMC-Reaction[$]-Activation-Energy_K` must be set. The reactants of each reaction are defined as follows: `DSMC-Reaction[$]-Reactants` contains the species number of up to three reactants and `DSMC-Reaction[$]-Products` contains the species number of up to four products. Nonreactive parts are mentioned in another list, the length of which is given by `DSMC-Reaction[$]-NumberOfNonReactives`. The list `DSMC-Reaction[$]-NonReactiveSpecies` is filled with the species numbers. If the reaction rate depends on the non-reaction partner (e.g. higher dissociation rate if the non-reactive partner is an atom), each reaction can be defined separately by simply defining the second reactant and product, respectively.

    DSMC-NumOfReactions = 1
    DSMC-Reaction1-ReactionModel         = TCE
    DSMC-Reaction1-Reactants             = (/1,0,0/)
    DSMC-Reaction1-Products              = (/2,0,3,0/)
    DSMC-Reaction1-Arrhenius-Prefactor   = 1.15E-08
    DSMC-Reaction1-Arrhenius-Powerfactor = -1.5
    DSMC-Reaction1-Activation-Energy_K   = 63280
    DSMC-Reaction1-NumberOfNonReactives  = 3
    DSMC-Reaction1-NonReactiveSpecies    = (/1,2,3/)

Therefore, in this example with one reaction and each of the three species as possible non-reactive partner as well as the corresponding backward reaction lead to a number of six reactions in total. For more information see Section {ref}`sec:DSMC-chemistry`.

### Running the code (Parallel computing)

In order to investigate the transient behavior, a longer simulation time was chosen. This results in comparatively long computing times, which is why the use of several computing cores is recommended. The number of cores may not exceed the number of cells. This results in a maximum of 4 cores for the described simulation. Another important note is that bash does not understand aliases which are not at the start of a line. Thus a copy of the **piclas** binary must be located in the current folder

    cp $PICLAS_PATH/build/bin/piclas .
    
or the whole path to the binary must be used instead. Assuming a run with 4 cores is desired and the **piclas** binary is located at the current directory, the command

    mpirun -np 4 piclas parameter.ini DSMC.ini | tee std.out

executes the code and dumps all output into the file *std.out*.

## Visualization, chemistry enabled (post-processing)

To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**, **VisIt** or any other
visualisation tool for which the program **piclas2vtk** is used.

The parameters for **piclas2vtk** are stored in the **parameter.ini** file under

    ! =============================================================================== !
    ! piclas2vtk
    ! =============================================================================== !
    NVisu         = 1
    VisuParticles = T

where `NVisu` is the polynomial visualization degree on which the field solution is interpolated. The flag `VisuParticles` activates
the output of particle position, velocity and species, which is disabled per default due to the usually large number of particles.

Run the command

    ./piclas2vtk parameter.ini dsmc_reservoir_chemisty_on_State_000.00000*

to convert the HDF5 file to the binary VTK format (`*.vtu`), which can then be opened with e.g. ParaView.

The `Projectname_visuPart_Timestamp.vtu` files contain simulation particle specific data like their position, temperature and species. The figure below illustrates the initiated species and the resulting change of species. Additionaly, the particle information can be used for e.g. the determination of the velocity and energy distribution function of each species. It should be mentioned that the particle movement was disabled and therefore each reaction product stays at the position of the given reactant.

```{figure} results/dsmc-reservoir-species.jpg
---
name: fig:dsmc-reservoir-species
width: 100%
---

Comparison of the present species (simulation particles) at start (left) and end (right) time of the simulation.
```

While the figure above is only capable of giving a general overview about the processes in the reservoir, an analysis of the `PartAnalyze.csv` shows, that nearly all carbon dioxide is dissociated at the end of the simulation. During this process the overall temperature and the reaction probabilities are decreasing.

```{figure} results/dsmc-reservoir-reaction.jpg
---
name: fig:dsmc-reservoir-reaction
width: 100%
---

Development of species composition (left) and translational temperature and dissociation rate (right) over time.
```
