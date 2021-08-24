# Direct Simulation Monte Carlo

To enable the simulation with DSMC, an appropriate time discretization method including the DSMC module should be chosen before the
code compilation. A stand-alone DSMC simulation can be enabled by compiling PICLas with the following parameter

    PICLAS_TIMEDISCMETHOD = DSMC

The DSMC method can then be enabled in the parameter file by

    UseDSMC = T

Additionally, the number of simulated physical models depending on the application can be controlled through

    Particles-DSMC-CollisMode = 1   ! Elastic collisions only
                                2   ! Internal energy exchange
                                3   ! Chemical reactions

`CollisMode = 1` can be utilized for the simulation of a non-reactive, cold atomic gas, where no chemical reactions or electronic
excitation is expected. `CollisMode = 2` should be chosen for non-reactive diatomic gas flows to include the internal energy
exchange (by default including the rotational and vibrational energy treatment). Finally, reactive gas flows can be simulated with
`CollisMode = 3`. The following sections describe the required definition of species parameter (Section {ref}`sec:DSMC-species`),
the parameters for the internal energy exchange (Section {ref}`sec:DSMC-relaxation`) and chemical reactions (Section
{ref}`sec:DSMC-chemistry`).

A fixed ("manual") simulation time step $\Delta t$ is defined by

    ManualTimeStep = 1.00E-7

(sec:DSMC-species)=
## Species Definition

For the DSMC simulation, additional species-specific parameters (collision model parameters, characteristic vibrational
temperature, etc.) are required. This file is also utilized for the definition of chemical reactions paths. To define a species,
its name as well as an `InteractionID` have to be defined

    Part-Species1-SpeciesName = CH4
    Part-Species1-InteractionID = 2

The name is at the moment only utilized to retrieve the electronic energy levels from an additional database. The interaction ID
determines the type of a species as follows

|   ID | Type                               |
| ---: | ---------------------------------- |
|    1 | Atom                               |
|    2 | Molecule (diatomic and polyatomic) |
|    4 | Electron                           |
|   10 | Atomic Ion                         |
|   20 | Molecular Ion                      |

Depending on the utilized collision model, different parameters have to be defined. As an example, the parameters for the Variable
Hard Sphere (VHS) collision cross-section model are be defined by the temperature exponent $\omega = [0,0.5]$, reference
temperature $T_{\mathrm{ref}}$ [K] and diameter $d_{\mathrm{ref}}$ [m]

    Part-Species1-omega = 0.24
    Part-Species1-Tref = 273
    Part-Species1-dref = 4.63E-10

More detail on the utilization of species-specific, collision-specific parameters and the utilization of the Variable Soft Sphere
(VSS) model are given in Section {ref}`sec:DSMC-collision`.

In case of chemical reactions, the input of a species-specific heat/enthalpy of formation [K] is required. A reliable source for
these reference values are the [Active Thermochemical Tables (ATcT)](https://atct.anl.gov/) by the Argonne National Laboratory
{cite}`Ruscic2004,Ruscic2005`.

    Part-Species1-HeatOfFormation_K = 0.0

In the case of ionization reactions the heat of formation is not required, as the ionization energy is read-in as the last
electronic level of the neutral (or previous) state. Therefore, an additional entry for each ionic species about its previous
state is required. This is done by providing the species index, an example is given below assuming that the first species is C,
whereas the second and thirds species are the first and second ionization levels, respectively.

    Part-Species2-PreviousState = 1
    Part-Species3-PreviousState = 2

If the previous state of the ionized species is not part of the simulation (e.g. a reservoir with N, N$_2^+$, and e), the previous
state of N$_2^+$ can be set to zero, however, in that case a heat/enthalpy of formation has to be provided, which includes the
ionization energy (as given in ATcT).

Diatomic molecular species require the definition of the characteristic temperature [K] and their dissociation energy [eV]
(which is at the moment utilized as a first guess for the upper bound of the temperature calculation as well as the threshold
energy for the dissociation by the Quantum-Kinetic chemistry model)

    Part-Species1-CharaTempVib = 4194.9
    Part-Species1-Ediss_eV = 4.53

Polyatomic molecular species require an additional flag, the input of the number of atoms  and whether the molecule is linear
(e.g. CO$_2$, $\xi_{\mathrm{rot}} = 2$) or non-linear (e.g. H$_2$O, CH$_4$, $\xi_{\mathrm{rot}} = 3$). The number of the
vibrational degrees of freedom is then given by

$$ \alpha = 3 N_{\mathrm{atom}} - 3 - \xi_{\mathrm{rot}} $$

As an example the parameters of CH$_3$ are given below. The molecule has four vibrational modes, with two of them having a
degeneracy of two. These values are simply given the according amount of times

    Part-Species1-NumOfAtoms = 4
    Part-Species1-LinearMolec = FALSE
    Part-Species1-CharaTempVib1 = 4320.6
    Part-Species1-CharaTempVib2 = 872.1
    Part-Species1-CharaTempVib3 = 4545.5
    Part-Species1-CharaTempVib4 = 4545.5
    Part-Species1-CharaTempVib5 = 2016.2
    Part-Species1-CharaTempVib6 = 2016.2

These parameters allow the simulation of non-reactive gases. Additional parameters required for the consideration of chemical
reaction are given in Section {ref}`sec:DSMC-chemistry`.

(sec:DSMC-collision)=
## Pairing & Collision Modelling

By default, a conventional statistical pairing algorithm randomly pairs particles within a cell. Here, the mesh should resolve
the mean free path to avoid numerical diffusion. To circumvent this requirement, an octree-based sorting and cell refinement
{cite}`Pfeiffer2013` can be enabled by

    Particles-DSMC-UseOctree        = T
    Particles-OctreePartNumNode     = 80        ! (3D default, 2D default: 40)
    Particles-OctreePartNumNodeMin  = 50        ! (3D default, 2D default: 28)

The algorithm refines a cell recursively as long as the mean free path is smaller than a characteristic length (approximated by
the cubic root of the cell volume) and the number of particles is greater than `Particles-OctreePartNumNode`. The latter condition
serves the purpose to accelerate the simulation by avoiding looking for the nearest neighbour in a cell with a large number of
particles. To avoid cells with a low particle number, the cell refinement is stopped when the particle number is below
`Particles-OctreePartNumNodeMin`. These two parameters have different default values for 2D/axisymmetric and 3D simulations.

To further reduce numerical diffusion, the nearest neighbour search for the particle pairing can be enabled

    Particles-DSMC-UseNearestNeighbour = T

An additional attempt to increase the quality of simulations results is to prohibit repeated collisions between particles
{cite}`Shevyrin2005,Akhlaghi2018`. This options is enabled by default in 2D/axisymmetric simulations, but disabled by default in
3D simulations.

    Particles-DSMC-ProhibitDoubleCollision = T

The Variable Hard Sphere (VHS) is utilized by default with collision-averaged parameters, which are given per species

    Part-Species1-omega = 0.24
    Part-Species1-Tref = 273
    Part-Species1-dref = 4.63E-10

To enable the Variable Soft Sphere (VSS) model, the additional $\alpha$ parameter is required

    Part-Species1-alphaVSS = 1.2

In order to enable the collision-specific definition of the VHS/VSS parameters, a different input is required

    ! Input in parameter.ini
    Particles-DSMC-averagedCollisionParameters = F
    ! Input in species.ini
    Part-Collision1 - partnerSpecies = (/1,1/)              ! Collision1: Parameters for the collision between equal species
    Part-Collision1 - Tref           = 273
    Part-Collision1 - dref           = 4.037e-10
    Part-Collision1 - omega          = .216
    Part-Collision1 - alphaVSS       = 1.448
    Part-Collision2 - partnerSpecies = (/2,1/)              ! Collision2: Parameters for the collision between species 2 and 1

The numbers in the `partnerSpecies` definition correspond to the species numbers and their order is irrelevant. Collision-specific
parameters can be obtained from e.g. {cite}`Swaminathan-Gopalan2016`.

### Cross-section based collision probabilities

Cross-section data to model collisional and relaxation probabilities (e.g. in case of electron-neutral collisions), analogous to
Monte Carlo Collisions, can be utilized and is described in Section {ref}`sec:tools-xsec-collision` and {ref}`sec:xsec-chemistry`.

(sec:DSMC-relaxation)=
## Inelastic Collisions \& Relaxation

To consider inelastic collisions and relaxation processes within PICLas, the chosen `CollisMode` has to be at least 2

    Particles-DSMC-CollisMode = 2

Two selection procedures are implemented, which differ whether only a single or as many as possible relaxation processes can occur
for a collision pair. The default model (`SelectionProcedure = 1`) allows the latter, so-called multi-relaxation method, whereas
`SelectionProcedure = 2` enables the prohibiting double-relaxation method {cite}`Haas1994b`

    Particles-DSMC-SelectionProcedure = 1    ! Multi-relaxation
                                        2    ! Prohibiting double-relaxation

Rotational, vibrational and electronic relaxation (not included by default, see Section {ref}`sec:DSMC-electronic-relaxation`
for details) processes are implemented in PICLas and their specific options to use either constant relaxation probabilities
(default) or variable, mostly temperature dependent, relaxation probabilities are discussed in the following sections.
To achieve consistency between continuum and particle-based relaxation modelling, the correction factor of Lumpkin
{cite}`Lumpkin1991` can be enabled (default = F):

    Particles-DSMC-useRelaxProbCorrFactor = T

### Rotational Relaxation

To adjust the rotational relaxation this variable has to be changed:

    Particles-DSMC-RotRelaxProb = 0.2   ! Value between 0 and 1 as a constant probability
                                    2   ! Model by Boyd
                                    3   ! Model by Zhang

If `RotRelaxProb` is between 0 and 1, it is set as a constant rotational relaxation probability (default = 0.2). `RotRelaxProb = 2`
activates the variable rotational relaxation model according to Boyd {cite}`Boyd1990a`. Consequently, for each molecular species
two additional parameters have to be defined, the rotational collision number and the rotational reference temperature. As an
example, nitrogen is used {cite}`Boyd1990b`.

    Part-Species1-CollNumRotInf = 23.3
    Part-Species1-TempRefRot = 91.5

It is not recommended to use this model with the prohibiting double-relaxation selection procedure (`Particles-DSMC-SelectionProcedure = 2`). Low collision energies result in high relaxation probabilities, which can lead to cumulative collision probabilities greater than 1.

If the relaxation probability is equal to 3, the relaxation model of Zhang et al. {cite}`Zhang2012` is used. However, it is only
implemented for nitrogen and not tested. It is not recommended for use.

### Vibrational Relaxation

Analogous to the rotational relaxation probability, the vibrational relaxation probability is implemented. This variable has to be
changed, if the vibrational relaxation probability should be adjusted:

    Particles-DSMC-VibRelaxProb = 0.004 ! Value between 0 and 1 as a constant probability
                                      2 ! Model by Boyd

If `VibRelaxProb` is between 0 and 1, it is used as a constant vibrational relaxation probability (default = 0.004). The variable
vibrational relaxation model of Boyd {cite}`Boyd1990b` can be activated with `VibRelaxProb = 2`. For each molecular species pair,
the constants A and B according to Millikan and White {cite}`MillikanWhite1963` (which will be used for the calculation of the
characteristic velocity and vibrational collision number according to Abe {cite}`Abe1994`) and the vibrational cross section have
to be defined. The given example below is a 2 species mixture of nitrogen and oxygen, using the values for A and B given by Farbar
{cite}`Farbar2010` and the vibrational cross section given by Boyd {cite}`Boyd1990b`:

    Part-Species1-MWConstA-1-1 = 220.00
    Part-Species1-MWConstA-1-2 = 115.10
    Part-Species1-MWConstB-1-1 = -12.27
    Part-Species1-MWConstB-1-2 = -6.92
    Part-Species1-VibCrossSection = 1e-19

    Part-Species2-MWConstA-2-2 = 129.00
    Part-Species2-MWConstA-2-1 = 115.10
    Part-Species2-MWConstB-2-2 = -9.76
    Part-Species2-MWConstB-2-1 = -6.92
    Part-Species2-VibCrossSection = 1e-19

It is not possible to calculate an instantaneous vibrational relaxation probability with this model {cite}`Boyd1992`. Thus, the
probability is calculated for every collision and is averaged. To avoid large errors in cells containing only a few particles,
a relaxation of this average probability is implemented. The relaxation factor $\alpha$ can be changed with the following parameter
in the ini file:

    Particles-DSMC-alpha = 0.99

The new probability is calculated with the vibrational relaxation probability of the $n^{\mathrm{th}}$ iteration $P^{n}_{\mathrm{v}}$,
the number of collision pairs $n_{\mathrm{pair}}$ and the average vibrational relaxation probability of the actual iteration
$P^{\mathrm{iter}}_{\mathrm{v}}$.

$$P^{n+1}_{\mathrm{v}}= P^{n}_{\mathrm{v}}  \cdot  \alpha^{2  \cdot  n_{\mathrm{pair}}} + (1-\alpha^{2  \cdot  n_{\mathrm{pair}}}) \cdot P^{\mathrm{iter}}_{\mathrm{v}} $$

This model is extended to more species by calculating a separate probability for each species. An initial vibrational relaxation
probability is set by calculating $\mathrm{INT}(1/(1-\alpha))$ vibrational relaxation probabilities for each species and cell by
using an instantaneous translational cell temperature.

(sec:DSMC-electronic-relaxation)=
### Electronic Relaxation

For the modelling of electronic relaxation, two models are available: the model by Liechty et al. {cite}`Liechty2011a`, where each
particle has a specific electronic state and the model by Burt and Eswar {cite}`Burt2015b`, where each particle has an electronic
distribution function attached. Both models utilize tabulated energy levels, which can be found in literature for a wide range of
species (e.g. for monatomic {cite}`NISTASD`, diatomic {cite}`Huber1979`, polyatomic {cite}`Herzberg1966` molecules). An example
database `DSMCSpecies_electronic_state_full_Data.h5` can be found in e.g.
`piclas/regressioncheck/checks/NIG_Reservoir/CHEM_EQUI_TCE_Air_5Spec`, where the energy levels are stored in containers and
accessed via the species name, e.g. `Part-Species1-SpeciesName=N2`. Each level is described by its degeneracy in the first column
and by the energy in [J] in the second column. To include electronic excitation in the simulation, the following parameters are
required

    Particles-DSMC-ElectronicModel  = 0     ! No electronic energy is considered (default)
                                    = 1     ! Model by Liechty
                                    = 2     ! Model by Burt
    Particles-DSMCElectronicDatabase = DSMCSpecies_electronic_state_full_Data.h5

In case of a large number of electronic levels, their number can be reduced by providing a relative merge tolerance.
Levels those relative differences are below this parameter will be merged:

    EpsMergeElectronicState = 1E-3

However, this option should be evaluated carefully based on the specific simulation case and tested against a zero/very
low merge tolerance. Finally, the default relaxation probability can be adjusted by

    Particles-DSMC-ElecRelaxProb = 0.01

An electronic state database can be created using a Fortran tool in `piclas/tools/electronic_data`. An alternative is to use the
Python-based script discussed in Section {ref}`sec:tools-xsec-collision` and to adapt it to electronic energy levels.

(sec:DSMC-chemistry)=
## Chemistry & Ionization

Three chemistry models are currently available in PICLas

  * Quantum Kinetic (QK)
  * Total Collision Energy (TCE)
  * Cross-section (XSec)

The model of each reaction can be chosen separately. If a collision pair has multiple reaction paths (e.g. CH3 + H, two possible
dissociation reactions and a recombination), the reaction paths of QK and TCE are treated together, meaning that it is decided
between those reaction paths. If a reaction path is selected, the reaction is performed and the following routines of the chemistry
module are not performed. It is recommended not to combine cross-section based reaction paths with other reaction paths using
QK/TCE for the same collision pair.

Chemical reactions and ionization processes require

    Particles-DSMC-CollisMode = 3

The reactions paths can then be defined in the species parameter file. First, the number of reactions to read-in has to be defined

    DSMC-NumOfReactions = 2

A reaction is then defined by

    DSMC-Reaction1-ReactionModel = TCE
    DSMC-Reaction1-Reactants     = (/1,1,0/)
    DSMC-Reaction1-Products      = (/2,1,2,0/)

where the reaction model can be defined as follows

| Model | Description                                       |
| ----: | ------------------------------------------------- |
|   TCE | Total Collision Energy: Arrhenius-based chemistry |
|    QK | Quantum Kinetic: Threshold-based chemistry        |
|  XSec | Cross-section based chemistry                     |
| phIon | Photo-ionization (e.g. N + ph -> N$^+$ + e)       |

The reactants (left-hand side) and products (right-hand side) are defined by their respective species index. The photo-ionization
reaction is a special case to model the ionization process within a defined volume by photon impact (see Section
{ref}`sec:particle-photo-ionization`). It should be noted that for the dissociation reaction, the first given species is the
molecule to be dissociated. The second given species is the non-reacting partner, which can either be defined specifically or set
to zero to define multiple possible collision partners. In the latter case, the number of non-reactive partners and their species
have to be given by

    DSMC-Reaction1-Reactants=(/1,0,0/)
    DSMC-Reaction1-Products=(/2,0,2,0/)
    DSMC-Reaction1-NumberOfNonReactives=3
    DSMC-Reaction1-NonReactiveSpecies=(/1,2,3/)

This allows to define a single reaction for an arbitrary number of collision partners. In the following, three possibilities to
model the reaction rates are presented.
### Total Collision Energy (TCE)

The Total Collision Energy (TCE) model {cite}`Bird1994` utilizes Arrhenius type reaction rates to reproduce the probabilities for
a chemical reaction. The extended Arrhenius equation is

$$k(T) = A T^b e^{-E_\mathrm{a}/T}$$

where $A$ is the prefactor ([1/s, m$^3$/s, m$^6$/s] depending on the reaction type), $b$ the power factor and $E_\mathrm{a}$ the
activation energy [K]. These parameters can be defined in PICLas as follows

    DSMC-Reaction1-Arrhenius-Prefactor=6.170E-9
    DSMC-Reaction1-Arrhenius-Powerfactor=-1.60
    DSMC-Reaction1-Activation-Energy_K=113200.0

An example initialization file for a TCE-based chemistry model can be found in the regression tests (e.g.
`regressioncheck/checks/NIG_Reservoir/CHEM_EQUI_TCE_Air_5Spec`).

### Quantum Kinetic Chemistry (QK)

The Quantum Kinetic (QK) model {cite}`Bird2011` chooses a different approach and models chemical reactions on the microscopic level.
Currently, the QK model is only available for ionization and dissociation reactions. It is possible to utilize TCE- and QK-based
reactions in the same simulation for different reactions paths for the same collision pair, such as the ionization and dissociation
reactions paths (e.g. N$_2$ + e can lead to a dissociation with the TCE model and to an ionization with the QK model).
An example setup can be found in the regression tests (e.g. `regressioncheck/checks/NIG_Reservoir/CHEM_QK_multi-ionization_C_to_C6+`).

Besides the reaction model, reactants and products definition no further parameter are required for the reaction. However,
the dissociation energy [eV] has to be defined on a species basis

    Part-Species1-Ediss_eV = 4.53

The ionization threshold is determined from the last level of the previous state of the ionized product and thus requires the
read-in of the electronic state database.
(sec:xsec-chemistry)=
### Cross-section Chemistry (XSec)

The cross-section based chemistry model utilizes experimentally measured or ab-initio calculated cross-sections (analogous to
the collision probability procedure described in Section {ref}`sec:tools-xsec-collision`). It requires the same database, where the
reaction paths are stored per particle pair, e.g. the `N2-electron` container contains the `REACTION` folder, which includes the
reactions named by their products, e.g. `N2Ion1-electron-electron`.

If the defined reaction cannot be found in the database, the code will abort. It should be noted that this model is not limited to
the utilization with MCC or a background gas and can be used with conventional DSMC as an alternative chemistry model. Here, the
probability will be added to the collision probability to reproduce the reaction rate. Examples of the utilization of this model
can be found in the regression tests (e.g. `regressioncheck/checks/NIG_Reservoir/CHEM_RATES_XSec_Chem_H2-e`).
### Backward Reaction Rates

Backward reaction rates can be calculated for any given forward reaction rate by using the equilibrium constant

$$K_\mathrm{equi} = \frac{k_\mathrm{f}}{k_\mathrm{b}}$$

where $K_\mathrm{equi}$ is calculated through partition functions. This option can be enabled for all given reaction paths in the
first parameter file by

    Particles-DSMC-BackwardReacRate = T

or it can be enabled for each reaction path separately in the species parameter file, e.g. to disable certain reaction paths or to
treat the backward reaction directly with a given Arrhenius rate

    DSMC-Reaction1-BackwardReac = T

It should be noted that if the backward reactions are enabled globally, a single backward reaction can be disabled by setting the
reaction-specific flag to false and vice versa, if the global backward rates are disabled then a single backward reaction can be
enabled. Since the partition functions are tabulated, a maximum temperature and the interval are required to define the temperature
range which is expected during the simulation.

    Particles-DSMC-PartitionMaxTemp = 120000.
    Particles-DSMC-PartitionInterval = 20.

 Should a collision temperature be outside of that range, the partition function will be calculated on the fly. Additional
 species-specific parameters are required in the species initialization file to calculate the rotational partition functions

    Part-Species1-SymmetryFactor=2
    ! Linear poly- and diatomic molecules
    Part-Species1-CharaTempRot=2.1
    ! Non-linear polyatomic molecules
    Part-Species1-CharaTempRot1=2.1
    Part-Species1-CharaTempRot2=2.1
    Part-Species1-CharaTempRot3=2.1

The rotational symmetry factor depends on the symmetry point group of the molecule and can be found in e.g. Table 2 in
{cite}`Fernandez-Ramos2007`. While linear polyatomic and diatomic molecules require a single characteristic rotational
temperature, three values have to be supplied for non-linear polyatomic molecules. Finally, electronic energy levels have to be
supplied to consider the electronic partition function. For this purpose, the user should provide an electronic state database as
presented in Section {ref}`sec:DSMC-electronic-relaxation`.

## Additional Features
### Deletion of Chemistry Products

Specified product species can be deleted immediately after the reaction occurs, e.g. if an ionization process with a background gas
is simulated and the neutral species as a result from dissociation are not of interest. To do so, the number of species to be
deleted and there indices have to be defined

    Particles-Chemistry-NumDeleteProducts = 2
    Particles-Chemistry-DeleteProductsList = (/2,3/)

### Ambipolar Diffusion

A simple ambipolar diffusion model in order to be able to ignore the self-induced electric fields, e.g. for the application in
hypersonic re-entry flows, where ionization reactions cannot be neglected, can be enabled by

    Particles-DSMC-AmbipolarDiffusion = T

Electrons are now attached to and move with the ions, although, they still have their own velocity vector and are part of the
pairing and collisional process (including chemical reactions). The velocity vector of the ion species is not additionally
corrected to account for the acceleration due to the self-induced electric fields. The restart from a state file without previously
enabled ambipolar diffusion is currently not supported. However, the simulation can be continued if a macroscopic output is
available with the macroscopic restart. In that case, the electrons are not inserted but paired with an ion and given the sampled
velocity from the macroscopic restart file.

(sec:DSMC-quality)=
## Ensuring Physical Simulation Results

To determine whether the DSMC related parameters are chosen correctly, so-called quality factors can be written out as part of the
regular DSMC state file output by

    Particles-DSMC-CalcQualityFactors = T

This flag writes out the spatial distribution of the mean and maximal collision probability (`DSMC_MeanCollProb` and
`DSMC_MaxCollProb`). On the one hand, maximal collision probabilities above unity indicate that the time step should be reduced.
On the other hand, very small collision probabilities mean that the time step can be further increased. Additionally, the ratio of
the mean collision separation distance to the mean free path is written out (`DSMC_MCSoverMFP`)

$$\frac{l_{\mathrm{mcs}}}{\lambda} < 1$$

The mean collision separation distance is determined during every collision and compared to the mean free path, where its ratio
should be less than unity. Values above unity indicate an insufficient particle discretization. In order to estimate the required
weighting factor $w$, the following equation can be utilized for a 3D simulation

$$w < \frac{1}{\left(\sqrt{2}\pi d_{\mathrm{ref}}^2 n^{2/3}\right)^3},$$

where $d_{\mathrm{ref}}$ is the reference diameter and $n$ the number density. Here, the largest number density within the
simulation domain should be used as the worst-case. For supersonic/hypersonic flows, the conditions behind a normal shock can be
utilized as a first guess. For a thruster/nozzle expansion simulation, the chamber or throat conditions are the limiting factor.

