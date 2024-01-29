(sec:background-gas)=
# Background Gas

A constant or spatially varying background gas (single species or mixture) can be utilized to enable efficient particle collisions between the
background gas and other particle species (represented by actual simulation particles). The assumption is that the density of the
background gas $n_{\mathrm{gas}}$ is much greater than the density of the particle species, e.g. the charged species in a plasma,
$n_{\mathrm{charged}}$

$$ n_{\mathrm{gas}} >> n_{\mathrm{charged}}.$$

Under this assumption, collisions within the particle species can be neglected and collisions between the background gas and
particle species do not alter the background gas conditions. It can be activated by using the regular particle insertion parameters
(as defined in Section {ref}`sec:particle-insertion`, the `-Init1` construct can be neglected since no other
initialization regions are allowed if the species is already defined a background species) and by defining the `SpaceIC` as
`background` as well as the number density [1/m$^3$] as shown below

    Part-Species1-SpaceIC     = background
    Part-Species1-PartDensity = 1E+22

Other species parameters such as mass, charge, temperature and velocity distribution for the background are also defined by the
regular read-in parameters. A mixture as a background gas can be simulated by simply defining multiple background species. Every
time step particles are generated from the background gas (for a mixture, the species of the generated particle is chosen
based on the species composition) and paired with the particle species. Subsequently, the collision probabilities are calculated
using the conventional DSMC routines and the VHS cross-section model. Afterwards, the collision process is performed (if the
probability is greater than a random number) and it is tested whether additional energy exchange and chemical reactions occur.
While the VHS model is sufficient to model collisions between neutral species, it cannot reproduce the phenomena of a
neutral-electron interaction. For this purpose, the cross-section based collision probabilities should be utilized, which are
discussed in the following section.

## Distribution from DSMC result

A spatially varying background gas distribution may be used by running a stand-alone DSMC simulation beforehand and using a time-averaged
DSMC state file (*PROJECT_DSMCState_\*.h5*) as input for the actual simulation by setting

    Particles-BGGas-UseDistribution              = T
    Particles-MacroscopicRestart-Filename        = DSMCResult.h5
    Part-SpeciesX-InitX-BGG-Distribution-SpeciesIndex = 1

where the first parameter activates the background gas distribution and the second parameter supplies the relative path to the file
from which the background gas density, velocity and temperature field is read (cell-constant values).
The third parameter defines which species index within the DSMC file is to be used as it may contain multiple species.

## Regions

Another possibility to define a non-constant background gas is available through the definition of regions. Multiple regions defined
by simple geometrical volumes (e.g. cylinder) can be mapped to different species. First, one or more regions are defined:

    Particles-BGGas-nRegions                  = 1
    Particles-BGGas-Region1-Type              = cylinder
    Particles-BGGas-Region1-RadiusIC          = 5E-6
    Particles-BGGas-Region1-CylinderHeightIC  = 5E-6
    Particles-BGGas-Region1-BasePointIC       = (/0.,0.,0./)
    Particles-BGGas-Region1-BaseVector1IC     = (/1.,0.,0./)
    Particles-BGGas-Region1-BaseVector2IC     = (/0.,1.,0./)

Here, a cylinder is defined by two base vectors (from which a normal is determined for the direction of the cylinder height),
basepoint, radius and cylinder height. The definition of the species is the same as described above, with the addition of an
additional parameter, defining in which region, these properties should be applied to:

    Part-Species1-Init1-BGG-Region            = 1

While a species can be part of different regions through multiple inits and multiple species can be part of a single region,
overlapping regions are not allowed. Whether an element is within a region is determined through the midpoint of the element and
thus it does not have to be fully enveloped.

## Trace species

If the number densities of the background gas species differ greatly and a specific background species is of interest (or the interaction with it) that has a lower number density compared to the other background species, it can be defined as a so-called trace species as shown below.

    Part-vMPF                        = T
    Part-Species1-Init1-TraceSpecies = T

The first flag enables the variables weighting factor feature in general (details about this feature can be found in Section {ref}`sec:variable-particle-weighting`). An additional flag defines the background gas species as a trace species, where multiple trace species can be defined. Finally, the weighting factors of the background species can be adopted to define the difference in the weighting factors.

    Part-Species1-MacroParticleFactor = 1E2
    Part-Species2-MacroParticleFactor = 1E4

Using the values above, each collision with the first background species will result in 100 collision tests using the simulation particle (ie. not the background species) and randomly generated background particles. Consequently, the number of samples for the trace species will be increased and simulation particles with the weighting factor of the trace background species will be introduced into the simulation.

## Cross-section based collision probability

For modelling of particle collisions with the Particle-in-Cell method, often the Monte Carlo Collision (MCC) algorithm is utilized.
Here, experimentally measured or ab-initio calculated cross-sections are typically utilized to determine the collision probability,
based on the cross-section [m$^2$], the timestep [s], the particle velocity [m/s] and the target gas number density [m$^{-3}$]

$$ P = 1 - e^{-\sigma v \Delta t n_{\mathrm{gas}}}.$$

In PICLas, the null collision method after {cite}`Birdsall1991`,{cite}`Vahedi1995` is available, where the number of collision
pairs is determined based a maximum collision frequency. Thus, the computational effort is reduced as not every particle has to be
checked for a collision, such as in the previously described DSMC-based background gas. To activate the MCC procedure, the
collision cross-sections have to be supplied via read-in from a database

    Particles-CollXSec-Database = MCC_Database.h5
    Particles-CollXSec-NullCollision = TRUE

Cross-section data can be retrieved from the [LXCat database](https://fr.lxcat.net/home/) {cite}`Pitchford2017` and converted with
a Python script provided in the tools folder: `piclas/tools/crosssection_database`. Details on how to create an own database with
custom cross-section data is given in Section {ref}`sec:tools-xsec-collision`. Finally, the input which species should be treated with the MCC
model is required

    Part-Species2-SpeciesName = electron
    Part-Species2-UseCollXSec = T

The read-in of the cross-section data is based on the provided species name and the species name of the background gas (e.g. if the
background species name is Ar, the code will look for a container named `Ar-electron` in the MCC database). Finally, the
cross-section based collision modelling (e.g. for neutral-charged collisions) and the VHS model (e.g. for neutral-neutral
collisions) can be utilized within a simulation for different species.

## Cross-section based vibrational relaxation probability

In the following, the utilization of cross-section data is extended to the determination of the vibrational relaxation probability.
When data is available, it will be read-in by the Python script described above. If different vibrational levels are available,
they will be summarized to a single relaxation probability. Afterwards the regular DSMC-based relaxation procedure will be
performed. To enable the utilization of these levels, the following flag shall be supplied

    Part-Species2-UseVibXSec = T

It should be noted that even if Species2 corresponds to an electron, the vibrational cross-section data will be read-in for any
molecule-electron pair. If both species are molecular, priority will be given to the species utilizing this flag.

(sec:background-gas-electronic-xsec)=
## Cross-section based electronic relaxation probability

In the following, the utilization of cross-section data is extended to the electronic excitation for neutral-electron collisions. When data is available, it will be read-in by the Python script described above. Each level will be handled separately, allowing the atom/molecule to be excited in each level. The cross-section data will be used to determine whether and which excitation will occur. During the excitation procedure the energy of the atom/molecule will be set to respective level. To enable this model, the following flags are required

    Particles-DSMC-ElectronicModel    = 3
    Part-Species1-UseElecXSec         = T

The species-specific flag `UseElecXSec` should be set to `TRUE` for the heavy-species and not the electron species. Currently, this model is only implemented for the background gas, however, an extension to regular DSMC simulations is envisioned. It can be used with cross-section based as well as Variable Hard Sphere (VHS) collision modelling. The output of the relaxation rates is provided through

    CalcRelaxProb = T

However, the electronic temperature is currently not added to the total temperature for the PartAnalyze.csv output. Additionally, the cell-local excitation rate [1/s] per electronic level of the respective species can be sampled and output as described in Section {ref}`sec:sampling-elec-excitation`.
