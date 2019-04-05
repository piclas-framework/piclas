\hypertarget{features_models}{}

# Features \& Models \label{chap:features_models}

The goal of PICLas is to enable to approximation of the complete Boltzmann equation:

$$ \frac{\partial f}{\partial t} + \mathbf{v}\cdot\frac{\partial f}{\partial \mathbf{x}} + \frac{\mathbf{F}}{m}\cdot\frac{\partial f}{\partial \mathbf{v}} = \left.\frac{\partial f}{\partial t}\right|_\mathrm{coll} $$

## Particle Tracking

### Linear

### Curved

## Boundary Conditions

### Field

#### Dielectric

### Particle

#### Specular/Reflective Wall

#### Porous Wall

The porous boundary condition uses a removal probability to determine whether a particle is deleted or reflected at the boundary.

## Particle Emission

### Surface Flux

#### Adaptive Boundaries

Multiple adaptive particle emission conditions can be defined.

## Particle in Cell

## Direct Simulation Monte Carlo

### Species Definition

### Relaxation

### Chemistry & Ionization

### Surface Chemistry

## Modelling of Continuum Gas Flows

Two methods are currently implemented to allow the simulation of gas flows in the continuum and transitional regime, where the DSMC method is computationally too expensive. The Fokker–Planck- and Bhatnagar-Gross-Krook-based approximation of the collision integral are compared in detail in paper to be published in Physics of Fluids.

### Fokker–Planck Collision Operator

The implementation of the FP-based collision operator is based on the publications by @Gorji2014 and @Pfeiffer2017. The collision integral is hereby approximated by a drift and diffusion process

$$  \left.\frac{\partial f}{\partial t}\right|_\mathrm{coll}\approx-\sum_{i=1}^3 {\frac{\partial }{\partial v_i}(A_i f)+\frac{1}{2}\sum_{i=1}^3 \sum_{j=1}^3\frac{\partial ^2 }{\partial v_i\partial v_j}(D_{ij}f)}, $$

where $\mathbf{A}$ is the drift vector and $\mathcal{D}$ the diffusion matrix.

The current implementation supports:

- 2 different methods: Cubic (only atomic species) and Ellipsoidal Statistical (ES)
- Single species, monoatomic and diatomic gases
- Thermal non-equilibrium with rotational and vibrational excitation (continuous or quantized treatment)

Relevant publications of the developers:

- Implementation of the cubic Fokker-Planck in PICLas (@Pfeiffer2017)
- Comparison of the cubic and ellipsoidal statistical Fokker-Planck (@Jun2019)

To enable the simulation with the FP module, the respective compiler setting has to be activated:

    PICLAS_TIMEDISCMETHOD = FP-Flow

A parameter file and species initialization file is required, analagous to the DSMC setup. To enable the simulation with the FP methods, select the Fokker-Planck method, cubic (`=1`) and ES (`=2`):

    Particles-FP-CollModel = 2

The **recommended method is ESFP**. The vibrational excitation can be controlled with the following flags, including the choice between continuous and quantized vibrational energy:

    Particles-FP-DoVibRelaxation = T
    Particles-FP-UseQuantVibEn   = T
    
An octree cell refinement until the given number of particles is reached can be utilized, which corresponds to an equal refinement in all three directions (x,y,z):

    Particles-FP-DoCellAdaptation = T
    Particles-FP-MinPartsPerCell  = 10

A coupled FP-DSMC simulation can be enabled, where the FP method will be utilized if the number density $[\text{m}^{-3}]$ is above a certain value:

    Particles-CoupledFPDSMC       = T
    Particles-FP-DSMC-SwitchDens  = 1E22

The flag `Particles-DSMC-CalcQualityFactors` controls the output of quality factors such as mean/maximal relaxation factor (mean: average over a cell, max: maximal value within the octree), max rotational relaxation factor, which are defined as

$$ \frac{\Delta t}{\tau} < 1,$$

where $\Delta t$ is the chosen time step and $1/\tau$ the relaxation frequency. The time step should be chosen as such that the relaxation factors are below unity. The `FP_DSMC_Ratio` gives the percentage of the sampled time during which the FP model was utilized. In a couple FP-DSMC simulation this variable indicates the boundary between FP and DSMC. However, a value below 1 can occur for pure FP simulations due to low particle numbers, when an element is skipped. Additionally, the Prandtl number utilized by the ESFP model is given.

### Bhatnagar-Gross-Krook Collision Operator

The implementation of the BGK-based collision operator is based on the publications by @Pfeiffer2018a and @Pfeiffer2018b. It allows the simulation of gas flows in the continuum and transitional regime, where the DSMC method is computationally too expensive. The collision integral is hereby approximated by a relaxation process:

$$ \left.\frac{\partial f}{\partial t}\right|_\mathrm{coll} \approx \nu(f^t-f), $$

where $f^t$ is the target distribution function and $\nu$ the relaxation frequency.

The current implementation supports:

- 4 different methods (i.e. different target distribution functions): Ellipsoidal Statistical, Shakov, standard BGK, and Unified
- Single species, monoatomic and diatomic gases
- Thermal non-equilibrium with rotational and vibrational excitation (continuous or quantized treatment)

Relevant publications of the developers:

- Implementation and comparison of the ESBGK, SBGK, and Unified models in PICLas for atomic species @Pfeiffer2018a
- Extension of the modelling to diatomic species including quantized vibrational energy treatment, validation of ESBGK with the Mach 20 hypersonic flow measurements of the heat flux on a $70^\circ$ cone @Pfeiffer2018b
- Simulation of a nozzle expansion (including the pressure chamber) with ESBGK, SBGK and coupled ESBGK-DSMC, comparison to experimental measurements (*to be published*)
- Simulation of the carbon dioxide hypersonic flow around a flat-faced cylinder, comparison of ESBGK, SBGK and DSMC regarding the shock structure and heat flux  (*to be published*)

To enable the simulation with the BGK module, the respective compiler setting has to be activated:

    PICLAS_TIMEDISCMETHOD = BGK

A parameter file and species initialization file is required, analagous to the DSMC setup. To enable the simulation with the BGK methods, select the BGK method, ES (`=1`), Shakov (`=2`), Standard BGK (`=3`), and Unified (`=4`):

    Particles-BGK-CollModel = 1

The **recommended method is ESBGK**. The vibrational excitation can be controlled with the following flags, including the choice between continuous and quantized vibrational energy:

    Particles-BGK-DoVibRelaxation = T
    Particles-BGK-UseQuantVibEn   = T
    
An octree cell refinement until the given number of particles is reached can be utilized, which corresponds to an equal refinement in all three directions (x,y,z):

    Particles-BGK-DoCellAdaptation = T
    Particles-BGK-MinPartsPerCell  = 10

It is recommended to utilize at least between 7 and 10 particles per (sub)cell. To enable the cell refinement above certain number density, the following option can be utilized

    Particles-BGK-SplittingDens = 1E23

A coupled BGK-DSMC simulation can be enabled, where the BGK method will be utilized if the number density $[\text{m}^{-3}]$ is above a certain value:

    Particles-CoupledBGKDSMC       = T
    Particles-BGK-DSMC-SwitchDens  = 1E22

The flag `Particles-DSMC-CalcQualityFactors` controls the output of quality factors such as mean/maximal relaxation factor (mean: average over a cell, max: maximal value within the octree), max rotational relaxation factor, which are defined as

$$ \frac{\Delta t}{\tau} < 1,$$

where $\Delta t$ is the chosen time step and $1/\tau$ the relaxation frequency. The time step should be chosen as such that the relaxation factors are below unity. The `BGK_DSMC_Ratio` gives the percentage of the sampled time during which the BGK model was utilized. In a couple BGK-DSMC simulation this variable indicates the boundary between BGK and DSMC. However, a value below 1 can occur for pure BGK simulations due to low particle numbers, when an element is skipped.

An option is available to utilize a moving average for the variables used in the calculation of the relaxation frequency:

    Particles-BGK-MovingAverage = T

The purpose is to increase the sample size for steady gas flows. An extension of this feature to account for unsteady flows is to limit the moving average to the last $N$ number of iterations, where the first value of the array is deleted and the most current value is added at the end of the list:

    Particles-BGK-MovingAverageLength = 100

Although this feature was tested with a hypersonic flow around a $70^\circ$ blunted cone and a nozzle expansion, a clear advantage could not be observed, however, it might reduce the statistical noise for other application cases.