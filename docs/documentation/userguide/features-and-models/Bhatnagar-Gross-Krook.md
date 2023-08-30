(sec:BGK-Flow)=
# Bhatnagar-Gross-Krook Collision Operator

The implementation of the BGK-based collision operator is based on the publications by {cite}`Pfeiffer2018a` and {cite}`Pfeiffer2018b`.
It allows the simulation of gas flows in the continuum and transitional regimes, where the DSMC method is computationally too expensive.
The collision integral is hereby approximated by a relaxation process:

$$ \left.\frac{\partial f}{\partial t}\right|_{\mathrm{coll}} \approx \nu(f^t-f), $$

where $f^t$ is the target distribution function and $\nu$ the relaxation frequency.

The current implementation supports:

- Three different BGK methods (i.e. different target distribution functions): Ellipsoidal Statistical, Shakov, and standard BGK
- Single species: monatomic, diatomic, and polyatomic gases
- Gas mixtures with an arbitrary number of monatomic, diatomic, and polyatomic species
- Thermal non-equilibrium with rotational and vibrational excitation (continuous or quantized treatment)
- 2D axisymmetric simulations
- Variable time step (adaption of the distribution according to the maximal relaxation factor and linear scaling)

Relevant publications of the developers:

- Implementation and comparison of the ESBGK, SBGK, and Unified models in PICLas for atomic species {cite}`Pfeiffer2018a`
- Extension of the modeling to diatomic species including quantized vibrational energy treatment, validation of ESBGK with the Mach 20 hypersonic flow measurements of the heat flux on a $70^\circ$ cone {cite}`Pfeiffer2018b`
- Simulation of a nozzle expansion (including the pressure chamber) with ESBGK, SBGK and coupled ESBGK-DSMC, comparison to experimental measurements {cite}`Pfeiffer2019a`,{cite}`Pfeiffer2019b`
- Extension to polyatomic molecules, simulation of the carbon dioxide hypersonic flow around a flat-faced cylinder, comparison of ESBGK, SBGK and DSMC regarding the shock structure and heat flux {cite}`Pfeiffer2019c`
- Implemention of Brull's multi-species modeling for monatomic gas mixtures using Wilke's mixture rules and collision integrals for the calculation of transport coefficients {cite}`Pfeiffer2021`
- Extension of the implementation of Brull's ESBGK multi-species model to diatomic gas mixtures using Wilke's mixture rules (under review)
- Extension of the ESBGK method to multi-species modeling of polyatomic molecules, based on the ESBGK models of Mathiaud, Mieussens, Pfeiffer, and Brull, and including internal energies with multiple vibrational degrees of freedom, using Wilke's mixture rules and collision integrals in comparison (under review)

To enable the simulation with the BGK module, the respective compiler setting has to be activated:

    PICLAS_TIMEDISCMETHOD = BGK-Flow

A parameter file and species initialization file is required, analogous to the DSMC setup. It is recommended to utilize a previous DSMC parameter file to ensure a complete simulation setup. To enable the simulation with the BGK methods, select the BGK method, ES (`= 1`), Shakov (`= 2`), and standard BGK (`= 3`):

    Particles-BGK-CollModel = 1

The **recommended method is ESBGK**. If the simulation contains a gas mixture, a choice for the determination of the transport coefficients is available. The first model uses Wilke's mixture rules (`= 1`) to calculate the gas mixture viscosity and thermal conductivity. The **recommended second model utilizes collision integrals** (derived for the VHS model, `= 2`) to calculate these mixture properties.

    Particles-BGK-MixtureModel    = 2

The vibrational excitation can be controlled with the following flags, including the choice between continuous and quantized vibrational energy.

    Particles-BGK-DoVibRelaxation = T
    Particles-BGK-UseQuantVibEn   = T

An octree cell refinement can be utilized until the given number of particles is reached, which corresponds to an equal refinement in all three directions (x,y,z):

    Particles-BGK-DoCellAdaptation = T
    Particles-BGK-MinPartsPerCell  = 10

It is recommended to utilize at least between 7 and 10 particles per (sub)cell. To enable the cell refinement above a certain number density, the following option can be utilized:

    Particles-BGK-SplittingDens = 1E23

A coupled BGK-DSMC simulation can be enabled, where the BGK method will be utilized if the number density $[\text{m}^{-3}]$ is above a certain value:

    Particles-CoupledBGKDSMC       = T
    Particles-BGK-DSMC-SwitchDens  = 1E22

The flag `Particles-DSMC-CalcQualityFactors` controls the output of quality factors such as mean/maximum relaxation factor (mean: average over a cell, max: maximal value within the octree), max. rotational relaxation factor, which are defined as

$$ \frac{\Delta t}{\tau} < 1,$$

where $\Delta t$ is the chosen time step and $1/\tau$ the relaxation frequency. The time step should be chosen as such that the relaxation factors are below unity. The `BGK_DSMC_Ratio` gives the percentage of the sampled time during which the BGK model was utilized. In a coupled BGK-DSMC simulation this variable indicates the boundary between BGK and DSMC. However, a value below 1 can occur for pure BGK simulations due to low particle numbers, when an element is skipped.

An option is available to utilize a moving average for the variables used in the calculation of the relaxation frequency:

    Particles-BGK-MovingAverage = T

The purpose is to increase the sample size and reduce the noise for steady gas flows. For this, the factor $f$

    Particles-BGK-MovingAverageFac = 0.01

between zero and one must be defined with which the old $M^n$ and newly sampled moments $M$ are weighted to define the moments for the next time step $M^{n+1}$:

$$ M^{n+1}=f M+(1-f) M^n.$$
