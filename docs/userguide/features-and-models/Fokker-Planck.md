# Fokker-Planck Collision Operator

The implementation of the FP-based collision operator is based on the publications by {cite}`Gorji2014` and {cite}`Pfeiffer2017`.
It is a method, which allows the simulation of gas flows in the continuum and transitional regime, where the DSMC method is
computationally too expensive. The collision integral is hereby approximated by a drift and diffusion process

$$  \left.\frac{\partial f}{\partial t}\right|_{\mathrm{coll}}\approx-\sum_{i=1}^3 {\frac{\partial }{\partial v_i}(A_i f)+\frac{1}{2}\sum_{i=1}^3 \sum_{j=1}^3\frac{\partial ^2 }{\partial v_i\partial v_j}(D_{ij}f)}, $$

where $\mathbf{A}$ is the drift vector and $\mathcal{D}$ the diffusion matrix.

The current implementation supports:

- 2 different methods: Cubic and Ellipsoidal Statistical (ES)
- Single species, monatomic and polyatomic gases
- Thermal non-equilibrium with rotational and vibrational excitation (continuous or quantized treatment)
- 2D/Axisymmetric simulations
- Variable time step (adaption of the distribution according to the maximal relaxation factor and linear scaling)

Relevant publications of the developers:

- Implementation of the cubic Fokker-Planck in PICLas {cite}`Pfeiffer2017`
- Comparison of the cubic and ellipsoidal statistical Fokker-Planck {cite}`Jun2019`
- Simulation of a nozzle expansion (including the pressure chamber) with ESBGK, ESFP and coupled ESBGK-DSMC, comparison to
experimental measurements {cite}`Pfeiffer2019a`

To enable the simulation with the FP module, the respective compiler setting has to be activated:

    PICLAS_TIMEDISCMETHOD = FP-Flow

A parameter file and species initialization file is required, analogous to the DSMC setup. It is recommended to utilize a previous
DSMC parameter file to ensure a complete simulation setup. To enable the simulation with the FP methods, select the Fokker-Planck
method, cubic (`=1`) and ES (`=2`):

    Particles-FP-CollModel = 2

The vibrational excitation can be controlled with the following flags, including the choice between continuous and quantized
vibrational energy:

    Particles-FP-DoVibRelaxation = T
    Particles-FP-UseQuantVibEn   = T

An octree cell refinement until the given number of particles is reached can be utilized, which corresponds to an equal
refinement in all three directions (x,y,z):

    Particles-FP-DoCellAdaptation = T
    Particles-FP-MinPartsPerCell  = 10

A coupled FP-DSMC simulation can be enabled, where the FP method will be utilized if the number density $[\text{m}^{-3}]$ is above a certain value:

    Particles-CoupledFPDSMC       = T
    Particles-FP-DSMC-SwitchDens  = 1E22

The flag `Particles-DSMC-CalcQualityFactors` controls the output of quality factors such as mean/maximal relaxation factor (mean:
average over a cell, max: maximal value within the octree), max rotational relaxation factor, which are defined as

$$ \frac{\Delta t}{\tau} < 1,$$

where $\Delta t$ is the chosen time step and $1/\tau$ the relaxation frequency. The time step should be chosen as such that the
relaxation factors are below unity. The `FP_DSMC_Ratio` gives the percentage of the sampled time during which the FP model was utilized.
In a couple FP-DSMC simulation this variable indicates the boundary between FP and DSMC. However, a value below 1 can occur for
pure FP simulations due to low particle numbers, when an element is skipped. Additionally, the Prandtl number utilized by the ESFP
model is given.

