\hypertarget{features_models}{}

# Features \& Models \label{chap:features_models}

The goal of PICLas is to enable to approximation of the complete Boltzmann equation:

$$ \frac{\partial f}{\partial t} + \mathbf{v}\cdot\frac{\partial f}{\partial \mathbf{x}} + \frac{\mathbf{F}}{m}\cdot\frac{\partial f}{\partial \mathbf{v}} = \left.\frac{\partial f}{\partial t}\right|_\mathrm{coll} $$

## Particle Tracking

Three different particle tracking methods are implemented in PICLas. For conventional computations on linear meshes, the following tracking algorithm is recommended:

    TriaTracking = T
    DoRefMapping = F

The option DoRefMapping should be disabled. The two alternative tracking routines and their options are described in the following.

### DoRefMapping (NEEDS UPDATING)

    TriaTracking = F
    DoRefMapping = T

This method is the slowest implemented method for linear grids and large particle movements. A particle is mapped into 
a element to compute the particle position
in the reference space. This test determines in which element a particle is located. Each element has a slightly larger
reference space due to tolerance. Starting from reference values >=1. the best element is found and used for the 
hosting element. In order to take boundary interactions into account, all BC faces in the halo vicinity of the element
are checked for boundary interactions and a boundary condition is performed accordingly. This algorithm has a 
inherent self check. If a boundary condition is not detected, the particle position is located outside of all elements.
A fall-back algorithm is used to recompute the position and boundary interaction. Periodic domains are only possible
for Cartesian meshes. The particle position is used for periodic displacements.

| Option                 | Values     |  Notes                                                  |
|:----------------------:|:----------:|:-------------------------------------------------------:|
| CartesianPeriodic      | T/F        | If a fully periodic box (all 6 sides) is computed, the  |
|                        |            | intersections do not have to be computed. Instead, each |
|                        |            | particle can be simply shifted by the periodic vector.  |
| FastPeriodic           | T/F        | Moves particle the whole periodic distance once, which  |
|                        |            | can be several times the mesh size in this direction.   |


### Tracing  (NEEDS UPDATING)

    TriaTracking = F
    DoRefMapping = F

This method traces the particles throughout the domain. The initial element is determined by computing the intersection
between the particle-element-origin vector and each element face. If non of the six element faces are hit, the particle is 
located inside of this element. Next, the particle trajectory is traced throughout the domain. Hence, each face is checked
for an intersection and a particle mapped accordingly into the neighbor element or perform a boundary condition. This 
algorithm has no inherent self-consistency check. For critical intersections (beginning,end of particle path or close to 
edges of faces) an additional safety check is performed by recomputing the element check and if it fails a re-localization of 
the particle. Particles traveling parallel to faces are in a undefined state and a currently removed. This prints a warning
message. Note, the tracing on periodic meshes works only for non-mpi computations. Periodic displacement requires 
additional coding.


| Option                 | Values     |  Notes                                                  |
|:----------------------:|:----------:|:-------------------------------------------------------:|
| CountNbOfLostParts     | T/F        | Count number of lost particles due to tolerance issues. |
|                        |            | This number is a global number, summed over the full t. |

### Parameters for DoRefMapping and Tracing  (NEEDS UPDATING)

Following parameters can be used for both schemes.

| Option                 | Values     |  Notes                                                  | 
|:----------------------:|:----------:|:-------------------------------------------------------:|
| MeasureTrackTime       | T/F        | Measure the time required for tracking and init local.  |
| RefMappingGuess        | 1-4        | Prediction of particle position in reference space:     |
|                        | 1          | Assumption of a linear element coord system.            |
|                        | 2          | Gauss point which is closest to the particle.           |
|                        | 3          | CL point which is closest to the particle.              |
|                        | 4          | Trival guess: element origin                            |
| RefMappingEps          | 1e-4       | Tolerance of the Newton algorithm for mapping in ref.   |
|                        |            | space. It is the L2 norm of the delta Xi in ref space.  |
| BezierElevation        | 0-50       | Increase polinomial degree of BezierControlPoints to    |
|                        |            | construct a thighter bounding box for each side.        |
| BezierSampleN          | NGeo       | Polynomial degree to sample sides for SurfaceFlux and   |
|                        |            | Sampling of DSMC surface data.                          |
| BezierNewtonAngle      | <PI/2      | Angle to switch between Clipping and a Newton algorithm.|
| BezierClipTolerance    | 1e-8       | Tolerance of Bezier-Clipping and Bezier-Newton          |
| BezierClipHit          | 1e-6       | Tolerance to increase sides and path during Bezier-Algo.|
| BezierSplitLimit       | 0.6        | Minimum degrees of side during clipping. A larger       |
|                        |            | surface is spit in two to increase convergence rate and |
|                        |            | predict several intersections.                          |
| BezierClipMaxIntersec  | 2*NGeo     | Maximum number of roots for curvilinear faces.          |
| epsilontol             | 100*epsM   | Tolerance for linear and bilinear algorithm.            |

### Possible outdated  (NEEDS UPDATING)

| Option                 | Values     |  Notes                                                  | 
|:----------------------:|:----------:|:-------------------------------------------------------:|
| BezierEpsilonBilinear  | T/F        | Tolerance for linear-bilinear side. Obsolet.            |

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