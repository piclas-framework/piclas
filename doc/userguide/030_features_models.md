\hypertarget{features_models}{}

# Features \& Models \label{chap:features_models}

The goal of PICLas is to enable to approximation of the complete Boltzmann equation:

$$ \frac{\partial f}{\partial t} + \mathbf{v}\cdot\frac{\partial f}{\partial \mathbf{x}} + \frac{\mathbf{F}}{m}\cdot\frac{\partial f}{\partial \mathbf{v}} = \left.\frac{\partial f}{\partial t}\right|_{\mathrm{coll}} $$

## Particle Tracking

Three different particle tracking methods are implemented in PICLas and are selected via

    TrackingMethod = triatracking ! Define Method that is used for tracking of
                                  ! particles:
                                  ! refmapping (1): reference mapping of particle
                                  ! position with (bi-)linear and bezier (curved)
                                  ! description of sides. 
                                  ! tracing (2): tracing of particle path with
                                  ! (bi-)linear and bezier (curved) description of
                                  ! sides.
                                  ! triatracking (3): tracing of particle path with
                                  ! triangle-aproximation of (bi-)linear sides.

For conventional computations on (bi-, tri-) linear meshes, the following tracking algorithm is recommended:

    TrackingMethod = triatracking

The two alternative tracking routines and their options are described in the following.

### DoRefMapping

    TrackingMethod = refmapping

This method is the slowest implemented method for linear grids and large particle displacements.
A particle is mapped into a element to compute the particle position in the reference space. 
This test determines in which element a particle is located. Each element has a slightly larger
reference space due to tolerance. Starting from reference values >=1. the best element is found and used for the
hosting element. In order to take boundary interactions into account, all BC faces in the halo vicinity of the element
are checked for boundary interactions and a boundary condition is performed accordingly. This algorithm has an
inherent self-check. If a boundary condition is not detected, the particle position is located outside of all elements.
A fall-back algorithm is then used to recompute the position and boundary interaction. Periodic domains are only possible
for Cartesian meshes. The particle position is used for periodic displacements.

|      Option       | Values |                          Notes                          |
| :---------------: | :----: | :-----------------------------------------------------: |
| CartesianPeriodic |  T/F   | If a fully periodic box (all 6 sides) is computed, the  |
|                   |        | intersections do not have to be computed. Instead, each |
|                   |        | particle can be simply shifted by the periodic vector.  |
|   FastPeriodic    |  T/F   | Moves particle the whole periodic distance once, which  |
|                   |        |  can be several times the mesh size in this direction.  |


### Tracing

    TrackingMethod = tracing

This method traces the particle trajectory throughout the domain. The initial element is determined by computing the intersection
between the particle-element-origin vector and each element face. If none of the six element faces are hit, the particle is
located inside of this element. Next, the particle trajectory is traced throughout the domain. Hence, each face is checked
for an intersection and a particle assigned accordingly to neighbor elements or the interaction with boundary conditions occur. This
algorithm has no inherent self-consistency check. For critical intersections (beginning or end of a particle path or if a particle is located close to
the edges of element faces) an additional safety check is performed by recomputing the element check and if it fails a re-localization of
the particle is required. Particles traveling parallel to element faces are in an undefined state and are currently removed from the computation. 
This leads to a warning message. Note that tracing on periodic meshes works only for non-mpi computations. Periodic displacement requires
additional coding.


|       Option       | Values |                          Notes                          |
| :----------------: | :----: | :-----------------------------------------------------: |
| CountNbOfLostParts |  T/F   | Count number of lost particles due to tolerance issues. |
|                    |        | This number is a global number, summed over the full t. |

### Parameters for DoRefMapping and Tracing  (NEEDS UPDATING)

Following parameters can be used for both schemes.

|        Option         |  Values  |                          Notes                           |
| :-------------------: | :------: | :------------------------------------------------------: |
|   MeasureTrackTime    |   T/F    |  Measure the time required for tracking and init local.  |
|    RefMappingGuess    |   1-4    |   Prediction of particle position in reference space:    |
|                       |    1     |       Assumption of a linear element coord system.       |
|                       |    2     |      Gauss point which is closest to the particle.       |
|                       |    3     |        CL point which is closest to the particle.        |
|                       |    4     |              Trivial guess: element origin               |
|     RefMappingEps     |   1e-4   |  Tolerance of the Newton algorithm for mapping in ref.   |
|                       |          |  space. It is the L2 norm of the delta Xi in ref space.  |
|    BezierElevation    |   0-50   |   Increase polynomial degree of BezierControlPoints to   |
|                       |          |     construct a tighter bounding box for each side.      |
|     BezierSampleN     |   NGeo   |  Polynomial degree to sample sides for SurfaceFlux and   |
|                       |          |              Sampling of DSMC surface data.              |
|   BezierNewtonAngle   |  $<PI/2$   | Angle to switch between Clipping and a Newton algorithm. |
|  BezierClipTolerance  |   1e-8   |      Tolerance of Bezier-Clipping and Bezier-Newton      |
|     BezierClipHit     |   1e-6   | Tolerance to increase sides and path during Bezier-Algo. |
|   BezierSplitLimit    |   0.6    |    Minimum degrees of side during clipping. A larger     |
|                       |          | surface is spit in two to increase convergence rate and  |
|                       |          |              predict several intersections.              |
| BezierClipMaxIntersec |  2*NGeo  |      Maximum number of roots for curvilinear faces.      |
|      epsilontol       | 100*epsM |       Tolerance for linear and bilinear algorithm.       |

### Possible outdated  (NEEDS UPDATING)

|        Option         | Values |                    Notes                     |
| :-------------------: | :----: | :------------------------------------------: |
| BezierEpsilonBilinear |  T/F   | Tolerance for linear-bilinear side. Obsolet. |

## Boundary Conditions - Field Solver

To-do: Modification of boundaries with the PICLas parameter file (order is of importance)

### Maxwell's Equations

To-do

Dielectric -> type 100?

### Poisson's Equation

To-do

### Dielectric Materials

Dielectric material properties can be considered by defining regions (or specific elements)
in the computational domain, where permittivity and permeability constants for linear isotropic 
non-lossy dielectrics are used. The interfaces between dielectrics and vacuum regions must be separated 
by element-element interfaces due to the DGSEM (Maxwell) and HDG (Poisson) solver requirements, but
can vary spatially within these elements.

The dielectric module is activated by setting

    DoDielectric = T

and specifying values for the permittivity and permeability constants

    DielectricEpsR = X
    DielectricMuR = X

Furthermore, the corresponding regions in which the dielectric materials are found must be defined, 
e.g., simple boxes via

    xyzDielectricMinMax  = (/0.0 , 1.0 , 0.0 , 1.0 , 0.0 , 1.0/)

for the actual dielectric region (vector with 6 entries yielding $x$-min/max, $y$-min/max and 
$z$-min/max) or the inverse (vacuum, define all elements which are NOT dielectric) by

    xyzPhysicalMinMaxDielectric = (/0.0 , 1.0 , 0.0 , 1.0 , 0.0 , 1.0/)

Spherical regions can be defined by setting a radius value

    DielectricRadiusValue = X

and special pre-defined regions (which also consider spatially varying material properties) may also be 
used, e.g., 

    DielectricTestCase = FishEyeLens

where the following pre-defined cases are available as given in table \ref{tab:dielectric_test_cases}.

Table: Dielectric Test Cases \label{tab:dielectric_test_cases}

| Option                       | Additional Parameters                           | Notes                                                                                                               |
| :-------------------------:  | :------------------------:                        | :-------------------------------------------------------:                                                           |
| `FishEyeLens`                | none                                            | function with radial dependence: $\varepsilon_{r}=n_{0}^{2}/(1 + (r/r_{max})^{2})^{2}$                              |
| `Circle`                     | `DielectricRadiusValue, DielectricRadiusValueB` | Circular dielectric in x-y-direction (constant in z-direction)  with optional cut-out radius DielectricRadiusValueB |
| `DielectricResonatorAntenna` | `DielectricRadiusValue`                         | Circular dielectric in x-y-direction (only elements with $z>0$)                                                     |
| `FH_lens`                    | none                                            | specific geometry (see `SUBROUTINE SetGeometry` for more information)                                               |

For the Maxwell solver (DGSEM), the interface fluxes between vacuum and dielectric regions can
either be conserving or non-conserving, which is selected by

    DielectricFluxNonConserving = T

which uses non-conserving fluxes. This is recommended for improved simulation results, as described in [@Copplestone2019b].
When particles are to be considered in a simulation, these are generally removed from dielectric
materials during the emission (inserting) stage, but may be allowed to exist within dielectrics by
setting

    DielectricNoParticles = F

which is set true by default, hence, removing the particles.


## Boundary Conditions - Particle Solver

Within the parameter file it is possible to define different particle boundary conditions. The number of boundaries is defined by

    Part-nBounds=2
    Part-Boundary1-SourceName=BC_OPEN
    Part-Boundary1-Condition=open
    Part-Boundary2-SourceName=BC_WALL
    Part-Boundary2-Condition=reflective
    Part-Boundary2-SurfaceModel=2

The `Part-Boundary1-SourceName=` corresponds to the name given during the preprocessing step with HOPR. The available conditions (`Part-Boundary1-Condition=`) are described in the table below.

| Condition    | Description                                                                                                                                                                                 |
| :----------: | :----------------------------------------------------------------------------------------------                                                                      |
| `open`       | Every particle crossing the boundary will be deleted.                                                                                                                                       |
| `reflective` | Allows the definition of specular and diffuse reflection. A perfect specular reflection is performed, if no other parameters are given (discussed in more detail in the following section). |
| `symmetric`  | A perfect specular reflection, without sampling of particle impacts.                                                                                                                        |

For `reflective` boundaries, an additional option `Part-Boundary2-SurfaceModel` is available, that
is used for heterogeneous reactions (reactions have reactants in two or more phases) or secondary electron emission models. These models are described in \ref{sec:chem_reac}.

### Diffuse Wall

Gas-surface interaction can be modelled with the extended Maxwellian model [@Padilla2009], using accommodation coefficients of the form

$$\alpha = \frac{E_i-E_r}{E_i - E_w}$$

where $i$, $r$ and $w$ denote the incident, reflected and wall energy, respectively.  The coefficient `MomentumACC` is utilized to decide whether a diffuse (`MomentumACC` $>R$) or specular reflection (`MomentumACC` $<R$) occurs upon particle impact, where $R=[0,1)$ is a random number. Separate accommodation coefficients can be defined for the translation (`TransACC`), rotational (`RotACC`), vibrational (`VibACC`) and electronic energy (`ElecACC`) accommodation at a constant wall temperature [K].

    Part-Boundary2-SourceName=BC_WALL
    Part-Boundary2-Condition=reflective
    Part-Boundary2-MomentumACC=1.
    Part-Boundary2-WallTemp=300.
    Part-Boundary2-TransACC=1.
    Part-Boundary2-VibACC=1.
    Part-Boundary2-RotACC=1.
    Part-Boundary2-ElecACC=1.

Additionally, a wall velocity [m/s] and voltage [V] can be given

    Part-Boundary2-WallVelo=(/0,0,100/)
    Part-Boundary2-Voltage=100

A linear temperature gradient across a boundary can be defined by supplying a second wall temperature and the start and end vector

    Part-Boundary2-WallTemp2=500.
    Part-Boundary2-TemperatureGradientStart=(/0.,0.,0./)
    Part-Boundary2-TemperatureGradientEnd=(/0.,0.,1./)

Between these two points the temperature will be interpolated, where the start vector corresponds to the first wall temperature, while the end vector to the second wall temperature. Beyond these position values, the first and second temperature will be used as the constant wall temperature, respectively.

### Porous Wall / Pump

The porous boundary condition uses a removal probability to determine whether a particle is deleted or reflected at the boundary. The main application of the implemented condition is to model a pump, according to [@Lei2017]. It is defined by giving the number of porous boundaries and the respective boundary number (`BC=2` corresponds to the `BC_WALL` boundary defined in the previous section) on which the porous condition is.

    Part-nPorousBC=1
    Part-PorousBC1-BC=2
    Part-PorousBC1-Pressure=5.
    Part-PorousBC1-Temperature=300.
    Part-PorousBC1-Type=pump
    Part-PorousBC1-PumpingSpeed=2e-9
    Part-PorousBC1-DeltaPumpingSpeed-Kp=0.1
    Part-PorousBC1-DeltaPumpingSpeed-Ki=0.0

The removal probability is determined through the given pressure [Pa] and temperature [K] at the boundary. A pumping speed can be given as a first guess, however, the pumping speed $S$ [$m^3/s$] will be adapted if the proportional factor ($K_{\mathrm{p}}$, `DeltaPumpingSpeed-Kp`) is greater than zero

$$ S^{n+1}(t) = S^{n}(t) + K_{\mathrm{p}} \Delta p(t) + K_{\mathrm{i}} \int_0^t \Delta p(t') dt',$$

where $\Delta p$ is the pressure difference between the given pressure and the actual pressure at the pump. An integral factor ($K_{\mathrm{i}}$, `DeltaPumpingSpeed-Ki`) can be utilized to mimic a PI controller. The proportional and integral factors are relative to the given pressure. However, the integral factor has not yet been thoroughly tested. The removal probability $\alpha$ is then calculated by

$$\alpha = \frac{S n \Delta t}{N_{\mathrm{pump}} w} $$

where $n$ is the sampled, cell-local number density and $N_{\mathrm{pump}}$ is the total number of impinged particle at the pump during the previous time step. $\Delta t$ is the time step and $w$ the weighting factor. The pumping speed $S$ is only adapted if the resulting removal probability $\alpha$ is between zero and unity. The removal probability is not species-specific.

To reduce the influence of statistical fluctuations, the relevant macroscopic values (pressure difference $\Delta p$ and number density $n$) can be sampled for $N$ iterations by defining (for all porous boundaries)

    Part-PorousBC-IterationMacroVal=10

A porous region on the specified boundary can be defined. At the moment, only the `circular` option is implemented. The origin of the circle/ring on the surface and the radius have to be given. In the case of a ring, a maximal and minimal radius is required (`-rmax` and `-rmin`, respectively), while for a circle only the input of maximal radius is sufficient.

    Part-PorousBC1-Region=circular
    Part-PorousBC1-normalDir=1
    Part-PorousBC1-origin=(/5e-6,5e-6/)
    Part-PorousBC1-rmax=2.5e-6

The absolute coordinates are defined as follows for the respective normal direction.

| Normal Direction | Coordinates |
| :--------------: | :---------: |
|      x (=1)      |    (y,z)    |
|      y (=2)      |    (z,x)    |
|      z (=3)      |    (x,y)    |

Using the regions, multiple pumps can be defined on a single boundary. Additionally, the BC can be used as a sensor by defining the respective type:

    Part-PorousBC1-BC=3
    Part-PorousBC1-Pressure=5.
    Part-PorousBC1-Temperature=300.
    Part-PorousBC1-Type=sensor

Together with a region definition, a pump as well as a sensor can be defined on a single and/or multiple boundaries, allowing e.g. to determine the pressure difference between the pump and a remote area of interest.

### Surface Chemistry \label{sec:chem_reac}

Modelling of reactive surfaces is enabled by setting `Part-BoundaryX-Condition=reflective` and an
appropriate particle boundary surface model `Part-BoundaryX-SurfaceModel`.
The available conditions (`Part-BoundaryX-SurfaceModel=`) are described in the table below.

| SurfaceModel | Description                                                                                                                                                                     |
| :----------: | :-----------------------------------------------------------------                                                                                                              |
| 0 (default)  | Standard extended Maxwellian scattering                                                                                                                                         |
| 2            | Simple recombination on surface collision, where an impinging particle as given by Ref. [@Reschke2019].                                                                         |
| 3            | Kinetic Monte Carlo surface: Replicates surfaces with a specified lattice structure, either fcc(100) or fcc(111) and models complete catalysis as given by Ref. [@Reschke2019]. |
| 5            | Secondary electron emission as given by Ref. [@Levko2015].                                                                                                                      |
| 101          | Evaporation from surfaces according to a Maxwellian velocity distribution.                                                                                                      |
| 102          | Evaporation according to MD-fitted velocity distributions.                                                                                                                      |

For surface sampling output, where the surface is split into, e.g., $3\times3$ sub-surfaces, the following parameters mus be set

    BezierSampleN = 3
    DSMC-nSurfSample = 3
    Part-WriteMacroSurfaceValues = T
    Particles-DSMC-CalcSurfaceVal = T
    Part-IterationForMacroVal = 200

where `BezierSampleN=DSMC-nSurfSample`. In this example, sampling is performed over 200 interations.

### Deposition of Charges on Dielectric Surfaces

Charged particles can be absorbed (or reflected and leave their charge behind) at dielectric surfaces
when using the deposition method `cell_volweight_mean`. The boundary can be used by specifying

    ```
    Part-Boundary1-Condition         = reflective
    Part-Boundary1-Dielectric        = T
    Part-Boundary1-NbrOfSpeciesSwaps = 3
    Part-Boundary1-SpeciesSwaps1     = (/1,0/) ! e-
    Part-Boundary1-SpeciesSwaps2     = (/2,2/) ! Ar
    Part-Boundary1-SpeciesSwaps3     = (/3,2/) ! Ar+
    ```

which sets the boundary dielectric and the given species swap parameters effectively remove
electrons ($e^{-}$) on impact, reflect $Ar$ atoms and neutralize $Ar^{+}$ ions by swapping these to $Ar$ atoms.
Note that currently only singly charged particles can be handled this way. When multiple charged
particles would be swapped, their complete charge mus be deposited at the moment.

The boundary must also be specified as an *inner* boundary via

    BoundaryName                     = BC_INNER
    BoundaryType                     = (/100,0/)

or directly in the *hopr.ini* file that is used for creating the mesh.

## Particle Initialization & Emission

The following section gives an overview of the available options regarding the definition of species and particle initialization and emission. Simulation particles can be inserted initially within the computational domain and/or emitted at every time step. First of all, the number of species is defined by

    Part-nSpecies=1
    Part-MaxParticleNumber=1000000

The maximum particle number is defined per core and should be chosen according to the number of simulation particles you expect, including a margin to account for imbalances due transient flow features and/or the occurrence of new particles due to chemical reactions. Example: A node of a HPC cluster has 2 CPUs, each has 12 cores. Thus, the node has 24 cores that share a total of 128GB RAM. Allocating 1000000 particles per core means, you can simulate up to 24 Million particles on a single node in this example (assuming an even particle distribution). The limiting factor here is the amount of RAM available per node.

Regardless whether a standalone PIC, DSMC, or a coupled simulation is performed, the atomic mass [kg], the charge [C] and the weighting factor $w$ [-] are required for each species.

    Part-Species1-MassIC=5.31352E-26
    Part-Species1-ChargeIC=0.0
    Part-Species1-MacroParticleFactor=5E2

Species that are not part of the initialization or emission but might occur as a result of e.g. chemical reactions should also be defined with these parameters.

**WORK IN PROGRESS**

Different velocity distributions are available for the initialization of particles.

| Distribution | Description                                             |
| ------------ | ------------------------------------------------------- |
| maxwell      | Maxwell-Boltzmann distribution                          |
| maxwell_lpn  | Maxwell-Boltzmann distribution for low particle numbers |
| ...          | many many more                                          |

### Initialization

### Surface Flux

A surface flux enables the emission of particles at a boundary in order to simulate, e.g. a free-stream. They are defined species-specific and can overlap. First, the number of surface fluxes has to be given

    Part-Species1-nSurfaceFluxBCs=1

The surface flux is mapped to a certain boundary by giving its boundary number (e.g. `BC=1` corresponds to the previously defined boundary `BC_OPEN`)

    Part-Species1-Surfaceflux1-BC=1
    Part-Species1-Surfaceflux1-VeloIC=1500
    Part-Species1-Surfaceflux1-VeloVecIC=(/-1.0,0.0,0.0/)
    Part-Species1-Surfaceflux1-velocityDistribution=maxwell_lpn
    Part-Species1-Surfaceflux1-MWTemperatureIC=300.
    Part-Species1-Surfaceflux1-PartDensity=1E20

The drift velocity is defined by the direction vector `VeloVecIC`, which is a unit vector, and a velocity magnitude [m/s]. The thermal velocity of particle is determined based on the defined velocity distribution and the given translation temperature `MWTemperatureIC` [K]. Finally, the 'real' number density is defined by `PartDensity` [1/m$^3$], from which the actual number of simulation particles will be determined (depending on the chosen weighting factor).

In the case of molecules, the rotational and vibrational temperature [K] have to be defined. If electronic excitation is considered, the electronic temperature [K] has to be defined

    Part-Species1-Surfaceflux1-TempRot=300.
    Part-Species1-Surfaceflux1-TempVib=300.
    Part-Species1-Surfaceflux1-TempElec=300.


#### Circular Inflow

The emission of particles from a surface flux can be limited to the area within a circle or a ring. The respective boundary has to coincide or be parallel to the xy-, xz, or yz-planes. This allows to define inflow boundaries without specifically meshing the geometrical feature, e.g. small orifices. The feature can be enabled per species and surface flux

    Part-Species1-Surfaceflux1-CircularInflow=TRUE

The normal direction of the respective boundary has to be defined by

    Part-Species1-Surfaceflux1-axialDir=1

Finally, the origin of the circle/ring on the surface and the radius have to be given. In the case of a ring, a maximal and minimal radius is required (`-rmax` and `-rmin`, respectively), while for a circle only the input of maximal radius is sufficient.

    Part-Species1-Surfaceflux1-origin=(/5e-6,5e-6/)
    Part-Species1-Surfaceflux1-rmax=2.5e-6
    Part-Species1-Surfaceflux1-rmin=1e-6

The absolute coordinates are defined as follows for the respective normal direction.

| Normal Direction | Coordinates |
| :--------------: | :---------: |
|      x (=1)      |    (y,z)    |
|      y (=2)      |    (z,x)    |
|      z (=3)      |    (x,y)    |

Multiple circular inflows can be defined on a single boundary through multiple surface fluxes, e.g. to enable the simulation of multiple inlets on a chamber wall.

#### Adaptive Boundaries

Different adaptive boundaries can be defined as a part of a surface flux to model subsonic in- and outflows, where the emission is adapted based on the prevalent conditions at the boundary. The modelling is based on the publications by [@Farbar2014] and [@Lei2017].

    Part-Species1-Surfaceflux1-Adaptive=TRUE
    Part-Species1-Surfaceflux1-Adaptive-Type=1

An overview over the available types is given below.

 * `Type=1`: Constant static pressure and temperature inlet, defined as Type 1 in Ref. [@Farbar2014]
 * `Type=2`: Constant static pressure outlet, defined as Type 1 in Ref. [@Farbar2014]
 * `Type=3`: Constant mass flow and temperature inlet, where the given mass flow and sampled velocity are used to determine the number of particles for the surface flux. It requires the BC to be defined as `open`. Defined as Type 2 in Ref. [@Farbar2014]
 * `Type=4`: Constant mass flow inlet and temperature inlet, where number of particles to be inserted is determined directly from the mass flow and the number of particles leaving the domain, $N_{\mathrm{in}}=N_{\dot{m}} + N_{\mathrm{out}}$. Defined as cf_3 in Ref. [@Lei2017]

Depending of the type of the chosen boundary type either the mass flow [kg/s] or the static pressure [Pa] have to be given

    Part-Species1-Surfaceflux1-Adaptive-Massflow=1.00E-14
    Part-Species1-Surfaceflux1-Adaptive-Pressure=10

The adaptive boundaries require the sampling of macroscopic properties such as flow velocity at the boundary. To compensate for the statistical fluctuations a relaxation factor $f_{\mathrm{relax}}$ is utilized and the current value of the sampled variable $v^{n}$ is updated according to

$$v^{n}= (1-f_{\mathrm{relax}})\,v^{n-1} + f_{\mathrm{relax}} v^{\mathrm{samp}} $$

The relaxation factor $f_{\mathrm{relax}}$ is defined by

    Part-AdaptiveWeightingFactor = 0.001

The adaptive particle emission can be combined with the circular inflow feature. In this context when the area of the actual emission circle/ring is very small, it is preferable to utilize the `Type=4` constant mass flow condition. `Type=3` assumes an open boundary and accounts for particles leaving the domain through that boundary already when determining the number of particles to be inserted. As a result, this method tends to overpredict the given mass flow, when the emission area is very small and large sample size would be required to have enough particles that leave the domain through the emission area. For the `Type=4` method, the actual number of particles leaving the domain through the circular inflow is counted and the mass flow adapted accordingly, thus the correct mass flow can be reproduced.

Additionally, the `Type=4` method can be utilized in combination with a reflective boundary condition to model diffusion and leakage (e.g. in vacuum tanks) based on a diffusion rate $Q$ [Pa m$^3$ s$^{-1}$]. The input mass flow [kg s$^{-1}$] for the simulation is then determined by

$$\dot{m} = \frac{QM}{1000RT},$$

where $R=8.314$ J mol$^{-1}$K$^{-1}$ is the gas constant, $M$ the molar mass in [g mol$^{-1}$] and $T$ is the gas temperature [K].

To verify the resulting mass flow rate of an adaptive surface flux, the following option can be enabled

    CalcMassflowRate = T

This will output a species-specific mass flow rate [kg s^$-1$] for each surface flux condition in the `PartAnalyze.csv`, which gives the current mass flow for the time step. Positive values correspond to a net mass flux into the domain and negative values vice versa. It should be noted that while multiple adaptive boundaries are possible, adjacent boundaries that share a mesh element should be avoided or treated carefully.

#### Missing descriptions

SimpleRadialVeloFit, ReduceNoise, DoForceFreeSurfaceFlux

DoPoissonRounding: [@Tysanner2004]

AcceptReject, ARM_DmaxSampleN: [@Garcia2006]

## Particle-In-Cell \label{sec:pic}

### Charge and Current Deposition \label{sec:pic_deposition}

Charge and current deposition can be performed using different methods, among others, shape
functions, B-splines or locally volume-weighted approaches.

#### Linear Distribution Over Cell Interfaces
A linear deposition method that also considers neighbouring elements can be selected by

    PIC-Deposition-Type = cell_volweight_mean

The method also considers the corner nodes of each element to which all neighbouring elements
contribute, hence, resulting in a non-local deposition scheme.

#### Shape Function

High-order field solvers require deposition methods that reduce the noise, e.g., shape functions [@Jacobs2006]. The standard 3D shape function is selected by

    PIC-Deposition-Type = shape_function

or

    PIC-Deposition-Type = shape_function_simple

where `shape_function_simple` is faster for small numbers of elements per processor (high parallelization).

The shape function sphere might be truncated at walls or open boundaries, which can be prevented by
using a local deposition method near boundaries. The deposition of particles in elements where the shape
function might be truncated is changed to *cell_volweight* for these elements via

    PIC-shapefunction-local-depo-BC = T

The following polynomial isotropic shape functions are all designed to be used in three dimensions, where reductions to 2D and 1D are applied.

##### Shape Function 1D
A one-dimensional shape function in $x$-direction is given by

$$
S_{1D}(r,R,\alpha)=\frac{\Gamma(\alpha+3/2)}{\sqrt{\pi}R\Gamma(\alpha+1)\Delta y \Delta z}\left( 1-\left( \frac{r}{R} \right)^{2} \right)^{\alpha}~,
$$

which is normalized to give $\int_{z_{1}}^{z_{2}}\int_{y_{1}}^{y_{2}}\int_{-R}^{R}S_{1D}(r,R,\alpha)dxdydz=1$,
where the radius ${r=|\boldsymbol{x}-\boldsymbol{x}_{n}|=|x-x_{n}|}$ is the distance between the position of the
grid point at position $\boldsymbol{x}$ and the $n$-th particle at position $\boldsymbol{x}_{n}$,
$R$ is the cut-off radius, $\Delta y=y_{2}-y_{1}$ and $\Delta z=z_{2}-z_{1}$ are the domain lengths in $y$- and $z$-direction,
respectively, and $\Gamma(z)$ is the gamma function given by

$$
  \Gamma(z)=\int_{0}^{\infty}x^{z-1}\exp(-x)dx~.
$$

The direction in which deposition is performed is chosen via

    PIC-shapefunction1d-direction = 1 ! for x-direction
                                    2 ! for y-direction
                                    3 ! for z-direction


##### Shape Function 2D
A two-dimensional shape function in $x$-$y$-direction is given by

$$
S_{2D}(r,R,\alpha)=\frac{\alpha+1}{\pi R^{2} \Delta z}\left( 1-\left( \frac{r}{R} \right)^{2} \right)^{\alpha}~,
$$

which is normalized to give $\int_{z_{1}}^{z_{2}}\int_{0}^{2\pi}\int_{0}^{R}S_{2D}(r,R,\alpha)rdr d\phi d\theta=1$,
where the radius ${r=|\boldsymbol{x}-\boldsymbol{x}_{n}|}$ is the distance between the position of the
grid point at position $\boldsymbol{x}$ and the $n$-th particle at position $\boldsymbol{x}_{n}$,
$R$ is the cut-off radius and $\Delta z=z_{2}-z_{1}$ is the domain length in $z$-direction.
The perpendicular direction to the two axes, in which deposition is performed is chosen via

    PIC-shapefunction1d-direction = 1 ! for const. depo in x-direction
                                    2 ! for const. depo in y-direction
                                    3 ! for const. depo in z-direction

when the charge is to be deposited const. along the $x$- or $y$- or $z$-direction.
If the charge is to be deposited over the area instead of the volume, the flag

    PIC-shapefunction-3D-deposition=F

must be set, which simply sets $\Delta z=1$ for the example described above.

##### Shape Function 3D
A three-dimensional shape function in $x$-$y$-direction is given by [@Stock2012]

$$
S_{3D}(r,R,\alpha)=\frac{\Gamma(\alpha+5/2)}{\pi^{3/2}R^{3}\Gamma(\alpha+1)}\left( 1-\left( \frac{r}{R} \right)^{2} \right)^{\alpha}~,
$$

which is normalized to give $\int_{0}^{\pi}\int_{0}^{2\pi}\int_{0}^{R}S_{2D}(r,R,\alpha)r^{2}\sin(\phi)dr d\phi d\theta=1$,
where the radius ${r=|\boldsymbol{x}-\boldsymbol{x}_{n}|}$ is the distance between the position of the
grid point at position $\boldsymbol{x}$ and the $n$-th particle at position $\boldsymbol{x}_{n}$ and
$R$ is the cut-off radius.

## Background Field \label{sec:superB}

Certain application cases allow the utilization of a constant magnetic background field. The magnetic field resulting from certain types of coils and permanent magnets can be calculated during the initialization within PICLas or with the standalone tool **superB** (see Section \ref{sec:compileroptions} for compilation).

The background field can be enabled by

    PIC-BG-Field = T

The first option is to use a previously calculated background field. It can be read-in with

    PIC-BGFileName = BField.h5
    PIC-NBG = 1
    PIC-BGFieldScaling = 1.

Additionally, the polynomial degree for the background field can be set by ``PIC-NBG`` and might differ from the actually read-in polynomial degree. Optionally, the read-in field can be scaled by the last of the three parameters above.

The second option is to calculate the magnetic field during the initialization, which will produce an output of the field. The calculation is enabled by

    PIC-CalcBField = T

For this purpose, different coil and permanent magnet geometries can be defined. For visualization purposes, the geometry of the respective coils and permanent magnets can be directly written out as a VTK with

    PIC-CalcBField-OutputVTK = T

In the following the parameters for different coils and permanent magnets based on the implementation by Hinsberger [@Hinsberger2017] are presented.

### Magnetic Field by Permanent Magnets

First, the total number of permanent magnets has to be defined and the type selected. Options are `cuboid`, `sphere`, `cylinder` and `conic`.

    NumOfPermanentMagnets = 1
    PermanentMagnet1-Type = cuboid
                            sphere
                            cylinder
                            conic

All options require the input of a base/origin vector, a number of discretization nodes (results in a different number of total points depending on the chosen geometry) and a magnetisation in [A/m]

    PermanentMagnet1-BasePoint = (/0.,0.,0./)
    PermanentMagnet1-NumNodes = 10
    PermanentMagnet1-Magnetisation = (/0.,0.,1./)

The geometries require different input parameters given below

    ! Three vectors spanning the cuboid
    PermanentMagnet1-BaseVector1 = (/1.,0.,0./)
    PermanentMagnet1-BaseVector2 = (/0.,1.,0./)
    PermanentMagnet1-BaseVector3 = (/0.,0.,1./)
    ! Radius required for a spherical, cylindrical and conical magnet
    PermanentMagnet1-Radius = 1.
    ! Height vector required for a cylindrical and conical magnet
    PermanentMagnet1-HeightVector = (/0.,0.,1./)
    ! Second radius only required for a conical magnet
    PermanentMagnet1-Radius2 = 1.

### Magnetic Field by Coils

The total number of coils and the respective type of the cross-section (`linear`,`circle`,`rectangle`,`custom`) is defined by

    NumOfCoils = 1
    Coil1-Type = linear
                 circle
                 rectangle
                 custom

All options require the input of a base/origin vector, a length vector (vector normal to the cross-section of the coil) and the current in [A]

    Coil1-BasePoint = (/0.0,0.0,0.0/)
    Coil1-LengthVector = (/0.0,1.0,0.0/)
    Coil1-Current = 1.

The first option `linear` represents a simple linear conductor (e.g. a straight wire) and requires only the input of a number of discretization points

    Coil1-NumNodes = 5

The other three types, which are actually coils, are described by the number of loops and the number of discretization points per loop

    Coil1-LoopNum = 10
    Coil1-PointsPerLoop = 10

The cross-section of the coil is defined in a plane normal to the `-LengthVector`. A circular coil cross-section requires simply the input of a radius while a rectangular coil cross-section is spanned by two vectors (`-RectVec1` and `-RectVec2`) and an additional vector, which must be orthogonal to the `-LengthVector` to define the orientation of the cross-section (`-AxisVec1`). In these two cases, the base/origin vector defines the middle point of the cross-section.

    ! Circular coil cross-section
    Coil1-Radius = 1.
    ! Rectangular coil cross-section
    Coil1-RectVec1 = (/1.0,0.0/)
    Coil1-RectVec2 = (/0.0,1.0/)
    Coil1-AxisVec1 = (/0.0,0.0,1.0/)

The last cross-section type `custom` allows the definition of a cross-section as a combination of multiple linear (`line`) and circular (`circle`) segments and also requires an additional vector to define the orientation of the cross-section (`-AxisVec1`)

    Coil1-AxisVec1 = (/0.0,0.0,1.0/)
    Coil1-NumOfSegments = 3
    ! Linear segment defined by 
    Coil1-Segment1-SegmentType = line
    Coil1-Segment1-NumOfPoints = 5
    Coil1-Segment1-LineVector = (/1.0,1.0/)
    ! Circular segment connected to the previous segment
    Coil1-Segment2-SegmentType = circle
    Coil1-Segment2-NumOfPoints = 5
    Coil1-Segment2-Radius = 1.
    Coil1-Segment2-Phi1 = 90.
    Coil1-Segment2-Phi2 = 0.
    ! Linear segment connected to the previous segment, closing the cross-section
    Coil1-Segment3-SegmentType = line
    Coil1-Segment3-NumOfPoints = 5
    Coil1-Segment3-LineVector = (/-2.0,0.0/)

The `-NumOfPoints` controls the number of discretization points per segment. A linear segment is simply described by a vector in the cross-section plane. The circular segment is defined by a radius and the initial as well as final angle of the segment. It should be noted that the base point defines the start of the first segment as opposed to the circular and rectangular cross-sections, where it is the middle point of the cross-section.

## Direct Simulation Monte Carlo

To enable the simulation with DSMC, an appropriate time discretization method including the DSMC module should be chosen before the code compilation. A stand-alone DSMC simulation can be enabled by compiling PICLas with the following parameter

    PICLAS_TIMEDISCMETHOD = DSMC

The DSMC method can then be enabled in the parameter file by

    UseDSMC = T

Additionally, the number of simulated physical models depending on the application can be controlled through

    Particles-DSMC-CollisMode = 1   ! Elastic collisions only
                                2   ! Internal energy exchange
                                3   ! Chemical reactions

`CollisMode = 1` can be utilized for the simulation of a non-reactive, cold atomic gas, where no chemical reactions or electronic excitation is expected. `CollisMode = 2` should be chosen for non-reactive diatomic gas flows to include the internal energy exchange (by default including the rotational and vibrational energy treatment). Finally, reactive gas flows can be simulated with `CollisMode = 3`. The following sections describe the required definition of species parameter (Section \ref{sec:dsmc_species}), the parameters for the internal energy exchange (Section \ref{sec:dsmc_relaxation}) and chemical reactions (Section \ref{sec:dsmc_chemistry}).

The simulation time step $\Delta t$ is defined by

    Particles-ManualTimeStep = 1.00E-7

### Macroscopic Restart \label{sec:macro_restart}

The so-called macroscopic restart, allows to restart the DSMC simulation by using a DSMC output file of a previous simulation run (the regular state file has still to be supplied). This enables to change the weighting factor, without beginning a new simulation.

    Particles-MacroscopicRestart = T
    Particles-MacroscopicRestart-Filename = Test_DSMCState.h5

The particle velocity distribution within the domain is then generated assuming a Maxwell-Boltzmann distribution, using the translational temperature per direction of each species per cell. The rotational and vibrational energy per species is initialized assuming an equilibrium distribution.

### Variable Time Step \label{sec:vartimestep}

A spatially variable time step (VTS) can be activated for steady-state DSMC, BGK and FP simulations, where two options are currently available and described in the following:

* Distribution: use a simulation result to adapt the time step in order to resolve physical parameters (e.g. collision frequency)
* Linear scaling: use a linearly increasing/decreasing time step along a given direction

#### Distribution

The first option is to adapt the time step during a simulation restart based on certain parameters of the DSMC/BGK/FP simulation such as maximal collision probability (DSMC), mean collision separation distance over mean free path (DSMC), maximal relaxation factor (BGK/FP) and particle number. This requires the read-in of a DSMC state file that includes DSMC quality factors (see Section \ref{sec:dsmc_quality}).

    Part-VariableTimeStep-Distribution = T
    Part-VariableTimeStep-Distribution-Adapt = T
    Part-VariableTimeStep-Distribution-MaxFactor = 1.0
    Part-VariableTimeStep-Distribution-MinFactor = 0.1
    Part-VariableTimeStep-Distribution-MinPartNum = 10          ! Optional
    ! DSMC only
    Part-VariableTimeStep-Distribution-TargetMCSoverMFP = 0.3   ! Default = 0.25
    Part-VariableTimeStep-Distribution-TargetMaxCollProb = 0.8  ! Default = 0.8
    ! BGK/FP only
    Part-VariableTimeStep-Distribution-TargetMaxRelaxFactor = 0.8
    ! Restart from a given DSMC state file (Disable if no adaptation is performed!)
    Particles-MacroscopicRestart = T
    Particles-MacroscopicRestart-Filename = Test_DSMCState.h5

The second flag allows to enable/disable the adaptation of the time step distribution. Typically, a simulation would be performed until a steady-state (or close to it, e.g. the particle number is not increasing significantly anymore) is reached with a uniform time step. Then a restart with the above options would be performed, where the time step distribution is adapted using the DSMC output of the last simulation. Now, the user can decide to continue adapting the time step with the subsequent DSMC outputs (Note: Do not forget to update the DSMCState file name!) or to disable the adaptation and to continue the simulation with the distribution from the last simulation (the adapted particle time step is saved within the regular state file). It should be noted that if after a successful restart at e.g. $t=2$, and the simulation fails during the runtime at $t=2.5$ before the next state file could be written out at $t=3$, an adaptation for the next simulation attempt shoud NOT be performed as the adapted time step is stored in the output of new restart file at the restart time $t=2$. Restart files from which the restart is performed are overwritten after a successful restart.

The `MaxFactor` and `MinFactor` allow to limit the adapted time step within a range of $f_{\mathrm{min}} \Delta t$ and $f_{\mathrm{max}} \Delta t$. The time step adaptation can be used to increase the number of particles by defining a minimum particle number (e.g `MinPartNum` = 10, optional). For DSMC, the parameters `TargetMCSoverMFP` (ratio of the mean collision separation distance over mean free path) and `TargetMaxCollProb` (maximum collision probability) allow to modify the target values for the adaptation. For the BGK and FP methods, the time step can be adapted according to a target maximal relaxation frequency.

The last two flags enable to initialize the particles distribution from the given DSMC state file, using the macroscopic properties such as flow velocity, number density and temperature (see Section \ref{sec:macro_restart}). Strictly speaking, the VTS procedure only requires the `Filename` for the read-in of the aforementioned parameters, however, it is recommended to perform a macroscopic restart to initialize the correct particle number per cells. Otherwise, cells with a decreased/increased time step will require some time until the additional particles have reached/left the cell.

The time step adaptation can also be utilized in coupled BGK-DSMC simulations, where the time step will be adapted in both regions according to the respective criteria as the BGK factors are zero in the DSMC region and vice versa. Attention should be payed in the transitional region between BGK and DSMC, where the factors are potentially calculated for both methods. Here, the time step required to fulfil the maximal collision probability criteria will be utilized as it is the more stringent one.

#### Linear scaling

The second option is to use a linearly increasing time step along a given direction. This option does not require a restart or a previous simulation result. Currently, only the increase of the time step along the **x-direction** is implemented. With the start point and end point, the region in which the linear increase should be performed can be defined. To define the domain border as the end point in maximal x-direction, the vector `(/-99999.,0.0,0.0/)` should be supplied. Finally, the `ScaleFactor` defines the maximum time step increase towards the end point $\Delta t (x_{\mathrm{end}})=f \Delta t$.

    Part-VariableTimeStep-LinearScaling = T
    Part-VariableTimeStep-ScaleFactor   = 2
    Part-VariableTimeStep-Direction     =      (/1.0,0.0,0.0/)
    Part-VariableTimeStep-StartPoint    =     (/-0.4,0.0,0.0/)
    Part-VariableTimeStep-EndPoint      =  (/-99999.,0.0,0.0/)

Besides DSMC, the linear scaling is available for the BGK and FP method. Finally, specific options for 2D/axisymmetric simulations are discussed in Section \ref{sec:2DAxi_vts}.

### 2D/Axisymmetric Simulation \label{sec:2DAxi}

For two-dimensional and axisymmetric cases, the computational effort can be greatly reduced. Two-dimensional and axisymmetric simulations require a mesh in the $xy$-plane, where the $x$-axis is the rotational axis and $y$ ranges from zero to a positive value. Additionally, the mesh shall be centered around zero in the $z$-direction with a single cell row, such as that $|z_{\mathrm{min}}|=|z_{\mathrm{max}}|$. The rotational symmetry axis shall be defined as a separate boundary with the `symmetric_axis` boundary condition

Part-Boundary4-SourceName=SYMAXIS
Part-Boundary4-Condition=symmetric_axis

The boundaries (or a single boundary definition for both boundary sides) in the $z$-direction should be defined as symmetry sides with the `symmetric` condition

Part-Boundary5-SourceName=SYM
Part-Boundary5-Condition=symmetric

To enable two-dimensional simulations, the following flag is required

    Particles-Symmetry2D=T

It should be noted that the two-dimensional mesh assumes a length of $\Delta z = 1$, regardless of the actual dimension in $z$. Therefore, the weighting factor should be adapted accordingly.

To enable axisymmetric simulations, the following flag is required

    Particles-Symmetry2DAxisymmetric=T

To fully exploit rotational symmetry, a radial weighting can be enabled, which will linearly increase the weighting factor $w$ towards $y_{\mathrm{max}}$ (i.e. the domain border in $y$-direction), depending on the current $y$-position of the particle.

    Particles-RadialWeighting=T
    Particles-RadialWeighting-PartScaleFactor=100

A radial weighting factor of 100 means that the weighting factor at $y_{\mathrm{max}}$ will be $100w$. Although greatly reducing the number of particles, this introduces the need to delete and create (in the following "clone") particles, which travel upwards and downwards in the $y$-direction, respectively. If the new weighting factor is smaller than the previous one, a cloning probability is calculated by

$$ P_{\mathrm{clone}} = \frac{w_{\mathrm{old}}}{w_{\mathrm{new}}} - \mathrm{INT}\left(\frac{w_{\mathrm{old}}}{w_{\mathrm{new}}}\right)\qquad \text{for}\quad w_{\mathrm{new}}<w_{\mathrm{old}}.$$

For the deletion process, a deletion probability is calculated, if the new weighting factor is greater than the previous

$$ P_{\mathrm{delete}} = 1 - P_{\mathrm{clone}}\qquad \text{for}\quad w_{\mathrm{old}}<w_{\mathrm{new}}.$$

If the ratio between the old and the new weighting factor is $w_{\mathrm{old}}/w_{\mathrm{new}}> 2$, the time step or the radial weighting factor should be reduced as the creation of more than one clone per particle per time step is not allowed. The same applies if the deletion probability is above $0.5$.

For the cloning procedure, two methods are implemented, where the information of the particles to be cloned are stored for a given number of iterations (`CloneDelay=10`) and inserted at the old position. The difference is whether the list is inserted chronologically (`CloneMode=1`) or randomly (`CloneMode=2`) after the first number of delay iterations.

    Particles-RadialWeighting-CloneMode=2
    Particles-RadialWeighting-CloneDelay=10

This serves the purpose to avoid the so-called particle avalanche phenomenon [@Galitzine2015], where clones travel on the exactly same path as the original in the direction of a decreasing weight. They have a zero relative velocity (due to the same velocity vector) and thus a collision probability of zero. Combined with the nearest neighbor pairing, this would lead to an ever-increasing number of identical particles travelling on the same path. An indicator how often identical particle pairs are encountered per time step during collisions is given as an output (`2D_IdenticalParticles`, to enable the output see Section \ref{sec:dsmc_quality}). Additionally, it should be noted that a large delay of the clone insertion might be problematic for time-accurate simulations. However, for the most cases, values for the clone delay between 2 and 10 should be sufficient to avoid the avalance phenomenon.

Another issue is the particle emission on large sides in $y$-dimension close to the rotational axis. As particles are inserted linearly along the $y$-direction of the side, a higher number density is inserted closer to the axis. This effect is directly visible in the free-stream in the cells downstream, when using mortar elements, or in the heatflux (unrealistic peak) close to the rotational axis. It can be avoided by splitting the surface flux emission side into multiple subsides with the following flag (default value is 20)

    Particles-RadialWeighting-SurfFluxSubSides = 20

An alternative to the particle position-based weighting is the cell-local radial weighting, which can be enabled by

    Particles-RadialWeighting-CellLocalWeighting = T

However, this method is not preferable if the cell dimensions in $y$-direction are large, resulting in numerical artifacts due to the clustered cloning processes at cell boundaries.

Besides DSMC, 2D/axisymmetric simulations are also possible the BGK/FP particle method with the same parameters as discussed above (for more informatino about the BGK and FP methods, see Section \ref{sec:continuum}).

#### Variable Time Step: Linear scaling \label{sec:2DAxi_vts}

The linear scaling of the variable time step is implemented slightly different to the 3D case. Here, a particle-based time step is used, where the time step of the particle is determined on its current position. The first scaling is applied in the radial direction, where the time step is increased towards the radial domain border. Thus, $\Delta t (y_{\mathrm{max}}) = f \Delta t$ and $\Delta t (y_{\mathrm{min}} = 0) = \Delta t$.

    Part-VariableTimeStep-LinearScaling = T
    Part-VariableTimeStep-ScaleFactor = 2

Additionally, the time step can be varied along the x-direction by defining a "stagnation" point, towards which the time step is decreased from the minimum x-coordinate ($\Delta t (x_{\mathrm{min}}) = f_{\mathrm{front}}\Delta t$) and away from which the time step is increased again towards the maximum x-coordinate ($\Delta t (x_{\mathrm{max}}) = f_{\mathrm{back}}\Delta t$). Therefore, only at the stagnation point, the time step defined during the initialization is used.

    Part-VariableTimeStep-Use2DFunction = T
    Part-VariableTimeStep-StagnationPoint = 0.0
    Part-VariableTimeStep-ScaleFactor2DFront = 2.0
    Part-VariableTimeStep-ScaleFactor2DBack = 2.0

### Species Definition \label{sec:dsmc_species}

For the DSMC simulation, additional species-specific parameters (collision model parameters, characteristic vibrational temperature, etc.) are required. This file is also utilized for the definition of chemical reactions paths. To define a species, its name as well as an `InteractionID` have to be defined

    Part-Species1-SpeciesName = CH4
    Part-Species1-InteractionID = 2

The name is at the moment only utilized to retrieve the electronic energy levels from an additional database. The interaction ID determines the type of a species as follows

|   ID | Type                               |
| ---: | ---------------------------------- |
|    1 | Atom                               |
|    2 | Molecule (diatomic and polyatomic) |
|    4 | Electron                           |
|   10 | Atomic Ion                         |
|   20 | Molecular Ion                      |

Depending on the utilized collision model, different parameters have to be defined. As an example, the parameters for the Variable Hard Sphere (VHS) collision cross-section model are be defined by the temperature exponent $\omega$, reference temperature $T_{\mathrm{ref}}$ and diameter $d_{\mathrm{ref}}$

    Part-Species1-omegaVHS = 0.24
    Part-Species1-VHSReferenceTemp = 273
    Part-Species1-VHSReferenceDiam = 4.63E-10

It should be noted that although species-specific $\omega$ values can be read-in, DSMC in PICLas should only be utilized with a single $\omega$ at the moment. Other collisional models and their respective parameters are given in Section \ref{sec:dsmc_collision}.

Diatomic molecular species require the definition of the characteristic temperature [K] and their dissociation energy [eV] (which is at the moment only utilized as a first guess for the upper bound of the temperature calculation)

    Part-Species1-CharaTempVib = 4194.9
    Part-Species1-Ediss_eV = 4.53

Polyatomic molecular species require an additional flag, the input of the number of atoms  and whether the molecule is linear (e.g. CO$_2$, $\xi_{\mathrm{rot}} = 2$) or non-linear (e.g. H$_2$O, CH$_4$, $\xi_{\mathrm{rot}} = 3$). The number of the vibrational degrees of freedom is then given by

$$ \alpha = 3 N_{\mathrm{atom}} - 3 - \xi_{\mathrm{rot}} $$

As an example the parameters of CH$_3$ are given below. The molecule has four vibrational modes, with two of them having a degeneracy of two. These values are simply given the according amount of times

    Part-Species1-NumOfAtoms = 4
    Part-Species1-LinearMolec = FALSE
    Part-Species1-CharaTempVib1 = 4320.6
    Part-Species1-CharaTempVib2 = 872.1
    Part-Species1-CharaTempVib3 = 4545.5
    Part-Species1-CharaTempVib4 = 4545.5
    Part-Species1-CharaTempVib5 = 2016.2
    Part-Species1-CharaTempVib6 = 2016.2

These parameters allow the simulation of non-reactive gases. Additional parameters required for the consideration of chemical reaction are given in Section \ref{sec:dsmc_chemistry}.

### Pairing & Collision Modelling \label{sec:dsmc_collision}

WIP: octree, nearest neighbor, VHS

Particles-DSMC-ProhibitDoubleCollision [@Shevyrin2005,@Akhlaghi2018]

### Inelastic Collisions \& Relaxation \label{sec:dsmc_relaxation}

To consider inelastic collisions and relaxation processes within PICLas, the chosen `CollisMode` has to be at least 2

    Particles-DSMC-CollisMode = 2

Two selection procedures are implemented, which differ whether only a single or as many as possible relaxation processes can occur for a collision pair. The default model (`SelectionProcedure = 1`) allows the latter, so-called multi-relaxation method, while `SelectionProcedure = 2` enables the prohibiting double-relaxation method [@Haas1994b]

    Particles-DSMC-SelectionProcedure = 1    ! Multi-relaxation
                                        2    ! Prohibiting double-relaxation

Rotational, vibrational and electronic relaxation (not included by default, see Section \ref{sec:dsmc_electronic_relaxation} for details) processes are implemented in PICLas and their specific options to use either constant relaxation probabilities (default) or variable, mostly temperature dependent, relaxation probabilities are discussed in the following sections. To achieve consistency between continuum and particle-based relaxation modelling, the correction factor of Lumpkin [@Lumpkin1991] can be enabled (default = F):

    Particles-DSMC-useRelaxProbCorrFactor = T

#### Rotational Relaxation \label{sec:dsmc_rotational_relaxation}

To adjust the rotational relaxation this variable has to be changed:

    Particles-DSMC-RotRelaxProb = 0.2   ! Value between 0 and 1 as a constant probability
                                    2   ! Model by Boyd
                                    3   ! Model by Zhang

If `RotRelaxProb` is between 0 and 1, it is set as a constant rotational relaxation probability (default = 0.2). `RotRelaxProb = 2` activates the variable rotational relaxation model according to Boyd [@Boyd1990a]. Consequently, for each molecular species two additional parameters have to be defined, the rotational collision number and the rotational reference temperature. As an example, nitrogen is used [@Boyd1990b].

    Part-Species1-CollNumRotInf = 23.3
    Part-Species1-TempRefRot = 91.5

It is not recommended to use this model with the prohibiting double-relaxation selection procedure (`Particles-DSMC-SelectionProcedure = 2`). Low collision energies result in high relaxation probabilities, which can lead to cumulative collision probabilities greater than 1.

If the relaxation probability is equal to 3, the relaxation model of Zhang et al. [@Zhang2012] is used. However, it is only implemented for nitrogen and not tested. It is not recommended for use.

#### Vibrational Relaxation \label{sec:dsmc_vibrational_relaxation}

Analogous to the rotational relaxation probability, the vibrational relaxation probability is implemented. This variable has to be changed, if the vibrational relaxation probability should be adjusted:

    Particles-DSMC-VibRelaxProb = 0.004 ! Value between 0 and 1 as a constant probability
                                      2 ! Model by Boyd

If `VibRelaxProb` is between 0 and 1, it is used as a constant vibrational relaxation probability (default = 0.004). The variable vibrational relaxation model of Boyd [@Boyd1990b] can be activated with `VibRelaxProb = 2`. For each molecular species pair, the constants A and B according to Millikan and White [@MillikanWhite1963] (which will be used for the calculation of the characteristic velocity and vibrational collision number according to Abe [@Abe1994]) and the vibrational cross section have to be defined. The given example below is a 2 species mixture of nitrogen and oxygen, using the values for A and B given by Farbar [@Farbar2010] and the vibrational cross section given by Boyd [@Boyd1990b]:

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

It is not possible to calculate an instantaneous vibrational relaxation probability with this model [@Boyd1992]. Thus, the probability is calculated for every colission and is averaged. To avoid large errors in cells containing only a few particles, a relaxation of this average probability is implemented. The relaxation factor $\alpha$ can be changed with the following parameter in the ini file:

    Particles-DSMC-alpha = 0.99

The new probability is calculated with the vibrational relaxation probability of the $n^{\mathrm{th}}$ iteration $P^{n}_{\mathrm{v}}$, the number of collision pairs $n_{\mathrm{pair}}$ and the average vibrational relaxation probability of the actual iteration $P^{\mathrm{iter}}_{\mathrm{v}}$.

$$P^{n+1}_{\mathrm{v}}= P^{n}_{\mathrm{v}}  \cdot  \alpha^{2  \cdot  n_{\mathrm{pair}}} + (1-\alpha^{2  \cdot  n_{\mathrm{pair}}}) \cdot P^{\mathrm{iter}}_{\mathrm{v}} $$

This model is extended to more species by calculating a separate probability for each species. An initial vibrational relaxation probability is set by calculating $\mathrm{INT}(1/(1-\alpha))$ vibrational relaxation probabilities for each species and cell by using an instantaneous translational cell temperature.

#### Electronic Relaxation \label{sec:dsmc_electronic_relaxation}

The modelling of electronic states is based on the work of [@Liechty2011a] and uses tabulated energy levels, which can be found in literature for a wide range of species (e.g. for monatomic [@NISTASD], diatomic [@Huber1979], polyatomic [@Herzberg1966] molecules). An example database `DSMCSpecies_electronic_state_full_Data.h5` can be found in e.g. `piclas/regressioncheck/checks/NIG_Reservoir/CHEM_EQUI_TCE_Air_5Spec`, where the energy levels are stored in containers and accessed via the species name, e.g. `Part-Species1-SpeciesName=N2`. Each level is described by its degeneracy in the first column and by the energy in [J] in the seconed column. To include electronic excitation in the simulation, the following parameters are required

    Particles-DSMC-ElectronicModel  = T
    Particles-DSMCElectronicDatabase = DSMCSpecies_electronic_state_full_Data.h5

In case of a large number of electronic levels, their number can be reduced by providing a relative merge tolerance. Levels those relative differences are below this parameter will be merged:

    EpsMergeElectronicState = 1E-3

However, this option should be evaluated carefully based on the specific simulation case and tested against a zero/very low merge tolerance. Finally, the default relaxation probability can be adjusted by

    Particles-DSMC-ElecRelaxProb = 0.01

An electronic state database can be created using a Fortran tool in `piclas/tools/electronic_data`. An alternative is to use the Python-based script discussed in Section \ref{sec:tools_mcc} and to adapt it to electronic energy levels.

### Chemistry & Ionization \label{sec:dsmc_chemistry}

WIP

### Ensuring Physical Simulation Results \label{sec:dsmc_quality}

To determine whether the DSMC related parameters are chosen correctly, so-called quality factors can be written out as part of the regular DSMC state file output by

    Particles-DSMC-CalcQualityFactors = T

This flag writes out the spatial distribution of the mean and maximal collision probability (`DSMC_MeanCollProb` and `DSMC_MaxCollProb`). On the one hand, maximal collision probabilities above unity indicate that the time step should be reduced. On the other hand, very small collision probabilities mean that the time step can be further increased. Additionally, the ratio of the mean collision separation distance to the mean free path is written out (`DSMC_MCSoverMFP`)

$$\frac{l_{\mathrm{mcs}}}{\lambda} < 1$$

The mean collision separation distance is determined during every collision and compared to the mean free path, where its ratio should be less than unity. Values above unity indicate an insufficient particle discretization. In order to estimate the required weighting factor $w$, the following equation can be utilized for a 3D simulation

$$w < \frac{1}{\left(\sqrt{2}\pi d_{\mathrm{ref}}^2 n^{2/3}\right)^3},$$

where $d_{\mathrm{ref}}$ is the reference diameter and $n$ the number density. Here, the largest number density within the simulation domain should be used as the worst-case. For supersonic/hypersonic flows, the conditions behind a normal shock can be utilized as a first guess. For a thruster/nozzle expansion simulation, the chamber or throat conditions are the limiting factor.

## Background Gas

A constant background gas can be utilized to enable efficient particle collisions between the background gas and other particle species (represented by actual simulation particles). The assumption is that the density of the background gas $n_{\mathrm{gas}}$ is much greater than the density of the particle species, e.g. the charged species in a plasma, $n_{\mathrm{charged}}$

$$ n_{\mathrm{gas}} >> n_{\mathrm{charged}}.$$

Under this assumption, collisions within the particle species can be neglected and collisions between the background gas and particle species do not change the conditions of the background gas. It can be activated by defining the species (as defined in Section \ref{sec:dsmc_species}) that should act as the background gas and the number density in m$^{-3}$.

    Particles-DSMCBackgroundGas        = 1
    Particles-DSMCBackgroundGasDensity = 9.64E+21

Other species parameters such as mass, charge, temperature and velocity distribution for the background are defined by the regular read-in parameters

    Part-Species1-SpaceIC              = cuboid
    Part-Species1-velocityDistribution = maxwell_lpn
    Part-Species1-MWTemperatureIC      = 300.0
    Part-Species1-ChargeIC             = 0
    Part-Species1-MassIC               = 6.6464764E-27
    Part-Species1-TempElec             = 300.0

Every time step, particles are generated from the background species and paired with the particle species. Consequently, the collision probabilities are calculated using the conventional DSMC routines and the VHS cross-section model. Aftwards, the collilsion process is performed (if the probability is greater than a random number) and it is tested whether additional energy exchange and chemical reactions occur. While the VHS model is sufficient to model collisions between neutral species, it cannot reproduce the phenomena of a neutral-electron interaction. For this purpose, the cross-section based collision probabilities should be utilized, which are discussed in the following.

### Cross-section based collision probability

For modelling of particle collisions with the Particle-in-Cell method, often the Monte Carlo Collision (MCC) algorithm is utilized. Here, experimentally measured or ab-initio calculated cross-sections are typically utilized to determine the collision probability. In PICLas, the null collision method after [@Birdsall1991],[@Vahedi1995] is implemented, where the number of collision pairs is determined based a maximum collision frequency. Thus, the computational effort is reduced as not every particle has to be checked for a collision, such as in the previously described DSMC-based background gas. To activate the MCC procedure, the collision cross-sections have to be supplied via read-in from a database

    Particles-CollXSec-Database = MCC_Database.h5

An example database, containing the effective collision cross-sections of Argon-electron and Helium-electron, is provided in the tools folder: `piclas/tools/crosssection_database`. Details on how to create an own database with custom cross-section data is given in Section \ref{sec:tools_mcc}. Finally, the input which species should be treated with the MCC model is required

    Part-Species2-SpeciesName = electron
    Part-Species2-UseCollXSec = T

The read-in of the cross-section data is based on the provided species name and the species name of the background gas (e.g. if the background species name is Ar, the code will look for a container named `Ar-electron` in the MCC database). Finally, the cross-section based collision modelling (e.g. for neutral-charged collisions) and the VHS model (e.g. for neutral-neutral collisions) can be utilized within a simulation for different species.

## Modelling of Continuum Gas Flows \label{sec:continuum}

Two methods are currently implemented to allow the simulation of gas flows in the continuum and transitional regime, where the DSMC method is computationally too expensive. The Fokker–Planck- and Bhatnagar-Gross-Krook-based approximation of the collision integral are compared in detail in paper to be published in Physics of Fluids. It is recommended to utilize a previous DSMC parameter file to ensure a complete simulation setup.

### Fokker–Planck Collision Operator \label{sec:fpflow}

The implementation of the FP-based collision operator is based on the publications by [@Gorji2014] and [@Pfeiffer2017]. The collision integral is hereby approximated by a drift and diffusion process

$$  \left.\frac{\partial f}{\partial t}\right|_{\mathrm{coll}}\approx-\sum_{i=1}^3 {\frac{\partial }{\partial v_i}(A_i f)+\frac{1}{2}\sum_{i=1}^3 \sum_{j=1}^3\frac{\partial ^2 }{\partial v_i\partial v_j}(D_{ij}f)}, $$

where $\mathbf{A}$ is the drift vector and $\mathcal{D}$ the diffusion matrix.

The current implementation supports:

- 2 different methods: Cubic (only atomic species) and Ellipsoidal Statistical (ES)
- Single species, monoatomic and polyatomic gases
- Thermal non-equilibrium with rotational and vibrational excitation (continuous or quantized treatment)
- 2D/Axisymmetric simulations
- Variable time step (adaption of the distribution according to the maximal relaxation factor and linear scaling)

Relevant publications of the developers:

- Implementation of the cubic Fokker-Planck in PICLas [@Pfeiffer2017]
- Comparison of the cubic and ellipsoidal statistical Fokker-Planck [@Jun2019]
- Simulation of a nozzle expansion (including the pressure chamber) with ESBGK, ESFP and coupled ESBGK-DSMC, comparison to experimental measurements [@Pfeiffer2019a]

To enable the simulation with the FP module, the respective compiler setting has to be activated:

    PICLAS_TIMEDISCMETHOD = FP-Flow

A parameter file and species initialization file is required, analogous to the DSMC setup. To enable the simulation with the FP methods, select the Fokker-Planck method, cubic (`=1`) and ES (`=2`):

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

### Bhatnagar-Gross-Krook Collision Operator \label{sec:bgk}

The implementation of the BGK-based collision operator is based on the publications by [@Pfeiffer2018a] and [@Pfeiffer2018b]. It allows the simulation of gas flows in the continuum and transitional regime, where the DSMC method is computationally too expensive. The collision integral is hereby approximated by a relaxation process:

$$ \left.\frac{\partial f}{\partial t}\right|_{\mathrm{coll}} \approx \nu(f^t-f), $$

where $f^t$ is the target distribution function and $\nu$ the relaxation frequency.

The current implementation supports:

- 4 different methods (i.e. different target distribution functions): Ellipsoidal Statistical, Shakov, standard BGK, and Unified
- Single species, monoatomic and polyatomic gases
- Thermal non-equilibrium with rotational and vibrational excitation (continuous or quantized treatment)
- 2D/Axisymmetric simulations
- Variable time step (adaption of the distribution according to the maximal relaxation factor and linear scaling)

Relevant publications of the developers:

- Implementation and comparison of the ESBGK, SBGK, and Unified models in PICLas for atomic species [@Pfeiffer2018a]
- Extension of the modelling to diatomic species including quantized vibrational energy treatment, validation of ESBGK with the Mach 20 hypersonic flow measurements of the heat flux on a $70^\circ$ cone [@Pfeiffer2018b]
- Simulation of a nozzle expansion (including the pressure chamber) with ESBGK, SBGK and coupled ESBGK-DSMC, comparison to experimental measurements [@Pfeiffer2019a],[@Pfeiffer2019b]
- Extension to polyatomic molecules, simulation of the carbon dioxide hypersonic flow around a flat-faced cylinder, comparison of ESBGK, SBGK and DSMC regarding the shock structure and heat flux [@Pfeiffer2019c]

To enable the simulation with the BGK module, the respective compiler setting has to be activated:

    PICLAS_TIMEDISCMETHOD = BGK-Flow

A parameter file and species initialization file is required, analogous to the DSMC setup. To enable the simulation with the BGK methods, select the BGK method, ES (`=1`), Shakov (`=2`), Standard BGK (`=3`), and Unified (`=4`):

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

## Macroscopic Bodies

In order to enable macroscopic bodies inside of computational domain, an appropriate time step
method has to be compiled first. A standalone executable is created by enabling the following parameter with cmake

    PICLAS_TIMEDISCMETHOD = MACROBODY

This time step method combines the DSMC method timestep with additionally simulating Macroscopic
Bodies. Here, the octree used for DSMC and sampling of macroscopic values considers macroscopic bodies that
may occupy some elements and the appropriate volume portions are calculated.

### Starting simulations with macroscopic bodies

In order to start a smulation, the number of macroscopic bodies is enabled by defining

    MacroBody-nMacroBody=1

and setting it greater than 0. For each defined macroscopic body properties are assigned with (Spheres)

    MacroBody1-center=(/ 10.0 , 0.1 , 0.5 /)
    MacroBody1-radius=0.5e-5

For surface interactions of DSMC particles with the bodies the appropriate accomodation coefficients and temperature can be defined

    MacroBody1-momentumAcc=1.
    MacroBody1-transAcc=1.
    MacroBody1-temp=300

If the the mesh cell sizes are by more than one order of magnitude bigger than the smallest
macroscopic body then the volume estimation must be controlled and adapted.
The volume portions are calculated with a Monte Carlo approach with nPointsMCVolumeEstimate number
of point per octree-subcell.
with:

    Particles-nPointsMCVolumeEstimate=1000.


|   octree level | Total MC points per Element     |
| -------------: | :------------------------------ |
|     lvl=0      |                  1000           |
|     lvl=1      |                  8000           |
|     lvl=2      |                  56000          |

### Macroscopic bodies, which move trough domain
Note: Currently, restart and nProcs>1 only works for non-moving and non-size-changing Macro bodies.
To let the macroscopic bodies move (MacroBody1-velocity $> 0$) inside the computational domain, enable

    MacroBody1-velocity=(/-1165.,0.,0./)

In order to enable influences of particle collisions on macroscopic bodies velocity (acceleration, brownian motion), the following parameter needs to be enabled

    MacroBody-AccelerationEnabled=T

For influences of fluxes onto the macroscopic bodies' temperature, size or mass activate:

    MacroBody-FluxesEnabled=T

### Output of MacroBodies

Information of macroscopic bodies is written into the state file (VarNamesmacroParticles,MacroPartData) at every analyzedt.
In order to add additional elemdata output (for debug/analyze purposes) of ElemHasMacroBody of the first defined macroscopic body, activate

    MacroBody-WriteElemData=T
