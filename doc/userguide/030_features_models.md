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

|      Option       | Values |                          Notes                          |
| :---------------: | :----: | :-----------------------------------------------------: |
| CartesianPeriodic |  T/F   | If a fully periodic box (all 6 sides) is computed, the  |
|                   |        | intersections do not have to be computed. Instead, each |
|                   |        | particle can be simply shifted by the periodic vector.  |
|   FastPeriodic    |  T/F   | Moves particle the whole periodic distance once, which  |
|                   |        |  can be several times the mesh size in this direction.  |


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
|   BezierNewtonAngle   |  <PI/2   | Angle to switch between Clipping and a Newton algorithm. |
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

## Boundary Conditions

To-do: Modification of boundaries with the PICLas parameter file (order is of importance)

### Field

Dielectric -> type 100?

### Particle

Within the parameter file it is possible to define different particle boundary conditions. The number of boundaries is defined by

    Part-nBounds=2
    Part-Boundary1-SourceName=BC_OPEN
    Part-Boundary1-Condition=open
    Part-Boundary2-SourceName=BC_WALL
    Part-Boundary2-Condition=reflective

The `Part-Boundary1-SourceName=` corresponds to the name given during the preprocessing step with HOPR. The available conditions (`Part-Boundary1-Condition=`) are described in the table below.

|  Condition   | Description                                                                                                                                                                                 |
| :----------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|    `open`    | Every particle crossing the boundary will be deleted.                                                                                                                                       |
| `reflective` | Allows the definition of specular and diffuse reflection. A perfect specular reflection is performed, if no other parameters are given (discussed in more detail in the following section). |
| `symmetric`  | A perfect specular reflection, without sampling of particle impacts.                                                                                                                        |

#### Diffuse Wall

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

#### Porous Wall / Pump

The porous boundary condition uses a removal probability to determine whether a particle is deleted or reflected at the boundary. The main application of the implemented condition is to model a pump, according to [@Lei2017]. It is defined by giving the number of porous boundaries and the respective boundary number (`BC=2` corresponds to the `BC_WALL` boundary defined in the previous section) on which the porous condition is. 

    Part-nPorousBC=1
    Part-PorousBC1-BC=2
    Part-PorousBC1-Pressure=5.
    Part-PorousBC1-Temperature=300.
    Part-PorousBC1-PumpingSpeed=2e-9
    Part-PorousBC1-DeltaPumpingSpeed-Kp=0.1
    Part-PorousBC1-DeltaPumpingSpeed-Ki=0.0

The removal probability is determined through the given pressure [Pa] and temperature [K] at the boundary. A pumping speed can be given as a first guess, however, the pumping speed $S$ [$m^3/s$] will be adapted if the proportional factor ($K_\mathrm{p}$, `DeltaPumpingSpeed-Kp`) is greater than zero

$$ S^{n+1}(t) = S^{n}(t) + K_\mathrm{p} \Delta p(t) + K_\mathrm{i} \int_0^t \Delta p(t') dt',$$

where $\Delta p$ is the pressure difference between the given pressure and the actual pressure at the pump. An integral factor ($K_\mathrm{i}$, `DeltaPumpingSpeed-Ki`) can be utilized to mimic a PI controller. The proportional and integral factors are relative to the given pressure. However, the integral factor has not yet been thoroughly tested. The removal probability $\alpha$ is then calculated by

$$\alpha = \frac{S n \Delta t}{N_\mathrm{pump} w} $$

where $n$ is the sampled, cell-local number density and $N_\mathrm{pump}$ is the total number of impinged particle at the pump during the previous time step. $\Delta t$ is the time step and $w$ the weighting factor. The pumping speed $S$ is only adapted if the resulting removal probability $\alpha$ is between zero and unity. The removal probability is not species-specific.

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

Using the regions, multiple pumps can be defined on a single boundary.

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

The adaptive particle emission can be combined with the circular inflow feature. In this context when the area of the actual emission circle/ring is very small, it is preferable to utilize the `Type=4` constant mass flow condition. `Type=3` assumes an open boundary and accounts for particles leaving the domain through that boundary already when determining the number of particles to be inserted. As a result, this method tends to overpredict the given mass flow, when the emission area is very small and large sample size would required to have enough particles that leave the method. For `Type=4` method, the actual number of particles leaving the domain through the circular inflow is counted, and thus the correct mass flow can be reproduced.

It should be noted that while multiple adaptive boundaries are possible, adjacent boundaries that share a mesh element should be avoided or treated carefully.

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

Depending on the utilized collision model, different parameters have to be defined. As an example, the parameters for the Variable Hard Sphere (VHS) collision cross-section model are be defined by the temperature exponent $\omega$, reference temperature $T_\mathrm{ref}$ and diameter $d_\mathrm{ref}$

    Part-Species1-omegaVHS = 0.24
    Part-Species1-VHSReferenceTemp = 273
    Part-Species1-VHSReferenceDiam = 4.63E-10

It should be noted that although species-specific $\omega$ values can be read-in, DSMC in PICLas should only be utilized with a single $\omega$ at the moment. Other collisional models and their respective parameters are given in Section \ref{sec:dsmc_collision}.

Diatomic molecular species require the definition of the characteristic temperature [K] and their dissociation energy [eV] (which is at the moment only utilized as a first guess for the upper bound of the temperature calculation)

    Part-Species1-CharaTempVib = 4194.9
    Part-Species1-Ediss_eV = 4.53

Polyatomic molecular species require an additional flag, the input of the number of atoms  and whether the molecule is linear (e.g. CO$_2$, $\xi_\mathrm{rot} = 2$) or non-linear (e.g. H$_2$O, CH$_4$, $\xi_\mathrm{rot} = 3$). The number of the vibrational degrees of freedom is then given by

$$ \alpha = 3 N_\mathrm{atom} - 3 - \xi_\mathrm{rot} $$

As an example the parameters of CH$_3$ are given below. The molecule has four vibrational modes, with two of them having a degeneracy of two. These values are simply given the according amount of times

    Part-Species2-NumOfAtoms = 4
    Part-Species2-LinearMolec = FALSE
    Part-Species2-CharaTempVib1 = 4320.6
    Part-Species2-CharaTempVib2 = 872.1
    Part-Species2-CharaTempVib3 = 4545.5
    Part-Species2-CharaTempVib4 = 4545.5
    Part-Species2-CharaTempVib5 = 2016.2
    Part-Species2-CharaTempVib6 = 2016.2

These parameters allow the simulation of non-reactive gases. Additional parameters required for the consideration of chemical reaction are given in Section \ref{sec:dsmc_chemistry}.

### Pairing & Collision Modelling \label{sec:dsmc_collision}

WIP: octree, VHS

### Relaxation \label{sec:dsmc_relaxation}

WIP

### Chemistry & Ionization \label{sec:dsmc_chemistry}

WIP

### Ensuring Physical Simulation Results \label{sec:dsmc_quality}

To determine whether the DSMC related parameters are chosen correctly, so-called quality factors can be written out as part of the regular DSMC state file output by

    Particles-DSMC-CalcQualityFactors = T

This flag writes out the spatial distribution of the mean and maximal collision probability (`DSMC_MeanCollProb` and `DSMC_MaxCollProb`). On the one hand, maximal collision probabilities above unity indicate that the time step should be reduced. On the other hand, very small collision probabilities mean that the time step can be further increased. Additionally, the ratio of the mean collision separation distance to the mean free path is written out (`DSMC_MCSoverMFP`)

$$\frac{l_\mathrm{mcs}}{\lambda} < 1$$

The mean collision separation distance is determined during every collision and compared to the mean free path, where its ratio should be less than unity. Values above unity indicate an insufficient particle discretization. In order to estimate the required weighting factor $w$, the following equation can be utilized for a 3D simulation

$$w < \frac{1}{\left(\sqrt{2}\pi d_\mathrm{ref}^2 n^{2/3}\right)^3},$$

where $d_\mathrm{ref}$ is the reference diameter and $n$ the number density. Here, the largest number density within the simulation domain should be used as the worst-case. For supersonic/hypersonic flows, the conditions behind a normal shock can be utilized as a first guess. For a thruster/nozzle expansion simulation, the chamber or throat conditions are the limiting factor.

## Surface Chemistry

WIP

## Modelling of Continuum Gas Flows

Two methods are currently implemented to allow the simulation of gas flows in the continuum and transitional regime, where the DSMC method is computationally too expensive. The Fokker–Planck- and Bhatnagar-Gross-Krook-based approximation of the collision integral are compared in detail in paper to be published in Physics of Fluids. It is recommended to utilize a previous DSMC parameter file to ensure a complete simulation setup.

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

## Output of Macroscopic Variables

In general, simulation results are either available spatially resolved based on the mesh (classic CFD results) and/or as integral values (e.g. for reservoir/heat bath simulations).

### Field Variables

WIP

### Flowfield and Surface Variables

A sampling over a certain number of iterations is performed to calculate the macroscopic values such as number density, bulk velocity and temperature from the microscopic particle information. Output and sampling on surfaces can be enabled by

    Particles-DSMC-CalcSurfaceVal = T

Parameters indicating the quality of the simulation (e.g. the maximal collision probability in case of DSMC) can be enabled by

    Particles-DSMC-CalcQualityFactors = T

Two variants are available in PICLas, allowing to sample a certain amount of the simulation duration or to sample continuously during the simulation and output the result after the given number of iterations.

The first variant is usually utilized to sample at the end of a simulation, when the steady condition is reached. The first parameter `Part-TimeFracForSampling` defines the percentage that shall be sampled relative to the simulation end time $T_\mathrm{end}$ (Parameter: `TEnd`)

    Part-TimeFracForSampling = 0.1
    Particles-NumberForDSMCOutputs = 2

`Particles-NumberForDSMCOutputs` defines the number of outputs during the sampling time. Example: The simulation end time is $T_\mathrm{end}=1$, thus sampling will begin at $T=0.9$ and the first output will be written at $T=0.95$. At this point the sample will NOT be resetted but continued. Therefore, the second and last output at $T=T_\mathrm{end}=1.0$ is not independent of the previous result but contains the sample of the complete sampling duration. It should be noted that if a simulation is continued at e.g. $T=0.95$, sampling with the given parameters will begin immediately.

The second variant can be used to produce outputs for unsteady simulations, while still to be able to sample for a number of iterations (Parameter: `Part-IterationForMacroVal`). The first two flags allow to enable the output of flowfield/volume and surface values, respectively.

    Part-WriteMacroVolumeValues = T
    Part-WriteMacroSurfaceValues = T
    Part-IterationForMacroVal = 100

Example: The simulation end time is $T_\mathrm{end}=1$ with a time step of $\Delta t = 0.001$. With the parameters given above, we would sample for 100 iterations up to $T = 0.1$ and get the first output. Afterwards, the sample is deleted and the sampling begins anew for the following output at $T=0.2$. This procedure is repeated until the simulation end, resulting in 10 outputs with independent samples.

### Integral Variables

PartAnalyze/FieldAnalyze
