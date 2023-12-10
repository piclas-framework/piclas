(sec:Radiation)=
# Radiation 

## Radiation Coupling

To account for high temperature gas radiation, the DSMC method can be coupled with a radiation module {cite}`Beyer2022JQSRT,Beyer2022RGD,Pfeiffer2024ICARUS`. For a deeper physical knowledge of the implemented options see Ref. {cite}`Beyer2022JQSRT`.
The radiation module calculates cell-local emission and absorption coefficients on the same computational mesh as the flow field solver using a line-by-line (LBL) method, and the radiative transfer solver simulates the radiative energy transfer through the computational domain using a photon Monte Carlo (pMC) solver. The LBL solver can also be used in a stand-alone version. The enitre radiation simulation can be run highly parallel using a MPI-shared memory concept (MPI 3.0). The LBL solver uses a domain decomposition and the pMC solver distributes the photon bundles across the different cores.
The radiation module can be enabled by compiling PICLas with the following parameter

    PICLAS_TIMEDISCMETHOD = Radiation

In addition, the following parameter must be set for larger computational meshes

    PICLAS_INKIND8 = ON

For the radiation module, the following options are available and must be set accordingly

    Radiation-RadType = 1   ! Full radiation module with radiative energy transfer
                        2   ! Blackbody radiation in entire computational domain
                        3   ! Emission and absorption coefficients only (stand-alone radiation solver)
                        4   ! Shock tube mode

To define the species which are used in the radiation simulation, the following parameters must be set in one single ini-file (as an example, atomic nitrogen is chosen)

    Part-Species1-MassIC                = 2.32600E-26 ! N Molecular Mass
    Part-Species1-MacroParticleFactor   = 2.5E10

    Part-Species1-SpeciesName           = N
    Part-Species1-InteractionID         = 1
    Part-Species1-Tref                  = 273     ! in K
    Part-Species1-dref                  = 3.0E-10 ! in m
    Part-Species1-omega                 = 0.24

    Part-Species1-RadiationIonizationEn = 117345  ! Ionization Energy, cm-1
    Part-Species1-RadiationRadius_A     = 0.7     ! Radius, A
    Part-Species1-Starkex               = 0.33    ! Stark Index
    Part-Species1-NuclCharge            = 1       ! Charge (1:neutral particles, 2:ions)
    Radiation-Species1-SpectraFileName  = Ni.dat

The `Radiation-Species[$]-SpectraFileName` contains information about level energy and radiative transition lines. An example is given in the regression-tests and can e.g. be built using information available in the NIST database. For convenience, the calculation of each species can be disabled with

    Part-Species[$]-DoRadiation = F

In a first step of the tool chain, the local emission and absorption coefficients are determined on the same computational mesh used by the flow field solver. Therefore, the same mesh has to be read in. If the emission in a single cell is to be calculated (`Radiation-RadType = 1`), a mesh file containing one single cell must be provided.

Different radiation mechanisms can be considered, depending on the energy state of the involved electron

    Radiation-bb-atoms     = T   ! atomic line radiation (bound-bound)
    Radiation-bb-molecules = T   ! molecular band radiation (bound-bound)
    Radiation-bf           = T   ! recombination radiation (bound-free)
    Radiation-ff           = T   ! Bremsstrahlung (free-free)

If the complete radiation module (1) is selected, the flow field information (DSMCState) previously calculated with PICLas is used as input by setting the following parameters

    Radiation-MacroRadInput = T
    Radiation-MacroInput-Filename = PROJECTNAME_DSMCState_000.001000000.h5

Often, radiation solvers use the translational temperature of the electrons as $T_\mathrm{E}$ to determine upper state densities. However, PICLas also offers the possibility to use the actual electronic excitation temperature

    Radiation-UseElectronicExcitation = T

For a single cell (`Radiation-RadType = 1`), temperatures and densities can be given by

    Part-Species[$]-RadiationTtrans  = 10000.
    Part-Species[$]-RadiationTelec   = 10000.
    Part-Species[$]-RadiationTvib    = 10000.
    Part-Species[$]-RadiationTrot    = 10000.
    Part-Species[$]-RadiationNumDens = 1E20.

The units are Kelvin and m^-3. For electrons, they are red in by

    Radiation-TElectrons             = 10000.
    Radiation-NumDensElectrons       = 1E20.

The wavelength range and the discretization for the simulation can be set with the following parameters. Radiative transitions outside of this range will not be considered in the simulation!

    Radiation-MinWaveLen   = 200.    ! minimum wavelength in nanometers
    Radiation-MaxWaveLen   = 1000.   ! maximum wavelength in nanometers
    Radiation-WaveLenDiscr = 500000  ! Number of spectral discretization points

To save computational memory, wavelengths can also be spectrally binned together.

    Radiation-WaveLenReductionFactor       = 10

<!-- Additionally, a time step needs to be set with $ManualTimeStep \ge t_\mathrm{end}$. -->

For `Radiation-RadType = 3`, the next step is to determine the radiative energy transfer through the computational domain (change in radiation intensity due to the emission and absorption of the surrounding gas). This is done with a photon Monte Carlo method, where photon bundles are traced through the computational domain by solving the radiative transfer equation. Due to the time scales involved, only the steady-state solution is used, furthermore, scattering effects within the gas are neglected. For shock tubes (`Radiation-RadType = 4`), a simplified version of the radiative energy transfer along a tanget slab is calculated. The required diameter of the shock tube can be set (in meters) with

    Radiation-ShockTubeDiameter = 0.16

For the radiative energy transfer of `Radiation-RadType = 3`, the number of photon bundles in each computational cell is set by

    Radiation-NumPhotonsPerCell = 200

Instead of placing the same number of photons in each computational cell and giving them different energies, it is also possible to give them the same energy and redistribute them across the computational domain according to the emitted energy.

    Radiation-AdaptivePhotonNumEmission = T   ! true:photons have the same energy, false:number of photons per cell is equal

The initial properties of the photon bundles are set randomly. However, if a more uniform distribution within a cell is desired, different options are available. For the position, they are placed in the cell using a 2,3-Halton sequence. Additionally, their directions can be distributed along a spiral configuration

    Radiation-PhotonPosModel = 2 ! 1:random 2:Halton
    Radiation-DirectionModel = 2 ! 1:random 2:spiral

The absorption along their paths can be calculated (1) analytically or (2) stochastically, however, it is strongly recommended to do it analytically

    Radiation-AbsorptionModel = 1

To determine the wavelength of each photon bundle, two different methods are implemented, (1) an acceptance-rejection method and (2) a bisection method

    Radiation-PhotonWaveLengthModel = 1   ! 1:Acceptance-Rejection 2:Bisection

The acceptance-rejection method is more computationally intensive than the bisection method. However, it should be more accurate, because it consideres individual wavelengths rather than integral values. Thus, the acceptance-rejection method can never select a wavelength that has no emission, while the bisection method has a very low probability of doing so.

For surface properties for photons on a reflecting wall, different options are possible. They can either be specularly reflected 

    Part-Boundary3-Condition=reflective
    Part-Boundary3-PhotonSpecularReflection = T

or diffuse. If they are reflected diffusely, the energy accommodation coefficient can bet set ([0,1]) depending on the surface properties.

    Part-Boundary3-Condition=reflective
    Part-Boundary3-PhotonSpecularReflection = F
    Part-Boundary3-PhotonEnACC = 0.5


In addition, PICLas offers the possibility to sample spectra on a virtual sensor and to compare them with measured data. Two different options are available: (1) with a viewing angle and (2) along a line-of-sight (LOS)

    Radiation-RadObservationPointMethod    = 2   !1:observation angle, 2:LOS

The location of the virtual sensor can be defined with

    Radiation-ObservationMidPoint          = (/-1.0,1E-5,0.0/)

and the view direction with

    Radiation-ObservationViewDirection     = (/1.0,0.0,0.0/)

Its diameter can be set with

    Radiation-ObservationDiameter          = 0.1

and the viewing angle with 

    Radiation-ObservationAngularAperture   = 0.404533

Along a LOS, it is also possible to use as many photons per cell as wavelength discretizations. These photons then receive the actual energy of the corresponding wavelength in order to obtain a perfect sampling of the corresponding energies for spectral comparisons.

    Radiation-ObservationCalcFullSpectra   = T

To simulate with a high resolution and to match the units of the radiance, the output can be expressed in different units, e.g. to have a spectral discretization of 1/A and to get the radiance in /nm, the following parameters must be set

    Radiation-WaveLenReductionFactorOutput = 10

To account for instrumental broadening, the radiance profile can be mathematically convolved with a trapezoid

    Radiation-ObservationDoConvolution = T

The trapezoid can be defined by a topwidth and a basewidth, both read-in in angstroms

    Radiation-ObservationSlitFunction      = (/1.7,3.42/)

The simulations can also be run on a two-dimensional rotationally symmetric mesh. To do this, the following options must be set. Different tracking routines are used than with an axisymmetric particle solver, therefore, the RadialWeighting-PartScaleFactor can have different values

    Particles-Symmetry2D                         = T
    Particles-Symmetry2DAxisymmetric             = T
    Particles-RadialWeighting                    = T
    Particles-RadialWeighting-PartScaleFactor    = 10000
    Particles-RadialWeighting-CloneMode          = 2
    Particles-RadialWeighting-CloneDelay         = 6
    Particles-RadialWeighting-CellLocalWeighting = F


## Raytracing

In addition to the radiation coupling, a ray tracing model is implemented. A boundary must be defined from which rays or photons are emitted in a preliminary step, which are tracked throughout the domain until they are absorbed at a boundary. The volumes and surface elements are sampled by passing photons and from this information the ionization within each volume element and secondary electron emission from each surface is calculated in the actual plasma simulation. The output of the sampling procedure can be viewed as irradiated volumes and surfaces and is written to the output files EUV_RadiationVolState.h5 and EUV_RadiationSurfState.h5, which can be converted to .vtk format with piclas2vtk.
Raytracing is activated with

    UseRayTracing = T

It only requires only one single section in the parameter.ini file. The user must specify a single rectangular and planar particle-boundary (here with index 5)

    RayTracing-PartBound = 5

from which the number of rays
    
    RayTracing-NumRays = 200000000

shall be emitted in the direction

    RayTracing-RayDirection = (/0.0, 0.0, -1.0/)

Currently, all coordinates of this boundary must have the same z-coordinate and it must extend into the complete domain into the x- and y-direction. The parameter

    RayTracing-PulseDuration = 1e-9

defines the pulse duration $\tau$ that defines the temporal shape of the light intensity function $I\propto\exp(-(t/\tau)^2)$ in [s].

    RayTracing-NbrOfPulses

defines the number of pulses that are performed and 

    RayTracing-WaistRadius

the waist radius $w_{b}$ that defines the spatial intensity via $I\propto\exp(-(r/w_b)^2)$ in [m].
The wavelength in [m] is given by 

    RayTracing-WaveLength = 50e-9

and the repetition rate of the pulses in [Hz] by

    RayTracing-RepetitionRate = 2e3

The time-averaged (over one pulse) power density of the pulse in [W/m2] is used in 

    RayTracing-PowerDensity = 1e3

which is converted to an average pulse energy considering the irradiated area, hence, the same parameter can be used for the quarter and the full mesh setups.

To account for the reflectivity of specific surfaces, the absorption rate $A_{\nu}=1-R$ ($R$ is the reflectivity) for photons must be supplied for each particle boundary. This is done by setting the parameter

    Part-Boundary1-PhotonEnACC = 1.0
    
to a value between zero and unity. Additionally, it is possibly to switch between perfect angular reflection and diffuse reflection for each boundary.

    Part-Boundary$-PhotonSpecularReflection = T

The parameter

    RayTracing-ForceAbsorption=T

activates sampling of photons on surfaces independent of what happens to them there. They might be reflected or absorbed. If this parameter is set to `false``, then only absorbed photons will be sampled on surfaces. By also sampling reflected photons, the statistic is improved, hence, it should always be activated.

The angle under which photons are emitted from the particle-boundary is calculated from the normal vector of the boundary and the parameter

    RayTracing-RayDirection

Because the interaction of every ray with every volume and surface element in three dimensions would lead to an unfeasible amount of memory usage if stored on the hard drive, the calculated volume and surface intersections need to be agglomerated in such a way that the details of the irradiated geometry are preserved. One goal is to have a clean cut between shadowed and illuminated regions where this interface cuts through surface of volume elements and without the need for a cut-cell method that splits the actual mesh elements. This is achieved by introducing a super-sampling method in the volume as well as on the surface elements. For the volumetric sampling, the originally cell-constant value is distributed among the volume sub-cells depending on a element-specific number of sub cells $n_{cells} = (N_{cell} + 1)^3$ , where
$N_{cell}$ is the polynomial degree in each element used for visualization of the super-sampled
ray tracing solution. The polynomial degree $N_{cell}$ is chosen between unity and a user-defined
value, which can be automatically selected depending on the different criteria

    RayTracing-VolRefineMode = 0 ! 0: do nothing (default)
                                 ! 1: refine below user-defined z-coordinate with NMax
                                 ! 2: scale N according to the mesh element volume between NMin>=1 and NMax>=2
                                 ! 3: refine below user-defined z-coordinate and scale N according to the mesh element volume between NMin>=1 and NMax>=2 (consider only elements below the user-defined z-coordinate for the scaling)

The maximum polynomial degree within refined volume elements for photon tracking (p-adaption) can hereby be set using

    RayTracing-NMax = 1

In contrast to the volume super-sampling, only one global parameter is used to refine the all surfaces for sampling. Each surface can be split into $n^2$ sub-surfaces on which the sampling is performed via the parameter

    RayTracing-nSurfSample = 2

The surfaces (quadrilaterals) are therefore equidistantly divided at the midpoint of each edge to create approximately equal-sized sub-surfaces.