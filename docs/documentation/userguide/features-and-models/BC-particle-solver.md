(sec:particle-boundary-conditions)=
# Boundary Conditions - Particle Solver

Within the parameter file it is possible to define different particle boundary conditions. The number of boundaries is defined by

    Part-nBounds = 2
    Part-Boundary1-SourceName   = BC_OPEN
    Part-Boundary1-Condition    = open
    Part-Boundary2-SourceName   = BC_WALL
    Part-Boundary2-Condition    = reflective

The `Part-Boundary1-SourceName=` corresponds to the name given during the preprocessing step with PyHOPE. The available conditions
(`Part-Boundary1-Condition=`) are described in the table below.

|         Condition          | Description                                                                                                                                                                    |
| :------------------------: | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|           `open`           | Every particle crossing the boundary will be deleted                                                                                                                           |
|        `symmetric`         | A perfect specular reflection, without sampling of particle impacts                                                                                                            |
|      `symmetric_axis`      | Definition of the axis of rotation in axisymmetric 2D simulations: Section {ref}`sec:2D-axisymmetric`                                                                          |
|      `symmetric_dim`       | Definition of symmetrical boundaries in 1D and 2D simulations: Section {ref}`sec:2D-axisymmetric`, Section {ref}`sec:1D-sym`                                                   |
|        `reflective`        | Definition of different surface models: Section {ref}`sec:particle-boundary-conditions-reflective`, Section {ref}`sec:surface-chemistry`, Section {ref}`sec:catalytic-surface` |
|       `rot_periodic`       | Definition of rotational periodicity: Section {ref}`sec:particle-boundary-conditions-rotBC`                                                                                    |
| `rot_periodic_inter_plane` | Extension of rotational periodicity, allowing non-conformal interfaces and varying periodicity                                                                                 |

(sec:particle-boundary-conditions-reflective)=
## Reflective Wall

A reflective boundary can be defined with

    Part-Boundary2-SourceName   = BC_WALL
    Part-Boundary2-Condition    = reflective

A perfect specular reflection is performed, if no other parameters are given. Gas-surface interactions can be modelled with the
extended Maxwellian model {cite}`Padilla2009`, using accommodation coefficients of the form

$$\alpha = \frac{E_i-E_r}{E_i - E_w}$$

where $i$, $r$ and $w$ denote the incident, reflected and wall energy, respectively.  The coefficient `MomentumACC` is utilized to
decide whether a diffuse (`MomentumACC` $>R$) or specular reflection (`MomentumACC` $<R$) occurs upon particle impact, where
$R=[0,1)$ is a random number. Separate accommodation coefficients can be defined for the translation (`TransACC`), rotational
(`RotACC`), vibrational (`VibACC`) and electronic energy (`ElecACC`) accommodation at a constant wall temperature [K].

    Part-Boundary2-MomentumACC  = 1.
    Part-Boundary2-WallTemp     = 300.
    Part-Boundary2-TransACC     = 1.
    Part-Boundary2-VibACC       = 1.
    Part-Boundary2-RotACC       = 1.
    Part-Boundary2-ElecACC      = 1.

An additional option `Part-Boundary2-SurfaceModel` is available, that is used for heterogeneous reactions (reactions that have reactants
in two or more phases) or secondary electron emission models. These models are described in detail in Section {ref}`sec:surface-chemistry`.

(sec:particle-boundary-conditions-reflective-wallvelo)=
### Wall movement (Linear & rotational)

Additionally, a linear wall velocity [m/s] can be given

    Part-Boundary2-WallVelo = (/0,0,100/)

In the case of rotating walls the `-RotVelo` flag, a rotation frequency [Hz], and the rotation axis (x=1, y=2, z=3) must be set.
Note that the definition of the rotational direction is defined by the sign of the frequency using the right-hand rule.

    Part-Boundary2-RotVelo = T
    Part-Boundary2-RotFreq = 100
    Part-Boundary2-RotAxis = 3

The wall velocity will then be superimposed onto the particle velocity.

### Linear temperature gradient

A linear temperature gradient across a boundary can be defined by supplying a second wall temperature and the start and end vector
as well as an optional direction to which the gradient shall be limited (default: 0, x = 1, y = 2, z = 3)

    Part-Boundary2-WallTemp2      = 500.
    Part-Boundary2-TempGradStart  = (/0.,0.,0./)
    Part-Boundary2-TempGradEnd    = (/1.,0.,1./)
    Part-Boundary2-TempGradDir    = 0

In the default case of the `TempGradDir = 0`, the temperature will be interpolated between the start and end vector, where the
start vector corresponds to the first wall temperature `WallTemp`, and the end vector to the second wall temperature `WallTemp2`.
Position values (which are projected onto the temperature gradient vector) beyond the gradient vector utilize the first (Start)
and second temperature (End) as the constant wall temperature, respectively. In the special case of `TempGradDir = 1/2/3`, the
temperature gradient will only be applied along the chosen the direction. As opposed to the default case, the positions of the
surfaces are not projected onto the gradient vector before checking whether they are inside the box spanned by `TempGradStart` and
`TempGradEnd`. The applied surface temperature is output in the `DSMCSurfState` as `Wall_Temperature` for verification.

### Radiative equilibrium

Another option is to adapt the wall temperature based on the heat flux assuming that the wall is in radiative equilibrium.
The temperature is then calculated from

$$ q_w = \varepsilon \sigma T_w^4,$$

where $\varepsilon$ is the radiative emissivity of the wall (default = 1) and
$\sigma = \pu{5.67E-8 Wm^{-2}K^{-4}}$ is the Stefan-Boltzmann constant. The adaptive boundary is enabled by

    Part-AdaptWallTemp = T
    Part-Boundary1-UseAdaptedWallTemp = T
    Part-Boundary1-RadiativeEmissivity = 0.8

If provided, the wall temperature will be adapted during the next output of macroscopic variables, where the heat flux calculated
during the preceding sampling period is utilized to determine the side-local temperature. The temperature is included in the
`State` file and thus available during a restart of the simulation. The surface output (in `DSMCSurfState`) will additionally
include the temperature distribution in the `Wall_Temperature` variable (see Section {ref}`sec:sampled-flow-field-and-surface-variables`).
To continue the simulation without further adapting the temperature, the first flag has to be disabled (`Part-AdaptWallTemp = F`).
It should be noted that the adaptation should be performed multiple times to achieve a converged temperature distribution.

(sec:particle-boundary-conditions-rotBC)=
## Rotational Periodicity

The rotational periodic boundary condition can be used in order to reduce the computational effort in case of an existing
rotational periodicity. In contrast to symmetric boundary conditions, a macroscopic flow velocity in azimuthal direction can be
simulated (e.g. circular flow around a rotating cylinder). Exactly two corresponding boundaries must be defined by setting
`rot_periodic` as the BC condition and the rotating angle for each BC. Multiple pairs of boundary conditions with different angles
can be defined.

    Part-Boundary1-SourceName       = BC_Rot_Peri_plus
    Part-Boundary1-Condition        = rot_periodic
    Part-Boundary1-RotPeriodicAngle = 90.

    Part-Boundary2-SourceName       = BC_Rot_Peri_minus
    Part-Boundary2-Condition        = rot_periodic
    Part-Boundary2-RotPeriodicAngle = -90.

CAUTION! The correct sign for the rotating angle must be determined. The position of particles that cross one rotational
periodic boundary is transformed according to this angle, which is defined by the right-hand rule and the rotation axis:

    Part-RotPeriodicAxi = 1    ! (x = 1, y = 2, z = 3)

The usage of rotational periodic boundary conditions is limited to cases, where the rotational periodic axis is one of the three
Cartesian coordinate axis (x, y, z) with its origin at (0, 0, 0).

### Intermediate Plane Definition
If several segments with different rotation angles are defined, exactly two corresponding boundary conditions must be defined for each segment.
In pyHope, multiple internal boundaries must be allowed by setting

    CheckInternalBoundaries = F

Since the plane between these segments with different rotational symmetry angles represents a non-conforming connection, additional
two boundary conditions must be defined as `rot_periodic_inter_plane` at this intermediate plane. Both boundary conditions must refer to each other in the
definition in order to ensure the connection.

    Part-Boundary40-SourceName       = BC_INT_R1_BOT
    Part-Boundary40-Condition        = rot_periodic_inter_plane
    Part-Boundary40-AssociatedPlane  = 41

    Part-Boundary41-SourceName       = BC_INT_S1_TOP
    Part-Boundary41-Condition        = rot_periodic_inter_plane
    Part-Boundary41-AssociatedPlane  = 40

Note that using the intermediate plane definition with two corresponding boundary conditions allows the user to mesh the segments independently,
creating a non-conforming interface at the intermediate plane. However, use of these non-conformal grids is so far only possible in standalone DSMC simulations.

## Porous Wall / Pump

The porous boundary condition uses a removal probability to determine whether a particle is deleted or reflected at the boundary.
The main application of the implemented condition is to model a pump, according to {cite}`Lei2017`. It is defined by giving the
number of porous boundaries and the respective boundary number (`BC=2` corresponds to the `BC_WALL` boundary defined in the
previous section) on which the porous condition is.

    Surf-nPorousBC=1
    Surf-PorousBC1-BC=2
    Surf-PorousBC1-Type=pump
    Surf-PorousBC1-Pressure=5.
    Surf-PorousBC1-PumpingSpeed=2e-9
    Surf-PorousBC1-DeltaPumpingSpeed-Kp=0.1
    Surf-PorousBC1-DeltaPumpingSpeed-Ki=0.0

Currently, two porous BC types are available, `pump` and `sensor`. For the former, the removal probability is determined through
the given pressure [Pa] at the boundary. A pumping speed can be given as a first guess, however, the pumping speed $S$ [$m^3/s$]
will be adapted if the proportional factor ($K_{\mathrm{p}}$, `DeltaPumpingSpeed-Kp`) is greater than zero

$$ S^{n+1}(t) = S^{n}(t) + K_{\mathrm{p}} \Delta p(t) + K_{\mathrm{i}} \int_0^t \Delta p(t') dt',$$

where $\Delta p$ is the pressure difference between the given pressure and the actual pressure at the pump. An integral factor
($K_{\mathrm{i}}$, `DeltaPumpingSpeed-Ki`) can be utilized to mimic a PI controller. The proportional and integral factors are
relative to the given pressure. However, the integral factor has not yet been thoroughly tested. The removal probability $\alpha$
is then calculated by

$$\alpha = \frac{S n \Delta t}{N_{\mathrm{pump}} w} $$

where $n$ is the sampled, cell-local number density and $N_{\mathrm{pump}}$ is the total number of impinged particle at the pump
during the previous time step. $\Delta t$ is the time step and $w$ the weighting factor. The pumping speed $S$ is only adapted if
the resulting removal probability $\alpha$ is between zero and unity. The removal probability is not species-specific.

To reduce the influence of statistical fluctuations, the relevant macroscopic values (pressure difference $\Delta p$ and number
density $n$) can be sampled for $N$ iterations by defining (for all porous boundaries)

    AdaptiveBC-SamplingIteration=10

A porous region on the specified boundary can be defined. At the moment, only the `circular` option is implemented. The origin of
the circle/ring on the surface and the radius have to be given. In the case of a ring, a maximal and minimal radius is required
(`-rmax` and `-rmin`, respectively), whereas for a circle only the input of maximal radius is sufficient.

    Surf-PorousBC1-Region=circular
    Surf-PorousBC1-normalDir=1
    Surf-PorousBC1-origin=(/5e-6,5e-6/)
    Surf-PorousBC1-rmax=2.5e-6

The absolute coordinates are defined as follows for the respective normal direction.

| Normal Direction | Coordinates |
| :--------------: | :---------: |
|      x (=1)      |    (y,z)    |
|      y (=2)      |    (z,x)    |
|      z (=3)      |    (x,y)    |

Using the regions, multiple pumps can be defined on a single boundary. Additionally, the BC can be used as a sensor by defining
the respective type:

    Surf-PorousBC1-BC=3
    Surf-PorousBC1-Type=sensor

Together with a region definition, a pump as well as a sensor can be defined on a single and/or multiple boundaries, allowing e.g.
to determine the pressure difference between the pump and a remote area of interest.

(sec:surface-chemistry)=
## Surface Chemistry

Modelling of reactive surfaces is enabled by setting `Part-BoundaryX-Condition=reflective` and an
appropriate particle boundary surface model `Part-BoundaryX-SurfaceModel`:

    Part-Boundary1-SurfaceModel = 0

The available conditions (`Part-BoundaryX-SurfaceModel=`) are described in the table below, ranging from simple empirical models and secondary electron/ion emission to finite-rate catalysis modelling including a surface treatment.

|    Model    | Description                                                                                                                                                                                  |
| :---------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 0 (default) | Standard extended Maxwellian scattering                                                                                                                                                      |
|      1      | Empirical modelling of sticking coefficient/probability                                                                                                                                      |
|      2      | Fixed probability surface chemistry                                                                                                                                                          |
|      3      | Secondary electron emission as a square-fit ($a + b*E + c*E^2$) for electron energies above the work function W                                                                              |
|      4      | Secondary electron emission as a power-fit ($a E^b + c$) for electron energies above the work function W                                                                                     |
|      5      | Secondary electron emission as given by Ref. {cite}`Levko2015`.                                                                                                                              |
|      7      | Secondary electron emission due to ion impact (SEE-I with $Ar^{+}$ on different metals) as used in Ref. {cite}`Pflug2014` and given by Ref. {cite}`Depla2009` with a default yield of 13 \%. |
|      8      | Secondary electron emission due to ion impact (SEE-E with $e^{-}$ on dielectric surfaces) as used in Ref. {cite}`Liu2010` and given by Ref. {cite}`Morozov2004`.                             |
|      9      | Secondary electron emission due to ion impact (SEE-I with $Ar^{+}$) with a constant yield of 1 \%. Emitted electrons have an energy of 6.8 eV upon emission.                                 |
|     10      | Secondary electron emission due to ion impact (SEE-I with $Ar^{+}$ on copper) as used in Ref. {cite}`Theis2021` originating from {cite}`Phelps1999`                                          |
|     11      | Secondary electron emission due to electron impact (SEE-E with $e^{-}$ on quartz (SiO$_{2}$)) as described in Ref. {cite}`Zeng2020` originating from {cite}`Dunaevsky2003`                   |
|     12      | Secondary electron emission due to electron impact as described in Ref. {cite}`Seiler1983`                                                                                                   |
|     13      | Secondary electron emission due to electron impact according to the Vaughan formula described in Ref. {cite}`Villemant2019`                                                                                                   |
|     20      | Finite-rate catalysis model, Section {ref}`sec:catalytic-surface`                                                                                                                            |

### Empirical model for a sticking coefficient

To model the sticking of gas particles on cold surfaces, an empirical model is available, which is based on experimental measurements. The sticking coefficient is modelled through the product of a non-bounce probability $B(\alpha)$ and a condensation probability $C(\alpha,T)$

$$ p_s (\alpha,T) = B(\alpha) C(\alpha,T) $$

The non-bounce probability introduces a linear dependency on the impact angle $\alpha$

$$
B(\alpha) = \begin{cases}
   1 , & |\alpha| < \alpha_{\mathrm{B}} \\
   \dfrac{90^{\circ}-|\alpha|}{90^{\circ}-|\alpha_{\mathrm{B}}|} , & \alpha_{\mathrm{B}} \leq |\alpha| \leq 90^{\circ} \\
\end{cases}
$$

$\alpha_{\mathrm{B}}$ is a model-dependent cut-off angle. The condensation probability introduces a linear dependency on the surface temperature $T$

$$
C(\alpha, T) = \begin{cases}
   1 , & T < T_1 \\
   \dfrac{T_2(\alpha)-T}{T_2(\alpha)-T_1(\alpha)} , & T_1 \leq T \leq T_2 \\
   0 , & T > T_2 \\
\end{cases}
$$

The temperature limits $T_1$ and $T_2$ are model parameters and can be given for different impact angle ranges defined by the maximum impact angle $\alpha_{\mathrm{max}}$. These model parameters are read-in through the species database and have to be provided in the `/Surface-Chemistry/StickingCoefficient` dataset in the following format (example values):

| $\alpha_{\mathrm{max}}$ [deg] | $\alpha_{\mathrm{B}}$ [deg] | $T_1$ [K] | $T_2$ [K] |
| ----------------------------: | --------------------------: | --------: | --------: |
|                            45 |                          80 |        50 |       100 |
|                            90 |                          70 |        20 |        50 |

In this example, within impact angles of $0°\leq\alpha\leq45°$, the model parameters of the first row will be used and for $45°<\alpha\leq90°$ the second row. The number of rows is not limited. The species database is read-in by

    Particles-Species-Database = SpeciesDatabase.h5

As additional output, the cell-local sticking coefficient will be added to the sampled surface output. A particle sticking to the surface will be deleted and its energy added to the heat flux sampling. This model can be combined with the linear temperature gradient and radiative equilibrium modelling as described in Section {ref}`sec:particle-boundary-conditions-reflective`.

### Fixed probability surface chemistry

This simple fixed-probability surface chemistry model allows the user to define arbitrary surface reactions, by defining the
impacting species, the products and a fixed event probability. The reaction is then assigned to the boundaries by specifying their number and index as defined previously.
This model corresponds to `Part-BoundaryX-SurfaceModel = 2`, which is set automatically when a reaction of this type is defined.

    Surface-NumOfReactions               = 1
    Surface-Reaction1-Type               = P
    Surface-Reaction1-Reactants          = (/1,0/)
    Surface-Reaction1-Products           = (/2,1,0/)
    Surface-Reaction1-EventProbability   = 0.25
    Surface-Reaction1-NumOfBoundaries    = 2
    Surface-Reaction1-Boundaries         = (/1,3/)

Optionally, a reaction-specific accommodation coefficient for the products can be defined, otherwise the surface-specific accommodation will be utilized for the product species:

    Surface-Reaction1-ProductAccommodation = 0.

In the case that the defined event does not occur, a regular interaction using the surface-specific accommodation coefficients is performed. Examples are provided as part of the regression tests: `regressioncheck/NIG_DSMC/SURF_PROB_DifferentProbs` and `regressioncheck/NIG_DSMC/SURF_PROB_MultiReac`.

(sec:BC-see)=
### Secondary Electron Emission (SEE)

Different models are implemented for secondary electron emission that are based on either electron or ion bombardment, depending on
the surface material. All models require the specification of the electron species that is emitted from the surface via

    Part-SpeciesA-PartBoundB-ResultSpec = C

where electrons of species `C` are emitted from boundary `B` on the impact of species `A`. Some of the models allow the choice of an angle and energy distribution function to define the velocity vector of the secondary. The available options are:

|                  Name | Description                                                                                                                                                                                     |                  Source                  |
| --------------------: | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------------: |
|     deltadistribution | Random velocity vector and complete remaining impact energy                                                                                                                                     |                    -                     |
|        uniform-energy | Random velocity vector and random uniform distribution of the remaining impact energy                                                                                                           |                    -                     |
| Chung-Everhart-cosine | Angle distribution according to $\cos$ in the normal direction, equally distributed in the tangential direction, and a Chung-Everhart distribution of the energy $f = \frac{E}{(E+W)^4}$        | {cite}`Chung1974`, {cite}`Greenwood2002` |
|                cosine | Angle distribution according to $\cos$ in the normal direction, equally distributed in the tangential direction                                                                                 |                                          |

For the `Chung-Everhart-cosine` distribution, in the case of 2 or more secondaries, we are currently sampling each energy independently, which can result in an energy
addition and thus energy conservation violation. An output to monitor the percentage of violations and energy addition as a percentage of the impact energy per SEE event can be enabled through `CalcEnergyViolationSEE = T`.

For a simulation using variable particle weights (`Part-vMPF = T`) as described in Section {ref}`sec:split-merge`, the models 3, 4, 12, and 13 support the emission of only a single secondary, weighted according to the calculated yield. This feature can be enabled per boundary:

    Part-Boundary1-SurfMod-vMPF = T

Additional parameters, which can be utilized per boundary for models 3, 4, 12, and 13 are

    Part-Boundary1-SurfMod-SubtractWorkFunction = .FALSE.       ! Default: .TRUE.
    Part-Boundary1-SurfMod-ReflectElectron      = .TRUE.        ! Default: .FALSE.

The first parameter allows to disable the subtraction of the work function from the incident electron energy (example in `piclas/regressioncheck/NIG_PIC_poisson_Leapfrog/BC_SEE_EnergyDistribution_Constant`). The second parameter can be used to utilize the yield function as a probability for a reflection if the incident energy is lower than the work function (example in `piclas/regressioncheck/NIG_PIC_poisson_Leapfrog/BC_SEE_SquareFit_ReflectBelowThreshold`).

#### Model 3/4

These models use a square- (= 3) or power-fit (= 4) and an additional threshold to model the secondary electron emission yield.
It is assumed that the impacting particle is absorbed.

$$\text{Square-fit: }\gamma = (a E + b E^2 + c)H(E-W),$$
$$\text{Power-fit: }\gamma = (a E^b + c)H(E-W),$$

where $a$, $b$, $c$ are the fitting coefficients and $W$ is the material-dependent work function [eV] above which the yield is calculated.
The parameters are read-in through:

    Part-BoundaryB-SurfModSEEFitCoeff   = (/0.1,0.5,0.25,9/)      ! (/a,b,c,W/)

Additionally, the energy distribution can be selected with

    Part-BoundaryB-SurfModEnergyDistribution = Chung-Everhart-cosine

It should be noted that per default the impact energy is reduced by the work function before the energy distribution. An example of the model usage is given in the regression test: `piclas/regressioncheck/NIG_PIC_poisson_Leapfrog/BC_SEE_PowerFit/`. Fit coefficients can be found for example in {cite}`Goebel2008`.

#### Model 5

The model by Levko {cite}`Levko2015` can be applied for copper electrodes for electron and ion bombardment and is activated via
`Part-BoundaryX-SurfaceModel=5`. For ions, a fixed emission yield of 0.02 is used and for electrons an energy-dependent function is
employed.

#### Model 7

The model by Depla {cite}`Depla2009` can be used for various metal surfaces and features a default emission yield of 13 \% and is
activated via `Part-BoundaryX-SurfaceModel=7` and is intended for the impact of $Ar^{+}$ ions. For more details, see the original
publication.

The emission yield and energy can be varied for this model by setting

    SurfModEmissionYield  = 1.45 ! ratio of emitted electron flux vs. impacting ion flux [-]
    SurfModEmissionEnergy = 6.8  ! [eV]

respectively.
The emission yield represents the ratio of emitted electrons vs. impacting ions and the emission energy is given in electronvolt.
If the energy is not set, the emitted electron will have the same velocity as the impacting ion.

Additionally, a uniform energy distribution function for the emitted electrons can be set via

    SurfModEnergyDistribution = uniform-energy

which will scale the energy of the emitted electron to fit a uniform distribution function.

#### Model 8

The model by Morozov {cite}`Morozov2004` can be applied for dielectric surfaces and is activated via
`Part-BoundaryX-SurfaceModel=8` and has an additional parameter for setting the reference electron temperature (see model for
details) via `Part-SurfaceModel-SEE-Te`, which takes the electron temperature in Kelvin as input (default is 50 eV, which
corresponds to 11604 K).
The emission yield is determined from an energy-dependent function.
The model can be switched to an automatic determination of the bulk electron temperature via

    Part-SurfaceModel-SEE-Te-automatic = T ! Activate automatic bulk temperature calculation
    Part-SurfaceModel-SEE-Te-Spec      = 2 ! Species ID used for automatic temperature calculation (must correspond to electrons)

where the species ID must be supplied, which corresponds to the electron species for which, during `Part-AnalyzeStep`, the global
translational temperature is determined and subsequently used to adjust the energy dependence of the SEE model. The global (bulk)
electron temperature is written to *PartAnalyze.csv* as *XXX-BulkElectronTemp-[K]*.

#### Model 10

An energy-dependent model of secondary electron emission due to $Ar^{+}$ ion impact on a copper cathode as used in
Ref. {cite}`Theis2021` originating from {cite}`Phelps1999` is
activated via `Part-BoundaryX-SurfaceModel=10`. For more details, see the original publications.

#### Model 11

An energy-dependent model (linear and power fit of measured SEE yields) of secondary electron emission due to $e^{-}$ impact on a
quartz (SiO$_{2}$) surface as described in Ref. {cite}`Zeng2020` originating from {cite}`Dunaevsky2003` is
activated via `Part-BoundaryX-SurfaceModel=11`. For more details, see the original publications.

#### Model 12

This model relies on a semi-empirical formulation by Seiler {cite}`Seiler1983`. It is assumed that the impacting particle is absorbed.

$$\gamma = a \cdot 1.11 \cdot \left(\frac{E}{b}\right)^{-0.35}\left(1-e^{-2.3\left(\frac{E}{b}\right)^{1.35}}\right)$$

where $a$ and $b$ [eV] are material-specific coefficients and $W$ is the work function [eV] above which the yield is calculated.
The parameters are read-in through:

    Part-BoundaryB-SurfModSEEFitCoeff   = (/1.0,700,0.0,9/)      ! (/a,b,c,W/)

Additionally, the energy distribution can be selected with

    Part-BoundaryB-SurfModEnergyDistribution = Chung-Everhart-cosine

It should be noted that per default the impact energy is reduced by the work function before the energy distribution. An example of the model usage is given in the regression test: `piclas/regressioncheck/NIG_DSMC/BC_SEE_Model_12/`.

#### Model 13

This model relies on the Vaughan formula given by Villeman {cite}`Villemant2019`. It is assumed that the impaction particle is absorbed.

$$\gamma = a \left( \frac{E}{b} \cdot e^{1-\frac{E}{b}} \right)^c$$

where $a$, $b$ [eV], and $c$ are material-specific coefficients and $W$ is the work function [eV] above which the yield is calculated.
The parameters are read-in through:

    Part-BoundaryB-SurfModSEEFitCoeff   = (/2.016,299,0.563,0/)      ! (/a,b,c,W/)

Additionally, the energy distribution can be selected with

    Part-BoundaryB-SurfModEnergyDistribution = cosine

Using the cosine energy distribution, the angle distribution is according to $\cos$ in the normal direction and equally distributed in the tangential direction. Using the SEE model 13 (and the cosine energy distribution), the energy of all secondary emitted electrons is set to 2 eV.

If a work function greater than zero is set, the impact energy is reduced by the work function before the energy distribution per default. An example of the model usage is given in the regression test: `piclas/regressioncheck/NIG_DSMC/BC_SEE_Model_13/`.

(sec:catalytic-surface)=
## Catalytic Surfaces

Catalytic reactions can be modeled in PICLas using a finite-rate reaction model with an implicit treatment of the reactive surface. For a better resolution of the parameters, the catalytic boundaries are discretized into a certain number of subsides. A definition of the boundary temperature in the parameter input file is required in all cases. Different types of surfaces can be defined by the lattice constant of the unit cell `Part-BoundaryX-LatticeVec` and the number of particles in the unit cell `Part-BoundaryX-NbrOfMol-UnitCell`. These parameters are used in the calculation of the number of active sites.

By default, the simulation is started with a clean surface, but an initial species-specific coverage can be specified by `Part-BoundaryX-SpeciesX-Coverage`, which represents the relative number of active sites that are occupied by adsorbate particles. Maximum values for the coverage values can be specified by:

    Part-Boundary1-Species1-MaxCoverage
    Part-Boundary1-MaxTotalCoverage

Multi-layer adsorption is enabled by a maximal total coverage greater than 1.

The reaction paths are defined in the input parameter file. First, the number of gas-surface reactions to be read in must be defined:

    Surface-NumOfReactions = 2

A catalytic reaction and the boundary on which it takes place is then defined by

    Surface-Reaction1-SurfName           = Adsorption
    Surface-Reaction1-Type               = A
    Surface-Reaction1-Reactants          = (/1,0/)
    Surface-Reaction1-Products           = (/2,1,0/)
    Surface-Reaction1-NumOfBoundaries    = 2
    Surface-Reaction1-Boundaries         = (/1,3/)

All reactants and products are defined by their respective species index. In the case of multiple reacting, the order does not influence the input. The following options are available for the catalytic reaction type:

| Model | Description                                                 |
| ----: | ----------------------------------------------------------- |
|     A | Adsorption: Kisliuk or Langmuir model                       |
|     D | Desorption: Polanyi-Wigner model                            |
|    ER | Eley-Rideal reaction: Arrhenius based chemistry             |
|    LH | Langmuir-Hinshelwood reaction: Arrhenius based chemistry    |
|   LHD | Langmuir-Hisnhelwood reaction with instantaneous desorption |

For the treatment of multiple reaction paths of the same species, a possible bias in the reaction rate is avoided by a randomized treatment. Bulk species can participate in the reaction. In this case, the bulk species is defined by `Surface-Species` and the corresponding species index. All reaction types allow for the definition of a reaction enthalpy. In addition, this value can be linearly increased (negative factor) or decreased (positive factor) by a scaling factor for the heat of reaction. Both values are given in [K].

    Surface-Reaction1-ReactHeat      = 17101.4
    Surface-Reaction1-HeatScaling    = 1202.9

Depending on the reaction type, different additional parameters have to be defined. More details on the specific cases are given in the following subsections. Example input files for the adsorption and desorption of CO and O2 on a palladium surface can be found in the regression tests `regressioncheck/WEK_DSMC/SurfChem_AdsorpDesorp_CO` and `regressioncheck/WEK_DSMC/SurfChem_AdsorpDesorp_O2`.

### Adsorption

For the modelling of the adsorption of a gas particle on the surface, two models are available: the simple Langmuir model, with a linear dependence of the adsorption probability on the surface coverage, and the precursor-based Kisliuk model:

$$ S = S_0 (1 + K (1/\theta^{\alpha} - 1))^{-1}$$

Here, $S_0$ is the binding coefficient for a clean surface, $\alpha$ is the dissociation constant (2 for dissociative adsorption) and $K$ is the equilibrium constant between adsorption and desorption from the precursor state. For $K = 1$, the model simplifies to the Langmuir case. The parameters can be defined in PICLas as follows:

    Surface-Reaction1-StickingCoefficient  = 0.2
    Surface-Reaction1-DissOrder            = 1
    Surface-Reaction1-EqConstant           = 0.6

A special case of adsorption is the dissociative adsorption (`Surface-ReactionX-DissociativeAdsorption = true`), where only half of the molecule binds to the surface, while the other half remains in the gas phase. The adsorbate half `Surface-ReactionX-AdsorptionProduct` and the gas phase product `Surface-ReactionX-GasPhaseProduct` are specified by their respective species indices. The adsorption probability is calculated analogously to the general case.

Lateral interactions between multiple adsorbate species, which can disfavor further adsorption can be taken into account by the command `Surface-ReactionX-Inhibition` and the species index of the inhibiting species.

### Desorption

The desorption of an adsorbate particle into the gas phase is modelled by the Polanyi-Wigner equation.

$$k(T) = A T^b \theta^{\alpha}_{A} e^{-E_\mathrm{a}/T}$$

where $A$ is the prefactor ([1/s, m$^2$/s] depending on the dissociation constant), $\alpha$ the dissociation constant and $E_\mathrm{a}$ the
activation energy [K]. These parameters can be defined in PICLas as follows:

    Surface-ReactionX-Prefactor
    Surface-ReactionX-Energy

### Catalytic Reaction

Catalytic reactions can be modelled in PICLas using a finite-rate reaction model with an implicit treatment of
the reactive surface. The adsorbate is not represented by particles but by a coverage value per species, which
is stored for every sub-surface element of a catalytic boundary. Gas phase species are treated as usual, so a
reaction either consumes an impacting particle, inserts new particles into the gas phase, or both.

Two mechanisms drive the chemistry:

* **Impact-driven reactions** (`A`, `ER`) are evaluated whenever a particle hits a catalytic boundary. They use
  the state of the sub-surface element that was hit.
* **Time-step-driven reactions** (`D`, `LH`, `LHD`) are evaluated once per time step for every sub-surface
  element of a catalytic boundary, independently of any particle impact.

A definition of the boundary temperature is required in all cases, since it enters every reaction rate.

#### Surface properties of the boundary

The number of active sites per area follows from the lattice constant of the unit cell and the number of
particles per unit cell:

    Part-Boundary1-LatticeVector    = 0.389E-9
    Part-Boundary1-NbrOfMol-UnitCell = 1

If `LatticeVector` is left at zero, a generic monolayer with $10^{19}$ sites per m$^2$ is assumed instead. Note
that the number of sites of a sub-surface element is the site density times its area, so a fine surface
subdivision on a fine mesh can leave only a few sites per sub-element. In that case the coverage becomes
strongly quantised, because a single adsorbing simulation particle already changes it by
`MacroParticleFactor` divided by the number of sites. It is worth checking this ratio before interpreting
coverage results.

By default the simulation starts with a clean surface. A species-specific initial coverage, the relative number
of occupied active sites, can be prescribed together with the maximum values:

    Part-Boundary1-Species1-Coverage    = 0.1
    Part-Boundary1-Species1-MaxCoverage = 0.333
    Part-Boundary1-MaxTotalCoverage     = 1.0

Multi-layer adsorption is enabled by a maximum total coverage greater than 1. The initial coverage must not
exceed the species maximum, and the sum of all initial coverages must not exceed the total maximum; both are
checked during the read-in.

Choosing an initial coverage exactly equal to `MaxCoverage` is a degenerate case: the free-site fraction is then
zero and the species can never adsorb until another reaction frees sites. This is legal but rarely intended.

#### Definition of a reaction

First the total number of gas-surface reactions is declared:

    Surface-NumOfReactions = 2

A reaction, its type and the boundaries it acts on are then defined by:

    Surface-Reaction1-SurfName        = Adsorption_CO
    Surface-Reaction1-Type            = A
    Surface-Reaction1-Reactants       = (/1,0/)
    Surface-Reaction1-Products        = (/1,0,0/)
    Surface-Reaction1-NumOfBoundaries = 2
    Surface-Reaction1-Boundaries      = (/1,3/)

All reactants and products are given by their species index. `Reactants` holds up to two, `Products` up to
three entries, unused slots are zero. A reaction is only evaluated on the boundaries listed in `Boundaries`.

| Model | Description                                                  |
| ----: | ------------------------------------------------------------ |
|     A | Adsorption: Kisliuk or Langmuir model                        |
|     D | Desorption: Polanyi-Wigner model                             |
|    ER | Eley-Rideal reaction: Arrhenius based chemistry              |
|    LH | Langmuir-Hinshelwood reaction: Arrhenius based chemistry     |
|   LHD | Langmuir-Hinshelwood reaction with instantaneous desorption  |

#### Meaning of reactants and products per type

The species indices are interpreted differently depending on the reaction type. This is the most common source
of input errors, so it is listed explicitly:

| Type | `Reactants`                                            | `Products`                                                     |
| ---: | ------------------------------------------------------ | -------------------------------------------------------------- |
|    A | impacting **gas phase** species                        | slot 1: adsorbate, slot 2: gas phase fragment, slot 3: unused   |
|    D | **adsorbed** species that desorbs                      | gas phase products                                             |
|   ER | impacting gas species **and** the consumed adsorbate   | gas phase products                                             |
|   LH | **adsorbed** species                                   | adsorbed products, they stay on the surface                     |
|  LHD | **adsorbed** species                                   | gas phase products, inserted directly                           |

For `D`, `LH` and `LHD` the entry in `Reactants` is therefore not a gas phase collision partner but the
adsorbate whose coverage drives the rate and which is consumed by the reaction. It cannot be derived from the
products: for an associative desorption $2\,\mathrm{O(ads)} \rightarrow \mathrm{O_2(gas)}$ the product carries
no information about the consumed adsorbate. The entry is mandatory.

Note that an adsorbate is generally a **separate species index** from its gas phase counterpart, so that it can
carry its own `MaxCoverage` and its own coverage entry. Sharing one index between the gas species and the
adsorbate is possible but is a modelling decision, not an assumption of the code.

#### Reaction enthalpy and energy accommodation

All reaction types allow the definition of a reaction enthalpy, which can be linearly scaled with the coverage.
Both values are given in [K] and converted internally:

    Surface-Reaction1-ReactHeat   = 17101.4
    Surface-Reaction1-HeatScaling = 1202.9

The released energy follows $E(\theta) = E_\mathrm{React} - \theta\,E_\mathrm{Scaling}$, i.e. a positive scaling
factor reduces the enthalpy at high coverage. This reproduces the weaker binding of an adsorbate when
neighbouring sites are already occupied. A scaling factor larger than the enthalpy itself would turn the
released energy negative above a certain coverage, which is unphysical.

The fraction of that energy which is transferred to the solid is set by

    Surface-Reaction1-EnergyAccommodation = 1.0

with values between 0 and 1, default 1. The remaining fraction stays with the desorbing product. Note that this
coefficient currently serves a second purpose as the translational accommodation of the product in the
post-reaction velocity sampling, so the two are not independent.

#### Restricting a reaction to a coverage or temperature window

Both windows are optional and switched on separately. They are checked before the rate is evaluated; a reaction
outside its window is skipped entirely.

A **coverage window** limits the reaction to a range of the total coverage and, in addition, to a range of every
individual species coverage:

    Surface-Reaction1-CoverageDependence          = TRUE
    Surface-Reaction1-MinimumTotalCoverage        = 0.1
    Surface-Reaction1-MaximumTotalCoverage        = 0.8
    Surface-Reaction1-Species1-MinimumCoverage    = 0.0
    Surface-Reaction1-Species1-MaximumCoverage    = 0.5
    Surface-Reaction1-Species2-MinimumCoverage    = 0.2
    Surface-Reaction1-Species2-MaximumCoverage    = 1.0

The species-specific bounds are read for **all** species when `CoverageDependence` is enabled, defaulting to
`0.0` and `1.0`. All of them have to be satisfied simultaneously. Typical uses are a reaction that requires a
co-adsorbate to be present, or one that is blocked once a poisoning species accumulates.

A **temperature window** restricts the reaction to a range of the wall temperature:

    Surface-Reaction1-TemperatureDependence = TRUE
    Surface-Reaction1-MinimumTemperature    = 300.
    Surface-Reaction1-MaximumTemperature    = 1200.

Defaults are `0.` and `10000.`, so an enabled dependence without further input has no effect. This is intended
for rate expressions that were fitted over a limited temperature range and should not be extrapolated.

#### Adsorption reactions

Two models are available for the adsorption of a gas particle: the Langmuir model with a linear dependence of
the adsorption probability on the coverage, and the precursor-based Kisliuk model:

$$ S = S_0 \left(1 + K \left(1/\theta_\mathrm{free}^{\alpha} - 1\right)\right)^{-1} $$

Here $S_0$ is the sticking coefficient of the clean surface, $K$ the equilibrium constant between adsorption and
desorption from the precursor state, and $\alpha$ the number of adjacent free sites required by the process.
For $K = 1$ the model reduces to the Langmuir case, where $S = S_0\,\theta_\mathrm{free}^{\alpha}$.

$\theta_\mathrm{free}$ is the **available** free-site fraction,
and is clamped to $[0,1]$, so that $S_0$ retains its meaning as the sticking coefficient of the clean surface.
The number of required sites enters through the exponent, the availability through the base.

    Surface-Reaction1-StickingCoefficient = 0.2
    Surface-Reaction1-EqConstant          = 0.6
    Surface-Reaction1-DissOrder           = 1

Adsorption is additionally suppressed when the resulting total coverage would exceed `MaxTotalCoverage`.

##### Dissociative adsorption

Dissociative adsorption, where one fragment binds to the surface while the other returns to the gas phase, is
expressed through the product slots. No separate switch is required:

| `Products` | Meaning                                                                    |
| ---------- | -------------------------------------------------------------------------- |
| `0,0,0`    | the impacting species itself is adsorbed                                   |
| `X,0,0`    | `X` is adsorbed                                                            |
| `X,Y,0`    | `X` is adsorbed, `Y` is released into the gas phase (dissociative)          |

The detection depends only on the number of occupied slots, never on the identity of a species, so an adsorbate
sharing the index of its gas phase counterpart is handled correctly. A non-zero second slot requires an explicit
first slot, since the split of the molecule between surface and gas phase would otherwise be undefined. A third
entry is rejected for adsorption reactions.

`DissOrder` defaults to **2 for dissociative** and **1 for non-dissociative** adsorption, matching the two
adjacent sites a dissociating molecule needs. An explicit value overrides this.

Dissociative adsorption in which *both* fragments remain on the surface is not covered by this convention.

##### Inhibition and promotion by co-adsorbates

Lateral interactions with other adsorbates modify the free-site fraction. Inhibitors block sites and reduce it,
promotors enhance the adsorption and increase it:

    Surface-Reaction1-Inhibition = TRUE
    Surface-Reaction1-Inhibitors = (/3,0,0/)
    Surface-Reaction1-Promotion  = TRUE
    Surface-Reaction1-Promotors  = (/4,5,0/)

Up to three species can be given for each. Every contribution enters normalised by the maximum coverage of that
species. Since the free-site fraction is clamped to $[0,1]$, a strong promotion cannot push the sticking
coefficient above $S_0$.

#### Desorption reactions

Desorption into the gas phase is modelled by the Polanyi-Wigner equation:

$$ k(T) = \nu\,\left(\theta_A N_\mathrm{s}\right)^{\alpha}\,e^{-E_\mathrm{a}/T} $$

with the prefactor $\nu$, the site density $N_\mathrm{s}$, the desorption order $\alpha$ and the activation
energy $E_\mathrm{a}$ in [K]:

    Surface-Reaction3-Prefactor = 1E13
    Surface-Reaction3-Energy    = 17688.8
    Surface-Reaction3-DissOrder = 1

`DissOrder` is the desorption order, and only the values 1 and 2 are supported. A value of 2 denotes associative
desorption, $2\,A_\mathrm{(ads)} \rightarrow A_2$, in which case two adsorbates are consumed per event and a
separate prefactor conversion is applied. The default is 1.

The activation energy can depend linearly on the coverage through the lateral interaction parameter,
$E_\mathrm{a} = E_0 + \theta\,W$, which lowers the barrier at high coverage for a negative $W$:

    Surface-Reaction3-LateralInteraction = -18410.8

If `Prefactor` is left at zero, the prefactor is instead computed from a coverage-dependent correlation,
$\nu = 10^{(C_a + C_b \theta)}$:

    Surface-Reaction3-Ca = 16
    Surface-Reaction3-Cb = -15

The number of desorption events of a time step follows from the rate, the time step, the area of the
sub-surface element and an exponentially distributed random variate, limited by the available adsorbate.
Fractional events are accumulated across time steps and released once they add up to a full simulation
particle, so the desorption rate is reproduced correctly even when much less than one particle per time step
desorbs.

#### Catalytic reactions

The Eley-Rideal and the Langmuir-Hinshelwood reactions use Arrhenius-type rates together with the coverage of
all surface-bound reactants. The activation energy and the prefactor are read as for the desorption:

    Surface-Reaction2-Prefactor = 1E-16
    Surface-Reaction2-Energy    = 7246.38

The prefactor is a bimolecular rate coefficient in [m$^3$/s] for the Eley-Rideal reaction and in [m$^2$/s] for
the Langmuir-Hinshelwood case. For a Langmuir-Hinshelwood reaction with two adsorbed reactants the coverage of
each of them enters the rate, so the units of the prefactor depend on the number of reactants.

##### Eley-Rideal

An Eley-Rideal `ER` reaction is a reaction which an incoming gas phase particle reacts directly with an already
adsorbed species upon impact, forming a product that is released into the gas phase. Since the reaction happens
in a single collision without the impacting particle first equilibrating with the surface, the activation barrier
is overcome by the translational energy the particle brings along rather than by the thermal energy of the wall.

##### Langmuir-Hinshelwood

For the type `LH` the product species stays adsorbed on the surface until a separate desorption reaction takes
place. The product coverage is limited by its `MaxCoverage`. For reactions combined with very high desorption
rates the type `LHD` is more appropriate: the products are inserted directly into the gas phase without an
intermediate desorption step.


#### Diffusion

Two variants of an instantaneous diffusion are available, both equivalent to an averaging of the coverage:

    Surface-Diffusion      = TRUE
    Surface-TotalDiffusion = TRUE

`Surface-Diffusion` averages the coverage over the sub-surface elements of each catalytic boundary separately,
`Surface-TotalDiffusion` over all catalytic boundaries together. The averaging is global, i.e. independent of
the number of processes and of the domain decomposition.

#### Surface sampling and output

With `Particles-DSMC-CalcSurfaceVal = TRUE` the surface output contains the coverage of every species,
`SpecXXX_Coverage`, as an instantaneous value per sub-surface element, and the catalytic heat flux,
`Catalytic_HeatFlux`, in W/m$^2$. The latter is a sampling quantity: the released and consumed energies are
accumulated only while the sampling is active and are normalised by the sampled area and the sampling duration,
exactly like the remaining wall quantities. The catalytic contribution is also included in `Total_HeatFlux`.

A positive catalytic heat flux denotes energy transferred into the wall. Adsorption and exothermic reactions
contribute positively, desorption negatively, so the sign of the net value indicates which process dominates.

#### Parameter read-in from the species database

All rate parameters of a catalytic reaction can be retrieved from the species database instead of the parameter
file. The reaction is looked up by the name given in `SurfName`, e.g. `Adsorption_CO_Pt`:

    Surface-Reaction1-SurfName = Adsorption_CO_Pt
    OverwriteCatParameters     = FALSE

Only setup-independent numbers are taken from the database: reaction enthalpy, scaling factor, accommodation
coefficient, prefactors, activation energies, the coverage and temperature windows and the desorption order.
Reactants, products, reaction type, boundaries and the reaction name always come from the parameter file, since
species indices are specific to a setup and a database is meant to be shared between them.

If the dataset of a reaction is not found in the database, that reaction falls back to the parameter file while
the others keep their database values. Setting `OverwriteCatParameters = TRUE` forces the parameter file for all
reactions.

## Deposition of Charges on resolved Dielectric Surfaces
This deposition of charges is designed for thick layers of dielectric materials, which are resolved by mesh
elements can directly be modelled as described in Section {ref}`sec:dielectric-materials`.
Charged particles can be absorbed (or reflected and leave their charge behind) at dielectric surfaces
when using the deposition method `cell_volweight_mean`. The boundary can be used by specifying

    Part-Boundary1-Condition         = reflective
    Part-Boundary1-Dielectric        = T
    Part-Boundary1-NbrOfSpeciesSwaps = 3
    Part-Boundary1-SpeciesSwaps1     = (/1,0/) ! e-
    Part-Boundary1-SpeciesSwaps2     = (/2,2/) ! Ar
    Part-Boundary1-SpeciesSwaps3     = (/3,2/) ! Ar+

which sets the boundary dielectric and the given species swap parameters effectively remove
electrons ($e^{-}$) on impact, reflect $Ar$ atoms and neutralize $Ar^{+}$ ions by swapping these to $Ar$ atoms.
Note that currently only singly charged particles can be handled this way. When multiple charged
particles would be swapped, their complete charge mus be deposited at the moment.

The boundary must also be specified as an *inner* boundary via

    BoundaryName                     = BC_INNER
    BoundaryType                     = (/100,0/)

or directly in the *hopr.ini* file that is used for creating the mesh.

(sec:distributed-capacitance-boundary-condition-for-particles)=
## Deposition of charges on distributed capacitance boundary condition (DCBC) surfaces
Charged particles impacting on distributed capacitance boundary condition (DCBC) surfaces are deposited on the surface element face
via linear weighting, which calculates the surface charge density $\sigma$ in the equation given in
{ref}`sec:distributed-capacitance-boundary-condition` and requires the following parameter settings

    Part-Boundary1-Condition        = reflective ! Surface charging requires the boundary condition "reflective"
    Part-Boundary1-UseSurfaceCharge = T          ! Activate surface charging
    Part-Boundary1-DC-BiasVoltage   = 1000.0     ! Electric potential "Phi_0" of the dielectric layer
    Part-Boundary1-DC-Permittivity  = 10.0       ! Relative permittivity "eps_r" of the dielectric layer
    Part-Boundary1-DC-Thickness     = 2.0e-3     ! Thickness "d" of the dielectric layer

and the last three parameters listed here are described in Section {ref}`sec:distributed-capacitance-boundary-condition`.
