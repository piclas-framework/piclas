(sec:particle-boundary-conditions)=
# Boundary Conditions - Particle Solver

Within the parameter file it is possible to define different particle boundary conditions. The number of boundaries is defined by

    Part-nBounds=2
    Part-Boundary1-SourceName=BC_OPEN
    Part-Boundary1-Condition=open
    Part-Boundary2-SourceName=BC_WALL
    Part-Boundary2-Condition=reflective
    Part-Boundary2-SurfaceModel=2

The `Part-Boundary1-SourceName=` corresponds to the name given during the preprocessing step with HOPR. The available conditions
(`Part-Boundary1-Condition=`) are described in the table below.

|   Condition    | Description                                                                                                          |
| :------------: | :------------------------------------------------------------------------------------------------------------------- |
|     `open`     | Every particle crossing the boundary will be deleted.                                                                |
|  `symmetric`   | A perfect specular reflection, without sampling of particle impacts.                                                 |
|  `reflective`  | Allows the definition of specular and diffuse reflection, Section {ref}`sec:particle-boundary-conditions-reflective` |
| `rot_periodic` | Allows the definition of rotational periodicity, Section {ref}`sec:particle-boundary-conditions-rotBC`               |

(sec:particle-boundary-conditions-reflective)=
## Reflective Wall

A reflective boundary can be defined with

    Part-Boundary2-SourceName=BC_WALL
    Part-Boundary2-Condition=reflective

A perfect specular reflection is performed, if no other parameters are given. Gas-surface interactions can be modelled with the
extended Maxwellian model {cite}`Padilla2009`, using accommodation coefficients of the form

$$\alpha = \frac{E_i-E_r}{E_i - E_w}$$

where $i$, $r$ and $w$ denote the incident, reflected and wall energy, respectively.  The coefficient `MomentumACC` is utilized to
decide whether a diffuse (`MomentumACC` $>R$) or specular reflection (`MomentumACC` $<R$) occurs upon particle impact, where
$R=[0,1)$ is a random number. Separate accommodation coefficients can be defined for the translation (`TransACC`), rotational
(`RotACC`), vibrational (`VibACC`) and electronic energy (`ElecACC`) accommodation at a constant wall temperature [K].

    Part-Boundary2-MomentumACC=1.
    Part-Boundary2-WallTemp=300.
    Part-Boundary2-TransACC=1.
    Part-Boundary2-VibACC=1.
    Part-Boundary2-RotACC=1.
    Part-Boundary2-ElecACC=1.

An additional option `Part-Boundary2-SurfaceModel` is available, that is used for heterogeneous reactions (reactions that have reactants
in two or more phases) or secondary electron emission models. These models are described in detail in Section {ref}`sec:surface-chemistry`.

(sec:particle-boundary-conditions-reflective-wallvelo)=
### Wall movement (Linear & rotational)

Additionally, a linear wall velocity [m/s] can be given

    Part-Boundary2-WallVelo=(/0,0,100/)

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
temperature gradient will only be applied along the chosen the direction. As oppposed to the default case, the positions of the
surfaces are not projected onto the gradient vector before checking wether they are inside the box spanned by `TempGradStart` and
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
It should be noted that the the adaptation should be performed multiple times to achieve a converged temperature distribution.

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
periodic boundary is tranformed according to this angle, which is defined by the right-hand rule and the rotation axis:

    Part-RotPeriodicAxi = 1    ! (x = 1, y = 2, z = 3)

The usage of rotational periodic boundary conditions is limited to cases, where the rotational periodic axis is one of the three
Cartesian coordinate axis (x, y, z) with its origin at (0, 0, 0).

### Intermediate Plane Definition
If several segments with different rotation angles are defined, exactly two corresponding BCs must be defined for each segment.
Since the plane between these segments with different rotational symmetry angles represents a non-conforming connection, additional 
two BCs must be defined as `rot_periodic_inter_plane` at this intermediate plane. Both BCs must refer to each other in the 
definition in order to ensure the connection.

    Part-Boundary40-SourceName       = BC_INT_R1_BOT
    Part-Boundary40-Condition        = rot_periodic_inter_plane
    Part-Boundary40-AssociatedPlane  = 41

    Part-Boundary41-SourceName       = BC_INT_S1_TOP
    Part-Boundary41-Condition        = rot_periodic_inter_plane
    Part-Boundary41-AssociatedPlane  = 40 

Note that using the intermediate plane definition with two corresponding BCs allows the user to mesh the segments independently, 
creating a non-conforming interface at the intermediate plane. However, use of these non-conformal grids has so far only been 
possible in pure DSMC simulations.

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
appropriate particle boundary surface model `Part-BoundaryX-SurfaceModel`.
The available conditions (`Part-BoundaryX-SurfaceModel=`) are described in the table below.

|    Model    | Description                                                                                                                                                                                  |
| :---------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 0 (default) | Standard extended Maxwellian scattering                                                                                                                                                      |
|      1      | Empirical modelling of sticking coefficient/probability                                                                                                                                                  |
|      5      | Secondary electron emission as given by Ref. {cite}`Levko2015`.                                                                                                                              |
|      7      | Secondary electron emission due to ion impact (SEE-I with $Ar^{+}$ on different metals) as used in Ref. {cite}`Pflug2014` and given by Ref. {cite}`Depla2009` with a default yield of 13 \%. |
|      8      | Secondary electron emission due to ion impact (SEE-E with $e^{-}$ on dielectric surfaces) as used in Ref. {cite}`Liu2010` and given by Ref. {cite}`Morozov2004`.                             |
|      9      | Secondary electron emission due to ion impact (SEE-I with $Ar^{+}$) with a constant yield of 1 \%. Emitted electrons have an energy of 6.8 eV upon emission.                                 |
|     10      | Secondary electron emission due to ion impact (SEE-I with $Ar^{+}$ on copper) as used in Ref. {cite}`Theis2021` originating from {cite}`Phelps1999`                                          |
|     11      | Secondary electron emission due to electron impact (SEE-E with $e^{-}$ on quartz (SiO$_{2}$)) as described in Ref. {cite}`Zeng2020` originating from {cite}`Dunaevsky2003`                   |

For surface sampling output, where the surface is split into, e.g., $3\times3$ sub-surfaces, the following parameters mus be set

    BezierSampleN                 = 3
    DSMC-nSurfSample              = 3
    Part-WriteMacroSurfaceValues  = T
    Particles-DSMC-CalcSurfaceVal = T
    Part-IterationForMacroVal     = 200

where `BezierSampleN=DSMC-nSurfSample`. In this example, sampling is performed over and every 200 iterations.

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

    Particles-Species-Database = Species_Database.h5

As additional output, the cell-local sticking coefficient will be added to the sampled surface output. A particle sticking to the surface will be deleted and its energy added to the heat flux sampling. This model can be combined with the linear temperature gradient and radiative equilibrium modelling as described in Section {ref}`sec:particle-boundary-conditions-reflective`.

### Secondary Electron Emission (SEE)

Different models are implemented for secondary electron emission that are based on either electron or ion bombardment, depending on
the surface material. All models require the specification of the electron species that is emitted from the surface via

    Part-SpeciesA-PartBoundB-ResultSpec = C

where electrons of species `C` are emitted from boundary `B` on the impact of species `A`.

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
The emission yield represents the ratio of emitted electrons vs. impacting ions and the emission energy is given in electronvolts.
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

## Deposition of Charges on Dielectric Surfaces

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

