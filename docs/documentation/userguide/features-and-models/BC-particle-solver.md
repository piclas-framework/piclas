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

|   Condition    | Description                                                                                                                                                                                 |
| :------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
|     `open`     | Every particle crossing the boundary will be deleted.                                                                                                                                       |
|  `reflective`  | Allows the definition of specular and diffuse reflection. A perfect specular reflection is performed, if no other parameters are given (discussed in more detail in the following section). |
|  `symmetric`   | A perfect specular reflection, without sampling of particle impacts.                                                                                                                        |
| `rot_periodic` | Allows the definition of rotational periodicity.                                                                                                                                            |

For `reflective` boundaries, an additional option `Part-Boundary2-SurfaceModel` is available, that
is used for heterogeneous reactions (reactions that have reactants in two or more phases) or secondary electron emission models.
These models are described in {ref}`sec:surface-chemistry`.

For `rot_periodic` exactly two corresponding boundaries must be defined. Every particle crossing one of these boundaries will be
inserted at the corresponding other boundary that is rotationally shifted.

## Diffuse Wall

Gas-surface interaction can be modelled with the extended Maxwellian model {cite}`Padilla2009`, using accommodation coefficients
of the form

$$\alpha = \frac{E_i-E_r}{E_i - E_w}$$

where $i$, $r$ and $w$ denote the incident, reflected and wall energy, respectively.  The coefficient `MomentumACC` is utilized to
decide whether a diffuse (`MomentumACC` $>R$) or specular reflection (`MomentumACC` $<R$) occurs upon particle impact, where
$R=[0,1)$ is a random number. Separate accommodation coefficients can be defined for the translation (`TransACC`), rotational
(`RotACC`), vibrational (`VibACC`) and electronic energy (`ElecACC`) accommodation at a constant wall temperature [K].

    Part-Boundary2-SourceName=BC_WALL
    Part-Boundary2-Condition=reflective
    Part-Boundary2-MomentumACC=1.
    Part-Boundary2-WallTemp=300.
    Part-Boundary2-TransACC=1.
    Part-Boundary2-VibACC=1.
    Part-Boundary2-RotACC=1.
    Part-Boundary2-ElecACC=1.

### Wall movement (Linear & rotational)

Additionally, a linear wall velocity [m/s] can be given

    Part-Boundary2-WallVelo=(/0,0,100/)

In the case of rotating walls the `-RotVelo` flag, a rotation frequency [Hz], a origin of rotation axis (x, y, z coordinates) and
the rotation axis vector must be set. Note that the definition of rotation direction is given by the rotation axis and the
right-hand rule.

    Part-Boundary2-RotVelo = T
    Part-Boundary2-RotFreq = 100
    Part-Boundary2-RotOrg = (/0.,0.,0./)
    Part-Boundary2-RotAxi = (/0.,0.,1./)

### Linear temperature gradient

A linear temperature gradient across a boundary can be defined by supplying a second wall temperature and the start and end vector

    Part-Boundary2-WallTemp2=500.
    Part-Boundary2-TemperatureGradientStart=(/0.,0.,0./)
    Part-Boundary2-TemperatureGradientEnd=(/0.,0.,1./)

Between these two points the temperature will be interpolated, where the start vector corresponds to the first wall temperature,
whereas the end vector to the second wall temperature. Beyond these position values, the first and second temperature will be used
as the constant wall temperature, respectively.

### Radiative equilibrium

Another option is to adapt the wall temperature based on the heat flux assuming that the wall is in radiative equilibrium.
The temperature is then calculated from

$$ q_w = \varepsilon \sigma T_w^4,$$

where $\varepsilon$ is the radiative emissivity of the wall (default = 1) and
$\sigma = \SI{5.67E-8}{\watt\per\square\meter\per\kelvin\tothe{4}}$ is the Stefan-Boltzmann constant. The adaptive boundary is
enabled by

    Part-AdaptWallTemp = T
    Part-Boundary1-UseAdaptedWallTemp = T
    Part-Boundary1-RadiativeEmissivity = 0.8

If provided, the wall temperature will be adapted during the next output of macroscopic variables, where the heat flux calculated
during the preceding sampling period is utilized to determine the side-local temperature. The temperature is included in the
`State` file and thus available during a restart of the simulation. The surface output (in `DSMCSurfState`) will additionally
include the temperature distribution in the `Wall_Temperature` variable (see Section {ref}`sec:sampled-flow-field-and-surface-variables`).
To continue the simulation without further adapting the temperature, the first flag has to be disabled (`Part-AdaptWallTemp = F`).
It should be noted that the the adaptation should be performed multiple times to achieve a converged temperature distribution.

## Rotational Periodicity

The rotational periodic boundary condition can be used in order to reduce the computational effort in case of an existing
rotational periodicity. In contrast to symmetric boundary conditions, a macroscopic flow velocity in azimuthal direction can be
simulated (e.g. circular flow around a rotating cylinder). Exactly two corresponding boundaries must be defined by setting
`rot_periodic` as BC condition and the rotation direction for each BCs.

    Part-Boundary1-SourceName=BC_Rot_Peri_plus
    Part-Boundary1-Condition=rot_periodic
    Part-Boundary1-RotPeriodicDir=1

    Part-Boundary2-SourceName=BC_Rot_Peri_minus
    Part-Boundary2-Condition=rot_periodic
    Part-Boundary2-RotPeriodicDir=-1

CAUTION! The correct sign for the direction must be determined. Here, the rotation direction is defined by the rotation axis
`Part-RotPeriodicAxi` that must be defined separately, and the right-hand rule.

    Part-RotPeriodicAxi=1    ! (x=1, y=2, z=3)

The usage of rotational periodic boundary conditions is limited to cases, where the rotational periodic axis is one of the three
Cartesian coordinate axis (x, y, z) with its origin at (0, 0, 0). Finally the rotation angle `Part-RotPeriodicAngle` [°] must be
defined.

    Part-RotPeriodicAngle=90

## Porous Wall / Pump

The porous boundary condition uses a removal probability to determine whether a particle is deleted or reflected at the boundary.
The main application of the implemented condition is to model a pump, according to {cite}`Lei2017`. It is defined by giving the
number of porous boundaries and the respective boundary number (`BC=2` corresponds to the `BC_WALL` boundary defined in the
previous section) on which the porous condition is.

    Surf-nPorousBC=1
    Surf-PorousBC1-BC=2
    Surf-PorousBC1-Pressure=5.
    Surf-PorousBC1-Temperature=300.
    Surf-PorousBC1-Type=pump
    Surf-PorousBC1-PumpingSpeed=2e-9
    Surf-PorousBC1-DeltaPumpingSpeed-Kp=0.1
    Surf-PorousBC1-DeltaPumpingSpeed-Ki=0.0

The removal probability is determined through the given pressure [Pa] and temperature [K] at the boundary. A pumping speed can be
given as a first guess, however, the pumping speed $S$ [$m^3/s$] will be adapted if the proportional factor ($K_{\mathrm{p}}$,
`DeltaPumpingSpeed-Kp`) is greater than zero

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
    Surf-PorousBC1-Pressure=5.
    Surf-PorousBC1-Temperature=300.
    Surf-PorousBC1-Type=sensor

Together with a region definition, a pump as well as a sensor can be defined on a single and/or multiple boundaries, allowing e.g.
to determine the pressure difference between the pump and a remote area of interest.

(sec:surface-chemistry)=
## Surface Chemistry

Modelling of reactive surfaces is enabled by setting `Part-BoundaryX-Condition=reflective` and an
appropriate particle boundary surface model `Part-BoundaryX-SurfaceModel`.
The available conditions (`Part-BoundaryX-SurfaceModel=`) are described in the table below.

|    Model    | Description                                                                                                                                                                         |
| :---------: | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 0 (default) | Standard extended Maxwellian scattering                                                                                                                                             |
|      2      | Simple recombination on surface collision, where an impinging particle as given by Ref. {cite}`Reschke2019`.                                                                             |
|      3      | Kinetic Monte Carlo surface: Replicates surfaces with a specified lattice structure, either fcc(100) or fcc(111) and models complete catalysis as given by Ref. {cite}`Reschke2019`.     |
|      5      | Secondary electron emission as given by Ref. {cite}`Levko2015`.                                                                                                                          |
|      7      | Secondary electron emission due to ion impact (SEE-I with $Ar^{+}$ on different metals) as used in Ref. {cite}`Pflug2014` and given by Ref. {cite}`Depla2009` with a constant yield of 13 \%. |
|     101     | Evaporation from surfaces according to a Maxwellian velocity distribution.                                                                                                          |
|     102     | Evaporation according to MD-fitted velocity distributions.                                                                                                                          |

For surface sampling output, where the surface is split into, e.g., $3\times3$ sub-surfaces, the following parameters mus be set

    BezierSampleN = 3
    DSMC-nSurfSample = 3
    Part-WriteMacroSurfaceValues = T
    Particles-DSMC-CalcSurfaceVal = T
    Part-IterationForMacroVal = 200

where `BezierSampleN=DSMC-nSurfSample`. In this example, sampling is performed over 200 iterations.

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

