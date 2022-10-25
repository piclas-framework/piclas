# Boundary Conditions - Field Solver

Boundary conditions are defined in the mesh creation step in the hopr.ini file and can be modified when running PICLas in the
corresponding *parameter.ini* file. In the *hopr.ini* file, which is read by the *hopr* executable, a boundary is defined by

    BoundaryName = BC_Inflow   ! BC index 1 (from  position in the parameter file)
    BoundaryType = (/4,0,0,0/) ! (/ Type, curveIndex, State, alpha /)

where the name of the boundary is directly followed by its type definition, which contains information on the BC type, curving,
state and periodicity. This can be modified in the *parameter.ini* file, which is read by the *piclas* executable via

    BoundaryName = BC_Inflow ! BC name in the mesh.h5 file
    BoundaryType = (/5,0/)   ! (/ Type, State /)

In this case the boundary type is changed from 4 (in the mesh file) to 5 in the simulation.

## Maxwell's Equations

The boundary conditions used for Maxwell's equations are defined by the first integer value in the *BoundaryType* vector (consisting of the *Type* and *State*) and
include, periodic, Dirichlet, Silver-Mueller, perfectly conducting, symmetry and reference state boundaries as detailed in the following table.

| BoundaryType |   Type    | State                                                                                                                     |
| :----------: | :-------: | :------------------------------------------------------------------------------------------------------------------------ |
|   (/1,1/)    | periodic  | 1: positive direction of the 1st periodicity vector                                                                       |
|   (/1,-1/)   | periodic  | -1: negative (opposite) direction of the 1st periodicity vector                                                           |
|              |           |                                                                                                                           |
|   (/2,2/)    | Dirichlet | 2: Coaxial waveguide                                                                                                      |
|   (/2,3/)    | Dirichlet | 3: Resonator                                                                                                              |
|   (/2,4/)    | Dirichlet | 4: Electromagnetic dipole (implemented via RHS source terms and shape function deposition)                                |
|   (/2,40/)   | Dirichlet | 40: Electromagnetic dipole without initial condition (implemented via RHS source terms and shape function deposition)     |
|   (/2,41/)   | Dirichlet | 41: Pulsed Electromagnetic dipole (implemented via RHS source terms and shape function deposition)                        |
|   (/2,5/)    | Dirichlet | 5: Transversal Electric (TE) plane wave in a circular waveguide                                                           |
|   (/2,7/)    | Dirichlet | 7: Special manufactured Solution                                                                                          |
|   (/2,10/)   | Dirichlet | 10: Issautier 3D test case with source (Stock et al., div. correction paper), domain [0;1]^3                              |
|   (/2,12/)   | Dirichlet | 12: Plane wave                                                                                                            |
|  (/2,121/)   | Dirichlet | 121: Pulsed plane wave (infinite spot size) and temporal Gaussian                                                         |
|   (/2,14/)   | Dirichlet | 14: Gaussian pulse is initialized inside the domain (usually used as initial condition and not BC)                        |
|   (/2,15/)   | Dirichlet | 15: Gaussian pulse with optional delay time *tDelayTime*                                                                  |
|   (/2,16/)   | Dirichlet | 16: Gaussian pulse which is initialized in the domain and used as a boundary condition for t>0                            |
|   (/2,50/)   | Dirichlet | 50: Initialization and BC Gyrotron - including derivatives                                                                |
|   (/2,51/)   | Dirichlet | 51: Initialization and BC Gyrotron - including derivatives (nothing is set for z>eps)                                     |
|              |           |                                                                                                                           |
|   (/3,0/)    |    SM     | 1st order absorbing BC (Silver-Mueller) - Munz et al. 2000 / Computer Physics Communication 130, 83-117 with fix          |
|              |           | of div. correction field for low B-fields that only set the correction fields when ABS(B)>1e-10                           |
|   (/5,0/)    |    SM     | 1st order absorbing BC (Silver-Mueller) - Munz et al. 2000 / Computer Physics Communication 130, 83-117                   |
|   (/6,0/)    |    SM     | 1st order absorbing BC (Silver-Mueller) - Munz et al. 2000 / Computer Physics Communication 130, 83-117 with fix          |
|              |           | of div. correction field for low B-fields that only set the correction fields when B is significantly large compared to E |
|              |           |                                                                                                                           |
|   (/4,0/)    |    PEC    | Perfectly conducting surface (Munz, Omnes, Schneider 2000, pp. 97-98)                                                     |
|              |           |                                                                                                                           |
|   (/10,0/)   | Symmetry  | Symmetry BC (perfect MAGNETIC conductor, PMC)                                                                             |
|              |           |                                                                                                                           |
|   (/20,0/)   |    Ref    | Use state that is read from .h5 file and interpolated to the BC                                                           |

Dielectric -> type 100?

## Poisson's Equation

The boundary conditions used for Maxwell's equations are defined by the first integer value in the *BoundaryType* vector (consisting of the *Type* and *State*) and
include, periodic, Dirichlet (via pre-defined function, zero-potential or *RefState*), Neumann and reference state boundaries
as detailed in the following table.

| BoundaryType |   Type    | State                                                                                                                      |
| :----------: | :-------: | :------------------------------------------------------------------------------------------------------------------------- |
|   (/1,1/)    | periodic  | 1: positive direction of the 1st periodicity vector                                                                        |
|   (/1,-1/)   | periodic  | -1: negative (opposite) direction of the 1st periodicity vector                                                            |
|              |           |                                                                                                                            |
|   (/2,0/)    | Dirichlet | 0: Phi=0                                                                                                                   |
|   (/2,2/)    | Dirichlet | 2: Automatic adjustment for Phi to meet const. input power, see {ref}`sec:fixed-coupled-power`                             |
|  (/2,1001/)  | Dirichlet | 1001: linear potential y-z via Phi = 2340y + 2340z                                                                         |
|  (/2,101/)   | Dirichlet | 101: linear in z-direction: z=-1: 0, z=1, 1000                                                                             |
|  (/2,103/)   | Dirichlet | 103: dipole                                                                                                                |
|  (/2,104/)   | Dirichlet | 104: solution to Laplace's equation: Phi_xx + Phi_yy + Phi_zz = 0                                                          |
|              |           | $\Phi=(COS(x)+SIN(x))(COS(y)+SIN(y))(COSH(SQRT(2.0)z)+SINH(SQRT(2.0)z))$                                                   |
|  (/2,200/)   | Dirichlet | 200: Dielectric Sphere of Radius R in constant electric field E_0 from book: John David Jackson, Classical Electrodynamics |
|  (/2,300/)   | Dirichlet | 300: Dielectric Slab in z-direction of half width R in constant electric field E_0:                                        |
|              |           | adjusted from CASE(200)                                                                                                    |
|  (/2,301/)   | Dirichlet | 301: like CASE=300, but only in positive z-direction the dielectric region is assumed                                      |
|  (/2,400/)   | Dirichlet | 400: Point Source in Dielectric Region with                                                                                |
|              |           | epsR_1  = 1  for x $<$ 0 (vacuum)                                                                                          |
|              |           | epsR_2 != 1 for x $>$ 0 (dielectric region)                                                                                |
|              |           |                                                                                                                            |
|   (/4,0/)    | Dirichlet | zero-potential (Phi=0)                                                                                                     |
|              |           |                                                                                                                            |
|   (/5,1/)    | Dirichlet | 1: use RefState Nbr 1, see details below                                                                                   |
|              |           |                                                                                                                            |
|   (/6,1/)    | Dirichlet | 1: use RefState Nbr 1, see details below                                                                                   |
|              |           |                                                                                                                            |
|   (/7,1/)    | Dirichlet | 1: use LinState Nbr 1, linear function for Phi, see {ref}`sec:linear-potential`                                            |
|              |           |                                                                                                                            |
|   (/10,0/)   |  Neumann  | zero-gradient (dPhi/dn=0)                                                                                                  |
|   (/11,0/)   |  Neumann  | q*n=1                                                                                                                      |

### RefState boundaries {-}

For each boundary of type *5* (reference state boundary *RefState*), e.g., by setting the boundary in the *parameter.ini* file

    BoundaryName = BC_WALL ! BC name in the mesh.h5 file
    BoundaryType = (/5,1/) ! (/ Type, State/)

the corresponding *RefState* number must also be supplied in the parameter.ini file (here 1) and is selected from its position
in the parameter file.
Each *RefState* is defined in the *parameter.ini* file by supplying a value for the voltage an alternating frequency for the cosine
function (a frequency of 0 results in a fixed potential over time) and phase shift

    RefState = (/-0.18011, 1.0, 0.0/) ! RefState Nbr 1: Voltage, Frequency and Phase shift

This yields the three parameters used in the cosine function

    Phi(t) = A*COS(2*pi*f*t + psi)

where *A=-0.18011* is the amplitude, *t* is the time, *f=1* is the frequency and *psi=0* is the phase shift.

Similar to boundary type *5* is type *6*, which simply uses a cosine function that always has the same sign, depending on the
amplitude *A*

    Phi(t) = (A/2) * (COS(2*pi*f*t + psi) + 1)

(sec:linear-potential)=
### Linear potential function

A linear function that ramps the electric potential from 0 V to a user-defined value can be applied to a boundary via

    BoundaryName = BC_WALL ! BC name in the mesh.h5 file
    BoundaryType = (/7,1/) ! 1: 1st LinState

Additionally, this specific boundary condition requires a starting position `LinPhiBasePoint` and
a direction along which the potential varies `LinPhiNormal`. The distance along which the potential varies as well as the final
value are defined by `LinPhiHeight` and `LinPhi`, respectively.
Coordinates below and above this distance are simply set to 0 V and the defined value, respectively.
The example below creates a linear ramp from 0 V to 1000 V starting at 1 mm in z-direction and ramps the value over 10 mm in the
same direction.

    LinPhiBasePoint = (/0. , 0. , 1e-3/)
    LinPhiNormal    = (/0. , 0. , 1.0/)
    LinPhiHeight    = 10e-3
    LinPhi          = 1000.

The linear potential uses the same functionality as RefState, hence, when two different functions are to be defined use the
following example

    BoundaryName    = BC_right
    BoundaryType    = (/7,1/)          ! 7: Dirichlet with linear ramp 1st LinState
    LinPhiBasePoint = (/0. , 0. , 0./) ! 1st LinState
    LinPhiNormal    = (/1. , 0. , 0./) ! 1st LinState
    LinPhiHeight    = 1.0              ! 1st LinState
    LinPhi          = 1e3              ! 1st LinState

    BoundaryName    = BC_left
    BoundaryType    = (/7,2/)          ! 7: Dirichlet with linear ramp 2nd LinState
    LinPhiBasePoint = (/0. , 0. , 0./) ! 2nd LinState
    LinPhiNormal    = (/1. , 0. , 0./) ! 2nd LinState
    LinPhiHeight    = 1.0              ! 2nd LinState
    LinPhi          = 0.0              ! 2nd LinState

(sec:fixed-coupled-power)=
### Fixed coupled power (const. input power)

An automatic adjustment of the electric potential to ensure that a fixed power input to the system is achieved requires the
following parameters

    BoundaryName = BC_WALL ! BC name in the mesh.h5 file
    BoundaryType = (/2,2/) ! all BCs with this type will be adjusted to the same electric potential that is adjusted over time

Additionally, a starting value for the potential, lower and upper boundaries and a relaxation factor are required
as well as the target input power, which is set via

    CoupledPowerPotential = (/10. , 1000. , 2000./) ! lower, starting and maximum values for the electric potential
    CoupledPowerRelaxFac  = 0.05  ! the new potential is updated by 5% in each time step
    CoupledPowerTarget    = 1e-10 ! target power of 1e-10 Watt

The values in `CoupledPowerPotential` correspond to the lower boundary, the starting value and the upper boundary, respectively.
When a simulation is restarted from a state file, the last known value of the BC will be used instead of the starting value, which
is only applied when starting a fresh simulation from scratch.

### Zero potential enforcement

It is important to note that when no Dirichlet boundary conditions are selected by the user, the code automatically enforces mixed
boundaries on either Neumann or periodic boundaries. Depending on the simulation domain, the direction with the largest extent is
selected and on those boundaries an additional Dirichlet boundary condition with $\phi=0$ is enforced to ensure convergence of the
HDG solver. The boundary conditions selected by the user are only altered at these locations and not removed.
The information regarding the direction that is selected for this purpose is printed to std.out with the following line

     |      Zero potential side activated in direction (1: x, 2: y, 3: z) |        1 |  OUTPUT |

To selected the direction by hand, simply supply the desired direction via

    HDGZeroPotentialDir = 1

with 1: x-, 2: y-, 3: z-direction.

## Dielectric Materials

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

where the following pre-defined cases are available as given in table {numref}`tab:dielectric_test_cases`

```{table} Dielectric Test Cases
---
name: tab:dielectric_test_cases
---
  |            Option            |   Additional Parameters   |                                         Notes                                          |
  | :--------------------------: | :-----------------------: | :------------------------------------------------------------------------------------: |
  |        `FishEyeLens`         |           none            | function with radial dependence: $\varepsilon_{r}=n_{0}^{2}/(1 + (r/r_{max})^{2})^{2}$ |
  |           `Circle`           | `DielectricRadiusValue`,  |             Circular dielectric in x-y-direction (constant in z-direction)             |
  |                              | `DielectricRadiusValueB`, |             with optional cut-out radius DielectricRadiusValueB along the              |
  |                              |  `DielectricCircleAxis`   |                           axis given by DielectricCircleAxis                           |
  | `DielectricResonatorAntenna` |  `DielectricRadiusValue`  |            Circular dielectric in x-y-direction (only elements with $z>0$)             |
  |          `FH_lens`           |           none            |          specific geometry (`SUBROUTINE SetGeometry` yields more information)          |
```

For the Maxwell solver (DGSEM), the interface fluxes between vacuum and dielectric regions can
either be conserving or non-conserving, which is selected by

    DielectricFluxNonConserving = T

which uses non-conserving fluxes. This is recommended for improved simulation results, as described in {cite}`Copplestone2019b`.
When particles are to be considered in a simulation, these are generally removed from dielectric
materials during the emission (inserting) stage, but may be allowed to exist within dielectrics by
setting

    DielectricNoParticles = F

which is set true by default, hence, removing the particles.


