# Features of the Particle Solver

This section describes common features, which are available to the particle-based methods such as PIC, DSMC, MCC, BGK and FP.

(sec:macroscopic-restart)=
## Macroscopic Restart

The so-called macroscopic restart, allows to restart the simulation by using an output file of a previous simulation run (the
regular state file has still to be supplied). This enables to change the weighting factor, without beginning a new simulation.

    Particles-MacroscopicRestart = T
    Particles-MacroscopicRestart-Filename = Test_DSMCState.h5

The particle velocity distribution within the domain is then generated assuming a Maxwell-Boltzmann distribution, using the
translational temperature per direction of each species per cell. The rotational and vibrational energy per species is initialized
assuming an equilibrium distribution.

(sec:variable-time-step)=
## Variable Time Step

A spatially variable or species-specific time step (VTS) can be activated for steady-state simulations, where three options are
currently available and described in the following:

* Distribution: use a simulation result to adapt the time step in order to resolve physical parameters (e.g. collision frequency)
* Linear scaling: use a linearly increasing/decreasing time step along a given direction
* Species-specific: every species can have its own time step

### Distribution

The first option is to adapt the time step during a simulation restart based on certain parameters of the simulation such as
maximal collision probability (DSMC), mean collision separation distance over mean free path (DSMC), maximal relaxation factor
(BGK/FP) and particle number. This requires the read-in of a DSMC state file that includes DSMC quality factors (see Section
{ref}`sec:DSMC-quality`).

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

The second flag allows to enable/disable the adaptation of the time step distribution. Typically, a simulation would be performed
until a steady-state (or close to it, e.g. the particle number is not increasing significantly anymore) is reached with a uniform
time step. Then a restart with the above options would be performed, where the time step distribution is adapted using the DSMC
output of the last simulation. Now, the user can decide to continue adapting the time step with the subsequent DSMC outputs (Note:
Do not forget to update the DSMCState file name!) or to disable the adaptation and to continue the simulation with the distribution
from the last simulation (the adapted particle time step is saved within the regular state file). The DSMC state file after a
successful simulation run will contain the adapted time step (as a factor of the `ManualTimeStep`) as an additional output.
The time step factor within the DSMC state is always given priority over the time step stored in the state file during the
adaptation step.

The `MaxFactor` and `MinFactor` allow to limit the adapted time step within a range of $f_{\mathrm{min}} \Delta t$ and
$f_{\mathrm{max}} \Delta t$. The time step adaptation can be used to increase the number of particles by defining a minimum
particle number (e.g `MinPartNum` = 10, optional). For DSMC, the parameters `TargetMCSoverMFP` (ratio of the mean collision
separation distance over mean free path) and `TargetMaxCollProb` (maximum collision probability) allow to modify the target values
for the adaptation. For the BGK and FP methods, the time step can be adapted according to a target maximal relaxation frequency.

The last two flags enable to initialize the particles distribution from the given DSMC state file, using the macroscopic properties
such as flow velocity, number density and temperature (see Section {ref}`sec:macroscopic-restart`). Strictly speaking, the VTS procedure
only requires the `Filename` for the read-in of the aforementioned parameters, however, it is recommended to perform a macroscopic
restart to initialize the correct particle number per cells. Otherwise, cells with a decreased/increased time step will require
some time until the additional particles have reached/left the cell.

The time step adaptation can also be utilized in coupled BGK-DSMC simulations, where the time step will be adapted in both regions
according to the respective criteria as the BGK factors are zero in the DSMC region and vice versa. Attention should be payed in
the transitional region between BGK and DSMC, where the factors are potentially calculated for both methods. Here, the time step
required to fulfil the maximal collision probability criteria will be utilized as it is the more stringent one.

### Linear scaling

The second option is to use a linearly increasing time step along a given direction. This option does not require a restart or a
previous simulation result. Currently, only a scaling along the three axis is implemented, where the sign additionally defines
the direction of the scaling. With the start point and end point, the region in which the linear increase should be performed can
be defined. To define the domain border as the end point in maximal x-direction, the vector `(/-99999.,0.0,0.0/)` should be supplied.
Finally, the `ScaleFactor` defines the maximum time step increase towards the end point $\Delta t (x_{\mathrm{end}})=f \Delta t$.

    Part-VariableTimeStep-LinearScaling = T
    Part-VariableTimeStep-ScaleFactor   = 2
    Part-VariableTimeStep-Direction     =      (/1.0,0.0,0.0/)
    Part-VariableTimeStep-StartPoint    =     (/-0.4,0.0,0.0/)
    Part-VariableTimeStep-EndPoint      =  (/-99999.,0.0,0.0/)

Besides DSMC, the linear scaling is available for the BGK and FP method. Finally, specific options for 2D/axisymmetric simulations
are discussed in Section {ref}`sec:2D-axisymmetric`.

### Species-specific time step

This option is decoupled from the other two time step options as the time step is not applied on a per-particle basis but for each species. Currently, its main application is for PIC-MCC simulations (only Poisson field solver with Euler, Leapfrog and Boris-Leapfrog time discretization methods), where there are large differences in the time scales (e.g. electron movement requires a time step of several orders of magnitude smaller than the ions). The species-specific time step is actvitated per species by setting a factor

    Part-Species1-TimeStepFactor = 0.01

that is multiplied with the provided time step. If no time step factor is provided, the default time step will be utilized. In this example, the species will be effectively simulated with a time step 100 smaller than the given time step. Since the number of iterations remains the same, the species will effectively progress slower in time and the defined `tend` is not applicable anymore. This model should not be confused with a sub-cycling method and only be used for steady-state cases.

To accelerate the convergence to steady-state, the following flag can be used to perform collisions and reactions at the regular time step:

    Part-VariableTimeStep-DisableForMCC = T

For species with a time step factor lower than 1, it is compared with a random number to decide whether the collision/reaction is performed for that species.

## Symmetric Simulations

For one-dimensional (e.g. shock-tubes), two-dimensional (e.g. cylinder) and axisymmetric (e.g. re-entry capsules) cases, the
computational effort can be greatly reduced.

(sec:1D-sym)=
### 1D Simulations

To enable one-dimensional simulations, the symmetry order has to be set

    Particles-Symmetry-Order=1

The calculation is performed along the $x$-axis. The $y$ and $z$ dimension should be centered to the $xz$-plane (i.e.
$|y_{\mathrm{min}}|=|y_{\mathrm{max}}|$). All sides of the hexahedrons must be parallel to the $xy$-, $xz$-, and $yz$-plane.
Boundaries in $y$ and $z$ direction shall be defined as `symmetric_dim`.

    Part-Boundary5-SourceName=SYM
    Part-Boundary5-Condition=symmetric_dim

The code assumes a length of $\Delta y = \Delta z = 1$, regardless of the actual dimension in $y$ and $z$. Therefore, the weighting factor should be adapted accordingly.

(sec:2D-axisymmetric)=
### 2D/Axisymmetric Simulations

To enable two-dimensional simulations, the symmetry order has to be set

    Particles-Symmetry-Order=2

Two-dimensional and axisymmetric simulations require a mesh in the $xy$-plane, where the $x$-axis is the rotational axis and $y$
ranges from zero to a positive value. Additionally, the mesh shall be centered around zero in the $z$-direction with a single cell
row, such as that $|z_{\mathrm{min}}|=|z_{\mathrm{max}}|$. It should be noted that when converting or creating the mesh using `HOPR`,
it is recommended to define the space filling curve by

    sfc_type = mortonZ    ! alternative: hilbertZ

in the `hopr.ini` parameter file. This utilizes a 2D space filling curve, which can reduce the simulation duration significantly as
it improves the distribution of the cells per processor. This applies even to certain 3D cases, such as a channel flow.

The rotational symmetry axis shall be defined as a separate boundary with the `symmetric_axis` boundary condition

    Part-Boundary4-SourceName=SYMAXIS
    Part-Boundary4-Condition=symmetric_axis

The boundaries (or a single boundary definition for both boundary sides) in the $z$-direction should be defined as symmetry sides
with the `symmetric_dim` condition

    Part-Boundary5-SourceName=SYM
    Part-Boundary5-Condition=symmetric_dim

It should be noted that the two-dimensional mesh assumes a length of $\Delta z = 1$, regardless of the actual dimension in $z$.
Therefore, the weighting factor should be adapted accordingly.

To enable axisymmetric simulations, the following flag is required

    Particles-Symmetry2DAxisymmetric=T

To fully exploit rotational symmetry, a radial weighting can be enabled, which will linearly increase the weighting factor $w$
towards $y_{\mathrm{max}}$ (i.e. the domain border in $y$-direction), depending on the current $y$-position of the particle.

    Part-Weight-Type = radial
    Part-Weight-Radial-ScaleFactor = 100

A radial weighting factor of 100 means that the weighting factor at $y_{\mathrm{max}}$ will be $100w$. Although greatly reducing
the number of particles, this introduces the need to delete and create (in the following "clone") particles, which travel upwards
and downwards in the $y$-direction, respectively. If the new weighting factor is smaller than the previous one, a cloning
probability is calculated by

$$ P_{\mathrm{clone}} = \frac{w_{\mathrm{old}}}{w_{\mathrm{new}}} - \mathrm{INT}\left(\frac{w_{\mathrm{old}}}{w_{\mathrm{new}}}\right)\qquad \text{for}\quad w_{\mathrm{new}}<w_{\mathrm{old}}.$$

For the deletion process, a deletion probability is calculated, if the new weighting factor is greater than the previous

$$ P_{\mathrm{delete}} = 1 - P_{\mathrm{clone}}\qquad \text{for}\quad w_{\mathrm{old}}<w_{\mathrm{new}}.$$

If the ratio between the old and the new weighting factor is $w_{\mathrm{old}}/w_{\mathrm{new}}> 2$, the time step or the radial
weighting factor should be reduced as the creation of more than one clone per particle per time step is not allowed. The same
applies if the deletion probability is above $0.5$.

For the cloning procedure, two methods are implemented, where the information of the particles to be cloned are stored for a
given number of iterations (`CloneDelay=10`) and inserted at the old position. The difference is whether the list is inserted
chronologically (`CloneMode=1`) or randomly (`CloneMode=2`) after the first number of delay iterations.

    Part-Weight-CloneMode=2
    Part-Weight-CloneDelay=10

This serves the purpose to avoid the so-called particle avalanche phenomenon {cite}`Galitzine2015`, where clones travel on the
exactly same path as the original in the direction of a decreasing weight. They have a zero relative velocity (due to the same
velocity vector) and thus a collision probability of zero. Combined with the nearest neighbor pairing, this would lead to an
ever-increasing number of identical particles travelling on the same path. An indicator how often identical particle pairs are
encountered per time step during collisions is given as an output (`2D_IdenticalParticles`, to enable the output see Section
{ref}`sec:DSMC-quality`). Additionally, it should be noted that a large delay of the clone insertion might be problematic for
time-accurate simulations. However, for the most cases, values for the clone delay between 2 and 10 should be sufficient to
avoid the avalance phenomenon.

Another issue is the particle emission on large sides in $y$-dimension close to the rotational axis. As particles are inserted
linearly along the $y$-direction of the side, a higher number density is inserted closer to the axis. This effect is directly
visible in the free-stream in the cells downstream, when using mortar elements, or in the heat flux (unrealistic peak) close to
the rotational axis. It can be avoided by splitting the surface flux emission side into multiple subsides with the following flag
(default value is 20)

    Part-Weight-SurfFluxSubSides = 20

An alternative to the particle position-based weighting is the cell-local radial weighting, which can be enabled by

    Part-Weight-UseCellAverage = T

However, this method is not preferable if the cell dimensions in $y$-direction are large, resulting in numerical artifacts due to
the clustered cloning processes at cell boundaries.
#### Variable Time Step: Linear scaling

The linear scaling of the variable time step is implemented slightly different to the 3D case. Here, a particle-based time step is
used, where the time step of the particle is determined on its current position. The first scaling is applied in the radial
direction, where the time step is increased towards the radial domain border. Thus, $\Delta t (y_{\mathrm{max}}) = f \Delta t$
and $\Delta t (y_{\mathrm{min}} = 0) = \Delta t$.

    Part-VariableTimeStep-LinearScaling = T
    Part-VariableTimeStep-ScaleFactor = 2

Additionally, the time step can be varied along the x-direction by defining a "stagnation" point, towards which the time step is
decreased from the minimum x-coordinate ($\Delta t (x_{\mathrm{min}}) = f_{\mathrm{front}}\Delta t$) and away from which the time
step is increased again towards the maximum x-coordinate ($\Delta t (x_{\mathrm{max}}) = f_{\mathrm{back}}\Delta t$). Therefore,
only at the stagnation point, the time step defined during the initialization is used.

    Part-VariableTimeStep-Use2DFunction = T
    Part-VariableTimeStep-StagnationPoint = 0.0
    Part-VariableTimeStep-ScaleFactor2DFront = 2.0
    Part-VariableTimeStep-ScaleFactor2DBack = 2.0

(sec:variable-particle-weighting)=
## Variable Particle Weighting

The particle weighting scheme for DSMC, BGK, and FP simulations can be defined by

    Part-Weight-Type = constant

The available options are described in the table below.

| Type       | Description                                                                                                                             |
| ---------- | :-------------------------------------------------------------------------------------------------------------------------------------- |
| constant   | Default case, constant weighting as set by `Part-Species$-MacroParticleFactor`                                                          |
| radial     | Radial weighting in y-direcition for the axisymmetric simulation based on the particle position, see Section {ref}`sec:2D-axisymmetric` |
| linear     | Linear weighting along an axis or user-defined vector, see Section {ref}`sec:linear-particle-weighting`                                 |
| cell_local | Cell-local weighting distribution, based on previous simulation results, see Section {ref}`sec:celllocal-particle-weighting`            |

An additional type of variable particle weighting scheme *split & merge*, currently only available for PIC simulations, is described in Section {ref}`sec:split-merge`.

(sec:linear-particle-weighting)=
### Linear Weighting in 3D

Linearly increasing particle weights in 3D can be defined similarly to the radial weighting in the 2D axisymmetric case.

    Part-Weight-Type = linear

As a first step, the scaling direction needs to be defined, either as one of the coordinate axes

    Part-Weight-Linear-CoordinateAxis = 2 ! y-axis

or by start and end coordinates from which a scaling vector is generated

    Part-Weight-Linear-StartPointForScaling = (/0,0,0/)
    Part-Weight-Linear-EndPointForScaling = (/1,1,1/)

In each case, scaling points with a given particle weight can be defined along the axis or the vector. The weight between two scaling points is then linearly interpolated based on the relative particle position (a value between 0 and 1). If the weight is to be increased along one of the coordinate axes, the absolute position [m] has to be given for the scaling points.

In the example below, the weights are increased from $10^8$ to $10^{10}$, over a length of 1 cm along the axes.

    Part-Weight-Linear-nScalePoints = 2

    Part-Weight-Linear-ScalePoint1-Coordinate = 0.0
    Part-Weight-Linear-ScalePoint1-Factor = 1E8

    Part-Weight-Linear-ScalePoint2-Coordinate = 0.01
    Part-Weight-Linear-ScalePoint2-Factor = 1E10

Analogously to the radial weighting, the particles will be cloned or deleted when they move to a new cell with a different weighting factor. It should be noted that while a change of the scaling parameters between simulation runs is possible, it is recommended to perform a macroscopic restart (see Section {ref}`sec:macroscopic-restart`) to avoid large jumps in the weighting factor distribution during the first iteration.

The cloning and deletion probability is again calculated by:

$$ P_{\mathrm{clone}} = \frac{w_{\mathrm{old}}}{w_{\mathrm{new}}} - \mathrm{INT}\left(\frac{w_{\mathrm{old}}}{w_{\mathrm{new}}}\right)\qquad \text{for}\quad w_{\mathrm{new}}<w_{\mathrm{old}}.$$

$$ P_{\mathrm{delete}} = 1 - P_{\mathrm{clone}}\qquad \text{for}\quad w_{\mathrm{old}}<w_{\mathrm{new}}.$$

The same parameters as for radial weighting can be utilized for the linear weighting:

    Part-Weight-CloneMode=2
    Part-Weight-CloneDelay=10
    Part-Weight-SurfFluxSubSides = 20
    Part-Weight-UseCellAverage = T

A detailed description of these parameters can be found in Section {ref}`sec:2D-axisymmetric`.

(sec:celllocal-particle-weighting)=
### Cell-local Weighting

An automatic determination of the optimal particle weights in each cell can be performed during a macroscopic restart (see Section {ref}`sec:macroscopic-restart`), starting from a simulation with constant, radial or linear weighting. The process is enabled by

    Part-Weight-Type = cell_local

The adaption is based on multiple criteria. If a quality factor {ref}`sec:DSMC-quality` is not resolved in a cell, the weighting factor is lowered. For a 3D case, the following equation is used to set the bound of the weight.

$$w < \frac{1}{\left(\sqrt{2}\pi d_{\mathrm{ref}}^2 n^{2/3}\right)^3}$$

The threshold for the adaption is set by

    Part-Weight-CellLocal-QualityFactor = 0.8

which in case of DSMC corresponds to the mean collision separation distance over mean free path (DSMC_MCS_over_MFP) and for BGK/FP to the maximal relaxation factor (BGK_MaxRelaxationFactor/FP_MaxRelaxationFactor). If the read-in cell values from the DSMCState are above this threshold the weighting factor will be adapted to resolve it. For DSMC, the scaling factor is put to the power of 3 for 3D and 2 for 2D to increase the number of particles and thus the resolution of the mean collision separation distance accordingly.

If all quality factors are resolved, the weight is adapted in such a way, that the simulation particle numbers stay in a predefined range.

    Part-Weight-CellLocal-MinParticleNumber = 10
    Part-Weight-CellLocal-MaxParticleNumber = 100

A further refinement, with a higher minimum particle number can be additionally given close to the symmetry axis in axisymmetric simulations (determined by using 5% of the maximum y-coordinate of the domain) or close to catalytic boundaries.

    Part-Weight-CellLocal-SymAxis-MinPartNum = 200
    Part-Weight-CellLocal-Cat-MinPartNum = 200

For coupled BGK/FP-DSMC simulations, additional parameters are available to scale the weight in the BGK/FP regions respectively, where the other criteria shown above are resolved. The following parameters are directly multiplied by read-in weighting factor:

    Part-Weight-CellLocal-RefineFactorBGK = 10
    Part-Weight-CellLocal-RefineFactorFP  = 10

These values would correspond to an increase of the weighting factor of 10 and thus a reduction of the particle number in the BGK/FP regions. In addition, the upper limit of the weights can be disabled for simulations in which the weight should only be lowered.

    Part-Weight-CellLocal-IncludeMaxPartNum = F

To avoid large jumps in the optimal weight between neighboring cells, a median filtering routine can be applied, together with the number of times that the algorithm should be called.

    Part-Weight-CellLocal-ApplyMedianFilter = T
    Part-Weight-CellLocal-RefinementNumber = 2

If a restart should be performed from an already adapted simulation, without further optimization of the weights, the following flag can be set in the input file.

    Part-Weight-CellLocal-SkipAdaption = T

In this case, the weight distribution from the given restart files is read in and used. As for the linear weighting in 3D, particles are cloned or deleted with a given probability when moving to a new cell.

(sec:split-merge)=
### Split-And-Merge

Variable particle weighting based on a split and merge algorithm is currently supported for PIC (with and without background gas) or a background gas (an additional trace species feature is described in Section {ref}`sec:background-gas`). The general functionality can be enabled with the following flag:

    Part-vMPF                           = T

The split and merge algorithm is called at the end of every time step. In order to manipulate the number of particles per species per cell, merge and split thresholds can be defined as is shown in the following.

    Part-Species2-vMPFMergeThreshold    = 100

The merge routine randomly deletes particles until the desired number of particles is reached and the weighting factor is adopted accordingly. Afterwards, the particle velocities $v_i$ are scaled in order to ensure momentum and energy conservation with

$$ \alpha = \frac{E^{\mathrm{old}}_{\mathrm{trans}}}{E^{\mathrm{new}}_{\mathrm{trans}}},$$
$$ v^{\mathrm{new}}_{i} = v_{\mathrm{bulk}} + \alpha (v_{i} - v^{\mathrm{new}}_{\mathrm{bulk}}).$$

Internal degrees of freedom are conserved analogously. In the case of quantized energy treatment (for vibrational and electronic excitation), energy is only conserved over time, where the energy difference (per species and energy type) in each time step due to quantized energy levels is stored and accounted for in the next merge process.

    Part-Species2-vMPFSplitThreshold    = 10

The split routine clones particles until the desired number of particles is reached.
The algorithm randomly chooses a particle, clones it and assigns the half weighting factor to both particles.
If the resulting weighting factor would drop below a specific limit that is defined by

    Part-vMPFSplitLimit = 1.0 ! default value is 1

the split is stopped and the desired number of particles will not be reached.
The basic functionality of both routines is verified during the nightly regression testing in `piclas/regressioncheck/NIG_code_analyze/vMPF_SplitAndMerge_Reservoir`.
It is important to note that the variable `Part-SpeciesX-vMPFMergeThreshold` might be initially ignored, when using cell-local particle emission and `Part-SpeciesX-vMPFSplitThreshold` is set, as
described in Section {ref}`sec:particle-cell-local`. In this case, the split threshold variable will be utilized as the target number of particles per cell during the insertion and the weighting factor will be determined from the given number density and cell volume.

## Virtual Cell Merge

In the case of very complex geometries, very small cells can sometimes be created to represent the geometry, in which, however, only very few particles are present. This results in a very large statistical noise in these cells and the collision process is only inadequately covered due to the poor statistics. In this case, cells can be virtually merged with neighbouring cells during the runtime using predefined parameters. The merged cells are then treated as one large cell. This means that DSMC collision pairings and the BGK relaxation process are performed with all particles within the merged cell, i.e. all particles of the individual cells. The averaging then also takes place for all particles in the merged cell and all individual cells then receive the values of the merged cell in the output.
The virtual cell merge can be enabled with the following flag:

    Part-DoVirtualCellMerge             = T

Currently, only merging based on the number of particles within the cell is implemented. The minimum number of particles below which the cell is merged is defined with the following parameter:

    Part-MinPartNumCellMerge            = 3

Furthermore, the spread or aggressiveness of the merge algorithm can be changed, i.e. how deep the merge extends into the mesh starting from each cell. 0 is the least aggressive merge, 3 the most aggressive merge.

    Part-CellMergeSpread                = 0

There is also the possibility to define a maximum number of cells that can be merged. In this way, a desired "resolution" of the virtual cells can be achieved.

    Part-MaxNumbCellsMerge              = 5

