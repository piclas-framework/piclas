(sec:discrete-velocity-method)=
# Discrete velocity method

The discrete velocity method (DVM) is a numerical approach to solve the BGK approximation of the Boltzmann equation by discretizing the velocity space into a finite set of velocities. This allows for the simulation of rarefied gas or plasma dynamics in a completely deterministic manner.
The DVM module can be enabled by compiling PICLas with the following parameters


    PICLAS_EQNSYSNAME = discrete_velocity
    PICLAS_PARTICLES = OFF

DVM for neutral gas should be used with the exponential differencing method (ED-DVM), consisting in a finite volume method with second-order accuracy in time and space {cite}`Garmirian2025`, by setting

    PICLAS_TIMEDISCMETHOD = ED-DVM

Simulating plasma, with a combination of DVM for particle flow and HDG for the electric potential, is also possible with

    PICLAS_TIMEDISCMETHOD = PLOESMA

The different species parameters can be retrieved from the species database (see {ref}`sec:unified-species-database`). For example:

    DVM-Species-Database = SpeciesDatabase.h5
    DVM-nSpecies = 2
    DVM-Species1-SpeciesName = Ar
    DVM-Species2-SpeciesName = He

For each species, the velocity discretization can be specified with

    DVM-Species1-VeloDiscretization

which should be set to 3 for near-equilibrium simulations (Gauss-Hermite quadrature) or to 2 for strong non-equilibrium (uniform grid). For the uniform grid, the limits of the grid (in m/s) can be set in each direction $(v_x,v_y,v_z)$ with

    DVM-Species1-VeloMin = (/-2000.,-1500.,-1500./)
    DVM-Species1-VeloMax = (/2000.,1500.,1500./)

For the Gauss-Hermite quadrature, only the reference translational temperature, which should be close to the actual flow temperature, has to be set in every direction:

    DVM-Species1-GaussHermiteTemp = (300.,300.,300.)

The number of discrete velocities in each direction can be set with

    DVM-Species1-nVeloX = 25
    DVM-Species1-nVeloY = 15
    DVM-Species1-nVeloZ = 15

If the flow has symmetry properties, the number of discrete velocities can be reduced by only specifying the number of velocities for directions where it needs to be discretized. In this case, a dimension parameter has to be set accordingly:

    DVM-Species1-Dimension = 1 or 2

Several BGK collision models are available through

    DVM-Collisions = T
    DVM-BGKCollModel = 1 to 6


## Finite volume solver

The DVM module uses a second-order finite volume solver based on flux reconstruction at the element interfaces.
In order to avoid nonphysical oscillations, gradient limiters can be applied by choosing the value of `Grad-LimiterType` (see {ref}`sec:drift-diffusion`).

## Initial and boundary conditions

Initial and boundary conditions are set using `RefState-FV`. Their content is, in order: density $n$, velocity in three directions $(u_x,u_y,u_z)$, temperature $T$, six components of traceless pressure tensor $(P_{xx},P_{yy},P_{zz},P_{xy},P_{xz},P_{yz})$, heat flux $(q_x,q_y,q_z)$. The number of conditions per species is given by `IniRefState-FV`. For example, for two species initialized uniformly at equilibrium in a domain with diffuse walls:

    IniExactFunc-FV  = 1  ! uniform
    IniRefState-FV   = 3  ! per species

    ! initial conditions (RefState number 1)
    RefState-FV =(/6.5E19, 0, 0, 0, 273., 0, 0, 0,0,0,0,0,0,0/)
    RefState-FV =(/13.E19, 0, 0, 0, 273., 0, 0, 0,0,0,0,0,0,0/)

    BoundaryName = BC_x+
    BoundaryType-FV = (/4,2/) ! diffuse wall using RefState number 2
    RefState-FV =(/1., 0, 350., 0, 273., 0, 0, 0,0,0,0,0,0,0/)
    RefState-FV =(/1., 0, 350., 0, 273., 0, 0, 0,0,0,0,0,0,0/)

    BoundaryName = BC_x-
    BoundaryType-FV = (/4,3/) ! diffuse wall using  RefState number 3
    RefState-FV =(/1., 0, -350., 0, 273., 0, 0, 0,0,0,0,0,0,0/)
    RefState-FV =(/1., 0, -350., 0, 273., 0, 0, 0,0,0,0,0,0,0/)

For diffuse walls, only velocity and temperature actually matter, the other values are dummy variables.

## Visualisation

Convert state file with *piclas2vtk* to view macroscopic flow values in the *DVM_Solution.vtu* file.
