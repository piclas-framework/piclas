(sec:drift-diffusion)=
# Drift-diffusion model

The drift-diffusion scheme {cite}`Dujko2020` can be used to model the plasma as a continuum instead of kinetic particles. In PICLas, it is utilized to model electrons as a fluid, whereas ions are treated kinetically with PIC, resulting in a hybrid model. The executable to use this model can be created by

    mkdir build_electron_fluid && cd build_electron_fluid
    cmake .. -DPICLAS_EQNSYSNAME=drift_diffusion -DPICLAS_TIMEDISCMETHOD=Explicit-FV -DLIBS_USE_PETSC=ON

A tutorial, where the model is applied to a streamer formation, can be found here: {ref}`sec:tutorial-streamer`.

For the fluid model, only the first conservation law, namely the continuity equation, is considered (first-order model) {cite}`Markosyan2015`. Accordingly, the change in electron density is determined as follows:

$$
\frac{\partial n}{\partial t} = \nabla \cdot \left( \mu(E) E n \right) + D(E) \cdot \nabla n + \nu_I(E, t) \tag{1}
$$

where the electron mobility $\mu$ and diffusion coefficient $D$ are functions of the local electric field.
The change in ion density $n_{\text{ion}}$ is determined consistently based on the ionization rate $\nu_I$, and a corresponding number of ions is introduced as particles:

$$
\frac{\partial n_{\text{ion}}}{\partial t} = \nu_I(E, t) \tag{2}
$$

The electric field is then calculated in PICLas — as in a classical PIC simulation — by solving the Poisson equation.
However, in this hybrid approach, the electron density used in the Poisson equation is taken from the solution of the above drift-diffusion equation.

## Finite volume solver

The drift-diffusion equation is solved using a second-order finite volume solver based on flux reconstruction at the element interfaces. For this linear reconstruction, gradients of the current solution (here the electron density) are computed at each time step from direct neighbour elements in a least-squares approach.

In order to avoid nonphysical oscillations, gradient limiters can be applied by choosing the value of `Grad-LimiterType`.


| **Grad-LimiterType** |                                             **Description**                                              |
| :------------------: | :------------------------------------------------------------------------------------------------------: |
|          0           |                       Gradients are set to 0. Emulates first-order finite volumes.                       |
|          1           |                        Barth-Jespersen (a.k.a. minmax) limiter {cite}`Barth1989`.                        |
|          4           | Venkatakrishnan limiter {cite}`Venkatakrishnan1995`. Additional parameter K can be set with `Grad-VktK`. |
|          9           |                                               No limiter.                                                |



## Boundary conditions

Overview of the boundary conditions for electron density

    BoundaryType-FV = (/2,1/) ! 2: Dirichlet, 1: number of the RefState (0 to use ExactFunc)
    BoundaryType-FV = (/3,0/) ! 3: Neumann

## Plasma chemistry

The diffusion coefficient is considered as a scalar in the model. However, in the case of streamer discharges, the diffusion tensor can be anisotropic. This means that diffusion can occur at different rates in different directions, which would require a tensor instead of a single scalar diffusion coefficient {cite}`Markosyan2015`.

Uses *Diffusion-Coefficients* in the SpeciesDatabase.h5 file.

Select different models via

    Part-Species1-SpeciesName = N2
    BGGas-DriftDiff-Database  = Phelps

Part-Species1-SpeciesName = N2: This line indicates that the species to be used in the simulation is nitrogen gas (N2).

BGGas-DriftDiff-Database = Phelps: This line indicates that the database used for diffusion and drift coefficients is selected according to the Phelps model.

## Visualisation

Convert state file with *piclas2vtk* to view charge $\rho$ and current density $j$ as well as the electric potential $\Phi$ and
field strengths $E$ in the *Solution.vtu* file.

Furthermore, the electron number density $n_{e}$ is written to *ElemData.vtu* and is labelled *ElectronDensityCell*.
