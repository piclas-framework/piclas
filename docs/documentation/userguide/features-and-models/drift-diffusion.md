(sec:drift-diffusion)=
# Drift-diffusion model

The drift-diffusion scheme {cite}`Dujko2020` can be used to model the plasma as a continuum instead of kinetic particles. In PICLas, it is utilized to model electrons as a fluid, whereas ions are treated kinetically, resulting in a hybrid model. The executable to use this model can be created by

    mkdir build_electron_fluid && cd build_electron_fluid
    cmake .. -DPICLAS_EQNSYSNAME=drift_diffusion -DPICLAS_TIMEDISCMETHOD=Explicit-FV -DLIBS_USE_PETSC=ON

A tutorial, where the model is applied to a streamer formation, can be found here: {ref}`sec:tutorial-streamer`.

TODO: add different time integration possibilities for field (FV solver for electrons) and charged heavy species (particle push +
HDG solver for electric fields)

## First-order fluid model

The first-order-fluid model is derived from the Boltzmann equation. Continuity and the balance of momentum equations are used and the set is truncated at the momentum balance equation {cite}`Dujko2020JPD_I`. The first-order-fluid model considers only the first two balance laws from the system. For electrons and ions, it reads as {cite}`Markosyan2015`:

$$
\frac{\partial n}{\partial t} = \nabla \cdot \left( \mu(E) E n \right) + D(E) \cdot \nabla n + \nu_I(E, t) \tag{1}
$$

$$
\frac{\partial n_{\text{ion}}}{\partial t} = \nu_I(E, t) \tag{2}
$$

Where $E = |E|$, $n_{\text{ion}}$ is the ion density, and where mobility $\mu$, diffusion $D$ and $\nu_I$ are functions of the local electric field.

Coupling these two equations with the Poisson equation results in the following form {cite}`Markosyan2015`:

$$
\frac{\partial E}{\partial x} = -\varepsilon_0 e n_{\text{ion}} \tag{3}
$$

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

TODO: Introduce new subroutine, similar to AddBRElectronFluidToPartSource() but for the drift-diffusion electrons.
See CalcSourceHDG() on how the drift-diffusion is added to the source terms.
