# Drift-diffusion model
Streamers are growing ionized fingers that appear when ionizable matter is suddenly exposed to high voltages. Streamers occur in lightning and sprites as well as in industrial applications such as lighting, treatment of polluted gases and water, disinfection plasma jets and bullets and plasma-assisted combustion. Further optimization and understanding of such applications depend on an accurate knowledge of the electron dynamics during streamer development \citep{Dujko2020JPD}.

Streamer discharge simulations are basically modeled in three ways. These are the Particle Model, Fluid Model, and Hybrid Model.
In this section, the Hybrid Model which combines the 1st order fluid model with particle model approach is explained. In the Hybrid Model Simulations, ions are represented by particles.

The first-order fluid model is the multiply used reaction–drift-diffusion model[1]. It is also called 'classical fluid model´. The drift-diffusion scheme can be used to model the electrons species in the plasma as a continuum instead of kinetic particles.

The piclas executable to use the model in created by

    mkdir build_electron_fluid && cd build_electron_fluid
    cmake .. -DPICLAS_EQNSYSNAME=drift_diffusion -DPICLAS_TIMEDISCMETHOD=Explicit-FV -D PICLAS_PETSC=ON

TODO: add different time integration possibilities for field (FV solver for electrons) and charged heavy species (particle push +
HDG solver for electric fields)

**The First-Order Fluid Model** 

The first-order-fluid model is derived from the Boltzmann equation. Continuity and the balance of momentum equations are used and the set is truncated at the momentum balance equation [3]. 

The first-order-fluid model (Classical Model) considers only the first two balance laws from the system. For electrons and ions, it reads as[2]:

$$
\frac{\partial n}{\partial t} = \nabla \cdot \left( \mu(E) E n \right) + D(E) \cdot \nabla n + \nu_I(E, t) \tag{1}
$$

$$
\frac{\partial n_{\text{ion}}}{\partial t} = \nu_I(E, t) \tag{2}
$$

\text{where } E = |E|, \, n_{\text{ion}} \text{ is the ion density, and where mobility } \mu, \text{ diffusion } D, \text{ and } \nu_I \text{ are functions of the local electric field.}




Coupling these two equations with the Poisson equation results in the following form[2]: 

$$
\frac{\partial E}{\partial x} = -\varepsilon_0 e n_{\text{ion}} \tag{3}
$$

## Boundary conditions

Overview of boundary conditions

    BoundaryType-FV = (/3,0/) ! 3: open BC
    BoundaryType-FV = (/4,0/) ! 4: Neumann, 0: electron density
    BoundaryType-FV = (/4,1/) ! 4: Neumann, 1: number of the RefState

TODO: implement RefState (the same as HDG solver uses for the electric potential BC) + add new reggie

To create steady propagation conditions for the negative front, the electric feld on the left boundary x = 0 is fxed to the time 
independent value E0.

$$
E(0, t) = E_0 > 0. \tag{4}
$$
 

The electric feld for x > 0 is calculated by integrating the first order fluid model equation (3) numerically over x, with the (4) as a boundary condiition. The system length _L_, the number of grid points, and the grid spacing are adjusted as boundary conditions.

## Plasma chemistry

The diffusion coefficient is considered as a scalar in the model. However, in the case of streamer discharges, the diffusion tensor can be anisotropic. This means that diffusion can occur at different rates in different directions, which would require a tensor instead of a single scalar diffusion coefficient [2].

Uses *Diffusion-Coefficients* in the SpeciesDatabase.h5 file.

Select different models via

    Part-Species1-SpeciesName = N2
    BGGas-DriftDiff-Database  = Phelps



Part-Species1-SpeciesName = N2: This line indicates that the species to be used in the simulation is nitrogen gas (N2).

BGGas-DriftDiff-Database = Phelps: This line indicates that the database used for diffusion and drift coefficients is selected according to the Phelps model.

### Visualisation 

Convert state file with `piclas2vtk` to view charge $\rho$ and current density $j$ as well as the electric potential $\Phi$ and
field strengths $E$ in the *Solution.vtu* file.

TODO: Introduce new subroutine, similar to AddBRElectronFluidToPartSource() but for the drift-diffusion electrons.
See CalcSourceHDG() on how the drift-diffusion is added to the source terms.

Furthermore, the electron number density $n_{e}$ is written to *ElemData.vtu* and is labelled *Electron\_Density*.

TODO: combine *Electron\_Density* from the drift-diffusion output with the *ElectronDensityCell* output that is used for PIC and
PIC-BR. See CalcSourceHDG() on how the drift-diffusion is added to the source terms.

### References 
[1] Dujko, S., Markosyan, A. H., White, R. D., & Ebert, U. (2020). High-order fluid model for streamer discharges: II. Numerical solution and investigation of planar fronts. Journal of Physics D: Applied Physics, 53(33), 335202. https://doi.org/10.1088/1361-6463/ab8b93

[2] Markosyan, A. H., Teunissen, J., Dujko, S., & Ebert, U. (2015). Comparing plasma fluid models of different order for 1D streamer ionization fronts. Journal of Physics D: Applied Physics. Received 23 March 2015, revised 6 August 2015, accepted for publication 7 September 2015, published 8 October 2015.

[3] Dujko, S., Markosyan, A. H., White, R. D., & Ebert, U. (2020). High-order fluid model for streamer discharges: I. Derivation of model and transport data. Journal of Physics D: Applied Physics, 53(33), 335201. https://doi.org/10.1088/1361-6463/ab8b92

