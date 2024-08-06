# Drift-diffusion model
Streamers are growing ionized fingers that appear when ionizable matter is suddenly exposed to high voltages [1]. Streamers occur in lightning and sprites as well as in industrial applications such as lighting, treatment of polluted gases and water, disinfection plasma jets and bullets and plasma-assisted combustion. Further optimization and understanding of such applications depend on an accurate knowledge of the electron dynamics during streamer development [2].

Streamer discharge simulations are basically modeled in three ways. These are the Particle Model, Fluid Model, and Hybrid Model.
In this section, the Hybrid Model which combines the 1st order fluid model with particle model approach is explained. In the Hybrid Model Simulations, ions are represented by particles.

The first-order fluid model is the multiply used reaction–drift-diffusion model[1]. It is also called 'cllassical fluid model´. The drift-diffusion scheme can be used to model the electrons species in the plasma as a continuum instead of kinetic particles.

The piclas executable to use the model in created by

    mkdir build_electron_fluid && cd build_electron_fluid
    cmake .. -DPICLAS_EQNSYSNAME=drift_diffusion -DPICLAS_TIMEDISCMETHOD=Explicit-FV -D PICLAS_PETSC=ON

TODO: add different time integration possibilities for field (FV solver for electrons) and charged heavy species (particle push +
HDG solver for electric fields)

**The First-Order Fluid Model** 

The first-order-fluid model is derived from the Boltzmann equation. 
Continuity and the balance of momentum equations are used and the set is truncated at the momentum balance equation [3].

 
$$
\frac{\partial n}{\partial t} + \nabla \cdot (nv) = C1  
$$

$$
\frac{\partial (nmv)}{\partial t} + \nabla \cdot (nm\langle cc \rangle) - neE = C_{mc}  
$$

After some implementation and the simplfications are applied to the these constitutive equations, the first-order-fluid model is given by the equations : 

The equation is given by:

$$
\frac{\partial n}{\partial t} = \nabla \cdot (D \cdot \nabla n) + \nabla \cdot (\mu n \mathbf{E}) + n (\nu_I - \nu_A)
$$


$$
\frac{\partial n_p}{\partial t} = n \nu_I  
$$

$$
\frac{\partial n_n}{\partial t} = n \nu_A  
$$

coupled to the Poisson equation for the electric field, 

$$
\nabla^2 \phi = \frac{e}{\epsilon_0} (n - n_p + n_n), \quad E = -\nabla \phi
$$

Here, $$ n_p $$ and $$ n_n $$ represent the positive and negative ion densities, respectively, while `$$ \nu_I $$ and $$ \nu_A $$ denote the ionization and attachment collision frequencies due to electron–molecule collisions. Additionally, $$ \phi $$ stands for the electric potential.

As a summary, for electrons and ions equations read as: 

$$
\frac{\partial n}{\partial t} = \nabla \cdot (\mu(E)n \mathbf{E} + D(E) \cdot \nabla n) + n \nu_I(E)
$$

$$
\frac{\partial n_{\text{ion}}}{\partial t} = n \nu_I(E)
$$

## Boundary conditions

Overview of boundary conditions

    BoundaryType-FV = (/3,0/) ! 3: open BC
    BoundaryType-FV = (/4,0/) ! 4: Neumann, 0: electron density
    BoundaryType-FV = (/4,1/) ! 4: Neumann, 1: number of the RefState

TODO: implement RefState (the same as HDG solver uses for the electric potential BC) + add new reggie

## Plasma chemistry
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

[2] 

[3] Dujko, S., Markosyan, A. H., White, R. D., & Ebert, U. (2020). High-order fluid model for streamer discharges: I. Derivation of model and transport data. Journal of Physics D: Applied Physics, 53(33), 335201. https://doi.org/10.1088/1361-6463/ab8b92