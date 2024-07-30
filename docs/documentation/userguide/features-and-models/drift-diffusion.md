# Drift-diffusion model

Streamers occur in nature (lightning and sprites, and technology (plasma-assisted combustion, plasma medicine, disinfection etc). 
Streamer models can be categorized into four types: 

**- Particle Models**

These models are typically of the particle-in-cell (PIC)type. With these models, a large number of particles is followed as they move through the simulation domain, so thatone has direct information about the particle distributionin phase space.

**- Kinetic Models**

The second type of models are the so-called kinetic models, that couple the full Boltzmann equation with the Poisson equation. Such models are computationally very costly, because they require a numerical grid that covers the full phase space. However, advances in computing power and algorithms have made some of these fully kinetic simulations possible.

**- Fluid Models**

The third type of models are the plasma fluid models,which describe the electron dynamics in plasma based on macroscopic quantities like electron density, average elec-tron velocity, average electron energy etc.

**- Hybrid Models**

Hybrid models combine the strengths of the fluidand particle approaches, by combining the fast speed of fluidsimulations with the accurate particle kinetics of particlemodels [1].

The focus of this section is on hybrid modelling approach. The main difference between the hybrid model and the 1st order fluid model is that, in the hybrid model, ions are represented by particles, not by constants.


The drift-diffusion scheme can be used to model the electrons species in the plasma as a continuum instead of kinetic particles.
The piclas executable to use the model in created by

    mkdir build_electron_fluid && cd build_electron_fluid
    cmake .. -DPICLAS_EQNSYSNAME=drift_diffusion -DPICLAS_TIMEDISCMETHOD=Explicit-FV -D PICLAS_PETSC=ON

TODO: add different time integration possibilities for field (FV solver for electrons) and charged heavy species (particle push +
HDG solver for electric fields)

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

### Visualisation

Convert state file with `piclas2vtk` to view charge $\rho$ and current density $j$ as well as the electric potential $\Phi$ and
field strengths $E$ in the *Solution.vtu* file.

TODO: Introduce new subroutine, similar to AddBRElectronFluidToPartSource() but for the drift-diffusion electrons.
See CalcSourceHDG() on how the drift-diffusion is added to the source terms.

Furthermore, the electron number density $n_{e}$ is written to *ElemData.vtu* and is labelled *Electron\_Density*.

TODO: combine *Electron\_Density* from the drift-diffusion output with the *ElectronDensityCell* output that is used for PIC and
PIC-BR. See CalcSourceHDG() on how the drift-diffusion is added to the source terms.
