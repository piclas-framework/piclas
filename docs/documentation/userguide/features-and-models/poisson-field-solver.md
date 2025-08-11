(sec:poisson)=
# Field Solver - Poisson Equation
To numerically solve electrostatic problems, PICLas offers a solver for Poisson's equation

$$\nabla \cdot \textbf{E} = \frac{\rho}{\varepsilon_{0}}$$

where the charge density $\rho$ is given by the charged particles within the system.
To enable a simulation using the Poisson field solver, the code must be compiled using the following setting

    PICLAS_EQNSYSNAME     = poisson
    PICLAS_TIMEDISCMETHOD = Euler-Explicit, Leapfrog, Boris-Leapfrog, Higuera-Cary, RK3, RK4

where one of the available time discretization methods for the particle evolution must be chosen for `PICLAS_TIMEDISCMETHOD`, which
are also referred to as particle pushers.
There are different options available that yield different functionalities.

| `PICLAS_TIMEDISCMETHOD` | Electric field | Magnetic field | Order of Convergence |
| ----------------------: | -------------- | -------------- | -------------------- |
|        `Euler-Explicit` | yes            | no             | 1                    |
|              `Leapfrog` | yes            | no             | 2                    |
|        `Boris-Leapfrog` | yes            | yes            | 2                    |
|          `Higuera-Cary` | yes            | yes            | 2                    |
|                   `RK3` | yes            | yes            | 3                    |
|                   `RK4` | yes            | yes            | 4                    |

Note that high-order time discretization methods in general allow for a larger time step and are usually more costly per time step.
To see the available parameter input file options, simply run

    ./bin/piclas --help HDG

## PETSc Solver
A multitude of different numerical methods to solve the resulting system of linear equations is given by the implemented PETSc
library {cite}`petsc-web-page`, {cite}`petsc-user-ref`, {cite}`petsc-efficient`. For detailed installation steps of PETSc within PICLas, see Section {ref}`sec:petsc-installation`.
To use PETSc, another flag must be set during the compilation of PICLas

    LIBS_USE_PETSC = ON

and the parameter input file for the simulation requires setting

    PrecondType = 1
    epsCG       = 1e-6      ! default = 1e-50
    MaxIterCG   = 10000     ! default

where the following options are possible:

| `PrecondType` | Preconditioner                                                                                                                                                                                                                                                                                             | Solver                                                                            |
| ------------: | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------- |
|           `1` | [Block Jacobi](https://petsc.org/main/manualpages/PC/PCBJACOBI/)                                                                                                                                                                                                                                           | [Conjugate gradient](https://petsc.org/main/manualpages/KSP/KSPCG/)               |
|           `2` | [Block Jacobi](https://petsc.org/main/manualpages/PC/PCBJACOBI/)                                                                                                                                                                                                                                           | [Pipelined conjugate gradient](https://petsc.org/main/manualpages/KSP/KSPPIPECG/) |
|           `3` | [BoomerAMG](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html) (requires [Hypre](https://hypre.readthedocs.io/en/latest/index.html), [configuration](https://mooseframework.inl.gov/application_usage/hypre.html)), Fallback: [Block Jacobi](https://petsc.org/main/manualpages/PC/PCBJACOBI/) | [Flexible GMRES](https://petsc.org/main/manualpages/KSP/KSPFGMRES/)               |
|          `10` | [Direct Cholesky](https://petsc.org/main/manualpages/PC/PCCHOLESKY/) (requires [MUMPS](https://mumps-solver.org/index.php))                                                                                                                                                                                | -                                                                                 |

The parameters `epsCG` and `MaxIterCG` correspond to the absolute tolerance (`atol`) and the maximum number of iterations (`maxits`) as described in the [PETSc manual](https://petsc.org/main/manual/ksp/#convergence-tests). In general, it is recommended to the use the direct solver `PrecondType = 10` (if the case fits into the memory) or the preconditioner iterative solver combination `PrecondType = 1`. When a large number of computational nodes is used, `PrecondType = 2` can offer a performance benefit, while `PrecondType = 3` might offer better convergence properties for large, complex cases.

Other potential combinations and/or specific preconditioner/solver options can be easily implemented in the routine [`PETScSetSolver`](https://github.com/piclas-framework/piclas/blob/master/src/hdg/hdg.f90). An overview over the available options is available in the [Summary of Sparse Linear Solvers Available In PETSc](https://petsc.org/main/overview/linear_solve_table/).

## Internal CG Solver
The formerly used numerical method for solving the resulting system of linear equations, is the Conjugate Gradient Method. The following
input parameter can be set to control the simulation

    epsCG     = 1e-6 ! default value for the abort residual
    MaxIterCG = 500  ! default value for the number of CG solver iterations

where `epsCG` is the residual of the CG solver and `MaxIterCG` are the maximum number of iteration performed in one time step to
solve the system.
Furthermore, the residual can be either set absolute (default) or relative via

    useRelativeAbortCrit = F ! default

Note that the same parameter setting for `epsCG` will result in a higher residual without the PETSC library as compared with PETSC solver.

## Symmetric Simulations

The Poisson solver is suitable for 3D, 2D/1D, and axially symmetric 2D simulations, but the field solver generally works in 3D Cartesian coordinates when no specific settings are selected.
For 2D axisymmetric simulations, the solver can be switched from $x,y$ coordinates to $x,r$ by setting the parameter `Particles-Symmetry2DAxisymmetric`.
Note that the mesh should be arranged in a way that $x$ is the axial direction, $y$ is the radial $r$ direction and, $z$ is the third symmetric or periodic dimension.

| `Particles-Symmetry2DAxisymmetric` | **Coordinates** | **Equation**                                                                                                          |
| ---------------------------------: | --------------- | --------------------------------------------------------------------------------------------------------------------- |
|                                `F` | $x,y$           | $\dfrac{d^2\Phi}{dx^2}+\dfrac{d^2\Phi}{dy^2} = \dfrac{\rho(x,y)}{\varepsilon_{0}}$                                    |
|                                `T` | $x,r$           | $\dfrac{d^2\Phi}{dx^2}+\dfrac{1}{r}\dfrac{d}{dr}\left(r\dfrac{d\Phi}{dr}\right) = \dfrac{\rho(x,r)}{\varepsilon_{0}}$ |

When particles are present in the simulation, an additional parameter `Particles-Symmetry-Order` (default is 3) sets the
dimensionality for particle behaviour

| `Particles-Symmetry-Order` | **Coordinates** | `Particles-Symmetry2DAxisymmetric` |
| -------------------------: | --------------- | ---------------------------------- |
|                          1 | 1D              | `F`                                |
|                          2 | 2D              | `F` ($x,y$), `T` ($x,r$)           |
|                          3 | 3D (default)    | `F`                                |

There are different deposition methods available and their compatibility with the parameters `Particles-Symmetry-Order` and `Particles-Symmetry2DAxisymmetric` is described in
Section {ref}`userguide/features-and-models/PIC:Compatibility with Field Solver Symmetry Settings`

Examples: Use the following parameter for these types of simulations

    ! 2D Cartesian
    Particles-Symmetry-Order = 2

    ! 2D Axisymmetrc
    Particles-Symmetry-Order = 2
    Particles-Symmetry2DAxisymmetric = T

    ! 1D
    Particles-Symmetry-Order = 1

