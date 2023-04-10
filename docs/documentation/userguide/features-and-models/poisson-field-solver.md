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

|  `PICLAS_TIMEDISCMETHOD` |   Electric field  |  Magnetic field  | Order of Convergence |
| -----------------------: | ----------------- | ---------------- |   -----------------  |
|     `Euler-Explicit`     |        yes        |        no        |           1          |
|        `Leapfrog`        |        yes        |        no        |           2          |
|     `Boris-Leapfrog`     |        yes        |        yes       |           2          |
|      `Higuera-Cary`      |        yes        |        yes       |           2          |
|           `RK3`          |        yes        |        yes       |           3          |
|           `RK4`          |        yes        |        yes       |           4          |

Note that high-order time discretization methods in general allow for a larger time step and are usually more costly per time step.
To see the available parameter input file options, simply run

    ./bin/piclas --help HDG

## CG Solver
The default numerical method for solving the resulting system of linear equations, is the Conjugate Gradient Method. The following
input parameter can be set to control the simulation

    epsCG     = 1e-6 ! default value for the abort residual
    MaxIterCG = 500  ! default value for the number of CG solver iterations

where `epsCG` is the residual of the CG solver and `MaxIterCG` are the maximum number of iteration performed in one time step to
solve the system.
Furthermore, the residual can be either set absolute (default) or relative via

    useRelativeAbortCrit = F ! default

## PETSc Solver
A multitude of different numerical methods to solve the resulting system of linear equations is given by the implemented PETSc
library {cite}`petsc-web-page`, {cite}`petsc-user-ref`, {cite}`petsc-efficient`. For detailed installation steps of PETSc within PICLas, see Section {ref}`sec:petsc-installation`.
To use PETSc, another flag must be set during the compilation of PICLas

    PICLAS_PETSC = ON

and the parameter input file for the simulation requires setting

    PrecondType = 2 ! default

where the following options are possible

|  `PrecondType`  | Iterative or Direct |      Method      |
| --------------: |  -----------------  | ---------------- |
|       `1`       |      iterative      |  Krylov subspace |
|       `2`       |      iterative      |  Krylov subspace |
|       `3`       |      iterative      |  Krylov subspace |
|       `10`      |        direct       |                  |

Note that the same parameter setting for `epsCG` will result in a smaller residual with PETSc as compared with the default CG solver
without using the PETSc library.


