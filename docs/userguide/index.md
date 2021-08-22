% ```{include} ../../README.md
% ---
% relative-docs: docs/
% relative-images:
% ---
% ```

# Welcome to the PICLas Documentation!

![alt](../logo.png)

[**PICLas**](https://github.com/piclas-framework/piclas)  is a three-dimensional simulation
framework for Particle-in-Cell, Direct Simulation Monte Carlo and other particle methods that can be coupled for
the simulation of collisional plasma flows.
It features a high-order discontinuous Galerkin (DG) simulation module for the solution of the time-dependent Maxwell equations on
unstructured hexahedral elements in three space dimensions.
The code was specifically designed for very high order accurate simulations on massively parallel systems.
It is licensed under GPLv3, written in Fortran and parallelized with MPI.

```{toctree}
---
maxdepth: 2
caption: User Guide
---
introduction.md
```

```{toctree}
---
maxdepth: 2
caption: Developer Guide
---
developerguide/index.md
```

```{toctree}
---
maxdepth: 1
caption: References
---
references.md
```

