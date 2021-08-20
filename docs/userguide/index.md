% ```{include} ../../README.md
% ---
% relative-docs: docs/
% relative-images:
% ---
% ```

# Welcome to PICLas Testing's documentation!

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
maxdepth: 1
caption: User Guide
numbered:
---
introduction.md
installation.md
meshing.md
workflow.md
features-and-models/index.md
visu_output.md
tools.md
tutorials/index.md
cluster_guide.md
appendix.md
references.md
```








