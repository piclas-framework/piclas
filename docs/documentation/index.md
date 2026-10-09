# Welcome to the PICLas Documentation!

[**PICLas**](https://github.com/piclas-framework/piclas) is a three-dimensional simulation
framework for Particle-in-Cell, Direct Simulation Monte Carlo and other particle methods that can be coupled for
the simulation of collisional plasma flows.
It features high-order discontinuous Galerkin (DG) and hybridizable discontinuous Galerkin (HDG) simulation modules for the solution of the time-dependent Maxwell equations and electrostatic Poisson equation on
unstructured hexahedral elements in three space dimensions.
The code was specifically designed for very high order accurate simulations on massively parallel systems.
It is licensed under GPLv3, written in Fortran and parallelized with MPI.

::::{grid} 5
:gutter: 2

:::{grid-item}
:class: sd-text-center
```{button-ref} code_of_conduct
:ref-type: doc
:color: primary
:shadow:
Code of Conduct
```
:::

:::{grid-item}
:class: sd-text-center
```{button-ref} acknowledge
:ref-type: doc
:color: primary
:shadow:
Acknowledge PICLas
```
:::

:::{grid-item}
:class: sd-text-center
```{button-ref} userguide/index
:ref-type: doc
:color: primary
:shadow:
User Guide
```
:::

:::{grid-item}
:class: sd-text-center
```{button-ref} developerguide/index
:ref-type: doc
:color: primary
:shadow:
Developer Guide
```
:::

:::{grid-item}
:class: sd-text-center
```{button-ref} references
:ref-type: doc
:color: primary
:shadow:
References
```
:::
::::

::::{grid} 3
:gutter: 1

:::{grid-item}
```{image} ../examples/Laser_pulse.jpg
:alt: Laser pulse
:width: 100%
```
:::

:::{grid-item}
```{image} ../examples/Gyrotron.jpg
:alt: Gyrotron
:width: 100%
```
:::

:::{grid-item}
```{image} ../examples/electric_propulsion_new.jpg
:alt: Electric propulsion
:width: 100%
```
:::

:::{grid-item}
```{image} ../examples/Streamer.jpg
:alt: Streamer
:width: 100%
```
:::

:::{grid-item}
```{image} ../logo.png
:alt: PICLas
:width: 100%
```
:::

:::{grid-item}
```{image} ../examples/NozzleExpansion.jpg
:alt: Nozzle expansion
:width: 100%
```
:::

:::{grid-item}
```{image} ../examples/Titan_Aerocapture.jpg
:alt: Titan aerocapture
:width: 100%
```
:::

:::{grid-item}
```{image} ../examples/Magnetron.jpg
:alt: Magnetron
:width: 100%
```
:::

:::{grid-item}
```{image} ../examples/VacuumPump.jpg
:alt: Vacuum pump
:width: 100%
```
:::

::::

```{toctree}
---
maxdepth: 1
hidden: true
---
code_of_conduct.md
```

```{toctree}
---
maxdepth: 1
hidden: true
---
acknowledge.md
```

```{toctree}
---
maxdepth: 1
hidden: true
caption: User Guide
---
Overview <userguide/index>
```

```{toctree}
---
maxdepth: 2
hidden: true
numbered:
---
userguide/installation.md
userguide/meshing.md
userguide/workflow.md
userguide/features-and-models/index.md
userguide/visu_output.md
userguide/tools.md
userguide/tutorials/index.md
userguide/cluster_guide.md
userguide/appendix.md
```

```{toctree}
---
maxdepth: 1
hidden: true
---
references.md
```

```{toctree}
---
maxdepth: 1
hidden: true
caption: Developer Guide
---
Overview <developerguide/index>
```

```{toctree}
---
maxdepth: 2
hidden: true
numbered:
---
developerguide/git_workflow.md
developerguide/styleguide.md
developerguide/linting.md
developerguide/code_extension.md
developerguide/documentation.md
developerguide/release.md
developerguide/mpi.md
developerguide/reggie.md
developerguide/unittest.md
developerguide/troubleshooting.md
developerguide/performance.md
```
