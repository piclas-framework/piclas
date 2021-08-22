# Features & Models

```{toctree}
---
maxdepth: 1
caption: Table of Contents
---
particle-tracking.md
BC-field-solver.md
BC-particle-solver.md
particle-initialization-and-emission.md
PIC.md
superB.md
DSMC.md
BGG.md
Fokker-Planck.md
Bhatnagar-Gross-Krook.md
features-particle-solver.md
```

The goal of PICLas is to enable to approximation of the complete Boltzmann equation:

$$ \frac{\partial f}{\partial t} + \mathbf{v}\cdot\frac{\partial f}{\partial \mathbf{x}} + \frac{\mathbf{F}}{m}\cdot\frac{\partial f}{\partial \mathbf{v}} = \left.\frac{\partial f}{\partial t}\right|_{\mathrm{coll}} $$
