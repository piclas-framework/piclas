# Tutorials

This chapter will give a detailed overview of simulations with **PICLas**.
It assumes that you are familiar with setting the compiler options and compiling the code.
The paths to the executables are omitted. It is assumed that you have either added aliases for **piclas**,
**hopr** and **piclas2vtk**, or that you added the binary directories to your `$PATH` variable as
described in {ref}`sec:directory-paths`.

Each tutorial is equipped with .ini files, *hopr.ini*, *parameter.ini*, *posti.ini* and for DSMC setups *DSMC.ini*,
as well as the mesh file *\*\_mesh.h5* in the HDF5 format (created with **hopr**).

    hopr.ini
    parameter.ini
    DSMC.ini
    posti.ini
    mesh.h5

It is suggested to copy each folder to a new directory, where you can run and modify the parameter (*.ini) files.



```{toctree}
---
maxdepth: 2
caption: Table of Contents
---
pic-poisson-plasma-wave/pic-poisson-plasma-wave.md
dsmc-reservoir/dsmc-reservoir.md
dsmc-cone-2D/dsmc-cone-2D.md
dsmc-cone-3D/dsmc-cone-3D.md
```
