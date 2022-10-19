(sec:tutorial-dsmc-cone-3D-gmsh)=
# Hypersonic Flow around the 70° Cone (DSMC) - 3D Mesh with Gmsh

With the validation case of a 70° blunted cone already used in the previous tutorial ({ref}`sec:tutorial-dsmc-cone-2D`), the 3D mesh generation using [Gmsh](https://gmsh.info/) is presented in greater detail in this tutorial.
Before starting, copy the `dsmc-cone-gmsh` directory from the tutorial folder in the top level directory to a separate location

    cp -r $PICLAS_PATH/tutorials/dsmc-cone-gmsh .
    cd dsmc-cone-gmsh

The general information needed to setup a DSMC simulation is given in the previous tutorials {ref}`sec:tutorial-dsmc-reservoir` and {ref}`sec:tutorial-dsmc-cone-2D`. The following focuses on the mesh generation with Gmsh and case-specific differences for the DSMC simulation.

## Mesh generation with Gmsh

First, create a new file in gmsh: `70DegCone_3D.geo`. In general, the mesh can be generated using the GUI or by using the `.geo` script environment. In the GUI, the script can be edited via `Edit script` and loaded with `Reload script`. This tutorial focuses on the scripting approach.

After opening the `.geo` script file, select the OpenCASCADE CAD kernel and open the provided `70DegCone_3D_model.step` file with the following commands:

    SetFactory("OpenCASCADE");
    (v) = 70DegCone_3D_model.step;

physical groups

mesh options

mesh + refinements

edit script

save, export

## Flow simulation with DSMC

changes compared to 2D simulation
{ref}`sec:tutorial-dsmc-cone-2D`