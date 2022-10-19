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
    v() = ShapeFromFile("70degCone_3D_model.step");

The simulation domain is created next by adding a cylindrical section and subtracting the volume of the cone.

    Cylinder(2) = {-50, 0, 0, 100, 0, 0, 50, Pi/6};
    BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; Delete; };

Physical groups are used to define the boundary conditions at all surfaces:

    Physical Surface("BC_INFLOW", 29) = {4, 1};
    Physical Surface("BC_SYMMETRY", 30) = {3, 5};
    Physical Surface("BC_OUTFLOW", 31) = {2};
    Physical Surface("BC_WALL", 32) = {7, 8, 9, 10, 11, 6};

The mesh options can be set with the following commands:

    Mesh.MeshSizeMin = 1;
    Mesh.MeshSizeMax = 10;
    Field[1] = MathEval;
    Field[1].F = "0.2";
    Field[2] = Restrict;
    Field[2].SurfacesList = {7, 8, 9};
    Background Field = 2;
    Mesh.Algorithm = 1;
    Mesh.Algorithm3D = 7;
    Mesh.SubdivisionAlgorithm = 2;
    Mesh.OptimizeNetgen = 1;

With the added `Field` options, the size of the mesh can be specified using an explicit mathematical function using `MathEval` and restriced to specific surfaces with `Restrict`. In this tutorial, a mesh refinement at the frontal wall of the cone is enabled with this.
Different meshing algorithms can be chosen within Gmsh

mesh

save, export

## Flow simulation with DSMC

changes compared to 2D simulation
{ref}`sec:tutorial-dsmc-cone-2D`