(sec:tutorial-dsmc-cone-3D-gmsh)=
# Hypersonic Flow around the 70° Cone (DSMC) - 3D Mesh with Gmsh

With the validation case of a 70° blunted cone already used in the previous tutorial ({ref}`sec:tutorial-dsmc-cone-2D`), the 3D mesh generation using [Gmsh](https://gmsh.info/) is presented in greater detail in this tutorial.
Before starting, copy the `dsmc-cone-gmsh` directory from the tutorial folder in the top level directory to a separate location

    cp -r $PICLAS_PATH/tutorials/dsmc-cone-3D .
    cd dsmc-cone-3D

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

    Physical Surface("IN", 29) = {4, 1};
    Physical Surface("SYM", 30) = {3, 5};
    Physical Surface("OUT", 31) = {2};
    Physical Surface("WALL", 32) = {7, 8, 9, 10, 11, 6};

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

The commands `Mesh.MeshSizeMin` and `Mesh.MeshSizeMax` define the minimum and maximum mesh element sizes. With the prescribed `Field` options, the size of the mesh can be specified using an explicit mathematical function using `MathEval` and restriced to specific surfaces with `Restrict`. In this tutorial, a mesh refinement at the frontal wall of the cone is enabled with this. `Background Field = 2` sets `Field[2]` as background field.
Different meshing algorithms for creating the 2D and 3D meshes can be chosen within Gmsh. The command `Mesh.SubdisionAlgorithm = 2` enables the generation of a fully hexahedral mesh by subdivision of cells.
`Mesh.OptimizeNetgen` improves the mesh quality additionally.

Next, the 3D mesh is created:

    Mesh 3;

The following commands are required to save all elements even if they are not part of a physical group and to use the ASCII format, before saving the mesh as `70degCone_3D.msh`:

    Mesh.SaveAll = 1;
    Mesh.Binary = 0;
    Mesh.MshFileVersion = 4.1;
    Save "70degCone_3D.msh";

The mesh file in the file format `.h5` used by **piclas** has to be converted using HOPR by supplying an input file `hopr.ini` using the corresponding mode:

    Mode = 5

As another possibility, the `SplitToHex` option can be enabled in the `hopr.ini` file instead of using the `SubdivionAlgorithm` command in Gmsh. The expected result for the 3D mesh is shown in {numref}`fig:dsmc-cone-gmsh-mesh`.

```{figure} mesh/dsmc-cone-gmsh-mesh.jpg
---
name: fig:dsmc-cone-gmsh-mesh
width: 600px
---

3D mesh of the 70° cone.
```

## Flow simulation with DSMC

For the general information on the setup of the DSMC simulation, please see the previous tutorial {ref}`sec:tutorial-dsmc-cone-2D`. In this tutorial, only the changes in the `parameter.ini` file for the 3D simulation compared to the 2D simulation are explained further.

First, the mesh file name is adapted. The number of boundaries is reduced from five to four, as only one symmetrical boundary is used in this 3D simulation. The maximum particle number per processor is increased due to the changed simulation domain. Additionally, the octree is adapted to appropriate values for a 3D simulation.

    MeshFile                        = 70degCone_3D_mesh.h5
    Part-nBounds                    = 4
    Part-maxParticleNumber          = 1500000
    Particles-OctreePartNumNode     = 80
    Particles-OctreePartNumNodeMin  = 60

Compared to the `parameter.ini` for the 2D simulation, the symmetrical boundaries, the commands for the 2D axisymmetric simulation and the radial weighting are deleted:

    Part-Boundary4-SourceName                   = SYMAXIS
    Part-Boundary4-Condition                    = symmetric_axis
    Part-Boundary5-SourceName                   = ROTSYM
    Part-Boundary5-Condition                    = symmetric
    Particles-Symmetry-Order                    = 2
    Particles-Symmetry2DAxisymmetric            = T
    Particles-RadialWeighting                   = T
    Particles-RadialWeighting-PartScaleFactor   = 60
    Particles-RadialWeighting-CloneMode         = 2
    Particles-RadialWeighting-CloneDelay        = 5

Instead, a new symmetrical boundary is added:

    Part-Boundary4-SourceName  = SYM
    Part-Boundary4-Condition   = symmetric

An exemplary simulation result using the 3D mesh generated with Gmsh is shown in {numref}`fig:dsmc-cone-gmsh-visu`.

```{figure} results/dsmc-cone-gmsh-visu.jpg
---
name: fig:dsmc-cone-gmsh-visu
width: 600px
---

Translational temperature around the 70° cone.
```