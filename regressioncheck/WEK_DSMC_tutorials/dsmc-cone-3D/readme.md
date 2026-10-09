# WEK_DSMC_tutorials: dsmc-cone-3D
- runs the tutorial in `piclas/tutorials/dsmc-cone-3D` by first executing `gmsh` to create the `.msh`, then creating the mesh using
`pyhope` and then executing `piclas`
  - First, run `gmsh 70degCone_3D_model.geo -parse_and_exit` to create the `70degCone_3D.msh` file, which is used in `pyhope.ini`
  - Second, build the `.h5` mesh file using `pyhope`
  - Third, run `piclas`