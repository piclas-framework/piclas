# Inner reflective boundaries (Cubit mesh)
* Mesh with 8 elements (4 inner BC sides) generated with Cubit, exported with the Gambit ASCII format
* Using inner boundaries of the Type 100 (in hopr.ini) for the field solver to model a dielectric (meshed) part
* Reflective boundaries for DSMC with surface sampling output
  * hopr.ini

            BoundaryName=BC_inner
            BoundaryType=(/100,0,0,0/)

  * parameter.ini

            Part-Boundary2-SourceName = BC_inner
            Part-Boundary2-Condition  = reflective