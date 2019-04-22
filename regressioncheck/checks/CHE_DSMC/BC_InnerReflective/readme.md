# Inner reflective boundaries
* Using inner boundaries of the Type 100 (in preproc.ini) for the field solver to model a dielectric (meshed) part
* Reflective boundaries for DSMC with surface sampling output
  * preproc.ini

            BoundaryName=BC_inner
            BoundaryType=(/100,0,0,0/)

  * parameter.ini

            Part-Boundary2-SourceName = BC_inner
            Part-Boundary2-Condition  = reflective
