# 2D dielectric surface charge deposition via inner BCs (=100)
- simulations
  - particle removal in dielectric regions during initial emission
  - inner BC (100) which are treated as normal inner faces in the HDG solver
  - particles impacting on the dielectric surface are removed (via swap species procedures) and their charge is deposited (via cell_volweight_mean) to the surface nodes (cell vertices)
  - Test restart functionality of surface charging module by using two restart files
    1. 2Dplasma_test_State_000.00000000000000000.h5
    1. 2Dplasma_test_State_000.00000005000000000.h5

- post-externals
  - post-piclas-DoDeposition-F/
    - run piclas and restart from 1e-7 with DoDeposition=F and a) Part-Boundary4-Dielectric=T and b) Part-Boundary4-Dielectric=F to test if this
    does not crash the program
  - post-vtk-statefile-conversion/
    - run piclas2vtk conversion with VisuParticles = T and NVisu = 1 applied to 2Dplasma_test_State_000.00000010000000000.h5
