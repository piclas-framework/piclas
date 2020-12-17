# 3D dielectric surface charge deposition via inner BCs (=100) no PartStateBoundary
NOTE: Because of TriaTracking=T, the decomposition of each inner face into two triangles depends on
the number of chosen procs, hence, the choice between a concave and convex side is arbitrary and
leads to different results for MPI=3,7,12 as compared with MPI=1,2 when checking the "PartStateBoundary" option. 
This test is therefore not conducted here, but in another reggie example. This test is conducted for 
MPI=1,2,3,7,12 and tests only DG,DG_Source,ElemData,DG_SourceExt (the higher number of procs correctly 
checks the partitioning of the dielectic<->vacuum interfaces
in MPI).

Tests
    * Dielectric sphere in sphere example (linear mesh NGeo=1)
    * Particle removal in dielectric regions during initial emission
    * Inner BC (100) which are treated as normal inner faces in the HDG solver
    * Particles impacting on the dielectric surface are swapped to no particles (removed)
    * Their charge is deposited (via cell_volweight_mean) to the surface nodes (cell vertices)
    * Checks dielectric interface detection for large number of MPI ranks
