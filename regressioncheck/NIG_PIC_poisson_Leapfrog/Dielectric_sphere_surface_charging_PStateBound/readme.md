# 3D dielectric surface charge deposition via inner BCs (=100)
NOTE: Because of TriaTracking=T, the decomposition of each inner face into two triangles depends on
the number of chosen procs, hence, the choice between a concave and convex side is arbitrary and
leads to different results for MPI=3,7,12 as compared with MPI=1,2. 
Therefore, MPI=3,7,12 are not used and a higher tolerance is chosen for the comparison of the output files for MPI=3.
Furthermore, the analysis and comparison with reference files for PartStateBoundary and SurfaceData is deactivated until a solution is found.

Tests
    * Dielectric sphere in sphere example (linear mesh NGeo=1)
    * Particle removal in dielectric regions during initial emission
    * Inner BC (100) which are treated as normal inner faces in the HDG solver
    * Particles impacting on the dielectric surface are swapped to no particles (removed)
    * Their charge is deposited (via cell_volweight_mean) to the surface nodes (cell vertices)
    * PartStateBoundary output and DSMC SurfState output
